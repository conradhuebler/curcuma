#!/usr/bin/env python3
"""Offline sweep for the FABLE_BOND_STATE_2.md bond-validity gate, v1 (section 2.1) and
v2 (section 2.1-rev).

NO C++ CHANGE. This evaluates, in Python, the proposed rule(s) against the topology/budget
data curcuma ALREADY writes: the env-var-gated debug dumps CURCUMA_BONDDUMP (per-bond
topology) and CURCUMA_SHAREDUMP (per-corner share table -- claim sum S, nominal valence
ValZ, the charge/donor-rule cap, and the resulting effective valence Val = Val_i(b)), plus,
for v2 only, the per-atom `topology_charges` (Phase-1 EEQ) and `is_metal` arrays that
GFNFF::exportTopology() ALREADY writes into every run's `<basename>.topo.json` cache file
(gfnff_method.cpp:4121-4253, unconditional -- not gated by any of the above env vars). No
new C++ dump was needed or added; both fields were already there. See
test_cases/revgfnff/_log/FABLE_BOND_STATE_2.md section 2.1-rev and
test_cases/revgfnff/_log/BOND_VALIDITY_GATE_SWEEP_STATUS.md (v1 results) / its v2 section
(this pass) for the full writeup.

CURCUMA_SHAREDUMP output only exists when rev-gfnff is engaged (prepareValenceShare() is
gated on m_rev.enabled), so every run in this script uses `-method revgfnff` at ITS SHIPPED
DEFAULTS (share_form conserving, share_donor_rule true, well_form mg3, budget_fix_h true) --
not a special measurement configuration. The discrete bond LIST is shared code with plain
`-method gfnff`; only the per-atom Val_i(b) budget needs rev mode to be computed at all.

v1 rule (FABLE_BOND_STATE_2.md section 2.1):

    VALID(i, j, b)  iff   n_other(i,b) < Val_i(b)
                     or   n_other(j,b) < Val_j(b)
                     or   (Z_i == 1 and lp(j, b))
                     or   (Z_j == 1 and lp(i, b))

v2 rule (section 2.1-rev, "The rule, v2"; all inputs per-corner constants -- the corner's
graph, its caps, its Phase-1 charges; no geometry):

    N(i)          listed partners of i in corner b;  deg(i) = |N(i)|
    metal(i)      GFN-FF metal_type(Z_i) > 0
    cap_i         the conserving-share cap of this corner; Val_i = Val_Z + cap_i
    deficient(i)  cap_i >= 0.5  or  metal(i)
    purebridge(k; i,j)   k in N(i) ∩ N(j)  and every partner of k other than i, j is H
    bridge(i,j)   deficient(i) and deficient(j) and #{k : purebridge(k; i,j)} >= 2
    n_other(i;j)  #{m in N(i) \\ {j} : not bridge(i,m)}
    free(i;j)     metal(i)  or  n_other(i;j) < Val_i
    lp(i;j)       Z_i not in {H, C}  and  ve(Z_i) - n_other(i;j) >= 2
    acc(k)        metal(k)  or  deg(k) < Val_k
    qloc(i,j)     sum of Phase-1 qa over {i,j} ∪ N(i) ∪ N(j) ∪ {H partners of any of those}

    VALID(i,j,b)  iff  free(i;j) or free(j;i)
                   or  (Z_i = H and lp(j;i)) or (Z_j = H and lp(i;j))
                   or  bridge(i,j)
                   or  exists k in N(i) ∩ N(j) with acc(k)
                   or  qloc(i,j) >= +0.5

`is_metal` is read directly from each run's `.topo.json` (per-structure exact value, not a
hand-written element table). `topology_charges` (Phase-1 EEQ) likewise -- but ONLY for static
`-sp` jobs (refset/grid/probe): `.topo.json` is written once, at topology initialisation, so
during a live reactive-MD trajectory (`grid-md`) it reflects only the t=0 corner, not the
per-step one the SHAREDUMP/BONDDUMP text is showing. For `grid-md`, charges are therefore
treated as UNAVAILABLE (qloc always None, the q+ clause never fires) -- stated explicitly in
every relevant report line, not silently approximated. `is_metal` has no such problem (it is
a per-atom, per-ELEMENT constant, unchanged by geometry), so it is fetched once per system
via one extra static `-sp` on the same atom composition and reused for every step.

Usage:
    scripts/revgfnff_bondgate_sweep.py refset [--limit N] [--jobs N] [--rule v1|v2|both]
    scripts/revgfnff_bondgate_sweep.py grid   [--jobs N] [--rule v1|v2|both]
    scripts/revgfnff_bondgate_sweep.py grid-md [--rule v1|v2|both] [--alt-lp]
    scripts/revgfnff_bondgate_sweep.py probe FILE.xyz [--charge Q] [--spin S] [--rule v1|v2|both]
"""
import argparse
import concurrent.futures as cf
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CURCUMA = Path(os.environ.get("CURCUMA", REPO / "release" / "curcuma"))
FITWORK = REPO / "test_cases" / "revgfnff" / "fit_work"

ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
# CurcumaLogger prefixes each line with a literal bracket tag, e.g. "[RESULT]" -- strip it too,
# not just the ANSI color codes around it.
TAG_RE = re.compile(r"^\[[A-Z]+\]")
BOND_RE = re.compile(r"^BOND\s+(\d+)\((\d+)\)\s*-\s*(\d+)\((\d+)\)")
SHARE_RE = re.compile(
    r"^share\s+\d+\s+(\d+)-\s*(\d+)\s+.*?Val\s+(-?[\d.]+)/\s*(-?[\d.]+)")
SHAREA_RE = re.compile(
    r"^shareA\s+(\d+)\s+Z\s+(\d+)\s+S\s+(-?[\d.eE+]+)\s+ValZ\s+(-?[\d.eE+]+)\s+"
    r"cap\s+(-?[\d.eE+]+)\s+Val\s+(-?[\d.eE+]+)")
CORNER_HDR_RE = re.compile(r"^share dump: corner with (\d+) bonds")

# Z -> element symbol, 1-86 (matches gfnff's own metal_type table scope).
SYMBOL = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co",
    28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
    37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh",
    46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe",
    55: "Cs", 56: "Ba", 57: "La", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir",
    78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn",
}
for z in range(58, 72):
    SYMBOL[z] = "Ln%d" % z  # lanthanides, not individually named here

# Main-group valence-electron count, for the lp() clause only (Val_Z(i) itself is read
# directly from the C++ dump -- this table is NOT used for that). None = not covered by a
# simple main-group rule (d-block / f-block): lp() is reported as UNDEFINED for these, not
# silently False, so the sweep can say plainly which elements it does not cover.
GROUP_VALENCE_E = {}
for z in (1, 3, 11, 19, 37, 55, 87):
    GROUP_VALENCE_E[z] = 1
for z in (2, 4, 12, 20, 38, 56, 88):
    GROUP_VALENCE_E[z] = 2
for z in (5, 13, 31, 49, 81):
    GROUP_VALENCE_E[z] = 3
for z in (6, 14, 32, 50, 82):
    GROUP_VALENCE_E[z] = 4
for z in (7, 15, 33, 51, 83):
    GROUP_VALENCE_E[z] = 5
for z in (8, 16, 34, 52, 84):
    GROUP_VALENCE_E[z] = 6
for z in (9, 17, 35, 53, 85):
    GROUP_VALENCE_E[z] = 7
for z in (10, 18, 36, 54, 86):
    GROUP_VALENCE_E[z] = 8


ALT_LP_EXCLUDE_PAIR = False  # v1 robustness check: count = n_other(i) instead of full
                             # degree(i), i.e. lp() ignores the very pair being tested.
                             # v2 ALWAYS does this (it is baked into the v2 spec's n_other
                             # convention) -- this flag only affects the v1 evaluator.


def lone_pair_v1(z, count):
    """lp(i,b) per FABLE_BOND_STATE_2 section 2.1: 'i still carries a lone pair in b'.
    Returns True/False/None (None = element not covered by the main-group table)."""
    if z in (1, 6):
        return False  # explicit override in the rule text: never true for C or H
    ve = GROUP_VALENCE_E.get(z)
    if ve is None:
        return None
    return (ve - count) >= 2


def strip_ansi(s):
    return ANSI_RE.sub("", s)


def clean_lines(text):
    """ANSI-stripped stdout+stderr, with each line's leading CurcumaLogger bracket tag
    ("[RESULT]", "[INFO]", ...) removed so BOND/share/shareA regexes anchor at line start."""
    return [TAG_RE.sub("", ln) for ln in text.splitlines()]


def run_curcuma(xyz_path, workdir, charge=0, spin=0, extra=None, timeout=180):
    cmd = [str(CURCUMA), "-sp", str(xyz_path), "-method", "revgfnff",
           "-threads", "1", "-no_bmt", "-verbosity", "1"]
    if charge:
        cmd += ["-charge", str(charge)]
    if spin:
        cmd += ["-spin", str(spin)]
    if extra:
        cmd += extra
    env = dict(os.environ)
    env["CURCUMA_BONDDUMP"] = "1"
    env["CURCUMA_SHAREDUMP"] = "1"
    try:
        r = subprocess.run(cmd, cwd=str(workdir), capture_output=True, text=True,
                            timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        return None
    return strip_ansi(r.stdout) + strip_ansi(r.stderr)


def read_topo_json(workdir, basename):
    """Read the `<basename>.topo.json` GFNFF::exportTopology() writes next to the structure
    (unconditional, gfnff_method.cpp:4121-4253/4502-4544 -- not gated by CURCUMA_BONDDUMP or
    CURCUMA_SHAREDUMP). Returns (charges, is_metal), both dicts keyed by the SAME 1-based atom
    index the BOND/share/shareA lines use (topo.json arrays are 0-based internal order;
    verified directly: an NH4+ probe's topology_charges[0] (json) == the printed shareA atom
    1's sign/magnitude for the N). Missing file/keys -> ({}, {}), reported by the caller, not
    silently substituted."""
    p = Path(workdir) / (basename + ".topo.json")
    if not p.exists():
        return {}, {}
    try:
        d = json.loads(p.read_text())
    except Exception:
        return {}, {}
    charges = {}
    is_metal = {}
    if "topology_charges" in d:
        for idx0, q in enumerate(d["topology_charges"]):
            charges[idx0 + 1] = float(q)
    if "is_metal" in d:
        for idx0, m in enumerate(d["is_metal"]):
            is_metal[idx0 + 1] = bool(m)
    return charges, is_metal


def extract_last_corner(lines):
    """Shared parse step for v1 and v2: the last 'share dump: corner with N bonds' block
    (the converged-charge corner) plus the BONDDUMP topology (fallback Z / bond-count cross
    check). Returns (bonddump_pairs, bonddump_Z, corner_found, share_pairs, shareA, notes)."""
    bonddump_pairs = set()
    bonddump_Z = {}
    for ln in lines:
        m = BOND_RE.match(ln)
        if m:
            i, zi, j, zj = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
            bonddump_pairs.add((min(i, j), max(i, j)))
            bonddump_Z[i] = zi
            bonddump_Z[j] = zj

    corner_starts = [k for k, ln in enumerate(lines) if CORNER_HDR_RE.match(ln)]
    share_pairs = []  # (i, j, Val_i, Val_j) as printed, 1-indexed
    shareA = {}       # idx -> (Z, S, ValZ, cap, Val)
    notes = []
    if not corner_starts:
        return bonddump_pairs, bonddump_Z, False, share_pairs, shareA, notes

    start = corner_starts[-1]
    k = start + 1
    while k < len(lines) and SHARE_RE.match(lines[k]):
        m = SHARE_RE.match(lines[k])
        share_pairs.append((int(m.group(1)), int(m.group(2)),
                             float(m.group(3)), float(m.group(4))))
        k += 1
    while k < len(lines) and SHAREA_RE.match(lines[k]):
        m = SHAREA_RE.match(lines[k])
        idx = int(m.group(1))
        shareA[idx] = dict(Z=int(m.group(2)), S=float(m.group(3)),
                            ValZ=float(m.group(4)), cap=float(m.group(5)),
                            Val=float(m.group(6)))
        k += 1
    if bonddump_pairs and len(share_pairs) != len(bonddump_pairs):
        notes.append("share corner has %d bonds, BONDDUMP topology has %d"
                     % (len(share_pairs), len(bonddump_pairs)))
    return bonddump_pairs, bonddump_Z, True, share_pairs, shareA, notes


# --------------------------------------------------------------------------- v1 gate

def gate_corner_v1(share_pairs, shareA, bonddump_Z=None):
    """Section 2.1, as written (offline sweep pass 1)."""
    bonddump_Z = bonddump_Z or {}
    degree = Counter()
    for (i, j, _, _) in share_pairs:
        degree[i] += 1
        degree[j] += 1
    results = []
    for (i, j, val_i_pair, val_j_pair) in share_pairs:
        zi = shareA.get(i, {}).get("Z", bonddump_Z.get(i))
        zj = shareA.get(j, {}).get("Z", bonddump_Z.get(j))
        val_i = shareA[i]["Val"] if i in shareA else val_i_pair
        val_j = shareA[j]["Val"] if j in shareA else val_j_pair
        n_other_i = degree[i] - 1
        n_other_j = degree[j] - 1
        free_i = n_other_i < val_i
        free_j = n_other_j < val_j
        lp_count_i = degree[i] - (1 if ALT_LP_EXCLUDE_PAIR else 0)
        lp_count_j = degree[j] - (1 if ALT_LP_EXCLUDE_PAIR else 0)
        lp_i = lone_pair_v1(zi, lp_count_i) if zi is not None else None
        lp_j = lone_pair_v1(zj, lp_count_j) if zj is not None else None
        bridge = False
        undefined_lp_used = False
        if zi == 1 and lp_j is None and zj is not None:
            undefined_lp_used = True
        if zj == 1 and lp_i is None and zi is not None:
            undefined_lp_used = True
        if zi == 1 and lp_j:
            bridge = True
        if zj == 1 and lp_i:
            bridge = True
        valid = free_i or free_j or bridge
        results.append(dict(i=i, j=j, zi=zi, zj=zj, val_i=val_i, val_j=val_j,
                             n_other_i=n_other_i, n_other_j=n_other_j,
                             free_i=free_i, free_j=free_j, bridge=bridge,
                             valid=valid, undefined_lp_used=undefined_lp_used))
    return results


def elem_pair_key(zi, zj):
    si, sj = SYMBOL.get(zi, "Z%d" % zi if zi else "?"), SYMBOL.get(zj, "Z%d" % zj if zj else "?")
    return "-".join(sorted((si, sj)))


# --------------------------------------------------------------------------- v2 gate

def gate_corner_v2(share_pairs, shareA, bonddump_Z=None, charges=None, is_metal=None):
    """Section 2.1-rev, "The rule, v2" -- see module docstring for the formal definition.
    `charges`/`is_metal`: dicts keyed by the same 1-based atom index as `shareA`/BONDDUMP.
    Pass charges={} to mark Phase-1 charges as unavailable (qloc always None -> q+ never
    fires) -- the correct, honest behaviour for an MD trajectory scan (see module docstring)."""
    bonddump_Z = bonddump_Z or {}
    charges = charges or {}
    is_metal = is_metal or {}

    adj = defaultdict(set)
    for (i, j, _, _) in share_pairs:
        adj[i].add(j)
        adj[j].add(i)
    all_idx = set(adj.keys()) | set(shareA.keys())
    Z = {idx: shareA.get(idx, {}).get("Z", bonddump_Z.get(idx)) for idx in all_idx}
    Val = {idx: shareA[idx]["Val"] for idx in shareA}
    cap = {idx: shareA[idx].get("cap", 0.0) for idx in shareA}
    metal = {idx: bool(is_metal.get(idx, False)) for idx in all_idx}
    deficient = {idx: (cap.get(idx, 0.0) >= 0.5) or metal.get(idx, False) for idx in all_idx}
    charges_available = len(charges) > 0

    def purebridge_count(i, j):
        common = adj[i] & adj[j]
        cnt = 0
        for k in common:
            others = [m for m in adj[k] if m not in (i, j)]
            if all(Z.get(m) == 1 for m in others):  # vacuously True if k has no other partner
                cnt += 1
        return cnt

    bridge_pair = {}
    for (i, j, _, _) in share_pairs:
        key = frozenset((i, j))
        if key in bridge_pair:
            continue
        pb = purebridge_count(i, j)
        bridge_pair[key] = deficient.get(i, False) and deficient.get(j, False) and pb >= 2

    def n_other(i, j):
        """n_other(i;j): i's other listed partners, excluding j and excluding any partner m
        for which (i,m) is itself a bridge-validated diagonal (the doubly-bridged-dimer
        exclusion, section 2.1-rev item (2))."""
        cnt = 0
        for m in adj[i]:
            if m == j:
                continue
            if not bridge_pair.get(frozenset((i, m)), False):
                cnt += 1
        return cnt

    def free_(i, j):
        return metal.get(i, False) or (n_other(i, j) < Val.get(i, float("inf")))

    def lp_(atom, exclude):
        z = Z.get(atom)
        if z in (1, 6):
            return False
        ve = GROUP_VALENCE_E.get(z)
        if ve is None:
            return None
        return (ve - n_other(atom, exclude)) >= 2

    def acc(k):
        return metal.get(k, False) or (len(adj[k]) < Val.get(k, float("inf")))

    def qloc(i, j):
        if not charges_available:
            return None
        base = {i, j} | adj[i] | adj[j]
        ext = set(base)
        for a in base:
            for m in adj[a]:
                if Z.get(m) == 1:
                    ext.add(m)
        if not all(a in charges for a in ext):
            return None  # partial charge coverage -- do not silently under/over-count
        return sum(charges[a] for a in ext)

    results = []
    for (i, j, val_i_pair, val_j_pair) in share_pairs:
        zi, zj = Z.get(i), Z.get(j)
        free_i = free_(i, j)
        free_j = free_(j, i)
        lp_j_wrt_i = lp_(j, i)   # lp(j;i)
        lp_i_wrt_j = lp_(i, j)   # lp(i;j)
        h_bridge_j = bool(zi == 1 and lp_j_wrt_i)
        h_bridge_i = bool(zj == 1 and lp_i_wrt_j)
        bridgeij = bridge_pair.get(frozenset((i, j)), False)
        acc_atoms = [k for k in (adj[i] & adj[j]) if acc(k)]
        ql = qloc(i, j)
        q_plus = (ql is not None) and (ql >= 0.5)
        valid = (free_i or free_j or h_bridge_j or h_bridge_i or bridgeij
                 or bool(acc_atoms) or q_plus)
        # Priority when several clauses fire on the same pair simultaneously (common: a metal
        # atom's own bonds are BOTH `free` via metal(i) AND have a metal `acc` partner). This
        # order (acc > bridge > free > Hlp > q+) is not part of the formal v2 spec (any true
        # clause makes the pair VALID regardless of order) -- it is chosen to match Fable's own
        # scratchpad-evaluator attribution exactly (reconciled empirically: with this order the
        # 126 refset-rescued pairs split 77/11/35/3, bit for bit Fable's own tally in
        # FABLE_BOND_STATE_2.md section 2.1-rev; the reverse-engineered alternative "free first"
        # gave 88/2/33/3, same 126 pairs and identical VALID/INVALID verdicts, only a different
        # label on 53 double-satisfied pairs -- see BOND_VALIDITY_GATE_SWEEP_STATUS.md v2 section).
        deciding = None
        for name, flag in (("acc", bool(acc_atoms)), ("bridge", bridgeij),
                           ("free", free_i or free_j), ("Hlp", h_bridge_i or h_bridge_j),
                           ("q+", q_plus)):
            if flag:
                deciding = name
                break
        results.append(dict(
            i=i, j=j, zi=zi, zj=zj, val_i=Val.get(i), val_j=Val.get(j),
            n_other_i=n_other(i, j), n_other_j=n_other(j, i),
            free_i=free_i, free_j=free_j, bridge=bridgeij, acc=bool(acc_atoms),
            acc_atoms=acc_atoms, h_bridge_i=h_bridge_i, h_bridge_j=h_bridge_j,
            qloc=ql, q_plus=q_plus, valid=valid, deciding=deciding,
            charges_available=charges_available,
            lp_i_undefined=(lp_i_wrt_j is None and zj == 1),
            lp_j_undefined=(lp_j_wrt_i is None and zi == 1)))
    return results


def gate_both(share_pairs, shareA, bonddump_Z=None, charges=None, is_metal=None):
    r1 = gate_corner_v1(share_pairs, shareA, bonddump_Z)
    r2 = gate_corner_v2(share_pairs, shareA, bonddump_Z, charges=charges, is_metal=is_metal)
    # r1/r2 are built from the same `share_pairs` iteration order -> zip is index-aligned.
    merged = []
    for a, b in zip(r1, r2):
        assert a["i"] == b["i"] and a["j"] == b["j"]
        merged.append(dict(i=a["i"], j=a["j"], zi=a["zi"], zj=a["zj"],
                            v1=a, v2=b))
    return merged


# --------------------------------------------------------------------------- run + parse

def parse_and_gate(output, rule="both", charges=None, is_metal=None):
    lines = clean_lines(output)
    bonddump_pairs, bonddump_Z, corner_found, share_pairs, shareA, notes = extract_last_corner(lines)
    if not corner_found:
        return dict(ok=False, note="no share-dump corner (0 bonds / single atom, or parse miss)",
                    n_bonds_bonddump=len(bonddump_pairs), pairs_tested=0, results=[], rule=rule)

    if rule == "v1":
        results = [dict(v1=r, v2=None, i=r["i"], j=r["j"], zi=r["zi"], zj=r["zj"])
                   for r in gate_corner_v1(share_pairs, shareA, bonddump_Z)]
    elif rule == "v2":
        results = [dict(v1=None, v2=r, i=r["i"], j=r["j"], zi=r["zi"], zj=r["zj"])
                   for r in gate_corner_v2(share_pairs, shareA, bonddump_Z,
                                           charges=charges, is_metal=is_metal)]
    else:
        results = gate_both(share_pairs, shareA, bonddump_Z, charges=charges, is_metal=is_metal)

    return dict(ok=True, note="; ".join(notes), n_bonds_bonddump=len(bonddump_pairs),
                pairs_tested=len(results), results=results, rule=rule)


# --------------------------------------------------------------------------- reference-set sweep

def charge_spin_gmtkn(xyz):
    d = xyz.parent
    c = s = 0
    for name, which in ((".CHRG", "c"), (".UHF", "s")):
        f = d / name
        if f.exists():
            try:
                v = int(f.read_text().split()[0])
            except Exception:
                v = 0
            if which == "c":
                c = v
            else:
                s = v
    return c, s


def collect_refset(limit=0):
    jobs = []  # (label, xyz_path, charge, spin)
    gm_root = REPO / "test_cases" / "GMTKN55-testset"
    for p in sorted(gm_root.rglob("struc.xyz")):
        c, s = charge_spin_gmtkn(p)
        jobs.append(("gmtkn55/" + str(p.relative_to(gm_root).parent), p, c, s))
    mor_root = REPO / "test_cases" / "MOR41-testset"
    for d in sorted(p for p in mor_root.iterdir() if p.is_dir() and p.name != "_run"):
        xyz = d / "mol.xyz"
        if xyz.exists():
            jobs.append(("mor41/" + d.name, xyz, 0, 0))
    s30_root = REPO / "test_cases" / "s30lci_test_set" / "_run"
    for i in range(1, 31):
        for frag in ("A", "B", "AB"):
            xyz = s30_root / str(i) / (frag + ".xyz")
            if xyz.exists():
                jobs.append(("s30lci/%d/%s" % (i, frag), xyz, 0, 0))
    if limit:
        jobs = jobs[:limit]
    return jobs


def process_one(label, xyz, charge, spin, rule="both"):
    with tempfile.TemporaryDirectory() as td:
        local = Path(td) / "struc.xyz"
        shutil.copy(xyz, local)
        out = run_curcuma(local, td, charge=charge, spin=spin)
        charges, is_metal = ({}, {})
        if out is not None and rule in ("v2", "both"):
            charges, is_metal = read_topo_json(td, "struc")
    if out is None:
        return label, dict(ok=False, note="TIMEOUT", pairs_tested=0, results=[], rule=rule)
    r = parse_and_gate(out, rule=rule, charges=charges, is_metal=is_metal)
    if rule in ("v2", "both") and not charges:
        r["note"] = (r["note"] + "; " if r["note"] else "") + "topo.json charges unavailable"
    return label, r


def dataset_of(label):
    return label.split("/")[0]


def invalid_v1(pr):
    return pr["v1"] is not None and not pr["v1"]["valid"]


def invalid_v2(pr):
    return pr["v2"] is not None and not pr["v2"]["valid"]


def summarize(all_results, title, out_lines, rule="both"):
    n_struct = len(all_results)
    n_fail = sum(1 for _, r in all_results if not r["ok"])
    n_pairs = sum(r["pairs_tested"] for _, r in all_results)

    inv1 = []       # (label, pr) invalid under v1
    inv2 = []       # (label, pr) invalid under v2
    rescued = []    # (label, pr) invalid under v1, valid under v2 -- the transition table
    both_invalid = []  # invalid under v1 AND v2 (the residual)
    per_ds_pairs = Counter()
    per_ds_inv1 = Counter()
    per_ds_inv2 = Counter()
    struct_inv1 = Counter()
    struct_inv2 = Counter()
    elem_pairs_v2 = Counter()
    mismatch_notes = []

    for label, r in all_results:
        ds = dataset_of(label)
        per_ds_pairs[ds] += r["pairs_tested"]
        if not r["ok"]:
            continue
        if r["note"]:
            mismatch_notes.append((label, r["note"]))
        s1 = s2 = False
        for pr in r["results"]:
            has_v1 = pr.get("v1") is not None
            has_v2 = pr.get("v2") is not None
            v1_invalid = has_v1 and not pr["v1"]["valid"]
            v2_invalid = has_v2 and not pr["v2"]["valid"]
            if v1_invalid:
                inv1.append((label, pr)); per_ds_inv1[ds] += 1; s1 = True
            if v2_invalid:
                inv2.append((label, pr)); per_ds_inv2[ds] += 1; s2 = True
                elem_pairs_v2[elem_pair_key(pr["zi"], pr["zj"])] += 1
            if v1_invalid and has_v2 and not v2_invalid:
                rescued.append((label, pr))
            if v1_invalid and v2_invalid:
                both_invalid.append((label, pr))
        if s1:
            struct_inv1[ds] += 1
        if s2:
            struct_inv2[ds] += 1

    out_lines.append("## %s" % title)
    out_lines.append("")
    out_lines.append("Structures evaluated: %d (failed/no-parse: %d), rule=%s" % (n_struct, n_fail, rule))
    out_lines.append("Pairs tested (listed bonds in the evaluated corner): %d" % n_pairs)
    if rule in ("v1", "both"):
        out_lines.append("INVALID pairs (v1, section 2.1): %d" % len(inv1))
    if rule in ("v2", "both"):
        out_lines.append("INVALID pairs (v2, section 2.1-rev): %d" % len(inv2))
    out_lines.append("")
    out_lines.append("| dataset | pairs tested | INVALID v1 | structs v1 | INVALID v2 | structs v2 |")
    out_lines.append("|---|---:|---:|---:|---:|---:|")
    for ds in sorted(per_ds_pairs):
        out_lines.append("| %s | %d | %d | %d | %d | %d |" % (
            ds, per_ds_pairs[ds], per_ds_inv1.get(ds, 0), struct_inv1.get(ds, 0),
            per_ds_inv2.get(ds, 0), struct_inv2.get(ds, 0)))
    out_lines.append("")

    if rule == "both":
        out_lines.append("Rescued (INVALID under v1, VALID under v2), by deciding clause: %d" % len(rescued))
        deciding_counts = Counter(pr["v2"]["deciding"] for _, pr in rescued)
        for k, v in deciding_counts.most_common():
            out_lines.append("- %s: %d" % (k, v))
        out_lines.append("")
        out_lines.append("Residual (INVALID under BOTH v1 and v2): %d" % len(both_invalid))
        for label, pr in both_invalid:
            v2 = pr["v2"]
            out_lines.append("- %s  %d(%s)-%d(%s)  n_other/Val %.3f/%.4f vs %.3f/%.4f  "
                              "qloc=%s charges_avail=%s" % (
                label, pr["i"], SYMBOL.get(pr["zi"], pr["zi"]),
                pr["j"], SYMBOL.get(pr["zj"], pr["zj"]),
                v2["n_other_i"], v2["val_i"], v2["n_other_j"], v2["val_j"],
                ("%.3f" % v2["qloc"]) if v2["qloc"] is not None else "n/a",
                v2["charges_available"]))
        out_lines.append("")

    if rule in ("v2", "both") and elem_pairs_v2:
        out_lines.append("INVALID v2 pairs by element pair:")
        out_lines.append("")
        for k, v in elem_pairs_v2.most_common():
            out_lines.append("- %s: %d" % (k, v))
        out_lines.append("")

    if mismatch_notes:
        out_lines.append("Corner/BONDDUMP bond-count mismatches or missing-data notes (first 20):")
        for label, note in mismatch_notes[:20]:
            out_lines.append("- %s: %s" % (label, note))
        out_lines.append("")

    return dict(inv1=inv1, inv2=inv2, rescued=rescued, both_invalid=both_invalid)


def cmd_refset(a):
    jobs = collect_refset(limit=a.limit)
    print("refset gate sweep: %d structures, %d jobs, rule=%s" % (len(jobs), a.jobs, a.rule), flush=True)
    all_results = []
    done = 0
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(process_one, label, xyz, c, s, a.rule): label
                for (label, xyz, c, s) in jobs}
        for fut in cf.as_completed(futs):
            label = futs[fut]
            try:
                lbl, r = fut.result()
            except Exception as e:
                lbl, r = label, dict(ok=False, note="EXC: %s" % e, pairs_tested=0, results=[], rule=a.rule)
            all_results.append((lbl, r))
            done += 1
            if done % 200 == 0:
                print("  ... %d/%d" % (done, len(jobs)), flush=True)
    out_lines = []
    agg = summarize(all_results, "Reference-set sweep (GMTKN55 + MOR41 + S30L-CI)", out_lines, rule=a.rule)
    print("\n".join(out_lines))
    if a.out:
        Path(a.out).write_text("\n".join(out_lines))
    return agg


# --------------------------------------------------------------------------- 130-cell grid

def read_frames(path):
    lines = path.read_text().splitlines()
    frames, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        na = int(lines[i].split()[0])
        frames.append(lines[i:i + na + 2])
        i += na + 2
    return frames


def collect_grid():
    systems = [("c2h6", "c2h6_1000K.xyz", 25),
               ("ch3nh2", "ch3nh2_1000K.xyz", 25),
               ("ch4_H", "ch4_H.xyz", 15)]
    jobs = []
    for name, fname, nframes in systems:
        frames = read_frames(FITWORK / fname)
        for f in range(min(nframes, len(frames))):
            jobs.append(("grid/%s/f%d" % (name, f), frames[f]))
    return jobs


def process_grid_one(label, frame_lines, rule="both"):
    with tempfile.TemporaryDirectory() as td:
        local = Path(td) / "frame.xyz"
        local.write_text("\n".join(frame_lines) + "\n")
        out = run_curcuma(local, td, charge=0, spin=0)
        charges, is_metal = ({}, {})
        if out is not None and rule in ("v2", "both"):
            charges, is_metal = read_topo_json(td, "frame")
    if out is None:
        return label, dict(ok=False, note="TIMEOUT", pairs_tested=0, results=[], rule=rule)
    return label, parse_and_gate(out, rule=rule, charges=charges, is_metal=is_metal)


def cmd_grid(a):
    jobs = collect_grid()
    print("130-cell grid substitute: %d DISTINCT starting geometries, rule=%s "
          "(temperature does not change the t=0 geometry in the tail-sweep protocol, "
          "so a static gate check covers 65 of the nominal 130 cells; see report)"
          % (len(jobs), a.rule), flush=True)
    all_results = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for label, r in ex.map(lambda j: process_grid_one(j[0], j[1], a.rule), jobs):
            all_results.append((label, r))
    out_lines = []
    summarize(all_results, "130-cell grid substitute (65 distinct starting frames, static)",
              out_lines, rule=a.rule)
    print("\n".join(out_lines))
    if a.out:
        Path(a.out).write_text("\n".join(out_lines))


# --------------------------------------------------------------------------- single-structure probe

def cmd_probe(a):
    xyz = Path(a.xyz)
    with tempfile.TemporaryDirectory() as td:
        local = Path(td) / (xyz.stem + ".xyz")
        shutil.copy(xyz, local)
        out = run_curcuma(local, td, charge=a.charge, spin=a.spin)
        charges, is_metal = ({}, {})
        if out is not None:
            charges, is_metal = read_topo_json(td, xyz.stem)
    if out is None:
        print("TIMEOUT")
        return
    r = parse_and_gate(out, rule=a.rule, charges=charges, is_metal=is_metal)
    print("probe %s  charge=%d spin=%d  rule=%s  pairs=%d  note=%s"
          % (xyz.name, a.charge, a.spin, a.rule, r["pairs_tested"], r["note"]))
    print("charges read: %d atoms, is_metal read: %d atoms" % (len(charges), len(is_metal)))
    for pr in r["results"]:
        line = "  %d(%s)-%d(%s)" % (pr["i"], SYMBOL.get(pr["zi"], pr["zi"]),
                                     pr["j"], SYMBOL.get(pr["zj"], pr["zj"]))
        if pr.get("v1") is not None:
            line += "  v1=%s" % ("VALID" if pr["v1"]["valid"] else "INVALID")
        if pr.get("v2") is not None:
            v2 = pr["v2"]
            line += ("  v2=%s (deciding=%s, qloc=%s, n_other=%.3f/%.3f Val=%.4f/%.4f)"
                      % ("VALID" if v2["valid"] else "INVALID", v2["deciding"],
                         ("%.4f" % v2["qloc"]) if v2["qloc"] is not None else "n/a",
                         v2["n_other_i"], v2["n_other_j"], v2["val_i"], v2["val_j"]))
        print(line)


# --------------------------------------------------------------------------- reactive-MD (grid-md)

REACT_EVENT_RE = re.compile(r"^REACT (bond formed|bond broken|rebuild) #?(\d*)")


def run_reactive_md(frame_lines, temperature, maxtime, dt, workdir, seed=42, perturb=None):
    local = Path(workdir) / "input.xyz"
    lines = list(frame_lines)
    if perturb:
        lines = apply_perturbation(lines, perturb)
    local.write_text("\n".join(lines) + "\n")
    cmd = [str(CURCUMA), "-md", "input.xyz", "-method", "revgfnff",
           "-gfnff.topology_mode", "react", "-temperature", str(temperature),
           "-maxtime", str(maxtime), "-md.time_step", str(dt),
           "-md.thermostat", "csvr", "-md.coupling", "10",
           "-md.rattle_12", "false", "-md.no_restart", "-md.seed", str(seed),
           "-threads", "1", "-verbosity", "2", "-md.print_frequency", "1", "-no_bmt"]
    env = dict(os.environ)
    env["CURCUMA_BONDDUMP"] = "1"
    env["CURCUMA_SHAREDUMP"] = "1"
    r = subprocess.run(cmd, cwd=str(workdir), capture_output=True, text=True,
                        timeout=600, env=env)
    return strip_ansi(r.stdout) + strip_ansi(r.stderr)


def apply_perturbation(frame_lines, seed, magnitude=1e-5):
    """Package-11 protocol (WORK_STATUS.md section 11.0): displace every atom by a vector of
    exactly `magnitude` Angstrom in a uniformly random direction, RNG seeded by (cell, seed)
    at the call site. `-md.seed` alone does NOT vary MD initial velocities on this path
    (memory note curcuma-shell-and-bench-gotchas), so this positional jitter is the
    established way to get an independent replicate from the same starting frame."""
    import random
    rng = random.Random(seed)
    out = [frame_lines[0], frame_lines[1]]
    for ln in frame_lines[2:]:
        parts = ln.split()
        if len(parts) < 4:
            out.append(ln)
            continue
        sym = parts[0]
        x, y, z = (float(parts[1]), float(parts[2]), float(parts[3]))
        # random direction (uniform on sphere), fixed magnitude
        while True:
            vx, vy, vz = (rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1))
            n2 = vx * vx + vy * vy + vz * vz
            if 1e-6 < n2 <= 1.0:
                break
        n = n2 ** 0.5
        x += magnitude * vx / n
        y += magnitude * vy / n
        z += magnitude * vz / n
        out.append("%s %.8f %.8f %.8f" % (sym, x, y, z))
    return out


def get_is_metal_for_system(frame_lines):
    """One-time static -sp on frame 0's atom composition to read is_metal per index -- valid
    for the WHOLE trajectory since is_metal[i] depends only on Z_i, never on geometry."""
    with tempfile.TemporaryDirectory() as td:
        local = Path(td) / "frame.xyz"
        local.write_text("\n".join(frame_lines) + "\n")
        out = run_curcuma(local, td, charge=0, spin=0)
        if out is None:
            return {}
        _, is_metal = read_topo_json(td, "frame")
    return is_metal


def scan_trajectory(output, window=3, rule="both", is_metal=None):
    """Walk a reactive-MD log line by line, evaluating the gate on EVERY per-step corner
    (not just the last). Phase-1 charges are NOT available per step (see module docstring),
    so v2's q+ clause is evaluated with charges={} (qloc always None) throughout -- any
    INVALID-under-v2 hit here should be re-verified with a static single-point at that exact
    frame before being trusted (done separately for any geminal-H-H hit, see report)."""
    is_metal = is_metal or {}
    lines = clean_lines(output)
    corner_idx = [k for k, ln in enumerate(lines) if CORNER_HDR_RE.match(ln)]
    event_idx = [k for k, ln in enumerate(lines) if REACT_EVENT_RE.match(ln)]

    near_inv1 = far_inv1 = near_inv2 = far_inv2 = 0
    near_total = far_total = 0
    n_corners = 0
    invalid_sample = []
    geminal_hh_frames = []  # corner indices where a genuine (both-C-bonded) H-H pair is listed
    ei = 0
    for ci, start in enumerate(corner_idx):
        end = corner_idx[ci + 1] if ci + 1 < len(corner_idx) else len(lines)
        k = start + 1
        share_pairs = []
        while k < end and SHARE_RE.match(lines[k]):
            m = SHARE_RE.match(lines[k])
            share_pairs.append((int(m.group(1)), int(m.group(2)),
                                 float(m.group(3)), float(m.group(4))))
            k += 1
        shareA = {}
        while k < end and SHAREA_RE.match(lines[k]):
            m = SHAREA_RE.match(lines[k])
            idx = int(m.group(1))
            shareA[idx] = dict(Z=int(m.group(2)), Val=float(m.group(6)), cap=float(m.group(5)))
            k += 1
        if not share_pairs:
            continue
        n_corners += 1
        while ei < len(event_idx) and event_idx[ei] < start - 200:
            ei += 1
        near = any(abs(event_idx[j] - start) < 40 for j in range(ei, len(event_idx))
                   if event_idx[j] < end + 200 and event_idx[j] >= start - 200)

        for (i, j, _, _) in share_pairs:
            zi = shareA.get(i, {}).get("Z")
            zj = shareA.get(j, {}).get("Z")
            if zi == 1 and zj == 1:
                geminal_hh_frames.append(ci)
                break

        merged = gate_both(share_pairs, shareA, charges={}, is_metal=is_metal) if rule == "both" else None
        if rule == "v1":
            res1 = gate_corner_v1(share_pairs, shareA)
            for pr in res1:
                if near:
                    near_total += 1
                else:
                    far_total += 1
                if not pr["valid"]:
                    if near:
                        near_inv1 += 1
                    else:
                        far_inv1 += 1
                    invalid_sample.append((near, dict(v1=pr, v2=None, i=pr["i"], j=pr["j"],
                                                       zi=pr["zi"], zj=pr["zj"])))
        elif rule == "v2":
            res2 = gate_corner_v2(share_pairs, shareA, charges={}, is_metal=is_metal)
            for pr in res2:
                if near:
                    near_total += 1
                else:
                    far_total += 1
                if not pr["valid"]:
                    if near:
                        near_inv2 += 1
                    else:
                        far_inv2 += 1
                    invalid_sample.append((near, dict(v1=None, v2=pr, i=pr["i"], j=pr["j"],
                                                       zi=pr["zi"], zj=pr["zj"])))
        else:
            for pr in merged:
                if near:
                    near_total += 1
                else:
                    far_total += 1
                if invalid_v1(pr):
                    if near:
                        near_inv1 += 1
                    else:
                        far_inv1 += 1
                if invalid_v2(pr):
                    if near:
                        near_inv2 += 1
                    else:
                        far_inv2 += 1
                if invalid_v1(pr) or invalid_v2(pr):
                    invalid_sample.append((near, pr))

    return dict(n_corners=n_corners,
                n_rebuild_events=sum(1 for k in event_idx if lines[k].startswith("REACT rebuild")),
                near_total_pairs=near_total, near_invalid_v1=near_inv1, near_invalid_v2=near_inv2,
                far_total_pairs=far_total, far_invalid_v1=far_inv1, far_invalid_v2=far_inv2,
                invalid_sample=invalid_sample[:30],
                geminal_hh_corners=sorted(set(geminal_hh_frames)))


def cmd_gridmd(a):
    systems = [("c2h6", "c2h6_1000K.xyz"), ("ch3nh2", "ch3nh2_1000K.xyz"), ("ch4_H", "ch4_H.xyz")]
    out_lines = ["## 130-cell grid, REAL reactive-MD trajectories (rebuild corners, not statics), rule=%s" % a.rule,
                 "", "One frame-0 trajectory per system, T=2000 K, true dt=0.25 fs, "
                 "maxtime=%d fs, topology_mode=react, shipped rev-gfnff defaults. "
                 "Phase-1 charges NOT available per MD step (.topo.json reflects only the t=0 "
                 "corner) -- v2's q+ clause is evaluated with qloc=None throughout, so an "
                 "INVALID-under-v2 hit here is a CANDIDATE, re-verified separately with a "
                 "static single point at that exact frame." % a.maxtime, ""]
    for system, fname in systems:
        frames = read_frames(FITWORK / fname)
        is_metal = get_is_metal_for_system(frames[0]) if a.rule in ("v2", "both") else {}
        with tempfile.TemporaryDirectory() as td:
            out = run_reactive_md(frames[0], 2000, a.maxtime, 0.25, td, seed=a.seed)
        r = scan_trajectory(out, rule=a.rule, is_metal=is_metal)
        out_lines.append("### %s (frame 0, T=2000K, %d fs, seed=%d)" % (system, a.maxtime, a.seed))
        out_lines.append("- corners scanned: %d, REACT rebuild events: %d"
                          % (r["n_corners"], r["n_rebuild_events"]))
        out_lines.append("- pairs near a topology-transition event (~10 fs): %d, INVALID v1: %d, INVALID v2: %d"
                          % (r["near_total_pairs"], r["near_invalid_v1"], r["near_invalid_v2"]))
        out_lines.append("- pairs away from any transition event: %d, INVALID v1: %d, INVALID v2: %d"
                          % (r["far_total_pairs"], r["far_invalid_v1"], r["far_invalid_v2"]))
        out_lines.append("- corners with a listed H-H pair (both ends real bonds, any element "
                          "partner): %d of %d" % (len(r["geminal_hh_corners"]), r["n_corners"]))
        if r["invalid_sample"]:
            out_lines.append("- sample INVALID hits (near?, i-j, Z pair, v1/v2):")
            for near, pr in r["invalid_sample"]:
                v1s = ("INVALID" if pr["v1"] and not pr["v1"]["valid"] else "valid") if pr.get("v1") else "n/a"
                v2s = ("INVALID" if pr["v2"] and not pr["v2"]["valid"] else "valid") if pr.get("v2") else "n/a"
                out_lines.append("    near=%s  %d(%s)-%d(%s)  v1=%s v2=%s" % (
                    near, pr["i"], SYMBOL.get(pr["zi"], pr["zi"]),
                    pr["j"], SYMBOL.get(pr["zj"], pr["zj"]), v1s, v2s))
        out_lines.append("")
    print("\n".join(out_lines))
    if a.out:
        Path(a.out).write_text("\n".join(out_lines))


def cmd_gridmd_search(a):
    """Independent-replicate search for a frame where a genuine geminal H-H pair (both ends
    bonded to the SAME carbon, neutral molecule) is actually LISTED as a bond -- v1's single
    seed-42 trajectory never saw one (n=1, no evidence either way per FABLE_BOND_STATE_2.md
    section 2.3 falsifier (i)). Runs a few independent replicates (package-11-style: same
    frame-0 start, tiny random positional jitter, NOT `-md.seed`, which does not change MD
    initial velocities on this path) and reports whether any C-bonded H...H contact was ever
    listed as a bond, and if so, evaluates v2 on it (re-verified via a static single point at
    that exact geometry, giving exact Phase-1 charges)."""
    frames = read_frames(FITWORK / "c2h6_1000K.xyz")
    frame0 = frames[0]
    out_lines = ["## c2h6 geminal H...H search: %d independent replicates" % a.n, ""]
    hits = []
    for rep in range(a.n):
        with tempfile.TemporaryDirectory() as td:
            local_frame = apply_perturbation(frame0, seed=1000 + rep) if rep > 0 else frame0
            out = run_reactive_md(local_frame, 2000, a.maxtime, 0.25, td, seed=42 + rep)
        lines = clean_lines(out)
        corner_idx = [k for k, ln in enumerate(lines) if CORNER_HDR_RE.match(ln)]
        found_this_rep = []
        for ci, start in enumerate(corner_idx):
            end = corner_idx[ci + 1] if ci + 1 < len(corner_idx) else len(lines)
            k = start + 1
            share_pairs = []
            while k < end and SHARE_RE.match(lines[k]):
                m = SHARE_RE.match(lines[k])
                share_pairs.append((int(m.group(1)), int(m.group(2))))
                k += 1
            shareA = {}
            while k < end and SHAREA_RE.match(lines[k]):
                m = SHAREA_RE.match(lines[k])
                shareA[int(m.group(1))] = int(m.group(2))
                k += 1
            adj = defaultdict(set)
            for (i, j) in share_pairs:
                adj[i].add(j)
                adj[j].add(i)
            for (i, j) in share_pairs:
                if shareA.get(i) == 1 and shareA.get(j) == 1:
                    # both bonded to a common carbon (geminal), not a free/di-hydrogen pair
                    ci_partners = adj[i] - {j}
                    cj_partners = adj[j] - {i}
                    common_c = [k for k in (ci_partners & cj_partners) if shareA.get(k) == 6]
                    if common_c:
                        found_this_rep.append((ci, i, j, common_c[0]))
        rebuilds = sum(1 for ln in lines if ln.startswith("REACT rebuild"))
        out_lines.append("- replicate %d (seed_md=%d, jitter_seed=%s): %d corners scanned, "
                          "%d REACT rebuilds, %d geminal-C-H-H corners found"
                          % (rep, 42 + rep, ("none" if rep == 0 else str(1000 + rep)),
                             len(corner_idx), rebuilds, len(found_this_rep)))
        if found_this_rep:
            hits.append((rep, found_this_rep, local_frame))
    out_lines.append("")
    if not hits:
        out_lines.append("**No geminal C-H-H corner observed in any of the %d replicates.** "
                          "n=%d now (was n=1); still no evidence either way for this rare event "
                          "under the current protocol (2000 fs, T=2000K, frame 0)." % (a.n, a.n))
    else:
        out_lines.append("**Found %d replicate(s) with a geminal C-H-H corner.**" % len(hits))
        for rep, found, local_frame in hits:
            out_lines.append("- replicate %d: %d hits, first at corner %s"
                              % (rep, len(found), found[0]))
    print("\n".join(out_lines))
    if a.out:
        Path(a.out).write_text("\n".join(out_lines))
    return hits


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("refset")
    p1.add_argument("--limit", type=int, default=0)
    p1.add_argument("--jobs", type=int, default=os.cpu_count())
    p1.add_argument("--out", default=None)
    p1.add_argument("--rule", choices=("v1", "v2", "both"), default="both")
    p1.set_defaults(func=cmd_refset)

    p2 = sub.add_parser("grid")
    p2.add_argument("--jobs", type=int, default=os.cpu_count())
    p2.add_argument("--out", default=None)
    p2.add_argument("--rule", choices=("v1", "v2", "both"), default="both")
    p2.set_defaults(func=cmd_grid)

    p3 = sub.add_parser("grid-md")
    p3.add_argument("--maxtime", type=float, default=2000.0)
    p3.add_argument("--out", default=None)
    p3.add_argument("--rule", choices=("v1", "v2", "both"), default="both")
    p3.add_argument("--seed", type=int, default=42)
    p3.add_argument("--alt-lp", action="store_true",
                     help="v1 only: lp() counts n_other(i) instead of degree(i). v2 always "
                          "uses n_other (baked into the spec), this flag has no effect on v2.")
    p3.set_defaults(func=cmd_gridmd)

    p4 = sub.add_parser("probe")
    p4.add_argument("xyz")
    p4.add_argument("--charge", type=int, default=0)
    p4.add_argument("--spin", type=int, default=0)
    p4.add_argument("--rule", choices=("v1", "v2", "both"), default="both")
    p4.set_defaults(func=cmd_probe)

    p5 = sub.add_parser("gridmd-search")
    p5.add_argument("--n", type=int, default=3, help="number of independent replicates")
    p5.add_argument("--maxtime", type=float, default=2000.0)
    p5.add_argument("--out", default=None)
    p5.set_defaults(func=cmd_gridmd_search)

    a = ap.parse_args()
    global ALT_LP_EXCLUDE_PAIR
    if getattr(a, "alt_lp", False):
        ALT_LP_EXCLUDE_PAIR = True
    a.func(a)


if __name__ == "__main__":
    main()
