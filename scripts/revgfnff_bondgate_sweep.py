#!/usr/bin/env python3
"""Offline sweep for the FABLE_BOND_STATE_2.md section 2.1 bond-validity gate.

NO C++ CHANGE. This evaluates, in Python, the proposed rule

    VALID(i, j, b)  iff   n_other(i,b) < Val_i(b)
                     or   n_other(j,b) < Val_j(b)
                     or   (Z_i == 1 and lp(j, b))
                     or   (Z_j == 1 and lp(i, b))

against the topology/budget data curcuma ALREADY prints via the env-var-gated debug dumps
CURCUMA_BONDDUMP (per-bond topology: element, index, cached bond list) and CURCUMA_SHAREDUMP
(per-corner share table: for EVERY atom, Z, the corner's claim sum S, nominal valence ValZ, the
charge/donor-rule cap `cap`, and the resulting effective valence `Val` = Val_i(b) -- i.e. the
`prepareConservingShare()` output the gate rule needs is read directly off the C++, not
re-derived). See test_cases/revgfnff/_log/FABLE_BOND_STATE_2.md section 4, measurement plan
step 1, and test_cases/revgfnff/_log/BOND_VALIDITY_GATE_SWEEP_STATUS.md for the results.

CURCUMA_SHAREDUMP output only exists when rev-gfnff is engaged (prepareValenceShare() is gated
on m_rev.enabled), so every run in this script uses `-method revgfnff` at ITS SHIPPED DEFAULTS
(share_form conserving, share_donor_rule true, well_form mg3, budget_fix_h true) -- not a
special measurement configuration. The discrete bond LIST (perceiveGeometricBonds()) is shared
code, independent of the rev flag, so this is the same topology plain `-method gfnff` would
perceive; only the per-atom Val_i(b) budget needs rev mode to be computed at all. This is a
genuine limitation of the existing dumps (there is no plain-gfnff equivalent of the budget
cap), stated here rather than worked around.

Usage:
    scripts/revgfnff_bondgate_sweep.py refset [--limit N] [--jobs N]      # GMTKN55+MOR41+S30L-CI
    scripts/revgfnff_bondgate_sweep.py grid   [--jobs N]                  # 130-cell grid frames
"""
import argparse
import concurrent.futures as cf
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
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


ALT_LP_EXCLUDE_PAIR = False  # robustness check: count = n_other(i) instead of full degree(i),
                             # i.e. lp() ignores the very pair being tested. See report section
                             # on the lp() self-reference ambiguity (ch3nh2 hot-MD finding).


def lone_pair(z, count):
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


def parse_and_gate(output):
    """Parse one curcuma run's stdout for the LAST share-dump corner and evaluate the gate on
    every pair it lists. Returns a dict with keys: n_bonds_share, n_bonds_bonddump, atoms
    (idx -> dict(Z, Val)), pairs (list of (i,j)), results (list of dict per pair), and any
    parse-failure notes."""
    lines = clean_lines(output)
    bonddump_pairs = set()
    bonddump_Z = {}
    for ln in lines:
        m = BOND_RE.match(ln)
        if m:
            i, zi, j, zj = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))
            bonddump_pairs.add((min(i, j), max(i, j)))
            bonddump_Z[i] = zi
            bonddump_Z[j] = zj

    # Take the LAST "share dump: corner with N bonds" block (the converged-charge corner).
    corner_starts = [k for k, ln in enumerate(lines) if CORNER_HDR_RE.match(ln)]
    share_pairs = []  # (i, j, Val_i, Val_j) as printed, 1-indexed
    shareA = {}       # idx -> (Z, S, ValZ, cap, Val)
    if corner_starts:
        start = corner_starts[-1]
        end = corner_starts[-1] + 1
        # consume "share" lines then "shareA" lines until a non-matching line
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

    notes = []
    if not corner_starts:
        return dict(ok=False, note="no share-dump corner (0 bonds / single atom, or parse miss)",
                    n_bonds_bonddump=len(bonddump_pairs), pairs_tested=0, results=[])
    if len(bonddump_pairs) and len(share_pairs) != len(bonddump_pairs):
        notes.append("share corner has %d bonds, BONDDUMP topology has %d"
                     % (len(share_pairs), len(bonddump_pairs)))

    # degree in the corner = number of listed partners per atom index
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
        lp_i = lone_pair(zi, lp_count_i) if zi is not None else None
        lp_j = lone_pair(zj, lp_count_j) if zj is not None else None
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
    return dict(ok=True, note="; ".join(notes), n_bonds_bonddump=len(bonddump_pairs),
                pairs_tested=len(results), results=results)


def elem_pair_key(zi, zj):
    si, sj = SYMBOL.get(zi, "Z%d" % zi if zi else "?"), SYMBOL.get(zj, "Z%d" % zj if zj else "?")
    return "-".join(sorted((si, sj)))


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


def process_one(label, xyz, charge, spin):
    with tempfile.TemporaryDirectory() as td:
        local = Path(td) / "struc.xyz"
        shutil.copy(xyz, local)
        out = run_curcuma(local, td, charge=charge, spin=spin)
    if out is None:
        return label, dict(ok=False, note="TIMEOUT", pairs_tested=0, results=[])
    return label, parse_and_gate(out)


def dataset_of(label):
    return label.split("/")[0]


def summarize(all_results, title, out_lines):
    n_struct = len(all_results)
    n_fail = sum(1 for _, r in all_results if not r["ok"])
    n_pairs = sum(r["pairs_tested"] for _, r in all_results)
    invalid = []
    per_ds_invalid = Counter()
    per_ds_pairs = Counter()
    per_ds_struct_with_invalid = Counter()
    elem_pair_counts = Counter()
    undefined_lp_elems = Counter()
    mismatch_notes = []
    for label, r in all_results:
        ds = dataset_of(label)
        per_ds_pairs[ds] += r["pairs_tested"]
        if not r["ok"]:
            continue
        if r["note"]:
            mismatch_notes.append((label, r["note"]))
        struct_has_invalid = False
        for pr in r["results"]:
            if pr["undefined_lp_used"]:
                z = pr["zi"] if pr["zi"] != 1 else pr["zj"]
                undefined_lp_elems[SYMBOL.get(z, "Z%d" % z)] += 1
            if not pr["valid"]:
                invalid.append((label, pr))
                per_ds_invalid[ds] += 1
                struct_has_invalid = True
                elem_pair_counts[elem_pair_key(pr["zi"], pr["zj"])] += 1
        if struct_has_invalid:
            per_ds_struct_with_invalid[ds] += 1

    out_lines.append("## %s" % title)
    out_lines.append("")
    out_lines.append("Structures evaluated: %d (failed/no-parse: %d)" % (n_struct, n_fail))
    out_lines.append("Pairs tested (listed bonds in the evaluated corner): %d" % n_pairs)
    out_lines.append("INVALID pairs total: %d" % len(invalid))
    out_lines.append("")
    out_lines.append("| dataset | pairs tested | INVALID pairs | structures with >=1 INVALID |")
    out_lines.append("|---|---:|---:|---:|")
    for ds in sorted(per_ds_pairs):
        out_lines.append("| %s | %d | %d | %d |" % (
            ds, per_ds_pairs[ds], per_ds_invalid.get(ds, 0), per_ds_struct_with_invalid.get(ds, 0)))
    out_lines.append("")
    if elem_pair_counts:
        out_lines.append("INVALID pairs by element pair:")
        out_lines.append("")
        out_lines.append("| element pair | count |")
        out_lines.append("|---|---:|")
        for k, v in elem_pair_counts.most_common():
            out_lines.append("| %s | %d |" % (k, v))
        out_lines.append("")
    else:
        out_lines.append("No INVALID pairs found.")
        out_lines.append("")
    if undefined_lp_elems:
        out_lines.append("H-bridge clause evaluated against an element with NO main-group "
                          "lp() rule (d/f-block; lp() treated as unknown, bridge clause could "
                          "not fire either way) -- counts of (H, X) pairs hitting this:")
        out_lines.append("")
        for k, v in undefined_lp_elems.most_common():
            out_lines.append("- %s: %d" % (k, v))
        out_lines.append("")
    if mismatch_notes:
        out_lines.append("Corner/BONDDUMP bond-count mismatches (first 20):")
        for label, note in mismatch_notes[:20]:
            out_lines.append("- %s: %s" % (label, note))
        out_lines.append("")
    out_lines.append("All INVALID pairs (structure, i-j, Z pair, n_other/Val each end):")
    out_lines.append("")
    for label, pr in invalid:
        out_lines.append("- %s  %d(%s)-%d(%s)  n_other/Val: %.3f/%.4f vs %.3f/%.4f" % (
            label, pr["i"], SYMBOL.get(pr["zi"], pr["zi"]),
            pr["j"], SYMBOL.get(pr["zj"], pr["zj"]),
            pr["n_other_i"], pr["val_i"], pr["n_other_j"], pr["val_j"]))
    out_lines.append("")
    return invalid


def cmd_refset(a):
    jobs = collect_refset(limit=a.limit)
    print("refset gate sweep: %d structures, %d jobs" % (len(jobs), a.jobs), flush=True)
    all_results = []
    done = 0
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(process_one, label, xyz, c, s): label
                for (label, xyz, c, s) in jobs}
        for fut in cf.as_completed(futs):
            label = futs[fut]
            try:
                lbl, r = fut.result()
            except Exception as e:
                lbl, r = label, dict(ok=False, note="EXC: %s" % e, pairs_tested=0, results=[])
            all_results.append((lbl, r))
            done += 1
            if done % 200 == 0:
                print("  ... %d/%d" % (done, len(jobs)), flush=True)
    out_lines = []
    summarize(all_results, "Reference-set sweep (GMTKN55 + MOR41 + S30L-CI)", out_lines)
    print("\n".join(out_lines))
    if a.out:
        Path(a.out).write_text("\n".join(out_lines))


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


def process_grid_one(label, frame_lines):
    with tempfile.TemporaryDirectory() as td:
        local = Path(td) / "frame.xyz"
        local.write_text("\n".join(frame_lines) + "\n")
        out = run_curcuma(local, td, charge=0, spin=0)
    if out is None:
        return label, dict(ok=False, note="TIMEOUT", pairs_tested=0, results=[])
    return label, parse_and_gate(out)


REACT_EVENT_RE = re.compile(r"^REACT (bond formed|bond broken|rebuild) #?(\d*)")


def gate_one_corner(share_pairs, shareA, bonddump_Z=None):
    """Evaluate the section-2.1 gate on one already-extracted corner. share_pairs: list of
    (i, j, val_i_pair, val_j_pair). shareA: idx -> dict(Z, S, ValZ, cap, Val)."""
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
        lp_i = lone_pair(zi, lp_count_i) if zi is not None else None
        lp_j = lone_pair(zj, lp_count_j) if zj is not None else None
        bridge = (zi == 1 and lp_j) or (zj == 1 and lp_i)
        valid = free_i or free_j or bridge
        results.append(dict(i=i, j=j, zi=zi, zj=zj, val_i=val_i, val_j=val_j,
                             n_other_i=n_other_i, n_other_j=n_other_j, valid=valid))
    return results


def run_reactive_md(system, frame_lines, temperature, maxtime, dt, workdir):
    local = Path(workdir) / "input.xyz"
    local.write_text("\n".join(frame_lines) + "\n")
    cmd = [str(CURCUMA), "-md", "input.xyz", "-method", "revgfnff",
           "-gfnff.topology_mode", "react", "-temperature", str(temperature),
           "-maxtime", str(maxtime), "-md.time_step", str(dt),
           "-md.thermostat", "csvr", "-md.coupling", "10",
           "-md.rattle_12", "false", "-md.no_restart", "-md.seed", "42",
           "-threads", "1", "-verbosity", "2", "-md.print_frequency", "1", "-no_bmt"]
    env = dict(os.environ)
    env["CURCUMA_BONDDUMP"] = "1"
    env["CURCUMA_SHAREDUMP"] = "1"
    r = subprocess.run(cmd, cwd=str(workdir), capture_output=True, text=True,
                        timeout=600, env=env)
    return strip_ansi(r.stdout) + strip_ansi(r.stderr)


def scan_trajectory(output, window=3):
    """Walk a reactive-MD log line by line, evaluating the gate on EVERY per-step corner
    (not just the last), and tag whether each corner falls within `window` corner-blocks of a
    REACT bond-formed/broken/rebuild line (a genuine topology TRANSITION corner) or not (an
    ordinary settled corner). Returns aggregate counts and a sample of INVALID hits."""
    lines = clean_lines(output)
    corner_idx = [k for k, ln in enumerate(lines) if CORNER_HDR_RE.match(ln)]
    event_idx = [k for k, ln in enumerate(lines) if REACT_EVENT_RE.match(ln)]

    all_invalid = []
    near_invalid = 0
    near_total_pairs = 0
    far_invalid = 0
    far_total_pairs = 0
    n_corners = 0
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
            shareA[idx] = dict(Z=int(m.group(2)), Val=float(m.group(6)))
            k += 1
        if not share_pairs:
            continue
        n_corners += 1
        # nearest REACT event line index (by line distance) to this corner block
        while ei < len(event_idx) and event_idx[ei] < start - 200:
            ei += 1
        near = any(abs(event_idx[j] - start) < 40 for j in range(ei, len(event_idx))
                   if event_idx[j] < end + 200 and event_idx[j] >= start - 200)
        res = gate_one_corner(share_pairs, shareA)
        for pr in res:
            if near:
                near_total_pairs += 1
            else:
                far_total_pairs += 1
            if not pr["valid"]:
                all_invalid.append((near, pr))
                if near:
                    near_invalid += 1
                else:
                    far_invalid += 1
    return dict(n_corners=n_corners, n_rebuild_events=sum(1 for k in event_idx
                                                           if lines[k].startswith("REACT rebuild")),
                near_total_pairs=near_total_pairs, near_invalid=near_invalid,
                far_total_pairs=far_total_pairs, far_invalid=far_invalid,
                invalid_sample=all_invalid[:20])


def cmd_gridmd(a):
    systems = [("c2h6", "c2h6_1000K.xyz"), ("ch3nh2", "ch3nh2_1000K.xyz"), ("ch4_H", "ch4_H.xyz")]
    out_lines = ["## 130-cell grid, REAL reactive-MD trajectories (rebuild corners, not statics)",
                 "", "One unperturbed frame-0 trajectory per system, T=2000 K, true dt=0.25 fs, "
                 "maxtime=%d fs, topology_mode=react, shipped rev-gfnff defaults "
                 "(mg3/conserving/donor_rule/budget_fix_h). NOT the full 130-cell x 6-replicate "
                 "protocol of package 11 -- 3 trajectories, chosen to actually pass through "
                 "REACT rebuild events (unlike the static substitute above) rather than "
                 "reproduce the full statistic." % a.maxtime, ""]
    for system, fname in systems:
        frames = read_frames(FITWORK / fname)
        with tempfile.TemporaryDirectory() as td:
            out = run_reactive_md(system, frames[0], 2000, a.maxtime, 0.25, td)
        r = scan_trajectory(out)
        out_lines.append("### %s (frame 0, T=2000K, %d fs)" % (system, a.maxtime))
        out_lines.append("- corners scanned: %d, REACT rebuild events: %d"
                          % (r["n_corners"], r["n_rebuild_events"]))
        out_lines.append("- pairs near a topology-transition event (within ~10 fs): %d, INVALID: %d"
                          % (r["near_total_pairs"], r["near_invalid"]))
        out_lines.append("- pairs away from any transition event: %d, INVALID: %d"
                          % (r["far_total_pairs"], r["far_invalid"]))
        if r["invalid_sample"]:
            out_lines.append("- sample INVALID hits (near?, i-j, Z pair):")
            for near, pr in r["invalid_sample"]:
                out_lines.append("    near=%s  %d(%s)-%d(%s)  n_other/Val %.2f/%.3f vs %.2f/%.3f" % (
                    near, pr["i"], SYMBOL.get(pr["zi"], pr["zi"]),
                    pr["j"], SYMBOL.get(pr["zj"], pr["zj"]),
                    pr["n_other_i"], pr["val_i"], pr["n_other_j"], pr["val_j"]))
        out_lines.append("")
    print("\n".join(out_lines))
    if a.out:
        Path(a.out).write_text("\n".join(out_lines))


def cmd_grid(a):
    jobs = collect_grid()
    print("130-cell grid substitute: %d DISTINCT starting geometries "
          "(temperature does not change the t=0 geometry in the tail-sweep protocol, "
          "so a static gate check covers 65 of the nominal 130 cells; see report)"
          % len(jobs), flush=True)
    all_results = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for label, r in ex.map(lambda j: process_grid_one(*j), jobs):
            all_results.append((label, r))
    out_lines = []
    summarize(all_results, "130-cell grid substitute (65 distinct starting frames, static)",
              out_lines)
    print("\n".join(out_lines))
    if a.out:
        Path(a.out).write_text("\n".join(out_lines))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("refset")
    p1.add_argument("--limit", type=int, default=0)
    p1.add_argument("--jobs", type=int, default=os.cpu_count())
    p1.add_argument("--out", default=None)
    p1.set_defaults(func=cmd_refset)

    p2 = sub.add_parser("grid")
    p2.add_argument("--jobs", type=int, default=os.cpu_count())
    p2.add_argument("--out", default=None)
    p2.set_defaults(func=cmd_grid)

    p3 = sub.add_parser("grid-md")
    p3.add_argument("--maxtime", type=float, default=2000.0)
    p3.add_argument("--out", default=None)
    p3.add_argument("--alt-lp", action="store_true",
                     help="lp() counts n_other(i) instead of degree(i) (excludes the pair "
                          "under test from its own lone-pair count)")
    p3.set_defaults(func=cmd_gridmd)

    a = ap.parse_args()
    global ALT_LP_EXCLUDE_PAIR
    if getattr(a, "alt_lp", False):
        ALT_LP_EXCLUDE_PAIR = True
    a.func(a)


if __name__ == "__main__":
    main()
