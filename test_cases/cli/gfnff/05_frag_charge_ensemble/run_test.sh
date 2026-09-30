#!/bin/bash

# Test: -gfnff.frag_charge_model ensemble - identity of the (now opt-out) reference rule, and the
#       properties the ensemble model exists for (Known Issue #31, docs/FRAG_CHARGE_MODEL.md)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 24, 2026) - AI-generated, machine-tested; human production testing pending
#
# Ported to master (Sep 27, 2026) from reactff2-llm. The original drove every evaluation through
# that branch's batch single-point mode (-batch/-batch_out JSONL with per-frame charges and
# -batch_reuse_topology), none of which exists on master. This version runs one plain
# `-sp -dump_gradient` per frame (12-digit energy + analytic gradient in Eh/A) instead; the
# thresholds are unchanged. Two checks had to change shape:
#   - PLACEMENT no longer reads per-atom charges (master's -sp does not export them). It checks
#     the stronger, index-free statement instead: the default (ensemble) energy of H3O+(H2O)2 with
#     a WATER listed first equals the reference-rule energy of the same structure re-ordered so
#     the HYDRONIUM is fragment 0 - i.e. the carrier no longer depends on the atom numbering.
#   - HISTORY (topology built at another geometry) needs a reused calculator instance, which the
#     master CLI cannot express; it lives in the C++ test gfnff_frag_charge_history
#     (test_cases/test_gfnff_frag_charge_history.cpp).
#
# Default is frag_charge_model=ensemble at frag_charge_s_max=1.0 (placement rule only, continuous
# window OFF). REF (-gfnff.frag_charge_model reference) and ENS (ensemble, s_max pinned EXPLICITLY
# to 1.1 so the window is engaged) are therefore passed explicitly wherever the check needs them.
#
#   1. IDENTITY  - explicit reference mode reproduces pinned 12-digit energies of GMTKN55
#                  G21EA/EA_25 (Cl2-) and WATER27/H3OpH2O2 (H3O+(H2O)2): the legacy rule is
#                  untouched by the port. Pinned from master's own binary (see the value comments).
#   2. LABEL     - Cl2- at the pass-1 split (2.6409 A) with a water H-bonded at either end: ENS
#                  gives the same interaction at both ends (|dE| < 1e-6 kcal/mol); REF must NOT
#                  (> 1 kcal/mol, liveness).
#   3. CONTINUITY- across the pass-1 split of Cl2- (r_thr -/+ 1e-4 A) ENS changes by
#                  < 0.01 kcal/mol; REF jumps by > 50 kcal/mol (liveness).
#   4. PLACEMENT - see above, plus E(default) - E(REF, water first) < -100 kcal/mol (liveness).
#   6. GRADIENT  - analytic vs central FD (h = 1e-4 A) at Cl2- 2.70 A + water, < 1e-5 Eh/A.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 05: frag_charge_model ensemble"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import math, os, subprocess, sys, tempfile, shutil
B = os.environ["CURCUMA"]
K = 627.509474
REF = ["-gfnff.frag_charge_model", "reference"]
ENS = ["-gfnff.frag_charge_model", "ensemble", "-gfnff.frag_charge_s_max", "1.1"]
DEF = []   # compiled-in default (ensemble, s_max 1.0)
fails = 0

def sp(frame, charge, extra):
    """One single point; returns (energy_eh, flat gradient in Eh/A)."""
    d = tempfile.mkdtemp()
    try:
        with open(os.path.join(d, "s.xyz"), "w") as f:
            f.write(f"{len(frame)}\n\n" + "".join(f"{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n" for a in frame))
        cmd = [B, "-sp", "s.xyz", "-method", "gfnff", "-gfnff.cache_topology", "false",
               "-charge", str(charge), "-threads", "1", "-verbosity", "0", "-no_bmt",
               "-dump_gradient", "g.txt"] + extra
        subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=600)
        e, g = None, []
        for l in open(os.path.join(d, "g.txt")):
            if l.startswith("# energy"):
                e = float(l.split()[2])
            elif not l.startswith("#") and l.strip():
                g += [float(x) for x in l.split()]
        return e, g
    finally:
        shutil.rmtree(d, ignore_errors=True)

def energies(frames, charge, extra):
    return [sp(fr, charge, extra)[0] for fr in frames]

def check(name, ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + f"{name}: {msg}")
    if not ok: fails += 1

def rd(p):
    L = open(p).read().split("\n"); n = int(L[0])
    return [(l.split()[0], *map(float, l.split()[1:4])) for l in L[2:2 + n]]

def water(end_z, direction, dXH=2.3):
    zH = end_z + direction * dXH; zO = zH + direction * 0.97; th = math.radians(104.5)
    return [("H", 0.0, 0.0, zH), ("O", 0.0, 0.0, zO), ("H", 0.97 * math.sin(th), 0.0, zO - direction * 0.97 * math.cos(th))]

ea25 = rd("cl2m_ea25.xyz"); h3o = rd("h3op_h2o2.xyz")
# 1. identity of the (now explicit) reference rule. Values pinned from master's own binary
#    (Sep 27, 2026); EA_25 equals reactff2-llm's pin, H3OpH2O2 equals reactff2-llm's pre-Sep-25
#    value (master keeps the bonded-triple ATM term that reactff2-llm switched off).
#    Re-pinned Sep 30, 2026 after merging origin/master's CODATA-2018 unit-constants
#    unification (Known Issue #35: gfnff shift <= 1.5e-4 kcal/mol ~= 2.4e-7 Eh) - both values
#    moved by ~2-4e-9 Eh, well inside that documented bound.
H3O_PIN = 0.466545742354   # master = reactff2-llm before its Sep-25 ATM default flip (after it: 0.466545743025)
e1 = sp(ea25, -1, REF)[0]; e2 = sp(h3o, 1, REF)[0]
check("identity EA_25 (reference)", abs(e1 - (-0.980160561416)) < 1e-10, f"{e1:.12f} vs -0.980160561416")
check("identity H3OpH2O2 (reference)", abs(e2 - H3O_PIN) < 1e-10, f"{e2:.12f} vs {H3O_PIN:.12f}")
# 2. label symmetry at the split
r = 2.6409
X = [("Cl", 0, 0, 0), ("Cl", 0, 0, r)]
def eint(extra):
    e = energies([X + water(0.0, -1), X + water(r, +1)], -1, extra)
    return (e[0] - e[1]) * K
g_ens, g_ref = eint(ENS), eint(REF)
check("label ensemble", abs(g_ens) < 1e-6, f"end A - end B = {g_ens:+.2e} kcal/mol")
check("label reference (liveness)", abs(g_ref) > 1.0, f"end A - end B = {g_ref:+.2f} kcal/mol")

# 3. continuity across the pass-1 split: locate the reference's jump by bisection, then compare
def ecl(rr, extra):
    return energies([[("Cl", 0, 0, 0), ("Cl", 0, 0, r0)] for r0 in rr], -1, extra)
lo, hi = 2.55, 2.66
elo, ehi = ecl([lo, hi], REF)
for _ in range(24):
    mid = 0.5 * (lo + hi); em = ecl([mid], REF)[0]
    if abs(em - elo) < abs(em - ehi): lo, elo = mid, em
    else: hi, ehi = mid, em
rs = 0.5 * (lo + hi)
er = ecl([rs - 1e-4, rs + 1e-4], REF); ee = ecl([rs - 1e-4, rs + 1e-4], ENS)
j_ref, j_ens = abs(er[1] - er[0]) * K, abs(ee[1] - ee[0]) * K
check("continuity ensemble", j_ens < 0.01, f"split at {rs:.6f} A, |dE| over 2e-4 A = {j_ens:.4f} kcal/mol")
check("continuity reference (liveness)", j_ref > 50.0, f"|dE| over 2e-4 A = {j_ref:.2f} kcal/mol")

# 4. placement by parity, index-free: hydronium (atoms 7-10) written first -> the reference rule
#    puts the +1 on it; the default must give the same energy on the original (water-first) file.
h3o_first = h3o[6:10] + h3o[0:6]
e_def = sp(h3o, 1, DEF)[0]; e_ref_first = sp(h3o_first, 1, REF)[0]
check("placement H3O+ (default == reference with hydronium first)", abs(e_def - e_ref_first) * K < 1e-6,
      f"{(e_def - e_ref_first) * K:+.2e} kcal/mol")
check("placement H3O+ (liveness vs water-first reference)", (e_def - e2) * K < -100.0,
      f"E(default) - E(reference, water first) = {(e_def - e2) * K:+.1f} kcal/mol")

# 6. gradient
geo = [("Cl", 0, 0, 0), ("Cl", 0.02, -0.01, 2.70)] + water(0.0, -1, 2.3)
h = 1e-4
ga = sp(geo, -1, ENS)[1]; dev = 0.0
for a in range(len(geo)):
    for c in range(3):
        ep, em = [], []
        for s_, out in ((1, ep), (-1, em)):
            g = [list(x) for x in geo]; g[a][1 + c] += s_ * h
            out.append(sp([tuple(x) for x in g], -1, ENS)[0])
        fd = (ep[0] - em[0]) / (2 * h)
        dev = max(dev, abs(ga[3 * a + c] - fd))
check("gradient ensemble", dev < 1e-5, f"max |g_an - g_fd| = {dev:.2e} Eh/A")
sys.exit(1 if fails else 0)
PYEOF
    local rc=$?
    if [ $rc -eq 0 ]; then TESTS_PASSED=$((TESTS_PASSED + 1)); else TESTS_FAILED=$((TESTS_FAILED + 1)); fi
    TESTS_RUN=$((TESTS_RUN + 1))
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
