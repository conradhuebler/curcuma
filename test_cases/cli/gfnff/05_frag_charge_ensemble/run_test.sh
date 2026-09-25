#!/bin/bash

# Test: -gfnff.frag_charge_model ensemble - identity of the (now opt-out) reference rule, and the
#       four properties the ensemble model exists for (test_cases/revgfnff/_log/FRAG_CHARGE_STATUS.md)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 24, 2026) - AI-generated, machine-tested; human production testing pending
# Updated (Sep 24, 2026, default flip - FRAG_CHARGE_STATUS.md section 17): the DEFAULT is now
# frag_charge_model=ensemble at frag_charge_s_max=1.0 (placement rule only, continuous window
# OFF - PARAM help text: "1.0 switches the window off"). Every check below therefore now passes
# REF (-gfnff.frag_charge_model reference) or ENS (ensemble, s_max pinned EXPLICITLY to 1.1, the
# value section 5 recommends for plain GFN-FF and section 5/9/10/11 measured against) instead of
# relying on whatever the compiled-in default happens to be - this keeps testing BOTH the
# placement fix (now the default, no window needed to see it) and the continuous-window behaviour
# (still opt-in) without weakening either assertion when the default changes again later.
#
#   1. IDENTITY  - explicit reference mode reproduces the pinned 12-digit energies of GMTKN55
#                  G21EA/EA_25 (Cl2-) and WATER27/H3OpH2O2 (H3O+(H2O)2) of the package-25 binary
#                  (i.e. the legacy rule is untouched by the default flip).
#   2. LABEL     - Cl2- at the pass-1 split (2.6409 A) with a water H-bonded at either end: the
#                  ensemble (explicit s_max 1.1) gives the same interaction at both ends
#                  (|dE| < 1e-6 kcal/mol); explicit reference must NOT (> 1 kcal/mol, liveness).
#   3. CONTINUITY- across the pass-1 split of Cl2- (r_thr -/+ 1e-4 A) the ensemble (explicit
#                  s_max 1.1, i.e. the window engaged) energy changes by < 0.01 kcal/mol; explicit
#                  reference jumps by > 50 kcal/mol (liveness). At the new default s_max=1.0 the
#                  window is off and the step remains by design (PARAM help text) - not what this
#                  check is testing.
#   4. PLACEMENT - H3O+(H2O)2 (atom 1 is a water O): explicit reference puts +1 on that water, the
#                  ensemble on the hydronium (its three H-bonded-O fragment sums to +1 +- 0.1).
#   5. HISTORY   - Cl2- built at 2.50 A and evaluated at 2.70 A with the topology kept equals the
#                  fresh evaluation to < 0.01 kcal/mol (explicit reference: ~99 kcal/mol apart).
#   6. GRADIENT  - analytic vs central FD (h = 1e-4 A) at Cl2- 2.70 A + water, < 1e-5 Eh/A.
# RE-PINNED Sep 25, 2026 (feature/multi-gpu merge): the H3OpH2O2 identity value moved
# 0.466545744259 -> 0.466545743025 (-1.23e-9 Eh) solely because the reference-less bonded-triple
# ATM term is now off by default (-gfnff.dispersion_atm true restores the old value exactly) -
# the same operator-accepted default change that re-pinned cli_gfnff_04. EA_25 is unaffected.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 05: frag_charge_model ensemble"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import json, math, os, subprocess, sys, tempfile, shutil
B = os.environ["CURCUMA"]
K = 627.509474
REF = ["-gfnff.frag_charge_model", "reference"]
ENS = ["-gfnff.frag_charge_model", "ensemble", "-gfnff.frag_charge_s_max", "1.1"]
fails = 0

def batch(frames, charge, extra, mode="fresh"):
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "s.xyz"), "w") as f:
        for fr in frames:
            f.write(f"{len(fr)}\n\n" + "".join(f"{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n" for a in fr))
    cmd = [B, "-sp", "s.xyz", "-method", "gfnff", "-batch", "true", "-batch_out", "o.jsonl",
           "-gfnff.cache_topology", "false", "-charge", str(charge), "-threads", "1", "-verbosity", "0",
           "-no_bmt", "-gradient", "true"]
    cmd += ["-batch_reuse_topology", "false"] if mode == "fresh" else ["-batch_reuse_topology", "true", "-gfnff.reuse_topology_check", "false"]
    subprocess.run(cmd + extra, cwd=d, capture_output=True, text=True, timeout=600)
    recs = [json.loads(l) for l in open(os.path.join(d, "o.jsonl"))]
    shutil.rmtree(d, ignore_errors=True)
    return recs

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
# 1. identity of the (now explicit) reference rule
e1 = batch([ea25], -1, REF)[0]["energy_eh"]; e2 = batch([h3o], 1, REF)[0]["energy_eh"]
check("identity EA_25 (reference)", abs(e1 - (-0.980160564980)) < 1e-10, f"{e1:.12f} vs -0.980160564980")
check("identity H3OpH2O2 (reference)", abs(e2 - 0.466545743025) < 1e-10, f"{e2:.12f} vs 0.466545743025")

# 2. label symmetry at the split
r = 2.6409
X = [("Cl", 0, 0, 0), ("Cl", 0, 0, r)]
def eint(extra):
    rec = batch([X, X + water(0.0, -1), X + water(r, +1)], -1, extra)
    return (rec[1]["energy_eh"] - rec[2]["energy_eh"]) * K
g_ens, g_ref = eint(ENS), eint(REF)
check("label ensemble", abs(g_ens) < 1e-6, f"end A - end B = {g_ens:+.2e} kcal/mol")
check("label reference (liveness)", abs(g_ref) > 1.0, f"end A - end B = {g_ref:+.2f} kcal/mol")

# 3. continuity across the pass-1 split: locate the reference's jump by bisection, then compare
def ecl(rr, extra):
    return [x["energy_eh"] for x in batch([[("Cl", 0, 0, 0), ("Cl", 0, 0, r0)] for r0 in rr], -1, extra)]
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

# 4. placement by parity
rec = batch([h3o], 1, ENS)[0]
q = rec["charges"]
qh3o = q[6] + q[7] + q[8] + q[9]      # O 7 and its three H (atoms 7-10)
check("placement H3O+", abs(qh3o - 1.0) < 0.1 and (rec["energy_eh"] - e2) * K < -100.0,
      f"q(H3O) = {qh3o:+.3f}, E - E(reference) = {(rec['energy_eh'] - e2) * K:+.1f} kcal/mol")

# 5. topology history
geo = [("Cl", 0, 0, 0), ("Cl", 0.02, -0.01, 2.70)]
ek = batch([[("Cl", 0, 0, 0), ("Cl", 0, 0, 2.50)], geo], -1, ENS, "kept")[1]["energy_eh"]
ef = batch([geo], -1, ENS)[0]["energy_eh"]
check("history ensemble", abs(ek - ef) * K < 0.01, f"kept - fresh = {(ek - ef) * K:+.4f} kcal/mol")

# 6. gradient
geo = [("Cl", 0, 0, 0), ("Cl", 0.02, -0.01, 2.70)] + water(0.0, -1, 2.3)
h = 1e-4
frames = [geo]
for a in range(len(geo)):
    for c in range(3):
        for s in (1, -1):
            g = [list(x) for x in geo]; g[a][1 + c] += s * h; frames.append([tuple(x) for x in g])
recs = batch(frames, -1, ENS)
ga = recs[0]["gradient_eh_ang"]; dev = 0.0; k = 1
ga = [v for row in ga for v in row] if ga and isinstance(ga[0], list) else ga
for a in range(len(geo)):
    for c in range(3):
        fd = (recs[k]["energy_eh"] - recs[k + 1]["energy_eh"]) / (2 * h); k += 2
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
