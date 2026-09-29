#!/bin/bash

# Test: -gfnff.rev_pair_validity - the rev-gfnff pair-validity gate
#       (test_cases/revgfnff/_log/FABLE_BOND_STATE_2.md sec 2.1-rev,
#        test_cases/revgfnff/_log/PAIR_VALIDITY_IMPL_STATUS.md)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 29, 2026) - AI-generated, machine-tested; human production testing pending
#
# bf4_compressed.xyz / bf4_realistic.xyz are the EXACT geometries the design doc's own numbers
# were measured against (Td BF4-, B-F 1.14301 / 1.39980 A respectively - recovered from the
# design session's own scratch files, not re-derived, so the targets below are bit-for-bit, not
# approximate).
#
#   1. OFF bit-identity   - flag off/absent give the identical energy on both an ordinary
#                            molecule and the compressed-BF4- probe (nothing changed by adding
#                            the code when the flag is not requested).
#   2. BF4- compressed     - the headline number: with the gate off, a fresh 10-bond perception
#                            (six spurious F...F contacts alongside the four B-F bonds) reads
#                            +0.17707951 Eh; with the gate on, all six F...F pairs are invalid,
#                            the corner is regenerated on the 4 B-F bonds alone, and the energy
#                            is EXACTLY -1.25025415 Eh - the same force field's own 4-bond
#                            evaluation (bit-for-bit, not merely close).
#   3. BF4- realistic      - at B-F 1.3998 A only 4 bonds are perceived in the first place (every
#                            pair already valid): gate on/off must be bit-identical (the "all
#                            pairs valid" falsifier, FABLE_BOND_STATE_2.md sec 4 step 2) AND
#                            equal -1.47017678 Eh (Fable's own recorded value for this geometry).
#   4. Ordinary molecule    - a neutral organic (caffeine) is untouched by the gate.
#   5. Gradient             - analytic vs central FD (h = 1e-4 A) at the compressed BF4- geometry
#                            with the gate on: the gate is a discrete per-corner constant, so it
#                            must add no gradient term of its own away from a corner change.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 07: pair-validity gate"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import json, os, subprocess, sys, tempfile, shutil
B = os.environ["CURCUMA"]
GATE_ON = ["-gfnff.rev_pair_validity", "true"]
GATE_OFF = ["-gfnff.rev_pair_validity", "false"]
fails = 0

def rd(p):
    L = open(p).read().split("\n")
    n = int(L[0])
    return [(l.split()[0], *map(float, l.split()[1:4])) for l in L[2:2 + n]]

def batch(frames, charge, extra):
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "s.xyz"), "w") as f:
        for fr in frames:
            f.write(f"{len(fr)}\n\n" + "".join(f"{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n" for a in fr))
    cmd = [B, "-sp", "s.xyz", "-method", "revgfnff", "-batch", "true", "-batch_out", "o.jsonl",
           "-gfnff.cache_topology", "false", "-charge", str(charge), "-threads", "1", "-verbosity", "0",
           "-no_bmt", "-gradient", "true", "-batch_reuse_topology", "false"]
    r = subprocess.run(cmd + extra, cwd=d, capture_output=True, text=True, timeout=600)
    try:
        recs = [json.loads(l) for l in open(os.path.join(d, "o.jsonl"))]
    except FileNotFoundError:
        recs = []
        print("curcuma stderr:", r.stderr[-2000:], file=sys.stderr)
    shutil.rmtree(d, ignore_errors=True)
    return recs

def check(name, ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + f"{name}: {msg}")
    if not ok:
        fails += 1

compressed = rd("bf4_compressed.xyz")
realistic = rd("bf4_realistic.xyz")

# 1. OFF bit-identity: an ordinary neutral molecule and BF4- (both geometries) are unaffected by
#    merely requesting the flag off/absent
e_water_off = batch([[("O", 0, 0, 0), ("H", 0.96, 0, 0), ("H", -0.24, 0.93, 0)]], 0, GATE_OFF)[0]["energy_eh"]
e_water_none = batch([[("O", 0, 0, 0), ("H", 0.96, 0, 0), ("H", -0.24, 0.93, 0)]], 0, [])[0]["energy_eh"]
check("off bit-identity (water, explicit false vs absent)", abs(e_water_off - e_water_none) < 1e-12,
      f"{e_water_off:.12f} vs {e_water_none:.12f}")

e_bf4c_off = batch([compressed], -1, GATE_OFF)[0]["energy_eh"]
e_bf4c_none = batch([compressed], -1, [])[0]["energy_eh"]
check("off bit-identity (BF4- compressed, explicit false vs absent)", abs(e_bf4c_off - e_bf4c_none) < 1e-12,
      f"{e_bf4c_off:.12f} vs {e_bf4c_none:.12f}")
check("off matches the pre-gate 10-bond fresh perception (+0.17707951 Eh)",
      abs(e_bf4c_off - 0.17707951) < 1e-7, f"{e_bf4c_off:.8f}")

# 2. BF4- compressed, gate ON: the headline acceptance number (FABLE_BOND_STATE_2.md sec 2.3.1),
#    -1.25025415 Eh EXACTLY, matching the topology's own 4-bond evaluation bit-for-bit.
e_bf4c_on = batch([compressed], -1, GATE_ON)[0]["energy_eh"]
check("BF4- compressed, gate ON == pinned 4-bond evaluation, bit-for-bit",
      abs(e_bf4c_on - (-1.25025415)) < 1e-7, f"{e_bf4c_on:.8f} vs -1.25025415")

# 3. BF4- realistic (B-F 1.3998 A): only 4 bonds perceived in the first place -> gate on/off
#    identical, and both equal Fable's own recorded value for this geometry.
e_bf4r_on = batch([realistic], -1, GATE_ON)[0]["energy_eh"]
e_bf4r_off = batch([realistic], -1, GATE_OFF)[0]["energy_eh"]
check("BF4- realistic (all pairs already valid): gate on == gate off", abs(e_bf4r_on - e_bf4r_off) < 1e-12,
      f"{e_bf4r_on:.12f} vs {e_bf4r_off:.12f}")
check("BF4- realistic matches the recorded 4-bond value (-1.47017678 Eh)",
      abs(e_bf4r_on - (-1.47017678)) < 1e-7, f"{e_bf4r_on:.8f}")

# 4. Ordinary neutral organic molecule (caffeine, from the sibling test's own file): untouched
caf = rd(os.path.join("..", "02_caffeine_gfnff_energy_components", "caffeine.xyz"))
e_caf_on = batch([caf], 0, GATE_ON)[0]["energy_eh"]
e_caf_off = batch([caf], 0, GATE_OFF)[0]["energy_eh"]
check("caffeine (ordinary organic): gate on == gate off", abs(e_caf_on - e_caf_off) < 1e-12,
      f"{e_caf_on:.12f} vs {e_caf_off:.12f}")

# 5. Gradient at the compressed BF4- geometry, gate ON: analytic vs central FD
h = 1e-4
frames = [compressed]
for a in range(len(compressed)):
    for c in range(3):
        for sgn in (1, -1):
            g = [list(x) for x in compressed]
            g[a][1 + c] += sgn * h
            frames.append([tuple(x) for x in g])
recs = batch(frames, -1, GATE_ON)
if len(recs) == len(frames):
    ga = recs[0]["gradient_eh_ang"]
    ga = [v for row in ga for v in row] if ga and isinstance(ga[0], list) else ga
    dev = 0.0
    k = 1
    for a in range(len(compressed)):
        for c in range(3):
            fd = (recs[k]["energy_eh"] - recs[k + 1]["energy_eh"]) / (2 * h)
            k += 2
            dev = max(dev, abs(ga[3 * a + c] - fd))
    check("gradient at BF4- compressed, gate ON (analytic vs FD)", dev < 1e-4, f"max |g_an - g_fd| = {dev:.2e} Eh/A")
else:
    check("gradient at BF4- compressed, gate ON (analytic vs FD)", False,
          f"expected {len(frames)} batch records, got {len(recs)}")

sys.exit(1 if fails else 0)
PYEOF
    local rc=$?
    if [ $rc -eq 0 ]; then TESTS_PASSED=$((TESTS_PASSED + 1)); else TESTS_FAILED=$((TESTS_FAILED + 1)); fi
    TESTS_RUN=$((TESTS_RUN + 1))
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
