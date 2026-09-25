#!/bin/bash

# Test: energy-only GFN-FF calls on a REUSED calculator must use the CN of their own geometry
#       (test_cases/revgfnff/_log/STALE_CN_STATUS.md, fix A)
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Sep 24, 2026) - AI-generated, machine-tested; human production testing pending
#
# Before the fix the Coulomb self-energy's chi(CN) = chi_base + cnf*sqrt(CN) term read a CN that
# only gradient calls refreshed, so an energy-only call on a reused calculator (batch reuse, FD
# tests, energy-based line searches / Hessians) evaluated it with the CN of another geometry.
#
#   1. FD      - Cl2- (q=-1) at 2.73 A: analytic gradient vs central FD whose energies are
#                energy-only calls on ONE reused calculator, < 1e-4 Eh/A.
#                (pre-fix binary: 1.77e-2; after fix A 2.41e-5, the remainder is the frozen D4
#                C6 of STALE_CN_STATUS.md section 3, not addressed by fix A)
#   2. E vs G  - on a reused calculator, an energy-only call and a gradient call at the same
#                geometry give the same energy (< 1e-10 Eh): Cl2- built at 2.73 A, evaluated at
#                2.65 A (pre-fix: 1.4e-3 Eh apart), and caffeine perturbed by up to 0.05 A.
#   3. COULOMB - the reused energy-only Coulomb term at 2.65 A equals the fresh one (< 1e-9 Eh;
#                pre-fix: 1.4e-3 Eh).

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 06: stale CN in energy-only calls"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import json, os, random, subprocess, sys, tempfile, shutil
B = os.environ["CURCUMA"]
fails = 0

def batch(frames, charge, grad, reuse):
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "s.xyz"), "w") as f:
        for fr in frames:
            f.write(f"{len(fr)}\n\n" + "".join(f"{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n" for a in fr))
    cmd = [B, "-sp", "s.xyz", "-method", "gfnff", "-batch", "true", "-batch_out", "o.jsonl",
           "-gfnff.cache_topology", "false", "-charge", str(charge), "-threads", "1", "-verbosity", "0",
           "-no_bmt", "-gradient", "true" if grad else "false"]
    cmd += ["-batch_reuse_topology", "true", "-gfnff.reuse_topology_check", "false"] if reuse else ["-batch_reuse_topology", "false"]
    subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=600)
    recs = [json.loads(l) for l in open(os.path.join(d, "o.jsonl"))]
    shutil.rmtree(d, ignore_errors=True)
    return recs

def check(name, ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + f"{name}: {msg}")
    if not ok: fails += 1

# 1. FD on a reused calculator, energy-only
geo = [("Cl", 0.0, 0.0, 0.0), ("Cl", 2.73, 0.0, 0.0)]
h = 1e-4
frames = [geo]
for a in range(len(geo)):
    for c in range(3):
        for s in (1, -1):
            g = [list(x) for x in geo]; g[a][1 + c] += s * h; frames.append([tuple(x) for x in g])
ga = batch([geo], -1, True, False)[0]["gradient_eh_ang"]
recs = batch(frames, -1, False, True)
dev = 0.0; k = 1
for a in range(len(geo)):
    for c in range(3):
        fd = (recs[k]["energy_eh"] - recs[k + 1]["energy_eh"]) / (2 * h); k += 2
        dev = max(dev, abs(ga[a][c] - fd))
check("FD Cl2- 2.73 A, reused energy-only", dev < 1e-4, f"max |g_an - g_fd| = {dev:.2e} Eh/A (pre-fix 1.77e-2)")

# 2./3. energy-only vs gradient call on a reused calculator, and Coulomb vs fresh
pair = [geo, [("Cl", 0.0, 0.0, 0.0), ("Cl", 2.65, 0.0, 0.0)]]
eE = batch(pair, -1, False, True)[1]; eG = batch(pair, -1, True, True)[1]; eF = batch(pair[1:], -1, False, False)[0]
check("E vs G Cl2- 2.73 -> 2.65 A", abs(eE["energy_eh"] - eG["energy_eh"]) < 1e-10,
      f"|E_energy-only - E_gradient| = {abs(eE['energy_eh'] - eG['energy_eh']):.2e} Eh (pre-fix 1.4e-3)")
dc = abs(eE["terms"]["Coulomb"] - eF["terms"]["Coulomb"])
check("Coulomb reused vs fresh Cl2- 2.65 A", dc < 1e-9, f"|dCoulomb| = {dc:.2e} Eh (pre-fix 1.4e-3)")

L = open("caffeine.xyz").read().split("\n"); n = int(L[0])
caf = [(l.split()[0], *map(float, l.split()[1:4])) for l in L[2:2 + n]]
rnd = random.Random(7)
frs = [caf] + [[(a[0], a[1] + rnd.uniform(-0.05, 0.05), a[2] + rnd.uniform(-0.05, 0.05), a[3] + rnd.uniform(-0.05, 0.05)) for a in caf] for _ in range(3)]
rE = batch(frs, 0, False, True); rG = batch(frs, 0, True, True)
d = max(abs(x["energy_eh"] - y["energy_eh"]) for x, y in zip(rE[1:], rG[1:]))
check("E vs G caffeine (3 perturbed frames)", d < 1e-10, f"max |E_energy-only - E_gradient| = {d:.2e} Eh")
sys.exit(1 if fails else 0)
PYEOF
    local rc=$?
    if [ $rc -eq 0 ]; then TESTS_PASSED=$((TESTS_PASSED + 1)); else TESTS_FAILED=$((TESTS_FAILED + 1)); fi
    TESTS_RUN=$((TESTS_RUN + 1))
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
