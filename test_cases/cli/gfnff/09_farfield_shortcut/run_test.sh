#!/bin/bash

# Test: exact far-field shortcut of the Coulomb term (-gfnff.farfield_shortcut), CPU
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026) - AI-generated, machine-tested; human production testing pending
#
# For gamma_ij*r >= 6.5 (gamma_ij = 1/sqrt(alpha_i+alpha_j)) erf(gamma_ij*r) is exactly 1.0 in FP64, so the pair
# energy/force of FFWorkspace::calcCoulomb and the EEQ matrix elements are evaluated without erf/exp. The shortcut
# must not change a single bit: with -gfnff.farfield_shortcut true (default) and false
#   1. SP       - energy + gradient dump (-dump_gradient, 15 digits) are byte-identical for the 1410-atom
#                 peo201-chain, caffeine, and the charged hydronium-water2 / cl2-anion
#                 (-gfnff.frag_charge_model reference), at 2 threads (CPU is bit-reproducible at a fixed thread count);
#   2. FORCED   - two neutral waters 60 Bohr apart (every inter-molecular pair is far field, so the shortcut
#                 certainly fires): byte-identical too, and the energy equals the sum of the two monomers
#                 (+ the point-charge interaction of the monomers) within 1e-3 Eh (liveness of the geometry);
#   3. MD       - 3 fs of MD on caffeine and on the two-water system: identical Epot series (md diagnostics: per-term energies and EEQ charges of every step).
# Docs: docs/GFNFF_STATUS.md, section "CPU far-field shortcut".

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 09: exact far-field shortcut of the Coulomb term"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import os, re, subprocess, sys, tempfile, shutil, json, math
B = os.environ["CURCUMA"]
REF = ["-gfnff.frag_charge_model", "reference"]
fails = 0
ANSI = re.compile(r"\x1b\[[0-9;]*m")

def check(name, ok, msg=""):
    global fails
    print(("PASS " if ok else "FAIL ") + name + ("   " + msg if msg else ""), flush=True)
    if not ok: fails += 1

def work(xyz_text=None, xyz_file=None):
    d = tempfile.mkdtemp()
    if xyz_file: shutil.copy(xyz_file, os.path.join(d, "m.xyz"))
    else: open(os.path.join(d, "m.xyz"), "w").write(xyz_text)
    return d

def sp(src, charge, flag, extra):
    d = work(*src)
    try:
        cmd = [B, "-sp", "m.xyz", "-method", "gfnff", "-charge", str(charge), "-threads", "2", "-verbosity", "0",
               "-no_bmt", "-gfnff.cache_topology", "false", "-gfnff.farfield_shortcut", flag, "-dump_gradient", "g.txt"] + extra
        p = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=600)
        g = open(os.path.join(d, "g.txt")).read() if os.path.exists(os.path.join(d, "g.txt")) else None
        e = None
        for l in (g or "").splitlines():
            if l.startswith("# energy"): e = float(l.split()[2])
        return p.returncode, e, g
    finally:
        shutil.rmtree(d, ignore_errors=True)

def md(src, charge, flag, extra):
    d = work(*src)
    try:
        cmd = [B, "-md", "m.xyz", "-method", "gfnff", "-charge", str(charge), "-threads", "2", "-no_bmt",
               "-gfnff.cache_topology", "false", "-gfnff.farfield_shortcut", flag, "-max_time", "3", "-seed", "1",
               "-md_diagnostics", "true", "-dump_frequency", "1"] + extra
        p = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=600)
        ep = []
        fn = os.path.join(d, "m.diag.jsonl")
        if os.path.exists(fn):
            for l in open(fn):
                if l.strip():
                    j = json.loads(l)
                    ep.append((j["energy"], j["charges"], j["gradient_norm"]))  # per-term energies, EEQ charges, per-atom gradient norms
        return p.returncode, ep
    finally:
        shutil.rmtree(d, ignore_errors=True)

# two neutral waters, 60 Bohr (31.75 A) apart
w = [("O", 0.0, 0.0, 0.0), ("H", 0.7570, 0.5860, 0.0), ("H", -0.7570, 0.5860, 0.0)]
shift = 60.0 * 0.529177210903
two = w + [(a, x + shift, y, z) for a, x, y, z in w]
two_txt = "6\n\n" + "".join(f"{a} {x:.10f} {y:.10f} {z:.10f}\n" for a, x, y, z in two)
one_txt = "3\n\n" + "".join(f"{a} {x:.10f} {y:.10f} {z:.10f}\n" for a, x, y, z in w)

for name, src, q, ex in (("peo201-chain (1410 atoms)", (None, "chain.xyz"), 0, []),
                         ("caffeine", (None, "caffeine.xyz"), 0, []),
                         ("hydronium-water2 +1", (None, "h3op_h2o2.xyz"), 1, REF),
                         ("cl2-anion -1", (None, "cl2m.xyz"), -1, REF),
                         ("two waters 60 Bohr apart", (two_txt, None), 0, [])):
    r1, e1, g1 = sp(src, q, "true", ex)
    r0, e0, g0 = sp(src, q, "false", ex)
    check(f"SP {name}: energy + gradient dump byte-identical", r1 == 0 and r0 == 0 and g1 is not None and g1 == g0,
          f"E true {e1!r} false {e0!r}")

# liveness of the forced case: dimer energy ~ 2 x monomer energy (neutral, no overlap, charges polar -> dipole-dipole < 1e-3 Eh)
_, e2, _ = sp((two_txt, None), 0, "true", [])
_, e1m, _ = sp((one_txt, None), 0, "true", [])
check("two waters 60 Bohr apart: E(dimer) = 2 E(monomer) within 1e-3 Eh", e2 is not None and e1m is not None and abs(e2 - 2 * e1m) < 1e-3,
      f"E dimer {e2}, 2 x monomer {None if e1m is None else 2 * e1m}")

for name, src, q in (("caffeine", (None, "caffeine.xyz"), 0), ("two waters 60 Bohr apart", (two_txt, None), 0)):
    r1, ep1 = md(src, q, "true", [])
    r0, ep0 = md(src, q, "false", [])
    check(f"MD 3 fs {name}: Epot series identical", r1 == 0 and r0 == 0 and len(ep1) >= 3 and ep1 == ep0,
          f"{len(ep1)} / {len(ep0)} frames")
sys.exit(1 if fails else 0)
PYEOF
    local rc=$?
    if [ $rc -eq 0 ]; then TESTS_PASSED=$((TESTS_PASSED + 1)); else TESTS_FAILED=$((TESTS_FAILED + 1)); fi
    TESTS_RUN=$((TESTS_RUN + 1))
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
