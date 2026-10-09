#!/bin/bash

# Test: EEQ solver control - the CLI chooses the solver, auto switches at the thresholds, a typo is an error
# Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026) - AI-generated, machine-tested; human production testing pending
#
#   1. AUTO      - default on 1410 atoms: "ppcg (auto)"; on caffeine (24 atoms): "cholesky (auto)"
#   2. THRESHOLD - -eeq_solver.eeq_ppcg_min_atoms 5000 turns auto on 1410 atoms into "cholesky (auto)"
#   3. EXPLICIT  - -eeq_solver.solve_method cholesky|ppcg|pcg|lu|ldlt runs exactly that solver, at a size where
#                  auto would pick another one (no size threshold overrides a named solver); energies agree
#                  to 1e-7 Eh; -gfnff.solve is an equivalent alias, -eeq_solver.solve_method wins over it
#   4. TYPO      - an unknown name ends with an error naming the valid solvers (it used to mean auto)
#   5. GPU       - only if a CUDA plugin and a device are present: lu has no device implementation and is an
#                  error with -gpu cuda, not a silent substitute
# The verbosity-1 line "EEQ solver (final charges): <name>[ (auto)] [CPU|GPU ...]" is what is read.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 06: EEQ solver control"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import os, re, shutil, subprocess, sys, tempfile
B = os.environ["CURCUMA"]
ANSI = re.compile(r"\x1b\[[0-9;]*m")
fails = 0

def check(name, ok, detail=""):
    global fails
    print(("PASS  " if ok else "FAIL  ") + name + ("   " + detail if detail else ""))
    if not ok:
        fails += 1

def sp(xyz, extra):
    """Fresh directory (no topology cache); returns (rc, solver line, energy, output)."""
    d = tempfile.mkdtemp()
    try:
        shutil.copy(xyz, os.path.join(d, "m.xyz"))
        p = subprocess.run([B, "-sp", "m.xyz", "-method", "gfnff", "-no_bmt", "-verbosity", "1", "-threads", "2"] + extra,
                           cwd=d, capture_output=True, text=True, timeout=600)
        out = ANSI.sub("", p.stdout + p.stderr)
        m = re.findall(r"EEQ solver \(final charges\): ([^\n]*)", out)
        e = re.search(r"Single Point Energy\s*=\s*(-?\d+\.\d+)", out)
        return p.returncode, (m[-1] if m else ""), (float(e.group(1)) if e else None), out
    finally:
        shutil.rmtree(d, ignore_errors=True)

CH, CAF = "chain.xyz", "caffeine.xyz"
S = "-eeq_solver.solve_method"

rc, line, e_auto, _ = sp(CH, [])
check("auto, 1410 atoms -> ppcg", rc == 0 and line.startswith("ppcg (auto) [CPU]"), line)
rc, line, _, _ = sp(CAF, [])
check("auto, 24 atoms -> cholesky", rc == 0 and line.startswith("cholesky (auto) [CPU]"), line)
rc, line, _, _ = sp(CH, ["-eeq_solver.eeq_ppcg_min_atoms", "5000"])
check("threshold argument moves the switch", rc == 0 and line.startswith("cholesky (auto) [CPU]"), line)

energies = {}
for name in ("cholesky", "ppcg", "pcg", "lu", "ldlt"):
    rc, line, e, _ = sp(CH, [S, name])
    energies[name] = e
    check(f"explicit {name} runs {name}", rc == 0 and line.startswith(name + " [CPU]"), line)
ref = energies["cholesky"]
def dev(names):
    return max(abs(energies[n] - ref) if energies[n] is not None and ref is not None else 1.0 for n in names)
# pcg stops at a relative residual of 1e-6 for large systems (about 4e-7 Eh here); the others are exact to the print precision
check("direct solvers and ppcg agree", dev(("ppcg", "lu", "ldlt")) < 1e-7, f"max |dE| to cholesky = {dev(('ppcg', 'lu', 'ldlt')):.2e} Eh")
check("pcg agrees within its tolerance", dev(("pcg",)) < 1e-5, f"|dE| = {dev(('pcg',)):.2e} Eh")
rc, line, _, _ = sp(CAF, [S, "ppcg"])
check("explicit ppcg below the threshold", rc == 0 and line.startswith("ppcg [CPU]"), line)
rc, line, _, _ = sp(CH, ["-gfnff.solve", "lu"])
check("-gfnff.solve alias", rc == 0 and line.startswith("lu [CPU]"), line)
rc, line, _, _ = sp(CH, ["-gfnff.solve", "lu", S, "ppcg"])
check("-eeq_solver.solve_method wins over -gfnff.solve", rc == 0 and line.startswith("ppcg [CPU]"), line)

rc, line, _, out = sp(CAF, [S, "choleski"])
check("typo is an error", rc != 0 and "unknown EEQ solver 'choleski'" in out and "valid:" in out, f"rc={rc}")

plugin = os.path.join(os.path.dirname(os.path.abspath(B)), "libcurcuma_cuda.so")
have_gpu = os.path.exists(plugin) and shutil.which("nvidia-smi") is not None \
    and subprocess.run(["nvidia-smi", "-L"], capture_output=True).returncode == 0
if have_gpu:
    rc, _, _, out = sp(CAF, ["-gpu", "cuda", S, "lu"])
    check("GPU: lu is an error, not a substitute", rc != 0 and "no device implementation" in out, f"rc={rc}")
else:
    print("SKIP  GPU check (no CUDA plugin or device)")
sys.exit(1 if fails else 0)
PYEOF
    local rc=$?
    if [ $rc -eq 0 ]; then TESTS_PASSED=$((TESTS_PASSED + 1)); else TESTS_FAILED=$((TESTS_FAILED + 1)); fi
    TESTS_RUN=$((TESTS_RUN + 1))
    print_test_summary
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}

main "$@"
