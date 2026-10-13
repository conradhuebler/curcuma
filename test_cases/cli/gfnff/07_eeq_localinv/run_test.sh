#!/bin/bash

# Test: local-inverse preconditioner of the projected PCG (-eeq_solver.eeq_ppcg_precond localinv)
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026) - AI-generated, machine-tested; human production testing pending
#
#   1. SAME RESULT  - single point of peo201-water500 (2910 atoms, 501 fragments): localinv gives the diag energy to
#                     1e-8 Eh, and the Phase-2 solve needs at most half the iterations (CURCUMA_PPCG_ITERS lines)
#   2. ONE FRAGMENT - the same on peo201-chain (1410 atoms, nfrag = 1)
#   3. FALLBACK     - CURCUMA_LOCALINV_TEST_REJECT=1 rejects every build: warning, every solve on diag, energy
#                     bit-identical to the diag run (never a wrong preconditioner)
#   4. REPRODUCIBLE - two MD runs (2 fs) with localinv at 4 threads give bit-identical charges in every step
#      MATVEC       - symmetric threaded dense matvec of the ppcg (upper block triangle, fixed tiles, slices summed in
#                     fixed order, no atomics): same SP
#                     energy and same iteration counts at -threads 1 and 4, for diag and localinv
#   5. REBUILD      - eeq_localinv_rebuild_bohr 100 builds once in a 3 fs MD, 0 builds every Phase-2 solve
#                     (default 1.0 Bohr)
#   6. TYPO         - an unknown preconditioner name is an error
#   7. GPU          - only with a CUDA plugin and a device: localinv runs in WP7-E (lines "precond=localinv gpu"),
#                     G built on the device ("build=device"), fewer iterations than diag, charges within 1e-8 e of
#                     the GPU diag run in every MD step; the host build (CURCUMA_LOCALINV_HOST_BUILD=1) gives the
#                     same charges; an unusable G (CURCUMA_LOCALINV_TEST_INDEFINITE=1: -G; =3: true-residual check
#                     forced to fail) is caught, the solve redone with Jacobi (charges of diag), and -gpu_strict
#                     true then exits with code 3

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../test_utils.sh"

TEST_NAME="gfnff - 07: EEQ local-inverse preconditioner"

main() {
    test_header "$TEST_NAME"
    cd "$SCRIPT_DIR"
    set +e
    CURCUMA="$CURCUMA" python3 - <<'PYEOF'
import json, os, re, shutil, subprocess, sys, tempfile
B = os.environ["CURCUMA"]
ANSI = re.compile(r"\x1b\[[0-9;]*m")
fails = 0

def check(name, ok, detail=""):
    global fails
    print(("PASS  " if ok else "FAIL  ") + name + ("   " + detail if detail else ""))
    if not ok:
        fails += 1

def run(mode, xyz, extra, env_add=None, threads=4, fs=2):
    """Fresh directory (no topology cache). Returns rc, output, PPCG lines, energy, per-step charges."""
    d = tempfile.mkdtemp()
    try:
        shutil.copy(xyz, os.path.join(d, "m.xyz"))
        env = dict(os.environ, CURCUMA_PPCG_ITERS="1", OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
        env.update(env_add or {})
        cmd = [B, "-" + mode, "m.xyz", "-method", "gfnff", "-no_bmt", "-threads", str(threads)]
        if mode == "md":
            cmd += ["-max_time", str(fs), "-md_diagnostics", "true", "-dump_frequency", "1", "-seed", "1"]  # explicit seed: the CSVR noise comes from the per-instance RNG seeded by -seed (default -1 = time), so two runs are only comparable with a fixed seed
        p = subprocess.run(cmd + extra, cwd=d, capture_output=True, text=True, timeout=600, env=env)
        out = ANSI.sub("", p.stdout + p.stderr)
        pp = [l for l in out.splitlines() if l.startswith("PPCG")]
        e = re.search(r"Single Point Energy\s*=\s*(-?\d+\.\d+)", out)
        q = []
        fn = os.path.join(d, "m.diag.jsonl")
        if os.path.exists(fn):
            q = [json.loads(l)["charges"] for l in open(fn) if l.strip()]
        return p.returncode, out, pp, (float(e.group(1)) if e else None), q
    finally:
        shutil.rmtree(d, ignore_errors=True)

def iters(lines, tag):
    return [int(re.search(r"iters=(\d+)", l).group(1)) for l in lines if tag in l]

LI = ["-eeq_solver.eeq_ppcg_precond", "localinv"]
DIAG = ["-eeq_solver.eeq_ppcg_precond", "diag"]
W, CH = "water500.xyz", "chain.xyz"

for name, xyz in (("501 fragments", W), ("one fragment", CH)):
    rc0, _, pp0, e0, _ = run("sp", xyz, DIAG)
    rc1, _, pp1, e1, _ = run("sp", xyz, LI)
    it0 = max(iters(pp0, "iters="), default=0)              # Phase 2 is the larger count with diag
    it1 = iters(pp1, "precond=localinv")
    ok = rc0 == 0 and rc1 == 0 and e0 is not None and e1 is not None and abs(e0 - e1) < 1e-8 and it1 \
        and max(it1) * 2 <= it0
    check(f"SP {name}: same energy, at most half the iterations", ok,
          f"E diag {e0} localinv {e1}, Phase-2 iterations diag {it0} localinv {it1}")

# Default since Oct 10, 2026: localinv is used without any flag; "diag" restores the diagonal Jacobi.
rcd0, _, ppdef, edef, _ = run("sp", W, [])
rcd1, _, ppdg, edg, _ = run("sp", W, DIAG)
check("default preconditioner is localinv; diag selectable, same energy",
      rcd0 == 0 and rcd1 == 0 and iters(ppdef, "precond=localinv") and not iters(ppdg, "precond=localinv")
      and edef is not None and edg is not None and abs(edef - edg) < 1e-8,
      f"E default {edef} diag {edg}")

rc2, out2, pp2, e2, _ = run("sp", W, LI, {"CURCUMA_LOCALINV_TEST_REJECT": "1"})
rc0, _, _, e0, _ = run("sp", W, DIAG)
check("rejected build falls back to diag, energy bit-identical", rc2 == 0 and e2 == e0 and e0 is not None
      and "preconditioner rejected" in out2 and not iters(pp2, "precond=localinv"),
      f"E {e2} vs {e0}")

rca, _, _, _, qa = run("md", W, LI)
rcb, _, _, _, qb = run("md", W, LI)
check("localinv MD reproducible at a fixed thread count", rca == 0 and rcb == 0 and len(qa) >= 2 and qa == qb,
      f"steps {len(qa)}/{len(qb)}")

for pname, pextra in (("diag", DIAG), ("localinv", LI)):
    # The matrix itself differs by ~1e-13 between thread counts (threaded assembly), so only energy and
    # iteration counts are compared here; the bitwise independence of the symmetric matvec was checked separately.
    s1 = run("sp", W, pextra, threads=1)
    s4 = run("sp", W, pextra, threads=4)
    check(f"threaded matvec, SP {pname}: energy and ppcg iterations identical at 1 and 4 threads",
          s1[0] == 0 and s4[0] == 0 and s1[3] is not None and s1[3] == s4[3]
          and iters(s1[2], "iters=") == iters(s4[2], "iters=") and iters(s1[2], "iters="),
          f"E {s1[3]} vs {s4[3]}, iterations {iters(s1[2], 'iters=')} vs {iters(s4[2], 'iters=')}")

rc, _, pp, _, _ = run("md", W, LI + ["-eeq_solver.eeq_localinv_rebuild_bohr", "100"], fs=3)
nb = sum(1 for l in pp if "precond=localinv" in l and "build_ms=0.0 " not in l)
nsolve = sum(1 for l in pp if "precond=localinv" in l)
check("rebuild threshold 100 Bohr: one build", rc == 0 and nb == 1 and nsolve >= 3, f"builds {nb}, localinv solves {nsolve}")
rc, _, pp, _, _ = run("md", W, LI + ["-eeq_solver.eeq_localinv_rebuild_bohr", "0"], fs=3)
nb = sum(1 for l in pp if "precond=localinv" in l and "build_ms=0.0 " not in l)
nsolve = sum(1 for l in pp if "precond=localinv" in l)
check("rebuild threshold 0: a build per solve", rc == 0 and nb == nsolve and nsolve >= 3, f"builds {nb}, solves {nsolve}")

rc, out, _, _, _ = run("sp", CH, ["-eeq_solver.eeq_ppcg_precond", "lokalinv"])
check("unknown preconditioner is an error", rc != 0 and "eeq_ppcg_precond 'lokalinv' unknown" in out, f"rc={rc}")

plugin = os.path.join(os.path.dirname(os.path.abspath(B)), "libcurcuma_cuda.so")
have_gpu = os.path.exists(plugin) and shutil.which("nvidia-smi") is not None \
    and subprocess.run(["nvidia-smi", "-L"], capture_output=True).returncode == 0
if have_gpu:
    G = ["-gpu", "cuda", "-gpu_strict", "true"]
    rcd, _, ppd, _, qd = run("md", W, G + DIAG, fs=3)
    rcl, outl, ppl, _, ql = run("md", W, G + LI, fs=3)
    itd = [i for i in iters(ppd, "iters=")]
    itl = iters(ppl, "precond=localinv gpu")
    n = min(len(qd), len(ql))
    dq = max((abs(a - b) for i in range(n) for a, b in zip(qd[i], ql[i])), default=1.0)
    check("GPU MD: localinv runs in WP7-E with fewer iterations", rcd == 0 and rcl == 0 and len(itl) >= 2
          and max(itl) < min(itd[-len(itl):]), f"diag {itd} localinv {itl}")
    check("GPU MD: the local inverse is built on the device", len(iters(ppl, "build=device")) == len(itl)
          and "LOCALINV GPU device build ok=1" in outl, f"device-built solves {len(iters(ppl, 'build=device'))} of {len(itl)}")
    check("GPU MD: charges agree with GPU diag", n >= 2 and dq < 1e-8, f"steps {n}, max |dq| = {dq:.2e} e")

    def dq_to_diag(q):
        m = min(len(qd), len(q))
        return m, max((abs(a - b) for i in range(m) for a, b in zip(qd[i], q[i])), default=1.0)
    # Host build (comparison switch): same charges, solves tagged build=host.
    rch, outh, pph, _, qh = run("md", W, G + LI, {"CURCUMA_LOCALINV_HOST_BUILD": "1"}, fs=3)
    m, dqh = dq_to_diag(qh)
    check("GPU MD: CURCUMA_LOCALINV_HOST_BUILD=1 builds on the host, same charges", rch == 0
          and iters(pph, "build=host") and not iters(pph, "build=device") and m >= 2 and dqh < 1e-8,
          f"rc={rch}, max |dq| = {dqh:.2e} e")
    # Unusable preconditioner on the device: G -> -G (test hook). Every solve breaks down (rz <= 0)
    # and is redone with Jacobi: same charges as diag; each re-solve is a counted GPU fallback, so
    # -gpu_strict true stops with exit code 3.
    for mode, tag, what in (("1", "solve repeated with diagonal Jacobi", "negative definite G"),
                            ("3", "true residual above 2 tol", "true-residual check forced to fail")):
        rci, outi, _, _, qi = run("md", W, ["-gpu", "cuda"] + LI, {"CURCUMA_LOCALINV_TEST_INDEFINITE": mode}, fs=3)
        m, dqi = dq_to_diag(qi)
        check(f"GPU MD: {what}: Jacobi re-solve, charges of diag", rci == 0 and tag in outi
              and "localinv re-solve with Jacobi" in outi and m >= 2 and dqi < 1e-8,
              f"rc={rci}, steps {m}, max |dq| = {dqi:.2e} e")
        rcs, outs, _, _, _ = run("md", W, G + LI, {"CURCUMA_LOCALINV_TEST_INDEFINITE": mode}, fs=3)
        check(f"GPU MD: {what} with -gpu_strict true: exit code 3", rcs == 3 and "-gpu_strict" in outs, f"rc={rcs}")
    rcr, outr, ppr, _, qr = run("md", W, G + LI, {"CURCUMA_LOCALINV_TEST_REJECT": "1"}, fs=3)
    check("GPU MD: rejected build runs diag, no GPU fallback", rcr == 0 and "preconditioner rejected" in outr
          and not iters(ppr, "precond=localinv gpu") and len(qr) == len(qd)
          and max((abs(a - b) for i in range(len(qr)) for a, b in zip(qd[i], qr[i])), default=1.0) < 1e-8,
          f"rc={rcr}")
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
