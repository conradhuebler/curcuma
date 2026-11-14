#!/usr/bin/env python3
"""Claude Generated (Oct 2026): guards of the CPU projected PCG (EEQ Phase 2) after two speed-ups.

Uses only the curcuma binary and test_cases/structures/clusters/urea-water-cluster197.unrelaxed.xyz (197 atoms,
many fragments), -eeq_solver.solve_method ppcg, default preconditioner localinv. Checks
  1. DEFERRED CONSTRAINT BLOCK - the PPCG path leaves the constraint rows/columns of A unwritten; a forced PPCG
     failure (eeq_ppcg_max_iter 2) fills them before the Schur-Cholesky fallback: warning, energy line and gradient
     dump bit-identical to -solve_method cholesky. A converged ppcg run agrees with cholesky to 1e-8 Eh.
  2. NO DEFINITENESS CHECK - a normal localinv run: no guard warning, check_ms=0.0 on every build line, the first
     solve of the build passed the true-residual check (true_res/tol <= 10).
  3. GUARDS - CURCUMA_LOCALINV_TEST_INDEFINITE: 1 (G -> -G) and 2 (G_ii -> -G_ii, every 10th atom) are rejected by the
     S = C G C^T factorisation, 3 by the forced true-residual check, 4 (G - s v v^T with C v = 0: S unchanged, G
     indefinite) only by the CG (rz <= 0 / pAp <= 0). Each: its warning, every solve on diag afterwards, energy and
     gradient dump bit-identical to the eeq_ppcg_precond diag run (never a wrong preconditioner).
  4. MD - 3 fs (6 steps of 0.5 fs) with mode 4 and eeq_localinv_rebuild_bohr 0 (a new, indefinite G every step): rc 0, no Schur
     fallback, per-step charges bit-identical to the diag MD.
Usage: check_eeq_ppcg_guards.py CURCUMA_BINARY SOURCE_DIR
"""
import hashlib, json, os, re, shutil, subprocess, sys, tempfile

curcuma, src = sys.argv[1], sys.argv[2]
xyz = os.path.join(src, "test_cases/structures/clusters/urea-water-cluster197.unrelaxed.xyz")
ANSI = re.compile(r"\x1b\[[0-9;]*m")
failed = 0


def check(ok, what):
    global failed
    print(("  PASS  " if ok else "  FAIL  ") + what)
    failed += 0 if ok else 1


def run(mode, extra, env_add=None):
    """Fresh directory per run (no topology cache). Returns rc, log, energy line, gradient md5, per-step charges."""
    d = tempfile.mkdtemp(prefix="curcuma_ppcg_guards_")
    try:
        shutil.copy(xyz, os.path.join(d, "m.xyz"))
        env = dict(os.environ, CURCUMA_PPCG_ITERS="1")
        env.update(env_add or {})
        cmd = [curcuma, "-" + mode, "m.xyz", "-method", "gfnff", "-no_bmt", "-threads", "4", "-verbosity", "1",
               "-gfnff.cache_topology", "false", "-eeq_solver.solve_method", "ppcg"]
        if mode == "sp":
            cmd += ["-dump_gradient", "grad.txt"]
        else:
            cmd += ["-md.max_time", "3", "-md.time_step", "0.5", "-md.md_diagnostics", "true", "-md.print_frequency", "1", "-md.dump_frequency", "1",
                    "-md.write_restart_frequency", "-1", "-md.seed", "7"]
        p = subprocess.run(cmd + extra, cwd=d, capture_output=True, text=True, timeout=600, env=env)
        log = ANSI.sub("", p.stdout + p.stderr)
        g = os.path.join(d, "grad.txt")
        md5 = hashlib.md5(open(g, "rb").read()).hexdigest() if os.path.exists(g) else None
        e = next((l for l in open(g) if l.startswith("# energy")), None) if md5 else None
        q = []
        fn = os.path.join(d, "m.diag.jsonl")
        if os.path.exists(fn):
            q = [json.loads(l)["charges"] for l in open(fn) if l.strip()]
        return p.returncode, log, e, md5, q
    finally:
        shutil.rmtree(d, ignore_errors=True)


def ppcg_lines(log):
    return [l for l in log.splitlines() if l.startswith("PPCG")]


def energy(e_line):
    return float(e_line.split()[2])


# 1. deferred constraint block
rc_c, log_c, e_c, g_c, _ = run("sp", ["-eeq_solver.solve_method", "cholesky"])
rc_f, log_f, e_f, g_f, _ = run("sp", ["-eeq_solver.eeq_ppcg_max_iter", "2"])
check(rc_c == 0 and rc_f == 0 and e_c is not None and e_f is not None, "cholesky and forced-failure ppcg runs finished")
check("did not converge, using the exact Schur-Cholesky" in log_f, "forced PPCG failure falls back to Schur-Cholesky")
check(e_f == e_c and g_f == g_c, "fallback after the deferred constraint block: energy and gradient bit-identical to cholesky")
rc_p, log_p, e_p, g_p, _ = run("sp", [])
check(rc_p == 0 and e_p is not None and abs(energy(e_p) - energy(e_c)) < 1e-8,
      f"converged ppcg energy = cholesky to 1e-8 Eh ({energy(e_p) - energy(e_c):+.2e})" if e_p else "converged ppcg run")

# 2. normal localinv run: no definiteness check, true-residual check passed
pp = ppcg_lines(log_p)
builds = [l for l in pp if re.search(r"build_ms=[1-9]|build_ms=0\.[1-9]", l)]
ratios = [float(m.group(1)) for l in pp for m in [re.search(r"true_res/tol=([-0-9.e+]+)", l)] if m and float(m.group(1)) >= 0]
check(len(builds) >= 1 and all("check_ms=0.0" in l for l in builds), f"{len(builds)} localinv builds, none with a definiteness check")
check(len(ratios) >= 1 and max(ratios) <= 10.0, f"true-residual check after the first solve passed (max ratio {max(ratios) if ratios else -1:.3g})")
check("unusable" not in log_p and "rejected" not in log_p and all("precond=localinv" in l for l in pp[-1:]),
      "no guard fired in the normal run, last solve on localinv")

# 3. guards against an unusable G
rc_d, log_d, e_d, g_d, _ = run("sp", ["-eeq_solver.eeq_ppcg_precond", "diag"])
expect = {"1": "rejected (S = C G C^T is not positive definite)", "2": "rejected (S = C G C^T is not positive definite)",
          "3": "true residual above the tolerance", "4": "indefinite in the CG"}
for mode, msg in expect.items():
    rc_m, log_m, e_m, g_m, _ = run("sp", [], {"CURCUMA_LOCALINV_TEST_INDEFINITE": mode})
    pm = ppcg_lines(log_m)
    check(rc_m == 0 and msg in log_m, f"mode {mode}: caught ({msg})")
    check(len(pm) > 0 and all("precond=diag" in l for l in pm), f"mode {mode}: every solve reported on diag")
    check(e_m == e_d and g_m == g_d, f"mode {mode}: energy and gradient bit-identical to the diag run")

# 4. MD with a new indefinite G every step
rc_md, log_md, _, _, q_d = run("md", ["-eeq_solver.eeq_ppcg_precond", "diag"])
rc_mi, log_mi, _, _, q_i = run("md", ["-eeq_solver.eeq_localinv_rebuild_bohr", "0"], {"CURCUMA_LOCALINV_TEST_INDEFINITE": "4"})
n_caught = log_mi.count("indefinite in the CG")
check(rc_md == 0 and rc_mi == 0 and len(q_d) >= 5, f"MD runs finished ({len(q_d)} diagnostic steps)")
check(n_caught >= len(q_d) and "Schur-Cholesky" not in log_mi, f"MD: indefinite G caught {n_caught} times, no Schur fallback")
check(q_i == q_d, "MD: per-step charges bit-identical to the diag MD")

print(f"{failed} check(s) failed" if failed else "all checks passed")
sys.exit(1 if failed else 0)
