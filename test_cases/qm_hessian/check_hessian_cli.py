#!/usr/bin/env python3
"""-hessian capability with -method hf: analytic (-hessian.analytic true) vs numerical frequencies.

Runs `curcuma -hessian <xyz> -method hf` twice (default seminumerical path and the analytic RHF Hessian) in a
temporary directory with -no_bmt, parses the printed vibrational frequencies (cm^-1) and compares the 3N-6 largest.
This checks the integration seam (Eh/Bohr^2 from the engine -> Eh/Angstrom^2 of the capability -> frequencies),
not the accuracy of either Hessian (ctest qm_hess_* does that). Tolerance: the numerical path uses
-hessian.finite_diff_step 0.0005 and a tight SCF; its truncation error is well below 0.5 cm^-1 for water.
Usage: check_hessian_cli.py --curcuma BIN --xyz FILE. Claude Generated (Nov 2026).
"""
import argparse, os, re, shutil, subprocess, sys, tempfile, time

TOL_CM = 0.5
NUM = re.compile(r"^-?\d+\.\d+")


def frequencies(text):
    lines = text.splitlines()
    out = []
    start = None
    for i, l in enumerate(lines):
        if "Vibrational Frequencies" in l:
            start = i + 1
    if start is None:
        return out
    for l in lines[start:]:
        toks = re.sub(r"\x1b\[[0-9;]*m", "", l).replace("[OK]", " ").split()
        if not toks:
            continue
        vals = []
        for t in toks:
            m = NUM.match(t)
            if not m:
                vals = None
                break
            vals.append((float(m.group(0)), "(i)" in t))
        if vals is None:
            break
        out += vals
    return out


def run(curcuma, xyz, extra, wd):
    cmd = [curcuma, "-hessian", os.path.basename(xyz), "-method", "hf", "-no_bmt", "-verbosity", "1",
           "-qm.scf_threshold", "1e-10", "-qm.scf_energy_threshold", "1e-12",
           "-qm.scf_gradient_threshold_tight", "1e-10"] + extra
    t0 = time.time()
    p = subprocess.run(cmd, cwd=wd, capture_output=True, text=True, timeout=1700)
    dt = time.time() - t0
    if p.returncode != 0:
        print(p.stdout[-3000:], p.stderr[-3000:])
        sys.exit(f"curcuma failed: {' '.join(cmd)}")
    return p.stdout + p.stderr, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--fallback", metavar="FLAGS",
                    help="refusal check, one string (--fallback=\"-qm.ri jk\"): -hessian.analytic true with these flags must warn and run the numerical path")
    a = ap.parse_args()
    wd = tempfile.mkdtemp(prefix="qm_hess_cli_")
    shutil.copy(a.xyz, wd)
    n = int(open(a.xyz).readline().split()[0])
    if a.fallback:
        out, _ = run(os.path.abspath(a.curcuma), a.xyz, ["-hessian.analytic", "true"] + a.fallback.split(), wd)
        shutil.rmtree(wd, ignore_errors=True)
        ok = "using the numerical Hessian" in out and "Analytic RHF Hessian" not in out
        ok &= len([f for f, im in frequencies(out) if f > 100.0]) == 3 * n - 6
        for l in out.splitlines():
            if "using the numerical Hessian" in l:
                print(re.sub(r"\x1b\[[0-9;]*m", "", l).strip())
        print("ok" if ok else "FAIL")
        sys.exit(0 if ok else 1)
    out_num, t_num = run(os.path.abspath(a.curcuma), a.xyz, ["-hessian.finite_diff_step", "0.0005"], wd)
    out_ana, t_ana = run(os.path.abspath(a.curcuma), a.xyz, ["-hessian.analytic", "true"], wd)
    shutil.rmtree(wd, ignore_errors=True)
    if "Analytic RHF Hessian" not in out_ana:
        print(out_ana[-3000:])
        sys.exit("analytic path did not run")
    fn = sorted(f for f, im in frequencies(out_num) if not im)[-(3 * n - 6):]
    fa = sorted(f for f, im in frequencies(out_ana) if not im)[-(3 * n - 6):]
    if len(fn) != 3 * n - 6 or len(fa) != 3 * n - 6:
        print(out_num[-2000:], out_ana[-2000:])
        sys.exit("could not parse the frequencies")
    d = max(abs(x - y) for x, y in zip(fn, fa))
    print("numerical  :", " ".join(f"{x:.2f}" for x in fn), f"({t_num:.1f} s)")
    print("analytic   :", " ".join(f"{x:.2f}" for x in fa), f"({t_ana:.1f} s)")
    print(f"max |diff| = {d:.3f} cm^-1 (tolerance {TOL_CM})")
    print("ok" if d <= TOL_CM else "FAIL")
    sys.exit(0 if d <= TOL_CM else 1)


if __name__ == "__main__":
    main()
