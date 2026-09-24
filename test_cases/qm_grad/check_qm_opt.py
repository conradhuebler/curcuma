#!/usr/bin/env python3
"""-opt with the native hf-3c must reach the independent reference minimum.

Runs `curcuma -opt <start> -method hf-3c` (tight convergence), re-evaluates the
optimised structure with `-sp`, and compares energy and (Kabsch-aligned) geometry
with the reference minimum made by scripts/hf3c_opt_reference.py (PySCF +
simple-dftd3, BFGS to |g| < 1e-7). Claude Generated (Sep 2026).

Usage: check_qm_opt.py --curcuma <bin> --start <xyz> --ref <xyz>
"""
import argparse, os, re, shutil, subprocess, sys, tempfile


def read_xyz(path):
    L = open(path).read().splitlines(); n = int(L[0])
    return [l.split()[0] for l in L[2:2 + n]], [list(map(float, l.split()[1:4])) for l in L[2:2 + n]]


def kabsch_rmsd(a, b):
    import numpy as np
    a = np.array(a); b = np.array(b); a -= a.mean(0); b -= b.mean(0)
    U, S, Vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return float(np.sqrt(((a @ R - b) ** 2).sum(1).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--start", required=True)
    ap.add_argument("--ref", required=True)
    a = ap.parse_args()
    cur = os.path.abspath(a.curcuma)
    e_ref = float(re.search(r"E=(-?\d+\.\d+)", open(a.ref).read().splitlines()[1]).group(1))
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(a.start, os.path.join(tmp, "start.xyz"))
        subprocess.run([cur, "-opt", "start.xyz", "-method", "hf-3c", "-no_bmt", "-verbosity", "0",
                        "-qm.scf_threshold", "1e-9", "-opt.convergence_preset", "tight"],
                       cwd=tmp, capture_output=True, text=True, timeout=1500, check=True)
        out = subprocess.run([cur, "-sp", "start.opt.xyz", "-method", "hf-3c", "-no_bmt", "-verbosity", "1",
                              "-qm.scf_threshold", "1e-9"], cwd=tmp, capture_output=True, text=True,
                             timeout=600, check=True).stdout
        e = float(re.search(r"E\(HF-3c\)\s*=\s*(-?\d+\.\d+)", out).group(1))
        _, x_opt = read_xyz(os.path.join(tmp, "start.opt.xyz"))
    _, x_ref = read_xyz(a.ref)
    rmsd = kabsch_rmsd(x_opt, x_ref)
    ok = abs(e - e_ref) <= 1e-9 and rmsd <= 1e-4   # observed: 3e-12 Eh, 1e-6 A
    print(f"E(opt) {e:.10f}  E(ref) {e_ref:.10f}  diff {e - e_ref:+.2e}  RMSD {rmsd:.2e} A  {'ok' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
