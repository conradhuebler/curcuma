#!/usr/bin/env python3
"""-opt with the native hf-3c must reach the independent reference minimum.

Runs `curcuma -opt <start> -method hf-3c` (tight convergence), re-evaluates the
optimised structure with `-sp`, and compares energy and (aligned with `curcuma -rmsd`) geometry
with the reference minimum made by scripts/hf3c_opt_reference.py (PySCF +
simple-dftd3, BFGS to |g| < 1e-7). Claude Generated (Sep 2026).

Usage: check_qm_opt.py --curcuma <bin> --start <xyz> --ref <xyz>
"""
import argparse, json, os, re, shutil, subprocess, sys, tempfile


def rmsd_curcuma(cur, tmp, a, b):
    """Optimal-rotation RMSD of two xyz files from curcuma's own -rmsd (no reordering, no numpy).

    Reads the full-precision "rmsd" from the <b>.rmsd.json that -rmsd writes next to the inputs.
    """
    subprocess.run([cur, "-rmsd", a, b, "-no_reorder", "true", "-no_bmt", "-verbosity", "0"],
                   cwd=tmp, capture_output=True, text=True, timeout=60, check=True)
    with open(os.path.join(tmp, os.path.splitext(b)[0] + ".rmsd.json")) as f:
        return float(json.load(f)["rmsd"])


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
        shutil.copy(a.ref, os.path.join(tmp, "ref.xyz"))
        rmsd = rmsd_curcuma(cur, tmp, "start.opt.xyz", "ref.xyz")
    ok = abs(e - e_ref) <= 1e-9 and rmsd <= 1e-4   # observed: 3e-12 Eh, 1e-6 A
    print(f"E(opt) {e:.10f}  E(ref) {e_ref:.10f}  diff {e - e_ref:+.2e}  RMSD {rmsd:.2e} A  {'ok' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
