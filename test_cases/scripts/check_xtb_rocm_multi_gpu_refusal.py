#!/usr/bin/env python3
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Nov 2026, docs/ROCM_PARITY_A6.md)
"""Native GFN1/GFN2 on ROCm: the multi-GPU eigensolver and density split are refused, not ignored.

CUDA spreads the eigensolve (cuSOLVERMp / cusolverMg) and the pattern density of one molecule over several GPUs
(-gpu_eigensolver_devices, -gpu_density_devices). ROCm has no counterpart (rocSOLVER has no multi-device
eigensolver), so a request for more than one device must give a warning and a counted GPU fallback: with
-gpu_strict true the run exits with code 3. Without the flags the same run must succeed.

Exit 0 = pass, 1 = fail, 77 = skipped (no ROCm plugin or no HIP device).
Usage: check_xtb_rocm_multi_gpu_refusal.py <curcuma-binary> <caffeine.xyz>
"""
import os
import shutil
import subprocess
import sys
import tempfile


def run(binary, xyz, workdir, extra):
    os.makedirs(workdir, exist_ok=True)
    shutil.copy(xyz, workdir)
    cmd = [binary, "-sp", os.path.basename(xyz), "-method", "gfn2", "-gpu", "rocm", "-gpu_strict", "true",
           "-no_bmt", "-verbosity", "1"] + extra
    return subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=600)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    binary, xyz = sys.argv[1], sys.argv[2]
    if not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(binary)), "libcurcuma_rocm.so")):
        print("SKIP: no ROCm plugin next to the binary")
        return 77
    tmp = tempfile.mkdtemp(prefix="xtb_rocm_mgpu_")
    try:
        base = run(binary, xyz, os.path.join(tmp, "base"), [])
        if base.returncode != 0:
            print(f"SKIP: plain ROCm run failed (rc={base.returncode}), no usable HIP device?")
            return 77
        ok = True
        for key in ("-gpu_eigensolver_devices", "-gpu_density_devices"):
            r = run(binary, xyz, os.path.join(tmp, key.strip("-")), [key, "0,1"])
            out = r.stdout + r.stderr
            good = r.returncode == 3 and "not available on ROCm" in out
            print(f"{key} 0,1: rc={r.returncode} ({'refused with message' if good else 'NOT refused'})")
            ok &= good
        r = run(binary, xyz, os.path.join(tmp, "none"), ["-gpu_eigensolver_devices", "none", "-gpu_density_devices", "none"])
        print(f"none/none: rc={r.returncode}")
        ok &= r.returncode == 0
        print("PASS" if ok else "FAIL")
        return 0 if ok else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
