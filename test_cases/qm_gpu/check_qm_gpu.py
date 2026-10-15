#!/usr/bin/env python3
"""CPU vs GPU single point of the native QM engine (Q-G1: J/K on the device).

Claude Generated (Oct 2026). Exit 77 (ctest SKIP) when no GPU plugin / device is usable.
The energy must agree with the CPU path to --tol Eh; the GPU run must not report a fallback.
"""
import argparse, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--xyz", required=True)
ap.add_argument("--method", default="hf")
ap.add_argument("--basis", default="def2-SVP")
ap.add_argument("--tol", type=float, default=1e-9)
a = ap.parse_args()


def run(gpu):
    cmd = [a.curcuma, "-sp", a.xyz, "-method", a.method, "-qm.basis", a.basis,
           "-gpu", gpu, "-qm.threads", "2", "-verbosity", "2", "-no_bmt"]
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)
    m = re.search(r"Single Point Energy = (-?\d+\.\d+)", txt)
    if not m:
        print(txt[-2000:]); sys.exit("no energy found for gpu=" + gpu)
    return float(m.group(1)), txt


e_cpu, _ = run("none")
e_gpu, txt = run("cuda")
if "J/K on cuda" not in txt:
    print("no usable CUDA device/plugin: skipped"); sys.exit(77)
d = abs(e_cpu - e_gpu)
print(f"CPU {e_cpu:.8f}  GPU {e_gpu:.8f}  diff {d:.2e}")
sys.exit(0 if d <= a.tol else 1)
