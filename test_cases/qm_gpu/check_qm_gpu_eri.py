#!/usr/bin/env python3
"""Device-built ERI tensor (Q-G2) vs the CPU buildERI, element by element, plus the energy.

Claude Generated (Oct 2026). Runs curcuma with CURCUMA_QM_GPU_ERI_CHECK=1, which makes the
engine download the device tensor and compare it with qmint::buildERI (same basis, same
screening) and print one line "QM GPU ERI check: ... max|dev-cpu|=X ...". Passes when
  max|dev-cpu| <= --tol-eri (default 1e-12, absolute, FP64) and
  |E_gpu - E_cpu| <= --tol-e (default 1e-9 Eh).
Exit 77 (ctest SKIP) when no CUDA plugin / device is usable.
"""
import argparse, os, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--xyz", required=True)
ap.add_argument("--method", default="hf")
ap.add_argument("--basis", default="def2-SVP")
ap.add_argument("--tol-eri", type=float, default=1e-12)
ap.add_argument("--tol-e", type=float, default=1e-9)
ap.add_argument("--extra", default="", help="extra curcuma flags, e.g. -qm.cartesian_d true (one string, use --extra=...)")
a = ap.parse_args()


def run(gpu, check):
    cmd = [a.curcuma, "-sp", a.xyz, "-method", a.method, "-qm.basis", a.basis,
           "-gpu", gpu, "-qm.threads", "2", "-verbosity", "2", "-no_bmt"] + a.extra.split()
    env = dict(os.environ)
    if check:
        env["CURCUMA_QM_GPU_ERI_CHECK"] = "1"
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=1200, env=env)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)
    m = re.search(r"Single Point Energy = (-?\d+\.\d+)", txt)
    if not m:
        print(txt[-2000:])
        sys.exit("no energy found for gpu=" + gpu)
    return float(m.group(1)), txt


e_gpu, txt = run("cuda", True)
if "J/K on cuda" not in txt:
    print("no usable CUDA device/plugin: skipped")
    sys.exit(77)
if "built 4-centre ERI on cuda" not in txt:
    print(txt[-2000:])
    sys.exit("the device ERI build (Q-G2) did not run (fallback?)")
m = re.search(r"QM GPU ERI check: n=(\d+) elements=(\d+) max\|dev-cpu\|=(\S+) max\|cpu\|=(\S+) "
              r"above_1e-12=(\d+)", txt)
if not m:
    print(txt[-2000:])
    sys.exit("no ERI check line found")
n, nel, dmax, cmax, nbad = int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4)), int(m.group(5))
e_cpu, _ = run("none", False)
de = abs(e_gpu - e_cpu)
ok = dmax <= a.tol_eri and de <= a.tol_e
print(f"n={n} elements={nel} max|ERI_gpu-ERI_cpu|={dmax:.3e} (max|ERI|={cmax:.3e}, >1e-12: {nbad})  "
      f"E_cpu={e_cpu:.8f} E_gpu={e_gpu:.8f} diff={de:.2e}  {'ok' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
