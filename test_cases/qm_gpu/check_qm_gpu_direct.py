#!/usr/bin/env python3
"""Integral-direct J/K on the device (Q-G3) vs the CPU, J/K matrices and energy.

Claude Generated (Oct 2026). Runs curcuma with -gpu cuda and CURCUMA_QM_GPU_JK_CHECK=1: the
first device integral-direct build (a full density) is compared with the CPU DirectJK (same
screen; also the number of computed shell quartets), with the stored device tensor (Q-G2) when
it fits and with the stored CPU tensor when the host may hold it. Then the SCF energy
(E_elec, 10 decimals) is compared with a CPU run. Passes when
  max|J_dev - J_cpu_direct|, max|K_dev - K_cpu_direct| <= --tol-jk (default 1e-10, absolute),
  the same against the stored device/CPU tensor when available,
  the quartet counts agree (device == CPU DirectJK),
  |E_gpu - E_cpu| <= --tol-e (default 1e-9 Eh),
  and the run reports no GPU fallback.
--extra selects the path, e.g. "-qm.gpu_direct on" (forced) or
"-qm.gpu_max_memory_mb 100" (the stored tensor does not fit: auto must go direct on the device).
Exit 77 (ctest SKIP) when no CUDA plugin / device is usable.
"""
import argparse, os, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--xyz", required=True)
ap.add_argument("--method", default="hf")
ap.add_argument("--basis", default="def2-SVP")
ap.add_argument("--tol-jk", type=float, default=1e-10)
ap.add_argument("--tol-e", type=float, default=1e-9)
ap.add_argument("--extra", default="", help="extra curcuma flags (one string, use --extra=...)")
ap.add_argument("--cpu-extra", default="", help="extra flags for the CPU reference run only")
a = ap.parse_args()


def run(gpu, check, extra):
    cmd = [a.curcuma, "-sp", a.xyz, "-method", a.method, "-qm.basis", a.basis,
           "-gpu", gpu, "-qm.threads", "4", "-verbosity", "2", "-no_bmt"] + extra.split()
    env = dict(os.environ)
    if check:
        env["CURCUMA_QM_GPU_JK_CHECK"] = "1"
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000, env=env)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)
    m = re.search(r"E_elec\s*[:=]?\s*(-?\d+\.\d+) Eh", txt)
    if not m:
        print(txt[-3000:])
        sys.exit("no E_elec found for gpu=" + gpu)
    return float(m.group(1)), txt


e_gpu, txt = run("cuda", True, a.extra)
if "no GPU plugin" in txt or "not loadable" in txt:
    print("no usable CUDA device/plugin: skipped")
    sys.exit(77)
if "integral-direct J/K on cuda" not in txt:
    print(txt[-3000:])
    sys.exit("the device integral-direct path (Q-G3) did not run")
if "GPU fallbacks during this run" in txt:
    print(txt[-3000:])
    sys.exit("the GPU run reported a fallback")
m = re.search(r"QM GPU JK check: (.*)", txt)
if not m:
    print(txt[-3000:])
    sys.exit("no JK check line found")
line = m.group(1)
print(line)
kv = dict(re.findall(r"(\w+)=(\S+)", line))
ok = True
if kv["quartets_dev"] != kv["quartets_cpu"]:
    print(f"quartet count differs: device {kv['quartets_dev']} cpu {kv['quartets_cpu']}")
    ok = False
for key in ("dJ_cpu_direct", "dK_cpu_direct", "dJ_dev_stored", "dK_dev_stored", "dJ_cpu_stored", "dK_cpu_stored"):
    v = kv.get(key)
    if v is None or v.startswith("n/a"):
        continue
    if v == "failed" or not (float(v) <= a.tol_jk):
        print(f"{key}={v} above {a.tol_jk}")
        ok = False
e_cpu, _ = run("none", False, a.cpu_extra)
de = abs(e_gpu - e_cpu)
ok = ok and de <= a.tol_e
print(f"E_elec cpu={e_cpu:.10f} gpu-direct={e_gpu:.10f} diff={de:.2e}  {'ok' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
