#!/usr/bin/env python3
"""Device-resident SCF loop (Q-G4) vs the host SCF loop and the CPU.

Claude Generated (Oct 2026). Three runs of the same single point:
  device loop : -gpu cuda -qm.gpu_scf on   (must report "device-resident loop", no fallback)
  host loop   : -gpu cuda -qm.gpu_scf off  (same device J/K builder, host DIIS/eigensolve)
  CPU         : -gpu none                  (--cpu-extra appended)
Passes when E_elec (10 decimals) of the device loop agrees with both within --tol-e (default
1e-9 Eh) and the device loop needs the same number of SCF iterations as the host loop.
--extra selects the J/K path, e.g. "-qm.gpu_direct on" (device integral-direct, incremental
Fock builds). Exit 77 (ctest SKIP) when no CUDA plugin / device is usable.
"""
import argparse, re, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--xyz", required=True)
ap.add_argument("--method", default="hf")
ap.add_argument("--basis", default="def2-SVP")
ap.add_argument("--tol-e", type=float, default=1e-9)
ap.add_argument("--extra", default="", help="extra flags for both -gpu cuda runs (one string, use --extra=...)")
ap.add_argument("--cpu-extra", default="", help="extra flags for the CPU reference run only")
a = ap.parse_args()


def run(gpu, extra):
    cmd = [a.curcuma, "-sp", a.xyz, "-method", a.method, "-qm.basis", a.basis,
           "-gpu", gpu, "-qm.threads", "4", "-verbosity", "2", "-no_bmt"] + extra.split()
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)
    m = re.search(r"E_elec\s*[:=]?\s*(-?\d+\.\d+) Eh", txt)
    it = re.search(r"SCF converged in (\d+) iterations", txt)
    if not m or not it:
        print(txt[-3000:])
        sys.exit(f"no E_elec / iteration count found for gpu={gpu} {extra}")
    return float(m.group(1)), int(it.group(1)), txt


e_dev, it_dev, txt = run("cuda", a.extra + " -qm.gpu_scf on")
if "no GPU plugin" in txt or "not loadable" in txt:
    print("no usable CUDA device/plugin: skipped")
    sys.exit(77)
if "device-resident loop" not in txt:
    print(txt[-3000:])
    sys.exit("the device-resident SCF loop (Q-G4) did not run")
if "GPU fallbacks during this run" in txt:
    print(txt[-3000:])
    sys.exit("the device-loop run reported a fallback")
for line in txt.splitlines():
    if "QM: device SCF:" in line or "QM: SCF loop" in line:
        print(line.strip())
e_host, it_host, txt_h = run("cuda", a.extra + " -qm.gpu_scf off")
if "device-resident loop" in txt_h:
    sys.exit("-qm.gpu_scf off still ran the device loop")
e_cpu, it_cpu, _ = run("none", a.cpu_extra)
d_host, d_cpu = abs(e_dev - e_host), abs(e_dev - e_cpu)
ok = d_host <= a.tol_e and d_cpu <= a.tol_e and it_dev == it_host
print(f"E_elec device-loop={e_dev:.10f} host-loop={e_host:.10f} cpu={e_cpu:.10f}  "
      f"diff host {d_host:.2e} cpu {d_cpu:.2e}  iterations {it_dev}/{it_host}/{it_cpu}  {'ok' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
