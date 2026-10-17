#!/usr/bin/env python3
"""Native canonical MP2 analytic gradient against PySCF (reference_grad.json) and against the finite
difference of curcuma's own MP2 energy.

Claude Generated (Oct 2026). Runs `dump_qm_gradient <xyz> --method <mp2|scs-mp2|sos-mp2> [--frozen_core] [--fd]`
(test_cases/qm_grad/dump_qm_gradient.cpp, Eh/Bohr) and checks
  (a) method mp2 only (PySCF has no SCS/SOS gradient), references from scripts/mp2_gradient_reference.py
      (PySCF 2.14, def2-SVP spherical, RHF conv_tol 1e-13 conv_tol_grad 1e-11):
      max |g - g_FD(PySCF)| <= --tol (default 1e-8 Eh/Bohr; four-point central difference of PySCF's MP2
      energy, h = 1e-3 Bohr), max |g - g_PySCF analytic| <= --tol-analytic (default 2.5e-7: PySCF's
      analytic MP2 gradient itself differs from its own FD by up to 1.1e-7 on H2O) and |E - E_PySCF| <= 1e-8 Eh;
  (b) with --fd: max |g - FD(own energy)| <= --tol-fd (central difference, h = 1e-4 Bohr);
  (c) |sum over atoms of g| <= 1e-9 per component (translational invariance).
--gpu (with --curcuma): runs `curcuma -sp -gradient -dump_gradient -gpu cuda` instead (Eh/Angstrom, converted)
with CURCUMA_QM_GPU_GRAD_CHECK=1; the 2e derivative part must have run on the device, no fallback may be
reported, and the printed device-vs-CPU difference of the 2e part must be <= --tol-dev (1e-9 Eh/Bohr); exit 77
(ctest SKIP) without a usable CUDA plugin/device.
Exit 0 = pass.
"""
import argparse, json, os, re, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--dump", required=True)
ap.add_argument("--ref", required=True)
ap.add_argument("--mol", required=True)
ap.add_argument("--method", default="mp2")
ap.add_argument("--fc", action="store_true")
ap.add_argument("--fd", action="store_true")
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--tol", type=float, default=1e-8)
ap.add_argument("--tol-analytic", type=float, default=2.5e-7)
ap.add_argument("--tol-fd", type=float, default=1e-7)
ap.add_argument("--extra", default="")
ap.add_argument("--gpu", action="store_true")
ap.add_argument("--curcuma", default="")
ap.add_argument("--tol-dev", type=float, default=1e-9)
a = ap.parse_args()
BOHR = 0.529177210903

R = json.load(open(a.ref))["cases"][a.mol]
xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", R["xyz"])
if a.gpu:
    with tempfile.TemporaryDirectory() as tmp:
        dump = os.path.join(tmp, "grad.txt")
        cmd = [os.path.abspath(a.curcuma), "-sp", xyz, "-method", a.method, "-gradient", "-dump_gradient", dump,
               "-qm.scf_threshold", "1e-10", "-qm.threads", str(a.threads), "-gpu", "cuda", "-no_bmt", "-verbosity", "2"]
        if a.fc:
            cmd += ["-qm.frozen_core", "true"]
        cmd += a.extra.split()
        env = dict(os.environ, CURCUMA_QM_GPU_GRAD_CHECK="1")
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=5000, cwd=tmp, env=env)
        txt = re.sub(r"\x1b\[[0-9;]*m", "", p.stdout + p.stderr)
        if "no GPU plugin" in txt or "not loadable" in txt:
            print("no usable CUDA device/plugin: skipped")
            sys.exit(77)
        if p.returncode != 0 or not os.path.exists(dump):
            print(" ".join(cmd))
            print(txt[-3000:])
            sys.exit("curcuma -sp -gradient -gpu cuda failed")
        L = open(dump).read().splitlines()
    m = re.search(r"QM GPU MP2 grad check: .*max\|dev-cpu\|=(\S+)", txt)
    for line in txt.splitlines():
        if "QM GPU MP2 grad check" in line or "MP2 gradient timing" in line:
            print(line.strip())
    if not m or "2e derivative integrals on the device" not in txt:
        sys.exit("the device MP2 2e gradient did not run")
    if "GPU fallbacks during this run" in txt:
        print(txt[-3000:])
        sys.exit("the -gpu run reported a fallback")
    ddev = float(m.group(1))
    print(f"device vs CPU (2e part): {ddev:.2e} {'ok' if ddev <= a.tol_dev else 'FAIL'}")
    got = {"energy": float(re.search(r"energy (-?\d+\.\d+)", L[1]).group(1)),
           "gradient": [float(x) * BOHR for l in L[2:] for x in l.split()]}
    got["natoms"] = len(got["gradient"]) // 3
    gpu_ok = ddev <= a.tol_dev
else:
    gpu_ok = True
    cmd = [os.path.abspath(a.dump), xyz, "--method", a.method, "--basis", "def2-SVP", "--threads", str(a.threads)]
    if a.fc:
        cmd.append("--frozen_core")
    if a.fd:
        cmd.append("--fd")
    cmd += a.extra.split()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=5000)
    if p.returncode != 0:
        print(" ".join(cmd))
        print(p.stdout[-2000:], p.stderr[-2000:])
        sys.exit("dump_qm_gradient failed")
    got = json.loads(p.stdout)
g = got["gradient"]
n = got["natoms"]
ok = gpu_ok
if a.method == "mp2":
    ref = R["fc" if a.fc else "full"]
    d = max(abs(x - y) for x, y in zip(g, ref["fd_gradient"]))
    da = max(abs(x - y) for x, y in zip(g, ref["gradient"]))
    de = abs(got["energy"] - ref["e_tot"])
    good = d <= a.tol and da <= a.tol_analytic and de <= 1e-8
    ok &= good
    print(f"{a.mol}{' fc' if a.fc else ''}: E diff {de:.2e}  max|g - g_FD(PySCF)| = {d:.2e}  "
          f"max|g - g_PySCF analytic| = {da:.2e} (max|g| {max(map(abs, g)):.3e}) {'ok' if good else 'FAIL'}")
if "fd_gradient" in got:
    d = max(abs(x - y) for x, y in zip(g, got["fd_gradient"]))
    good = d <= a.tol_fd
    ok &= good
    print(f"{a.mol}{' fc' if a.fc else ''} {a.method}: max|g - FD(own energy)| = {d:.2e} {'ok' if good else 'FAIL'}")
s = [sum(g[3 * i + k] for i in range(n)) for k in range(3)]
good = max(map(abs, s)) <= 1e-9
ok &= good
print(f"sum over atoms = ({s[0]:.2e}, {s[1]:.2e}, {s[2]:.2e}) {'ok' if good else 'FAIL'}")
sys.exit(0 if ok else 1)
