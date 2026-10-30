#!/usr/bin/env python3
"""Kohn-Sham DFT on the GPU (cuda/qm_gpu_xc.cu, device SCF loop, device KS gradient) vs the CPU of the same binary.

Claude Generated (Oct 2026, GPU KS-DFT G1-G4). Modes:
  xc    : -gpu cuda -qm.gpu_scf off with CURCUMA_QM_GPU_XC_CHECK=1; every V_xc build also runs on the CPU from the SAME
          density and prints "QM GPU XC check: ...". Passes when max rel |dV| <= --tol-v and rel |dE_xc| <= --tol-exc for
          every build, the electron count and the number of points above the threshold agree.
  scf   : -gpu cuda (device run, --dev-extra) vs -gpu none (--cpu-extra): |E(dev) - E(cpu)| <= --tol-e and the SCF
          iteration counts differ by at most --tol-it; --loop device requires the device-resident loop
          ("device-resident loop" in the output), --loop host forbids it.
  grad  : -dump_gradient of a device run (CURCUMA_QM_GPU_XC_CHECK=1, CURCUMA_QM_GPU_GRAD_CHECK=1) vs a CPU run:
          XC gradient parts vs the CPU from the same density (<= --tol-xcg Eh/Bohr), device 2e part vs the CPU 2e part
          (<= --tol-2e), total gradient of the two runs (<= --tol-total Eh/Bohr).
  opt   : -opt (tight) with -gpu cuda and with -gpu none from the same start; the two minima must agree in energy
          (<= --tol-opt-e Eh, read from the -sp of each optimised structure on the CPU) and geometry (largest change
          of an interatomic distance <= --tol-rmsd Angstrom); the device run must have used the device XC gradient (no fallback).
  fallback : the device run must report a counted GPU fallback "qm_xc" and still give the CPU energy (--tol-e);
          with --strict, -gpu_strict true must end the run with exit code 3.
Common: --extra for both runs. A run reporting any GPU fallback fails (except mode fallback).
Exit 77 (ctest SKIP) when no CUDA plugin / device is usable.
"""
import argparse, os, re, subprocess, sys, tempfile

BOHR_IN_ANGSTROM = 0.529177210903  # CODATA 2018, as src/core/units.h

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--xyz", required=True)
ap.add_argument("--method", default="pbe")
ap.add_argument("--basis", default="def2-SVP")
ap.add_argument("--mode", choices=["xc", "scf", "grad", "opt", "fallback"], required=True)
ap.add_argument("--loop", choices=["device", "host", "any"], default="any")
ap.add_argument("--extra", default="")
ap.add_argument("--dev-extra", default="")
ap.add_argument("--cpu-extra", default="")
ap.add_argument("--tol-v", type=float, default=1e-11)
ap.add_argument("--tol-exc", type=float, default=1e-11)
ap.add_argument("--tol-e", type=float, default=1e-9)
ap.add_argument("--tol-it", type=int, default=1)
ap.add_argument("--tol-xcg", type=float, default=1e-9)
ap.add_argument("--tol-2e", type=float, default=1e-9)
ap.add_argument("--tol-total", type=float, default=1e-8)
ap.add_argument("--strict", action="store_true")
ap.add_argument("--tol-opt-e", type=float, default=1e-8)
ap.add_argument("--tol-rmsd", type=float, default=1e-4)
a = ap.parse_args()
a.curcuma, a.xyz = os.path.abspath(a.curcuma), os.path.abspath(a.xyz)
work = tempfile.mkdtemp(prefix="qm_gpu_dft_")


def run(gpu, extra, env_extra=None, dump=False, threshold="1e-10"):
    cmd = [a.curcuma, "-sp", a.xyz, "-method", a.method, "-qm.basis", a.basis, "-gpu", gpu,
           "-qm.scf_threshold", threshold, "-qm.threads", "4", "-verbosity", "2", "-no_bmt"]
    out = os.path.join(work, f"g_{gpu}.txt")
    if dump:
        cmd += ["-dump_gradient", out]
    cmd += (a.extra + " " + extra).split()
    env = dict(os.environ)
    env.update(env_extra or {})
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3000, env=env, cwd=work)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    m = list(re.finditer(r"E_elec:\s*(-?\d+\.\d+) Eh", txt))
    it = re.search(r"SCF converged in (\d+) iterations", txt)
    g = []
    if dump and os.path.exists(out):
        for line in open(out):
            if not line.startswith("#") and line.strip():
                g += [float(x) * BOHR_IN_ANGSTROM for x in line.split()]  # Eh/Angstrom -> Eh/Bohr
        os.remove(out)
    return (float(m[-1].group(1)) if m else None), (int(it.group(1)) if it else None), g, txt, r.returncode


def skip_check(txt):
    if "not loadable" in txt or "no GPU plugin" in txt:
        print("no usable CUDA device/plugin: skipped")
        sys.exit(77)


def no_fallback(txt, what):
    if "GPU fallbacks during this run" in txt:
        print(txt[-3000:])
        sys.exit(f"the {what} run reported a GPU fallback")


ok = True
if a.mode == "xc":
    e, it, _, txt, _ = run("cuda", a.dev_extra + " -qm.gpu_scf off", {"CURCUMA_QM_GPU_XC_CHECK": "1"})
    skip_check(txt)
    no_fallback(txt, "device")
    rows = re.findall(r"QM GPU XC check: n=(\d+) points=(\d+) max\|dV\|=(\S+) max\|V\|=(\S+) rel_dV=(\S+) dExc=(\S+) "
                      r"rel_dExc=(\S+) dEx=(\S+) dEc=(\S+) dN=(\S+) used dev=(\d+) cpu=(\d+) dev_ms=(\S+) cpu_ms=(\S+)", txt)
    if not rows:
        print(txt[-3000:])
        sys.exit("no 'QM GPU XC check' line: V_xc did not run on the device")
    worst_v = max(float(r[4]) for r in rows)
    worst_e = max(float(r[6]) for r in rows)
    worst_n = max(abs(float(r[9])) for r in rows)
    same_used = all(r[10] == r[11] for r in rows)
    ok = worst_v <= a.tol_v and worst_e <= a.tol_exc and worst_n <= 1e-10 and same_used
    print(f"{len(rows)} V_xc builds, n={rows[0][0]}, {rows[0][1]} points: max rel|dV| {worst_v:.2e} (tol {a.tol_v:.0e}), "
          f"max rel|dExc| {worst_e:.2e} (tol {a.tol_exc:.0e}), max|dN| {worst_n:.2e}, points used equal: {same_used}; "
          f"last build dev {rows[-1][12]} ms / cpu {rows[-1][13]} ms  {'ok' if ok else 'FAIL'}")
elif a.mode == "scf":
    e_d, it_d, _, txt, _ = run("cuda", a.dev_extra)
    skip_check(txt)
    no_fallback(txt, "device")
    if e_d is None or it_d is None:
        print(txt[-3000:])
        sys.exit("device run: no converged energy")
    dev_loop = "device-resident loop" in txt
    if a.loop == "device" and not dev_loop:
        print(txt[-3000:])
        sys.exit("the device-resident SCF loop did not run")
    if a.loop == "host" and dev_loop:
        sys.exit("the device-resident SCF loop ran although the host loop was expected")
    if "V_xc/E_xc on cuda" not in txt:
        print(txt[-3000:])
        sys.exit("V_xc did not run on the device")
    e_c, it_c, _, txt_c, _ = run("none", a.cpu_extra)
    d = abs(e_d - e_c)
    ok = d <= a.tol_e and abs(it_d - it_c) <= a.tol_it
    print(f"E_elec device={e_d:.10f} cpu={e_c:.10f} diff {d:.2e} (tol {a.tol_e:.0e}); iterations {it_d}/{it_c}; "
          f"loop {'device' if dev_loop else 'host'}  {'ok' if ok else 'FAIL'}")
    for line in txt.splitlines():
        if "QM: device SCF:" in line or "QM: SCF loop" in line:
            print(line.strip())
elif a.mode == "grad":
    env = {"CURCUMA_QM_GPU_XC_CHECK": "1", "CURCUMA_QM_GPU_GRAD_CHECK": "1"}
    e_d, it_d, g_d, txt, _ = run("cuda", a.dev_extra, env, dump=True)
    skip_check(txt)
    no_fallback(txt, "device")
    m = re.search(r"QM GPU XC grad check: natoms=(\d+) max\|dbasis\|=(\S+) max\|dpoints\|=(\S+) max\|dweights\|=(\S+) "
                  r"max\|basis\|=(\S+) points_used dev=(\d+) cpu=(\d+) dev_ms=(\S+) cpu_ms=(\S+)", txt)
    m2 = re.search(r"QM GPU grad check: natoms=(\d+) max\|dev-cpu\|=(\S+)", txt)
    if not m or not m2:
        print(txt[-3000:])
        sys.exit("XC gradient or 2e gradient did not run on the device")
    dxc = max(float(m.group(2)), float(m.group(3)), float(m.group(4)))
    d2e = float(m2.group(2))
    e_c, it_c, g_c, _, _ = run("none", a.cpu_extra, dump=True)
    if not g_d or len(g_d) != len(g_c):
        sys.exit("no gradient dump")
    dt = max(abs(x - y) for x, y in zip(g_d, g_c))
    ok = dxc <= a.tol_xcg and d2e <= a.tol_2e and dt <= a.tol_total and m.group(6) == m.group(7)
    print(f"XC gradient dev vs cpu (same density): basis {float(m.group(2)):.2e}, points {float(m.group(3)):.2e}, "
          f"weights {float(m.group(4)):.2e} (tol {a.tol_xcg:.0e}; dev {m.group(8)} ms, cpu {m.group(9)} ms); "
          f"2e part {d2e:.2e} (tol {a.tol_2e:.0e}); total gradient dev run vs cpu run {dt:.2e} Eh/Bohr "
          f"(tol {a.tol_total:.0e})  {'ok' if ok else 'FAIL'}")
elif a.mode == "opt":
    import shutil, time

    def read_xyz(path):
        L = open(path).read().splitlines()
        n = int(L[0])
        return [list(map(float, l.split()[1:4])) for l in L[2:2 + n]]

    def dist_dev(x, y):
        """max |r_ij(x) - r_ij(y)| over all atom pairs (rotation/translation invariant, no numpy needed)."""
        import math
        m = 0.0
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                m = max(m, abs(math.dist(x[i], x[j]) - math.dist(y[i], y[j])))
        return m

    res = {}
    for gpu in ["cuda", "none"]:
        d = os.path.join(work, gpu)
        os.makedirs(d)
        shutil.copy(a.xyz, os.path.join(d, "start.xyz"))
        t0 = time.time()
        r = subprocess.run([a.curcuma, "-opt", "start.xyz", "-method", a.method, "-qm.basis", a.basis, "-gpu", gpu,
                            "-no_bmt", "-verbosity", "1", "-qm.threads", "4", "-qm.scf_threshold", "1e-9",
                            "-opt.convergence_preset", "tight"] + (a.extra + " " + (a.dev_extra if gpu == "cuda" else a.cpu_extra)).split(),
                           cwd=d, capture_output=True, text=True, timeout=3000,
                           env=dict(os.environ, CURCUMA_QM_SCF_TRACE="1"))
        wall = time.time() - t0
        txt = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
        if gpu == "cuda":
            skip_check(txt)
            no_fallback(txt, "device -opt")
            if "2e=device" not in txt:  # CURCUMA_QM_SCF_TRACE line of every gradient
                print(txt[-3000:])
                sys.exit("the -gpu cuda -opt did not run the gradient on the device")
        opt = os.path.join(d, "start.opt.xyz")
        if not os.path.exists(opt):
            print(txt[-3000:])
            sys.exit(f"no optimised structure for -gpu {gpu}")
        sp = subprocess.run([a.curcuma, "-sp", "start.opt.xyz", "-method", a.method, "-qm.basis", a.basis, "-gpu", "none",
                             "-no_bmt", "-verbosity", "2", "-qm.scf_threshold", "1e-10"] + a.extra.split(),
                            cwd=d, capture_output=True, text=True, timeout=3000)
        st = re.sub(r"\x1b\[[0-9;]*m", "", sp.stdout + sp.stderr)
        e = float(list(re.finditer(r"E_elec:\s*(-?\d+\.\d+) Eh", st))[-1].group(1))
        res[gpu] = (e, read_xyz(opt), wall)
    de = abs(res["cuda"][0] - res["none"][0])
    rmsd = dist_dev(res["cuda"][1], res["none"][1])
    ok = de <= a.tol_opt_e and rmsd <= a.tol_rmsd
    print(f"-opt minimum: E(gpu run) {res['cuda'][0]:.10f} E(cpu run) {res['none'][0]:.10f} diff {de:.2e} Eh "
          f"(tol {a.tol_opt_e:.0e}); max |dr_ij| {rmsd:.2e} A (tol {a.tol_rmsd:.0e}); wall gpu {res['cuda'][2]:.1f} s "
          f"cpu {res['none'][2]:.1f} s  {'ok' if ok else 'FAIL'}")
else:  # fallback
    e_d, it_d, _, txt, rc = run("cuda", a.dev_extra)
    skip_check(txt)
    if a.strict:
        _, _, _, txt_s, rc_s = run("cuda", a.dev_extra + " -gpu_strict true")
        ok = rc_s == 3 and "qm_xc" in txt_s
        print(f"-gpu_strict true: exit code {rc_s} (expected 3), qm_xc named: {'qm_xc' in txt_s}  {'ok' if ok else 'FAIL'}")
        sys.exit(0 if ok else 1)
    counted = "GPU fallbacks during this run" in txt and "qm_xc" in txt
    e_c, it_c, _, _, _ = run("none", a.cpu_extra)
    d = abs(e_d - e_c) if e_d is not None and e_c is not None else 1.0
    ok = counted and d <= a.tol_e
    print(f"qm_xc fallback counted: {counted}; E_elec device-run={e_d} cpu={e_c} diff {d:.2e}  {'ok' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
