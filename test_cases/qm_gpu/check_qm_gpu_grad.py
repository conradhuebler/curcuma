#!/usr/bin/env python3
"""Two-electron (Q-G7) and one-electron (Q-G8) gradient on the device vs the CPU.

Claude Generated (Oct 2026). Two `curcuma -sp ... -dump_gradient` runs of the same structure:
  device : -gpu cuda with CURCUMA_QM_GPU_GRAD_CHECK=1 (the engine also evaluates the CPU 2e
           gradient from the same density and prints "QM GPU grad check: ...")
  CPU    : -gpu none
Passes when
  * the device path ran (check line present, no counted GPU fallback),
  * max|g2e(device) - g2e(CPU)| <= --tol-2e (default 1e-9 Eh/Bohr, same density),
  * if the 1e part ran on the device ("QM GPU grad1e check: ..." line, required with --require-1e):
    max|g1e(device) - g1e(CPU)| <= --tol-1e (default 1e-10 Eh/Bohr, same P and W) and its
    translational invariance |sum| <= --tol-sum (--no-2e: the 2e part is expected on the CPU),
  * |sum over atoms of g2e(device)| <= --tol-sum per direction (translational invariance),
  * max|g(device run) - g(CPU run)| <= --tol-total (default 1e-8 Eh/Bohr; two independent SCFs,
    scf_threshold 1e-10),
  * with --fd: the device-run gradient agrees with the central finite difference of the device-run
    energy (Angstrom displacements --fd-step, energies from the 15-decimal dump header, both in
    Eh/Angstrom as the getGradient() contract) within --tol-fd.
--extra adds flags to both runs, --dev-extra to the device run only (e.g. "-qm.gpu_direct on").
Exit 77 (ctest SKIP) when no CUDA plugin / device is usable.
"""
import argparse, os, re, subprocess, sys, tempfile

BOHR_IN_ANGSTROM = 0.529177210903  # CODATA 2018, as src/core/units.h

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--xyz", required=True)
ap.add_argument("--method", default="hf")
ap.add_argument("--basis", default="def2-SVP")
ap.add_argument("--extra", default="", help="flags for both runs (one string, use --extra=...)")
ap.add_argument("--dev-extra", default="", help="flags for the device run only")
ap.add_argument("--tol-2e", type=float, default=1e-9)
ap.add_argument("--tol-1e", type=float, default=1e-10)
ap.add_argument("--require-1e", action="store_true", help="the 1e part must have run on the device (Q-G8)")
ap.add_argument("--no-2e", action="store_true", help="the 2e part is expected on the CPU (e.g. -qm.gpu_grad off)")
ap.add_argument("--tol-sum", type=float, default=1e-10)
ap.add_argument("--tol-total", type=float, default=1e-8)
ap.add_argument("--fd", action="store_true")
ap.add_argument("--fd-step", type=float, default=1e-4, help="Angstrom")
ap.add_argument("--tol-fd", type=float, default=1e-7, help="Eh/Angstrom (h = 1e-4: 2.1e-8 measured for H2O)")
a = ap.parse_args()
work = tempfile.mkdtemp(prefix="qm_gpu_grad_")


def run(gpu, extra, xyz, check):
    out = os.path.join(work, f"g_{gpu}.txt")
    cmd = [a.curcuma, "-sp", xyz, "-method", a.method, "-qm.basis", a.basis, "-gpu", gpu,
           "-qm.scf_threshold", "1e-10", "-qm.threads", "4", "-verbosity", "1", "-no_bmt",
           "-dump_gradient", out] + extra.split()
    env = dict(os.environ)
    env["CURCUMA_QM_GPU_GRAD_CHECK"] = "1" if check else "0"
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3000, env=env, cwd=work)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout + r.stderr)
    if not os.path.exists(out):
        print(txt[-3000:])
        sys.exit(f"no gradient dump for gpu={gpu}")
    e = None
    g = []
    for line in open(out):
        if line.startswith("# energy"):
            e = float(line.split()[2])
        elif not line.startswith("#") and line.strip():
            g += [float(x) for x in line.split()]  # Eh/Angstrom
    os.remove(out)
    return e, g, txt


e_dev, g_dev, txt = run("cuda", a.extra + " " + a.dev_extra, a.xyz, True)
if "not loadable" in txt or "no GPU plugin" in txt:
    print("no usable CUDA device/plugin: skipped")
    sys.exit(77)
if "GPU fallbacks during this run" in txt:
    print(txt[-3000:])
    sys.exit("the device run reported a GPU fallback")
m = re.search(r"QM GPU grad check: natoms=(\d+) max\|dev-cpu\|=(\S+) max\|cpu\|=(\S+) "
              r"sum_dev=\((\S+),(\S+),(\S+)\) sum_cpu=\((\S+),(\S+),(\S+)\) quartets=(\d+) dev_ms=(\S+) cpu_ms=(\S+)", txt)
d2e, sum_dev = 0.0, 0.0
if m:
    print(m.group(0))
    d2e = float(m.group(2))
    sum_dev = max(abs(float(m.group(k))) for k in (4, 5, 6))
elif not a.no_2e:
    print(txt[-3000:])
    sys.exit("no 'QM GPU grad check' line: the 2e gradient did not run on the device")
m1 = re.search(r"QM GPU grad1e check: natoms=(\d+) max\|dev-cpu\|=(\S+) max\|cpu\|=(\S+) "
               r"sum_dev=\((\S+),(\S+),(\S+)\) sum_cpu=\((\S+),(\S+),(\S+)\) dev_ms=(\S+) cpu_ms=(\S+)", txt)
d1e, sum1 = 0.0, 0.0
if m1:
    print(m1.group(0))
    d1e = float(m1.group(2))
    sum1 = max(abs(float(m1.group(k))) for k in (4, 5, 6))
elif a.require_1e:
    print(txt[-3000:])
    sys.exit("no 'QM GPU grad1e check' line: the 1e gradient did not run on the device")
e_cpu, g_cpu, _ = run("none", a.extra, a.xyz, False)
d_tot = max(abs(x - y) for x, y in zip(g_dev, g_cpu)) * BOHR_IN_ANGSTROM  # Eh/Bohr
ok = d2e <= a.tol_2e and sum_dev <= a.tol_sum and d_tot <= a.tol_total and d1e <= a.tol_1e and sum1 <= a.tol_sum
if m1:
    print(f"1e: max|dev-cpu| {d1e:.2e} Eh/Bohr, |sum_dev| {sum1:.1e}")
print(f"2e: max|dev-cpu| {d2e:.2e} Eh/Bohr, |sum_dev| {sum_dev:.1e}; total gradient device run vs CPU run "
      f"{d_tot:.2e} Eh/Bohr; E diff {abs(e_dev - e_cpu):.1e} Eh")

if a.fd:
    lines = open(a.xyz).read().splitlines()
    nat = int(lines[0])
    atoms = [l.split() for l in lines[2:2 + nat]]
    h = a.fd_step
    dmax = 0.0
    for i in range(nat):
        for k in range(3):
            es = []
            for sgn in (1, -1):
                at = [r[:] for r in atoms]
                at[i][k + 1] = repr(float(at[i][k + 1]) + sgn * h)
                p = os.path.join(work, "fd.xyz")
                with open(p, "w") as f:
                    f.write(f"{nat}\n\n" + "\n".join(" ".join(r[:4]) for r in at) + "\n")
                es.append(run("cuda", a.extra + " " + a.dev_extra + " -qm.scf_warm_start false", p, False)[0])
            fd = (es[0] - es[1]) / (2 * h)  # Eh/Angstrom
            dmax = max(dmax, abs(fd - g_dev[3 * i + k]))
    print(f"device run: max|g - FD(own energy, h = {h} Angstrom)| = {dmax:.2e} Eh/Angstrom")
    ok &= dmax <= a.tol_fd
print("ok" if ok else "FAIL")
sys.exit(0 if ok else 1)
