#!/usr/bin/env python3
"""WP10 density fitting (RI-J / RI-JK) checks for the native HF engine. Claude Generated (Oct 2026).

Modes (--check):
  energy          curcuma -sp -method hf -qm.ri <ri> vs the PySCF total energy in reference.json
                  (--tol, default 1e-8 Eh; same method and auxiliary basis, so near rounding) and vs ORCA 6.1
                  (--tol-orca, default 1e-7 Eh; ORCA uses its own thresholds).
  gradient        -sp -gradient -dump_gradient with -qm.ri <ri> vs the PySCF analytic DF gradient in reference_grad.json
                  (--tol-grad, default 5e-10 Eh/Bohr; scripts/ri_gradient_reference.py), plus sum of the components ~ 0.
  fd              dump_qm_gradient --ri <ri> --fd (--dump): analytic RI gradient vs the central finite difference of
                  curcuma's own RI energy (1e-7 Eh/Bohr, as ctest -L qm_grad; observed ~3e-9 = FD truncation).
  opt             -opt -qm.ri jk (tight) must reach the PySCF RI-JK minimum in reference_grad.json ("opt"):
                  energy 1e-9 Eh, interatomic distances 1e-4 Angstrom.
  gradient_error  a gradient that cannot be formed (3-index densities above -qm.ri_max_memory_mb) must fail loud
                  (non-zero exit, clear message), never return a number. (Before Oct 2026 every RI gradient
                  request was this error; the analytic RI gradient replaced it.)
  hf3c            -method hf-3c ignores -qm.ri (warning) and returns the exact-integral energy.
  gpu_grad        -gpu cuda -gradient with CURCUMA_QM_GPU_RI_GRAD_CHECK=1: device 3-index RI gradient vs the CPU one
                  (<= 1e-10 Eh/Bohr), full gradient vs PySCF (--tol-grad), no fallback. Exit 77 without plugin/device.
  gpu             -gpu cuda with CURCUMA_QM_GPU_RI_CHECK=1: device J/K of the first build vs the CPU contraction
                  (<= 1e-10), energy vs PySCF, no fallback. Exit 77 (ctest SKIP) without plugin/device.
Plain python3 (no numpy): parses E_elec and E_nn (verbosity 2, 10 decimals).
"""
import argparse, json, os, re, shutil, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--curcuma", required=True)
ap.add_argument("--ref", required=True)
ap.add_argument("--mol", required=True)
ap.add_argument("--ri", default="jk")
ap.add_argument("--check", default="energy")
ap.add_argument("--set", default="pyscf", help="reference set: pyscf | pyscf_cartesian")
ap.add_argument("--tol", type=float, default=1e-8)
ap.add_argument("--tol-orca", type=float, default=1e-7)
ap.add_argument("--extra", default="", help="extra curcuma flags (one string, use --extra=...)")
ap.add_argument("--dump", default="", help="dump_qm_gradient binary (--check fd)")
ap.add_argument("--tol-grad", type=float, default=5e-10)
a = ap.parse_args()

ref = json.load(open(a.ref))[a.mol]
xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), ref["xyz"])


def run(args, env=None):
    cmd = [a.curcuma, "-sp", xyz, "-no_bmt"] + args + a.extra.split()
    e = dict(os.environ)
    if env:
        e.update(env)
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000, env=e)
    return out.returncode, re.sub(r"\x1b\[[0-9;]*m", "", out.stdout + out.stderr)


def energy(txt):
    ee = re.search(r"E_elec: (-?\d+\.\d+) Eh", txt)
    en = re.search(r"E_nn\s*: (-?\d+\.\d+) Eh", txt)
    if not ee or not en:
        print(txt[-3000:])
        sys.exit("no E_elec / E_nn in the output")
    return float(ee.group(1)) + float(en.group(1))


hf = ["-method", "hf", "-qm.scf_threshold", "1e-10", "-verbosity", "2"]

if a.check == "energy":
    rc, txt = run(hf + ["-qm.ri", a.ri])
    if rc != 0:
        print(txt[-3000:])
        sys.exit(f"curcuma exited with {rc}")
    if a.ri != "off" and f"RI-{a.ri.upper()} with" not in txt:
        print(txt[-3000:])
        sys.exit("the RI path did not run")
    e = energy(txt)
    d = e - ref[a.set][a.ri]
    ok = abs(d) <= a.tol
    line = f"{a.mol} ri={a.ri} set={a.set} E={e:.10f} ref={ref[a.set][a.ri]:.10f} diff={d:+.2e} (tol {a.tol:.0e})"
    if a.set == "pyscf":
        do = e - ref["orca"][a.ri]
        ok = ok and abs(do) <= a.tol_orca
        line += f"  ORCA diff={do:+.2e} (tol {a.tol_orca:.0e})"
    print(line + ("  ok" if ok else "  FAIL"))
    sys.exit(0 if ok else 1)

BOHR = 0.529177210903  # curcuma's CODATA 2018 value; -dump_gradient writes Eh/Angstrom


def read_dump(path):
    """-dump_gradient file -> (energy, gradient Eh/Bohr atom-major)."""
    e, g = None, []
    for line in open(path):
        if line.startswith("# energy"):
            e = float(line.split()[2])
        elif not line.startswith("#") and line.strip():
            g += [float(x) * BOHR for x in line.split()]
    return e, g


def gradient_ref():
    return json.load(open(os.path.join(os.path.dirname(os.path.abspath(a.ref)), "reference_grad.json")))


if a.check == "gradient":
    gref = gradient_ref()["cases"][a.mol + ("_cartesian" if a.set == "pyscf_cartesian" else "")]
    xyz = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", gref["xyz"])
    with tempfile.TemporaryDirectory() as tmp:
        dump = os.path.join(tmp, "g.txt")
        cmd = [a.curcuma, "-sp", xyz, "-no_bmt", "-method", "hf", "-qm.ri", a.ri, "-qm.scf_threshold", "1e-10",
               "-gradient", "-dump_gradient", dump, "-verbosity", "1"] + a.extra.split()
        if gref.get("cartesian"):
            cmd += ["-qm.cartesian_d", "true"]
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=3000, cwd=tmp)
        if out.returncode != 0 or not os.path.exists(dump):
            print((out.stdout + out.stderr)[-3000:])
            sys.exit(f"curcuma exited with {out.returncode} / no gradient file")
        e, g = read_dump(dump)
    r = gref[a.ri]
    d = max(abs(x - y) for x, y in zip(g, r["gradient"]))
    sums = [abs(sum(g[k::3])) for k in range(3)]
    de = e - r["energy"]
    ok = d <= a.tol_grad and max(sums) <= 1e-10 and abs(de) <= 1e-8 and len(g) == len(r["gradient"])
    print(f"{a.mol} ri={a.ri} max|g - PySCF| = {d:.2e} Eh/Bohr (tol {a.tol_grad:.0e}, max|g| {r['max_abs']:.3e}), "
          f"max|sum| {max(sums):.1e}, dE {de:+.1e}  {'ok' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)

if a.check == "fd":
    cmd = [a.dump, xyz, "--ri", a.ri, "--fd", "--threads", "4"] + a.extra.split()
    got = json.loads(subprocess.run(cmd, capture_output=True, text=True, timeout=3000, check=True).stdout)
    d = max(abs(x - y) for x, y in zip(got["gradient"], got["fd_gradient"]))
    ok = d <= 1e-7
    print(f"{a.mol} ri={a.ri} max|g - FD(own RI energy)| = {d:.2e} Eh/Bohr (h = {got['fd_step_bohr']} Bohr, tol 1e-7)  "
          f"{'ok' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)

if a.check == "opt":
    o = gradient_ref()["opt"][a.mol]
    start = os.path.join(os.path.dirname(os.path.abspath(a.ref)), "..", o["start"])
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(start, os.path.join(tmp, "start.xyz"))
        subprocess.run([a.curcuma, "-opt", "start.xyz", "-method", "hf", "-qm.ri", "jk", "-no_bmt", "-verbosity", "0",
                        "-qm.scf_threshold", "1e-9", "-opt.convergence_preset", "tight"] + a.extra.split(),
                       cwd=tmp, capture_output=True, text=True, timeout=3000, check=True)
        subprocess.run([a.curcuma, "-sp", "start.opt.xyz", "-method", "hf", "-qm.ri", "jk", "-no_bmt", "-verbosity", "0",
                        "-qm.scf_threshold", "1e-10", "-gradient", "-dump_gradient", "g.txt"] + a.extra.split(),
                       cwd=tmp, capture_output=True, text=True, timeout=3000, check=True)
        e, g = read_dump(os.path.join(tmp, "g.txt"))
        L = open(os.path.join(tmp, "start.opt.xyz")).read().splitlines()
        nat = int(L[0])
        x = [list(map(float, l.split()[1:4])) for l in L[2:2 + nat]]
    xr = [o["geometry_angstrom"][3 * i:3 * i + 3] for i in range(nat)]
    dist = lambda X, i, j: sum((X[i][k] - X[j][k]) ** 2 for k in range(3)) ** 0.5
    dd = max(abs(dist(x, i, j) - dist(xr, i, j)) for i in range(nat) for j in range(i + 1, nat))
    ok = abs(e - o["energy"]) <= 1e-9 and dd <= 1e-4
    print(f"{a.mol} -opt -qm.ri jk: E {e:.10f} ref {o['energy']:.10f} diff {e - o['energy']:+.1e}, max|d r_ij| {dd:.1e} A, "
          f"max|g| at the curcuma minimum {max(map(abs, g)):.1e} Eh/Bohr  {'ok' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)

if a.check == "gradient_error":
    # H2O/def2-SVP with def2-universal-JKFIT: tensor 0.52 MB (naux 113, n 24), the gradient needs as much again
    rc, txt = run(["-method", "hf", "-qm.ri", a.ri, "-gradient", "-qm.ri_max_memory_mb", "0.8", "-verbosity", "1"])
    msg = "RI gradient not possible" in txt
    print(f"exit code {rc}, message {'found' if msg else 'missing'}")
    sys.exit(0 if (rc != 0 and msg) else 1)

if a.check == "hf3c":
    rc1, t1 = run(["-method", "hf-3c", "-qm.ri", a.ri, "-verbosity", "1"])
    rc0, t0 = run(["-method", "hf-3c", "-verbosity", "1"])
    m1 = re.search(r"Single Point Energy = (-?\d+\.\d+)", t1)
    m0 = re.search(r"Single Point Energy = (-?\d+\.\d+)", t0)
    warned = "hf-3c: -qm.ri" in t1 and "ignored" in t1
    ok = rc0 == 0 and rc1 == 0 and m0 and m1 and m0.group(1) == m1.group(1) and warned
    print(f"hf-3c with -qm.ri {a.ri}: {m1.group(1) if m1 else 'NA'}, without: {m0.group(1) if m0 else 'NA'}, "
          f"warning {'found' if warned else 'missing'}")
    sys.exit(0 if ok else 1)

if a.check == "gpu_grad":
    gref = gradient_ref()["cases"][a.mol][a.ri]
    with tempfile.TemporaryDirectory() as tmp:
        dump = os.path.join(tmp, "g.txt")
        rc, txt = run(hf + ["-qm.ri", a.ri, "-gpu", "cuda", "-qm.threads", "4", "-gradient", "-dump_gradient", dump],
                      env={"CURCUMA_QM_GPU_RI_GRAD_CHECK": "1"})
        if "no GPU plugin" in txt or "not loadable" in txt or "cudaSetDevice failed" in txt:
            print("no usable CUDA device/plugin: skipped")
            sys.exit(77)
        m = re.search(r"QM GPU RI grad check: natoms=(\d+) max\|dev-cpu\|=(\S+)", txt)
        if rc != 0 or not m or not os.path.exists(dump):
            print(txt[-3000:])
            sys.exit("device RI gradient did not run (no check line)")
        if "GPU fallbacks during this run" in txt:
            print(txt[-3000:])
            sys.exit("the run reported a GPU fallback")
        e, g = read_dump(dump)
    dd = float(m.group(2))
    d = max(abs(x - y) for x, y in zip(g, gref["gradient"]))
    ok = dd <= 1e-10 and d <= a.tol_grad
    print(f"{a.mol} ri={a.ri} device vs CPU 3-index gradient: max|d|={dd:.2e}; full gradient vs PySCF {d:.2e} Eh/Bohr  "
          f"{'ok' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)

if a.check == "gpu":
    rc, txt = run(hf + ["-qm.ri", a.ri, "-gpu", "cuda", "-qm.threads", "4"], env={"CURCUMA_QM_GPU_RI_CHECK": "1"})
    if "no GPU plugin" in txt or "not loadable" in txt or "cudaSetDevice failed" in txt:
        print("no usable CUDA device/plugin: skipped")
        sys.exit(77)
    m = re.search(r"QM GPU RI check: n=(\d+) naux=(\d+) max\|dJ\|=(\S+) max\|dK\|=(\S+)", txt)
    if rc != 0 or not m:
        print(txt[-3000:])
        sys.exit("device RI path did not run (no check line)")
    if "GPU fallbacks during this run" in txt:
        print(txt[-3000:])
        sys.exit("the run reported a GPU fallback")
    dJ, dK = float(m.group(3)), float(m.group(4))
    e = energy(txt)
    d = e - ref["pyscf"][a.ri]
    ok = dJ <= 1e-10 and dK <= 1e-10 and abs(d) <= a.tol
    print(f"{a.mol} ri={a.ri} device vs CPU: max|dJ|={dJ:.2e} max|dK|={dK:.2e}; E={e:.10f} diff PySCF={d:+.2e}  "
          f"{'ok' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)

sys.exit("unknown --check " + a.check)
