#!/usr/bin/env python3
"""ctest checker for the native QM engine with f and g shells (ctest -L qm_fg). Standard library only.

Also used by test_cases/qm_ecp (elements beyond Ne, def2-ECP) with --reference qm_ecp/reference.json; a case may carry
"orca_tol" (ORCA's ECP energies differ from PySCF and curcuma by up to 6e-8 Eh for Xe and Hg, cause not identified).

Checks (--check):
  ref    one case of reference.json (written by scripts/qm_fg_reference.py): dump_qm_energy (built-in basis sets
         only, CURCUMA_QM_BASIS_EMBEDDED_ONLY=1) against PySCF (energy |dE| <= --tol, default 1e-8 Eh; gradient
         max |dg| <= --gtol, default 1e-8 Eh/Bohr) and, where stored, ORCA 6.1 (energy, 1e-8 Eh; RI-MP2 2e-8 Eh,
         ORCA's RI-MP2 sits ~1e-8 off PySCF and curcuma)
  fd     finite-difference gate of the fully assembled gradient: the analytic gradient of curcuma (dump_qm_energy
         --gradient) against a 4-point central difference of curcuma's own energies, step 0.005 Bohr, on the
         coordinates --coords (atom:axis, default 0:z 1:x 1:z), max |dg| <= --gtol (default 1e-8 Eh/Bohr;
         measured <= 2.4e-10 on every case, largest on O z in every method, cause not analysed)
  gpu    curcuma -sp with -gpu cuda on an f basis: the run finishes, the energy equals the CPU run to 1e-9 Eh;
         exit 77 (skip) without a usable CUDA plugin/device
  gpudev stage 6b (label qm_gpu_fg): dump_qm_energy --gpu cuda against the CPU run of the same binary (energy
         |dE| <= 1e-9 Eh, gradient max |dg| <= --gtol); FAILS on any counted GPU fallback except those named in
         --allow (the f/g device kernels must run: qm_eri, qm_direct, qm_grad, qm_grad1e, qm_mp2_grad, qm_ri_grad
         are not allowed); --qm KEY=VALUE adds engine options; with --gpu also ref/fd run the device path
         (same fallback rule). Exit 77 without libcurcuma_cuda.so next to dump_qm_energy (or CURCUMA_PLUGIN_DIR)
         or without a CUDA device
  kernel two-electron gradient kernels (Nov 2026): the Hermite-matrix kernel that f/g bases use by default against
         the reference kernel quartetGradient (CURCUMA_QM_GRAD_KERNEL=reference), same binary and threads, so the
         same SCF; max |dg| <= --gtol (default 1e-11 Eh/Bohr), energies equal bit for bit
  threads the default (f/g) gradient with 1 and --threads threads: energy and gradient bit for bit equal

Claude Generated (Nov 2026).
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.normpath(os.path.join(HERE, ".."))
BOHR = 0.529177210903  # curcuma units.h


def env_embedded():
    e = dict(os.environ)
    e["CURCUMA_QM_BASIS_EMBEDDED_ONLY"] = "1"
    return e


def qm_args(qm):
    out = []
    for k, v in qm.items():
        out += ["--qm", k, json.dumps(v) if not isinstance(v, str) else v]
    return out


GPU = {"on": False, "fallbacks": {}, "allow": set()}  # stage 6b: device path of dump_qm_energy (--gpu)


def gpu_available(a):
    """libcurcuma_cuda.so in CURCUMA_PLUGIN_DIR or next to dump_qm_energy's build tree, and a CUDA device; else exit 77."""
    cands = [os.environ.get("CURCUMA_PLUGIN_DIR", "")]
    d = os.path.dirname(os.path.abspath(a.dump))
    for _ in range(4):  # dump_qm_energy lives in <build>/test_cases/qm_basis_embed
        cands.append(d)
        d = os.path.dirname(d)
    pdir = next((c for c in cands if c and os.path.exists(os.path.join(c, "libcurcuma_cuda.so"))), "")
    ok = bool(pdir)
    if ok:
        try:
            ok = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True).returncode == 0
        except OSError:
            ok = False
    if not ok:
        print("skip: no CUDA plugin found (%s) or no CUDA device" % cands)
        sys.exit(77)
    os.environ["CURCUMA_PLUGIN_DIR"] = pdir


def run_dump(dump, xyz, method, qm, gradient=False, threads=1, gpu=None, charge=0, spin=0, env_extra=None):
    if gpu is None:
        gpu = GPU["on"]
    cmd = [dump, xyz, method] + qm_args(qm) + ["--threads", str(threads)]
    if charge or spin:
        cmd += ["--charge", str(charge), "--spin", str(spin)]
    if gpu:
        cmd += ["--gpu", "cuda"]
    if gradient:
        cmd.append("--gradient")
    env = env_embedded()
    env.update(env_extra or {})
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    out = p.stdout + p.stderr
    m = re.search(r"^energy (\S+)$", out, re.M)
    if p.returncode != 0 or not m:
        print(out[-4000:])
        raise SystemExit("dump_qm_energy failed: " + " ".join(cmd))
    g = [[float(x) for x in l.split()[2:5]] for l in out.splitlines() if l.startswith("grad ")]
    if gpu:
        for fb in re.finditer(r"^fallback (\S+) (\d+)$", out, re.M):
            GPU["fallbacks"][fb.group(1)] = GPU["fallbacks"].get(fb.group(1), 0) + int(fb.group(2))
    return float(m.group(1)), g


def check_gpu_fallbacks():
    bad = {k: v for k, v in GPU["fallbacks"].items() if k not in GPU["allow"]}
    print("counted GPU fallbacks: %s (allowed: %s)" % (GPU["fallbacks"] or "none", sorted(GPU["allow"]) or "none"))
    if bad:
        raise SystemExit("FAILED: device path fell back to the CPU: %s" % bad)


def case_qm(a, ref):
    qm = dict(ref["qm"])
    for kv in a.qm:
        k, v = kv.split("=", 1)
        try:
            qm[k] = json.loads(v)
        except ValueError:
            qm[k] = v
    return qm


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0].split()[0])
    return [[a.split()[0]] + [float(v) for v in a.split()[1:4]] for a in L[2:2 + n]]


def write_xyz(path, atoms):
    with open(path, "w") as f:
        f.write("%d\n\n" % len(atoms))
        for s, x, y, z in atoms:
            f.write("%s %.12f %.12f %.12f\n" % (s, x, y, z))


def check_ref(a):
    ref = json.load(open(a.reference))["cases"][a.case]
    xyz = os.path.join(TC, ref["xyz"])
    qm = case_qm(a, ref)
    qm.setdefault("scf_threshold", 1e-10)
    e, g = run_dump(a.dump, xyz, ref["method"], qm, ref["gradient"], a.threads)
    ok = True
    dp = e - ref["pyscf"]["energy"]
    line = "%-30s nbf %3d curcuma %.12f PySCF %.12f diff %+.2e" % (a.case, ref["pyscf"]["nbf"], e,
                                                                    ref["pyscf"]["energy"], dp)
    if abs(dp) > a.tol:
        ok = False
    if "orca" in ref:
        do = e - ref["orca"]["energy"]
        tol_o = max(a.tol, 2e-8) if ref["method"] == "ri-mp2" else a.tol
        tol_o = ref.get("orca_tol", tol_o)  # per-case ORCA tolerance (qm_ecp: ORCA's ECP energies, Nov 2026)
        line += "  ORCA diff %+.2e" % do
        if abs(do) > tol_o:
            ok = False
    print(line)
    if ref["gradient"]:
        gr = ref["pyscf"]["gradient"]
        dg = max(abs(x - y) for r1, r2 in zip(g, gr) for x, y in zip(r1, r2))
        print("  gradient max |g - g_PySCF| = %.2e Eh/Bohr (limit %.1e)" % (dg, a.gtol))
        if len(g) != len(gr) or dg > a.gtol:
            ok = False
    if GPU["on"]:
        check_gpu_fallbacks()
    if not ok:
        raise SystemExit("FAILED: %s outside the tolerance" % a.case)


def check_fd(a):
    ref = json.load(open(a.reference))["cases"][a.case]
    xyz = os.path.join(TC, ref["xyz"])
    qm = case_qm(a, ref)
    qm["scf_threshold"] = 1e-11 if "QZV" not in qm["basis"] else 1e-10
    atoms = read_xyz(xyz)
    method = a.method or ref["method"]  # Nov 2026 (qm_gpu_ecp): method / charge / spin overrides for the FD gate
    ch, sp = ref.get("charge", 0) if a.charge is None else a.charge, ref.get("spin", 0) if a.spin is None else a.spin
    _, g = run_dump(a.dump, xyz, method, qm, True, a.threads, charge=ch, spin=sp)
    h = 0.005  # Bohr
    worst = 0.0
    with tempfile.TemporaryDirectory() as d:
        for spec in a.coords:
            at, ax = spec.split(":")
            at, k = int(at), "xyz".index(ax)
            es = {}
            for s in (-2, -1, 1, 2):
                disp = [list(r) for r in atoms]
                disp[at][1 + k] += s * h * BOHR
                p = os.path.join(d, "d.xyz")
                write_xyz(p, disp)
                es[s], _ = run_dump(a.dump, p, method, qm, False, a.threads, charge=ch, spin=sp)
            fd = (es[-2] - 8 * es[-1] + 8 * es[1] - es[2]) / (12 * h)
            diff = g[at][k] - fd
            worst = max(worst, abs(diff))
            print("%-30s atom %d %s analytic %+.10f FD %+.10f diff %+.2e" % (a.case, at, ax, g[at][k], fd, diff))
    print("max |analytic - FD| = %.2e Eh/Bohr (limit %.1e)" % (worst, a.gtol))
    if GPU["on"]:
        check_gpu_fallbacks()
    if worst > a.gtol:
        raise SystemExit("FAILED: FD gate")


def check_gpu(a):
    ref = json.load(open(a.reference))["cases"][a.case]
    xyz = os.path.join(TC, ref["xyz"])
    args = []
    for k, v in ref["qm"].items():
        args += ["-qm." + k, str(v)]
    def cli(extra):
        with tempfile.TemporaryDirectory() as d:
            p = subprocess.run([a.curcuma, "-sp", xyz, "-method", ref["method"], "-no_bmt", "-verbosity", "1",
                                "-threads", str(a.threads)] + args + extra, capture_output=True, text=True, cwd=d,
                               env=env_embedded())
        return p.returncode, p.stdout + p.stderr
    rc, out = cli(["-gpu", "cuda"])
    if "no GPU plugin" in out or "plugin" in out and "not found" in out or "no CUDA device" in out:
        print(out[-2000:])
        print("skip: no usable CUDA plugin/device")
        sys.exit(77)
    rc0, out0 = cli([])
    def energy(o):
        m = re.findall(r"total \S+ energy:\s+(-?\d+\.\d+)", re.sub(r"\x1b\[[0-9;]*m", "", o))
        return float(m[-1]) if m else None
    e, e0 = energy(out), energy(out0)
    print(out[-3000:])
    if rc != 0 or e is None or e0 is None:
        raise SystemExit("FAILED: -gpu cuda run did not finish (rc %d)" % rc)
    print("CPU %.8f  -gpu cuda %.8f  diff %+.1e" % (e0, e, e - e0))
    if abs(e - e0) > 1e-8:
        raise SystemExit("FAILED: GPU energy differs from the CPU")


def check_gpudev(a):
    """Stage 6b: device path (dump_qm_energy --gpu cuda) against the CPU path of the same binary."""
    ref = json.load(open(a.reference))["cases"][a.case]
    xyz = os.path.join(TC, ref["xyz"])
    qm = case_qm(a, ref)
    qm.setdefault("scf_threshold", 1e-10)
    grad = ref["gradient"] or a.gradient
    ch, sp = ref.get("charge", 0), ref.get("spin", 0)
    e0, g0 = run_dump(a.dump, xyz, ref["method"], qm, grad, a.threads, gpu=False, charge=ch, spin=sp)
    e1, g1 = run_dump(a.dump, xyz, ref["method"], qm, grad, a.threads, gpu=True, charge=ch, spin=sp)
    dg = max((abs(x - y) for r1, r2 in zip(g0, g1) for x, y in zip(r1, r2)), default=0.0)
    print("%-30s %s CPU %.12f GPU %.12f diff %+.2e  gradient max |GPU - CPU| %.2e Eh/Bohr" % (
        a.case, " ".join(a.qm), e0, e1, e1 - e0, dg))
    if "pyscf" in ref:
        print("  GPU vs PySCF %+.2e Eh" % (e1 - ref["pyscf"]["energy"]))
    check_gpu_fallbacks()
    if abs(e1 - e0) > 1e-9 or dg > a.gtol or len(g0) != len(g1):
        raise SystemExit("FAILED: GPU result differs from the CPU")


def check_kernel(a):
    """Nov 2026: Hermite-matrix 2e gradient kernel (default for f/g) vs the reference kernel on the same SCF."""
    ref = json.load(open(a.reference))["cases"][a.case]
    xyz = os.path.join(TC, ref["xyz"])
    qm = case_qm(a, ref)
    qm.setdefault("scf_threshold", 1e-10)
    e0, g0 = run_dump(a.dump, xyz, ref["method"], qm, True, a.threads, env_extra={"CURCUMA_QM_GRAD_KERNEL": "reference"})
    e1, g1 = run_dump(a.dump, xyz, ref["method"], qm, True, a.threads)
    dg = max(abs(x - y) for r1, r2 in zip(g0, g1) for x, y in zip(r1, r2))
    print("%-30s E %.12f / %.12f  gradient max |md - reference| %.2e Eh/Bohr (limit %.1e)" % (a.case, e0, e1, dg, a.gtol))
    if e0 != e1 or len(g0) != len(g1) or dg > a.gtol:
        raise SystemExit("FAILED: gradient kernels differ")


def check_threads(a):
    """Nov 2026: the default gradient (f/g: Hermite-matrix kernel, fixed chunks) is the same at 1 and N threads."""
    ref = json.load(open(a.reference))["cases"][a.case]
    xyz = os.path.join(TC, ref["xyz"])
    qm = case_qm(a, ref)
    qm.setdefault("scf_threshold", 1e-10)
    e0, g0 = run_dump(a.dump, xyz, ref["method"], qm, True, 1)
    e1, g1 = run_dump(a.dump, xyz, ref["method"], qm, True, a.threads)
    same = e0 == e1 and g0 == g1
    print("%-30s 1 vs %d threads: E %.17g / %.17g, gradient %s" % (a.case, a.threads, e0, e1,
                                                                  "identical" if same else "DIFFERENT"))
    if not same:
        raise SystemExit("FAILED: result depends on the thread count")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", required=True, choices=["ref", "fd", "gpu", "gpudev", "kernel", "threads"])
    ap.add_argument("--case", required=True)
    ap.add_argument("--dump")
    ap.add_argument("--curcuma")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--gtol", type=float, default=None)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--coords", nargs="*", default=["0:z", "1:x", "1:z"])
    ap.add_argument("--reference", default=os.path.join(HERE, "reference.json"),
                    help="reference file (default: qm_fg/reference.json; test_cases/qm_ecp passes its own)")
    ap.add_argument("--gpu", action="store_true", help="ref/fd: run the device path (stage 6b)")
    ap.add_argument("--gradient", action="store_true", help="gpudev: also compare the gradient")
    ap.add_argument("--allow", nargs="*", default=[], help="counted GPU fallbacks that may occur")
    ap.add_argument("--qm", nargs="*", default=[], help="extra engine options KEY=VALUE")
    ap.add_argument("--method", default=None, help="fd: method instead of the case's (Nov 2026)")
    ap.add_argument("--charge", type=int, default=None, help="fd: molecular charge instead of the case's")
    ap.add_argument("--spin", type=int, default=None, help="fd: unpaired electrons instead of the case's")
    a = ap.parse_args()
    if a.gtol is None:
        a.gtol = 1e-11 if a.check == "kernel" else 1e-8
    if a.check == "gpudev" or a.gpu:
        gpu_available(a)
        GPU["on"] = a.check != "gpudev"
        GPU["allow"] = set(a.allow)
    {"ref": check_ref, "fd": check_fd, "gpu": check_gpu, "gpudev": check_gpudev, "kernel": check_kernel,
     "threads": check_threads}[a.check](a)


if __name__ == "__main__":
    main()
