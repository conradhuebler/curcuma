#!/usr/bin/env python3
"""ctest checker for the open-shell native QM engine (stage 5b: UHF/UKS; ctest -L "qm_uhf|qm_uks"). Standard library only.

Checks (--check):
  kernel   test_qm_uks_kernel against libxc (reference.json "kernel", scripts/qm_uhf_reference.py): eps, v_rho, v_sigma
           per point, |d| <= 1e-10 (1 + |ref|); a spin channel with rho_s < 1e-10 is not compared (libxc's density
           clamp, see check_kernel) and the other values of such a point to 1e-7; plus the program's own finite-difference check of v (fd_max <= 1e-5,
           noise-weighted, step 1e-4 x) and the closed-shell limit rho_a = rho_b against the RKS kernel (cs_max <= 1e-12 relative)
  ref      one case: dump_qm_energy (built-in basis sets only) against PySCF (energy |dE| <= --tol, default 1e-8 Eh;
           <S^2> |d| <= 1e-6; gradient max |dg| <= --gtol, default 1e-8 Eh/Bohr) and ORCA 6.1 where stored (energy,
           --tol; <S^2> 1e-5)
  fd       FD gate of the fully assembled gradient: analytic gradient vs a 4-point central difference of curcuma's own
           energies (step 0.005 Bohr) on --coords, max |dg| <= --gtol
  rhf      closed-shell limit: -qm.reference uhf on a singlet (water) equals the RHF/RKS energy to 1e-9 Eh and gradient
           to 1e-8 Eh/Bohr, <S^2> = 0 to 1e-9 (--case hf|pbe|b3lyp)
  direct   -qm.scf_direct on equals the stored-tensor energy of a case to 1e-9 Eh (and gradient to 1e-8 Eh/Bohr)
  reject   open-shell input to a closed-shell-only method or a contradicting spin state is an error with a message
           (--case mp2|ri-mp2|ccsd|lmp2|b2plyp|tpss|scan|parity|rohf)
  gpu      stage 6a: dump_qm_energy --gpu cuda against the CPU run of the same case: energy |dE| <= 1e-9 Eh, <S^2>
           1e-8, gradient max |dg| <= --gtol; the only counted fallback allowed is qm_uhf_scf (host SCF loop), so J/K,
           V_xc and every gradient part ran on the device; --extra-qm adds -qm settings (e.g. scf_direct on)
  --gpu    (with ref / fd) run the device path (dump_qm_energy --gpu cuda) and require the same fallback rule;
           every GPU check exits 77 (skip) without a CUDA plugin next to the binary or without a CUDA device
  opt      curcuma -opt of the distorted CH3 radical with hf: exit 0, energy below the start, last gradient norm <= 1e-4
           Eh/Bohr, planar minimum (C within 0.01 A of the H3 plane)

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


GPU = {"on": False, "fallbacks": {}}
ALLOWED_GPU_FALLBACKS = {"qm_uhf_scf"}  # the UHF/UKS SCF loop itself stays on the host (stage 6a)


def gpu_available(a):
    """CUDA plugin next to curcuma (CURCUMA_PLUGIN_DIR) and a device visible to nvidia-smi; else exit 77."""
    pdir = a.plugin_dir or (os.path.dirname(a.curcuma) if a.curcuma else "")
    ok = bool(pdir) and os.path.exists(os.path.join(pdir, "libcurcuma_cuda.so"))
    if ok:
        try:
            ok = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True).returncode == 0
        except OSError:
            ok = False
    if not ok:
        print("skip: no CUDA plugin in %r or no CUDA device" % pdir)
        sys.exit(77)
    os.environ["CURCUMA_PLUGIN_DIR"] = pdir


def run_dump(dump, xyz, method, qm, charge=0, spin=0, gradient=False, threads=1, expect_fail=False, gpu=None):
    cmd = [dump, xyz, method, "--charge", str(charge), "--spin", str(spin), "--threads", str(threads)]
    if gpu is None:
        gpu = GPU["on"]
    if gpu:
        cmd += ["--gpu", "cuda"]
    for k, v in qm.items():
        cmd += ["--qm", k, json.dumps(v) if not isinstance(v, str) else v]
    if gradient:
        cmd.append("--gradient")
    p = subprocess.run(cmd, capture_output=True, text=True, env=env_embedded())
    out = p.stdout + p.stderr
    if expect_fail:
        return p.returncode, out
    m = re.search(r"^energy (\S+)$", out, re.M)
    if p.returncode != 0 or not m:
        print(out[-4000:])
        raise SystemExit("dump_qm_energy failed: " + " ".join(cmd))
    g = [[float(x) for x in l.split()[2:5]] for l in out.splitlines() if l.startswith("grad ")]
    if gpu:
        for fb in re.finditer(r"^fallback (\S+) (\d+)$", out, re.M):
            GPU["fallbacks"][fb.group(1)] = GPU["fallbacks"].get(fb.group(1), 0) + int(fb.group(2))
    s2 = re.search(r"^s2 (\S+)$", out, re.M)
    return float(m.group(1)), g, (float(s2.group(1)) if s2 else None)


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0].split()[0])
    return [[a.split()[0]] + [float(v) for v in a.split()[1:4]] for a in L[2:2 + n]]


def write_xyz(path, atoms):
    with open(path, "w") as f:
        f.write("%d\n\n" % len(atoms))
        for s, x, y, z in atoms:
            f.write("%s %.12f %.12f %.12f\n" % (s, x, y, z))


def load_case(a):
    return json.load(open(a.reference))["cases"][a.case]


def maxdiff(g1, g2):
    return max(abs(x - y) for r1, r2 in zip(g1, g2) for x, y in zip(r1, r2))


def check_kernel(a):
    ref = json.load(open(a.reference))["kernel"]
    p = subprocess.run([a.exe, a.reference], capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stdout[-2000:], p.stderr[-2000:])
        raise SystemExit("test_qm_uks_kernel failed")
    got = json.loads(p.stdout)
    ok = True
    for name, r in ref["kernels"].items():
        g = got[name]
        worst = 0.0
        for key in ("eps", "vrho", "vsigma"):
            if key not in r:
                continue
            for i, (x, y) in enumerate(zip(g[key], r[key])):
                xs = x if isinstance(x, list) else [x]
                ys = y if isinstance(y, list) else [y]
                # A spin channel below 1e-10: libxc raises each spin density to the functional's dens_threshold (e.g.
                # 1e-12 for GGA_C_PBE, measured) and curcuma does not (as for the closed-shell kernels), so the
                # derivatives of that channel differ (ill-defined at rho_s = 0) and eps differs by up to 4e-10 Eh.
                small = [ref["points"][i][0] < 1e-10, ref["points"][i][1] < 1e-10]
                skip = {"eps": [False], "vrho": small, "vsigma": [small[0], small[0] or small[1], small[1]]}[key]
                tol = 1e-7 if (small[0] or small[1]) else 1e-10
                for c, (u, v) in enumerate(zip(xs, ys)):
                    if skip[c]:
                        continue
                    d = abs(u - v) / (1.0 + abs(v))
                    worst = max(worst, d)
                    if d > tol:
                        ok = False
                        print("  %s %s[%d] point %d: curcuma %.15e libxc %.15e" % (name, key, c, i, u, v))
        print("%-14s max |d|/(1+|ref|) vs libxc %.2e   fd_max %.2e %s   closed-shell limit %.2e %s" %
              (name, worst, g["fd_max"], g.get("fd_worst"), g["cs_max"], g.get("cs_worst")))
        if g["fd_max"] > 1e-5 or g["cs_max"] > 1e-12:
            ok = False
    if not ok:
        raise SystemExit("FAILED: kernel check")


def check_ref(a):
    ref = load_case(a)
    xyz = os.path.join(TC, ref["xyz"])
    qm = dict(ref["qm"])
    qm.setdefault("scf_threshold", 1e-10)
    e, g, s2 = run_dump(a.dump, xyz, ref["method"], qm, ref["charge"], ref["spin"], ref["gradient"], a.threads)
    ok = True
    P = ref["pyscf"]
    dp = e - P["energy"]
    line = "%-30s nbf %3d curcuma %.12f PySCF %.12f diff %+.2e  <S^2> %.8f (PySCF %.8f, %+.1e)" % (
        a.case, P["nbf"], e, P["energy"], dp, s2 if s2 is not None else float("nan"), P["s2"],
        (s2 - P["s2"]) if s2 is not None else float("nan"))
    if abs(dp) > a.tol or s2 is None or abs(s2 - P["s2"]) > 1e-6:
        ok = False
    if "orca" in ref:
        O = ref["orca"]
        do = e - O["energy"]
        line += "\n  ORCA %.12f diff %+.2e" % (O["energy"], do)
        if O.get("s2") is not None:
            line += "  <S^2> ORCA %.6f (%+.1e)" % (O["s2"], s2 - O["s2"])
            if abs(s2 - O["s2"]) > 1e-5:
                ok = False
        if abs(do) > ref.get("orca_tol", a.tol):
            ok = False
    print(line)
    if ref["gradient"]:
        dg = maxdiff(g, P["gradient"])
        print("  gradient max |g - g_PySCF| = %.2e Eh/Bohr (limit %.1e)" % (dg, a.gtol))
        if len(g) != len(P["gradient"]) or dg > a.gtol:
            ok = False
    if GPU["on"]:
        check_gpu_fallbacks()
    if not ok:
        raise SystemExit("FAILED: %s outside the tolerance" % a.case)


def check_fd(a):
    ref = load_case(a)
    xyz = os.path.join(TC, ref["xyz"])
    qm = dict(ref["qm"])
    qm["scf_threshold"] = 1e-11
    if ref["method"] == "hf-3c":
        qm.pop("basis", None)
    atoms = read_xyz(xyz)
    _, g, _ = run_dump(a.dump, xyz, ref["method"], qm, ref["charge"], ref["spin"], True, a.threads)
    h = 0.005
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
                es[s], _, _ = run_dump(a.dump, p, ref["method"], qm, ref["charge"], ref["spin"], False, a.threads)
            fd = (es[-2] - 8 * es[-1] + 8 * es[1] - es[2]) / (12 * h)
            diff = g[at][k] - fd
            worst = max(worst, abs(diff))
            print("%-30s atom %d %s analytic %+.10f FD %+.10f diff %+.2e" % (a.case, at, ax, g[at][k], fd, diff))
    print("max |analytic - FD| = %.2e Eh/Bohr (limit %.1e)%s" % (worst, a.gtol, " (device path)" if GPU["on"] else ""))
    if GPU["on"]:
        check_gpu_fallbacks()
    if worst > a.gtol:
        raise SystemExit("FAILED: FD gate")


def check_rhf(a):
    xyz = os.path.join(TC, "structures/small/water.xyz")
    qm = {"basis": "def2-SVP", "scf_threshold": 1e-10}
    if a.case != "hf":
        qm["grid"] = "sg1"
    e0, g0, _ = run_dump(a.dump, xyz, a.case, qm, 0, 0, True, a.threads)
    e1, g1, s2 = run_dump(a.dump, xyz, a.case, dict(qm, reference="uhf"), 0, 0, True, a.threads)
    dg = maxdiff(g0, g1)
    print("%s water: RHF %.12f UHF %.12f diff %+.2e; gradient max diff %.2e; <S^2> %.2e" % (a.case, e0, e1, e1 - e0, dg, s2))
    if abs(e1 - e0) > 1e-9 or dg > 1e-8 or s2 is None or abs(s2) > 1e-9:
        raise SystemExit("FAILED: closed-shell limit of the unrestricted path")


def check_direct(a):
    ref = load_case(a)
    xyz = os.path.join(TC, ref["xyz"])
    qm = dict(ref["qm"], scf_threshold=1e-10)
    e0, g0, _ = run_dump(a.dump, xyz, ref["method"], dict(qm, scf_direct="off"), ref["charge"], ref["spin"], ref["gradient"],
                         a.threads)
    e1, g1, _ = run_dump(a.dump, xyz, ref["method"], dict(qm, scf_direct="on"), ref["charge"], ref["spin"], ref["gradient"],
                         a.threads)
    dg = maxdiff(g0, g1) if ref["gradient"] else 0.0
    print("%s: stored %.12f direct %.12f diff %+.2e; gradient max diff %.2e" % (a.case, e0, e1, e1 - e0, dg))
    if abs(e1 - e0) > 1e-9 or dg > 1e-8:
        raise SystemExit("FAILED: integral-direct UHF differs from the stored tensor")


REJECT = {
    "mp2": ("mp2", {}, 0, 1, "closed shell only"),
    "ri-mp2": ("ri-mp2", {}, 0, 1, "closed shell only"),
    "ccsd": ("ccsd", {}, 0, 1, "closed shell only"),
    "lmp2": ("lmp2", {}, 0, 1, "closed shell only"),
    "b2plyp": ("b2plyp", {"grid": "sg1"}, 0, 1, "closed shell only"),
    "tpss": ("tpss", {"grid": "sg1"}, 0, 1, "not implemented for tpss"),
    "scan": ("scan", {"grid": "sg1"}, 0, 1, "not implemented for scan"),
    "parity": ("hf", {}, 0, 0, "pass -spin 1"),
    "rohf": ("hf", {"reference": "rhf"}, 0, 1, "ROHF is not implemented"),
}


def check_reject(a):
    method, qm, charge, spin, text = REJECT[a.case]
    xyz = os.path.join(TC, "structures/small/hydroxyl.xyz")
    rc, out = run_dump(a.dump, xyz, method, dict(qm, basis="def2-SVP"), charge, spin, False, a.threads, expect_fail=True)
    clean = re.sub(r"\x1b\[[0-9;]*m", "", out)
    hit = [l for l in clean.splitlines() if text in l]
    print("%s on OH (charge %d, -spin %d): exit %d; message: %s" % (method, charge, spin, rc, hit[0].strip() if hit else "none"))
    if rc == 0 or not hit or re.search(r"^energy ", clean, re.M):
        print(clean[-3000:])
        raise SystemExit("FAILED: no clear error for open-shell input to %s" % method)


def cli(curcuma, args, cwd):
    p = subprocess.run([curcuma] + args, capture_output=True, text=True, cwd=cwd, env=env_embedded())
    return p.returncode, re.sub(r"\x1b\[[0-9;]*m", "", p.stdout + p.stderr)


def check_gpu_fallbacks():
    bad = {k: v for k, v in GPU["fallbacks"].items() if k not in ALLOWED_GPU_FALLBACKS}
    print("counted GPU fallbacks: %s" % (GPU["fallbacks"] or "none"))
    if bad:
        raise SystemExit("FAILED: device path fell back to the CPU: %s" % bad)


def check_gpu(a):
    ref = load_case(a)
    xyz = os.path.join(TC, ref["xyz"])
    qm = dict(ref["qm"], scf_threshold=1e-10)
    for kv in a.extra_qm:
        k, v = kv.split("=", 1)
        qm[k] = json.loads(v) if v[:1].isdigit() or v in ("true", "false") else v
    e0, g0, s0 = run_dump(a.dump, xyz, ref["method"], qm, ref["charge"], ref["spin"], ref["gradient"], a.threads, gpu=False)
    e1, g1, s1 = run_dump(a.dump, xyz, ref["method"], qm, ref["charge"], ref["spin"], ref["gradient"], a.threads, gpu=True)
    dg = maxdiff(g0, g1) if ref["gradient"] else 0.0
    ds = abs(s1 - s0) if s0 is not None and s1 is not None else float("nan")
    print("%-30s %s CPU %.12f GPU %.12f diff %+.2e  <S^2> diff %.1e  gradient max |GPU - CPU| %.2e Eh/Bohr" % (
        a.case, " ".join(a.extra_qm), e0, e1, e1 - e0, ds, dg))
    if "pyscf" in ref:
        print("  GPU vs PySCF %+.2e Eh%s" % (e1 - ref["pyscf"]["energy"],
              ("  gradient max |GPU - PySCF| %.2e" % maxdiff(g1, ref["pyscf"]["gradient"])) if ref["gradient"] else ""))
    check_gpu_fallbacks()
    if abs(e1 - e0) > 1e-9 or not ds <= 1e-8 or dg > a.gtol or len(g0) != len(g1):
        raise SystemExit("FAILED: open-shell GPU result differs from the CPU")


def check_opt(a):
    xyz = os.path.join(TC, "structures/small/methyl.pyramidal.xyz")
    with tempfile.TemporaryDirectory() as d:
        rc, out = cli(a.curcuma, ["-opt", xyz, "-method", a.case, "-spin", "1", "-no_bmt", "-verbosity", "1",
                                  "-qm.basis", "def2-SVP", "-threads", str(a.threads)], d)
        res = [os.path.join(d, f) for f in os.listdir(d) if f.endswith(".opt.xyz")]
        atoms = read_xyz(res[0]) if res else []
    # optimiser table: step, energy (Eh), energy change, RMSD change, gradient norm (Eh/Bohr), time
    rows = [(float(m.group(2)), float(m.group(5))) for m in
            re.finditer(r"^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+\d+\.\d+\s*$", out, re.M)]
    # CH3 radical minimum (UHF/def2-SVP) is planar: distance of C from the plane of the three H
    oop = None
    if len(atoms) == 4:
        import math
        h = [a_[1:] for a_ in atoms[1:]]
        u = [h[1][k] - h[0][k] for k in range(3)]
        v = [h[2][k] - h[0][k] for k in range(3)]
        n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        nn = math.sqrt(sum(x * x for x in n))
        oop = abs(sum((atoms[0][1 + k] - h[0][k]) * n[k] for k in range(3))) / nn
    print("%s -opt CH3 (distorted, -spin 1): exit %d, %d steps, E %s -> %s Eh, last |g| %s Eh/Bohr, C out of plane %s A" %
          (a.case, rc, len(rows), rows[0][0] if rows else None, rows[-1][0] if rows else None,
           rows[-1][1] if rows else None, oop))
    if rc != 0 or len(rows) < 2 or not rows[-1][0] < rows[0][0] - 1e-3 or rows[-1][1] > 1e-4 or oop is None or oop > 1e-2:
        print(out[-4000:])
        raise SystemExit("FAILED: open-shell optimisation")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", required=True, choices=["kernel", "ref", "fd", "rhf", "direct", "reject", "gpu", "opt"])
    ap.add_argument("--case", default="")
    ap.add_argument("--dump")
    ap.add_argument("--exe")
    ap.add_argument("--curcuma")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--gtol", type=float, default=1e-8)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--coords", nargs="*", default=["0:x", "1:y", "2:z"])
    ap.add_argument("--reference", default=os.path.join(HERE, "reference.json"))
    ap.add_argument("--gpu", action="store_true", help="ref/fd: run the CUDA path (stage 6a)")
    ap.add_argument("--plugin-dir", default="", help="directory of libcurcuma_cuda.so (default: that of --curcuma)")
    ap.add_argument("--extra-qm", nargs="*", default=[], help="gpu: extra -qm settings KEY=VALUE")
    a = ap.parse_args()
    if a.gpu or a.check == "gpu":
        gpu_available(a)
        GPU["on"] = a.gpu
    {"kernel": check_kernel, "ref": check_ref, "fd": check_fd, "rhf": check_rhf, "direct": check_direct,
     "reject": check_reject, "gpu": check_gpu, "opt": check_opt}[a.check](a)


if __name__ == "__main__":
    main()
