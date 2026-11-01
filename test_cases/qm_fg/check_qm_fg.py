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


def run_dump(dump, xyz, method, qm, gradient=False, threads=1):
    cmd = [dump, xyz, method] + qm_args(qm) + ["--threads", str(threads)]
    if gradient:
        cmd.append("--gradient")
    p = subprocess.run(cmd, capture_output=True, text=True, env=env_embedded())
    out = p.stdout + p.stderr
    m = re.search(r"^energy (\S+)$", out, re.M)
    if p.returncode != 0 or not m:
        print(out[-4000:])
        raise SystemExit("dump_qm_energy failed: " + " ".join(cmd))
    g = [[float(x) for x in l.split()[2:5]] for l in out.splitlines() if l.startswith("grad ")]
    return float(m.group(1)), g


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
    qm = dict(ref["qm"])
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
    if not ok:
        raise SystemExit("FAILED: %s outside the tolerance" % a.case)


def check_fd(a):
    ref = json.load(open(a.reference))["cases"][a.case]
    xyz = os.path.join(TC, ref["xyz"])
    qm = dict(ref["qm"])
    qm["scf_threshold"] = 1e-11 if "QZV" not in qm["basis"] else 1e-10
    atoms = read_xyz(xyz)
    _, g = run_dump(a.dump, xyz, ref["method"], qm, True, a.threads)
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
                es[s], _ = run_dump(a.dump, p, ref["method"], qm, False, a.threads)
            fd = (es[-2] - 8 * es[-1] + 8 * es[1] - es[2]) / (12 * h)
            diff = g[at][k] - fd
            worst = max(worst, abs(diff))
            print("%-30s atom %d %s analytic %+.10f FD %+.10f diff %+.2e" % (a.case, at, ax, g[at][k], fd, diff))
    print("max |analytic - FD| = %.2e Eh/Bohr (limit %.1e)" % (worst, a.gtol))
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", required=True, choices=["ref", "fd", "gpu"])
    ap.add_argument("--case", required=True)
    ap.add_argument("--dump")
    ap.add_argument("--curcuma")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--gtol", type=float, default=None)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--coords", nargs="*", default=["0:z", "1:x", "1:z"])
    ap.add_argument("--reference", default=os.path.join(HERE, "reference.json"),
                    help="reference file (default: qm_fg/reference.json; test_cases/qm_ecp passes its own)")
    a = ap.parse_args()
    if a.gtol is None:
        a.gtol = 1e-8
    {"ref": check_ref, "fd": check_fd, "gpu": check_gpu}[a.check](a)


if __name__ == "__main__":
    main()
