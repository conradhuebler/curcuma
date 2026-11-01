#!/usr/bin/env python3
"""Native HF-3c beyond Ne against ORCA 6.1 (ctest -L qm_hf3c_z). Standard library only.

reference_orca_z.json (scripts/hf3c_orca_reference.py): ORCA 6.1.0 `! HF-3c VeryTightSCF` per case: E(HF/MINIX),
E(D3BJ), E(gCP+bas) (= gCP + SRB, 9 printed decimals), the total, and for some cases the EnGrad gradient (Eh/Bohr).

Checks (--check):
  ref   curcuma -sp -method hf-3c: every term and the total against ORCA, |dE| <= --tol (default 1e-7 Eh; the
        gCP+SRB term is compared at 1e-7 too, ORCA prints it with 9 decimals); with a stored gradient also
        -dump_gradient against ORCA's EnGrad, max |dg| <= --gtol (default 1e-6 Eh/Bohr)
  fd    FD gate of the fully assembled hf-3c gradient: -dump_gradient against a 4-point central difference of
        curcuma's own -sp energies (step 0.005 Bohr) on --coords (atom:axis), max |dg| <= --gtol (default 1e-8)
  opt   curcuma -opt -method hf-3c from the library geometry: the final energy against ORCA's `! HF-3c Opt
        TightOpt` minimum from the same start (|dE| <= --tol, default 1e-7 Eh; curcuma L-BFGS, preset verytight)

Claude Generated (Nov 2026). AI-generated, machine-tested.
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.normpath(os.path.join(HERE, ".."))
BOHR = 0.529177210903  # curcuma units.h
PATTERN = {"hf": r"E\(HF/MINIX\)", "d3": r"E\(D3BJ\)", "gcp": r"E\(gCP\)", "srb": r"E\(SRB\)", "total": r"E\(HF-3c\)"}
SCF = ["-qm.scf_threshold", "1e-9"]


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    return [(t.split()[0], [float(v) for v in t.split()[1:4]]) for t in lines[2:2 + n]]


def write_xyz(path, atoms):
    with open(path, "w") as f:
        f.write("%d\n\n" % len(atoms))
        for s, x in atoms:
            f.write("%s %.12f %.12f %.12f\n" % (s, x[0], x[1], x[2]))


def run_sp(curcuma, xyz, tmp, gradient=False):
    cmd = [curcuma, "-sp", xyz, "-method", "hf-3c", "-verbosity", "1", "-no_bmt"] + SCF
    gfile = os.path.join(tmp, "grad.txt")
    if gradient:
        cmd += ["-dump_gradient", gfile]
    p = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, timeout=1200)
    got = {}
    for key, pat in PATTERN.items():
        m = re.findall(pat + r"\s*=\s*(-?\d+\.\d+)", p.stdout)
        if not m:
            print(p.stdout[-3000:], p.stderr[-2000:])
            sys.exit("no '%s' line in the curcuma output (exit %d)" % (key, p.returncode))
        got[key] = float(m[-1])
    if gradient:  # -dump_gradient writes Eh/Angstrom and the energy with 15 decimals in its header
        lines = open(gfile).read().splitlines()
        got["gradient"] = [[float(v) * BOHR for v in l.split()] for l in lines if l.strip() and not l.startswith("#")]
        m = re.search(r"# energy\s+(-?\d+\.\d+)", "\n".join(lines))
        if m:
            got["total_full"] = float(m.group(1))
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curcuma", required=True)
    ap.add_argument("--ref", default=os.path.join(HERE, "reference_orca_z.json"))
    ap.add_argument("--case", required=True)
    ap.add_argument("--check", choices=["ref", "fd", "opt"], default="ref")
    ap.add_argument("--tol", type=float, default=1e-7)
    ap.add_argument("--gtol", type=float, default=None)
    ap.add_argument("--coords", nargs="*", default=["0:z", "1:x", "1:z"])
    a = ap.parse_args()
    curcuma = os.path.abspath(a.curcuma)
    ref = json.load(open(a.ref))
    fail = False
    with tempfile.TemporaryDirectory() as tmp:
        if a.check == "opt":
            r = ref["opt"][a.case]
            xyz = os.path.join(tmp, "start.xyz")
            write_xyz(xyz, read_xyz(os.path.join(TC, r["xyz"])))
            p = subprocess.run([curcuma, "-opt", xyz, "-method", "hf-3c", "-verbosity", "1", "-no_bmt"] + SCF
                               + ["-opt.optimizer", "lbfgs", "-opt.convergence_preset", "verytight"],
                               cwd=tmp, capture_output=True, text=True, timeout=1800)
            outs = [f for f in glob.glob(os.path.join(tmp, "*.opt.xyz")) + glob.glob(os.path.join(tmp, "*opt*.xyz"))]
            if p.returncode != 0 or not outs:
                print(p.stdout[-3000:])
                sys.exit("opt failed (exit %d)" % p.returncode)
            got = run_sp(curcuma, sorted(outs, key=os.path.getmtime)[-1], tmp)
            d = got["total"] - r["total"]
            ok = abs(d) <= a.tol
            print("%-6s opt  curcuma %.12f  ORCA %.12f  diff %+.2e  %s" % (a.case, got["total"], r["total"], d,
                                                                         "ok" if ok else "FAIL (tol %.0e)" % a.tol))
            sys.exit(0 if ok else 1)

        r = ref["cases"][a.case]
        xyz = os.path.join(TC, r["xyz"])
        if a.check == "ref":
            gtol = a.gtol if a.gtol is not None else 1e-6
            got = run_sp(curcuma, xyz, tmp, gradient="gradient" in r)
            got["gcp_srb"] = got["gcp"] + got["srb"]
            for key in ("hf", "d3", "gcp_srb", "total"):
                d = got[key] - r[key]
                ok = abs(d) <= a.tol
                fail |= not ok
                print("%-6s %-7s curcuma %18.12f  ORCA %18.12f  diff %+.2e  %s" % (
                    a.case, key, got[key], r[key], d, "ok" if ok else "FAIL (tol %.0e)" % a.tol))
            if "gradient" in r:
                dg = max(abs(x - y) for gr, go in zip(got["gradient"], r["gradient"]) for x, y in zip(gr, go))
                ok = dg <= gtol
                fail |= not ok
                print("%-6s gradient max |dg| vs ORCA EnGrad %.2e Eh/Bohr  %s" % (a.case, dg, "ok" if ok else
                                                                                   "FAIL (tol %.0e)" % gtol))
        else:  # fd
            gtol = a.gtol if a.gtol is not None else 1e-8
            atoms = read_xyz(xyz)
            g = run_sp(curcuma, xyz, tmp, gradient=True)["gradient"]
            h = 0.005  # Bohr
            for c in a.coords:
                i, ax = c.split(":")
                i, k = int(i), "xyz".index(ax)
                e = {}
                for s in (-2, -1, 1, 2):
                    moved = [(sym, list(x)) for sym, x in atoms]
                    moved[i][1][k] += s * h * BOHR
                    f = os.path.join(tmp, "d.xyz")
                    write_xyz(f, moved)
                    e[s] = run_sp(curcuma, f, tmp, gradient=True)["total_full"]  # 15-decimal energy
                fd = (8.0 * (e[1] - e[-1]) - (e[2] - e[-2])) / (12.0 * h)
                d = g[i][k] - fd
                ok = abs(d) <= gtol
                fail |= not ok
                print("%-6s FD %s analytic %+.10f  FD %+.10f  diff %+.2e  %s" % (a.case, c, g[i][k], fd, d,
                                                                                "ok" if ok else "FAIL (tol %.0e)" % gtol))
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
