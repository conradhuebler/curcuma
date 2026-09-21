#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""Arbitrate a curcuma-vs-xtb gradient disagreement with finite differences.

`gradient_compare.py` reports WHERE the two analytic gradients differ; it cannot say
WHICH one is right. This does, using the only referee that needs no third code: the
central finite difference of each program's OWN total energy. A correct analytic
gradient must reproduce the finite difference of the energy it belongs to.

For the component with the largest |g_curcuma - g_xtb| it prints
    curcuma analytic | xtb analytic | FD of curcuma energy | FD of xtb energy
all in Eh/Bohr, at two step sizes.

Usage:
    python scripts/gradient_arbitrate.py <xyz> <method> [charge] [uhf]
"""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CURCUMA = str(ROOT / "release" / "curcuma")
XTB = "/opt/bin/xtb"
AU = 0.52917721092
XTB_FLAG = {"gfn1": ["--gfn", "1"], "gfn2": ["--gfn", "2"], "gfnff": ["--gfnff"]}


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        els.append(p[0])
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return els, xyz


def write_xyz(path, els, xyz):
    with open(path, "w") as f:
        f.write(f"{len(els)}\n\n")
        for e, r in zip(els, xyz):
            f.write(f"{e} {r[0]:.12f} {r[1]:.12f} {r[2]:.12f}\n")


def clean(w):
    for s in ("gradient", "energy", "xtbrestart", "charges", "wbo", "xtbtopo.mol"):
        p = os.path.join(w, s)
        if os.path.exists(p):
            os.remove(p)


def e_cur(path, method, chrg, uhf, w):
    cmd = [CURCUMA, "-sp", str(path), "-method", method, "-no_bmt", "-verbosity", "0"]
    if chrg:
        cmd += ["-charge", str(chrg)]
    if uhf and method != "gfnff":
        cmd += ["-spin", str(uhf)]
    o = subprocess.run(cmd, capture_output=True, text=True, cwd=w).stdout
    m = re.search(r"Single Point Energy = (-?\d+\.\d+)", o)
    return float(m.group(1)) if m else None


def e_xtb(path, method, chrg, uhf, w):
    clean(w)
    cmd = [XTB, str(path)] + XTB_FLAG[method] + ["--sp", "--acc", "0.0001"]
    if chrg:
        cmd += ["--chrg", str(chrg)]
    if uhf:
        cmd += ["--uhf", str(uhf)]
    o = subprocess.run(cmd, capture_output=True, text=True, cwd=w).stdout
    m = re.search(r":: total energy\s+(-?\d+\.\d+)", o)
    return float(m.group(1)) if m else None


def g_cur(path, method, chrg, uhf, w):
    d = os.path.join(w, "c.grad")
    if os.path.exists(d):
        os.remove(d)
    cmd = [CURCUMA, "-sp", str(path), "-method", method, "-gradient",
           "-dump_gradient", d, "-no_bmt", "-verbosity", "0"]
    if chrg:
        cmd += ["-charge", str(chrg)]
    if uhf and method != "gfnff":
        cmd += ["-spin", str(uhf)]
    subprocess.run(cmd, capture_output=True, text=True, cwd=w)
    if not os.path.exists(d):
        return None
    return [[float(v) * AU for v in l.split()]
            for l in open(d) if not l.startswith("#") and len(l.split()) == 3]


def g_xtb(path, method, chrg, uhf, w):
    clean(w)
    cmd = [XTB, str(path)] + XTB_FLAG[method] + ["--grad", "--acc", "0.0001"]
    if chrg:
        cmd += ["--chrg", str(chrg)]
    if uhf:
        cmd += ["--uhf", str(uhf)]
    subprocess.run(cmd, capture_output=True, text=True, cwd=w)
    gf = os.path.join(w, "gradient")
    if not os.path.exists(gf):
        return None
    return [[float(v.replace("D", "E")) for v in l.split()]
            for l in open(gf) if len(l.split()) == 3]


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: gradient_arbitrate.py <xyz> <method> [charge] [uhf]")
    xyz_in = Path(sys.argv[1]).resolve()
    method = sys.argv[2]
    chrg = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    uhf = int(sys.argv[4]) if len(sys.argv) > 4 else 0

    w = tempfile.mkdtemp(prefix="arb_")
    els, xyz = read_xyz(xyz_in)
    gc = g_cur(xyz_in, method, chrg, uhf, w)
    gx = g_xtb(xyz_in, method, chrg, uhf, w)
    if gc is None or gx is None or len(gc) != len(gx):
        sys.exit("could not obtain both analytic gradients")

    diffs = sorted(((abs(gc[i][k] - gx[i][k]), i, k)
                    for i in range(len(gc)) for k in range(3)), reverse=True)
    dmax, i, k = diffs[0]
    print(f"{xyz_in.name}  method={method}  charge={chrg} uhf={uhf}  nat={len(gc)}")
    print(f"largest disagreement: atom {i} ({els[i]}) {'xyz'[k]}   |dG| = {dmax:.3e} Eh/Bohr")

    tmp = os.path.join(w, "t.xyz")
    print(f"{'h/Ang':>7s} {'FD(curcuma E)':>15s} {'FD(xtb E)':>15s}")
    for h in (0.002, 0.005):
        plus = [r[:] for r in xyz]
        plus[i][k] += h
        write_xyz(tmp, els, plus)
        cp, xp = e_cur(tmp, method, chrg, uhf, w), e_xtb(tmp, method, chrg, uhf, w)
        minus = [r[:] for r in xyz]
        minus[i][k] -= h
        write_xyz(tmp, els, minus)
        cm, xm = e_cur(tmp, method, chrg, uhf, w), e_xtb(tmp, method, chrg, uhf, w)
        fc = (cp - cm) / (2 * h) * AU if None not in (cp, cm) else float("nan")
        fx = (xp - xm) / (2 * h) * AU if None not in (xp, xm) else float("nan")
        print(f"{h:7} {fc:15.6f} {fx:15.6f}")
    print(f"{'analytic':>7s} {gc[i][k]:15.6f} {gx[i][k]:15.6f}   <- curcuma | xtb")


if __name__ == "__main__":
    main()
