#!/usr/bin/env python3
"""Idealised test geometries for the QM wall-time comparison (Claude Generated, Oct 2026).

Generates the structures in code so that no structure file lives outside
test_cases/structures/ (repository rule). All geometries are idealised, not optimised:
the timing comparison only needs the same geometry in every program.

  alkaneN     all-trans n-alkane CnH(2n+2), r(CC)=1.53 A, r(CH)=1.09 A, tetrahedral angles
  benzene     test_cases/structures/small/benzene.ideal-d6h.xyz (copied)
  anthracene  idealised D2h, r(CC)=1.40 A, r(CH)=1.09 A
  waterN      N water molecules (r(OH)=0.9572 A, 104.52 deg) on a cubic 2.9 A grid,
              orientation alternating by index

Usage: geometries.py OUTDIR name [name ...]   e.g. alkane20 benzene water4
"""
import math
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _norm(v):
    n = math.sqrt(sum(x * x for x in v))
    return [x / n for x in v]


def _add(a, b, s=1.0):
    return [a[i] + s * b[i] for i in range(3)]


def alkane(n, rcc=1.53, rch=1.09):
    theta = math.radians(109.4712)
    dx, dy = rcc * math.sin(theta / 2), rcc * math.cos(theta / 2)
    # carbons on a zigzag, plus one virtual carbon at each end for the terminal H
    pos = [[i * dx, (i % 2) * dy, 0.0] for i in range(-1, n + 1)]
    atoms = [("C", p) for p in pos[1:-1]]
    alpha = math.radians(109.4712 / 2)
    z = [0.0, 0.0, 1.0]
    hs = []
    for k in range(1, n + 1):
        c, a, b = pos[k], pos[k - 1], pos[k + 1]
        u1, u2 = _norm(_add(a, c, -1)), _norm(_add(b, c, -1))
        bis = _norm([-(u1[i] + u2[i]) for i in range(3)])
        for s in (+1, -1):
            d = [bis[i] * math.cos(alpha) + s * z[i] * math.sin(alpha) for i in range(3)]
            hs.append(("H", _add(c, d, rch)))
        if k == 1:
            hs.append(("H", _add(c, u1, rch)))
        if k == n:
            hs.append(("H", _add(c, u2, rch)))
    return atoms + hs


def anthracene(rcc=1.40, rch=1.09):
    # three hexagons with vertical shared edges; centres at x = 0, +-sqrt(3) rcc
    cs = []
    for cx in (-math.sqrt(3) * rcc, 0.0, math.sqrt(3) * rcc):
        for k in range(6):
            ang = math.radians(30 + 60 * k)
            p = [cx + rcc * math.cos(ang), rcc * math.sin(ang), 0.0]
            if not any(math.dist(p, q) < 0.1 for q in cs):
                cs.append(p)
    atoms = [("C", p) for p in cs]
    for p in cs:
        nb = [q for q in cs if 0.1 < math.dist(p, q) < rcc + 0.1]
        if len(nb) == 2:
            d = _norm([-sum(_norm(_add(q, p, -1))[i] for q in nb) for i in range(3)])
            atoms.append(("H", _add(p, d, rch)))
    return atoms


def water_cluster(n, a=2.9, roh=0.9572, ang=104.52):
    h = math.radians(ang / 2)
    m = math.ceil(n ** (1.0 / 3.0))
    atoms = []
    for k in range(n):
        ix, iy, iz = k % m, (k // m) % m, k // (m * m)
        o = [ix * a, iy * a, iz * a]
        # alternate the molecular plane between xy, yz, xz
        e1, e2 = [([1, 0, 0], [0, 1, 0]), ([0, 1, 0], [0, 0, 1]), ([0, 0, 1], [1, 0, 0])][k % 3]
        atoms.append(("O", o))
        for s in (+1, -1):
            d = [e1[i] * math.cos(h) * (-1) ** k + s * e2[i] * math.sin(h) for i in range(3)]
            atoms.append(("H", _add(o, d, roh)))
    return atoms


def read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0])
    out = []
    for ln in lines[2:2 + n]:
        p = ln.split()
        out.append((p[0], [float(x) for x in p[1:4]]))
    return out


def build(name):
    if name.startswith("alkane"):
        return alkane(int(name[6:]))
    if name == "benzene":
        return read_xyz(os.path.join(REPO, "test_cases/structures/small/benzene.ideal-d6h.xyz"))
    if name == "anthracene":
        return anthracene()
    if name.startswith("water"):
        return water_cluster(int(name[5:]))
    raise SystemExit(f"unknown structure {name}")


def write(outdir, name):
    atoms = build(name)
    # sanity: no two atoms closer than 0.9 A
    for i in range(len(atoms)):
        for j in range(i):
            if math.dist(atoms[i][1], atoms[j][1]) < 0.9:
                raise SystemExit(f"{name}: atoms {i},{j} too close")
    path = os.path.join(outdir, name + ".xyz")
    with open(path, "w") as f:
        f.write(f"{len(atoms)}\n{name} (scripts/qm_timing/geometries.py)\n")
        for el, p in atoms:
            f.write(f"{el} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}\n")
    return path


if __name__ == "__main__":
    os.makedirs(sys.argv[1], exist_ok=True)
    for nm in sys.argv[2:]:
        print(write(sys.argv[1], nm))
