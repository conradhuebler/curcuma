#!/usr/bin/env python3
"""Configuration of a methyl hexopyranoside (C7H14O6) found as a fragment in an xyz file.

Claude Generated (Oct 2026). Answers: which sugar (gluco, galacto, manno, ...), D or L, alpha or beta, and the ring
conformation, from the 3D coordinates alone.

Method (conformation independent, Haworth convention):
  1. Bonds from covalent radii; fragments; the fragment with formula C7H14O6 is the guest.
  2. The pyranose ring is the 6-ring with one O and five C. C1 is the ring carbon bonded to the ring O and to an
     exocyclic O, C5 the other ring carbon bonded to the ring O (exocyclic substituent CH2OH), numbering C1..C5 away
     from the ring O.
  3. The ring normal is the best-fit plane normal, oriented so that the C5 substituent (CH2OH) points up.
  4. Each heavy substituent of C1..C4 is up (same face as CH2OH) or down. In the D series seen from above the
     numbering runs clockwise; the sign of sum(r_i x r_i+1) along the normal gives D (clockwise) or L.
  5. With CH2OH up, the faces (C2, C3, C4) identify the sugar (D-aldohexoses, Haworth): gluco (down, up, down),
     galacto (down, up, up), manno (up, up, down), allo (down, down, down), altro (up, down, down),
     gulo (down, down, up), ido (up, down, up), talo (up, up, up). The anomer is beta if the C1 substituent is on the
     CH2OH face, else alpha (mirrored for the L series).
  6. Axial or equatorial is the angle of each ring-substituent bond to the normal (axial above 45 degrees from the plane).

Usage:  pyranoside_stereo.py FILE.xyz [--frames]   (prints one line per frame of the guest)
Standard library plus numpy.
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import structlib as sl  # noqa: E402

RADII = {"H": .31, "C": .76, "N": .71, "O": .66, "S": 1.05, "F": .57, "Cl": 1.02, "Br": 1.2, "I": 1.39}
SUGARS = {("down", "up", "down"): "gluco", ("down", "up", "up"): "galacto", ("up", "up", "down"): "manno",
          ("down", "down", "down"): "allo", ("up", "down", "down"): "altro", ("down", "down", "up"): "gulo",
          ("up", "down", "up"): "ido", ("up", "up", "up"): "talo"}


def bonds(sym, xyz):
    n = len(sym)
    nb = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if sl.math.dist(xyz[i], xyz[j]) < 1.25 * (RADII.get(sym[i], 1.4) + RADII.get(sym[j], 1.4)) + 0.1:
                nb[i].append(j)
                nb[j].append(i)
    return nb


def fragments(nb):
    seen, out = set(), []
    for s in range(len(nb)):
        if s in seen:
            continue
        stack, comp = [s], []
        while stack:
            a = stack.pop()
            if a in seen:
                continue
            seen.add(a)
            comp.append(a)
            stack.extend(nb[a])
        out.append(sorted(comp))
    return out


def find_ring(sym, nb, atoms):
    """Return [ring O, then the five ring carbons in path order] for the 6-ring with one O inside `atoms`."""
    aset = set(atoms)
    for o in [a for a in atoms if sym[a] == "O" and len([x for x in nb[a] if sym[x] == "C"]) == 2]:
        ca, cb = [x for x in nb[o] if sym[x] == "C"]

        def dfs(path):
            if len(path) == 5:
                return path if path[-1] == cb else None
            for x in nb[path[-1]]:
                if x in aset and sym[x] == "C" and x not in path:
                    r = dfs(path + [x])
                    if r:
                        return r
            return None
        r = dfs([ca])
        if r:
            return [o] + r
    return None


def analyse(sym, xyz):
    nb = bonds(sym, xyz)
    out = []
    for frag in fragments(nb):
        if sl.hill_formula([sym[i] for i in frag]) != "C7H14O6":
            continue
        ring = find_ring(sym, nb, frag)
        if not ring:
            out.append({"error": "no pyranose ring"})
            continue
        o5, a, b, c, d, e = ring     # a..e in path order, which end is C1 is decided below
        path = [a, b, c, d, e]

        def exo_o(ci):
            return [x for x in nb[ci] if sym[x] == "O" and x not in ring]

        def exo_c(ci):
            return [x for x in nb[ci] if sym[x] == "C" and x not in ring]
        # C1 carries the exocyclic O (acetal), C5 the exocyclic C
        if exo_o(path[0]) and exo_c(path[-1]):
            c1, c2, c3, c4, c5 = path
        elif exo_o(path[-1]) and exo_c(path[0]):
            c5, c4, c3, c2, c1 = path
        else:
            out.append({"error": "could not assign C1/C5"})
            continue
        P = np.array(xyz)
        rp = P[[c1, c2, c3, c4, c5, o5]]
        cen = rp.mean(axis=0)
        u, s, vt = np.linalg.svd(rp - cen)
        n = vt[2]
        sub = {1: exo_o(c1)[0], 2: exo_o(c2)[0], 3: exo_o(c3)[0], 4: exo_o(c4)[0], 5: exo_c(c5)[0]}
        ch2oh = P[sub[5]] - P[c5]
        if np.dot(ch2oh, n) < 0:
            n = -n
        # rotation sense of C1->C2->C3->C4->C5->O5 seen from +n
        seq = P[[c1, c2, c3, c4, c5, o5]] - cen
        m = sum(np.cross(seq[i], seq[(i + 1) % 6]) for i in range(6))
        series = "D" if np.dot(m, n) < 0 else "L"      # clockwise from the CH2OH side = D
        face = {}
        axial = {}
        ringc = {1: c1, 2: c2, 3: c3, 4: c4, 5: c5}
        for k in range(1, 6):
            v = P[sub[k]] - P[ringc[k]]
            v = v / np.linalg.norm(v)
            face[k] = "up" if np.dot(v, n) > 0 else "down"
            axial[k] = "ax" if abs(np.dot(v, n)) > 0.7071 else "eq"
        sugar = SUGARS.get((face[2], face[3], face[4]), "?")
        anomer = "beta" if face[1] == face[5] else "alpha"
        out.append({"series": series, "sugar": sugar, "anomer": anomer,
                    "faces": "".join("U" if face[k] == "up" else "d" for k in range(1, 6)),
                    "axial_eq": "".join("a" if axial[k] == "ax" else "e" for k in range(1, 6)),
                    "atoms": [int(i) for i in frag]})
    return out


def analyse_file(path):
    """Analyse every frame of an xyz file. Returns a list (one entry per frame) of lists of guest results."""
    with open(path, encoding="utf-8") as f:
        lines = f.read().split("\n")
    i, results = 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        sym, xyz = [], []
        for k in range(n):
            p = lines[i + 2 + k].split()
            sym.append(p[0])
            xyz.append((float(p[1]), float(p[2]), float(p[3])))
        results.append(analyse(sym, xyz))
        i += 2 + n
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--frames", action="store_true", help="one line per frame instead of a summary")
    ap.add_argument("--dump-guest", metavar="FILE", help="write the guest fragment of the first frame as xyz")
    args = ap.parse_args()
    with open(args.file, encoding="utf-8") as f:
        lines = f.read().split("\n")
    i, frame, results = 0, 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        sym, xyz = [], []
        for k in range(n):
            p = lines[i + 2 + k].split()
            sym.append(p[0])
            xyz.append((float(p[1]), float(p[2]), float(p[3])))
        res = analyse(sym, xyz)
        results.append(res)
        if args.dump_guest and frame == 0 and res and "atoms" in res[0]:
            with open(args.dump_guest, "w", encoding="utf-8") as g:
                idx = res[0]["atoms"]
                g.write(f"{len(idx)}\nguest fragment of {os.path.basename(args.file)}\n")
                for k in idx:
                    g.write(f"{sym[k]} {xyz[k][0]:.8f} {xyz[k][1]:.8f} {xyz[k][2]:.8f}\n")
        if args.frames:
            for r in res:
                print(frame, {k: v for k, v in r.items() if k != "atoms"})
        frame += 1
        i += 2 + n
    from collections import Counter
    c = Counter()
    for res in results:
        for r in res:
            c[tuple(sorted((k, v) for k, v in r.items() if k in ("series", "sugar", "anomer", "faces", "axial_eq", "error")))] += 1
    print(f"{os.path.basename(args.file)}: {len(results)} frame(s), guest fragments found {sum(len(r) for r in results)}")
    for k, v in c.most_common():
        print("  ", v, "x", dict(k))


if __name__ == "__main__":
    main()
