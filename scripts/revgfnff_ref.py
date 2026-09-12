#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP2
"""r2SCAN-3c reference campaign for rev-gfnff (ORCA, rigid scans and single points).

Every result lands immediately in test_cases/revgfnff/ref/<class>/<system>/ as
energies.json + gradients.json + points.xyz + meta.json, so the fitter can use the
campaign while it is still running. Jobs that already have a complete energies.json
are skipped, which makes the campaign restartable and extendable.

Classes (roadmap WP2):
  A  dissociation curves: one bond of a minimal molecule stretched rigidly (the fragment
     on one side moves along the bond axis), 20 points from 0.75 to 3.5 x r_eq; RKS from
     the inside out and UKS broken-symmetry from the outside in (both kept)
  C  hyper-coordination: a radical (H, CH3) or a proton approaches a saturated centre,
     15 points from 1.0 to 3.0 A (UKS doublets, RKS for the cations)
  D  off-equilibrium: GFN-FF MD snapshots at 1000 / 2000 K, single points
  E  charged cases for the charge model: Cl2-/F2- curves, proton-transfer transits,
     with Hirshfeld charges
  B  NEB-TS paths (not yet in this driver)

Reference geometries: the hand-written start geometry is relaxed with curcuma GFN2 and then
with r2SCAN-3c (ORCA Opt); the optimised molecule is the origin of every rigid scan.

Units in the files: energies Eh, gradients Eh/Bohr (raw ORCA, unit recorded in the file),
coordinates Angstrom.

Usage:
    python scripts/revgfnff_ref.py plan
    python scripts/revgfnff_ref.py run --classes A C --jobs 6 --nprocs 4
    python scripts/revgfnff_ref.py run --classes A --only h2 hf --jobs 2
    python scripts/revgfnff_ref.py status
"""
import argparse
import datetime
import gzip
import itertools
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
CURCUMA = REPO / "release" / "curcuma"
REF = REPO / "test_cases" / "revgfnff" / "ref"
GMTKN = REPO / "test_cases" / "GMTKN55-testset"
ORCA = os.environ.get("ORCA_PATH", "/opt/orca/orca")
KEYWORDS = "r2SCAN-3c TightSCF"
BOHR = 0.529177210903

# ------------------------------------------------------------------ molecules (start geometries)
# name: (charge, mult, [(sym, x, y, z), ...])  -- rough, relaxed by GFN2 and r2SCAN-3c before use
MOL = {
    "h2":     (0, 1, [("H", 0, 0, 0), ("H", 0, 0, 0.74)]),
    "hf":     (0, 1, [("H", 0, 0, 0), ("F", 0, 0, 0.92)]),
    "hcl":    (0, 1, [("H", 0, 0, 0), ("Cl", 0, 0, 1.27)]),
    "f2":     (0, 1, [("F", 0, 0, 0), ("F", 0, 0, 1.41)]),
    "cl2":    (0, 1, [("Cl", 0, 0, 0), ("Cl", 0, 0, 1.99)]),
    "n2":     (0, 1, [("N", 0, 0, 0), ("N", 0, 0, 1.10)]),
    "co":     (0, 1, [("C", 0, 0, 0), ("O", 0, 0, 1.13)]),
    "o2":     (0, 3, [("O", 0, 0, 0), ("O", 0, 0, 1.21)]),
    "h2o":    (0, 1, [("O", 0, 0, 0), ("H", 0.76, 0.59, 0), ("H", -0.76, 0.59, 0)]),
    "nh3":    (0, 1, [("N", 0, 0, 0), ("H", 0.94, 0, -0.33), ("H", -0.47, 0.81, -0.33), ("H", -0.47, -0.81, -0.33)]),
    "ch4":    (0, 1, [("C", 0, 0, 0), ("H", 0.63, 0.63, 0.63), ("H", -0.63, -0.63, 0.63), ("H", -0.63, 0.63, -0.63), ("H", 0.63, -0.63, -0.63)]),
    "c2h6":   (0, 1, [("C", 0, 0, 0.77), ("C", 0, 0, -0.77), ("H", 1.02, 0, 1.16), ("H", -0.51, 0.88, 1.16), ("H", -0.51, -0.88, 1.16),
                      ("H", -1.02, 0, -1.16), ("H", 0.51, 0.88, -1.16), ("H", 0.51, -0.88, -1.16)]),
    "c2h4":   (0, 1, [("C", 0, 0, 0.67), ("C", 0, 0, -0.67), ("H", 0.92, 0, 1.23), ("H", -0.92, 0, 1.23), ("H", 0.92, 0, -1.23), ("H", -0.92, 0, -1.23)]),
    "c2h2":   (0, 1, [("C", 0, 0, 0.60), ("C", 0, 0, -0.60), ("H", 0, 0, 1.66), ("H", 0, 0, -1.66)]),
    "hcn":    (0, 1, [("C", 0, 0, 0), ("N", 0, 0, 1.16), ("H", 0, 0, -1.07)]),
    "h2co":   (0, 1, [("C", 0, 0, 0), ("O", 0, 0, 1.21), ("H", 0.94, 0, -0.54), ("H", -0.94, 0, -0.54)]),
    "ch3oh":  (0, 1, [("C", 0, 0, 0), ("O", 1.43, 0, 0), ("H", 1.75, 0.90, 0), ("H", -0.36, 1.03, 0), ("H", -0.36, -0.51, 0.89), ("H", -0.36, -0.51, -0.89)]),
    "ch3nh2": (0, 1, [("C", 0, 0, 0), ("N", 1.47, 0, 0), ("H", 1.85, 0.80, 0.45), ("H", 1.85, -0.80, 0.45), ("H", -0.36, 1.03, 0), ("H", -0.36, -0.51, 0.89), ("H", -0.36, -0.51, -0.89)]),
    "ch2nh":  (0, 1, [("C", 0, 0, 0), ("N", 1.27, 0, 0), ("H", 1.75, 0.90, 0), ("H", -0.55, 0.93, 0), ("H", -0.55, -0.93, 0)]),
    "n2h4":   (0, 1, [("N", 0, 0, 0.72), ("N", 0, 0, -0.72), ("H", 0.95, 0, 1.10), ("H", -0.50, 0.80, 1.10), ("H", -0.95, 0, -1.10), ("H", 0.50, -0.80, -1.10)]),
    "n2h2":   (0, 1, [("N", 0, 0, 0.62), ("N", 0, 0, -0.62), ("H", 0.95, 0, 1.00), ("H", -0.95, 0, -1.00)]),
    "h2o2":   (0, 1, [("O", 0, 0, 0.73), ("O", 0, 0, -0.73), ("H", 0.95, 0, 1.00), ("H", -0.60, 0.75, -1.00)]),
    "ch3f":   (0, 1, [("C", 0, 0, 0), ("F", 1.38, 0, 0), ("H", -0.36, 1.03, 0), ("H", -0.36, -0.51, 0.89), ("H", -0.36, -0.51, -0.89)]),
    "ch3cl":  (0, 1, [("C", 0, 0, 0), ("Cl", 1.78, 0, 0), ("H", -0.36, 1.03, 0), ("H", -0.36, -0.51, 0.89), ("H", -0.36, -0.51, -0.89)]),
    "nh2oh":  (0, 1, [("N", 0, 0, 0), ("O", 1.45, 0, 0), ("H", 1.80, 0.90, 0), ("H", -0.40, 0.90, 0.40), ("H", -0.40, -0.90, 0.40)]),
    "hocl":   (0, 1, [("O", 0, 0, 0), ("H", 0.97, 0, 0), ("Cl", -0.60, 1.58, 0)]),
    # class C partners / class E species
    "ch3":    (0, 2, [("C", 0, 0, 0), ("H", 1.08, 0, 0), ("H", -0.54, 0.94, 0), ("H", -0.54, -0.94, 0)]),
    "cl2m":   (-1, 2, [("Cl", 0, 0, 0), ("Cl", 0, 0, 2.60)]),
    "f2m":    (-1, 2, [("F", 0, 0, 0), ("F", 0, 0, 1.95)]),
    # class B radical partner (BH76 H-transfer, RKT14)
    "oh":     (0, 2, [("O", 0, 0, 0), ("H", 0, 0, 0.97)]),
    # class L (lost scans): formamidine H2N-CH=NH, planar hand-built starting guess,
    # relaxed like every other MOL entry; atom order C, N(imine), N(amine), H(vinyl),
    # H(imine, the one whose C-N-H angle TODO#1 scans), H(amine) x2
    "formamidine": (0, 1, [("C", 0.0, 0.0, 0.0), ("N", 1.279, 0.0, 0.0), ("N", -0.6995, 1.2117, 0.0),
                           ("H", -0.546, -0.9455, 0.0), ("H", 1.7825, -0.8721, 0.0),
                           ("H", -1.7019, 1.1767, 0.0), ("H", -0.1681, 2.0623, 0.0)]),
}

# class A: (molecule, bond i-j, label); the j side moves
CURVES = [
    ("h2", 0, 1, "H-H"), ("hf", 0, 1, "H-F"), ("hcl", 0, 1, "H-Cl"), ("f2", 0, 1, "F-F"), ("cl2", 0, 1, "Cl-Cl"),
    ("n2", 0, 1, "N#N"), ("co", 0, 1, "C#O"), ("o2", 0, 1, "O=O"),
    ("h2o", 0, 1, "O-H"), ("nh3", 0, 1, "N-H"), ("ch4", 0, 1, "C-H"), ("c2h6", 0, 1, "C-C"), ("c2h4", 0, 1, "C=C"),
    ("c2h2", 0, 1, "C#C"), ("hcn", 0, 1, "C#N"), ("hcn", 0, 2, "HC-H"), ("h2co", 0, 1, "C=O"), ("ch3oh", 0, 1, "C-O"),
    ("ch3oh", 1, 2, "HO-H"), ("ch3nh2", 0, 1, "C-N"), ("ch2nh", 0, 1, "C=N"), ("n2h4", 0, 1, "N-N"), ("n2h2", 0, 1, "N=N"),
    ("h2o2", 0, 1, "O-O"), ("ch3f", 0, 1, "C-F"), ("ch3cl", 0, 1, "C-Cl"), ("nh2oh", 0, 1, "N-O"), ("hocl", 0, 2, "O-Cl"),
]
CURVE_GRID = [0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.80, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50]

# class C: (name, host, attacker, host anchor atom, direction rule, charge, mult, distances)
# direction rule: ("opposite", k): attacker approaches anchor along -(anchor->k) (backside of the anchor-k bond)
#                 ("lonepair",): along the negative mean of the anchor's bond vectors (lone-pair side)
APPROACH_GRID = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.75, 1.9, 2.1, 2.3, 2.5, 2.7, 2.85, 3.0]
APPROACH = [
    ("ch4_H", "ch4", "H", 0, ("opposite", 1), 0, 2),
    ("nh3_H", "nh3", "H", 0, ("lonepair",), 0, 2),
    ("h2o_H", "h2o", "H", 0, ("lonepair",), 0, 2),
    ("hf_H", "hf", "H", 1, ("opposite", 0), 0, 2),
    ("h2_H", "h2", "H", 1, ("opposite", 0), 0, 2),
    ("ch4_CH3", "ch4", "ch3", 0, ("opposite", 1), 0, 2),
    ("n2h4_H", "n2h4", "H", 0, ("lonepair",), 0, 2),
    ("ch4_Hp", "ch4", "H", 0, ("opposite", 1), 1, 1),
    ("nh3_Hp", "nh3", "H", 0, ("lonepair",), 1, 1),
    ("h2o_Hp", "h2o", "H", 0, ("lonepair",), 1, 1),
    ("h2_Hp", "h2", "H", 1, ("opposite", 0), 1, 1),
]

# class D: molecules sampled by GFN-FF MD
OFFEQ = ["ch4", "nh3", "h2o", "c2h6", "c2h4", "ch3oh", "ch3nh2", "h2co", "ch3f", "ch3cl"]
OFFEQ_T = [1000, 2000]
OFFEQ_FRAMES = 25   # per temperature -> 50 per molecule

# class E: charged rigid curves (same grid as A), Hirshfeld charges printed
ECURVES = [("cl2m", 0, 1, "Cl-Cl(-)"), ("f2m", 0, 1, "F-F(-)")]

# Species whose unconstrained r2SCAN-3c Opt (via reference_geometry()) diverges instead of
# finding the bound anion minimum -- skip the Opt for these and scan around the given fixed
# bond length (Angstrom) instead. Claude Generated (Sep 2026).
ECURVES_FIXED_R = {
    "f2m": 1.92,  # F2- radical anion; Opt from the 1.95 A hand-written guess ran away to
                  # F-F > 4e5 A (coordinator-flagged, test_cases/revgfnff/ref/E/f2m_F-F-)
}


# ------------------------------------------------------------------ helpers


def write_xyz(path, atoms, comment=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{len(atoms)}\n{comment}\n" + "".join(f"{s} {x:.6f} {y:.6f} {z:.6f}\n" for s, x, y, z in atoms))


def read_xyz(path):
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    return [(t[0], float(t[1]), float(t[2]), float(t[3])) for t in (ln.split() for ln in lines[2:2 + n])]


def read_xyz_frames(path):
    lines = path.read_text().splitlines()
    frames, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        frames.append([(t[0], float(t[1]), float(t[2]), float(t[3])) for t in (ln.split() for ln in lines[i + 2:i + 2 + n])])
        i += n + 2
    return frames


def git_head():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=REPO).stdout.strip()
    except Exception:
        return "?"


def meta(extra):
    m = {"date": datetime.datetime.now().isoformat(timespec="seconds"), "host": platform.node(),
         "curcuma_head": git_head(), "orca": ORCA, "script": "scripts/revgfnff_ref.py", "argv": sys.argv[1:]}
    m.update(extra)
    return m


def bonded_side(atoms, i, j):
    """Atoms on j's side of bond i-j (BFS over a 1.3*rcov connectivity, i excluded)."""
    rcov = {"H": 0.32, "C": 0.75, "N": 0.71, "O": 0.63, "F": 0.64, "Cl": 0.99}
    n = len(atoms)
    adj = [[] for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            d = math.dist(atoms[a][1:], atoms[b][1:])
            if d < 1.3 * (rcov[atoms[a][0]] + rcov[atoms[b][0]]):
                adj[a].append(b)
                adj[b].append(a)
    side, stack = set(), [j]
    while stack:
        k = stack.pop()
        if k in side or k == i:
            continue
        side.add(k)
        stack.extend(adj[k])
    if i in side:
        raise RuntimeError("bond is in a ring")
    return sorted(side)


def stretched(atoms, i, j, r):
    """Rigidly move the j side so that the i-j distance becomes r."""
    side = bonded_side(atoms, i, j)
    pi, pj = atoms[i][1:], atoms[j][1:]
    d = math.dist(pi, pj)
    u = [(pj[k] - pi[k]) / d for k in range(3)]
    shift = [(r - d) * u[k] for k in range(3)]
    out = []
    for a, (s, x, y, z) in enumerate(atoms):
        if a in side:
            out.append((s, x + shift[0], y + shift[1], z + shift[2]))
        else:
            out.append((s, x, y, z))
    return out


def approach(host, attacker_atoms, anchor, rule, dist):
    """Place the attacker (its first atom) at `dist` from the host anchor along the rule's direction."""
    pa = host[anchor][1:]
    if rule[0] == "opposite":
        pk = host[rule[1]][1:]
        v = [pa[c] - pk[c] for c in range(3)]
    else:  # lone pair: opposite to the mean bond vector
        rcov = {"H": 0.32, "C": 0.75, "N": 0.71, "O": 0.63, "F": 0.64, "Cl": 0.99}
        v = [0.0, 0.0, 0.0]
        for b, (s, x, y, z) in enumerate(host):
            if b == anchor:
                continue
            if math.dist(pa, (x, y, z)) < 1.3 * (rcov[host[anchor][0]] + rcov[s]):
                for c, q in enumerate((x, y, z)):
                    v[c] -= q - pa[c]
        if all(abs(c) < 1e-6 for c in v):  # perfectly symmetric centre (e.g. planar): use z
            v = [0.0, 0.0, 1.0]
    nv = math.sqrt(sum(c * c for c in v))
    u = [c / nv for c in v]
    # orient the attacker: its first atom at the tip, the rest pointing away from the host
    first = attacker_atoms[0][1:]
    placed = []
    for s, x, y, z in attacker_atoms:
        rel = [x - first[0], y - first[1], z - first[2]]
        # attacker's local +x axis is mapped onto u (only used for CH3 -- planar, its x axis is in-plane;
        # rotate so the C3 axis (z) points along u instead)
        if len(attacker_atoms) > 1:
            axis = [0.0, 0.0, 1.0]
            rel = rotate_onto(rel, axis, u)
        placed.append((s, pa[0] + dist * u[0] + rel[0], pa[1] + dist * u[1] + rel[1], pa[2] + dist * u[2] + rel[2]))
    return host + placed


def rotate_onto(vec, a, b):
    """Rotate vec by the rotation that maps unit vector a onto unit vector b."""
    cross = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    dot = sum(a[i] * b[i] for i in range(3))
    s = math.sqrt(sum(c * c for c in cross))
    if s < 1e-9:
        return vec if dot > 0 else [-v for v in vec]
    k = [c / s for c in cross]
    ang = math.atan2(s, dot)
    cs, sn = math.cos(ang), math.sin(ang)
    kxv = [k[1] * vec[2] - k[2] * vec[1], k[2] * vec[0] - k[0] * vec[2], k[0] * vec[1] - k[1] * vec[0]]
    kdv = sum(k[i] * vec[i] for i in range(3))
    return [vec[i] * cs + kxv[i] * sn + k[i] * kdv * (1 - cs) for i in range(3)]


# ------------------------------------------------------------------ geometry helpers (E/L/B, Sep 2026)
# Claude Generated (Sep 2026) - WP2 class E remainder, lost scans (class L), NEB-TS paths (class B)


def rotate_about_axis(vec, axis, angle_rad):
    """Rodrigues rotation of a 3-vector about a UNIT vector axis by angle_rad."""
    cs, sn = math.cos(angle_rad), math.sin(angle_rad)
    kxv = [axis[1] * vec[2] - axis[2] * vec[1], axis[2] * vec[0] - axis[0] * vec[2], axis[0] * vec[1] - axis[1] * vec[0]]
    kdv = sum(axis[i] * vec[i] for i in range(3))
    return [vec[i] * cs + kxv[i] * sn + axis[i] * kdv * (1 - cs) for i in range(3)]


def bend_angle(atoms, i, j, k, target_deg):
    """Rigidly rotate the k-side of atom j (bonded_side(atoms, j, k)) about j so that the
    angle i-j-k becomes target_deg. Used for the CH2- / formamidine / FCH3F- rigid bends."""
    pi, pj, pk = atoms[i][1:], atoms[j][1:], atoms[k][1:]
    v1 = [pi[c] - pj[c] for c in range(3)]
    v2 = [pk[c] - pj[c] for c in range(3)]
    n1, n2 = math.sqrt(sum(c * c for c in v1)), math.sqrt(sum(c * c for c in v2))
    u1, u2 = [c / n1 for c in v1], [c / n2 for c in v2]
    cur = math.degrees(math.acos(max(-1.0, min(1.0, sum(u1[c] * u2[c] for c in range(3))))))
    axis = [u1[1] * u2[2] - u1[2] * u2[1], u1[2] * u2[0] - u1[0] * u2[2], u1[0] * u2[1] - u1[1] * u2[0]]
    naxis = math.sqrt(sum(c * c for c in axis))
    axis = [0.0, 0.0, 1.0] if naxis < 1e-9 else [c / naxis for c in axis]
    delta = math.radians(target_deg - cur)
    side = bonded_side(atoms, j, k)
    out = []
    for a, (s, x, y, z) in enumerate(atoms):
        if a in side:
            rel = rotate_about_axis([x - pj[0], y - pj[1], z - pj[2]], axis, delta)
            out.append((s, pj[0] + rel[0], pj[1] + rel[1], pj[2] + rel[2]))
        else:
            out.append((s, x, y, z))
    return out


def set_distance(p_from, p_move, r):
    """Return a point r Angstrom from p_from, along the current p_from->p_move direction."""
    v = [p_move[c] - p_from[c] for c in range(3)]
    n = math.sqrt(sum(c * c for c in v))
    return [p_from[c] + r * v[c] / n for c in range(3)]


def build_3atom_endpoint(ts_atoms, bonded_pair, bond_len, free_idx, far_dist=4.5):
    """One NEB endpoint for a 3-atom A...H...B abstraction TS: set the `bonded_pair`
    distance to `bond_len` (moving the second atom of the pair along the current direction
    from the first) and push `free_idx` out to far_dist from whichever bonded_pair atom it
    currently sits closest to. Reused for all 3-atom BH76 H-transfer endpoints (class B)."""
    atoms = [list(a) for a in ts_atoms]
    a, b = bonded_pair
    new_b = set_distance(atoms[a][1:], atoms[b][1:], bond_len)
    atoms[b] = [atoms[b][0]] + new_b
    da = math.dist(atoms[free_idx][1:], atoms[a][1:])
    db = math.dist(atoms[free_idx][1:], atoms[b][1:])
    anchor = a if da < db else b
    new_free = set_distance(atoms[anchor][1:], atoms[free_idx][1:], far_dist)
    atoms[free_idx] = [atoms[free_idx][0]] + new_free
    return [tuple(a) for a in atoms]


def stack_along_z(frag_a, frag_b, gap):
    """Concatenate two fragments, translating frag_b along +z so its nearest edge sits
    `gap` Angstrom beyond frag_a's farthest edge. No rotation -- a starting guess for an
    unconstrained Opt, not a relaxed complex."""
    za_max = max(a[3] for a in frag_a)
    zb_min = min(a[3] for a in frag_b)
    dz = (za_max + gap) - zb_min
    return list(frag_a) + [(s, x, y, z + dz) for s, x, y, z in frag_b]


def _kabsch(P, Q):
    """Rotate+translate P (n,3) onto Q (n,3); returns (P_aligned, rmsd)."""
    Pc, Qc = P - P.mean(axis=0), Q - Q.mean(axis=0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    P_aligned = (R @ Pc.T).T + Q.mean(axis=0)
    rmsd = math.sqrt(float(((P_aligned - Q) ** 2).sum()) / len(Q))
    return P_aligned, rmsd


def match_transit_endpoints(reactant, ts):
    """Reorder+rigidly align `reactant` onto `ts`'s atom order (same element multiset),
    brute-forcing same-element permutations (systems here have <=8 atoms) and keeping the
    lowest-RMSD correspondence. Returns (aligned_reactant, rmsd)."""
    elems = [a[0] for a in ts]
    groups = {}
    for idx, e in enumerate(elems):
        groups.setdefault(e, []).append(idx)
    r_by_elem = {}
    for e in groups:
        r_by_elem[e] = [i for i, a in enumerate(reactant) if a[0] == e]
    elem_list = sorted(groups)
    perm_lists = [list(itertools.permutations(r_by_elem[e])) for e in elem_list]
    Q = np.array([a[1:] for a in ts], dtype=float)
    best = None
    for combo in itertools.product(*perm_lists):
        order = [None] * len(ts)
        for e, perm in zip(elem_list, combo):
            for slot, ridx in zip(groups[e], perm):
                order[slot] = ridx
        P = np.array([reactant[i][1:] for i in order], dtype=float)
        P_aligned, rmsd = _kabsch(P, Q)
        if best is None or rmsd < best[0]:
            best = (rmsd, P_aligned)
    rmsd, P_aligned = best
    return [(elems[k], float(P_aligned[k][0]), float(P_aligned[k][1]), float(P_aligned[k][2])) for k in range(len(ts))], rmsd


def transit_points(reactant, ts, fracs):
    """Linear-synchronous-transit points between `reactant` (aligned onto ts's atom order)
    and `ts`, at each fraction in `fracs` (0=reactant, 1=ts; >1 extrapolates past the TS,
    used to build the degenerate mirror 'product' at frac=2 for class B)."""
    aligned, rmsd = match_transit_endpoints(reactant, ts)
    pts = []
    for f in fracs:
        pts.append([(aligned[k][0],
                     aligned[k][1] + f * (ts[k][1] - aligned[k][1]),
                     aligned[k][2] + f * (ts[k][2] - aligned[k][2]),
                     aligned[k][3] + f * (ts[k][3] - aligned[k][3])) for k in range(len(ts))])
    return pts, rmsd


# ------------------------------------------------------------------ ORCA


UKS_INSIDE_OUT = False  # --uks-inside-out: start the broken-symmetry series at the compressed end
SLOWCONV = False  # --slowconv: damped SCF for hard broken-symmetry series


def orca_input(points, charge, mult, nprocs, uks=False, broken_sym=False, hirshfeld=False, extra=""):
    """Multi-job input: one EnGrad single point per geometry ($new_job chain, MOs carried over)."""
    blocks = []
    for k, atoms in enumerate(points):
        kw = f"! {KEYWORDS} EnGrad" + (" UKS" if uks else "") + (" SlowConv" if SLOWCONV else "")
        b = [kw, f"%pal nprocs {nprocs} end", "%maxcore 2000"]
        if broken_sym and k == 0:
            b.append("%scf BrokenSym 1,1" + (" MaxIter 500" if SLOWCONV else "") + " end")
        elif SLOWCONV:
            b.append("%scf MaxIter 500 end")
        if hirshfeld:
            b.append("%output Print[P_Hirshfeld] 1 end")
        if extra:
            b.append(extra)
        b.append(f"* xyz {charge} {mult}")
        b += [f"{s} {x:.8f} {y:.8f} {z:.8f}" for s, x, y, z in atoms]
        b.append("*")
        blocks.append("\n".join(b))
    return "\n$new_job\n".join(blocks) + "\n"


def orca_opt_input(atoms, charge, mult, nprocs, uks=False, broken_sym=False):
    kw = f"! {KEYWORDS} Opt" + (" UKS" if uks else "")
    lines = [kw, f"%pal nprocs {nprocs} end", "%maxcore 2000"]
    if broken_sym:
        lines.append("%scf BrokenSym 1,1 end")
    lines.append(f"* xyz {charge} {mult}")
    return "\n".join(lines + [f"{s} {x:.8f} {y:.8f} {z:.8f}" for s, x, y, z in atoms] + ["*", ""])


def run_orca(workdir, inp_name, extra_keep=None):
    """Run ORCA in workdir. Deletes ORCA scratch files afterwards, keeping only what the repo
    needs; extra_keep additionally keeps an explicit set of exact filenames (used for NEB
    campaigns to retain the per-image geometries and the converged TS/CI structures without
    keeping the much larger .gbw/.bas/.tmp scratch ORCA writes per image)."""
    t0 = time.time()
    with (workdir / inp_name.replace(".inp", ".out")).open("w") as out:
        proc = subprocess.run([ORCA, inp_name], stdout=out, stderr=subprocess.STDOUT, cwd=workdir)
    # keep only what the repo needs: input, output (gzipped by the caller), the optimised geometry
    keep = {inp_name, inp_name.replace(".inp", ".out"), inp_name.replace(".inp", ".xyz"), inp_name.replace(".inp", ".out.gz"),
            "energies.json", "gradients.json", "points.xyz", "meta.json", "ref.xyz", "start.xyz", "start.opt.xyz",
            "frames.xyz", "input.xyz", "product.xyz", "ts_guess.xyz", "reactant.xyz"}
    if extra_keep:
        keep |= set(extra_keep)
    for f in workdir.glob("*"):
        if f.is_file() and f.name not in keep:
            try:
                f.unlink()
            except OSError:
                pass
    return proc.returncode, time.time() - t0


RE_E = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
RE_S2 = re.compile(r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)")


def parse_orca(out_text):
    """Per job: coordinates (A), energy (Eh), gradient (Eh/Bohr), last <S**2>, Hirshfeld charges."""
    jobs = re.split(r"\$+\s+JOB NUMBER\s+\d+\s+\$+", out_text)
    if len(jobs) == 1:  # single job: no marker
        jobs = [out_text]
    else:
        jobs = jobs[1:]
    results = []
    for chunk in jobs:
        rec = {"energy_eh": None, "gradient_eh_bohr": None, "coords_ang": None, "s2": None, "hirshfeld": None}
        m = re.search(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n((?:\s*\S+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\n)+)", chunk)
        if m:
            rec["coords_ang"] = [(t[0], float(t[1]), float(t[2]), float(t[3])) for t in (ln.split() for ln in m.group(1).strip().splitlines())]
        es = RE_E.findall(chunk)
        if es:
            rec["energy_eh"] = float(es[-1])
        m = re.search(r"CARTESIAN GRADIENT\n-+\n\n((?:\s*\d+\s+\S+\s+:\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\n)+)", chunk)
        if m:
            rec["gradient_eh_bohr"] = [[float(t[3]), float(t[4]), float(t[5])] for t in (ln.split() for ln in m.group(1).strip().splitlines())]
        s2 = RE_S2.findall(chunk)
        if s2:
            rec["s2"] = float(s2[-1])
        m = re.search(r"HIRSHFELD ANALYSIS.*?\n\s*ATOM\s+CHARGE\s+SPIN\s*\n((?:\s*\d+\s+\S+\s+-?\d+\.\d+\s+-?\d+\.\d+\n)+)", chunk, re.S)
        if m:
            rec["hirshfeld"] = [float(ln.split()[2]) for ln in m.group(1).strip().splitlines()]
        results.append(rec)
    return results


def orca_version(out_text):
    m = re.search(r"Program Version\s+(\S+)", out_text)
    return m.group(1) if m else "?"


# ------------------------------------------------------------------ jobs


class Job:
    def __init__(self, cls, system, points, charge, mult, uks=False, broken_sym=False, hirshfeld=False,
                 labels=None, tag="", note=""):
        self.cls, self.system, self.points = cls, system, points
        self.charge, self.mult, self.uks, self.broken_sym, self.hirshfeld = charge, mult, uks, broken_sym, hirshfeld
        self.labels = labels or [str(k) for k in range(len(points))]
        self.tag = tag
        self.note = note  # Claude Generated (Sep 2026) -- recorded verbatim into meta.json when set

    @property
    def dir(self):
        return REF / self.cls / self.system

    def done(self):
        e = self.dir / "energies.json"
        if not e.exists():
            return False
        try:
            d = json.loads(e.read_text())
            return d.get("n_ok", 0) == len(self.points)
        except Exception:
            return False

    def run(self, nprocs, log=None):
        d = self.dir
        d.mkdir(parents=True, exist_ok=True)
        inp = orca_input(self.points, self.charge, self.mult, nprocs, self.uks, self.broken_sym, self.hirshfeld)
        (d / "job.inp").write_text(inp)
        rc, wall = run_orca(d, "job.inp")
        out = (d / "job.out").read_text(errors="replace")
        with gzip.open(d / "job.out.gz", "wt") as gz:
            gz.write(out)
        (d / "job.out").unlink()
        recs = parse_orca(out)
        energies, ok = [], 0
        with (d / "points.xyz").open("w") as fx:
            for k, atoms in enumerate(self.points):
                rec = recs[k] if k < len(recs) else {"energy_eh": None, "gradient_eh_bohr": None, "s2": None, "hirshfeld": None}
                e = rec["energy_eh"] if rec["gradient_eh_bohr"] is not None else None
                if e is not None:
                    ok += 1
                energies.append({"point": k, "label": self.labels[k], "energy_eh": e, "s2": rec["s2"],
                                 "hirshfeld": rec["hirshfeld"], "gradient_eh_bohr": rec["gradient_eh_bohr"] if e is not None else None})
                fx.write(f"{len(atoms)}\nE={e if e is not None else 'nan'} charge={self.charge} mult={self.mult} point={k} label={self.labels[k]}\n")
                fx.write("".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, x, y, z in atoms))
        (d / "gradients.json").write_text(json.dumps({"unit": "Eh/Bohr", "gradients": [e["gradient_eh_bohr"] for e in energies]}))
        for e in energies:
            e.pop("gradient_eh_bohr")
        (d / "energies.json").write_text(json.dumps({"class": self.cls, "system": self.system, "tag": self.tag,
                                                     "charge": self.charge, "mult": self.mult, "uks": self.uks,
                                                     "broken_sym": self.broken_sym, "n_points": len(self.points), "n_ok": ok,
                                                     "orca_rc": rc, "wall_s": round(wall, 1), "points": energies}, indent=1))
        meta_extra = {"orca_version": orca_version(out), "keywords": KEYWORDS + " EnGrad",
                      "nprocs": nprocs, "uks": self.uks, "broken_sym": self.broken_sym}
        if self.note:
            meta_extra["note"] = self.note
        (d / "meta.json").write_text(json.dumps(meta(meta_extra), indent=1))
        return ok, len(self.points), wall


# ------------------------------------------------------------------ reference geometries


def reference_geometry(name, nprocs, log):
    """GFN2 (curcuma) then r2SCAN-3c (ORCA) optimisation; cached in ref/_geom/<name>/ref.xyz."""
    charge, mult, atoms = MOL[name]
    d = REF / "_geom" / name
    ref = d / "ref.xyz"
    if ref.exists():
        return read_xyz(ref), charge, mult
    d.mkdir(parents=True, exist_ok=True)
    start = d / "start.xyz"
    write_xyz(start, atoms, f"hand-written start geometry {name}")
    if len(atoms) > 1:
        subprocess.run([str(CURCUMA), "-opt", "start.xyz", "-method", "gfn2", "-charge", str(charge), "-spin", str(mult - 1),
                        "-no_bmt", "-verbosity", "0", "-threads", "1"], capture_output=True, text=True, cwd=d)
        g = d / "start.opt.xyz"
        atoms = read_xyz(g) if g.exists() else atoms
        (d / "opt.inp").write_text(orca_opt_input(atoms, charge, mult, nprocs, uks=(mult > 1)))
        rc, wall = run_orca(d, "opt.inp")
        out = (d / "opt.out").read_text(errors="replace")
        with gzip.open(d / "opt.out.gz", "wt") as gz:
            gz.write(out)
        (d / "opt.out").unlink()
        final = d / "opt.xyz"
        if final.exists() and "OPTIMIZATION RUN DONE" in out:
            atoms = read_xyz(final)
            log(f"  {name}: r2SCAN-3c geometry optimised ({wall:.0f} s)")
        else:
            log(f"  {name}: ORCA optimisation FAILED (rc={rc}), keeping the GFN2 geometry")
    write_xyz(ref, atoms, f"r2SCAN-3c optimised reference geometry ({name}, charge {charge}, mult {mult})")
    (d / "meta.json").write_text(json.dumps(meta({"keywords": KEYWORDS + " Opt"}), indent=1))
    return atoms, charge, mult


def md_snapshots(name, T, nframes, log):
    """GFN-FF MD snapshots for class D, cached in ref/_md/<name>_<T>/frames.xyz."""
    d = REF / "_md" / f"{name}_{T}"
    frames_file = d / "frames.xyz"
    if frames_file.exists():
        return read_xyz_frames(frames_file)
    atoms, charge, mult = reference_geometry(name, 4, log)
    d.mkdir(parents=True, exist_ok=True)
    write_xyz(d / "input.xyz", atoms, f"{name} reference geometry")
    steps_fs = nframes * 100 * 0.5  # dump every 100 steps of 0.5 fs
    subprocess.run([str(CURCUMA), "-md", "input.xyz", "-method", "gfnff", "-temperature", str(T), "-maxtime", str(int(steps_fs)),
                    "-md.time_step", "0.5", "-md.dump_frequency", "100", "-md.seed", "7", "-md.no_restart", "-md.thermostat", "csvr",
                    "-md.coupling", "100", "-no_bmt", "-verbosity", "0", "-threads", "1"], capture_output=True, text=True, cwd=d)
    trj = d / "input.trj.xyz"
    if not trj.exists():
        trj = d / "input.snapshots" / "input.trj.xyz"
    frames = read_xyz_frames(trj)[1:nframes + 1] if trj.exists() else []
    with frames_file.open("w") as f:
        for k, fr in enumerate(frames):
            f.write(f"{len(fr)}\n{name} GFN-FF MD {T} K frame {k}\n" + "".join(f"{s} {x:.6f} {y:.6f} {z:.6f}\n" for s, x, y, z in fr))
    shutil.rmtree(d / "input.snapshots", ignore_errors=True)
    return frames


def translate(atoms, dx, dy, dz):
    return [(s, x + dx, y + dy, z + dz) for s, x, y, z in atoms]


def bond_len(mol, nprocs, log):
    """Equilibrium bond length (Angstrom) of a cached 2-atom reference geometry."""
    atoms, _, _ = reference_geometry(mol, nprocs, log)
    return math.dist(atoms[0][1:], atoms[1][1:])


def read_gmtkn(subset, name):
    """Read a GMTKN55 benchmark structure: (atoms, charge, mult) from its .CHRG/.UHF files."""
    d = GMTKN / subset / name
    atoms = read_xyz(d / "struc.xyz")
    charge = int((d / ".CHRG").read_text().strip()) if (d / ".CHRG").exists() else 0
    uhf = int((d / ".UHF").read_text().strip()) if (d / ".UHF").exists() else 0
    return atoms, charge, uhf + 1


def optimize_endpoint(atoms, charge, mult, nprocs, log, workdir, tag, uks=False, broken_sym=False):
    """r2SCAN-3c unconstrained Opt of a hand-built NEB endpoint guess; cached as
    workdir/<tag>.xyz so a restarted campaign does not redo it."""
    workdir.mkdir(parents=True, exist_ok=True)
    cached = workdir / f"{tag}.xyz"
    if cached.exists():
        return read_xyz(cached)
    inp_name = f"{tag}_opt.inp"
    (workdir / inp_name).write_text(orca_opt_input(atoms, charge, mult, nprocs, uks=uks, broken_sym=broken_sym))
    rc, wall = run_orca(workdir, inp_name)
    out = (workdir / inp_name.replace(".inp", ".out")).read_text(errors="replace")
    with gzip.open(workdir / inp_name.replace(".inp", ".out.gz"), "wt") as gz:
        gz.write(out)
    (workdir / inp_name.replace(".inp", ".out")).unlink()
    final = workdir / inp_name.replace(".inp", ".xyz")
    if final.exists() and "OPTIMIZATION RUN DONE" in out:
        out_atoms = read_xyz(final)
        log(f"    {tag}: r2SCAN-3c endpoint optimised ({wall:.0f} s)")
    else:
        out_atoms = atoms
        log(f"    {tag}: ORCA endpoint optimisation FAILED (rc={rc}), keeping the starting guess")
    write_xyz(cached, out_atoms, f"r2SCAN-3c optimised NEB endpoint ({tag}, charge {charge}, mult {mult})")
    return out_atoms


# ------------------------------------------------------------------ class B: NEB-TS paths
# Claude Generated (Sep 2026)


def orca_neb_input(reactant, charge, mult, nprocs, nimages, uks=False, broken_sym=False,
                    has_ts_guess=False, maxiter=None):
    kw = f"! NEB-TS {KEYWORDS}" + (" UKS" if uks else "")
    lines = [kw, f"%pal nprocs {nprocs} end", "%maxcore 2000", "%neb",
              '  NEB_END_XYZFILE "product.xyz"', f"  NImages {nimages}", "  PREOPT_ENDS false"]
    if has_ts_guess:
        lines.append('  NEB_TS_XYZFILE "ts_guess.xyz"')
    if maxiter:
        lines.append(f"  MAXITER {maxiter}")
    lines.append("end")
    if broken_sym:
        lines.append("%scf BrokenSym 1,1 end")
    lines.append(f"* xyz {charge} {mult}")
    lines += [f"{s} {x:.8f} {y:.8f} {z:.8f}" for s, x, y, z in reactant]
    lines.append("*")
    return "\n".join(lines) + "\n"


class NebJob:
    """One NEB-TS path: builds reactant/product/TS-guess via `build_fn(nprocs, log)` ->
    (reactant, product, ts_guess_or_None, charge, mult, note), runs ORCA NEB-TS (retrying
    once at 12 images / MAXITER 300 on failure per the WP2 hard rule), then a normal EnGrad
    multi-job pass (reusing orca_input/parse_orca) on every converged image + the TS +
    both endpoints, written in the same energies.json/gradients.json/points.xyz/meta.json
    layout as classes A/C/D/E with an added per-point "role" field."""

    def __init__(self, system, build_fn, mult=1, uks=False, broken_sym=False, nimages=8, tag=""):
        self.cls = "B"
        self.system, self.build_fn = system, build_fn
        self.uks, self.broken_sym, self.nimages, self.tag = uks, broken_sym, nimages, tag
        self.charge, self.mult = 0, mult  # all 15 roadmap paths are neutral; mult known upfront for `plan`
        self.points = [None] * (nimages + 3)  # +2 endpoints +1 TS, for plan-time point counts

    @property
    def dir(self):
        return REF / self.cls / self.system

    def done(self):
        e = self.dir / "energies.json"
        if not e.exists():
            return False
        try:
            return json.loads(e.read_text()).get("status") in ("ok", "failed")
        except Exception:
            return False

    def _attempt(self, d, reactant, charge, mult, nprocs, nimages, has_ts_guess, maxiter, log):
        inp = orca_neb_input(reactant, charge, mult, nprocs, nimages, self.uks, self.broken_sym, has_ts_guess, maxiter)
        (d / "job.inp").write_text(inp)
        # ORCA writes the converged band as ONE multi-frame trajectory (job_MEP_trj.xyz, ordered
        # reactant..product, nimages+2 frames, each frame's comment line carries its own energy)
        # plus the separately-refined saddle point (job_NEB-TS_converged.xyz) -- not per-image
        # job_im<k>.xyz files (those exist only as transient per-iteration scratch, confirmed
        # against a live test run, 2026-09-11).
        keep = {"job_MEP_trj.xyz", "job_NEB-TS_converged.xyz", "job_NEB-CI_converged.xyz", "job.NEB.log"}
        rc, wall = run_orca(d, "job.inp", extra_keep=keep)
        out = (d / "job.out").read_text(errors="replace")
        with gzip.open(d / "job.out.gz", "wt") as gz:
            gz.write(out)
        (d / "job.out").unlink()
        ok = "THE TS OPTIMIZATION HAS CONVERGED" in out and (d / "job_NEB-TS_converged.xyz").exists() \
            and (d / "job_MEP_trj.xyz").exists()
        return ok, wall, rc, out, nimages

    def run(self, nprocs, log):
        d = self.dir
        d.mkdir(parents=True, exist_ok=True)
        try:
            reactant, product, ts_guess, charge, mult, note = self.build_fn(nprocs, log)
        except Exception as exc:
            (d / "energies.json").write_text(json.dumps({"class": "B", "system": self.system, "tag": self.tag,
                                                          "charge": self.charge, "mult": self.mult, "status": "blocked",
                                                          "n_points": len(self.points), "n_ok": 0, "wall_s": 0.0,
                                                          "error": f"endpoint build failed: {exc}"}, indent=1))
            (d / "meta.json").write_text(json.dumps(meta({"status": "blocked", "error": str(exc)}), indent=1))
            log(f"  B {self.system}: BLOCKED (endpoint build) {exc}")
            return "blocked", len(self.points), 0.0
        self.charge, self.mult = charge, mult
        write_xyz(d / "reactant.xyz", reactant, f"{self.system} reactant")
        write_xyz(d / "product.xyz", product, f"{self.system} product")
        if ts_guess is not None:
            write_xyz(d / "ts_guess.xyz", ts_guess, f"{self.system} TS guess")
        t0 = time.time()
        ok, wall, rc, out, nim = self._attempt(d, reactant, charge, mult, nprocs, self.nimages,
                                                ts_guess is not None, None, log)
        attempts = 1
        if not ok:
            log(f"  B {self.system}: attempt 1 did not converge (rc={rc}), retrying with 12 images / MAXITER 300")
            ok, wall2, rc, out, nim = self._attempt(d, reactant, charge, mult, nprocs, 12,
                                                     ts_guess is not None, 300, log)
            wall += wall2
            attempts = 2
        if not ok:
            (d / "energies.json").write_text(json.dumps({"class": "B", "system": self.system, "tag": self.tag,
                                                          "charge": charge, "mult": mult, "status": "failed",
                                                          "n_points": len(self.points), "n_ok": 0,
                                                          "attempts": attempts, "wall_s": round(wall, 1)}, indent=1))
            (d / "meta.json").write_text(json.dumps(meta({"status": "failed", "attempts": attempts, "note": note,
                                                           "keywords": KEYWORDS + " NEB-TS", "nprocs": nprocs,
                                                           "uks": self.uks, "charge": charge, "mult": mult}), indent=1))
            log(f"  B {self.system}: FAILED after {attempts} attempt(s), {wall:.0f} s")
            return "failed", len(self.points), wall
        # gather the converged path (single multi-frame trajectory, reactant..product) + the
        # separately-refined TS, then a normal EnGrad pass on all of them
        images = read_xyz_frames(d / "job_MEP_trj.xyz")
        ts_conv = read_xyz(d / "job_NEB-TS_converged.xyz")
        points = images + [ts_conv]
        roles = ["reactant"] + ["image"] * (len(images) - 2) + ["product", "ts"]
        labels = [f"role={r}_{k}" for k, r in enumerate(roles)]
        eg_inp = orca_input(points, charge, mult, nprocs, uks=self.uks, broken_sym=self.broken_sym, hirshfeld=False)
        (d / "job_engrad.inp").write_text(eg_inp)
        rc2, wall2 = run_orca(d, "job_engrad.inp")
        eg_out = (d / "job_engrad.out").read_text(errors="replace")
        with gzip.open(d / "job_engrad.out.gz", "wt") as gz:
            gz.write(eg_out)
        (d / "job_engrad.out").unlink()
        recs = parse_orca(eg_out)
        energies, ok_n = [], 0
        with (d / "points.xyz").open("w") as fx:
            for k, atoms in enumerate(points):
                rec = recs[k] if k < len(recs) else {"energy_eh": None, "gradient_eh_bohr": None, "s2": None}
                e = rec["energy_eh"] if rec.get("gradient_eh_bohr") is not None else None
                if e is not None:
                    ok_n += 1
                energies.append({"point": k, "label": labels[k], "role": roles[k], "energy_eh": e,
                                 "s2": rec.get("s2"), "gradient_eh_bohr": rec.get("gradient_eh_bohr") if e is not None else None})
                fx.write(f"{len(atoms)}\nE={e if e is not None else 'nan'} charge={charge} mult={mult} point={k} role={roles[k]}\n")
                fx.write("".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, x, y, z in atoms))
        (d / "gradients.json").write_text(json.dumps({"unit": "Eh/Bohr", "gradients": [e["gradient_eh_bohr"] for e in energies]}))
        for e in energies:
            e.pop("gradient_eh_bohr")
        total_wall = wall + wall2
        (d / "energies.json").write_text(json.dumps({"class": "B", "system": self.system, "tag": self.tag,
                                                      "charge": charge, "mult": mult, "uks": self.uks,
                                                      "broken_sym": self.broken_sym, "n_points": len(points), "n_ok": ok_n,
                                                      "orca_rc": rc2, "status": "ok", "attempts": attempts,
                                                      "wall_s": round(total_wall, 1), "points": energies}, indent=1))
        (d / "meta.json").write_text(json.dumps(meta({"status": "ok", "attempts": attempts, "note": note,
                                                       "keywords": KEYWORDS + " NEB-TS then EnGrad", "nprocs": nprocs,
                                                       "uks": self.uks, "broken_sym": self.broken_sym,
                                                       "charge": charge, "mult": mult, "nimages": nim}), indent=1))
        log(f"  B {self.system}: {ok_n}/{len(points)} points ok, {total_wall:.0f} s ({attempts} NEB attempt(s))")
        return "ok", len(points), total_wall


def _bh76_3atom_builder(ts_subdir, reactant_spec, product_spec, mult):
    a, b, molA, freeA = reactant_spec
    c, d_, molB, freeB = product_spec

    def build(nprocs, log):
        ts_atoms, _, _ = read_gmtkn("BH76", ts_subdir)
        reactant = build_3atom_endpoint(ts_atoms, (a, b), bond_len(molA, nprocs, log), freeA)
        product = build_3atom_endpoint(ts_atoms, (c, d_), bond_len(molB, nprocs, log), freeB)
        return reactant, product, ts_atoms, 0, mult, f"3-atom endpoints from BH76/{ts_subdir}, bonded pairs set to {molA}/{molB} eq. length"
    return build


def _n2_h_builder():
    def build(nprocs, log):
        ts_atoms, _, _ = read_gmtkn("BH76", "hn2ts")
        reactant = build_3atom_endpoint(ts_atoms, (0, 1), bond_len("n2", nprocs, log), 2)
        product, _, _ = read_gmtkn("BH76", "hn2")
        return reactant, product, ts_atoms, 0, 2, "reactant: N2 bonded pair set to eq length, H pushed to 4.5 A; product: BH76 hn2 (bound N2H)"
    return build


def _n2h_h_builder():
    def build(nprocs, log):
        hn2_atoms, _, _ = read_gmtkn("BH76", "hn2")
        h_atoms, _, _ = read_gmtkn("BH76", "h")
        guess = stack_along_z(hn2_atoms, h_atoms, 3.5)
        reactant = optimize_endpoint(guess, 0, 1, nprocs, log, REF / "B" / "n2h_h_n2h2", "reactant", uks=True, broken_sym=True)
        product, _, _ = reference_geometry("n2h2", nprocs, log)
        return reactant, product, None, 0, 1, "reactant: BH76 hn2 + h stacked 3.5 A apart, r2SCAN-3c Opt (UKS BS); product: class A n2h2 reference geometry"
    return build


def _chain_builder(mol_a, mol_b, mol_product, tag):
    def build(nprocs, log):
        a_atoms, _, _ = reference_geometry(mol_a, nprocs, log)
        b_atoms, _, _ = reference_geometry(mol_b, nprocs, log)
        guess = stack_along_z(a_atoms, b_atoms, 3.5)
        reactant = optimize_endpoint(guess, 0, 1, nprocs, log, REF / "B" / f"{mol_a}_{mol_b}_{mol_product}", "reactant")
        product, _, _ = reference_geometry(mol_product, nprocs, log)
        return reactant, product, None, 0, 1, tag
    return build


def _n2h4_2nh3_builder():
    def build(nprocs, log):
        n2h4_atoms, _, _ = reference_geometry("n2h4", nprocs, log)
        h2_atoms, _, _ = reference_geometry("h2", nprocs, log)
        guess = stack_along_z(n2h4_atoms, h2_atoms, 3.5)
        reactant = optimize_endpoint(guess, 0, 1, nprocs, log, REF / "B" / "n2h4_h2_2nh3", "reactant")
        nh3_atoms, _, _ = reference_geometry("nh3", nprocs, log)
        sep = 4.0
        mol_a = translate(nh3_atoms, 0.0, 0.0, -sep)
        mol_b = translate(nh3_atoms, 0.0, 0.0, sep)
        # reactant atom order: N,N,H1,H2,H3,H4(orig N2H4),H5,H6(orig H2) -- product keeps that
        # slot order, H1,H2,H5 assigned to NH3_a (built around N1), H3,H4,H6 to NH3_b (N2)
        product = [mol_a[0], mol_b[0], mol_a[1], mol_a[2], mol_b[1], mol_b[2], mol_a[3], mol_b[3]]
        return reactant, product, None, 0, 1, "reactant: class A n2h4 + h2 stacked 3.5 A apart, r2SCAN-3c Opt; product: 2 separate optimised NH3 (class A geometry) 8 A apart"
    return build


def _px13_transit_builder(subset_name):
    def build(nprocs, log):
        reactant, _, _ = read_gmtkn("PX13", subset_name)
        ts_atoms, _, _ = read_gmtkn("PX13", f"{subset_name}_ts")
        pts, rmsd = transit_points(reactant, ts_atoms, [0.0, 2.0])
        return pts[0], pts[1], ts_atoms, 0, 1, f"reactant aligned onto PX13/{subset_name}_ts (Kabsch rmsd {rmsd:.4f} A); product = 2*TS - reactant (degenerate mirror)"
    return build


# (name, builder, mult, uks, broken_sym)  -- the 15 roadmap NEB-TS paths (all neutral, charge 0)
BPATHS = [
    ("rkt06_h_h2", _bh76_3atom_builder("RKT06", (0, 1, "h2", 2), (0, 2, "h2", 1), 2), 2, True, False),
    ("hclhts_h_hcl", _bh76_3atom_builder("hclhts", (1, 0, "hcl", 2), (1, 2, "hcl", 0), 2), 2, True, False),
    ("hfhts_h_hf", _bh76_3atom_builder("hfhts", (1, 0, "hf", 2), (1, 2, "hf", 0), 2), 2, True, False),
    ("rkt01_h_hcl_h2_cl", _bh76_3atom_builder("RKT01", (0, 1, "hcl", 2), (0, 2, "h2", 1), 2), 2, True, False),
    ("rkt10_f_h2_hf_h", _bh76_3atom_builder("RKT10", (0, 2, "h2", 1), (0, 1, "hf", 2), 2), 2, True, False),
    ("hf2ts_h_f2_hf_f", _bh76_3atom_builder("hf2ts", (1, 2, "f2", 0), (0, 1, "hf", 2), 2), 2, True, False),
    ("rkt14_h_oh_h2_o", _bh76_3atom_builder("RKT14", (1, 0, "oh", 2), (0, 2, "h2", 1), 3), 3, True, False),
    ("n2_h_n2h", _n2_h_builder(), 2, True, False),
    ("n2h_h_n2h2", _n2h_h_builder(), 1, True, True),
    ("n2_h2_n2h2", _chain_builder("n2", "h2", "n2h2", "N2+H2 -> N2H2 chain start"), 1, False, False),
    ("n2h2_h2_n2h4", _chain_builder("n2h2", "h2", "n2h4", "N2H2+H2 -> N2H4 chain step 2"), 1, False, False),
    ("n2h4_h2_2nh3", _n2h4_2nh3_builder(), 1, False, False),
    ("px13_hf_2", _px13_transit_builder("hf_2"), 1, False, False),
    ("px13_h2o_2", _px13_transit_builder("h2o_2"), 1, False, False),
    ("px13_nh3_2", _px13_transit_builder("nh3_2"), 1, False, False),
]


# ------------------------------------------------------------------ class E remainder + class L (lost scans)
# Claude Generated (Sep 2026)

FCH3F_UMBRELLA_GRID = [75, 80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5, 100, 105]  # F-C-H angle, deg
NH4_NH3_FRAC_GRID = [-1.0, -0.85, -0.7, -0.55, -0.4, -0.2, 0.0, 0.2, 0.4, 0.55, 0.7, 0.85, 1.0]
TRANSIT_FRAC_GRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # reactant(0) -> TS(1)
CH2M_ANGLES = [90, 100, 110, 120, 130, 145]  # TODO#3b grid, C-H fixed at 1.129 A
FORMAMIDINE_ANGLES = [119, 130, 140, 150, 160, 170, 179]  # TODO#1 grid, imine C-N-H


def build_e_extra_jobs(only, nprocs, log):
    """Class E remainder (WP2 item 1): SN2 umbrella, (HF)2/(H2O)2 proton-transfer transits,
    NH4+ + NH3 proton transfer, HCOO-...HF stretch. All with Hirshfeld charges."""
    jobs = []

    def want(name):
        return not only or name in only

    if want("fch3f_umbrella"):
        ts_atoms, charge, mult = read_gmtkn("BH76", "fch3fts")  # F,C,H,H,H,F ; charge -1, mult 1
        pts, labels = [], []
        for ang in FCH3F_UMBRELLA_GRID:
            a = ts_atoms
            for h_idx in (2, 3, 4):
                a = bend_angle(a, 0, 1, h_idx, ang)
            pts.append(a)
            labels.append(f"FCH={ang}")
        jobs.append(Job("E", "fch3f_umbrella", pts, charge, mult, hirshfeld=True, labels=labels,
                        tag="[F...CH3...F]- umbrella scan from BH76/fch3fts (TODO#6)"))

    if want("hf2_transit") or want("h2o2_transit") or want("nh3_2_transit"):
        for tagname, subset in (("hf2_transit", "hf_2"), ("h2o2_transit", "h2o_2"), ("nh3_2_transit", "nh3_2")):
            if not want(tagname):
                continue
            reactant, _, _ = read_gmtkn("PX13", subset)
            ts_atoms, _, _ = read_gmtkn("PX13", f"{subset}_ts")
            pts, rmsd = transit_points(reactant, ts_atoms, TRANSIT_FRAC_GRID)
            labels = [f"frac={f:.2f}" for f in TRANSIT_FRAC_GRID]
            log(f"    {tagname}: reactant/TS atom-matching Kabsch rmsd {rmsd:.4f} A")
            jobs.append(Job("E", tagname, pts, 0, 1, hirshfeld=True, labels=labels,
                            tag=f"PX13/{subset} proton-transfer transit, reactant->TS (Kabsch-matched)"))

    if want("nh4_nh3_pt"):
        d_nn = 2.75  # N...N literature range for [H3N...H...NH3]+, ~2.70-2.80 A (not fitted, a reasonable choice)
        rterm = 1.021  # N-Hterm (NH4+-like)
        costheta, sintheta = -1.0 / 3.0, math.sqrt(1.0 - 1.0 / 9.0)  # tetrahedral

        def term_h(nz, sign, phi_deg):
            phi = math.radians(phi_deg)
            return ("H", sintheta * math.cos(phi) * rterm, sintheta * math.sin(phi) * rterm, nz + sign * costheta * rterm)

        n1z, n2z = -d_nn / 2, d_nn / 2
        frame = [("N", 0.0, 0.0, n1z), ("N", 0.0, 0.0, n2z)]
        frame += [term_h(n1z, +1, phi) for phi in (0, 120, 240)]
        frame += [term_h(n2z, -1, phi) for phi in (60, 180, 300)]
        rcov_offset = d_nn / 2 - rterm
        pts = [frame + [("H", 0.0, 0.0, f * rcov_offset)] for f in NH4_NH3_FRAC_GRID]
        labels = [f"zbridge={f * rcov_offset:.3f}" for f in NH4_NH3_FRAC_GRID]
        jobs.append(Job("E", "nh4_nh3_pt", pts, 1, 1, hirshfeld=True, labels=labels,
                        tag=f"hand-built [H3N...H...NH3]+ bridging-H scan (N...N={d_nn} A)"))

    if want("ahb21_21_stretch"):
        atoms, charge, mult = read_gmtkn("AHB21", "21")  # O,O,C,H,F,H ; charge -1
        req = math.dist(atoms[4][1:], atoms[5][1:])
        pts = [stretched(atoms, 4, 5, f * req) for f in CURVE_GRID]
        labels = [f"r={f * req:.4f}" for f in CURVE_GRID]
        jobs.append(Job("E", "ahb21_21_stretch", pts, charge, mult, hirshfeld=True, labels=labels,
                        tag="HCOO-...HF F-H stretch from GMTKN55 AHB21/21"))
    return jobs


def build_lost_scan_jobs(only, nprocs, log):
    """Class L (WP2 item 2): recompute the lost scans of REV_GFNFF_TODO.md to +-0.1 kcal/mol."""
    jobs = []

    def want(name):
        return not only or name in only

    if want("ch2m_bend"):
        r = 1.129
        pts, labels = [], []
        for ang in CH2M_ANGLES:
            half = math.radians(ang / 2.0)
            pts.append([("C", 0.0, 0.0, 0.0), ("H", r * math.cos(half), r * math.sin(half), 0.0),
                        ("H", r * math.cos(half), -r * math.sin(half), 0.0)])
            labels.append(f"HCH={ang}")
        # CH2- has 9 electrons (6 C + 1 + 1 H + 1 extra) -- odd, so it is a doublet radical
        # anion, not the singlet originally coded here (ORCA rejects mult=1/9 electrons
        # outright: "multiplicity (1) is odd and number of electrons (9) is odd -> impossible").
        # Claude Generated (Sep 2026)
        jobs.append(Job("L", "ch2m_bend", pts, -1, 2, uks=True, hirshfeld=True, labels=labels,
                        tag="CH2- bending scan, C-H fixed 1.129 A, UKS doublet (TODO#3b)"))

    if want("formamidine_imine"):
        atoms, charge, mult = reference_geometry("formamidine", nprocs, log)
        pts = [bend_angle(atoms, 0, 1, 4, ang) for ang in FORMAMIDINE_ANGLES]  # C=0, N(imine)=1, H(imine)=4
        labels = [f"CNH={ang}" for ang in FORMAMIDINE_ANGLES]
        jobs.append(Job("L", "formamidine_imine", pts, charge, mult, labels=labels,
                        tag="formamidine imine C-N-H angle scan (TODO#1)"))

    if want("cl2m_at_ea25"):
        atoms, charge, mult = reference_geometry("cl2m", nprocs, log)
        pts = [stretched(atoms, 0, 1, 2.7300)]
        jobs.append(Job("L", "cl2m_at_ea25", pts, charge, mult, uks=True, hirshfeld=True, labels=["r=2.7300"],
                        tag="Cl2- at the GMTKN55 G21EA/EA_25 geometry (TODO#3)"))
    if want("cl_minus"):
        jobs.append(Job("L", "cl_minus", [[("Cl", 0.0, 0.0, 0.0)]], -1, 1, hirshfeld=True, labels=["Cl-"],
                        tag="Cl- single point (TODO#3 dissociation reference)"))
    if want("cl_radical"):
        jobs.append(Job("L", "cl_radical", [[("Cl", 0.0, 0.0, 0.0)]], 0, 2, uks=True, hirshfeld=True, labels=["Cl"],
                        tag="Cl radical single point, UKS doublet (TODO#3 dissociation reference)"))
    return jobs


# ------------------------------------------------------------------ plan


def build_jobs(classes, only, nprocs, log):
    jobs = []
    if "A" in classes:
        for mol, i, j, label in CURVES:
            if only and mol not in only and f"{mol}_{label}" not in only:
                continue
            atoms, charge, mult = reference_geometry(mol, nprocs, log)
            req = math.dist(atoms[i][1:], atoms[j][1:])
            pts = [stretched(atoms, i, j, f * req) for f in CURVE_GRID]
            labels = [f"r={f * req:.4f}" for f in CURVE_GRID]
            sysname = f"{mol}_{label.replace('#', 'T').replace('=', 'D')}"
            if mult == 1:
                jobs.append(Job("A", sysname + "_rks", pts, charge, 1, labels=labels, tag=f"{label} rigid stretch, RKS"))
                if UKS_INSIDE_OUT:  # retry strategy: flip the spins near r_eq and follow the curve outward
                    jobs.append(Job("A", sysname + "_uks", pts, charge, 1, uks=True, broken_sym=True,
                                    labels=labels, tag=f"{label} rigid stretch, UKS broken symmetry, inside out"))
                else:
                    jobs.append(Job("A", sysname + "_uks", list(reversed(pts)), charge, 1, uks=True, broken_sym=True,
                                    labels=list(reversed(labels)), tag=f"{label} rigid stretch, UKS broken symmetry, outside in"))
            else:
                jobs.append(Job("A", sysname + "_uks", pts, charge, mult, uks=True, labels=labels, tag=f"{label} rigid stretch, UKS mult {mult}"))
    if "C" in classes:
        for name, host, att, anchor, rule, charge, mult in APPROACH:
            if only and name not in only:
                continue
            hatoms, _, _ = reference_geometry(host, nprocs, log)
            aatoms = [("H", 0.0, 0.0, 0.0)] if att == "H" else reference_geometry(att, nprocs, log)[0]
            pts = [approach(hatoms, aatoms, anchor, rule, dd) for dd in APPROACH_GRID]
            labels = [f"d={dd:.2f}" for dd in APPROACH_GRID]
            jobs.append(Job("C", name, pts, charge, mult, uks=(mult > 1), labels=labels, tag=f"{att} approaching {host} atom {anchor}"))
    if "E" in classes:
        for mol, i, j, label in ECURVES:
            if only and mol not in only:
                continue
            note = ""
            if mol in ECURVES_FIXED_R:
                # r2SCAN-3c unconstrained Opt diverges for this species (SCF wanders into a
                # dissociative channel instead of the bound anion minimum) -- skip Opt
                # entirely and scan around a fixed literature bond length instead.
                # Claude Generated (Sep 2026), coordinator-flagged: f2m_F-F- originally came
                # out with F-F from 4e5 to 2e6 A (the diverged Opt geometry fed straight into
                # the rigid stretch grid).
                charge, mult, atoms = MOL[mol]
                req = ECURVES_FIXED_R[mol]
                note = f"r2SCAN-3c unconstrained Opt diverges for {mol}; using fixed r={req} A (no Opt run), Sep 2026 fix"
            else:
                atoms, charge, mult = reference_geometry(mol, nprocs, log)
                req = math.dist(atoms[i][1:], atoms[j][1:])
            pts = [stretched(atoms, i, j, f * req) for f in CURVE_GRID]
            labels = [f"r={f * req:.4f}" for f in CURVE_GRID]
            jobs.append(Job("E", f"{mol}_{label.replace('(', '').replace(')', '')}", pts, charge, mult, uks=True, hirshfeld=True,
                            labels=labels, tag=f"{label} rigid stretch with Hirshfeld charges", note=note))
        jobs += build_e_extra_jobs(only, nprocs, log)
    if "L" in classes:
        jobs += build_lost_scan_jobs(only, nprocs, log)
    if "B" in classes:
        for name, build_fn, mult, uks, broken_sym in BPATHS:
            if only and name not in only:
                continue
            jobs.append(NebJob(name, build_fn, mult=mult, uks=uks, broken_sym=broken_sym, nimages=8, tag=f"NEB-TS path {name}"))
    if "D" in classes:
        for mol in OFFEQ:
            if only and mol not in only:
                continue
            _, charge, mult = reference_geometry(mol, nprocs, log)
            for T in OFFEQ_T:
                frames = md_snapshots(mol, T, OFFEQ_FRAMES, log)
                if frames:
                    jobs.append(Job("D", f"{mol}_{T}K", frames, charge, mult, labels=[f"frame{k}" for k in range(len(frames))],
                                    tag=f"GFN-FF MD snapshots at {T} K"))
    return jobs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["plan", "run", "status"])
    ap.add_argument("--classes", nargs="+", default=["A", "C"],
                     help="A dissociation curves, C hyper-coordination, D off-eq MD snapshots, "
                          "E charged/proton-transfer (Hirshfeld), L lost-scan reproductions, B NEB-TS paths")
    ap.add_argument("--only", nargs="*", help="system names (molecule for A/D/E, approach name for C, "
                                              "job name for E-extra/L/B)")
    ap.add_argument("--jobs", type=int, default=6, help="concurrent ORCA jobs")
    ap.add_argument("--nprocs", type=int, default=4, help="cores per ORCA job")
    ap.add_argument("--uks-inside-out", action="store_true", help="retry strategy for UKS curves that failed outside-in")
    ap.add_argument("--slowconv", action="store_true", help="damped SCF (SlowConv, MaxIter 500) for hard series")
    args = ap.parse_args()
    global UKS_INSIDE_OUT, SLOWCONV
    UKS_INSIDE_OUT, SLOWCONV = args.uks_inside_out, args.slowconv

    REF.mkdir(parents=True, exist_ok=True)
    logdir = REF / "_log"
    logdir.mkdir(exist_ok=True)
    logfile = logdir / f"campaign_{datetime.datetime.now():%Y%m%d_%H%M%S}.log"

    def log(msg):
        line = f"[{datetime.datetime.now():%H:%M:%S}] {msg}"
        print(line, flush=True)
        with logfile.open("a") as f:
            f.write(line + "\n")

    if args.cmd == "status":
        sys.path.insert(0, str(REPO / "scripts"))
        from revgfnff_data import summary
        summary()
        return

    jobs = build_jobs(args.classes, args.only, args.nprocs, log)
    todo = [j for j in jobs if not j.done()]
    log(f"{len(jobs)} jobs planned, {len(todo)} to run ({sum(len(j.points) for j in todo)} points), "
        f"{args.jobs} x {args.nprocs} cores")
    if args.cmd == "plan":
        for j in jobs:
            print(f"  {'done' if j.done() else 'todo'} {j.cls} {j.system:28s} {len(j.points):3d} pts charge {j.charge} mult {j.mult} {'UKS' if j.uks else 'RKS'}{' BS' if j.broken_sym else ''}  {j.tag}")
        return

    def work(j):
        try:
            ok, n, wall = j.run(args.nprocs, log)
            if isinstance(ok, int):  # Job: NebJob already logs its own detailed status internally
                log(f"  {j.cls} {j.system}: {ok}/{n} points ok, {wall:.0f} s")
        except Exception as exc:
            log(f"  {j.cls} {j.system}: FAILED {exc}")

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        list(ex.map(work, todo))
    log("campaign finished")


if __name__ == "__main__":
    main()
