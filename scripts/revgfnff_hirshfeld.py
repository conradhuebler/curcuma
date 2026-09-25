#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP2, class H reporting
"""Per-bond-type Hirshfeld charge table from the class-H reference scans.

Class H (scripts/revgfnff_ref.py) re-runs every class-A rigid stretch with
`%output Print[P_Hirshfeld] 1`; the points are identical to class A's (same cached reference
geometry, same CURVE_GRID), so point k of H/<system> is point k of A/<system>.

Why (FABLE_ROADMAP_REVIEW.md item 7): GFN-FF's Coulomb term drifts by -20..-35 kcal/mol
between r_eq and 3.5 r_eq for polar heavy-heavy bonds (C=O, C=N, N-O, O-O...). If the
reference charges move by a comparable amount along the stretch, the drift is physics and the
bond term's depth D has to carry it; if they barely move, it is curcuma's EEQ (chi(CN)) and
belongs to the charge model. Falsifier: Hirshfeld far-point charges agree with EEQ within
0.1 e -> the drift is physics.

    python scripts/revgfnff_hirshfeld.py            # per-curve table
    python scripts/revgfnff_hirshfeld.py --pair     # aggregated per element pair
    python scripts/revgfnff_hirshfeld.py --missing  # what is still absent
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import revgfnff_ref as R

NEAR, FAR = 1.00, 3.50  # CURVE_GRID multipliers: the equilibrium and the farthest point


import math


def pick(info, req, f):
    """Point whose label r matches f * r_eq (the labels are r values in Angstrom)."""
    want = f * req
    best, bd = None, None
    for p in info["points"]:
        if p.get("hirshfeld") is None:
            continue
        try:
            r = float(p["label"].split("=")[1])
        except (IndexError, ValueError):
            continue
        d = abs(r - want)
        if bd is None or d < bd:
            best, bd = p, d
    return best


def series():
    """Yield (sysname, tag, info, mol, i, j, label, req, near, far) for every usable class-H run."""
    for mol, i, j, label in R.CURVES:
        sysname = f"{mol}_{label.replace('#', 'T').replace('=', 'D')}"
        atoms = R.read_xyz(R.REF / "_geom" / mol / "ref.xyz")
        req = math.dist(atoms[i][1:], atoms[j][1:])
        for tag in ("_rks", "_uks"):
            d = R.REF / "H" / (sysname + tag)
            if not (d / "energies.json").exists():
                continue
            info = json.loads((d / "energies.json").read_text())
            near, far = pick(info, req, NEAR), pick(info, req, FAR)
            if near is None or far is None:
                continue
            yield sysname, tag, info, mol, i, j, label, req, near, far


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", action="store_true", help="aggregate per element pair")
    ap.add_argument("--missing", action="store_true", help="list class-A curves without a class-H run")
    args = ap.parse_args()

    rows = list(series())
    if args.missing:
        have = {r[0] for r in rows}
        miss = [f"{m}_{l}" for m, _, _, l in R.CURVES if f"{m}_{l.replace('#', 'T').replace('=', 'D')}" not in have]
        print(f"class H: {len(have)} curve-series present, {len(miss)} missing")
        for m in miss:
            print("   missing:", m)
        return
    def side_sums(mol, i, j):
        atoms = R.read_xyz(R.REF / "_geom" / mol / "ref.xyz")
        return R.bonded_side(atoms, i, j)  # atoms on j's side

    if args.pair:
        # RKS and broken-symmetry UKS are different electronic states and are aggregated
        # separately (the BS solution is a diradical at r_eq but the physical homolysis limit
        # at the far point -- mixing them would average two unrelated numbers).
        agg = {}
        for sysname, tag, info, mol, i, j, label, req, near, far in rows:
            key = (tuple(sorted((R.MOL[mol][2][i][0], R.MOL[mol][2][j][0]))), tag)
            agg.setdefault(key, []).append((sysname + tag, near, far, i, j))
        print("pair  state  n | mean dq_i  mean dq_j | max|dq|  largest")
        for key in sorted(agg):
            vals = [(name, far["hirshfeld"][i] - near["hirshfeld"][i],
                     far["hirshfeld"][j] - near["hirshfeld"][j]) for name, near, far, i, j in agg[key]]
            di = [v[1] for v in vals]
            dj = [v[2] for v in vals]
            worst = max(vals, key=lambda v: max(abs(v[1]), abs(v[2])))
            print(f"{key[0][0]:>2s}-{key[0][1]:<2s} {key[1].strip('_'):4s} {len(vals):2d} | "
                  f"{sum(di) / len(di):+.3f}   {sum(dj) / len(dj):+.3f}   | "
                  f"{max(max(abs(x) for x in di), max(abs(x) for x in dj)):.3f}     {worst[0]}")
        return

    print("system                       bond  r_eq   q_i(r_eq) q_i(far) dq_i | q_j(r_eq) q_j(far) dq_j | frag_j dq  n_ok")
    for sysname, tag, info, mol, i, j, label, req, near, far in rows:
        qn, qf = near["hirshfeld"], far["hirshfeld"]
        side = side_sums(mol, i, j)
        qs_n = sum(qn[k] for k in side)
        qs_f = sum(qf[k] for k in side)
        print(f"{sysname + tag:26s} {label:6s} {req:.3f}  {qn[i]:+.3f} {qf[i]:+.3f} {qf[i] - qn[i]:+.3f} |"
              f" {qn[j]:+.3f} {qf[j]:+.3f} {qf[j] - qn[j]:+.3f} | {qs_f - qs_n:+.3f} {info['n_ok']:3d}/{info['n_points']}")


if __name__ == "__main__":
    main()
