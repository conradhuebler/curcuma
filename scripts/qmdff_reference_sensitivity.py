#!/usr/bin/env python3
"""
How much does the choice of the QMDFF training conformer matter?

QMDFF is fitted to ONE Hessian in ONE basin. A conformer search cannot know the
global minimum in advance, so in practice the force field is fitted wherever the
search happens to be. This script quantifies the consequence: it repeats the fit
at several different reference conformers and measures, for each, how well the
resulting QMDFF@QM surface reproduces the reference ranking of the whole ensemble.

It also reports how the quality depends on the *distance* from the training basin,
by splitting the ensemble into conformers close to and far from the reference
(heavy-atom RMSD), which is the direct test of the "degrades away from the fitting
point" caveat.

Usage:
    scripts/qmdff_reference_sensitivity.py ensemble.xyz --refs 0,5,17,30 [options]

Reuses the machinery of qmdff_conformer_evaluation.py.

Claude Generated (August 2026)
Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qmdff_conformer_evaluation import (  # noqa: E402
    HARTREE_TO_KJ, pearson, read_xyz_frames, relative, run_single_point,
    spearman, write_xyz)


def heavy_atom_rmsd(frame_a, frame_b):
    """Plain (unaligned-agnostic) heavy-atom RMSD after centroid removal.

    No Kabsch rotation: the conformers come from the same optimisation pipeline
    and are already in a common frame. This is only used to bin conformers by
    distance from the training basin, so a monotone proxy is sufficient.
    """
    def coords(frame):
        out = []
        for line in frame[2:]:
            parts = line.split()
            if len(parts) >= 4 and parts[0].upper() not in ("H", "1"):
                out.append([float(parts[1]), float(parts[2]), float(parts[3])])
        return out

    a, b = coords(frame_a), coords(frame_b)
    if len(a) != len(b) or not a:
        return float("nan")
    for pts in (a, b):
        for axis in range(3):
            mean = sum(p[axis] for p in pts) / len(pts)
            for p in pts:
                p[axis] -= mean
    return math.sqrt(sum((p[k] - q[k]) ** 2 for p, q in zip(a, b) for k in range(3)) / len(a))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ensemble")
    parser.add_argument("--refs", default="", help="comma-separated conformer indices")
    parser.add_argument("--n-refs", type=int, default=4,
                        help="if --refs is empty: spread this many references over the "
                             "GFN2 energy range (minimum, maximum and evenly in between)")
    parser.add_argument("--curcuma", default="release/curcuma")
    parser.add_argument("--workdir", default="qmdff_refscan")
    parser.add_argument("--frames", type=int, default=0)
    parser.add_argument("--method", default="gfn2")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    curcuma = os.path.abspath(args.curcuma)
    workdir = os.path.abspath(args.workdir)
    os.makedirs(workdir, exist_ok=True)

    frames = read_xyz_frames(args.ensemble)
    if args.frames > 0:
        frames = frames[:args.frames]
    natoms = int(frames[0][0].split()[0])
    paths = []
    for index, frame in enumerate(frames):
        path = os.path.join(workdir, f"conf_{index:04d}.xyz")
        write_xyz(path, frame)
        paths.append(path)
    print(f"ensemble: {args.ensemble} — {len(frames)} conformers, {natoms} atoms\n")

    # Reference ranking
    print(f"[1] {args.method} single points ...", flush=True)
    e_ref = []
    for index, path in enumerate(paths):
        energy, _ = run_single_point(curcuma, path, args.method, args.threads)
        e_ref.append(energy)
        print(f"\r  {index + 1}/{len(paths)}", end="", flush=True)
    print()
    rel_ref = relative(e_ref)

    if args.refs:
        refs = [int(v) for v in args.refs.split(",")]
    else:
        order = sorted(range(len(e_ref)), key=lambda i: e_ref[i])
        picks = [round(i * (len(order) - 1) / max(1, args.n_refs - 1))
                 for i in range(args.n_refs)]
        refs = [order[p] for p in picks]
    print(f"training conformers: {refs}\n")

    rows = []
    for ref in refs:
        tag = f"ref{ref:04d}"
        fit_dir = os.path.join(workdir, tag)
        os.makedirs(fit_dir, exist_ok=True)
        fit_xyz = os.path.join(fit_dir, f"{tag}.xyz")
        fitted = os.path.join(fit_dir, f"{tag}.param.json")
        shutil.copy(paths[ref], fit_xyz)

        if not os.path.exists(fitted):
            print(f"[fit] conformer {ref} ({args.method} Hessian, "
                  f"{6 * natoms} gradient calls) ...", flush=True)
            subprocess.run([curcuma, "-qmdfffit", fit_xyz, "-method", args.method,
                            "-threads", str(args.threads), "-v", "1"],
                           capture_output=True, text=True, cwd=fit_dir, timeout=86400)
        if not os.path.exists(fitted):
            print(f"  fit at {ref} FAILED", file=sys.stderr)
            continue

        energies, rmsds = [], []
        for index, path in enumerate(paths):
            shutil.copy(fitted, path.replace(".xyz", ".param.json"))
            energy, _ = run_single_point(curcuma, path, "qmdff", args.threads)
            energies.append(energy if energy is not None else float("nan"))
            rmsds.append(heavy_atom_rmsd(frames[ref], frames[index]))

        valid = [i for i, v in enumerate(energies) if v == v]
        x = [rel_ref[i] for i in valid]
        y = relative([energies[i] for i in valid])
        rho = spearman(x, y)
        r = pearson(x, y)
        mad = sum(abs(a - b) for a, b in zip(x, y)) / len(x)

        true_min = e_ref.index(min(e_ref))
        order = sorted(valid, key=lambda i: energies[i])
        rank = order.index(true_min) + 1

        # near / far split by heavy-atom RMSD from the training conformer
        finite = [i for i in valid if rmsds[i] == rmsds[i]]
        median_rmsd = sorted(rmsds[i] for i in finite)[len(finite) // 2] if finite else 0.0
        near = [i for i in finite if rmsds[i] <= median_rmsd]
        far = [i for i in finite if rmsds[i] > median_rmsd]

        def subset_mad(subset):
            if len(subset) < 3:
                return float("nan")
            xs = [rel_ref[i] for i in subset]
            ys = [energies[i] * HARTREE_TO_KJ for i in subset]
            shift = sum(a - b for a, b in zip(xs, ys)) / len(subset)
            return sum(abs(a - (b + shift)) for a, b in zip(xs, ys)) / len(subset)

        rows.append({"ref": ref, "ref_rank_in_reference": sorted(rel_ref).index(rel_ref[ref]) + 1,
                     "spearman": rho, "pearson": r, "mad_kj": mad,
                     "spread_kj": max(y), "rank_of_reference_minimum": rank,
                     "mad_near_kj": subset_mad(near), "mad_far_kj": subset_mad(far),
                     "median_rmsd": median_rmsd})
        print(f"  ref {ref:>4}: rho={rho:6.3f}  MAD={mad:6.1f}  rank(min)={rank}/{len(order)}")

    print("\n" + "=" * 92)
    print("SENSITIVITY OF QMDFF@%s TO THE TRAINING CONFORMER" % args.method.upper())
    print("=" * 92)
    print(f"{'train':>6}{'its own rank':>14}{'Spearman':>11}{'Pearson':>10}"
          f"{'MAD':>8}{'spread':>9}{'rank of min':>13}{'MAD near':>11}{'MAD far':>10}")
    print("-" * 92)
    for row in rows:
        print(f"{row['ref']:>6}{row['ref_rank_in_reference']:>10}/{len(frames):<4}"
              f"{row['spearman']:>11.3f}{row['pearson']:>10.3f}{row['mad_kj']:>8.1f}"
              f"{row['spread_kj']:>9.1f}{row['rank_of_reference_minimum']:>8}/{len(frames):<5}"
              f"{row['mad_near_kj']:>11.1f}{row['mad_far_kj']:>10.1f}")
    print("-" * 92)
    if rows:
        rhos = [row["spearman"] for row in rows]
        print(f"Spearman across training choices: {min(rhos):.3f} .. {max(rhos):.3f}"
              f"   (spread {max(rhos) - min(rhos):.3f})")
        near = [r["mad_near_kj"] for r in rows if r["mad_near_kj"] == r["mad_near_kj"]]
        far = [r["mad_far_kj"] for r in rows if r["mad_far_kj"] == r["mad_far_kj"]]
        if near and far:
            print(f"MAD near the training basin {sum(near)/len(near):.1f} kJ/mol  vs  "
                  f"far {sum(far)/len(far):.1f} kJ/mol   "
                  f"(split at the median heavy-atom RMSD)")
    print(f"\n{args.method} spread: {max(rel_ref):.1f} kJ/mol")

    with open(os.path.join(workdir, "reference_sensitivity.json"), "w") as handle:
        json.dump({"ensemble": args.ensemble, "method": args.method, "rows": rows},
                  handle, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
