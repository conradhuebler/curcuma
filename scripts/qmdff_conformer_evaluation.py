#!/usr/bin/env python3
"""
QMDFF@GFN2 use-case evaluation on a conformer ensemble.

Background (Obsidian: "Beschreibung und Rangfolge von Konformeren"): a conformer
search explores on a cheap surface and ranks on an accurate one. Measured on a
107-atom peptide, GFN-FF and GFN2 correlate at r = -0.32 .. -0.46 — the exploration
surface actively misleads the ranking. QMDFF@GFN2, a force field fitted to the
ranking method itself, is the structural fix: exploration and ranking then live on
the same surface, and the force-field term decomposition finally describes the
*right* surface.

This script measures whether that holds, on a real ensemble:

  1. GFN2 single points for every conformer          -> the reference ranking
  2. GFN-FF single points                            -> the status-quo exploration surface
  3. one QMDFF fit at a reference conformer (GFN2 Hessian + charges)
  4. QMDFF single points for every conformer, with that ONE fixed parameter set
  5. correlation / ranking statistics, and the exact variance decomposition
     Var(E) = sum_i Cov(E_i, E) over the force-field terms

Point 4 deliberately uses a single parameter set for the whole ensemble, because
that is what a conformer search would do — and it is exactly where QMDFF is
expected to degrade, since it is fitted in one minimum (the note's first caveat).

Usage:
    scripts/qmdff_conformer_evaluation.py ensemble.xyz [options]

    --curcuma PATH     curcuma binary (default: release/curcuma)
    --workdir DIR      scratch directory (default: ./qmdff_eval)
    --frames N         use only the first N conformers (default: all)
    --ref-frame I      reference conformer for the fit (default: the GFN2 minimum)
    --threads N        threads per curcuma call (default: 4)
    --skip-fit         reuse an existing fit in the workdir

Claude Generated (August 2026)
Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys

HARTREE_TO_KJ = 2625.4996394799

# Terms as reported by FFWorkspace at verbosity 3 ("=== CPU ENERGY TERMS ===")
TERMS = ["bond", "angle", "dihedral", "inversion", "stors", "batm", "atm",
         "disp", "brep", "nbrep", "coulomb", "hbond", "xbond"]

TERM_LABEL = {
    "bond": "bond stretch", "angle": "angle bend", "dihedral": "torsion",
    "inversion": "out-of-plane", "stors": "triple-bond torsion",
    "batm": "bonded ATM", "atm": "ATM", "disp": "dispersion",
    "brep": "bonded repulsion", "nbrep": "repulsion",
    "coulomb": "electrostatics", "hbond": "hydrogen bonds", "xbond": "halogen bonds",
}


# --------------------------------------------------------------------------- io

def read_xyz_frames(path):
    frames = []
    with open(path) as handle:
        lines = handle.read().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        natoms = int(lines[i].split()[0])
        frames.append(lines[i:i + natoms + 2])
        i += natoms + 2
    return frames


def write_xyz(path, frame):
    with open(path, "w") as handle:
        handle.write("\n".join(frame) + "\n")


# ---------------------------------------------------------------------- running

ENERGY_RE = re.compile(r"Single Point Energy\s*=\s*(-?\d+\.?\d*(?:[eE][-+]?\d+)?)")
TERM_RE = re.compile(r"^\s*(\w+)\s*=\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)\s*$")


def run_single_point(curcuma, xyz, method, threads, want_terms=False, extra=None):
    """Return (energy_hartree, {term: energy}) or (None, {}) on failure."""
    cmd = [curcuma, "-sp", xyz, "-method", method, "-no_bmt",
           "-threads", str(threads), "-v", "3" if want_terms else "1"]
    if extra:
        cmd += extra
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=3600).stdout
    except subprocess.TimeoutExpired:
        return None, {}

    energy = None
    for match in ENERGY_RE.finditer(out):
        energy = float(match.group(1))

    terms = {}
    if want_terms and "CPU ENERGY TERMS" in out:
        block = out.split("CPU ENERGY TERMS")[-1].split("CPU ENERGY END")[0]
        for line in block.splitlines():
            match = TERM_RE.match(line)
            if match and match.group(1) in TERMS:
                terms[match.group(1)] = float(match.group(2))
    return energy, terms


# ------------------------------------------------------------------- statistics

def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            mean_rank = (i + j) / 2.0
            for k in range(i, j + 1):
                out[order[k]] = mean_rank
            i = j + 1
        return out
    return pearson(ranks(x), ranks(y))


def relative(values):
    lo = min(values)
    return [(v - lo) * HARTREE_TO_KJ for v in values]


def variance_decomposition(term_series, total):
    """Exact decomposition Var(E) = sum_i Cov(E_i, E); returns share per term."""
    n = len(total)
    mean_total = sum(total) / n
    var_total = sum((t - mean_total) ** 2 for t in total) / n
    shares = {}
    for term, series in term_series.items():
        mean_term = sum(series) / n
        cov = sum((a - mean_term) * (b - mean_total)
                  for a, b in zip(series, total)) / n
        shares[term] = cov / var_total if var_total > 0 else float("nan")
    return shares, var_total


# ------------------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ensemble")
    parser.add_argument("--curcuma", default="release/curcuma")
    parser.add_argument("--workdir", default="qmdff_eval")
    parser.add_argument("--frames", type=int, default=0)
    parser.add_argument("--ref-frame", type=int, default=-1)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--skip-fit", action="store_true")
    args = parser.parse_args()

    curcuma = os.path.abspath(args.curcuma)
    workdir = os.path.abspath(args.workdir)
    os.makedirs(workdir, exist_ok=True)

    frames = read_xyz_frames(args.ensemble)
    if args.frames > 0:
        frames = frames[:args.frames]
    natoms = int(frames[0][0].split()[0])
    print(f"ensemble: {args.ensemble} — {len(frames)} conformers, {natoms} atoms\n")

    paths = []
    for index, frame in enumerate(frames):
        path = os.path.join(workdir, f"conf_{index:04d}.xyz")
        write_xyz(path, frame)
        paths.append(path)

    # --- 1) GFN2 reference ranking -----------------------------------------
    print("[1/4] GFN2 single points (reference ranking) ...", flush=True)
    e_gfn2 = []
    for index, path in enumerate(paths):
        energy, _ = run_single_point(curcuma, path, "gfn2", args.threads)
        if energy is None:
            print(f"  conformer {index}: GFN2 FAILED", file=sys.stderr)
            return 1
        e_gfn2.append(energy)
        print(f"\r  {index + 1}/{len(paths)}", end="", flush=True)
    print()

    ref = args.ref_frame if args.ref_frame >= 0 else e_gfn2.index(min(e_gfn2))
    print(f"  reference conformer for the fit: {ref} (GFN2 minimum)\n")

    # --- 2) GFN-FF ----------------------------------------------------------
    print("[2/4] GFN-FF single points (status-quo exploration surface) ...", flush=True)
    e_gfnff, terms_gfnff = [], {t: [] for t in TERMS}
    for index, path in enumerate(paths):
        energy, terms = run_single_point(curcuma, path, "gfnff", args.threads, want_terms=True)
        e_gfnff.append(energy if energy is not None else float("nan"))
        for term in TERMS:
            terms_gfnff[term].append(terms.get(term, 0.0))
        print(f"\r  {index + 1}/{len(paths)}", end="", flush=True)
    print()

    # --- 3) the QMDFF fit at the reference conformer -------------------------
    fit_source = os.path.join(workdir, "qmdff_reference.xyz")
    fitted_param = fit_source.replace(".xyz", ".param.json")
    provenance = os.path.join(workdir, "fit_provenance.json")
    if not (args.skip_fit and os.path.exists(fitted_param)):
        print(f"\n[3/4] QMDFF fit at conformer {ref} (GFN2 Hessian, {6 * natoms} gradient calls) ...",
              flush=True)
        shutil.copy(paths[ref], fit_source)
        for stale in (fitted_param, "scf.json", "hessian.json"):
            if os.path.exists(stale):
                os.remove(stale)
        result = subprocess.run(
            [curcuma, "-qmdfffit", fit_source, "-method", "gfn2",
             "-threads", str(args.threads), "-v", "2"],
            capture_output=True, text=True, cwd=workdir, timeout=86400)
        for line in result.stdout.splitlines():
            if any(key in line for key in ("qmdff_fit", "verified", "Verification", "WARN")):
                print("   ", line.strip())
        if not os.path.exists(fitted_param):
            print("QMDFF fit produced no parameter file", file=sys.stderr)
            return 1
        with open(provenance, "w") as handle:
            json.dump({"fitted_at_frame": ref, "n_frames_when_fitted": len(frames)}, handle)
    else:
        # --skip-fit reuses a parameter set that may have been fitted at a DIFFERENT
        # conformer (e.g. from a shorter previous run). Report what was actually used —
        # a force field fitted away from the ensemble minimum is the realistic case, but
        # it must not be silently mislabelled.
        if os.path.exists(provenance):
            with open(provenance) as handle:
                ref = json.load(handle).get("fitted_at_frame", ref)
        print(f"\n[3/4] reusing existing fit {fitted_param} (fitted at conformer {ref})")

    # --- 4) QMDFF@GFN2 with that ONE parameter set --------------------------
    # The auto-parameter cache keys on <basename>.param.json, so copying the fitted
    # set next to each conformer applies the same force field to the whole ensemble —
    # which is what a conformer search would do, and where QMDFF is expected to degrade.
    print("\n[4/4] QMDFF@GFN2 single points with the fixed fitted parameter set ...", flush=True)
    e_qmdff, terms_qmdff = [], {t: [] for t in TERMS}
    for index, path in enumerate(paths):
        shutil.copy(fitted_param, path.replace(".xyz", ".param.json"))
        energy, terms = run_single_point(curcuma, path, "qmdff", args.threads, want_terms=True)
        e_qmdff.append(energy if energy is not None else float("nan"))
        for term in TERMS:
            terms_qmdff[term].append(terms.get(term, 0.0))
        print(f"\r  {index + 1}/{len(paths)}", end="", flush=True)
    print("\n")

    # --- statistics ----------------------------------------------------------
    rel_gfn2 = relative(e_gfn2)
    report = {"ensemble": args.ensemble, "n_conformers": len(frames),
              "n_atoms": natoms, "fitted_at_frame": ref,
              "gfn2_minimum_frame": e_gfn2.index(min(e_gfn2))}

    print("=" * 78)
    print("AGREEMENT WITH THE GFN2 RANKING SURFACE")
    print("=" * 78)
    print(f"{'method':<14}{'Pearson r':>12}{'Spearman rho':>14}"
          f"{'MAD [kJ/mol]':>14}{'spread':>10}{'rank of min':>13}")
    print("-" * 78)

    for name, energies in (("GFN-FF", e_gfnff), ("QMDFF@GFN2", e_qmdff)):
        valid = [i for i, v in enumerate(energies) if v == v]
        if len(valid) < 3:
            print(f"{name:<14} insufficient data")
            continue
        x = [rel_gfn2[i] for i in valid]
        y = relative([energies[i] for i in valid])
        r = pearson(x, y)
        rho = spearman(x, y)
        mad = sum(abs(a - b) for a, b in zip(x, y)) / len(x)
        # Where does the true GFN2 minimum land in this method's own ranking?
        order = sorted(valid, key=lambda i: energies[i])
        true_min = e_gfn2.index(min(e_gfn2))
        rank = order.index(true_min) + 1 if true_min in order else -1
        print(f"{name:<14}{r:>12.3f}{rho:>14.3f}{mad:>14.1f}"
              f"{max(y):>10.1f}{rank:>7d}/{len(order):<5d}")
        report[name] = {"pearson": r, "spearman": rho, "mad_kj": mad,
                        "spread_kj": max(y), "rank_of_gfn2_minimum": rank}

    print(f"\nGFN2 conformer spread: {max(rel_gfn2):.1f} kJ/mol")
    report["gfn2_spread_kj"] = max(rel_gfn2)

    # --- term interpretability -----------------------------------------------
    print("\n" + "=" * 78)
    print("TERM VARIANCE DECOMPOSITION   Var(E) = sum_i Cov(E_i, E)")
    print("Which force-field terms actually carry the conformational energy spread?")
    print("=" * 78)
    print(f"{'term':<22}{'GFN-FF':>14}{'QMDFF@GFN2':>16}")
    print("-" * 78)

    shares_gfnff, var_gfnff = variance_decomposition(terms_gfnff, e_gfnff)
    shares_qmdff, var_qmdff = variance_decomposition(terms_qmdff, e_qmdff)
    report["variance_shares"] = {"gfnff": shares_gfnff, "qmdff": shares_qmdff}

    ordering = sorted(TERMS, key=lambda t: -abs(shares_qmdff.get(t, 0.0)))
    for term in ordering:
        a, b = shares_gfnff.get(term, 0.0), shares_qmdff.get(term, 0.0)
        if abs(a) < 5e-4 and abs(b) < 5e-4:
            continue
        print(f"{TERM_LABEL[term]:<22}{100 * a:>13.1f}%{100 * b:>15.1f}%")
    print("-" * 78)
    print(f"{'sum (must be 100%)':<22}{100 * sum(shares_gfnff.values()):>13.1f}%"
          f"{100 * sum(shares_qmdff.values()):>15.1f}%")
    print(f"{'std. dev. of E':<22}{math.sqrt(var_gfnff) * HARTREE_TO_KJ:>12.1f} "
          f"{math.sqrt(var_qmdff) * HARTREE_TO_KJ:>15.1f}   kJ/mol")

    with open(os.path.join(workdir, "evaluation.json"), "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"\nwritten: {os.path.join(workdir, 'evaluation.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
