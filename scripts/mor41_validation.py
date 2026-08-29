#!/usr/bin/env python3
# Claude Generated (Jul 2026)
"""MOR41 benchmark: compare Curcuma native GFN-FF/GFN1/GFN2 vs xtb 6.6.1.

MOR41 (Dohm, Hansen, Steinmetz, Grimme, Checinski, JCTC 2018) is a set of 41
realistic closed-shell metal-organic reactions. Every structure is a neutral
closed-shell singlet (verified via the SI: Table S14 = large S-T gaps, and an
even electron count for all 95 mol.xyz).

For each of the 95 structure folders this runs a single point with both engines
on the provided cartesian mol.xyz (Angstrom, charge 0, multiplicity 1), for the
requested methods:
    curcuma  -sp mol.xyz -method {gfnff|gfn1|gfn2}
    xtb      mol.xyz --sp {--gfnff|--gfn 1|--gfn 2}
Energies are cached in _run/energies.json keyed by (structure, engine, method),
so re-runs and --only subsets never recompute.

The 41 reaction energies dE = sum(coeff*E) [kcal/mol] are then assembled from
reactions.dat (Table S1) per method and compared:
  - Primary:   curcuma vs xtb          (does native reproduce the reference impl?)
  - Secondary: each engine vs DLPNO ref (method accuracy context)

Read-only w.r.t. the test set; writes only under _run/.

Usage:
    python scripts/mor41_validation.py                     # all methods, all reactions
    python scripts/mor41_validation.py --method gfnff      # one method
    python scripts/mor41_validation.py --method gfn2 --only 5 12 41
"""
import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TESTSET = REPO / "test_cases" / "MOR41-testset"
RUNDIR = TESTSET / "_run"
REACTIONS = TESTSET / "reactions.dat"
CURCUMA = REPO / "release" / "curcuma"
XTB = Path.home() / "Downloads" / "xtb-6.6.1" / "bin" / "xtb"

AU2KCAL = 627.509474           # Hartree -> kcal/mol

# curcuma method name -> xtb command-line flag(s) for the same method.
METHODS = {
    "gfnff": ["--gfnff"],
    "gfn1": ["--gfn", "1"],
    "gfn2": ["--gfn", "2"],
}

# ------------------------------------------------------------------ parsing


def parse_curcuma_energy(stdout):
    m = re.search(r"Single Point Energy\s*=\s*(-?\d+\.\d+)\s*Eh", stdout)
    return float(m.group(1)) if m else None


def parse_xtb_energy(stdout):
    # xtb 6.6.1 prints "| TOTAL ENERGY  ...  Eh |" plus a plain
    # "total energy   :  -123.456789 Eh" summary line. Try both.
    m = re.search(r"TOTAL ENERGY[^\n]*?(-?\d+\.\d+)\s*Eh", stdout)
    if not m:
        m = re.search(r"total energy\s*:\s*(-?\d+\.\d+)\s*Eh", stdout)
    if not m:
        for line in stdout.splitlines():
            if "energy" in line.lower():
                mm = re.search(r"(-?\d+\.\d+)\s*Eh", line)
                if mm:
                    return float(mm.group(1))
    return float(m.group(1)) if m else None


# ------------------------------------------------------------------ engines


def run_curcuma(xyz_path, method):
    cmd = [str(CURCUMA), "-sp", str(xyz_path), "-method", method,
           "-charge", "0", "-verbosity", "0", "-no_bmt"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT (curcuma)\n"
    return parse_curcuma_energy(proc.stdout), proc.stdout + proc.stderr


def run_xtb(xyz_path, method, workdir):
    """Copy xyz into a clean workdir and run xtb --sp there (charge 0, singlet)."""
    workdir.mkdir(parents=True, exist_ok=True)
    local_xyz = workdir / xyz_path.name
    shutil.copy(xyz_path, local_xyz)
    cmd = [str(XTB), str(local_xyz), "--sp"] + METHODS[method]
    env = dict(os.environ)
    env["XTBPATH"] = ""   # avoid picking up stray param files
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "4")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800,
                              cwd=str(workdir), env=env)
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT (xtb)\n"
    # clean xtb scratch
    for scratch in ("xtbrestart", "charges", "wbo", "gfnff_topo",
                    "gfnff_charges", "gfnff_adjacency", "xtbtopo.mol",
                    "energy", "gradient", ".xtboptok"):
        sf = workdir / scratch
        if sf.exists():
            try:
                sf.unlink()
            except OSError:
                pass
    return parse_xtb_energy(proc.stdout), proc.stdout + proc.stderr


# ------------------------------------------------------------------ cache


def load_cache():
    if (RUNDIR / "energies.json").exists():
        return json.loads((RUNDIR / "energies.json").read_text())
    return {}


def save_cache(cache):
    RUNDIR.mkdir(parents=True, exist_ok=True)
    (RUNDIR / "energies.json").write_text(json.dumps(cache, indent=1, sort_keys=True))


def cache_key(structure, engine, method):
    return f"{structure}|{engine}|{method}"


def energy_for(structure, engine, method, cache, logdir, recompute=False):
    """Return cached energy or compute it (curcuma or xtb) and cache."""
    key = cache_key(structure, engine, method)
    if not recompute and key in cache:
        return cache[key]
    xyz = TESTSET / structure / "mol.xyz"
    if not xyz.exists():
        print(f"  ! {structure}: no mol.xyz", flush=True)
        cache[key] = None
        return None
    logdir.mkdir(parents=True, exist_ok=True)
    if engine == "cur":
        e, out = run_curcuma(xyz, method)
        (logdir / f"{structure}.{method}.cur.log").write_text(out)
    else:
        e, out = run_xtb(xyz, method, logdir / f"{structure}_{method}_xtb")
        (logdir / f"{structure}.{method}.xtb.log").write_text(out)
    cache[key] = e
    tag = "cur" if engine == "cur" else "xtb"
    print(f"  {structure:8s} {method:5s} {tag}={e}", flush=True)
    return e


# ------------------------------------------------------------------ reactions


def parse_reactions():
    """Parse reactions.dat -> list of dicts {id, ref, terms:{name:coeff}}.

    Products carry +coeff, reactants -coeff. Line form:
        id  ref   c name c name | c name c name
    """
    reactions = []
    for raw in REACTIONS.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        head, _, rest = line.partition("|")
        hp = head.split()
        rid = int(hp[0])
        ref = float(hp[1])
        terms = {}

        def add(tokens, sign):
            i = 0
            while i < len(tokens):
                coeff = int(tokens[i])
                name = tokens[i + 1]
                terms[name] = terms.get(name, 0) + sign * coeff
                i += 2

        add(hp[2:], +1)          # products (left of '|')
        add(rest.split(), -1)    # reactants (right of '|')
        reactions.append({"id": rid, "ref": ref, "terms": terms})
    return reactions


def reaction_energy(terms, engine, method, cache, logdir, recompute):
    """dE in kcal/mol, or None if any term energy is missing."""
    total = 0.0
    for name, coeff in terms.items():
        e = energy_for(name, engine, method, cache, logdir, recompute)
        if e is None:
            return None
        total += coeff * e
    return total * AU2KCAL


# ------------------------------------------------------------------ stats/output


def stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    n = len(vals)
    md = sum(vals) / n                       # signed mean deviation
    mad = sum(abs(v) for v in vals) / n
    mx = max(abs(v) for v in vals)
    rms = math.sqrt(sum(v * v for v in vals) / n)
    return n, md, mad, mx, rms


def write_method_csv(method, rows):
    path = RUNDIR / f"mor41_results_{method}.csv"
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "de_cur", "de_xtb", "de_ref",
                    "d_cur_xtb", "d_cur_ref", "d_xtb_ref", "status"])
        for r in rows:
            def fmt(x):
                return f"{x:.4f}" if x is not None else ""
            w.writerow([r["id"], fmt(r["de_cur"]), fmt(r["de_xtb"]),
                        f"{r['ref']:.2f}", fmt(r["d_cur_xtb"]),
                        fmt(r["d_cur_ref"]), fmt(r["d_xtb_ref"]), r["status"]])
    return path


def write_structure_report(methods, reactions, cache):
    """Per-structure curcuma-vs-xtb energy table (cleanest implementation check).

    Reaction energies mix several structures and are corrupted by any single xtb
    SCF pathology; the per-structure delta E(cur)-E(xtb) [kcal/mol] isolates
    exactly which molecules native curcuma fails to reproduce."""
    structures = sorted({name for rx in reactions for name in rx["terms"]})
    path = RUNDIR / "mor41_per_structure.md"
    with path.open("w") as f:
        f.write("# MOR41 per-structure single-point agreement: curcuma - xtb\n\n")
        f.write("dE = (E_curcuma - E_xtb) * 627.509474 [kcal/mol], same method, "
                "identical mol.xyz, charge 0, singlet.\n\n")
        for method in methods:
            deltas = []
            f.write(f"## {method}\n\n")
            f.write("| structure | E_cur/Eh | E_xtb/Eh | dE(cur-xtb)/kcal |\n")
            f.write("|---|------:|------:|------:|\n")
            recs = []
            for s in structures:
                ec = cache.get(cache_key(s, "cur", method))
                ex = cache.get(cache_key(s, "xtb", method))
                d = (ec - ex) * AU2KCAL if (ec is not None and ex is not None) else None
                recs.append((s, ec, ex, d))
                if d is not None:
                    deltas.append(d)
            for s, ec, ex, d in sorted(recs, key=lambda r: -abs(r[3]) if r[3] is not None else 1):
                f.write("| {s} | {c} | {x} | {d} |\n".format(
                    s=s,
                    c=f"{ec:.6f}" if ec is not None else "-",
                    x=f"{ex:.6f}" if ex is not None else "-",
                    d=f"{d:.2f}" if d is not None else "-"))
            st = stats(deltas)
            if st:
                n, md, mad, mx, rms = st
                f.write(f"\n- n={n}  MD={md:+.3f}  MAD={mad:.3f}  max={mx:.3f}  "
                        f"RMSD={rms:.3f} kcal/mol\n")
            within = sum(1 for d in deltas if abs(d) < 1.0)
            f.write(f"- structures within 1.0 kcal/mol of xtb: {within}/{len(deltas)}\n\n")
    return path


def write_markdown(all_rows):
    path = RUNDIR / "mor41_results.md"
    with path.open("w") as f:
        f.write("# MOR41: Curcuma native vs xtb 6.6.1 (gfnff / gfn1 / gfn2)\n\n")
        f.write("Reaction energies dE = sum(coeff*E) in kcal/mol.\n")
        f.write("d_cur_xtb = curcuma - xtb (primary implementation check); ")
        f.write("d_*_ref = engine - DLPNO-CCSD(T) reference (context).\n\n")
        for method, rows in all_rows.items():
            f.write(f"## {method}\n\n")
            f.write("| id | dE_cur | dE_xtb | dE_ref | d_cur_xtb | d_cur_ref | d_xtb_ref | status |\n")
            f.write("|---|------:|------:|------:|--------:|--------:|--------:|:------|\n")
            for r in rows:
                def g(x):
                    return f"{x:.2f}" if x is not None else "-"
                f.write("| {i} | {c} | {x} | {ref:.2f} | {cx} | {cr} | {xr} | {s} |\n".format(
                    i=r["id"], c=g(r["de_cur"]), x=g(r["de_xtb"]), ref=r["ref"],
                    cx=g(r["d_cur_xtb"]), cr=g(r["d_cur_ref"]), xr=g(r["d_xtb_ref"]),
                    s=r["status"]))
            f.write("\n### Summary (kcal/mol)\n\n")
            for label, key in [("curcuma vs xtb", "d_cur_xtb"),
                               ("curcuma vs ref", "d_cur_ref"),
                               ("xtb vs ref", "d_xtb_ref")]:
                s = stats([r[key] for r in rows])
                if s is None:
                    f.write(f"- {label}: no data\n")
                else:
                    n, md, mad, mx, rms = s
                    f.write(f"- {label}: n={n}  MD={md:+.3f}  MAD={mad:.3f}  "
                            f"max={mx:.3f}  RMSD={rms:.3f}\n")
            outliers = [r for r in rows if r["d_cur_xtb"] is not None
                        and abs(r["d_cur_xtb"]) > 1.0]
            f.write(f"\n### curcuma-vs-xtb outliers (|d| > 1.0): {len(outliers)}\n\n")
            for r in outliers:
                f.write(f"- reaction {r['id']}: d_cur_xtb={r['d_cur_xtb']:.3f} "
                        f"(dE_cur={r['de_cur']:.2f}, dE_xtb={r['de_xtb']:.2f})\n")
            missing = [r["id"] for r in rows if r["status"] != "ok"]
            if missing:
                f.write(f"\n### incomplete reactions (missing energies): {missing}\n")
            f.write("\n")
    return path


# ------------------------------------------------------------------ main


def main():
    ap = argparse.ArgumentParser(description="MOR41 curcuma-vs-xtb validation")
    ap.add_argument("--method", choices=list(METHODS) + ["all"], default="all")
    ap.add_argument("--only", type=int, nargs="*", default=None,
                    help="reaction ids to run (default all 41)")
    ap.add_argument("--recompute", action="store_true",
                    help="ignore cache and recompute energies")
    args = ap.parse_args()

    methods = list(METHODS) if args.method == "all" else [args.method]
    reactions = parse_reactions()
    if args.only:
        want = set(args.only)
        reactions = [r for r in reactions if r["id"] in want]
    RUNDIR.mkdir(parents=True, exist_ok=True)
    cache = load_cache()
    logdir = RUNDIR / "logs"

    all_rows = {}
    for method in methods:
        print(f"\n=== method {method} ===", flush=True)
        rows = []
        for rx in reactions:
            de_cur = reaction_energy(rx["terms"], "cur", method, cache, logdir, args.recompute)
            de_xtb = reaction_energy(rx["terms"], "xtb", method, cache, logdir, args.recompute)
            save_cache(cache)   # persist incrementally (crash-safe)
            ref = rx["ref"]
            d_cur_xtb = (de_cur - de_xtb) if (de_cur is not None and de_xtb is not None) else None
            d_cur_ref = (de_cur - ref) if de_cur is not None else None
            d_xtb_ref = (de_xtb - ref) if de_xtb is not None else None
            status = "ok" if (de_cur is not None and de_xtb is not None) else "FAIL"
            rows.append({"id": rx["id"], "ref": ref, "de_cur": de_cur, "de_xtb": de_xtb,
                         "d_cur_xtb": d_cur_xtb, "d_cur_ref": d_cur_ref,
                         "d_xtb_ref": d_xtb_ref, "status": status})
            print(f"[{rx['id']:2d}] dE_cur={de_cur}  dE_xtb={de_xtb}  "
                  f"d_cur_xtb={d_cur_xtb}", flush=True)
        csvp = write_method_csv(method, rows)
        all_rows[method] = rows
        s = stats([r["d_cur_xtb"] for r in rows])
        if s:
            n, md, mad, mx, rms = s
            print(f"--- {method} curcuma vs xtb (n={n}): "
                  f"MD={md:+.3f} MAD={mad:.3f} max={mx:.3f} RMSD={rms:.3f} kcal/mol")
        print(f"    wrote {csvp}")

    mdp = write_markdown(all_rows)
    sp = write_structure_report(methods, reactions, cache)
    save_cache(cache)
    print(f"\nWrote {mdp}\nWrote {sp}")


if __name__ == "__main__":
    main()
