#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""GMTKN55 benchmark: compare Curcuma native GFN-FF/GFN1/GFN2 vs xtb.

Runs a single point with both engines on every struc.xyz of the fetched
GMTKN55 set (scripts/fetch_testset.py fetch gmtkn55; test_cases/GMTKN55-testset,
54 subsets / 2462 structures) and compares curcuma vs xtb energies directly
(same philosophy as scripts/mor41_validation.py / scripts/s30l_gfnff_compare.py:
this is a reproduces-the-reference-implementation check, not a GMTKN55 WTMAD-2
accuracy run - the .res files are tmer2++ shell scripts, a separate and much
larger undertaking, out of scope here).

KNOWN LIMITATION (verified Sep 2026): curcuma's `-sp` command has no working
CLI path to request open-shell (UHF) occupation for the native gfn1/gfn2
solver - `-spin N` only sets inert Molecule metadata (never read by any
energy_calculators/ code) and `-multi N` is not routed to the top-level
`controller["multi"]` EnergyCalculator actually reads for `-sp`. Verified on
RC21/me (CH3 radical, GFN2): curcuma gives the IDENTICAL energy for
`-spin 0` and `-spin 1` (-3.56397536 Eh, always closed-shell), vs xtb's true
UHF doublet energy -3.56265832 Eh (0.83 kcal/mol off, small but silently
wrong here and can be much larger for other radicals). Because of this,
gfn1/gfn2 SKIP every structure with a nonzero .UHF file rather than silently
reporting a closed-shell number as if it were the requested open-shell one;
gfnff runs everything since GFN-FF (force field, no explicit electronic
open-shell term) does not depend on this path.

Energies are cached in _run/energies.json keyed by (subset, name, engine,
method), so re-runs and --subset/--limit subsets never recompute.

Read-only w.r.t. the test set; writes only under _run/.

Usage:
    python scripts/gmtkn55_compare.py --subset ACONF S22       # smoke test
    python scripts/gmtkn55_compare.py --method gfnff           # one method, all subsets
    python scripts/gmtkn55_compare.py --limit 200              # first 200 structures
    python scripts/gmtkn55_compare.py                          # full sweep (slow, hours)
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
TESTSET = REPO / "test_cases" / "GMTKN55-testset"
RUNDIR = TESTSET / "_run"
CURCUMA = REPO / "release" / "curcuma"

AU2KCAL = 627.509474           # Hartree -> kcal/mol
NON_SUBSET_DIRS = {"_utils", "_results", ".git"}

METHODS = {
    "gfnff": ["--gfnff"],
    "gfn1": ["--gfn", "1"],
    "gfn2": ["--gfn", "2"],
}
OPEN_SHELL_UNSUPPORTED = {"gfn1", "gfn2"}   # see module docstring


def find_xtb():
    env = os.environ.get("XTB_BIN")
    if env and Path(env).exists():
        return Path(env)
    which = shutil.which("xtb")
    if which:
        return Path(which)
    for cand in ("/opt/xtb/bin/xtb",
                 Path.home() / "Downloads" / "xtb-dist" / "bin" / "xtb",
                 Path.home() / "Downloads" / "xtb-6.6.1" / "bin" / "xtb"):
        p = Path(cand)
        if p.exists():
            return p
    return None


# ------------------------------------------------------------------ discovery


def discover_structures(subset_filter=None):
    """List (subset, name, xyz_path, charge, uhf) for every struc.xyz found."""
    structures = []
    subsets = sorted(p.name for p in TESTSET.iterdir()
                      if p.is_dir() and p.name not in NON_SUBSET_DIRS)
    if subset_filter:
        want = set(subset_filter)
        subsets = [s for s in subsets if s in want]
    for subset in subsets:
        for mol_dir in sorted((TESTSET / subset).iterdir()):
            xyz = mol_dir / "struc.xyz"
            if not xyz.exists():
                continue
            charge = 0
            chrg_f = mol_dir / ".CHRG"
            if chrg_f.exists():
                charge = int(chrg_f.read_text().strip())
            uhf = 0
            uhf_f = mol_dir / ".UHF"
            if uhf_f.exists():
                text = uhf_f.read_text().strip()
                uhf = int(text) if text else 0
            structures.append((subset, mol_dir.name, xyz, charge, uhf))
    return structures


# ------------------------------------------------------------------ parsing


def parse_curcuma_energy(stdout):
    m = re.search(r"Single Point Energy\s*=\s*(-?\d+\.\d+)\s*Eh", stdout)
    return float(m.group(1)) if m else None


def parse_xtb_energy(stdout):
    """Total energy from an xtb run, or None if xtb did not produce one.

    NaN handling (Sep 2026): GFN-FF's bond charge factor overflows to NaN for strongly
    ionic bonds, and xtb then prints "TOTAL ENERGY NaN Eh" while still printing finite
    numbers for the individual terms above it. The old loose fallback scanned every line
    containing "energy" and took the first number followed by "Eh", which for such a run
    silently returned the ANGLE energy as the total — that is how PX13/hf_4_ts entered
    this comparison as -0.000479853103 Eh instead of as a failure. Detect NaN first, and
    keep the fallback anchored to a total-energy line.
    """
    if re.search(r"(?i)total energy\s*:?\s*N[Aa]N", stdout):
        return None
    m = re.search(r"TOTAL ENERGY[^\n]*?(-?\d+\.\d+)\s*Eh", stdout)
    if not m:
        m = re.search(r"total energy\s*:?\s*(-?\d+\.\d+)\s*Eh", stdout)
    return float(m.group(1)) if m else None


# ------------------------------------------------------------------ engines


def run_curcuma(xyz_path, method, charge, timeout):
    cmd = [str(CURCUMA), "-sp", str(xyz_path), "-method", method,
           "-charge", str(charge), "-verbosity", "0", "-no_bmt"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT (curcuma)\n"
    return parse_curcuma_energy(proc.stdout), proc.stdout + proc.stderr


def run_xtb(xtb_bin, xyz_path, method, charge, uhf, workdir, timeout):
    workdir.mkdir(parents=True, exist_ok=True)
    local_xyz = workdir / "mol.xyz"
    shutil.copy(xyz_path, local_xyz)
    (workdir / ".CHRG").write_text(f"{charge:+d}\n")
    if uhf:
        (workdir / ".UHF").write_text(f"{uhf:d}\n")
    cmd = [str(xtb_bin), str(local_xyz), "--sp"] + METHODS[method]
    env = dict(os.environ)
    env["XTBPATH"] = ""
    env["OMP_NUM_THREADS"] = env.get("OMP_NUM_THREADS", "4")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                              cwd=str(workdir), env=env)
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT (xtb)\n"
    for scratch in ("xtbrestart", "charges", "wbo", "gfnff_topo", "gfnff_charges",
                    "gfnff_adjacency", "xtbtopo.mol", "energy", "gradient", ".xtboptok"):
        sf = workdir / scratch
        if sf.exists():
            try:
                sf.unlink()
            except OSError:
                pass
    return parse_xtb_energy(proc.stdout), proc.stdout + proc.stderr


# ------------------------------------------------------------------ cache


def load_cache():
    p = RUNDIR / "energies.json"
    return json.loads(p.read_text()) if p.exists() else {}


def save_cache(cache):
    RUNDIR.mkdir(parents=True, exist_ok=True)
    (RUNDIR / "energies.json").write_text(json.dumps(cache, indent=1, sort_keys=True))


def cache_key(subset, name, engine, method):
    return f"{subset}/{name}|{engine}|{method}"


# ------------------------------------------------------------------ stats


def stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None
    n = len(vals)
    md = sum(vals) / n
    mad = sum(abs(v) for v in vals) / n
    mx = max(abs(v) for v in vals)
    rms = math.sqrt(sum(v * v for v in vals) / n)
    return n, md, mad, mx, rms


# ------------------------------------------------------------------ main


def main():
    ap = argparse.ArgumentParser(description="GMTKN55 curcuma-vs-xtb validation",
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", choices=list(METHODS) + ["all"], default="all")
    ap.add_argument("--subset", nargs="*", default=None,
                     help="restrict to these subset names (default: all 54)")
    ap.add_argument("--limit", type=int, default=None,
                     help="only the first N discovered structures (smoke test)")
    ap.add_argument("--recompute", action="store_true")
    ap.add_argument("--xtb", type=Path, default=None, help="override xtb binary path")
    ap.add_argument("--timeout", type=int, default=600)
    args = ap.parse_args()

    if not CURCUMA.exists():
        raise SystemExit(f"curcuma binary not found at {CURCUMA} - build release/ first")
    xtb_bin = args.xtb or find_xtb()
    if not xtb_bin or not xtb_bin.exists():
        raise SystemExit("xtb binary not found - set XTB_BIN or pass --xtb PATH")
    print(f"xtb: {xtb_bin}", flush=True)

    methods = list(METHODS) if args.method == "all" else [args.method]
    structures = discover_structures(args.subset)
    if not structures:
        raise SystemExit(f"no struc.xyz found under {TESTSET} - "
                          f"fetch it first: python scripts/fetch_testset.py fetch gmtkn55")
    if args.limit:
        structures = structures[:args.limit]
    print(f"{len(structures)} structures across "
          f"{len({s for s, *_ in structures})} subset(s)", flush=True)

    RUNDIR.mkdir(parents=True, exist_ok=True)
    cache = load_cache()
    logdir = RUNDIR / "logs"

    for method in methods:
        print(f"\n=== method {method} ===", flush=True)
        rows = []
        n_skipped = 0
        by_subset = {}
        for i, (subset, name, xyz, charge, uhf) in enumerate(structures):
            skip = uhf and method in OPEN_SHELL_UNSUPPORTED
            key_cur = cache_key(subset, name, "cur", method)
            key_xtb = cache_key(subset, name, "xtb", method)

            if skip:
                n_skipped += 1
                e_cur, e_xtb = None, None
                status = "skip-open-shell"
            else:
                if not args.recompute and key_cur in cache:
                    e_cur = cache[key_cur]
                else:
                    e_cur, out = run_curcuma(xyz, method, charge, args.timeout)
                    cache[key_cur] = e_cur
                    (logdir / subset).mkdir(parents=True, exist_ok=True)
                    (logdir / subset / f"{name}.{method}.cur.log").write_text(out)
                if not args.recompute and key_xtb in cache:
                    e_xtb = cache[key_xtb]
                else:
                    e_xtb, out = run_xtb(xtb_bin, xyz, method, charge, uhf,
                                          logdir / subset / f"{name}_{method}_xtb",
                                          args.timeout)
                    cache[key_xtb] = e_xtb
                status = "ok" if (e_cur is not None and e_xtb is not None) else "FAIL"

            d = ((e_cur - e_xtb) * AU2KCAL) if (e_cur is not None and e_xtb is not None) else None
            rows.append({"subset": subset, "name": name, "atoms": None, "charge": charge,
                         "uhf": uhf, "e_cur": e_cur, "e_xtb": e_xtb, "d_kcal": d,
                         "status": status})
            by_subset.setdefault(subset, []).append(d)

            if (i + 1) % 50 == 0 or i == len(structures) - 1:
                save_cache(cache)
                print(f"  [{i + 1}/{len(structures)}] {subset}/{name} d={d} "
                      f"(skipped so far: {n_skipped})", flush=True)

        save_cache(cache)
        write_csv(method, rows)
        write_subset_summary(method, rows, by_subset)
        s = stats([r["d_kcal"] for r in rows])
        if s:
            n, md, mad, mx, rms = s
            print(f"--- {method} curcuma vs xtb (n={n}, skipped={n_skipped}): "
                  f"MD={md:+.3f} MAD={mad:.3f} max={mx:.3f} RMSD={rms:.3f} kcal/mol")
        else:
            print(f"--- {method}: no comparable structures (skipped={n_skipped})")


def write_csv(method, rows):
    path = RUNDIR / f"gmtkn55_results_{method}.csv"
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["subset", "name", "charge", "uhf", "e_cur_Eh", "e_xtb_Eh",
                    "d_kcal", "status"])
        for r in rows:
            def fmt(x, nd=6):
                return f"{x:.{nd}f}" if x is not None else ""
            w.writerow([r["subset"], r["name"], r["charge"], r["uhf"],
                        fmt(r["e_cur"]), fmt(r["e_xtb"]), fmt(r["d_kcal"], 3), r["status"]])
    print(f"    wrote {path}")
    return path


def write_subset_summary(method, rows, by_subset):
    path = RUNDIR / f"gmtkn55_summary_{method}.md"
    with path.open("w") as f:
        f.write(f"# GMTKN55: curcuma native vs xtb ({method})\n\n")
        f.write("dE = (E_curcuma - E_xtb) * 627.509474 [kcal/mol], per-structure single "
                "point on the published GMTKN55 geometry, identical charge/UHF. This is a "
                "reproduces-the-reference-implementation check, not a GMTKN55 WTMAD-2 "
                "accuracy statistic.\n\n")
        if method in OPEN_SHELL_UNSUPPORTED:
            f.write("**Open-shell structures (.UHF != 0) are skipped** - curcuma's `-sp` "
                    "has no working path to request UHF occupation for gfn1/gfn2 (see "
                    "script docstring); reporting a closed-shell number for those would be "
                    "silently wrong.\n\n")
        f.write("| subset | n | skipped | MD | MAD | max | RMSD |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|\n")
        for subset in sorted(by_subset):
            deltas = by_subset[subset]
            skipped = sum(1 for d in deltas if d is None)
            s = stats(deltas)
            if s:
                n, md, mad, mx, rms = s
                f.write(f"| {subset} | {n} | {skipped} | {md:+.3f} | {mad:.3f} | "
                        f"{mx:.3f} | {rms:.3f} |\n")
            else:
                f.write(f"| {subset} | 0 | {skipped} | - | - | - | - |\n")
        all_deltas = [d for ds in by_subset.values() for d in ds]
        s = stats(all_deltas)
        if s:
            n, md, mad, mx, rms = s
            f.write(f"\n**Overall**: n={n}  MD={md:+.3f}  MAD={mad:.3f}  "
                    f"max={mx:.3f}  RMSD={rms:.3f} kcal/mol\n")
        outliers = [r for r in rows if r["d_kcal"] is not None and abs(r["d_kcal"]) > 5.0]
        outliers.sort(key=lambda r: -abs(r["d_kcal"]))
        f.write(f"\n## Outliers (|d| > 5.0 kcal/mol): {len(outliers)}\n\n")
        for r in outliers[:50]:
            f.write(f"- {r['subset']}/{r['name']}: d={r['d_kcal']:+.2f} kcal/mol "
                    f"(charge={r['charge']}, uhf={r['uhf']})\n")
    print(f"    wrote {path}")
    return path


if __name__ == "__main__":
    main()
