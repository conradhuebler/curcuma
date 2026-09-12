#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP0b
"""Per-term GFN-FF decomposition of GMTKN55 barrier heights (and any other subset).

For every reaction of the selected subsets the script runs a curcuma single point
with the verbosity-2 term table on each involved structure (cached per structure
in _run/terms_<method>.json), forms the reaction sum of every term, and reports
per subset:

  * MAD / MSE of the total against the published reference (plus the cached gfn2
    total as a diagnostic column, no term table exists for it),
  * the mean signed contribution of every GFN-FF term to the reaction energy,
  * which term's contribution best explains the error across the subset's
    reactions (Pearson r and the least-squares slope of error vs contribution).

There is no reference decomposition, so "the term that carries the error" is a
correlation statement, not a proof; it tells where to look first.

Every structure is evaluated in its own scratch directory: all GMTKN55 files are
named struc.xyz and GFN-FF writes struc.topo.json next to the file it reads, so
a shared directory would silently reuse another molecule's topology.

Usage:
    python scripts/revgfnff_barrier_terms.py                      # the 7 barrier subsets
    python scripts/revgfnff_barrier_terms.py --subset PX13 BH76
    python scripts/revgfnff_barrier_terms.py --method gfnff --recompute
"""
import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gmtkn55_reactions import (AU2KCAL, CHEM_CLASS, RESULTS, RUNDIR, TESTSET,  # noqa: E402
                               evaluate, load_energies, load_reactions, reaction_flags,
                               structure_dir, structure_meta)

REPO = Path(__file__).resolve().parents[1]
CURCUMA = REPO / "release" / "curcuma"
BARRIER_SUBSETS = CHEM_CLASS["barriers"]

# curcuma verbosity-2 [RESULT] row name -> short term key (order = report order)
TERM_CANON = OrderedDict([
    ("Bond", "bond"), ("Angle", "angle"), ("Dihedral", "torsion"), ("Inversion", "inversion"),
    ("sTors", "stors"), ("Repulsion (bonded)", "rep_b"), ("Repulsion (nonbond)", "rep_nb"),
    ("Repulsion", "rep"), ("Coulomb", "coulomb"), ("Dispersion", "disp"), ("H-bonds", "hb"),
    ("X-bonds", "xb"), ("ATM (3-body)", "atm"), ("BATM", "batm"), ("OverCoord", "over"),
])
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def curcuma_terms(xyz, charge, uhf, method, timeout=300):
    """{'total': Eh, terms: {key: Eh}} from a verbosity-2 single point in a private scratch dir."""
    work = Path(tempfile.mkdtemp(prefix="revgfnff_terms_"))
    try:
        local = work / "struc.xyz"
        shutil.copy(xyz, local)
        cmd = [str(CURCUMA), "-sp", "struc.xyz", "-method", method, "-charge", str(charge),
               "-spin", str(uhf), "-verbosity", "2", "-no_bmt", "-threads", "1"]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=work)
        except subprocess.TimeoutExpired:
            return {"total": None, "terms": {}, "status": "TIMEOUT"}
        out = ANSI.sub("", proc.stdout)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    terms, total = {}, None
    for line in out.splitlines():
        m = re.match(r"\[RESULT\]\s*(\w[\w ()\-]*?)\s+([+-]?\d+\.\d+)(\s|$)", line)
        if not m or "ms" in line or "wall" in line or "%" in line:
            continue
        name, val = m.group(1).strip(), float(m.group(2))
        if name.lower() == "total":
            total = val
        elif name in TERM_CANON:
            terms[TERM_CANON[name]] = val
    m = re.search(r"Single Point Energy\s*=\s*(-?\d+\.\d+)\s*Eh", out)
    if m:
        total = float(m.group(1))
    status = "ok" if total is not None and terms else "NO_TERMS"
    if total is not None and (math.isnan(total) or any(math.isnan(v) for v in terms.values())):
        status = "NAN"
    return {"total": total, "terms": terms, "status": status}


def load_term_cache(method):
    p = RUNDIR / f"terms_{method}.json"
    return json.loads(p.read_text()) if p.exists() else {}


def save_term_cache(method, cache):
    RUNDIR.mkdir(exist_ok=True)
    (RUNDIR / f"terms_{method}.json").write_text(json.dumps(cache, indent=1, sort_keys=True))


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None, None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return None, None
    return sxy / math.sqrt(sxx * syy), sxy / sxx


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subset", nargs="*", default=BARRIER_SUBSETS)
    ap.add_argument("--method", default="gfnff")
    ap.add_argument("--recompute", action="store_true")
    ap.add_argument("--timeout", type=int, default=300)
    args = ap.parse_args()

    reactions, _ = load_reactions(args.subset)
    cache = load_term_cache(args.method)
    energies_gfn2 = load_energies("cur", "gfn2")
    needed = OrderedDict()
    for rx in reactions:
        for s in rx.species:
            needed[f"{structure_dir(rx.subset)}/{s}"] = (rx.subset, s)
    todo = [k for k in needed if args.recompute or k not in cache or cache[k].get("status") != "ok"]
    print(f"{len(reactions)} reactions, {len(needed)} structures, {len(todo)} to compute ({args.method})")
    for i, key in enumerate(todo, 1):
        subset, name = needed[key]
        charge, uhf = structure_meta(subset, name)
        xyz = TESTSET / structure_dir(subset) / name / "struc.xyz"
        cache[key] = curcuma_terms(xyz, charge, uhf, args.method, args.timeout)
        if i % 20 == 0 or i == len(todo):
            save_term_cache(args.method, cache)
            print(f"  {i}/{len(todo)} ({key}: {cache[key]['status']})")
    save_term_cache(args.method, cache)

    term_keys = [k for k in TERM_CANON.values()
                 if any(k in cache.get(key, {}).get("terms", {}) for key in needed)]
    rows, by_subset = [], defaultdict(list)
    for rx in reactions:
        keys = [f"{structure_dir(rx.subset)}/{s}" for s in rx.species]
        ok = all(cache.get(k, {}).get("status") == "ok" for k in keys)
        contrib = {}
        total = None
        if ok:
            total = sum(c * cache[k]["total"] for k, c in zip(keys, rx.coeffs)) * AU2KCAL
            for t in term_keys:
                contrib[t] = sum(c * cache[k]["terms"].get(t, 0.0) for k, c in zip(keys, rx.coeffs)) * AU2KCAL
        e_gfn2 = evaluate(rx, energies_gfn2)
        charged, open_shell = reaction_flags(rx)
        row = {"subset": rx.subset, "reaction": rx.label, "charged": int(charged), "open_shell": int(open_shell),
               "ref": rx.ref, "total": total, "error": None if total is None else total - rx.ref,
               "gfn2_error": None if e_gfn2 is None else e_gfn2 - rx.ref, "pbeh3c_error": None if rx.pbeh3c is None else rx.pbeh3c - rx.ref,
               "contrib": contrib}
        rows.append(row)
        by_subset[rx.subset].append(row)

    RESULTS.mkdir(exist_ok=True)
    out_csv = RESULTS / f"barrier_terms_{args.method}.csv"
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["subset", "reaction", "charged", "open_shell", "ref_kcal", f"{args.method}_kcal", "error_kcal",
                    "gfn2_error_kcal", "pbeh3c_error_kcal"] + [f"{t}_kcal" for t in term_keys])
        for r in rows:
            w.writerow([r["subset"], r["reaction"], r["charged"], r["open_shell"], f"{r['ref']:.3f}",
                        "" if r["total"] is None else f"{r['total']:.3f}",
                        "" if r["error"] is None else f"{r['error']:.3f}",
                        "" if r["gfn2_error"] is None else f"{r['gfn2_error']:.3f}",
                        "" if r["pbeh3c_error"] is None else f"{r['pbeh3c_error']:.3f}"]
                       + ["" if r["total"] is None else f"{r['contrib'][t]:.3f}" for t in term_keys])

    lines = [f"# GFN-FF ({args.method}) per-term decomposition of GMTKN55 reaction energies", "",
             "AI-generated (scripts/revgfnff_barrier_terms.py), machine-evaluated. kcal/mol. Errors = method - published reference.",
             "gfn2 column: native gfn2 total error from the energy cache (diagnostic, no term table).",
             "Term columns: mean signed contribution of the term to the reaction energy over the subset;",
             "r/slope: Pearson correlation and least-squares slope of the total error against that term's",
             "contribution across the subset's reactions (a correlation, not a decomposition of the error).", ""]
    lines += ["## Per subset", "", "| subset | n | n ok | ref mean | MAD | MSE | max | gfn2 MAD | pbeh3c MAD | " +
              " | ".join(f"{t} mean" for t in term_keys) + " |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---:|" + "---:|" * len(term_keys)]
    explain = []
    for subset in args.subset:
        rs = by_subset.get(subset, [])
        ok = [r for r in rs if r["error"] is not None]
        if not rs:
            continue
        errs = [r["error"] for r in ok]
        mad = sum(abs(e) for e in errs) / len(errs) if errs else float("nan")
        mse = sum(errs) / len(errs) if errs else float("nan")
        mx = max(errs, key=abs) if errs else float("nan")
        g2 = [r["gfn2_error"] for r in rs if r["gfn2_error"] is not None]
        pb = [r["pbeh3c_error"] for r in rs if r["pbeh3c_error"] is not None]
        g2mad = sum(abs(e) for e in g2) / len(g2) if g2 else float("nan")
        pbmad = sum(abs(e) for e in pb) / len(pb) if pb else float("nan")
        refmean = sum(r["ref"] for r in rs) / len(rs)
        means = {t: (sum(r["contrib"][t] for r in ok) / len(ok) if ok else float("nan")) for t in term_keys}
        lines.append(f"| {subset} | {len(rs)} | {len(ok)} | {refmean:.2f} | {mad:.2f} | {mse:+.2f} | {mx:+.2f} | {g2mad:.2f} | {pbmad:.2f} | "
                     + " | ".join(f"{means[t]:+.2f}" for t in term_keys) + " |")
        corr = []
        for t in term_keys:
            r_, slope = pearson([r["contrib"][t] for r in ok], errs)
            if r_ is not None:
                corr.append((abs(r_), t, r_, slope))
        corr.sort(reverse=True)
        explain.append((subset, corr[:3]))
    lines += ["", "## Which term tracks the error (top 3 per subset)", "", "| subset | term | r | slope |", "|---|---|---:|---:|"]
    for subset, top in explain:
        for _, t, r_, slope in top:
            lines.append(f"| {subset} | {t} | {r_:+.2f} | {slope:+.2f} |")
    lines += ["", "## Reactions above 10 kcal/mol", "", "| subset | reaction | ref | error | gfn2 err | largest term contributions |", "|---|---|---:|---:|---:|---|"]
    for r in sorted((r for r in rows if r["error"] is not None and abs(r["error"]) > 10), key=lambda r: -abs(r["error"])):
        big = sorted(r["contrib"].items(), key=lambda kv: -abs(kv[1]))[:3]
        lines.append(f"| {r['subset']} | {r['reaction']} | {r['ref']:.1f} | {r['error']:+.1f} | "
                     f"{'-' if r['gfn2_error'] is None else f'{r['gfn2_error']:+.1f}'} | "
                     + ", ".join(f"{t} {v:+.1f}" for t, v in big) + " |")
    out_md = RESULTS / f"barrier_terms_{args.method}.md"
    out_md.write_text("\n".join(lines) + "\n")
    n_ok = sum(r["error"] is not None for r in rows)
    print(f"{n_ok}/{len(rows)} reactions scored; wrote {out_csv} and {out_md}")


if __name__ == "__main__":
    main()
