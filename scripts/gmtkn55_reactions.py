#!/usr/bin/env python3
# Claude Generated (Sep 2026) - rev-gfnff WP0a
"""GMTKN55 reaction-level scoring from cached single points.

Reads the tmer2++ ``.res`` scripts of every subset (stoichiometry + published
reference energy), cross-checks them against the upstream
``_results/PBEh-3c_reactions.csv`` (which also carries the 30 BH76RC reactions
that have no ``.res`` of their own and the PBEh-3c value per reaction), and
scores any engine/method whose single points are cached in
``_run/energies.json`` (written by scripts/gmtkn55_compare.py).

Nothing is computed here: a structure without a cached energy makes the
reaction "unavailable" and is counted, never silently dropped.

Outputs (under test_cases/GMTKN55-testset/_results/):
  reactions_<engine>_<method>.csv   one row per reaction (all 1505)
  reactions_summary.md              MAD/RMSD/max per subset, per chemistry
                                    class and per official GMTKN55 category,
                                    plus WTMAD-2, for every method requested

Also importable: load_reactions(), load_energies(), structure_meta(),
evaluate(). scripts/revgfnff_barrier_terms.py builds on it.

Usage:
    python scripts/gmtkn55_reactions.py                        # cur gfnff/gfn1/gfn2 + PBEh-3c
    python scripts/gmtkn55_reactions.py --method gfnff --engine cur xtb
    python scripts/gmtkn55_reactions.py --subset AL2X6 --print  # show every reaction
"""
import argparse
import ast
import csv
import json
import math
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TESTSET = REPO / "test_cases" / "GMTKN55-testset"
RUNDIR = TESTSET / "_run"
RESULTS = TESTSET / "_results"
CSV_UPSTREAM = RESULTS / "PBEh-3c_reactions.csv"

AU2KCAL = 627.509474

# Chemistry classes used for the rev-gfnff work (docs/REV_GFNFF_TODO.md section 9
# had a similar split; that tooling was never committed, this map is explicit).
CHEM_CLASS = OrderedDict([
    ("conformers", ["ACONF", "Amino20x4", "BUT14DIOL", "ICONF", "MCONF", "PCONF21", "SCONF", "UPU23"]),
    ("intramolecular_nci", ["IDISP"]),
    ("nci", ["ADIM6", "CARBHB12", "HAL59", "HEAVY28", "PNICO23", "RG18", "S22", "S66", "WATER27"]),
    ("charged_nci", ["AHB21", "CHB6", "IL16"]),
    ("isomerisation", ["ISO34", "ISOL24", "C60ISO", "TAUT15", "PArel"]),
    ("barriers", ["BH76", "BHPERI", "BHDIV10", "BHROT27", "INV24", "PX13", "WCPT18"]),
    ("reactions_closed_shell", ["G2RC", "FH51", "DC13", "BSR36", "DARC", "CDIE20", "NBPRC", "AL2X6", "ALK8",
                                "HEAVYSB11", "YBDE18", "ALKBDE10", "MB16-43", "BH76RC", "DIPCS10", "PA26"]),
    ("reactions_open_shell", ["W4-11", "G21EA", "G21IP", "SIE4x4", "RC21", "RSE43"]),
])
SUBSET_CLASS = {s: c for c, subs in CHEM_CLASS.items() for s in subs}

# Official GMTKN55 categories (upstream _utils/constants.py); imported from there
# when available so the two never drift, hard copy as fallback.
def _load_upstream(name):
    """Import one upstream _utils module by file path (the package __init__ pulls in tqdm)."""
    import importlib.util
    path = TESTSET / "_utils" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"gmtkn55_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


try:
    _const = _load_upstream("constants")
    extract_species_from_path = _load_upstream("res_file").extract_species_from_path
    SMALL_REACTION_DIRS, LARGE_REACTION_DIRS, BARRIER_DIRS = (
        _const.SMALL_REACTION_DIRS, _const.LARGE_REACTION_DIRS, _const.BARRIER_DIRS)
    INTERMOL_NCI_DIRS, INTRAMOL_NCI_DIRS = _const.INTERMOL_NCI_DIRS, _const.INTRAMOL_NCI_DIRS
except Exception as exc:  # pragma: no cover
    sys.exit(f"cannot load the upstream GMTKN55 _utils modules from {TESTSET}: {exc}")

OFFICIAL_CATEGORY = OrderedDict([
    ("smallreactions", SMALL_REACTION_DIRS),
    ("largereactions", LARGE_REACTION_DIRS),
    ("barrierheights", BARRIER_DIRS),
    ("intermolecular", INTERMOL_NCI_DIRS),
    ("intramolecular", INTRAMOL_NCI_DIRS),
])
SUBSET_OFFICIAL = {s: c for c, subs in OFFICIAL_CATEGORY.items() for s in subs}


def structure_dir(subset):
    """Directory holding a subset's structures (BH76RC reuses BH76's)."""
    return "BH76" if subset == "BH76RC" else subset


class Reaction:
    __slots__ = ("subset", "species", "coeffs", "ref", "pbeh3c", "index")

    def __init__(self, subset, species, coeffs, ref, pbeh3c=None, index=0):
        self.subset = subset
        self.species = list(species)
        self.coeffs = list(coeffs)
        self.ref = float(ref)
        self.pbeh3c = pbeh3c
        self.index = index

    @property
    def key(self):
        return (self.subset, tuple(self.species), tuple(self.coeffs))

    @property
    def label(self):
        return " ".join(f"{c:+d}*{s}" for s, c in zip(self.species, self.coeffs))


# ------------------------------------------------------------------ .res parsing


def _species_tokens(tokens):
    out = []
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        path = tok.split("/$f")[0] if "/$f" in tok else tok.split("/")[0]
        if "{" in path:
            out.extend(extract_species_from_path(path))
        else:
            out.append(path)
    return out


def parse_res(subset):
    """All $tmer reactions of one subset's .res, in file order (duplicates kept)."""
    res = TESTSET / subset / ".res"
    reactions = []
    if not res.exists():
        return reactions
    for line in res.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#") or not s.startswith("$tmer"):
            continue
        toks = s.split()
        if "x" not in toks or "$w" not in toks:
            raise ValueError(f"{res}: malformed line: {line}")
        ix, iw = toks.index("x"), toks.index("$w")
        species = _species_tokens(toks[1:ix])
        coeffs = [int(t) for t in toks[ix + 1:iw]]
        ref = float(toks[iw + 1])
        if len(species) != len(coeffs):
            raise ValueError(f"{res}: {len(species)} species vs {len(coeffs)} coefficients: {line}")
        reactions.append(Reaction(subset, species, coeffs, ref, index=len(reactions)))
    return reactions


def load_csv():
    """Upstream reaction list: {(subset, species, coeffs): (ref, pbeh3c)}, in file order."""
    rows = []
    with CSV_UPSTREAM.open() as fh:
        for r in csv.DictReader(fh):
            species = tuple(ast.literal_eval(r["Reaction"]))
            coeffs = tuple(int(c) for c in ast.literal_eval(r["Stochiometry"]))
            rows.append((r["Subset"], species, coeffs, float(r["ReferenceValue"]), float(r["MethodValue"])))
    return rows


def load_reactions(subsets=None, strict=True):
    """Merge .res (primary) with the upstream CSV (cross-check + BH76RC + PBEh-3c)."""
    csv_rows = load_csv()
    csv_by_subset = defaultdict(list)
    for row in csv_rows:
        csv_by_subset[row[0]].append(row)
    wanted = subsets or sorted(csv_by_subset)
    reactions, problems = [], []
    for subset in wanted:
        res_rx = parse_res(subset)
        rows = list(csv_by_subset.get(subset, []))
        if not res_rx:
            if not rows:
                problems.append(f"{subset}: neither .res nor CSV rows")
                continue
            for i, (s, sp, co, ref, pb) in enumerate(rows):
                reactions.append(Reaction(s, sp, co, ref, pb, i))
            continue
        # match .res reactions to CSV rows by (species, coeffs), consuming duplicates in order
        pool = defaultdict(list)
        for row in rows:
            pool[(row[1], row[2])].append(row)
        for rx in res_rx:
            k = (tuple(rx.species), tuple(rx.coeffs))
            if pool.get(k):
                row = pool[k].pop(0)
                rx.pbeh3c = row[4]
                if abs(row[3] - rx.ref) > 0.006:  # CSV refs are rounded to 2 decimals; .res wins
                    problems.append(f"{subset}/{rx.label}: ref .res {rx.ref} vs CSV {row[3]}")
            else:
                problems.append(f"{subset}/{rx.label}: not in CSV")
            reactions.append(rx)
        leftover = sum(len(v) for v in pool.values())
        if leftover:
            problems.append(f"{subset}: {leftover} CSV rows without .res counterpart")
    if problems and strict:
        sys.exit("reaction list inconsistent:\n  " + "\n  ".join(problems))
    return reactions, problems


# ------------------------------------------------------------------ structures & energies

_META = {}


def structure_meta(subset, name):
    """(charge, uhf) of one structure from its .CHRG/.UHF files (cached)."""
    k = (subset, name)
    if k not in _META:
        d = TESTSET / structure_dir(subset) / name
        charge = int((d / ".CHRG").read_text().strip()) if (d / ".CHRG").exists() else 0
        uhf = 0
        if (d / ".UHF").exists():
            t = (d / ".UHF").read_text().strip()
            uhf = int(t) if t else 0
        _META[k] = (charge, uhf)
    return _META[k]


def load_energies(engine, method, cache=None):
    """{'subset/name': Eh or None} for one engine/method from _run/energies.json."""
    if cache is None:
        cache = json.loads((RUNDIR / "energies.json").read_text())
    suffix = f"|{engine}|{method}"
    out = {}
    for k, v in cache.items():
        if k.endswith(suffix):
            key = k[: -len(suffix)]
            out[key] = None if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)
    return out


def evaluate(rx, energies):
    """Reaction energy in kcal/mol, or None if any single point is missing."""
    total = 0.0
    for s, c in zip(rx.species, rx.coeffs):
        e = energies.get(f"{structure_dir(rx.subset)}/{s}")
        if e is None:
            return None
        total += c * e
    return total * AU2KCAL


def reaction_flags(rx):
    charged = any(structure_meta(rx.subset, s)[0] != 0 for s in rx.species)
    open_shell = any(structure_meta(rx.subset, s)[1] != 0 for s in rx.species)
    return charged, open_shell


# ------------------------------------------------------------------ statistics


def stat(errors):
    errs = [e for e in errors if e is not None]
    if not errs:
        return None
    n = len(errs)
    md = sum(errs) / n
    mad = sum(abs(e) for e in errs) / n
    rmsd = math.sqrt(sum(e * e for e in errs) / n)
    mx = max(errs, key=abs)
    return {"n": n, "MD": md, "MAD": mad, "RMSD": rmsd, "max": mx}


def wtmad2(per_subset, subsets_in_category=None):
    """WTMAD-2 as in the upstream statistics.py: sum_i N_i*(<|dE|>/|dE|_i)*MAD_i / sum_i N_i,
    with <|dE|> the mean over ALL subsets of the mean absolute reference energy."""
    all_avg = [v["mean_abs_ref"] for v in per_subset.values()]
    if not all_avg:
        return None
    mean_abs = sum(all_avg) / len(all_avg)
    subs = list(per_subset) if subsets_in_category is None else subsets_in_category
    if not subs:
        return None
    num = den = 0.0
    for s in subs:
        v = per_subset.get(s)
        if not v or v["stat"] is None:
            continue
        num += v["n_ref"] * (mean_abs / v["mean_abs_ref"]) * v["stat"]["MAD"]
        den += v["n_ref"]
    return num / den if den else None


def score(reactions, values):
    """values: {rx.index_key: kcal or None}. Returns per-subset/class/category tables."""
    by_subset = defaultdict(list)
    for rx in reactions:
        by_subset[rx.subset].append(rx)
    per_subset = OrderedDict()
    for s in sorted(by_subset):
        rxs = by_subset[s]
        errs = [None if values[id(rx)] is None else values[id(rx)] - rx.ref for rx in rxs]
        errs_cs = [e for rx, e in zip(rxs, errs) if not reaction_flags(rx)[1]]
        per_subset[s] = {
            "n_ref": len(rxs),
            "n_avail": sum(e is not None for e in errs),
            "mean_abs_ref": sum(abs(rx.ref) for rx in rxs) / len(rxs),
            "stat": stat(errs),
            "stat_closed": stat(errs_cs),
            "n_closed": len(errs_cs),
        }

    def group(mapping):
        out = OrderedDict()
        for cls, subs in mapping.items():
            errs, errs_cs = [], []
            for s in subs:
                for rx in by_subset.get(s, []):
                    e = None if values[id(rx)] is None else values[id(rx)] - rx.ref
                    errs.append(e)
                    if not reaction_flags(rx)[1]:
                        errs_cs.append(e)
            out[cls] = {"stat": stat(errs), "stat_closed": stat(errs_cs),
                        "n_ref": sum(1 for s in subs for _ in by_subset.get(s, [])),
                        "wtmad2": wtmad2(per_subset, [s for s in subs if s in per_subset])}
        return out

    return {"subset": per_subset, "chem": group(CHEM_CLASS), "official": group(OFFICIAL_CATEGORY),
            "wtmad2_total": wtmad2(per_subset),
            "total": stat([None if values[id(rx)] is None else values[id(rx)] - rx.ref for rx in reactions])}


# ------------------------------------------------------------------ output


def fmt(st, key="MAD"):
    return "-" if st is None else f"{st[key]:.3f}"


def write_csv(path, reactions, values, method_label):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["subset", "chem_class", "official", "reaction", "coeffs", "charged", "open_shell",
                    "ref_kcal", f"{method_label}_kcal", "error_kcal", "pbeh3c_kcal"])
        for rx in reactions:
            ch, os_ = reaction_flags(rx)
            v = values[id(rx)]
            w.writerow([rx.subset, SUBSET_CLASS.get(rx.subset, "?"), SUBSET_OFFICIAL.get(rx.subset, "?"),
                        " ".join(rx.species), " ".join(str(c) for c in rx.coeffs), int(ch), int(os_),
                        f"{rx.ref:.4f}", "" if v is None else f"{v:.4f}",
                        "" if v is None else f"{v - rx.ref:.4f}",
                        "" if rx.pbeh3c is None else f"{rx.pbeh3c:.4f}"])


def write_summary(path, reactions, tables, labels):
    lines = ["# GMTKN55 reaction-level scoring (kcal/mol)", "",
             "AI-generated (scripts/gmtkn55_reactions.py), machine-evaluated from cached single points.",
             f"Reactions: {len(reactions)} (from .res + upstream CSV). Errors are method - published reference.",
             "'closed' columns exclude every reaction that involves a structure with a nonzero .UHF.", ""]
    lines += ["## Totals", "", "| method | n | MAD | RMSD | max | MAD closed | n closed | WTMAD-2 |", "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for lab in labels:
        t = tables[lab]
        cs = stat([None if t["_values"][id(rx)] is None else t["_values"][id(rx)] - rx.ref
                   for rx in reactions if not reaction_flags(rx)[1]])
        tot = t["total"]
        lines.append(f"| {lab} | {tot['n'] if tot else 0} | {fmt(tot)} | {fmt(tot, 'RMSD')} | {fmt(tot, 'max')} | "
                     f"{fmt(cs)} | {cs['n'] if cs else 0} | {t['wtmad2_total']:.3f} |" if t['wtmad2_total'] is not None else
                     f"| {lab} | {tot['n'] if tot else 0} | {fmt(tot)} | {fmt(tot, 'RMSD')} | {fmt(tot, 'max')} | {fmt(cs)} | {cs['n'] if cs else 0} | - |")
    for title, key, mapping in (("Chemistry classes (rev-gfnff)", "chem", CHEM_CLASS),
                                ("Official GMTKN55 categories", "official", OFFICIAL_CATEGORY)):
        lines += ["", f"## {title}", "", "| class | n | " + " | ".join(f"{l} MAD (closed) | {l} max" for l in labels) + " |",
                  "|---|---:|" + "---:|---:|" * len(labels)]
        for cls in mapping:
            row = [cls, str(tables[labels[0]][key][cls]["n_ref"])]
            for lab in labels:
                g = tables[lab][key][cls]
                row += [f"{fmt(g['stat'])} ({fmt(g['stat_closed'])})", fmt(g["stat"], "max")]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("| WTMAD-2 | | " + " | ".join(
            (f"{tables[l][key][cls]['wtmad2']:.3f}" if False else "") for l in labels for cls in []) + "")
        lines.pop()  # (no per-class WTMAD-2 row in this table; see per-category line below)
        if key == "official":
            for lab in labels:
                lines.append("")
                lines.append(f"WTMAD-2 {lab}: total {tables[lab]['wtmad2_total']:.3f}; " + "; ".join(
                    f"{cls} {tables[lab]['official'][cls]['wtmad2']:.3f}" for cls in mapping
                    if tables[lab]['official'][cls]['wtmad2'] is not None))
    lines += ["", "## Per subset", "", "| subset | class | n | avail | " + " | ".join(f"{l} MAD | {l} max" for l in labels) + " |",
              "|---|---|---:|---:|" + "---:|---:|" * len(labels)]
    for s in tables[labels[0]]["subset"]:
        row = [s, SUBSET_CLASS.get(s, "?"), str(tables[labels[0]]["subset"][s]["n_ref"]),
               str(tables[labels[0]]["subset"][s]["n_avail"])]
        for lab in labels:
            st = tables[lab]["subset"][s]["stat"]
            row += [fmt(st), fmt(st, "max")]
        lines.append("| " + " | ".join(row) + " |")
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", nargs="+", default=["gfnff", "gfn1", "gfn2"])
    ap.add_argument("--engine", nargs="+", default=["cur"], help="cur and/or xtb (cache engines)")
    ap.add_argument("--no-pbeh3c", action="store_true", help="omit the upstream PBEh-3c column")
    ap.add_argument("--subset", nargs="*", help="restrict to these subsets")
    ap.add_argument("--print", action="store_true", help="print every reaction with all values")
    ap.add_argument("--lenient", action="store_true", help="report .res/CSV mismatches instead of aborting")
    args = ap.parse_args()

    reactions, problems = load_reactions(args.subset, strict=not args.lenient)
    for p in problems:
        print("WARN", p)
    n_cs = sum(1 for rx in reactions if not reaction_flags(rx)[1])
    print(f"{len(reactions)} reactions from {len({rx.subset for rx in reactions})} subsets "
          f"({n_cs} closed-shell, {len(reactions) - n_cs} involve an open-shell structure)")

    cache = json.loads((RUNDIR / "energies.json").read_text())
    tables, labels = OrderedDict(), []
    RESULTS.mkdir(exist_ok=True)
    for engine in args.engine:
        for method in args.method:
            energies = load_energies(engine, method, cache)
            if not energies:
                print(f"no cached energies for {engine}/{method}, skipped")
                continue
            values = {id(rx): evaluate(rx, energies) for rx in reactions}
            lab = f"{engine}_{method}"
            t = score(reactions, values)
            t["_values"] = values
            tables[lab] = t
            labels.append(lab)
            write_csv(RESULTS / f"reactions_{lab}.csv", reactions, values, lab)
            miss = sum(v is None for v in values.values())
            print(f"{lab}: MAD {fmt(t['total'])} RMSD {fmt(t['total'], 'RMSD')} max {fmt(t['total'], 'max')} "
                  f"WTMAD-2 {t['wtmad2_total']:.3f}  ({miss} reactions without energy)")
    if not args.no_pbeh3c:
        values = {id(rx): rx.pbeh3c for rx in reactions}
        t = score(reactions, values)
        t["_values"] = values
        tables["pbeh3c"] = t
        labels.append("pbeh3c")
        print(f"pbeh3c: MAD {fmt(t['total'])} WTMAD-2 {t['wtmad2_total']:.3f} (upstream CSV; upstream file says total 11.129)")
    if args.print:
        for rx in reactions:
            vals = " ".join(f"{lab}={tables[lab]['_values'][id(rx)]:.2f}" if tables[lab]['_values'][id(rx)] is not None
                            else f"{lab}=NA" for lab in labels)
            print(f"{rx.subset:10s} {rx.label:60s} ref={rx.ref:9.2f} {vals}")
    if labels:
        out = RESULTS / "reactions_summary.md"
        write_summary(out, reactions, tables, labels)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
