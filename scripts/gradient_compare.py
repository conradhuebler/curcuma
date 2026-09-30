#!/usr/bin/env python3
# Claude Generated (Sep 2026)
"""Set-wide analytic-gradient comparison: curcuma native vs the xtb binary.

Until Sep 2026 curcuma's gradients were only spot-checked (three small molecules
in ctest plus ad-hoc finite differences), while the energies were validated over
thousands of structures. This closes that gap with the same method as
`gmtkn55_compare.py`: run both codes on identical geometry/charge/multiplicity and
diff the full 3N gradient vector.

UNITS -- the one thing to get right here:
  * curcuma's `-dump_gradient` writes EnergyCalculator::Gradient(), whose contract
    is **Eh/Angstrom** (native xTB converts explicitly in xtb_native.cpp; GFN-FF
    since Sep 2026 in gfnff_method.cpp).
  * xtb's Turbomole `gradient` file is **Eh/Bohr**.
  So curcuma is multiplied by `au` (Bohr->Angstrom) to reach Eh/Bohr before the diff.
  Everything below is reported in Eh/Bohr.

xtb APPENDS a cycle to an existing `gradient` file, so every run happens in a fresh
scratch directory.

Usage:
    python scripts/gradient_compare.py --set gmtkn55 --method gfn1
    python scripts/gradient_compare.py --set mor41 --method gfnff --limit 20
    python scripts/gradient_compare.py --set gmtkn55 --method all --recompute
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CURCUMA = ROOT / "release" / "curcuma"
AU = 0.52917721092  # Angstrom per Bohr (curcuma src/core/global.h)

SETS = {
    # name: (directory, per-structure xyz filename)
    "gmtkn55": (ROOT / "test_cases" / "GMTKN55-testset", "struc.xyz"),
    "mor41": (ROOT / "test_cases" / "MOR41-testset", "mol.xyz"),
    "s30lci": (ROOT / "test_cases" / "s30lci_test_set", "struc.xyz"),
}
XTB_FLAG = {"gfn1": ["--gfn", "1"], "gfn2": ["--gfn", "2"], "gfnff": ["--gfnff"]}


def discover(setdir, xyzname):
    """Yield (label, xyz_path, charge, uhf) for every structure in the set."""
    out = []
    for sub in sorted(p for p in setdir.iterdir() if p.is_dir() and not p.name.startswith("_")):
        # two layouts: <set>/<subset>/<name>/x.xyz and <set>/<name>/x.xyz
        direct = sub / xyzname
        if direct.exists():
            out.append((sub.name, direct))
            continue
        for st in sorted(p for p in sub.iterdir() if p.is_dir()):
            f = st / xyzname
            if f.exists():
                out.append((f"{sub.name}/{st.name}", f))
    res = []
    for label, f in out:
        d = f.parent
        chrg = int((d / ".CHRG").read_text().split()[0]) if (d / ".CHRG").exists() else 0
        uhf = int((d / ".UHF").read_text().split()[0]) if (d / ".UHF").exists() else 0
        res.append((label, f, chrg, uhf))
    return res


def curcuma_gradient(xyz, method, chrg, uhf, workdir):
    """Analytic gradient from curcuma, converted to Eh/Bohr."""
    dump = os.path.join(workdir, "cur.grad")
    if os.path.exists(dump):
        os.remove(dump)
    cmd = [str(CURCUMA), "-sp", str(xyz), "-method", method, "-gradient",
           "-dump_gradient", dump, "-no_bmt", "-verbosity", "0"]
    if chrg:
        cmd += ["-charge", str(chrg)]
    if uhf and method != "gfnff":
        cmd += ["-spin", str(uhf)]
    subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
    if not os.path.exists(dump):
        return None
    g = []
    for line in open(dump):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) == 3:
            try:
                g.append([float(x) * AU for x in p])   # Eh/Ang -> Eh/Bohr
            except ValueError:
                return None
    return g or None


def xtb_gradient(xyz, method, chrg, uhf, workdir):
    """Analytic gradient from the xtb binary (Turbomole `gradient` file, Eh/Bohr)."""
    for stale in ("gradient", "energy", "xtbrestart", "charges", "wbo", "xtbtopo.mol"):
        p = os.path.join(workdir, stale)
        if os.path.exists(p):
            os.remove(p)
    cmd = ["/opt/bin/xtb", str(xyz)] + XTB_FLAG[method] + ["--grad", "--acc", "0.0001"]
    if chrg:
        cmd += ["--chrg", str(chrg)]
    if uhf:
        cmd += ["--uhf", str(uhf)]
    subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
    gf = os.path.join(workdir, "gradient")
    if not os.path.exists(gf):
        return None
    rows = []
    for line in open(gf):
        p = line.split()
        if len(p) == 3:
            try:
                rows.append([float(x.replace("D", "E")) for x in p])
            except ValueError:
                pass          # the coordinate block carries a 4th (element) column
    return rows or None


def compare(method, structures, cache, recompute, timeout_note):
    stats, worst, skipped = [], [], 0
    for i, (label, xyz, chrg, uhf) in enumerate(structures, 1):
        key = f"{label}|{method}"
        if not recompute and key in cache:
            v = cache[key]
            if v is None:
                skipped += 1
                continue
            stats.append(v)
            worst.append((v, label))
            continue
        # A FRESH directory per structure. Every GMTKN55 structure file is named
        # `struc.xyz`, and GFN-FF caches its perceived topology next to the basename
        # (`struc.topo.json`), so a shared scratch directory silently applies one
        # molecule's topology to the next — that produced 1e+4 Eh/Bohr nonsense on the
        # first run of this script. xtb's `gradient`/`xtbrestart` have the same hazard.
        workdir = tempfile.mkdtemp(prefix="gradcmp_")
        c = curcuma_gradient(xyz, method, chrg, uhf, workdir)
        x = xtb_gradient(xyz, method, chrg, uhf, workdir)
        shutil.rmtree(workdir, ignore_errors=True)
        if c is None or x is None or len(c) != len(x):
            cache[key] = None
            skipped += 1
            continue
        dmax = max(abs(c[a][k] - x[a][k]) for a in range(len(c)) for k in range(3))
        cache[key] = dmax
        stats.append(dmax)
        worst.append((dmax, label))
        if i % 100 == 0:
            print(f"  [{i}/{len(structures)}] {label} dmax={dmax:.2e} (skipped {skipped})",
                  flush=True)
    return stats, worst, skipped


def main():
    ap = argparse.ArgumentParser(description="curcuma vs xtb analytic gradients")
    ap.add_argument("--set", default="gmtkn55", choices=sorted(SETS))
    ap.add_argument("--method", default="all", choices=["gfn1", "gfn2", "gfnff", "all"])
    ap.add_argument("--limit", type=int, default=0, help="only the first N structures")
    ap.add_argument("--recompute", action="store_true")
    args = ap.parse_args()

    setdir, xyzname = SETS[args.set]
    if not setdir.exists():
        sys.exit(f"test set not found: {setdir} (see docs/TESTSET_RETRIEVAL.md)")
    structures = discover(setdir, xyzname)
    if args.limit:
        structures = structures[: args.limit]
    print(f"{len(structures)} structures in {args.set}")

    rundir = setdir / "_run"
    rundir.mkdir(exist_ok=True)
    cachefile = rundir / "gradients.json"
    cache = json.loads(cachefile.read_text()) if cachefile.exists() else {}

    methods = ["gfn1", "gfn2", "gfnff"] if args.method == "all" else [args.method]
    for m in methods:
        print(f"\n=== method {m} ===", flush=True)
        stats, worst, skipped = compare(m, structures, cache, args.recompute, None)
        cachefile.write_text(json.dumps(cache))
        if not stats:
            print("  no comparable structures")
            continue
        worst.sort(reverse=True)
        mad = sum(stats) / len(stats)
        print(f"--- {m}: max-component |dG| curcuma-vs-xtb over n={len(stats)} "
              f"(skipped {skipped}): mean {mad:.3e}  worst {max(stats):.3e} Eh/Bohr")
        for thr in (1e-6, 1e-5, 1e-4, 1e-3):
            print(f"      above {thr:.0e}: {sum(1 for v in stats if v > thr)}")
        print("      worst 5: " + ", ".join(f"{n} {v:.2e}" for v, n in worst[:5]))


if __name__ == "__main__":
    main()
