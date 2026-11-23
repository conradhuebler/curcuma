#!/usr/bin/env python3
"""Atom-order (row permutation) scan for GFN-FF single points.

Claude Generated (Nov 2026). Writes N random row permutations of an XYZ file (coordinate
lines copied textually, only their order changes; seed derived from the file name so each
structure has its own reproducible set), runs one program per order in its own directory,
parses the energy terms and the gradient, maps the gradient back to the original atom order
and prints the spread per term (std and max-min over all orders, identity included).

Programs:
  curcuma  -sp <xyz> -method gfnff -gfnff.cache_topology false -verbosity 3 -dump_gradient
  pprcht   pprcht/gfnff app (external/gfnff, cmake -Dbuild_exe=ON); terms from its result print
           (dihedral/inversion split only from an instrumented build, docs/GFNFF_ORDER_DEPENDENCE.md)
  xtb      xtb --gfnff (terms from the "SUMMARY" block)

Usage:
  gfnff_order_scan.py --xyz FILE --n 30 --program curcuma --bin release/curcuma --out DIR
                      [--threads 4] [--extra "-gfnff.some_switch true"] [--keep]
  gfnff_order_scan.py --compare curcuma_dir/<structure>/scan.json pprcht_dir/<structure>/scan.json
Output: DIR/<structure>/scan.json (per order: terms, gradient max diff vs identity) and a
one-line-per-term summary on stdout. Only standard library + numpy.
"""
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

import numpy as np


def read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    nat = int(lines[0].split()[0])
    return lines[1], lines[2:2 + nat]


def write_xyz(path, comment, rows):
    with open(path, "w") as f:
        f.write(f"{len(rows)}\n{comment}\n" + "\n".join(rows) + "\n")


def permutations(nat, n, name):
    seed = int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    perms = [list(range(nat))]
    for _ in range(n):
        perms.append([int(x) for x in rng.permutation(nat)])
    return perms  # perms[o][r] = original index of the atom on row r


# --- parsers --------------------------------------------------------------------------
CURCUMA_TERMS = {  # "=== CPU ENERGY TERMS ===" block of -verbosity 3 (17 digits)
    k: r"^\s+" + k + r"\s+=\s+([-+]\d\.\d+e[-+]\d+)"
    for k in ("bond", "angle", "dihedral", "inversion", "stors", "batm", "atm", "disp", "brep", "nbrep",
              "coulomb", "hbond", "xbond")
}


def parse_terms(text, patterns):
    out = {}
    for key, pat in patterns.items():
        m = re.findall(pat, text, re.M)
        if m:
            out[key] = float(m[-1])
    return out


def run_curcuma(binary, xyz, wdir, threads, extra):
    cmd = [binary, "-sp", os.path.basename(xyz), "-method", "gfnff", "-gfnff.cache_topology", "false",
           "-threads", str(threads), "-verbosity", "3", "-dump_gradient", "grad.txt", "-no_bmt"] + extra
    with open(os.path.join(wdir, "run.log"), "w") as log:
        rc = subprocess.call(cmd, cwd=wdir, stdout=log, stderr=subprocess.STDOUT)
    text = re.sub(r"\x1b\[[0-9;]*m", "", open(os.path.join(wdir, "run.log"), errors="replace").read())
    terms = parse_terms(text, CURCUMA_TERMS)
    gpath = os.path.join(wdir, "grad.txt")
    grad = None
    if os.path.exists(gpath):
        lines = open(gpath).read().splitlines()
        m = re.search(r"energy\s+(-?\d+\.\d+)", lines[1])
        if m:
            terms["total"] = float(m.group(1))
        grad = np.array([[float(x) for x in l.split()] for l in lines[2:] if l.strip()])
    return rc, terms, grad


PPRCHT_TERMS = {
    "total": r"total energy\s+(-?\d+\.\d+)",
    "bond": r"bond energy\s+(-?\d+\.\d+)",
    "angle": r"angle energy\s+(-?\d+\.\d+)",
    "torsion_block": r"torsion energy\s+(-?\d+\.\d+)",
    "repulsion": r"repulsion energy\s+(-?\d+\.\d+)",
    "coulomb": r"electrostat energy\s+(-?\d+\.\d+)",
    "dispersion": r"dispersion energy\s+(-?\d+\.\d+)",
    "hbond": r"HB energy\s+(-?\d+\.\d+)",
    "xbond": r"XB energy\s+(-?\d+\.\d+)",
    "batm": r"bonded atm energy\s+(-?\d+\.\d+)",
}


def run_pprcht(binary, xyz, wdir, threads, extra):
    env = dict(os.environ, OMP_NUM_THREADS="1", PPRCHT_INVDBG="1")
    with open(os.path.join(wdir, "run.log"), "w") as log:
        rc = subprocess.call([binary, os.path.basename(xyz)] + extra, cwd=wdir, stdout=log,
                             stderr=subprocess.STDOUT, env=env)
    text = open(os.path.join(wdir, "run.log"), errors="replace").read()
    terms = parse_terms(text, PPRCHT_TERMS)
    m = re.findall(r"PPTORS_DIH_INV\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", text)
    if m:  # instrumented build only: dihedral and inversion separately
        terms["torsion"] = float(m[-1][0])
        terms["inversion"] = float(m[-1][1])
    g = re.findall(r"Atom\s+\d+ :\s+(\S+)\s+(\S+)\s+(\S+)", text)
    grad = np.array([[float(x) for x in r] for r in g]) / 0.529177210544 if g else None  # Eh/a0 -> Eh/A
    return rc, terms, grad


XTB_TERMS = {
    "total": r":: total energy\s+(-?\d+\.\d+)",
    "bond": r":: bond energy\s+(-?\d+\.\d+)",
    "angle": r":: angle energy\s+(-?\d+\.\d+)",
    "torsion_block": r":: torsion energy\s+(-?\d+\.\d+)",
    "repulsion": r":: repulsion energy\s+(-?\d+\.\d+)",
    "coulomb": r":: electrostat energy\s+(-?\d+\.\d+)",
    "dispersion": r":: dispersion energy\s+(-?\d+\.\d+)",
    "hbond": r":: HB energy\s+(-?\d+\.\d+)",
    "xbond": r":: XB energy\s+(-?\d+\.\d+)",
    "batm": r":: bonded atm energy\s+(-?\d+\.\d+)",
}


def run_xtb(binary, xyz, wdir, threads, extra):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    with open(os.path.join(wdir, "run.log"), "w") as log:
        rc = subprocess.call([binary, os.path.basename(xyz), "--gfnff", "--norestart"] + extra, cwd=wdir,
                             stdout=log, stderr=subprocess.STDOUT, env=env)
    text = open(os.path.join(wdir, "run.log"), errors="replace").read()
    return rc, parse_terms(text, XTB_TERMS), None


RUNNERS = {"curcuma": run_curcuma, "pprcht": run_pprcht, "xtb": run_xtb}

# curcuma term -> pprcht/xtb term (repulsion is brep + nbrep in curcuma)
COMMON = {"total": "total", "bond": "bond", "angle": "angle", "dihedral": "torsion", "inversion": "inversion",
          "coulomb": "coulomb", "disp": "dispersion", "hbond": "hbond", "xbond": "xbond", "batm": "batm",
          "repulsion": "repulsion"}


def compare(path_a, path_b):
    """Per order and term: |curcuma - reference| for two scan.json files of the same structure."""
    ra, rb = json.load(open(path_a)), json.load(open(path_b))
    assert [r["perm"] for r in ra] == [r["perm"] for r in rb], "different permutation sets"
    worst = {}
    for x, y in zip(ra, rb):
        tx, ty = dict(x["terms"]), y["terms"]
        if "brep" in tx:
            tx["repulsion"] = tx["brep"] + tx.get("nbrep", 0.0)
        for kc, kr in COMMON.items():
            if kc in tx and kr in ty:
                worst[kc] = max(worst.get(kc, 0.0), abs(tx[kc] - ty[kr]))
    print(f"{path_a} vs {path_b}: max over {len(ra)} orders of |a - b| (Eh)")
    for k, v in worst.items():
        print(f"  {k:10s} {v:.3e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compare", nargs=2, metavar=("A", "B"), help="compare two scan.json files and exit")
    ap.add_argument("--xyz")
    ap.add_argument("--n", type=int, default=30)
    ap.add_argument("--program", choices=RUNNERS, default="curcuma")
    ap.add_argument("--bin")
    ap.add_argument("--out")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--extra", default="", help="extra arguments, split on whitespace")
    ap.add_argument("--keep", action="store_true", help="keep the per-order run directories")
    a = ap.parse_args()
    if a.compare:
        return compare(*a.compare)
    if not (a.xyz and a.bin and a.out):
        ap.error("--xyz, --bin and --out are required for a scan")

    name = os.path.splitext(os.path.basename(a.xyz))[0]
    comment, rows = read_xyz(a.xyz)
    perms = permutations(len(rows), a.n, name)
    base = os.path.join(a.out, name)
    os.makedirs(base, exist_ok=True)
    binary = os.path.abspath(a.bin) if os.path.exists(a.bin) else a.bin
    results = []
    gref = None
    for o, p in enumerate(perms):
        wdir = os.path.join(base, f"o{o:02d}")
        os.makedirs(wdir, exist_ok=True)
        xyz = os.path.join(wdir, f"{name}.xyz")
        write_xyz(xyz, comment, [rows[i] for i in p])
        rc, terms, grad = RUNNERS[a.program](binary, xyz, wdir, a.threads, a.extra.split())
        gdiff = None
        if grad is not None and len(grad) == len(rows):
            gorig = np.zeros_like(grad)
            gorig[p] = grad  # row r holds original atom p[r]
            if gref is None:
                gref = gorig
            gdiff = float(np.abs(gorig - gref).max())
        results.append({"order": o, "perm": p, "rc": rc, "terms": terms, "grad_maxdiff": gdiff})
        if not a.keep and o > 0:
            shutil.rmtree(wdir)
    json.dump(results, open(os.path.join(base, "scan.json"), "w"), indent=1)

    keys = sorted({k for r in results for k in r["terms"]})
    bad = [r["order"] for r in results if r["rc"] != 0]
    print(f"{name} {a.program} orders={len(results)} failed={bad}")
    for k in keys:
        v = np.array([r["terms"][k] for r in results if k in r["terms"]])
        if len(v) < 2:
            continue
        print(f"  {k:14s} n={len(v):3d} std={np.std(v) * 1e9:14.1f} nEh  span={(v.max() - v.min()) * 1e9:14.1f} nEh"
              f"  distinct={len(np.unique(np.round(v, 10)))}")
    gd = [r["grad_maxdiff"] for r in results if r["grad_maxdiff"] is not None]
    if gd:
        print(f"  gradient       max |g(order)-g(identity)| = {max(gd):.3e} Eh/A over {len(gd)} orders")


if __name__ == "__main__":
    sys.exit(main())
