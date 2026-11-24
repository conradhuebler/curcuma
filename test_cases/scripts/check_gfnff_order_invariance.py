#!/usr/bin/env python3
"""GFN-FF energy terms under row permutations of the input (atom renumbering).

Claude Generated (Nov 2026), docs/GFNFF_ORDER_DEPENDENCE.md.
For three small structures and 8 seeded permutations each (plus the identity), one -sp per order:
  1. default (inversion_sort_xtb true since Nov 24, 2026, dgam_pi_membership false): every term,
     inversion included, must agree over all orders to TOL, except the Coulomb spread of caffeine
     (dgam piadr index cutoff, pprcht and xtb). It is printed and checked to be PRESENT, so a silent
     change of the default shows up here.
  2. -gfnff.inversion_sort_xtb false (pprcht ssort, the default before Nov 24, 2026): the inversion
     spread of C6H5COOH and the acetic acid dimer must still be PRESENT (regression of the old variant).
  3. -gfnff.inversion_sort_xtb true -gfnff.dgam_pi_membership true: every term must agree to TOL.
Usage: check_gfnff_order_invariance.py <curcuma> <source dir>
"""
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

TOL = 1e-8          # Eh; summation-order rounding of permuted lists, measured <= 6e-9 (see doc)
EXPECT_MIN = 1e-7   # Eh; the reference-caused spreads are 1e-6 .. 1e-3 on these structures
NPERM = 8
TERMS = ("bond", "angle", "dihedral", "inversion", "stors", "batm", "disp", "brep", "nbrep",
         "coulomb", "hbond", "xbond")
# variant -> (extra flags, {structure: terms whose spread is expected})
VARIANTS = {
    "default": ([], {"caffeine": ("coulomb",)}),
    "pprcht": (["-gfnff.inversion_sort_xtb", "false"],
               {"C6H5COOH": ("inversion",), "acetic_acid_dimer": ("inversion",), "caffeine": ("coulomb",)}),
    "switches": (["-gfnff.inversion_sort_xtb", "true", "-gfnff.dgam_pi_membership", "true"], {}),
}
CASES = ("C6H5COOH", "acetic_acid_dimer", "caffeine")


def run(curcuma, rows, comment, wdir, extra):
    os.makedirs(wdir, exist_ok=True)
    with open(os.path.join(wdir, "s.xyz"), "w") as f:
        f.write(f"{len(rows)}\n{comment}\n" + "\n".join(rows) + "\n")
    cmd = [curcuma, "-sp", "s.xyz", "-method", "gfnff", "-gfnff.cache_topology", "false", "-threads", "1",
           "-verbosity", "3", "-no_bmt"] + extra
    out = subprocess.run(cmd, cwd=wdir, capture_output=True, text=True, errors="replace").stdout
    out = re.sub(r"\x1b\[[0-9;]*m", "", out)
    terms = {}
    for k in TERMS:
        m = re.findall(r"^\s+" + k + r"\s+=\s+([-+]\d\.\d+e[-+]\d+)", out, re.M)
        if not m:
            raise RuntimeError(f"term {k} missing in output ({wdir})")
        terms[k] = float(m[-1])
    return terms


def main():
    curcuma, src = sys.argv[1], sys.argv[2]
    failures = []
    with tempfile.TemporaryDirectory() as tmp:
        for name in CASES:
            lines = open(os.path.join(src, "test_cases/molecules/larger", name + ".xyz")).read().splitlines()
            nat = int(lines[0].split()[0])
            rows = lines[2:2 + nat]
            rng = np.random.default_rng(12345 + nat)
            perms = [list(range(nat))] + [list(rng.permutation(nat)) for _ in range(NPERM)]
            for label, (extra, expect_map) in VARIANTS.items():
                expected = expect_map.get(name, ())
                res = [run(curcuma, [rows[i] for i in p], lines[1], os.path.join(tmp, f"{name}_{label}_{o}"), extra)
                       for o, p in enumerate(perms)]
                for k in TERMS:
                    v = np.array([r[k] for r in res])
                    spread = v.max() - v.min()
                    if k in expected:
                        ok = spread > EXPECT_MIN
                        print(f"{name:18s} {label:8s} {k:10s} spread {spread:.3e} Eh  "
                              f"(reference behaviour, expected > {EXPECT_MIN:.0e}) {'ok' if ok else 'FAIL'}")
                    else:
                        ok = spread <= TOL
                        if not ok or k in ("inversion", "coulomb"):
                            print(f"{name:18s} {label:8s} {k:10s} spread {spread:.3e} Eh  (<= {TOL:.0e}) "
                                  f"{'ok' if ok else 'FAIL'}")
                    if not ok:
                        failures.append(f"{name}/{label}/{k}")
    print("FAILED: " + ", ".join(failures) if failures else "all checks passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
