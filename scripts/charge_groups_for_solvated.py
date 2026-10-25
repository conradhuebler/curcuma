#!/usr/bin/env python3
# Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated (Oct 2026): EEQ charge groups for a solvated polymer, for -gfnff.charge_groups FILE with "rest": "per_fragment".
"""Write a charge-group file for a polymer (polymerbuild unit groups) in water with isolated counter-ions.

Groups written:
  - every unit group of polymer_A_final.groups.json, atom indices shifted by --polymer-offset;
  - the polymer atoms outside all unit groups (neutral units, end groups) as ONE more group, charge --polymer-rest-charge;
    the polymer fragment must be covered completely, because with rest per_fragment a partly covered fragment is rejected;
  - every isolated anion (atom of --anion without a bonded neighbour by covalent radii) as a group, charge --anion-charge.
All other fragments (water) are left out: GFN-FF makes each of them its own group with charge 0. The sum of the group charges
is checked against --charge (the total charge passed to curcuma); a mismatch is an error and nothing is written.

Usage:
  charge_groups_for_solvated.py SYSTEM.xyz --polymer-groups polymer_A_final.groups.json --polymer-offset 0 \
      --anion Cl --charge 0 --out groups.json
Python standard library and numpy only. Atom indices are 0-based.
"""
import argparse
import json
import sys

import numpy as np

# covalent radii in Angstrom (Pyykko & Atsumi 2009, single bond); bond if d < 1.2 (r_i + r_j)
RCOV = {"H": 0.32, "B": 0.85, "C": 0.75, "N": 0.71, "O": 0.63, "F": 0.64, "Na": 1.55, "Mg": 1.39, "Si": 1.16, "P": 1.11,
        "S": 1.03, "Cl": 0.99, "K": 1.96, "Ca": 1.71, "Br": 1.14, "I": 1.33}
BOND_SCALE = 1.2


def read_xyz(path):
    with open(path) as f:
        lines = f.read().split("\n")
    n = int(lines[0])
    el, xyz = [], []
    for line in lines[2:2 + n]:
        p = line.split()
        el.append(p[0][0].upper() + p[0][1:].lower())
        xyz.append([float(v) for v in p[1:4]])
    return el, np.array(xyz)


def fragments(el, xyz):
    """connected components by covalent radii; returns (fragment id per atom (ordered by smallest atom), neighbour count)"""
    n = len(el)
    r = np.array([RCOV.get(e, 1.5) for e in el])
    parent = list(range(n))

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    nb = np.zeros(n, dtype=int)
    for i in range(n):
        d = np.linalg.norm(xyz[i + 1:] - xyz[i], axis=1)
        for k in np.nonzero(d < BOND_SCALE * (r[i + 1:] + r[i]))[0]:
            j = i + 1 + int(k)
            nb[i] += 1
            nb[j] += 1
            parent[root(i)] = root(j)
    label, frag = {}, []
    for i in range(n):
        frag.append(label.setdefault(root(i), len(label)))
    return frag, nb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("system", help="solvated system (xyz)")
    ap.add_argument("--polymer-groups", required=True, help="unit groups from polymerbuild (polymer_A_final.groups.json)")
    ap.add_argument("--polymer-offset", type=int, default=0, help="index of the first polymer atom in SYSTEM (0-based)")
    ap.add_argument("--polymer-rest-charge", type=float, default=0.0, help="charge of the polymer atoms outside the unit groups")
    ap.add_argument("--anion", action="append", default=[], help="element of isolated anions (repeatable), e.g. Cl")
    ap.add_argument("--anion-charge", type=float, default=-1.0, help="charge of each isolated anion")
    ap.add_argument("--charge", type=float, default=0.0, help="total charge of SYSTEM (checked against the group sum)")
    ap.add_argument("--out", required=True, help="output JSON")
    a = ap.parse_args()

    el, xyz = read_xyz(a.system)
    n = len(el)
    frag, nb = fragments(el, xyz)
    with open(a.polymer_groups) as f:
        units = json.load(f)["groups"]
    off = a.polymer_offset
    if not 0 <= off < n:
        sys.exit(f"error: polymer offset {off} outside 0..{n - 1}")
    fpoly = frag[off]
    poly_atoms = [i for i in range(n) if frag[i] == fpoly]
    poly_set = set(poly_atoms)
    groups, owner = [], {}
    for g in units:
        atoms = sorted(int(x) + off for x in g["atoms"])
        for i in atoms:
            if i not in poly_set:
                sys.exit(f"error: unit atom {i} (offset {off}) is not in the polymer fragment of atom {off}")
            if i in owner:
                sys.exit(f"error: atom {i} is in two unit groups")
            owner[i] = len(groups)
        groups.append({"atoms": atoms, "charge": float(g["charge"])})
    n_units = len(groups)
    rest = [i for i in poly_atoms if i not in owner]
    if rest:
        groups.append({"atoms": rest, "charge": a.polymer_rest_charge})
    anions = [i for i in range(n) if el[i] in a.anion and nb[i] == 0]
    for i in anions:
        groups.append({"atoms": [i], "charge": a.anion_charge})
    total = sum(g["charge"] for g in groups)
    n_frag = len(set(frag))
    print(f"system {n} atoms, {n_frag} fragments (covalent radii x {BOND_SCALE}); polymer fragment {len(poly_atoms)} atoms from atom {off}: "
          f"{n_units} unit groups, rest group {len(rest)} atoms (charge {a.polymer_rest_charge:+.3f}); {len(anions)} isolated anions "
          f"({','.join(a.anion) or '-'}, {a.anion_charge:+.3f} each); {n_frag - 1 - len(anions)} uncovered fragments (charge 0 in GFN-FF)")
    if abs(total - a.charge) > 1e-6:
        sys.exit(f"error: group charges sum to {total:.6f}, total charge is {a.charge}; nothing written")
    with open(a.out, "w") as f:
        json.dump({"rest": "per_fragment", "groups": groups}, f, indent=1)
    print(f"written {a.out}: {len(groups)} groups, sum {total:+.6f}")
    # Claude Generated (Oct 2026): recommended flags (docs/POLYMERBUILD_CHARGE.md, "Recommended run recipe").
    # cache_topology false and topology_mode constant are what curcuma uses with a group file anyway (a group file
    # bypasses the topology cache; constant is the default); they are written out so the run does not depend on that.
    print(f"recommended: -method gfnff -charge {a.charge:g} -gfnff.charge_groups {a.out} -gfnff.cn_zero_guard true "
          f"-gfnff.cache_topology false -gfnff.topology_mode constant")


if __name__ == "__main__":
    main()
