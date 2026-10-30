#!/usr/bin/env python3
"""PySCF references for the native DFT grid test (WP4, test_cases/qm_grid/reference.json).

For every case the molecular grid is built with PySCF's own routines and the curcuma settings:
    gen_grid.gen_atomic_grids(mol, atom_grid=(n_rad, n_ang), radi_method=..., prune=None)
    gen_grid.get_partition(mol, tab, radii_adjust, radi.BRAGG_RADII, becke_scheme)
(radi.treutler, gen_grid.original_becke | gen_grid.stratmann, radi.treutler_atomic_radii_adjust |
radi.becke_atomic_radii_adjust | None). PySCF stores the points of an atom in chunks of 12 radial shells
(angular outer, radial inner); they are re-indexed here to curcuma's order (radial outer, angular inner), so the
sampled points can be compared index by index.

The euler_maclaurin case has no PySCF radial scheme: the Murray-Handy-Laming formula is written out below and
passed as radi_method (checks the partition and assembly, not the radial formula itself).

Geometries: Angstrom from the xyz files, converted to Bohr with curcuma's CODATA 2018 factor 0.529177210903 and
passed to PySCF with unit='Bohr' (PySCF's own factor differs in the 11th digit).
Density cases: RHF/def2-SVP (spherical), conv_tol 1e-12, conv_tol_grad 1e-9; rho and grad rho from
numint.eval_ao(deriv=1) + eval_rho(xctype='GGA').

Run:  ~/.venvs/qmref/bin/python scripts/qm_grid_pyscf_reference.py
Claude Generated (Oct 2026).
"""
import json
import os
import sys

import numpy as np
import pyscf
from pyscf import dft, gto, scf
from pyscf.data import radii
from pyscf.dft import gen_grid, numint, radi

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.join(HERE, "..", "test_cases")
OUT = os.path.join(TC, "qm_grid", "reference.json")
BOHR = 0.529177210903  # curcuma units.h
STRIDE = 211

XYZ = {
    "H2O": "qm_1e/H2O.xyz",
    "CH4": "qm_1e/CH4.xyz",
    "LiH": "qm_1e/LiH.xyz",
    "Ne": "qm_1e/Ne.xyz",
    "benzene": "qm_hf3c/benzene.xyz",
}
LEVELS = {"sg0": (23, 170), "sg1": (50, 194), "sg2": (75, 302), "sg3": (99, 590)}

# case: (molecule, level, radial, partition, radii_adjust, density)
CASES = {
    "H2O_sg1": ("H2O", "sg1", "treutler", "becke", "treutler", True),
    "H2O_sg3": ("H2O", "sg3", "treutler", "becke", "treutler", False),
    "H2O_stratmann": ("H2O", "sg1", "treutler", "stratmann", "treutler", False),
    "H2O_adjust_becke": ("H2O", "sg1", "treutler", "becke", "becke", False),
    "H2O_adjust_none": ("H2O", "sg1", "treutler", "becke", "none", False),
    "H2O_euler": ("H2O", "sg1", "euler_maclaurin", "becke", "treutler", False),
    "CH4_sg1": ("CH4", "sg1", "treutler", "becke", "treutler", True),
    "benzene_sg1": ("benzene", "sg1", "treutler", "becke", "treutler", True),
    "LiH_sg2": ("LiH", "sg2", "treutler", "becke", "treutler", True),
    "Ne_sg0": ("Ne", "sg0", "treutler", "becke", "treutler", True),
}


def euler_maclaurin(n, chg, *args, **kwargs):
    """Murray-Handy-Laming m_R = 2 with R = Bragg radius (same formula as qm_grid.cpp)."""
    R = radii.BRAGG[chg]
    i = np.arange(1, n + 1, dtype=float)
    q = n + 1 - i
    r = R * (i / q) ** 2
    r2dr = 2.0 * R ** 3 * (n + 1) * i ** 5 / q ** 7
    return r, r2dr / r ** 2


def read_mol(name):
    path = os.path.join(TC, XYZ[name])
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        atoms.append([s, [float(x) / BOHR, float(y) / BOHR, float(z) / BOHR]])
    return gto.M(atom=atoms, unit="Bohr", basis="def2-svp", cart=False, verbose=0), XYZ[name]


def build_grid(mol, level, radial, partition, adjust):
    nr, na = LEVELS[level]
    rm = radi.treutler if radial == "treutler" else euler_maclaurin
    tab = gen_grid.gen_atomic_grids(mol, atom_grid=(nr, na), radi_method=rm, prune=None)
    adj = {"treutler": radi.treutler_atomic_radii_adjust, "becke": radi.becke_atomic_radii_adjust, "none": None}[adjust]
    scheme = gen_grid.original_becke if partition == "becke" else gen_grid.stratmann
    coords, weights = gen_grid.get_partition(mol, tab, adj, radi.BRAGG_RADII, scheme)
    # curcuma order: per atom, radial i outer, angular j inner
    per = nr * na
    perm = np.empty(mol.natm * per, dtype=int)
    for A in range(mol.natm):
        for i in range(nr):
            c = i // 12
            cs = min(12, nr - 12 * c)
            for j in range(na):
                perm[A * per + i * na + j] = A * per + 12 * c * na + j * cs + (i - 12 * c)
    return coords[perm], weights[perm], nr, na


def grid_ref(mol, coords, weights, nr, na):
    atm = mol.atom_coords()
    gauss = pz2 = 0.0
    for A in range(mol.natm):
        d = coords - atm[A]
        r2 = np.einsum("ij,ij->i", d, d)
        gauss += np.dot(weights, np.exp(-r2))
        pz2 += np.dot(weights, d[:, 2] ** 2 * np.exp(-2.0 * r2))
    per = nr * na
    atom_w = [float(weights[A * per:(A + 1) * per].sum()) for A in range(mol.natm)]
    idx = np.arange(0, len(weights), STRIDE)
    pts = [[int(p), float(coords[p, 0]), float(coords[p, 1]), float(coords[p, 2]), float(weights[p])] for p in idx]
    return {"npts": int(len(weights)), "n_radial": nr, "n_angular": na, "weight_sum": float(weights.sum()),
            "gauss_integral": float(gauss), "pz2_integral": float(pz2), "atom_weight_sums": atom_w,
            "stride": STRIDE, "points": pts}


def density_ref(mol, coords, weights):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.kernel()
    assert mf.converged
    dm = mf.make_rdm1()
    nel = 0.0
    for p0 in range(0, len(weights), 20000):
        ao = numint.eval_ao(mol, coords[p0:p0 + 20000], deriv=0)
        rho = numint.eval_rho(mol, ao, dm, xctype="LDA")
        nel += float(np.dot(weights[p0:p0 + 20000], rho))
    idx = np.arange(0, len(weights), STRIDE)
    ao = numint.eval_ao(mol, coords[idx], deriv=1)
    rho = numint.eval_rho(mol, ao, dm, xctype="GGA")
    pts = [[int(p), float(rho[0, k]), float(rho[1, k]), float(rho[2, k]), float(rho[3, k])] for k, p in enumerate(idx)]
    return {"energy": float(mf.e_tot), "num_electrons": int(mol.nelectron), "nelec_grid": nel,
            "deviation": nel - mol.nelectron, "points": pts}


def main():
    out = {"_comment": [
        "Native DFT grid (WP4) references from PySCF {} (venv ~/.venvs/qmref), scripts/qm_grid_pyscf_reference.py."
        .format(pyscf.__version__),
        "Grid: gen_atomic_grids(atom_grid=(n_rad, n_ang), prune=None) + get_partition(radii_adjust, BRAGG_RADII, "
        "becke_scheme); points re-indexed to curcuma's order (radial outer, angular inner); every 211th point stored.",
        "Density: RHF/def2-SVP spherical, conv_tol 1e-12, conv_tol_grad 1e-9; points: [index, rho, drho/dx, drho/dy, "
        "drho/dz] in a.u.; nelec_grid = sum w rho over the whole grid.",
        "Geometries in Bohr via 0.529177210903 (curcuma units.h). Claude Generated (Oct 2026)."],
        "cases": {}}
    for case, (mname, level, radial, part, adj, dens) in CASES.items():
        mol, xyz = read_mol(mname)
        coords, weights, nr, na = build_grid(mol, level, radial, part, adj)
        entry = {"xyz": xyz, "grid": level, "radial": radial, "partition": part, "radii_adjust": adj,
                 "grid_ref": grid_ref(mol, coords, weights, nr, na)}
        if dens:
            entry["density_ref"] = density_ref(mol, coords, weights)
        out["cases"][case] = entry
        print(case, entry["grid_ref"]["npts"], entry.get("density_ref", {}).get("deviation"), flush=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", os.path.normpath(OUT))


if __name__ == "__main__":
    sys.exit(main())
