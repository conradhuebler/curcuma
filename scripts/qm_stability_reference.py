#!/usr/bin/env python3
"""References for the native QM wavefunction stability analysis (-qm.scf_stability), test_cases/qm_stability/reference.json.

For every case PySCF converges the same SCF as curcuma (PySCF 2.14 / libxc 7.0.0, spherical basis by PySCF name,
geometry Angstrom -> Bohr with curcuma's CODATA 2018 factor, init_guess atom, conv_tol 1e-12, conv_tol_grad 1e-9;
KS on curcuma's grid as scripts/qm_uhf_reference.py) and forms the stability operator of curcuma's convention
DENSE, column by column, from PySCF's UHF response function (exact f_xc, no finite difference):

  (A + B) x = F_vv x - x F_oo + C_v^T G C_o,   D_s = C_v x_s C_o^T + transpose,
  G_a, G_b = UHF response (J[D_a + D_b] - a_x K[D_s] + f_xc^s [D_a, D_b]),  mf.gen_response(hermi=1) of the UHF object
  RHF -> RHF (singlet): D_a = D_b = D;  RHF -> UHF (triplet): D_a = D, D_b = -D;  UHF -> UHF: both spins.

Stored: energy, the lowest eigenvalue of each operator (numpy eigh of the symmetrised dense matrix), PySCF's own
stability() verdict, and for an unstable case the energy PySCF reaches by following (stability() rotation + SCF, at
most three times; RHF -> UHF continues as UHF from the rotated orbitals).

Usage: <venv with pyscf>/bin/python scripts/qm_stability_reference.py [--only CASE ...]

Claude Generated (Nov 2026). Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse
import json
import os

import numpy as np
import pyscf
from pyscf import dft, gto, scf
from pyscf.dft import gen_grid, radi

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.normpath(os.path.join(HERE, "..", "test_cases"))
OUT = os.path.join(TC, "qm_stability", "reference.json")
BOHR = 0.529177210903  # curcuma units.h
XC = {"b3lyp": "HYB_GGA_XC_B3LYP5", "pbe": "GGA_X_PBE,GGA_C_PBE"}

# id: (structure, method, basis, reference, charge, spin)
CASES = {
    "h2_r2.0_hf": ("small/dihydrogen.r2.0.xyz", "hf", "def2-SVP", "rhf", 0, 0),
    "h2_r2.0_uhf": ("small/dihydrogen.r2.0.xyz", "hf", "def2-SVP", "uhf", 0, 0),
    "h2o_hf": ("small/water.xyz", "hf", "def2-SVP", "rhf", 0, 0),
    "h2o_r1.8_hf": ("small/water.r1.8.xyz", "hf", "def2-SVP", "rhf", 0, 0),
    "o3_hf": ("small/ozone.xyz", "hf", "def2-SVP", "rhf", 0, 0),
    "n2_r2.0_hf": ("small/dinitrogen.r2.0.xyz", "hf", "def2-SVP", "rhf", 0, 0),
    "h2_r2.0_b3lyp": ("small/dihydrogen.r2.0.xyz", "b3lyp", "def2-SVP", "rhf", 0, 0),
    "o3_b3lyp": ("small/ozone.xyz", "b3lyp", "def2-SVP", "rhf", 0, 0),
    "h2_r2.0_ub3lyp": ("small/dihydrogen.r2.0.xyz", "b3lyp", "def2-SVP", "uhf", 0, 0),
    "ch3_uhf": ("small/methyl.xyz", "hf", "def2-SVP", "uhf", 0, 1),
}


def read_xyz(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0])
    atoms = []
    for ln in lines[2:2 + n]:
        p = ln.split()
        atoms.append((p[0], [float(v) / BOHR for v in p[1:4]]))
    return atoms


def grid(mf):
    mf.grids.atom_grid = (50, 194)  # curcuma sg1
    mf.grids.prune = None
    mf.grids.radi_method = radi.treutler
    mf.grids.becke_scheme = gen_grid.original_becke
    mf.grids.radii_adjust = radi.treutler_atomic_radii_adjust
    mf.small_rho_cutoff = 0.0


def make_mf(mol, method, ref):
    if method == "hf":
        mf = scf.RHF(mol) if ref == "rhf" else scf.UHF(mol)
    else:
        mf = dft.RKS(mol) if ref == "rhf" else dft.UKS(mol)
        mf.xc = XC[method]
        grid(mf)
    mf.init_guess = "atom"
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.max_cycle = 300
    mf.verbose = 0
    return mf


def dense_lowest(mfu, Ca, Cb, Fa, Fb, na, nb, mode):
    """Lowest eigenvalue of curcuma's (A + B), mode singlet | triplet | uhf, from the UHF response of mfu."""
    vind = mfu.gen_response(mo_coeff=(Ca, Cb), mo_occ=mfu.mo_occ, hermi=1)
    blocks = []
    for C, F, no in ((Ca, Fa, na), (Cb, Fb, nb)):
        Co, Cv = C[:, :no], C[:, no:]
        blocks.append((Co, Cv, Co.T @ F @ Co, Cv.T @ F @ Cv, Cv.shape[1], no))
    spins = [0, 1] if mode == "uhf" else [0]
    dims = [blocks[s][4] * blocks[s][5] for s in spins]
    dim = sum(dims)
    M = np.zeros((dim, dim))
    for k in range(dim):
        e = np.zeros(dim)
        e[k] = 1.0
        if mode == "uhf":
            xa = e[:dims[0]].reshape(blocks[0][4], blocks[0][5], order="F")
            xb = e[dims[0]:].reshape(blocks[1][4], blocks[1][5], order="F")
        else:
            xa = e.reshape(blocks[0][4], blocks[0][5], order="F")
            xb = xa if mode == "singlet" else -xa
        ds = []
        for s, x in ((0, xa), (1, xb)):
            Co, Cv = blocks[s][0], blocks[s][1]
            d = Cv @ x @ Co.T
            ds.append(d + d.T)
        v = vind(np.array(ds))
        out = []
        for s, x in ((0, xa), (1, xb)):
            Co, Cv, Foo, Fvv = blocks[s][:4]
            out.append(Fvv @ x - x @ Foo + Cv.T @ v[s] @ Co)
        if mode == "uhf":
            M[:, k] = np.concatenate([out[0].ravel(order="F"), out[1].ravel(order="F")])
        else:
            M[:, k] = out[0].ravel(order="F")
    M = 0.5 * (M + M.T)
    return float(np.linalg.eigvalsh(M)[0])


def run_case(cid):
    xyz, method, basis, ref, charge, spin = CASES[cid]
    mol = gto.M(atom=read_xyz(os.path.join(TC, "structures", xyz)), unit="Bohr", basis=basis, charge=charge, spin=spin,
                cart=False, verbose=0)
    mf = make_mf(mol, method, ref)
    e = mf.kernel()
    assert mf.converged, cid
    res = {"xyz": "structures/" + xyz, "method": method, "basis": basis, "reference": ref, "charge": charge, "spin": spin,
           "energy": float(e)}
    if ref == "rhf":
        mfu = scf.addons.convert_to_uhf(mf)
        C, F = mf.mo_coeff, mf.get_fock()
        no = mol.nelectron // 2
        res["rhf_rhf"] = dense_lowest(mfu, C, C, F, F, no, no, "singlet")
        res["rhf_uhf"] = dense_lowest(mfu, C, C, F, F, no, no, "triplet")
        st = mf.stability(external=True, return_status=True)
        res["pyscf_stable_internal"] = bool(st[2])
        res["pyscf_stable_external"] = bool(st[3])
        if not st[3]:  # follow RHF -> UHF as UHF
            mu = make_mf(mol, method, "uhf")
            mo = st[1]
            for _ in range(3):
                dma, dmb = mu.make_rdm1(mo, mfu.mo_occ)
                mu.kernel(dm0=(dma, dmb))
                r = mu.stability(return_status=True)
                mo, stable = r[0], r[2]
                if stable:
                    break
            res["followed_energy"] = float(mu.e_tot)
            res["followed_s2"] = float(mu.spin_square()[0])
    else:
        Fa, Fb = mf.get_fock()
        Ca, Cb = mf.mo_coeff
        na, nb = mol.nelec
        res["uhf_uhf"] = dense_lowest(mf, Ca, Cb, Fa, Fb, na, nb, "uhf")
        r = mf.stability(return_status=True)
        mo, stable = r[0], r[2]
        res["pyscf_stable_internal"] = bool(stable)
        if not stable:
            for _ in range(3):
                dma, dmb = mf.make_rdm1(mo, mf.mo_occ)
                mf.kernel(dm0=(dma, dmb))
                r = mf.stability(return_status=True)
                mo, stable = r[0], r[2]
                if stable:
                    break
            res["followed_energy"] = float(mf.e_tot)
            res["followed_s2"] = float(mf.spin_square()[0])
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    data = {"generator": "scripts/qm_stability_reference.py", "pyscf": pyscf.__version__, "cases": {}}
    if os.path.exists(OUT):
        with open(OUT) as f:
            data = json.load(f)
        data["pyscf"] = pyscf.__version__
    for cid in CASES:
        if a.only and cid not in a.only:
            continue
        data["cases"][cid] = run_case(cid)
        print(cid, json.dumps(data["cases"][cid]))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(data, f, indent=1, sort_keys=True)
        f.write("\n")


if __name__ == "__main__":
    main()
