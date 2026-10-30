#!/usr/bin/env python3
"""PySCF/libxc references for the native LDA test (WP5, test_cases/qm_lda/reference.json).

kernel:   libxc (through pyscf.dft.libxc.eval_xc, spin 0) on 61 densities rho = 1e-12 ... 1e2 (log-spaced):
          eps and v = d(rho eps)/d rho for LDA_X, LDA_C_VWN (VWN5) and LDA_C_VWN_RPA.
energies: RKS/def2-SVP (spherical) with the SAME grid as curcuma (qm_grid.cpp, WP4):
              grids.atom_grid = (n_rad, n_ang), prune = None, radi_method = radi.treutler,
              becke_scheme = gen_grid.original_becke, radii_adjust = radi.treutler_atomic_radii_adjust
          xc = 'LDA_X,LDA_C_VWN' (vwn5) or 'LDA_X,LDA_C_VWN_RPA' (vwn_rpa); small_rho_cutoff = 0 (no grid pruning by
          density), conv_tol 1e-12, conv_tol_grad 1e-9. Stored: E_total, E_xc (veff.exc), N on the grid.
          RI cases: mf.density_fit(auxbasis=<curcuma's def2-universal-JFIT/-JKFIT.turbomole>, parsed by
          rimp2_reference.parse_turbomole); for a pure functional PySCF fits J only, as curcuma's -qm.ri j|jk for lda.
Geometries: Angstrom xyz files converted to Bohr with curcuma's CODATA 2018 factor (as qm_grid_pyscf_reference.py).

Run:  ~/.venvs/qmref/bin/python scripts/qm_lda_pyscf_reference.py
Claude Generated (Oct 2026, WP5).
"""
import json
import os

import numpy as np
import importlib.util

import pyscf
from pyscf import dft, gto
from pyscf.dft import gen_grid, libxc, radi

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.join(HERE, "..", "test_cases")
OUT = os.path.join(TC, "qm_lda", "reference.json")
BOHR = 0.529177210903  # curcuma units.h

XYZ = {
    "H2O": "qm_1e/H2O.xyz",
    "CH4": "qm_1e/CH4.xyz",
    "NH3": "qm_1e/NH3.xyz",
    "HF": "qm_1e/HF.xyz",
    "LiH": "qm_1e/LiH.xyz",
    "Ne": "qm_1e/Ne.xyz",
    "benzene": "qm_hf3c/benzene.xyz",
}
LEVELS = {"sg0": (23, 170), "sg1": (50, 194), "sg2": (75, 302), "sg3": (99, 590)}
XC = {"vwn5": "LDA_X,LDA_C_VWN", "vwn_rpa": "LDA_X,LDA_C_VWN_RPA"}

QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")
AUX = {"j": "def2-universal-JFIT.turbomole", "jk": "def2-universal-JKFIT.turbomole"}

# case: (molecule, grid level, vwn variant, ri mode)
CASES = {
    "H2O_sg1": ("H2O", "sg1", "vwn5", "off"),
    "CH4_sg1": ("CH4", "sg1", "vwn5", "off"),
    "NH3_sg1": ("NH3", "sg1", "vwn5", "off"),
    "HF_sg1": ("HF", "sg1", "vwn5", "off"),
    "LiH_sg1": ("LiH", "sg1", "vwn5", "off"),
    "Ne_sg1": ("Ne", "sg1", "vwn5", "off"),
    "benzene_sg1": ("benzene", "sg1", "vwn5", "off"),
    "H2O_sg3": ("H2O", "sg3", "vwn5", "off"),
    "H2O_rpa_sg1": ("H2O", "sg1", "vwn_rpa", "off"),
    "NH3_rpa_sg2": ("NH3", "sg2", "vwn_rpa", "off"),
    "H2O_rij_sg1": ("H2O", "sg1", "vwn5", "j"),
    "benzene_rijk_sg1": ("benzene", "sg1", "vwn5", "jk"),
}


def parse_turbomole(path):
    spec = importlib.util.spec_from_file_location("rimp2_reference", os.path.join(HERE, "rimp2_reference.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.parse_turbomole(path)


def read_mol(name):
    path = os.path.join(TC, XYZ[name])
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        atoms.append([s, [float(x) / BOHR, float(y) / BOHR, float(z) / BOHR]])
    return gto.M(atom=atoms, unit="Bohr", basis="def2-svp", cart=False, verbose=0)


def kernel_ref():
    rho = np.logspace(-12, 2, 61)
    out = {"rho": rho.tolist()}
    for key, name in (("slater", "LDA_X"), ("vwn5", "LDA_C_VWN"), ("vwn_rpa", "LDA_C_VWN_RPA")):
        exc, vxc = libxc.eval_xc(name, rho, spin=0, deriv=1)[:2]
        out[key] = {"libxc": name, "eps": exc.tolist(), "v": vxc[0].tolist()}
    return out


def energy_ref(mol_name, level, variant, ri):
    mol = read_mol(mol_name)
    mf = dft.RKS(mol)
    if ri != "off":
        mf = mf.density_fit(auxbasis=parse_turbomole(os.path.join(QMDIR, AUX[ri])))
    mf.xc = XC[variant]
    mf.grids.atom_grid = LEVELS[level]
    mf.grids.prune = None
    mf.grids.radi_method = radi.treutler
    mf.grids.becke_scheme = gen_grid.original_becke
    mf.grids.radii_adjust = radi.treutler_atomic_radii_adjust
    mf.small_rho_cutoff = 0.0
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.max_cycle = 200
    mf.kernel()
    assert mf.converged, (mol_name, level, variant)
    dm = mf.make_rdm1()
    n, exc, _ = mf._numint.nr_rks(mol, mf.grids, mf.xc, dm)
    return {"molecule": mol_name, "xyz": XYZ[mol_name], "grid": level, "vwn_variant": variant, "ri": ri, "xc": mf.xc,
            "npts": int(np.count_nonzero(mf.grids.weights)), "energy": float(mf.e_tot), "exc": float(exc),
            "nelec_grid": float(n), "num_electrons": int(mol.nelectron)}


def main():
    ref = {"generator": "scripts/qm_lda_pyscf_reference.py", "pyscf": pyscf.__version__,
           "libxc": libxc.__version__ if hasattr(libxc, "__version__") else "pyscf-bundled",
           "kernel": kernel_ref(), "energies": {}}
    for case, (mol, level, var, ri) in CASES.items():
        ref["energies"][case] = energy_ref(mol, level, var, ri)
        print(case, ref["energies"][case]["energy"], ref["energies"][case]["exc"], flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
