#!/usr/bin/env python3
"""PySCF references for the native Kohn-Sham gradient test (WP8 DFT, test_cases/qm_dft_grad/reference.json).

gradients: RKS/def2-SVP (spherical) with the SAME grid as curcuma (qm_grid.cpp, WP4; as the WP5-WP7 scripts):
              grids.atom_grid = (n_rad, n_ang), prune = None, radi_method = radi.treutler,
              becke_scheme = gen_grid.original_becke, radii_adjust = radi.treutler_atomic_radii_adjust,
              small_rho_cutoff = 0, conv_tol 1e-12, conv_tol_grad 1e-10
          xc = 'LDA_X,LDA_C_VWN' (lda vwn5) | 'LDA_X,LDA_C_VWN_RPA' (lda vwn_rpa) | 'GGA_X_PBE,GGA_C_PBE' (pbe)
               | 'HYB_GGA_XC_B3LYP5' (b3lyp vwn5) | 'HYB_GGA_XC_B3LYP' (b3lyp vwn_rpa)
          mf.nuc_grad_method() with grid_response False (PySCF default: basis-function part only) and True
          (Becke weight derivative + moving points, grids_response_becke). Stored in Eh/Bohr, atom-major.
          RI: rij = density_fit(JFIT, only_dfj=True) (J fitted, K exact, = curcuma -qm.ri j); rijk = density_fit(JKFIT)
          (= curcuma -qm.ri jk; for a pure functional only J is used). auxbasis_response True (PySCF default).
optimisations: the start structures test_cases/qm_grad/{H2O,formaldehyde}_distorted.xyz minimised with
          scipy.optimize.minimize(BFGS, gtol 1e-7 Eh/Bohr) on the PySCF energy and its grid_response=True gradient
          (sg1 grid as above; the grid moves with the atoms, the gradient is the exact derivative of that energy).
          Stored: minimum energy and geometry (Angstrom) per functional.
Geometries: Angstrom xyz files converted to Bohr with curcuma's CODATA 2018 factor.

--set mgga (WP11 stage 2): the meta-GGAs tpss ('MGGA_X_TPSS,MGGA_C_TPSS') and scan ('MGGA_X_SCAN,MGGA_C_SCAN'),
          cases CASES_MGGA / OPT_MGGA, written to test_cases/qm_mgga_grad/reference.json (same settings as above;
          grids "NxM" = n_rad x n_ang, as curcuma -qm.grid NxM).

Run:  ~/.venvs/qmref/bin/python scripts/qm_dft_grad_pyscf_reference.py [--only grad|opt] [--set ks|mgga]
Claude Generated (Oct 2026, WP8; --set mgga WP11 stage 2).
"""
import argparse
import importlib.util
import json
import os

import numpy as np
import pyscf
import scipy
import scipy.optimize
from pyscf import dft, gto
from pyscf.dft import gen_grid, radi

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.join(HERE, "..", "test_cases")
OUT = os.path.join(TC, "qm_dft_grad", "reference.json")
BOHR = 0.529177210903  # curcuma units.h
QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")
AUX = {"j": "def2-universal-JFIT.turbomole", "jk": "def2-universal-JKFIT.turbomole"}

XYZ = {
    "H2O": "qm_1e/H2O.xyz",
    "CH4": "qm_1e/CH4.xyz",
    "NH3": "qm_1e/NH3.xyz",
    "HF": "qm_1e/HF.xyz",
    "benzene": "qm_hf3c/benzene.xyz",
    "H2O_distorted": "qm_grad/H2O_distorted.xyz",
    "formaldehyde_distorted": "qm_grad/formaldehyde_distorted.xyz",
}
LEVELS = {"sg0": (23, 170), "sg1": (50, 194), "sg2": (75, 302), "sg3": (99, 590)}
XC = {("lda", "vwn5"): "LDA_X,LDA_C_VWN", ("lda", "vwn_rpa"): "LDA_X,LDA_C_VWN_RPA",
      ("pbe", "vwn5"): "GGA_X_PBE,GGA_C_PBE",
      ("b3lyp", "vwn5"): "HYB_GGA_XC_B3LYP5", ("b3lyp", "vwn_rpa"): "HYB_GGA_XC_B3LYP",
      ("tpss", "vwn5"): "MGGA_X_TPSS,MGGA_C_TPSS", ("scan", "vwn5"): "MGGA_X_SCAN,MGGA_C_SCAN"}

# case: (molecule, functional, vwn, grid, ri)
CASES = {}
for _m in ("H2O", "NH3", "CH4", "HF", "benzene"):
    for _f, _v in (("lda", "vwn5"), ("pbe", "vwn5"), ("b3lyp", "vwn5"), ("b3lyp", "vwn_rpa")):
        CASES[f"{_m}_{_f}_{_v}_sg1"] = (_m, _f, _v, "sg1", "off")
CASES["H2O_lda_vwn_rpa_sg1"] = ("H2O", "lda", "vwn_rpa", "sg1", "off")
CASES["H2O_b3lyp_vwn5_sg3"] = ("H2O", "b3lyp", "vwn5", "sg3", "off")
CASES["H2O_b3lyp_vwn5_rij_sg1"] = ("H2O", "b3lyp", "vwn5", "sg1", "j")
CASES["H2O_b3lyp_vwn5_rijk_sg1"] = ("H2O", "b3lyp", "vwn5", "sg1", "jk")
CASES["NH3_pbe_vwn5_rijk_sg1"] = ("NH3", "pbe", "vwn5", "sg1", "jk")
CASES["NH3_lda_vwn5_rij_sg1"] = ("NH3", "lda", "vwn5", "sg1", "j")
CASES["benzene_b3lyp_vwn_rpa_rijk_sg1"] = ("benzene", "b3lyp", "vwn_rpa", "sg1", "jk")

OPT = {f"{m}_{f}": (m, f) for m in ("H2O_distorted", "formaldehyde_distorted") for f in ("lda", "pbe", "b3lyp")}

# WP11 stage 2 (Claude Generated, Oct 2026): meta-GGA gradient references (vwn is unused, kept for the case format)
CASES_MGGA = {}
for _f in ("tpss", "scan"):
    for _m in ("H2O", "NH3", "CH4", "HF", "benzene"):
        CASES_MGGA[f"{_m}_{_f}_sg1"] = (_m, _f, "vwn5", "sg1", "off")
    CASES_MGGA[f"H2O_{_f}_sg3"] = ("H2O", _f, "vwn5", "sg3", "off")
    CASES_MGGA[f"H2O_{_f}_rij_sg1"] = ("H2O", _f, "vwn5", "sg1", "j")
    CASES_MGGA[f"NH3_{_f}_150x974"] = ("NH3", _f, "vwn5", "150x974", "off")
CASES_MGGA["benzene_tpss_rijk_sg1"] = ("benzene", "tpss", "vwn5", "sg1", "jk")
OPT_MGGA = {"H2O_distorted_tpss": ("H2O_distorted", "tpss"), "NH3_tpss": ("NH3", "tpss"),
            "formaldehyde_distorted_tpss": ("formaldehyde_distorted", "tpss")}


def parse_turbomole(path):
    spec = importlib.util.spec_from_file_location("rimp2_reference", os.path.join(HERE, "rimp2_reference.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.parse_turbomole(path)


def read_xyz(name):
    lines = open(os.path.join(TC, XYZ[name])).read().splitlines()
    n = int(lines[0].split()[0])
    sym, xyz = [], []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        sym.append(s)
        xyz.append([float(x), float(y), float(z)])
    return sym, np.array(xyz)


def make_mol(sym, xyz_ang):
    atoms = [[s, (np.array(x) / BOHR).tolist()] for s, x in zip(sym, xyz_ang)]
    return gto.M(atom=atoms, unit="Bohr", basis="def2-svp", cart=False, verbose=0)


def make_mf(mol, functional, vwn, level, ri):
    mf = dft.RKS(mol)
    if ri != "off":
        mf = mf.density_fit(auxbasis=parse_turbomole(os.path.join(QMDIR, AUX[ri])), only_dfj=(ri == "j"))
    mf.xc = XC[(functional, vwn)]
    mf.grids.atom_grid = LEVELS[level] if level in LEVELS else tuple(int(v) for v in level.split("x"))
    mf.grids.prune = None
    mf.grids.radi_method = radi.treutler
    mf.grids.becke_scheme = gen_grid.original_becke
    mf.grids.radii_adjust = radi.treutler_atomic_radii_adjust
    mf.small_rho_cutoff = 0.0
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 200
    return mf


def grad_ref(mol_name, functional, vwn, level, ri):
    sym, xyz = read_xyz(mol_name)
    mol = make_mol(sym, xyz)
    mf = make_mf(mol, functional, vwn, level, ri)
    mf.kernel()
    assert mf.converged, (mol_name, functional, vwn, level, ri)
    out = {"molecule": mol_name, "xyz": XYZ[mol_name], "functional": functional, "vwn": vwn, "grid": level,
           "ri": ri, "xc": mf.xc, "energy": float(mf.e_tot)}
    for resp in (False, True):
        g = mf.nuc_grad_method()
        g.grid_response = resp
        out["gradient_response" if resp else "gradient"] = g.kernel().reshape(-1).tolist()
    return out


def opt_ref(mol_name, functional):
    sym, xyz0 = read_xyz(mol_name)
    vwn = "vwn5"
    cache = {}

    def eg(x_bohr):
        key = tuple(np.round(x_bohr, 14))
        if key not in cache:
            mol = make_mol(sym, (x_bohr.reshape(-1, 3) * BOHR))
            mf = make_mf(mol, functional, vwn, "sg1", "off")
            mf.kernel()
            assert mf.converged
            g = mf.nuc_grad_method()
            g.grid_response = True
            cache.clear()
            cache[key] = (float(mf.e_tot), g.kernel().reshape(-1))
        return cache[key]

    res = scipy.optimize.minimize(lambda x: eg(x)[0], (xyz0 / BOHR).reshape(-1), jac=lambda x: eg(x)[1],
                                  method="BFGS", options={"gtol": 1e-7, "maxiter": 400})
    e, g = eg(res.x)
    return {"molecule": mol_name, "start": XYZ[mol_name], "functional": functional, "vwn": vwn, "grid": "sg1",
            "energy": e, "max_gradient": float(np.abs(g).max()), "iterations": int(res.nit),
            "symbols": sym, "geometry_angstrom": (res.x.reshape(-1, 3) * BOHR).tolist()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["grad", "opt"])
    ap.add_argument("--set", choices=["ks", "mgga"], default="ks")
    a = ap.parse_args()
    global OUT
    cases, opts = CASES, OPT
    if a.set == "mgga":
        OUT = os.path.join(TC, "qm_mgga_grad", "reference.json")
        cases, opts = CASES_MGGA, OPT_MGGA
    ref = {}
    if os.path.exists(OUT):
        ref = json.load(open(OUT))
    ref.update({"generator": "scripts/qm_dft_grad_pyscf_reference.py", "pyscf": pyscf.__version__,
                "scipy": scipy.__version__})
    if a.only in (None, "grad"):
        ref["gradients"] = {}
        for case, args in cases.items():
            ref["gradients"][case] = grad_ref(*args)
            r = ref["gradients"][case]
            d = np.abs(np.array(r["gradient"]) - np.array(r["gradient_response"])).max()
            print(f"{case:34s} E {r['energy']:.10f}  max|g_resp - g| {d:.2e}", flush=True)
    if a.only in (None, "opt"):
        ref["optimisations"] = {}
        for case, args in opts.items():
            ref["optimisations"][case] = opt_ref(*args)
            r = ref["optimisations"][case]
            print(f"opt {case:30s} E {r['energy']:.10f} max|g| {r['max_gradient']:.1e} ({r['iterations']} it)",
                  flush=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=1)
    print("written", OUT)


if __name__ == "__main__":
    main()
