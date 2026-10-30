#!/usr/bin/env python3
"""PySCF/libxc references for the native B3LYP test (WP7, test_cases/qm_b3lyp/reference.json).

kernel:   libxc (through pyscf.dft.libxc.eval_xc, spin 0, deriv 1) on the (rho, sigma) table of the PBE test:
          29 densities rho = 1e-12 ... 1e2 (log-spaced) times 16 sigma values: 0, absolute 1e-30 ... 1e10, and sigma
          from the reduced gradient s = 0.01 ... 1e3 (sigma = 4 (3 pi^2)^(2/3) rho^(8/3) s^2), plus 8 pairs with
          sigma = 8e12 (largest sigma asked for, very large B88 x). Stored per kernel: eps, v_rho, v_sigma for
          GGA_X_B88 and GGA_C_LYP. libxc raises sigma to a floor before evaluating ("sigma_floor") and
          drops densities below its dens_threshold ("dens_threshold", measured here; floor = dens_threshold^(8/3)); the C++ driver evaluates
          curcuma at max(sigma, floor) for the comparison and also reports the raw sigma.
energies: RKS/def2-SVP (spherical) with the SAME grid as curcuma (qm_grid.cpp, WP4):
              grids.atom_grid = (n_rad, n_ang), prune = None, radi_method = radi.treutler,
              becke_scheme = gen_grid.original_becke, radii_adjust = radi.treutler_atomic_radii_adjust
          xc = 'HYB_GGA_XC_B3LYP5' (libxc id 475, VWN5; = curcuma -qm.vwn_variant vwn5, ORCA "B3LYP") or
               'HYB_GGA_XC_B3LYP'  (libxc id 402, VWN-RPA; = curcuma vwn_rpa, Gaussian, ORCA "B3LYP/G");
          the libxc ids are spelled out because PySCF's alias 'B3LYP' has changed meaning between versions.
          small_rho_cutoff = 0, conv_tol 1e-12, conv_tol_grad 1e-9. Stored: E_total, E_xc (nr_rks on the converged
          density, density-functional part only, as curcuma's E_xc), N on the grid, and the exact-exchange energy
          E_x^HF part = -(a_x/4) Tr(P K[P]) with PySCF's K (fitted K for rijk).
          RI cases: rij  = mf.density_fit(JFIT file, only_dfj=True): J fitted, K exact (= curcuma -qm.ri j);
                    rijk = mf.density_fit(JKFIT file): J and K fitted (= curcuma -qm.ri jk).
Geometries: Angstrom xyz files converted to Bohr with curcuma's CODATA 2018 factor (as qm_pbe_pyscf_reference.py).

Run:  ~/.venvs/qmref/bin/python scripts/qm_b3lyp_pyscf_reference.py
Claude Generated (Oct 2026, WP7).
"""
import importlib.util
import json
import os

import numpy as np
import pyscf
from pyscf import dft, gto
from pyscf.dft import gen_grid, libxc, radi

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.join(HERE, "..", "test_cases")
OUT = os.path.join(TC, "qm_b3lyp", "reference.json")
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
XC = {"vwn5": "HYB_GGA_XC_B3LYP5", "vwn_rpa": "HYB_GGA_XC_B3LYP"}
QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")
AUX = {"j": "def2-universal-JFIT.turbomole", "jk": "def2-universal-JKFIT.turbomole"}

# case: (molecule, grid level, ri mode, vwn variant)
CASES = {}
for _v in ("vwn5", "vwn_rpa"):
    for _m in ("H2O", "CH4", "NH3", "HF", "LiH", "Ne", "benzene"):
        CASES[f"{_m}_{_v}_sg1"] = (_m, "sg1", "off", _v)
CASES["H2O_vwn5_sg3"] = ("H2O", "sg3", "off", "vwn5")
CASES["NH3_vwn_rpa_sg2"] = ("NH3", "sg2", "off", "vwn_rpa")
CASES["H2O_vwn5_rij_sg1"] = ("H2O", "sg1", "j", "vwn5")
CASES["benzene_vwn5_rijk_sg1"] = ("benzene", "sg1", "jk", "vwn5")
CASES["H2O_vwn_rpa_rijk_sg1"] = ("H2O", "sg1", "jk", "vwn_rpa")


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


def eval_gga(name, rho, sigma):
    r = np.array([rho, np.zeros_like(rho), np.zeros_like(rho), np.sqrt(sigma)])
    exc, vxc = libxc.eval_xc(name, r, spin=0, deriv=1)[:2]
    return exc, vxc[0], vxc[1]


def dens_threshold(name):
    """Largest rho = 10^k at which libxc returns eps = 0 (density dropped)."""
    thr = 0.0
    for k in range(-30, -5):
        if eval_gga(name, np.array([10.0 ** k]), np.array([0.0]))[0][0] == 0.0:
            thr = 10.0 ** k
    return thr


def kernel_ref():
    rhos = np.logspace(-12, 2, 29)
    abs_sigma = [0.0, 1e-30, 1e-20, 1e-10, 1e-4, 1.0, 1e4, 1e10]
    s_values = [0.01, 0.1, 0.5, 1.0, 3.0, 10.0, 100.0, 1000.0]
    rho, sigma = [], []
    for r in rhos:
        for sg in abs_sigma:
            rho.append(r)
            sigma.append(sg)
        for s in s_values:
            rho.append(r)
            sigma.append(4.0 * (3.0 * np.pi ** 2) ** (2.0 / 3.0) * r ** (8.0 / 3.0) * s * s)
    for r in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1.0, 1e2):
        rho.append(r)
        sigma.append(8e12)
    rho, sigma = np.array(rho), np.array(sigma)
    out = {"rho": rho.tolist(), "sigma": sigma.tolist(),
           "dens_threshold": {"b88": dens_threshold("GGA_X_B88"), "lyp": dens_threshold("GGA_C_LYP")}}
    # libxc raises sigma to sigma_threshold^2, sigma_threshold = dens_threshold^(4/3) (as measured for PBE in WP6:
    # 1e-15 -> 1e-40). A direct scan does not work for LYP: its sigma term is ~exp(-c rho^(-1/3)) ~ 0 at small rho.
    out["sigma_floor"] = {k: v ** (8.0 / 3.0) for k, v in out["dens_threshold"].items()}
    for key, name in (("b88", "GGA_X_B88"), ("lyp", "GGA_C_LYP")):
        e, vr, vs = eval_gga(name, rho, sigma)
        out[key] = {"libxc": name, "eps": e.tolist(), "vrho": vr.tolist(), "vsigma": vs.tolist()}
    return out


def energy_ref(mol_name, level, ri, vwn):
    mol = read_mol(mol_name)
    mf = dft.RKS(mol)
    if ri != "off":
        mf = mf.density_fit(auxbasis=parse_turbomole(os.path.join(QMDIR, AUX[ri])), only_dfj=(ri == "j"))
    mf.xc = XC[vwn]
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
    assert mf.converged, (mol_name, level, ri, vwn)
    dm = mf.make_rdm1()
    n, exc, _ = mf._numint.nr_rks(mol, mf.grids, mf.xc, dm)
    hyb = float(libxc.hybrid_coeff(mf.xc))
    exx = -0.25 * hyb * float(np.einsum("ij,ji->", dm, mf.get_k(mol, dm)))
    return {"molecule": mol_name, "xyz": XYZ[mol_name], "grid": level, "ri": ri, "vwn": vwn, "xc": mf.xc,
            "hybrid_coeff": hyb, "exx": exx,
            "npts": int(np.count_nonzero(mf.grids.weights)), "energy": float(mf.e_tot), "exc": float(exc),
            "nelec_grid": float(n), "num_electrons": int(mol.nelectron)}


def main():
    ref = {"generator": "scripts/qm_b3lyp_pyscf_reference.py", "pyscf": pyscf.__version__,
           "libxc": libxc.__version__ if hasattr(libxc, "__version__") else "pyscf-bundled",
           "kernel": kernel_ref(), "energies": {}}
    print("sigma floors", ref["kernel"]["sigma_floor"], "density thresholds", ref["kernel"]["dens_threshold"],
          flush=True)
    for case, (mol, level, ri, vwn) in CASES.items():
        ref["energies"][case] = energy_ref(mol, level, ri, vwn)
        print(case, ref["energies"][case]["energy"], ref["energies"][case]["exc"], flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
