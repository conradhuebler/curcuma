#!/usr/bin/env python3
"""PySCF/libxc references for the native PBE test (WP6, test_cases/qm_pbe/reference.json).

kernel:   libxc (through pyscf.dft.libxc.eval_xc, spin 0, deriv 1) on a (rho, sigma) table: 29 densities
          rho = 1e-12 ... 1e2 (log-spaced) times 16 sigma values: 0, absolute 1e-30 ... 1e10, and sigma from the
          reduced gradient s = 0.01 ... 1e3 (sigma = 4 (3 pi^2)^(2/3) rho^(8/3) s^2; physical s is ~0 ... 3).
          Stored per kernel: eps, v_rho, v_sigma for GGA_X_PBE, GGA_C_PBE, and eps, v_rho for LDA_C_PW_MOD.
          libxc raises sigma to a floor before evaluating (sigma_threshold^2 with sigma_threshold =
          dens_threshold^(4/3): 1e-40 for GGA_X_PBE, 1e-32 for GGA_C_PBE, measured here and stored as "sigma_floor");
          the C++ driver evaluates curcuma at max(sigma, floor) for the comparison and also reports the raw sigma = 0.
energies: RKS/def2-SVP (spherical) with the SAME grid as curcuma (qm_grid.cpp, WP4):
              grids.atom_grid = (n_rad, n_ang), prune = None, radi_method = radi.treutler,
              becke_scheme = gen_grid.original_becke, radii_adjust = radi.treutler_atomic_radii_adjust
          xc = 'GGA_X_PBE,GGA_C_PBE' (= PySCF 'PBE'); small_rho_cutoff = 0, conv_tol 1e-12, conv_tol_grad 1e-9.
          Stored: E_total, E_xc (nr_rks on the converged density), N on the grid.
          RI cases: mf.density_fit(auxbasis=<curcuma's def2-universal-JFIT/-JKFIT.turbomole>) (J only for a GGA).
Geometries: Angstrom xyz files converted to Bohr with curcuma's CODATA 2018 factor (as qm_lda_pyscf_reference.py).

Run:  ~/.venvs/qmref/bin/python scripts/qm_pbe_pyscf_reference.py
Claude Generated (Oct 2026, WP6).
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
OUT = os.path.join(TC, "qm_pbe", "reference.json")
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
XC = "GGA_X_PBE,GGA_C_PBE"
QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")
AUX = {"j": "def2-universal-JFIT.turbomole", "jk": "def2-universal-JKFIT.turbomole"}

# case: (molecule, grid level, ri mode)
CASES = {
    "H2O_sg1": ("H2O", "sg1", "off"),
    "CH4_sg1": ("CH4", "sg1", "off"),
    "NH3_sg1": ("NH3", "sg1", "off"),
    "HF_sg1": ("HF", "sg1", "off"),
    "LiH_sg1": ("LiH", "sg1", "off"),
    "Ne_sg1": ("Ne", "sg1", "off"),
    "benzene_sg1": ("benzene", "sg1", "off"),
    "H2O_sg3": ("H2O", "sg3", "off"),
    "NH3_sg2": ("NH3", "sg2", "off"),
    "H2O_rij_sg1": ("H2O", "sg1", "j"),
    "benzene_rijk_sg1": ("benzene", "sg1", "jk"),
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


def eval_gga(name, rho, sigma):
    r = np.array([rho, np.zeros_like(rho), np.zeros_like(rho), np.sqrt(sigma)])
    exc, vxc = libxc.eval_xc(name, r, spin=0, deriv=1)[:2]
    return exc, vxc[0], vxc[1]


def sigma_floor(name):
    """Largest sigma at which libxc still returns the sigma = 0 value (floor = 10^k, k found by scan)."""
    rho = np.array([1e-12])
    e0 = eval_gga(name, rho, np.array([0.0]))[0][0]
    floor = 0.0
    for k in range(-60, -10):
        if eval_gga(name, rho, np.array([10.0 ** k]))[0][0] == e0:
            floor = 10.0 ** k
    return floor


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
    rho, sigma = np.array(rho), np.array(sigma)
    out = {"rho": rho.tolist(), "sigma": sigma.tolist(),
           "sigma_floor": {"pbe_x": sigma_floor("GGA_X_PBE"), "pbe_c": sigma_floor("GGA_C_PBE")}}
    for key, name in (("pbe_x", "GGA_X_PBE"), ("pbe_c", "GGA_C_PBE")):
        e, vr, vs = eval_gga(name, rho, sigma)
        out[key] = {"libxc": name, "eps": e.tolist(), "vrho": vr.tolist(), "vsigma": vs.tolist()}
    exc, vxc = libxc.eval_xc("LDA_C_PW_MOD", rho, spin=0, deriv=1)[:2]
    out["pw92"] = {"libxc": "LDA_C_PW_MOD", "eps": exc.tolist(), "vrho": vxc[0].tolist()}
    return out


def energy_ref(mol_name, level, ri):
    mol = read_mol(mol_name)
    mf = dft.RKS(mol)
    if ri != "off":
        mf = mf.density_fit(auxbasis=parse_turbomole(os.path.join(QMDIR, AUX[ri])))
    mf.xc = XC
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
    assert mf.converged, (mol_name, level, ri)
    dm = mf.make_rdm1()
    n, exc, _ = mf._numint.nr_rks(mol, mf.grids, mf.xc, dm)
    return {"molecule": mol_name, "xyz": XYZ[mol_name], "grid": level, "ri": ri, "xc": mf.xc,
            "npts": int(np.count_nonzero(mf.grids.weights)), "energy": float(mf.e_tot), "exc": float(exc),
            "nelec_grid": float(n), "num_electrons": int(mol.nelectron)}


def main():
    ref = {"generator": "scripts/qm_pbe_pyscf_reference.py", "pyscf": pyscf.__version__,
           "libxc": libxc.__version__ if hasattr(libxc, "__version__") else "pyscf-bundled",
           "kernel": kernel_ref(), "energies": {}}
    print("sigma floors", ref["kernel"]["sigma_floor"], flush=True)
    for case, (mol, level, ri) in CASES.items():
        ref["energies"][case] = energy_ref(mol, level, ri)
        print(case, ref["energies"][case]["energy"], ref["energies"][case]["exc"], flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
