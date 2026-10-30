#!/usr/bin/env python3
"""PySCF/libxc references for the native meta-GGAs TPSS and SCAN (WP11 stage 1, test_cases/qm_mgga/reference.json).

kernel:   libxc (pyscf.dft.libxc.eval_xc, spin 0, deriv 1, input rows [rho, grad_x, grad_y, grad_z, tau]) for
          MGGA_X_TPSS, MGGA_C_TPSS, MGGA_X_SCAN, MGGA_C_SCAN on a (rho, sigma, tau) table built from
            rho   = 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100
            s     = 0, 0.01, 0.3, 1, 3, 30          (sigma = 4 (3 pi^2)^(2/3) rho^(8/3) s^2)
            alpha = 0, 1e-3, 0.1, 0.5, 0.9, 0.99, 0.999, 1, 1.001, 1.01, 1.1, 1.5, 2, 5, 10, 100, 1e3 and both sides
                    (1 -+ 1e-9 relative) of the four SCAN interpolation cut-offs (libxc mgga_x_scan.mpl)
            tau   = tau_W + alpha tau_unif,  tau_W = sigma / (8 rho),  tau_unif = (3/10) (3 pi^2)^(2/3) rho^(5/3)
          (physically admissible: tau >= tau_W), plus a "clamp" set with tau = 0.5 tau_W and 0.9 tau_W (not
          admissible; libxc lowers sigma to 8 rho tau, curcuma does the same). libxc raises sigma to 1e-40
          (sigma_threshold^2, dens_threshold 1e-15 of these functionals, libxc 7 functionals.c); stored as
          "sigma_floor", the C++ driver evaluates curcuma at max(sigma, floor). Tau = 0 (s = 0, alpha = 0) is raised to
          1e-20 by both (tau_threshold).
          Stored per kernel: eps, v_rho, v_sigma, v_tau. Also eps and v_rho of LDA_C_PW_MOD (scale of the
          correlation errors, as in the PBE test).
energies: RKS/def2-SVP (spherical) with the SAME grid as curcuma (qm_grid.cpp, WP4):
              grids.atom_grid = (n_rad, n_ang), prune = None, radi_method = radi.treutler,
              becke_scheme = gen_grid.original_becke, radii_adjust = radi.treutler_atomic_radii_adjust
          xc = 'MGGA_X_TPSS,MGGA_C_TPSS' or 'MGGA_X_SCAN,MGGA_C_SCAN'; small_rho_cutoff = 0, conv_tol 1e-12,
          conv_tol_grad 1e-9, max_cycle 300. Stored: E_total, E_xc (nr_rks on the converged density), N on the grid.
          RI cases: mf.density_fit(auxbasis=<curcuma's def2-universal-JFIT/-JKFIT.turbomole>) (J only, pure DFT).
Geometries: Angstrom xyz files converted to Bohr with curcuma's CODATA 2018 factor (as qm_pbe_pyscf_reference.py).

Run:  ~/.venvs/qmref/bin/python scripts/qm_mgga_pyscf_reference.py
Claude Generated (Oct 2026, WP11 stage 1).
"""
import importlib.util
import json
import math
import os

import numpy as np
import pyscf
from pyscf import dft, gto
from pyscf.dft import gen_grid, libxc, radi

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.join(HERE, "..", "test_cases")
OUT = os.path.join(TC, "qm_mgga", "reference.json")
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
XC = {"tpss": "MGGA_X_TPSS,MGGA_C_TPSS", "scan": "MGGA_X_SCAN,MGGA_C_SCAN"}
QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")
AUX = {"j": "def2-universal-JFIT.turbomole", "jk": "def2-universal-JKFIT.turbomole"}
KERNELS = (("tpss_x", "MGGA_X_TPSS"), ("tpss_c", "MGGA_C_TPSS"), ("scan_x", "MGGA_X_SCAN"), ("scan_c", "MGGA_C_SCAN"))
SIGMA_FLOOR = 1e-40  # (dens_threshold^(4/3))^2, dens_threshold = 1e-15 (libxc 7.0.0 functionals.c, work_mgga_inc.c)

# case: (functional, molecule, grid level, ri mode)
CASES = {}
for fn in ("tpss", "scan"):
    for mol in XYZ:
        CASES[f"{mol}_{fn}_sg1"] = (fn, mol, "sg1", "off")
    for mol in ("H2O", "NH3"):
        for lv in ("sg2", "sg3"):
            CASES[f"{mol}_{fn}_{lv}"] = (fn, mol, lv, "off")
    CASES[f"benzene_{fn}_sg3"] = (fn, "benzene", "sg3", "off")
    CASES[f"H2O_{fn}_rij_sg1"] = (fn, "H2O", "sg1", "j")
    CASES[f"benzene_{fn}_rijk_sg1"] = (fn, "benzene", "sg1", "jk")


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


def scan_cutoffs():
    """The four alpha cut-offs of libxc's SCAN interpolation functions (mgga_x_scan.mpl)."""
    L = -math.log(np.finfo(float).eps)
    out = []
    for c1, c2, d in ((0.667, 0.8, 1.24), (0.64, 1.5, 0.7)):
        Ld = -math.log(np.finfo(float).eps / d)
        out += [L / (L + c1), (Ld + c2) / Ld]
    return out


def eval_mgga(name, rho, sigma, tau):
    r = np.array([rho, np.sqrt(sigma), np.zeros_like(rho), np.zeros_like(rho), tau])
    exc, vxc = libxc.eval_xc(name, r, spin=0, deriv=1)[:2]
    return exc, vxc[0], vxc[1], vxc[3]


def kernel_ref():
    rhos = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0]
    s_values = [0.0, 0.01, 0.3, 1.0, 3.0, 30.0]
    alphas = [0.0, 1e-3, 0.1, 0.5, 0.9, 0.99, 0.999, 1.0, 1.001, 1.01, 1.1, 1.5, 2.0, 5.0, 10.0, 100.0, 1e3]
    for c in scan_cutoffs():
        alphas += [c * (1 - 1e-9), c * (1 + 1e-9)]
    cs = (3.0 * math.pi ** 2) ** (2.0 / 3.0)
    rho, sigma, tau, alpha, kind = [], [], [], [], []
    for r in rhos:
        for s in s_values:
            sg = 4.0 * cs * r ** (8.0 / 3.0) * s * s
            tw, tu = sg / (8.0 * r), 0.3 * cs * r ** (5.0 / 3.0)
            for a in alphas:
                rho.append(r)
                sigma.append(sg)
                tau.append(tw + a * tu)
                alpha.append(a)
                kind.append("admissible")
            if s > 0.0:
                for f in (0.5, 0.9):
                    rho.append(r)
                    sigma.append(sg)
                    tau.append(f * tw)
                    alpha.append(-1.0)
                    kind.append("clamp")
    rho, sigma, tau = np.array(rho), np.array(sigma), np.array(tau)
    out = {"rho": rho.tolist(), "sigma": sigma.tolist(), "tau": tau.tolist(), "alpha": alpha, "kind": kind,
           "sigma_floor": SIGMA_FLOOR, "scan_cutoffs": scan_cutoffs()}
    sg_eval = np.maximum(sigma, SIGMA_FLOOR)
    for key, name in KERNELS:
        e, vr, vs, vt = eval_mgga(name, rho, sg_eval, tau)
        out[key] = {"libxc": name, "eps": e.tolist(), "vrho": vr.tolist(), "vsigma": vs.tolist(), "vtau": vt.tolist()}
    exc, vxc = libxc.eval_xc("LDA_C_PW_MOD", rho, spin=0, deriv=1)[:2]
    out["pw92"] = {"libxc": "LDA_C_PW_MOD", "eps": exc.tolist(), "vrho": vxc[0].tolist()}
    return out


def energy_ref(fn, mol_name, level, ri):
    mol = read_mol(mol_name)
    mf = dft.RKS(mol)
    if ri != "off":
        mf = mf.density_fit(auxbasis=parse_turbomole(os.path.join(QMDIR, AUX[ri])))
    mf.xc = XC[fn]
    mf.grids.atom_grid = LEVELS[level]
    mf.grids.prune = None
    mf.grids.radi_method = radi.treutler
    mf.grids.becke_scheme = gen_grid.original_becke
    mf.grids.radii_adjust = radi.treutler_atomic_radii_adjust
    mf.small_rho_cutoff = 0.0
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.max_cycle = 300
    mf.kernel()
    assert mf.converged, (fn, mol_name, level, ri)
    dm = mf.make_rdm1()
    n, exc, _ = mf._numint.nr_rks(mol, mf.grids, mf.xc, dm)
    return {"functional": fn, "molecule": mol_name, "xyz": XYZ[mol_name], "grid": level, "ri": ri, "xc": mf.xc,
            "npts": int(np.count_nonzero(mf.grids.weights)), "energy": float(mf.e_tot), "exc": float(exc),
            "nelec_grid": float(n), "num_electrons": int(mol.nelectron)}


def main():
    ref = {"generator": "scripts/qm_mgga_pyscf_reference.py", "pyscf": pyscf.__version__,
           "libxc": libxc.__version__ if hasattr(libxc, "__version__") else "pyscf-bundled",
           "kernel": kernel_ref(), "energies": {}}
    print("kernel points", len(ref["kernel"]["rho"]), flush=True)
    for case, (fn, mol, level, ri) in CASES.items():
        ref["energies"][case] = energy_ref(fn, mol, level, ri)
        print(case, ref["energies"][case]["energy"], ref["energies"][case]["exc"], flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=0)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
