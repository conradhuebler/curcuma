#!/usr/bin/env python3
"""PySCF/libxc references for the native double hybrids (WP12 stage 1, test_cases/qm_dh/reference.json).

PySCF 2.14 / libxc 7.0.0 have no double-hybrid definitions (no B2PLYP/DSD in XC_CODES or PySCF's aliases), so the
reference is assembled from components, each with the same definition as curcuma's table
(src/core/energy_calculators/qm_methods/double_hybrid_method.cpp):

kernel:   libxc (eval_xc, spin 0, deriv 1) of GGA_C_P86 (PZ81 local part), GGA_C_P86VWN (VWN5 local part) on
            rho = 1e-12 .. 1e2 (11 values) x s = 0, 0.01, 0.3, 1, 3, 30   (sigma = 4 (3 pi^2)^(2/3) rho^(8/3) s^2)
          and MGGA_C_BC95 on the same (rho, s) times alpha = 0, 1e-3, 0.1, 0.5, 1, 2, 10, 1e3
          (tau = tau_W + alpha tau_unif) plus two points with tau < tau_W per (rho, s > 0) (libxc and curcuma lower
          sigma to 8 rho tau). sigma raised to libxc's floor (dens_threshold^(4/3))^2 (1e-40 for P86, 1e-14^(8/3) for
          BC95), stored as sigma_floor per kernel.
energies: RKS/def2-SVP (spherical) on curcuma's grid (Treutler x Lebedev, Becke, no pruning, Treutler radii adjust),
          small_rho_cutoff 0, conv_tol 1e-12, conv_tol_grad 1e-9, xc = 'a_x*HF + w_x*GGA_X_..., w_c*<corr>' (the KS part
          of the row), then with the converged KS orbitals and eigenvalues
            RI-PT2: B = L^-1 (P|mn) with curcuma's def2-universal-JKFIT.turbomole (as scripts/rimp2_reference.py),
                    E_os = sum (ia|jb)^2 / D, E_ss = sum (ia|jb) [(ia|jb) - (ib|ja)] / D, full and frozen core
            canonical PT2: pyscf.mp.MP2(mf) (e_corr_os, e_corr_ss) for H2O and benzene
          Stored unscaled; curcuma's E = E_KS + c_os E_os + c_ss E_ss is checked component by component.
Geometries: Angstrom xyz -> Bohr with curcuma's CODATA 2018 factor (as qm_mgga_pyscf_reference.py).

Run:  ~/.venvs/qmref/bin/python scripts/qm_dh_pyscf_reference.py
Claude Generated (Oct 2026, WP12 stage 1).
"""
import importlib.util
import json
import math
import os

import numpy as np
import pyscf
from pyscf import dft, gto, mp
from pyscf.dft import gen_grid, libxc, radi

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.join(HERE, "..", "test_cases")
OUT = os.path.join(TC, "qm_dh", "reference.json")
QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")
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
LEVELS = {"sg1": (50, 194), "sg3": (99, 590)}
# KS part of every table row (curcuma method name -> PySCF xc string); the -d3bj rows of b2plyp/b2gp-plyp share
# the KS part of the plain names.
XC = {
    "b2plyp": "0.53*HF + 0.47*GGA_X_B88, 0.73*GGA_C_LYP",
    "b2gp-plyp": "0.65*HF + 0.35*GGA_X_B88, 0.64*GGA_C_LYP",
    "dsd-blyp": "0.75*HF + 0.25*GGA_X_B88, 0.53*GGA_C_LYP",
    "dsd-blyp-d3bj": "0.71*HF + 0.29*GGA_X_B88, 0.54*GGA_C_LYP",
    "dsd-pbep86": "0.72*HF + 0.28*GGA_X_PBE, 0.44*GGA_C_P86",
    "dsd-pbep86-d3bj": "0.69*HF + 0.31*GGA_X_PBE, 0.44*GGA_C_P86",
    "dsd-pbeb95-d3bj": "0.66*HF + 0.34*GGA_X_PBE, 0.55*MGGA_C_BC95",
    "revdsd-pbep86": "0.69*HF + 0.31*GGA_X_PBE, 0.4224*GGA_C_P86",
    # -qm.p86_lda vwn5 (ORCA's P86)
    "dsd-pbep86@vwn5": "0.72*HF + 0.28*GGA_X_PBE, 0.44*GGA_C_P86VWN",
}

CASES = {}
for fn in ("b2plyp", "b2gp-plyp", "dsd-blyp", "dsd-pbep86", "dsd-pbeb95-d3bj", "revdsd-pbep86"):
    for mol in XYZ:
        if fn == "dsd-pbeb95-d3bj" and mol == "LiH":
            continue  # PySCF's SCF stalls at |g| ~ 6e-5 with MGGA_C_BC95 (300 cycles, also newton()): no reference
        CASES[f"{mol}_{fn}_sg1"] = (fn, mol, "sg1")
for fn in ("dsd-blyp-d3bj", "dsd-pbep86-d3bj", "dsd-pbep86@vwn5"):
    for mol in ("H2O", "NH3", "benzene"):
        CASES[f"{mol}_{fn}_sg1"] = (fn, mol, "sg1")
for fn in ("b2plyp", "dsd-pbeb95-d3bj"):
    CASES[f"H2O_{fn}_sg3"] = (fn, "H2O", "sg3")
CANONICAL = ("H2O", "benzene")


def load_rimp2():
    spec = importlib.util.spec_from_file_location("rimp2_reference", os.path.join(HERE, "rimp2_reference.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_mol(name):
    lines = open(os.path.join(TC, XYZ[name])).read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        atoms.append([s, [float(x) / BOHR, float(y) / BOHR, float(z) / BOHR]])
    return gto.M(atom=atoms, unit="Bohr", basis="def2-svp", cart=False, verbose=0), atoms


def kernel_ref():
    rhos = [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0]
    s_values = [0.0, 0.01, 0.3, 1.0, 3.0, 30.0]
    alphas = [0.0, 1e-3, 0.1, 0.5, 1.0, 2.0, 10.0, 1e3]
    cs = (3.0 * math.pi ** 2) ** (2.0 / 3.0)
    g_rho, g_sig = [], []
    m_rho, m_sig, m_tau, m_kind = [], [], [], []
    for r in rhos:
        for s in s_values:
            sg = 4.0 * cs * r ** (8.0 / 3.0) * s * s
            g_rho.append(r)
            g_sig.append(sg)
            tw, tu = sg / (8.0 * r), 0.3 * cs * r ** (5.0 / 3.0)
            for a in alphas:
                m_rho.append(r)
                m_sig.append(sg)
                m_tau.append(tw + a * tu)
                m_kind.append("admissible")
            if s > 0.0:
                for f in (0.5, 0.9):
                    m_rho.append(r)
                    m_sig.append(sg)
                    m_tau.append(f * tw)
                    m_kind.append("clamp")
    out = {"gga": {"rho": g_rho, "sigma": g_sig}, "mgga": {"rho": m_rho, "sigma": m_sig, "tau": m_tau, "kind": m_kind}}
    gr, gs = np.array(g_rho), np.array(g_sig)
    for key, name, floor in (("p86", "GGA_C_P86", 1e-40), ("p86vwn", "GGA_C_P86VWN", 1e-40)):
        sg = np.maximum(gs, floor)
        r = np.array([gr, np.sqrt(sg), np.zeros_like(gr), np.zeros_like(gr)])
        exc, vxc = libxc.eval_xc(name, r, spin=0, deriv=1)[:2]
        out["gga"][key] = {"libxc": name, "sigma_floor": floor, "eps": exc.tolist(), "vrho": vxc[0].tolist(),
                           "vsigma": vxc[1].tolist()}
    mr, ms, mt = np.array(m_rho), np.array(m_sig), np.array(m_tau)
    floor = (1e-14 ** (4.0 / 3.0)) ** 2  # MGGA_C_BC95 dens_threshold 1e-14 (libxc mgga_c_bc95.c)
    sg = np.maximum(ms, floor)
    r = np.array([mr, np.sqrt(sg), np.zeros_like(mr), np.zeros_like(mr), mt])
    exc, vxc = libxc.eval_xc("MGGA_C_BC95", r, spin=0, deriv=1)[:2]
    out["mgga"]["b95"] = {"libxc": "MGGA_C_BC95", "sigma_floor": floor, "eps": exc.tolist(), "vrho": vxc[0].tolist(),
                          "vsigma": vxc[1].tolist(), "vtau": vxc[3].tolist()}
    return out


def rks(mol, xc, level):
    mf = dft.RKS(mol)
    mf.xc = xc
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
    return mf


def energy_ref(rim, aux, fn, mol_name, level):
    mol, atoms = read_mol(mol_name)
    mf = rks(mol, XC[fn], level)
    assert mf.converged, (fn, mol_name, level)
    dm = mf.make_rdm1()
    _, exc, _ = mf._numint.nr_rks(mol, mf.grids, mf.xc, dm)
    nfc = rim.frozen_core([(a[0], None) for a in atoms])
    out = {"functional": fn, "molecule": mol_name, "xyz": XYZ[mol_name], "grid": level, "xc": mf.xc,
           "e_ks": float(mf.e_tot), "exc": float(exc), "nbf": mol.nao_nr(), "n_frozen_core": nfc}
    for tag, frz in (("full", 0), ("fc", nfc)):
        eos, ess, naux = rim.manual_rimp2(mol, mf, aux, frz)
        out["ri_" + tag] = {"e_os": float(eos), "e_ss": float(ess), "naux": int(naux)}
    if mol_name in CANONICAL and level == "sg1":
        for tag, frz in (("full", None), ("fc", nfc)):
            pt = mp.MP2(mf, frozen=frz)
            pt.kernel()
            out["canonical_" + tag] = {"e_os": float(pt.e_corr_os), "e_ss": float(pt.e_corr_ss)}
    return out


def main():
    rim = load_rimp2()
    aux = rim.parse_turbomole(os.path.join(QMDIR, "def2-universal-JKFIT.turbomole"))
    ref = {"generator": "scripts/qm_dh_pyscf_reference.py", "pyscf": pyscf.__version__,
           "libxc": libxc.__version__ if hasattr(libxc, "__version__") else "pyscf-bundled",
           "aux": "def2-universal-JKFIT.turbomole", "kernel": kernel_ref(), "energies": {}}
    for case, (fn, mol, level) in CASES.items():
        e = energy_ref(rim, aux, fn, mol, level)
        ref["energies"][case] = e
        print(f"{case:32s} E_KS {e['e_ks']:.10f} OS {e['ri_full']['e_os']:.10f} SS {e['ri_full']['e_ss']:.10f}",
              flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, separators=(",", ":"))
    print("wrote", OUT, os.path.getsize(OUT), "bytes")


if __name__ == "__main__":
    main()
