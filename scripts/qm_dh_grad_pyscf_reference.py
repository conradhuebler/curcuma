#!/usr/bin/env python3
"""PySCF/libxc references for the double-hybrid analytic gradient (WP12 stage 2, test_cases/qm_dh_grad/reference.json).

kernel:     libxc 7.0.0 (PySCF eval_xc, spin 0, deriv 2) second derivatives v2rho2, v2rhosigma, v2sigma2 (GGA) and
            additionally v2rhotau, v2sigmatau, v2tau2 (MGGA_C_BC95) of GGA_X_B88, GGA_C_LYP, GGA_X_PBE, GGA_C_P86,
            GGA_C_P86VWN, MGGA_C_BC95 on rho = 1e-6 .. 1e2 x s = 0.01 .. 10 x alpha = 0.1, 1, 5 (tau = tau_W + alpha
            tau_unif; GGA kernels ignore tau); sigma > 0 everywhere (P86's d2/dsigma2 diverges at sigma = 0).
gradients:  PySCF has no double-hybrid gradient (no DH in PySCF/libxc). Reference = central 4-point finite difference
            (h = 1e-3 Bohr) of the PySCF double-hybrid ENERGY E = E_KS + c_os E_os + c_ss E_ss assembled from
            components exactly as scripts/qm_dh_pyscf_reference.py (RKS on curcuma's grid, which PySCF rebuilds at every
            displaced geometry = moving grid; RI-PT2 with def2-universal-JKFIT or pyscf.mp.MP2), i.e. an independent
            program's energy surface. curcuma's analytic gradient WITH grid response is compared with it.
Run:  ~/.venvs/qmref/bin/python scripts/qm_dh_grad_pyscf_reference.py [case ...]   (cases: update only those)
Claude Generated (Oct 2026, WP12 stage 2).
"""
import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
import pyscf
from pyscf import gto, mp
from pyscf.dft import libxc

HERE = os.path.dirname(os.path.abspath(__file__))
TC = os.path.join(HERE, "..", "test_cases")
OUT = os.path.join(TC, "qm_dh_grad", "reference.json")
QMDIR = os.path.join(HERE, "..", "src", "core", "energy_calculators", "qm_methods")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


DH = load("qm_dh_pyscf_reference", os.path.join(HERE, "qm_dh_pyscf_reference.py"))
RIM = load("rimp2_reference", os.path.join(HERE, "rimp2_reference.py"))
# c_os, c_ss of the table rows (double_hybrid_method.cpp)
SCALE = {"b2plyp": (0.27, 0.27), "b2gp-plyp": (0.36, 0.36), "dsd-blyp": (0.46, 0.60), "dsd-pbep86": (0.51, 0.36),
         "dsd-pbeb95-d3bj": (0.46, 0.09), "revdsd-pbep86": (0.5935, 0.0566)}
# (case, functional, molecule, pt2 ri|canonical, frozen core)
GRAD_CASES = [
    ("H2O_b2plyp_ri", "b2plyp", "H2O", "ri", False),
    ("H2O_b2plyp_ri_fc", "b2plyp", "H2O", "ri", True),
    ("H2O_b2plyp_canonical", "b2plyp", "H2O", "canonical", False),
    ("NH3_dsd-blyp_ri", "dsd-blyp", "NH3", "ri", False),
    ("H2O_dsd-pbep86_ri", "dsd-pbep86", "H2O", "ri", False),
    ("HF_revdsd-pbep86_ri_fc", "revdsd-pbep86", "HF", "ri", True),
    ("H2O_dsd-pbeb95-d3bj_ri", "dsd-pbeb95-d3bj", "H2O", "ri", False),
    ("CH4_b2gp-plyp_canonical_fc", "b2gp-plyp", "CH4", "canonical", True),
    # P86 on VWN5 (smooth): PZ81 (libxc LDA_C_PZ, used by GGA_C_P86) is discontinuous at r_s = 1, so a finite
    # difference of a PZ81-P86 energy picks up jumps of ~1e-9 Eh per grid point crossing the seam (WP12 doc, stage 2)
    ("H2O_dsd-pbep86_vwn5_ri", "dsd-pbep86@vwn5", "H2O", "ri", False),
    ("HF_revdsd-pbep86_vwn5_ri_fc", "revdsd-pbep86@vwn5", "HF", "ri", True),
]
DH.XC["revdsd-pbep86@vwn5"] = "0.69*HF + 0.31*GGA_X_PBE, 0.4224*GGA_C_P86VWN"


def kernel_ref():
    rhos = [1e-6, 1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0]
    s_values = [0.01, 0.3, 1.0, 3.0, 10.0]
    alphas = [0.1, 1.0, 5.0]
    cs = (3.0 * math.pi ** 2) ** (2.0 / 3.0)
    R, S, T = [], [], []
    for r in rhos:
        for s in s_values:
            sg = 4.0 * cs * r ** (8.0 / 3.0) * s * s
            tw, tu = sg / (8.0 * r), 0.3 * cs * r ** (5.0 / 3.0)
            for a in alphas:
                R.append(r)
                S.append(sg)
                T.append(tw + a * tu)
    out = {"rho": R, "sigma": S, "tau": T}
    r, s, t = np.array(R), np.array(S), np.array(T)
    z = np.zeros_like(r)
    for key, name in (("b88", "GGA_X_B88"), ("lyp", "GGA_C_LYP"), ("pbex", "GGA_X_PBE"), ("p86", "GGA_C_P86"),
                      ("p86vwn", "GGA_C_P86VWN")):
        exc, vxc, fxc = libxc.eval_xc(name, np.array([r, np.sqrt(s), z, z]), spin=0, deriv=2)[:3]
        out[key] = {"libxc": name, "v": [vxc[0].tolist(), vxc[1].tolist()],
                    "h": {"rr": fxc[0].tolist(), "rs": fxc[1].tolist(), "ss": fxc[2].tolist()}}
    exc, vxc, fxc = libxc.eval_xc("MGGA_C_BC95", np.array([r, np.sqrt(s), z, z, z, t]), spin=0, deriv=2)[:3]
    out["b95"] = {"libxc": "MGGA_C_BC95", "v": [vxc[0].tolist(), vxc[1].tolist(), vxc[3].tolist()],
                  "h": {"rr": fxc[0].tolist(), "rs": fxc[1].tolist(), "ss": fxc[2].tolist(), "tt": fxc[4].tolist(),
                        "rt": fxc[6].tolist(), "st": fxc[9].tolist()}}
    return out


def dh_energy(atoms, fn, pt2, fc, aux):
    mol = gto.M(atom=atoms, unit="Bohr", basis="def2-svp", cart=False, verbose=0)
    mf = DH.rks(mol, DH.XC[fn], "sg1")
    assert mf.converged, (fn, atoms)
    nfc = RIM.frozen_core([(a[0], None) for a in atoms]) if fc else 0
    if pt2 == "ri":
        eos, ess, _ = RIM.manual_rimp2(mol, mf, aux, nfc)
    else:
        pt = mp.MP2(mf, frozen=nfc if nfc else None)
        pt.kernel()
        eos, ess = pt.e_corr_os, pt.e_corr_ss
    cos, css = SCALE[fn.split("@")[0]]
    return float(mf.e_tot + cos * eos + css * ess)


def fd_gradient(name, fn, mol_name, pt2, fc, aux, h=1e-3):
    _, atoms = DH.read_mol(mol_name)
    e0 = dh_energy(atoms, fn, pt2, fc, aux)
    g = []
    for i in range(len(atoms)):
        for k in range(3):
            def e(d):
                at = [[a[0], list(a[1])] for a in atoms]
                at[i][1][k] += d
                return dh_energy(at, fn, pt2, fc, aux)
            g.append((-e(2 * h) + 8 * e(h) - 8 * e(-h) + e(-2 * h)) / (12 * h))
    return {"functional": fn, "molecule": mol_name, "xyz": DH.XYZ[mol_name], "pt2": pt2, "frozen_core": fc,
            "grid": "sg1", "energy": e0, "fd_step_bohr": h, "gradient": g}


def main():
    aux = RIM.parse_turbomole(os.path.join(QMDIR, "def2-universal-JKFIT.turbomole"))
    ref = {"generator": "scripts/qm_dh_grad_pyscf_reference.py", "pyscf": pyscf.__version__,
           "libxc": getattr(libxc, "__version__", "pyscf-bundled"), "kernel": kernel_ref(), "gradients": {}}
    only = sys.argv[1:]
    if only and os.path.exists(OUT):  # update the named cases of an existing file (kernel part regenerated)
        ref["gradients"] = json.load(open(OUT))["gradients"]
    for case, fn, mol, pt2, fc in GRAD_CASES:
        if only and case not in only:
            continue
        t0 = time.time()
        ref["gradients"][case] = fd_gradient(case, fn, mol, pt2, fc, aux)
        print(f"{case:32s} E {ref['gradients'][case]['energy']:.10f} max|g| "
              f"{max(map(abs, ref['gradients'][case]['gradient'])):.3e} ({time.time() - t0:.0f} s)", flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, separators=(",", ":"))
    print("wrote", OUT, os.path.getsize(OUT), "bytes")


if __name__ == "__main__":
    main()
