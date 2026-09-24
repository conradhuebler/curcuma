#!/usr/bin/env python3
"""Reference gradients for the native HF / HF-3c gradient (test_cases/qm_grad/reference.json).

  hf     : PySCF analytic RHF gradient (mf.nuc_grad_method()), basis def2-SVP from
           PySCF's own library (spherical 5d, = curcuma's def2-SVP.dat: benzene
           HF/def2-SVP agrees to 8 decimals) or MINIX parsed from MINIX.dat
  hf-3c  : PySCF RHF/MINIX gradient + simple-dftd3 D3(BJ, method="hf3c") gradient
           + simple-dftd3 GeometricCounterpoise(method="hf3c") gradient
All in Eh/Bohr, atom-major (x0,y0,z0,x1,...). Requirements: pip install pyscf dftd3.
Claude Generated (Sep 2026).
"""
import json, os, sys
import numpy as np
from pyscf import gto, scf
import pyscf
from dftd3.interface import GeometricCounterpoise, RationalDampingParam, DispersionModel
from dftd3.library import get_api_version

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
from hf3c_reference import load_gamess_basis, read_xyz, MINIX, AATOAU  # noqa: E402

OUT = os.path.join(ROOT, "test_cases/qm_grad/reference.json")
CASES = [  # (name, xyz relative to test_cases/, method)
    ("H2O_svp", "qm_1e/H2O.xyz", "hf"),
    ("NH3_svp", "qm_1e/NH3.xyz", "hf"),
    ("CH4_svp", "qm_1e/CH4.xyz", "hf"),
    ("HF_svp", "qm_1e/HF.xyz", "hf"),
    ("LiH_svp", "qm_1e/LiH.xyz", "hf"),
    ("H2O_distorted_svp", "qm_grad/H2O_distorted.xyz", "hf"),
    ("formaldehyde_distorted_svp", "qm_grad/formaldehyde_distorted.xyz", "hf"),
    ("HCN_svp", "qm_hf3c/HCN.xyz", "hf"),
    ("H2O_hf3c", "qm_1e/H2O.xyz", "hf-3c"),
    ("NH3_hf3c", "qm_1e/NH3.xyz", "hf-3c"),
    ("H2O_distorted_hf3c", "qm_grad/H2O_distorted.xyz", "hf-3c"),
    ("formaldehyde_distorted_hf3c", "qm_grad/formaldehyde_distorted.xyz", "hf-3c"),
    ("water_dimer_hf3c", "qm_hf3c/water_dimer.xyz", "hf-3c"),
    ("BF3_hf3c", "qm_hf3c/BF3.xyz", "hf-3c"),
    ("LiF_hf3c", "qm_hf3c/LiF.xyz", "hf-3c"),
    ("He_CH4_hf3c", "qm_hf3c/He_CH4.xyz", "hf-3c"),
    ("benzene_hf3c", "qm_hf3c/benzene.xyz", "hf-3c"),
]


def main():
    minix = load_gamess_basis(MINIX)
    out = {"_about": "Reference gradients (Eh/Bohr) for ctest -L qm_grad",
           "_sources": {"pyscf": pyscf.__version__, "simple-dftd3": get_api_version()},
           "cases": {}}
    for name, rel, method in CASES:
        atoms = read_xyz(os.path.join(ROOT, "test_cases", rel))
        basis = {s: minix[s] for s, _ in atoms} if method == "hf-3c" else "def2-svp"
        mol = gto.M(atom=[[s, x] for s, x in atoms], basis=basis, unit='Angstrom',
                    cart=False, verbose=0)
        mf = scf.RHF(mol); mf.conv_tol = 1e-12; mf.conv_tol_grad = 1e-9; mf.max_cycle = 200
        e = mf.kernel()
        if not mf.converged:
            sys.exit(f"{name}: PySCF RHF did not converge")
        g = mf.nuc_grad_method().kernel()
        entry = {"xyz": rel, "method": method, "hf_energy": e, "hf_gradient": g.ravel().tolist()}
        if method == "hf-3c":
            num = np.array([gto.charge(s) for s, _ in atoms])
            pos = np.array([x for _, x in atoms]) * AATOAU
            d3 = DispersionModel(num, pos).get_dispersion(RationalDampingParam(method="hf3c"), grad=True)
            gc = GeometricCounterpoise(num, pos, method="hf3c").get_counterpoise(grad=True)
            entry["d3_gradient"] = np.asarray(d3["gradient"]).ravel().tolist()
            entry["gcp_gradient"] = np.asarray(gc["gradient"]).ravel().tolist()
            g = g + np.asarray(d3["gradient"]) + np.asarray(gc["gradient"])
            e = e + float(d3["energy"]) + float(gc["energy"])
        entry["energy"] = e
        entry["gradient"] = g.ravel().tolist()
        out["cases"][name] = entry
        print(f"{name:28s} E={e:.10f} max|g|={np.abs(g).max():.4e}")
    json.dump(out, open(OUT, "w"), indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
