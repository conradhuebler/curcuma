#!/usr/bin/env python3
"""Reference RI-JK HF/def2-SVP minima for the -opt test with -qm.ri jk (test_cases/qm_ri/reference_grad.json, key "opt").

PySCF RHF/def2-SVP with mf.density_fit(auxbasis=<curcuma's def2-universal-JKFIT.turbomole>), analytic DF gradient
(auxbasis_response), minimised with scipy BFGS (gtol 1e-7 Eh/Bohr) in Cartesian Bohr coordinates.
geomeTRIC/pyberny are not installed in the venv, hence scipy, as scripts/hf3c_opt_reference.py.
Usage: ~/.venvs/qmref/bin/python scripts/ri_opt_reference.py      (adds/updates "opt" in the json)
Claude Generated (Oct 2026, RI gradient).
"""
import json, os, sys
import numpy as np
from pyscf import gto, scf
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rimp2_reference import parse_turbomole, read_xyz, ROOT, QMDIR  # noqa: E402

BOHR = 0.529177210903  # curcuma's CODATA 2018 value (src/core/units.h)
OUT = os.path.join(ROOT, "test_cases/qm_ri/reference_grad.json")
CASES = [("H2O", "qm_grad/H2O_distorted.xyz"), ("formaldehyde", "qm_grad/formaldehyde_distorted.xyz")]


def minimum(path, aux):
    atoms = read_xyz(path)
    sym = [s for s, _ in atoms]
    x0 = np.array([x for _, x in atoms]).ravel() / BOHR
    dm = [None]

    def f(x):
        mol = gto.M(atom=[[s, p] for s, p in zip(sym, x.reshape(-1, 3))], basis="def2-svp", unit="Bohr", verbose=0)
        mf = scf.RHF(mol).density_fit(auxbasis=aux)
        mf.conv_tol = 1e-13
        mf.conv_tol_grad = 1e-10
        e = mf.kernel(dm0=dm[0])
        dm[0] = mf.make_rdm1()
        return e, mf.nuc_grad_method().kernel().ravel()

    r = minimize(f, x0, jac=True, method="BFGS", options={"gtol": 1e-7, "maxiter": 500})
    print(f"{path}: E_min {r.fun:.10f}  |g|max {np.abs(r.jac).max():.2e}  iterations {r.nit}")
    return {"start": os.path.relpath(path, os.path.join(ROOT, "test_cases")), "energy": float(r.fun),
            "max_abs_gradient": float(np.abs(r.jac).max()), "symbols": sym,
            "geometry_angstrom": [float(v) for v in (r.x * BOHR)]}


def main():
    aux = parse_turbomole(os.path.join(QMDIR, "def2-universal-JKFIT.turbomole"))
    ref = json.load(open(OUT))
    ref["opt"] = {"_comment": "RI-JK HF/def2-SVP minima (PySCF DF + scipy BFGS gtol 1e-7), scripts/ri_opt_reference.py"}
    for name, rel in CASES:
        ref["opt"][name] = minimum(os.path.join(ROOT, "test_cases", rel), aux)
    json.dump(ref, open(OUT, "w"), indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
