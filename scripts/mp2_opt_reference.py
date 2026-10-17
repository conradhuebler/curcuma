#!/usr/bin/env python3
"""Reference MP2/def2-SVP minima for the -opt -method mp2 tests (test_cases/qm_mp2/reference_grad.json, key "opt").

PySCF RHF (conv_tol 1e-12, conv_tol_grad 1e-11) + mp.MP2(mf).nuc_grad_method() (full MP2, def2-SVP spherical),
CPHF tol 1e-13, minimised with scipy BFGS (gtol 1e-7 Eh/Bohr) from the distorted start structures. Stored per case: the start
xyz, the minimum energy and all interatomic distances of the minimum (Angstrom, pairs i<j in atom order), so
that the test needs no structure file and no numpy. Run with the PySCF venv after mp2_gradient_reference.py:
    ~/.venvs/qmref/bin/python scripts/mp2_opt_reference.py
Claude Generated (Oct 2026, MP2 gradient).
"""
import json, os, sys
import numpy as np
from pyscf import gto, scf, mp
from pyscf.grad import mp2 as grad_mp2
from scipy.optimize import minimize

_cphf_solve = grad_mp2.cphf.solve


def _tight_solve(*args, **kw):  # as in mp2_gradient_reference.py: tight Z-vector (CPHF) solution
    kw["tol"] = 1e-13
    kw["max_cycle"] = 200
    return _cphf_solve(*args, **kw)


grad_mp2.cphf.solve = _tight_solve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "test_cases/qm_mp2/reference_grad.json")
BOHR = 0.529177210903
CASES = [("H2O", "qm_grad/H2O_distorted.xyz"), ("formaldehyde", "qm_grad/formaldehyde_distorted.xyz")]


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    return [(a.split()[0], tuple(map(float, a.split()[1:4]))) for a in L[2:2 + n]]


def main():
    data = json.load(open(OUT))
    data["opt"] = {"_about": "PySCF MP2/def2-SVP minima (scipy BFGS gtol 1e-7 Eh/Bohr, scripts/mp2_opt_reference.py); "
                             "distances in Angstrom for the atom pairs i<j", "cases": {}}
    for name, rel in CASES:
        atoms = read_xyz(os.path.join(ROOT, "test_cases", rel))
        sym = [s for s, _ in atoms]
        x0 = np.array([x for _, x in atoms]).ravel() / BOHR
        dm = [None]

        def f(x):
            mol = gto.M(atom=[[s, p] for s, p in zip(sym, x.reshape(-1, 3))], basis="def2-svp", unit="Bohr",
                        cart=False, verbose=0)
            mf = scf.RHF(mol)
            mf.conv_tol = 1e-12
            mf.conv_tol_grad = 1e-11
            mf.kernel(dm0=dm[0])
            dm[0] = mf.make_rdm1()
            p = mp.MP2(mf)
            p.kernel()
            g = p.nuc_grad_method().kernel()
            return p.e_tot, g.ravel()

        r = minimize(f, x0, jac=True, method="BFGS", options={"gtol": 1e-7, "maxiter": 500})
        pos = r.x.reshape(-1, 3) * BOHR
        dist = [float(np.linalg.norm(pos[i] - pos[j])) for i in range(len(sym)) for j in range(i + 1, len(sym))]
        data["opt"]["cases"][name] = {"start": rel, "e_min": float(r.fun), "distances": dist,
                                      "max_grad": float(np.abs(r.jac).max())}
        print(f"{name}: E_min {r.fun:.12f} |g|max {np.abs(r.jac).max():.2e} iters {r.nit} distances {dist}")
    json.dump(data, open(OUT, "w"), indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
