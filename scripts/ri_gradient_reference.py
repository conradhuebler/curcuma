#!/usr/bin/env python3
"""Reference analytic gradients of the density-fitted RHF (RI-J, RI-JK) for test_cases/qm_ri/reference_grad.json.

PySCF RHF/def2-SVP (PySCF library basis = curcuma's def2-SVP.dat), with the SAME auxiliary
basis data as curcuma: src/core/energy_calculators/qm_methods/def2-universal-{JFIT,JKFIT}.turbomole are parsed
here (parse_turbomole from scripts/rimp2_reference.py) and handed to PySCF as a basis dict.
  jk: mf = scf.RHF(mol).density_fit(auxbasis=JKFIT)              (Coulomb and exchange fitted)
  j : mf = scf.RHF(mol).density_fit(auxbasis=JFIT, only_dfj=True) (Coulomb fitted, exchange exact = ORCA RIJONX)
conv_tol 1e-14, conv_tol_grad 1e-11 (with conv_tol_grad 1e-8 the PySCF gradient moved by ~2e-9 Eh/Bohr).
Gradient: mf.nuc_grad_method().kernel() (pyscf.df.grad.rhf, auxbasis_response=True: the metric and the
auxiliary-centre derivatives are included, i.e. the derivative of the RI energy). Units Eh/Bohr, atom-major.
Run with the PySCF venv:
    ~/.venvs/qmref/bin/python scripts/ri_gradient_reference.py              # writes the json
    ~/.venvs/qmref/bin/python scripts/ri_gradient_reference.py name=x.xyz   # print only (absolute path or rel. to test_cases)
Claude Generated (Oct 2026, RI gradient).
"""
import json, os, sys
import numpy as np
import pyscf
from pyscf import gto, scf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from rimp2_reference import parse_turbomole, read_xyz, ROOT, QMDIR  # noqa: E402

OUT = os.path.join(ROOT, "test_cases/qm_ri/reference_grad.json")
CASES = [  # (name, xyz relative to test_cases/, cartesian)
    ("H2O", "qm_grad/H2O_distorted.xyz", False),
    ("NH3", "qm_1e/NH3.xyz", False),
    ("formaldehyde", "qm_grad/formaldehyde_distorted.xyz", False),
    ("benzene", "qm_hf3c/benzene.xyz", False),
    ("H2O_cartesian", "qm_grad/H2O_distorted.xyz", True),
]


def run(path, cart, aux):
    atoms = read_xyz(path)
    mol = gto.M(atom=[[s, x] for s, x in atoms], basis="def2-svp", unit="Angstrom", cart=cart, verbose=0)
    res = {}
    for mode in ("j", "jk"):
        mf = scf.RHF(mol)
        if mode == "jk":
            mf = mf.density_fit(auxbasis=aux["JKFIT"])
        else:
            mf = mf.density_fit(auxbasis=aux["JFIT"], only_dfj=True)
        mf.conv_tol = 1e-14
        mf.conv_tol_grad = 1e-11  # gradient error is first order in the orbital gradient
        mf.max_cycle = 200
        e = mf.kernel()
        if not mf.converged:
            sys.exit(f"{path} {mode}: PySCF SCF did not converge")
        gobj = mf.nuc_grad_method()
        assert getattr(gobj, "auxbasis_response", True)
        g = gobj.kernel()
        res[mode] = {"energy": e, "gradient": [float(x) for x in np.asarray(g).ravel()],
                     "max_abs": float(np.abs(g).max()), "sum": [float(x) for x in np.asarray(g).sum(axis=0)]}
    return res


def main():
    aux = {"JKFIT": parse_turbomole(os.path.join(QMDIR, "def2-universal-JKFIT.turbomole")),
           "JFIT": parse_turbomole(os.path.join(QMDIR, "def2-universal-JFIT.turbomole"))}
    extra = sys.argv[1:]
    cases = CASES if not extra else [(x.split("=", 1)[0], x.split("=", 1)[1], False) for x in extra]
    out = {"_comment": ["RI-J / RI-JK HF/def2-SVP analytic gradients (Eh/Bohr, atom-major) from PySCF, produced by",
                        "  ~/.venvs/qmref/bin/python scripts/ri_gradient_reference.py",
                        "  (same auxiliary basis files as curcuma; jk: density_fit(JKFIT); j: density_fit(JFIT, only_dfj=True);",
                        "  mf.nuc_grad_method() with auxbasis_response; conv_tol 1e-14, conv_tol_grad 1e-11).",
                        f"PySCF {pyscf.__version__}. Claude Generated (Oct 2026)."],
           "cases": {}}
    for name, rel, cart in cases:
        path = rel if os.path.isabs(rel) else os.path.join(ROOT, "test_cases", rel)
        r = run(path, cart, aux)
        for m in ("j", "jk"):
            print(f"{name:16s} {m:2s} E={r[m]['energy']:.10f} max|g|={r[m]['max_abs']:.6e} sum={['%.1e' % s for s in r[m]['sum']]}")
        out["cases"][name] = {"xyz": rel, "cartesian": cart, **r}
    if not extra:
        if os.path.exists(OUT) and "opt" in json.load(open(OUT)):
            out["opt"] = json.load(open(OUT))["opt"]  # kept: written by scripts/ri_opt_reference.py
        json.dump(out, open(OUT, "w"), indent=1)
        print("wrote", OUT)
    else:
        print(json.dumps(out["cases"]))


if __name__ == "__main__":
    main()
