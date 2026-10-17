#!/usr/bin/env python3
"""Reference analytic MP2 gradients for the native canonical MP2 gradient (test_cases/qm_mp2/reference_grad.json).

PySCF RHF (spherical def2-SVP from PySCF's own library, = curcuma's def2-SVP.dat; conv_tol 1e-12,
conv_tol_grad 1e-11, CPHF (Z-vector) tol 1e-13 -- a looser gradient criterion leaves 1e-8 Eh/Bohr of orbital error in the gradient)
followed by mp.MP2(mf).nuc_grad_method() (canonical, no density fitting, relaxed density via the
CPHF/Z-vector equations) with and without the chemical frozen core (frozen = number of 1s orbitals of
Li-Ne atoms). Stored per case: E_HF, E_MP2 (total) and the gradient in Eh/Bohr (atom-major, x y z) for
both settings. Run with the PySCF venv:
    ~/.venvs/qmref/bin/python scripts/mp2_gradient_reference.py [name=path.xyz ...]
Extra name=path pairs are printed only (not written to the json).

Measured (H2O, Oct 2026): PySCF's analytic MP2 gradient differs from the central finite difference of PySCF's
own MP2 energy by up to 1.1e-7 Eh/Bohr (cause not investigated; independent of the CPHF tolerance). Therefore
the script also stores "fd_gradient": the four-point central difference
    g = [E(-2h) - 8 E(-h) + 8 E(+h) - E(+2h)] / (12 h),  h = 1e-3 Bohr  (truncation O(h^4)),
of PySCF's MP2 energy (RHF conv_tol 1e-13); this is the tight reference of the ctest.
Claude Generated (Oct 2026, MP2 gradient).
"""
import json, os, sys
from pyscf import gto, scf, mp
from pyscf.grad import mp2 as grad_mp2
import pyscf

# PySCF's MP2 gradient solves the Z-vector (CPHF) with cphf.solve defaults (tol 1e-9, max_cycle 30); that left
# up to 1e-7 Eh/Bohr in the NH3 frozen-core gradient. Tighter here (tol 1e-13, 200 cycles).
_cphf_solve = grad_mp2.cphf.solve


def _tight_solve(*args, **kw):
    kw["tol"] = 1e-13
    kw["max_cycle"] = 200
    return _cphf_solve(*args, **kw)


grad_mp2.cphf.solve = _tight_solve

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "test_cases/qm_mp2/reference_grad.json")
CASES = [  # (name, xyz relative to test_cases/)
    ("H2O", "qm_1e/H2O.xyz"),
    ("NH3", "qm_1e/NH3.xyz"),
    ("CH4", "qm_1e/CH4.xyz"),
    ("HF", "qm_1e/HF.xyz"),
    ("H2O_distorted", "qm_grad/H2O_distorted.xyz"),
    ("formaldehyde_distorted", "qm_grad/formaldehyde_distorted.xyz"),
    ("formaldehyde", "qm_hf3c/formaldehyde.xyz"),
    ("benzene", "qm_hf3c/benzene.xyz"),
]


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    return [(a.split()[0], tuple(map(float, a.split()[1:4]))) for a in L[2:2 + n]]


def frozen_core(atoms):
    """Chemical frozen core as in qmint::frozenCoreOrbitals (H-He 0, Li-Ne 1)."""
    return sum(0 if gto.charge(s) <= 2 else 1 for s, _ in atoms)


FD_H = 1e-3
BOHR = 0.529177210903


def mp2_energy(sym, pos_bohr, frz):
    mol = gto.M(atom=[[s, p] for s, p in zip(sym, pos_bohr)], basis="def2-svp", unit="Bohr", cart=False, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-11
    mf.max_cycle = 200
    mf.kernel()
    p = mp.MP2(mf, frozen=frz)
    p.kernel()
    return p.e_tot


def fd_gradient(atoms, frz):
    """Four-point central difference of the PySCF MP2 energy (Eh/Bohr, atom-major)."""
    sym = [s for s, _ in atoms]
    x0 = [[c / BOHR for c in x] for _, x in atoms]
    g = []
    for i in range(len(sym)):
        for k in range(3):
            e = {}
            for m in (-2, -1, 1, 2):
                x = [r[:] for r in x0]
                x[i][k] += m * FD_H
                e[m] = mp2_energy(sym, x, frz)
            g.append((e[-2] - 8 * e[-1] + 8 * e[1] - e[2]) / (12 * FD_H))
    return g


def main():
    extra = sys.argv[1:]
    cases = CASES if not extra else [tuple(x.split("=", 1)) for x in extra]
    out = {"_about": "PySCF canonical RMP2/def2-SVP analytic gradients (Eh/Bohr) for ctest -L qm_mp2 "
                     "(scripts/mp2_gradient_reference.py)",
           "_sources": {"pyscf": pyscf.__version__, "basis": "def2-svp (PySCF library), spherical",
                        "scf": "RHF conv_tol 1e-12 conv_tol_grad 1e-11",
                        "mp2": "mp.MP2(mf).nuc_grad_method(), frozen=None and frozen=n_core, CPHF tol 1e-13"},
           "_units": "Hartree, Hartree/Bohr", "cases": {}}
    for name, rel in cases:
        path = rel if os.path.isabs(rel) else os.path.join(ROOT, "test_cases", rel)
        atoms = read_xyz(path)
        mol = gto.M(atom=[[s, x] for s, x in atoms], basis="def2-svp", unit="Angstrom", cart=False, verbose=0)
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-13
        mf.conv_tol_grad = 1e-11
        mf.max_cycle = 200
        e_hf = mf.kernel()
        if not mf.converged:
            sys.exit(f"{name}: PySCF RHF did not converge")
        entry = {"xyz": rel, "nbf": mol.nao_nr(), "hf": e_hf}
        nfc = frozen_core(atoms)
        for tag, frz in (("full", None), ("fc", nfc if nfc > 0 else None)):
            p = mp.MP2(mf, frozen=frz)
            p.kernel()
            g = p.nuc_grad_method().kernel()
            fd = fd_gradient(atoms, frz)
            entry[tag] = {"frozen": frz or 0, "e_tot": p.e_tot, "e_corr": p.e_corr,
                          "gradient": [float(x) for x in g.reshape(-1)], "fd_gradient": fd,
                          "fd_step_bohr": FD_H}
            print(f"{name:24s} {tag:4s} max|analytic - FD| (PySCF) = {max(abs(a - b) for a, b in zip(g.reshape(-1), fd)):.2e}")
            print(f"{name:24s} {tag:4s} E={p.e_tot:.12f} max|g|={abs(g).max():.6e} sum={g.sum(axis=0)}")
        out["cases"][name] = entry
    if not extra:
        if os.path.exists(OUT) and "opt" in json.load(open(OUT)):
            out["opt"] = json.load(open(OUT))["opt"]  # written by mp2_opt_reference.py, kept
        json.dump(out, open(OUT, "w"), indent=1)
        print("wrote", OUT)


if __name__ == "__main__":
    main()
