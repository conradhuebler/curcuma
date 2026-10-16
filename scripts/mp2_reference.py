#!/usr/bin/env python3
"""Reference MP2 energies for the native canonical MP2 (test_cases/qm_mp2/reference.json).

PySCF RHF (spherical def2-SVP from PySCF's own library, = curcuma's def2-SVP.dat; conv_tol 1e-12)
followed by mp.MP2(mf) (canonical, no density fitting) with and without the chemical frozen core
(frozen = number of 1s orbitals of Li-Ne atoms). Stored per case: E_HF, E_corr, E_OS, E_SS for
both settings (Hartree). Run with the PySCF venv:
    ~/.venvs/qmref/bin/python scripts/mp2_reference.py
Claude Generated (Oct 2026, WP11).
"""
import json, os, sys
from pyscf import gto, scf, mp
import pyscf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    return [(a.split()[0], tuple(map(float, a.split()[1:4]))) for a in L[2:2 + n]]

OUT = os.path.join(ROOT, "test_cases/qm_mp2/reference.json")
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


def frozen_core(atoms):
    """Chemical frozen core as in qmint::frozenCoreOrbitals (H-He 0, Li-Ne 1)."""
    return sum(0 if gto.charge(s) <= 2 else 1 for s, _ in atoms)


def main():
    extra = sys.argv[1:]  # optional "name=path.xyz" pairs (not written to the json)
    cases = CASES if not extra else [tuple(x.split("=", 1)) for x in extra]
    out = {"_about": "PySCF canonical RMP2/def2-SVP references for ctest -L qm_mp2 (scripts/mp2_reference.py)",
           "_sources": {"pyscf": pyscf.__version__, "basis": "def2-svp (PySCF library), spherical",
                        "scf": "RHF conv_tol 1e-12", "mp2": "mp.MP2(mf), frozen=None and frozen=n_core"},
           "_units": "Hartree", "cases": {}}
    for name, rel in cases:
        path = rel if os.path.isabs(rel) else os.path.join(ROOT, "test_cases", rel)
        atoms = read_xyz(path)
        mol = gto.M(atom=[[s, x] for s, x in atoms], basis="def2-svp", unit="Angstrom", cart=False, verbose=0)
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.max_cycle = 200
        e_hf = mf.kernel()
        if not mf.converged:
            sys.exit(f"{name}: PySCF RHF did not converge")
        entry = {"xyz": rel, "nbf": mol.nao_nr(), "hf": e_hf}
        nfc = frozen_core(atoms)
        for tag, frz in (("full", None), ("fc", nfc if nfc > 0 else None)):
            p = mp.MP2(mf, frozen=frz)
            p.kernel()
            entry[tag] = {"frozen": frz or 0, "e_corr": p.e_corr, "e_os": p.e_corr_os, "e_ss": p.e_corr_ss}
        out["cases"][name] = entry
        print(f"{name:24s} nbf={entry['nbf']:4d} E_HF={e_hf:.10f} Ecorr={entry['full']['e_corr']:.12f} "
              f"Ecorr(fc={nfc})={entry['fc']['e_corr']:.12f}")
    if not extra:
        json.dump(out, open(OUT, "w"), indent=1)
        print("wrote", OUT)


if __name__ == "__main__":
    main()
