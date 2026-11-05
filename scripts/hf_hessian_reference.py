#!/usr/bin/env python3
"""Reference Hessians for the native analytic RHF Hessian (test_cases/qm_hessian/reference.json).

PySCF analytic RHF Hessian (pyscf.hessian.rhf, mf.Hessian()), basis from PySCF's own library (spherical, the same
def2 sets as curcuma's built-in ones), tight SCF and CPHF. Hessian in Eh/Bohr^2, row/column index 3*atom + xyz;
harmonic frequencies (cm^-1) from pyscf.hessian.thermo.harmonic_analysis (translations and rotations projected).
Also stored: "hessian_fd", the 4-point central finite difference (h = 5e-4 Bohr) of PySCF's analytic gradient. It is
the tighter reference: PySCF's analytic Hessian differs from its own FD by up to 9.0e-8 Eh/Bohr^2 (distorted
formaldehyde, def2-SVP; not changed by conv_tol_cpscf 1e-12 or direct_scf_tol 1e-16; cause not investigated).

  hf_hessian_reference.py                      write reference.json (def2-SVP cases of ctest -L qm_hess)
  hf_hessian_reference.py --xyz F --basis B    print the Hessian and the frequencies of one structure as JSON

Requirements: pip install pyscf. Claude Generated (Nov 2026).
"""
import argparse, json, os, sys
import numpy as np
import pyscf
from pyscf import gto, scf
from pyscf.hessian import thermo

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "test_cases/qm_hessian/reference.json")
CASES = [  # (name, xyz relative to test_cases/, basis)
    ("H2O_svp", "qm_1e/H2O.xyz", "def2-svp"),
    ("NH3_svp", "qm_1e/NH3.xyz", "def2-svp"),
    ("CH4_svp", "qm_1e/CH4.xyz", "def2-svp"),
    ("HF_svp", "qm_1e/HF.xyz", "def2-svp"),
    ("H2O_distorted_svp", "qm_grad/H2O_distorted.xyz", "def2-svp"),
    ("formaldehyde_distorted_svp", "qm_grad/formaldehyde_distorted.xyz", "def2-svp"),
]


def read_xyz(path):
    lines = open(path).read().split("\n")
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        p = line.split()
        atoms.append([p[0], [float(p[1]), float(p[2]), float(p[3])]])
    return atoms


def scf_of(atoms, basis, unit="Angstrom"):
    mol = gto.M(atom=atoms, basis=basis, unit=unit, cart=False, verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 300
    mf.conv_tol_cpscf = 1e-12  # CPHF tolerance of pyscf.hessian.rhf (default 1e-9)
    e = mf.kernel()
    if not mf.converged:
        sys.exit("PySCF RHF did not converge")
    return mol, mf, e


def fd_hessian(atoms, basis, h=5e-4):
    """4-point central difference of PySCF's analytic gradient, displacements in Bohr."""
    bohr = 0.529177210903
    c0 = np.array([x for _, x in atoms]) / bohr
    n = len(atoms)
    H = np.zeros((3 * n, 3 * n))
    for i in range(n):
        for k in range(3):
            g = []
            for s in (2, 1, -1, -2):
                c = c0.copy()
                c[i, k] += s * h
                _, mf, _ = scf_of([[a[0], list(x)] for a, x in zip(atoms, c)], basis, unit="Bohr")
                g.append(mf.nuc_grad_method().kernel().ravel())
            H[3 * i + k] = (-g[0] + 8 * g[1] - 8 * g[2] + g[3]) / (12 * h)
    return H


def run(atoms, basis):
    mol, mf, e = scf_of(atoms, basis)
    hobj = mf.Hessian()
    hobj.max_cycle = 200
    h = hobj.kernel()  # (natm, natm, 3, 3)
    n = mol.natm
    H = h.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
    freq = thermo.harmonic_analysis(mol, h)["freq_wavenumber"]
    freq = [float(np.real(f)) if abs(np.imag(f)) < 1e-12 else -float(np.abs(np.imag(f))) for f in freq]
    g = mf.nuc_grad_method().kernel()
    return e, g, H, freq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz")
    ap.add_argument("--basis", default="def2-svp")
    a = ap.parse_args()
    if a.xyz:
        e, g, H, freq = run(read_xyz(a.xyz), a.basis)
        print(json.dumps({"energy": e, "gradient": g.ravel().tolist(), "hessian": H.ravel().tolist(),
                          "freq_wavenumber": freq, "pyscf": pyscf.__version__}))
        return
    out = {"_about": "PySCF analytic RHF Hessians (Eh/Bohr^2, index 3*atom+xyz) for ctest -L qm_hess",
           "_sources": {"pyscf": pyscf.__version__, "settings": "conv_tol 1e-13, conv_tol_grad 1e-10, conv_tol_cpscf 1e-12; hessian_fd: 4-point FD of the analytic gradient, h = 5e-4 Bohr"},
           "cases": {}}
    for name, rel, basis in CASES:
        e, g, H, freq = run(read_xyz(os.path.join(ROOT, "test_cases", rel)), basis)
        atoms = read_xyz(os.path.join(ROOT, "test_cases", rel))
        Hfd = fd_hessian(atoms, basis)
        out["cases"][name] = {"xyz": rel, "basis": basis, "energy": e, "gradient": g.ravel().tolist(),
                              "hessian": H.ravel().tolist(), "hessian_fd": Hfd.ravel().tolist(),
                              "freq_wavenumber": freq}
        print(f"{name:28s} E={e:.10f} max|H|={np.abs(H).max():.4e} max|H - FD(PySCF grad)|={np.abs(H - Hfd).max():.2e}")
    json.dump(out, open(OUT, "w"), indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
