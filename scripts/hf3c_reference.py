#!/usr/bin/env python3
"""Reference data for the native HF-3c (test_cases/qm_hf3c/reference.json).

HF-3c = HF/MINIX + D3(BJ) + gCP + SRB (Sure & Grimme, J. Comput. Chem. 34, 1672
(2013)). Every term comes from an independent program, none shares code with
curcuma:

  HF/MINIX      PySCF (RHF, spherical, conv_tol 1e-12) with the basis parsed from
                src/core/energy_calculators/qm_methods/MINIX.dat (ORCA export)
  D3(BJ)        simple-dftd3 Python API, RationalDampingParam(method="hf3c"),
                i.e. s6=1, s8=0.8777, a1=0.4171, a2=2.9149, s9=0 (two-body only)
  gCP + SRB     simple-dftd3 GeometricCounterpoise(method="hf3c"); the gCP part
                alone is GeometricCounterpoise(method="hf", basis="minix")

Cross-check against ORCA 6.1 (`! HF-3c`, H2O, from docs/NATIVE_QM_IMPLEMENTATION.md):
this chain gives -75.5014894618 Eh vs ORCA -75.501489461360 (4e-10).

Requirements:  pip install pyscf dftd3
Usage:         hf3c_reference.py            (rewrites reference.json)
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
MINIX = os.path.join(ROOT, "src/core/energy_calculators/qm_methods/MINIX.dat")
OUT = os.path.join(ROOT, "test_cases/qm_hf3c/reference.json")
AATOAU = 1.0 / 0.529177210903

SYM = {'HYDROGEN': 'H', 'HELIUM': 'He', 'LITHIUM': 'Li', 'BERYLLIUM': 'Be', 'BORON': 'B',
       'CARBON': 'C', 'NITROGEN': 'N', 'OXYGEN': 'O', 'FLUORINE': 'F', 'NEON': 'Ne'}

# name -> xyz path relative to test_cases/
MOLECULES = {m: f"qm_1e/{m}.xyz" for m in
             ["H2", "He", "LiH", "BeH2", "BH", "CH4", "NH3", "H2O", "HF", "Ne"]}
MOLECULES.update({m: f"qm_hf3c/{m}.xyz" for m in
                  ["benzene", "water_dimer", "formaldehyde", "HCN", "BF3", "LiF", "He_CH4"]})


def load_gamess_basis(path):
    """Parse the GAMESS-style $DATA block ORCA exports (one element per block)."""
    basis, cur = {}, None
    lines = open(path).read().split('$DATA')[1].splitlines()
    i = 0
    while i < len(lines):
        l = lines[i].strip(); i += 1
        if not l or l.startswith('!') or l.startswith('$'):
            continue
        if l.upper() in SYM:
            cur = SYM[l.upper()]; basis[cur] = []; continue
        t = l.split(); L = {'S': 0, 'P': 1, 'D': 2}[t[0]]; n = int(t[1]); prims = []
        for _ in range(n):
            p = lines[i].split(); i += 1
            prims.append([float(p[1]), float(p[2])])
        basis[cur].append([L] + prims)
    return basis


def read_xyz(path):
    L = open(path).read().splitlines(); n = int(L[0])
    return [(a.split()[0], tuple(map(float, a.split()[1:4]))) for a in L[2:2 + n]]


def main():
    B = load_gamess_basis(MINIX)
    out = {"_about": __doc__.strip().splitlines()[0],
           "_sources": {"pyscf": pyscf.__version__, "simple-dftd3": get_api_version(),
                        "basis": "MINIX.dat (ORCA 6.1 orca_exportbasis -b minix)"},
           "_units": "Hartree",
           "orca_h2o": {"hf": -75.48881664772478, "d3": -0.002646402235,
                        "gcp_srb": -0.010026411, "total": -75.501489461360,
                        "note": "ORCA 6.1 ! HF-3c, printed decomposition"},
           "molecules": {}}
    for name, rel in MOLECULES.items():
        atoms = read_xyz(os.path.join(ROOT, "test_cases", rel))
        mol = gto.M(atom=[[s, x] for s, x in atoms], basis={s: B[s] for s, _ in atoms},
                    unit='Angstrom', cart=False, verbose=0)
        mf = scf.RHF(mol); mf.conv_tol = 1e-12; mf.max_cycle = 200
        e_hf = mf.kernel()
        if not mf.converged:
            sys.exit(f"{name}: PySCF RHF did not converge")
        num = np.array([gto.charge(s) for s, _ in atoms])
        pos = np.array([x for _, x in atoms]) * AATOAU   # same conversion as curcuma
        e_d3 = float(DispersionModel(num, pos).get_dispersion(
            RationalDampingParam(method="hf3c"), grad=False)['energy'])
        e_tot_gcp = float(GeometricCounterpoise(num, pos, method="hf3c")
                          .get_counterpoise(grad=False)['energy'])
        e_gcp = float(GeometricCounterpoise(num, pos, method="hf", basis="minix")
                      .get_counterpoise(grad=False)['energy'])
        out["molecules"][name] = {
            "xyz": rel, "nbf": int(mol.nao), "hf": e_hf, "d3": e_d3,
            "gcp": e_gcp, "srb": e_tot_gcp - e_gcp, "total": e_hf + e_d3 + e_tot_gcp}
        print(f"{name:13s} nbf={mol.nao:3d} HF={e_hf:.10f} D3={e_d3:.10f} "
              f"gCP={e_gcp:.10f} SRB={e_tot_gcp - e_gcp:.10f} total={e_hf + e_d3 + e_tot_gcp:.10f}")
    json.dump(out, open(OUT, "w"), indent=1)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
