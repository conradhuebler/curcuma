#!/usr/bin/env python3
"""DFT-D4 reference values from the dftd4 Python API (WP12 stage 4), written to test_cases/qm_dh/reference_d4.json.

For every molecule x functional: D4(EEQ)-ATM energy (model "d4", dftd4 default cutoffs), its two-body and ATM parts
(get_pairwise_dispersion: additive / non-additive pairwise energies summed), EEQ charges and CN (get_properties) and
the analytic gradient (Eh/Bohr). Parameters: dftd4's own table (DampingParam(method=...), bj-eeq-atm) for pbe,
b3lyp, tpss, scan, b2plyp, b2gpplyp; explicit s6/s8/a1/a2/s9 for revdsd-pbep86 (Mehta et al. 2023 Table 1, PZ81 row;
dftd4's own revdsdpbep86 entry is the 2019 set and is NOT used).

Geometries: Angstrom from the xyz files, converted with 0.529177210903 Angstrom/Bohr (CODATA 2018, curcuma's
CurcumaUnit::Length::BOHR_TO_ANGSTROM).

Environment: dftd4 4.3.1 from PyPI in a venv (not the system python), e.g.
    python3 -m venv <dir>/d4venv && <dir>/d4venv/bin/pip install dftd4==4.3.1 tomli
    <dir>/d4venv/bin/python scripts/qm_d4_dftd4_reference.py
Claude Generated (Oct 2026).
"""
import json
import os
import sys

import numpy as np
import dftd4
from dftd4.interface import DampingParam, DispersionModel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TC = os.path.join(ROOT, "test_cases")
BOHR = 0.529177210903
SYM = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10}
MOLECULES = {
    "H2O": "qm_1e/H2O.xyz",
    "NH3": "qm_1e/NH3.xyz",
    "CH4": "qm_1e/CH4.xyz",
    "H2CO": "qm_grad/formaldehyde_hf3c_min_reference.xyz",
    "benzene": "qm_hf3c/benzene.xyz",
}
# curcuma key -> dftd4 parameters
FUNCTIONALS = {
    "pbe": {"method": "pbe"},
    "b3lyp": {"method": "b3lyp"},
    "tpss": {"method": "tpss"},
    "scan": {"method": "scan"},
    "b2plyp": {"method": "b2plyp"},
    "b2gp-plyp": {"method": "b2gpplyp"},
    "revdsd-pbep86": {"s6": 0.5917, "s8": 0.0, "a1": 0.371, "a2": 4.201, "s9": 1.0},
}


def read_xyz(path):
    lines = open(path).read().split("\n")
    n = int(lines[0])
    z, x = [], []
    for line in lines[2:2 + n]:
        t = line.split()
        z.append(SYM[t[0]])
        x.append([float(v) for v in t[1:4]])
    return np.array(z), np.array(x) / BOHR


def main():
    out = {"_about": (f"dftd4 {dftd4.library.get_api_version()} Python API, model d4, default cutoffs, "
                      "bj-eeq-atm; scripts/qm_d4_dftd4_reference.py"), "cases": {}}
    for mol, rel in MOLECULES.items():
        z, x = read_xyz(os.path.join(TC, rel))
        model = DispersionModel(z, x)
        props = model.get_properties()
        for fn, spec in FUNCTIONALS.items():
            param = DampingParam(**spec)
            r = model.get_dispersion(param, grad=True)
            pw = model.get_pairwise_dispersion(param)
            out["cases"][f"{mol}_{fn}"] = {
                "molecule": mol, "xyz": rel, "functional": fn, "dftd4_param": spec,
                "energy": float(r["energy"]),
                "two_body": float(np.sum(pw["additive pairwise energy"])),
                "three_body": float(np.sum(pw["non-additive pairwise energy"])),
                "gradient": [float(v) for v in np.asarray(r["gradient"]).reshape(-1)],
                "charges": [float(v) for v in props["partial charges"]],
                "cn": [float(v) for v in props["coordination numbers"]],
            }
    path = os.path.join(TC, "qm_dh", "reference_d4.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", path, len(out["cases"]), "cases")


if __name__ == "__main__":
    sys.exit(main())
