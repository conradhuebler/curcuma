#!/usr/bin/env python3
"""ORCA 6.1 reference energies for the built-in basis sets (test_cases/qm_basis_embed/reference_orca.json).

For H2O and NH3 (test_cases/qm_1e/*.xyz) every case below is run with ORCA ('! RHF <basis> NoRI|RIJONX|RIJK
VeryTightSCF', RI-MP2 with '%mp2 RI true' and no frozen core) and the total energy is written to the json. The
ctest qm_basis_orca_* compares curcuma with the built-in sets (CURCUMA_QM_BASIS_EMBEDDED_ONLY=1) against it.

Cases: HF with def2-SVP, def2-SV(P), def2-SVPD, def2-mTZVP, def2-mTZVPP (the orbital sets with at most d shells
for H-Ne, the engine's limit); on BeH2 HF with def2-SVP, -TZVP, -TZVPP, -TZVPD, -TZVPPD (d at most for H-Be); RI-JK with
def2/JK and RI-J with def2-mTZVPP/J on def2-SVP; RI-MP2 on def2-SVP with def2-SVP/C, def2-TZVP/C, def2-TZVPP/C,
def2-QZVPP/C and on def2-SVPD with def2-SVPD/C, def2-TZVPD/C, def2-TZVPPD/C, def2-QZVPPD/C (the ORCA /C sets are
def2-<X>-RIFIT for H-Ne, checked by scripts/qm_basis_export_tables.py --compare-orca).

With --pyscf PYTHON the RI-MP2 cases are also run with PySCF (RHF conv_tol 1e-13, DFMP2 with PySCF's own
def2-<x>-ri auxiliary data, cart=False) and stored as "pyscf": ORCA's RI-MP2 sits ~1e-8 Eh off curcuma and PySCF,
which agree to ~1e-11, so the ctest checks RI-MP2 against PySCF tightly and against ORCA at 2e-8.

Usage: python3 scripts/qm_basis_orca_reference.py --orca /path/to/orca [--pyscf /path/to/python] [--work DIR]

Claude Generated (Nov 2026). Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse
import json
import os
import re
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.path.join(SRC, "test_cases", "qm_basis_embed", "reference_orca.json")
MOLS = {"H2O": "test_cases/qm_1e/H2O.xyz", "NH3": "test_cases/qm_1e/NH3.xyz", "BeH2": "test_cases/qm_1e/BeH2.xyz"}

# case id -> (curcuma method, curcuma qm settings, ORCA keywords, ORCA %basis aux, ORCA %mp2)
CASES = {
    "hf_def2-SVP": ("hf", {"basis": "def2-SVP"}, "RHF def2-SVP NoRI", "", ""),
    "hf_def2-SV(P)": ("hf", {"basis": "def2-SV(P)"}, "RHF def2-SV(P) NoRI", "", ""),
    "hf_def2-SVPD": ("hf", {"basis": "def2-SVPD"}, "RHF def2-SVPD NoRI", "", ""),
    "hf_def2-mTZVP": ("hf", {"basis": "def2-mTZVP"}, "RHF def2-mTZVP NoRI", "", ""),
    "hf_def2-mTZVPP": ("hf", {"basis": "def2-mTZVPP"}, "RHF def2-mTZVPP NoRI", "", ""),
    "rijk_def2-SVP_def2-universal-JKFIT": ("hf", {"basis": "def2-SVP", "ri": "jk", "ri_aux": "def2-universal-JKFIT"},
                                           "RHF def2-SVP RIJK def2/JK", "", ""),
    "rij_def2-SVP_def2-mTZVPP-RIJ": ("hf", {"basis": "def2-SVP", "ri": "j", "ri_aux": "def2-mTZVPP-RIJ"},
                                     "RHF def2-SVP RIJONX", 'AuxJ "def2-mTZVPP/J"', ""),
}
for orb, auxes in (("def2-SVP", ("SVP", "TZVP", "TZVPP", "QZVPP")), ("def2-SVPD", ("SVPD", "TZVPD", "TZVPPD", "QZVPPD"))):
    for x in auxes:
        CASES["rimp2_%s_def2-%s-RIFIT" % (orb, x)] = (
            "ri-mp2", {"basis": orb, "ri_mp2_aux": "def2-%s-RIFIT" % x}, "RHF MP2 %s NoRI NoFrozenCore" % orb,
            'AuxC "def2-%s/C"' % x, "RI true")


# BeH2 only: the TZVP family has at most d shells for H-Be (chosen when the engine had no f shells; f/g cases since
# Nov 2026: scripts/qm_fg_reference.py, test_cases/qm_fg)
BEH2_CASES = {"hf_" + b: ("hf", {"basis": b}, "RHF %s NoRI" % b, "", "")
              for b in ("def2-SVP", "def2-TZVP", "def2-TZVPP", "def2-TZVPD", "def2-TZVPPD")}


def run_orca(orca, work, xyz, kw, aux, mp2):
    with open(os.path.join(SRC, xyz)) as f:
        atoms = f.read().splitlines()[2:]
    inp = "! %s VeryTightSCF\n" % kw
    if aux:
        inp += "%%basis %s end\n" % aux
    if mp2:
        inp += "%%mp2 %s end\n" % mp2
    inp += "* xyz 0 1\n" + "\n".join(a for a in atoms if a.strip()) + "\n*\n"
    d = tempfile.mkdtemp(dir=work)
    with open(os.path.join(d, "job.inp"), "w") as f:
        f.write(inp)
    env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(orca))
    out = subprocess.run([orca, "job.inp"], cwd=d, env=env, capture_output=True, text=True).stdout
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", out)
    if not m:
        raise RuntimeError("ORCA failed in " + d)
    return float(m[-1]), inp


PYSCF_SNIPPET = r"""
import sys
from pyscf import gto, scf, df
from pyscf.mp import dfmp2
xyz, basis, aux = sys.argv[1:4]
atoms = "; ".join(l for l in open(xyz).read().splitlines()[2:] if l.strip())
mol = gto.M(atom=atoms, basis=basis, cart=False, verbose=0)
mf = scf.RHF(mol); mf.conv_tol = 1e-13; mf.conv_tol_grad = 1e-10; mf.kernel()
pt = dfmp2.DFMP2(mf); pt.with_df = df.DF(mol, auxbasis=aux); pt.kernel()
print("%.12f" % (mf.e_tot + pt.e_corr))
"""


def run_pyscf(python, xyz, qm):
    basis = qm["basis"].lower()
    aux = qm["ri_mp2_aux"].lower().replace("-rifit", "-ri")
    out = subprocess.run([python, "-c", PYSCF_SNIPPET, os.path.join(SRC, xyz), basis, aux],
                         capture_output=True, text=True, check=True).stdout
    return float(out.split()[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", required=True)
    ap.add_argument("--work", default=None)
    ap.add_argument("--pyscf", default=None, help="python with PySCF for the RI-MP2 cross-check")
    a = ap.parse_args()
    work = a.work or tempfile.mkdtemp(prefix="qm_basis_orca_")
    os.makedirs(work, exist_ok=True)
    ref = {"_comment": [
        "ORCA 6.1.0 total energies (Eh) for the built-in basis-set tests (ctest qm_basis_orca_*). Claude Generated (Nov 2026).",
        "Written by scripts/qm_basis_orca_reference.py: '! RHF <basis> NoRI|RIJK def2/JK|RIJONX VeryTightSCF', RI-MP2 =",
        "'! RHF MP2 <basis> NoRI NoFrozenCore' + '%mp2 RI true' + '%basis AuxC \"def2-<X>/C\"'.",
        "curcuma side: dump_qm_energy with the settings in 'curcuma', CURCUMA_QM_BASIS_EMBEDDED_ONLY=1.",
        "pyscf (RI-MP2 only): PySCF 2.x DFMP2, own def2-<x>-ri data, RHF conv_tol 1e-13; ORCA's RI-MP2 is ~1e-8 Eh off both."]}
    for mol, xyz in MOLS.items():
        ref[mol] = {"xyz": "../../" + xyz}
        for case, (method, qm, kw, aux, mp2) in (BEH2_CASES if mol == "BeH2" else CASES).items():
            e, inp = run_orca(a.orca, work, xyz, kw, aux, mp2)
            ref[mol][case] = {"method": method, "curcuma": qm, "orca": e, "orca_input": inp.splitlines()[:3]}
            line = "%-4s %-40s %.12f" % (mol, case, e)
            if a.pyscf and method == "ri-mp2":
                ref[mol][case]["pyscf"] = run_pyscf(a.pyscf, xyz, qm)
                line += "  PySCF %.12f" % ref[mol][case]["pyscf"]
            print(line)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=1)
        f.write("\n")
    print("wrote", os.path.relpath(OUT, SRC))


if __name__ == "__main__":
    main()
