#!/usr/bin/env python3
"""Reference RI-MP2 energies for the native RI-MP2 (test_cases/qm_mp2/reference_ri.json).

PySCF RHF (spherical def2-SVP from PySCF's library, = curcuma's def2-SVP.dat; conv_tol 1e-12), then
density-fitted MP2 with the SAME auxiliary basis data as curcuma: the Turbomole files
src/core/energy_calculators/qm_methods/def2-universal-{JKFIT,JFIT}.turbomole are parsed here and handed
to PySCF as a basis dict (contraction coefficients as in the file; the fit depends only on the span).
Two independent PySCF evaluations per case, which must agree:
  - mp.dfmp2.DFMP2 with with_df = df.DF(mol, auxbasis=<dict>)
  - a manual evaluation: B = df.incore.cholesky_eri(mol, auxmol) (L^-1 (P|mn)), B_ia = C_o^T B C_v,
    E_OS/E_SS from (ia|jb) = sum_P B_ia^P B_jb^P (this one gives the OS/SS split stored below).
Optional SCF variant "rijk": RHF with mf.density_fit(auxbasis=JKFIT dict) (RI-JK SCF) before the
RI-MP2 step, as curcuma -qm.ri jk.
Run with the PySCF venv:
    ~/.venvs/qmref/bin/python scripts/rimp2_reference.py            # writes the json
    ~/.venvs/qmref/bin/python scripts/rimp2_reference.py name=x.xyz  # print only
The per-case "orca" entries and _sources.orca in the json were added by hand from ORCA 6.1 runs (see the
WP10 doc); re-running this script drops them.
Claude Generated (Oct 2026, RI-MP2).
"""
import json, os, sys
import numpy as np
import pyscf
from pyscf import gto, scf, mp, df

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
QMDIR = os.path.join(ROOT, "src/core/energy_calculators/qm_methods")
OUT = os.path.join(ROOT, "test_cases/qm_mp2/reference_ri.json")
LVAL = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5}

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


def parse_turbomole(path):
    """Turbomole basis file -> {Symbol: [[L, [exp, coef], ...], ...]} (file order)."""
    out, elem, lines = {}, None, open(path).read().splitlines()
    k = 0
    while k < len(lines):
        s = lines[k].strip()
        k += 1
        if not s or s.startswith("#") or s.startswith("$"):
            continue
        if s == "*":
            continue
        tok = s.split()
        if tok[0].isalpha():  # "h def2-universal-JKFIT": element header
            elem = tok[0].capitalize()
            out[elem] = []
            continue
        if len(tok) == 2 and tok[0].isdigit() and tok[1] in LVAL:  # "3 s": shell header
            nprim, L = int(tok[0]), LVAL[tok[1]]
            sh = [L]
            for _ in range(nprim):
                e, c = lines[k].replace("D", "E").replace("d", "e").split()[:2]
                k += 1
                sh.append([float(e), float(c)])
            out[elem].append(sh)
    return out


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    return [(a.split()[0], tuple(map(float, a.split()[1:4]))) for a in L[2:2 + n]]


def frozen_core(atoms):
    return sum(0 if gto.charge(s) <= 2 else 1 for s, _ in atoms)


def manual_rimp2(mol, mf, auxdict, nfc):
    auxmol = df.make_auxmol(mol, auxdict)
    Bp = df.incore.cholesky_eri(mol, auxmol=auxmol, aosym="s1")  # (naux, n*n), L^-1 (P|mn)
    n = mol.nao_nr()
    nocc = mol.nelectron // 2
    C, e = mf.mo_coeff, mf.mo_energy
    Co, Cv = C[:, nfc:nocc], C[:, nocc:]
    B = np.einsum("Pmn,mi,na->iaP", Bp.reshape(-1, n, n), Co, Cv, optimize=True)
    eo, ev = e[nfc:nocc], e[nocc:]
    eos = ess = 0.0
    for i in range(B.shape[0]):
        K = np.einsum("aP,jbP->ajb", B[i], B, optimize=True)  # (ia|jb)
        D = eo[i] + eo[None, :, None] - ev[:, None, None] - ev[None, None, :]
        Kx = K.transpose(2, 1, 0)  # (ib|ja)
        eos += np.sum(K * K / D)
        ess += np.sum(K * (K - Kx) / D)
    return eos, ess, auxmol.nao_nr()


def main():
    extra = sys.argv[1:]
    cases = CASES if not extra else [tuple(x.split("=", 1)) for x in extra]
    aux = {"def2-universal-JKFIT": parse_turbomole(os.path.join(QMDIR, "def2-universal-JKFIT.turbomole")),
           "def2-universal-JFIT": parse_turbomole(os.path.join(QMDIR, "def2-universal-JFIT.turbomole"))}
    out = {"_about": "PySCF RI-MP2/def2-SVP references for ctest -L qm_mp2 (scripts/rimp2_reference.py)",
           "_sources": {"pyscf": pyscf.__version__, "basis": "def2-svp (PySCF library), spherical",
                        "aux": "def2-universal-{JKFIT,JFIT}.turbomole of curcuma, parsed by the script",
                        "scf": "RHF conv_tol 1e-12 (exact) or density_fit(JKFIT) (rijk)",
                        "mp2": "manual B = cholesky_eri L^-1 (P|mn); cross-checked against mp.dfmp2.DFMP2"},
           "_units": "Hartree", "cases": {}}
    for name, rel in cases:
        path = rel if os.path.isabs(rel) else os.path.join(ROOT, "test_cases", rel)
        atoms = read_xyz(path)
        mol = gto.M(atom=[[s, x] for s, x in atoms], basis="def2-svp", unit="Angstrom", cart=False, verbose=0)
        nfc = frozen_core(atoms)
        entry = {"xyz": rel, "nbf": mol.nao_nr()}
        for scfmode in ("exact", "rijk"):
            if scfmode == "rijk" and name not in ("H2O", "formaldehyde", "benzene"):
                continue
            mf = scf.RHF(mol)
            if scfmode == "rijk":
                mf = mf.density_fit(auxbasis=aux["def2-universal-JKFIT"])
            mf.conv_tol = 1e-12
            mf.max_cycle = 200
            e_hf = mf.kernel()
            if not mf.converged:
                sys.exit(f"{name}: PySCF RHF did not converge")
            sub = {"hf": e_hf}
            for auxname, auxd in aux.items():
                if scfmode == "rijk" and auxname != "def2-universal-JKFIT":
                    continue
                res = {}
                for tag, frz in (("full", 0), ("fc", nfc)):
                    eos, ess, naux = manual_rimp2(mol, mf, auxd, frz)
                    pt = mp.dfmp2.DFMP2(mf, frozen=frz if frz else None)
                    pt.with_df = df.DF(mol, auxbasis=auxd)
                    pt.kernel()
                    dd = abs(pt.e_corr - (eos + ess))
                    if dd > 1e-10:
                        sys.exit(f"{name}/{auxname}/{tag}: manual {eos + ess:.12f} vs DFMP2 {pt.e_corr:.12f}")
                    res[tag] = {"frozen": frz, "naux": naux, "e_corr": eos + ess, "e_os": eos, "e_ss": ess,
                                "dfmp2_e_corr": pt.e_corr}
                sub[auxname] = res
                print(f"{name:24s} {scfmode:5s} {auxname:22s} nbf={mol.nao_nr():4d} naux={res['full']['naux']:5d} "
                      f"E_HF={e_hf:.10f} Ecorr={res['full']['e_corr']:.12f} fc={res['fc']['e_corr']:.12f} "
                      f"|man-DFMP2|={abs(res['full']['e_corr'] - res['full']['dfmp2_e_corr']):.1e}")
            entry[scfmode] = sub
        out["cases"][name] = entry
    if not extra:
        json.dump(out, open(OUT, "w"), indent=1)
        print("wrote", OUT)


if __name__ == "__main__":
    main()
