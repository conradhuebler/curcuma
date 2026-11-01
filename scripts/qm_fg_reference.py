#!/usr/bin/env python3
"""References for the native QM engine with f and g shells (test_cases/qm_fg/reference.json, ctest -L qm_fg).

Every case runs one curcuma method on one molecule of test_cases/qm_1e with a def2 basis that has f (def2-TZVP*) or
f and g (def2-QZVP*) shells, and stores the reference values written by this script:

  pyscf  PySCF 2.14 / libxc 7.0.0 (venv ~/.venvs/qmref), spherical basis (cart=False), the basis by PySCF name
         (def2-tzvp ...; same numbers as curcuma's built-in BSE copies for H-Ne); cases with cartesian_d use cart=True. Geometry: Angstrom -> Bohr with
         curcuma's CODATA 2018 factor.
           hf      scf.RHF, conv_tol 1e-13, conv_tol_grad 1e-10 (gradient cases
                   continued at 1e-11)
           ks      dft.RKS on curcuma's grid (Treutler x Lebedev, original Becke, no pruning, Treutler radii adjust,
                   small_rho_cutoff 0), conv_tol 1e-12, conv_tol_grad 1e-10; lda = LDA_X,LDA_C_VWN (vwn5),
                   pbe = GGA_X_PBE,GGA_C_PBE, b3lyp = HYB_GGA_XC_B3LYP5 (curcuma default vwn5), tpss
                   = MGGA_X_TPSS,MGGA_C_TPSS; gradients with grid_response = True (curcuma default)
           ri      density_fit(auxbasis = curcuma's def2-universal-JFIT/-JKFIT.turbomole parsed as in
                   scripts/rimp2_reference.py); pure functionals with ri j fit J only (as curcuma); RI-JK gradient
                   with auxbasis_response (PySCF default)
           mp2     mp.MP2 (all electrons), ri-mp2 = dfmp2.DFMP2 with PySCF's def2-<x>-ri (= def2-<X>-RIFIT), ccsd =
                   cc.CCSD (all electrons, conv_tol 1e-12)
           b2plyp  RKS '0.53*HF + 0.47*GGA_X_B88, 0.73*GGA_C_LYP' + 0.27 (E_os + E_ss) of RI-PT2 on the KS orbitals
                   with def2-universal-JKFIT (manual_rimp2 of scripts/rimp2_reference.py), no frozen core
  orca   ORCA 6.1.0 total energy: '! RHF <basis> NoRI VeryTightSCF' (HF), '! RHF <basis> RIJK def2/JK' (RI-JK),
         '! RHF MP2 <basis> NoRI NoFrozenCore' + '%mp2 RI true' + AuxC def2-<X>/C (RI-MP2). Only with --orca.

Usage: ~/.venvs/qmref/bin/python scripts/qm_fg_reference.py [--orca /path/to/orca] [--only CASE ...]

Claude Generated (Nov 2026). Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse
import importlib.util
import json
import os
import re
import subprocess
import tempfile

import numpy as np
import pyscf
from pyscf import cc, df, dft, gto, mp, scf
from pyscf.dft import gen_grid, libxc, radi
from pyscf.mp import dfmp2

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, ".."))
TC = os.path.join(SRC, "test_cases")
OUT = os.path.join(TC, "qm_fg", "reference.json")
QMDIR = os.path.join(SRC, "src", "core", "energy_calculators", "qm_methods")
BOHR = 0.529177210903  # curcuma units.h
LEVELS = {"sg1": (50, 194)}
XC = {"lda": "LDA_X,LDA_C_VWN", "pbe": "GGA_X_PBE,GGA_C_PBE", "b3lyp": "HYB_GGA_XC_B3LYP5",
      "tpss": "MGGA_X_TPSS,MGGA_C_TPSS", "b2plyp": "0.53*HF + 0.47*GGA_X_B88, 0.73*GGA_C_LYP"}
AUXFILE = {"def2-universal-JFIT": "def2-universal-JFIT.turbomole",
           "def2-universal-JKFIT": "def2-universal-JKFIT.turbomole"}

# case id -> dict(molecule, method, qm, gradient, orca keywords/aux/mp2 or None)
CASES = {}


def add(cid, mol, method, qm, gradient=False, orca=None):
    CASES[cid] = {"molecule": mol, "xyz": "qm_1e/%s.xyz" % mol, "method": method, "qm": qm, "gradient": gradient,
                  "orca": orca}


for mol, b in (("H2O", "def2-TZVP"), ("H2O", "def2-TZVPP"), ("H2O", "def2-TZVPD"), ("H2O", "def2-QZVP"),
               ("H2O", "def2-QZVPP"), ("NH3", "def2-TZVPP"), ("CH4", "def2-TZVPPD"), ("HF", "def2-QZVPPD")):
    add("hf_%s_%s" % (mol, b), mol, "hf", {"basis": b}, orca=("RHF %s NoRI" % b, "", ""))
for b in ("def2-TZVP", "def2-QZVP"):
    for fn in ("lda", "pbe", "b3lyp", "tpss"):
        add("%s_H2O_%s_sg1" % (fn, b), "H2O", fn, {"basis": b, "grid": "sg1"})
add("hfcart_H2O_def2-TZVP", "H2O", "hf", {"basis": "def2-TZVP", "cartesian_d": True})  # cartesian 6d/10f, PySCF cart
add("hfcart_H2O_def2-QZVP", "H2O", "hf", {"basis": "def2-QZVP", "cartesian_d": True})  # 15g
add("hfgrad_H2O_def2-TZVP", "H2O", "hf", {"basis": "def2-TZVP"}, gradient=True)
add("hfgrad_H2O_def2-QZVPP", "H2O", "hf", {"basis": "def2-QZVPP"}, gradient=True)
add("hfgrad_NH3_def2-TZVPD", "NH3", "hf", {"basis": "def2-TZVPD"}, gradient=True)
for fn in ("lda", "pbe", "b3lyp", "tpss"):
    add("%sgrad_H2O_def2-TZVP_sg1" % fn, "H2O", fn, {"basis": "def2-TZVP", "grid": "sg1"}, gradient=True)
add("b3lypgrad_H2O_def2-QZVP_sg1", "H2O", "b3lyp", {"basis": "def2-QZVP", "grid": "sg1"}, gradient=True)
add("rijk_H2O_def2-TZVP", "H2O", "hf", {"basis": "def2-TZVP", "ri": "jk", "ri_aux": "def2-universal-JKFIT"},
    gradient=True, orca=("RHF def2-TZVP RIJK def2/JK", "", ""))
add("rijk_H2O_def2-QZVP", "H2O", "hf", {"basis": "def2-QZVP", "ri": "jk", "ri_aux": "def2-universal-JKFIT"},
    orca=("RHF def2-QZVP RIJK def2/JK", "", ""))
add("pberij_H2O_def2-TZVP_sg1", "H2O", "pbe", {"basis": "def2-TZVP", "grid": "sg1", "ri": "j",
                                                "ri_aux": "def2-universal-JFIT"}, gradient=True)
add("mp2_H2O_def2-TZVP", "H2O", "mp2", {"basis": "def2-TZVP", "frozen_core": False})
add("rimp2_H2O_def2-TZVP", "H2O", "ri-mp2", {"basis": "def2-TZVP", "ri_mp2_aux": "def2-TZVP-RIFIT",
                                             "frozen_core": False},
    orca=("RHF MP2 def2-TZVP NoRI NoFrozenCore", 'AuxC "def2-TZVP/C"', "RI true"))
add("rimp2_H2O_def2-QZVPP", "H2O", "ri-mp2", {"basis": "def2-QZVPP", "ri_mp2_aux": "def2-QZVPP-RIFIT",
                                              "frozen_core": False},
    orca=("RHF MP2 def2-QZVPP NoRI NoFrozenCore", 'AuxC "def2-QZVPP/C"', "RI true"))
add("ccsd_H2O_def2-TZVP", "H2O", "ccsd", {"basis": "def2-TZVP", "frozen_core": False})
add("b2plyp_H2O_def2-TZVP_sg1", "H2O", "b2plyp", {"basis": "def2-TZVP", "grid": "sg1", "frozen_core": False})


def rimp2_module():
    spec = importlib.util.spec_from_file_location("rimp2_reference", os.path.join(HERE, "rimp2_reference.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_mol(name, basis, cart=False):
    lines = open(os.path.join(TC, "qm_1e", name + ".xyz")).read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        atoms.append([s, [float(x) / BOHR, float(y) / BOHR, float(z) / BOHR]])
    return gto.M(atom=atoms, unit="Bohr", basis=basis.lower(), cart=cart, verbose=0)


def aux_dict(name):
    return rimp2_module().parse_turbomole(os.path.join(QMDIR, AUXFILE[name]))


def make_rks(mol, fn, level):
    mf = dft.RKS(mol)
    mf.xc = XC[fn]
    mf.grids.atom_grid = LEVELS[level]
    mf.grids.prune = None
    mf.grids.radi_method = radi.treutler
    mf.grids.becke_scheme = gen_grid.original_becke
    mf.grids.radii_adjust = radi.treutler_atomic_radii_adjust
    mf.small_rho_cutoff = 0.0
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 300
    return mf


def pyscf_case(c):
    q = c["qm"]
    mol = read_mol(c["molecule"], q["basis"], bool(q.get("cartesian_d", False)))
    m = c["method"]
    out = {"nbf": mol.nao_nr()}
    if m in ("hf", "mp2", "ri-mp2", "ccsd"):
        mf = scf.RHF(mol)
        if q.get("ri") == "jk":
            mf = mf.density_fit(auxbasis=aux_dict(q["ri_aux"]))
        mf.conv_tol = 1e-13
        mf.conv_tol_grad = 1e-10
        mf.max_cycle = 300
        mf.kernel()
    else:
        mf = make_rks(mol, m, q["grid"])
        if q.get("ri") == "j":
            mf = mf.density_fit(auxbasis=aux_dict(q["ri_aux"]))
        mf.kernel()
    assert mf.converged, c
    if c["gradient"]:
        # gradients: one more pass at conv_tol_grad 1e-11 from the converged density (PySCF's flag may stay
        # False at this level because of rounding; the energy change is reported)
        e0 = mf.e_tot
        mf.conv_tol_grad = 1e-11
        mf.kernel(dm0=mf.make_rdm1())
        out["refine_dE"] = float(mf.e_tot - e0)
    e = float(mf.e_tot)
    if m == "mp2":
        e += float(mp.MP2(mf).kernel()[0])
    elif m == "ri-mp2":
        pt = dfmp2.DFMP2(mf)
        pt.with_df = df.DF(mol, auxbasis=q["ri_mp2_aux"].lower().replace("-rifit", "-ri"))
        pt.kernel()
        e += float(pt.e_corr)
    elif m == "ccsd":
        cs = cc.CCSD(mf)
        cs.conv_tol = 1e-12
        cs.conv_tol_normt = 1e-10
        cs.kernel()
        e += float(cs.e_corr)
    elif m == "b2plyp":
        eos, ess, naux = rimp2_module().manual_rimp2(mol, mf, aux_dict("def2-universal-JKFIT"), 0)
        out["e_ks"] = e
        out["e_os"], out["e_ss"] = float(eos), float(ess)
        e += 0.27 * (eos + ess)
    out["energy"] = e
    if c["gradient"]:
        g = mf.nuc_grad_method()
        if m not in ("hf",):
            g.grid_response = True
        out["gradient"] = np.asarray(g.kernel()).tolist()
    return out


def run_orca(orca, work, mol, kw, aux, mp2):
    with open(os.path.join(TC, "qm_1e", mol + ".xyz")) as f:
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
    res = subprocess.run([orca, "job.inp"], cwd=d, env=env, capture_output=True, text=True).stdout
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", res)
    if not m:
        raise RuntimeError("ORCA failed in " + d)
    return float(m[-1]), inp.splitlines()[:3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", default=None)
    ap.add_argument("--only", nargs="*", default=None)
    a = ap.parse_args()
    ref = {"generator": "scripts/qm_fg_reference.py", "pyscf": pyscf.__version__, "libxc": libxc.__version__,
           "cases": {}}
    if os.path.exists(OUT):
        old = json.load(open(OUT))
        ref["cases"] = old.get("cases", {})
    work = tempfile.mkdtemp(prefix="qm_fg_orca_")
    for cid, c in CASES.items():
        if a.only and cid not in a.only:
            continue
        r = {k: v for k, v in c.items() if k != "orca"}
        r["pyscf"] = pyscf_case(c)
        if a.orca and c["orca"]:
            e, inp = run_orca(a.orca, work, c["molecule"], *c["orca"])
            r["orca"] = {"energy": e, "input": inp}
        elif cid in ref["cases"] and "orca" in ref["cases"][cid]:
            r["orca"] = ref["cases"][cid]["orca"]
        ref["cases"][cid] = r
        print("%-32s pyscf %.12f%s" % (cid, r["pyscf"]["energy"],
                                       ("  orca %.12f" % r["orca"]["energy"]) if "orca" in r else ""), flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=1)
        f.write("\n")
    print("wrote", os.path.relpath(OUT, SRC))


if __name__ == "__main__":
    main()
