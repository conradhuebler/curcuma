#!/usr/bin/env python3
"""References for the native QM engine beyond Ne: all-electron Na-Kr and def2-ECP elements (test_cases/qm_ecp/reference.json,
ctest -L "qm_z|qm_ecp").

Every case runs one curcuma method on one structure of test_cases/structures with a def2 basis; the elements Rb-Rn get
the def2-ECP (curcuma -qm.ecp auto). Reference values written by this script:

  pyscf  PySCF 2.14 / libxc 7.0.0 (venv ~/.venvs/qmref), spherical basis; basis and ECP by PySCF name (def2-svp ...,
         ecp = the same name), except for elements PySCF's def2 library lacks (the lanthanides, here Ce): then basis
         and ECP are converted from curcuma's built-in BSE text (qm_basis_data_def2_*.h), flagged "bse_text".
         Geometry: Angstrom -> Bohr with curcuma's CODATA 2018 factor.
           hf      scf.RHF, init_guess atom, conv_tol 1e-13, conv_tol_grad 1e-10 (gradient cases continued at 1e-11)
           ks      dft.RKS on curcuma's grid (Treutler x Lebedev, original Becke, no pruning, Treutler radii adjust,
                   small_rho_cutoff 0), b3lyp = HYB_GGA_XC_B3LYP5, pbe = GGA_X_PBE,GGA_C_PBE; gradient with
                   grid_response = True
           ri      density_fit(auxbasis = def2-universal-jkfit) (RI-JK, gradient with auxbasis_response)
           mp2     mp.MP2 with frozen = pyscf.data.elements.chemcore(mol) when frozen_core (ECP electrons subtracted,
                   the convention curcuma follows), ri-mp2 = dfmp2.DFMP2 with auxbasis def2-universal-jkfit
  orca   ORCA 6.1.0 total energy: '! RHF <basis> NoRI VeryTightSCF' (def2-ECP automatically for Rb-Rn). Only with --orca.
         ORCA's ECP energies sit up to 6e-8 Eh from curcuma and PySCF for Xe and Hg (curcuma = PySCF to 1e-11 for Hg;
         for XeF2 curcuma is 2.7e-9 from PySCF and 2.1e-8 from ORCA), so the case stores "orca_tol".

Usage: ~/.venvs/qmref/bin/python scripts/qm_ecp_reference.py [--orca /path/to/orca] [--only CASE ...]

Claude Generated (Nov 2026). Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>. GPL-3.0.
"""
import argparse
import json
import os
import re
import subprocess
import tempfile

import numpy as np
import pyscf
from pyscf import df, dft, gto, mp, scf
from pyscf.data.elements import charge as zcharge, chemcore
from pyscf.dft import gen_grid, libxc, radi
from pyscf.mp import dfmp2

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, ".."))
TC = os.path.join(SRC, "test_cases")
OUT = os.path.join(TC, "qm_ecp", "reference.json")
QMDIR = os.path.join(SRC, "src", "core", "energy_calculators", "qm_methods")
BOHR = 0.529177210903  # curcuma units.h
XC = {"pbe": "GGA_X_PBE,GGA_C_PBE", "b3lyp": "HYB_GGA_XC_B3LYP5"}
DATAFILE = {"def2-SVP": "qm_basis_data_def2_sv.h", "def2-TZVP": "qm_basis_data_def2_tzv.h",
            "def2-ECP": "qm_basis_data_def2_ecp.h"}

CASES = {}


def add(cid, xyz, method, qm, gradient=False, orca=None, label="qm_ecp", orca_tol=None):
    CASES[cid] = {"xyz": "structures/" + xyz, "method": method, "qm": qm, "gradient": gradient, "orca": orca,
                  "label": label}
    if orca_tol:
        CASES[cid]["orca_tol"] = orca_tol


# --- all-electron Na-Kr (label qm_z) ---
add("hf_nacl_def2-SVP", "small/sodium-chloride.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI", label="qm_z")
add("hf_zncl2_def2-SVP", "metals/zinc-dichloride.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI", label="qm_z")
add("hf_hbr_def2-TZVP", "small/hydrogen-bromide.xyz", "hf", {"basis": "def2-TZVP"}, orca="RHF def2-TZVP NoRI", label="qm_z")
add("b3lyp_so2_def2-SVP_sg1", "small/sulfur-dioxide.xyz", "b3lyp", {"basis": "def2-SVP", "grid": "sg1"}, label="qm_z")
add("pbe_kcl_def2-TZVP_sg1", "small/potassium-chloride.xyz", "pbe", {"basis": "def2-TZVP", "grid": "sg1"}, label="qm_z")
add("hfgrad_h2s_def2-SVP", "small/hydrogen-sulfide.d-shell-validation.xyz", "hf", {"basis": "def2-SVP"}, gradient=True,
    label="qm_z")
add("b3lypgrad_so2_def2-SVP_sg1", "small/sulfur-dioxide.xyz", "b3lyp", {"basis": "def2-SVP", "grid": "sg1"},
    gradient=True, label="qm_z")
add("mp2fc_nacl_def2-SVP", "small/sodium-chloride.xyz", "mp2", {"basis": "def2-SVP", "frozen_core": True}, label="qm_z")
# --- def2-ECP elements (label qm_ecp) ---
add("hf_hi_def2-SVP", "small/hydrogen-iodide.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI")
add("hf_hi_def2-TZVP", "small/hydrogen-iodide.xyz", "hf", {"basis": "def2-TZVP"}, orca="RHF def2-TZVP NoRI")
add("hf_ch3i_def2-SVP", "small/iodomethane.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI")
add("hf_csi_def2-SVP", "small/caesium-iodide.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI")
add("hf_xef2_def2-SVP", "small/xenon-difluoride.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI",
    orca_tol=5e-8)
add("hf_hgcl2_def2-SVP", "metals/mercury-dichloride.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI",
    orca_tol=1e-7)
add("hf_pbh4_def2-SVP", "small/plumbane.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI")
add("hf_cef4_def2-SVP", "metals/cerium-tetrafluoride.xyz", "hf", {"basis": "def2-SVP"}, orca="RHF def2-SVP NoRI")
add("b3lyp_hi_def2-SVP_sg1", "small/hydrogen-iodide.xyz", "b3lyp", {"basis": "def2-SVP", "grid": "sg1"})
add("hfgrad_hit_def2-SVP", "small/hydrogen-iodide.tilted.xyz", "hf", {"basis": "def2-SVP"}, gradient=True)
add("hfgrad_hit_def2-TZVP", "small/hydrogen-iodide.tilted.xyz", "hf", {"basis": "def2-TZVP"}, gradient=True)
add("b3lypgrad_hit_def2-SVP_sg1", "small/hydrogen-iodide.tilted.xyz", "b3lyp", {"basis": "def2-SVP", "grid": "sg1"},
    gradient=True)
add("rijkgrad_hit_def2-SVP", "small/hydrogen-iodide.tilted.xyz", "hf",
    {"basis": "def2-SVP", "ri": "jk", "ri_aux": "def2-universal-JKFIT"}, gradient=True)
add("mp2fc_hi_def2-SVP", "small/hydrogen-iodide.xyz", "mp2", {"basis": "def2-SVP", "frozen_core": True})
add("rimp2fc_hi_def2-SVP", "small/hydrogen-iodide.xyz", "ri-mp2",
    {"basis": "def2-SVP", "frozen_core": True, "ri_mp2_aux": "def2-universal-JKFIT"})
# FD-gate-only cases (no PySCF number needed beyond the energy; the checker's fd mode uses curcuma alone)
add("mp2grad_hit_def2-SVP", "small/hydrogen-iodide.tilted.xyz", "mp2", {"basis": "def2-SVP", "frozen_core": True},
    gradient=True)
add("b2plypgrad_hit_def2-SVP_sg1", "small/hydrogen-iodide.tilted.xyz", "b2plyp",
    {"basis": "def2-SVP", "grid": "sg1", "frozen_core": False}, gradient=True)


def element_block(name, Z):
    txt = open(os.path.join(QMDIR, DATAFILE[name])).read()
    key = "kBlocks_" + re.sub(r"[^A-Za-z0-9]", "_", name) + "["
    i = txt.index(key)
    body = txt[i:txt.index("};", i)]
    for m in re.finditer(r"// (\d+) (\w+)\n\s*R\"BASIS\((.*?)\)BASIS\"", body, re.S):
        if int(m.group(1)) == Z:
            return m.group(3)
    raise KeyError((name, Z))


def bse_basis(name, Z):
    lines = [x for x in element_block(name, Z).splitlines() if x.strip()]
    shells, i = [], 1
    while i < len(lines):
        a = lines[i].split()
        n = int(a[1])
        prim = [[float(b.split()[1]), float(b.split()[2])] for b in lines[i + 1:i + 1 + n]]
        shells.append(["SPDFGHI".index(a[0])] + prim)
        i += n + 1
    return shells


def bse_ecp(Z):
    nc, ch, cur = None, {}, None
    for x in element_block("def2-ECP", Z).splitlines():
        x = x.strip()
        if not x or x == "*":
            continue
        if x.startswith("ncore"):
            nc = int(x.replace("=", " ").split()[1])
            continue
        if nc is None:
            continue
        p = x.split()
        if len(p) == 1:
            cur = -1 if len(p[0]) == 1 else "spdfghi".index(p[0][0])
            ch.setdefault(cur, [])
            continue
        ch[cur].append((float(p[2]), float(p[0]), int(p[1])))
    out = []
    for l, terms in sorted(ch.items()):
        r = [[] for _ in range(7)]
        for z, c, n in terms:
            r[n].append([z, c])
        out.append([l, r])
    return [nc, out]


def read_mol(path, basis):
    lines = open(os.path.join(TC, path)).read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        atoms.append([s, [float(x) / BOHR, float(y) / BOHR, float(z) / BOHR]])
    els = sorted(set(a[0] for a in atoms))
    bse = [e for e in els if 57 <= zcharge(e) <= 71]  # lanthanides: not in PySCF's def2 library
    if bse:
        b = {e: bse_basis(basis, zcharge(e)) for e in els}
        ecp = {e: bse_ecp(zcharge(e)) for e in els if zcharge(e) > 36}
        return gto.M(atom=atoms, unit="Bohr", basis=b, ecp=ecp, verbose=0), True
    return gto.M(atom=atoms, unit="Bohr", basis=basis.lower(), ecp=basis.lower(), verbose=0), False


def make_rks(mol, fn):
    mf = dft.RKS(mol)
    mf.xc = XC[fn]
    mf.grids.atom_grid = (50, 194)
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
    mol, bse = read_mol(c["xyz"], q["basis"])
    m = c["method"]
    out = {"nbf": mol.nao_nr(), "bse_text": bse}
    if m in ("hf", "mp2", "ri-mp2"):
        mf = scf.RHF(mol)
        if q.get("ri") == "jk":
            mf = mf.density_fit(auxbasis="def2-universal-jkfit")
        mf.init_guess = "atom"
        mf.conv_tol = 1e-13
        mf.conv_tol_grad = 1e-10
        mf.max_cycle = 300
        mf.kernel()
    elif m in XC:
        mf = make_rks(mol, m)
        mf.kernel()
    else:
        return None  # FD-only case
    assert mf.converged, c
    if c["gradient"]:
        e0 = mf.e_tot
        mf.conv_tol_grad = 1e-11
        mf.kernel(dm0=mf.make_rdm1())
        out["refine_dE"] = float(mf.e_tot - e0)
    e = float(mf.e_tot)
    nf = chemcore(mol) if q.get("frozen_core") else 0
    out["frozen"] = nf
    if m == "mp2":
        if c["gradient"]:
            return None
        e += float(mp.MP2(mf, frozen=nf).kernel()[0])
    elif m == "ri-mp2":
        pt = dfmp2.DFMP2(mf, frozen=nf)
        pt.with_df = df.DF(mol, auxbasis="def2-universal-jkfit")
        pt.kernel()
        e += float(pt.e_corr)
    out["energy"] = e
    if c["gradient"]:
        g = mf.nuc_grad_method()
        if m != "hf":
            g.grid_response = True
        out["gradient"] = np.asarray(g.kernel()).tolist()
    return out


def run_orca(orca, work, path, kw):
    with open(os.path.join(TC, path)) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    inp = "! %s VeryTightSCF\n* xyz 0 1\n%s\n*\n" % (kw, "\n".join(lines[2:2 + n]))
    d = tempfile.mkdtemp(dir=work)
    with open(os.path.join(d, "job.inp"), "w") as f:
        f.write(inp)
    env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(orca))
    res = subprocess.run([orca, "job.inp"], cwd=d, env=env, capture_output=True, text=True).stdout
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", res)
    if not m:
        raise RuntimeError("ORCA failed in " + d)
    return float(m[-1]), inp.splitlines()[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", default=None)
    ap.add_argument("--only", nargs="*", default=None)
    a = ap.parse_args()
    ref = {"generator": "scripts/qm_ecp_reference.py", "pyscf": pyscf.__version__, "libxc": libxc.__version__,
           "cases": {}}
    if os.path.exists(OUT):
        ref["cases"] = json.load(open(OUT)).get("cases", {})
    work = tempfile.mkdtemp(prefix="qm_ecp_orca_")
    for cid, c in CASES.items():
        if a.only and cid not in a.only:
            continue
        r = {k: v for k, v in c.items() if k != "orca"}
        p = pyscf_case(c)
        if p is not None:
            r["pyscf"] = p
        if a.orca and c["orca"]:
            e, inp = run_orca(a.orca, work, c["xyz"], c["orca"])
            r["orca"] = {"energy": e, "input": inp}
        elif cid in ref["cases"] and "orca" in ref["cases"][cid]:
            r["orca"] = ref["cases"][cid]["orca"]
        ref["cases"][cid] = r
        print("%-30s %s%s" % (cid, ("pyscf %.12f" % p["energy"]) if p else "(FD only)",
                              ("  orca %.12f" % r["orca"]["energy"]) if "orca" in r else ""), flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=1)
        f.write("\n")
    print("wrote", os.path.relpath(OUT, SRC))


if __name__ == "__main__":
    main()
