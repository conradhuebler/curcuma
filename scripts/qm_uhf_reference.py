#!/usr/bin/env python3
"""References for the open-shell native QM engine (stage 5b: UHF/UKS), test_cases/qm_uhf/reference.json,
ctest -L "qm_uhf|qm_uks".

Every case runs one curcuma method on one structure of test_cases/structures with charge and -spin (unpaired electrons
= multiplicity - 1). Reference values written by this script:

  pyscf   PySCF 2.14 / libxc 7.0.0 (venv ~/.venvs/qmref), spherical basis by PySCF name (def2-ECP for Rb-Rn);
          geometry Angstrom -> Bohr with curcuma's CODATA 2018 factor; mol.charge, mol.spin = 2 S.
            hf    scf.UHF, init_guess atom, conv_tol 1e-13, conv_tol_grad 1e-10 (gradient cases continued at 1e-11),
                  energy, <S^2> (mf.spin_square()), gradient
            ks    dft.UKS on curcuma's grid (Treutler x Lebedev, original Becke, no pruning, Treutler radii adjust,
                  small_rho_cutoff 0): lda = LDA_X,LDA_C_VWN | LDA_X,LDA_C_VWN_RPA, pbe = GGA_X_PBE,GGA_C_PBE,
                  b3lyp = HYB_GGA_XC_B3LYP5 (vwn5) | HYB_GGA_XC_B3LYP (vwn_rpa); gradient with grid_response = True
            ri    density_fit(auxbasis = def2-universal-jkfit) (RI-JK; with only_dfj = True RI-J), gradient with
                  auxbasis_response
  orca    ORCA 6.1.0 total energy and <S^2>: '! UHF <basis> NoRI VeryTightSCF', '* xyz <charge> <multiplicity>'
          (only with --orca, HF cases that name an ORCA input).
  kernel  libxc spin-polarised point values (eps, v_rho[2], v_sigma[3]) of LDA_C_VWN, LDA_C_VWN_RPA, LDA_C_PW_MOD,
          GGA_C_PBE, GGA_C_LYP, LDA_X, GGA_X_PBE, GGA_X_B88 on a fixed set of (rho_a, rho_b, sigma_aa, sigma_ab,
          sigma_bb) for ctest qm_uks_kernel.

Usage: ~/.venvs/qmref/bin/python scripts/qm_uhf_reference.py [--orca /path/to/orca] [--only CASE ...]

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
from pyscf import dft, gto, scf
from pyscf.dft import gen_grid, libxc, radi

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.normpath(os.path.join(HERE, ".."))
TC = os.path.join(SRC, "test_cases")
OUT = os.path.join(TC, "qm_uhf", "reference.json")
BOHR = 0.529177210903  # curcuma units.h
XC = {("lda", "vwn5"): "LDA_X,LDA_C_VWN", ("lda", "vwn_rpa"): "LDA_X,LDA_C_VWN_RPA", ("pbe", None): "GGA_X_PBE,GGA_C_PBE",
      ("b3lyp", "vwn5"): "HYB_GGA_XC_B3LYP5", ("b3lyp", "vwn_rpa"): "HYB_GGA_XC_B3LYP"}

CASES = {}


def add(cid, xyz, method, qm, charge=0, spin=1, gradient=False, orca=None, label="qm_uhf", pyscf=True):
    CASES[cid] = {"xyz": "structures/" + xyz, "method": method, "qm": qm, "charge": charge, "spin": spin,
                  "gradient": gradient, "orca": orca, "label": label}
    if not pyscf:
        CASES[cid]["fd_only"] = True  # no PySCF value (the checker's fd mode uses curcuma alone)


SVP = {"basis": "def2-SVP"}
TZVP = {"basis": "def2-TZVP"}
# --- UHF energies (label qm_uhf) ---
add("uhf_oh_svp", "small/hydroxyl.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_oh_tzvp", "small/hydroxyl.xyz", "hf", TZVP, orca="UHF def2-TZVP NoRI")
add("uhf_nh2_svp", "small/amidogen.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_ch3_svp", "small/methyl.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_o2_svp", "small/dioxygen.xyz", "hf", SVP, spin=2, orca="UHF def2-SVP NoRI")
add("uhf_o2_tzvp", "small/dioxygen.xyz", "hf", TZVP, spin=2, orca="UHF def2-TZVP NoRI")
add("uhf_no_svp", "small/nitric-oxide.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_cn_svp", "small/cyano-radical.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_ch2_svp", "small/methylene-triplet.xyz", "hf", SVP, spin=2, orca="UHF def2-SVP NoRI")
add("uhf_h2op_svp", "small/water-cation.xyz", "hf", SVP, charge=1, orca="UHF def2-SVP NoRI")
add("uhf_li_svp", "atoms/lithium-atom.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_na_svp", "atoms/sodium-atom.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_k_svp", "atoms/potassium-atom.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_cl_tzvp", "atoms/chlorine-atom.xyz", "hf", TZVP, orca="UHF def2-TZVP NoRI")
add("uhf_i_svp", "atoms/iodine-atom.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_cu_svp", "atoms/copper-atom.xyz", "hf", SVP, orca="UHF def2-SVP NoRI")
add("uhf_oh_svp_rij", "small/hydroxyl.xyz", "hf", {"basis": "def2-SVP", "ri": "j", "ri_aux": "def2-universal-JKFIT"})
add("uhf_oh_svp_rijk", "small/hydroxyl.xyz", "hf", {"basis": "def2-SVP", "ri": "jk", "ri_aux": "def2-universal-JKFIT"})
# --- UHF gradients ---
add("uhfgrad_ch3p_svp", "small/methyl.pyramidal.xyz", "hf", SVP, gradient=True)
add("uhfgrad_ch2d_svp", "small/methylene-triplet.distorted.xyz", "hf", SVP, spin=2, gradient=True)
add("uhfgrad_ch2d_svp_rijk", "small/methylene-triplet.distorted.xyz", "hf",
    {"basis": "def2-SVP", "ri": "jk", "ri_aux": "def2-universal-JKFIT"}, spin=2, gradient=True)
add("uhfgrad_ch2d_svp_rij", "small/methylene-triplet.distorted.xyz", "hf",
    {"basis": "def2-SVP", "ri": "j", "ri_aux": "def2-universal-JKFIT"}, spin=2, gradient=True)
add("uhfgrad_h2op_svp", "small/water-cation.xyz", "hf", SVP, charge=1, gradient=True)
# hf-3c open shell (UHF/MINIX + D3(BJ) + gCP + SRB): FD gate only; ORCA comparison in docs (stage 5b)
add("hf3cgrad_ch3p", "small/methyl.pyramidal.xyz", "hf-3c", {}, gradient=True, pyscf=False)
# --- UKS (label qm_uks), grid sg1 ---
G1 = {"grid": "sg1"}
# Not used for UKS: NO (PySCF's UKS/LDA SCF did not reach conv_tol 1e-13 on this grid) and OH (pi^3, a degenerate pair:
# on a finite grid the energy depends on the orientation of the hole, curcuma and PySCF landed 1.2e-7 Eh apart for LDA;
# Nov 2026). NH2 (2B1) and CH3 (2A2'') are non-degenerate; O2 (pi^2 triplet) is the symmetric pi_x pi_y occupation.
for mol, xyz, sp in (("nh2", "small/amidogen.xyz", 1), ("o2", "small/dioxygen.xyz", 2), ("ch3", "small/methyl.xyz", 1)):
    add("uks_lda_vwn5_%s_svp" % mol, xyz, "lda", dict(SVP, vwn_variant="vwn5", **G1), spin=sp, label="qm_uks")
    add("uks_pbe_%s_svp" % mol, xyz, "pbe", dict(SVP, **G1), spin=sp, label="qm_uks")
    add("uks_b3lyp_vwn5_%s_svp" % mol, xyz, "b3lyp", dict(SVP, vwn_variant="vwn5", **G1), spin=sp, label="qm_uks")
add("uks_lda_vwnrpa_nh2_svp", "small/amidogen.xyz", "lda", dict(SVP, vwn_variant="vwn_rpa", **G1), label="qm_uks")
add("uks_b3lyp_vwnrpa_o2_svp", "small/dioxygen.xyz", "b3lyp", dict(SVP, vwn_variant="vwn_rpa", **G1), spin=2,
    label="qm_uks")
# Cl and I atoms (degenerate p^5) are not used for UKS: PySCF's UKS SCF did not converge to 1e-13 on it (Nov 2026)
add("uks_b3lyp_vwn5_na_tzvp", "atoms/sodium-atom.xyz", "b3lyp", dict(TZVP, vwn_variant="vwn5", **G1), label="qm_uks")
add("uks_pbe_li_svp", "atoms/lithium-atom.xyz", "pbe", dict(SVP, **G1), label="qm_uks")
add("uks_b3lyp_vwn5_nh2_svp_rijk", "small/amidogen.xyz", "b3lyp",
    dict(SVP, vwn_variant="vwn5", ri="jk", ri_aux="def2-universal-JKFIT", **G1), label="qm_uks")
# --- UKS gradients ---
add("uksgrad_pbe_ch2d_svp", "small/methylene-triplet.distorted.xyz", "pbe", dict(SVP, **G1), spin=2, gradient=True,
    label="qm_uks")
add("uksgrad_b3lyp_ch3p_svp", "small/methyl.pyramidal.xyz", "b3lyp", dict(SVP, vwn_variant="vwn5", **G1),
    gradient=True, label="qm_uks")
add("uksgrad_lda_h2op_svp", "small/water-cation.xyz", "lda", dict(SVP, vwn_variant="vwn5", **G1), charge=1,
    gradient=True, label="qm_uks")


def read_mol(c):
    lines = open(os.path.join(TC, c["xyz"])).read().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        atoms.append([s, [float(x) / BOHR, float(y) / BOHR, float(z) / BOHR]])
    b = c["qm"]["basis"].lower()
    from pyscf.data.elements import charge as zcharge
    ecp = {a[0]: b for a in atoms if zcharge(a[0]) > 36}  # def2-ECP for Rb-Rn (curcuma -qm.ecp auto)
    return gto.M(atom=atoms, unit="Bohr", basis=b, ecp=ecp, charge=c["charge"], spin=c["spin"], verbose=0)


def grid(mf):
    mf.grids.atom_grid = (50, 194)
    mf.grids.prune = None
    mf.grids.radi_method = radi.treutler
    mf.grids.becke_scheme = gen_grid.original_becke
    mf.grids.radii_adjust = radi.treutler_atomic_radii_adjust
    mf.small_rho_cutoff = 0.0


def pyscf_case(c):
    q = c["qm"]
    mol = read_mol(c)
    m = c["method"]
    if m == "hf":
        mf = scf.UHF(mol)
    else:
        mf = dft.UKS(mol)
        mf.xc = XC[(m, q.get("vwn_variant"))]
        grid(mf)
    if q.get("ri") in ("j", "jk"):
        mf = mf.density_fit(auxbasis="def2-universal-jkfit", only_dfj=(q["ri"] == "j"))
    mf.init_guess = "atom"
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 500
    mf.kernel()
    if not mf.converged:  # second try from the result (Cu, K)
        mf.kernel(dm0=mf.make_rdm1())
    assert mf.converged, c
    out = {"nbf": mol.nao_nr()}
    if c["gradient"]:
        mf.conv_tol_grad = 1e-11
        mf.kernel(dm0=mf.make_rdm1())
    out["energy"] = float(mf.e_tot)
    out["s2"] = float(mf.spin_square()[0])
    if c["gradient"]:
        g = mf.nuc_grad_method()
        if m != "hf":
            g.grid_response = True
        if q.get("ri"):
            g.auxbasis_response = True
        out["gradient"] = np.asarray(g.kernel()).tolist()
    return out


def run_orca(orca, work, c):
    with open(os.path.join(TC, c["xyz"])) as f:
        lines = f.read().splitlines()
    n = int(lines[0].split()[0])
    inp = "! %s VeryTightSCF\n* xyz %d %d\n%s\n*\n" % (c["orca"], c["charge"], c["spin"] + 1, "\n".join(lines[2:2 + n]))
    d = tempfile.mkdtemp(dir=work)
    with open(os.path.join(d, "job.inp"), "w") as f:
        f.write(inp)
    env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(orca))
    res = subprocess.run([orca, "job.inp"], cwd=d, env=env, capture_output=True, text=True).stdout
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", res)
    s2 = re.findall(r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)", res)
    if not m:
        raise RuntimeError("ORCA failed in " + d)
    return {"energy": float(m[-1]), "s2": float(s2[-1]) if s2 else None, "input": inp.splitlines()[0]}


def kernel_points():
    """Point set of the kernel test: (rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb), incl. rho_b = 0 and equal spins."""
    pts = []
    for ra, rb in ((0.3, 0.1), (0.01, 0.004), (1.2, 1e-14), (1e-3, 1e-4), (5.0, 2.0), (0.2, 0.2), (0.05, 0.0),
                   (2e-6, 1e-6), (30.0, 29.0), (0.7, 0.01)):
        for saa, sbb, cosab in ((0.0, 0.0, 0.0), (0.05, 0.01, 0.6), (1e-4, 2e-5, -0.3), (2.0, 0.5, 1.0)):
            sab = cosab * np.sqrt(saa * sbb)
            if rb < 1e-10:  # a (nearly) empty spin channel has no gradient (sigma_bb > 0 at rho_b = 1e-14 is unphysical)
                sbb, sab = 0.0, 0.0
            pts.append([ra, rb, saa, sab, sbb])
    return pts


def kernel_reference():
    pts = np.array(kernel_points())
    ra, rb, saa, sab, sbb = pts.T
    ga = np.sqrt(saa)
    gb = np.sqrt(sbb)
    cos = np.where(ga * gb > 0, sab / np.maximum(ga * gb, 1e-300), 0.0)
    R = np.zeros((2, 4, len(pts)))
    R[0, 0], R[1, 0] = ra, rb
    R[0, 1] = ga
    R[1, 1] = gb * cos
    R[1, 2] = gb * np.sqrt(np.maximum(1 - cos ** 2, 0))
    out = {"points": pts.tolist(), "kernels": {}}
    for name in ("LDA_X", "LDA_C_VWN", "LDA_C_VWN_RPA", "LDA_C_PW_MOD", "GGA_C_PBE", "GGA_C_LYP", "GGA_X_PBE", "GGA_X_B88"):
        gga = name.startswith("GGA")
        exc, vxc = libxc.eval_xc(name, R if gga else R[:, 0], spin=1, deriv=1)[:2]
        k = {"eps": np.asarray(exc).tolist(), "vrho": np.asarray(vxc[0]).tolist()}
        if gga:
            k["vsigma"] = np.asarray(vxc[1]).tolist()
        out["kernels"][name] = k
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orca", default=None)
    ap.add_argument("--only", nargs="*", default=None)
    a = ap.parse_args()
    ref = {"generator": "scripts/qm_uhf_reference.py", "pyscf": pyscf.__version__, "libxc": libxc.__version__,
           "cases": {}}
    if os.path.exists(OUT):
        ref["cases"] = json.load(open(OUT)).get("cases", {})
    ref["kernel"] = kernel_reference()
    work = tempfile.mkdtemp(prefix="qm_uhf_orca_")
    for cid, c in CASES.items():
        if a.only and cid not in a.only:
            continue
        r = {k: v for k, v in c.items() if k != "orca"}
        if not c.get("fd_only"):
            r["pyscf"] = pyscf_case(c)
        if a.orca and c["orca"]:
            r["orca"] = run_orca(a.orca, work, c)
        elif cid in ref["cases"] and "orca" in ref["cases"][cid]:
            r["orca"] = ref["cases"][cid]["orca"]
        ref["cases"][cid] = r
        with open(OUT, "w") as f:  # after every case: a failing case does not lose the others
            json.dump(ref, f, indent=1)
            f.write("\n")
        if "pyscf" not in r:
            print("%-30s (FD only)" % cid)
            continue
        print("%-30s pyscf %.12f <S2> %.6f%s" % (cid, r["pyscf"]["energy"], r["pyscf"]["s2"],
                                                ("  orca %.12f <S2> %s" % (r["orca"]["energy"], r["orca"]["s2"]))
                                                if "orca" in r else ""), flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(ref, f, indent=1)
        f.write("\n")
    print("wrote", os.path.relpath(OUT, SRC))


if __name__ == "__main__":
    main()
