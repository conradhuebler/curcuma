#!/usr/bin/env python3
"""Correctness checks of the local-MP2 prototype before any truncation is measured.

For each small system:
  1. manual frozen-core RI-MP2 (lmp2_common.rimp2_canonical) vs PySCF mp.dfmp2.DFMP2 (same aux, same orbitals)
  2. canonical (non-RI) frozen-core MP2 with exact integrals on the same orbitals (size of the RI error)
  3. LMP2 equations with complete pair spaces built from all PAOs (force_pao) in localised (PM) orbitals:
     iterative result must equal RI-MP2; semicanonical (F_ik, k != i, dropped) must NOT (shows the coupling matters)
  4. semicanonical solve with canonical orbitals (loc=can) equals RI-MP2
Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/check_solver.py [systems...]
Claude Generated (Oct 2026).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp2_common as lc  # first: sets the thread environment before numpy loads
import numpy as np
from pyscf import df, mp, scf



def main():
    names = sys.argv[1:] or ["alkane4", "water4", "acene1"]
    rows = []
    for name in names:
        s = lc.setup(name, loc="pm")
        Bao = lc.b_local_ao(s)
        Bc = lc.b_canonical(s, Bao)
        e_ri = lc.rimp2_canonical(s, Bc)
        # 1. PySCF DFMP2 on the same orbitals
        mf = scf.RHF(s.mol)
        mf.mo_coeff, mf.mo_energy = s.C, s.eps
        mf.mo_occ = np.array([2.0] * s.nocc + [0.0] * s.nvir)
        mf.e_tot = s.ehf
        mf.converged = True  # PySCF then takes mo_energy as the (diagonal) Fock: canonical MP2 on these orbitals
        pt = mp.dfmp2.DFMP2(mf, frozen=s.nfc)
        pt.with_df = df.DF(s.mol, auxbasis=lc.AUX_JKFIT)
        pt.kernel()
        # 2. canonical MP2, exact integrals
        pc = mp.mp2.MP2(mf, frozen=s.nfc)
        pc.kernel()
        # 3. full PAO pair spaces, localised orbitals
        pairs = lc.pair_list(s)
        allat = list(range(s.natm))
        sp = lc.pair_space(s, allat, force_pao=True)
        spaces = [sp] * len(pairs)
        Ks = [lc.pair_ints(Bao, sp[0], i, j) for i, j in pairs]
        r_it = lc.solve_lmp2(s, pairs, spaces, Ks, coupled=True, e_tol=1e-11, r_tol=1e-9)
        r_sc = lc.solve_lmp2(s, pairs, spaces, Ks, coupled=False)
        # 4. canonical orbitals, semicanonical = exact
        lc.localise(s, "can")
        Bao_c = lc.b_local_ao(s)
        sp_c = lc.pair_space(s, allat)
        Ks_c = [lc.pair_ints(Bao_c, sp_c[0], i, j) for i, j in pairs]
        r_can = lc.solve_lmp2(s, pairs, [sp_c] * len(pairs), Ks_c, coupled=False)
        row = dict(system=name, nbf=s.n, naux=s.naux, nval=s.nval, m_pao=sp[0].shape[1], nvir=s.nvir,
                   e_rimp2=e_ri, d_dfmp2=e_ri - pt.e_corr, e_canon=pc.e_corr, ri_err=e_ri - pc.e_corr,
                   d_lmp2_iter=r_it["energy"] - e_ri, it=r_it["iterations"], d_semican_loc=r_sc["energy"] - e_ri,
                   d_can_semi=r_can["energy"] - e_ri)
        rows.append(row)
        print(row, flush=True)
    print("\n| system | nbf | naux | o_val | PAO space | E_corr RI-MP2 (fc) | - PySCF DFMP2 | RI - canonical MP2 |"
          " LMP2 iter (full PAO) - RI | iterations | semican. local - RI | semican. canonical - RI |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {lc.label(r['system'])} | {r['nbf']} | {r['naux']} | {r['nval']} | {r['m_pao']} (nvir {r['nvir']}) |"
              f" {r['e_rimp2']:.10f} | {r['d_dfmp2']:.1e} | {r['ri_err']:.2e} | {r['d_lmp2_iter']:.1e} | {r['it']} |"
              f" {r['d_semican_loc']:.2e} | {r['d_can_semi']:.1e} |")
    lc.save_json("check_solver.json", rows)


if __name__ == "__main__":
    main()
