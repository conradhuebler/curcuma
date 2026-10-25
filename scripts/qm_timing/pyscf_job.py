#!/usr/bin/env python3
"""One PySCF job for the QM wall-time comparison (Claude Generated, Oct 2026).

Run with the PySCF venv; prints one JSON line (prefix 'RESULT_JSON ') with energies and
sub-timings. Thread count comes from OMP_NUM_THREADS (set by qm_timing.py) and lib.num_threads.

  pyscf_job.py XYZ METHOD VARIANT THREADS
    METHOD hf    VARIANT exact (scf.RHF, PySCF default: in-core ERIs if they fit max_memory,
                               else integral-direct) | df (density_fit, def2-universal-jkfit)
    METHOD rimp2 VARIANT exact (exact RHF, then DF-MP2 with def2-universal-jkfit, frozen core)
    METHOD ccsd  VARIANT exact (exact RHF, then canonical RCCSD, frozen core)

Convergence: RHF conv_tol 1e-9 (energy), conv_tol_grad 1e-5; CCSD conv_tol 1e-8,
conv_tol_normt 1e-6. Basis def2-SVP (spherical, PySCF library = curcuma's def2-SVP.dat).
"""
import json
import sys
import time

from pyscf import gto, lib, scf, df


def frozen_core(mol):
    """Chemical frozen core as qmint::frozenCoreOrbitals (H-He 0, Li-Ne 1)."""
    return sum(0 if mol.atom_charge(i) <= 2 else 1 for i in range(mol.natm))


def main():
    xyz, method, variant, threads = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    lib.num_threads(threads)
    t0 = time.perf_counter()
    mol = gto.M(atom=xyz, basis="def2-svp", cart=False, verbose=4, max_memory=64000)
    out = {"nbf": mol.nao_nr(), "natoms": mol.natm, "threads": lib.num_threads()}
    if method == "hf" and variant == "df":
        mf = scf.RHF(mol).density_fit(auxbasis="def2-universal-jkfit")
    else:
        mf = scf.RHF(mol)
    mf.conv_tol = 1e-9
    mf.conv_tol_grad = 1e-5
    mf.max_cycle = 100
    t1 = time.perf_counter()
    e_hf = mf.kernel()
    t2 = time.perf_counter()
    out.update({"e_hf": e_hf, "scf_converged": bool(mf.converged), "t_scf": t2 - t1,
                "eri_incore": getattr(mf, "_eri", None) is not None})
    nfc = frozen_core(mol)
    out["frozen"] = nfc
    if method == "rimp2":
        from pyscf.mp import dfmp2
        pt = dfmp2.DFMP2(mf, frozen=nfc)
        pt.with_df = df.DF(mol, auxbasis="def2-universal-jkfit")
        t3 = time.perf_counter()
        ec = pt.kernel()[0]
        out.update({"e_corr": ec, "t_corr": time.perf_counter() - t3})
    elif method == "ccsd":
        from pyscf import cc
        c = cc.CCSD(mf, frozen=nfc)
        c.conv_tol = 1e-8
        c.conv_tol_normt = 1e-6
        c.max_cycle = 100
        t3 = time.perf_counter()
        ec = c.kernel()[0]
        out.update({"e_corr": ec, "t_corr": time.perf_counter() - t3, "cc_converged": bool(c.converged)})
    out["t_total_inside"] = time.perf_counter() - t0
    print("RESULT_JSON " + json.dumps(out), flush=True)


if __name__ == "__main__":
    main()
