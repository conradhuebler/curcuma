#!/usr/bin/env python3
"""Q4 cost part: Cholesky of the global RI metric vs all atom-centred local metric blocks (fit-domain radius R).

Only the geometry and the aux basis are needed (no SCF), so the sizes go beyond the energy systems:
alkanes C_nH_(2n+2) and water clusters (same generators as lmp2_common). Timings: scipy cho_factor (OpenBLAS of the
venv, LMP2_THREADS threads), global = min of 3, local = sum over one block per atom, min of 3 repeats; n=1 machine.
Flop model: global naux^3/3, local sum_A ndom_A^3/3.
Run it alone on the machine (timings).

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/q4_timing.py
Claude Generated (Oct 2026).
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp2_common as lc  # first: sets the thread environment before numpy loads
import numpy as np
import scipy.linalg
from pyscf import df


CASES = [("alkane", n) for n in (8, 16, 32, 64, 128)] + [("water", n) for n in (8, 16, 32, 64, 128)]
RADII = [6, 8, 10]


def tmin(f, rep=3):
    best = 1e9
    for _ in range(rep):
        t = time.perf_counter()
        f()
        best = min(best, time.perf_counter() - t)
    return best


def main():
    out = []
    for kind, n in CASES:
        atoms = lc.alkane(n) if kind == "alkane" else lc.water_cluster(n)
        mol = lc.build_mol(atoms)
        aux = df.make_auxmol(mol, lc.AUX_JKFIT)
        t0 = time.perf_counter()
        J = aux.intor("int2c2e")
        t_int = time.perf_counter() - t0
        naux = J.shape[0]
        A = lc._ao_atom(aux)
        X = mol.atom_coords()
        Dat = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        tg = tmin(lambda: scipy.linalg.cho_factor(J, lower=True))
        rec = {"system": f"{kind}{n}", "natm": mol.natm, "nbf": mol.nao_nr(), "naux": naux, "t_global_s": tg,
               "t_int2c_s": t_int, "flop_global": naux ** 3 / 3.0, "local": {}}
        for R in RADII:
            blocks = [np.where(Dat[a][A] <= R)[0] for a in range(mol.natm)]
            subs = [np.ascontiguousarray(J[np.ix_(b, b)]) for b in blocks]  # extraction not timed
            tl = tmin(lambda: [scipy.linalg.cho_factor(M, lower=True) for M in subs])
            nd = np.array([len(b) for b in blocks], dtype=float)
            rec["local"][str(R)] = {"t_local_s": tl, "ndom_mean": float(nd.mean()), "ndom_max": int(nd.max()),
                                    "flop_local": float(np.sum(nd ** 3) / 3.0)}
        out.append(rec)
        l8 = rec["local"]["8"]
        print(f"{kind}{n:<4d} natm={mol.natm:4d} naux={naux:6d} global {tg * 1e3:9.1f} ms | R=8: ndom={l8['ndom_mean']:.0f}"
              f" local {l8['t_local_s'] * 1e3:9.1f} ms  flop ratio local/global={l8['flop_local'] / rec['flop_global']:.2f}",
              flush=True)
        lc.save_json("q4_timing.json", out)
        del J


if __name__ == "__main__":
    main()
