#!/usr/bin/env python3
"""How much rank does any global low-rank G ~ J^-1 need for chemically useful fits? (capacity bound)

For the relative J-norm loss the best global rank-k operator on a data set is closed form (reduced-rank
regression, sl_common.global_lowrank_ops 'data-optimal'): G_k = J^-1/2 P_k P_k^T J^-1/2, P_k = top-k left singular
vectors of W = J^-1/2 V diag(w)^1/2. Two data sets: (a) 'ia': P_k fitted on 80 % of the (P|ia) columns, loss on the
held-out 20 %; (b) 'ao': P_k fitted on all AO-pair columns (P|mu nu), mu <= nu (the general rhs space, rho is a
linear combination of them), robust RI-J error for v = (P|rho). (P|rho) is an occupied-occupied quantity and is
outside the span of the dominant (P|ia) directions, so (a) is not used for the RI-J number.
Every learned or structured operator whose rank is bounded by k (the Schwarz structure has rank <= N r_l + r0) is
limited at least by this curve on the same data, up to the train/test difference.
Target: |dE_J| < 1e-6 Eh/atom (the vault falsifier for local fits).

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/schwarz_learn/rank_needed.py [systems...]
Writes rank_needed.json. Claude Generated (Oct 2026).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sl_common as sc
import lmp2_common as lc
import numpy as np

DEFAULT = ["alkane4", "alkane8", "alkane12", "water4", "water8", "acene2", "acene3"]
FRACS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 1.0]


def main():
    names = sys.argv[1:] or DEFAULT
    out = []
    for name in names:
        d = sc.load(name, keep_ao=True)
        A = d.Araw.reshape(-1, d.naux).T
        rng = np.random.default_rng(0)
        perm = rng.permutation(A.shape[1])
        nval = len(perm) // 5
        Vtr, Vva = A[:, perm[nval:]], A[:, perm[:nval]]
        w, U = np.linalg.eigh(d.J)
        Jmh = (U / np.sqrt(w)) @ U.T
        Wtr = Jmh @ Vtr
        Wtr = Wtr / np.linalg.norm(Wtr, axis=0)  # same as weighting with 1/(v^T J^-1 v)
        P, sv, _ = np.linalg.svd(Wtr, full_matrices=False)
        Wva = Jmh @ Vva
        Wva = Wva / np.linalg.norm(Wva, axis=0)
        wr = Jmh @ d.rho
        Wao = Jmh @ d.Vao
        Wao = Wao / np.maximum(np.linalg.norm(Wao, axis=0), 1e-300)
        Pa, _, _ = np.linalg.svd(Wao, full_matrices=False)
        del Wao, d.Vao
        rec = {"system": name, "natm": d.natm, "naux": d.naux, "E_J": d.E_J, "ntrain": Vtr.shape[1], "curve": []}
        kneed = None
        for f in FRACS:
            k = min(P.shape[1], int(round(f * d.naux)))
            Pk = P[:, :k]
            Rva = Wva - Pk @ (Pk.T @ Wva)
            val = float(np.mean(np.sum(Rva ** 2, axis=0)))
            Pak = Pa[:, :k]
            rr = wr - Pak @ (Pak.T @ wr)
            ej = float(-0.5 * rr @ rr)
            rec["curve"].append({"frac": f, "k": k, "val_rel_jnorm": val, "EJ_err": ej, "EJ_err_atom": ej / d.natm})
        # smallest k meeting the RI-J target, scanning all k
        Pr = Pa.T @ wr
        resid = wr @ wr - np.cumsum(Pr ** 2)
        ok = np.where(0.5 * np.maximum(resid, 0) / d.natm < 1e-6)[0]
        kneed = int(ok[0] + 1) if len(ok) else None
        rec["k_needed_EJ_1e-6_per_atom"] = kneed
        rec["rank_W"] = int(P.shape[1])
        out.append(rec)
        print(f"{name:9s} naux={d.naux} rank(data)={P.shape[1]} k needed for |dE_J|<1e-6/atom: {kneed} "
              f"({(kneed or 0) / d.naux:.2f} naux) | val@0.5={rec['curve'][4]['val_rel_jnorm']:.1e} "
              f"val@0.9={rec['curve'][8]['val_rel_jnorm']:.1e}", flush=True)
        lc.save_json("rank_needed.json", out)


if __name__ == "__main__":
    main()
