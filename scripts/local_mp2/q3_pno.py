#!/usr/bin/env python3
"""Q3 of the local-MP2 feasibility study: PNO truncation.

Per pair ij (all pairs, Pipek-Mezey orbitals) a base virtual space is chosen:
  base=full   : all canonical virtuals (isolates the PNO error)
  base=bp<t>[+r<R>] : PAOs of the pair domain as in q2_domains.py (domain + PNO error together)
Semicanonical amplitudes in the base space give the pair density D_ij (lmp2_common.make_pnos), PNOs with occupation
< T_CutPNO are dropped, the coupled LMP2 equations are iterated in the PNO spaces (|dE| < 1e-9, max|R| < 1e-7).
Reported: error without and with the semicanonical truncation correction dE_PNO = sum_ij [e_ij^sc(base) - e_ij^sc(PNO)],
PNOs per pair, and a flop proxy (definitions in the result document) against the canonical RI-MP2 energy step.

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/q3_pno.py [--base full|bp0.001...]
       [--tcut 1e-6,1e-7,...] [systems...]
Claude Generated (Oct 2026).
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp2_common as lc  # first: sets the thread environment before numpy loads
import numpy as np


DEFAULT = ["alkane4", "alkane8", "alkane12", "alkane16", "alkane20", "water4", "water8", "water12", "water16",
           "acene1", "acene2", "acene3", "acene4", "acene5"]
TCUTS = [1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 1e-9]
FSCREEN = 1e-5  # only for the coupling-count part of the cost proxy, not applied to the equations


def cost_proxy(s, pairs, m, d):
    """Flop proxies (counts, not timings). m[p] = PNOs, d[p] = base-space dimension of pair p."""
    o, v, n, naux = s.nval, s.nvir, s.n, s.naux
    pidx = {pq: k for k, pq in enumerate(pairs)}
    can = o * (o + 1) / 2 * 2.0 * v * v * naux  # RI-MP2 energy step, i<=j
    kint = float(np.sum(2.0 * m * m * naux))  # K_ij in the PNO basis from PNO-transformed B
    btr = float(np.sum(2 * 2.0 * m * d * naux))  # B_i, B_j -> PNO basis per pair
    it_all = it_scr = ovl_all = 0.0
    F = np.abs(s.Focc)
    for p, (i, j) in enumerate(pairs):
        for k in range(o):
            for (a, b, f) in (((k, j), k != i, F[i, k]), ((i, k), k != j, F[k, j])):
                if not b:
                    continue
                q = pidx[tuple(sorted(a))]
                c = 2.0 * (m[p] * m[p] * m[q] + m[p] * m[q] * m[q])  # S T S^T with stored S^{pq}
                so = 2.0 * m[p] * n * m[q]  # building S^{pq} = V_p^T (S V_q)
                it_all += c
                ovl_all += so
                if f > FSCREEN:
                    it_scr += c
    return {"canonical": can, "k_pno": kint, "b_transform": btr, "iter_all_k": it_all, "iter_screened": it_scr,
            "overlaps_all_k": ovl_all}


def main():
    args = sys.argv[1:]
    base, tcuts = "full", TCUTS
    if "--base" in args:
        k = args.index("--base")
        base = args[k + 1]
        args = args[:k] + args[k + 2:]
    if "--tcut" in args:
        k = args.index("--tcut")
        tcuts = [float(x) for x in args[k + 1].split(",")]
        args = args[:k] + args[k + 2:]
    names = args or DEFAULT
    out = []
    fn = f"q3_pno_{base}.json"
    for name in names:
        s = lc.setup(name, loc="pm")
        Bao = lc.b_local_ao(s)
        e_ref = lc.rimp2_canonical(s, lc.b_canonical(s, Bao))
        pairs = lc.pair_list(s)
        t0 = time.time()
        if base == "full":
            sp_full = (s.Cvir, s.evir)
            Bv = np.einsum("imQ,ma->iaQ", Bao, s.Cvir, optimize=True)
            bases = [sp_full] * len(pairs)
            Kb = [Bv[i] @ Bv[j].T for i, j in pairs]
        else:
            thr, _, r = base[2:].partition("+r")
            dom = [lc.extend_domain(s, lc.bp_domain(s, i, float(thr)), float(r) if r else 0.0) for i in range(s.nval)]
            bases, Kb = [], []
            for i, j in pairs:
                sp = lc.pair_space(s, sorted(set(dom[i]) | set(dom[j])))
                bases.append(sp)
                Kb.append(lc.pair_ints(Bao, sp[0], i, j))
        d = np.array([b[0].shape[1] for b in bases], dtype=float)
        e_base_sc = lc.solve_lmp2(s, pairs, bases, Kb, coupled=False)["energy"]
        tb = time.time() - t0
        dist = np.array([np.linalg.norm(s.centroid[i] - s.centroid[j]) for i, j in pairs])
        Fd = np.diag(s.Focc)
        for tc in tcuts:
            t1 = time.time()
            spaces, Ks, dE, m = [], [], 0.0, []
            for p, (i, j) in enumerate(pairs):
                V, e, K, de, npno = lc.make_pnos(bases[p], Kb[p], Fd[i], Fd[j], i, j, tc)
                spaces.append((V, e))
                Ks.append(K)
                dE += de
                m.append(npno)
            m = np.array(m, dtype=float)
            res = lc.solve_lmp2(s, pairs, spaces, Ks, coupled=True)
            rec = {"system": name, "class": lc.SYSTEMS[name][0], "base": base, "tcut": tc, "natm": s.natm,
                   "nbf": s.n, "nval": s.nval, "nvir": s.nvir, "naux": s.naux, "npairs": len(pairs),
                   "e_rimp2": e_ref, "e_pno": res["energy"], "dE_pno": dE, "e_base_semican": e_base_sc,
                   "err": res["energy"] - e_ref, "err_corr": res["energy"] + dE - e_ref,
                   "err_corr_per_atom": (res["energy"] + dE - e_ref) / s.natm,
                   "recovery_pct": 100 * res["energy"] / e_ref, "recovery_corr_pct": 100 * (res["energy"] + dE) / e_ref,
                   "iterations": res["iterations"], "pno_mean": float(m.mean()), "pno_max": int(m.max()),
                   "pno_mean_close": float(m[dist < 10].mean()), "n_close": int(np.sum(dist < 10)),
                   "base_dim_mean": float(d.mean()), "zero_pno_pairs": int(np.sum(m == 0)),
                   "cost": cost_proxy(s, pairs, m, d), "t_base_s": tb, "t_s": time.time() - t1}
            out.append(rec)
            print(f"{name:9s} {base:8s} T={tc:.0e} err={rec['err']:.3e} corr={rec['err_corr']:.3e}"
                  f" ({rec['recovery_corr_pct']:.4f}%) it={rec['iterations']} PNO mean={rec['pno_mean']:.1f}"
                  f" close={rec['pno_mean_close']:.1f} max={rec['pno_max']} base={rec['base_dim_mean']:.0f}"
                  f" cost(iter/can)={rec['cost']['iter_all_k'] / rec['cost']['canonical']:.2e} {rec['t_s']:.0f}s",
                  flush=True)
            lc.save_json(fn, out)
        del bases, Kb


if __name__ == "__main__":
    main()
