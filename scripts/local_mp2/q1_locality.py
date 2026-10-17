#!/usr/bin/env python3
"""Q1 of the local-MP2 feasibility study: locality of the localised occupied orbitals and decay of the exact pair
energies.

Exact pair energies in a localised basis: canonical RI-MP2 amplitudes T_{I'J'} (full virtual space) are rotated to the
localised occupied orbitals, T_ij = sum U_{I'i} U_{J'j} T_{I'J'}, same for K; e_ij = (2-delta_ij) sum_ab T_ij^ab
(2K_ij^ab - K_ij^ba) for i<=j. sum e_ij = E(RI-MP2) exactly (invariance), so the cutoff curves measure only the pair
truncation, nothing else. Pair distance = distance between orbital centroids <i|r|i>.

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/q1_locality.py [--loc pm,boys,ibo] [systems...]
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
RCUT = [6, 8, 10, 12, 15, 20, 25]  # Bohr
ETHR = [1e-4, 1e-5, 1e-6, 1e-7]
BINS = [0, 2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 40, 60]


def exact_pair_energies(s, Bc):
    o, v = s.nval, s.nvir
    eo, ev = s.eval_, s.evir
    Tc = np.empty((o, o, v, v))
    for i in range(o):
        K = np.einsum("aQ,jbQ->jab", Bc[i], Bc, optimize=True)
        Tc[i] = -K / (ev[None, :, None] + ev[None, None, :] - eo[i] - eo[:, None, None])
    U = s.U
    Tl = (U.T @ Tc.reshape(o, -1)).reshape(o, o, v, v).transpose(1, 0, 2, 3).reshape(o, -1)
    del Tc
    Tl = (U.T @ Tl).reshape(o, o, v, v).transpose(1, 0, 2, 3)
    Bl = np.einsum("Ii,IaQ->iaQ", U, Bc, optimize=True)
    E = np.zeros((o, o))
    for i in range(o):
        K = np.einsum("aQ,jbQ->jab", Bl[i], Bl, optimize=True)
        E[i] = np.einsum("jab,jab->j", Tl[i], 2 * K - K.transpose(0, 2, 1))
    del Tl
    return E  # E[i,j] + E[j,i] = pair energy for i != j


def bp_sizes(s, thr=0.02):
    return [len(lc.bp_domain(s, i, thr)) for i in range(s.nval)]


def main():
    args = sys.argv[1:]
    locs = ["pm", "boys", "ibo"]
    if "--loc" in args:
        k = args.index("--loc")
        locs = args[k + 1].split(",")
        args = args[:k] + args[k + 2:]
    names = args or DEFAULT
    out = []
    for name in names:
        s = lc.setup(name, loc=None)
        t0 = time.time()
        lc.localise(s, "can")
        Bc = lc.b_canonical(s, lc.b_local_ao(s))  # canonical B (independent of the localisation)
        e_ref = lc.rimp2_canonical(s, Bc)
        for loc in locs:
            t1 = time.time()
            lc.localise(s, loc)
            E = exact_pair_energies(s, Bc)
            o = s.nval
            iu = np.triu_indices(o)
            pe = np.where(iu[0] == iu[1], E[iu], E[iu] + E.T[iu])
            dist = np.linalg.norm(s.centroid[iu[0]] - s.centroid[iu[1]], axis=1)
            tot = pe.sum()
            rec = {"system": name, "class": lc.SYSTEMS[name][0], "loc": loc, "natm": s.natm, "nbf": s.n,
                   "nval": o, "npairs": len(pe), "e_rimp2": e_ref, "sum_pairs": float(tot),
                   "mean_spread_bohr": float(np.mean(s.spread)), "max_spread_bohr": float(np.max(s.spread)),
                   "bp002_mean": float(np.mean(bp_sizes(s))), "bp002_max": int(np.max(bp_sizes(s)))}
            rec["cut_err"] = {str(r): float(pe[dist > r].sum()) for r in RCUT}
            rec["cut_npairs"] = {str(r): int(np.sum(dist <= r)) for r in RCUT}
            rec["nsig"] = {f"{t:.0e}": int(np.sum(np.abs(pe) > t)) for t in ETHR}
            # error if all pairs below |e| < t are dropped (oracle selection, lower bound for any screening)
            rec["sig_err"] = {f"{t:.0e}": float(pe[np.abs(pe) <= t].sum()) for t in ETHR}
            dec = []
            for a, b in zip(BINS[:-1], BINS[1:]):
                m = (dist >= a) & (dist < b)
                if m.any():
                    dec.append([a, b, int(m.sum()), float(np.max(np.abs(pe[m]))), float(np.sum(pe[m]))])
            rec["decay"] = dec
            far = (dist > 8) & (np.abs(pe) > 1e-12)
            if far.sum() > 5:
                rec["slope_logR"] = float(np.polyfit(np.log(dist[far]), np.log(np.abs(pe[far])), 1)[0])
            rec["time_s"] = time.time() - t1
            out.append(rec)
            print(f"{name:9s} {loc:4s} o={o:3d} pairs={len(pe):5d} sum-ref={tot - e_ref:.1e} spread={rec['mean_spread_bohr']:.2f}"
                  f" BP0.02={rec['bp002_mean']:.1f} nsig(1e-6)={rec['nsig']['1e-06']} err(R>10)={rec['cut_err']['10']:.2e}"
                  f" slope={rec.get('slope_logR', float('nan')):.2f} {rec['time_s']:.1f}s", flush=True)
        del Bc
    lc.save_json("q1_locality.json", out)


if __name__ == "__main__":
    main()
