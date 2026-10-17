#!/usr/bin/env python3
"""Q2 of the local-MP2 feasibility study: Boughton-Pulay orbital domains + PAO pair domains (Pulay-Saeboe LMP2).

Per localised (Pipek-Mezey) orbital i the BP domain [i] (lmp2_common.bp_domain, criterion thr), optionally extended by
all atoms within r_ext Bohr; pair domain [ij] = [i] u [j] (pair domains overlap between pairs by construction). The
pair virtual space is spanned by the PAOs of [ij] (redundancy threshold 1e-6 on the normalised PAO overlap), the
coupled LMP2 equations are iterated to |dE| < 1e-9 Eh and max|R| < 1e-7 (DIIS). All pairs are kept (no pair
screening), so the error is the domain error alone. Reference: full RI-MP2 (frozen core) with the same aux basis.

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/q2_domains.py [--cfg bp0.02,bp0.02+r3] [systems...]
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
CFGS = ["bp0.02", "bp0.01", "bp0.001", "bp0.02+r3", "bp0.02+r5", "bp0.02+r7"]


def parse_cfg(c):
    thr, _, r = c[2:].partition("+r")
    return float(thr), float(r) if r else 0.0


def main():
    args = sys.argv[1:]
    cfgs = CFGS
    if "--cfg" in args:
        k = args.index("--cfg")
        cfgs = args[k + 1].split(",")
        args = args[:k] + args[k + 2:]
    names = args or DEFAULT
    out = []
    for name in names:
        s = lc.setup(name, loc="pm")
        Bao = lc.b_local_ao(s)
        e_ref = lc.rimp2_canonical(s, lc.b_canonical(s, Bao))
        pairs = lc.pair_list(s)
        for cfg in cfgs:
            t0 = time.time()
            thr, rext = parse_cfg(cfg)
            dom = [lc.extend_domain(s, lc.bp_domain(s, i, thr), rext / 1.0) for i in range(s.nval)]
            spaces, Ks, natp = [], [], []
            for i, j in pairs:
                at = sorted(set(dom[i]) | set(dom[j]))
                natp.append(len(at))
                sp = lc.pair_space(s, at)
                spaces.append(sp)
                Ks.append(lc.pair_ints(Bao, sp[0], i, j))
            t1 = time.time()
            semi = lc.solve_lmp2(s, pairs, spaces, Ks, coupled=False)["energy"]
            res = lc.solve_lmp2(s, pairs, spaces, Ks, coupled=True)
            m = np.array([sp[0].shape[1] for sp in spaces])
            rec = {"system": name, "class": lc.SYSTEMS[name][0], "cfg": cfg, "natm": s.natm, "nbf": s.n,
                   "nvir": s.nvir, "nval": s.nval, "e_rimp2": e_ref, "e_lmp2": res["energy"], "e_semican": semi,
                   "err": res["energy"] - e_ref, "err_per_atom": (res["energy"] - e_ref) / s.natm,
                   "recovery_pct": 100 * res["energy"] / e_ref, "iterations": res["iterations"],
                   "dom_atoms_mean": float(np.mean([len(d) for d in dom])),
                   "dom_atoms_max": int(max(len(d) for d in dom)),
                   "pair_atoms_mean": float(np.mean(natp)), "m_mean": float(m.mean()), "m_max": int(m.max()),
                   "t_build_s": t1 - t0, "t_solve_s": time.time() - t1}
            out.append(rec)
            print(f"{name:9s} {cfg:10s} E={res['energy']:.8f} err={rec['err']:.3e} ({rec['err_per_atom']:.1e}/atom,"
                  f" {rec['recovery_pct']:.3f}%) it={res['iterations']} [i]={rec['dom_atoms_mean']:.1f}"
                  f" [ij]={rec['pair_atoms_mean']:.1f} m={rec['m_mean']:.0f}/{s.nvir}"
                  f" build {rec['t_build_s']:.0f}s solve {rec['t_solve_s']:.0f}s", flush=True)
            lc.save_json("q2_domains.json", out)
            del spaces, Ks, res


if __name__ == "__main__":
    main()
