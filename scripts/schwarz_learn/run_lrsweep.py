#!/usr/bin/env python3
"""Learning-rate sweep (analogue of Koehler/Rheinbach Table 1) for G ~ J^-1 of C4H10 / def2-universal-JKFIT:
learned two-level Schwarz (free per subdomain, R_sd = 4 Bohr, r_l = 8, learned coarse block of rank n0 = natm)
vs learned global Q K^T at the Schwarz rank bound (natm*r_l + n0, 'rank-matched') and at the same parameter
count ('param-matched'). 2000 Adam steps (as in the paper), lr in {1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2},
3 seeds each. Reported per run: final running training loss, validation loss, normalised operator error, and the
closed-form data-optimal loss of a global operator of the same rank (the optimum the global model could reach).

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/schwarz_learn/run_lrsweep.py [system] [steps]
Writes lrsweep_<system>.json. Claude Generated (Oct 2026).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sl_common as sc
import sl_learn as sl
import lmp2_common as lc
import numpy as np

LRS = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "alkane4"
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    d = sc.load(name)
    rng = np.random.default_rng(0)
    V = sc.rhs_pool(d, rng)
    m = V.shape[1]
    perm = rng.permutation(np.arange(1 + d.natm, m))
    nval = len(perm) // 5
    tr = np.concatenate([np.arange(1 + d.natm), perm[nval:]])
    va = perm[:nval]
    nd = sl.NormData(d, V)
    dc = sc.decomposition(d, 4.0)
    cls, _, sizes = sl.subdomain_classes([dc], "free")
    ctx = (dc, cls[0], sl.galerkin_n(d, dc))
    r_l = 8
    p_s = sl.SchwarzModel(sizes, r_l, coarse="learned", n0=dc.n0).n_params()
    ranks = {"rank-matched": d.natm * r_l + dc.n0, "param-matched": max(1, round(p_s / (2 * d.naux)))}
    out = {"system": name, "naux": d.naux, "steps": steps, "schwarz_params": p_s, "ranks": ranks, "runs": [],
           "optimum": {}}
    Jni = d.Jinv * d.sc[:, None] * d.sc[None, :]
    for tag, k in ranks.items():
        op = sc.global_lowrank_ops(d, k, V[:, tr])[2]
        out["optimum"][tag] = {"rank": k, "train": sc.rel_jnorm_err(d, op.G, V[:, tr]),
                               "val": sc.rel_jnorm_err(d, op.G, V[:, va])}
    for lr in LRS:
        for seed in range(3):
            models = [("schwarz", sl.SchwarzModel(sizes, r_l, coarse="learned", n0=dc.n0, seed=seed), ctx)]
            for tag, k in ranks.items():
                models.append((f"global {tag}", sl.GlobalModel(d.naux, k, seed=seed), None))
            for label, mdl, cx in models:
                h = sl.train(mdl, [(nd, cx)], [tr], steps=steps, lr=lr, seed=seed, val_cols=[va], log_every=100)
                Gn = mdl.dense(cx, d.naux)
                rec = {"model": label, "lr": lr, "seed": seed, "params": mdl.n_params(),
                       "train": h["loss"][-1] if h["loss"] else float("nan"),
                       "val": h["val"][-1] if h["val"] else float("nan"),
                       "diverged": "diverged_step" in h, "time_s": h["time_s"],
                       "frob_rel_norm": float(np.linalg.norm(Gn - Jni) / np.linalg.norm(Jni)), "hist": h}
                out["runs"].append(rec)
                print(f"lr={lr:.0e} seed={seed} {label:22s} p={rec['params']:7d} train={rec['train']:.3e} "
                      f"val={rec['val']:.3e} div={rec['diverged']} ({h['time_s']:.0f}s)", flush=True)
                lc.save_json(f"lrsweep_{name}.json", out)


if __name__ == "__main__":
    main()
