#!/usr/bin/env python3
"""Per-molecule comparison (no transfer): learned two-level Schwarz operator (free parameters per subdomain, as in
Koehler/Rheinbach) vs non-learned comparators for G ~ J^-1 (def2-universal-JKFIT), on one molecule at a time.

Comparators (sl_common / sl_learn):
  exact J^-1; global rank k: truncated eigen J^-1, largest-eigen J, data-optimal (reduced-rank regression);
  classical Schwarz with exact local inverses (asm, asm-pu, ras-core; coarse none/add/bal);
  closed-form local reduced-rank blocks of rank r_l + Galerkin coarse (same structure and budget as the learned one);
  learned Schwarz (coarse learned or Galerkin), learned global Q K^T (rank-matched and parameter-matched);
  identity / Jacobi as preconditioner references.
Metrics: operator errors, validation loss, robust/plain RI-J and RI-MP2 energy errors with c = G v,
PCG (symmetrised G) and GMRES iterations for J x = (P|rho).
Data: pool = (P|rho), (P|rho_A), all (P|ia); 20 % of the (P|ia) columns held out for validation.

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/schwarz_learn/run_free.py <system> [steps] [seeds]
Writes free_<system>.json. Claude Generated (Oct 2026).
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sl_common as sc
import sl_learn as sl
import lmp2_common as lc
import numpy as np

CONFIGS = [(4.0, 8), (4.0, 32), (6.0, 8), (6.0, 32)]  # (R_sd in Bohr, local rank r_l)
LR = 1e-3  # the paper's default; the lr dependence is in run_lrsweep.py


def split(d, V, rng):
    m = V.shape[1]
    perm = rng.permutation(np.arange(1 + d.natm, m))
    nval = len(perm) // 5
    return np.concatenate([np.arange(1 + d.natm), perm[nval:]]), perm[:nval]


def main():
    name = sys.argv[1]
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    nseed = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    d = sc.load(name)
    rng = np.random.default_rng(0)
    V = sc.rhs_pool(d, rng)
    tr, va = split(d, V, rng)
    nd = sl.NormData(d, V)
    Vval = V[:, va]
    out = {"system": name, "natm": d.natm, "naux": d.naux, "E_J": d.E_J, "E_mp2": d.E_mp2, "ntrain": len(tr),
           "nval": len(va), "steps": steps, "lr": LR, "rows": []}
    t = time.perf_counter()
    for _ in range(3):
        t = time.perf_counter()
        sc.scipy.linalg.cho_factor(d.J, lower=True)
    out["t_cholesky_s"] = time.perf_counter() - t

    def add(op, extra=None, mp2=True):
        r = sc.evaluate(d, op, Vval=Vval, mp2=mp2)
        r.update(extra or {})
        out["rows"].append(r)
        print(f"{r['label'][:60]:60s} p={r['params']:8d} val={r.get('val_rel_jnorm', 0):.2e} "
              f"EJ={r['EJ_robust_err']:+.2e} MP2={r.get('Emp2_robust_err', 0):+.2e} frobN={r['frob_rel_norm']:.2e} "
              f"pcg={r['pcg_iter']} gmres={r['gmres_iter']}", flush=True)
        lc.save_json(f"free_{name}.json", out)

    add(sc.exact_op(d), {"kind": "exact"})
    for op in sc.reference_preconditioners(d):
        add(op, {"kind": "reference"}, mp2=False)
    out["local_density_fit_R8_EJ_err"] = sc.local_density_fit_error(d, 8.0)
    ranks = sorted({d.natm * r + d.natm for _, r in CONFIGS})
    for k in ranks:
        for op in sc.global_lowrank_ops(d, k, V[:, tr]):
            add(op, {"kind": "global-closed", "rank": k})
    for R in sorted({c[0] for c in CONFIGS}):
        dc = sc.decomposition(d, R)
        out.setdefault("subdomain_sizes", {})[str(R)] = [len(i) for i in dc.idx]
        for var in ("asm", "asm-pu", "ras-core"):
            for co in (None, "add", "bal"):
                add(sc.SchwarzExactOp(d, dc, var, co), {"kind": "schwarz-exact", "R": R, "variant": var,
                                                       "coarse": co})
    for R, r_l in CONFIGS:
        dc = sc.decomposition(d, R)
        cls, _, sizes = sl.subdomain_classes([dc], "free")
        P0n = sl.galerkin_n(d, dc)
        ctx = (dc, cls[0], P0n)
        # closed-form local reduced-rank blocks + Galerkin coarse (non-learned, same structure/budget)
        m = sl.SchwarzModel(sizes, r_l, coarse="galerkin", seed=0)
        blocks = sl.local_rrr_blocks(d, dc, V[:, tr], r_l)
        m.Q = [b[0] for b in blocks]
        m.K = [b[1] for b in blocks]
        Gn = m.dense(ctx, d.naux)
        add(sc.DenseOp(Gn / d.sc[:, None] / d.sc[None, :], m.n_params(), 0, f"local RRR r={r_l} + Galerkin (R={R})"),
            {"kind": "schwarz-rrr", "R": R, "r_l": r_l})
        for coarse in ("learned", "galerkin"):
            for seed in range(nseed):
                m = sl.SchwarzModel(sizes, r_l, coarse=coarse, n0=dc.n0, seed=seed)
                h = sl.train(m, [(nd, ctx)], [tr], steps=steps, lr=LR, seed=seed, val_cols=[va], log_every=250)
                Gn = m.dense(ctx, d.naux)
                flops = int(sum(4 * n * r_l for n in sizes) + (4 * d.naux + 4 * dc.n0 * dc.n0))
                add(sc.DenseOp(Gn / d.sc[:, None] / d.sc[None, :], m.n_params(), flops,
                               f"learned Schwarz r={r_l} coarse={coarse} (R={R}) seed {seed}"),
                    {"kind": "schwarz-learned", "R": R, "r_l": r_l, "coarse": coarse, "seed": seed, "hist": h,
                     "train_time_s": h["time_s"]})
        # global learned: rank-matched (rank bound natm*r_l + n0) and parameter-matched
        p_s = m.n_params()
        for tag, rg in (("rank-matched", d.natm * r_l + dc.n0), ("param-matched", max(1, round(p_s / (2 * d.naux))))):
            done = out.setdefault("_done_global", [])
            if [tag, rg] in done:
                continue
            done.append([tag, rg])
            for seed in range(nseed):
                g = sl.GlobalModel(d.naux, rg, seed=seed)
                h = sl.train(g, [(nd, None)], [tr], steps=steps, lr=LR, seed=seed, val_cols=[va], log_every=250)
                Gn = g.dense(None, d.naux)
                add(sc.DenseOp(Gn / d.sc[:, None] / d.sc[None, :], g.n_params(), 4 * d.naux * rg,
                               f"learned global rank {rg} ({tag}) seed {seed}"),
                    {"kind": "global-learned", "rank": rg, "match": tag, "seed": seed, "hist": h,
                     "train_time_s": h["time_s"]})
            if tag == "param-matched":  # closed-form optimum at the same rank for reference
                for op in sc.global_lowrank_ops(d, rg, V[:, tr]):
                    add(op, {"kind": "global-closed", "rank": rg})


if __name__ == "__main__":
    main()
