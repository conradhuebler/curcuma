#!/usr/bin/env python3
"""Transfer test: learned Schwarz operator with parameters shared across translation-equivalent subdomains,
trained on small molecules, applied unchanged to larger unseen molecules of the same family.

Sharing rule (sl_learn.subdomain_classes, share='shared'): two atom-centred subdomains share Q, K when their key
(element and relative position of every atom within R_sd, rounded to 1e-3 Bohr) is identical; then their local
metric blocks are identical (checked, printed as max |Jn_ii - Jn_jj|). Coarse level: Galerkin (fixed, recomputed
on each molecule), because the paper's learned coarse block Q0 K0^T is tied to the number of atoms of the
training molecule and cannot be transferred. Subdomains of a test molecule whose key was never seen get a zero
local block (count reported as 'unseen').
Comparators on each test molecule (computed on the test molecule itself, i.e. no transfer needed): exact,
classical Schwarz (asm, asm-pu + bal, ras-core + bal), global data-optimal at the Schwarz rank bound (an oracle:
it is fitted on the test molecule's own data), Jacobi.

Families: 'alkanes' (train C6H14 + C8H18, test C12, C16, C20), 'acenes' (train naphthalene + anthracene, test
tetracene, pentacene). Water clusters have no equivalent subdomains (every key unique), so no transfer test.

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/schwarz_learn/run_transfer.py <family> [steps] [seeds]
Writes transfer_<family>.json. Claude Generated (Oct 2026).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sl_common as sc
import sl_learn as sl
import lmp2_common as lc
import numpy as np

FAMILIES = {"alkanes": (["alkane6", "alkane8"], ["alkane12", "alkane16", "alkane20"]),
            "acenes": (["acene2", "acene3"], ["acene4", "acene5"])}
CONFIGS = [(4.0, 8), (4.0, 32), (6.0, 8), (6.0, 32)]
LR = 1e-3


def main():
    fam = sys.argv[1]
    steps = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    nseed = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    train_names, test_names = FAMILIES[fam]
    trd = [sc.load(n) for n in train_names]
    rng = np.random.default_rng(0)
    pools, tr_cols, va_cols = [], [], []
    for d in trd:
        V = sc.rhs_pool(d, rng)
        m = V.shape[1]
        perm = rng.permutation(np.arange(1 + d.natm, m))
        nval = len(perm) // 5
        pools.append(V)
        tr_cols.append(np.concatenate([np.arange(1 + d.natm), perm[nval:]]))
        va_cols.append(perm[:nval])
    nds = [sl.NormData(d, V) for d, V in zip(trd, pools)]
    out = {"family": fam, "train": train_names, "test": test_names, "steps": steps, "lr": LR, "rows": [],
           "classes": {}}

    def add(d, op, extra, mp2=True):
        r = sc.evaluate(d, op, Vval=d.Veval, mp2=mp2)
        r.update(extra)
        r["system"] = d.name
        r["natm"] = d.natm
        r["naux"] = d.naux
        out["rows"].append(r)
        print(f"{d.name:9s} {r['label'][:58]:58s} p={r['params']:8d} val={r['val_rel_jnorm']:.2e} "
              f"EJ={r['EJ_robust_err']:+.2e} "
              f"MP2={r.get('Emp2_robust_err', 0):+.2e} jw={r['jw_identity_err']:.2e} gmres={r['gmres_iter']}",
              flush=True)
        lc.save_json(f"transfer_{fam}.json", out)

    models = {}
    for R, r_l in CONFIGS:
        dcs = [sc.decomposition(d, R) for d in trd]
        cls, keymap, sizes = sl.subdomain_classes(dcs, "shared")
        dev, refblocks = sl.check_equivalence(trd, dcs, cls)
        out["classes"][f"{R}_{r_l}"] = {"n_classes": len(sizes), "n_subdomains": sum(len(c) for c in cls),
                                        "max_dev_equal_keys": dev}
        print(f"# R={R} r_l={r_l}: {len(sizes)} classes for {sum(len(c) for c in cls)} subdomains, "
              f"max |Jn_ii - Jn_jj| within classes {dev:.1e}", flush=True)
        items = [(nd, (dc, c, sl.galerkin_n(d, dc))) for nd, dc, c, d in zip(nds, dcs, cls, trd)]
        for seed in range(nseed):
            mdl = sl.SchwarzModel(sizes, r_l, coarse="galerkin", seed=seed)
            h = sl.train(mdl, items, tr_cols, steps=steps, lr=LR, seed=seed, val_cols=va_cols, log_every=250)
            models[(R, r_l, seed)] = (mdl, keymap, h, refblocks)
            print(f"#   seed {seed}: final loss {h['loss'][-1]:.3e} val {h['val'][-1]:.3e} ({h['time_s']:.0f}s)",
                  flush=True)
    # evaluation on training and test molecules
    for name in train_names + test_names:
        d = trd[train_names.index(name)] if name in train_names else sc.load(name)
        role = "train" if name in train_names else "test"
        # evaluation pool: (P|rho), (P|rho_A), 2000 random (P|ia) (same seed for every operator)
        d.Veval = sc.rhs_pool(d, np.random.default_rng(7), n_mp2=2000)
        add(d, sc.exact_op(d), {"kind": "exact", "role": role}, mp2=False)
        add(d, sc.reference_preconditioners(d)[1], {"kind": "reference", "role": role}, mp2=False)
        for R in sorted({c[0] for c in CONFIGS}):
            dc = sc.decomposition(d, R)
            for var, co in (("asm", None), ("asm-pu", "bal"), ("ras-core", "bal")):
                add(d, sc.SchwarzExactOp(d, dc, var, co), {"kind": "schwarz-exact", "R": R, "variant": var,
                                                          "coarse": co, "role": role})
        for R, r_l in CONFIGS:
            dc = sc.decomposition(d, R)
            P0n = sl.galerkin_n(d, dc)
            k = min(d.naux, d.natm * r_l + dc.n0)
            V = sc.rhs_pool(d, np.random.default_rng(1))
            ops = sc.global_lowrank_ops(d, k, V)
            add(d, ops[2], {"kind": "global-oracle", "rank": k, "R": R, "r_l": r_l, "role": role})
            for seed in range(nseed):
                mdl, keymap, h, refblocks = models[(R, r_l, seed)]
                cls = np.array([keymap.get(kk, -1) for kk in dc.keys])
                unseen = int(np.sum(cls < 0))
                # equivalence check against the training blocks of the same class (copy: do not extend refs)
                dev_t, _ = sl.check_equivalence([d], [dc], [cls], dict(refblocks))
                Gn = mdl.dense((dc, cls, P0n), d.naux)
                flops = int(sum(4 * len(idx) * r_l for idx in dc.idx) + 4 * d.naux + 4 * dc.n0 ** 2)
                add(d, sc.DenseOp(Gn / d.sc[:, None] / d.sc[None, :], mdl.n_params(), flops,
                                  f"learned shared Schwarz r={r_l} + Galerkin (R={R}) seed {seed}"),
                    {"kind": "schwarz-learned-shared", "R": R, "r_l": r_l, "seed": seed, "unseen": unseen,
                     "max_dev_vs_train_blocks": dev_t,
                     "role": role, "train_time_s": h["time_s"], "final_loss": h["loss"][-1],
                     "final_val": h["val"][-1]}, mp2=(seed == 0))
        if role == "test":
            d.Araw = None


if __name__ == "__main__":
    main()
