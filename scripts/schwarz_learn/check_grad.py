#!/usr/bin/env python3
"""Finite-difference check of the analytic gradients in sl_learn (Schwarz free/shared, learned/Galerkin coarse,
global model) on C4H10 / JKFIT. Prints max relative deviation per model; expected <~1e-6 at h=1e-3 (smaller h is round-off dominated: 1e-6 gives ~1e-4).

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/schwarz_learn/check_grad.py
Claude Generated (Oct 2026).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sl_common as sc
import sl_learn as sl
import numpy as np


def fd_check(model, items, cols, nprobe=12, h=1e-3, seed=0):
    rng = np.random.default_rng(seed)
    L0, g = sl.loss_and_grad(model, items, cols)
    params = model.params()
    worst = 0.0
    for _ in range(nprobe):
        k = rng.integers(len(params))
        p = params[k]
        j = rng.integers(p.size)
        old = p.flat[j]
        p.flat[j] = old + h
        Lp = sl.loss_and_grad(model, items, cols, need_grad=False)[0]
        p.flat[j] = old - h
        Lm = sl.loss_and_grad(model, items, cols, need_grad=False)[0]
        p.flat[j] = old
        fd = (Lp - Lm) / (2 * h)
        an = g[k].flat[j]
        worst = max(worst, abs(fd - an) / max(abs(fd), abs(an), 1e-12))
    return L0, worst


def main():
    d = sc.load("alkane4")
    rng = np.random.default_rng(0)
    V = sc.rhs_pool(d, rng, n_mp2=200)
    nd = sl.NormData(d, V)
    cols = [np.arange(64)]
    for R in (4.0,):
        dc = sc.decomposition(d, R)
        for share in ("free", "shared"):
            cls, _, sizes = sl.subdomain_classes([dc], share)
            for coarse in ("learned", "galerkin"):
                m = sl.SchwarzModel(sizes, 4, coarse=coarse, n0=dc.n0, seed=1)
                ctx = (dc, cls[0], sl.galerkin_n(d, dc))
                L0, w = fd_check(m, [(nd, ctx)], cols)
                print(f"Schwarz {share:6s} coarse={coarse:8s} classes={len(sizes):3d} L={L0:.4e} max rel FD dev {w:.2e}")
    g = sl.GlobalModel(d.naux, 8, seed=1)
    L0, w = fd_check(g, [(nd, None)], cols)
    print(f"Global rank 8 L={L0:.4e} max rel FD dev {w:.2e}")


if __name__ == "__main__":
    main()
