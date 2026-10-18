#!/usr/bin/env python3
"""Learned operators for the Schwarz experiment (docs/SCHWARZ_LEARNED_OPERATOR_2026-10.md): models, loss, Adam.

Models (all act in Coulomb-normalised variables, Gn ~ Jn^-1, see sl_common):
  SchwarzModel  Koehler/Rheinbach eq. 17: Gn = Phi Q0 K0^T Phi^T + sum_i R_i^T D_i^1/2 Q_i K_i^T D_i^1/2 R_i
                subdomains/PU/coarse space from sl_common.decomposition (atom-centred, multipole coarse basis).
                share='free'   : own Q_i, K_i per subdomain (as in the paper, not transferable)
                share='shared' : one Q, K per class of translation-equivalent subdomains (identical key, i.e.
                                 identical element pattern and relative positions within R_sd); transferable to
                                 another molecule for every subdomain whose key was seen in training
                coarse='learned' (Q0 K0^T, rank r0, as in the paper) | 'galerkin' (fixed Phi (Phi^T Jn Phi)^-1
                Phi^T, not trained, recomputed on every molecule -> transferable) | 'none'
  GlobalModel   Gn = Q K^T, Q, K in R^(naux x r) (the paper's global low-rank attention baseline)

Loss (analogue of the paper's sample-weighted wMSE; here the quantity the robust RI energy depends on):
  L = mean_b  e_b^T Jn e_b / (vn_b^T cn_b),   e_b = Gn vn_b - cn_b,   cn_b = Jn^-1 vn_b
  = mean relative J-norm fitting error (2 x relative robust RI-J energy error per right-hand side).
Optimiser: Adam (beta 0.9/0.999, eps 1e-8, no weight decay, as AdamW(wd=0) in the paper) with
ReduceLROnPlateau (factor 0.5, patience 200 steps on the 50-step running mean of the minibatch loss, relative
threshold 1e-4), minibatch 256 columns drawn from the training pool. Hand-written (no torch in the venv).
Initialisation: N(0,1) * alpha with alpha^2 sqrt(r) = 0.01/sqrt(n_block) (deviation from the paper's h-based
scale, which has no analogue here), or a warm start from the closed-form local reduced-rank blocks.

Claude Generated (Oct 2026). AI-generated, machine-tested only.
"""
import time

import numpy as np

import sl_common as sc


class NormData:
    """Normalised quantities of one system for training: Jn, vn = v/sc, cn = Jn^-1 vn, weights."""

    def __init__(self, d, V):
        self.d = d
        self.Jn = d.Jn
        self.Vn = V / d.sc[:, None]
        Jni = d.Jinv * d.sc[:, None] * d.sc[None, :]
        self.Cn = Jni @ self.Vn
        self.w = 1.0 / np.einsum("pm,pm->m", self.Vn, self.Cn)


def subdomain_classes(dcs, share):
    """dcs: list of decompositions (one per training molecule). Returns (list of class-id arrays per molecule,
    dict key -> class id, block sizes per class)."""
    keymap, sizes, cls = {}, [], []
    for dc in dcs:
        c = []
        for i, k in enumerate(dc.keys):
            kk = k if share == "shared" else (id(dc), i)
            if kk not in keymap:
                keymap[kk] = len(sizes)
                sizes.append(len(dc.idx[i]))
            c.append(keymap[kk])
        cls.append(np.array(c))
    return cls, keymap, sizes


def check_equivalence(ds, dcs, cls, ref=None):
    """Max |J_ii - J_jj| over subdomains in the same class (validates the translation-equivalence keys).
    ref: dict class -> reference block (filled on the way; pass the training dict to check a test molecule)."""
    ref = {} if ref is None else ref
    worst = 0.0
    for d, dc, c in zip(ds, dcs, cls):
        for i, idx in enumerate(dc.idx):
            Ji = d.Jn[np.ix_(idx, idx)]
            if c[i] in ref:
                worst = max(worst, float(np.max(np.abs(Ji - ref[c[i]]))))
            elif c[i] >= 0:
                ref[c[i]] = Ji
    return worst, ref


class SchwarzModel:
    def __init__(self, sizes, r_l, coarse="galerkin", r0=None, n0=None, seed=0, init="random"):
        rng = np.random.default_rng(seed)
        self.r_l, self.coarse = r_l, coarse
        self.Q, self.K = [], []
        for n in sizes:
            a = np.sqrt(0.01 / np.sqrt(n) / np.sqrt(r_l))
            self.Q.append(a * rng.standard_normal((n, r_l)))
            self.K.append(a * rng.standard_normal((n, r_l)))
        self.Q0 = self.K0 = None
        if coarse == "learned":
            self.r0 = r0 or n0
            a = np.sqrt(0.01 / np.sqrt(n0) / np.sqrt(self.r0))
            self.Q0 = a * rng.standard_normal((n0, self.r0))
            self.K0 = a * rng.standard_normal((n0, self.r0))

    def params(self):
        p = self.Q + self.K
        if self.coarse == "learned":
            p = p + [self.Q0, self.K0]
        return p

    def n_params(self):
        return int(sum(x.size for x in self.params()))

    # ---- one molecule: ctx = (dc, cls, P0n) ----
    def forward(self, ctx, Vn):
        dc, cls, P0n = ctx
        Y = np.zeros_like(Vn)
        for i, idx in enumerate(dc.idx):
            c = cls[i]
            if c < 0:
                continue
            Vi = dc.sq[i][:, None] * Vn[idx]
            Y[idx] += dc.sq[i][:, None] * (self.Q[c] @ (self.K[c].T @ Vi))
        if self.coarse == "learned":
            Y += dc.Phi @ (self.Q0 @ (self.K0.T @ (dc.Phi.T @ Vn)))
        elif self.coarse == "galerkin":
            Y += P0n @ Vn
        return Y

    def backward(self, ctx, Vn, Gy, grads):
        """Accumulate dL/dparams into grads (same list layout as params()) given dL/dY = Gy."""
        dc, cls, _ = ctx
        nb = len(self.Q)
        for i, idx in enumerate(dc.idx):
            c = cls[i]
            if c < 0:
                continue
            Vi = dc.sq[i][:, None] * Vn[idx]
            Gi = dc.sq[i][:, None] * Gy[idx]
            grads[c] += Gi @ (Vi.T @ self.K[c])
            grads[nb + c] += Vi @ (Gi.T @ self.Q[c])
        if self.coarse == "learned":
            V0, G0 = dc.Phi.T @ Vn, dc.Phi.T @ Gy
            grads[2 * nb] += G0 @ (V0.T @ self.K0)
            grads[2 * nb + 1] += V0 @ (G0.T @ self.Q0)

    def dense(self, ctx, naux, fallback=None):
        """Dense Gn on one molecule. Subdomains with class -1 (key unseen in training) use the fallback block
        function fallback(i) (e.g. the exact local inverse with PU weighting) if given, else zero."""
        dc, cls, P0n = ctx
        G = np.zeros((naux, naux))
        for i, idx in enumerate(dc.idx):
            c = cls[i]
            if c < 0:
                B = fallback(i) if fallback else 0.0
            else:
                B = self.Q[c] @ self.K[c].T
            G[np.ix_(idx, idx)] += dc.sq[i][:, None] * B * dc.sq[i][None, :]
        if self.coarse == "learned":
            G += dc.Phi @ self.Q0 @ self.K0.T @ dc.Phi.T
        elif self.coarse == "galerkin":
            G += P0n
        return G


class GlobalModel:
    def __init__(self, naux, r, seed=0):
        rng = np.random.default_rng(seed)
        a = np.sqrt(0.01 / np.sqrt(naux) / np.sqrt(r))
        self.Q = a * rng.standard_normal((naux, r))
        self.K = a * rng.standard_normal((naux, r))

    def params(self):
        return [self.Q, self.K]

    def n_params(self):
        return self.Q.size + self.K.size

    def forward(self, ctx, Vn):
        return self.Q @ (self.K.T @ Vn)

    def backward(self, ctx, Vn, Gy, grads):
        grads[0] += Gy @ (Vn.T @ self.K)
        grads[1] += Vn @ (Gy.T @ self.Q)

    def dense(self, ctx, naux, fallback=None):
        return self.Q @ self.K.T


def galerkin_n(d, dc):
    """Galerkin coarse operator in normalised variables: Phi (Phi^T Jn Phi)^-1 Phi^T."""
    A0 = dc.Phi.T @ d.Jn @ dc.Phi
    return dc.Phi @ np.linalg.solve(A0, dc.Phi.T)


def loss_and_grad(model, items, idx_per_item, need_grad=True):
    """items: list of (NormData, ctx); idx_per_item: column indices per item for this batch."""
    params = model.params()
    grads = [np.zeros_like(p) for p in params] if need_grad else None
    tot, cnt = 0.0, 0
    for (nd, ctx), cols in zip(items, idx_per_item):
        if len(cols) == 0:
            continue
        Vn, Cn, w = nd.Vn[:, cols], nd.Cn[:, cols], nd.w[cols]
        E = model.forward(ctx, Vn) - Cn
        JE = nd.Jn @ E
        tot += float(np.sum(w * np.einsum("pm,pm->m", E, JE)))
        cnt += len(cols)
        if need_grad:
            model.backward(ctx, Vn, 2.0 * JE * w, grads)
    if need_grad:
        for g in grads:
            g /= cnt
    return tot / cnt, grads


def train(model, items, train_cols, steps=2000, lr=1e-3, batch=256, seed=0, log_every=100, patience=200,
          val_cols=None, verbose=False):
    """Adam + ReduceLROnPlateau. train_cols: list (per item) of training column indices. Returns history dict."""
    rng = np.random.default_rng(1000 + seed)
    params = model.params()
    m = [np.zeros_like(p) for p in params]
    v = [np.zeros_like(p) for p in params]
    b1, b2, eps = 0.9, 0.999, 1e-8
    sizes = np.array([len(c) for c in train_cols], dtype=float)
    hist = {"step": [], "loss": [], "lr": [], "val": []}
    window, best, bad = [], np.inf, 0
    t0 = time.time()
    for t in range(1, steps + 1):
        # draw the batch across molecules proportional to their pool sizes
        pick = rng.choice(len(train_cols), batch, p=sizes / sizes.sum())
        cols = [rng.choice(train_cols[k], np.sum(pick == k), replace=False) if np.sum(pick == k) else []
                for k in range(len(train_cols))]
        L, g = loss_and_grad(model, items, cols)
        if not np.isfinite(L):
            hist["diverged_step"] = t
            break
        for p, gg, mm, vv in zip(params, g, m, v):
            mm *= b1
            mm += (1 - b1) * gg
            vv *= b2
            vv += (1 - b2) * gg * gg
            p -= lr * (mm / (1 - b1 ** t)) / (np.sqrt(vv / (1 - b2 ** t)) + eps)
        window.append(L)
        if len(window) > 50:
            window.pop(0)
        run = np.mean(window)
        if run < best * (1 - 1e-4):
            best, bad = run, 0
        else:
            bad += 1
            if bad > patience:
                lr *= 0.5
                bad = 0
        if t % log_every == 0 or t == steps:
            hist["step"].append(t)
            hist["loss"].append(float(run))
            hist["lr"].append(lr)
            if val_cols is not None:
                hist["val"].append(loss_and_grad(model, items, val_cols, need_grad=False)[0])
            if verbose:
                print(f"   step {t:5d} loss {run:.3e} lr {lr:.1e}" +
                      (f" val {hist['val'][-1]:.3e}" if val_cols is not None else ""), flush=True)
    hist["time_s"] = time.time() - t0
    return hist


def local_rrr_blocks(d, dc, V, r):
    """Closed-form (non-learned) local blocks of rank r: for subdomain i, B_i = X X^T with X = Jn_ii^-1/2 P_r,
    P_r = top-r left singular vectors of Jn_ii^-1/2 (D_i^1/2 R_i Vn) diag(w)^1/2 (local reduced-rank fit of the
    local inverse to the local data). Returns list of (Q, K) = (X, X)."""
    nd = NormData(d, V)
    out = []
    for i, idx in enumerate(dc.idx):
        Ji = d.Jn[np.ix_(idx, idx)]
        w, U = np.linalg.eigh(Ji)
        Jmh = (U / np.sqrt(w)) @ U.T
        Wi = Jmh @ (dc.sq[i][:, None] * nd.Vn[idx]) * np.sqrt(nd.w)
        P, _, _ = np.linalg.svd(Wi, full_matrices=False)
        X = Jmh @ P[:, :r]
        out.append((X, X.copy()))
    return out
