#!/usr/bin/env python3
"""Shared code of the learned two-level Schwarz operator experiment (docs/SCHWARZ_LEARNED_OPERATOR_2026-10.md).

Question: does the learned, softmax-free two-level overlapping Schwarz operator of Koehler and Rheinbach
(arXiv:2606.18525, eq. 17)
    M^-1 = Phi Q0 K0^T Phi^T + sum_i R_i^T D_i^1/2 Q_i K_i^T D_i^1/2 R_i
help to approximate the inverse Coulomb metric J^-1 (def2-universal-JKFIT) of RI, compared with non-learned
operators of the same form?

How an approximate inverse G is used (decision, see the doc): only as the map rhs -> fitting coefficients,
c = G v with v = (P|rho) (RI-J) or v = (P|ia) (RI-MP2), followed by the robust (Dunlap) energy expression
    E_J  = (rho|c) - 1/2 c^T J c                       error = -1/2 dc^T J dc,  dc = (G - J^-1) v
    (ia|jb) = A_ia.c_jb + c_ia.A_jb - c_ia^T J c_jb    error = -dc_ia^T J dc_jb
so the energy error is second order in the fit error and no symmetric factor J^-1/2 is needed. The non-robust
expressions (1/2 v^T G v, A^T G_sym A) are reported for comparison.

Units: atomic units. Learning is done in Coulomb-normalised variables (Jn = S^-1 J S^-1, S = diag(J)^1/2),
Gn ~ Jn^-1, G = S^-1 Gn S^-1. This fixed diagonal scaling only conditions the optimisation.

Reuses scripts/local_mp2/lmp2_common.py (geometries, JKFIT parser, cached RI-JK RHF/def2-SVP).
Claude Generated (Oct 2026). AI-generated, machine-tested only.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "local_mp2"))
import lmp2_common as lc  # sets thread environment before numpy loads

import numpy as np
import scipy.linalg

# extra alkanes for training (alkane4/8/12/16/20 exist in lmp2_common)
for _n in (6, 10):
    lc.SYSTEMS[f"alkane{_n}"] = ("alkanes", lambda n=_n: lc.alkane(n))
    lc.LABEL[f"alkane{_n}"] = f"C{_n}H{2 * _n + 2}"


# ----------------------------------------------------------------------------------------------------------
# Per-system data: metric, scaling, right-hand sides
# ----------------------------------------------------------------------------------------------------------
class Data:
    pass


def load(name, mp2=True, verbose=True, keep_ao=False):
    """System setup + RI quantities. Fields: J, Jinv (dense exact), sc (= diag(J)^1/2), Jn, rho, rho_A, Araw
    (canonical valence i, all virtuals a; (o, v, naux)), E_J, E_mp2 (exact-fit references); keep_ao: Vao =
    (P|mu nu) for mu <= nu as columns (naux x npair)."""
    t0 = time.time()
    s = lc.setup(name, loc="can" if mp2 else None, verbose=False)
    d = Data()
    d.name, d.s, d.natm, d.naux = name, s, s.natm, s.naux
    d.J = lc.ri_metric(s)
    d.cf = scipy.linalg.cho_factor(d.J, lower=True)
    d.Jinv = scipy.linalg.cho_solve(d.cf, np.eye(d.naux))
    d.Jinv = 0.5 * (d.Jinv + d.Jinv.T)
    d.sc = np.sqrt(np.diag(d.J))
    d.Jn = d.J / d.sc[:, None] / d.sc[None, :]
    d.aux_atom = s.aux_atom
    d.xyz = s.xyz
    d.Z = np.array([s.mol.atom_charge(a) for a in range(s.natm)])
    I3 = lc.int3c(s)
    Docc = 2 * s.Cocc @ s.Cocc.T
    d.rho_A = np.zeros((s.natm, d.naux))
    for a in range(s.natm):
        mu = np.where(s.ao_atom == a)[0]
        d.rho_A[a] = np.einsum("mn,mnP->P", Docc[mu], I3[mu], optimize=True)
    d.rho = d.rho_A.sum(axis=0)
    d.c_rho = scipy.linalg.cho_solve(d.cf, d.rho)
    d.E_J = 0.5 * d.rho @ d.c_rho
    if mp2:
        d.Araw = np.einsum("mi,na,mnP->iaP", s.Cval, s.Cvir, I3, optimize=True)
        d.eo, d.ev = s.eval_, s.evir
        Bglob = scipy.linalg.solve_triangular(s.Lj, d.Araw.reshape(-1, d.naux).T, lower=True).T.reshape(d.Araw.shape)
        d.E_mp2 = lc.rimp2_canonical(s, Bglob)
        del Bglob
    if keep_ao:
        iu = np.triu_indices(s.n)
        d.Vao = np.ascontiguousarray(I3[iu].T)
    s._int3c = None
    del I3
    if verbose:
        print(f"# {name}: natm={d.natm} naux={d.naux} cond(J)={np.linalg.cond(d.J):.1e} E_J={d.E_J:.8f}"
              + (f" E_mp2={d.E_mp2:.10f}" if mp2 else "") + f" ({time.time() - t0:.1f}s)", flush=True)
    return d


def rhs_pool(d, rng, n_mp2=None):
    """Training / validation right-hand sides in the dual space (columns): (P|rho), (P|rho_A), (P|ia).
    Returns V (naux x m). n_mp2 limits the number of (P|ia) columns (random subset)."""
    A = d.Araw.reshape(-1, d.naux)
    idx = np.arange(A.shape[0])
    if n_mp2 is not None and n_mp2 < len(idx):
        idx = rng.choice(idx, n_mp2, replace=False)
    return np.concatenate([d.rho[:, None], d.rho_A.T, A[idx].T], axis=1)


# ----------------------------------------------------------------------------------------------------------
# Operators: every approximate inverse is an object with apply(V) (V: naux x m), dense() and counts
# ----------------------------------------------------------------------------------------------------------
class DenseOp:
    """G given as a dense matrix (used for global low-rank forms and for evaluation)."""

    def __init__(self, G, params, flops, label):
        self.G, self.params, self.flops, self.label = G, params, flops, label

    def apply(self, V):
        return self.G @ V

    def dense(self):
        return self.G


def exact_op(d):
    return DenseOp(d.Jinv, d.naux * (d.naux + 1) // 2, 2 * d.naux ** 2, "exact J^-1 (Cholesky)")


def global_lowrank_ops(d, k, Vtrain=None):
    """Global rank-k forms (non-learned): (a) truncated eigen-decomposition of J^-1 (largest eigenvalues of J^-1 =
    smallest of J; Frobenius-optimal for J^-1), (b) largest eigenvalues of J, (c) data-optimal rank k for the
    relative J-norm loss on Vtrain (reduced-rank regression: H = J^-1/2 P_k J^-1/2, P_k = top-k left singular
    vectors of W = J^-1/2 V diag(w)^1/2). Parameter count = naux*k (symmetric factor)."""
    w, U = np.linalg.eigh(d.J)
    ops = []
    G = (U[:, :k] / w[:k]) @ U[:, :k].T
    ops.append(DenseOp(G, d.naux * k + k, 2 * 2 * d.naux * k, f"global trunc-eig J^-1, rank {k}"))
    G = (U[:, -k:] / w[-k:]) @ U[:, -k:].T
    ops.append(DenseOp(G, d.naux * k + k, 2 * 2 * d.naux * k, f"global large-eig J, rank {k}"))
    if Vtrain is not None:
        Jmh = (U / np.sqrt(w)) @ U.T  # J^-1/2
        C = d.Jinv @ Vtrain
        wt = 1.0 / np.einsum("pm,pm->m", Vtrain, C)
        W = (Jmh @ Vtrain) * np.sqrt(wt)
        P, _, _ = np.linalg.svd(W, full_matrices=False)
        Pk = P[:, :k]
        X = Jmh @ Pk
        ops.append(DenseOp(X @ X.T, d.naux * k, 2 * 2 * d.naux * k, f"global data-optimal, rank {k}"))
    return ops


# ----------------------------------------------------------------------------------------------------------
# Atom-centred overlapping subdomains, partition of unity, coarse space
# ----------------------------------------------------------------------------------------------------------
class Decomp:
    pass


def decomposition(d, R_sd, coarse="charge"):
    """Subdomain i = all aux functions on atoms within R_sd (Bohr) of atom i (the 'core' = functions on atom i).
    Functions inside a subdomain are ordered by (relative position of their atom to atom i, element), then file
    order, so translation-equivalent subdomains get identical local matrices (checked in shared_classes).
    m_P = number of subdomains containing P, D_i = diag(1/m_P).
    Coarse space Phi (naux x n0): 'charge' = one column per atom, the charge integrals q_P = int chi_P of the
    functions on that atom (normalised); 'charge+dipole' adds the three dipole-moment columns per atom."""
    dc = Decomp()
    dc.R = R_sd
    natm = d.natm
    Dat = np.linalg.norm(d.xyz[:, None, :] - d.xyz[None, :, :], axis=2)
    fn_of_atom = [np.where(d.aux_atom == a)[0] for a in range(natm)]
    dc.idx, dc.keys, dc.core = [], [], []
    for i in range(natm):
        members = [j for j in range(natm) if Dat[i, j] <= R_sd + 1e-9]
        rel = {j: tuple(np.round(d.xyz[j] - d.xyz[i], 3) + 0.0) for j in members}
        members.sort(key=lambda j: (rel[j], d.Z[j]))
        idx = np.concatenate([fn_of_atom[j] for j in members])
        dc.idx.append(idx)
        dc.keys.append(tuple((int(d.Z[j]),) + rel[j] for j in members))
        dc.core.append(np.isin(idx, fn_of_atom[i]))
    m = np.zeros(d.naux)
    for idx in dc.idx:
        m[idx] += 1
    dc.m = m
    dc.sq = [1.0 / np.sqrt(m[idx]) for idx in dc.idx]  # diag of D_i^1/2
    dc.Phi = coarse_space(d, coarse)
    dc.n0 = dc.Phi.shape[1]
    return dc


def coarse_space(d, kind):
    """Coarse basis in coefficient space from the multipole integrals of the aux functions (Becke grid)."""
    from pyscf import dft
    aux = d.s.auxmol
    g = dft.gen_grid.Grids(aux)
    g.level = 3
    g.build()
    ao = dft.numint.eval_ao(aux, g.coords)  # (ngrid, naux)
    q = g.weights @ ao
    cols = []
    for a in range(d.natm):
        on = d.aux_atom == a
        v = np.where(on, q, 0.0)
        cols.append(v / np.linalg.norm(v))
        if kind == "charge+dipole":
            for c in range(3):
                mu = (g.weights * (g.coords[:, c] - d.xyz[a, c])) @ ao
                v = np.where(on, mu, 0.0)
                if np.linalg.norm(v) > 1e-12:
                    cols.append(v / np.linalg.norm(v))
    return np.array(cols).T


def galerkin_coarse(d, dc):
    """P0 = Phi (Phi^T J Phi)^-1 Phi^T (exact inverse on the coarse space)."""
    A0 = dc.Phi.T @ d.J @ dc.Phi
    return dc.Phi @ np.linalg.solve(A0, dc.Phi.T)


class SchwarzExactOp:
    """Classical (non-learned) Schwarz operators with exact local inverses J_ii^-1 = (R_i J R_i^T)^-1.
    variant: 'asm'      sum R_i^T J_ii^-1 R_i                        (plain additive Schwarz, symmetric)
             'asm-pu'   sum R_i^T D^1/2 J_ii^-1 D^1/2 R_i            (paper's weighting, symmetric)
             'ras-core' rows of J_ii^-1 belonging to the core atom i  (restricted, non-symmetric local inverse)
    coarse:  None | 'add' (+ Phi (Phi^T J Phi)^-1 Phi^T) | 'bal' (P0 + (I - P0 J) M (I - J P0), balancing/deflation).
    """

    def __init__(self, d, dc, variant="asm-pu", coarse="bal"):
        n = d.naux
        M = np.zeros((n, n))
        params = 0
        flops = 0
        for i, idx in enumerate(dc.idx):
            Li = np.linalg.inv(d.J[np.ix_(idx, idx)])
            Li = 0.5 * (Li + Li.T)
            ni = len(idx)
            params += ni * (ni + 1) // 2
            if variant == "asm":
                M[np.ix_(idx, idx)] += Li
                flops += 2 * ni * ni
            elif variant == "asm-pu":
                M[np.ix_(idx, idx)] += dc.sq[i][:, None] * Li * dc.sq[i][None, :]
                flops += 2 * ni * ni
            elif variant == "ras-core":
                core = dc.core[i]
                M[np.ix_(idx[core], idx)] += Li[core]
                flops += 2 * core.sum() * ni
            else:
                raise ValueError(variant)
        self.M1 = M
        if coarse is None:
            G = M
        else:
            P0 = galerkin_coarse(d, dc)
            params += dc.n0 * (dc.n0 + 1) // 2 + n  # coarse factor + Phi (one nonzero block per atom column)
            flops += 4 * n + 2 * dc.n0 ** 2
            if coarse == "add":
                G = M + P0
            elif coarse == "bal":
                IP = np.eye(n) - P0 @ d.J
                G = P0 + IP @ M @ IP.T
                # two extra coarse solves and two products with the stored J Phi (n x n0)
                flops += 2 * (4 * n + 2 * dc.n0 ** 2) + 2 * 2 * n * dc.n0
                params += n * dc.n0
            else:
                raise ValueError(coarse)
        self.G = G
        self.params, self.flops = params, flops
        self.label = f"Schwarz exact {variant}" + (f" + coarse {coarse}" if coarse else "") + f" (R={dc.R})"

    def apply(self, V):
        return self.G @ V

    def dense(self):
        return self.G


# ----------------------------------------------------------------------------------------------------------
# Evaluation
# ----------------------------------------------------------------------------------------------------------
def rel_jnorm_err(d, G, V):
    """mean over columns of dc^T J dc / (v^T J^-1 v) (= 2 x relative robust RI-J error per rhs)."""
    C = d.Jinv @ V
    E = G @ V - C
    return float(np.mean(np.einsum("pm,pm->m", E, d.J @ E) / np.einsum("pm,pm->m", V, C)))


def evaluate(d, op, Vval=None, mp2=True, pcg=True):
    """Operator and downstream errors of an approximate inverse (op.dense() is used; sizes here are small)."""
    G = op.dense()
    r = {"label": op.label, "params": int(op.params), "flops_apply": int(op.flops),
         "params_dense": d.naux * (d.naux + 1) // 2, "flops_dense": 2 * d.naux ** 2}
    r["frob_rel"] = float(np.linalg.norm(G - d.Jinv) / np.linalg.norm(d.Jinv))
    # Coulomb-normalised Frobenius error (insensitive to the near-null space of J)
    Gn = G * d.sc[:, None] * d.sc[None, :]
    Jni = d.Jinv * d.sc[:, None] * d.sc[None, :]
    r["frob_rel_norm"] = float(np.linalg.norm(Gn - Jni) / np.linalg.norm(Jni))
    # J-weighted identity error ||J^1/2 G J^1/2 - I||_F / sqrt(naux)
    w, U = np.linalg.eigh(d.J)
    Jh = (U * np.sqrt(w)) @ U.T
    H = Jh @ G @ Jh
    r["jw_identity_err"] = float(np.linalg.norm(H - np.eye(d.naux)) / np.sqrt(d.naux))
    Hs = 0.5 * (H + H.T)
    ev = np.linalg.eigvalsh(Hs)
    r["eig_JG_min"], r["eig_JG_max"] = float(ev[0]), float(ev[-1])
    # RI-J
    c = G @ d.rho
    e_rob = d.rho @ c - 0.5 * c @ d.J @ c
    r["EJ_robust_err"] = float(e_rob - d.E_J)
    r["EJ_robust_err_atom"] = float((e_rob - d.E_J) / d.natm)
    Gs = 0.5 * (G + G.T)
    r["EJ_plain_err"] = float(0.5 * d.rho @ Gs @ d.rho - d.E_J)
    if Vval is not None:
        r["val_rel_jnorm"] = rel_jnorm_err(d, G, Vval)
    if mp2:
        r["Emp2_robust_err"], r["Emp2_plain_err"] = mp2_errors(d, G)
        r["Emp2_robust_err_atom"] = r["Emp2_robust_err"] / d.natm
    if pcg:
        r["pcg_iter"] = pcg_iterations(d, Gs, d.rho)
        r["gmres_iter"] = gmres_iterations(d, G, d.rho)
    return r


def gmres_iterations(d, M, b, tol=1e-10, maxit=400):
    """Right-preconditioned full GMRES (no restart) for J x = b with x = M y; iterations until the residual
    estimate ||r||/||b|| < tol. Works for non-symmetric M. Returns maxit if not converged."""
    n = len(b)
    beta = np.linalg.norm(b)
    Vk = np.zeros((maxit + 1, n))
    H = np.zeros((maxit + 1, maxit))
    Vk[0] = b / beta
    cs, sn = np.zeros(maxit), np.zeros(maxit)
    g = np.zeros(maxit + 1)
    g[0] = beta
    for k in range(maxit):
        w = d.J @ (M @ Vk[k])
        for _ in range(2):  # classical Gram-Schmidt, twice
            h = Vk[:k + 1] @ w
            w -= h @ Vk[:k + 1]
            H[:k + 1, k] += h
        H[k + 1, k] = np.linalg.norm(w)
        if H[k + 1, k] > 1e-300:
            Vk[k + 1] = w / H[k + 1, k]
        for j in range(k):  # apply previous Givens rotations
            t = cs[j] * H[j, k] + sn[j] * H[j + 1, k]
            H[j + 1, k] = -sn[j] * H[j, k] + cs[j] * H[j + 1, k]
            H[j, k] = t
        rr = np.hypot(H[k, k], H[k + 1, k])
        cs[k], sn[k] = H[k, k] / rr, H[k + 1, k] / rr
        H[k, k], H[k + 1, k] = rr, 0.0
        g[k + 1] = -sn[k] * g[k]
        g[k] = cs[k] * g[k]
        if abs(g[k + 1]) / beta < tol:
            return k + 1
    return maxit


def reference_preconditioners(d):
    """Identity and Jacobi (diag J)^-1 as preconditioner references."""
    return [DenseOp(np.eye(d.naux), 0, 0, "none (identity)"),
            DenseOp(np.diag(1.0 / np.diag(d.J)), d.naux, d.naux, "Jacobi diag(J)^-1")]


def local_density_fit_error(d, R):
    """Reference outside the operator class (q4_metric.py (b)): fit each atomic density rho_A with the aux
    functions within R of atom A (exact local inverse), robust RI-J energy. Not a linear map of v = (P|rho)."""
    Dat = np.linalg.norm(d.xyz[:, None, :] - d.xyz[None, :, :], axis=2)
    c = np.zeros(d.naux)
    for a in range(d.natm):
        idx = np.where(Dat[a][d.aux_atom] <= R)[0]
        c[idx] += np.linalg.solve(d.J[np.ix_(idx, idx)], d.rho_A[a][idx])
    return float(d.rho @ c - 0.5 * c @ d.J @ c - d.E_J)


def mp2_errors(d, G):
    """RI-MP2 (frozen core, canonical) with c_ia = G (P|ia): robust and plain (A^T G_sym A) energy errors."""
    o, v, n = d.Araw.shape
    A = d.Araw
    C = (A.reshape(-1, n) @ G.T).reshape(o, v, n)
    CJ = (C.reshape(-1, n) @ d.J).reshape(o, v, n)
    Gs = 0.5 * (G + G.T)
    AG = (A.reshape(-1, n) @ Gs).reshape(o, v, n)
    e_r = e_p = 0.0
    for i in range(o):
        K = (np.einsum("aP,jbP->jab", A[i], C, optimize=True) + np.einsum("aP,jbP->jab", C[i], A, optimize=True)
             - np.einsum("aP,jbP->jab", CJ[i], C, optimize=True))
        Kp = np.einsum("aP,jbP->jab", AG[i], A, optimize=True)
        D = d.eo[i] + d.eo[:, None, None] - d.ev[None, :, None] - d.ev[None, None, :]
        e_r += np.sum(K * (2 * K - K.transpose(0, 2, 1)) / D)
        e_p += np.sum(Kp * (2 * Kp - Kp.transpose(0, 2, 1)) / D)
    return float(e_r - d.E_mp2), float(e_p - d.E_mp2)


def pcg_iterations(d, M, b, tol=1e-10, maxit=2000):
    """Iterations of preconditioned CG for J x = b, ||r||/||b|| < tol (M must be SPD; if not, returns -1 when
    the iteration breaks down)."""
    x = np.zeros_like(b)
    r = b.copy()
    z = M @ r
    p = z.copy()
    rz = r @ z
    nb = np.linalg.norm(b)
    for it in range(1, maxit + 1):
        Ap = d.J @ p
        pAp = p @ Ap
        if pAp <= 0 or rz <= 0:
            return -1
        a = rz / pAp
        x += a * p
        r -= a * Ap
        if np.linalg.norm(r) / nb < tol:
            return it
        z = M @ r
        rz_new = r @ z
        p = z + (rz_new / rz) * p
        rz = rz_new
    return maxit
