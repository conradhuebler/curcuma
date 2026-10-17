#!/usr/bin/env python3
"""Q5 of the local-MP2 feasibility study: preconditioners for the coupled PNO-LMP2 amplitude equations.

The LMP2 residual (lmp2_common docstring) is linear in T:  R(T) = K + A T.  In the pseudo-canonical PNO basis of
pair p=(i,j) the pair's own block of A is diagonal (eps_a + eps_b - F_ii - F_jj); pairs couple only through F_ik
(k != i) and F_kj (k != j) times PNO overlap matrices S^{pq} = V_p^T S V_q. Here A is applied in the PNO spaces with
the S^{pq} stored for every coupling with |F| > 1e-10 (energy checked against the AO-embedded solver).

Preconditioners (all used inside the same outer iteration T <- T + M^-1(-R), optionally with the same DIIS):
  jacobi       : M = diagonal (pseudo-canonical denominators) -- the standard
  block-pair   : M = block diagonal per pair = identical to jacobi in this basis (run as a check)
  ras-pair:tau : restricted additive Schwarz, one subdomain per pair p = {p} + every pair q coupled to p with
                 |F| > tau; the subdomain system A_s x = -R_s is solved iteratively (inner Jacobi+DIIS to a relative
                 residual 1e-4, at most 50 steps), only the update of p is kept (RAS, cf. localinv of the EEQ solver)
  ras-group:g:tau : subdomains from groups of <= g localised orbitals (recursive coordinate bisection of the
                 centroids); owned pairs = pairs whose first orbital is in the group, overlap = pairs coupled to an
                 owned pair with |F| > tau; tau = inf gives non-overlapping block Jacobi over groups
Convergence: |dE| < 1e-8 Eh and max|R| < 1e-6 Eh. Cost: S T S^T block products counted as flops
(2 m_p m_q (m_p + m_q) each) for the outer residual and for all inner subdomain solves, plus wall time (n=1).

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/q5_schwarz.py [--tcut 1e-7] [systems...]
Claude Generated (Oct 2026).
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp2_common as lc  # first: sets the thread environment before numpy loads
import numpy as np


DEFAULT = ["alkane8", "alkane12", "alkane16", "water8", "water12", "water16", "acene2", "acene3", "acene4"]
VARIANTS = ["jacobi", "jacobi+diis", "block-pair+diis", "ras-pair:1e-2+diis", "ras-pair:1e-3+diis",
            "ras-pair:1e-2", "ras-group:6:inf+diis", "ras-group:6:1e-2+diis", "ras-group:6:1e-3+diis"]
E_TOL, R_TOL, MAXIT = 1e-8, 1e-6, 200


class PNOSystem:
    """Pairs, PNO spaces, couplings with stored overlaps; linear operator A restricted to any pair subset."""

    def __init__(self, s, tcut):
        self.s = s
        Bao = lc.b_local_ao(s)
        self.e_ref = lc.rimp2_canonical(s, lc.b_canonical(s, Bao))
        Bv = np.einsum("imQ,ma->iaQ", Bao, s.Cvir, optimize=True)
        self.pairs = lc.pair_list(s)
        self.pidx = {pq: k for k, pq in enumerate(self.pairs)}
        Fd = np.diag(s.Focc)
        self.V, self.K, self.den, self.eps, self.dE = [], [], [], [], 0.0
        for i, j in self.pairs:
            V, e, K, de, _ = lc.make_pnos((s.Cvir, s.evir), Bv[i] @ Bv[j].T, Fd[i], Fd[j], i, j, tcut)
            self.V.append(V)
            self.eps.append(e)
            self.K.append(K)
            self.den.append(e[:, None] + e[None, :] - Fd[i] - Fd[j])
            self.dE += de
        self.m = np.array([v.shape[1] for v in self.V])
        SV = [s.S @ v for v in self.V]
        # couplings: (p, q, coef, transpose?) ; T'_{kj} = T_q or T_q^T
        ent = []
        for p, (i, j) in enumerate(self.pairs):
            for k in range(s.nval):
                if k != i and abs(s.Focc[i, k]) > 1e-10:
                    ent.append((p, self.pidx[(min(k, j), max(k, j))], -s.Focc[i, k], k > j))
                if k != j and abs(s.Focc[k, j]) > 1e-10:
                    ent.append((p, self.pidx[(min(i, k), max(i, k))], -s.Focc[k, j], i > k))
        self.ent = ent
        self.S = [self.V[p].T @ SV[q] for p, q, _, _ in ent]
        self.Fabs = np.array([abs(c) for _, _, c, _ in ent])
        self.flop = np.array([2.0 * self.m[p] * self.m[q] * (self.m[p] + self.m[q]) for p, q, _, _ in ent])
        self.by_p = [[] for _ in self.pairs]
        for e_, (p, q, c, tr) in enumerate(ent):
            self.by_p[p].append(e_)

    def apply_A(self, T, members=None, ent_sel=None):
        """(A T)_p for p in members (default all); couplings restricted to ent_sel (indices). Returns dict/list."""
        if members is None:
            out = [self.den[p] * T[p] for p in range(len(self.pairs))]
            sel = range(len(self.ent))
        else:
            out = {p: self.den[p] * T[p] for p in members}
            sel = ent_sel
        fl = 0.0
        for e_ in sel:
            p, q, c, tr = self.ent[e_]
            if self.m[p] == 0 or self.m[q] == 0:
                continue
            Tq = T[q].T if tr else T[q]
            out[p] = out[p] + c * (self.S[e_] @ Tq @ self.S[e_].T)
            fl += self.flop[e_]
        return out, fl

    def energy(self, T):
        return sum(lc.pair_energy(T[p], self.K[p], i, j) for p, (i, j) in enumerate(self.pairs))


def groups_rcb(cen, g):
    """Recursive coordinate bisection of orbital centroids into groups of <= g orbitals (deterministic)."""
    def rec(idx):
        if len(idx) <= g:
            return [sorted(idx)]
        X = cen[idx]
        ax = int(np.argmax(np.ptp(X, axis=0)))
        o = [idx[k] for k in np.argsort(X[:, ax], kind="stable")]
        h = len(o) // 2
        return rec(o[:h]) + rec(o[h:])
    return rec(list(range(len(cen))))


def make_subdomains(P, kind, tau, g=6):
    """List of (members, owned) pair-index lists."""
    npair = len(P.pairs)
    nb = [set() for _ in range(npair)]
    for e_, (p, q, c, tr) in enumerate(P.ent):
        if abs(c) > tau:
            nb[p].add(q)
            nb[q].add(p)
    subs = []
    if kind == "pair":
        for p in range(npair):
            subs.append((sorted({p} | nb[p]), [p]))
    else:
        grp = groups_rcb(P.s.centroid, g)
        gid = {i: k for k, G in enumerate(grp) for i in G}
        for k in range(len(grp)):
            owned = [p for p, (i, j) in enumerate(P.pairs) if gid[i] == k]
            mem = set(owned)
            for p in owned:
                mem |= nb[p]
            subs.append((sorted(mem), owned))
    for k, (mem, own) in enumerate(subs):
        ms = set(mem)
        ent_sel = [e_ for p in mem for e_ in P.by_p[p] if P.ent[e_][1] in ms]
        subs[k] = (mem, own, ent_sel)
    return subs


def subdomain_solve(P, sub, rhs, tol_rel=1e-4, maxit=50):
    """Solve A_s x = rhs on the subdomain by Jacobi+DIIS; returns x (dict), inner iterations, flops."""
    mem, own, ent_sel = sub
    x = {p: rhs[p] / P.den[p] for p in mem}
    r0 = max((np.max(np.abs(rhs[p])) if rhs[p].size else 0.0) for p in mem)
    if r0 == 0.0:
        return x, 0, 0.0
    shapes = [x[p].shape for p in mem]
    offs = np.concatenate([[0], np.cumsum([np.prod(sh) for sh in shapes])]).astype(int)
    dd = lc.DIIS(6)
    fl = 0.0
    it = 0
    for it in range(1, maxit + 1):
        Ax, f = P.apply_A(x, members=mem, ent_sel=ent_sel)
        fl += f
        r = {p: rhs[p] - Ax[p] for p in mem}
        rm = max((np.max(np.abs(r[p])) if r[p].size else 0.0) for p in mem)
        if rm < tol_rel * r0:
            break
        xn = np.concatenate([(x[p] + r[p] / P.den[p]).ravel() for p in mem])
        er = np.concatenate([(r[p] / P.den[p]).ravel() for p in mem])
        xn = dd.update(xn, er)
        x = {p: xn[offs[k]:offs[k + 1]].reshape(shapes[k]) for k, p in enumerate(mem)}
    return x, it, fl


def run(P, variant):
    use_diis = variant.endswith("+diis")
    v = variant.replace("+diis", "")
    subs = None
    if v.startswith("ras-pair"):
        subs = make_subdomains(P, "pair", float(v.split(":")[1]))
    elif v.startswith("ras-group"):
        _, g, tau = v.split(":")
        subs = make_subdomains(P, "group", float(tau), int(g))
    T = [-P.K[p] / P.den[p] for p in range(len(P.pairs))]
    sizes = [t.size for t in T]
    offs = np.concatenate([[0], np.cumsum(sizes)]).astype(int)
    unflat = lambda x: [x[offs[p]:offs[p + 1]].reshape(T[p].shape) for p in range(len(P.pairs))]
    dd = lc.DIIS(8) if use_diis else None
    e_old = P.energy(T)
    fl_out = fl_in = 0.0
    inner = []
    t0 = time.time()
    hist = []
    it = 0
    for it in range(1, MAXIT + 1):
        AT, f = P.apply_A(T)
        fl_out += f
        R = [P.K[p] + AT[p] for p in range(len(P.pairs))]
        rmax = max((np.max(np.abs(r)) if r.size else 0.0) for r in R)
        if subs is None:
            dT = [-R[p] / P.den[p] for p in range(len(P.pairs))]
        else:
            dT = [None] * len(P.pairs)
            negR = [-r for r in R]
            for sub in subs:
                x, nin, f = subdomain_solve(P, sub, negR)
                fl_in += f
                inner.append(nin)
                for p in sub[1]:
                    dT[p] = x[p]
        Tn = [T[p] + dT[p] for p in range(len(P.pairs))]
        if dd is not None:
            Tn = unflat(dd.update(np.concatenate([t.ravel() for t in Tn]), np.concatenate([d.ravel() for d in dT])))
        T = Tn
        e = P.energy(T)
        hist.append([e, rmax])
        if abs(e - e_old) < E_TOL and rmax < R_TOL:
            break
        e_old = e
    res = {"variant": variant, "iterations": it, "energy": e, "flop_outer": fl_out, "flop_inner": fl_in,
           "inner_mean": float(np.mean(inner)) if inner else 0.0, "t_s": time.time() - t0,
           "history": hist, "converged": it < MAXIT}
    if subs is not None:
        res["sub_mean"] = float(np.mean([len(sb[0]) for sb in subs]))
        res["sub_max"] = int(max(len(sb[0]) for sb in subs))
        res["nsub"] = len(subs)
    return res


def main():
    args = sys.argv[1:]
    tcut = 1e-7
    if "--tcut" in args:
        k = args.index("--tcut")
        tcut = float(args[k + 1])
        args = args[:k] + args[k + 2:]
    variants = VARIANTS
    if "--variants" in args:
        k = args.index("--variants")
        variants = args[k + 1].split(",")
        args = args[:k] + args[k + 2:]
    names = args or DEFAULT
    out = []
    for name in names:
        s = lc.setup(name, loc="pm")
        t0 = time.time()
        P = PNOSystem(s, tcut)
        # cross-check of the PNO-space operator against the AO-embedded solver of lmp2_common
        chk = lc.solve_lmp2(s, P.pairs, list(zip(P.V, P.eps)), P.K, coupled=True, e_tol=1e-10, r_tol=1e-8)
        offd = np.abs(s.Focc - np.diag(np.diag(s.Focc)))
        print(f"# {name}: pairs={len(P.pairs)} couplings={len(P.ent)} PNO mean={P.m.mean():.1f} max={P.m.max()}"
              f" max|F_ik|={offd.max():.3f} n(|F|>1e-2)={int(np.sum(offd > 1e-2)) // 2}"
              f" n(|F|>1e-3)={int(np.sum(offd > 1e-3)) // 2} setup {time.time() - t0:.0f}s", flush=True)
        for var in variants:
            r = run(P, var)
            r.update({"system": name, "class": lc.SYSTEMS[name][0], "tcut": tcut, "npairs": len(P.pairs),
                      "ncoupl": len(P.ent), "pno_mean": float(P.m.mean()), "e_rimp2": P.e_ref,
                      "err_corr": r["energy"] + P.dE - P.e_ref, "flop_per_matvec": float(P.flop.sum()),
                      "max_offdiag_F": float(offd.max()), "e_ao_solver": chk["energy"],
                      "d_vs_ao_solver": r["energy"] - chk["energy"]})
            errs = [abs(h[0] - chk["energy"]) for h in r["history"]]
            r["it_to_tight_1e-8"] = next((k + 1 for k in range(len(errs)) if max(errs[k:]) < 1e-8), None)
            out.append(r)
            print(f"{name:9s} {var:24s} it={r['iterations']:3d} (tight {r['it_to_tight_1e-8']}) E={r['energy']:.10f} err(+dE)={r['err_corr']:.2e}"
                  f" matvec-equiv={(r['flop_outer'] + r['flop_inner']) / P.flop.sum():7.1f}"
                  f" dAO={r['energy'] - chk['energy']:.1e} inner={r['inner_mean']:.1f} sub={r.get('sub_mean', 0):.1f} {r['t_s']:.0f}s", flush=True)
            lc.save_json(f"q5_schwarz_{tcut:.0e}.json", out)


if __name__ == "__main__":
    main()
