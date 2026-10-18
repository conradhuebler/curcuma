#!/usr/bin/env python3
"""Stage B of the local-CCSD prototype: CCSD restricted to overlapping pair domains ("Schwarz domains").

docs/LOCAL_CCSD_PROTOTYPE_2026-10.md has the results. Python/PySCF only, nothing in the curcuma engine.

Setup: exact-integral RHF/def2-SVP, frozen core (ccsd_canonical.rhf), Pipek-Mezey valence orbitals
(lmp2_common.localise), all pairs i <= j kept (no pair screening, no weak-pair treatment).
Amplitude space (the truncation that is measured):
  - pair p = (i,j): t_ij^ab = [V_p T_p V_p^T]_ab with V_p (v x m_p) an orthonormal, pseudo-canonical basis of the pair
    virtual space, written in the canonical virtual MOs; t_ji = t_ij^T
  - singles of orbital i: t_i = V_ii t~_i (the pair space of the diagonal pair ii, as in DLPNO's singles PNOs)
  - spaces: "bp<thr>[+r<R>]" = PAOs of the Boughton-Pulay domain [i] u [j], each extended by all atoms within R Bohr
    (lmp2_common.bp_domain/extend_domain/pair_space, 1e-6 redundancy cut); "pno<tcut>" = PNOs from semicanonical
    MP2 amplitudes in all virtuals (lmp2_common.make_pnos); "<pao>|pno<tcut>" = PNOs inside the PAO pair space
Equations (projected / Galerkin local CCSD): the full CCSD residual R(t) of ccsd_canonical.residual is evaluated
densely for the embedded amplitudes and projected, R_p = V_p^T R_ij V_p and R_i = V_ii^T R_i. All couplings between
pairs (Fock, ladder, ring, T1 terms) are therefore exact for the truncated amplitudes; the only approximation is the
amplitude space itself. The dense evaluation means the wall time is that of canonical CCSD; the cost of a real
local implementation is given as counted flops (proxies below), not measured.
Preconditioners (both with Pulay DIIS on the pair-space vectors, error vector = the update):
  jacobi : T_p <- T_p - R_p / (eps_a + eps_b - F_ii - F_jj), t~_i <- t~_i - R_i / (eps_a - F_ii)
  ras-pair:<tau>, ras-group:<g>:<tau> : restricted additive Schwarz on overlapping pair subdomains (q5_schwarz of
    the local-MP2 study, reused unchanged): each subdomain solves A0_s x = -R_s, A0 = the Fock part of the doubles
    Jacobian (pseudo-canonical diagonal + occupied coupling -F_ik S^pq T_q S^pq^T), inner Jacobi+DIIS to a relative
    residual 1e-4; singles stay Jacobi. A0 is the LMP2 operator: the ladder, ring and T1 parts of the CCSD Jacobian
    are not in the preconditioner.
Flop proxies (counts, per outer iteration, pairs i <= j): ladder = sum_p 2 m_p^4 (particle-particle ladder in the pair
basis, a lower bound for a local CCSD residual); fock_all = sum over all couplings (p, q=(kj)|(ik), every k) of the
block products S T S^T, 2 m_p m_q (m_p + m_q); Schwarz inner work = the same block products inside all subdomain solves.
Canonical reference per iteration: o(o+1)/2 * 2 v^4.

Usage: LCCSD_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_ccsd/ccsd_local.py accuracy [systems...]
       ... ccsd_local.py schwarz [--space pno1e-7] [systems...]
       ... ccsd_local.py check [systems...]     (complete spaces reproduce canonical CCSD)
Claude Generated (Oct 2026). AI-generated, machine-tested only.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccsd_canonical as cc0  # first: sets the thread environment before numpy loads
import lmp2_common as lc
import q5_schwarz as q5

import json
import time

import numpy as np

ACC_SYSTEMS = ["alkane4", "alkane6", "alkane8", "alkane10", "water2", "water3", "water4", "water5", "water6",
               "acene1", "acene2", "acene3"]
ACC_SPACES = ["full", "bp0.02", "bp0.02+r3", "bp0.02+r5", "bp0.02+r7", "pno1e-6", "pno1e-7", "pno1e-8",
              "bp0.02+r5|pno1e-7", "bp0.02+r7|pno1e-8"]
SCH_SYSTEMS = ["alkane6", "alkane8", "alkane10", "water4", "water6", "acene2", "acene3"]
SCH_VARIANTS = ["jacobi+diis", "jacobi", "ras-pair:1e-2+diis", "ras-pair:1e-3+diis", "ras-group:6:inf+diis",
                "ras-group:6:1e-2+diis", "ras-group:6:1e-3+diis"]


# ----------------------------------------------------------------------------------------------------------
# System, canonical reference (cached)
# ----------------------------------------------------------------------------------------------------------
def setup(name):
    t0 = time.time()
    s = cc0.rhf(name)
    lc.localise(s, "pm")
    s.Ican = cc0.MOInts(s.mol, s.Cval, s.Cvir, s.F)
    s.I = s.Ican.rotate_occ(s.U)
    cache = os.path.join(lc.work_dir(), f"canonical_{name}.json")
    if os.path.exists(cache):
        ref = json.load(open(cache))
    else:
        r = cc0.solve_ccsd(s.Ican, e_tol=1e-11, r_tol=1e-9)
        rm = cc0.solve_ccsd(s.Ican, e_tol=1e-12, r_tol=1e-10, mp2_only=True)
        ref = {"e_ccsd": r["e_corr"], "e_mp2": r["e_iter0"], "e_mp2_iter": rm["e_corr"],
               "iterations": r["iterations"], "t_s": r["t_s"], "nbf": s.n, "o": s.nval, "v": s.nvir,
               "natm": s.natm, "e_hf": s.ehf}
        json.dump(ref, open(cache, "w"), indent=1)
    s.ref = ref
    print(f"# {name}: {cc0.label(name)} natm={s.natm} nbf={s.n} o={s.nval} v={s.nvir} E_HF={s.ehf:.10f}"
          f" E_CCSD(can)={ref['e_ccsd']:.10f} ({ref['iterations']} it) E_MP2={ref['e_mp2']:.10f}"
          f" setup {time.time() - t0:.0f}s", flush=True)
    return s


# ----------------------------------------------------------------------------------------------------------
# Pair spaces
# ----------------------------------------------------------------------------------------------------------
def parse_pao(spec):
    """'bp0.02+r5' -> (0.02, 5.0)."""
    part = spec.split("+r")
    return float(part[0][2:]), (float(part[1]) if len(part) > 1 else 0.0)


def build_spaces(s, spec):
    """Returns dict(pairs, V (list v x m), eps, dE_pno, info). V in canonical virtual coordinates."""
    pairs = lc.pair_list(s)
    o = s.nval
    Fd = np.diag(s.Focc)
    pao, tcut = None, None
    for part in spec.split("|"):
        if part.startswith("bp"):
            pao = part
        elif part.startswith("pno"):
            tcut = float(part[3:])
        elif part != "full":
            raise ValueError(spec)
    CvS = s.Cvir.T @ s.S
    if pao is not None:
        thr, rext = parse_pao(pao)
        dom = [lc.extend_domain(s, lc.bp_domain(s, i, thr), rext) for i in range(o)]
        base = []
        for i, j in pairs:
            X, e = lc.pair_space(s, sorted(set(dom[i]) | set(dom[j])))
            base.append((CvS @ X, e))
        dom_atoms = [len(d) for d in dom]
    else:
        base = [(np.eye(s.nvir), s.evir)] * len(pairs)
        dom_atoms = [s.natm] * o
    dE = 0.0
    if tcut is not None:
        V, eps = [], []
        ovov = s.I.ovov
        for p, (i, j) in enumerate(pairs):
            Vb, eb = base[p]
            Kb = Vb.T @ ovov[i, :, j, :] @ Vb
            Vp, ep, _, de, _ = lc.make_pnos((Vb, eb), Kb, Fd[i], Fd[j], i, j, tcut)
            V.append(Vp)
            eps.append(ep)
            dE += de
    else:
        V, eps = [b[0] for b in base], [b[1] for b in base]
    return {"pairs": pairs, "V": V, "eps": eps, "dE_pno": dE, "dom_atoms": dom_atoms}


class LocalAmplitudes:
    """Pair-space amplitude containers, embedding into dense (o,o,v,v) tensors and projection of residuals."""

    def __init__(self, s, sp):
        self.s, self.o, self.v = s, s.nval, s.nvir
        self.pairs, self.V, self.eps = sp["pairs"], sp["V"], sp["eps"]
        self.pidx = {pq: k for k, pq in enumerate(self.pairs)}
        Fd = np.diag(s.Focc)
        self.den = [e[:, None] + e[None, :] - Fd[i] - Fd[j] for e, (i, j) in zip(self.eps, self.pairs)]
        self.sing = [self.pidx[(i, i)] for i in range(self.o)]
        self.den1 = [self.eps[p] - Fd[i] for i, p in enumerate(self.sing)]
        self.m = np.array([V.shape[1] for V in self.V])
        self.shapes = [(m, m) for m in self.m] + [(self.m[p],) for p in self.sing]
        self.offs = np.concatenate([[0], np.cumsum([int(np.prod(sh)) for sh in self.shapes])]).astype(int)

    def embed(self, T, t):
        o, v = self.o, self.v
        t2 = np.zeros((o, o, v, v))
        for p, (i, j) in enumerate(self.pairs):
            V = self.V[p]
            x = V @ T[p] @ V.T
            t2[i, j] = x
            if i != j:
                t2[j, i] = x.T
        t1 = np.zeros((o, v))
        for i, p in enumerate(self.sing):
            t1[i] = self.V[p] @ t[i]
        return t1, t2

    def project(self, R1, R2):
        R = [V.T @ R2[i, j] @ V for V, (i, j) in zip(self.V, self.pairs)]
        r = [self.V[p].T @ R1[i] for i, p in enumerate(self.sing)]
        return R, r

    def flat(self, T, t):
        return np.concatenate([x.ravel() for x in T] + [x.ravel() for x in t])

    def unflat(self, x):
        n2 = len(self.pairs)
        parts = [x[self.offs[k]:self.offs[k + 1]].reshape(self.shapes[k]) for k in range(len(self.shapes))]
        return parts[:n2], parts[n2:]

    def n_t2(self):
        return int(sum((1 if i == j else 2) * m * m for m, (i, j) in zip(self.m, self.pairs)))


# ----------------------------------------------------------------------------------------------------------
# Fock-coupling operator A0 for the Schwarz preconditioner (duck-typed to q5_schwarz.PNOSystem)
# ----------------------------------------------------------------------------------------------------------
class FockOperator:
    apply_A = q5.PNOSystem.apply_A
    energy = None

    def __init__(self, s, A):
        self.s, self.pairs, self.pidx, self.den, self.m = s, A.pairs, A.pidx, A.den, A.m
        ent = []
        for p, (i, j) in enumerate(self.pairs):
            for k in range(s.nval):
                if k != i and abs(s.Focc[i, k]) > 1e-10:
                    ent.append((p, self.pidx[(min(k, j), max(k, j))], -s.Focc[i, k], k > j))
                if k != j and abs(s.Focc[k, j]) > 1e-10:
                    ent.append((p, self.pidx[(min(i, k), max(i, k))], -s.Focc[k, j], i > k))
        self.ent = ent
        self.S = [A.V[p].T @ A.V[q] for p, q, _, _ in ent]
        self.flop = np.array([2.0 * self.m[p] * self.m[q] * (self.m[p] + self.m[q]) for p, q, _, _ in ent])
        self.by_p = [[] for _ in self.pairs]
        for e_, (p, q, c, tr) in enumerate(ent):
            self.by_p[p].append(e_)


def flop_proxies(s, A):
    o, v = s.nval, s.nvir
    m = A.m.astype(float)
    ladder = float(np.sum(2 * m ** 4))
    fock_all = 0.0
    for p, (i, j) in enumerate(A.pairs):
        for k in range(o):
            for q in ((min(k, j), max(k, j)), (min(i, k), max(i, k))):
                mq = m[A.pidx[q]]
                fock_all += 2 * m[p] * mq * (m[p] + mq)
    return {"ladder": ladder, "fock_all": fock_all, "canonical_ladder": o * (o + 1) / 2 * 2.0 * v ** 4}


# ----------------------------------------------------------------------------------------------------------
# Local solver
# ----------------------------------------------------------------------------------------------------------
def solve_local(s, A, mp2_only=False, variant="jacobi+diis", e_tol=1e-9, r_tol=1e-7, maxit=150, verbose=False):
    I = s.I
    use_diis = variant.endswith("+diis")
    v = variant.replace("+diis", "")
    P = subs = None
    if v.startswith("ras-"):
        P = FockOperator(s, A)
        if v.startswith("ras-pair"):
            subs = q5.make_subdomains(P, "pair", float(v.split(":")[1]))
        else:
            _, g, tau = v.split(":")
            subs = q5.make_subdomains(P, "group", float(tau), int(g))
    o = s.nval
    T = [-(V.T @ I.ovov[i, :, j, :] @ V) / d for V, (i, j), d in zip(A.V, A.pairs, A.den)]
    t = [np.zeros(A.m[p]) for p in A.sing]
    dd = lc.DIIS(8) if use_diis else None
    t1, t2 = A.embed(T, t)
    e_old = cc0.energy(t1, t2, I)
    hist, inner, fl_in = [e_old], [], 0.0
    t0 = time.time()
    it = 0
    for it in range(1, maxit + 1):
        R1, R2 = cc0.residual(t1, t2, I, mp2_only=mp2_only)
        R, r = A.project(R1, R2)
        rmax = max(max((np.max(np.abs(x)) if x.size else 0.0) for x in R),
                   max((np.max(np.abs(x)) if x.size else 0.0) for x in r))
        if subs is None:
            dT = [-R[p] / A.den[p] for p in range(len(A.pairs))]
        else:
            dT = [None] * len(A.pairs)
            negR = [-x for x in R]
            for sub in subs:
                x, nin, f = q5.subdomain_solve(P, sub, negR)
                fl_in += f
                inner.append(nin)
                for p in sub[1]:
                    dT[p] = x[p]
        dt = [-r[i] / A.den1[i] for i in range(o)]
        Tn = [T[p] + dT[p] for p in range(len(A.pairs))]
        tn = [t[i] + dt[i] for i in range(o)]
        if dd is not None:
            Tn, tn = A.unflat(dd.update(A.flat(Tn, tn), A.flat(dT, dt)))
        T, t = Tn, tn
        t1, t2 = A.embed(T, t)
        e = cc0.energy(t1, t2, I)
        hist.append(e)
        if verbose:
            print(f"   it {it:3d} E={e:.12f} dE={e - e_old:.2e} max|R|={rmax:.2e} {time.time() - t0:.0f}s", flush=True)
        if abs(e - e_old) < e_tol and rmax < r_tol:
            break
        e_old = e
    res = {"variant": variant, "energy": float(e), "iterations": it, "converged": it < maxit, "history": hist,
           "t_s": time.time() - t0, "flop_inner": fl_in, "inner_mean": float(np.mean(inner)) if inner else 0.0}
    if subs is not None:
        res["sub_mean"] = float(np.mean([len(sb[0]) for sb in subs]))
        res["nsub"] = len(subs)
        res["flop_fock_matvec"] = float(P.flop.sum())
    return res


# ----------------------------------------------------------------------------------------------------------
# Drivers
# ----------------------------------------------------------------------------------------------------------
def save(fname, obj):
    lc.save_json(fname, obj)


def run_check(names):
    """Complete pair spaces (all virtuals) must reproduce canonical CCSD and MP2 in the localised basis."""
    out = []
    for name in names:
        s = setup(name)
        A = LocalAmplitudes(s, build_spaces(s, "full"))
        r = solve_local(s, A, e_tol=1e-11, r_tol=1e-9)
        rm = solve_local(s, A, mp2_only=True, e_tol=1e-12, r_tol=1e-10)
        rec = {"system": name, "e_local_full": r["energy"], "d_ccsd": r["energy"] - s.ref["e_ccsd"],
               "iterations": r["iterations"], "d_mp2": rm["energy"] - s.ref["e_mp2"], "iterations_mp2": rm["iterations"]}
        out.append(rec)
        print(f"{name}: full-space local CCSD - canonical = {rec['d_ccsd']:.1e} ({r['iterations']} it, can "
              f"{s.ref['iterations']}); LMP2 - MP2 = {rec['d_mp2']:.1e} ({rm['iterations']} it)", flush=True)
        save("check_full.json", out)


def run_accuracy(names, spaces):
    out = []
    for name in names:
        s = setup(name)
        natm = s.natm
        for spec in spaces:
            t0 = time.time()
            sp = build_spaces(s, spec)
            A = LocalAmplitudes(s, sp)
            r = solve_local(s, A)
            rm = solve_local(s, A, mp2_only=True, e_tol=1e-10, r_tol=1e-8)
            fl = flop_proxies(s, A)
            ref = s.ref
            rec = {"system": name, "label": cc0.label(name), "natm": natm, "nbf": s.n, "o": s.nval, "v": s.nvir,
                   "space": spec, "e_ccsd_can": ref["e_ccsd"], "e_mp2_can": ref["e_mp2"],
                   "e_ccsd_loc": r["energy"], "e_mp2_loc": rm["energy"], "dE_pno": sp["dE_pno"],
                   "err_ccsd": r["energy"] - ref["e_ccsd"], "err_mp2": rm["energy"] - ref["e_mp2"],
                   "err_ccsd_mp2corr": r["energy"] + (ref["e_mp2"] - rm["energy"]) - ref["e_ccsd"],
                   "err_ccsd_pnocorr": r["energy"] + sp["dE_pno"] - ref["e_ccsd"],
                   "iterations": r["iterations"], "iterations_mp2": rm["iterations"], "converged": r["converged"],
                   "m_mean": float(A.m.mean()), "m_max": int(A.m.max()), "dom_atoms_mean": float(np.mean(sp["dom_atoms"])),
                   "n_t2": A.n_t2(), "n_t2_canonical": s.nval ** 2 * s.nvir ** 2, "flops": fl,
                   "t_s": time.time() - t0, "it_can": ref["iterations"]}
            rec["frac_t2"] = rec["n_t2"] / rec["n_t2_canonical"]
            out.append(rec)
            print(f"{name:9s} {spec:20s} m={A.m.mean():6.1f} it={r['iterations']:3d}/{ref['iterations']:3d}"
                  f" err CCSD={rec['err_ccsd']:.2e} ({rec['err_ccsd'] / natm:.1e}/atom) MP2={rec['err_mp2']:.2e}"
                  f" CCSD+MP2corr={rec['err_ccsd_mp2corr']:.2e} CCSD+dEpno={rec['err_ccsd_pnocorr']:.2e}"
                  f" t2 kept {100 * rec['frac_t2']:.1f}% ladder/can={fl['ladder'] / fl['canonical_ladder']:.3f}"
                  f" {rec['t_s']:.0f}s", flush=True)
            save("accuracy.json", out)


def run_schwarz(names, space, variants):
    out = []
    for name in names:
        s = setup(name)
        A = LocalAmplitudes(s, build_spaces(s, space))
        fl = flop_proxies(s, A)
        tight = solve_local(s, A, e_tol=1e-12, r_tol=1e-10)
        for var in variants:
            # run every variant to a tight criterion; the comparison uses the iterations until the energy stays
            # within 1e-8 Eh of the tight solution (a loose stop can end early while still 1e-7 off)
            r = solve_local(s, A, variant=var, e_tol=1e-10, r_tol=1e-8)
            errs = [abs(h - tight["energy"]) for h in r["history"]]
            for key, thr in (("it_to_tight_1e-6", 1e-6), ("it_to_tight_1e-8", 1e-8)):
                r[key] = next((k for k in range(len(errs)) if max(errs[k:]) < thr), None)
            r.update({"system": name, "space": space, "npairs": len(A.pairs), "m_mean": float(A.m.mean()),
                      "flops": fl, "e_tight": tight["energy"], "it_tight": tight["iterations"],
                      "d_tight": r["energy"] - tight["energy"]})
            out.append(r)
            print(f"{name:9s} {var:24s} it={r['iterations']:3d} (to 1e-8: {r['it_to_tight_1e-8']}) dE={r['d_tight']:.1e}"
                  f" inner/ladder={r['flop_inner'] / fl['ladder']:.2f} inner/fock_all={r['flop_inner'] / fl['fock_all']:.2f}"
                  f" inner mean={r['inner_mean']:.1f} sub={r.get('sub_mean', 0):.1f} {r['t_s']:.0f}s", flush=True)
            save(f"schwarz_{space.replace('|', '_')}.json", out)


def main():
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    mode, args = args[0], args[1:]
    space = "pno1e-7"
    if "--space" in args:
        k = args.index("--space")
        space = args[k + 1]
        args = args[:k] + args[k + 2:]
    spaces = ACC_SPACES
    if "--spaces" in args:
        k = args.index("--spaces")
        spaces = args[k + 1].split(",")
        args = args[:k] + args[k + 2:]
    variants = SCH_VARIANTS
    if "--variants" in args:
        k = args.index("--variants")
        variants = args[k + 1].split(",")
        args = args[:k] + args[k + 2:]
    if mode == "check":
        run_check(args or ["alkane4", "water3", "acene1"])
    elif mode == "accuracy":
        run_accuracy(args or ACC_SYSTEMS, spaces)
    elif mode == "schwarz":
        run_schwarz(args or SCH_SYSTEMS, space, variants)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
