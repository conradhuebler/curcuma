#!/usr/bin/env python3
"""Q4 of the local-MP2 feasibility study (WP14 prerequisite): locality of the inverse RI metric and local fit domains.

Aux basis def2-universal-JKFIT (curcuma's file). Three measurements per system:
  (a) decay of |J^-1_PQ| and of the Coulomb-normalised inverse (Jn = D^-1/2 J D^-1/2, D = diag J; elements of Jn^-1)
      with the distance of the atoms carrying P and Q (max per distance bin, relative to the largest element)
  (b) RI-J Coulomb energy E_J = 1/2 (rho|rho) of the SCF density with local fits: rho = sum_A rho_A,
      rho_A = sum_{mu on A, nu} D_mu,nu mu nu; rho_A is fitted with the aux functions on atoms within R of A,
      exact inverse of that local metric block; robust (Dunlap) energy (rho|c) - 1/2 c J c and plain 1/2 (rho|c)
  (c) RI-MP2 (frozen core) with local fits: (ia|P) for Pipek-Mezey orbital i and canonical virtual a is fitted with
      the aux functions on atoms within R of the centroid of i; robust (ia|jb) = (ia|P)c_jb + c_ia(P|jb) - c_ia J c_jb;
      the factors are rotated back to canonical occupied orbitals (linear) and the canonical energy formula is used,
      so only the fit changes. Reference: global fit with the same aux basis.
  (d) timing (numpy/scipy, this machine, n=1): Cholesky of the global metric vs all local blocks of (b) and (c).
Falsifier (vault/WP14): |error| < 1e-6 Eh/atom at R <= 8 Bohr over >= 3 systems and 3 sizes.

Usage: LMP2_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_mp2/q4_metric.py [systems...]
Claude Generated (Oct 2026).
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp2_common as lc  # first: sets the thread environment before numpy loads
import numpy as np
import scipy.linalg


DEFAULT = ["alkane4", "alkane8", "alkane12", "alkane16", "alkane20", "water4", "water8", "water12", "water16",
           "acene1", "acene2", "acene3", "acene4", "acene5"]
RADII = [4, 6, 8, 10, 12, 16, 20]  # Bohr
BINS = [0, 0.01, 3, 5, 8, 12, 16, 20, 30, 50]


def local_solve(J, idx, rhs):
    t = time.perf_counter()
    cf = scipy.linalg.cho_factor(J[np.ix_(idx, idx)], lower=True)
    t_f = time.perf_counter() - t
    return scipy.linalg.cho_solve(cf, rhs), t_f


def main():
    names = sys.argv[1:] or DEFAULT
    out = []
    for name in names:
        s = lc.setup(name, loc="pm")
        J = lc.ri_metric(s)
        naux = s.naux
        t_glob = 1e9
        for _ in range(3):  # min of 3 (the first call pays thread start-up)
            t = time.perf_counter()
            cfJ = scipy.linalg.cho_factor(J, lower=True)
            t_glob = min(t_glob, time.perf_counter() - t)
        rec = {"system": name, "class": lc.SYSTEMS[name][0], "natm": s.natm, "nbf": s.n, "naux": naux,
               "t_global_chol_s": t_glob}
        ev = np.linalg.eigvalsh(J)
        rec["cond_J"] = float(ev[-1] / ev[0])
        # (a) inverse decay
        Jinv = scipy.linalg.cho_solve(cfJ, np.eye(naux))
        dJ = np.sqrt(np.diag(J))
        Jn_inv = Jinv * dJ[:, None] * dJ[None, :]  # inverse of the Coulomb-normalised metric
        A = s.aux_atom
        Dat = np.linalg.norm(s.xyz[:, None, :] - s.xyz[None, :, :], axis=2)
        dPQ = Dat[np.ix_(A, A)]
        dec = []
        for lo_, hi in zip(BINS[:-1], BINS[1:]):
            m = (dPQ >= lo_) & (dPQ < hi)
            if m.any():
                dec.append([lo_, hi, float(np.max(np.abs(Jinv[m])) / np.max(np.abs(Jinv))),
                            float(np.max(np.abs(Jn_inv[m])) / np.max(np.abs(Jn_inv))),
                            float(np.max(np.abs(J[m])) / np.max(np.abs(J)))])
        rec["inv_decay"] = dec
        del Jinv, Jn_inv, dPQ
        # (b) RI-J energy
        I3 = lc.int3c(s)
        Docc = 2 * s.Cocc @ s.Cocc.T
        rho_A = np.zeros((s.natm, naux))
        for a in range(s.natm):
            mu = np.where(s.ao_atom == a)[0]
            rho_A[a] = np.einsum("mn,mnP->P", Docc[mu], I3[mu], optimize=True)
        rho = rho_A.sum(axis=0)
        c_glob = scipy.linalg.cho_solve(cfJ, rho)
        EJ = 0.5 * rho @ c_glob
        rec["E_J_global"] = float(EJ)
        rec["rij"] = {}
        for R in RADII:
            c = np.zeros(naux)
            tl, nd = 0.0, []
            for a in range(s.natm):
                idx = np.where(Dat[a][A] <= R)[0]
                ca, tf = local_solve(J, idx, rho_A[a][idx])
                c[idx] += ca
                tl += tf
                nd.append(len(idx))
            e_rob = rho @ c - 0.5 * c @ J @ c
            e_pl = 0.5 * rho @ c
            rec["rij"][str(R)] = {"err_robust": float(e_rob - EJ), "err_plain": float(e_pl - EJ),
                                  "err_robust_per_atom": float((e_rob - EJ) / s.natm), "ndom_mean": float(np.mean(nd)),
                                  "ndom_max": int(max(nd)), "t_local_chol_s": tl}
        # (c) RI-MP2 with local fits
        Araw_l = np.einsum("mi,na,mnP->iaP", s.L, s.Cvir, I3, optimize=True)  # (ia|P), local i
        Araw = np.einsum("Ii,iaP->IaP", s.U, Araw_l, optimize=True)  # canonical i'
        Bglob = scipy.linalg.solve_triangular(s.Lj, Araw.reshape(-1, naux).T, lower=True).T.reshape(Araw.shape)
        e_glob = lc.rimp2_canonical(s, Bglob)
        del Bglob
        rec["E_mp2_global_fit"] = float(e_glob)
        rec["rimp2"] = {}
        eo, evv = s.eval_, s.evir
        for R in RADII:
            Cl = np.zeros_like(Araw_l)
            tl, nd = 0.0, []
            for i in range(s.nval):
                idx = np.where(np.linalg.norm(s.xyz - s.centroid[i], axis=1)[A] <= R)[0]
                ci, tf = local_solve(J, idx, Araw_l[i][:, idx].T)
                Cl[i][:, idx] = ci.T
                tl += tf
                nd.append(len(idx))
            Cc = np.einsum("Ii,iaP->IaP", s.U, Cl, optimize=True)
            G = Cc @ J  # (o, v, naux)
            AmG = Araw - G
            e = 0.0
            for i in range(s.nval):
                K = np.einsum("aP,jbP->jab", AmG[i], Cc, optimize=True) + np.einsum("aP,jbP->jab", Cc[i], Araw,
                                                                                     optimize=True)
                D = eo[i] + eo[:, None, None] - evv[None, :, None] - evv[None, None, :]
                e += np.sum(K * (2 * K - K.transpose(0, 2, 1)) / D)
            rec["rimp2"][str(R)] = {"err": float(e - e_glob), "err_per_atom": float((e - e_glob) / s.natm),
                                    "ndom_mean": float(np.mean(nd)), "ndom_max": int(max(nd)), "t_local_chol_s": tl}
            del Cl, Cc, G, AmG
        out.append(rec)
        r8, m8 = rec["rij"]["8"], rec["rimp2"]["8"]
        print(f"{name:9s} naux={naux} cond={rec['cond_J']:.1e} chol {t_glob * 1e3:.1f} ms | RI-J R=8: "
              f"{r8['err_robust_per_atom']:.1e}/atom ndom={r8['ndom_mean']:.0f} | RI-MP2 R=8: {m8['err_per_atom']:.1e}/atom"
              f" ndom={m8['ndom_mean']:.0f}", flush=True)
        lc.save_json("q4_metric.json", out)
        s._int3c = None
        del I3, Araw, Araw_l


if __name__ == "__main__":
    main()
