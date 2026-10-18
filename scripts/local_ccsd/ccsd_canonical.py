#!/usr/bin/env python3
"""Stage A of the local-CCSD prototype: a simple closed-shell CCSD (RCCSD), dense numpy, plus its validation.

docs/LOCAL_CCSD_PROTOTYPE_2026-10.md has the results. Nothing here touches the curcuma engine.

Reference and integrals: PySCF RHF/def2-SVP (spherical) with exact four-index integrals (no RI), conv_tol 1e-12,
the same path as the local-MP2 study (scripts/local_mp2/). curcuma's own SCF orbitals are not exported by the binary,
so the bridge to curcuma is a comparison of E(HF) and E(MP2) with `curcuma -sp -method mp2 -qm.frozen_core true` on
the identical geometry (column "curcuma" of the validation table). Frozen core = chemical core (1s of Li-Ne).

Equations: spin-adapted RCCSD in the formulation of Hirata et al., J. Chem. Phys. 120, 2581 (2004), Eqs. 35-45,
transcribed from PySCF's pyscf/cc/rintermediates.py and rccsd.update_amps (chemists' notation integrals,
ovov[i,a,j,b] = (ia|jb), ...). Two changes against that transcription:
  - the code returns the residual R = N - e * t (N = PySCF's numerator, e_ia = f_ii - f_aa, e_ijab = e_ia + e_jb),
    so the same function serves non-canonical (localised) occupied orbitals: only diag(f) is moved to the
    denominator, the off-diagonal f_ik stays in Foo/Loo. R = 0 is the CCSD equation in any occupied basis.
  - the t1-dressed ladder W_abcd tau_ijcd is split into the bare (ac|bd) GEMM plus an o^3 v^3 term, so no second
    v^4 array is built.
Energy E = 2 sum f_ia t_ia + sum_ijab tau_ijab [2 (ia|jb) - (ib|ja)], tau = t2 + t1 t1.
Solver: Jacobi update t <- t - R/D (D = -e > 0) with Pulay DIIS (8 vectors) on (t1, t2), start t1 = 0,
t2 = MP2 amplitudes, so the iteration-0 energy is the MP2 energy (the MP2-limit check).

Validation (main): for each molecule E_corr(this code) against PySCF cc.CCSD (conv_tol 1e-12), ORCA 6.1 CCSD
(`! RHF CCSD def2-SVP ExtremeSCF`, %mdci STol 1e-10), the iteration-0 energy against PySCF mp.MP2 and curcuma's
native MP2, one residual against PySCF rccsd.update_amps with random amplitudes, and the occupied-rotation
invariance (Pipek-Mezey orbitals give the same E_corr).

Usage: LCCSD_WORK=<dir> ~/.venvs/qmref/bin/python -I scripts/local_ccsd/ccsd_canonical.py [names...]
       env ORCA=<orca binary> (default the 6.1 install), CURCUMA=<binary> (default release/curcuma)
Claude Generated (Oct 2026). AI-generated, machine-tested only.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "local_mp2"))
os.environ.setdefault("LMP2_WORK", os.environ.get("LCCSD_WORK", os.path.join(os.getcwd(), "local_ccsd_work")))
import lmp2_common as lc  # first: sets the thread environment before numpy loads

import json
import re
import subprocess
import time

import numpy as np
from pyscf import ao2mo, cc, gto, lib, mp, scf

ROOT = lc.ROOT
ORCA = os.environ.get("ORCA", "/home/conrad/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2/orca")
CURCUMA = os.environ.get("CURCUMA", os.path.join(ROOT, "release", "curcuma"))


def work_dir():
    return lc.work_dir()


# ----------------------------------------------------------------------------------------------------------
# Molecules
# ----------------------------------------------------------------------------------------------------------
def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    return [(a.split()[0], tuple(map(float, a.split()[1:4]))) for a in L[2:2 + n]]


SMALL = {  # validation set of stage A (test_cases files are read, not copied)
    "H2O": lambda: read_xyz(os.path.join(ROOT, "test_cases/qm_1e/H2O.xyz")),
    "NH3": lambda: read_xyz(os.path.join(ROOT, "test_cases/qm_1e/NH3.xyz")),
    "CH4": lambda: read_xyz(os.path.join(ROOT, "test_cases/qm_1e/CH4.xyz")),
    "formaldehyde": lambda: read_xyz(os.path.join(ROOT, "test_cases/qm_hf3c/formaldehyde.xyz")),
    "ethane": lambda: lc.alkane(2),
    "water_dimer": lambda: read_xyz(os.path.join(ROOT, "test_cases/qm_hf3c/water_dimer.xyz")),
}


def geometry(name):
    """Atoms (symbol, (x,y,z) Angstrom) for a stage-A name or a stage-B name (alkane<n>, water<n>, acene<n>)."""
    if name in SMALL:
        return SMALL[name]()
    for pre, fn in (("alkane", lc.alkane), ("water", lc.water_cluster), ("acene", lc.acene)):
        if name.startswith(pre):
            return fn(int(name[len(pre):]))
    raise KeyError(name)


def label(name):
    if name.startswith("alkane"):
        n = int(name[6:])
        return f"C{n}H{2 * n + 2}"
    if name.startswith("water"):
        return f"(H2O){name[5:]}"
    return lc.LABEL.get(name, name)


# ----------------------------------------------------------------------------------------------------------
# RHF and MO integrals
# ----------------------------------------------------------------------------------------------------------
def rhf(name, frozen_core=True):
    """Exact-integral RHF/def2-SVP, cached in the work directory. Returns a System-like object (lmp2_common names)."""
    atoms = geometry(name)
    mol = lc.build_mol(atoms)
    cache = os.path.join(work_dir(), f"rhf_exact_{name}.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        C, eps, ehf = z["C"], z["eps"], float(z["ehf"])
    else:
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.conv_tol_grad = 1e-8
        mf.max_cycle = 300
        ehf = mf.kernel()
        if not mf.converged:  # e.g. anthracene stalls at conv_tol_grad 1e-8 with DIIS: finish with second-order SCF
            mf2 = scf.RHF(mol).newton()
            mf2.conv_tol, mf2.conv_tol_grad, mf2.max_cycle = 1e-12, 1e-8, 100
            ehf = mf2.kernel(mf.mo_coeff, mf.mo_occ)
            mf = mf2
        if not mf.converged:
            sys.exit(f"{name}: RHF not converged")
        C, eps = mf.mo_coeff, mf.mo_energy
        np.savez(cache, C=C, eps=eps, ehf=ehf)
    s = lc.System()
    s.name, s.atoms, s.mol, s.ehf, s.C, s.eps = "exact_" + name, atoms, mol, ehf, C, eps
    s.n = mol.nao_nr()
    s.nocc = mol.nelectron // 2
    s.nfc = lc.nfrozen(mol) if frozen_core else 0
    s.nval, s.nvir = s.nocc - s.nfc, s.n - s.nocc
    s.S = mol.intor("int1e_ovlp")
    SC = s.S @ C
    s.F = (SC * eps) @ SC.T
    s.Cocc, s.Cval, s.Cvir = C[:, :s.nocc], C[:, s.nfc:s.nocc], C[:, s.nocc:]
    s.eval_, s.evir = eps[s.nfc:s.nocc], eps[s.nocc:]
    s.ao_atom = lc._ao_atom(mol)
    s.xyz = mol.atom_coords()
    s.natm = mol.natm
    return s


class MOInts:
    """Chemists'-notation MO integral blocks over active occupied (columns of Co) and virtual (Cv) orbitals.
    vvvv is stored in physicists' order V[a,b,c,d] = (ac|bd) as a (v^2, v^2) matrix for the ladder GEMM."""

    def __init__(self, mol, Co, Cv, Fao, eri_ao=None):
        t0 = time.time()
        o, v = Co.shape[1], Cv.shape[1]
        self.o, self.v = o, v
        if eri_ao is None:
            eri_ao = mol.intor("int2e", aosym="s8")
        g = lambda *cs: ao2mo.incore.general(eri_ao, cs, compact=False).reshape([c.shape[1] for c in cs])
        self.oooo = g(Co, Co, Co, Co)
        self.ovoo = g(Co, Cv, Co, Co)
        self.ovov = g(Co, Cv, Co, Cv)
        self.oovv = g(Co, Co, Cv, Cv)
        self.ovvo = g(Co, Cv, Cv, Co)
        self.ovvv = g(Co, Cv, Cv, Cv)
        X = ao2mo.incore.general(eri_ao, (Cv, Cv, Cv, Cv), compact=False).reshape(v, v, v, v)  # (ac|bd) [a,c,b,d]
        V = np.empty((v, v, v, v))
        for a in range(v):
            V[a] = X[a].transpose(1, 0, 2)  # [b,c,d]
        del X
        self.vvvv = V.reshape(v * v, v * v)
        self.foo = Co.T @ Fao @ Co
        self.fvv = Cv.T @ Fao @ Cv
        self.fov = Co.T @ Fao @ Cv
        self.t_build = time.time() - t0

    def rotate_occ(self, U):
        """New object with occupied orbitals Co U (U orthogonal o x o); virtual blocks shared."""
        r = MOInts.__new__(MOInts)
        r.o, r.v, r.vvvv, r.fvv, r.t_build = self.o, self.v, self.vvvv, self.fvv, 0.0
        r.oooo = lib.einsum("pi,qj,rk,sl,pqrs->ijkl", U, U, U, U, self.oooo)
        r.ovoo = lib.einsum("pi,rk,sl,pars->iakl", U, U, U, self.ovoo)
        r.ovov = lib.einsum("pi,rj,parb->iajb", U, U, self.ovov)
        r.oovv = lib.einsum("pi,qj,pqab->ijab", U, U, self.oovv)
        r.ovvo = lib.einsum("pi,sj,pabs->iabj", U, U, self.ovvo)
        r.ovvv = lib.einsum("pi,pabc->iabc", U, self.ovvv)
        r.foo = U.T @ self.foo @ U
        r.fov = U.T @ self.fov
        return r


# ----------------------------------------------------------------------------------------------------------
# RCCSD residual (Hirata et al. 2004 via PySCF rintermediates), energy, solver
# ----------------------------------------------------------------------------------------------------------
def denominators(I):
    eo, ev = np.diag(I.foo), np.diag(I.fvv)
    eia = eo[:, None] - ev[None, :]
    return eia, eia[:, None, :, None] + eia[None, :, None, :]


def energy(t1, t2, I):
    tau = t2 + np.einsum("ia,jb->ijab", t1, t1)
    return (2 * np.einsum("ia,ia", I.fov, t1) + 2 * np.einsum("ijab,iajb", tau, I.ovov)
            - np.einsum("ijab,ibja", tau, I.ovov))


def residual(t1, t2, I, mp2_only=False):
    """CCSD residuals (R1, R2) = (N1 - e_ia t1, N2 - e_ijab t2); zero at the solution in any occupied basis.
    mp2_only: the (non-canonical) MP2 residual K + f_vv t - f_oo t terms only (for local MP2 at the same truncation)."""
    o, v = I.o, I.v
    eia, eijab = denominators(I)
    fov, foo, fvv = I.fov, I.foo, I.fvv
    ovov, oovv, ovvo, ovoo, ovvv = I.ovov, I.oovv, I.ovvo, I.ovoo, I.ovvv
    if mp2_only:
        Loo = foo - np.diag(np.diag(foo))
        Lvv = fvv - np.diag(np.diag(fvv))
        N2 = ovov.transpose(0, 2, 1, 3).copy()
        tmp = lib.einsum("ac,ijcb->ijab", Lvv, t2) - lib.einsum("ki,kjab->ijab", Loo, t2)
        N2 += tmp + tmp.transpose(1, 0, 3, 2)
        return np.zeros_like(t1), N2 - eijab * t2

    # kappa intermediates (Eqs. 37-39)
    Foo = (2 * lib.einsum("kcld,ilcd->ki", ovov, t2) - lib.einsum("kdlc,ilcd->ki", ovov, t2)
           + 2 * lib.einsum("kcld,ic,ld->ki", ovov, t1, t1) - lib.einsum("kdlc,ic,ld->ki", ovov, t1, t1) + foo)
    Fvv = (-2 * lib.einsum("kcld,klad->ac", ovov, t2) + lib.einsum("kdlc,klad->ac", ovov, t2)
           - 2 * lib.einsum("kcld,ka,ld->ac", ovov, t1, t1) + lib.einsum("kdlc,ka,ld->ac", ovov, t1, t1) + fvv)
    Fov = 2 * np.einsum("kcld,ld->kc", ovov, t1) - np.einsum("kdlc,ld->kc", ovov, t1) + fov
    # lambda intermediates (Eqs. 40-41), built before the diagonal is removed from Foo/Fvv
    Loo = Foo + np.einsum("kc,ic->ki", fov, t1) + 2 * np.einsum("lcki,lc->ki", ovoo, t1) - np.einsum("kcli,lc->ki", ovoo, t1)
    Lvv = Fvv - np.einsum("kc,ka->ac", fov, t1) + 2 * np.einsum("kdac,kd->ac", ovvv, t1) - np.einsum("kcad,kd->ac", ovvv, t1)
    Foo -= np.diag(np.diag(foo))
    Fvv -= np.diag(np.diag(fvv))
    Loo -= np.diag(np.diag(foo))
    Lvv -= np.diag(np.diag(fvv))

    # T1 numerator
    N1 = -2 * np.einsum("kc,ka,ic->ia", fov, t1, t1)
    N1 += np.einsum("ac,ic->ia", Fvv, t1)
    N1 -= np.einsum("ki,ka->ia", Foo, t1)
    N1 += 2 * np.einsum("kc,kica->ia", Fov, t2) - np.einsum("kc,ikca->ia", Fov, t2)
    N1 += np.einsum("kc,ic,ka->ia", Fov, t1, t1)
    N1 += fov
    N1 += 2 * np.einsum("kcai,kc->ia", ovvo, t1) - np.einsum("kiac,kc->ia", oovv, t1)
    N1 += 2 * lib.einsum("kdac,ikcd->ia", ovvv, t2) - lib.einsum("kcad,ikcd->ia", ovvv, t2)
    N1 += 2 * lib.einsum("kdac,kd,ic->ia", ovvv, t1, t1) - lib.einsum("kcad,kd,ic->ia", ovvv, t1, t1)
    N1 += -2 * lib.einsum("lcki,klac->ia", ovoo, t2) + lib.einsum("kcli,klac->ia", ovoo, t2)
    N1 += -2 * lib.einsum("lcki,lc,ka->ia", ovoo, t1, t1) + lib.einsum("kcli,lc,ka->ia", ovoo, t1, t1)

    # T2 numerator
    tmp2 = lib.einsum("kibc,ka->abic", oovv, -t1) + ovvv.transpose(1, 3, 0, 2)
    tmp = lib.einsum("abic,jc->ijab", tmp2, t1)
    N2 = tmp + tmp.transpose(1, 0, 3, 2)
    tmp2 = lib.einsum("kcai,jc->akij", ovvo, t1) + ovoo.transpose(1, 3, 0, 2)
    tmp = lib.einsum("akij,kb->ijab", tmp2, t1)
    N2 -= tmp + tmp.transpose(1, 0, 3, 2)
    del tmp2
    N2 += ovov.transpose(0, 2, 1, 3)

    # chi intermediates (Eqs. 42-45)
    Woooo = (lib.einsum("lcki,jc->klij", ovoo, t1) + lib.einsum("kclj,ic->klij", ovoo, t1)
             + lib.einsum("kcld,ijcd->klij", ovov, t2) + lib.einsum("kcld,ic,jd->klij", ovov, t1, t1)
             + I.oooo.transpose(0, 2, 1, 3))
    Wvoov = (lib.einsum("kcad,id->akic", ovvv, t1) - lib.einsum("kcli,la->akic", ovoo, t1)
             + ovvo.transpose(2, 0, 3, 1) - 0.5 * lib.einsum("ldkc,ilda->akic", ovov, t2)
             - 0.5 * lib.einsum("lckd,ilad->akic", ovov, t2) - lib.einsum("ldkc,id,la->akic", ovov, t1, t1)
             + lib.einsum("ldkc,ilad->akic", ovov, t2))
    Wvovo = (lib.einsum("kdac,id->akci", ovvv, t1) - lib.einsum("lcki,la->akci", ovoo, t1)
             + oovv.transpose(2, 0, 3, 1) - 0.5 * lib.einsum("lckd,ilda->akci", ovov, t2)
             - lib.einsum("lckd,id,la->akci", ovov, t1, t1))

    tau = t2 + np.einsum("ia,jb->ijab", t1, t1)
    N2 += lib.einsum("klij,klab->ijab", Woooo, tau)
    # ladder: sum_cd W_abcd tau_ijcd with W_abcd = (ac|bd) - sum_k (kd|ac) t_kb - sum_k (kc|bd) t_ka
    N2 += (tau.reshape(o * o, v * v) @ I.vvvv).reshape(o, o, v, v)
    Y = lib.einsum("kdac,ijcd->ijka", ovvv, tau)
    tmp = -lib.einsum("kb,ijka->ijab", t1, Y)
    N2 += tmp + tmp.transpose(1, 0, 3, 2)
    del Y
    tmp = lib.einsum("ac,ijcb->ijab", Lvv, t2)
    N2 += tmp + tmp.transpose(1, 0, 3, 2)
    tmp = lib.einsum("ki,kjab->ijab", Loo, t2)
    N2 -= tmp + tmp.transpose(1, 0, 3, 2)
    tmp = 2 * lib.einsum("akic,kjcb->ijab", Wvoov, t2) - lib.einsum("akci,kjcb->ijab", Wvovo, t2)
    N2 += tmp + tmp.transpose(1, 0, 3, 2)
    tmp = lib.einsum("akic,kjbc->ijab", Wvoov, t2)
    N2 -= tmp + tmp.transpose(1, 0, 3, 2)
    tmp = lib.einsum("bkci,kjac->ijab", Wvovo, t2)
    N2 -= tmp + tmp.transpose(1, 0, 3, 2)
    return N1 - eia * t1, N2 - eijab * t2


def mp2_amplitudes(I):
    eia, eijab = denominators(I)
    return np.zeros((I.o, I.v)), I.ovov.transpose(0, 2, 1, 3) / eijab


def solve_ccsd(I, e_tol=1e-11, r_tol=1e-9, maxit=100, verbose=False, mp2_only=False):
    """Dense RCCSD with Jacobi update and DIIS. In a non-canonical occupied basis the start amplitudes are not
    the MP2 solution, but the converged result is the same. Returns dict(e_corr, e_iter0, iterations, t1, t2)."""
    eia, eijab = denominators(I)
    t1, t2 = mp2_amplitudes(I)
    e0 = energy(t1, t2, I)
    e_old, hist = e0, [e0]
    dd = lc.DIIS(8)
    o, v = I.o, I.v
    t0 = time.time()
    for it in range(1, maxit + 1):
        R1, R2 = residual(t1, t2, I, mp2_only=mp2_only)
        rmax = max(np.max(np.abs(R1)), np.max(np.abs(R2)))
        d1, d2 = R1 / eia, R2 / eijab  # = -R/D: the Jacobi step
        x = dd.update(np.concatenate([(t1 + d1).ravel(), (t2 + d2).ravel()]),
                      np.concatenate([d1.ravel(), d2.ravel()]))
        t1, t2 = x[:o * v].reshape(o, v), x[o * v:].reshape(o, o, v, v)
        e = energy(t1, t2, I)
        hist.append(e)
        if verbose:
            print(f"   it {it:3d} E={e:.12f} dE={e - e_old:.2e} max|R|={rmax:.2e} {time.time() - t0:.1f}s", flush=True)
        if abs(e - e_old) < e_tol and rmax < r_tol:
            break
        e_old = e
    return {"e_corr": float(e), "e_iter0": float(e0), "iterations": it, "t1": t1, "t2": t2, "history": hist,
            "t_s": time.time() - t0}


# ----------------------------------------------------------------------------------------------------------
# External references
# ----------------------------------------------------------------------------------------------------------
def pyscf_refs(s):
    mf = scf.RHF(s.mol)
    mf.mo_coeff, mf.mo_energy, mf.e_tot = s.C, s.eps, s.ehf
    mf.mo_occ = np.array([2.0] * s.nocc + [0.0] * s.nvir)
    mf.converged = True
    frozen = s.nfc if s.nfc else None
    e_mp2 = mp.MP2(mf, frozen=frozen).kernel()[0]
    mycc = cc.CCSD(mf, frozen=frozen)
    mycc.conv_tol, mycc.conv_tol_normt, mycc.max_cycle, mycc.verbose = 1e-12, 1e-10, 200, 0
    mycc.kernel()
    return {"e_mp2": float(e_mp2), "e_ccsd": float(mycc.e_corr), "ccsd_converged": bool(mycc.converged),
            "ccsd_iterations": None}, mf


def check_update_amps(s, mf, I):
    """Max deviation of one residual against PySCF rccsd.update_amps on random amplitudes (canonical orbitals)."""
    from pyscf.cc import rccsd
    mycc = rccsd.RCCSD(mf, frozen=s.nfc if s.nfc else None)
    eris = mycc.ao2mo()
    rng = np.random.default_rng(7)
    t1 = 0.02 * rng.standard_normal((I.o, I.v))
    t2 = 0.02 * rng.standard_normal((I.o, I.o, I.v, I.v))
    t2 = 0.5 * (t2 + t2.transpose(1, 0, 3, 2))
    p1, p2 = rccsd.update_amps(mycc, t1, t2, eris)
    # same Fock matrix as PySCF (eris.fock is rebuilt from the density, differs from diag(eps) by ~1e-10)
    I2 = MOInts.__new__(MOInts)
    I2.__dict__.update(I.__dict__)
    o = I.o
    I2.foo, I2.fvv, I2.fov = eris.fock[:o, :o], eris.fock[o:, o:], eris.fock[:o, o:]
    eia = eris.mo_energy[:o, None] - eris.mo_energy[None, o:]
    eijab = eia[:, None, :, None] + eia[None, :, None, :]
    R1, R2 = residual(t1, t2, I2)
    # PySCF: t_new = N / e  ->  R = e (t_new - t)
    return float(max(np.max(np.abs(eia * (p1 - t1) - R1)), np.max(np.abs(eijab * (p2 - t2) - R2))))


def orca_ccsd(name, atoms):
    d = os.path.join(work_dir(), "orca_ccsd", name)
    os.makedirs(d, exist_ok=True)
    outp = os.path.join(d, "job.out")
    if not (os.path.exists(outp) and "ORCA TERMINATED NORMALLY" in open(outp).read()):
        with open(os.path.join(d, "job.inp"), "w") as f:
            f.write("! RHF CCSD def2-SVP ExtremeSCF\n%mdci\n  STol 1e-10\n  MaxIter 200\nend\n%maxcore 8000\n"
                    "* xyz 0 1\n")
            for sym, x in atoms:
                f.write(f"{sym} {x[0]:.10f} {x[1]:.10f} {x[2]:.10f}\n")
            f.write("*\n")
        env = dict(os.environ, LD_LIBRARY_PATH=os.path.dirname(ORCA))
        with open(outp, "w") as fo:
            subprocess.run([ORCA, "job.inp"], cwd=d, stdout=fo, stderr=subprocess.STDOUT, env=env, check=False)
    txt = open(outp).read()
    if "ORCA TERMINATED NORMALLY" not in txt:
        return None
    ehf, ecorr = orca_property_energies(os.path.join(d, "job.property.txt"))
    block = txt.split("COUPLED CLUSTER ITERATIONS")[-1].split("have converged")[0]
    nit = len(re.findall(r"^\s*\d+\s+-\d+\.\d+\s+-\d+\.\d+", block, re.M)) - 1  # row 0 = MP2 guess
    return {"e_hf": ehf, "e_ccsd": ecorr, "iterations": nit, "out": outp}


def orca_property_energies(path):
    """(reference energy, correlation energy) at full precision from ORCA's job.property.txt (job.out prints 9 dp)."""
    txt = open(path).read()
    val = lambda key: float(re.findall(key + r"[^\n]*\n[^\n]*\n\s*\n\s*0\s+(-?\d\.\d+e[-+]\d+)", txt)[-1])
    return val("&refEnergy"), val("&corrEnergy")


def curcuma_mp2(name, atoms):
    d = os.path.join(work_dir(), "curcuma_mp2", name)
    os.makedirs(d, exist_ok=True)
    xyz = os.path.join(d, f"{name}.xyz")
    with open(xyz, "w") as f:
        f.write(f"{len(atoms)}\n{name}\n")
        for sym, x in atoms:
            f.write(f"{sym} {x[0]:.10f} {x[1]:.10f} {x[2]:.10f}\n")
    r = subprocess.run([CURCUMA, "-sp", xyz, "-method", "mp2", "-qm.frozen_core", "true", "-no_bmt"], cwd=d,
                       capture_output=True, text=True)
    txt = re.sub(r"\x1b\[[0-9;]*m", "", r.stdout)
    ehf = re.findall(r"E\(HF\)\s+=\s+(-?\d+\.\d+)", txt)
    ec = re.findall(r"E\(corr, mp2\s*\)\s+=\s+(-?\d+\.\d+)", txt)
    if not ehf or not ec:
        return None
    return {"e_hf": float(ehf[-1]), "e_mp2": float(ec[-1])}


def main():
    names = sys.argv[1:] or ["H2O", "NH3", "CH4", "formaldehyde", "ethane", "water_dimer"]
    out = []
    for name in names:
        s = rhf(name)
        I = MOInts(s.mol, s.Cval, s.Cvir, s.F)
        r = solve_ccsd(I)
        ref, mf = pyscf_refs(s)
        dupd = check_update_amps(s, mf, I)
        lc.localise(s, "pm")
        rl = solve_ccsd(I.rotate_occ(s.U))
        orc = orca_ccsd(name, s.atoms)
        cur = curcuma_mp2(name, s.atoms)
        rec = {"system": name, "nbf": s.n, "o": s.nval, "v": s.nvir, "nfc": s.nfc, "e_hf": s.ehf,
               "e_ccsd": r["e_corr"], "iterations": r["iterations"], "e_iter0": r["e_iter0"],
               "pyscf_ccsd": ref["e_ccsd"], "pyscf_mp2": ref["e_mp2"], "d_update_amps": dupd,
               "e_ccsd_pm": rl["e_corr"], "iterations_pm": rl["iterations"], "orca": orc, "curcuma": cur}
        out.append(rec)
        print(f"{name:13s} nbf={s.n:3d} o={s.nval} v={s.nvir} E_corr={r['e_corr']:.12f} it={r['iterations']}"
              f"  -PySCF={r['e_corr'] - ref['e_ccsd']:.1e}"
              f"  -ORCA={(r['e_corr'] - orc['e_ccsd']) if orc else float('nan'):.1e}"
              f"  E0-MP2(PySCF)={r['e_iter0'] - ref['e_mp2']:.1e}"
              f"  E0-MP2(curcuma)={(r['e_iter0'] - cur['e_mp2']) if cur else float('nan'):.1e}"
              f"  HF-curcuma={(s.ehf - cur['e_hf']) if cur else float('nan'):.1e}"
              f"  HF-ORCA={(s.ehf - orc['e_hf']) if orc else float('nan'):.1e}"
              f"  dUpdAmps={dupd:.1e}  PM-can={rl['e_corr'] - r['e_corr']:.1e} (it {rl['iterations']})", flush=True)
        json.dump(out, open(os.path.join(work_dir(), "stageA_validation.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
