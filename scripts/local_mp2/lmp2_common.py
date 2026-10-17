#!/usr/bin/env python3
"""Shared code of the local-MP2 feasibility prototype (WP15/WP14, docs/LOCAL_MP2_FEASIBILITY_2026-10.md).

Python/PySCF only, nothing here touches the curcuma engine. Pieces:
  - geometries generated in code (alkane zig-zag, linear acenes) or cut from a library structure (water clusters)
  - RI-JK RHF/def2-SVP with curcuma's def2-universal-JKFIT data (the same Turbomole file curcuma reads),
    cached per system in the work directory
  - frozen-core valence occupied orbitals, localised (Pipek-Mezey, Boys or IBO)
  - RI factors B_{i mu}^Q = sum_P (i mu|P) [L^-T]_{PQ}, J = L L^T (Cholesky of the aux metric)
  - Boughton-Pulay orbital domains, projected atomic orbitals (PAOs), pair-specific pseudo-canonical virtuals
  - the LMP2 residual equations in pair-specific virtual spaces (Pulay-Saeboe form), solved in an
    AO-embedded way: tau_ij = V_ij T_ij V_ij^T, the inter-pair coupling sum_k F_ik tau_kj is one GEMM
  - PNO construction from semicanonical pair amplitudes

Equations (closed shell, localised occupied i,j,k; pair virtuals a,b of pair ij orthonormal and
pseudo-canonical, eps_a; F_ik the occupied Fock matrix in the localised basis):
  R_ij = K_ij + (eps_a + eps_b - F_ii - F_jj) T_ij
         - V_ij^T S [ sum_{k!=i} F_ik tau_kj + sum_{k!=j} F_kj tau_ik ] S V_ij        (= 0 at convergence)
  E    = sum_{i<=j} (2 - delta_ij) sum_ab T_ij^ab (2 K_ij^ab - K_ij^ba),  K_ij^ab = (i a|j b) (RI)
With complete pair spaces this is canonical RI-MP2 in a rotated occupied basis (checked numerically by the
scripts); with diagonal F_ik (k != i) dropped it is the semicanonical (non-iterative) local MP2.

Claude Generated (Oct 2026). AI-generated, machine-tested only.
"""
import os

# Thread counts must be fixed before numpy/OpenBLAS load, otherwise every process spawns one BLAS thread per core
# (several concurrent runs then oversubscribe the machine). LMP2_THREADS (default 24) sets OMP and BLAS threads.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, os.environ.get("LMP2_THREADS", "24"))

import hashlib
import json
import sys
import time

import numpy as np
import scipy.linalg
from pyscf import df, gto, lib, lo, scf

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
QMDIR = os.path.join(ROOT, "src/core/energy_calculators/qm_methods")
LIB = os.path.join(ROOT, "test_cases/structures")
LVAL = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4, "h": 5}
BOHR = 0.52917721092  # Angstrom per Bohr (PySCF's value)

lib.num_threads(int(os.environ.get("LMP2_THREADS", "24")))


def work_dir():
    """Directory for caches and json results (env LMP2_WORK, default ./local_mp2_work)."""
    d = os.environ.get("LMP2_WORK", os.path.join(os.getcwd(), "local_mp2_work"))
    os.makedirs(d, exist_ok=True)
    return d


# ----------------------------------------------------------------------------------------------------------
# Auxiliary basis: curcuma's Turbomole file (same parser as scripts/rimp2_reference.py)
# ----------------------------------------------------------------------------------------------------------
def parse_turbomole(path):
    """Turbomole basis file -> {Symbol: [[L, [exp, coef], ...], ...]} (file order)."""
    out, elem, lines = {}, None, open(path).read().splitlines()
    k = 0
    while k < len(lines):
        s = lines[k].strip()
        k += 1
        if not s or s.startswith("#") or s.startswith("$") or s == "*":
            continue
        tok = s.split()
        if tok[0].isalpha():
            elem = tok[0].capitalize()
            out[elem] = []
            continue
        if len(tok) == 2 and tok[0].isdigit() and tok[1] in LVAL:
            nprim, L = int(tok[0]), LVAL[tok[1]]
            sh = [L]
            for _ in range(nprim):
                e, c = lines[k].replace("D", "E").replace("d", "e").split()[:2]
                k += 1
                sh.append([float(e), float(c)])
            out[elem].append(sh)
    return out


AUX_JKFIT = parse_turbomole(os.path.join(QMDIR, "def2-universal-JKFIT.turbomole"))


# ----------------------------------------------------------------------------------------------------------
# Geometries (Angstrom). All generated deterministically; nothing is written to the structure library.
# ----------------------------------------------------------------------------------------------------------
def alkane(n, rcc=1.53, rch=1.09):
    """All-anti n-alkane C_nH_(2n+2): planar zig-zag backbone (CCC 109.47 deg), tetrahedral H, staggered."""
    th = np.deg2rad(109.4712206)
    dx, dy = rcc * np.sin(th / 2), rcc * np.cos(th / 2)
    C = [np.array([k * dx, (k % 2) * dy, 0.0]) for k in range(-1, n + 1)]  # incl. two virtual end neighbours
    atoms = [("C", C[k]) for k in range(1, n + 1)]
    ca, sa = np.cos(np.pi - th / 2), np.sin(np.pi - th / 2)  # H directions relative to the CCC bisector
    for k in range(1, n + 1):
        u1 = (C[k - 1] - C[k]) / rcc
        u2 = (C[k + 1] - C[k]) / rcc
        b = (u1 + u2) / np.linalg.norm(u1 + u2)
        for s in (1, -1):  # two out-of-plane H per carbon
            d = ca * b + s * sa * np.array([0, 0, 1.0])
            atoms.append(("H", C[k] + rch * d / np.linalg.norm(d)))
        if k == 1:
            atoms.append(("H", C[k] + rch * u1))  # H in place of the missing chain continuation
        if k == n:
            atoms.append(("H", C[k] + rch * u2))
    return [(s, tuple(map(float, x))) for s, x in atoms]


def acene(nring, rcc=1.40, rch=1.08):
    """Linear acene with nring fused benzene rings (1 benzene, 2 naphthalene, ... 5 pentacene), D2h, ideal."""
    pts = []
    for r in range(nring):
        cx = r * 2 * rcc * np.cos(np.pi / 6)
        for k in range(6):
            a = np.pi / 6 + k * np.pi / 3
            p = np.array([cx + rcc * np.cos(a), rcc * np.sin(a), 0.0])
            if not any(np.linalg.norm(p - q) < 0.1 for q in pts):
                pts.append(p)
    atoms = [("C", p) for p in pts]
    for p in pts:
        nb = [q for q in pts if 0.1 < np.linalg.norm(p - q) < rcc + 0.1]
        if len(nb) == 2:
            d = (p - nb[0]) + (p - nb[1])
            atoms.append(("H", p + rch * d / np.linalg.norm(d)))
    return [(s, tuple(map(float, x))) for s, x in atoms]


WATER_SOURCE = "bulk/peo201x2-water1500.gfnff-min.xyz"
WATER_SOURCE_SHA256 = "889ee11848a489754b6394d13e88c80563be18d87679ba03a640c63ed312f10d"


def water_cluster(n):
    """(H2O)_n cut from the library structure peo201x2-water1500.gfnff-min (GFN-FF minimum of a PEO/water mixture).

    Waters = O atoms with exactly two H within 1.2 A and no C within 1.6 A. Seed = the water whose 16 nearest
    waters (O-O) have the smallest maximum distance (most compact environment); the n-cluster is the seed plus
    its n-1 nearest waters by O-O distance, so the clusters are nested. Geometry of each water unchanged."""
    path = os.path.join(LIB, WATER_SOURCE)
    raw = open(path, "rb").read()
    if hashlib.sha256(raw).hexdigest() != WATER_SOURCE_SHA256:
        sys.exit("water source structure changed (sha256)")
    L = raw.decode().splitlines()
    nat = int(L[0])
    sym = np.array([l.split()[0] for l in L[2:2 + nat]])
    X = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + nat]])
    O, H, Cc = np.where(sym == "O")[0], np.where(sym == "H")[0], np.where(sym == "C")[0]
    waters = []
    for o in O:
        dh = np.linalg.norm(X[H] - X[o], axis=1)
        hs = H[dh < 1.2]
        if len(hs) != 2:
            continue
        if len(Cc) and np.min(np.linalg.norm(X[Cc] - X[o], axis=1)) < 1.6:
            continue
        waters.append((o, hs[0], hs[1]))
    W = np.array(waters)
    XO = X[W[:, 0]]
    d = np.linalg.norm(XO[:, None, :] - XO[None, :, :], axis=2)
    srt = np.sort(d, axis=1)[:, 16]
    seed = int(np.argmin(srt))
    order = np.argsort(d[seed], kind="stable")[:n]
    atoms = []
    for w in order:
        for a in W[w]:
            atoms.append((str(sym[a]), tuple(map(float, X[a]))))
    return atoms


SYSTEMS = {}
for _n in (4, 8, 12, 16, 20):
    SYSTEMS[f"alkane{_n}"] = ("alkanes", lambda n=_n: alkane(n))
for _n in (4, 8, 12, 16):
    SYSTEMS[f"water{_n}"] = ("water clusters", lambda n=_n: water_cluster(n))
for _n in (1, 2, 3, 4, 5):
    SYSTEMS[f"acene{_n}"] = ("acenes", lambda n=_n: acene(n))
LABEL = {"acene1": "benzene", "acene2": "naphthalene", "acene3": "anthracene", "acene4": "tetracene",
         "acene5": "pentacene"}
for _n in (4, 8, 12, 16, 20):
    LABEL[f"alkane{_n}"] = f"C{_n}H{2 * _n + 2}"
for _n in (4, 8, 12, 16):
    LABEL[f"water{_n}"] = f"(H2O){_n}"


def label(name):
    return LABEL.get(name, name)


# ----------------------------------------------------------------------------------------------------------
# SCF and system setup
# ----------------------------------------------------------------------------------------------------------
class System:
    """Everything a local-MP2 step needs for one molecule (atomic units inside)."""
    pass


def build_mol(atoms):
    return gto.M(atom=[[s, x] for s, x in atoms], basis="def2-svp", unit="Angstrom", cart=False, verbose=0)


def nfrozen(mol):
    return sum(0 if mol.atom_charge(a) <= 2 else 1 for a in range(mol.natm))


def setup(name, loc="pm", verbose=True):
    """RI-JK RHF/def2-SVP (JKFIT), frozen core, localised valence occupied orbitals. SCF cached in work_dir()."""
    t0 = time.time()
    atoms = SYSTEMS[name][1]()
    mol = build_mol(atoms)
    auxmol = df.make_auxmol(mol, AUX_JKFIT)
    cache = os.path.join(work_dir(), f"scf_{name}.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        C, eps, ehf = z["C"], z["eps"], float(z["ehf"])
    else:
        mf = scf.RHF(mol).density_fit(auxbasis=AUX_JKFIT)
        mf.conv_tol = 1e-10
        mf.conv_tol_grad = 1e-6
        mf.max_cycle = 200
        ehf = mf.kernel()
        if not mf.converged:
            sys.exit(f"{name}: SCF not converged")
        C, eps = mf.mo_coeff, mf.mo_energy
        np.savez(cache, C=C, eps=eps, ehf=ehf)
    s = System()
    s.name, s.atoms, s.mol, s.auxmol, s.ehf = name, atoms, mol, auxmol, ehf
    s.C, s.eps = C, eps
    s.n, s.naux = mol.nao_nr(), auxmol.nao_nr()
    s.nocc = mol.nelectron // 2
    s.nfc = nfrozen(mol)
    s.nval = s.nocc - s.nfc
    s.nvir = s.n - s.nocc
    s.S = mol.intor("int1e_ovlp")
    SC = s.S @ C
    s.F = (SC * eps) @ SC.T  # AO Fock matrix consistent with the converged orbitals
    s.Cocc, s.Cval, s.Cvir = C[:, :s.nocc], C[:, s.nfc:s.nocc], C[:, s.nocc:]
    s.eval_, s.evir = eps[s.nfc:s.nocc], eps[s.nocc:]
    s.ao_atom = _ao_atom(mol)
    s.aux_atom = _ao_atom(auxmol)
    s.xyz = mol.atom_coords()  # Bohr
    s.natm = mol.natm
    s.localise = lambda method: localise(s, method)
    if loc:
        localise(s, loc)
    if verbose:
        print(f"# {name}: {label(name)} natm={s.natm} nbf={s.n} naux={s.naux} nocc={s.nocc} nfc={s.nfc} "
              f"nvir={s.nvir} E_HF(RI-JK)={ehf:.10f} setup {time.time() - t0:.1f}s", flush=True)
    return s


def _ao_atom(mol):
    out = []
    for b in range(mol.nbas):
        nf = (2 * mol.bas_angular(b) + 1) * mol.bas_nctr(b)
        out += [mol.bas_atom(b)] * nf
    return np.array(out)


def localise(s, method):
    """Localise the frozen-core valence occupied space. Sets s.L (AO coeffs), s.U (canonical->local), s.Focc,
    s.centroid (Bohr). method: pm (Pipek-Mezey, PySCF default meta-Lowdin populations), boys, ibo, can."""
    cache = os.path.join(work_dir(), f"loc_{s.name}_{method}.npy")
    if method == "can":
        L = s.Cval.copy()
    elif os.path.exists(cache):
        L = np.load(cache)
    else:
        if method == "pm":
            m = lo.PM(s.mol, s.Cval)
            m.conv_tol = 1e-10
            L = m.kernel()
        elif method == "boys":
            m = lo.Boys(s.mol, s.Cval)
            m.conv_tol = 1e-10
            L = m.kernel()
        elif method == "ibo":
            L = lo.ibo.ibo(s.mol, s.Cval, locmethod="IBO")
        else:
            raise ValueError(method)
        np.save(cache, L)
    # deterministic order: by centroid x, then y, z; deterministic sign: largest coefficient positive
    s.locmethod = method
    r = s.mol.intor("int1e_r")
    cen = np.array([[L[:, i] @ r[c] @ L[:, i] for c in range(3)] for i in range(L.shape[1])])
    order = np.lexsort((np.round(cen[:, 2], 6), np.round(cen[:, 1], 6), np.round(cen[:, 0], 6)))
    L = L[:, order]
    cen = cen[order]
    for i in range(L.shape[1]):
        k = np.argmax(np.abs(L[:, i]))
        if L[k, i] < 0:
            L[:, i] *= -1
    s.L = L
    s.U = s.Cval.T @ s.S @ L  # L = Cval U
    s.Focc = (s.U.T * s.eval_) @ s.U
    s.centroid = cen
    r2 = s.mol.intor("int1e_r2")
    s.spread = np.sqrt(np.maximum(np.array([L[:, i] @ r2 @ L[:, i] for i in range(L.shape[1])])
                                  - np.sum(cen ** 2, axis=1), 0))
    return L


# ----------------------------------------------------------------------------------------------------------
# RI integrals
# ----------------------------------------------------------------------------------------------------------
def ri_metric(s):
    if not hasattr(s, "J"):
        s.J = s.auxmol.intor("int2c2e")
        s.Lj = np.linalg.cholesky(s.J)
    return s.J


def int3c(s):
    """(mu nu|P) as array (n, n, naux)."""
    if not hasattr(s, "_int3c"):
        s._int3c = df.incore.aux_e2(s.mol, s.auxmol, intor="int3c2e", aosym="s1").reshape(s.n, s.n, s.naux)
    return s._int3c


def b_local_ao(s):
    """Fitted factors B[i, mu, Q] with localised valence i and AO mu: (i mu|j nu) = sum_Q B B."""
    ri_metric(s)
    I3 = int3c(s)
    X = np.einsum("mi,mnP->inP", s.L, I3, optimize=True)  # (i mu|P)
    o = s.nval
    Y = scipy.linalg.solve_triangular(s.Lj, X.reshape(-1, s.naux).T, lower=True)  # L^-1 (P|i mu)
    return np.ascontiguousarray(Y.T).reshape(o, s.n, s.naux)


def b_canonical(s, B_loc_ao=None):
    """B[i', a, Q] for canonical valence i' and canonical virtual a (full RI-MP2)."""
    if B_loc_ao is None:
        B_loc_ao = b_local_ao(s)
    Bl = np.einsum("imQ,ma->iaQ", B_loc_ao, s.Cvir, optimize=True)  # local i, canonical a
    return np.einsum("ji,iaQ->jaQ", s.U, Bl, optimize=True)  # back to canonical i'


def rimp2_canonical(s, Bc):
    """Full frozen-core RI-MP2 correlation energy (canonical formula)."""
    eo, ev = s.eval_, s.evir
    e = 0.0
    for i in range(Bc.shape[0]):
        K = np.einsum("aQ,jbQ->jab", Bc[i], Bc, optimize=True)
        D = eo[i] + eo[:, None, None] - ev[None, :, None] - ev[None, None, :]
        e += np.sum(K * (2 * K - K.transpose(0, 2, 1)) / D)
    return e


# ----------------------------------------------------------------------------------------------------------
# Domains and pair virtual spaces
# ----------------------------------------------------------------------------------------------------------
def bp_domain(s, i, thr):
    """Boughton-Pulay domain of local orbital i: add atoms by decreasing Mulliken gross population until the
    least-squares fit of L_i with the AOs of the selected atoms misses less than thr: 1 - <L_i|P_A L_i> < thr."""
    c = s.L[:, i]
    Sc = s.S @ c
    pop = np.array([np.sum((c * Sc)[s.ao_atom == A]) for A in range(s.natm)])
    order = np.argsort(-pop, kind="stable")
    sel = []
    for A in order:
        sel.append(A)
        idx = np.where(np.isin(s.ao_atom, sel))[0]
        b = Sc[idx]
        f = 1.0 - b @ np.linalg.solve(s.S[np.ix_(idx, idx)], b)
        if f < thr:
            break
    return sorted(int(a) for a in sel)


def extend_domain(s, atoms, rext_bohr):
    """Add every atom within rext_bohr of any atom of the domain."""
    if rext_bohr <= 0:
        return atoms
    d = np.linalg.norm(s.xyz[:, None, :] - s.xyz[None, atoms, :], axis=2).min(axis=1)
    return sorted(set(atoms) | set(int(a) for a in np.where(d <= rext_bohr)[0]))


def paos(s):
    """Normalised projected atomic orbitals: columns of (1 - D_occ S), D_occ incl. core."""
    if not hasattr(s, "Cpao"):
        P = np.eye(s.n) - s.Cocc @ s.Cocc.T @ s.S
        nrm = np.sqrt(np.einsum("mi,mn,ni->i", P, s.S, P))
        s.Cpao = P / nrm
    return s.Cpao


def pair_space(s, atoms, lindep=1e-6, force_pao=False):
    """Orthonormal pseudo-canonical virtuals spanned by the PAOs of the atom set. Returns (V (n x m), eps (m)).
    For the full atom set the canonical virtuals are returned (same span) unless force_pao."""
    if len(atoms) == s.natm and not force_pao:
        return s.Cvir, s.evir
    P = paos(s)[:, np.isin(s.ao_atom, atoms)]
    Sd = P.T @ s.S @ P
    w, v = np.linalg.eigh(Sd)
    keep = w > lindep
    X = P @ (v[:, keep] / np.sqrt(w[keep]))
    e, q = np.linalg.eigh(X.T @ s.F @ X)
    return X @ q, e


def pair_list(s):
    o = s.nval
    return [(i, j) for i in range(o) for j in range(i, o)]


# ----------------------------------------------------------------------------------------------------------
# LMP2 solver (AO-embedded coupling), PNOs
# ----------------------------------------------------------------------------------------------------------
def pair_energy(T, K, i, j):
    return (1.0 if i == j else 2.0) * np.sum(T * (2 * K - K.T))


def pair_ints(Bao, V, i, j):
    """K_ij^ab = (i a|j b) in the pair virtual basis V."""
    Bi = V.T @ Bao[i]
    Bj = Bi if i == j else V.T @ Bao[j]
    return Bi @ Bj.T


class DIIS:
    """Pulay DIIS on flat amplitude vectors with the residual (preconditioned update) as error vector."""

    def __init__(self, nmax=8):
        self.nmax, self.x, self.e = nmax, [], []

    def update(self, x, err):
        self.x.append(x.copy())
        self.e.append(err.copy())
        if len(self.x) > self.nmax:
            self.x.pop(0)
            self.e.pop(0)
        n = len(self.x)
        if n < 2:
            return x
        B = -np.ones((n + 1, n + 1))
        B[n, n] = 0
        for a in range(n):
            for b in range(a, n):
                B[a, b] = B[b, a] = self.e[a] @ self.e[b]
        rhs = np.zeros(n + 1)
        rhs[n] = -1
        try:
            c = np.linalg.solve(B, rhs)[:n]
        except np.linalg.LinAlgError:
            return x
        return sum(ci * xi for ci, xi in zip(c, self.x))


def solve_lmp2(s, pairs, spaces, Ks, coupled=True, diis=True, e_tol=1e-9, r_tol=1e-7, maxit=100,
               verbose=False):
    """Iterate the LMP2 residual equations (module docstring). pairs: list of (i,j), i<=j; spaces[p] = (V, eps);
    Ks[p] = K_ij in that space. coupled=False gives the semicanonical (one-step) result.
    Returns dict(energy, iterations, T, pair_energies, history)."""
    o, n = s.nval, s.n
    Fd = np.diag(s.Focc)
    Foff = s.Focc - np.diag(Fd)
    den = [spaces[p][1][:, None] + spaces[p][1][None, :] - Fd[i] - Fd[j] for p, (i, j) in enumerate(pairs)]
    T = [-Ks[p] / den[p] for p in range(len(pairs))]
    SV = [s.S @ spaces[p][0] for p in range(len(pairs))]
    sizes = [t.size for t in T]
    offs = np.concatenate([[0], np.cumsum(sizes)])
    flat = lambda TT: np.concatenate([t.ravel() for t in TT])
    unflat = lambda x: [x[offs[p]:offs[p + 1]].reshape(T[p].shape) for p in range(len(pairs))]
    energy = lambda TT: sum(pair_energy(TT[p], Ks[p], i, j) for p, (i, j) in enumerate(pairs))
    e_old = energy(T)
    hist = [e_old]
    if not coupled:
        return {"energy": e_old, "iterations": 0, "T": T, "history": hist,
                "pair_energies": [pair_energy(T[p], Ks[p], i, j) for p, (i, j) in enumerate(pairs)]}
    dd = DIIS() if diis else None
    it = 0
    for it in range(1, maxit + 1):
        tau = np.zeros((o, o, n, n))
        for p, (i, j) in enumerate(pairs):
            V = spaces[p][0]
            t = V @ T[p] @ V.T
            tau[i, j] = t
            if i != j:
                tau[j, i] = t.T
        Y = (Foff @ tau.reshape(o, -1)).reshape(o, o, n, n)
        del tau
        R = []
        for p, (i, j) in enumerate(pairs):
            Cij = Y[i, j] + Y[j, i].T
            R.append(Ks[p] + den[p] * T[p] - SV[p].T @ Cij @ SV[p])
        del Y
        rmax = max((np.max(np.abs(r)) if r.size else 0.0) for r in R)
        Tn = [T[p] - R[p] / den[p] for p in range(len(pairs))]
        if dd is not None:
            x = dd.update(flat(Tn), flat([R[p] / den[p] for p in range(len(pairs))]))
            Tn = unflat(x)
        T = Tn
        e = energy(T)
        hist.append(e)
        if verbose:
            print(f"   it {it:3d} E={e:.12f} dE={e - e_old:.2e} max|R|={rmax:.2e}", flush=True)
        if abs(e - e_old) < e_tol and rmax < r_tol:
            break
        e_old = e
    return {"energy": e, "iterations": it, "T": T, "history": hist,
            "pair_energies": [pair_energy(T[p], Ks[p], i, j) for p, (i, j) in enumerate(pairs)]}


def make_pnos(space, K, Fii, Fjj, i, j, tcut):
    """PNOs of one pair from semicanonical amplitudes in the given space.
    T = -K/(e_a+e_b-F_ii-F_jj); Ttilde = (4T - 2T^T)/(1+delta_ij); D = Ttilde^T T + Ttilde T^T (symmetric);
    keep eigenvectors with occupation >= tcut, re-canonicalise. Returns (V_pno, eps_pno, K_pno, dE_trunc, npno)."""
    V, e = space
    den = e[:, None] + e[None, :] - Fii - Fjj
    T = -K / den
    Tt = (4 * T - 2 * T.T) / (2.0 if i == j else 1.0)
    D = Tt.T @ T + Tt @ T.T
    w, Q = np.linalg.eigh(0.5 * (D + D.T))
    keep = w >= tcut
    Q = Q[:, keep]
    if Q.shape[1] == 0:
        return V[:, :0], e[:0], K[:0, :0], pair_energy(T, K, i, j), 0
    ep, W = np.linalg.eigh((Q.T * e) @ Q)
    QW = Q @ W
    Kp = QW.T @ K @ QW
    Tp = -Kp / (ep[:, None] + ep[None, :] - Fii - Fjj)
    dE = pair_energy(T, K, i, j) - pair_energy(Tp, Kp, i, j)
    return V @ QW, ep, Kp, dE, Q.shape[1]


def save_json(fname, obj):
    """Write results to work_dir(); env LMP2_TAG (optional) is appended to the name, e.g. q2_domains_alkanes.json."""
    tag = os.environ.get("LMP2_TAG", "")
    if tag:
        fname = fname.replace(".json", f"_{tag}.json")
    path = os.path.join(work_dir(), fname)
    json.dump(obj, open(path, "w"), indent=1)
    print("# wrote", path)
