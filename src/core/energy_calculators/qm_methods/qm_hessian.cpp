/*
 * <Native QM: analytic closed-shell RHF Hessian (second nuclear derivatives)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested only; see docs/QM_HF_HESSIAN.md for what was
 * compared with what.
 *
 * ---------------------------------------------------------------------------------------------------------
 * Theory (closed-shell RHF, real orbitals, spin-summed density P = 2 C_o C_o^T)
 * ---------------------------------------------------------------------------------------------------------
 *
 *   E = sum P_mn h_mn + 1/2 sum D_mnls (mn|ls) + E_nn,   D_mnls = P_mn P_ls - 1/4 (P_ml P_ns + P_ms P_nl)
 *
 * Gradient (WP8, qm_integrals.cpp):   E^x = sum P h^x + 1/2 sum D (mn|ls)^x - sum W S^x + E_nn^x,
 *   W = 2 sum_i eps_i C_i C_i^T = 1/2 P F P   (energy-weighted density).
 *
 * Differentiating the gradient once more (Pople, Krishnan, Schlegel, Binkley, Int. J. Quantum Chem. Symp. 13,
 * 225 (1979); Yamaguchi, Osamura, Goddard, Schaefer, "A New Dimension to Quantum Chemistry: Analytic Derivative
 * Methods in Ab Initio Molecular Electronic Structure Theory", Oxford 1994, ch. 4 and 10; Helgaker, Jorgensen,
 * Olsen, Molecular Electronic-Structure Theory, ch. 9 for the integrals):
 *
 *   H_xy = [ sum P h^xy + 1/2 sum D (mn|ls)^xy - sum W S^xy + E_nn^xy ]          (fixed density, "skeleton")
 *        + Tr( P^y F^x ) - Tr( W^y S^x )                                         (response)
 *
 *   F^x  = h^x + G^x[P],  G^x[P]_mn = sum_ls [ (mn|ls)^x - 1/2 (ml|ns)^x ] P_ls  (derivative Fock, P fixed)
 *   P^y  = dP/dy = -1/2 P S^y P + 2 ( C_v U^y C_o^T + C_o U^y^T C_v^T )
 *   W^y  = 1/2 ( P^y F P + P (F^y + G[P^y]) P + P F P^y ),  G[D] = J[D] - 1/2 K[D]
 *
 * The orbital response U^y (virtual x occupied) solves the coupled-perturbed HF equations, obtained by
 * differentiating the canonical condition (C^T F C)_ai = 0 with C^y = C U^y and U_ia + U_ai + S^y_ai = 0
 * (the occupied-occupied and virtual-virtual blocks are fixed at U_pq = -1/2 S^y_pq, which is what the first
 * term of P^y represents):
 *
 *   (eps_a - eps_i) U_ai + [ C_v^T G[ D(U) ] C_o ]_ai = -B^y_ai,     D(U) = 2 ( C_v U C_o^T + C_o U^T C_v^T )
 *   B^y_ai = F^y_ai - eps_i S^y_ai + [ C_v^T G[ -1/2 P S^y P ] C_o ]_ai         (MO blocks of AO matrices)
 *
 * The operator on the left is the RHF orbital Hessian A_ai,bj = (eps_a - eps_i) delta + 4 (ai|bj) - (ab|ij)
 * - (aj|bi), symmetric and positive definite at a stable RHF minimum, so it is solved column by column with a
 * preconditioned conjugate gradient (preconditioner 1/(eps_a - eps_i)); every iteration needs one J/K build of
 * a symmetric response density, done by the engine's own exact J/K builder (stored tensor or integral-direct).
 * No orbital-energy derivatives are needed because W is written in the invariant form 1/2 P F P.
 *
 * Derivative integrals: every derivative of a cartesian Gaussian with respect to its centre is a combination of
 * Gaussians with shifted powers, d/dA_x [x_A^l e^{-a r_A^2}] = 2a x_A^(l+1) e^.. - l x_A^(l-1) e^.., applied once
 * or twice per centre. All integrals are McMurchie-Davidson (Helgaker 9.5, 9.9): Hermite coefficients E^{ij}_t
 * of each Gaussian pair, Hermite Coulomb integrals R_tuv from the Boys function. The centres of an integral are
 * independent variables; the derivative with respect to the last centre (the nucleus C of V, the fourth shell D
 * of an ERI) follows from translational invariance, d/dC = -(d/dA + d/dB), applied once or twice. An atom's
 * derivative is the sum over the centres that sit on it.
 *
 * Basis: all derivative integrals are formed in the cartesian basis; densities go in as Q P Q^T and derivative
 * matrices come back as Q^T M Q, Q the geometry-independent cartesian -> spherical map (as the WP8 gradient).
 *
 * Units: everything in atomic units; the Hessian is Eh/Bohr^2 (mass independent, cartesian). The capability
 * (src/capabilities/hessian.cpp) converts to its internal Eh/Angstrom^2.
 *
 * Not implemented: open shell (UHF), Kohn-Sham, ECP, density fitting, g and higher shells, GPU, screening of
 * the derivative quartets (every canonical quartet with more than one centre is computed).
 *
 * This program is free software under GPL-3.0
 */

#include "qm_hessian.h"
#include "qm_engine.h"

#include "src/core/curcuma_logger.h"
#include "src/core/units.h"

#include <fmt/format.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace qmhess {

namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;

// ---------------------------------------------------------------------------------------------------------
// Shells (consecutive cartesian components with the same centre, exponents and L; as qmint::groupShells)
// ---------------------------------------------------------------------------------------------------------
struct HShell {
    int first = 0, L = 0, atom = 0;
    double x = 0, y = 0, z = 0;
    std::vector<double> exps;
    std::vector<std::array<int, 3>> lmn;
    std::vector<std::vector<double>> coef;  ///< [component][primitive], pre-normalised (basissetparser)
};

std::vector<HShell> groupShells(const std::vector<GTO::Orbital>& basis)
{
    std::vector<HShell> shells;
    for (size_t i = 0; i < basis.size(); ++i) {
        int l, m, k;
        GTO::orbitalTypeToComponents(basis[i].type, l, m, k);
        const int L = l + m + k;
        bool join = false;
        if (!shells.empty()) {
            const HShell& s = shells.back();
            const GTO::Orbital& prev = basis[i - 1];
            join = s.L == L && prev.atom == basis[i].atom && prev.exponents == basis[i].exponents && prev.x == basis[i].x
                && prev.y == basis[i].y && prev.z == basis[i].z && (int)s.lmn.size() < (L + 1) * (L + 2) / 2;
        }
        if (!join) {
            HShell s;
            s.first = (int)i;
            s.L = L;
            s.atom = basis[i].atom;
            s.x = basis[i].x;
            s.y = basis[i].y;
            s.z = basis[i].z;
            s.exps = basis[i].exponents;
            shells.push_back(std::move(s));
        }
        shells.back().lmn.push_back({ l, m, k });
        shells.back().coef.push_back(basis[i].coefficients);
    }
    return shells;
}

// ---------------------------------------------------------------------------------------------------------
// Hermite simplex: index of (t,u,v), t+u+v <= L
// ---------------------------------------------------------------------------------------------------------
struct Simplex {
    int L = 0;
    std::vector<std::array<int, 3>> tri;
};

const Simplex& simplex(int L)
{
    static const std::vector<Simplex> table = [] {
        std::vector<Simplex> t(25);
        for (int l = 0; l < 25; ++l) {
            t[l].L = l;
            for (int a = 0; a <= l; ++a)
                for (int b = 0; b <= l - a; ++b)
                    for (int c = 0; c <= l - a - b; ++c) t[l].tri.push_back({ a, b, c });
        }
        return t;
    }();
    return table[L];
}

// ---------------------------------------------------------------------------------------------------------
// 1D Hermite expansion coefficients E^{ij}_t (Helgaker 9.5.6/9.5.7), E^{00}_0 = 1 (K factor kept apart)
// ---------------------------------------------------------------------------------------------------------
struct Herm1D {
    int imax = 0, jmax = 0, tdim = 0;
    std::vector<double> e;
    double operator()(int t, int i, int j) const
    {
        if (t > i + j) return 0.0;
        return e[((size_t)t * (imax + 1) + i) * (jmax + 1) + j];
    }
    double& at(int t, int i, int j) { return e[((size_t)t * (imax + 1) + i) * (jmax + 1) + j]; }
};

void hermite1D(int imax, int jmax, double XPA, double XPB, double p, Herm1D& E)
{
    E.imax = imax;
    E.jmax = jmax;
    E.tdim = imax + jmax + 1;
    E.e.assign((size_t)E.tdim * (imax + 1) * (jmax + 1), 0.0);
    E.at(0, 0, 0) = 1.0;
    const double oo2p = 0.5 / p;
    for (int i = 0; i <= imax; ++i)
        for (int j = 0; j <= jmax; ++j) {
            if (i == 0 && j == 0) continue;
            // raise i from (i-1, j) when i > 0, else raise j from (0, j-1)
            const int pi = i > 0 ? i - 1 : i, pj = i > 0 ? j : j - 1;
            const double X = i > 0 ? XPA : XPB;
            for (int t = 0; t <= i + j; ++t) {
                double v = X * E(t, pi, pj);
                if (t > 0) v += oo2p * E(t - 1, pi, pj);
                if (t + 1 <= pi + pj) v += (t + 1) * E(t + 1, pi, pj);
                E.at(t, i, j) = v;
            }
        }
}

/// Hermite Coulomb integrals R^0_tuv(alpha, X, Y, Z) for t+u+v <= L (Helgaker 9.9.18-9.9.20), dense (L+1)^3.
void hermiteR(int L, double alpha, double X, double Y, double Z, double scale, std::vector<double>& R0)
{
    thread_local std::vector<double> F, Rn;
    const int D = L + 1;
    F.resize(D + 1);
    qmint::boysFunction(L, alpha * (X * X + Y * Y + Z * Z), F);
    Rn.assign((size_t)D * D * D * D, 0.0);
    auto idx = [D](int n, int t, int u, int v) { return (((size_t)n * D + t) * D + u) * D + v; };
    double m2a = 1.0;
    for (int n = 0; n <= L; ++n) {
        Rn[idx(n, 0, 0, 0)] = m2a * F[n];
        m2a *= -2.0 * alpha;
    }
    for (int n = L - 1; n >= 0; --n)
        for (int N = 1; N <= L - n; ++N)
            for (int t = 0; t <= N; ++t)
                for (int u = 0; u <= N - t; ++u) {
                    const int v = N - t - u;
                    double val;
                    if (t > 0)
                        val = X * Rn[idx(n + 1, t - 1, u, v)] + (t > 1 ? (t - 1) * Rn[idx(n + 1, t - 2, u, v)] : 0.0);
                    else if (u > 0)
                        val = Y * Rn[idx(n + 1, t, u - 1, v)] + (u > 1 ? (u - 1) * Rn[idx(n + 1, t, u - 2, v)] : 0.0);
                    else
                        val = Z * Rn[idx(n + 1, t, u, v - 1)] + (v > 1 ? (v - 1) * Rn[idx(n + 1, t, u, v - 2)] : 0.0);
                    Rn[idx(n, t, u, v)] = val;
                }
    R0.assign((size_t)D * D * D, 0.0);
    for (int t = 0; t <= L; ++t)
        for (int u = 0; u <= L - t; ++u)
            for (int v = 0; v <= L - t - u; ++v) R0[((size_t)t * D + u) * D + v] = scale * Rn[idx(0, t, u, v)];
}

// ---------------------------------------------------------------------------------------------------------
// Centre derivatives of one 1D Gaussian factor: (d/dA)^n x_A^l e^{-a x_A^2} = sum_k c_k x_A^{l_k} e^{..}
// ---------------------------------------------------------------------------------------------------------
struct Terms1D {
    int n = 0;
    double c[16];
    int l[16];
};

void derivTerms1D(int l, double alpha, int nder, Terms1D& out)
{
    out.n = 1;
    out.c[0] = 1.0;
    out.l[0] = l;
    for (int d = 0; d < nder; ++d) {
        Terms1D nw;
        nw.n = 0;
        for (int k = 0; k < out.n; ++k) {
            // merge equal powers to keep the list short
            auto add = [&nw](double c, int p) {
                for (int q = 0; q < nw.n; ++q)
                    if (nw.l[q] == p) {
                        nw.c[q] += c;
                        return;
                    }
                nw.c[nw.n] = c;
                nw.l[nw.n] = p;
                ++nw.n;
            };
            add(2.0 * alpha * out.c[k], out.l[k] + 1);
            if (out.l[k] > 0) add(-out.l[k] * out.c[k], out.l[k] - 1);
        }
        out = nw;
    }
}

/// e[t] = sum over derivative terms of the pair (nA derivatives on the A factor, nB on the B factor) of E_t.
void pair1D(const Herm1D& E, int la, int lb, double alpha, double beta, int nA, int nB, int tmax, double* e)
{
    Terms1D ta, tb;
    derivTerms1D(la, alpha, nA, ta);
    derivTerms1D(lb, beta, nB, tb);
    for (int t = 0; t <= tmax; ++t) {
        double s = 0.0;
        for (int i = 0; i < ta.n; ++i)
            for (int j = 0; j < tb.n; ++j) s += ta.c[i] * tb.c[j] * E(t, ta.l[i], tb.l[j]);
        e[t] = s;
    }
}

/// Derivative variants of a Gaussian pair (bra (A,B) or ket (C,D)): derivative counts per direction on the
/// first (A) and second (B) centre. 0 none; 1-3 A_k; 4-6 B_k; 7-12 A_k A_l (k <= l); 13-21 A_k B_l; 22-27 B_k B_l.
struct Variant {
    int nA[3], nB[3];
};
constexpr int kSym[6][2] = { { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 1 }, { 1, 2 }, { 2, 2 } };
inline int symIndex(int k, int l)
{
    if (k > l) std::swap(k, l);
    static const int m[3][3] = { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } };
    return m[k][l];
}
const std::vector<Variant>& variants()
{
    static const std::vector<Variant> v = [] {
        std::vector<Variant> out;
        auto mk = [](std::array<int, 3> a, std::array<int, 3> b) {
            Variant x;
            for (int d = 0; d < 3; ++d) {
                x.nA[d] = a[d];
                x.nB[d] = b[d];
            }
            return x;
        };
        auto unit = [](int k) {
            std::array<int, 3> e = { 0, 0, 0 };
            e[k] = 1;
            return e;
        };
        auto two = [](int k, int l) {
            std::array<int, 3> e = { 0, 0, 0 };
            e[k] += 1;
            e[l] += 1;
            return e;
        };
        out.push_back(mk({ 0, 0, 0 }, { 0, 0, 0 }));
        for (int k = 0; k < 3; ++k) out.push_back(mk(unit(k), { 0, 0, 0 }));
        for (int k = 0; k < 3; ++k) out.push_back(mk({ 0, 0, 0 }, unit(k)));
        for (int s = 0; s < 6; ++s) out.push_back(mk(two(kSym[s][0], kSym[s][1]), { 0, 0, 0 }));
        for (int k = 0; k < 3; ++k)
            for (int l = 0; l < 3; ++l) out.push_back(mk(unit(k), unit(l)));
        for (int s = 0; s < 6; ++s) out.push_back(mk({ 0, 0, 0 }, two(kSym[s][0], kSym[s][1])));
        return out;
    }();
    return v;
}
constexpr int vA(int k) { return 1 + k; }
constexpr int vB(int k) { return 4 + k; }
inline int vAA(int k, int l) { return 7 + symIndex(k, l); }
constexpr int vAB(int k, int l) { return 13 + 3 * k + l; }
inline int vBB(int k, int l) { return 22 + symIndex(k, l); }
/// ket tables keep variants 0, A_k, A_k A_l (the ket's first shell is centre C): 10 of them
constexpr int kKetVariants[10] = { 0, 1, 2, 3, 7, 8, 9, 10, 11, 12 };

/// Hermite vector of a component pair for one variant over simplex(Lh): out[h] = scale ex[t] ey[u] ez[v].
void pairVector(const Herm1D& Ex, const Herm1D& Ey, const Herm1D& Ez, const std::array<int, 3>& la,
                const std::array<int, 3>& lb, double alpha, double beta, const Variant& v, int Lh, double scale,
                double* out)
{
    double ex[32], ey[32], ez[32];
    pair1D(Ex, la[0], lb[0], alpha, beta, v.nA[0], v.nB[0], Lh, ex);
    pair1D(Ey, la[1], lb[1], alpha, beta, v.nA[1], v.nB[1], Lh, ey);
    pair1D(Ez, la[2], lb[2], alpha, beta, v.nA[2], v.nB[2], Lh, ez);
    const Simplex& S = simplex(Lh);
    for (size_t h = 0; h < S.tri.size(); ++h) out[h] = scale * ex[S.tri[h][0]] * ey[S.tri[h][1]] * ez[S.tri[h][2]];
}

/// Overlap-type value of a pair (Hermite index 0 only): ex[0] ey[0] ez[0] with derivative counts nA, nB.
double pairE0(const Herm1D& Ex, const Herm1D& Ey, const Herm1D& Ez, const std::array<int, 3>& la,
              const std::array<int, 3>& lb, double alpha, double beta, const int nA[3], const int nB[3])
{
    double e[3];
    pair1D(Ex, la[0], lb[0], alpha, beta, nA[0], nB[0], 0, &e[0]);
    pair1D(Ey, la[1], lb[1], alpha, beta, nA[1], nB[1], 0, &e[1]);
    pair1D(Ez, la[2], lb[2], alpha, beta, nA[2], nB[2], 0, &e[2]);
    return e[0] * e[1] * e[2];
}

// ---------------------------------------------------------------------------------------------------------
// Centre-pair second-derivative blocks and translational invariance
// ---------------------------------------------------------------------------------------------------------
using Mat3 = Eigen::Matrix3d;

/// Hc[X][Y] for X, Y < nind are given (independent centres); fill the dependent centre nind:
/// d/dZ = -sum_X d/dX  =>  H[X][Z] = -sum_Y H[X][Y],  H[Z][Y] = -sum_X H[X][Y],  H[Z][Z] = sum_XY H[X][Y].
void extendTranslational(int nind, Mat3 Hc[4][4])
{
    const int Z = nind;
    for (int X = 0; X < nind; ++X) {
        Mat3 s = Mat3::Zero();
        for (int Y = 0; Y < nind; ++Y) s += Hc[X][Y];
        Hc[X][Z] = -s;
    }
    for (int Y = 0; Y < nind; ++Y) {
        Mat3 s = Mat3::Zero();
        for (int X = 0; X < nind; ++X) s += Hc[X][Y];
        Hc[Z][Y] = -s;
    }
    Mat3 s = Mat3::Zero();
    for (int X = 0; X < nind; ++X)
        for (int Y = 0; Y < nind; ++Y) s += Hc[X][Y];
    Hc[Z][Z] = s;
}

void distribute(int ncent, const int atom[4], const Mat3 Hc[4][4], double w, Matrix& hess)
{
    for (int X = 0; X < ncent; ++X)
        for (int Y = 0; Y < ncent; ++Y)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l) hess(3 * atom[X] + k, 3 * atom[Y] + l) += w * Hc[X][Y](k, l);
}

/// Fixed work chunks (independent of the thread count, so the sums are bitwise reproducible).
constexpr int kChunks = 16;

/// Called first inside every parallel loop body: the small Eigen products of one worker must not open nested OpenMP
/// teams. With an outer team of one thread Eigen's GEMM otherwise forks omp_get_max_threads() threads per product
/// (measured: dimethyl ether def2-TZVP, -threads 1, > 45 min instead of minutes). Only this task's ICV changes.
inline void noNestedThreads()
{
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
}

}  // namespace

int maxAngularMomentum(const std::vector<GTO::Orbital>& basis)
{
    int L = 0;
    for (const auto& o : basis) {
        int l, m, n;
        GTO::orbitalTypeToComponents(o.type, l, m, n);
        L = std::max(L, l + m + n);
    }
    return L;
}

// =========================================================================================================
// One-electron part: S^x, (T+V)^x matrices and sum P h^xy - sum W S^xy
// =========================================================================================================
void oneElectronDerivatives(const std::vector<GTO::Orbital>& basis, const std::vector<int>& atomZ,
                            const Matrix& atomPosBohr, const Matrix& P, const Matrix& W, std::vector<Matrix>& S1,
                            std::vector<Matrix>& H1, Matrix& hess, int threads)
{
    const int n = (int)basis.size();
    const int natoms = (int)atomZ.size();
    const int n3 = 3 * natoms;
    const std::vector<HShell> sh = groupShells(basis);
    const int nsh = (int)sh.size();
    S1.assign(n3, Matrix::Zero(n, n));
    H1.assign(n3, Matrix::Zero(n, n));
    std::vector<Matrix> hchunk(kChunks, Matrix::Zero(n3, n3));
    const std::vector<Variant>& V = variants();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
#endif
    for (int ch = 0; ch < kChunks; ++ch) {
        noNestedThreads();
        Matrix& hloc = hchunk[ch];
        std::vector<double> vec, RC;
        // rows of S1/H1 that belong to shell A are written by the chunk owning A only (disjoint rows)
        for (int iA = ch; iA < nsh; iA += kChunks) {
            const HShell& A = sh[iA];
            for (int iB = 0; iB < nsh; ++iB) {
                const HShell& B = sh[iB];
                const bool same_atom = A.atom == B.atom;
                const int Lh = A.L + B.L + 2;
                const int nh = (int)simplex(Lh).tri.size();
                const int D = Lh + 1;
                const double ABx = A.x - B.x, ABy = A.y - B.y, ABz = A.z - B.z;
                const double rab2 = ABx * ABx + ABy * ABy + ABz * ABz;
                for (size_t ia = 0; ia < A.exps.size(); ++ia)
                    for (size_t ib = 0; ib < B.exps.size(); ++ib) {
                        const double al = A.exps[ia], be = B.exps[ib], p = al + be;
                        const double K = std::exp(-al * be / p * rab2);
                        if (K == 0.0) continue;
                        const double Px = (al * A.x + be * B.x) / p, Py = (al * A.y + be * B.y) / p,
                                     Pz = (al * A.z + be * B.z) / p;
                        Herm1D Ex, Ey, Ez;
                        hermite1D(A.L + 2, B.L + 4, Px - A.x, Px - B.x, p, Ex);
                        hermite1D(A.L + 2, B.L + 4, Py - A.y, Py - B.y, p, Ey);
                        hermite1D(A.L + 2, B.L + 4, Pz - A.z, Pz - B.z, p, Ez);
                        const double prefS = std::pow(kPi / p, 1.5) * K;
                        // R tensors of every nucleus for this primitive pair
                        std::vector<std::vector<double>> Rnuc(natoms);
                        for (int C = 0; C < natoms; ++C) {
                            if (same_atom && A.atom == C) continue;  // one-centre: total derivative is zero
                            hermiteR(Lh, p, Px - atomPosBohr(C, 0), Py - atomPosBohr(C, 1), Pz - atomPosBohr(C, 2),
                                     -atomZ[C] * 2.0 * kPi / p * K, Rnuc[C]);
                        }
                        vec.resize((size_t)28 * nh);
                        for (size_t ca = 0; ca < A.lmn.size(); ++ca)
                            for (size_t cb = 0; cb < B.lmn.size(); ++cb) {
                                const int mu = A.first + (int)ca, nu = B.first + (int)cb;
                                const double cc = A.coef[ca][ia] * B.coef[cb][ib];
                                if (cc == 0.0) continue;
                                const auto& la = A.lmn[ca];
                                const auto& lb = B.lmn[cb];
                                const double pmn = P(mu, nu), wmn = W(mu, nu);
                                // --- overlap and kinetic (two centres; skipped on one atom: zero total) ---
                                if (!same_atom) {
                                    auto S = [&](const Variant& v) { return prefS * pairE0(Ex, Ey, Ez, la, lb, al, be, v.nA, v.nB); };
                                    // T = -1/2 <a| nabla^2 |b> = -1/2 sum_d <a| (d/dB_d)^2 |b>
                                    auto T = [&](const Variant& v) {
                                        double s = 0.0;
                                        for (int d = 0; d < 3; ++d) {
                                            int nb[3] = { v.nB[0], v.nB[1], v.nB[2] };
                                            nb[d] += 2;
                                            s += pairE0(Ex, Ey, Ez, la, lb, al, be, v.nA, nb);
                                        }
                                        return -0.5 * prefS * s;
                                    };
                                    for (int k = 0; k < 3; ++k) {
                                        S1[3 * A.atom + k](mu, nu) += cc * S(V[vA(k)]);
                                        S1[3 * B.atom + k](mu, nu) += cc * S(V[vB(k)]);
                                        H1[3 * A.atom + k](mu, nu) += cc * T(V[vA(k)]);
                                        H1[3 * B.atom + k](mu, nu) += cc * T(V[vB(k)]);
                                    }
                                    Mat3 Hs[4][4], Ht[4][4];
                                    for (int k = 0; k < 3; ++k)
                                        for (int l = 0; l < 3; ++l) {
                                            Hs[0][0](k, l) = S(V[vAA(k, l)]);
                                            Hs[0][1](k, l) = S(V[vAB(k, l)]);
                                            Hs[1][1](k, l) = S(V[vBB(k, l)]);
                                            Ht[0][0](k, l) = T(V[vAA(k, l)]);
                                            Ht[0][1](k, l) = T(V[vAB(k, l)]);
                                            Ht[1][1](k, l) = T(V[vBB(k, l)]);
                                        }
                                    Hs[1][0] = Hs[0][1].transpose();
                                    Ht[1][0] = Ht[0][1].transpose();
                                    const int at2[4] = { A.atom, B.atom, 0, 0 };
                                    distribute(2, at2, Hs, -cc * wmn, hloc);
                                    distribute(2, at2, Ht, cc * pmn, hloc);
                                }
                                // --- nuclear attraction, nucleus C (third centre, d/dC by translational invariance) ---
                                for (int v = 0; v < 28; ++v)
                                    pairVector(Ex, Ey, Ez, la, lb, al, be, V[v], Lh, 1.0, &vec[(size_t)v * nh]);
                                const Simplex& Sx = simplex(Lh);
                                for (int C = 0; C < natoms; ++C) {
                                    if (Rnuc[C].empty()) continue;
                                    const std::vector<double>& R = Rnuc[C];
                                    auto val = [&](int v) {
                                        const double* e = &vec[(size_t)v * nh];
                                        double s = 0.0;
                                        for (int h = 0; h < nh; ++h) {
                                            const auto& t = Sx.tri[h];
                                            s += e[h] * R[((size_t)t[0] * D + t[1]) * D + t[2]];
                                        }
                                        return s;
                                    };
                                    for (int k = 0; k < 3; ++k) {
                                        const double dA = val(vA(k)), dB = val(vB(k));
                                        H1[3 * A.atom + k](mu, nu) += cc * dA;
                                        H1[3 * B.atom + k](mu, nu) += cc * dB;
                                        H1[3 * C + k](mu, nu) -= cc * (dA + dB);
                                    }
                                    Mat3 Hv[4][4];
                                    for (int k = 0; k < 3; ++k)
                                        for (int l = 0; l < 3; ++l) {
                                            Hv[0][0](k, l) = val(vAA(k, l));
                                            Hv[0][1](k, l) = val(vAB(k, l));
                                            Hv[1][1](k, l) = val(vBB(k, l));
                                        }
                                    Hv[1][0] = Hv[0][1].transpose();
                                    extendTranslational(2, Hv);
                                    const int at3[4] = { A.atom, B.atom, C, 0 };
                                    distribute(3, at3, Hv, cc * pmn, hloc);
                                }
                            }
                    }
            }
        }
    }
    hess = Matrix::Zero(n3, n3);
    for (int ch = 0; ch < kChunks; ++ch) hess += hchunk[ch];
}

// =========================================================================================================
// Two-electron part: derivative Fock matrices G^x[P] and 1/2 sum D (mn|ls)^xy, integral-direct
// =========================================================================================================
namespace {

/// One primitive pair of a shell pair: Gaussian product data, 1D Hermite tables raised by two on both shells,
/// and the ket-variant Hermite matrices EK[w] (ncomp_pair x simplex(LA+LB+2), contraction coefficients and K in).
struct HPrim {
    int ia = 0, ib = 0;
    double alpha = 0, beta = 0, p = 0, Px = 0, Py = 0, Pz = 0, K = 0;
    Herm1D Ex, Ey, Ez;
    std::vector<Eigen::MatrixXd> EK;  ///< 10 ket variants (kKetVariants)
};
struct HPair {
    int A = 0, B = 0, Lh = 0, nh = 0, ncomp = 0;
    std::vector<HPrim> prims;
};

HPair makePair(const std::vector<HShell>& sh, int iA, int iB, bool ket_tables)
{
    HPair pr;
    pr.A = iA;
    pr.B = iB;
    const HShell& A = sh[iA];
    const HShell& B = sh[iB];
    pr.Lh = A.L + B.L + 2;
    pr.nh = (int)simplex(pr.Lh).tri.size();
    pr.ncomp = (int)(A.lmn.size() * B.lmn.size());
    const double ABx = A.x - B.x, ABy = A.y - B.y, ABz = A.z - B.z;
    const double rab2 = ABx * ABx + ABy * ABy + ABz * ABz;
    const std::vector<Variant>& V = variants();
    for (size_t ia = 0; ia < A.exps.size(); ++ia)
        for (size_t ib = 0; ib < B.exps.size(); ++ib) {
            HPrim pp;
            pp.ia = (int)ia;
            pp.ib = (int)ib;
            pp.alpha = A.exps[ia];
            pp.beta = B.exps[ib];
            pp.p = pp.alpha + pp.beta;
            pp.K = std::exp(-pp.alpha * pp.beta / pp.p * rab2);
            if (pp.K == 0.0) continue;
            pp.Px = (pp.alpha * A.x + pp.beta * B.x) / pp.p;
            pp.Py = (pp.alpha * A.y + pp.beta * B.y) / pp.p;
            pp.Pz = (pp.alpha * A.z + pp.beta * B.z) / pp.p;
            hermite1D(A.L + 2, B.L + 2, pp.Px - A.x, pp.Px - B.x, pp.p, pp.Ex);
            hermite1D(A.L + 2, B.L + 2, pp.Py - A.y, pp.Py - B.y, pp.p, pp.Ey);
            hermite1D(A.L + 2, B.L + 2, pp.Pz - A.z, pp.Pz - B.z, pp.p, pp.Ez);
            if (ket_tables) {
                pp.EK.assign(10, Eigen::MatrixXd::Zero(pr.ncomp, pr.nh));
                std::vector<double> buf(pr.nh);
                for (int w = 0; w < 10; ++w)
                    for (size_t ca = 0; ca < A.lmn.size(); ++ca)
                        for (size_t cb = 0; cb < B.lmn.size(); ++cb) {
                            const double c = A.coef[ca][ia] * B.coef[cb][ib] * pp.K;
                            pairVector(pp.Ex, pp.Ey, pp.Ez, A.lmn[ca], B.lmn[cb], pp.alpha, pp.beta,
                                       V[kKetVariants[w]], pr.Lh, c, buf.data());
                            const int row = (int)(ca * B.lmn.size() + cb);
                            for (int h = 0; h < pr.nh; ++h) pp.EK[w](row, h) = buf[h];
                        }
            }
            pr.prims.push_back(std::move(pp));
        }
    return pr;
}

/// Bra-variant Hermite matrices (28 variants) of every primitive pair of a pair.
void braTables(const std::vector<HShell>& sh, const HPair& pr, std::vector<std::vector<Eigen::MatrixXd>>& EB)
{
    const HShell& A = sh[pr.A];
    const HShell& B = sh[pr.B];
    const std::vector<Variant>& V = variants();
    EB.resize(pr.prims.size());
    std::vector<double> buf(pr.nh);
    for (size_t q = 0; q < pr.prims.size(); ++q) {
        const HPrim& pp = pr.prims[q];
        EB[q].assign(28, Eigen::MatrixXd::Zero(pr.ncomp, pr.nh));
        for (int v = 0; v < 28; ++v)
            for (size_t ca = 0; ca < A.lmn.size(); ++ca)
                for (size_t cb = 0; cb < B.lmn.size(); ++cb) {
                    const double c = A.coef[ca][pp.ia] * B.coef[cb][pp.ib] * pp.K;
                    pairVector(pp.Ex, pp.Ey, pp.Ez, A.lmn[ca], B.lmn[cb], pp.alpha, pp.beta, V[v], pr.Lh, c,
                               buf.data());
                    const int row = (int)(ca * B.lmn.size() + cb);
                    for (int h = 0; h < pr.nh; ++h) EB[q][v](row, h) = buf[h];
                }
    }
}

/// Digest one block I[ab][cd] of a canonical shell quartet into G = J - 1/2 K of density P (one ordering, degeneracy
/// weights deg/4 for J and deg/8 for K); the caller symmetrises G <- G + G^T at the end (the other orderings).
void digest(const HShell& A, const HShell& B, const HShell& C, const HShell& Dd, const Eigen::MatrixXd& I, double deg,
            const Matrix& P, Matrix& G)
{
    const double fJ = deg / 4.0, fK = 0.5 * deg / 8.0;
    const int nB = (int)B.lmn.size(), nD = (int)Dd.lmn.size();
    for (int a = 0; a < (int)A.lmn.size(); ++a)
        for (int b = 0; b < nB; ++b) {
            const int m = A.first + a, nn = B.first + b;
            for (int c = 0; c < (int)C.lmn.size(); ++c)
                for (int d = 0; d < nD; ++d) {
                    const double v = I(a * nB + b, c * nD + d);
                    if (v == 0.0) continue;
                    const int l = C.first + c, s = Dd.first + d;
                    G(m, nn) += fJ * v * P(l, s);
                    G(l, s) += fJ * v * P(m, nn);
                    G(m, l) -= fK * v * P(nn, s);
                    G(nn, s) -= fK * v * P(m, l);
                    G(m, s) -= fK * v * P(nn, l);
                    G(nn, l) -= fK * v * P(m, s);
                }
        }
}

/// Shared driver: canonical shell quartets, per primitive quartet the Hermite products. mode 0 = derivative
/// Fock + second derivatives (the Hessian), mode 1 = undifferentiated J - K/2 (validation of the machinery).
void twoElectronDriver(const std::vector<GTO::Orbital>& basis, const Matrix& P, int natoms, int threads, int mode,
                       std::vector<Matrix>& G1, Matrix& hess, Matrix& G0)
{
    const int n = (int)basis.size();
    const int n3 = 3 * natoms;
    const std::vector<HShell> sh = groupShells(basis);
    const int nsh = (int)sh.size();
    std::vector<HPair> pairs;
    for (int i = 0; i < nsh; ++i)
        for (int j = 0; j <= i; ++j) pairs.push_back(makePair(sh, i, j, true));
    const int npair = (int)pairs.size();

    const int nmat = mode == 0 ? n3 : 1;
    std::vector<std::vector<Matrix>> gchunk(kChunks, std::vector<Matrix>(nmat, Matrix::Zero(n, n)));
    std::vector<Matrix> hchunk(kChunks, Matrix::Zero(mode == 0 ? n3 : 1, mode == 0 ? n3 : 1));
    const double pref0 = 2.0 * std::pow(kPi, 2.5);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
#endif
    for (int ch = 0; ch < kChunks; ++ch) {
        noNestedThreads();
        std::vector<Matrix>& gl = gchunk[ch];
        Matrix& hl = hchunk[ch];
        std::vector<std::vector<Eigen::MatrixXd>> EB;
        std::vector<double> R0;
        // pairs are assigned round-robin from the end (the later bra pairs carry more ket pairs)
        for (int ij = npair - 1 - ch; ij >= 0; ij -= kChunks) {
            const HPair& bra = pairs[ij];
            const HShell& A = sh[bra.A];
            const HShell& B = sh[bra.B];
            if (mode == 0) braTables(sh, bra, EB);
            for (int kl = 0; kl <= ij; ++kl) {
                const HPair& ket = pairs[kl];
                const HShell& C = sh[ket.A];
                const HShell& Dd = sh[ket.B];
                const int atoms[4] = { A.atom, B.atom, C.atom, Dd.atom };
                if (mode == 0 && atoms[0] == atoms[1] && atoms[0] == atoms[2] && atoms[0] == atoms[3]) continue;
                const double deg = (bra.A == bra.B ? 1.0 : 2.0) * (ket.A == ket.B ? 1.0 : 2.0) * (ij == kl ? 1.0 : 2.0);
                const int nab = bra.ncomp, ncd = ket.ncomp, nh = bra.nh, nk = ket.nh;
                const int Lmax = A.L + B.L + C.L + Dd.L + 2;  // highest Hermite order any variant pairs to
                const int DR = Lmax + 1;
                const Simplex& Sh = simplex(bra.Lh);
                const Simplex& Sk = simplex(ket.Lh);

                // Density weights of the second-derivative term: 1/2 deg D_abcd
                Eigen::MatrixXd Dq;
                if (mode == 0) {
                    Dq.resize(nab, ncd);
                    const int nB = (int)B.lmn.size(), nD = (int)Dd.lmn.size();
                    for (int a = 0; a < (int)A.lmn.size(); ++a)
                        for (int b = 0; b < nB; ++b)
                            for (int c = 0; c < (int)C.lmn.size(); ++c)
                                for (int d = 0; d < nD; ++d) {
                                    const int m = A.first + a, nn = B.first + b, l = C.first + c, s = Dd.first + d;
                                    Dq(a * nB + b, c * nD + d)
                                        = 0.5 * deg * (P(m, nn) * P(l, s) - 0.25 * (P(m, l) * P(nn, s) + P(m, s) * P(nn, l)));
                                }
                }
                // Contracted first-derivative blocks for centres A (0), B (1), C (2), per direction
                Eigen::MatrixXd I1[3][3], I0;
                if (mode == 0)
                    for (int X = 0; X < 3; ++X)
                        for (int k = 0; k < 3; ++k) I1[X][k] = Eigen::MatrixXd::Zero(nab, ncd);
                else
                    I0 = Eigen::MatrixXd::Zero(nab, ncd);
                Mat3 Hc[4][4];
                for (int X = 0; X < 4; ++X)
                    for (int Y = 0; Y < 4; ++Y) Hc[X][Y].setZero();

                Eigen::MatrixXd Rmat(nh, nk), M[10], Z(nab, nh);
                for (size_t qb = 0; qb < bra.prims.size(); ++qb) {
                    const HPrim& pb = bra.prims[qb];
                    for (const HPrim& pk : ket.prims) {
                        const double p = pb.p, q = pk.p;
                        const double rho = p * q / (p + q);
                        hermiteR(Lmax, rho, pb.Px - pk.Px, pb.Py - pk.Py, pb.Pz - pk.Pz,
                                 pref0 / (p * q * std::sqrt(p + q)), R0);
                        for (int h = 0; h < nh; ++h) {
                            const auto& th = Sh.tri[h];
                            const int Lhh = th[0] + th[1] + th[2];
                            for (int k = 0; k < nk; ++k) {
                                const auto& tk = Sk.tri[k];
                                const int Lkk = tk[0] + tk[1] + tk[2];
                                if (Lhh + Lkk > Lmax) {
                                    Rmat(h, k) = 0.0;
                                    continue;
                                }
                                const double r = R0[((size_t)(th[0] + tk[0]) * DR + (th[1] + tk[1])) * DR + (th[2] + tk[2])];
                                Rmat(h, k) = (Lkk & 1) ? -r : r;
                            }
                        }
                        if (mode == 1) {
                            // undifferentiated (ab|cd): bra variant 0 (from the ket-type tables of the bra pair)
                            I0.noalias() += pb.EK[0] * (Rmat * pk.EK[0].transpose());
                            continue;
                        }
                        // M_w = Rmat EK_w^T  (nh x ncd), ket variants: 0 none, 1-3 C_k, 4-9 C_k C_l
                        for (int w = 0; w < 10; ++w) M[w].noalias() = Rmat * pk.EK[w].transpose();
                        const std::vector<Eigen::MatrixXd>& eb = EB[qb];
                        for (int k = 0; k < 3; ++k) {
                            I1[0][k].noalias() += eb[vA(k)] * M[0];
                            I1[1][k].noalias() += eb[vB(k)] * M[0];
                            I1[2][k].noalias() += eb[0] * M[1 + k];
                        }
                        // second derivatives contracted with Dq: sum_{ab,cd} Dq (EB_v M_w) = <EB_v, Dq M_w^T>
                        Z.noalias() = Dq * M[0].transpose();
                        for (int k = 0; k < 3; ++k)
                            for (int l = 0; l < 3; ++l) {
                                if (l >= k) {
                                    const double aa = eb[vAA(k, l)].cwiseProduct(Z).sum();
                                    const double bb = eb[vBB(k, l)].cwiseProduct(Z).sum();
                                    Hc[0][0](k, l) += aa;
                                    Hc[1][1](k, l) += bb;
                                    if (l != k) {
                                        Hc[0][0](l, k) += aa;
                                        Hc[1][1](l, k) += bb;
                                    }
                                }
                                Hc[0][1](k, l) += eb[vAB(k, l)].cwiseProduct(Z).sum();
                            }
                        for (int l = 0; l < 3; ++l) {
                            Z.noalias() = Dq * M[1 + l].transpose();
                            for (int k = 0; k < 3; ++k) {
                                Hc[0][2](k, l) += eb[vA(k)].cwiseProduct(Z).sum();
                                Hc[1][2](k, l) += eb[vB(k)].cwiseProduct(Z).sum();
                            }
                        }
                        for (int s = 0; s < 6; ++s) {
                            Z.noalias() = Dq * M[4 + s].transpose();
                            const double cc = eb[0].cwiseProduct(Z).sum();
                            Hc[2][2](kSym[s][0], kSym[s][1]) += cc;
                            if (kSym[s][0] != kSym[s][1]) Hc[2][2](kSym[s][1], kSym[s][0]) += cc;
                        }
                    }
                }
                if (mode == 1) {
                    digest(A, B, C, Dd, I0, deg, P, gl[0]);
                    continue;
                }
                // Derivative Fock: centre D from translational invariance, then per atom (centres on one atom add up)
                Eigen::MatrixXd ID[3];
                for (int k = 0; k < 3; ++k) ID[k] = -(I1[0][k] + I1[1][k] + I1[2][k]);
                int uatoms[4], nu = 0;
                for (int X = 0; X < 4; ++X) {
                    bool seen = false;
                    for (int u = 0; u < nu; ++u) seen = seen || uatoms[u] == atoms[X];
                    if (!seen) uatoms[nu++] = atoms[X];
                }
                for (int u = 0; u < nu; ++u)
                    for (int k = 0; k < 3; ++k) {
                        Eigen::MatrixXd blk = Eigen::MatrixXd::Zero(nab, ncd);
                        for (int X = 0; X < 3; ++X)
                            if (atoms[X] == uatoms[u]) blk += I1[X][k];
                        if (atoms[3] == uatoms[u]) blk += ID[k];
                        digest(A, B, C, Dd, blk, deg, P, gl[3 * uatoms[u] + k]);
                    }
                // Second derivatives: H[B][A], H[C][A], H[C][B] by symmetry, centre D by translational invariance
                Hc[1][0] = Hc[0][1].transpose();
                Hc[2][0] = Hc[0][2].transpose();
                Hc[2][1] = Hc[1][2].transpose();
                extendTranslational(3, Hc);
                distribute(4, atoms, Hc, 1.0, hl);
            }
        }
    }
    if (mode == 1) {
        G0 = Matrix::Zero(n, n);
        for (int ch = 0; ch < kChunks; ++ch) G0 += gchunk[ch][0];
        G0 = Matrix(G0 + G0.transpose());
        return;
    }
    G1.assign(n3, Matrix::Zero(n, n));
    hess = Matrix::Zero(n3, n3);
    for (int ch = 0; ch < kChunks; ++ch) {
        for (int x = 0; x < n3; ++x) G1[x] += gchunk[ch][x];
        hess += hchunk[ch];
    }
    for (int x = 0; x < n3; ++x) G1[x] = Matrix(G1[x] + G1[x].transpose());
}

}  // namespace

void twoElectronDerivatives(const std::vector<GTO::Orbital>& basis, const Matrix& P, int natoms,
                            std::vector<Matrix>& G1, Matrix& hess, int threads)
{
    Matrix unused;
    twoElectronDriver(basis, P, natoms, threads, 0, G1, hess, unused);
}

Matrix fockTwoElectronReference(const std::vector<GTO::Orbital>& basis, const Matrix& P, int threads)
{
    std::vector<Matrix> unusedG;
    Matrix unusedH, G0;
    twoElectronDriver(basis, P, 0, threads, 1, unusedG, unusedH, G0);
    return G0;
}

Matrix nuclearRepulsionHessian(const std::vector<int>& Z, const Matrix& pos)
{
    const int N = (int)Z.size();
    Matrix H = Matrix::Zero(3 * N, 3 * N);
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j) {
            const Eigen::Vector3d d(pos(i, 0) - pos(j, 0), pos(i, 1) - pos(j, 1), pos(i, 2) - pos(j, 2));
            const double r2 = d.squaredNorm();
            if (r2 < 1.0e-24) continue;
            const double r = std::sqrt(r2);
            // d2(1/r)/dd_k dd_l = (3 d_k d_l - r^2 delta_kl) / r^5, d = R_i - R_j
            const Mat3 T = (3.0 * d * d.transpose() - r2 * Mat3::Identity()) / (r2 * r2 * r);
            const double zz = (double)Z[i] * Z[j];
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l) {
                    H(3 * i + k, 3 * i + l) += zz * T(k, l);
                    H(3 * j + k, 3 * j + l) += zz * T(k, l);
                    H(3 * i + k, 3 * j + l) -= zz * T(k, l);
                    H(3 * j + k, 3 * i + l) -= zz * T(k, l);
                }
        }
    return H;
}

}  // namespace qmhess

// =========================================================================================================
// QMEngine::computeHessian: assembly, CPHF, response terms
// =========================================================================================================
bool QMEngine::computeHessian(const qmhess::HessianOptions& opt, Matrix& hessian_bohr, qmhess::HessianInfo& info)
{
    using clk = std::chrono::steady_clock;
    auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    const auto t0 = clk::now();
    info = qmhess::HessianInfo();
    auto refuse = [&](const std::string& why) {
        info.error = "analytic Hessian: " + why;
        return false;
    };
    // Scope (Nov 2026): closed-shell RHF with exact integrals, all-electron, s-f shells, CPU.
    if (m_functional != QMFunctional::HF) return refuse("only -method hf (closed-shell RHF) is implemented");
    if (m_open_shell) return refuse("open shell (UHF) is not implemented");
    if (hasECP()) return refuse("effective core potentials are not implemented (use an all-electron basis)");
    if (riActive()) return refuse("density fitting (-qm.ri) is not implemented, use -qm.ri off");
    if (!m_scf_ready || !m_scf_converged) return refuse("needs a converged SCF of this geometry");
    const int Lmax = qmhess::maxAngularMomentum(m_gto_basis);
    if (Lmax > 3) return refuse("g and higher shells are not validated (s, p, d, f only)");

    const int natoms = m_atomcount;
    const int n3 = 3 * natoms;
    const int nocc = m_num_electrons / 2;
    const int nbf = m_nbf;
    const int nvir = nbf - nocc;
    if (nocc <= 0 || nvir <= 0) return refuse("no occupied or no virtual orbitals");

    const Matrix& C = m_mo;
    const Vector& eps = m_energies;
    const Matrix& P = m_density;
    const Matrix& F = m_fock;
    Matrix W = Matrix::Zero(nbf, nbf);
    for (int i = 0; i < nocc; ++i) W += 2.0 * eps(i) * (C.col(i) * C.col(i).transpose());
    const bool sph = m_Q.size() != 0;
    const Matrix Pc = sph ? Matrix(m_Q * P * m_Q.transpose()) : P;
    const Matrix Wc = sph ? Matrix(m_Q * W * m_Q.transpose()) : W;
    auto toActive = [&](const Matrix& Mc) { return sph ? Matrix(m_Q.transpose() * Mc * m_Q) : Mc; };
    Matrix posBohr(natoms, 3);
    for (int i = 0; i < natoms; ++i)
        for (int k = 0; k < 3; ++k) posBohr(i, k) = CurcumaUnit::Length::angstrom_to_bohr(m_geometry(i, k));

    // --- skeleton (fixed-density) terms and derivative matrices ---
    std::vector<Matrix> S1c, H1c, G1c;
    Matrix h1e, h2e;
    qmhess::oneElectronDerivatives(m_gto_basis, m_core_charges, posBohr, Pc, Wc, S1c, H1c, h1e, m_threads);
    const auto t1 = clk::now();
    qmhess::twoElectronDerivatives(m_gto_basis, Pc, natoms, G1c, h2e, m_threads);
    const auto t2 = clk::now();
    const Matrix hnn = qmhess::nuclearRepulsionHessian(m_core_charges, posBohr);
    Matrix H = h1e + h2e + hnn;

    std::vector<Matrix> S1(n3), F1(n3);
    for (int x = 0; x < n3; ++x) {
        S1[x] = toActive(S1c[x]);
        F1[x] = toActive(H1c[x] + G1c[x]);
    }
    S1c.clear();
    H1c.clear();
    G1c.clear();

    // --- CPHF: (eps_a - eps_i) U_ai + [Cv^T G[D(U)] Co]_ai = -B_ai, preconditioned CG, all 3N columns in lockstep ---
    const Matrix Co = C.leftCols(nocc), Cv = C.rightCols(nvir);
    Matrix delta(nvir, nocc);
    for (int a = 0; a < nvir; ++a)
        for (int i = 0; i < nocc; ++i) delta(a, i) = eps(nocc + a) - eps(i);
    if (delta.minCoeff() <= 1.0e-8) return refuse("HOMO-LUMO gap is zero or negative, CPHF not applicable");
    // G[D] = J[D] - 1/2 K[D] for a batch of symmetric response densities. With the stored active-basis tensor
    // (the default up to -qm.eri_max_memory_mb) all densities of one CG iteration go through one pass as GEMMs:
    //   J_x = T vec(D_x),  T((m,n),(l,s)) = (mn|ls);   K_x(m, .) = C_m vec(D_x),  C_m(n,(l,s)) = (ml|ns).
    // Otherwise (integral-direct, -gpu) the engine's buildJK runs once per density.
    const bool batched = !useDirectSCF() && !gpuRequested();
    auto Gblock = [&](const std::vector<Matrix>& Ds) {
        const int M = (int)Ds.size();
        std::vector<Matrix> out(M);
        info.jk_builds += M;
        if (M == 0) return out;
        if (!batched) {
            for (int x = 0; x < M; ++x) {
                Matrix J, K;
                buildJK(Ds[x], J, K);
                out[x] = J - 0.5 * K;
            }
            return out;
        }
        const qmint::ERITensor& T = eriActive();
        const Eigen::Index n2 = (Eigen::Index)nbf * nbf;
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> Tm(T.data(), n2, n2);
        Eigen::MatrixXd Dm(n2, M);
        for (int x = 0; x < M; ++x)
            for (int i = 0; i < nbf; ++i)
                for (int j = 0; j < nbf; ++j) Dm(i * nbf + j, x) = Ds[x](i, j);
        const Eigen::MatrixXd Jm = Tm * Dm;
        // K: the row (m,l) of T, read as the n x n matrix R_ml(n,s) = (ml|ns), times DL[l](s,x) = D_x(l,s)
        std::vector<Eigen::MatrixXd> DL(nbf, Eigen::MatrixXd(nbf, M));
        for (int x = 0; x < M; ++x)
            for (int l = 0; l < nbf; ++l)
                for (int s = 0; s < nbf; ++s) DL[l](s, x) = Ds[x](l, s);
        Eigen::MatrixXd Km(n2, M);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(m_threads > 0 ? m_threads : 1)
#endif
        for (int m = 0; m < nbf; ++m) {
            qmhess::noNestedThreads();
            Eigen::MatrixXd acc = Eigen::MatrixXd::Zero(nbf, M);
            for (int l = 0; l < nbf; ++l) {
                Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> R(
                    T.data() + ((size_t)m * nbf + l) * (size_t)n2, nbf, nbf);
                acc.noalias() += R * DL[l];
            }
            Km.middleRows((Eigen::Index)m * nbf, nbf) = acc;
        }
        for (int x = 0; x < M; ++x) {
            out[x].resize(nbf, nbf);
            for (int i = 0; i < nbf; ++i)
                for (int j = 0; j < nbf; ++j) out[x](i, j) = Jm(i * nbf + j, x) - 0.5 * Km(i * nbf + j, x);
        }
        return out;
    };
    auto responseDensity = [&](const Matrix& U) {
        const Matrix X = Cv * U * Co.transpose();
        return Matrix(2.0 * (X + X.transpose()));
    };
    std::vector<Matrix> Dsx(n3);
    for (int x = 0; x < n3; ++x) Dsx[x] = -0.5 * P * S1[x] * P;  // occupied-occupied (metric) part of P^x
    const std::vector<Matrix> GDs = Gblock(Dsx);
    std::vector<Matrix> U(n3), r(n3), z(n3), pdir(n3), Gacc(n3);
    std::vector<double> rz(n3, 0.0), rmax(n3, 0.0);
    std::vector<int> iters(n3, 0);
    std::vector<char> active(n3, 0);
    for (int x = 0; x < n3; ++x) {
        Matrix B = Cv.transpose() * (F1[x] + GDs[x]) * Co;
        const Matrix Svo = Cv.transpose() * S1[x] * Co;
        for (int i = 0; i < nocc; ++i) B.col(i) -= eps(i) * Svo.col(i);
        r[x] = -B;
        U[x] = Matrix::Zero(nvir, nocc);
        Gacc[x] = Matrix::Zero(nbf, nbf);  // G[D(U)] accumulated along the CG steps (linear in U)
        rmax[x] = r[x].cwiseAbs().maxCoeff();
        if (rmax[x] > opt.cphf_threshold) {
            active[x] = 1;
            z[x] = r[x].cwiseQuotient(delta);
            pdir[x] = z[x];
            rz[x] = r[x].cwiseProduct(z[x]).sum();
        }
    }
    for (int it = 1; it <= opt.cphf_max_iterations; ++it) {
        std::vector<int> cols;
        std::vector<Matrix> Dp;
        for (int x = 0; x < n3; ++x)
            if (active[x]) {
                cols.push_back(x);
                Dp.push_back(responseDensity(pdir[x]));
            }
        if (cols.empty()) break;
        const std::vector<Matrix> Gp = Gblock(Dp);
        for (size_t c = 0; c < cols.size(); ++c) {
            const int x = cols[c];
            const Matrix Ap = delta.cwiseProduct(pdir[x]) + Cv.transpose() * Gp[c] * Co;
            const double pAp = pdir[x].cwiseProduct(Ap).sum();
            if (!(pAp > 0.0))
                return refuse(fmt::format("CPHF operator not positive definite (column {}, p^T A p = {:.3e}); "
                                          "the RHF solution is not a stable minimum", x, pAp));
            const double alpha = rz[x] / pAp;
            U[x] += alpha * pdir[x];
            Gacc[x] += alpha * Gp[c];
            r[x] -= alpha * Ap;
            rmax[x] = r[x].cwiseAbs().maxCoeff();
            iters[x] = it;
            if (rmax[x] <= opt.cphf_threshold) {
                active[x] = 0;
                continue;
            }
            z[x] = r[x].cwiseQuotient(delta);
            const double rz_new = r[x].cwiseProduct(z[x]).sum();
            pdir[x] = z[x] + (rz_new / rz[x]) * pdir[x];
            rz[x] = rz_new;
        }
    }
    std::vector<Matrix> P1(n3), Gresp(n3);
    for (int x = 0; x < n3; ++x) {
        if (rmax[x] > opt.cphf_threshold)
            return refuse(fmt::format("CPHF did not converge in {} iterations (column {}, max|r| = {:.3e} > {:.1e})",
                                      opt.cphf_max_iterations, x, rmax[x], opt.cphf_threshold));
        info.cphf_iterations = std::max(info.cphf_iterations, iters[x]);
        info.cphf_residual = std::max(info.cphf_residual, rmax[x]);
        P1[x] = Dsx[x] + responseDensity(U[x]);
        Gresp[x] = GDs[x] + Gacc[x];  // G[P^x]
    }
    const auto t3 = clk::now();

    // --- response terms: H_xy += Tr(P^y F^x) - Tr(W^y S^x), W^y = 1/2 (P^y F P + P (F^y + G[P^y]) P + P F P^y) ---
    std::vector<Matrix> W1(n3);
    for (int y = 0; y < n3; ++y) {
        const Matrix PyFP = P1[y] * F * P;
        W1[y] = 0.5 * (PyFP + PyFP.transpose() + P * (F1[y] + Gresp[y]) * P);
    }
    for (int x = 0; x < n3; ++x)
        for (int y = 0; y < n3; ++y) H(x, y) += P1[y].cwiseProduct(F1[x]).sum() - W1[y].cwiseProduct(S1[x]).sum();

    double asym = 0.0;
    for (int x = 0; x < n3; ++x)
        for (int y = 0; y < x; ++y) asym = std::max(asym, std::abs(H(x, y) - H(y, x)));
    info.max_asymmetry = asym;
    hessian_bohr = 0.5 * (H + H.transpose());
    const auto t4 = clk::now();
    info.t_1e_ms = ms(t0, t1);
    info.t_2e_ms = ms(t1, t2);
    info.t_cphf_ms = ms(t2, t3);
    info.t_total_ms = ms(t0, t4);
    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "QM: analytic RHF Hessian ({} atoms, {} functions): 1e {:.1f} ms, 2e derivative integrals {:.1f} ms, CPHF "
            "{:.1f} ms ({} J/K builds, <= {} CG iterations, max|r| {:.2e}), total {:.1f} ms, asymmetry {:.2e} Eh/Bohr^2, "
            "{} threads",
            natoms, nbf, info.t_1e_ms, info.t_2e_ms, info.t_cphf_ms, info.jk_builds, info.cphf_iterations,
            info.cphf_residual, info.t_total_ms, asym, m_threads));
    return true;
}
