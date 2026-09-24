/*
 * <Native KS-DFT 1-Electron GTO Integrals -- implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (WP1): Obara-Saika overlap, gradient-form kinetic,
 * McMurchie-Davidson nuclear attraction, and the cartesian->spherical d
 * transform. See qm_integrals.hpp for conventions and literature.
 *
 * This program is free software under GPL-3.0
 */

#include "qm_integrals.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

namespace qmint {

// Full-precision pi (the global.h `pi` is only 12 digits).
static const double PI = 3.14159265358979323846264338327950288;

// ---------------------------------------------------------------------------
// Small integer helpers
// ---------------------------------------------------------------------------

static inline double dfact(int n)
{
    // (2n-1)!! with (-1)!! = 1, (-2)!! = 1 (used for l=0 normalization).
    if (n <= 0) return 1.0;
    double r = 1.0;
    for (int k = 1; k <= n; ++k) r *= (2.0 * k - 1);
    return r;
}

double doubleFactorial(int n) { return dfact(n); }

double primitiveNorm(double alpha, int l, int m, int n)
{
    // N(alpha,L) = (2*alpha/pi)^(3/4) * (4*alpha)^(L/2) / sqrt((2l-1)!!(2m-1)!!(2n-1)!!)
    const int L = l + m + n;
    double denom = std::sqrt(dfact(l) * dfact(m) * dfact(n));
    if (denom == 0.0) denom = 1.0;
    return std::pow(2.0 * alpha / PI, 0.75) * std::pow(4.0 * alpha, L / 2.0) / denom;
}

// ---------------------------------------------------------------------------
// OS overlap 1D auxiliary: S(i,j) along one cartesian axis.
//   S(0,0) = 1
//   S(i+1,j) = PA * S(i,j) + (1/(2g)) * (i * S(i-1,j) + j * S(i,j-1))
//   S(i,j+1) = PB * S(i,j) + (1/(2g)) * (i * S(i-1,j) + j * S(i,j-1))
//   S with a negative index = 0.
// Returns a table S[i][j] for i in 0..imax, j in 0..jmax (inclusive).
// ---------------------------------------------------------------------------

static std::vector<std::vector<double>> overlap1DTable(int imax, int jmax,
                                                       double PA, double PB, double gamma)
{
    const double g2 = 0.5 / gamma;
    std::vector<std::vector<double>> S(imax + 1, std::vector<double>(jmax + 1, 0.0));
    S[0][0] = 1.0;
    // Raise i (A side), j held at 0.
    for (int i = 1; i <= imax; ++i) {
        double prev = (i >= 2) ? S[i - 2][0] : 0.0;
        S[i][0] = PA * S[i - 1][0] + g2 * ((i - 1) * prev);
    }
    // Raise j (B side).
    for (int i = 0; i <= imax; ++i) {
        for (int j = 1; j <= jmax; ++j) {
            double si1 = (i >= 1) ? S[i - 1][j - 1] : 0.0;  // i*S(i-1,j-1)
            double sj1 = (j >= 2) ? S[i][j - 2] : 0.0;      // (j-1)*S(i,j-2)
            S[i][j] = PB * S[i][j - 1] + g2 * (i * si1 + (j - 1) * sj1);
        }
    }
    return S;
}

// Overlap of two primitives: (PI/g)^(3/2) * K_AB * Sx[l1,l2] * Sy[m1,m2] * Sz[n1,n2].
// Centers A,B and exponents alpha1,alpha2 in atomic units.
static inline double gaussianProductK(double alpha1, double alpha2,
                                       double Ax, double Ay, double Az,
                                       double Bx, double By, double Bz,
                                       double& gamma, double& Px, double& Py, double& Pz)
{
    gamma = alpha1 + alpha2;
    Px = (alpha1 * Ax + alpha2 * Bx) / gamma;
    Py = (alpha1 * Ay + alpha2 * By) / gamma;
    Pz = (alpha1 * Az + alpha2 * Bz) / gamma;
    const double zeta = alpha1 * alpha2 / gamma;
    const double R2 = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);
    return std::exp(-zeta * R2);  // K_AB (the (PI/g)^{3/2} prefactor is added by the caller)
}

// ---------------------------------------------------------------------------
// Self-overlap of a contracted orbital (for renormalization): uses the same
// OS overlap primitive. Both centers coincide (PA = PB = 0) but the recurrence
// still works; we just call the general primitive overlap with A == B.
// ---------------------------------------------------------------------------

static double primitiveOverlap(int l1, int m1, int n1, int l2, int m2, int n2,
                               double alpha1, double alpha2,
                               double Ax, double Ay, double Az,
                               double Bx, double By, double Bz)
{
    double gamma, Px, Py, Pz;
    double K = gaussianProductK(alpha1, alpha2, Ax, Ay, Az, Bx, By, Bz, gamma, Px, Py, Pz);
    if (K == 0.0) return 0.0;
    const double pref = std::pow(PI / gamma, 1.5) * K;
    auto Sx = overlap1DTable(std::max(l1, l2), std::max(l1, l2), Px - Ax, Px - Bx, gamma);
    auto Sy = overlap1DTable(std::max(m1, m2), std::max(m1, m2), Py - Ay, Py - By, gamma);
    auto Sz = overlap1DTable(std::max(n1, n2), std::max(n1, n2), Pz - Az, Pz - Bz, gamma);
    return pref * Sx[l1][l2] * Sy[m1][m2] * Sz[n1][n2];
}

void normalizeOrbitalSelfOverlap(GTO::Orbital& orb)
{
    if (orb.exponents.empty()) return;
    int l, m, n;
    GTO::orbitalTypeToComponents(orb.type, l, m, n);
    const size_t K = orb.exponents.size();
    double sii = 0.0;
    for (size_t a = 0; a < K; ++a) {
        for (size_t b = 0; b < K; ++b) {
            sii += orb.coefficients[a] * orb.coefficients[b] *
                   primitiveOverlap(l, m, n, l, m, n,
                                    orb.exponents[a], orb.exponents[b],
                                    orb.x, orb.y, orb.z, orb.x, orb.y, orb.z);
        }
    }
    if (sii > 1e-15) {
        const double inv = 1.0 / std::sqrt(sii);
        for (size_t a = 0; a < K; ++a) orb.coefficients[a] *= inv;
    }
}

// ---------------------------------------------------------------------------
// Contracted overlap matrix S
// ---------------------------------------------------------------------------

static double contractedOverlap(const GTO::Orbital& a, const GTO::Orbital& b)
{
    int l1, m1, n1, l2, m2, n2;
    GTO::orbitalTypeToComponents(a.type, l1, m1, n1);
    GTO::orbitalTypeToComponents(b.type, l2, m2, n2);
    double sum = 0.0;
    for (size_t ia = 0; ia < a.exponents.size(); ++ia) {
        for (size_t ib = 0; ib < b.exponents.size(); ++ib) {
            double gamma, Px, Py, Pz;
            double K = gaussianProductK(a.exponents[ia], b.exponents[ib],
                                        a.x, a.y, a.z, b.x, b.y, b.z, gamma, Px, Py, Pz);
            if (K == 0.0) continue;
            const double pref = std::pow(PI / gamma, 1.5) * K;
            auto Sx = overlap1DTable(std::max(l1, l2), std::max(l1, l2), Px - a.x, Px - b.x, gamma);
            auto Sy = overlap1DTable(std::max(m1, m2), std::max(m1, m2), Py - a.y, Py - b.y, gamma);
            auto Sz = overlap1DTable(std::max(n1, n2), std::max(n1, n2), Pz - a.z, Pz - b.z, gamma);
            sum += a.coefficients[ia] * b.coefficients[ib] * pref *
                   Sx[l1][l2] * Sy[m1][m2] * Sz[n1][n2];
        }
    }
    return sum;
}

Matrix buildOverlap(const std::vector<GTO::Orbital>& basis)
{
    const int n = (int)basis.size();
    Matrix S = Matrix::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double v = contractedOverlap(basis[i], basis[j]);
            S(i, j) = v;
            S(j, i) = v;
        }
    }
    return S;
}

// ---------------------------------------------------------------------------
// Contracted kinetic matrix T  via  T_ab = (1/2) <grad g_a | grad g_b>
//
//   <d/dx g_a | d/dx g_b> = l1*l2 <g_a^{l1-1}|g_b^{l2-1}>
//                        - 2*beta*l1 <g_a^{l1-1}|g_b^{l2+1}>
//                        - 2*alpha*l2 <g_a^{l1+1}|g_b^{l2-1}>
//                        + 4*alpha*beta <g_a^{l1+1}|g_b^{l2+1}>
//   where <g_a^p|g_b^q> is the overlap primitive with the x-angular-momentum
//   shifted to (p,q) (y,z unchanged). alpha = a exponent, beta = b exponent.
// ---------------------------------------------------------------------------

// General overlap primitive with arbitrary angular momenta (la,ma,na) on center A
// and (lb,mb,nb) on center B. Used by the kinetic kernel (which shifts each axis
// independently via the gradient identity) -- geometric displacements are tied
// to the correct axis, unlike a single "shifted x" helper.
static double primitiveOverlapGeneral(int la, int ma, int na,
                                              int lb, int mb, int nb,
                                              double gamma, double K,
                                              double PAx, double PBx,
                                              double PAy, double PBy,
                                              double PAz, double PBz)
{
    if (K == 0.0) return 0.0;
    const double pref = std::pow(PI / gamma, 1.5) * K;
    auto Sx = overlap1DTable(std::max(la, lb), std::max(la, lb), PAx, PBx, gamma);
    auto Sy = overlap1DTable(std::max(ma, mb), std::max(ma, mb), PAy, PBy, gamma);
    auto Sz = overlap1DTable(std::max(na, nb), std::max(na, nb), PAz, PBz, gamma);
    return pref * Sx[la][lb] * Sy[ma][mb] * Sz[na][nb];
}

// Primitive kinetic energy <g_a| -1/2 nabla^2 |g_b> for arbitrary cartesian powers,
// via the gradient identity documented above. Split out of contractedKineticPair
// (Sep 2026) so the analytic gradient can evaluate it at shifted powers.
static double primitiveKinetic(int l1, int m1, int n1, int l2, int m2, int n2,
                               double alpha, double beta,
                               double Ax, double Ay, double Az,
                               double Bx, double By, double Bz)
{
    double gamma, Px, Py, Pz;
    double K = gaussianProductK(alpha, beta, Ax, Ay, Az, Bx, By, Bz, gamma, Px, Py, Pz);
    if (K == 0.0) return 0.0;
    const double PAx = Px - Ax, PBx = Px - Bx;
    const double PAy = Py - Ay, PBy = Py - By;
    const double PAz = Pz - Az, PBz = Pz - Bz;

    // x contribution: <d/dx g_a | d/dx g_b>
    double tx = 0.0;
    if (l1 > 0 && l2 > 0) tx += l1 * l2 * primitiveOverlapGeneral(l1 - 1, m1, n1, l2 - 1, m2, n2, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
    if (l1 > 0)           tx -= 2.0 * beta * l1 * primitiveOverlapGeneral(l1 - 1, m1, n1, l2 + 1, m2, n2, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
    if (l2 > 0)           tx -= 2.0 * alpha * l2 * primitiveOverlapGeneral(l1 + 1, m1, n1, l2 - 1, m2, n2, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
                           tx += 4.0 * alpha * beta * primitiveOverlapGeneral(l1 + 1, m1, n1, l2 + 1, m2, n2, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
    // y contribution: <d/dy g_a | d/dy g_b>
    double ty = 0.0;
    if (m1 > 0 && m2 > 0) ty += m1 * m2 * primitiveOverlapGeneral(l1, m1 - 1, n1, l2, m2 - 1, n2, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
    if (m1 > 0)           ty -= 2.0 * beta * m1 * primitiveOverlapGeneral(l1, m1 - 1, n1, l2, m2 + 1, n2, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
    if (m2 > 0)           ty -= 2.0 * alpha * m2 * primitiveOverlapGeneral(l1, m1 + 1, n1, l2, m2 - 1, n2, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
                           ty += 4.0 * alpha * beta * primitiveOverlapGeneral(l1, m1 + 1, n1, l2, m2 + 1, n2, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
    // z contribution: <d/dz g_a | d/dz g_b>
    double tz = 0.0;
    if (n1 > 0 && n2 > 0) tz += n1 * n2 * primitiveOverlapGeneral(l1, m1, n1 - 1, l2, m2, n2 - 1, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
    if (n1 > 0)           tz -= 2.0 * beta * n1 * primitiveOverlapGeneral(l1, m1, n1 - 1, l2, m2, n2 + 1, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
    if (n2 > 0)           tz -= 2.0 * alpha * n2 * primitiveOverlapGeneral(l1, m1, n1 + 1, l2, m2, n2 - 1, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);
                           tz += 4.0 * alpha * beta * primitiveOverlapGeneral(l1, m1, n1 + 1, l2, m2, n2 + 1, gamma, K, PAx, PBx, PAy, PBy, PAz, PBz);

    return 0.5 * (tx + ty + tz);
}

static double contractedKineticPair(const GTO::Orbital& a, const GTO::Orbital& b)
{
    int l1, m1, n1, l2, m2, n2;
    GTO::orbitalTypeToComponents(a.type, l1, m1, n1);
    GTO::orbitalTypeToComponents(b.type, l2, m2, n2);
    double sum = 0.0;
    for (size_t ia = 0; ia < a.exponents.size(); ++ia)
        for (size_t ib = 0; ib < b.exponents.size(); ++ib)
            sum += a.coefficients[ia] * b.coefficients[ib] *
                   primitiveKinetic(l1, m1, n1, l2, m2, n2, a.exponents[ia], b.exponents[ib],
                                    a.x, a.y, a.z, b.x, b.y, b.z);
    return sum;
}

Matrix buildKinetic(const std::vector<GTO::Orbital>& basis)
{
    const int n = (int)basis.size();
    Matrix T = Matrix::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double v = contractedKineticPair(basis[i], basis[j]);
            T(i, j) = v;
            T(j, i) = v;
        }
    }
    return T;
}

// ---------------------------------------------------------------------------
// McMurchie-Davidson nuclear attraction.
//   Hermite expansion coeffs E_t^x(i,j) (Helgaker 9.9.2-9.9.8), two-pass:
//     E_0(0,0) = K_AB
//     E_0(i+1,0) = PA_x*E_0(i,0) + 1/(2g)*(i)*E_0(i-1,0)
//     E_0(i,j+1) = PB_x*E_0(i,j) + 1/(2g)*(i*E_0(i-1,j) + j*E_0(i,j-1))
//     E_{t+1}(i,j) = PA_x*E_t(i,j) + 1/(2g)*(i*E_t(i-1,j) + j*E_t(i,j-1) + t*E_{t-1}(i,j))
//   (Negative indices contribute 0.) E_0(0,0) carries K_AB, so the prefactor
//   (2*pi/g) sits only in R_000^(N). Identical recurrences hold for y, z.
//
//   Boys F_N(T): downward recurrence from a zero asymptote (stable for all T),
//   F_N(0) = 1/(2N+1).
//
//   Coulomb auxiliary (Helgaker 9.9.10-9.9.14):
//     R_000^(N) = (2*pi/g) * F_N(T),  T = g*|P-C|^2
//     R_{p,u,v}^(N) = (p-1)*R_{p-2,u,v}^(N+1) - PC_x*R_{p-1,u,v}^(N+1)
//     R_{t,q,v}^(N) = (q-1)*R_{t,q-2,v}^(N+1) - PC_y*R_{t,q-1,v}^(N+1)
//     R_{t,u,r}^(N) = (r-1)*R_{t,u,r-2}^(N+1) - PC_z*R_{t,u,r-1}^(N+1)
//
//   Primitive: <g_a|1/r_C|g_b> = sum_{t,u,v} E_t^x(la,lb) E_u^y(ma,mb) E_v^z(na,nb)
//                               * R_{t,u,v}^(0).
//   V_ij = - sum_C Z_C * <g_i|1/r_C|g_j>  (Hamiltonian sign).
// ---------------------------------------------------------------------------

// Hermite coefficients E[t][i][j] along one axis. exponents alpha,beta, center
// displacements PA = P-A, PB = P-B, gamma. Returns a (tmax+1) x (iA+1) x (iB+1)
// table; E[0][0][0] = K (the Gaussian product K_AB).
static std::vector<std::vector<std::vector<double>>> hermiteCoeffs(int iA, int iB,
                                                                    double PA, double PB,
                                                                    double gamma, double K)
{
    const int tmax = iA + iB;
    const double g2 = 0.5 / gamma;
    std::vector<std::vector<std::vector<double>>> E(
        tmax + 1, std::vector<std::vector<double>>(iA + 1, std::vector<double>(iB + 1, 0.0)));
    if (K == 0.0) return E;
    E[0][0][0] = K;
    // Standard McMurchie-Davidson forward recursion (Helgaker 9.5.5/9.5.6). It
    // raises the polynomial powers i and j and fills every Hermite order t at
    // once:
    //   E^{i+1,j}_t = (1/(2p)) E^{i,j}_{t-1} + PA E^{i,j}_t + (t+1) E^{i,j}_{t+1}
    //   E^{i,j+1}_t = (1/(2p)) E^{i,j}_{t-1} + PB E^{i,j}_t + (t+1) E^{i,j}_{t+1}
    // Only rows (i-1,j) and (i,j-1) are read, both already complete, so no entry
    // is read before it is written.
    auto at = [&](int t, int i, int j) -> double {
        if (t < 0 || t > tmax || i < 0 || j < 0 || i > iA || j > iB) return 0.0;
        return E[t][i][j];
    };
    for (int i = 1; i <= iA; ++i)
        for (int t = 0; t <= tmax; ++t)
            E[t][i][0] = g2 * at(t - 1, i - 1, 0) + PA * at(t, i - 1, 0)
                       + (t + 1) * at(t + 1, i - 1, 0);
    for (int j = 1; j <= iB; ++j)
        for (int i = 0; i <= iA; ++i)
            for (int t = 0; t <= tmax; ++t)
                E[t][i][j] = g2 * at(t - 1, i, j - 1) + PB * at(t, i, j - 1)
                           + (t + 1) * at(t + 1, i, j - 1);
    return E;
}

// Boys function array F_0..F_maxN via downward recurrence from a zero asymptote.
// Fills F[0..maxN] (resized, capacity reused) -- the blocked kernels call this once
// per primitive quartet, so it must not allocate (Claude Generated, Sep 2026).
static void boysArrayInto(int maxN, double T, std::vector<double>& F)
{
    F.assign(maxN + 1, 0.0);
    if (T < 1e-14) {
        for (int n = 0; n <= maxN; ++n) F[n] = 1.0 / (2.0 * n + 1.0);
        return;
    }
    // Large T: the downward recurrence below needs a starting index far above T,
    // and a fixed maxN+25 start is NOT enough -- F_0(37) came out 50x too small and
    // F_0(T >= 50) as exactly 0, which silently wrecked every integral whose
    // Gaussian pair sits away from the nucleus (i.e. every real molecule; the
    // single-atom T == 0 cases were the only ones this survived). For T >= 1 use
    // the closed form F_0(T) = 0.5 sqrt(pi/T) erf(sqrt(T)) and recur UPWARD,
    // F_{n+1} = ((2n+1) F_n - e^{-T}) / (2T), which is the stable direction there
    // (the subtraction cancels only at small T, hence the split at T = 1).
    // Worst relative error over T in [1e-14, 800] and n <= 5: 4.7e-14.
    if (T >= 1.0) {
        const long double Tl = (long double)T;
        const long double eT = expl(-Tl);
        F[0] = (double)(0.5L * sqrtl(PI / Tl) * erfl(sqrtl(Tl)));
        for (int n = 0; n < maxN; ++n)
            F[n + 1] = (double)(((2.0L * n + 1.0L) * (long double)F[n] - eT) / (2.0L * Tl));
        return;
    }
    const int M = maxN + 25;
    thread_local std::vector<long double> G;
    G.assign(M + 1, 0.0L);
    const long double eT = expl(-(long double)T);
    for (int n = M - 1; n >= 0; --n)
        G[n] = (2.0L * T * G[n + 1] + eT) / (2.0L * n + 1.0L);
    for (int n = 0; n <= maxN; ++n) F[n] = (double)G[n];
    return;
}

static std::vector<double> boysArray(int maxN, double T)
{
    std::vector<double> F;
    boysArrayInto(maxN, T, F);
    return F;
}

// R_{t,u,v}^(N) block: t in [0,tmax], u in [0,umax], v in [0,vmax], N in [0,Nmax].
// Nmax must be >= tmax+umax+vmax for the recursion to reach N=0 at the top.
static void buildRblock(int tmax, int umax, int vmax, int Nmax,
                        double PCx, double PCy, double PCz, double gamma,
                        const std::vector<double>& boys,
                        std::vector<double>& R)  // size (tmax+1)(umax+1)(vmax+1)(Nmax+1)
{
    auto idx = [&](int t, int u, int v, int N) {
        return ((t * (umax + 1) + u) * (vmax + 1) + v) * (Nmax + 1) + N;
    };
    std::fill(R.begin(), R.end(), 0.0);
    const double pref = 2.0 * PI / gamma;
    // base R_000^(N) = (2 pi / gamma) (-2 gamma)^N F_N(T). The (-2 gamma)^N factor is
    // required by the t >= 2 recurrence below (Helgaker 9.9.13: R^n_{t+1} =
    // X_PC R^n_t + t R^{n+1}_{t-1}, whose R^{n+1} carries one more (-2 gamma)).
    // Without it R^0_{200} came out +(2 pi/gamma) F_1 instead of -(2 pi/gamma) 2 gamma
    // F_1 -- wrong magnitude AND sign -- so every nuclear-attraction element with
    // total angular momentum >= 2 was wrong (V_pp exactly 2x too large on He) and
    // the HF SCF had no fixed point. Verified against the exact pi/6 for <px|1/r|px>
    // and against ORCA's He core Hamiltonian (h_pp = +0.372308 Eh).
    //
    // The displacement term carries a PLUS sign: this is the bra recurrence
    // (differentiation w.r.t. P), identical in form to buildRblockERI's. The former
    // minus reproduced every on-centre integral (P-C == 0 kills the term) but gave
    // off-centre elements the right magnitude with an inverted sign -- e.g.
    // <pz_A|1/r_B|s_A> came out -0.238696 instead of +0.238695.
    {
        double neg2p = 1.0;  // (-2 gamma)^N
        for (int N = 0; N <= Nmax; ++N) {
            R[idx(0, 0, 0, N)] = pref * neg2p * boys[N];
            neg2p *= -2.0 * gamma;
        }
    }
    for (int t = 0; t <= tmax; ++t) {
        // R[t][0][0]: base if t==0 else raise t
        if (t >= 1) {
            for (int N = 0; N <= Nmax - 1; ++N) {
                double low = (t >= 2) ? (t - 1) * R[idx(t - 2, 0, 0, N + 1)] : 0.0;
                R[idx(t, 0, 0, N)] = low + PCx * R[idx(t - 1, 0, 0, N + 1)];
            }
        }
        // raise u (v = 0)
        for (int u = 1; u <= umax; ++u) {
            for (int N = 0; N <= Nmax - 1; ++N) {
                double low = (u >= 2) ? (u - 1) * R[idx(t, u - 2, 0, N + 1)] : 0.0;
                R[idx(t, u, 0, N)] = low + PCy * R[idx(t, u - 1, 0, N + 1)];
            }
        }
        // raise v for all u
        for (int u = 0; u <= umax; ++u) {
            for (int v = 1; v <= vmax; ++v) {
                for (int N = 0; N <= Nmax - 1; ++N) {
                    double low = (v >= 2) ? (v - 1) * R[idx(t, u, v - 2, N + 1)] : 0.0;
                    R[idx(t, u, v, N)] = low + PCz * R[idx(t, u, v - 1, N + 1)];
                }
            }
        }
    }
}

static double primitiveNuclearAttraction(int la, int ma, int na, int lb, int mb, int nb,
                                         double alpha, double beta,
                                         double Ax, double Ay, double Az,
                                         double Bx, double By, double Bz,
                                         double Cx, double Cy, double Cz)
{
    double gamma, Px, Py, Pz;
    double K = gaussianProductK(alpha, beta, Ax, Ay, Az, Bx, By, Bz, gamma, Px, Py, Pz);
    if (K == 0.0) return 0.0;

    auto Ex = hermiteCoeffs(la, lb, Px - Ax, Px - Bx, gamma, K);
    auto Eu = hermiteCoeffs(ma, mb, Py - Ay, Py - By, gamma, 1.0);
    auto Ev = hermiteCoeffs(na, nb, Pz - Az, Pz - Bz, gamma, 1.0);
    // NOTE: K is folded into Ex only; Eu/Ev are started from 1.0 so their tables
    // carry the pure Hermite coefficients without re-introducing K.

    const int tmax = la + lb, umax = ma + mb, vmax = na + nb;
    const int Nmax = tmax + umax + vmax;
    const double T = gamma * ((Px - Cx) * (Px - Cx) + (Py - Cy) * (Py - Cy) + (Pz - Cz) * (Pz - Cz));
    auto F = boysArray(Nmax, T);

    std::vector<double> R((tmax + 1) * (umax + 1) * (vmax + 1) * (Nmax + 1), 0.0);
    buildRblock(tmax, umax, vmax, Nmax, Px - Cx, Py - Cy, Pz - Cz, gamma, F, R);

    const int strideN = (umax + 1) * (vmax + 1) * (Nmax + 1);
    const int strideU = (vmax + 1) * (Nmax + 1);
    const int strideV = (Nmax + 1);
    double val = 0.0;
    for (int t = 0; t <= tmax; ++t) {
        const double et = Ex[t][la][lb];
        if (et == 0.0) continue;
        for (int u = 0; u <= umax; ++u) {
            const double eu = Eu[u][ma][mb];
            if (eu == 0.0) continue;
            for (int v = 0; v <= vmax; ++v) {
                const double ev = Ev[v][na][nb];
                if (ev == 0.0) continue;
                const double r = R[t * strideN + u * strideU + v * strideV + 0];
                val += et * eu * ev * r;
            }
        }
    }
    return val;
}

static double contractedNuclearPair(const GTO::Orbital& a, const GTO::Orbital& b,
                                    double Cx, double Cy, double Cz)
{
    int l1, m1, n1, l2, m2, n2;
    GTO::orbitalTypeToComponents(a.type, l1, m1, n1);
    GTO::orbitalTypeToComponents(b.type, l2, m2, n2);
    double sum = 0.0;
    for (size_t ia = 0; ia < a.exponents.size(); ++ia) {
        for (size_t ib = 0; ib < b.exponents.size(); ++ib) {
            sum += a.coefficients[ia] * b.coefficients[ib] *
                   primitiveNuclearAttraction(l1, m1, n1, l2, m2, n2,
                                             a.exponents[ia], b.exponents[ib],
                                             a.x, a.y, a.z, b.x, b.y, b.z,
                                             Cx, Cy, Cz);
        }
    }
    return sum;
}

Matrix buildNuclearAttraction(const std::vector<GTO::Orbital>& basis,
                              const std::vector<int>& atomZ,
                              const Matrix& atomPosBohr)
{
    const int n = (int)basis.size();
    Matrix V = Matrix::Zero(n, n);
    const int natoms = (int)atomZ.size();
    for (int C = 0; C < natoms; ++C) {
        const double Z = atomZ[C];
        if (Z == 0.0) continue;
        const double Cx = atomPosBohr(C, 0);
        const double Cy = atomPosBohr(C, 1);
        const double Cz = atomPosBohr(C, 2);
        for (int i = 0; i < n; ++i) {
            for (int j = i; j < n; ++j) {
                double v = -Z * contractedNuclearPair(basis[i], basis[j], Cx, Cy, Cz);
                V(i, j) += v;
                if (i != j) V(j, i) += v;
            }
        }
    }
    return V;
}

Matrix buildCoreHamiltonian(const std::vector<GTO::Orbital>& basis,
                            const std::vector<int>& atomZ,
                            const Matrix& atomPosBohr)
{
    Matrix T = buildKinetic(basis);
    Matrix V = buildNuclearAttraction(basis, atomZ, atomPosBohr);
    return T + V;
}

// ---------------------------------------------------------------------------
// Cartesian 6d -> spherical 5d transform.
//
// The 5 real spherical harmonics (in the convention used by most programs,
// matching ORCA's def2-SVP pure-d ordering d_{-2..+2}) expressed as vectors over
// the cartesian d block ordered [DXX, DYY, DZZ, DXY, DXZ, DYZ]:
//   d(z2)      ~ ( 2*DZZ - DXX - DYY )           -> [-1, -1, 2, 0, 0, 0]
//   d(xz)      ~ DXZ                            -> [ 0,  0, 0, 0, 1, 0]
//   d(yz)      ~ DYZ                            -> [ 0,  0, 0, 0, 0, 1]
//   d(x2-y2)   ~ ( DXX - DYY )                   -> [ 1, -1, 0, 0, 0, 0]
//   d(xy)      ~ DXY                             -> [ 0,  0, 0, 1, 0, 0]
// Each raw vector is normalized against the actual cartesian d-block overlap
// (v^T S_d v), so the resulting spherical set is orthonormal regardless of the
// primitive/shell normalization convention. s and p blocks pass through
// unchanged (identity). M_sph = Q^T M_cart Q.
// ---------------------------------------------------------------------------

static bool isDShell(const GTO::Orbital& o)
{
    return o.type == GTO::DXX || o.type == GTO::DYY || o.type == GTO::DZZ ||
           o.type == GTO::DXY || o.type == GTO::DXZ || o.type == GTO::DYZ;
}

Matrix buildSphericalTransform(const std::vector<GTO::Orbital>& basis,
                               const Matrix& Scart)
{
    const int n = (int)basis.size();
    // Locate d shells: contiguous runs of 6 cartesian d on the same atom.
    // Order within a shell (as emitted by the parser): DXX, DYY, DZZ, DXY, DXZ, DYZ.
    // Raw spherical d vectors in that cartesian order.
    static const double raw[5][6] = {
        { -1, -1,  2,  0,  0,  0 },  // z2
        {  0,  0,  0,  0,  1,  0 },  // xz
        {  0,  0,  0,  0,  0,  1 },  // yz
        {  1, -1,  0,  0,  0,  0 },  // x2-y2
        {  0,  0,  0,  1,  0,  0 },  // xy
    };

    bool anyD = false;
    for (const auto& o : basis) if (isDShell(o)) { anyD = true; break; }
    if (!anyD) return Matrix();  // empty -> caller keeps cartesian

    // nbf_sph = nbf_cart - (#d shells)  (each 6d shell collapses to 5d)
    int nDshells = 0;
    for (int i = 0; i < n; ) {
        if (isDShell(basis[i])) {
            ++nDshells;
            i += 6;  // assume contiguous block of 6
        } else {
            ++i;
        }
    }
    const int nsph = n - nDshells;
    Matrix Q = Matrix::Zero(n, nsph);

    int sphCol = 0;
    int i = 0;
    while (i < n) {
        if (!isDShell(basis[i])) {
            Q(i, sphCol) = 1.0;  // 1:1 s/p passthrough
            ++i;
            ++sphCol;
            continue;
        }
        // d shell: cartesian rows i..i+5, build 5 orthonormal spherical columns.
        const int base = i;
        // S_d block (6x6) from the cartesian overlap.
        Matrix Sd(6, 6);
        for (int a = 0; a < 6; ++a)
            for (int b = 0; b < 6; ++b)
                Sd(a, b) = Scart(base + a, base + b);
        // Orthonormalize the 5 raw vectors against S_d (modified Gram-Schmidt
        // in the S_d metric), producing 5 columns Q[base..base+5][sphCol..sphCol+4].
        for (int k = 0; k < 5; ++k) {
            Eigen::Matrix<double, 6, 1> v;
            for (int a = 0; a < 6; ++a) v(a) = raw[k][a];
            // remove projections onto already-built spherical columns
            for (int p = 0; p < k; ++p) {
                Eigen::Matrix<double, 6, 1> uprev;
                for (int a = 0; a < 6; ++a) uprev(a) = Q(base + a, sphCol + p);
                double sprod = uprev.transpose() * Sd * v;
                double norm2 = uprev.transpose() * Sd * uprev;
                if (std::abs(norm2) > 1e-15) v -= (sprod / norm2) * uprev;
            }
            double norm2 = v.transpose() * Sd * v;
            if (norm2 < 1e-18) {
                // degenerate -- shouldn't happen for a proper d block
                for (int a = 0; a < 6; ++a) Q(base + a, sphCol + k) = 0.0;
            } else {
                double inv = 1.0 / std::sqrt(norm2);
                for (int a = 0; a < 6; ++a) Q(base + a, sphCol + k) = v(a) * inv;
            }
        }
        i += 6;
        sphCol += 5;
    }
    return Q;
}

Matrix applySphericalTransform(const Matrix& Mcart, const Matrix& Q)
{
    if (Q.size() == 0) return Mcart;
    return Q.transpose() * Mcart * Q;
}

// ===========================================================================
// WP2 -- 4-centre ERI (McMurchie-Davidson)
//
//   (ab|cd) = (2*pi^(5/2))/(p q sqrt(p+q))
//             * sum_{t,u,v} sum_{tau,ups,om}
//               E^ab_{t,u,v} E^cd_{tau,ups,om} (-1)^(tau+ups+om)
//               * R^0_{t+tau, u+ups, v+om}
//
//   p = a+b (bra composite), q = c+d (ket composite), rho = pq/(p+q),
//   Boys arg T = rho*|P-Q|^2. R auxiliary (Helgaker 9.9.13/9.9.17):
//     R^n_{0,0,0} = (-2*rho)^n * F_n(T)
//     R^n_{t+1,u,v} = t*R^{n+1}_{t-1,u,v} + X_{PQ}*R^{n+1}_{t,u,v}      (x)
//     (and u/v analogously), X_{PQ} = P - Q.  NOTE the PLUS sign and the
//   (P-Q) displacement: this is the BRA recurrence (differentiation w.r.t. P).
//   The ket Hermite orders enter with the (-1)^(tau+ups+om) sign because the
//   ket raise is differentiation w.r.t. Q (opposite sign). The recurrence is
//   gamma-independent -- rho appears only in the Boys argument and the base
//   prefactor (-2*rho)^n, never as a recurrence coefficient.
//
//   This is a DIFFERENT auxiliary from the WP1 nuclear buildRblock (which uses
//   base (2*pi/g)*(-2g)^n*F_n and displacement (C-P), i.e. differentiation w.r.t.
//   the source C); the two are not interchangeable, so the ERI gets its own R
//   build. Both carry the (-2*rho)^n / (-2*gamma)^n factor in the base -- the WP1
//   build was missing it until Jul 2026, which made every nuclear-attraction
//   element with total angular momentum >= 2 wrong (see the note in buildRblock).
//   The Hermite expansion coefficients hermiteCoeffs (Helgaker 9.9.2-9.9.8) and
//   the Boys function boysArray are reused unchanged from WP1.
// ===========================================================================

// ERI Coulomb auxiliary R^n_{t,u,v} for the combined (bra+ket) Hermite orders.
// Filled for t in [0,tmax], u in [0,umax], v in [0,vmax], N in [0,Nmax], with
// Nmax >= tmax+umax+vmax. Displacement W = P - Q (bra centre minus ket centre).
// Base R^n_{0,0,0} = (-2*rho)^n * F_n(T); recurrence with +W displacement.
static void buildRblockERI(int tmax, int umax, int vmax, int Nmax,
                           double Wx, double Wy, double Wz, double rho,
                           const std::vector<double>& boys,
                           std::vector<double>& R)
{
    auto idx = [&](int t, int u, int v, int N) {
        return ((t * (umax + 1) + u) * (vmax + 1) + v) * (Nmax + 1) + N;
    };
    std::fill(R.begin(), R.end(), 0.0);
    // base R^n_{0,0,0} = (-2*rho)^n * F_n(T)
    double b = -2.0 * rho;
    double bp = 1.0;
    for (int N = 0; N <= Nmax; ++N) {
        R[idx(0, 0, 0, N)] = bp * boys[N];
        bp *= b;
    }
    for (int t = 0; t <= tmax; ++t) {
        if (t >= 1) {
            for (int N = 0; N <= Nmax - 1; ++N) {
                double low = (t >= 2) ? (t - 1) * R[idx(t - 2, 0, 0, N + 1)] : 0.0;
                R[idx(t, 0, 0, N)] = low + Wx * R[idx(t - 1, 0, 0, N + 1)];
            }
        }
        for (int u = 1; u <= umax; ++u) {
            for (int N = 0; N <= Nmax - 1; ++N) {
                double low = (u >= 2) ? (u - 1) * R[idx(t, u - 2, 0, N + 1)] : 0.0;
                R[idx(t, u, 0, N)] = low + Wy * R[idx(t, u - 1, 0, N + 1)];
            }
        }
        for (int u = 0; u <= umax; ++u) {
            for (int v = 1; v <= vmax; ++v) {
                for (int N = 0; N <= Nmax - 1; ++N) {
                    double low = (v >= 2) ? (v - 1) * R[idx(t, u, v - 2, N + 1)] : 0.0;
                    R[idx(t, u, v, N)] = low + Wz * R[idx(t, u, v - 1, N + 1)];
                }
            }
        }
    }
}

// R^0_{tuv} for t+u+v <= L only (Claude Generated, Sep 2026). Same recursion as
// buildRblockERI -- base R^n_{000} = (-2 rho)^n F_n(T), then z before y before x --
// but built level by level over the simplex t+u+v <= L-n instead of the full
// (L+1)^4 box, and only level n = 0 is kept. For the blocked ERI and gradient
// kernels this was the dominant cost (the box holds ~24x more entries than the
// simplex at L = 5). Output layout: R0[(t*(L+1) + u)*(L+1) + v].
static void buildR0Simplex(int L, double Wx, double Wy, double Wz, double rho,
                           const std::vector<double>& boys, std::vector<double>& R0)
{
    const int s = L + 1;
    const size_t sz = (size_t)s * s * s;
    thread_local std::vector<double> lvl;
    R0.assign(sz, 0.0);
    lvl.assign(sz, 0.0);
    auto at = [s](int t, int u, int v) { return ((size_t)t * s + u) * s + v; };
    const double b = -2.0 * rho;
    double bpow[64];  // L = LA+LB+LC+LD (+1 for the gradient); 63 is far beyond any basis here
    bpow[0] = 1.0;
    for (int n = 1; n <= L; ++n) bpow[n] = bpow[n - 1] * b;

    // cur = level n+1 (starts at n = L: only R_000), nxt = level n
    std::vector<double>* cur = &R0;
    std::vector<double>* nxt = &lvl;
    (*cur)[0] = bpow[L] * boys[L];
    for (int n = L - 1; n >= 0; --n) {
        std::vector<double>& c = *cur;
        std::vector<double>& x = *nxt;
        x[0] = bpow[n] * boys[n];
        const int top = L - n;
        for (int t = 0; t <= top; ++t)
            for (int u = 0; u <= top - t; ++u)
                for (int v = 0; v <= top - t - u; ++v) {
                    if (t + u + v == 0) continue;
                    double r;
                    if (v > 0)
                        r = (v >= 2 ? (v - 1) * c[at(t, u, v - 2)] : 0.0) + Wz * c[at(t, u, v - 1)];
                    else if (u > 0)
                        r = (u >= 2 ? (u - 1) * c[at(t, u - 2, v)] : 0.0) + Wy * c[at(t, u - 1, v)];
                    else
                        r = (t >= 2 ? (t - 1) * c[at(t - 2, u, v)] : 0.0) + Wx * c[at(t - 1, u, v)];
                    x[at(t, u, v)] = r;
                }
        std::swap(cur, nxt);
    }
    if (cur != &R0) R0.swap(*cur);
}

// Primitive (ga b gb | gc gd) for cartesian angular momenta (la,ma,na) on A,
// (lb,mb,nb) on B, (lc,mc,nc) on C, (ld,md,nd) on D. Exponents alpha..delta,
// centres A..D in Bohr. Returns the chemists' integral (ab|cd).
static double primitiveERI(int la, int ma, int na, int lb, int mb, int nb,
                           int lc, int mc, int nc, int ld, int md, int nd,
                           double alpha, double beta, double gamma, double delta,
                           double Ax, double Ay, double Az,
                           double Bx, double By, double Bz,
                           double Cx, double Cy, double Cz,
                           double Dx, double Dy, double Dz)
{
    // Bra pair (ab): product centre P, composite p, K_AB.
    double p, Px, Py, Pz;
    double Kab = gaussianProductK(alpha, beta, Ax, Ay, Az, Bx, By, Bz, p, Px, Py, Pz);
    if (Kab == 0.0) return 0.0;
    // Ket pair (cd): product centre Q, composite q, K_CD.
    double q, Qx, Qy, Qz;
    double Kcd = gaussianProductK(gamma, delta, Cx, Cy, Cz, Dx, Dy, Dz, q, Qx, Qy, Qz);
    if (Kcd == 0.0) return 0.0;

    // Hermite expansion coefficients. K is folded into the x-axis table only
    // (Ex[0][0][0] = K_AB, Eu/Ev start from 1.0), mirroring the WP1 nuclear
    // kernel convention. Same for the ket pair (K_CD into Ex2).
    auto Ex = hermiteCoeffs(la, lb, Px - Ax, Px - Bx, p, Kab);
    auto Eu = hermiteCoeffs(ma, mb, Py - Ay, Py - By, p, 1.0);
    auto Ev = hermiteCoeffs(na, nb, Pz - Az, Pz - Bz, p, 1.0);
    auto Ex2 = hermiteCoeffs(lc, ld, Qx - Cx, Qx - Dx, q, Kcd);
    auto Eu2 = hermiteCoeffs(mc, md, Qy - Cy, Qy - Dy, q, 1.0);
    auto Ev2 = hermiteCoeffs(nc, nd, Qz - Cz, Qz - Dz, q, 1.0);

    const int tmax = la + lb, umax = ma + mb, vmax = na + nb;
    const int taumax = lc + ld, upsmax = mc + md, ommax = nc + nd;
    // Combined Hermite orders and the maximum Boys order needed.
    const int Ttot = tmax + taumax, Utot = umax + upsmax, Vtot = vmax + ommax;
    const int Nmax = Ttot + Utot + Vtot;
    const double rho = p * q / (p + q);
    const double T = rho * ((Px - Qx) * (Px - Qx) + (Py - Qy) * (Py - Qy) + (Pz - Qz) * (Pz - Qz));
    auto F = boysArray(Nmax, T);

    std::vector<double> R((Ttot + 1) * (Utot + 1) * (Vtot + 1) * (Nmax + 1), 0.0);
    buildRblockERI(Ttot, Utot, Vtot, Nmax, Px - Qx, Py - Qy, Pz - Qz, rho, F, R);

    const int strideN = (Utot + 1) * (Vtot + 1) * (Nmax + 1);
    const int strideU = (Vtot + 1) * (Nmax + 1);
    const int strideV = (Nmax + 1);

    // Overall primitive prefactor (outside the R sum).
    const double pref = 2.0 * std::pow(PI, 2.5) / (p * q * std::sqrt(p + q));

    double val = 0.0;
    for (int t = 0; t <= tmax; ++t) {
        const double et = Ex[t][la][lb];
        if (et == 0.0) continue;
        for (int u = 0; u <= umax; ++u) {
            const double eu = Eu[u][ma][mb];
            if (eu == 0.0) continue;
            for (int v = 0; v <= vmax; ++v) {
                const double ev = Ev[v][na][nb];
                if (ev == 0.0) continue;
                const double bra = et * eu * ev;
                for (int tau = 0; tau <= taumax; ++tau) {
                    const double et2 = Ex2[tau][lc][ld];
                    if (et2 == 0.0) continue;
                    for (int ups = 0; ups <= upsmax; ++ups) {
                        const double eu2 = Eu2[ups][mc][md];
                        if (eu2 == 0.0) continue;
                        for (int om = 0; om <= ommax; ++om) {
                            const double ev2 = Ev2[om][nc][nd];
                            if (ev2 == 0.0) continue;
                            // (-1)^(tau+ups+om): ket Hermite sign.
                            const double sgn = ((tau + ups + om) & 1) ? -1.0 : 1.0;
                            const int Tidx = t + tau, Uidx = u + ups, Vidx = v + om;
                            const double r = R[Tidx * strideN + Uidx * strideU + Vidx * strideV + 0];
                            val += bra * et2 * eu2 * ev2 * sgn * r;
                        }
                    }
                }
            }
        }
    }
    return pref * val;
}

// Contracted (a b | c d) over the primitive contractions of the four AOs -- the
// textbook one-integral-at-a-time form. buildERI() below computes the same numbers
// shell-quartet-blocked; this function stays as the readable reference.
double contractedERI(const GTO::Orbital& a, const GTO::Orbital& b,
                                const GTO::Orbital& c, const GTO::Orbital& d)
{
    int la, ma, na, lb, mb, nb, lc, mc, nc, ld, md, nd;
    GTO::orbitalTypeToComponents(a.type, la, ma, na);
    GTO::orbitalTypeToComponents(b.type, lb, mb, nb);
    GTO::orbitalTypeToComponents(c.type, lc, mc, nc);
    GTO::orbitalTypeToComponents(d.type, ld, md, nd);
    double sum = 0.0;
    for (size_t ia = 0; ia < a.exponents.size(); ++ia) {
        const double ca = a.coefficients[ia];
        for (size_t ib = 0; ib < b.exponents.size(); ++ib) {
            const double cb = b.coefficients[ib];
            for (size_t ic = 0; ic < c.exponents.size(); ++ic) {
                const double cc = c.coefficients[ic];
                for (size_t id = 0; id < d.exponents.size(); ++id) {
                    sum += ca * cb * cc * d.coefficients[id] *
                           primitiveERI(la, ma, na, lb, mb, nb, lc, mc, nc, ld, md, nd,
                                        a.exponents[ia], b.exponents[ib],
                                        c.exponents[ic], d.exponents[id],
                                        a.x, a.y, a.z, b.x, b.y, b.z,
                                        c.x, c.y, c.z, d.x, d.y, d.z);
                }
            }
        }
    }
    return sum;
}

void ERITensor::set8(int mu, int nu, int lam, int sig, double v)
{
    // 8-fold chemists' symmetry: swap within bra {mu,nu}, within ket {lam,sig},
    // and swap the two pairs. Diagonal/identical index tuples just overwrite.
    at(mu, nu, lam, sig) = v;
    at(nu, mu, lam, sig) = v;
    at(mu, nu, sig, lam) = v;
    at(nu, mu, sig, lam) = v;
    at(lam, sig, mu, nu) = v;
    at(sig, lam, mu, nu) = v;
    at(lam, sig, nu, mu) = v;
    at(sig, lam, nu, mu) = v;
}

// ---------------------------------------------------------------------------
// Shell-quartet-blocked ERI build (Claude Generated, Sep 2026)
//
// The one-integral-at-a-time form above (contractedERI) recomputes, for every
// AO *component* quartet, the Gaussian products, all six Hermite tables, the
// Boys function and the whole R block -- although none of these depends on the
// cartesian powers, only on the shells. A d-shell quartet thus repeats the same
// primitive work 6^4 = 1296 times. This is the same pattern the native xTB
// overlap/multipole kernels had (see docs/SQM_PERFORMANCE.md, "shell-pair-blocked
// integral kernels"), and it is fixed the same way:
//
//   1. group the flat AO list into shells (same centre, exponents, L);
//   2. precompute, ONCE per shell pair and primitive pair, the product exponent,
//      product centre and the three Hermite tables up to the shell maxima
//      (a table built for (LA, LB) holds every lower (i, j) as well);
//   3. per shell quartet and primitive quartet, evaluate Boys + R once and
//      contract it with the tables for every component quartet.
//
// On top: Schwarz screening |(ab|cd)| <= sqrt((ab|ab)) sqrt((cd|cd)) at shell-pair
// level, and threads over bra shell pairs. Distinct canonical quartets write
// disjoint tensor entries, so the threads need no synchronisation. The layout is
// what a GPU port needs too: flat pair tables, one work item per quartet.
// ---------------------------------------------------------------------------

namespace {

struct EriShell {
    int first = 0;                              ///< first AO index in the flat basis
    int L = 0;                                  ///< total angular momentum
    double x = 0, y = 0, z = 0;                 ///< centre (Bohr)
    std::vector<double> exps;                   ///< primitive exponents
    std::vector<std::array<int, 3>> lmn;        ///< cartesian powers per component
    std::vector<std::vector<double>> coef;      ///< [component][primitive], pre-normalised
};

std::vector<EriShell> groupShells(const std::vector<GTO::Orbital>& basis)
{
    std::vector<EriShell> shells;
    const int n = static_cast<int>(basis.size());
    for (int i = 0; i < n; ++i) {
        int l, m, k;
        GTO::orbitalTypeToComponents(basis[i].type, l, m, k);
        const int L = l + m + k;
        bool join = false;
        if (!shells.empty()) {
            const EriShell& s = shells.back();
            const GTO::Orbital& prev = basis[i - 1];
            join = s.L == L && prev.atom == basis[i].atom && prev.exponents == basis[i].exponents
                && prev.x == basis[i].x && prev.y == basis[i].y && prev.z == basis[i].z
                && static_cast<int>(s.lmn.size()) < (L + 1) * (L + 2) / 2;
        }
        if (!join) {
            EriShell s;
            s.first = i;
            s.L = L;
            s.x = basis[i].x; s.y = basis[i].y; s.z = basis[i].z;
            s.exps = basis[i].exponents;
            shells.push_back(std::move(s));
        }
        shells.back().lmn.push_back({ l, m, k });
        shells.back().coef.push_back(basis[i].coefficients);
    }
    return shells;
}

/// One primitive pair of a shell pair: everything the quartet loop needs.
struct PrimPair {
    int ia = 0, ib = 0;
    double p = 0, Px = 0, Py = 0, Pz = 0;
    std::vector<double> Ex, Ey, Ez;  ///< flat [t][i][j], dims (LA+LB+1)(LA+1)(LB+1); K folded into Ex
};

struct ShellPair {
    int A = 0, B = 0;
    std::vector<PrimPair> prims;
    double schwarz = 0.0;  ///< max over components of sqrt(|(ab|ab)|)
};

std::vector<double> flattenHermite(const std::vector<std::vector<std::vector<double>>>& E)
{
    std::vector<double> f;
    for (const auto& t : E)
        for (const auto& i : t)
            f.insert(f.end(), i.begin(), i.end());
    return f;
}

// raiseA/raiseB = 1 build the Hermite tables one angular-momentum step higher on
// that shell -- what the gradient kernel needs for the (l+1) half of d/dA g_l.
ShellPair makeShellPair(const std::vector<EriShell>& sh, int A, int B, int raiseA = 0, int raiseB = 0)
{
    ShellPair sp;
    sp.A = A;
    sp.B = B;
    const EriShell& a = sh[A];
    const EriShell& b = sh[B];
    for (size_t ia = 0; ia < a.exps.size(); ++ia)
        for (size_t ib = 0; ib < b.exps.size(); ++ib) {
            PrimPair pp;
            pp.ia = static_cast<int>(ia);
            pp.ib = static_cast<int>(ib);
            const double K = gaussianProductK(a.exps[ia], b.exps[ib], a.x, a.y, a.z, b.x, b.y, b.z,
                                              pp.p, pp.Px, pp.Py, pp.Pz);
            if (K == 0.0) continue;  // underflow: contributes exactly 0, as in primitiveERI
            const int LA = a.L + raiseA, LB = b.L + raiseB;
            pp.Ex = flattenHermite(hermiteCoeffs(LA, LB, pp.Px - a.x, pp.Px - b.x, pp.p, K));
            pp.Ey = flattenHermite(hermiteCoeffs(LA, LB, pp.Py - a.y, pp.Py - b.y, pp.p, 1.0));
            pp.Ez = flattenHermite(hermiteCoeffs(LA, LB, pp.Pz - a.z, pp.Pz - b.z, pp.p, 1.0));
            sp.prims.push_back(std::move(pp));
        }
    return sp;
}

/// Contracted (AB|CD) for every component quartet of four shells, written to
/// out[((ca*ncB + cb)*ncC + cc)*ncD + cd]. Per primitive quartet the Boys
/// function and R block are built once for all components.
void shellQuartet(const std::vector<EriShell>& sh, const ShellPair& bra, const ShellPair& ket,
                  std::vector<double>& out, std::vector<double>& R)
{
    const EriShell& A = sh[bra.A];
    const EriShell& B = sh[bra.B];
    const EriShell& C = sh[ket.A];
    const EriShell& D = sh[ket.B];
    const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
    out.assign((size_t)nA * nB * nC * nD, 0.0);

    const int LAB = A.L + B.L, LCD = C.L + D.L;
    const int Ltot = LAB + LCD;
    // Hermite table strides: E[t][i][j] with i <= LA, j <= LB.
    const int sAB_t = (A.L + 1) * (B.L + 1), sAB_i = B.L + 1;
    const int sCD_t = (C.L + 1) * (D.L + 1), sCD_i = D.L + 1;
    // R block strides (every Cartesian direction up to Ltot).
    // R^0 on the simplex, layout [(t*(Ltot+1) + u)*(Ltot+1) + v]
    const int strideV = 1;
    const int strideU = Ltot + 1;
    const int strideN = (Ltot + 1) * strideU;

    for (const PrimPair& pb : bra.prims) {
        for (const PrimPair& pk : ket.prims) {
            const double p = pb.p, q = pk.p;
            const double rho = p * q / (p + q);
            const double Wx = pb.Px - pk.Px, Wy = pb.Py - pk.Py, Wz = pb.Pz - pk.Pz;
            const double T = rho * (Wx * Wx + Wy * Wy + Wz * Wz);
            thread_local std::vector<double> F;
            boysArrayInto(Ltot, T, F);
            buildR0Simplex(Ltot, Wx, Wy, Wz, rho, F, R);
            const double pref = 2.0 * std::pow(PI, 2.5) / (p * q * std::sqrt(p + q));

            for (int ca = 0; ca < nA; ++ca) {
                const double c_a = A.coef[ca][pb.ia];
                const auto& la = A.lmn[ca];
                for (int cb = 0; cb < nB; ++cb) {
                    const double c_ab = c_a * B.coef[cb][pb.ib];
                    const auto& lb = B.lmn[cb];
                    const int tmax = la[0] + lb[0], umax = la[1] + lb[1], vmax = la[2] + lb[2];
                    for (int cc = 0; cc < nC; ++cc) {
                        const auto& lc = C.lmn[cc];
                        for (int cd = 0; cd < nD; ++cd) {
                            const auto& ld = D.lmn[cd];
                            const double cprod = c_ab * C.coef[cc][pk.ia] * D.coef[cd][pk.ib];
                            const int taumax = lc[0] + ld[0], upsmax = lc[1] + ld[1], ommax = lc[2] + ld[2];
                            double val = 0.0;
                            for (int t = 0; t <= tmax; ++t) {
                                const double et = pb.Ex[t * sAB_t + la[0] * sAB_i + lb[0]];
                                if (et == 0.0) continue;
                                for (int u = 0; u <= umax; ++u) {
                                    const double eu = pb.Ey[u * sAB_t + la[1] * sAB_i + lb[1]];
                                    if (eu == 0.0) continue;
                                    for (int v = 0; v <= vmax; ++v) {
                                        const double ev = pb.Ez[v * sAB_t + la[2] * sAB_i + lb[2]];
                                        if (ev == 0.0) continue;
                                        const double braE = et * eu * ev;
                                        for (int tau = 0; tau <= taumax; ++tau) {
                                            const double et2 = pk.Ex[tau * sCD_t + lc[0] * sCD_i + ld[0]];
                                            if (et2 == 0.0) continue;
                                            for (int ups = 0; ups <= upsmax; ++ups) {
                                                const double eu2 = pk.Ey[ups * sCD_t + lc[1] * sCD_i + ld[1]];
                                                if (eu2 == 0.0) continue;
                                                for (int om = 0; om <= ommax; ++om) {
                                                    const double ev2 = pk.Ez[om * sCD_t + lc[2] * sCD_i + ld[2]];
                                                    if (ev2 == 0.0) continue;
                                                    const double sgn = ((tau + ups + om) & 1) ? -1.0 : 1.0;
                                                    val += braE * et2 * eu2 * ev2 * sgn
                                                         * R[(t + tau) * strideN + (u + ups) * strideU + (v + om) * strideV];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            out[((size_t)(ca * nB + cb) * nC + cc) * nD + cd] += cprod * pref * val;
                        }
                    }
                }
            }
        }
    }
}

/// Where every shell lands in the ACTIVE basis: first output index, width, and
/// (spherical only) the cartesian -> output transform T (ncomp x width) cut out
/// of Q. Q is block-diagonal per shell, so the transform can be applied to each
/// shell quartet right after it is computed -- the full cartesian tensor is never
/// formed. For s/p shells T is the identity.
struct ShellOutput {
    int n = 0;                          ///< active basis dimension
    bool spherical = false;
    std::vector<int> first, width;
    std::vector<Eigen::MatrixXd> T;
    std::vector<bool> identity;
    ~ShellOutput();  // out of line: the implicit inline one trips -Winline
};
ShellOutput::~ShellOutput() = default;

ShellOutput makeShellOutput(const std::vector<EriShell>& sh, int ncart, const Matrix* Q)
{
    ShellOutput o;
    o.spherical = (Q != nullptr && Q->size() != 0);
    o.n = o.spherical ? (int)Q->cols() : ncart;
    const int ns = (int)sh.size();
    o.first.assign(ns, 0);
    o.width.assign(ns, 0);
    o.T.assign(ns, Eigen::MatrixXd());
    o.identity.assign(ns, true);
    for (int A = 0; A < ns; ++A) {
        const int nc = (int)sh[A].lmn.size();
        if (!o.spherical) { o.first[A] = sh[A].first; o.width[A] = nc; continue; }
        int lo = o.n, hi = -1;
        for (int j = 0; j < o.n; ++j)
            for (int a = 0; a < nc; ++a)
                if ((*Q)(sh[A].first + a, j) != 0.0) { lo = std::min(lo, j); hi = std::max(hi, j); }
        o.first[A] = lo;
        o.width[A] = hi - lo + 1;
        o.T[A] = Eigen::MatrixXd::Zero(nc, o.width[A]);
        for (int a = 0; a < nc; ++a)
            for (int j = 0; j < o.width[A]; ++j)
                o.T[A](a, j) = (*Q)(sh[A].first + a, lo + j);
        o.identity[A] = (nc == o.width[A]) && o.T[A].isIdentity(0.0);
    }
    return o;
}

/// Contract one index of a 4-index block (dims d[0..3]) with Tk, in place.
void transformAxis(std::vector<double>& blk, std::vector<double>& tmp, int d[4], int k,
                   const Eigen::MatrixXd& Tk)
{
    size_t outer = 1, inner = 1;
    for (int m = 0; m < k; ++m) outer *= d[m];
    for (int m = k + 1; m < 4; ++m) inner *= d[m];
    const int nin = d[k], nout = (int)Tk.cols();
    tmp.assign(outer * nout * inner, 0.0);
    for (size_t o = 0; o < outer; ++o)
        for (int a = 0; a < nin; ++a)
            for (int j = 0; j < nout; ++j) {
                const double t = Tk(a, j);
                if (t == 0.0) continue;
                const double* x = &blk[(o * nin + a) * inner];
                double* y = &tmp[(o * nout + j) * inner];
                for (size_t r = 0; r < inner; ++r) y[r] += t * x[r];
            }
    blk.swap(tmp);
    d[k] = nout;
}

/// One shell quartet in the ACTIVE basis: buf holds the block with dims d[0..3].
void activeQuartet(const std::vector<EriShell>& sh, const ShellOutput& o, const ShellPair& bra,
                   const ShellPair& ket, std::vector<double>& buf, std::vector<double>& R,
                   std::vector<double>& tmp, int d[4])
{
    shellQuartet(sh, bra, ket, buf, R);
    d[0] = (int)sh[bra.A].lmn.size();
    d[1] = (int)sh[bra.B].lmn.size();
    d[2] = (int)sh[ket.A].lmn.size();
    d[3] = (int)sh[ket.B].lmn.size();
    if (!o.spherical) return;
    const int s4[4] = { bra.A, bra.B, ket.A, ket.B };
    for (int k = 0; k < 4; ++k)
        if (!o.identity[s4[k]])
            transformAxis(buf, tmp, d, k, o.T[s4[k]]);
}

/// Shell pairs A <= B with their primitive-pair tables and Schwarz factors
/// sqrt(max |(ab|ab)|) from the diagonal quartets.
std::vector<ShellPair> makeScreeningPairs(const std::vector<EriShell>& sh)
{
    const int ns = (int)sh.size();
    std::vector<ShellPair> pairs;
    pairs.reserve((size_t)ns * (ns + 1) / 2);
    for (int A = 0; A < ns; ++A)
        for (int B = A; B < ns; ++B)
            pairs.push_back(makeShellPair(sh, A, B));
    std::vector<double> buf, R;
    for (ShellPair& sp : pairs) {
        shellQuartet(sh, sp, sp, buf, R);
        const int nA = (int)sh[sp.A].lmn.size(), nB = (int)sh[sp.B].lmn.size();
        double mx = 0.0;
        for (int a = 0; a < nA; ++a)
            for (int b = 0; b < nB; ++b)
                mx = std::max(mx, std::abs(buf[((size_t)(a * nB + b) * nA + a) * nB + b]));
        sp.schwarz = std::sqrt(mx);
    }
    return pairs;
}

}  // namespace

ERITensor buildERI(const std::vector<GTO::Orbital>& basis, int threads, double screening,
                   const Matrix* Q)
{
    const std::vector<EriShell> sh = groupShells(basis);
    const ShellOutput o = makeShellOutput(sh, (int)basis.size(), Q);
    ERITensor eri(o.n);
    if (o.n == 0) return eri;

    const std::vector<ShellPair> pairs = makeScreeningPairs(sh);
    const int np = (int)pairs.size();

    // Canonical shell quartets: pair(AB) <= pair(CD). Every AO quartet of a shell
    // quartet is written with set8, so A==B or AB==CD blocks just rewrite the same
    // value into already-written entries.
#ifdef _OPENMP
#pragma omp parallel num_threads(threads > 0 ? threads : 1)
#endif
    {
        std::vector<double> buf, R, tmp;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int i = 0; i < np; ++i) {
            const ShellPair& bra = pairs[i];
            for (int j = i; j < np; ++j) {
                const ShellPair& ket = pairs[j];
                if (bra.schwarz * ket.schwarz < screening) continue;
                int d[4];
                activeQuartet(sh, o, bra, ket, buf, R, tmp, d);
                const int oA = o.first[bra.A], oB = o.first[bra.B], oC = o.first[ket.A], oD = o.first[ket.B];
                for (int a = 0; a < d[0]; ++a)
                    for (int b = 0; b < d[1]; ++b)
                        for (int c = 0; c < d[2]; ++c)
                            for (int e = 0; e < d[3]; ++e)
                                eri.set8(oA + a, oB + b, oC + c, oD + e,
                                         buf[((size_t)(a * d[1] + b) * d[2] + c) * d[3] + e]);
            }
        }
    }
    (void)threads;
    return eri;
}

// ---------------------------------------------------------------------------
// Integral-direct J/K (Claude Generated, Sep 2026)
//
// Instead of storing (mu nu|lam sig) (n^4 doubles), every screened canonical shell
// quartet is recomputed per Fock build and digested straight into J and K
// (Almloef, Faegri, Korsell, J. Comput. Chem. 3, 385 (1982)). Each canonical
// quartet stands for up to 8 index orderings; it is applied once with the weight
// w = v * deg/8, deg = (A!=B ? 2:1)(C!=D ? 2:1)(AB!=CD ? 2:1), and the 8
// orderings are recovered by
//     Jt_ab += 4 w P_cd,  Jt_cd += 4 w P_ab,
//     Kt_ac += 2 w P_bd,  Kt_bd += 2 w P_ac,  Kt_ad += 2 w P_bc,  Kt_bc += 2 w P_ad,
//     J = (Jt + Jt^T)/2,  K = (Kt + Kt^T)/2
// (the same scheme as the libint Hartree-Fock example). Screening is
// density-weighted (Haeser, Ahlrichs, J. Comput. Chem. 10, 104 (1989)): a
// quartet is skipped when Q_AB Q_CD max(|P| over the six shell blocks it touches)
// is below the threshold -- so in an incremental build on dP most quartets drop
// out as the SCF converges.
// ---------------------------------------------------------------------------

struct DirectJK::Impl {
    std::vector<EriShell> sh;
    ShellOutput out;
    std::vector<ShellPair> pairs;
    double max_schwarz = 0.0;
};

DirectJK::DirectJK(const std::vector<GTO::Orbital>& basis, double screening, const Matrix* Q)
    : m_screening(screening)
{
    auto impl = std::make_shared<Impl>();
    impl->sh = groupShells(basis);
    impl->out = makeShellOutput(impl->sh, (int)basis.size(), Q);
    impl->pairs = makeScreeningPairs(impl->sh);
    for (const ShellPair& sp : impl->pairs)
        impl->max_schwarz = std::max(impl->max_schwarz, sp.schwarz);
    m_impl = impl;
}

int DirectJK::n() const { return m_impl ? m_impl->out.n : 0; }

long DirectJK::build(const Matrix& P, Matrix& J, Matrix& K, int threads) const
{
    const int n = this->n();
    J = Matrix::Zero(n, n);
    K = Matrix::Zero(n, n);
    if (!m_impl || n == 0) return 0;
    const Impl& I = *m_impl;
    const int ns = (int)I.sh.size();
    const int np = (int)I.pairs.size();
    const ShellOutput& o = I.out;

    // Shell-block maxima of |P| for the density-weighted screening.
    Eigen::MatrixXd Pmax = Eigen::MatrixXd::Zero(ns, ns);
    for (int A = 0; A < ns; ++A)
        for (int B = A; B < ns; ++B) {
            const double m = P.block(o.first[A], o.first[B], o.width[A], o.width[B]).cwiseAbs().maxCoeff();
            Pmax(A, B) = Pmax(B, A) = m;
        }
    const double Pall = Pmax.maxCoeff();
    const double thr = m_screening;
    long computed = 0;

#ifdef _OPENMP
#pragma omp parallel num_threads(threads > 0 ? threads : 1) reduction(+ : computed)
#endif
    {
        Matrix Jt = Matrix::Zero(n, n), Kt = Matrix::Zero(n, n);
        std::vector<double> buf, R, tmp;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int i = 0; i < np; ++i) {
            const ShellPair& bra = I.pairs[i];
            if (bra.schwarz * I.max_schwarz * Pall < thr) continue;
            const int A = bra.A, B = bra.B;
            for (int j = i; j < np; ++j) {
                const ShellPair& ket = I.pairs[j];
                const int C = ket.A, D = ket.B;
                const double pq = std::max({ Pmax(A, B), Pmax(C, D), Pmax(A, C), Pmax(A, D), Pmax(B, C), Pmax(B, D) });
                if (bra.schwarz * ket.schwarz * pq < thr) continue;
                int d[4];
                activeQuartet(I.sh, o, bra, ket, buf, R, tmp, d);
                ++computed;
                const double deg = (A != B ? 2.0 : 1.0) * (C != D ? 2.0 : 1.0) * (i != j ? 2.0 : 1.0);
                const double scale = deg / 8.0;
                const int oA = o.first[A], oB = o.first[B], oC = o.first[C], oD = o.first[D];
                for (int a = 0; a < d[0]; ++a) {
                    const int ia = oA + a;
                    for (int b = 0; b < d[1]; ++b) {
                        const int ib = oB + b;
                        const double Pab = P(ia, ib);
                        double jab = 0.0;
                        for (int c = 0; c < d[2]; ++c) {
                            const int ic = oC + c;
                            const double Pac = P(ia, ic), Pbc = P(ib, ic);
                            const double* row = &buf[((size_t)(a * d[1] + b) * d[2] + c) * d[3]];
                            double kac = 0.0, kbc = 0.0;
                            for (int e = 0; e < d[3]; ++e) {
                                const double w = scale * row[e];
                                if (w == 0.0) continue;
                                const int id = oD + e;
                                jab += w * P(ic, id);
                                Jt(ic, id) += 4.0 * w * Pab;
                                kac += w * P(ib, id);
                                kbc += w * P(ia, id);
                                Kt(ib, id) += 2.0 * w * Pac;
                                Kt(ia, id) += 2.0 * w * Pbc;
                            }
                            Kt(ia, ic) += 2.0 * kac;
                            Kt(ib, ic) += 2.0 * kbc;
                        }
                        Jt(ia, ib) += 4.0 * jab;
                    }
                }
            }
        }
#ifdef _OPENMP
#pragma omp critical(qm_direct_jk)
#endif
        {
            J += Jt;
            K += Kt;
        }
    }
    (void)threads;
    const Matrix Jsym = 0.5 * (J + J.transpose());
    const Matrix Ksym = 0.5 * (K + K.transpose());
    J = Jsym;
    K = Ksym;
    return computed;
}

// J_munu = sum_{lam,sig} (mu nu | lam sig) P_lamsig: with the tensor stored as an
// (n^2 x n^2) row-major matrix this is one matrix-vector product, J = ERI * vec(P).
// Written as one dot product per (mu, nu) row so it threads trivially (and maps
// 1:1 onto a cuBLAS GEMV for a GPU port).
Matrix buildCoulomb(const ERITensor& eri, const Matrix& P, int threads)
{
    const int n = eri.n();
    // Below ~32 functions the OpenMP start-up costs more than the O(n^4) work
    // (H2O/def2-SVP: 0.6 ms serial vs 3.1 ms on 4 threads).
    if (n < 32) threads = 1;
    Matrix J = Matrix::Zero(n, n);
    const Matrix Pc = P;  // row-major, contiguous
    const Eigen::Map<const Eigen::VectorXd> pvec(Pc.data(), (Eigen::Index)n * n);
    const double* base = eri.data();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads > 0 ? threads : 1) schedule(static)
#endif
    for (int mu = 0; mu < n; ++mu)
        for (int nu = 0; nu < n; ++nu) {
            const Eigen::Map<const Eigen::VectorXd> row(base + ((size_t)mu * n + nu) * n * n, (Eigen::Index)n * n);
            J(mu, nu) = row.dot(pvec);
        }
    (void)threads;
    return J;
}

// K_munu = sum_{lam,sig} (mu lam | nu sig) P_lamsig. For fixed (mu, lam) the block
// B(nu, sig) = (mu lam | nu sig) is a contiguous n x n row-major matrix, so
// K.row(mu) += (B * P.row(lam)^T)^T -- one GEMV per (mu, lam).
Matrix buildExchange(const ERITensor& eri, const Matrix& P, int threads)
{
    const int n = eri.n();
    // Below ~32 functions the OpenMP start-up costs more than the O(n^4) work
    // (H2O/def2-SVP: 0.6 ms serial vs 3.1 ms on 4 threads).
    if (n < 32) threads = 1;
    Matrix K = Matrix::Zero(n, n);
    const Matrix Pc = P;
    const double* base = eri.data();
    using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads > 0 ? threads : 1) schedule(static)
#endif
    for (int mu = 0; mu < n; ++mu) {
        Eigen::VectorXd acc = Eigen::VectorXd::Zero(n);
        for (int lam = 0; lam < n; ++lam) {
            const Eigen::Map<const RowMat> B(base + ((size_t)mu * n + lam) * n * n, n, n);
            acc.noalias() += B * Pc.row(lam).transpose();
        }
        K.row(mu) = acc.transpose();
    }
    (void)threads;
    return K;
}

// 4-index transform ERI_sph[i,j,k,l] = sum Q[a,i] Q[b,j] Q[c,k] Q[d,l] ERI[a,b,c,d]
// as four one-index transforms (O(n^5) instead of the O(n^8) direct sum), each
// using only the nonzeros of Q. Q is block-diagonal per shell -- s and p columns
// have one entry, a 5d column mixes the six cartesian d components of its own
// shell -- so a pass costs about n^4 multiply-adds, not n^5.
// Axis k of a tensor with dims (n0,n1,n2,n3) is contracted as
//   out[o][i][in] = sum_{(a,q) in column i} q * X[o][a][in],
// with o the indices before k and in the ones after (a contiguous run of length
// `inner`, so the update is a vectorisable axpy).
ERITensor applySphericalTransformERI(const ERITensor& eriCart, const Matrix& Q)
{
    const int ncart = eriCart.n();
    if (Q.size() == 0 || ncart == 0) return eriCart;  // no d -> keep cartesian
    const int nsph = (int)Q.cols();

    // Nonzero pattern of Q per spherical column.
    std::vector<std::vector<std::pair<int, double>>> col(nsph);
    for (int i = 0; i < nsph; ++i)
        for (int a = 0; a < ncart; ++a)
            if (Q(a, i) != 0.0) col[i].push_back({ a, Q(a, i) });

    std::vector<double> cur(eriCart.data(), eriCart.data() + (size_t)ncart * ncart * ncart * ncart);
    int dims[4] = { ncart, ncart, ncart, ncart };
    for (int k = 0; k < 4; ++k) {
        size_t outer = 1, inner = 1;
        for (int m = 0; m < k; ++m) outer *= dims[m];
        for (int m = k + 1; m < 4; ++m) inner *= dims[m];
        const int nin = dims[k];
        std::vector<double> next(outer * nsph * inner, 0.0);
        for (size_t o = 0; o < outer; ++o) {
            const double* X = cur.data() + o * nin * inner;
            double* Y = next.data() + o * nsph * inner;
            for (int i = 0; i < nsph; ++i) {
                double* yi = Y + (size_t)i * inner;
                for (const auto& aq : col[i]) {
                    const double* xa = X + (size_t)aq.first * inner;
                    const double q = aq.second;
                    for (size_t r = 0; r < inner; ++r) yi[r] += q * xa[r];
                }
            }
        }
        cur.swap(next);
        dims[k] = nsph;
    }
    ERITensor eriSph(nsph);
    std::copy(cur.begin(), cur.end(), eriSph.data());
    return eriSph;
}

// ===========================================================================
// Analytic nuclear gradient of the closed-shell RHF energy (WP8, Sep 2026)
// Claude Generated.
//
//   E = sum P_mn H_mn + 1/2 sum D_mnls (mn|ls) + E_nn,
//   D_mnls = P_mn P_ls - 1/4 (P_ml P_ns + P_ms P_nl)   (spin-summed P)
//
//   dE/dA = sum P_mn dH_mn/dA - sum W_mn dS_mn/dA + 1/2 sum D (mn|ls)^A + dE_nn/dA,
//   W_mn  = 2 sum_i^occ eps_i C_mi C_ni  (energy-weighted density; the -W dS term
//           is the orbital-orthonormality ("Pulay") contribution).
//   Pople, Krishnan, Schlegel, Binkley, Int. J. Quantum Chem. S13, 225 (1979);
//   Helgaker, Jorgensen, Olsen, Molecular Electronic-Structure Theory, ch. 9.
//
// Every derivative integral comes from the centre derivative of a cartesian
// Gaussian,  d/dA_x [x_A^l e^{-a r_A^2}] = 2a x_A^{l+1} e^{..} - l x_A^{l-1} e^{..},
// i.e. from ordinary integrals with the power on A shifted by +-1 (the
// normalisation constant of the original function stays attached). Only the
// derivative with respect to the centre of the FIRST function is ever needed:
//  - S, T, V (basis part): the matrices are symmetric, so d/dA of <m|O|n> summed
//    against a symmetric P gives 2 sum_{m on A} P_mn <dm|O|n>;
//  - V (operator part, the nucleus C itself moving): translational invariance,
//    d/dC <m|1/r_C|n> = -(d/dA_m + d/dB_n) <m|1/r_C|n>;
//  - ERI: D and (mn|ls) share the 8-fold symmetry, so the four centre
//    derivatives contribute equally: 1/2 sum D (mn|ls)^A = 2 sum_{m on A} D (dm n|ls).
// ===========================================================================

// d/dA_k of one primitive quantity f(l,m,n) evaluated with the powers of the
// function on A: 2 alpha f(l+1_k) - l_k f(l-1_k).
template <typename F>
static inline void centreDerivative(int l, int m, int n, double alpha, F f, double out[3])
{
    out[0] = 2.0 * alpha * f(l + 1, m, n) - (l > 0 ? l * f(l - 1, m, n) : 0.0);
    out[1] = 2.0 * alpha * f(l, m + 1, n) - (m > 0 ? m * f(l, m - 1, n) : 0.0);
    out[2] = 2.0 * alpha * f(l, m, n + 1) - (n > 0 ? n * f(l, m, n - 1) : 0.0);
}

Matrix gradientOneElectron(const std::vector<GTO::Orbital>& basis,
                           const std::vector<int>& atomZ, const Matrix& atomPosBohr,
                           const Matrix& P, const Matrix& W, int threads)
{
    const int n = (int)basis.size();
    const int natoms = (int)atomZ.size();
    Matrix grad = Matrix::Zero(natoms, 3);
#ifdef _OPENMP
#pragma omp parallel num_threads(threads > 0 ? threads : 1)
#endif
    {
        Matrix g = Matrix::Zero(natoms, 3);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int mu = 0; mu < n; ++mu) {
            const GTO::Orbital& a = basis[mu];
            int l1, m1, n1;
            GTO::orbitalTypeToComponents(a.type, l1, m1, n1);
            for (int nu = 0; nu < n; ++nu) {
                const GTO::Orbital& b = basis[nu];
                int l2, m2, n2;
                GTO::orbitalTypeToComponents(b.type, l2, m2, n2);
                const double pmn = P(mu, nu), wmn = W(mu, nu);
                double dS[3] = { 0, 0, 0 }, dT[3] = { 0, 0, 0 };
                std::vector<double> dV(3 * natoms, 0.0);  // <dm| -Z_C/r_C |n> per nucleus C
                for (size_t ia = 0; ia < a.exponents.size(); ++ia) {
                    const double alpha = a.exponents[ia];
                    for (size_t ib = 0; ib < b.exponents.size(); ++ib) {
                        const double beta = b.exponents[ib];
                        const double cc = a.coefficients[ia] * b.coefficients[ib];
                        double d[3];
                        centreDerivative(l1, m1, n1, alpha, [&](int x, int y, int z) {
                            return primitiveOverlap(x, y, z, l2, m2, n2, alpha, beta,
                                                    a.x, a.y, a.z, b.x, b.y, b.z); }, d);
                        for (int k = 0; k < 3; ++k) dS[k] += cc * d[k];
                        centreDerivative(l1, m1, n1, alpha, [&](int x, int y, int z) {
                            return primitiveKinetic(x, y, z, l2, m2, n2, alpha, beta,
                                                    a.x, a.y, a.z, b.x, b.y, b.z); }, d);
                        for (int k = 0; k < 3; ++k) dT[k] += cc * d[k];
                        for (int C = 0; C < natoms; ++C) {
                            const double Cx = atomPosBohr(C, 0), Cy = atomPosBohr(C, 1), Cz = atomPosBohr(C, 2);
                            centreDerivative(l1, m1, n1, alpha, [&](int x, int y, int z) {
                                return primitiveNuclearAttraction(x, y, z, l2, m2, n2, alpha, beta,
                                                                  a.x, a.y, a.z, b.x, b.y, b.z, Cx, Cy, Cz); }, d);
                            for (int k = 0; k < 3; ++k) dV[3 * C + k] += -atomZ[C] * cc * d[k];
                        }
                    }
                }
                const int A = a.atom;
                for (int k = 0; k < 3; ++k) {
                    double dVtot = 0.0;
                    for (int C = 0; C < natoms; ++C) {
                        dVtot += dV[3 * C + k];
                        // operator part: nucleus C moves (translational invariance)
                        g(C, k) -= 2.0 * pmn * dV[3 * C + k];
                    }
                    g(A, k) += 2.0 * pmn * (dT[k] + dVtot) - 2.0 * wmn * dS[k];
                }
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        grad += g;
    }
    (void)threads;
    return grad;
}

namespace {

/// Cartesian power triples of total angular momentum L (any order; used as a set).
std::vector<std::array<int, 3>> cartTriples(int L)
{
    std::vector<std::array<int, 3>> t;
    if (L < 0) return t;
    for (int x = L; x >= 0; --x)
        for (int y = L - x; y >= 0; --y)
            t.push_back({ x, y, L - x - y });
    return t;
}

/// For one canonical shell quartet (AB|CD): sum over components of
///   D_abcd * d(ab|cd)/dA,  .../dB,  .../dC     (3 components each),
/// the derivative with respect to D following from translational invariance.
/// `bra`/`ket` pairs carry Hermite tables raised by one on BOTH shells
/// (makeShellPair(..., 1, 1)); D weights are Dq[((a*nB+b)*nC+c)*nD+d].
/// Claude Generated (Sep 2026).
void quartetGradient(const std::vector<EriShell>& sh, const ShellPair& bra, const ShellPair& ket,
                     const std::vector<double>& Dq, double gA[3], double gB[3], double gC[3])
{
    const EriShell& A = sh[bra.A];
    const EriShell& B = sh[bra.B];
    const EriShell& C = sh[ket.A];
    const EriShell& D = sh[ket.B];
    const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nD = (int)D.lmn.size();
    for (int k = 0; k < 3; ++k) gA[k] = gB[k] = gC[k] = 0.0;

    // Tables were built with (L+1) on both shells of each pair.
    const int sAB_t = (A.L + 2) * (B.L + 2), sAB_i = B.L + 2;
    const int sCD_t = (C.L + 2) * (D.L + 2), sCD_i = D.L + 2;
    const int Ltot = A.L + B.L + C.L + D.L + 1;     // one centre differentiated at a time
    const int Lb = A.L + B.L + 1;                   // highest bra Hermite order needed
    const int hs = Lb + 1, hsz = hs * hs * hs;
    const int strideU = Ltot + 1, strideN = (Ltot + 1) * strideU;  // R^0 simplex layout

    // Ket powers on C that the C derivative needs: totals LC-1, LC, LC+1.
    const int cdim = C.L + 2;
    auto ccode = [cdim](int x, int y, int z) { return (x * cdim + y) * cdim + z; };
    std::vector<std::array<int, 3>> cpow;
    for (int L = C.L - 1; L <= C.L + 1; ++L)
        for (const auto& t : cartTriples(L)) cpow.push_back(t);

    thread_local std::vector<double> R, F, hbuf;
    hbuf.resize((size_t)cdim * cdim * cdim * nD * hsz);
    auto hptr = [&](int x, int y, int z, int d) { return &hbuf[((size_t)ccode(x, y, z) * nD + d) * hsz]; };

    for (const PrimPair& pb : bra.prims) {
        const double alpha = A.exps[pb.ia], beta = B.exps[pb.ib];
        for (const PrimPair& pk : ket.prims) {
            const double gamma = C.exps[pk.ia];
            const double p = pb.p, q = pk.p;
            const double rho = p * q / (p + q);
            const double Wx = pb.Px - pk.Px, Wy = pb.Py - pk.Py, Wz = pb.Pz - pk.Pz;
            boysArrayInto(Ltot, rho * (Wx * Wx + Wy * Wy + Wz * Wz), F);
            buildR0Simplex(Ltot, Wx, Wy, Wz, rho, F, R);
            const double pref = 2.0 * std::pow(PI, 2.5) / (p * q * std::sqrt(p + q));

            // h(t,u,v) = sum over ket Hermite terms of (-1)^(..) E E E R, for every ket
            // power combination (c', d) any derivative can reach.
            for (const auto& c : cpow)
                for (int d = 0; d < nD; ++d) {
                    const auto& ld = D.lmn[d];
                    double* h = hptr(c[0], c[1], c[2], d);
                    const int top = Ltot - (c[0] + c[1] + c[2] + ld[0] + ld[1] + ld[2]);
                    const int taumax = c[0] + ld[0], upsmax = c[1] + ld[1], ommax = c[2] + ld[2];
                    for (int t = 0; t <= std::min(top, Lb); ++t)
                        for (int u = 0; u <= std::min(top, Lb) - t; ++u)
                            for (int v = 0; v <= std::min(top, Lb) - t - u; ++v) {
                                double acc = 0.0;
                                for (int tau = 0; tau <= taumax; ++tau) {
                                    const double e1 = pk.Ex[tau * sCD_t + c[0] * sCD_i + ld[0]];
                                    if (e1 == 0.0) continue;
                                    for (int ups = 0; ups <= upsmax; ++ups) {
                                        const double e2 = pk.Ey[ups * sCD_t + c[1] * sCD_i + ld[1]];
                                        if (e2 == 0.0) continue;
                                        for (int om = 0; om <= ommax; ++om) {
                                            const double e3 = pk.Ez[om * sCD_t + c[2] * sCD_i + ld[2]];
                                            if (e3 == 0.0) continue;
                                            const double sgn = ((tau + ups + om) & 1) ? -1.0 : 1.0;
                                            acc += e1 * e2 * e3 * sgn
                                                 * R[(t + tau) * strideN + (u + ups) * strideU + (v + om)];
                                        }
                                    }
                                }
                                h[((size_t)t * hs + u) * hs + v] = acc;
                            }
                }

            // (a'b'|c'd) for given bra powers and a ket h block.
            auto val = [&](const int a[3], const int b[3], const double* h) {
                double s = 0.0;
                for (int t = 0; t <= a[0] + b[0]; ++t) {
                    const double et = pb.Ex[t * sAB_t + a[0] * sAB_i + b[0]];
                    if (et == 0.0) continue;
                    for (int u = 0; u <= a[1] + b[1]; ++u) {
                        const double eu = pb.Ey[u * sAB_t + a[1] * sAB_i + b[1]];
                        if (eu == 0.0) continue;
                        const double etu = et * eu;
                        const double* hrow = h + ((size_t)t * hs + u) * hs;
                        for (int v = 0; v <= a[2] + b[2]; ++v)
                            s += etu * pb.Ez[v * sAB_t + a[2] * sAB_i + b[2]] * hrow[v];
                    }
                }
                return s;
            };

            for (int cc = 0; cc < (int)C.lmn.size(); ++cc) {
                const auto& lc = C.lmn[cc];
                for (int cd = 0; cd < nD; ++cd) {
                    const double* h0 = hptr(lc[0], lc[1], lc[2], cd);
                    for (int ca = 0; ca < nA; ++ca) {
                        const auto& la = A.lmn[ca];
                        for (int cb = 0; cb < nB; ++cb) {
                            const double w = Dq[((size_t)(ca * nB + cb) * C.lmn.size() + cc) * nD + cd];
                            if (w == 0.0) continue;
                            const auto& lb = B.lmn[cb];
                            const double f = w * pref * A.coef[ca][pb.ia] * B.coef[cb][pb.ib]
                                           * C.coef[cc][pk.ia] * D.coef[cd][pk.ib];
                            int a[3] = { la[0], la[1], la[2] }, b[3] = { lb[0], lb[1], lb[2] };
                            for (int k = 0; k < 3; ++k) {
                                // d/dA: shift the power on A
                                a[k] += 1;
                                double dA = 2.0 * alpha * val(a, b, h0);
                                a[k] -= 2;
                                if (la[k] > 0) dA -= la[k] * val(a, b, h0);
                                a[k] += 1;
                                // d/dB: shift the power on B
                                b[k] += 1;
                                double dB = 2.0 * beta * val(a, b, h0);
                                b[k] -= 2;
                                if (lb[k] > 0) dB -= lb[k] * val(a, b, h0);
                                b[k] += 1;
                                // d/dC: shift the power on C (different ket h block)
                                int cp[3] = { lc[0], lc[1], lc[2] };
                                cp[k] += 1;
                                double dC = 2.0 * gamma * val(a, b, hptr(cp[0], cp[1], cp[2], cd));
                                if (lc[k] > 0) {
                                    cp[k] -= 2;
                                    dC -= lc[k] * val(a, b, hptr(cp[0], cp[1], cp[2], cd));
                                }
                                gA[k] += f * dA;
                                gB[k] += f * dB;
                                gC[k] += f * dC;
                            }
                        }
                    }
                }
            }
        }
    }
}

}  // namespace

// 1/2 sum D (mn|ls)^X over canonical shell quartets only (as in buildERI), each
// weighted by its permutational degeneracy f = (A!=B?2:1)(C!=D?2:1)(AB!=CD?2:1).
// Per quartet the derivatives on the centres of A, B and C are computed and the
// one on D follows from translational invariance, d/dD = -(d/dA + d/dB + d/dC);
// a quartet whose four shells share one atom contributes nothing and is skipped.
// Compared with the Sep-2026 first version (every ORDERED bra pair against the
// canonical kets, derivative on A only) this does ~4x fewer quartets.
Matrix gradientTwoElectron(const std::vector<GTO::Orbital>& basis, const Matrix& P,
                           int natoms, int threads)
{
    Matrix grad = Matrix::Zero(natoms, 3);
    if (basis.empty()) return grad;
    const std::vector<EriShell> sh = groupShells(basis);
    const int ns = (int)sh.size();
    std::vector<int> atomOf(ns);
    for (int s = 0; s < ns; ++s) atomOf[s] = basis[sh[s].first].atom;

    std::vector<ShellPair> pairs;
    pairs.reserve((size_t)ns * (ns + 1) / 2);
    for (int A = 0; A < ns; ++A)
        for (int B = A; B < ns; ++B)
            pairs.push_back(makeShellPair(sh, A, B, 1, 1));
    const int np = (int)pairs.size();

#ifdef _OPENMP
#pragma omp parallel num_threads(threads > 0 ? threads : 1)
#endif
    {
        Matrix g = Matrix::Zero(natoms, 3);
        std::vector<double> Dq;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int i = 0; i < np; ++i) {
            const ShellPair& bra = pairs[i];
            const EriShell& A = sh[bra.A];
            const EriShell& B = sh[bra.B];
            for (int j = i; j < np; ++j) {
                const ShellPair& ket = pairs[j];
                const int at[4] = { atomOf[bra.A], atomOf[bra.B], atomOf[ket.A], atomOf[ket.B] };
                if (at[0] == at[1] && at[0] == at[2] && at[0] == at[3]) continue;
                const EriShell& C = sh[ket.A];
                const EriShell& D = sh[ket.B];
                const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
                Dq.resize((size_t)nA * nB * nC * nD);
                for (int a = 0; a < nA; ++a)
                    for (int b = 0; b < nB; ++b)
                        for (int c = 0; c < nC; ++c)
                            for (int d = 0; d < nD; ++d) {
                                const int m = A.first + a, n = B.first + b, l = C.first + c, s = D.first + d;
                                Dq[((size_t)(a * nB + b) * nC + c) * nD + d] =
                                    P(m, n) * P(l, s) - 0.25 * (P(m, l) * P(n, s) + P(m, s) * P(n, l));
                            }
                const double f = 0.5 * (bra.A != bra.B ? 2.0 : 1.0) * (ket.A != ket.B ? 2.0 : 1.0)
                               * (i != j ? 2.0 : 1.0);
                double gA[3], gB[3], gC[3];
                quartetGradient(sh, bra, ket, Dq, gA, gB, gC);
                for (int k = 0; k < 3; ++k) {
                    g(at[0], k) += f * gA[k];
                    g(at[1], k) += f * gB[k];
                    g(at[2], k) += f * gC[k];
                    g(at[3], k) -= f * (gA[k] + gB[k] + gC[k]);
                }
            }
        }
#ifdef _OPENMP
#pragma omp critical
#endif
        grad += g;
    }
    (void)threads;
    return grad;
}

}  // namespace qmint