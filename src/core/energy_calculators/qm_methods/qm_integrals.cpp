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
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <mutex>
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
    // Rows in parallel (Claude Generated, Oct 2026): every element is computed by one thread exactly as
    // before, so the matrix is bit-identical for any thread count (OpenMP team size of the caller).
#pragma omp parallel for schedule(dynamic)
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
#pragma omp parallel for schedule(dynamic)  // bit-identical, see buildOverlap
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
// Reference Boys function F_n(T), n = 0..maxN, in long double (expl/erfl). Accurate
// but slow (x87); since Sep 2026 it only fills the interpolation table below and
// serves as the reference in test_qm_boys.
void boysFunctionReference(int maxN, double T, std::vector<double>& F)
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

// Fast Boys function (Claude Generated, Sep 2026): the standard grid + Taylor scheme.
// The ERI kernels evaluate F_n(T) once per primitive quartet, millions of times per
// Fock build; the long-double reference was a measurable share of that.
//   T <  36: F_N(T) from the nearest grid point T_i (spacing 0.05, |d| <= 0.025) by
//            F_N(T) = sum_k F_{N+k}(T_i) (T_i - T)^k / k!, k = 0..6 (truncation
//            ~ d^7/7! = 1.2e-15), then the downward recurrence
//            F_n = (2T F_{n+1} + e^{-T}) / (2n+1), which does not amplify errors.
//   T >= 36: F_0 = 0.5 sqrt(pi/T) (erfc(6) = 2e-17) and the upward recurrence
//            F_{n+1} = ((2n+1) F_n - e^{-T}) / (2T), stable for large T.
// Table from boysFunctionReference-quality long-double downward recurrence started
// far above n and T. Accuracy vs the reference: see test_qm_boys.
namespace {
struct BoysTable {
    static constexpr double dT = 0.05;
    static constexpr double Tmax = 36.0;
    static constexpr int NTAB = 40;   // n = 0..NTAB tabulated
    static constexpr int NTAYLOR = 7; // Taylor terms k = 0..6
    static constexpr int NMAX = NTAB - NTAYLOR + 1; // largest maxN served from the table
    int npts = 0;
    std::vector<double> f;  // [i * (NTAB + 1) + n]
    BoysTable()
    {
        npts = (int)std::lround(Tmax / dT) + 2;
        f.assign((size_t)npts * (NTAB + 1), 0.0);
        const int M = NTAB + 200;  // start of the downward recurrence (>> n and T)
        std::vector<long double> G(M + 2, 0.0L);
        for (int i = 0; i < npts; ++i) {
            const long double T = (long double)i * dT;
            const long double eT = expl(-T);
            G[M + 1] = 0.0L;
            for (int n = M; n >= 0; --n) G[n] = (2.0L * T * G[n + 1] + eT) / (2.0L * n + 1.0L);
            for (int n = 0; n <= NTAB; ++n) f[(size_t)i * (NTAB + 1) + n] = (double)G[n];
        }
    }
};
const BoysTable& boysTable()
{
    static const BoysTable table;  // thread-safe one-time init (C++11 magic static)
    return table;
}
}  // namespace

static void boysArrayInto(int maxN, double T, std::vector<double>& F)
{
    if (maxN > BoysTable::NMAX || T < 0.0) {
        boysFunctionReference(maxN, T, F);
        return;
    }
    F.resize(maxN + 1);
    if (T >= BoysTable::Tmax) {
        const double eT = std::exp(-T);
        F[0] = 0.5 * std::sqrt(PI / T);
        for (int n = 0; n < maxN; ++n) F[n + 1] = ((2.0 * n + 1.0) * F[n] - eT) / (2.0 * T);
        return;
    }
    const BoysTable& tab = boysTable();
    const int i = (int)(T / BoysTable::dT + 0.5);
    const double d = i * BoysTable::dT - T;  // T = T_i - d
    const double* fi = &tab.f[(size_t)i * (BoysTable::NTAB + 1) + maxN];
    // sum_k F_{maxN+k}(T_i) d^k / k!  (Horner)
    constexpr double inv[BoysTable::NTAYLOR] = { 1.0, 1.0, 1.0 / 2, 1.0 / 3, 1.0 / 4, 1.0 / 5, 1.0 / 6 };
    double acc = fi[BoysTable::NTAYLOR - 1];
    for (int k = BoysTable::NTAYLOR - 2; k >= 0; --k) acc = fi[k] + d * inv[k + 1] * acc;
    F[maxN] = acc;
    if (maxN > 0) {
        const double eT = std::exp(-T), T2 = 2.0 * T;
        for (int n = maxN - 1; n >= 0; --n) F[n] = (T2 * F[n + 1] + eT) / (2.0 * n + 1.0);
    }
}

void boysFunction(int maxN, double T, std::vector<double>& F) { boysArrayInto(maxN, T, F); }

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

// V_ij summed over all centres, in the order of the former loop (centres outer, primitive pairs
// inner, -Z * contracted sum per centre), so the result is bit-identical to
//   sum_C -Z_C sum_prim c_a c_b primitiveNuclearAttraction(..., C)  (the former per-centre loop).
// The Gaussian product and the three Hermite tables of a primitive pair do not depend on the centre;
// they are built once per primitive pair here instead of once per (primitive pair, centre), and R/F
// reuse thread-local buffers (Claude Generated, Oct 2026: the former per-call allocations dominated the
// O(n^2 natoms) build, 19 s single-threaded at 970 functions / 122 atoms).
static double nuclearPairAllCentres(const GTO::Orbital& a, const GTO::Orbital& b, const std::vector<int>& atomZ,
                                    const Matrix& pos)
{
    int la, ma, na, lb, mb, nb;
    GTO::orbitalTypeToComponents(a.type, la, ma, na);
    GTO::orbitalTypeToComponents(b.type, lb, mb, nb);
    const int tmax = la + lb, umax = ma + mb, vmax = na + nb, Nmax = tmax + umax + vmax;
    struct PrimData {
        double coef, gamma, Px, Py, Pz;
        std::vector<double> et, eu, ev;  // E^{la lb}_t, E^{ma mb}_u, E^{na nb}_v
    };
    thread_local std::vector<PrimData> prims;
    thread_local std::vector<double> R, F;
    prims.clear();
    for (size_t ia = 0; ia < a.exponents.size(); ++ia)
        for (size_t ib = 0; ib < b.exponents.size(); ++ib) {
            PrimData pd;
            pd.coef = a.coefficients[ia] * b.coefficients[ib];
            const double K = gaussianProductK(a.exponents[ia], b.exponents[ib], a.x, a.y, a.z, b.x, b.y, b.z, pd.gamma,
                                              pd.Px, pd.Py, pd.Pz);
            if (K != 0.0) {
                const auto Ex = hermiteCoeffs(la, lb, pd.Px - a.x, pd.Px - b.x, pd.gamma, K);
                const auto Eu = hermiteCoeffs(ma, mb, pd.Py - a.y, pd.Py - b.y, pd.gamma, 1.0);
                const auto Ev = hermiteCoeffs(na, nb, pd.Pz - a.z, pd.Pz - b.z, pd.gamma, 1.0);
                for (int t = 0; t <= tmax; ++t) pd.et.push_back(Ex[t][la][lb]);
                for (int u = 0; u <= umax; ++u) pd.eu.push_back(Eu[u][ma][mb]);
                for (int v = 0; v <= vmax; ++v) pd.ev.push_back(Ev[v][na][nb]);
            }
            prims.push_back(std::move(pd));  // K == 0: empty tables, contributes coef * 0.0
        }
    R.resize((size_t)(tmax + 1) * (umax + 1) * (vmax + 1) * (Nmax + 1));
    const int strideN = (umax + 1) * (vmax + 1) * (Nmax + 1);
    const int strideU = (vmax + 1) * (Nmax + 1);
    const int strideV = (Nmax + 1);
    double total = 0.0;
    for (int C = 0; C < (int)atomZ.size(); ++C) {
        const double Z = atomZ[C];
        if (Z == 0.0) continue;
        const double Cx = pos(C, 0), Cy = pos(C, 1), Cz = pos(C, 2);
        double sum = 0.0;
        for (const PrimData& pd : prims) {
            double val = 0.0;
            if (!pd.et.empty()) {
                const double T = pd.gamma * ((pd.Px - Cx) * (pd.Px - Cx) + (pd.Py - Cy) * (pd.Py - Cy) + (pd.Pz - Cz) * (pd.Pz - Cz));
                boysArrayInto(Nmax, T, F);
                buildRblock(tmax, umax, vmax, Nmax, pd.Px - Cx, pd.Py - Cy, pd.Pz - Cz, pd.gamma, F, R);
                for (int t = 0; t <= tmax; ++t) {
                    const double et = pd.et[t];
                    if (et == 0.0) continue;
                    for (int u = 0; u <= umax; ++u) {
                        const double eu = pd.eu[u];
                        if (eu == 0.0) continue;
                        for (int v = 0; v <= vmax; ++v) {
                            const double ev = pd.ev[v];
                            if (ev == 0.0) continue;
                            val += et * eu * ev * R[t * strideN + u * strideU + v * strideV + 0];
                        }
                    }
                }
            }
            sum += pd.coef * val;
        }
        total += -Z * sum;
    }
    return total;
}

Matrix buildNuclearAttraction(const std::vector<GTO::Orbital>& basis,
                              const std::vector<int>& atomZ,
                              const Matrix& atomPosBohr)
{
    const int n = (int)basis.size();
    Matrix V = Matrix::Zero(n, n);
    // Rows in parallel, centres C innermost (Claude Generated, Oct 2026; was serial with C outermost and
    // ~3.5 s at 560 cartesian functions). Each element still sums its centres in ascending order starting
    // from 0, so V is bit-identical to the serial loop for any thread count.
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            const double sum = nuclearPairAllCentres(basis[i], basis[j], atomZ, atomPosBohr);
            V(i, j) = sum;
            V(j, i) = sum;
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

/// Total angular momentum of one cartesian AO.
static int orbitalL(const GTO::Orbital& o)
{
    int l, m, n;
    GTO::orbitalTypeToComponents(o.type, l, m, n);
    return l + m + n;
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

    // Claude Generated (Nov 2026): f and g shells (l >= 3, (l+1)(l+2)/2 contiguous components in lexicographic
    // order) collapse to 2l+1 real solid harmonics as well. The d branch below is the Sep 2026 code unchanged.
    bool anyD = false;
    for (const auto& o : basis) if (isDShell(o) || orbitalL(o) >= 3) { anyD = true; break; }
    if (!anyD) return Matrix();  // empty -> caller keeps cartesian

    // nbf_sph = nbf_cart - sum over shells of (ncart - (2l+1))  (6d -> 5d, 10f -> 7f, 15g -> 9g)
    int nsph = 0;
    for (int i = 0; i < n; ) {
        const int L = orbitalL(basis[i]);
        if (isDShell(basis[i])) {
            nsph += 5;
            i += 6;  // assume contiguous block of 6
        } else if (L >= 3) {
            nsph += 2 * L + 1;
            i += (L + 1) * (L + 2) / 2;
        } else {
            ++nsph;
            ++i;
        }
    }
    Matrix Q = Matrix::Zero(n, nsph);

    int sphCol = 0;
    int i = 0;
    while (i < n) {
        const int Lsh = orbitalL(basis[i]);
        if (Lsh >= 3) {
            // f/g shell (Claude Generated, Nov 2026): 2l+1 real solid harmonics, qm_spherical_fg.cpp
            const int nc = (Lsh + 1) * (Lsh + 2) / 2;
            higherShellSphericalColumns(basis, Scart, i, Lsh, Q, sphCol);
            i += nc;
            sphCol += 2 * Lsh + 1;
            continue;
        }
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
    // No zero fill: every read below stays inside the previous level's simplex, and
    // the callers only read R0 on the simplex t+u+v <= L.
    R0.resize(sz);
    lvl.resize(sz);
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

// Raw-pointer twins of boysArrayInto and buildR0Simplex for the Hermite-matrix ERI kernel
// (Claude Generated, Oct 2026): the same operations in the same order (bitwise the same
// values), without the std::vector bookkeeping that cost ~30 % of a low-L quartet.
// F must hold maxN+1 values; maxN above the table falls back to the reference.
static void boysRaw(int maxN, double T, double* F)
{
    if (maxN > BoysTable::NMAX || T < 0.0) {
        std::vector<double> Fv;
        boysFunctionReference(maxN, T, Fv);
        for (int n = 0; n <= maxN; ++n) F[n] = Fv[n];
        return;
    }
    if (T >= BoysTable::Tmax) {
        F[0] = 0.5 * std::sqrt(PI / T);
        if (maxN > 0) {
            const double eT = std::exp(-T);
            for (int n = 0; n < maxN; ++n) F[n + 1] = ((2.0 * n + 1.0) * F[n] - eT) / (2.0 * T);
        }
        return;
    }
    const BoysTable& tab = boysTable();
    const int i = (int)(T / BoysTable::dT + 0.5);
    const double d = i * BoysTable::dT - T;
    const double* fi = &tab.f[(size_t)i * (BoysTable::NTAB + 1) + maxN];
    constexpr double inv[BoysTable::NTAYLOR] = { 1.0, 1.0, 1.0 / 2, 1.0 / 3, 1.0 / 4, 1.0 / 5, 1.0 / 6 };
    double acc = fi[BoysTable::NTAYLOR - 1];
    for (int k = BoysTable::NTAYLOR - 2; k >= 0; --k) acc = fi[k] + d * inv[k + 1] * acc;
    F[maxN] = acc;
    if (maxN > 0) {
        const double eT = std::exp(-T), T2 = 2.0 * T;
        for (int n = maxN - 1; n >= 0; --n) F[n] = (T2 * F[n + 1] + eT) / (2.0 * n + 1.0);
    }
}

/// boysRaw with the divisions replaced by multiplications with reciprocals (Claude
/// Generated, Oct 2026): 1/dT, 1/(2n+1) from a table, one 1/(2T); equal to boysRaw up to
/// rounding (a few ulp). Used by the templated ERI kernel, where the Boys function of
/// the low classes (ss|ss), (ss|sp) is a large share of the time.
static void boysFast(int maxN, double T, double* F)
{
    static constexpr double inv_odd[18] = { 1.0,        1.0 / 3,  1.0 / 5,  1.0 / 7,  1.0 / 9,  1.0 / 11,
                                            1.0 / 13,  1.0 / 15, 1.0 / 17, 1.0 / 19, 1.0 / 21, 1.0 / 23,
                                            1.0 / 25,  1.0 / 27, 1.0 / 29, 1.0 / 31, 1.0 / 33, 1.0 / 35 };
    if (maxN > 17 || maxN > BoysTable::NMAX || T < 0.0) {
        boysRaw(maxN, T, F);
        return;
    }
    if (T >= BoysTable::Tmax) {
        const double iT = 1.0 / T;
        F[0] = 0.5 * std::sqrt(PI * iT);
        if (maxN > 0) {
            const double eT = std::exp(-T), h = 0.5 * iT;
            for (int n = 0; n < maxN; ++n) F[n + 1] = ((2.0 * n + 1.0) * F[n] - eT) * h;
        }
        return;
    }
    const BoysTable& tab = boysTable();
    const int i = (int)(T * (1.0 / BoysTable::dT) + 0.5);
    const double d = i * BoysTable::dT - T;
    const double* fi = &tab.f[(size_t)i * (BoysTable::NTAB + 1) + maxN];
    constexpr double inv[BoysTable::NTAYLOR] = { 1.0, 1.0, 1.0 / 2, 1.0 / 3, 1.0 / 4, 1.0 / 5, 1.0 / 6 };
    double acc = fi[BoysTable::NTAYLOR - 1];
    for (int k = BoysTable::NTAYLOR - 2; k >= 0; --k) acc = fi[k] + d * inv[k + 1] * acc;
    F[maxN] = acc;
    if (maxN > 0) {
        const double eT = std::exp(-T), T2 = 2.0 * T;
        for (int n = maxN - 1; n >= 0; --n) F[n] = (T2 * F[n + 1] + eT) * inv_odd[n];
    }
}

/// R^0_tuv on the simplex t+u+v <= L into R0 (box layout (t*(L+1)+u)*(L+1)+v), scratch lvl
/// of the same size; recursion of buildR0Simplex (z before y before x), level by level.
static void buildR0Raw(int L, double Wx, double Wy, double Wz, double rho, const double* F, double* R0,
                              double* lvl)
{
    const int s = L + 1;
    const double b = -2.0 * rho;
    double bpow[64];
    bpow[0] = 1.0;
    for (int n = 1; n <= L; ++n) bpow[n] = bpow[n - 1] * b;
    // level n is written into R0 when (L - n) is even: level 0 ends in R0
    double* c = (L % 2 == 0) ? R0 : lvl;  // level L
    double* x = (L % 2 == 0) ? lvl : R0;
    c[0] = bpow[L] * F[L];
    for (int n = L - 1; n >= 0; --n) {
        const int top = L - n;
        x[0] = bpow[n] * F[n];
        for (int t = 0; t <= top; ++t)
            for (int u = 0; u <= top - t; ++u) {
                const int base = (t * s + u) * s;
                const int vmax = top - t - u;
                if (t + u > 0) {
                    if (u > 0)
                        x[base] = (u >= 2 ? (u - 1) * c[base - 2 * s] : 0.0) + Wy * c[base - s];
                    else
                        x[base] = (t >= 2 ? (t - 1) * c[base - 2 * s * s] : 0.0) + Wx * c[base - s * s];
                }
                if (vmax >= 1) x[base + 1] = 0.0 + Wz * c[base];
                for (int v = 2; v <= vmax; ++v) x[base + v] = (v - 1) * c[base + v - 2] + Wz * c[base + v - 1];
            }
        std::swap(c, x);
    }
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
    double schwarz = 0.0;  ///< sqrt(max |(ab|ab)|) of this primitive pair alone, coefficients included
    // Special members out of line: the implicit inline ones trip -Winline (sorting moves them).
    PrimPair();
    PrimPair(const PrimPair&);
    PrimPair(PrimPair&&) noexcept;
    PrimPair& operator=(const PrimPair&);
    PrimPair& operator=(PrimPair&&) noexcept;
    ~PrimPair();
};
PrimPair::PrimPair() = default;
PrimPair::PrimPair(const PrimPair&) = default;
PrimPair::PrimPair(PrimPair&&) noexcept = default;
PrimPair& PrimPair::operator=(const PrimPair&) = default;
PrimPair& PrimPair::operator=(PrimPair&&) noexcept = default;
PrimPair::~PrimPair() = default;

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
/// out[((ca*ncB + cb)*ncC + cc)*ncD + cd].
///
/// McMurchie-Davidson (Helgaker 9.9.33) in three steps (Claude Generated, Sep 2026;
/// the first blocked version did the six Hermite sums of every component quartet in
/// one nested loop, which was 78 % of the ERI time):
///
///  1. per primitive quartet: Boys function and R^0_tuv on the simplex t+u+v <= L;
///  2. ket side, per ket component pair (c,d): the ket Hermite sum
///       h_cd(t,u,v) = sum_{tau,nu,phi} (-1)^(tau+nu+phi) E^cd_tau E^cd_nu E^cd_phi R_{t+tau,u+nu,v+phi}
///     done one Cartesian direction at a time (z, then y, then x -- three short
///     sums instead of one triple sum), weighted with the ket contraction
///     coefficients and the prefactor and ACCUMULATED over the ket primitive pairs
///     into G_cd(t,u,v) -- the ket primitives are contracted before the bra is touched;
///  3. bra side, once per bra primitive pair: (ab|cd) += c_a c_b sum_tuv E^ab_t E^ab_u E^ab_v G_cd(t,u,v).
///
/// Same arithmetic as before, only regrouped, so the values agree to rounding.
void shellQuartet(const std::vector<EriShell>& sh, const ShellPair& bra, const ShellPair& ket,
                  std::vector<double>& out, std::vector<double>& R, double prim_thr = 0.0)
{
    const EriShell& A = sh[bra.A];
    const EriShell& B = sh[bra.B];
    const EriShell& C = sh[ket.A];
    const EriShell& D = sh[ket.B];
    const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
    const int nCD = nC * nD;
    out.assign((size_t)nA * nB * nC * nD, 0.0);

    const int LAB = A.L + B.L, LCD = C.L + D.L;
    const int Ltot = LAB + LCD;
    // Hermite table strides: E[t][i][j] with i <= LA, j <= LB.
    const int sAB_t = (A.L + 1) * (B.L + 1), sAB_i = B.L + 1;
    const int sCD_t = (C.L + 1) * (D.L + 1), sCD_i = D.L + 1;
    // R^0 (and the intermediates) live in a (Ltot+1)^3 box, index (t*s + u)*s + v.
    const int s = Ltot + 1;
    const size_t box = (size_t)s * s * s;
    auto at = [s](int t, int u, int v) { return ((size_t)t * s + u) * s + v; };
    // G_cd(t,u,v) for t+u+v <= LAB, in a (LAB+1)^3 box per ket component pair.
    const int g = LAB + 1;
    const size_t gsz = (size_t)g * g * g;

    thread_local std::vector<double> F, G, A1, A2;
    G.resize((size_t)nCD * gsz);
    A1.resize(box);
    A2.resize(box);

    // Primitive screening (Sep 2026): the primitive pairs are sorted by their own
    // Schwarz factor (makeScreeningPairs), and |(ab|cd)| of one primitive quartet is
    // bounded by pb.schwarz * pk.schwarz (Cauchy-Schwarz per primitive product), so
    // both loops stop at the first quartet below prim_thr. prim_thr = 0: every quartet.
    const double ket_max = ket.prims.empty() ? 0.0 : ket.prims.front().schwarz;
    for (const PrimPair& pb : bra.prims) {
        if (pb.schwarz * ket_max < prim_thr) break;
        std::fill(G.begin(), G.end(), 0.0);

        for (const PrimPair& pk : ket.prims) {
            if (pb.schwarz * pk.schwarz < prim_thr) break;
            const double p = pb.p, q = pk.p;
            const double rho = p * q / (p + q);
            const double Wx = pb.Px - pk.Px, Wy = pb.Py - pk.Py, Wz = pb.Pz - pk.Pz;
            boysArrayInto(Ltot, rho * (Wx * Wx + Wy * Wy + Wz * Wz), F);
            buildR0Simplex(Ltot, Wx, Wy, Wz, rho, F, R);
            const double pref = 2.0 * std::pow(PI, 2.5) / (p * q * std::sqrt(p + q));

            for (int cc = 0; cc < nC; ++cc) {
                const auto& lc = C.lmn[cc];
                const double c_c = pref * C.coef[cc][pk.ia];
                for (int cd = 0; cd < nD; ++cd) {
                    const double w = c_c * D.coef[cd][pk.ib];
                    if (w == 0.0) continue;
                    const auto& ld = D.lmn[cd];
                    const int tx = lc[0] + ld[0], ty = lc[1] + ld[1], tz = lc[2] + ld[2];
                    // signed ket Hermite coefficients (-1)^k E_k per direction
                    double ex[16], ey[16], ez[16];
                    for (int k = 0; k <= tx; ++k) ex[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ex[k * sCD_t + lc[0] * sCD_i + ld[0]];
                    for (int k = 0; k <= ty; ++k) ey[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ey[k * sCD_t + lc[1] * sCD_i + ld[1]];
                    for (int k = 0; k <= tz; ++k) ez[k] = ((k & 1) ? -1.0 : 1.0) * pk.Ez[k * sCD_t + lc[2] * sCD_i + ld[2]];

                    double* Gcd = &G[(size_t)(cc * nD + cd) * gsz];
                    if (tx + ty + tz == 0) {
                        // s-type ket pair: the three sums collapse to one scaled R.
                        const double wk = w * ex[0] * ey[0] * ez[0];
                        for (int t = 0; t <= LAB; ++t)
                            for (int u = 0; u <= LAB - t; ++u) {
                                const double* Rr = &R[at(t, u, 0)];
                                double* gr = Gcd + ((size_t)t * g + u) * g;
                                for (int v = 0; v <= LAB - t - u; ++v) gr[v] += wk * Rr[v];
                            }
                        continue;
                    }
                    // Only the index ranges that lead back to t+u+v <= LAB are built:
                    // x-index t' <= LAB+tx, y-index u' with u'+v <= LAB+ty, z-index v <= LAB.
                    // z: A1(t,u,v) = sum_phi ez_phi R(t,u,v+phi)
                    const int m1 = LAB + tx + ty, m2 = LAB + tx;
                    for (int t = 0; t <= m2; ++t)
                        for (int u = 0; u <= std::min(m1 - t, LAB + ty); ++u) {
                            const double* Rr = &R[at(t, u, 0)];
                            double* a1 = &A1[at(t, u, 0)];
                            const int vmax = std::min(m1 - t - u, LAB);
                            for (int v = 0; v <= vmax; ++v) {
                                double acc = 0.0;
                                for (int k = 0; k <= tz; ++k) acc += ez[k] * Rr[v + k];
                                a1[v] = acc;
                            }
                        }
                    // y: A2(t,u,v) = sum_nu ey_nu A1(t,u+nu,v),   u+v <= LAB
                    for (int t = 0; t <= m2; ++t)
                        for (int u = 0; u <= std::min(m2 - t, LAB); ++u) {
                            const int vmax = std::min(m2 - t - u, LAB - u);
                            for (int v = 0; v <= vmax; ++v) {
                                double acc = 0.0;
                                for (int k = 0; k <= ty; ++k) acc += ey[k] * A1[at(t, u + k, v)];
                                A2[at(t, u, v)] = acc;
                            }
                        }
                    // x: G_cd(t,u,v) += w * sum_tau ex_tau A2(t+tau,u,v),   t+u+v <= LAB
                    for (int t = 0; t <= LAB; ++t)
                        for (int u = 0; u <= LAB - t; ++u)
                            for (int v = 0; v <= LAB - t - u; ++v) {
                                double acc = 0.0;
                                for (int k = 0; k <= tx; ++k) acc += ex[k] * A2[at(t + k, u, v)];
                                Gcd[((size_t)t * g + u) * g + v] += w * acc;
                            }
                }
            }
        }

        // Bra contraction, once per bra primitive pair.
        for (int ca = 0; ca < nA; ++ca) {
            const auto& la = A.lmn[ca];
            const double c_a = A.coef[ca][pb.ia];
            for (int cb = 0; cb < nB; ++cb) {
                const double c_ab = c_a * B.coef[cb][pb.ib];
                if (c_ab == 0.0) continue;
                const auto& lb = B.lmn[cb];
                const int tx = la[0] + lb[0], ty = la[1] + lb[1], tz = la[2] + lb[2];
                double ex[16], ey[16], ez[16];
                for (int k = 0; k <= tx; ++k) ex[k] = pb.Ex[k * sAB_t + la[0] * sAB_i + lb[0]];
                for (int k = 0; k <= ty; ++k) ey[k] = pb.Ey[k * sAB_t + la[1] * sAB_i + lb[1]];
                for (int k = 0; k <= tz; ++k) ez[k] = pb.Ez[k * sAB_t + la[2] * sAB_i + lb[2]];
                double* o = &out[(size_t)(ca * nB + cb) * nCD];
                for (int k = 0; k < nCD; ++k) {
                    const double* Gcd = &G[(size_t)k * gsz];
                    double val = 0.0;
                    for (int t = 0; t <= tx; ++t) {
                        if (ex[t] == 0.0) continue;
                        for (int u = 0; u <= ty; ++u) {
                            const double etu = ex[t] * ey[u];
                            if (etu == 0.0) continue;
                            const double* row = Gcd + ((size_t)t * g + u) * g;
                            double sv = 0.0;
                            for (int v = 0; v <= tz; ++v) sv += ez[v] * row[v];
                            val += etu * sv;
                        }
                    }
                    o[k] += c_ab * val;
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
                   std::vector<double>& tmp, int d[4], double prim_thr)
{
    // shellQuartet does its expensive ket-side work once per PRIMITIVE quartet and
    // the bra side once per bra primitive, so the pair with more components goes
    // on the bra side: (ab|cd) = (cd|ab), computed swapped and transposed back.
    const size_t nab = sh[bra.A].lmn.size() * sh[bra.B].lmn.size();
    const size_t ncd = sh[ket.A].lmn.size() * sh[ket.B].lmn.size();
    if (ncd > nab) {
        shellQuartet(sh, ket, bra, tmp, R, prim_thr);
        buf.resize(nab * ncd);
        for (size_t x = 0; x < nab; ++x)
            for (size_t y = 0; y < ncd; ++y) buf[x * ncd + y] = tmp[y * nab + x];
    } else {
        shellQuartet(sh, bra, ket, buf, R, prim_thr);
    }
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
/// sqrt(max |(ab|ab)|) from the diagonal quartets -- per shell pair and, for the
/// primitive screening, per primitive pair (then sorted, largest first).
std::vector<ShellPair> makeScreeningPairs(const std::vector<EriShell>& sh, int threads = 1)
{
    const int ns = (int)sh.size();
    std::vector<ShellPair> pairs;
    pairs.reserve((size_t)ns * (ns + 1) / 2);
    for (int A = 0; A < ns; ++A)
        for (int B = A; B < ns; ++B) {
            pairs.emplace_back();
            pairs.back().A = A;
            pairs.back().B = B;
        }
    // Every pair is independent (its own primitive tables and diagonal quartets), so the
    // loop threads without changing any value (Claude Generated, Oct 2026: OpenMP).
    const int np = (int)pairs.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads > 0 ? threads : 1)
#endif
    for (int ip = 0; ip < np; ++ip) {
        std::vector<double> buf, R;
        ShellPair& sp = pairs[ip];
        sp = makeShellPair(sh, sp.A, sp.B);
        shellQuartet(sh, sp, sp, buf, R);
        const int nA = (int)sh[sp.A].lmn.size(), nB = (int)sh[sp.B].lmn.size();
        double mx = 0.0;
        for (int a = 0; a < nA; ++a)
            for (int b = 0; b < nB; ++b)
                mx = std::max(mx, std::abs(buf[((size_t)(a * nB + b) * nA + a) * nB + b]));
        sp.schwarz = std::sqrt(mx);

        ShellPair one;
        one.A = sp.A;
        one.B = sp.B;
        for (PrimPair& pp : sp.prims) {
            one.prims.assign(1, pp);
            shellQuartet(sh, one, one, buf, R);
            double m1 = 0.0;
            for (int a = 0; a < nA; ++a)
                for (int b = 0; b < nB; ++b)
                    m1 = std::max(m1, std::abs(buf[((size_t)(a * nB + b) * nA + a) * nB + b]));
            pp.schwarz = std::sqrt(m1);
        }
        std::stable_sort(sp.prims.begin(), sp.prims.end(),
                         [](const PrimPair& x, const PrimPair& y) { return x.schwarz > y.schwarz; });
    }
    (void)threads;
    return pairs;
}

// ---------------------------------------------------------------------------
// Hermite-matrix form of the McMurchie-Davidson quartet (Claude Generated, Oct 2026)
//
// shellQuartet above contracts the Hermite tables one cartesian direction at a time
// and transforms d shells to the spherical basis afterwards (transformAxis). That is
// the readable reference, but its per-primitive work is short loops with variable
// bounds, and it was 90 % of an integral-direct SCF. The same equation (Helgaker
// 9.9.33)
//
//   (ab|cd) = sum_{h in simplex(LAB)} sum_{k in simplex(LCD)} E^ab_h (-1)^|k| E^cd_k R_{h+k}
//
// is evaluated here as two small matrix products over precomputed per-primitive-pair
// matrices, the same idea as the device MD kernels (cuda/qm_gpu_direct_md_impl.cuh):
//
//   E^ab[ab'][h] = sum_{a,b} T_A[a][a'] T_B[b][b'] c_a c_b E^x_t E^y_u E^z_v ,  h = (t,u,v)
//
// with the contraction coefficients and the cartesian -> spherical transform T folded
// in, built ONCE per primitive pair (makeMDPair). Per primitive quartet the Boys
// values carry the prefactor 2 pi^(5/2) / (p q sqrt(p+q)), so R is already scaled, and
//
//   ket step  G[cd'][h] += sum_k (-1)^|k| E^cd[cd'][k] R_{h+k}     (per primitive quartet)
//   bra step  out[ab'][cd'] += sum_h E^ab[ab'][h] G[cd'][h]       (per bra primitive pair)
//
// R lives in the (L+1)^3 box of buildR0Simplex, whose index is linear in (t,u,v), so
// R_{h+k} = R[off(h) + off(k)]: one gather per k, then a vectorisable axpy over h.
// Same primitive screening (sorted primitive Schwarz factors, prim_thr) and the same
// sums as shellQuartet + transformAxis, only regrouped: values agree to rounding.
// ---------------------------------------------------------------------------

/// Hermite index triples (t,u,v) with t+u+v <= L, in loop order t, u, v (cached).
const std::vector<std::array<int, 3>>& hermiteSimplex(int L)
{
    static const std::vector<std::vector<std::array<int, 3>>> table = [] {
        std::vector<std::vector<std::array<int, 3>>> tab(17);
        for (int l = 0; l <= 16; ++l)
            for (int t = 0; t <= l; ++t)
                for (int u = 0; u <= l - t; ++u)
                    for (int v = 0; v <= l - t - u; ++v) tab[l].push_back({ t, u, v });
        return tab;
    }();
    return table[L];
}

/// One shell pair in the Hermite-matrix form: per primitive pair (same order as
/// ShellPair::prims, i.e. sorted by the primitive Schwarz factor) the composite
/// exponent, centre, Schwarz factor and the matrix E[ab'][h] (active basis).
struct MDPair {
    int LAB = 0;    ///< LA + LB
    int nh = 0;     ///< Hermite simplex size (LAB+1)(LAB+2)(LAB+3)/6
    int nout = 0;   ///< output components wA * wB in the active basis
    bool tabled = false;
    std::vector<double> p, ip, Px, Py, Pz, schwarz;  ///< ip = 1/p
    std::vector<double> E;  ///< [prim][ab'][h]
    ~MDPair();  // out of line: the implicit inline one trips -Winline
};
MDPair::~MDPair() = default;

MDPair makeMDPair(const std::vector<EriShell>& sh, const ShellOutput& o, const ShellPair& sp)
{
    MDPair m;
    const EriShell& A = sh[sp.A];
    const EriShell& B = sh[sp.B];
    const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size();
    const int wA = o.width[sp.A], wB = o.width[sp.B];
    m.LAB = A.L + B.L;
    const auto& hs = hermiteSimplex(m.LAB);
    m.nh = (int)hs.size();
    m.nout = wA * wB;
    m.tabled = true;
    const int sT = (A.L + 1) * (B.L + 1), sI = B.L + 1;  // Hermite table strides [t][i][j]
    const bool idA = !o.spherical || o.identity[sp.A];
    const bool idB = !o.spherical || o.identity[sp.B];
    std::vector<double> Ec((size_t)nA * nB * m.nh), Eb((size_t)wA * nB * m.nh);
    m.E.reserve(sp.prims.size() * (size_t)m.nout * m.nh);
    for (const PrimPair& pp : sp.prims) {
        m.p.push_back(pp.p);
        m.ip.push_back(1.0 / pp.p);
        m.Px.push_back(pp.Px);
        m.Py.push_back(pp.Py);
        m.Pz.push_back(pp.Pz);
        m.schwarz.push_back(pp.schwarz);
        // cartesian E[a][b][h] with the contraction coefficients
        for (int a = 0; a < nA; ++a) {
            const auto& la = A.lmn[a];
            for (int b = 0; b < nB; ++b) {
                const auto& lb = B.lmn[b];
                const double cab = A.coef[a][pp.ia] * B.coef[b][pp.ib];
                double* e = &Ec[((size_t)a * nB + b) * m.nh];
                for (int h = 0; h < m.nh; ++h) {
                    const int t = hs[h][0], u = hs[h][1], v = hs[h][2];
                    e[h] = (t > la[0] + lb[0] || u > la[1] + lb[1] || v > la[2] + lb[2])
                        ? 0.0
                        : cab * pp.Ex[t * sT + la[0] * sI + lb[0]] * pp.Ey[u * sT + la[1] * sI + lb[1]]
                            * pp.Ez[v * sT + la[2] * sI + lb[2]];
                }
            }
        }
        // transform index a, then b (identity for s/p and cartesian output)
        for (int a2 = 0; a2 < wA; ++a2)
            for (int b = 0; b < nB; ++b) {
                double* y = &Eb[((size_t)a2 * nB + b) * m.nh];
                for (int h = 0; h < m.nh; ++h) y[h] = 0.0;
                for (int a = 0; a < nA; ++a) {
                    const double t = idA ? (a == a2 ? 1.0 : 0.0) : o.T[sp.A](a, a2);
                    if (t == 0.0) continue;
                    const double* x = &Ec[((size_t)a * nB + b) * m.nh];
                    for (int h = 0; h < m.nh; ++h) y[h] += t * x[h];
                }
            }
        const size_t base = m.E.size();
        m.E.resize(base + (size_t)m.nout * m.nh, 0.0);
        for (int a2 = 0; a2 < wA; ++a2)
            for (int b2 = 0; b2 < wB; ++b2) {
                double* y = &m.E[base + ((size_t)a2 * wB + b2) * m.nh];
                for (int b = 0; b < nB; ++b) {
                    const double t = idB ? (b == b2 ? 1.0 : 0.0) : o.T[sp.B](b, b2);
                    if (t == 0.0) continue;
                    const double* x = &Eb[((size_t)a2 * nB + b) * m.nh];
                    for (int h = 0; h < m.nh; ++h) y[h] += t * x[h];
                }
            }
    }
    return m;
}

/// Per-thread scratch of mdQuartet (no allocation in the quartet loop).
struct MDScratch {
    std::vector<double> F, R, tmpR, G, r, tmp;
    std::vector<int> offa, offc;
    std::vector<double> sgnc;
    ~MDScratch();
};
MDScratch::~MDScratch() = default;

/// Contracted (ab|cd) of one shell quartet in the active basis, written to
/// out[ab' * ket.nout + cd'] (ab' = a' wB + b'). See the block comment above.
void mdQuartet(const MDPair& bra, const MDPair& ket, double prim_thr, double* out, MDScratch& w)
{
    const int Lt = bra.LAB + ket.LAB;
    const int s = Lt + 1;
    const int nha = bra.nh, nhc = ket.nh, nab = bra.nout, ncd = ket.nout;
    const auto& ha = hermiteSimplex(bra.LAB);
    const auto& hc = hermiteSimplex(ket.LAB);
    w.offa.resize(nha);
    w.offc.resize(nhc);
    w.sgnc.resize(nhc);
    for (int h = 0; h < nha; ++h) w.offa[h] = (ha[h][0] * s + ha[h][1]) * s + ha[h][2];
    for (int k = 0; k < nhc; ++k) {
        w.offc[k] = (hc[k][0] * s + hc[k][1]) * s + hc[k][2];
        w.sgnc[k] = ((hc[k][0] + hc[k][1] + hc[k][2]) & 1) ? -1.0 : 1.0;  // (-1)^(tau+nu+phi)
    }
    w.G.resize((size_t)ncd * nha);
    w.r.resize(nha);
    w.F.resize(Lt + 1);
    w.R.resize((size_t)s * s * s);
    w.tmpR.resize((size_t)s * s * s);
    for (int x = 0; x < nab * ncd; ++x) out[x] = 0.0;
    const int* offa = w.offa.data();
    double* G = w.G.data();
    double* r = w.r.data();
    static const double pref0 = 2.0 * std::pow(PI, 2.5);

    const int npb = (int)bra.p.size(), npk = (int)ket.p.size();
    const double ket_max = npk ? ket.schwarz[0] : 0.0;
    for (int ib = 0; ib < npb; ++ib) {
        if (bra.schwarz[ib] * ket_max < prim_thr) break;
        std::fill(w.G.begin(), w.G.end(), 0.0);
        const double p = bra.p[ib];
        for (int ik = 0; ik < npk; ++ik) {
            if (bra.schwarz[ib] * ket.schwarz[ik] < prim_thr) break;
            const double q = ket.p[ik];
            const double rho = p * q / (p + q);
            const double Wx = bra.Px[ib] - ket.Px[ik], Wy = bra.Py[ib] - ket.Py[ik], Wz = bra.Pz[ib] - ket.Pz[ik];
            double* F = w.F.data();
            boysRaw(Lt, rho * (Wx * Wx + Wy * Wy + Wz * Wz), F);
            const double pref = pref0 / (p * q * std::sqrt(p + q));
            for (int n = 0; n <= Lt; ++n) F[n] *= pref;  // R is linear in F: R carries the prefactor
            buildR0Raw(Lt, Wx, Wy, Wz, rho, F, w.R.data(), w.tmpR.data());
            const double* R = w.R.data();
            const double* Ek = &ket.E[(size_t)ik * ncd * nhc];
            // ket step: G[cd][h] += sum_k (-1)^|k| E^cd[cd][k] R[off(h) + off(k)]
            for (int k = 0; k < nhc; ++k) {
                const double* Rk = R + w.offc[k];
                const double sg = w.sgnc[k];
                for (int h = 0; h < nha; ++h) r[h] = sg * Rk[offa[h]];
                for (int cd = 0; cd < ncd; ++cd) {
                    const double e = Ek[(size_t)cd * nhc + k];
                    if (e == 0.0) continue;
                    double* g = G + (size_t)cd * nha;
                    for (int h = 0; h < nha; ++h) g[h] += e * r[h];
                }
            }
        }
        // bra step: out[ab][cd] += sum_h E^ab[ab][h] G[cd][h]
        const double* Eb = &bra.E[(size_t)ib * nab * nha];
        for (int ab = 0; ab < nab; ++ab) {
            const double* e = Eb + (size_t)ab * nha;
            double* o = out + (size_t)ab * ncd;
            for (int cd = 0; cd < ncd; ++cd) {
                const double* g = G + (size_t)cd * nha;
                // four partial sums: independent chains the compiler can vectorise
                double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
                int h = 0;
                for (; h + 3 < nha; h += 4) {
                    a0 += e[h] * g[h];
                    a1 += e[h + 1] * g[h + 1];
                    a2 += e[h + 2] * g[h + 2];
                    a3 += e[h + 3] * g[h + 3];
                }
                for (; h < nha; ++h) a0 += e[h] * g[h];
                o[cd] += (a0 + a1) + (a2 + a3);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Compile-time sizes for the common classes (Claude Generated, Oct 2026)
//
// mdQuartet above is the readable form. For the classes of s/p/d bases (LAB, LCD <= 4)
// the same code is instantiated with LAB and LCD as template parameters: the Hermite
// simplex sizes, the R box offsets and the loop trip counts become constants, so the
// short inner loops (1, 4, 10, 20 or 35 Hermite terms) are unrolled and vectorised.
// The arithmetic and its order are those of mdQuartet; mdQuartetDispatch picks the
// instantiation and falls back to mdQuartet for anything larger.
// ---------------------------------------------------------------------------

/// Hermite simplex t+u+v <= L in the loop order of hermiteSimplex: offsets in a box of
/// side S and the ket sign (-1)^(t+u+v).
#if defined(__GNUC__) && !defined(__clang__)
#define QM_UNROLL _Pragma("GCC unroll 16")
#elif defined(__clang__)
#define QM_UNROLL _Pragma("unroll")
#else
#define QM_UNROLL
#endif
// The R recursion is called once per primitive quartet; the compiler's unit-growth limit
// otherwise keeps it out of line (a call per primitive quartet).
#if defined(__GNUC__) || defined(__clang__)
#define QM_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#define QM_ALWAYS_INLINE inline
#endif

template <int L, int S>
struct HermiteBox {
    static constexpr int n = (L + 1) * (L + 2) * (L + 3) / 6;
    std::array<int, n> off{};
    std::array<double, n> sgn{};
    constexpr HermiteBox()
    {
        int k = 0;
        for (int t = 0; t <= L; ++t)
            for (int u = 0; u <= L - t; ++u)
                for (int v = 0; v <= L - t - u; ++v) {
                    off[k] = (t * S + u) * S + v;
                    sgn[k] = ((t + u + v) & 1) ? -1.0 : 1.0;
                    ++k;
                }
    }
};

/// buildR0Raw with L known at compile time (same recursion, same order).
template <int L>
QM_ALWAYS_INLINE void buildR0T(double Wx, double Wy, double Wz, double rho, const double* F, double* R0, double* lvl)
{
    constexpr int s = L + 1;
    const double b = -2.0 * rho;
    double bpow[L + 1];
    bpow[0] = 1.0;
    for (int n = 1; n <= L; ++n) bpow[n] = bpow[n - 1] * b;
    double* c = (L % 2 == 0) ? R0 : lvl;
    double* x = (L % 2 == 0) ? lvl : R0;
    c[0] = bpow[L] * F[L];
    // All trip counts are compile-time constants once the outer loops are unrolled; the
    // pragmas ask GCC for the complete unroll (the recursion is ~L^4/24 FMAs, the loop
    // control otherwise costs more than the arithmetic).
    QM_UNROLL
    for (int n = L - 1; n >= 0; --n) {
        const int top = L - n;
        x[0] = bpow[n] * F[n];
        QM_UNROLL
        for (int t = 0; t <= top; ++t) {
            QM_UNROLL
            for (int u = 0; u <= top - t; ++u) {
                const int base = (t * s + u) * s;
                const int vmax = top - t - u;
                if (t + u > 0) {
                    if (u > 0)
                        x[base] = (u >= 2 ? (u - 1) * c[base - 2 * s] : 0.0) + Wy * c[base - s];
                    else
                        x[base] = (t >= 2 ? (t - 1) * c[base - 2 * s * s] : 0.0) + Wx * c[base - s * s];
                }
                if (vmax >= 1) x[base + 1] = 0.0 + Wz * c[base];
                QM_UNROLL
                for (int v = 2; v <= vmax; ++v) x[base + v] = (v - 1) * c[base + v - 2] + Wz * c[base + v - 1];
            }
        }
        std::swap(c, x);
    }
}

template <int LAB, int LCD>
void mdQuartetT(const MDPair& bra, const MDPair& ket, double prim_thr, double* out, MDScratch& w)
{
    constexpr int Lt = LAB + LCD, s = Lt + 1;
    static constexpr HermiteBox<LAB, s> HA{};
    static constexpr HermiteBox<LCD, s> HC{};
    constexpr int nha = HermiteBox<LAB, s>::n, nhc = HermiteBox<LCD, s>::n;
    const int nab = bra.nout, ncd = ket.nout;
    w.G.resize((size_t)ncd * nha);
    double* G = w.G.data();
    double F[Lt + 1];
    double R[s * s * s], Rl[s * s * s];
    for (int x = 0; x < nab * ncd; ++x) out[x] = 0.0;
    static const double pref0 = 2.0 * std::pow(PI, 2.5);

    const int npb = (int)bra.p.size(), npk = (int)ket.p.size();
    const double ket_max = npk ? ket.schwarz[0] : 0.0;
    for (int ib = 0; ib < npb; ++ib) {
        if (bra.schwarz[ib] * ket_max < prim_thr) break;
        std::fill(G, G + (size_t)ncd * nha, 0.0);
        const double p = bra.p[ib], ipb = bra.ip[ib];
        for (int ik = 0; ik < npk; ++ik) {
            if (bra.schwarz[ib] * ket.schwarz[ik] < prim_thr) break;
            const double q = ket.p[ik];
            const double inv = 1.0 / (p + q);
            const double rho = p * q * inv;
            const double Wx = bra.Px[ib] - ket.Px[ik], Wy = bra.Py[ib] - ket.Py[ik], Wz = bra.Pz[ib] - ket.Pz[ik];
            boysFast(Lt, rho * (Wx * Wx + Wy * Wy + Wz * Wz), F);
            const double pref = pref0 * ipb * ket.ip[ik] * std::sqrt(inv);  // 2 pi^(5/2) / (p q sqrt(p+q))
            for (int n = 0; n <= Lt; ++n) F[n] *= pref;
            buildR0T<Lt>(Wx, Wy, Wz, rho, F, R, Rl);
            const double* Ek = &ket.E[(size_t)ik * ncd * nhc];
            QM_UNROLL
            for (int k = 0; k < nhc; ++k) {
                const double* Rk = R + HC.off[k];
                double r[nha];
                QM_UNROLL
                for (int h = 0; h < nha; ++h) r[h] = HC.sgn[k] * Rk[HA.off[h]];
                for (int cd = 0; cd < ncd; ++cd) {
                    const double e = Ek[(size_t)cd * nhc + k];
                    if (e == 0.0) continue;
                    double* g = G + (size_t)cd * nha;
                    for (int h = 0; h < nha; ++h) g[h] += e * r[h];
                }
            }
        }
        const double* Eb = &bra.E[(size_t)ib * nab * nha];
        for (int ab = 0; ab < nab; ++ab) {
            const double* e = Eb + (size_t)ab * nha;
            double* o = out + (size_t)ab * ncd;
            for (int cd = 0; cd < ncd; ++cd) {
                const double* g = G + (size_t)cd * nha;
                double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
                int h = 0;
                for (; h + 3 < nha; h += 4) {
                    a0 += e[h] * g[h];
                    a1 += e[h + 1] * g[h + 1];
                    a2 += e[h + 2] * g[h + 2];
                    a3 += e[h + 3] * g[h + 3];
                }
                for (; h < nha; ++h) a0 += e[h] * g[h];
                o[cd] += (a0 + a1) + (a2 + a3);
            }
        }
    }
}

/// mdQuartetT for LAB, LCD <= 4 (s, p, d shells), mdQuartet otherwise.
void mdQuartetDispatch(const MDPair& bra, const MDPair& ket, double prim_thr, double* out, MDScratch& w)
{
    using Fn = void (*)(const MDPair&, const MDPair&, double, double*, MDScratch&);
    static const Fn table[5][5] = {
        { mdQuartetT<0, 0>, mdQuartetT<0, 1>, mdQuartetT<0, 2>, mdQuartetT<0, 3>, mdQuartetT<0, 4> },
        { mdQuartetT<1, 0>, mdQuartetT<1, 1>, mdQuartetT<1, 2>, mdQuartetT<1, 3>, mdQuartetT<1, 4> },
        { mdQuartetT<2, 0>, mdQuartetT<2, 1>, mdQuartetT<2, 2>, mdQuartetT<2, 3>, mdQuartetT<2, 4> },
        { mdQuartetT<3, 0>, mdQuartetT<3, 1>, mdQuartetT<3, 2>, mdQuartetT<3, 3>, mdQuartetT<3, 4> },
        { mdQuartetT<4, 0>, mdQuartetT<4, 1>, mdQuartetT<4, 2>, mdQuartetT<4, 3>, mdQuartetT<4, 4> },
    };
    if (bra.LAB <= 4 && ket.LAB <= 4)
        table[bra.LAB][ket.LAB](bra, ket, prim_thr, out, w);
    else
        mdQuartet(bra, ket, prim_thr, out, w);
}

/// activeQuartet on the Hermite-matrix kernel: same output (buf, dims d) and the same
/// side choice (the pair with more output components on the bra side).
void activeQuartetMD(const ShellOutput& o, const ShellPair& braS, const ShellPair& ketS, const MDPair& bra,
                     const MDPair& ket, std::vector<double>& buf, MDScratch& w, int d[4], double prim_thr)
{
    d[0] = o.width[braS.A];
    d[1] = o.width[braS.B];
    d[2] = o.width[ketS.A];
    d[3] = o.width[ketS.B];
    const size_t nab = bra.nout, ncd = ket.nout;
    buf.resize(nab * ncd);
    if (ncd > nab) {
        w.tmp.resize(nab * ncd);
        mdQuartetDispatch(ket, bra, prim_thr, w.tmp.data(), w);
        for (size_t x = 0; x < nab; ++x)
            for (size_t y = 0; y < ncd; ++y) buf[x * ncd + y] = w.tmp[y * nab + x];
    } else {
        mdQuartetDispatch(bra, ket, prim_thr, buf.data(), w);
    }
}

/// MD tables for every shell pair whose Schwarz factor can still pass `screening`
/// (Q_AB * max Q * pq >= screening needs Q_AB >= screening / (max Q * pq)); the
/// margin 1e-6 covers any |P| up to 1e6. Pairs without tables take the reference
/// path (activeQuartet), so the result never depends on this cut.
std::vector<MDPair> makeMDPairs(const std::vector<EriShell>& sh, const ShellOutput& o,
                                const std::vector<ShellPair>& pairs, double screening, int threads)
{
    // CURCUMA_QM_ERI_REFERENCE=1: no tables, every quartet through the reference kernel
    // (shellQuartet + transformAxis), for comparing the two kernels (diagnostic only).
    static const bool reference = [] {
        const char* e = std::getenv("CURCUMA_QM_ERI_REFERENCE");
        return e != nullptr && e[0] != '\0' && e[0] != '0';
    }();
    if (reference) return std::vector<MDPair>(pairs.size());
    double qmax = 0.0;
    for (const ShellPair& sp : pairs) qmax = std::max(qmax, sp.schwarz);
    std::vector<MDPair> md(pairs.size());
    const int np = (int)pairs.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads > 0 ? threads : 1)
#endif
    for (int i = 0; i < np; ++i) {
        if (screening > 0.0 && pairs[i].schwarz * qmax < screening * 1.0e-6) continue;
        md[i] = makeMDPair(sh, o, pairs[i]);
    }
    (void)threads;
    return md;
}

}  // namespace

ERITensor buildERI(const std::vector<GTO::Orbital>& basis, int threads, double screening,
                   const Matrix* Q)
{
    const std::vector<EriShell> sh = groupShells(basis);
    const ShellOutput o = makeShellOutput(sh, (int)basis.size(), Q);
    ERITensor eri(o.n);
    if (o.n == 0) return eri;

    const std::vector<ShellPair> pairs = makeScreeningPairs(sh, threads);
    const std::vector<MDPair> md = makeMDPairs(sh, o, pairs, screening, threads);
    const int np = (int)pairs.size();
    // Primitive quartets 1000x below the shell-quartet threshold are skipped too;
    // screening = 0 keeps every quartet (exact).
    const double prim_thr = screening * 1.0e-3;

    // Canonical shell quartets: pair(AB) <= pair(CD). Every AO quartet of a shell
    // quartet is written with set8, so A==B or AB==CD blocks just rewrite the same
    // value into already-written entries.
#ifdef _OPENMP
#pragma omp parallel num_threads(threads > 0 ? threads : 1)
#endif
    {
        std::vector<double> buf, R, tmp;
        MDScratch mw;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int i = 0; i < np; ++i) {
            const ShellPair& bra = pairs[i];
            for (int j = i; j < np; ++j) {
                const ShellPair& ket = pairs[j];
                if (bra.schwarz * ket.schwarz < screening) continue;
                int d[4];
                if (md[i].tabled && md[j].tabled)
                    activeQuartetMD(o, bra, ket, md[i], md[j], buf, mw, d, prim_thr);
                else
                    activeQuartet(sh, o, bra, ket, buf, R, tmp, d, prim_thr);
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
// Flat tables for the device ERI kernel (Claude Generated, Oct 2026, Q-G2).
// Everything here is copied out of the SAME structures buildERI uses, so the
// device sees the same shells, primitive-pair order, Hermite tables and Schwarz
// factors and screens exactly like the CPU. The CPU path itself is untouched.
// ---------------------------------------------------------------------------
bool buildEriGpuTables(const std::vector<GTO::Orbital>& basis, double screening, const Matrix* Q,
                       EriGpuTables& t, std::string& why, bool list_quartets)
{
    t = EriGpuTables();
    const std::vector<EriShell> sh = groupShells(basis);
    const ShellOutput o = makeShellOutput(sh, (int)basis.size(), Q);
    t.n_out = o.n;
    if (o.n == 0) { why = "empty basis"; return false; }
    const int ns = (int)sh.size();
    for (const EriShell& s : sh) t.maxL = std::max(t.maxL, s.L);
    if (t.maxL > 2) {
        why = "angular momentum L=" + std::to_string(t.maxL) + " (device kernel handles s, p, d)";
        return false;
    }
    if (BoysTable::NMAX < 4 * t.maxL) { why = "Boys table too short"; return false; }

    for (int A = 0; A < ns; ++A) {
        const EriShell& s = sh[A];
        const int nc = (int)s.lmn.size(), np = (int)s.exps.size();
        t.sh_L.push_back(s.L);
        t.sh_ncomp.push_back(nc);
        t.sh_nprim.push_back(np);
        t.sh_lmn_off.push_back((int)t.lmn.size());
        for (const auto& c : s.lmn) t.lmn.insert(t.lmn.end(), c.begin(), c.end());
        t.sh_coef_off.push_back((int)t.coef.size());
        for (int c = 0; c < nc; ++c)
            for (int p = 0; p < np; ++p) t.coef.push_back(s.coef[c][p]);
        t.sh_out_first.push_back(o.first[A]);
        t.sh_out_width.push_back(o.width[A]);
        if (!o.spherical || o.identity[A]) {
            t.sh_T_off.push_back(-1);
        } else {
            t.sh_T_off.push_back((int)t.T.size());
            for (int a = 0; a < nc; ++a)
                for (int j = 0; j < o.width[A]; ++j) t.T.push_back(o.T[A](a, j));
        }
    }

    const std::vector<ShellPair> pairs = makeScreeningPairs(sh);
    const int np = (int)pairs.size();
    for (const ShellPair& sp : pairs) {
        t.pr_A.push_back(sp.A);
        t.pr_B.push_back(sp.B);
        t.pr_schwarz.push_back(sp.schwarz);
        t.pr_prim_off.push_back((int)t.pp_ia.size());
        t.pr_nprim.push_back((int)sp.prims.size());
        for (const PrimPair& pp : sp.prims) {
            t.pp_ia.push_back(pp.ia);
            t.pp_ib.push_back(pp.ib);
            t.pp_p.push_back(pp.p);
            t.pp_P.push_back(pp.Px);
            t.pp_P.push_back(pp.Py);
            t.pp_P.push_back(pp.Pz);
            t.pp_schwarz.push_back(pp.schwarz);
            t.pp_E_off.push_back((long long)t.E.size());
            t.E.insert(t.E.end(), pp.Ex.begin(), pp.Ex.end());
            t.E.insert(t.E.end(), pp.Ey.begin(), pp.Ey.end());
            t.E.insert(t.E.end(), pp.Ez.begin(), pp.Ez.end());
        }
    }

    t.prim_thr = screening * 1.0e-3;  // as in buildERI
    t.pref = 2.0 * std::pow(PI, 2.5);
    t.screening = screening;

    // Q-G3 (Claude Generated, Oct 2026): shell pairs bucketed by pair class LA*3 + LB, CPU
    // order inside a bucket. Integral-direct enumerates bucket x bucket on the device.
    {
        std::vector<int> cnt(9, 0);
        for (int i = 0; i < np; ++i) ++cnt[sh[pairs[i].A].L * 3 + sh[pairs[i].B].L];
        t.pc_off.assign(10, 0);
        for (int k = 0; k < 9; ++k) t.pc_off[k + 1] = t.pc_off[k] + cnt[k];
        std::vector<int> fillp(t.pc_off.begin(), t.pc_off.end() - 1);
        t.pc_pairs.resize(np);
        for (int i = 0; i < np; ++i) t.pc_pairs[fillp[sh[pairs[i].A].L * 3 + sh[pairs[i].B].L]++] = i;
    }
    if (!list_quartets) return true;

    // Canonical quartets (i <= j) passing the shell-level Schwarz test, bucketed by
    // class (LA, LB, LC, LD); the order inside a bucket is the CPU loop order.
    const int nkey = 81;  // L <= 2 per shell
    auto key = [&](int i, int j) {
        return ((sh[pairs[i].A].L * 3 + sh[pairs[i].B].L) * 3 + sh[pairs[j].A].L) * 3 + sh[pairs[j].B].L;
    };
    std::vector<long long> count(nkey, 0);
    for (int i = 0; i < np; ++i)
        for (int j = i; j < np; ++j)
            if (!(pairs[i].schwarz * pairs[j].schwarz < screening)) ++count[key(i, j)];
    std::vector<long long> start(nkey, 0), fill(nkey, 0);
    long long total = 0;
    for (int k = 0; k < nkey; ++k) {
        start[k] = total;
        total += count[k];
    }
    t.q_bra.resize(total);
    t.q_ket.resize(total);
    for (int i = 0; i < np; ++i)
        for (int j = i; j < np; ++j) {
            if (pairs[i].schwarz * pairs[j].schwarz < screening) continue;
            const int k = key(i, j);
            const long long at = start[k] + fill[k]++;
            t.q_bra[at] = i;
            t.q_ket[at] = j;
        }
    for (int k = 0; k < nkey; ++k) {
        if (count[k] == 0) continue;
        t.class_off.push_back(start[k]);
        t.class_L.push_back(k / 27);
        t.class_L.push_back((k / 9) % 3);
        t.class_L.push_back((k / 3) % 3);
        t.class_L.push_back(k % 3);
    }
    t.class_off.push_back(total);
    return true;
}

// Q-G7 (Claude Generated, Oct 2026): tables for the device 2e gradient. The shell/pair/primitive
// layout is the integral-direct one above, in the cartesian basis (the gradient is formed with the
// cartesian density, as gradientTwoElectron); only the Hermite tables are replaced by the raised
// ones (L+1 on both shells, the same hermiteCoeffs calls as makeShellPair(..., 1, 1)), in the same
// primitive-pair order, so pp_E_off keeps pointing at the matching primitive pair.
bool buildEriGradGpuTables(const std::vector<GTO::Orbital>& basis, double screening, EriGpuTables& t,
                           std::string& why)
{
    if (!buildEriGpuTables(basis, screening, nullptr, t, why, /*list_quartets=*/false)) return false;
    if (BoysTable::NMAX < 4 * t.maxL + 1) { why = "Boys table too short for the gradient"; return false; }
    const std::vector<EriShell> sh = groupShells(basis);
    const int ns = (int)sh.size();
    t.sh_atom.resize(ns);
    for (int s = 0; s < ns; ++s) t.sh_atom[s] = basis[sh[s].first].atom;
    const int np = (int)t.pr_A.size();
    t.E.clear();
    t.pp_ea.assign(t.pp_ia.size(), 0.0);
    t.pp_eb.assign(t.pp_ia.size(), 0.0);
    for (int i = 0; i < np; ++i) {
        const EriShell& a = sh[t.pr_A[i]];
        const EriShell& b = sh[t.pr_B[i]];
        for (int k = t.pr_prim_off[i]; k < t.pr_prim_off[i] + t.pr_nprim[i]; ++k) {
            const double ea = a.exps[t.pp_ia[k]], eb = b.exps[t.pp_ib[k]];
            double p, Px, Py, Pz;
            const double K = gaussianProductK(ea, eb, a.x, a.y, a.z, b.x, b.y, b.z, p, Px, Py, Pz);
            t.pp_ea[k] = ea;
            t.pp_eb[k] = eb;
            t.pp_E_off[k] = (long long)t.E.size();
            const auto ex = flattenHermite(hermiteCoeffs(a.L + 1, b.L + 1, Px - a.x, Px - b.x, p, K));
            const auto ey = flattenHermite(hermiteCoeffs(a.L + 1, b.L + 1, Py - a.y, Py - b.y, p, 1.0));
            const auto ez = flattenHermite(hermiteCoeffs(a.L + 1, b.L + 1, Pz - a.z, Pz - b.z, p, 1.0));
            t.E.insert(t.E.end(), ex.begin(), ex.end());
            t.E.insert(t.E.end(), ey.begin(), ey.end());
            t.E.insert(t.E.end(), ez.begin(), ez.end());
        }
    }
    return true;
}

EriGpuTables::EriGpuTables() = default;
EriGpuTables::EriGpuTables(const EriGpuTables&) = default;
EriGpuTables& EriGpuTables::operator=(const EriGpuTables&) = default;
EriGpuTables& EriGpuTables::operator=(EriGpuTables&&) noexcept = default;
EriGpuTables::~EriGpuTables() = default;

QMGpuEriTables EriGpuTables::view() const
{
    QMGpuEriTables v;
    v.n_out = n_out;
    v.nshell = (int)sh_L.size();
    v.npair = (int)pr_A.size();
    v.nprimpair = (int)pp_ia.size();
    v.prim_thr = prim_thr;
    v.pref = pref;
    v.sh_L = sh_L.data();
    v.sh_ncomp = sh_ncomp.data();
    v.sh_nprim = sh_nprim.data();
    v.sh_lmn_off = sh_lmn_off.data();
    v.sh_coef_off = sh_coef_off.data();
    v.sh_out_first = sh_out_first.data();
    v.sh_out_width = sh_out_width.data();
    v.sh_T_off = sh_T_off.data();
    v.lmn = lmn.data();
    v.n_lmn = (int)lmn.size();
    v.coef = coef.data();
    v.n_coef = (int)coef.size();
    v.T = T.data();
    v.n_T = (int)T.size();
    v.pr_A = pr_A.data();
    v.pr_B = pr_B.data();
    v.pr_prim_off = pr_prim_off.data();
    v.pr_nprim = pr_nprim.data();
    v.pp_ia = pp_ia.data();
    v.pp_ib = pp_ib.data();
    v.pp_p = pp_p.data();
    v.pp_P = pp_P.data();
    v.pp_schwarz = pp_schwarz.data();
    v.pp_E_off = pp_E_off.data();
    v.E = E.data();
    v.n_E = (long long)E.size();
    v.q_bra = q_bra.data();
    v.q_ket = q_ket.data();
    v.nquartet = (long long)q_bra.size();
    v.nclass = (int)class_L.size() / 4;
    v.class_off = class_off.data();
    v.class_L = class_L.data();
    const BoysTable& b = boysTable();
    v.boys_f = b.f.data();
    v.boys_npts = b.npts;
    v.boys_ntab = BoysTable::NTAB;
    v.boys_ntaylor = BoysTable::NTAYLOR;
    v.boys_dT = BoysTable::dT;
    v.boys_Tmax = BoysTable::Tmax;
    v.device_mem_limit_mb = device_mem_limit_mb;
    v.screening = screening;
    v.pr_schwarz = pr_schwarz.data();
    v.pc_off = pc_off.empty() ? nullptr : pc_off.data();
    v.pc_pairs = pc_pairs.data();
    return v;
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

namespace {
/// Index of the lowest set bit of a non-zero word.
inline int lowestBit(uint64_t x)
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int k = 0;
    while (!(x & 1)) { x >>= 1; ++k; }
    return k;
#endif
}
}  // namespace

struct DirectJK::Impl {
    std::vector<EriShell> sh;
    ShellOutput out;
    double screening = 0.0;
    // Pair tables, built once by the first build() with its thread count (Oct 2026; the
    // constructor has no thread count, and at 2000 basis functions the tables take seconds).
    mutable std::once_flag once;
    mutable std::vector<ShellPair> pairs;
    mutable std::vector<MDPair> md;  ///< Hermite-matrix tables (mdQuartet), same order as pairs
    mutable std::vector<int> order;  ///< pair indices by decreasing Schwarz factor (quartet enumeration)
    mutable std::vector<int> pos;    ///< inverse of order: sorted position of each pair
    mutable double max_schwarz = 0.0;
    void setup(int threads) const
    {
        std::call_once(once, [&] {
            pairs = makeScreeningPairs(sh, threads);
            for (const ShellPair& sp : pairs) max_schwarz = std::max(max_schwarz, sp.schwarz);
            md = makeMDPairs(sh, out, pairs, screening, threads);
            order.resize(pairs.size());
            for (size_t i = 0; i < order.size(); ++i) order[i] = (int)i;
            std::stable_sort(order.begin(), order.end(),
                             [&](int x, int y) { return pairs[x].schwarz > pairs[y].schwarz; });
            pos.resize(order.size());
            for (size_t k = 0; k < order.size(); ++k) pos[order[k]] = (int)k;
        });
    }
};

DirectJK::DirectJK(const std::vector<GTO::Orbital>& basis, double screening, const Matrix* Q)
    : m_screening(screening)
{
    auto impl = std::make_shared<Impl>();
    impl->sh = groupShells(basis);
    impl->out = makeShellOutput(impl->sh, (int)basis.size(), Q);
    impl->screening = screening;
    m_impl = impl;
}

int DirectJK::n() const { return m_impl ? m_impl->out.n : 0; }

long DirectJK::build(const Matrix& P, Matrix& J, Matrix& K, int threads, double thr_override, bool prim_density) const
{
    const int n = this->n();
    J = Matrix::Zero(n, n);
    K = Matrix::Zero(n, n);
    if (!m_impl || n == 0) return 0;
    const Impl& I = *m_impl;
    I.setup(threads);
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
    const double thr = thr_override >= 0.0 ? thr_override : m_screening;
    const double prim_thr = thr * 1.0e-3;  // primitive quartets, as in buildERI
    long computed = 0;

#ifdef _OPENMP
#pragma omp parallel num_threads(threads > 0 ? threads : 1) reduction(+ : computed)
#endif
    {
        Matrix Jt = Matrix::Zero(n, n), Kt = Matrix::Zero(n, n);
        std::vector<double> buf, R, tmp;
        MDScratch mw;
        // Canonical quartets {i, j} found over the pairs sorted by Schwarz factor
        // (Claude Generated, Oct 2026): for bra pair i at sorted position a, the kets are the
        // positions b <= a (Q_b >= Q_a), and since Q_a Q_b pq <= Q_a Q_b max|P| the scan stops
        // at the first b with Q_a Q_b max|P| < thr. The old i <= j double loop tested all
        // npair^2/2 pairs (1e11 tests per build at 2000 basis functions). The kets found are
        // then processed in pair-index order (shells close together in J/K, cache locality),
        // and the bra rows in pair-index order too; the quartets computed and the
        // per-quartet screen are the same as before, only their order differs.
        std::vector<int> kets;
        std::vector<uint64_t> kbits((size_t)np / 64 + 1, 0);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (int i = 0; i < np; ++i) {
            const ShellPair& bra = I.pairs[i];
            if (bra.schwarz * I.max_schwarz * Pall < thr) continue;
            const int A = bra.A, B = bra.B;
            const int a = I.pos[i];
            // ket indices into a bitmap and read back in increasing order: O(row + npair/64),
            // cheaper than sorting long rows
            int wlo = INT32_MAX, whi = -1;
            for (int b = 0; b <= a; ++b) {
                const int j = I.order[b];
                if (bra.schwarz * I.pairs[j].schwarz * Pall < thr) break;
                kbits[j >> 6] |= (uint64_t)1 << (j & 63);
                wlo = std::min(wlo, j >> 6);
                whi = std::max(whi, j >> 6);
            }
            kets.clear();
            for (int w = wlo; w <= whi; ++w) {
                uint64_t word = kbits[w];
                kbits[w] = 0;
                while (word) {
                    kets.push_back((w << 6) + lowestBit(word));
                    word &= word - 1;
                }
            }
            for (const int j : kets) {
                const ShellPair& ket = I.pairs[j];
                const int C = ket.A, D = ket.B;
                const double pq = std::max({ Pmax(A, B), Pmax(C, D), Pmax(A, C), Pmax(A, D), Pmax(B, C), Pmax(B, D) });
                if (bra.schwarz * ket.schwarz * pq < thr) continue;
                int d[4];
                // same expression as the device kernels (quartetOf in cuda/qm_gpu_eri_kernels.cuh)
                const double pthr = (prim_density && pq > 0.0) ? prim_thr / pq : prim_thr;
                if (I.md[i].tabled && I.md[j].tabled)
                    activeQuartetMD(o, bra, ket, I.md[i], I.md[j], buf, mw, d, pthr);
                else
                    activeQuartet(I.sh, o, bra, ket, buf, R, tmp, d, pthr);
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
    // Oct 2026 (WP8): fixed round-robin chunks of rows, each with its own accumulator, added in chunk order, so the
    // result is bitwise the same for every thread count (before: per-thread sums under omp critical).
    const int nchunk = std::max(1, std::min(64, n));
    std::vector<Matrix> gchunk(nchunk, Matrix::Zero(natoms, 3));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
#endif
    for (int ch = 0; ch < nchunk; ++ch) {
        Matrix& g = gchunk[ch];
        for (int mu = ch; mu < n; mu += nchunk) {
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
    }
    for (int ch = 0; ch < nchunk; ++ch) grad += gchunk[ch];
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
namespace {
// Shared driver of gradientTwoElectron, gradientTwoElectronExchange and gradientTwoElectronMP2: the
// canonical-quartet sum 1/2 sum D (mn|ls)^X with the quartet weights D filled per shell quartet by
// `fill(braA, braB, ketA, ketB, Dq)` (component layout Dq[((a*nB+b)*nC+c)*nD+d], cartesian basis). The
// fill must describe a D with the full 8-fold permutational symmetry of (mn|ls).
// Templated driver: Claude Generated (Oct 2026, MP2 gradient); the loop itself is the WP8 one unchanged.
template <class Fill>
Matrix gradientTwoElectronDriver(const std::vector<GTO::Orbital>& basis, const std::vector<EriShell>& sh,
                                 int natoms, int threads, Fill fill)
{
    Matrix grad = Matrix::Zero(natoms, 3);
    if (basis.empty()) return grad;
    const int ns = (int)sh.size();
    std::vector<int> atomOf(ns);
    for (int s = 0; s < ns; ++s) atomOf[s] = basis[sh[s].first].atom;

    std::vector<ShellPair> pairs;
    pairs.reserve((size_t)ns * (ns + 1) / 2);
    for (int A = 0; A < ns; ++A)
        for (int B = A; B < ns; ++B)
            pairs.push_back(makeShellPair(sh, A, B, 1, 1));
    const int np = (int)pairs.size();

    // Oct 2026 (WP8): fixed round-robin chunks of bra pairs with their own accumulators, added in chunk order:
    // bitwise the same result for every thread count (before: per-thread sums under omp critical).
    const int nchunk = std::max(1, std::min(256, np));
    std::vector<Matrix> gchunk(nchunk, Matrix::Zero(natoms, 3));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
#endif
    for (int ch = 0; ch < nchunk; ++ch) {
        Matrix& g = gchunk[ch];
        std::vector<double> Dq;
        for (int i = ch; i < np; i += nchunk) {
            const ShellPair& bra = pairs[i];
            for (int j = i; j < np; ++j) {
                const ShellPair& ket = pairs[j];
                const int at[4] = { atomOf[bra.A], atomOf[bra.B], atomOf[ket.A], atomOf[ket.B] };
                if (at[0] == at[1] && at[0] == at[2] && at[0] == at[3]) continue;
                fill(bra.A, bra.B, ket.A, ket.B, Dq);
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
    }
    for (int ch = 0; ch < nchunk; ++ch) grad += gchunk[ch];
    (void)threads;
    return grad;
}

// gradientTwoElectron (exchange_only = false, the WP8 expression unchanged) and gradientTwoElectronExchange
// (Claude Generated, Oct 2026, RI-J gradient: only the exchange part D = -1/4 (P_ml P_ns + P_ms P_nl)).
// WP8 KS (Oct 2026): kx scales the exchange part (a_x of a hybrid, 0 for LDA/GGA); kx = 1 keeps the HF expression
// bit for bit.
Matrix gradientTwoElectronImpl(const std::vector<GTO::Orbital>& basis, const Matrix& P,
                               int natoms, int threads, bool exchange_only, double kx = 1.0)
{
    if (basis.empty()) return Matrix::Zero(natoms, 3);
    const std::vector<EriShell> sh = groupShells(basis);
    auto fill = [&](int sA, int sB, int sC, int sD, std::vector<double>& Dq) {
        const EriShell& A = sh[sA];
        const EriShell& B = sh[sB];
        const EriShell& C = sh[sC];
        const EriShell& D = sh[sD];
        const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
        Dq.resize((size_t)nA * nB * nC * nD);
        for (int a = 0; a < nA; ++a)
            for (int b = 0; b < nB; ++b)
                for (int c = 0; c < nC; ++c)
                    for (int d = 0; d < nD; ++d) {
                        const int m = A.first + a, n = B.first + b, l = C.first + c, s = D.first + d;
                        if (exchange_only)
                            Dq[((size_t)(a * nB + b) * nC + c) * nD + d] =
                                -0.25 * (P(m, l) * P(n, s) + P(m, s) * P(n, l));
                        else if (kx == 1.0)
                            Dq[((size_t)(a * nB + b) * nC + c) * nD + d] =
                                P(m, n) * P(l, s) - 0.25 * (P(m, l) * P(n, s) + P(m, s) * P(n, l));
                        else if (kx == 0.0)
                            Dq[((size_t)(a * nB + b) * nC + c) * nD + d] = P(m, n) * P(l, s);
                        else
                            Dq[((size_t)(a * nB + b) * nC + c) * nD + d] =
                                P(m, n) * P(l, s) - 0.25 * kx * (P(m, l) * P(n, s) + P(m, s) * P(n, l));
                    }
    };
    return gradientTwoElectronDriver(basis, sh, natoms, threads, fill);
}

}  // namespace

Matrix gradientTwoElectron(const std::vector<GTO::Orbital>& basis, const Matrix& P,
                           int natoms, int threads)
{
    return gradientTwoElectronImpl(basis, P, natoms, threads, false);
}

// Declared in qm_ri.h (qm_integrals.hpp is included by the CUDA kernel units; a change there would
// rebuild them). Claude Generated (Oct 2026, RI-J gradient).
Matrix gradientTwoElectronExchange(const std::vector<GTO::Orbital>& basis, const Matrix& P,
                                   int natoms, int threads)
{
    return gradientTwoElectronImpl(basis, P, natoms, threads, true);
}

// Declared in qm_ri.h (same reason). Claude Generated (Oct 2026, WP8 Kohn-Sham gradient).
Matrix gradientTwoElectronScaled(const std::vector<GTO::Orbital>& basis, const Matrix& P, double kx,
                                 int natoms, int threads)
{
    return gradientTwoElectronImpl(basis, P, natoms, threads, false, kx);
}

// Declared in qm_ri.h (same reason). Claude Generated (Nov 2026, stage 5b): UHF/UKS two-electron gradient,
//   E_2 = 1/2 sum Pt_mn Pt_ls (mn|ls) - (kx/2) sum_s sum Ps_ml Ps_ns (mn|ls),   Pt = Pa + Pb,
//   D   = Pt_mn Pt_ls - (kx/2) sum_s (Ps_ml Ps_ns + Ps_ms Ps_nl)
// (Pa = Pb = P/2 gives the RHF weight of gradientTwoElectronScaled). Same canonical-quartet driver.
Matrix gradientTwoElectronUHF(const std::vector<GTO::Orbital>& basis, const Matrix& Pa, const Matrix& Pb, double kx,
                              int natoms, int threads)
{
    if (basis.empty()) return Matrix::Zero(natoms, 3);
    const std::vector<EriShell> sh = groupShells(basis);
    const Matrix Pt = Pa + Pb;
    auto fill = [&](int sA, int sB, int sC, int sD, std::vector<double>& Dq) {
        const EriShell& A = sh[sA];
        const EriShell& B = sh[sB];
        const EriShell& C = sh[sC];
        const EriShell& D = sh[sD];
        const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
        Dq.resize((size_t)nA * nB * nC * nD);
        for (int a = 0; a < nA; ++a)
            for (int b = 0; b < nB; ++b)
                for (int c = 0; c < nC; ++c)
                    for (int d = 0; d < nD; ++d) {
                        const int m = A.first + a, n = B.first + b, l = C.first + c, s = D.first + d;
                        const double kab = Pa(m, l) * Pa(n, s) + Pa(m, s) * Pa(n, l) + Pb(m, l) * Pb(n, s)
                            + Pb(m, s) * Pb(n, l);
                        Dq[((size_t)(a * nB + b) * nC + c) * nD + d] = Pt(m, n) * Pt(l, s) - 0.5 * kx * kab;
                    }
    };
    return gradientTwoElectronDriver(basis, sh, natoms, threads, fill);
}

// Declared in qm_mp2.h (same reason as above). Claude Generated (Oct 2026, MP2 gradient).
// D = P P - 1/4 (P P + P P)                                  (RHF, as gradientTwoElectron)
//   + (P2_mn P_ls + P_mn P2_ls) - 1/4 (P2_ml P_ns + P_ml P2_ns + P2_ms P_nl + P_ms P2_nl)   (P2 x P_HF)
//   + Dns, the non-separable MP2 part, given in the ACTIVE basis and taken to the cartesian one per
//     shell quartet with the transposed blocks of Q (Gamma_cart = Q Gamma_act Q^T on every index).
Matrix gradientTwoElectronMP2(const std::vector<GTO::Orbital>& basis, const Matrix& P, const Matrix& P2,
                              const double* Dns, const Matrix* Q, int natoms, int threads)
{
    if (basis.empty()) return Matrix::Zero(natoms, 3);
    const std::vector<EriShell> sh = groupShells(basis);
    const ShellOutput so = makeShellOutput(sh, (int)basis.size(), Q);
    const int na = so.n;
    std::vector<Eigen::MatrixXd> Tt(sh.size());
    for (size_t s = 0; s < sh.size(); ++s)
        if (so.spherical && !so.identity[s]) Tt[s] = so.T[s].transpose();  // (width x ncomp): active -> cart
    auto fill = [&](int sA, int sB, int sC, int sD, std::vector<double>& Dq) {
        const EriShell& A = sh[sA];
        const EriShell& B = sh[sB];
        const EriShell& C = sh[sC];
        const EriShell& D = sh[sD];
        const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
        // non-separable block: active-basis read, then active -> cartesian on each index
        thread_local std::vector<double> blk, tmp;
        int d[4] = { so.width[sA], so.width[sB], so.width[sC], so.width[sD] };
        const int f0 = so.first[sA], f1 = so.first[sB], f2 = so.first[sC], f3 = so.first[sD];
        blk.resize((size_t)d[0] * d[1] * d[2] * d[3]);
        for (int a = 0; a < d[0]; ++a)
            for (int b = 0; b < d[1]; ++b)
                for (int c = 0; c < d[2]; ++c) {
                    const double* src = Dns + ((((size_t)(f0 + a) * na + (f1 + b)) * na + (f2 + c)) * na + f3);
                    double* dst = &blk[(((size_t)a * d[1] + b) * d[2] + c) * d[3]];
                    for (int e = 0; e < d[3]; ++e) dst[e] = src[e];
                }
        const int s4[4] = { sA, sB, sC, sD };
        for (int k = 0; k < 4; ++k)
            if (so.spherical && !so.identity[s4[k]]) transformAxis(blk, tmp, d, k, Tt[s4[k]]);
        Dq.resize((size_t)nA * nB * nC * nD);
        for (int a = 0; a < nA; ++a)
            for (int b = 0; b < nB; ++b)
                for (int c = 0; c < nC; ++c)
                    for (int e = 0; e < nD; ++e) {
                        const int m = A.first + a, n = B.first + b, l = C.first + c, s = D.first + e;
                        const size_t idx = ((size_t)(a * nB + b) * nC + c) * nD + e;
                        Dq[idx] = P(m, n) * P(l, s) - 0.25 * (P(m, l) * P(n, s) + P(m, s) * P(n, l))
                                + (P2(m, n) * P(l, s) + P(m, n) * P2(l, s))
                                - 0.25 * (P2(m, l) * P(n, s) + P(m, l) * P2(n, s) + P2(m, s) * P(n, l) + P(m, s) * P2(n, l))
                                + blk[idx];
                    }
    };
    return gradientTwoElectronDriver(basis, sh, natoms, threads, fill);
}

// Declared in qm_ri.h (same reason as above). Claude Generated (Oct 2026, WP12 stage 2: double-hybrid gradient).
// The PT2 part (Dns optional); with include_ks also the KS part P x P of gradientTwoElectronScaled, so that
// computeDoubleHybridGradient needs one pass over the derivative integrals (Oct 31 2026, saves the second pass).
Matrix gradientTwoElectronDHPT2(const std::vector<GTO::Orbital>& basis, const Matrix& P, const Matrix& P2,
                                const double* Dns, const Matrix* Q, double kx, int natoms, int threads, bool include_ks)
{
    if (basis.empty()) return Matrix::Zero(natoms, 3);
    const std::vector<EriShell> sh = groupShells(basis);
    const ShellOutput so = makeShellOutput(sh, (int)basis.size(), Q);
    const int na = so.n;
    std::vector<Eigen::MatrixXd> Tt(sh.size());
    for (size_t s = 0; s < sh.size(); ++s)
        if (so.spherical && !so.identity[s]) Tt[s] = so.T[s].transpose();
    auto fill = [&](int sA, int sB, int sC, int sD, std::vector<double>& Dq) {
        const EriShell& A = sh[sA];
        const EriShell& B = sh[sB];
        const EriShell& C = sh[sC];
        const EriShell& D = sh[sD];
        const int nA = (int)A.lmn.size(), nB = (int)B.lmn.size(), nC = (int)C.lmn.size(), nD = (int)D.lmn.size();
        thread_local std::vector<double> blk, tmp;
        if (Dns) {
            int d[4] = { so.width[sA], so.width[sB], so.width[sC], so.width[sD] };
            const int f0 = so.first[sA], f1 = so.first[sB], f2 = so.first[sC], f3 = so.first[sD];
            blk.resize((size_t)d[0] * d[1] * d[2] * d[3]);
            for (int a = 0; a < d[0]; ++a)
                for (int b = 0; b < d[1]; ++b)
                    for (int c = 0; c < d[2]; ++c) {
                        const double* src = Dns + ((((size_t)(f0 + a) * na + (f1 + b)) * na + (f2 + c)) * na + f3);
                        double* dst = &blk[(((size_t)a * d[1] + b) * d[2] + c) * d[3]];
                        for (int e = 0; e < d[3]; ++e) dst[e] = src[e];
                    }
            const int s4[4] = { sA, sB, sC, sD };
            for (int k = 0; k < 4; ++k)
                if (so.spherical && !so.identity[s4[k]]) transformAxis(blk, tmp, d, k, Tt[s4[k]]);
        }
        Dq.resize((size_t)nA * nB * nC * nD);
        for (int a = 0; a < nA; ++a)
            for (int b = 0; b < nB; ++b)
                for (int c = 0; c < nC; ++c)
                    for (int e = 0; e < nD; ++e) {
                        const int m = A.first + a, n = B.first + b, l = C.first + c, s = D.first + e;
                        const size_t idx = ((size_t)(a * nB + b) * nC + c) * nD + e;
                        const double pt2 = (P2(m, n) * P(l, s) + P(m, n) * P2(l, s))
                                - 0.25 * kx * (P2(m, l) * P(n, s) + P(m, l) * P2(n, s) + P2(m, s) * P(n, l) + P(m, s) * P2(n, l))
                                + (Dns ? blk[idx] : 0.0);
                        // KS part: the expression of gradientTwoElectronImpl (general kx branch)
                        Dq[idx] = include_ks
                            ? (P(m, n) * P(l, s) - 0.25 * kx * (P(m, l) * P(n, s) + P(m, s) * P(n, l))) + pt2
                            : pt2;
                    }
    };
    return gradientTwoElectronDriver(basis, sh, natoms, threads, fill);
}

}  // namespace qmint