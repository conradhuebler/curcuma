/*
 * <Native KS-DFT 1-Electron GTO Integrals -- implementation>
 * Copyright (C) 2019 - 2026 Conrad Hbler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (WP1): Obara-Saika overlap, gradient-form kinetic,
 * McMurchie-Davidson nuclear attraction, and the cartesian->spherical d
 * transform. See dft_integrals.hpp for conventions and literature.
 *
 * This program is free software under GPL-3.0
 */

#include "dft_integrals.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

namespace dft1e {

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
static inline double primitiveOverlapGeneral(int la, int ma, int na,
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

static double contractedKineticPair(const GTO::Orbital& a, const GTO::Orbital& b)
{
    int l1, m1, n1, l2, m2, n2;
    GTO::orbitalTypeToComponents(a.type, l1, m1, n1);
    GTO::orbitalTypeToComponents(b.type, l2, m2, n2);
    double sum = 0.0;
    for (size_t ia = 0; ia < a.exponents.size(); ++ia) {
        const double alpha = a.exponents[ia];
        for (size_t ib = 0; ib < b.exponents.size(); ++ib) {
            const double beta = b.exponents[ib];
            double gamma, Px, Py, Pz;
            double K = gaussianProductK(alpha, beta, a.x, a.y, a.z, b.x, b.y, b.z,
                                        gamma, Px, Py, Pz);
            if (K == 0.0) continue;
            const double ca = a.coefficients[ia];
            const double cb = b.coefficients[ib];
            const double PAx = Px - a.x, PBx = Px - b.x;
            const double PAy = Py - a.y, PBy = Py - b.y;
            const double PAz = Pz - a.z, PBz = Pz - b.z;

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

            sum += 0.5 * ca * cb * (tx + ty + tz);
        }
    }
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
static std::vector<double> boysArray(int maxN, double T)
{
    std::vector<double> F(maxN + 1, 0.0);
    if (T < 1e-14) {
        for (int n = 0; n <= maxN; ++n) F[n] = 1.0 / (2.0 * n + 1.0);
        return F;
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
        return F;
    }
    const int M = maxN + 25;
    std::vector<long double> G(M + 1, 0.0L);
    const long double eT = expl(-(long double)T);
    for (int n = M - 1; n >= 0; --n)
        G[n] = (2.0L * T * G[n + 1] + eT) / (2.0L * n + 1.0L);
    for (int n = 0; n <= maxN; ++n) F[n] = (double)G[n];
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

// Contracted (a b | c d) over the primitive contractions of the four AOs.
static double contractedERIPair(const GTO::Orbital& a, const GTO::Orbital& b,
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

ERITensor buildERI(const std::vector<GTO::Orbital>& basis)
{
    const int n = (int)basis.size();
    ERITensor eri(n);
    if (n == 0) return eri;
    // Canonical quartet loop: mu<=nu, lam<=sig, and pair(mu,nu) <= pair(lam,sig).
    // Each canonical value is computed once and written to all 8 permutations.
    for (int mu = 0; mu < n; ++mu) {
        for (int nu = mu; nu < n; ++nu) {
            const long p1 = (long)mu * n + nu;  // pair index, mu<=nu
            for (int lam = 0; lam < n; ++lam) {
                for (int sig = lam; sig < n; ++sig) {
                    const long p2 = (long)lam * n + sig;  // pair index, lam<=sig
                    if (p1 > p2) continue;  // enforce pair(mu,nu) <= pair(lam,sig)
                    const double v = contractedERIPair(basis[mu], basis[nu],
                                                       basis[lam], basis[sig]);
                    eri.set8(mu, nu, lam, sig, v);
                }
            }
        }
    }
    return eri;
}

Matrix buildCoulomb(const ERITensor& eri, const Matrix& P)
{
    const int n = eri.n();
    Matrix J = Matrix::Zero(n, n);
    for (int mu = 0; mu < n; ++mu) {
        for (int nu = 0; nu < n; ++nu) {
            double s = 0.0;
            for (int lam = 0; lam < n; ++lam) {
                for (int sig = 0; sig < n; ++sig) {
                    s += P(lam, sig) * eri(mu, nu, lam, sig);  // (mu nu | lam sig)
                }
            }
            J(mu, nu) = s;
        }
    }
    return J;
}

Matrix buildExchange(const ERITensor& eri, const Matrix& P)
{
    const int n = eri.n();
    Matrix K = Matrix::Zero(n, n);
    for (int mu = 0; mu < n; ++mu) {
        for (int nu = 0; nu < n; ++nu) {
            double s = 0.0;
            for (int lam = 0; lam < n; ++lam) {
                for (int sig = 0; sig < n; ++sig) {
                    s += P(lam, sig) * eri(mu, lam, nu, sig);  // (mu lam | nu sig)
                }
            }
            K(mu, nu) = s;
        }
    }
    return K;
}

ERITensor applySphericalTransformERI(const ERITensor& eriCart, const Matrix& Q)
{
    const int ncart = eriCart.n();
    if (Q.size() == 0 || ncart == 0) return eriCart;  // no d -> keep cartesian
    const int nsph = (int)Q.cols();
    ERITensor eriSph(nsph);
    // ERI_sph[i,j,k,l] = sum_{a,b,c,d} Q[a,i] Q[b,j] Q[c,k] Q[d,l] ERI_cart[a,b,c,d]
    for (int i = 0; i < nsph; ++i)
        for (int j = 0; j < nsph; ++j)
            for (int k = 0; k < nsph; ++k)
                for (int l = 0; l < nsph; ++l) {
                    double s = 0.0;
                    for (int a = 0; a < ncart; ++a) {
                        const double qa = Q(a, i);
                        if (qa == 0.0) continue;
                        for (int b = 0; b < ncart; ++b) {
                            const double qb = Q(b, j);
                            if (qb == 0.0) continue;
                            for (int c = 0; c < ncart; ++c) {
                                const double qc = Q(c, k);
                                if (qc == 0.0) continue;
                                for (int d = 0; d < ncart; ++d) {
                                    const double qd = Q(d, l);
                                    if (qd == 0.0) continue;
                                    s += qa * qb * qc * qd * eriCart(a, b, c, d);
                                }
                            }
                        }
                    }
                    eriSph.at(i, j, k, l) = s;
                }
    return eriSph;
}

}  // namespace dft1e