/*
 * <Exact STO Overlap Integrals via Lofthus/Pople Method>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Exact two-center overlap integrals for Slater-Type Orbitals using the
 * Lofthus/Pople method with Ak/Bk auxiliary integrals in confocal
 * elliptical coordinates.
 *
 * References:
 *   E. Lofthus, Mol. Phys. 5, 105 (1962) — confocal elliptical formulation
 *   J. A. Pople, D. L. Beveridge, "Approximate Molecular Orbital Theory" (1970)
 *   R. S. Mulliken et al., J. Chem. Phys. 17, 1248 (1949)
 *
 * Claude Generated: Exact STO overlap for NDDO semi-empirical methods
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace STO {
namespace Lofthus {

// =========================================================================
// Factorial and combinatorial helpers
// =========================================================================

static inline double factorial(int n)
{
    double result = 1.0;
    for (int i = 2; i <= n; ++i)
        result *= i;
    return result;
}

/// a!/b! computed without overflow for moderate arguments
static inline double factorialRatio(int a, int b)
{
    if (a == 0) a = 1;
    if (b == 0) b = 1;
    int lo = std::min(a, b) + 1;
    int hi = std::max(a, b) + 1;
    double result = 1.0;
    for (int i = lo; i < hi; ++i)
        result *= i;
    return (a >= b) ? result : 1.0 / result;
}

/// Binomial coefficients C(n,k) for k=0..n, returned as vector
static inline std::vector<int> binomialCoefficients(int n)
{
    std::vector<int> result(n + 1, 0);
    for (int k = 0; k <= n; ++k) {
        double c = 1.0;
        int lo = std::min(k, n - k);
        for (int i = 1; i <= lo; ++i) {
            c *= (n - i + 1);
            c /= i;
        }
        result[k] = static_cast<int>(std::round(c));
    }
    return result;
}

// =========================================================================
// Polynomial representation for overlap integrals
//
// A bivariate polynomial in (A_k, B_k) is stored as three parallel vectors:
//   coeffs[i] * A_{ak[i]} * B_{bk[i]}
// =========================================================================

struct OverlapPolynomial {
    std::vector<int> coeffs;
    std::vector<int> ak;  // Ak index
    std::vector<int> bk;  // Bk index
};

/// Multiply polynomial p1 by polynomial p2, result stored in p1
static inline void contractPolynomials(
    std::vector<int>& c1, std::vector<int>& a1, std::vector<int>& b1,
    const std::vector<int>& c2, const std::vector<int>& a2, const std::vector<int>& b2)
{
    std::vector<int> fc, fa, fb;
    int sz1 = static_cast<int>(c1.size());
    int sz2 = static_cast<int>(c2.size());

    for (int i = 0; i < sz1; ++i) {
        for (int j = 0; j < sz2; ++j) {
            int cc = c1[i] * c2[j];
            int aa = a1[i] + a2[j];
            int bb = b1[i] + b2[j];

            // Check if this (aa,bb) pair already exists
            bool found = false;
            for (size_t k = 0; k < fc.size(); ++k) {
                if (fa[k] == aa && fb[k] == bb) {
                    fc[k] += cc;
                    if (fc[k] == 0) {
                        fc.erase(fc.begin() + k);
                        fa.erase(fa.begin() + k);
                        fb.erase(fb.begin() + k);
                    }
                    found = true;
                    break;
                }
            }
            if (!found) {
                fc.push_back(cc);
                fa.push_back(aa);
                fb.push_back(bb);
            }
        }
    }
    c1 = fc;
    a1 = fa;
    b1 = fb;
}

// =========================================================================
// Polynomial index functions for overlap integral types
// These encode the integrand structure in confocal elliptical coordinates
// =========================================================================

/// SS overlap polynomial indices for principal quantum numbers n, m
/// Represents the integrand in confocal elliptical coordinates:
///   S_ss ∝ Σ coeffs[i] * A_{ak[i]} * B_{bk[i]}
static inline OverlapPolynomial ssOverlapIndices(int n, int m)
{
    OverlapPolynomial result;

    if (n > m) {
        // (ξ² - η²)^m * (ξ + η)^(n-m)
        auto coeffs = binomialCoefficients(m);
        for (int i = 0; i <= m; ++i) {
            if ((m - i) % 2 != 0) coeffs[i] *= -1;  // alternating sign for η²
            result.ak.push_back(i * 2);
            result.bk.push_back((m - i) * 2);
        }
        result.coeffs = std::vector<int>(coeffs.begin(), coeffs.end());

        // Multiply by (ξ + η)^(n-m)
        auto c2 = binomialCoefficients(n - m);
        std::vector<int> a2(n - m + 1), b2(n - m + 1);
        for (int i = 0; i <= n - m; ++i) {
            a2[i] = i;
            b2[i] = n - m - i;
        }
        contractPolynomials(result.coeffs, result.ak, result.bk,
                           std::vector<int>(c2.begin(), c2.end()), a2, b2);
    }
    else if (n < m) {
        // (ξ² - η²)^n * (ξ - η)^(m-n)
        auto coeffs = binomialCoefficients(n);
        for (int i = 0; i <= n; ++i) {
            if ((n - i) % 2 != 0) coeffs[i] *= -1;
            result.ak.push_back(i * 2);
            result.bk.push_back((n - i) * 2);
        }
        result.coeffs = std::vector<int>(coeffs.begin(), coeffs.end());

        // Multiply by (ξ - η)^(m-n)
        auto c2 = binomialCoefficients(m - n);
        std::vector<int> a2(m - n + 1), b2(m - n + 1);
        for (int i = 0; i <= m - n; ++i) {
            if ((m - n - i) % 2 != 0) c2[i] *= -1;  // alternating sign for -η
            a2[i] = i;
            b2[i] = m - n - i;
        }
        contractPolynomials(result.coeffs, result.ak, result.bk,
                           std::vector<int>(c2.begin(), c2.end()), a2, b2);
    }
    else {
        // n == m: just (ξ² - η²)^n
        auto coeffs = binomialCoefficients(n);
        for (int i = 0; i <= n; ++i) {
            if ((n - i) % 2 != 0) coeffs[i] *= -1;
            result.ak.push_back(i * 2);
            result.bk.push_back((n - i) * 2);
        }
        result.coeffs = std::vector<int>(coeffs.begin(), coeffs.end());
    }

    return result;
}

/// SP sigma overlap polynomial indices (S orbital on bra, P orbital on ket)
/// Following Ulysses SPOvIntIndex: base=SS(n_s, n_p-1), factor [1,-1]
/// @param n_s Principal quantum number of S orbital (bra)
/// @param n_p Principal quantum number of P orbital (ket)
static inline OverlapPolynomial spSigmaOverlapIndices(int n_s, int n_p)
{
    auto result = ssOverlapIndices(n_s, n_p - 1);

    // Ulysses SPOvIntIndex: xk=[0,1], yk=[0,1], coeffs=[1,-1]
    // Factor: A_0*B_0 - A_1*B_1
    std::vector<int> c2 = {1, -1};
    std::vector<int> a2 = {0, 1};
    std::vector<int> b2 = {0, 1};
    contractPolynomials(result.coeffs, result.ak, result.bk, c2, a2, b2);

    return result;
}

/// PS sigma overlap polynomial indices (P orbital on bra, S orbital on ket)
/// Following Ulysses PSOvIntIndex: base=SS(n_p-1, n_s), factor [1,+1]
/// @param n_p Principal quantum number of P orbital (bra)
/// @param n_s Principal quantum number of S orbital (ket)
static inline OverlapPolynomial psSigmaOverlapIndices(int n_p, int n_s)
{
    auto result = ssOverlapIndices(n_p - 1, n_s);

    // Ulysses PSOvIntIndex: xk=[0,1], yk=[0,1], coeffs=[1,1]
    // Factor: A_0*B_0 + A_1*B_1  (note: both POSITIVE, unlike SP)
    std::vector<int> c2 = {1, 1};
    std::vector<int> a2 = {0, 1};
    std::vector<int> b2 = {0, 1};
    contractPolynomials(result.coeffs, result.ak, result.bk, c2, a2, b2);

    return result;
}

/// PP sigma overlap polynomial indices
/// Uses SS(n-1, m-1) base, multiplied by (ξ² - η²) factor
static inline OverlapPolynomial ppSigmaOverlapIndices(int n1, int n2)
{
    auto result = ssOverlapIndices(n1 - 1, n2 - 1);

    // Multiply by factor representing σ-σ component: (1 - A_2*B_2^{-1}...)
    // Following Ulysses: coeffs=[1,-1], xk=[0,2], yk=[0,2]
    std::vector<int> c2 = {1, -1};
    std::vector<int> a2 = {0, 2};
    std::vector<int> b2 = {0, 2};
    contractPolynomials(result.coeffs, result.ak, result.bk, c2, a2, b2);

    return result;
}

/// PP pi overlap polynomial indices
/// Uses SS(n-1, m-1) base with different angular factor
static inline OverlapPolynomial ppPiOverlapIndices(int n1, int n2)
{
    auto result = ssOverlapIndices(n1 - 1, n2 - 1);

    // Following Ulysses PPpiOvIntIndex:
    // coeffs=[1, 1, -1, -1], xk=[2, 0, 2, 0], yk=[0, 2, 2, 0]
    // Wait, let me re-read. Ulysses has:
    // xk = [0, 0, 0, 0] initially, then xk[0]=2, xk[2]=2 → [2, 0, 2, 0]
    // yk = [0, 0, 0, 0] initially, then yk[1]=2, yk[2]=2 → [0, 2, 2, 0]
    // coeffs = [1, 1, 1, 1] initially, then coeffs[2]=-1, coeffs[3]=-1 → [1, 1, -1, -1]
    std::vector<int> c2 = {1, 1, -1, -1};
    std::vector<int> a2 = {2, 0, 2, 0};
    std::vector<int> b2 = {0, 2, 2, 0};
    contractPolynomials(result.coeffs, result.ak, result.bk, c2, a2, b2);

    return result;
}

// =========================================================================
// Ak and Bk auxiliary integrals
// =========================================================================

/// Ak integrals: A_k(p) = ∫₁^∞ ξᵏ e^{-pξ} dξ
/// Recurrence: A_k = (e^{-p} + k·A_{k-1}) / p, starting with A_0 = e^{-p}/p
static inline std::vector<double> computeAk(double p, int kmax)
{
    std::vector<double> A(kmax + 1, 0.0);
    if (std::abs(p) < 1e-10) {
        // Limit p→0: A_k → ∞, but this shouldn't happen in practice
        p = 1e-10;
    }
    double expp = std::exp(-p);
    A[0] = expp / p;
    for (int k = 1; k <= kmax; ++k) {
        A[k] = (expp + k * A[k - 1]) / p;
    }
    return A;
}

/// Bk integrals: B_k(q) = ∫₋₁^₁ ξᵏ e^{-qξ} dξ
/// Uses Pople's adaptive method for numerical stability
static inline std::vector<double> computeBk(double q, int kmax)
{
    std::vector<double> B(kmax + 1, 0.0);
    double absq = std::abs(q);

    if (absq <= 1e-5) {
        // Small q limit: B_k = 2/(k+1) for even k, 0 for odd k
        for (int k = 0; k <= kmax; ++k) {
            if (k % 2 == 0) B[k] = 2.0 / (k + 1.0);
        }
    }
    else if (absq <= 0.5) {
        // Power series expansion (6 terms)
        for (int k = 0; k <= kmax; ++k) {
            double y = 0.0;
            for (int j = 0; j < 6; ++j) {
                double term = std::pow(-q, j) * (1.0 - std::pow(-1.0, k + j + 1))
                             / (factorial(j) * (k + j + 1));
                y += term;
            }
            B[k] = y;
        }
    }
    else if (absq <= 1.0 && kmax <= 5) {
        // Recurrence (stable for small k)
        double expq = std::exp(q);
        double expmq = 1.0 / expq;
        B[0] = (expq - expmq) / q;
        for (int k = 0; k < kmax; ++k) {
            B[k + 1] = ((k + 1) * B[k] + std::pow(-1.0, k + 1) * expq - expmq) / q;
        }
    }
    else if (absq <= 1.0) {
        // Power series (7 terms) for higher k
        for (int k = 0; k <= kmax; ++k) {
            double y = 0.0;
            for (int j = 0; j < 7; ++j) {
                y += std::pow(-q, j) * (1.0 - std::pow(-1.0, k + j + 1))
                    / (factorial(j) * (k + j + 1));
            }
            B[k] = y;
        }
    }
    else if (absq <= 2.0 && kmax <= 7) {
        double expq = std::exp(q);
        double expmq = 1.0 / expq;
        B[0] = (expq - expmq) / q;
        for (int k = 0; k < kmax; ++k) {
            B[k + 1] = ((k + 1) * B[k] + std::pow(-1.0, k + 1) * expq - expmq) / q;
        }
    }
    else if (absq <= 2.0) {
        for (int k = 0; k <= kmax; ++k) {
            double y = 0.0;
            for (int j = 0; j < 12; ++j) {
                y += std::pow(-q, j) * (1.0 - std::pow(-1.0, k + j + 1))
                    / (factorial(j) * (k + j + 1));
            }
            B[k] = y;
        }
    }
    else if (absq <= 3.0 && kmax <= 10) {
        double expq = std::exp(q);
        double expmq = 1.0 / expq;
        B[0] = (expq - expmq) / q;
        for (int k = 0; k < kmax; ++k) {
            B[k + 1] = ((k + 1) * B[k] + std::pow(-1.0, k + 1) * expq - expmq) / q;
        }
    }
    else if (absq <= 3.0) {
        for (int k = 0; k <= kmax; ++k) {
            double y = 0.0;
            for (int j = 0; j < 15; ++j) {
                y += std::pow(-q, j) * (1.0 - std::pow(-1.0, k + j + 1))
                    / (factorial(j) * (k + j + 1));
            }
            B[k] = y;
        }
    }
    else {
        // |q| > 3: recurrence is stable
        double expq = std::exp(q);
        double expmq = 1.0 / expq;
        B[0] = (expq - expmq) / q;
        for (int k = 0; k < kmax; ++k) {
            B[k + 1] = ((k + 1) * B[k] + std::pow(-1.0, k + 1) * expq - expmq) / q;
        }
    }

    return B;
}

// =========================================================================
// Overlap factor (normalization)
// =========================================================================

/// OvFactor: normalization constant for STO overlap integrals
/// Combines STO normalization with confocal elliptical coordinate scaling
///
/// @param n1, l1, m1 Quantum numbers for orbital 1
/// @param zeta1 Slater exponent for orbital 1
/// @param n2, l2, m2 Quantum numbers for orbital 2
/// @param zeta2 Slater exponent for orbital 2
/// @param R Distance in Bohr
static inline double overlapFactor(int n1, int l1, int m1, double zeta1,
                                    int n2, int l2, int m2, double zeta2, double R)
{
    // Lofthus angular factors
    double x = 1.0, y = 1.0, w = 1.0;

    if (l1 == 0 && l2 == 0) {
        w = 2.0;
    }
    else if (l1 == 1 && l2 == 1) {
        if (m1 == 0 && m2 == 0) {           // pp sigma
            x = 3.0; w = 2.0;
        }
        else if (std::abs(m1) == 1 && std::abs(m2) == 1) {  // pp pi
            x = 3.0; w = 4.0;
        }
    }
    else if ((l1 == 0 && l2 == 1) || (l1 == 1 && l2 == 0)) {
        if (std::abs(m1) + std::abs(m2) == 0) {  // sp sigma
            y = 3.0; w = 2.0;
        }
    }

    int nmax = std::max(n1, n2);
    int nmin = std::min(n1, n2);
    double t = factorialRatio(2 * nmax, 2 * nmin);

    double ovfactor = (x / w) * std::sqrt(y / t)
                     * std::pow(R, n1 + n2 + 1)
                     * std::sqrt(std::pow(zeta1, 2 * n1 + 1) * std::pow(zeta2, 2 * n2 + 1))
                     / factorial(2 * nmin);

    return ovfactor;
}

// =========================================================================
// Exact overlap integral computation
// =========================================================================

/// Evaluate overlap polynomial: Σ coeffs[i] * Ak[ak[i]] * Bk[bk[i]]
static inline double evaluatePolynomial(const OverlapPolynomial& poly,
                                         const std::vector<double>& Ak,
                                         const std::vector<double>& Bk)
{
    double sum = 0.0;
    for (size_t i = 0; i < poly.coeffs.size(); ++i) {
        sum += poly.coeffs[i] * Ak[poly.ak[i]] * Bk[poly.bk[i]];
    }
    return sum;
}

/// Exact 2-center SS overlap integral
/// @param n1 Principal quantum number for orbital 1
/// @param zeta1 Slater exponent for orbital 1
/// @param n2 Principal quantum number for orbital 2
/// @param zeta2 Slater exponent for orbital 2
/// @param R_bohr Distance in Bohr
static inline double exactSSOverlap(int n1, double zeta1, int n2, double zeta2, double R_bohr)
{
    if (R_bohr < 1e-10) return 1.0;

    double argA = 0.5 * R_bohr * (zeta1 + zeta2);
    double argB = 0.5 * R_bohr * (zeta1 - zeta2);

    int kmax = n1 + n2 + 2;
    auto Ak = computeAk(argA, kmax);
    auto Bk = computeBk(argB, kmax);

    double ovfac = overlapFactor(n1, 0, 0, zeta1, n2, 0, 0, zeta2, R_bohr);
    auto poly = ssOverlapIndices(n1, n2);

    return ovfac * evaluatePolynomial(poly, Ak, Bk);
}

/// Exact 2-center SP sigma overlap integral (diatomic frame)
/// S orbital on bra (atom A), P orbital on ket (atom B).
/// Multiply result by -direction_cosine for lab frame (Ulysses uses -SProt*sigma for SP).
///
/// @param n_s Principal quantum number of S orbital (bra)
/// @param zeta_s Slater exponent of S orbital
/// @param n_p Principal quantum number of P orbital (ket)
/// @param zeta_p Slater exponent of P orbital
/// @param R_bohr Distance in Bohr
static inline double exactSPSigmaOverlap(int n_s, double zeta_s,
                                           int n_p, double zeta_p, double R_bohr)
{
    if (R_bohr < 1e-10) return 0.0;

    // argA/argB follow Lofthus convention: bra exponent first
    double argA = 0.5 * R_bohr * (zeta_s + zeta_p);
    double argB = 0.5 * R_bohr * (zeta_s - zeta_p);

    int kmax = n_s + n_p + 2;
    auto Ak = computeAk(argA, kmax);
    auto Bk = computeBk(argB, kmax);

    double ovfac = overlapFactor(n_s, 0, 0, zeta_s, n_p, 1, 0, zeta_p, R_bohr);
    auto poly = spSigmaOverlapIndices(n_s, n_p);

    return ovfac * evaluatePolynomial(poly, Ak, Bk);
}

/// Exact 2-center PS sigma overlap integral (diatomic frame)
/// P orbital on bra (atom A), S orbital on ket (atom B).
/// Multiply result by +direction_cosine for lab frame (Ulysses uses +SProt*sigma for PS).
///
/// @param n_p Principal quantum number of P orbital (bra)
/// @param zeta_p Slater exponent of P orbital
/// @param n_s Principal quantum number of S orbital (ket)
/// @param zeta_s Slater exponent of S orbital
/// @param R_bohr Distance in Bohr
static inline double exactPSSigmaOverlap(int n_p, double zeta_p,
                                           int n_s, double zeta_s, double R_bohr)
{
    if (R_bohr < 1e-10) return 0.0;

    // argA/argB: bra (P) exponent first
    double argA = 0.5 * R_bohr * (zeta_p + zeta_s);
    double argB = 0.5 * R_bohr * (zeta_p - zeta_s);

    int kmax = n_p + n_s + 2;
    auto Ak = computeAk(argA, kmax);
    auto Bk = computeBk(argB, kmax);

    double ovfac = overlapFactor(n_p, 1, 0, zeta_p, n_s, 0, 0, zeta_s, R_bohr);
    auto poly = psSigmaOverlapIndices(n_p, n_s);

    return ovfac * evaluatePolynomial(poly, Ak, Bk);
}

/// Exact 2-center PP sigma overlap integral (diatomic frame)
static inline double exactPPSigmaOverlap(int n1, double zeta1,
                                           int n2, double zeta2, double R_bohr)
{
    if (R_bohr < 1e-10) return 1.0;

    double argA = 0.5 * R_bohr * (zeta1 + zeta2);
    double argB = 0.5 * R_bohr * (zeta1 - zeta2);

    int kmax = n1 + n2 + 2;
    auto Ak = computeAk(argA, kmax);
    auto Bk = computeBk(argB, kmax);

    double ovfac = overlapFactor(n1, 1, 0, zeta1, n2, 1, 0, zeta2, R_bohr);
    auto poly = ppSigmaOverlapIndices(n1, n2);

    return ovfac * evaluatePolynomial(poly, Ak, Bk);
}

/// Exact 2-center PP pi overlap integral (diatomic frame)
static inline double exactPPPiOverlap(int n1, double zeta1,
                                       int n2, double zeta2, double R_bohr)
{
    if (R_bohr < 1e-10) return 1.0;

    double argA = 0.5 * R_bohr * (zeta1 + zeta2);
    double argB = 0.5 * R_bohr * (zeta1 - zeta2);

    int kmax = n1 + n2 + 2;
    auto Ak = computeAk(argA, kmax);
    auto Bk = computeBk(argB, kmax);

    double ovfac = overlapFactor(n1, 1, 1, zeta1, n2, 1, 1, zeta2, R_bohr);
    auto poly = ppPiOverlapIndices(n1, n2);

    return ovfac * evaluatePolynomial(poly, Ak, Bk);
}

} // namespace Lofthus
} // namespace STO
