/*
 * <Native QM: cartesian -> spherical transform for f and g shells>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * The columns of the spherical-transform matrix Q (qmint::buildSphericalTransform) that belong to one f (l = 3) or
 * g (l = 4) shell: 2l+1 real solid harmonics over the (l+1)(l+2)/2 cartesian components, orthonormalised in the
 * overlap metric of the shell block (as the d shells in qm_integrals.cpp).
 *
 * Own translation unit: qm_integrals.cpp already hits GCC's inline-unit-growth limit (9 failed inlines of
 * GTO::orbitalTypeToComponents before this change), so every line added there shifts the inlining of its hot kernels.
 *
 * Literature: H. B. Schlegel, M. J. Frisch, Int. J. Quantum Chem. 54, 83 (1995), eq. 15;
 *             T. Helgaker, P. Jorgensen, J. Olsen, Molecular Electronic-Structure Theory (2000), eq. 6.4.47-6.4.50.
 *
 * Claude Generated (Nov 2026). AI-generated, machine-tested. GPL-3.0.
 */

#include "qm_integrals.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace qmint {

namespace {

double binom(int n, int k)
{
    if (k < 0 || k > n) return 0.0;
    double r = 1.0;
    for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

/// @brief Real solid harmonic S_lm as a polynomial over the cartesian monomials x^a y^b z^c of degree l,
/// coefficients in the lexicographic component order (GTO::lexicographicComponent). Unnormalised (the caller
/// normalises in the overlap metric):
///   S_lm = sum_{t=0}^{(l-|m|)/2} sum_{u=0}^{t} sum_{k} C_tuk x^(2t+|m|-2u-k) y^(2u+k) z^(l-2t-|m|),
///   C_tuk = (-1)^(t+(k-k0)/2) 4^-t binom(l,t) binom(l-t,|m|+t) binom(t,u) binom(|m|,k),
/// with k = 2v over the even (m >= 0, k0 = 0) or odd (m < 0, k0 = 1) integers 0 <= k <= |m|.
std::vector<double> solidHarmonic(int l, int m)
{
    const int ncomp = (l + 1) * (l + 2) / 2;
    std::vector<double> v(ncomp, 0.0);
    auto index = [l](int a, int b) {  // lexicographic: a descending, then b descending
        int idx = 0;
        for (int x = l; x > a; --x) idx += l - x + 1;
        return idx + (l - a - b);
    };
    const int am = std::abs(m), k0 = m < 0 ? 1 : 0;
    for (int t = 0; t <= (l - am) / 2; ++t)
        for (int u = 0; u <= t; ++u)
            for (int k = k0; k <= am; k += 2) {
                const int a = 2 * t + am - 2 * u - k, b = 2 * u + k, c = l - 2 * t - am;
                if (a < 0 || c < 0) continue;
                const double sign = (((t + (k - k0) / 2) & 1) ? -1.0 : 1.0);
                v[index(a, b)] += sign * std::pow(0.25, t) * binom(l, t) * binom(l - t, am + t) * binom(t, u)
                    * binom(am, k);
            }
    return v;
}

}  // namespace

void higherShellSphericalColumns(const std::vector<GTO::Orbital>& basis, const Matrix& Scart, int base, int L,
                                 Matrix& Q, int col)
{
    // A solid harmonic h = sum_k v_k x^a y^b z^c R(r); AO component k is c_k x^a y^b z^c R(r) with its own constant
    // c_k (primitive norm and S_ii = 1, same radial factor R for all components), so the coefficient of AO k is
    // v_k / c_k; the ratios c_k / c_0 are those of the first contraction coefficients. Then Gram-Schmidt in the S
    // metric of the block, as for d (different m are orthogonal already, so this only normalises).
    const int n = static_cast<int>(basis.size());
    const int nc = (L + 1) * (L + 2) / 2;
    if (base + nc > n) throw std::runtime_error("buildSphericalTransform: truncated f/g shell");
    for (int a = 0; a < nc; ++a) {
        int l, m, k;
        GTO::orbitalTypeToComponents(basis[base + a].type, l, m, k);
        int la, ma, ka;
        GTO::lexicographicComponent(L, a, la, ma, ka);
        if (l != la || m != ma || k != ka || basis[base + a].exponents != basis[base].exponents
            || basis[base + a].atom != basis[base].atom)
            throw std::runtime_error("buildSphericalTransform: f/g shell components not in the parser order");
    }
    Eigen::MatrixXd Sd(nc, nc);
    for (int a = 0; a < nc; ++a)
        for (int b = 0; b < nc; ++b) Sd(a, b) = Scart(base + a, base + b);
    const double c0 = basis[base].coefficients.empty() ? 1.0 : basis[base].coefficients[0];
    // ORCA order of the spherical components: m = 0, +1, -1, +2, -2, ... (as the d order of buildSphericalTransform)
    std::vector<int> mlist = { 0 };
    for (int m = 1; m <= L; ++m) {
        mlist.push_back(m);
        mlist.push_back(-m);
    }
    for (int kk = 0; kk < 2 * L + 1; ++kk) {
        const std::vector<double> h = solidHarmonic(L, mlist[kk]);
        Eigen::VectorXd v(nc);
        for (int a = 0; a < nc; ++a) {
            const double ca = basis[base + a].coefficients.empty() ? 1.0 : basis[base + a].coefficients[0];
            v(a) = h[a] * (c0 / ca);
        }
        for (int p = 0; p < kk; ++p) {
            Eigen::VectorXd uprev(nc);
            for (int a = 0; a < nc; ++a) uprev(a) = Q(base + a, col + p);
            const double sprod = uprev.dot(Sd * v);
            const double norm2 = uprev.dot(Sd * uprev);
            if (std::abs(norm2) > 1e-15) v -= (sprod / norm2) * uprev;
        }
        const double norm2 = v.dot(Sd * v);
        const double inv = norm2 < 1e-18 ? 0.0 : 1.0 / std::sqrt(norm2);
        for (int a = 0; a < nc; ++a) Q(base + a, col + kk) = v(a) * inv;
    }
}

}  // namespace qmint
