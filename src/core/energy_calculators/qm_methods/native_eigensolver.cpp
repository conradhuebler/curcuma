/*
 * <Native symmetric eigensolver — implementation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated. GPL-3.0.
 *
 * Householder tridiagonalization (tred2) + implicit-shift QL with eigenvector
 * accumulation (tql2). Classic, dependency-free symmetric eigensolver (EISPACK /
 * Numerical Recipes lineage), 0-indexed and expressed on Eigen storage. Used as the
 * `-eigensolver native` alternative to MKL dsyevd; see native_eigensolver.h. The
 * tridiagonal QL step is the drop-in point for a future Cuppen divide-and-conquer
 * solver (same d/e in, eigenpairs out).
 */

#include "native_eigensolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace curcuma::eigsolver {

namespace {

// sqrt(a^2 + b^2) without destructive under/overflow (Numerical Recipes pythag).
inline double pythag(double a, double b)
{
    const double aa = std::fabs(a), ab = std::fabs(b);
    if (aa > ab) { const double r = ab / aa; return aa * std::sqrt(1.0 + r * r); }
    if (ab == 0.0) return 0.0;
    const double r = aa / ab; return ab * std::sqrt(1.0 + r * r);
}

inline double sign(double a, double b) { return (b >= 0.0) ? std::fabs(a) : -std::fabs(a); }

/* ------------------------------------------------------------------ *
 *  tred2 — Householder reduction of a real symmetric matrix to       *
 *  tridiagonal form. On entry `z` holds the symmetric matrix; on     *
 *  exit `z` is the orthogonal transform Q (A = Q·T·Qᵀ), `d` the       *
 *  diagonal of T, `e` its sub-diagonal (e[0]=0).                     *
 * ------------------------------------------------------------------ */
void tred2(Eigen::MatrixXd& z, Eigen::VectorXd& d, Eigen::VectorXd& e)
{
    const int n = static_cast<int>(z.rows());
    for (int i = n - 1; i >= 1; --i) {
        const int l = i - 1;
        double h = 0.0, scale = 0.0;
        if (l > 0) {
            for (int k = 0; k <= l; ++k) scale += std::fabs(z(i, k));
            if (scale == 0.0) {
                e(i) = z(i, l);
            } else {
                for (int k = 0; k <= l; ++k) { z(i, k) /= scale; h += z(i, k) * z(i, k); }
                double f = z(i, l);
                const double g = (f >= 0.0) ? -std::sqrt(h) : std::sqrt(h);
                e(i) = scale * g;
                h -= f * g;
                z(i, l) = f - g;
                f = 0.0;
                for (int j = 0; j <= l; ++j) {
                    z(j, i) = z(i, j) / h;
                    double gg = 0.0;
                    for (int k = 0; k <= j; ++k)     gg += z(j, k) * z(i, k);
                    for (int k = j + 1; k <= l; ++k) gg += z(k, j) * z(i, k);
                    e(j) = gg / h;
                    f += e(j) * z(i, j);
                }
                const double hh = f / (h + h);
                for (int j = 0; j <= l; ++j) {
                    f = z(i, j);
                    const double gg = e(j) - hh * f;
                    e(j) = gg;
                    for (int k = 0; k <= j; ++k)
                        z(j, k) -= (f * e(k) + gg * z(i, k));
                }
            }
        } else {
            e(i) = z(i, l);
        }
        d(i) = h;
    }
    d(0) = 0.0;
    e(0) = 0.0;
    // Accumulate the transformation matrices.
    for (int i = 0; i < n; ++i) {
        const int l = i - 1;
        if (d(i) != 0.0) {
            for (int j = 0; j <= l; ++j) {
                double g = 0.0;
                for (int k = 0; k <= l; ++k) g += z(i, k) * z(k, j);
                for (int k = 0; k <= l; ++k) z(k, j) -= g * z(k, i);
            }
        }
        d(i) = z(i, i);
        z(i, i) = 1.0;
        for (int j = 0; j <= l; ++j) { z(j, i) = 0.0; z(i, j) = 0.0; }
    }
}

/* ------------------------------------------------------------------ *
 *  tql2 — eigenvalues/vectors of a symmetric tridiagonal matrix by   *
 *  the QL algorithm with implicit shifts. `d` holds the diagonal     *
 *  (→ eigenvalues), `e` the sub-diagonal, `z` the accumulated        *
 *  transform from tred2 (→ eigenvectors, columns). Returns false if  *
 *  an eigenvalue needs > 50 iterations to isolate.                   *
 * ------------------------------------------------------------------ */
bool tql2(Eigen::MatrixXd& z, Eigen::VectorXd& d, Eigen::VectorXd& e)
{
    const int n = static_cast<int>(z.rows());
    for (int i = 1; i < n; ++i) e(i - 1) = e(i);
    e(n - 1) = 0.0;

    for (int l = 0; l < n; ++l) {
        int iter = 0;
        int m = 0;
        do {
            for (m = l; m < n - 1; ++m) {
                const double dd = std::fabs(d(m)) + std::fabs(d(m + 1));
                if (std::fabs(e(m)) <= std::numeric_limits<double>::epsilon() * dd) break;
            }
            if (m != l) {
                if (iter++ == 50) return false;
                double g = (d(l + 1) - d(l)) / (2.0 * e(l));
                double r = pythag(g, 1.0);
                g = d(m) - d(l) + e(l) / (g + sign(r, g));
                double s = 1.0, c = 1.0, p = 0.0;
                int i = m - 1;
                for (; i >= l; --i) {
                    double f = s * e(i);
                    const double b = c * e(i);
                    r = pythag(f, g);
                    e(i + 1) = r;
                    if (r == 0.0) { d(i + 1) -= p; e(m) = 0.0; break; }
                    s = f / r;
                    c = g / r;
                    g = d(i + 1) - p;
                    r = (d(i) - g) * s + 2.0 * c * b;
                    p = s * r;
                    d(i + 1) = g + p;
                    g = c * r - b;
                    for (int k = 0; k < n; ++k) {
                        f = z(k, i + 1);
                        z(k, i + 1) = s * z(k, i) + c * f;
                        z(k, i)     = c * z(k, i) - s * f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d(l) -= p;
                e(l) = g;
                e(m) = 0.0;
            }
        } while (m != l);
    }
    return true;
}

} // namespace

bool solveSymmetric(const Eigen::MatrixXd& A,
                    Eigen::VectorXd& evals,
                    Eigen::MatrixXd& evecs)
{
    const int n = static_cast<int>(A.rows());
    if (n == 0 || A.cols() != n) return false;

    // Work on a symmetric copy (tred2 reads/writes in place; symmetrise to be safe
    // against tiny asymmetries in the input).
    Eigen::MatrixXd z = 0.5 * (A + A.transpose());
    Eigen::VectorXd d(n), e(n);

    tred2(z, d, e);
    if (!tql2(z, d, e)) return false;

    // Sort eigenpairs ascending (QL leaves them unordered).
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return d(a) < d(b); });

    evals.resize(n);
    evecs.resize(n, n);
    for (int j = 0; j < n; ++j) {
        evals(j) = d(idx[j]);
        evecs.col(j) = z.col(idx[j]);
    }
    return true;
}

} // namespace curcuma::eigsolver
