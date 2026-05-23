/*
 * <Shared Rational Function Optimization solver>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * See rf_solver.h for algorithm references and API documentation.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#ifndef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS
#endif
// EIGEN_USE_LAPACKE: enables LAPACK dsyevd (divide-and-conquer) inside
// SelfAdjointEigenSolver. Safe here — no variable named 'I' in this file.
#ifndef EIGEN_USE_LAPACKE
#define EIGEN_USE_LAPACKE
#endif

#include "rf_solver.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace RFSolver {

bool lanczosLowestEigenpair(const RFMatrix& hessian,
                            const RFVector& gradient,
                            const RFVector& start_vector,
                            RFVector& out_eigenvector,
                            double& out_eigenvalue,
                            int m_max,
                            double tol)
{
    // Symmetric Lanczos on the implicit augmented matrix
    //   A_aug = [[H, g], [g^T, 0]]
    // with full reorthogonalization. Never materializes A_aug.
    //
    // Implicit matvec:
    //   w.head(nvar) = H * v.head(nvar) + g * v(nvar)
    //   w(nvar)      = g . v.head(nvar)

    const int nvar = static_cast<int>(gradient.size());
    const int N = nvar + 1;
    if (hessian.rows() != nvar || hessian.cols() != nvar) return false;
    if (start_vector.size() != N || N == 0) return false;

    double v0_norm = start_vector.norm();
    if (v0_norm < 1e-14) return false;

    auto apply_A = [&](const Eigen::Ref<const RFVector>& v, RFVector& w) {
        w.resize(N);
        w.head(nvar).noalias() = hessian * v.head(nvar);
        w.head(nvar).noalias() += gradient * v(nvar);
        w(nvar) = gradient.dot(v.head(nvar));
    };

    m_max = std::min(m_max, N);

    // Scale for residual threshold — computed once, O(nvar^2).
    const double h_norm_sq = hessian.squaredNorm();
    const double g_norm_sq = gradient.squaredNorm();
    const double matrix_scale = std::sqrt(h_norm_sq + 2.0 * g_norm_sq) + 1e-12;

    RFMatrix V(N, m_max);
    std::vector<double> alpha; alpha.reserve(m_max);
    std::vector<double> beta;  beta.reserve(m_max);

    V.col(0) = start_vector / v0_norm;

    RFVector w;
    RFVector v_prev = RFVector::Zero(N);
    double beta_prev = 0.0;
    double prev_ritz = std::numeric_limits<double>::max();

    for (int j = 0; j < m_max; ++j) {
        apply_A(V.col(j), w);
        if (j > 0) w.noalias() -= beta_prev * v_prev;

        double a_j = V.col(j).dot(w);
        alpha.push_back(a_j);
        w.noalias() -= a_j * V.col(j);

        // One-pass Gram-Schmidt reorthogonalization (Parlett & Scott: one pass
        // is sufficient for double precision at these Lanczos dimensions).
        if (j > 0) {
            RFVector coeffs = V.leftCols(j + 1).transpose() * w;
            w.noalias() -= V.leftCols(j + 1) * coeffs;
        }

        double b_j = w.norm();
        beta.push_back(b_j);

        const bool breakdown = (b_j < 1e-12);

        // Check convergence via the tridiagonal Ritz value (O(j^3), tiny for j < 120).
        if (j + 1 >= 3 && (j % 3 == 2 || breakdown || j == m_max - 1)) {
            const int m = j + 1;
            RFMatrix T = RFMatrix::Zero(m, m);
            for (int k = 0; k < m; ++k) {
                T(k, k) = alpha[k];
                if (k + 1 < m) {
                    T(k, k + 1) = beta[k];
                    T(k + 1, k) = beta[k];
                }
            }
            Eigen::SelfAdjointEigenSolver<RFMatrix> tri_solver(T);
            if (tri_solver.info() != Eigen::Success) return false;

            double ritz = tri_solver.eigenvalues()(0);
            RFVector y = tri_solver.eigenvectors().col(0);

            double residual = std::abs(b_j * y(m - 1));
            const double scaled_residual_thr = std::max(tol * matrix_scale, tol);

            // Guard against false early convergence on sparse Krylov basis.
            const bool ritz_stable = (j + 1 >= 5) &&
                                     std::abs(ritz - prev_ritz) < tol;
            prev_ritz = ritz;

            if (residual < scaled_residual_thr || breakdown || ritz_stable) {
                out_eigenvector = V.leftCols(m) * y;
                out_eigenvector.normalize();
                out_eigenvalue = ritz;
                return true;
            }
        }

        if (breakdown) return false;

        v_prev = V.col(j);
        beta_prev = b_j;
        if (j + 1 < m_max) V.col(j + 1) = w / b_j;
    }

    return false;
}

RFVector calculateRFStep(const RFVector& gradient,
                         const RFMatrix& hessian,
                         RFVector& warm_start)
{
    const int nvar = static_cast<int>(gradient.size());
    const int nvar1 = nvar + 1;

    RFVector eigenvec;
    double eigenvalue = 0.0;
    bool lanczos_ok = false;

    // Lanczos path for nvar1 >= 50, matching XTB optimizer.f90:687.
    if (nvar1 >= 50) {
        const bool have_warm = (warm_start.size() == nvar1);
        RFVector start(nvar1);

        if (have_warm) {
            start = warm_start;
        } else {
            // Steepest-descent initial guess (XTB:691-694).
            start.head(nvar) = -gradient;
            start(nvar) = 1.0;
            start.normalize();
        }

        lanczos_ok = lanczosLowestEigenpair(hessian, gradient, start,
                                            eigenvec, eigenvalue);

        // Retry from steepest-descent if warm-start failed.
        if (!lanczos_ok && have_warm) {
            RFVector sd(nvar1);
            sd.head(nvar) = -gradient;
            sd(nvar) = 1.0;
            sd.normalize();
            lanczos_ok = lanczosLowestEigenpair(hessian, gradient, sd,
                                                eigenvec, eigenvalue);
        }
    }

    // Fallback: full SelfAdjointEigenSolver — always exact, O(N^3).
    if (!lanczos_ok) {
        RFMatrix A_aug = RFMatrix::Zero(nvar1, nvar1);
        A_aug.topLeftCorner(nvar, nvar) = hessian;
        A_aug.block(nvar, 0, 1, nvar) = gradient.transpose();
        A_aug.block(0, nvar, nvar, 1) = gradient;

        Eigen::SelfAdjointEigenSolver<RFMatrix> solver(A_aug);
        if (solver.info() != Eigen::Success) {
            // Hard fallback: steepest descent with small step.
            warm_start.resize(0);
            return -gradient * 0.1;
        }
        eigenvec = solver.eigenvectors().col(0);
        eigenvalue = solver.eigenvalues()(0);
    }

    // Ensure consistent sign convention: last component must be positive
    // so that dividing by it yields a minimizing direction.
    if (eigenvec(nvar) < 0.0) eigenvec = -eigenvec;

    if (std::abs(eigenvec(nvar)) < 1e-10) {
        warm_start.resize(0);
        return -gradient * 0.1;
    }

    // Cache Ritz vector as warm-start for the next call.
    warm_start = eigenvec;

    return eigenvec.head(nvar) / eigenvec(nvar);
}

} // namespace RFSolver
