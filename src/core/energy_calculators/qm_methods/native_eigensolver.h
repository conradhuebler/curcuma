/*
 * <Native symmetric eigensolver — self-contained alternative to MKL dsyevd>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated. GPL-3.0.
 *
 * A from-scratch, dependency-free (no MKL/LAPACK) eigensolver for real symmetric
 * matrices, selectable as an alternative to the MKL `dsyevd` path in the native xTB
 * SCF (`-eigensolver mkl|native`). MKL stays the default; this exists so the SCF can
 * run without a proprietary BLAS and as the CPU foundation for a future GPU port.
 *
 * Pipeline (the standard symmetric eigensolve):
 *   1. Householder tridiagonalization  A = Q·T·Qᵀ        (tridiagonalizeHouseholder)
 *   2. Symmetric tridiagonal eigensolve  T = V·Λ·Vᵀ      (tridiagonalEigen)
 *   3. Back-transform  eigenvectors = Q·V
 *
 * Step 2 currently uses the implicit-shift QL algorithm (tql2), which is correct and
 * self-contained. The Cuppen *divide-and-conquer* tridiagonal solver — the performance/
 * GPU-parallel core — is the planned drop-in replacement for tridiagonalEigen (same
 * interface); see docs/SQM_THREADING_WP.md WP4. Until it lands and is validated to MKL
 * accuracy, `native` is QL-based.
 *
 * Accuracy target: eigenpairs match LAPACK dsyevd to ~1e-10 (validated by
 * test_native_eigensolver). Eigenvalues are returned in ascending order with
 * Cᵀ·C = I (orthonormal eigenvectors).
 */

#pragma once

#include <Eigen/Dense>

namespace curcuma::eigsolver {

/**
 * @brief Eigenvalues (ascending) and eigenvectors of a real symmetric matrix.
 *
 * Self-contained (no MKL/LAPACK). On success fills `evals` (length n, ascending) and
 * `evecs` (n×n, columns = orthonormal eigenvectors) so that A ≈ evecs·diag(evals)·evecsᵀ.
 * Only the lower triangle of A is read (A is assumed symmetric). Returns false on a
 * non-square/empty matrix or a convergence failure.
 *
 * @param A     symmetric input matrix (n×n), read-only
 * @param evals output eigenvalues, ascending
 * @param evecs output eigenvectors (columns)
 */
bool solveSymmetric(const Eigen::MatrixXd& A,
                    Eigen::VectorXd& evals,
                    Eigen::MatrixXd& evecs);

} // namespace curcuma::eigsolver
