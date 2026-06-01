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
 *   1. Householder tridiagonalization  A = Q·T·Qᵀ        (tred2)
 *   2. Symmetric tridiagonal eigensolve  T = V·Λ·Vᵀ      (Cuppen divide-and-conquer)
 *   3. Back-transform  eigenvectors = Q·V
 *
 * Step 2 is the Cuppen divide-and-conquer (Cuppen 1981; Dongarra–Sorensen 1987;
 * Gu–Eisenstat 1994): recursively tear the tridiagonal into two halves plus a rank-1
 * coupling, solve each half, then merge by solving the rank-1-updated eigenproblem —
 * roots of the secular equation found by shifted bisection (dlaed4-style, accurate
 * denominators), with deflation (negligible weights + degenerate diagonals via Givens)
 * and Löwner/Gu–Eisenstat-reconstructed weights for orthogonal eigenvectors. Small blocks
 * (n ≤ 32) fall back to implicit-shift QL (tql2). The recursion's sub-problems are
 * independent (a natural parallel/GPU split); the current implementation is serial but
 * already on par with MKL dsyevd on the GFN workload (the surrounding BLAS threads).
 *
 * Accuracy: eigenpairs match LAPACK/Eigen to ~1e-12 (eigenvalues) and ~1e-14
 * (reconstruction, orthonormality) over random, degenerate, identity-like and tight-cluster
 * spectra up to n=558 (test_native_eigensolver). Eigenvalues ascending, Cᵀ·C = I.
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
