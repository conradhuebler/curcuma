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

#include "src/core/global.h"   // project Matrix (RowMajor) / Vector

#include <Eigen/Dense>

#include <functional>

namespace curcuma::eigsolver {

/**
 * @brief Eigenvalues (ascending) and eigenvectors of a real symmetric matrix.
 *
 * Self-contained (no MKL/LAPACK). On success fills `evals` (length n, ascending) and
 * `evecs` (n×n, columns = orthonormal eigenvectors) so that A ≈ evecs·diag(evals)·evecsᵀ.
 * Only the lower triangle of A is read (A is assumed symmetric). Returns false on a
 * non-square/empty matrix or a convergence failure.
 *
 * @param A        symmetric input matrix (n×n), read-only
 * @param evals    output eigenvalues, ascending
 * @param evecs    output eigenvectors (columns)
 * @param nThreads intra-eigensolver thread budget (default 1 = serial). When >1 and the
 *                 matrix is large enough, the divide-and-conquer's independent subtrees are
 *                 solved in parallel (the merges' BLAS gemms thread separately). Caller must
 *                 gate this (the native xTB passes its effectiveIntraThreads).
 */
bool solveSymmetric(const Eigen::MatrixXd& A,
                    Eigen::VectorXd& evals,
                    Eigen::MatrixXd& evecs,
                    int nThreads = 1);

/**
 * @brief 0 K density-matrix purification (TC2) — GEMM/trace only, no diagonalization.
 *
 * Builds the idempotent density projector `Ptil` onto the `nocc` lowest eigenstates of the
 * symmetric (orthonormal-basis) matrix `Atil`, plus its energy-weighted form `Wtil = Ptil·Atil·Ptil`,
 * using only matrix products and traces (Niklasson trace-correcting 2nd-order purification,
 * 2002). This is the GPU-portable density path: no eigenvalues/eigenvectors are formed.
 *
 * Intended for a closed-shell system with a HOMO–LUMO gap and integer occupation (`nocc` = number
 * of doubly-occupied orbitals). `Ptil` has eigenvalues 0/1 with `Tr(Ptil)=nocc`. Returns false if
 * it fails to reach an idempotent fixed point (e.g. gapless/metallic), so the caller can fall back
 * to a diagonalization.
 *
 * @param Atil    symmetric input (orthonormal-basis Fock Ã = L⁻¹·F·L⁻ᵀ), n×n
 * @param nocc    number of occupied (doubly-occupied) orbitals, 0 < nocc < n
 * @param Ptil    output idempotent density projector (n×n, Tr = nocc)
 * @param Wtil    output energy-weighted density Ptil·Atil·Ptil (n×n)
 * @param maxIter purification iteration cap (default 200)
 * @param tol     idempotency tolerance on Tr(P) − Tr(P²) (default 1e-10)
 */
bool purifyDensity(const Eigen::MatrixXd& Atil, int nocc,
                   Eigen::MatrixXd& Ptil, Eigen::MatrixXd& Wtil,
                   int maxIter = 200, double tol = 1e-10);

/**
 * @brief Lowest-`k` eigenpairs of a symmetric matrix by seeded block LOBPCG (iterative, GEMM-based).
 *
 * Locally Optimal Block Preconditioned Conjugate Gradient (Knyazev 2001) for the `k` smallest
 * eigenvalues/vectors of symmetric `A`, with a diagonal preconditioner and a robust generalized
 * Rayleigh-Ritz (Gram-matrix eigen-filtering of the [X|W|P] subspace, so near-dependent search
 * directions are dropped rather than corrupting the solve). The heavy work is GEMMs (A·X, the
 * subspace projections); only the small 3k×3k Ritz problems use a dense solver — the GPU-portable
 * iterative path, the partner of a future sparse/fragmented A.
 *
 * **Seeding:** if `X` enters as a valid n×k block (e.g. the previous SCF iteration's lowest-k
 * vectors) it is used as the initial guess; otherwise a deterministic random block is used. On
 * success `X` holds the k lowest eigenvectors (orthonormal columns, ascending) and `evals` the k
 * lowest eigenvalues — so the caller can store `X` to seed the next solve.
 *
 * **Caveat:** LOBPCG only pays when k ≪ n or A is sparse; for a dense matrix at k ≈ n/2 (a GFN
 * minimal basis at ~50% occupancy) it costs more than a full dsyevd — this is the opt-in
 * `-eigensolver lobpcg` research/experimental path, not a CPU win. Returns false if it fails to
 * converge to `tol` within `maxIter` (so the caller can fall back to a dense solve).
 *
 * **Guard vectors:** the block size is `k`, but only the lowest `nConverge` states must reach
 * `tol`; the extra `k − nConverge` columns are guard vectors that accelerate the wanted states
 * (the boundary eigenvector of a gapless block converges slowly, so request a few guards). Default
 * `nConverge = -1` means all `k`.
 *
 * @param A         symmetric input (orthonormal-basis Fock Ã), n×n
 * @param k         block size = number of eigenpairs returned, 0 < k ≤ n
 * @param X         in: optional n×k seed (else random); out: k lowest eigenvectors (columns)
 * @param evals     out: k lowest eigenvalues, ascending
 * @param maxIter   LOBPCG iteration cap (default 50)
 * @param tol       relative residual convergence tolerance (default 1e-8)
 * @param nConverge number of lowest states that must converge (default -1 = all k); the rest guard
 */
bool lobpcgLowest(const Eigen::MatrixXd& A, int k,
                  Eigen::MatrixXd& X, Eigen::VectorXd& evals,
                  int maxIter = 50, double tol = 1e-8, int nConverge = -1);

/* ===================================================================== *
 *  Raw-pointer building blocks — the pieces a GPU backend drives itself.
 *
 *  The Vulkan xTB engine (qm_methods/vulkan/xtb_vulkan_context.cpp) runs the
 *  tridiagonalization and the divide-and-conquer merges on the device, but the
 *  HOST half of each step is exactly the code above. These two entry points
 *  expose it on raw column-major arrays (what a device buffer maps to), so the
 *  backend calls the canonical implementation instead of keeping a copy.
 *  Claude Generated (Sep 2026).
 * ===================================================================== */

/**
 * @brief Symmetric tridiagonal eigensolve (implicit-shift QL) on raw arrays.
 *
 * The raw-pointer face of the D&C base case: `diag[n]` / `off[n-1]` (off[k] connects rows
 * k and k+1) in, ascending `eval[n]` and the matching eigenvectors `evec` (n×n,
 * COLUMN-major, column j ↔ eval[j]) out. Full rotation accumulation, so degenerate
 * spectra are handled without inverse-iteration clustering. Returns false if an
 * eigenvalue needs > 50 QL iterations.
 *
 * @param n     matrix dimension (> 0)
 * @param diag  diagonal, length n
 * @param off   sub-diagonal, length n-1 (unread for n == 1)
 * @param eval  out: eigenvalues, ascending, length n
 * @param evec  out: eigenvectors, column-major n*n
 */
bool solveTridiagonalQL(int n, const double* diag, const double* off,
                        double* eval, double* evec);

/**
 * @brief Secular-equation solver hook for rank1EigenDeflate.
 *
 * Called once per merge with the DEFLATED secular set: `D[k]` (strictly ascending),
 * `z[k]` (weights), `rho` > 0. Must write the k ascending roots to `lam[k]` and the
 * secular eigenvectors to `Wsec` (k×k, column-major, stride k). Return false to abort
 * the merge. An empty hook means "use the built-in host solver".
 */
using SecularSolveFn = std::function<bool(const double* D, const double* z, int k,
                                          double rho, double* lam, double* Wsec)>;

/**
 * @brief Eigenproblem of diag(D) + rho·z·zᵀ with deflation, on raw arrays.
 *
 * The rank-1 merge step of the Cuppen divide-and-conquer: normalise z, fold the sign of
 * rho, sort D, deflate (Givens rotations for degenerate diagonals + negligible weights),
 * solve the secular equation on what is left, assemble, undo the rotations and scatter
 * back to the input basis. Returns `eval[n]` and `W` (n×n COLUMN-major, eval[j] ↔ column j).
 *
 * Pass `secular` to run only the secular solve elsewhere (e.g. a GPU kernel); the
 * deflation bookkeeping around it — the part that is easy to get subtly wrong — stays
 * shared. With an empty hook this is the plain host solver.
 *
 * @param D_in     diagonal, length n
 * @param z_in     rank-1 update vector, length n
 * @param n        dimension (> 0)
 * @param rho_in   rank-1 coefficient (either sign)
 * @param eval     out: eigenvalues, length n (input order, NOT sorted)
 * @param W        out: eigenvectors, column-major n*n
 * @param secular  optional secular-solve hook (default: the built-in host solver)
 * @param nThreads thread budget for the built-in secular solve (ignored when hooked)
 */
bool rank1EigenDeflate(const double* D_in, const double* z_in, int n, double rho_in,
                       double* eval, double* W,
                       const SecularSolveFn& secular = SecularSolveFn(),
                       int nThreads = 1);

} // namespace curcuma::eigsolver
