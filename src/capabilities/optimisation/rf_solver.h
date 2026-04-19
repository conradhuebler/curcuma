/*
 * <Shared Rational Function Optimization solver>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Extracted from ANCOpt (port of XTB optimizer.f90 by Stefan Grimme).
 * Used by both ANCOptimizer and the native LBFGS/RFO optimizer.
 *
 * Algorithm references:
 *  - Banerjee et al. (1985) J. Phys. Chem. 89, 52  — RFO augmented eigenproblem
 *  - Lanczos (1950) J. Res. Natl. Bur. Stand. 45, 255
 *  - Parlett & Scott (1979) — "twice is enough" reorthogonalization
 *  - XTB: external/xtb/src/optimizer.f90, solver_sdavidson (line 696)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#pragma once

#ifndef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS
#endif

#include <Eigen/Dense>

using RFVector = Eigen::VectorXd;
using RFMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

namespace RFSolver {

/**
 * @brief Lanczos iteration for the lowest eigenpair of the RFO augmented matrix.
 *
 * Computes the lowest eigenpair of the implicit augmented matrix
 *   A_aug = [[H, g], [g^T, 0]]
 * without materializing A_aug. Equivalent to XTB's solver_sdavidson.
 *
 * Full (one-pass) Gram-Schmidt reorthogonalization at each step keeps
 * Ritz values accurate at moderate Lanczos dimensions (m_max <= 120).
 *
 * @param hessian        nvar×nvar symmetric Hessian (ANC or Cartesian)
 * @param gradient       nvar gradient vector
 * @param start_vector   size nvar+1 starting vector (warm-start or steepest-descent)
 * @param out_eigenvector output Ritz vector (size nvar+1), normalized
 * @param out_eigenvalue  output Ritz eigenvalue
 * @param m_max          maximum Lanczos steps (default 80)
 * @param tol            convergence tolerance relative to ||A_aug||_F (default 1e-6)
 * @return true if converged, false otherwise
 */
bool lanczosLowestEigenpair(const RFMatrix& hessian,
                            const RFVector& gradient,
                            const RFVector& start_vector,
                            RFVector& out_eigenvector,
                            double& out_eigenvalue,
                            int m_max = 80,
                            double tol = 1e-6);

/**
 * @brief Rational Function Optimization step via augmented eigenvalue problem.
 *
 * Solves the RFO augmented system:
 *   | H  g | * | p | = lambda * | p |
 *   | g' 0 |   | s |            | s |
 * and returns p/s as the displacement vector.
 *
 * For nvar+1 >= 50: Lanczos with warm-start from warm_start (in/out).
 * Fallback for small systems or Lanczos failure: full SelfAdjointEigenSolver.
 *
 * @param gradient   nvar gradient vector
 * @param hessian    nvar×nvar symmetric Hessian
 * @param warm_start In/out warm-start vector (size nvar+1); resized to zero if stale
 * @return RFO displacement vector (size nvar)
 */
RFVector calculateRFStep(const RFVector& gradient,
                         const RFMatrix& hessian,
                         RFVector& warm_start);

} // namespace RFSolver
