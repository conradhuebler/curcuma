/*
 * <Native QM: analytic closed-shell RHF Hessian (second nuclear derivatives)>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Nov 2026): analytic RHF Hessian with the coupled-perturbed HF (CPHF) equations.
 * Implementation in qm_hessian.cpp; the engine entry point is QMEngine::computeHessian (declared in
 * qm_engine.h). Derivation and validation: docs/QM_HF_HESSIAN.md.
 *
 * This program is free software under GPL-3.0
 */

#pragma once

#include "src/core/global.h"  // Matrix (Eigen, row-major)

#include "GTOIntegrals.hpp"  // GTO::Orbital

#include <string>
#include <vector>

namespace qmhess {

/// Options of the analytic Hessian (filled from the hessian module: -hessian.cphf_threshold, -hessian.cphf_max_iterations).
struct HessianOptions {
    double cphf_threshold = 1.0e-10;  ///< largest |residual| of every CPHF column (Eh) at convergence
    int cphf_max_iterations = 100;    ///< preconditioned conjugate-gradient iterations per column
};

/// What the analytic Hessian did (timings in ms, residuals, J/K builds); `error` is set when it was refused or failed.
struct HessianInfo {
    int cphf_iterations = 0;      ///< largest number of CG iterations over the 3N columns
    double cphf_residual = 0.0;   ///< largest final max|residual| over the columns
    double max_asymmetry = 0.0;   ///< max |H_xy - H_yx| before symmetrisation (Eh/Bohr^2)
    long jk_builds = 0;           ///< J/K builds of response densities (CPHF right-hand sides and iterations)
    double t_1e_ms = 0.0, t_2e_ms = 0.0, t_cphf_ms = 0.0, t_total_ms = 0.0;
    std::string error;
};

/**
 * @brief One-electron derivative integrals of the RHF Hessian (cartesian basis, atomic units).
 *
 * S1[3A+k] = dS/dA_k and H1[3A+k] = d(T+V)/dA_k (basis-function and nuclear (Hellmann-Feynman) parts) as full
 * n x n matrices; hess = sum_mn P_mn d2(T+V)_mn/dXdY - sum_mn W_mn d2S_mn/dXdY (3N x 3N, Eh/Bohr^2), i.e. the
 * one-electron second-derivative integrals contracted with the density P and the energy-weighted density W.
 * Claude Generated (Nov 2026).
 */
void oneElectronDerivatives(const std::vector<GTO::Orbital>& basis, const std::vector<int>& atomZ,
                            const Matrix& atomPosBohr, const Matrix& P, const Matrix& W, std::vector<Matrix>& S1,
                            std::vector<Matrix>& H1, Matrix& hess, int threads = 1);

/**
 * @brief Two-electron derivative integrals of the RHF Hessian (cartesian basis), integral-direct.
 *
 * G1[3A+k]_mn = sum_ls [ (mn|ls)^{A_k} - 1/2 (ml|ns)^{A_k} ] P_ls (the skeleton derivative Fock matrix,
 * density fixed) and hess = 1/2 sum D_mnls (mn|ls)^{XY} with D_mnls = P_mn P_ls - 1/4 (P_ml P_ns + P_ms P_nl).
 * Every canonical shell quartet is computed once (no n^4 tensor is stored); one-centre quartets are skipped.
 * Claude Generated (Nov 2026).
 */
void twoElectronDerivatives(const std::vector<GTO::Orbital>& basis, const Matrix& P, int natoms,
                            std::vector<Matrix>& G1, Matrix& hess, int threads = 1);

/// Validation helper: J(P) - 1/2 K(P) (cartesian) from the same shell-quartet digest as twoElectronDerivatives,
/// with undifferentiated integrals. Must equal the engine's J - K/2. Claude Generated (Nov 2026).
Matrix fockTwoElectronReference(const std::vector<GTO::Orbital>& basis, const Matrix& P, int threads = 1);

/// Second derivatives of E_nn = sum_{i<j} Z_i Z_j / R_ij (3N x 3N, Eh/Bohr^2). Claude Generated (Nov 2026).
Matrix nuclearRepulsionHessian(const std::vector<int>& Z, const Matrix& atomPosBohr);

/// Highest angular momentum of a flat cartesian basis (0 = s, 1 = p, ...).
int maxAngularMomentum(const std::vector<GTO::Orbital>& basis);

}  // namespace qmhess
