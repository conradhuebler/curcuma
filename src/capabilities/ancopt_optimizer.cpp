/*
 * <AncOpt Optimizer Implementation>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Port of AncOpt from XTB by Stefan Grimme
 * Original: external/xtb/src/optimizer.f90
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

// Route Eigen dense matrix ops (GEMV/GEMM) to BLAS if available. This is the
// same pattern used in eeq_solver.cpp. Must be defined before any Eigen header.
// For nvar~700 the H*v matvec in the Lanczos inner loop dominates the RF step;
// BLAS GEMV is roughly 5x faster than Eigen's generic path on typical hardware.
#ifndef EIGEN_USE_BLAS
#define EIGEN_USE_BLAS
#endif

#include "ancopt_optimizer.h"
#include "optimisation/rf_solver.h"
#include "src/core/curcuma_logger.h"
#include "src/core/molecule.h"
#include "src/core/parameter_registry.h"
#include "src/tools/geometry.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace Optimization {

// ============================================================================
// ANCCoordinates Implementation
// ============================================================================

void ANCCoordinates::allocate(int natoms, int num_vars, double h_low, double h_max) {
    n_atoms = natoms;
    n3 = 3 * natoms;
    nvar = num_vars;
    hlow = h_low;
    hmax = h_max;

    B = Matrix::Zero(n3, nvar);
    coord = Vector::Zero(nvar);
    coord_reference = Vector::Zero(nvar); // Claude Nov 2025: Initialize coord reference
    xyz_flat = Vector::Zero(n3);
    xyz_reference = Vector::Zero(n3); // Claude Nov 2025: Initialize xyz reference
    hess = Matrix::Zero(nvar, nvar);

    initialized = true;
}

void ANCCoordinates::deallocate() {
    initialized = false;
}

bool ANCCoordinates::generateANC(const Matrix& cartesian_hessian, const Vector& xyz, bool is_linear) {
    if (!initialized) {
        CurcumaLogger::error("ANC not initialized before generation");
        return false;
    }

    // Tier-L+ dispatch: for n3 > 1800 (~600 atoms), the O(n3^3) full eigendecomp
    // becomes the dominant runtime cost. Switch to a truncated Lanczos path that
    // computes only the top-nvar eigenpairs in O(n3^2 * m) time.
    if (n3 > 1800) {
        return generateANCLanczos(cartesian_hessian, xyz, is_linear, nvar);
    }

    // Diagonalize Hessian to get eigenvectors (approximate normal modes)
    Eigen::SelfAdjointEigenSolver<Matrix> solver(cartesian_hessian);
    if (solver.info() != Eigen::Success) {
        CurcumaLogger::error("Failed to diagonalize Hessian for ANC generation");
        return false;
    }

    // Bug 4 Fix: Apply XTB eigenvalue damping before selecting eigenvectors
    // Shifts all eigenvalues so the minimum is at least hlow
    // Reference: XTB type/anc.f90 generate_anc_blowup() Z. 94-108
    Vector evals = solver.eigenvalues();
    const Matrix& eigenvectors = solver.eigenvectors();

    double elow = std::numeric_limits<double>::max();
    for (int i = 0; i < evals.size(); ++i)
        if (std::abs(evals(i)) > 1e-10)
            elow = std::min(elow, evals(i));

    if (elow < std::numeric_limits<double>::max()) {
        double damp = std::max(hlow - elow, 0.0);
        if (damp > 0.0) {
            for (int i = 0; i < evals.size(); ++i)
                if (std::abs(evals(i)) > 1e-11)
                    evals(i) += damp;
        }
    }

    // Bug 4 Fix (Apr 2026): Apply detrotra8 — physics-based T/R mode identification.
    // XTB calls detrotra8() AFTER diagonalization (type/anc.f90 Zeile 167).
    // It distorts the geometry along each low-frequency mode and measures how much
    // interatomic distances change. True T/R modes don't change distances → eigenvalue
    // is zeroed. This is robust against imperfect T/R projection before diagonalization.
    // Reference: external/xtb/src/detrotra.f90 detrotra8()
    {
        struct ModeMetric { double metric; int index; };
        std::vector<ModeMetric> tr_metrics;
        tr_metrics.reserve(12);

        for (int ii = 0; ii < n3; ++ii) {
            if (evals(ii) > 0.05) continue; // only low-frequency modes (detrotra.f90 line 103)

            double c0 = 0.0;
            for (int i = 0; i < n_atoms; ++i) {
                for (int j = 0; j < i; ++j) {
                    // original distance
                    double dx = xyz(3*i+0) - xyz(3*j+0);
                    double dy = xyz(3*i+1) - xyz(3*j+1);
                    double dz = xyz(3*i+2) - xyz(3*j+2);
                    double a0 = std::sqrt(dx*dx + dy*dy + dz*dz);
                    // distorted position
                    double di0 = xyz(3*i+0) + eigenvectors(3*i+0, ii);
                    double di1 = xyz(3*i+1) + eigenvectors(3*i+1, ii);
                    double di2 = xyz(3*i+2) + eigenvectors(3*i+2, ii);
                    double dj0 = xyz(3*j+0) + eigenvectors(3*j+0, ii);
                    double dj1 = xyz(3*j+1) + eigenvectors(3*j+1, ii);
                    double dj2 = xyz(3*j+2) + eigenvectors(3*j+2, ii);
                    double ex = di0-dj0, ey = di1-dj1, ez = di2-dj2;
                    double b0 = std::sqrt(ex*ex + ey*ey + ez*ez);
                    c0 += (a0-b0)*(a0-b0);
                }
            }
            // Weight by |eigenvalue| (detrotra.f90 line 122)
            tr_metrics.push_back({std::sqrt(c0 / n_atoms) * std::abs(evals(ii)), ii});
        }

        std::sort(tr_metrics.begin(), tr_metrics.end(),
                  [](const ModeMetric& a, const ModeMetric& b){ return a.metric < b.metric; });

        int nend = is_linear ? 5 : 6;
        for (int i = 0; i < std::min(nend, (int)tr_metrics.size()); ++i)
            evals(tr_metrics[i].index) = 0.0; // mark as T/R
    }

    // Bug 5 Fix (Apr 2026): Select nvar modes by threshold (XTB-style) instead of fixed skip.
    // XTB: iterate from largest eigenvalue downward, include only |eigv| > thr.
    // 4-stage retry with falling threshold (generate_anc_blowup() Zeile 191-213).
    // After detrotra8, T/R modes have eigv=0 and are excluded by the threshold.
    // Sort selected modes ascending afterward (XTB calls sort() Zeile 220).
    B = Matrix::Zero(n3, nvar);
    hess = Matrix::Zero(nvar, nvar);
    {
        double thr = 1e-11;
        int nfound = 0;
        for (int itry = 0; itry < 4 && nfound < nvar; ++itry) {
            nfound = 0;
            B.setZero();
            hess.setZero();
            for (int i = n3 - 1; i >= 0 && nfound < nvar; --i) {
                if (std::abs(evals(i)) > thr) {
                    B.col(nfound) = eigenvectors.col(i);
                    hess(nfound, nfound) = std::min(std::max(evals(i), hlow), hmax);
                    ++nfound;
                }
            }
            if (nfound < nvar) thr *= 0.1;
        }
        if (nfound < nvar) {
            CurcumaLogger::warn_fmt("ANC: only {} of {} modes found — using available", nfound, nvar);
        }
        // Modes are in descending eigenvalue order; reverse to ascending (XTB sort())
        // Reverse columns of B and diagonal of hess
        for (int i = 0; i < nfound / 2; ++i) {
            int j = nfound - 1 - i;
            B.col(i).swap(B.col(j));
            std::swap(hess(i, i), hess(j, j));
        }
    }

    // Initialize internal coordinates from Cartesian
    setCartesian(xyz);

    CurcumaLogger::success("ANC generated with " + std::to_string(nvar) + " internal coordinates");

    return true;
}

bool ANCCoordinates::generateANCLanczos(const Matrix& H, const Vector& xyz,
                                         bool is_linear, int k_target) {
    // Truncated Lanczos eigendecomposition for large Cartesian Hessians.
    // Finds the nvar largest eigenpairs — the stiff (chemically meaningful)
    // modes that define the ANC basis. T/R modes have eigenvalue near zero
    // and are filtered out by the |eigv| > thr selection, so detrotra8
    // (which runs after full diagonalization) is not needed here.
    //
    // Cost: O(n3^2 * m) with m = k_target + oversampling, vs O(n3^3) for full.
    // For n3 = 4230, k_target = 500: ~8 Gflops vs ~80 Gflops — ~10x speedup.
    //
    // Algorithm: standard Lanczos with full reorthogonalization (two-pass
    // Gram-Schmidt per Parlett's "twice is enough"), then diagonalize the
    // resulting tridiagonal T and lift Ritz vectors back via V * Y.

    using clk = std::chrono::high_resolution_clock;
    auto t_start = clk::now();

    const int n = static_cast<int>(H.rows());
    if (n != n3) {
        CurcumaLogger::error("ANC Lanczos: Hessian dimension mismatch");
        return false;
    }

    // Lanczos subspace dimension: oversample to improve edge convergence.
    const int oversample = std::max(60, k_target / 5);
    const int m = std::min(n, k_target + oversample);

    Matrix V = Matrix::Zero(n, m);
    std::vector<double> alpha, beta;
    alpha.reserve(m);
    beta.reserve(m);

    // Starting vector: deterministic mild perturbation of the all-ones vector.
    // A pure ones vector can miss antisymmetric modes; the sinusoidal offset
    // seeds all symmetries while keeping results reproducible across runs.
    Vector v0(n);
    for (int i = 0; i < n; ++i) v0(i) = 1.0 + 0.01 * std::sin(0.1 * i + 1.0);
    v0.normalize();
    V.col(0) = v0;

    Vector w(n), v_prev(n);
    double beta_prev = 0.0;
    int m_actual = 0;

    for (int j = 0; j < m; ++j) {
        // w = H * v_j - beta_{j-1} * v_{j-1}
        w.noalias() = H * V.col(j);
        if (j > 0) w.noalias() -= beta_prev * v_prev;

        double a_j = V.col(j).dot(w);
        alpha.push_back(a_j);
        w.noalias() -= a_j * V.col(j);

        // Full reorthogonalization against V[:,0..j] — essential when seeking
        // many eigenpairs. Single GS pass can leave O(1) contamination; run
        // twice for orthogonality near machine precision.
        if (j > 0) {
            Vector c = V.leftCols(j + 1).transpose() * w;
            w.noalias() -= V.leftCols(j + 1) * c;
            c = V.leftCols(j + 1).transpose() * w;
            w.noalias() -= V.leftCols(j + 1) * c;
        }

        double b_j = w.norm();
        beta.push_back(b_j);
        m_actual = j + 1;

        if (b_j < 1e-12) break;  // Invariant subspace: Krylov space is complete.
        if (j + 1 < m) {
            v_prev = V.col(j);
            beta_prev = b_j;
            V.col(j + 1) = w / b_j;
        }
    }

    // Diagonalize the (small, cheap) tridiagonal Ritz matrix T.
    Matrix T = Matrix::Zero(m_actual, m_actual);
    for (int j = 0; j < m_actual; ++j) {
        T(j, j) = alpha[j];
        if (j + 1 < m_actual) {
            T(j, j + 1) = beta[j];
            T(j + 1, j) = beta[j];
        }
    }
    Eigen::SelfAdjointEigenSolver<Matrix> tri_solver(T);
    if (tri_solver.info() != Eigen::Success) {
        CurcumaLogger::error("ANC Lanczos: tridiagonal diagonalization failed");
        return false;
    }

    Vector ritz_vals = tri_solver.eigenvalues();   // ascending
    Matrix ritz_small = tri_solver.eigenvectors(); // m_actual x m_actual

    // XTB-style eigenvalue damping: shift so minimum |eigv| >= hlow. Mirrors the
    // full path at lines ~86-93 of generateANC.
    double elow = std::numeric_limits<double>::max();
    for (int i = 0; i < m_actual; ++i)
        if (std::abs(ritz_vals(i)) > 1e-10)
            elow = std::min(elow, ritz_vals(i));
    if (elow < std::numeric_limits<double>::max()) {
        double damp = std::max(hlow - elow, 0.0);
        if (damp > 0.0)
            for (int i = 0; i < m_actual; ++i)
                if (std::abs(ritz_vals(i)) > 1e-11)
                    ritz_vals(i) += damp;
    }

    // Lift all Ritz vectors to full space in one GEMM (BLAS-optimised).
    Matrix ritz_full = V.leftCols(m_actual) * ritz_small;  // n3 x m_actual

    // Threshold selection identical to the full path. Ritz values are in
    // ascending order, so iterate from the top (largest) downward. T/R modes
    // — which have ritz ~ 0 — are excluded by the |eigv| > thr filter.
    B = Matrix::Zero(n3, nvar);
    hess = Matrix::Zero(nvar, nvar);

    double thr = 1e-11;
    int nfound = 0;
    for (int itry = 0; itry < 4 && nfound < nvar; ++itry) {
        nfound = 0;
        B.setZero();
        hess.setZero();
        for (int i = m_actual - 1; i >= 0 && nfound < nvar; --i) {
            if (std::abs(ritz_vals(i)) > thr) {
                B.col(nfound) = ritz_full.col(i);
                hess(nfound, nfound) = std::min(std::max(ritz_vals(i), hlow), hmax);
                ++nfound;
            }
        }
        if (nfound < nvar) thr *= 0.1;
    }

    if (nfound < nvar) {
        CurcumaLogger::warn_fmt(
            "ANC Lanczos: captured {} of {} requested modes. Increase oversampling or "
            "lower anc_lanczos_k if convergence suffers.", nfound, nvar);
    }

    // Reverse collected (descending) order to ascending, matching XTB's sort().
    for (int i = 0; i < nfound / 2; ++i) {
        int j = nfound - 1 - i;
        B.col(i).swap(B.col(j));
        std::swap(hess(i, i), hess(j, j));
    }

    setCartesian(xyz);

    auto t_end = clk::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();
    CurcumaLogger::success_fmt(
        "ANC generated (Lanczos, n3={}, k={}, m={}) in {:.2f}s — top {} modes used",
        n3, k_target, m_actual, elapsed, nfound);
    return true;
}

void ANCCoordinates::getCartesian(Vector& xyz_out) const {
    // Claude Nov 2025: FIX - Transform displacement in internal coords back to Cartesian
    // The correct transformation is:
    //   delta_coord = coord - coord_reference  (accumulated displacement in internal space)
    //   delta_xyz = B * delta_coord  (transform displacement to Cartesian)
    //   xyz_new = xyz_reference + delta_xyz  (add to reference geometry)
    //
    // This fixes the bug where we were doing xyz = B * coord (which is completely wrong!)

    Vector delta_coord = coord - coord_reference;
    Vector delta_xyz = B * delta_coord;
    xyz_out = xyz_reference + delta_xyz;
}

void ANCCoordinates::setCartesian(const Vector& xyz_in) {
    // Claude Nov 2025: FIX - Set reference point and calculate internal coords
    // Store the reference geometry
    xyz_reference = xyz_in;
    xyz_flat = xyz_in;

    // Calculate reference internal coordinates
    // coord_ref = B^T * xyz (project Cartesian onto internal space)
    coord_reference = B.transpose() * xyz_in;
    coord = coord_reference;  // Initialize coord to reference
}

Vector ANCCoordinates::transformGradientToInternal(const Vector& cartesian_gradient) const {
    // Transform Cartesian gradient to internal gradient
    // g_internal = B^T * g_cartesian
    return B.transpose() * cartesian_gradient;
}

// ============================================================================
// ModelHessianParameters Implementation
// ============================================================================

ModelHessianParameters ModelHessianParameters::fromJson(const json& config) {
    ModelHessianParameters params;

    if (config.contains("model_hessian")) {
        std::string model_str = config["model_hessian"].get<std::string>();
        if (model_str == "lindh" || model_str == "lindh2007") {
            params.model = LINDH_2007;
        } else if (model_str == "lindh1995") {
            params.model = LINDH_1995;
        } else if (model_str == "lindh_d2") {
            params.model = LINDH_D2;
        } else if (model_str == "swart") {
            params.model = SWART;
        } else if (model_str == "read") {
            params.model = READ_FILE;
        }
    }

    if (config.contains("s6")) {
        params.s6 = config["s6"].get<double>();
    }

    return params;
}

// ============================================================================
// ANCOptimizer Implementation
// ============================================================================

ANCOptimizer::ANCOptimizer()
    : OptimizerDriver() {
    m_anc = std::make_unique<ANCCoordinates>();
}

bool ANCOptimizer::InitializeOptimizerInternal() {
    CurcumaLogger::info("Initializing AncOpt optimizer");

    // Claude Nov 2025: BUG FIX - Check molecule is valid before proceeding
    int atom_count = m_molecule.AtomCount();

    if (atom_count == 0) {
        CurcumaLogger::error_fmt("Cannot initialize ANC optimizer with empty molecule (AtomCount={})", atom_count);
        return false;
    }
    if (atom_count < 2) {
        CurcumaLogger::error_fmt("ANC optimizer requires at least 2 atoms (AtomCount={})", atom_count);
        return false;
    }

    CurcumaLogger::info_fmt("Molecule: {} atoms, charge {}", atom_count, m_molecule.Charge());

    // Load ANC-specific parameters from configuration
    loadANCParameters(m_configuration);

    // Set optimization level (affects thresholds)
    int opt_level = m_configuration.value("optimization_level", 0); // 0 = normal
    setOptimizationLevel(opt_level);

    // Check if molecule is linear
    bool is_linear = checkLinearMolecule(m_molecule);

    // Calculate number of internal coordinates
    int nvar = 3 * m_molecule.AtomCount() - (is_linear ? 5 : 6);

    // Claude Nov 2025: BUG FIX - Ensure nvar is positive
    if (nvar <= 0) {
        CurcumaLogger::error_fmt("Invalid number of internal coordinates: nvar={}", nvar);
        return false;
    }

    CurcumaLogger::info_fmt("Internal coordinates: {} (is_linear={})", nvar, is_linear);

    // Size-tier selection for the ANC basis. Thresholds on n3 = 3*N_atoms.
    //   Tier S/M (n3 <= 1800): full eigendecomposition, nvar = 3N-6.
    //   Tier L  (n3 > 1800):   truncated Lanczos, nvar capped at anc_lanczos_k.
    // The cap matters because dense BFGS still uses an nvar x nvar Hessian;
    // capping nvar bounds both setup cost (Lanczos dim m ≈ k + 60) and per-step
    // cost (BFGS rank-2 update is O(nvar^2)).
    const int n3_threshold = m_configuration.value("anc_lanczos_threshold", 1800);
    const int k_max = m_configuration.value("anc_lanczos_k", 500);
    const int n3_total = 3 * atom_count;
    const int nvar_full = nvar;

    const int lbfgs_xl_thr = m_configuration.value("anc_lbfgs_threshold", m_anc_lbfgs_threshold);
    // Use fmt::print directly — CurcumaLogger global verbosity may be suppressed by
    // EnergyCalculator initialization at this point. m_context.verbosity is the
    // user-requested level (set in LoadConfiguration) and is always reliable.
    if (n3_total > n3_threshold && nvar > k_max) {
        nvar = k_max;
        const char* step_method = (k_max > lbfgs_xl_thr) ? "L-BFGS (Tier XL)" : "BFGS+RF (Tier L)";
        if (m_context.verbosity >= 1)
            fmt::print("ANCOpt: {} atoms, n3={} (Tier L) — Lanczos ANC, nvar capped at {}/{}, step: {}\n"
                       "  For large-system minima, -optimizer lbfgspp is typically faster.\n",
                       atom_count, n3_total, k_max, nvar_full, step_method);
    } else if (nvar > lbfgs_xl_thr) {
        if (m_context.verbosity >= 1)
            fmt::print("ANCOpt: {} atoms, nvar={} (Tier XL) — L-BFGS in ANC subspace, {} history pairs.\n"
                       "  RF mode-following disabled. TS searches not supported at this size.\n",
                       atom_count, nvar, m_anc_lbfgs_history);
    } else if (nvar > 1000) {
        double mb_est = (static_cast<double>(nvar) * nvar * 8.0) / (1024.0 * 1024.0);
        if (m_context.verbosity >= 1)
            fmt::print("ANCOpt: {} atoms, nvar={} (Tier M) — full eigendecomp, ~{:.0f} MB, may take ~1 min setup.\n",
                       atom_count, nvar, mb_est);
    } else {
        if (m_context.verbosity >= 2)
            fmt::print("ANCOpt: {} atoms, nvar={} (Tier S) — full eigendecomp + RF.\n",
                       atom_count, nvar);
    }

    // Allocate ANC structure
    m_anc->allocate(m_molecule.AtomCount(), nvar, m_hlow, m_hmax);

    // Generate initial model Hessian
    Matrix cart_hess = generateModelHessian(m_molecule);

    // Project out translations and rotations
    projectTranslationsRotations(cart_hess, m_molecule);

    // Generate ANC from Hessian
    if (!generateANCFromHessian(cart_hess, m_molecule)) {
        CurcumaLogger::error("Failed to generate initial ANC");
        return false;
    }

    // Initialize optimization state
    m_gint = Vector::Zero(m_anc->nvar);
    m_gint_old = Vector::Zero(m_anc->nvar);
    m_displ = Vector::Zero(m_anc->nvar);
    m_micro_current = 0;
    m_needs_anc_regeneration = false;

    CurcumaLogger::param("Max displacement (ANC)", std::to_string(m_maxdispl) + " Bohr");
    CurcumaLogger::param("Frequency cutoffs", std::to_string(m_hlow) + " - " + std::to_string(m_hmax));
    CurcumaLogger::param("Micro-iterations", std::to_string(m_maxmicro));
    CurcumaLogger::param("Hessian update", m_hessian_update == BFGS ? "BFGS" : "Powell");

    return true;
}

Vector ANCOptimizer::CalculateOptimizationStep(const Vector& current_coordinates,
                                                const Vector& gradient) {
    using clk = std::chrono::high_resolution_clock;
    auto t0 = clk::now();

    // Transform Cartesian gradient to internal coordinates.
    // GFN-FF (and all XTB-interface methods) return the gradient in atomic units
    // (Eh/Bohr), but Curcuma geometry is stored in Angstrom.  The ANC model
    // Hessian is therefore built in Eh/Ang².  To keep the RF step and BFGS update
    // dimensionally consistent (Newton step [Ang], dg·dx [Eh]) we convert here:
    //
    //   dE/d(r_Ang) = dE/d(r_Bohr) × d(r_Bohr)/d(r_Ang)
    //               = gradient_Bohr  × ANGSTROM_TO_BOHR    [Eh/Ang]
    //
    // Note: ANGSTROM_TO_BOHR ≈ 1.88973, NOT BOHR_TO_ANGSTROM ≈ 0.52918.
    // The B matrix is dimensionless (eigenvectors), so g_int = B^T · g_cart_ang
    // inherits the Eh/Ang unit.  With H [Eh/Ang²] and g [Eh/Ang], the Newton
    // step p/s has units Ang — consistent with ANC coordinates stored in Ang.
    const Vector gradient_ang = gradient * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
    m_gint_old = m_gint;
    m_gint = m_anc->transformGradientToInternal(gradient_ang);
    m_gnorm_old = m_gnorm;
    m_gnorm = m_gint.norm();

    auto t1 = clk::now();
    m_time_transform += std::chrono::duration<double>(t1 - t0).count();

    // Check if we need to regenerate ANC
    if (m_needs_anc_regeneration || m_micro_current >= m_maxmicro) {
        auto tr0 = clk::now();
        if (m_anc->n3 > 1800)
            CurcumaLogger::info_fmt("Regenerating ANC (Tier L, Lanczos, n3={}, k={})",
                                     m_anc->n3, m_anc->nvar);
        else
            CurcumaLogger::info("Regenerating ANC (micro-iteration limit reached)");

        // Generate new model Hessian
        Matrix cart_hess = generateModelHessian(m_molecule);
        projectTranslationsRotations(cart_hess, m_molecule);

        // Regenerate ANC
        Vector current_xyz(3 * m_molecule.AtomCount());
        for (int i = 0; i < m_molecule.AtomCount(); ++i) {
            current_xyz.segment<3>(3 * i) = m_molecule.getGeometry().row(i);
        }
        m_anc->generateANC(cart_hess, current_xyz, checkLinearMolecule(m_molecule));

        m_micro_current = 0;
        m_needs_anc_regeneration = false;
        // ANC basis changed — L-BFGS history is in the old basis and must be discarded.
        m_lbfgs_s_int.clear();
        m_lbfgs_y_int.clear();
        m_lbfgs_rho_int.clear();

        // XTB optimizer.f90:407 — grow maxmicro by 10% each macro cycle, capped at 2x initial.
        // Rationale: as optimization progresses the Hessian model becomes less critical,
        // so ANC regeneration can happen less frequently, saving 3x O(N^3) work per skipped regen.
        if (m_current_iteration > 0) {
            m_maxmicro = std::min(static_cast<int>(m_maxmicro * 1.1),
                                  2 * m_maxmicro_initial);
        }

        // Critical: re-project gradient onto NEW ANC basis after regeneration.
        // m_gint was computed above with the OLD B matrix; m_anc->hess now uses the
        // NEW eigenvectors, so mixing them in the RF step would give a garbage direction.
        // Re-projecting here ensures m_gint and m_anc->hess are in the same basis.
        // Use gradient_ang (Eh/Ang) to stay consistent with the model Hessian (Eh/Ang²).
        m_gint = m_anc->transformGradientToInternal(gradient_ang);
        m_gnorm = m_gint.norm();
        m_gint_old = m_gint;                               // no valid old basis to compare
        m_displ = Vector::Zero(m_anc->nvar);               // old displacement is in old basis → invalid
        m_last_rf_eigenvector.resize(0);                   // basis change invalidates Lanczos warm-start
        m_time_anc_regen += std::chrono::duration<double>(clk::now() - tr0).count();
    }

    const int anc_lbfgs_thr = m_configuration.value("anc_lbfgs_threshold", m_anc_lbfgs_threshold);
    const bool use_lbfgs_xl = (m_anc->nvar > anc_lbfgs_thr);

    if (use_lbfgs_xl) {
        // ── Tier XL: L-BFGS in ANC subspace ────────────────────────────────
        // Dense nvar×nvar hess would be > 32 MB and O(nvar²) BFGS update too slow.
        // Store (s_k, y_k) pairs and apply H^-1 via two-loop recursion (O(nvar·m)).
        // RF mode-following not available; use backtracking trust-radius step instead.
        auto tlbfgs0 = clk::now();

        if (m_current_iteration > 0 && m_displ.norm() > 1e-10) {
            Vector dg = m_gint - m_gint_old;
            double sy = m_displ.dot(dg);
            if (sy > 1e-14) { // curvature condition
                if (static_cast<int>(m_lbfgs_s_int.size()) >= m_anc_lbfgs_history) {
                    m_lbfgs_s_int.erase(m_lbfgs_s_int.begin());
                    m_lbfgs_y_int.erase(m_lbfgs_y_int.begin());
                    m_lbfgs_rho_int.erase(m_lbfgs_rho_int.begin());
                }
                m_lbfgs_s_int.push_back(m_displ);
                m_lbfgs_y_int.push_back(dg);
                m_lbfgs_rho_int.push_back(1.0 / sy);
            }
        }

        Vector lbfgs_dir = -calculateLBFGSStepInternal(m_gint); // descent direction
        m_displ = lbfgs_dir;
        ++m_lbfgs_step_calls;
        m_time_lbfgs_step += std::chrono::duration<double>(clk::now() - tlbfgs0).count();

        // No dense hess to update; predict energy linearly (g·dx).
        m_depred = m_gint.dot(m_displ);

    } else {
        // ── Tier S/M/L: BFGS/Powell + RF step (original path) ──────────────
        auto tbfgs0 = clk::now();
        // Guard m_displ.norm() > 1e-10 skips the first step after ANC regeneration
        // (m_displ was zeroed above) and the very first iteration.
        if (m_current_iteration > 0 && m_displ.norm() > 1e-10) {
            Vector dg = m_gint - m_gint_old;

            if (m_hessian_update == BFGS) {
                updateHessianBFGS(m_anc->hess, m_displ, dg);
            } else {
                updateHessianPowell(m_anc->hess, m_displ, dg);
            }
        }
        m_time_bfgs += std::chrono::duration<double>(clk::now() - tbfgs0).count();

        // Calculate Rational Function step
        auto trf0 = clk::now();
        m_displ = calculateRationalFunctionStep(m_gint, m_anc->hess);
        m_time_rf_step += std::chrono::duration<double>(clk::now() - trf0).count();
    }

    // XTB optimizer.f90:707-712 — only per-component clipping, no global norm limit.
    // Earlier Curcuma versions imposed an additional global trust-radius on ||dx||,
    // which is NOT in the reference. For systems with many internal coordinates the
    // sqrt(nvar) scaling of the true step norm was shrunk to a tiny per-step radius,
    // producing "slurring" optimization with many small steps.
    double maxdispl_ang = m_maxdispl * 0.529177; // Bohr → Å
    for (int i = 0; i < m_displ.size(); ++i) {
        if (std::abs(m_displ(i)) > maxdispl_ang) {
            m_displ(i) = std::copysign(maxdispl_ang, m_displ(i));
        }
    }

    // XTB: exit micro-loop (= regenerate ANC) if step norm > 2.0 Bohr after iter 2.
    // This prevents the Hessian from being updated with a garbage step.
    // Reference: XTB optimizer.f90 relax() Zeile 733-736
    if (m_micro_current > 2 && m_displ.norm() > 2.0 * 0.529177) {
        CurcumaLogger::warn_fmt("Large step detected ({:.3f} A), forcing ANC regeneration",
                                m_displ.norm());
        m_needs_anc_regeneration = true;
    }

    // Predict energy change — XL tier already set m_depred (linear g·dx); others use 2nd order.
    if (!use_lbfgs_xl)
        m_depred = predictEnergyChange(m_gint, m_displ, m_anc->hess);

    // Adaptive step size scaling (XTB optimizer.f90:648).
    // XTB thresholds are in Eh/Bohr; m_gnorm is now in Eh/Ang after the unit fix.
    // Convert XTB values: X Eh/Bohr × ANGSTROM_TO_BOHR (Bohr/Ang) = X Eh/Ang equivalent.
    // dE/d(r_Ang) = dE/d(r_Bohr) × ANGSTROM_TO_BOHR  → multiply thresholds by A2B.
    // Note: XTB's elseif chain has dead branches (0.0006 and 0.0003 unreachable
    // after 0.002 triggers first) — Curcuma's correct descending order is intentional.
    double alp = 1.0;
    constexpr double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;
    if      (m_gnorm < 0.0003 * A2B) alp = 3.0;
    else if (m_gnorm < 0.0006 * A2B) alp = 2.0;
    else if (m_gnorm < 0.002  * A2B) alp = 1.5;

    // Update internal coordinates with scaled displacement
    m_anc->coord += m_displ * alp;

    // Transform back to Cartesian (returns absolute new coordinates)
    Vector new_xyz;
    m_anc->getCartesian(new_xyz);

    m_micro_current++;

    if (m_context.verbosity >= 2) {
        const char* step_type = use_lbfgs_xl ? "L-BFGS(XL)" : "BFGS+RF";
        const bool anc_regen = (m_micro_current == 0); // just regenerated
        fmt::print("  ANCOpt micro={}/{} type={} ||g_int||={:.3e} ||dx||={:.3e} depred={:.3e} Eh{}\n",
                   m_micro_current + 1, m_maxmicro,
                   step_type, m_gnorm, m_displ.norm(), m_depred,
                   anc_regen ? " [ANC regen]" : "");
    }
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("Displacement norm (ANC)", std::to_string(m_displ.norm()));
        CurcumaLogger::param("Predicted energy change", std::to_string(m_depred) + " Eh");
        CurcumaLogger::param("Step scaling (alp)", std::to_string(alp));
    }

    // Bug 1 Fix: Return displacement (new_xyz - current_xyz) NOT absolute coordinates
    // OptimizerDriver::Optimize() adds step to current coords: new = current + step
    return new_xyz - current_coordinates;
}

bool ANCOptimizer::CheckMethodSpecificConvergence() const {
    // XTB convergence: ALL THREE conditions simultaneously.
    // Reference: XTB optimizer.f90 relax() Zeile 742-747:
    //   if(abs(echng).lt.ethr .and. gnorm.lt.gthr .and. echng.lt.0) converged=.true.
    //
    // Old Curcuma code accepted (|dE| < ethr OR dE < 0) — too loose, caused premature
    // convergence when energy was still decreasing by large amounts.
    bool grad_converged = m_gnorm < m_gthr;
    bool energy_small   = std::abs(m_energy_change) < m_ethr;
    bool energy_lowered = m_energy_change < 0.0;

    // First iteration: no previous energy → skip energy condition
    bool energy_ok = (m_current_iteration <= 1) || (energy_small && energy_lowered);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info_fmt("Convergence: gnorm={:.3e} (thr={:.3e}), dE={:.3e} (thr={:.3e}), lowered={}",
            m_gnorm, m_gthr, m_energy_change, m_ethr, energy_lowered);
    }

    return grad_converged && energy_ok;
}

void ANCOptimizer::UpdateOptimizerState(const Vector& new_coordinates,
                                        const Vector& new_gradient,
                                        double new_energy) {
    // Bug 3 Fix: Track energy change before updating m_current_energy
    // CheckMethodSpecificConvergence() needs this to test echng < 0 condition
    m_energy_change = new_energy - m_current_energy;
    m_previous_energy = m_current_energy;

    // XTB reference (optimizer.f90 relax()): energy rises are NOT treated specially.
    // XTB simply continues the micro-loop; BFGS corrects the Hessian and the next
    // RF step will be better. ANC regeneration happens only after maxmicro steps or
    // large displacement (dsnrm > 2.0 Bohr) — NOT on every energy rise.
    //
    // Old Curcuma code set m_needs_anc_regeneration = true on any energy rise and
    // reduced trust radius by 0.25x. This caused expensive ANC regeneration (3x O(N³))
    // on every early-optimization step where energy briefly rises, which is why
    // the optimizer was slow: ANC was rebuilt every 1-2 steps instead of every 25.
    //
    // Fix: only track for convergence check. Large-step exit in CalculateOptimizationStep
    // handles the case where the step is genuinely too large.

    // Store current state
    m_current_energy = new_energy;
    m_current_gradient = new_gradient;
}

void ANCOptimizer::FinalizeOptimizationInternal() {
    // Phase timing summary — shown at verbosity >= 2 for profiling/diagnostics.
    // Uses fmt::print because CurcumaLogger global verbosity is set to 0
    // during optimization to suppress EnergyCalculator line-search spam.
    if (m_context.verbosity >= 2) {
        const int n_rf = m_rf_lanczos_calls + m_rf_fallback_calls;
        fmt::print("\nAncOpt phase timings (cumulative seconds):\n");
        fmt::print("  ANC regeneration    : {:.3f}\n", m_time_anc_regen);
        fmt::print("  BFGS/Powell update  : {:.3f}\n", m_time_bfgs);
        fmt::print("  RF step total       : {:.3f}\n", m_time_rf_step);
        fmt::print("    Lanczos-only      : {:.3f}\n", m_time_rf_lanczos);
        fmt::print("    Fallback-only     : {:.3f}\n", m_time_rf_fallback);
        fmt::print("  Gradient transform  : {:.3f}\n", m_time_transform);
        fmt::print("  RF calls Lanczos    : {}\n", m_rf_lanczos_calls);
        fmt::print("  RF calls fallback   : {}\n", m_rf_fallback_calls);
        if (n_rf > 0) {
            fmt::print("  RF avg time/call    : {:.2f} ms\n",
                       1000.0 * m_time_rf_step / n_rf);
        }
        if (m_rf_lanczos_calls > 0) {
            fmt::print("  Lanczos avg iters   : {:.1f}\n",
                       double(m_rf_lanczos_iters_total) / m_rf_lanczos_calls);
        }
        if (m_lbfgs_step_calls > 0) {
            fmt::print("  L-BFGS XL steps     : {}\n", m_lbfgs_step_calls);
            fmt::print("  L-BFGS XL time      : {:.3f} s\n", m_time_lbfgs_step);
        }
    }

    // Explicitly release ANC matrices before the optimizer object is destroyed.
    // This avoids heap-corruption crashes (signal 11 in __libc_free) that occur
    // when large Eigen allocations are freed in the wrong order during destruction.
    m_anc.reset();
    m_displ.resize(0);
    m_gint.resize(0);
    m_gint_old.resize(0);

    CurcumaLogger::success("AncOpt optimization finalized");
}

// ============================================================================
// Model Hessian Generation
// ============================================================================

Matrix ANCOptimizer::generateModelHessian(const Molecule& mol) {
    // Claude Generated (Mar 2026): Lindh 1995 model Hessian.
    // Reference: R. Lindh, A. Bernhardsson, G. Karlström, P.-A. Malmqvist,
    //            Chem. Phys. Lett. 241 (1995) 423-428
    //
    // The model Hessian is built from exponentially-damped bond-direction and
    // angle-bending contributions. It approximates the true Hessian well enough
    // to produce a good ANC basis, without needing an actual computation.
    //
    // STRETCHING: For each atom pair (i,j):
    //   e_ij = unit vector from i to j
    //   rho_ij = exp(-alpha_ij * (r_ij - r0_ij)^2)   [r in Bohr]
    //   H += k_str * rho_ij * (e_ij ⊗ e_ij)  (converted to Eh/Å²)
    //
    // BENDING: For each triplet (i,j,k) with j as central atom:
    //   Wilson B-matrix vectors di, dj, dk give d(angle)/d(Cartesian)
    //   rho = rho_ij * rho_jk
    //   H += k_bend * rho * (di ⊗ di + dj ⊗ dj + dk ⊗ dk + cross terms)

    CurcumaLogger::info("Generating model Hessian (Lindh 1995)");

    const int N = mol.AtomCount();
    const int n3 = 3 * N;
    Matrix H = Matrix::Zero(n3, n3);

    // Conversion: 1 Bohr = 0.529177 Å
    // Geometry stored in Å, parameters in Bohr.
    // Hessian units: Eh/Å² (consistent with gradient in Eh/Å).
    const double A2B = 1.0 / 0.529177; // Å → Bohr
    const double A2B2 = A2B * A2B;     // for converting Eh/Bohr² → Eh/Å²

    // Lindh 1995 reference distances r0 [Bohr] (Table 1), indexed by period:
    //   row 0 = H, row 1 = period 2 (Li-Ne), row 2 = period 3+ (Na+)
    // r0[ri][rj] = reference distance for pair in periods ri, rj
    const double r0[3][3] = {
        {1.35, 2.10, 2.53}, // H-H,     H-Li..Ne,  H-Na+
        {2.10, 2.87, 3.40}, // Li..Ne-H, Li..Ne pair, Li..Ne - Na+
        {2.53, 3.40, 3.40}  // Na+-H,   Na+-Li..Ne, Na+-Na+
    };

    // Lindh 1995 exponential decay parameters alpha [Bohr^-2] (Table 1)
    // alpha_HH=1.00, all others ~0.316 = sqrt(0.1)
    const double alpha_str[3][3] = {
        {1.000, 0.316, 0.316},
        {0.316, 0.316, 0.316},
        {0.316, 0.316, 0.316}
    };

    // For angle bending, Lindh uses softer alpha values (~0.12)
    const double alpha_bend[3][3] = {
        {0.120, 0.120, 0.120},
        {0.120, 0.120, 0.120},
        {0.120, 0.120, 0.120}
    };

    // Reference force constants [Eh/Bohr²]
    // Stretching and bending values from Lindh 1995
    const double k_str  = 0.45; // Eh/Bohr² (bond stretching)
    const double k_bend = 0.15; // Eh/rad² (angle bending, rad because B-matrix is dangle/dcoord)

    // Cutoff: only include pairs within 8 Bohr (~4.2 Å)
    const double r_cutoff_bohr = 8.0;

    // Assign each atom to a period row (for parameter table lookup)
    auto getPeriod = [](int Z) -> int {
        if (Z == 1) return 0;   // H
        if (Z <= 10) return 1;  // Li-Ne
        return 2;               // Na and beyond (simplified to period 3 parameters)
    };

    // Convert geometry to Bohr for all distance/vector calculations
    const auto& geom_ang = mol.getGeometry(); // Rows: atoms, Cols: x,y,z in Å

    // === BOND STRETCHING CONTRIBUTIONS ===
    // For each pair (i,j), add k_str*rho * (e_ij ⊗ e_ij) to the appropriate blocks.
    // This correctly represents the second derivative of a harmonic bond potential.
    for (int i = 0; i < N; ++i) {
        int pi = getPeriod(mol.Atom(i).first);
        for (int j = i + 1; j < N; ++j) {
            int pj = getPeriod(mol.Atom(j).first);

            // Distance in Bohr
            Eigen::Vector3d dr_ang = geom_ang.row(j) - geom_ang.row(i);
            double r_ang = dr_ang.norm();
            double r_bohr = r_ang * A2B;

            if (r_bohr > r_cutoff_bohr)
                continue;

            // Exponential damping: rho = exp(-alpha * (r - r0)^2)
            double dr = r_bohr - r0[pi][pj];
            double rho = std::exp(-alpha_str[pi][pj] * dr * dr);

            // Effective force constant in Eh/Bohr² → convert to Eh/Å²
            double k_eff_ang = k_str * rho * A2B2;

            // Bond unit vector (same direction in Å or Bohr space)
            Eigen::Vector3d e_ij = dr_ang / r_ang;

            // Outer product: e_ij ⊗ e_ij (dimensionless)
            // Second derivative of harmonic potential E=k/2*(r-r0)^2:
            //   d²E/dx_i dx_j(bond) = k * e_ij * e_ij^T  (i,j same atom → positive)
            //   d²E/dx_i dx_k(bond) = -k * e_ij * e_ij^T (i,k different atoms → negative)
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    double outer_ab = e_ij(a) * e_ij(b) * k_eff_ang;
                    H(3*i+a, 3*i+b) += outer_ab;  // atom i diagonal block: +
                    H(3*j+a, 3*j+b) += outer_ab;  // atom j diagonal block: +
                    H(3*i+a, 3*j+b) -= outer_ab;  // i-j off-diagonal block: -
                    H(3*j+a, 3*i+b) -= outer_ab;  // j-i off-diagonal block: -
                }
            }
        }
    }

    // === ANGLE BENDING CONTRIBUTIONS ===
    // For each triplet (i,j,k) with j as central atom, add Wilson B-matrix contributions.
    // Wilson B-vector: di = d(angle)/d(coord_i), units = 1/Å
    //   di = (cos(θ)*e_ji - e_jk) / (r_ji * sin(θ))  [in 1/Å if r in Å]
    //   dk = (cos(θ)*e_jk - e_ji) / (r_jk * sin(θ))
    //   dj = -(di + dk)
    // H contribution: k_bend * rho_ij * rho_jk * (da ⊗ db) for each atom pair (a,b)
    for (int j = 0; j < N; ++j) {
        int pj = getPeriod(mol.Atom(j).first);
        for (int i = 0; i < N; ++i) {
            if (i == j) continue;
            int pi = getPeriod(mol.Atom(i).first);

            Eigen::Vector3d dr_ji_ang = geom_ang.row(i) - geom_ang.row(j);
            double r_ji_ang = dr_ji_ang.norm();
            double r_ji_bohr = r_ji_ang * A2B;
            if (r_ji_bohr > r_cutoff_bohr) continue;

            double dri = r_ji_bohr - r0[pj][pi];
            double rho_ij = std::exp(-alpha_bend[pj][pi] * dri * dri);

            for (int k = i + 1; k < N; ++k) {
                if (k == j) continue;
                int pk = getPeriod(mol.Atom(k).first);

                Eigen::Vector3d dr_jk_ang = geom_ang.row(k) - geom_ang.row(j);
                double r_jk_ang = dr_jk_ang.norm();
                double r_jk_bohr = r_jk_ang * A2B;
                if (r_jk_bohr > r_cutoff_bohr) continue;

                double drk = r_jk_bohr - r0[pj][pk];
                double rho_jk = std::exp(-alpha_bend[pj][pk] * drk * drk);

                // Combined damping for angle i-j-k
                double rho = rho_ij * rho_jk;

                // Compute angle i-j-k (in radians)
                Eigen::Vector3d e_ji = dr_ji_ang / r_ji_ang;
                Eigen::Vector3d e_jk = dr_jk_ang / r_jk_ang;
                double cos_theta = e_ji.dot(e_jk);
                cos_theta = std::max(-1.0 + 1e-8, std::min(1.0 - 1e-8, cos_theta));
                double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
                if (sin_theta < 1e-6) continue; // Skip (nearly) linear angles

                // Wilson B-matrix vectors in 1/Å
                // These are the derivatives of the angle (in rad) with respect to Å
                Eigen::Vector3d di = (cos_theta * e_ji - e_jk) / (r_ji_ang * sin_theta);
                Eigen::Vector3d dk = (cos_theta * e_jk - e_ji) / (r_jk_ang * sin_theta);
                Eigen::Vector3d dj = -(di + dk);

                // Effective force constant for bending in Eh/Å²
                // k_bend [Eh/rad²] * di [rad/Å] * dj [rad/Å] → [Eh/Å²] ✓
                double k_eff = k_bend * rho;

                // Atoms involved: i (index i), j (index j), k (index k)
                const Eigen::Vector3d* dvec[3] = {&di, &dj, &dk};
                const int idx[3] = {i, j, k};

                // Add all 9 block contributions (symmetric)
                for (int a = 0; a < 3; ++a) {
                    for (int b = a; b < 3; ++b) {
                        for (int p = 0; p < 3; ++p) {
                            for (int q = 0; q < 3; ++q) {
                                double contrib = k_eff * (*dvec[a])(p) * (*dvec[b])(q);
                                H(3*idx[a]+p, 3*idx[b]+q) += contrib;
                                if (a != b) {
                                    H(3*idx[b]+q, 3*idx[a]+p) += contrib; // symmetry
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Symmetrize (should already be symmetric, but ensure numerical symmetry)
    for (int i = 0; i < n3; ++i)
        for (int j = i + 1; j < n3; ++j) {
            double avg = 0.5 * (H(i, j) + H(j, i));
            H(i, j) = avg;
            H(j, i) = avg;
        }

    CurcumaLogger::success("Lindh 1995 model Hessian generated (" + std::to_string(N) + " atoms)");
    return H;
}



bool ANCOptimizer::generateANCFromHessian(const Matrix& cart_hess, const Molecule& mol) {
    // Convert Cartesian coordinates to flat vector
    Vector xyz_flat(3 * mol.AtomCount());
    const auto& geom = mol.getGeometry();
    for (int i = 0; i < mol.AtomCount(); ++i) {
        xyz_flat.segment<3>(3 * i) = geom.row(i);
    }

    // Generate ANC
    bool is_linear = checkLinearMolecule(mol);
    return m_anc->generateANC(cart_hess, xyz_flat, is_linear);
}

// ============================================================================
// Rational Function Step Calculation
// ============================================================================

bool ANCOptimizer::lanczosLowestEigenpair(const Matrix& hessian,
                                          const Vector& gradient,
                                          const Vector& start_vector,
                                          Vector& out_eigenvector,
                                          double& out_eigenvalue,
                                          int m_max,
                                          double tol)
{
    // Delegate to shared RFSolver — keeps ANCOpt counter for timing reports.
    bool ok = RFSolver::lanczosLowestEigenpair(hessian, gradient, start_vector,
                                               out_eigenvector, out_eigenvalue,
                                               m_max, tol);
    if (ok) m_rf_lanczos_iters_total += m_max; // approximate; exact count not returned
    return ok;
}

Vector ANCOptimizer::calculateRationalFunctionStep(const Vector& gradient_internal,
                                                    const Matrix& hessian_internal) {
    /*
     * Rational Function Optimization (RFO) step — delegates to RFSolver::calculateRFStep.
     * Keeps per-call timing and Lanczos/fallback counters for the ANCOpt timing report.
     *
     * Algorithm: augmented eigenvalue problem [[H,g],[g',0]], Banerjee et al. 1985.
     * Lanczos path for nvar+1 >= 50, full SelfAdjointEigenSolver fallback.
     * See rf_solver.cpp for implementation details.
     */
    using rf_clk = std::chrono::high_resolution_clock;
    auto t0 = rf_clk::now();

    const int nvar1 = static_cast<int>(gradient_internal.size()) + 1;
    const bool had_warm = (m_last_rf_eigenvector.size() == nvar1);

    Vector result = RFSolver::calculateRFStep(gradient_internal, hessian_internal,
                                              m_last_rf_eigenvector);

    const bool used_lanczos = (m_last_rf_eigenvector.size() == nvar1);
    if (used_lanczos) {
        ++m_rf_lanczos_calls;
        m_time_rf_lanczos += std::chrono::duration<double>(rf_clk::now() - t0).count();
    } else {
        ++m_rf_fallback_calls;
        m_time_rf_fallback += std::chrono::duration<double>(rf_clk::now() - t0).count();
    }
    (void)had_warm;

    return result;
}

// ============================================================================
// L-BFGS Two-Loop Recursion in ANC Space (Tier XL)
// ============================================================================

Vector ANCOptimizer::calculateLBFGSStepInternal(const Vector& gradient)
{
    /*
     * L-BFGS two-loop recursion in ANC internal coordinates — Nocedal & Wright §7.2, Alg. 7.4.
     *
     * Used for Tier XL (nvar > m_anc_lbfgs_threshold) where the dense nvar×nvar Hessian
     * would exceed ~32 MB and the O(nvar²) BFGS rank-2 update becomes expensive.
     * Instead we store the last m_anc_lbfgs_history (s_k, y_k) pairs and apply the
     * implicit H^-1 product via the two-loop recursion (O(nvar * m_history) per step).
     *
     * Initial Hessian H₀ = γI with γ = (s_{k-1}·y_{k-1})/(y_{k-1}·y_{k-1}) (N&W eq. 7.20),
     * which scales the identity to approximate the local curvature.
     *
     * On empty history (first step after regeneration), returns gradient — caller negates.
     */
    const int m = static_cast<int>(m_lbfgs_s_int.size());
    if (m == 0) return gradient; // caller negates to steepest descent

    // First loop (backward).
    Vector q = gradient;
    std::vector<double> alpha_list(m);
    for (int i = m - 1; i >= 0; --i) {
        alpha_list[i] = m_lbfgs_rho_int[i] * m_lbfgs_s_int[i].dot(q);
        q.noalias() -= alpha_list[i] * m_lbfgs_y_int[i];
    }

    // Scale by γ from most recent pair (N&W eq. 7.20).
    const double sy = m_lbfgs_s_int.back().dot(m_lbfgs_y_int.back());
    const double yy = m_lbfgs_y_int.back().squaredNorm();
    const double gamma = (yy > 1e-14) ? sy / yy : 1.0;
    Vector z = gamma * q;

    // Second loop (forward).
    for (int i = 0; i < m; ++i) {
        double beta = m_lbfgs_rho_int[i] * m_lbfgs_y_int[i].dot(z);
        z.noalias() += m_lbfgs_s_int[i] * (alpha_list[i] - beta);
    }

    return z; // This is H^-1 * g. Caller negates to get descent direction.
}

// ============================================================================
// Hessian Update Methods
// ============================================================================

void ANCOptimizer::updateHessianBFGS(Matrix& hessian, const Vector& dx, const Vector& dg) {
    /*
     * BFGS Hessian update
     * H_new = H + (dg*dg^T)/(dg^T*dx) - (H*dx)*(H*dx)^T/(dx^T*H*dx)
     *
     * Ported from XTB bfgs.f90
     */

    double dg_dot_dx = dg.dot(dx);
    // Claude Generated (Mar 2026): BFGS requires positive curvature (dg^T*dx > 0).
    // Using abs() was wrong — a negative value means the curvature condition is violated
    // and the update would make the Hessian indefinite. Skip in that case.
    if (dg_dot_dx < 1e-10) {
        return; // Skip update if curvature condition fails (must be strictly positive)
    }

    Vector H_dx = hessian * dx;
    double dx_H_dx = dx.dot(H_dx);

    if (dx_H_dx < 1e-10) {
        return; // Avoid division by zero or indefinite update
    }

    // BFGS update
    hessian += (dg * dg.transpose()) / dg_dot_dx;
    hessian -= (H_dx * H_dx.transpose()) / dx_H_dx;

    // XTB: clamp diagonal elements to minimum hlow after update.
    // Prevents the Hessian from becoming too soft, which would cause huge RF steps.
    // Reference: XTB bfgs.f90 Zeile 96-101: if(abs(hess(ij)).lt.thr) hess(ij)=thr  (thr=0.01)
    for (int i = 0; i < hessian.rows(); ++i) {
        if (hessian(i, i) < m_anc->hlow)
            hessian(i, i) = m_anc->hlow;
    }
}

void ANCOptimizer::updateHessianPowell(Matrix& hessian, const Vector& dx, const Vector& dg) {
    /*
     * Powell (Symmetric Broyden) Hessian update
     *
     * Ported from XTB broyden.f90
     */

    Vector y = dg - hessian * dx;
    double dx_dot_y = dx.dot(y);

    if (std::abs(dx_dot_y) < 1e-10) {
        return; // Skip if denominator too small
    }

    // Powell update: H += (y * dx^T + dx * y^T) / (dx^T * dx) - (dx^T * y) / (dx^T * dx)^2 * dx * dx^T
    double dx_norm_sq = dx.squaredNorm();
    if (dx_norm_sq < 1e-10) {
        return;
    }

    Matrix update = (y * dx.transpose() + dx * y.transpose()) / dx_norm_sq;
    update -= (dx_dot_y / (dx_norm_sq * dx_norm_sq)) * (dx * dx.transpose());

    hessian += update;
}

// ============================================================================
// Helper Methods
// ============================================================================

double ANCOptimizer::predictEnergyChange(const Vector& gradient, const Vector& displacement,
                                         const Matrix& hessian) {
    /*
     * Predict energy change from 2nd order model:
     * ΔE = g·dx + 0.5·dx·H·dx
     *
     * Ported from XTB prdechng()
     */

    double linear_term = gradient.dot(displacement);
    double quadratic_term = 0.5 * displacement.dot(hessian * displacement);

    return linear_term + quadratic_term;
}

void ANCOptimizer::projectTranslationsRotations(Matrix& hessian, const Molecule& mol) {
    /*
     * Bug 5 Fix: Eckart T/R projection from Hessian
     * Builds 6 (5 for linear) translation/rotation vectors, orthonormalizes them
     * via Gram-Schmidt, and projects them out: H_proj = P * H * P^T
     * where P = I - V * V^T (V columns = orthonormal T/R vectors)
     *
     * Reference: XTB trproj() and detrotra8 subroutines
     */

    int n3 = 3 * mol.AtomCount();
    const auto& geom = mol.getGeometry();

    // Compute center of mass
    Eigen::Vector3d com = Eigen::Vector3d::Zero();
    double total_mass = 0.0;
    for (int i = 0; i < mol.AtomCount(); ++i) {
        double m = Elements::AtomicMass[mol.Atom(i).first];
        com += m * geom.row(i).transpose();
        total_mass += m;
    }
    com /= total_mass;

    // Build 6 T/R vectors in 3N space
    std::vector<Vector> tr_vecs;
    tr_vecs.reserve(6);

    // 3 Translations: T_d has component 1 at each atom's d-th coordinate
    for (int d = 0; d < 3; ++d) {
        Vector v = Vector::Zero(n3);
        for (int i = 0; i < mol.AtomCount(); ++i)
            v(3 * i + d) = 1.0;
        v.normalize();
        tr_vecs.push_back(v);
    }

    // 3 Rotations about x, y, z axes (cross product of axis with position)
    for (int d = 0; d < 3; ++d) {
        Vector v = Vector::Zero(n3);
        for (int i = 0; i < mol.AtomCount(); ++i) {
            Eigen::Vector3d r = geom.row(i).transpose() - com;
            // Rotation d: each atom contributes e_d x r (cross product)
            if (d == 0) { v(3 * i + 1) = -r(2); v(3 * i + 2) = r(1); }  // Rx: [0,-z,y]
            if (d == 1) { v(3 * i + 0) = r(2);  v(3 * i + 2) = -r(0); } // Ry: [z,0,-x]
            if (d == 2) { v(3 * i + 0) = -r(1); v(3 * i + 1) = r(0); }  // Rz: [-y,x,0]
        }
        if (v.norm() > 1e-10) {
            v.normalize();
            tr_vecs.push_back(v);
        }
    }

    // Gram-Schmidt orthonormalization
    std::vector<Vector> ortho_vecs;
    ortho_vecs.reserve(6);
    for (auto& v : tr_vecs) {
        Vector u = v;
        for (const auto& w : ortho_vecs)
            u -= w.dot(u) * w;
        double norm = u.norm();
        if (norm > 1e-10) {
            u /= norm;
            ortho_vecs.push_back(u);
        }
    }

    // Apply H_proj = P*H*P^T implicitly where P = I - V*V^T.
    // Expansion: P*H*P^T = H - V*(V^T*H) - (H*V)*V^T + V*(V^T*H*V)*V^T.
    // All products involve the thin V (n3 x k, k <= 6), so cost is O(n3^2 * k)
    // instead of O(n3^3) for the dense P*H*P^T multiply.
    const int k = static_cast<int>(ortho_vecs.size());
    if (k > 0) {
        Matrix V(n3, k);
        for (int j = 0; j < k; ++j) V.col(j) = ortho_vecs[j];

        Matrix HV = hessian * V;            // n3 x k
        Matrix VtH = V.transpose() * hessian; // k x n3
        Matrix VtHV = V.transpose() * HV;   // k x k

        hessian.noalias() -= V * VtH;
        hessian.noalias() -= HV * V.transpose();
        hessian.noalias() += V * (VtHV * V.transpose());
    }
}

void ANCOptimizer::setOptimizationLevel(int level) {
    /*
     * Set optimization thresholds based on level
     * Ported from XTB get_optthr()
     */

    // Gradient thresholds are stored in Eh/Ang (Curcuma's internal unit after the
    // Bohr→Ang gradient scaling in CalculateOptimizationStep).
    // Conversion: X Eh/Bohr = X × ANGSTROM_TO_BOHR Eh/Ang  (A2B ≈ 1.88973).
    // Reference values from XTB get_optthr() in Eh/Bohr; multiply by A2B for Eh/Ang.
    const double A2B = CurcumaUnit::Length::ANGSTROM_TO_BOHR;
    switch (level) {
    case -3: // crude
        m_ethr = 5e-4;
        m_gthr = 1e-2  * A2B;  // 1.890e-2 Eh/Ang
        m_acc = 3.0;
        break;
    case -2: // sloppy
        m_ethr = 1e-4;
        m_gthr = 6e-3  * A2B;  // 1.134e-2 Eh/Ang
        m_acc = 3.0;
        break;
    case -1: // loose
        m_ethr = 5e-5;
        m_gthr = 4e-3  * A2B;  // 7.559e-3 Eh/Ang
        m_acc = 2.0;
        break;
    case 0: // normal (default)
        m_ethr = 5e-6;
        m_gthr = 1e-3  * A2B;  // 1.890e-3 Eh/Ang
        m_acc = 1.0;
        break;
    case 1: // tight
        m_ethr = 1e-6;
        m_gthr = 8e-4  * A2B;  // 1.512e-3 Eh/Ang
        m_acc = 0.2;
        break;
    case 2: // vtight
        m_ethr = 1e-7;
        m_gthr = 2e-4  * A2B;  // 3.779e-4 Eh/Ang
        m_acc = 0.05;
        break;
    case 3: // extreme
        m_ethr = 5e-8;
        m_gthr = 5e-5  * A2B;  // 9.449e-5 Eh/Ang
        m_acc = 0.01;
        break;
    default:
        m_ethr = 5e-6;
        m_gthr = 1e-3  * A2B;  // 1.890e-3 Eh/Ang
        m_acc = 1.0;
    }

    CurcumaLogger::param("Optimization level", level);
    CurcumaLogger::param("Energy threshold", std::to_string(m_ethr) + " Eh");
    CurcumaLogger::param("Gradient threshold", std::to_string(m_gthr) + " Eh/Ang");
}

void ANCOptimizer::loadANCParameters(const json& config) {
    if (config.contains("maxdispl")) {
        m_maxdispl = config["maxdispl"].get<double>();
    }
    if (config.contains("hlow")) {
        m_hlow = config["hlow"].get<double>();
    }
    if (config.contains("hmax")) {
        m_hmax = config["hmax"].get<double>();
    }
    if (config.contains("maxmicro")) {
        m_maxmicro = config["maxmicro"].get<int>();
        m_maxmicro_initial = m_maxmicro; // base for ramp-up cap (XTB optimizer.f90:407)
    }
    if (config.contains("hessian_update")) {
        std::string update_str = config["hessian_update"].get<std::string>();
        if (update_str == "powell") {
            m_hessian_update = POWELL;
        } else {
            m_hessian_update = BFGS;
        }
    }

    m_model_hess_params = ModelHessianParameters::fromJson(config);
}

bool ANCOptimizer::checkLinearMolecule(const Molecule& mol) const {
    // Check if molecule is linear by computing moment of inertia
    // If smallest moment is < 1e-10, molecule is linear

    if (mol.AtomCount() < 3) {
        return true; // Atoms and diatomics are considered linear
    }

    // Compute center of mass
    Eigen::Vector3d com = Eigen::Vector3d::Zero();
    double total_mass = 0.0;
    const auto& geom = mol.getGeometry();

    for (int i = 0; i < mol.AtomCount(); ++i) {
        double mass = Elements::AtomicMass[mol.Atom(i).first];
        com += mass * geom.row(i).transpose();
        total_mass += mass;
    }
    com /= total_mass;

    // Compute inertia tensor (renamed from I to avoid conflict with C99 complex macro _Complex_I).
    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();
    for (int i = 0; i < mol.AtomCount(); ++i) {
        double mass = Elements::AtomicMass[mol.Atom(i).first];
        Eigen::Vector3d r = geom.row(i).transpose() - com;

        inertia(0, 0) += mass * (r(1) * r(1) + r(2) * r(2));
        inertia(1, 1) += mass * (r(0) * r(0) + r(2) * r(2));
        inertia(2, 2) += mass * (r(0) * r(0) + r(1) * r(1));
        inertia(0, 1) -= mass * r(0) * r(1);
        inertia(0, 2) -= mass * r(0) * r(2);
        inertia(1, 2) -= mass * r(1) * r(2);
    }
    inertia(1, 0) = inertia(0, 1);
    inertia(2, 0) = inertia(0, 2);
    inertia(2, 1) = inertia(1, 2);

    // Diagonalize
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertia);
    const Eigen::Vector3d& eigenvalues = solver.eigenvalues();

    // Check smallest eigenvalue
    double min_eigenvalue = eigenvalues.minCoeff();
    bool is_linear = (min_eigenvalue < 1e-10);

    if (is_linear) {
        CurcumaLogger::info("Molecule is linear");
    }

    return is_linear;
}

json ANCOptimizer::GetDefaultConfiguration() const {
    // Merge OptimizerDriver base config with ANCOpt-specific defaults from ParameterRegistry.
    // PARAM macros in ancopt_optimizer.h register all parameters under module "ancopt".
    json config = OptimizerDriver::GetDefaultConfiguration();
    json ancopt_defaults = ParameterRegistry::getInstance().getDefaultJson("ancopt");
    for (auto& [key, val] : ancopt_defaults.items()) {
        if (!config.contains(key))
            config[key] = val;
    }

    // Legacy string keys not in PARAM registry (kept for backward compatibility)
    if (!config.contains("model_hessian")) config["model_hessian"] = "lindh1995";
    if (!config.contains("s6"))           config["s6"] = 30.0;
    if (!config.contains("hessian_update")) config["hessian_update"] = "bfgs";

    return config;
}

} // namespace Optimization
