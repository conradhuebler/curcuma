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

#include "ancopt_optimizer.h"
#include "src/core/curcuma_logger.h"
#include "src/core/molecule.h"
#include "src/tools/geometry.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <iostream>

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

    // Select nvar eigenvectors as basis (skip rotations/translations)
    int skip = is_linear ? 5 : 6; // 3 translations + 2/3 rotations
    int start_idx = skip;

    // Check if we have enough eigenvectors
    if (n3 - skip < nvar) {
        CurcumaLogger::warn("Not enough eigenvectors after projection, using all available");
        start_idx = std::max(0, n3 - nvar);
    }

    // Copy selected eigenvectors to B matrix (raw eigenvectors, no scaling)
    // XTB uses eigenvectors directly as ANC basis without frequency-based scaling
    for (int i = 0; i < nvar && (start_idx + i) < n3; ++i) {
        B.col(i) = eigenvectors.col(start_idx + i);
    }

    // Set internal Hessian as diagonal of (damped) eigenvalues - XTB approach
    // This gives a physically meaningful starting Hessian for the RF step
    hess = Matrix::Zero(nvar, nvar);
    for (int i = 0; i < nvar; ++i) {
        int idx = start_idx + i;
        if (idx < evals.size()) {
            double ev = std::max(evals(idx), hlow); // Clamp to hlow minimum
            hess(i, i) = ev;
        }
    }

    // Initialize internal coordinates from Cartesian
    setCartesian(xyz);

    CurcumaLogger::success("ANC generated with " + std::to_string(nvar) + " internal coordinates");

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
    // Transform Cartesian gradient to internal coordinates
    m_gint_old = m_gint;
    m_gint = m_anc->transformGradientToInternal(gradient);
    m_gnorm_old = m_gnorm;
    m_gnorm = m_gint.norm();

    // Check if we need to regenerate ANC
    if (m_needs_anc_regeneration || m_micro_current >= m_maxmicro) {
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

        // Critical: re-project gradient onto NEW ANC basis after regeneration.
        // m_gint was computed above with the OLD B matrix; m_anc->hess now uses the
        // NEW eigenvectors, so mixing them in the RF step would give a garbage direction.
        // Re-projecting here ensures m_gint and m_anc->hess are in the same basis.
        m_gint = m_anc->transformGradientToInternal(gradient);
        m_gnorm = m_gint.norm();
        m_gint_old = m_gint;                               // no valid old basis to compare
        m_displ = Vector::Zero(m_anc->nvar);               // old displacement is in old basis → invalid
    }

    // Update Hessian with BFGS/Powell if we have previous step
    // Guard m_displ.norm() > 1e-10 naturally skips the first step after ANC regeneration
    // (m_displ was zeroed above) and the very first iteration.
    if (m_current_iteration > 0 && m_displ.norm() > 1e-10) {
        Vector dg = m_gint - m_gint_old;

        if (m_hessian_update == BFGS) {
            updateHessianBFGS(m_anc->hess, m_displ, dg);
        } else {
            updateHessianPowell(m_anc->hess, m_displ, dg);
        }
    }

    // Calculate Rational Function step
    m_displ = calculateRationalFunctionStep(m_gint, m_anc->hess);

    // Claude Generated (Mar 2026): Trust radius step limiting.
    // Limits the total RF step norm to m_trust_radius (Å).
    // This prevents overshooting when the model Hessian is too soft.
    // Each component is also clamped to m_maxdispl (per-component limit from XTB).
    double displ_norm = m_displ.norm();
    if (displ_norm > m_trust_radius) {
        m_displ *= m_trust_radius / displ_norm;  // Scale to trust radius
    }
    // Per-component limit (XTB-style maxdispl, in Å)
    double maxdispl_ang = m_maxdispl * 0.529177; // convert Bohr → Å
    for (int i = 0; i < m_displ.size(); ++i) {
        if (std::abs(m_displ(i)) > maxdispl_ang) {
            m_displ(i) = std::copysign(maxdispl_ang, m_displ(i));
        }
    }

    // Predict energy change
    m_depred = predictEnergyChange(m_gint, m_displ, m_anc->hess);

    // Bug 2 Fix: Adaptive step size scaling (XTB optimizer.f90 Z. 720-728)
    // XTB scales the step when gradient norm is small (near convergence)
    double alp = 1.0;
    if (m_gnorm < 0.0003) alp = 3.0;
    else if (m_gnorm < 0.0006) alp = 2.0;
    else if (m_gnorm < 0.002) alp = 1.5;

    // Update internal coordinates with scaled displacement
    m_anc->coord += m_displ * alp;

    // Transform back to Cartesian (returns absolute new coordinates)
    Vector new_xyz;
    m_anc->getCartesian(new_xyz);

    m_micro_current++;

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
    // XTB-style convergence: gradient AND energy conditions.
    // Reference: XTB optimizer.f90: converged = gconverged AND econverged
    bool grad_converged = m_gnorm < m_gthr;

    bool energy_ok;
    if (m_current_iteration <= 1) {
        energy_ok = true;  // First iteration: no previous energy to compare
    } else if (std::abs(m_energy_change) < m_ethr) {
        energy_ok = true;  // Within energy threshold (allows tiny fluctuations near minimum)
    } else {
        energy_ok = (m_energy_change < 0.0);  // Must be decreasing if not within threshold
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info_fmt("Convergence: gnorm={:.3e} (thr={:.3e}), dE={:.3e} (thr={:.3e})",
            m_gnorm, m_gthr, m_energy_change, m_ethr);
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

    // Claude Generated (Mar 2026): Trust radius adaptation.
    // Standard Fletcher-Reeves-style update based on actual/predicted ratio.
    // Reference: Nocedal & Wright "Numerical Optimization", trust region methods.
    if (m_energy_change > 0.0) {
        // Energy rose: model Hessian too soft, step too large.
        // Reduce trust radius aggressively and regenerate ANC from current geometry.
        m_trust_radius *= 0.25;
        m_trust_radius = std::max(m_trust_radius, 0.005 * 0.529177); // min ~0.003 Å
        m_needs_anc_regeneration = true;
        CurcumaLogger::warn_fmt("Energy rose, trust radius reduced to {:.4f} Å", m_trust_radius);
    } else if (m_depred < 0.0) {
        // Energy decreased. Update trust radius based on quality of prediction.
        double ratio = m_energy_change / m_depred; // ratio > 0 (both negative)
        double maxtr = m_maxdispl * 0.529177;       // max trust radius in Å
        if (ratio > 0.75) {
            // Excellent prediction: increase trust radius
            m_trust_radius = std::min(m_trust_radius * 2.0, maxtr);
        } else if (ratio < 0.25) {
            // Poor prediction (large overshoot): decrease trust radius
            m_trust_radius *= 0.5;
            m_trust_radius = std::max(m_trust_radius, 0.005 * 0.529177);
        }
        // else (0.25 ≤ ratio ≤ 0.75): keep trust radius unchanged
    }

    CurcumaLogger::param("Trust radius (Å)", fmt::format("{:.4f}", m_trust_radius));

    // Store current state
    m_current_energy = new_energy;
    m_current_gradient = new_gradient;
}

void ANCOptimizer::FinalizeOptimizationInternal() {
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

Vector ANCOptimizer::calculateRationalFunctionStep(const Vector& gradient_internal,
                                                    const Matrix& hessian_internal) {
    /*
     * Rational Function Optimization (RFO) step
     *
     * Solves augmented eigenvalue problem:
     * | H  g | * | dx |     = λ * | dx |
     * | g  0 |   | 1  |           | 1  |
     *
     * Ported from XTB optimizer.f90:676-729
     */

    int nvar = gradient_internal.size();
    int nvar1 = nvar + 1;

    // Construct augmented Hessian
    Matrix A_aug = Matrix::Zero(nvar1, nvar1);
    A_aug.topLeftCorner(nvar, nvar) = hessian_internal;
    A_aug.block(nvar, 0, 1, nvar) = gradient_internal.transpose();
    A_aug.block(0, nvar, nvar, 1) = gradient_internal;
    A_aug(nvar, nvar) = 0.0;

    // Solve eigenvalue problem
    Eigen::SelfAdjointEigenSolver<Matrix> solver(A_aug);
    if (solver.info() != Eigen::Success) {
        CurcumaLogger::error("RF solver failed - using steepest descent");
        return -gradient_internal * 0.1; // Fallback
    }

    // Get lowest eigenvalue and corresponding eigenvector
    // XTB takes the first (lowest/most negative) eigenvector for the downhill step
    // Reference: XTB optimizer.f90 - lambda = eigv(1)
    const Vector& eigenvalues = solver.eigenvalues();
    const Matrix& eigenvectors = solver.eigenvectors();

    // Use eigenvector with lowest eigenvalue (index 0 from SelfAdjointEigenSolver)
    Vector eigenvec = eigenvectors.col(0);

    // Extract displacement: eigenvec = [p; s], displacement = p/s
    // Ensure s > 0 (positive normalization) to get correct direction
    double last_elem = eigenvec(nvar);
    if (std::abs(last_elem) < 1e-10) {
        CurcumaLogger::warn("RF normalization failed - using steepest descent");
        return -gradient_internal * 0.1;
    }
    // Flip sign if s < 0 to ensure consistent direction
    if (last_elem < 0.0)
        eigenvec = -eigenvec;
    last_elem = eigenvec(nvar);

    Vector displacement = eigenvec.head(nvar) / last_elem;

    CurcumaLogger::param("RF eigenvalue", std::to_string(eigenvalues(0)));

    return displacement;
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

    // Build projection: P = I - V*V^T, apply H_proj = P * H * P^T
    Matrix P = Matrix::Identity(n3, n3);
    for (const auto& v : ortho_vecs)
        P -= v * v.transpose();

    hessian = P * hessian * P.transpose();
}

void ANCOptimizer::setOptimizationLevel(int level) {
    /*
     * Set optimization thresholds based on level
     * Ported from XTB get_optthr()
     */

    switch (level) {
    case -3: // crude
        m_ethr = 5e-4;
        m_gthr = 1e-2;
        m_acc = 3.0;
        break;
    case -2: // sloppy
        m_ethr = 1e-4;
        m_gthr = 6e-3;
        m_acc = 3.0;
        break;
    case -1: // loose
        m_ethr = 5e-5;
        m_gthr = 4e-3;
        m_acc = 2.0;
        break;
    case 0: // normal (default)
        m_ethr = 5e-6;
        m_gthr = 1e-3;
        m_acc = 1.0;
        break;
    case 1: // tight
        m_ethr = 1e-6;
        m_gthr = 8e-4;
        m_acc = 0.2;
        break;
    case 2: // vtight
        m_ethr = 1e-7;
        m_gthr = 2e-4;
        m_acc = 0.05;
        break;
    case 3: // extreme
        m_ethr = 5e-8;
        m_gthr = 5e-5;
        m_acc = 0.01;
        break;
    default:
        m_ethr = 5e-6;
        m_gthr = 1e-3;
        m_acc = 1.0;
    }

    CurcumaLogger::param("Optimization level", level);
    CurcumaLogger::param("Energy threshold", std::to_string(m_ethr) + " Eh");
    CurcumaLogger::param("Gradient threshold", std::to_string(m_gthr) + " Eh/Bohr");
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

    // Compute inertia tensor
    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    for (int i = 0; i < mol.AtomCount(); ++i) {
        double mass = Elements::AtomicMass[mol.Atom(i).first];
        Eigen::Vector3d r = geom.row(i).transpose() - com;

        I(0, 0) += mass * (r(1) * r(1) + r(2) * r(2));
        I(1, 1) += mass * (r(0) * r(0) + r(2) * r(2));
        I(2, 2) += mass * (r(0) * r(0) + r(1) * r(1));
        I(0, 1) -= mass * r(0) * r(1);
        I(0, 2) -= mass * r(0) * r(2);
        I(1, 2) -= mass * r(1) * r(2);
    }
    I(1, 0) = I(0, 1);
    I(2, 0) = I(0, 2);
    I(2, 1) = I(1, 2);

    // Diagonalize
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I);
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
    json config = OptimizerDriver::GetDefaultConfiguration();

    config["maxdispl"] = 0.3;
    config["hlow"] = 0.01;
    config["hmax"] = 5.0;
    config["maxmicro"] = 25;
    config["model_hessian"] = "lindh2007";
    config["s6"] = 30.0;
    config["hessian_update"] = "bfgs";
    config["optimization_level"] = 0; // normal

    return config;
}

} // namespace Optimization
