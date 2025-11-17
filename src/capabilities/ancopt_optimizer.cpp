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

    const Vector& eigenvalues = solver.eigenvalues();
    const Matrix& eigenvectors = solver.eigenvectors();

    // Select nvar eigenvectors as basis (skip rotations/translations)
    int skip = is_linear ? 5 : 6; // 3 translations + 2/3 rotations
    int start_idx = skip;

    // Check if we have enough eigenvectors
    if (n3 - skip < nvar) {
        CurcumaLogger::warn("Not enough eigenvectors after projection, using all available");
        start_idx = std::max(0, n3 - nvar);
    }

    // Copy selected eigenvectors to B matrix (transformation matrix)
    for (int i = 0; i < nvar && (start_idx + i) < n3; ++i) {
        B.col(i) = eigenvectors.col(start_idx + i);

        // Scale by frequency (optional, following XTB approach)
        double freq = eigenvalues(start_idx + i);
        if (freq > hlow && freq < hmax) {
            // Keep as is - within desired frequency range
        } else if (freq < hlow && freq > 0.0) {
            // Scale down low frequencies
            B.col(i) *= std::sqrt(freq / hlow);
        } else if (freq > hmax) {
            // Scale down high frequencies
            B.col(i) *= std::sqrt(hmax / freq);
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
    std::cerr << "[DEBUG ANCOptimizer::InitializeOptimizerInternal] m_molecule.AtomCount() = " << atom_count << std::endl;

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
    std::cerr << "[DEBUG] About to loadANCParameters" << std::endl;
    loadANCParameters(m_configuration);
    std::cerr << "[DEBUG] Done loadANCParameters" << std::endl;

    // Set optimization level (affects thresholds)
    int opt_level = m_configuration.value("optimization_level", 0); // 0 = normal
    setOptimizationLevel(opt_level);

    // Check if molecule is linear
    std::cerr << "[DEBUG] About to checkLinearMolecule" << std::endl;
    bool is_linear = checkLinearMolecule(m_molecule);
    std::cerr << "[DEBUG] is_linear = " << is_linear << std::endl;

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
    }

    // Update Hessian with BFGS/Powell if we have previous step
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

    // Limit step size
    for (int i = 0; i < m_displ.size(); ++i) {
        if (std::abs(m_displ(i)) > m_maxdispl) {
            m_displ(i) = std::copysign(m_maxdispl, m_displ(i));
        }
    }

    // Predict energy change
    m_depred = predictEnergyChange(m_gint, m_displ, m_anc->hess);

    // Update internal coordinates
    m_anc->coord += m_displ;

    // Transform back to Cartesian
    Vector new_xyz;
    m_anc->getCartesian(new_xyz);

    m_micro_current++;

    CurcumaLogger::param("Displacement norm (ANC)", std::to_string(m_displ.norm()));
    CurcumaLogger::param("Predicted energy change", std::to_string(m_depred) + " Eh");

    return new_xyz;
}

bool ANCOptimizer::CheckMethodSpecificConvergence() const {
    // AncOpt convergence: gradient norm below threshold
    bool grad_converged = m_gnorm < m_gthr;

    CurcumaLogger::info("Gradient norm: " + std::to_string(m_gnorm) + " Eh/Bohr (threshold: " + std::to_string(m_gthr) + ")");

    return grad_converged;
}

void ANCOptimizer::UpdateOptimizerState(const Vector& new_coordinates,
                                        const Vector& new_gradient,
                                        double new_energy) {
    // Check if step was too large (displacement norm > 2.0 in ANC space)
    if (m_displ.norm() > 2.0) {
        CurcumaLogger::warn("Large displacement detected - will regenerate ANC");
        m_needs_anc_regeneration = true;
    }

    // Store current state
    m_current_energy = new_energy;
    m_current_gradient = new_gradient;
}

void ANCOptimizer::FinalizeOptimizationInternal() {
    // Cleanup
    if (m_anc) {
        m_anc->deallocate();
    }

    CurcumaLogger::success("AncOpt optimization finalized");
}

// ============================================================================
// Model Hessian Generation
// ============================================================================

Matrix ANCOptimizer::generateModelHessian(const Molecule& mol) {
    CurcumaLogger::info("Generating model Hessian (Lindh 2007)");

    int n3 = 3 * mol.AtomCount();
    Matrix hessian = Matrix::Zero(n3, n3);

    // Simplified Lindh model Hessian
    // For now, use simple diagonal approximation based on atomic masses
    // TODO: Implement full Lindh (2007) model with bond/angle/dihedral terms

    const auto& geometry = mol.getGeometry();

    // Diagonal elements: simple harmonic approximation
    for (int i = 0; i < mol.AtomCount(); ++i) {
        double mass = Elements::AtomicMass[mol.Atom(i).first];
        double k_diagonal = m_model_hess_params.s6 / std::sqrt(mass);

        for (int xyz = 0; xyz < 3; ++xyz) {
            hessian(3 * i + xyz, 3 * i + xyz) = k_diagonal;
        }
    }

    // Off-diagonal elements: bond connectivity
    auto bonds = mol.Bonds();
    for (const auto& bond : bonds) {
        int i = bond.first;
        int j = bond.second;

        if (i >= mol.AtomCount() || j >= mol.AtomCount())
            continue;

        double rij = (geometry.row(i) - geometry.row(j)).norm();
        double mass_i = Elements::AtomicMass[mol.Atom(i).first];
        double mass_j = Elements::AtomicMass[mol.Atom(j).first];

        double k_bond = m_model_hess_params.s6 / (rij * std::sqrt(mass_i * mass_j));

        // Simplified: only diagonal blocks
        for (int xyz = 0; xyz < 3; ++xyz) {
            int idx_i = 3 * i + xyz;
            int idx_j = 3 * j + xyz;

            hessian(idx_i, idx_j) -= k_bond * 0.5;
            hessian(idx_j, idx_i) -= k_bond * 0.5;
        }
    }

    CurcumaLogger::success("Model Hessian generated");

    return hessian;
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
    const Vector& eigenvalues = solver.eigenvalues();
    const Matrix& eigenvectors = solver.eigenvectors();

    // Find eigenvalue closest to zero (but not exactly zero)
    int min_idx = 0;
    double min_abs = std::abs(eigenvalues(0));
    for (int i = 1; i < nvar1; ++i) {
        if (std::abs(eigenvalues(i)) < min_abs && std::abs(eigenvalues(i)) > 1e-10) {
            min_abs = std::abs(eigenvalues(i));
            min_idx = i;
        }
    }

    Vector eigenvec = eigenvectors.col(min_idx);

    // Extract displacement (divide by last element to normalize)
    double last_elem = eigenvec(nvar);
    if (std::abs(last_elem) < 1e-10) {
        CurcumaLogger::warn("RF normalization failed - using steepest descent");
        return -gradient_internal * 0.1;
    }

    Vector displacement = eigenvec.head(nvar) / last_elem;

    CurcumaLogger::param("RF eigenvalue", std::to_string(eigenvalues(min_idx)));

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
    if (std::abs(dg_dot_dx) < 1e-10) {
        return; // Skip update if curvature condition fails
    }

    Vector H_dx = hessian * dx;
    double dx_H_dx = dx.dot(H_dx);

    if (std::abs(dx_H_dx) < 1e-10) {
        return; // Avoid division by zero
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
     * Project out translations and rotations from Hessian
     * This is a simplified version - full implementation needs proper projection matrix
     *
     * Ported from XTB trproj()
     */

    // TODO: Implement full Eckart projection
    // For now, this is a placeholder that does nothing
    // The ANC generation will handle projection through eigenvector selection
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
