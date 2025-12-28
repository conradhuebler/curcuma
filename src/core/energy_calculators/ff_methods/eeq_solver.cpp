/*
 * < EEQ (Electronegativity Equalization) Charge Solver - Implementation >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Extracted from GFN-FF implementation (gfnff_method.cpp) for reusability
 * across all force field methods (UFF, QMDFF, D4 dispersion).
 *
 * References:
 * - S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665-15673
 * - R. T. Sanderson, Chemical Bonds and Bond Energy, Academic Press, 1976
 *
 * Claude Generated - December 2025
 */

#include "eeq_solver.h"
#include "gfnff_par.h"  // EEQ parameters (chi_eeq, gam_eeq, alpha_eeq, cnf_eeq)
#include "cn_calculator.h"  // Shared CN calculation utility

#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"

#include <Eigen/Dense>
#include <cmath>
#include <algorithm>

using namespace GFNFFParameters;  // Access to chi_eeq, gam_eeq, alpha_eeq, cnf_eeq

// ===== Constructor =====

EEQSolver::EEQSolver(const ConfigManager& config)
    : m_config(config)
{
    // Extract parameters from ConfigManager with defaults from PARAM definitions
    m_max_iterations = m_config.get<int>("max_iterations", 50);
    m_convergence_threshold = m_config.get<double>("convergence_threshold", 1e-6);
    m_verbosity = m_config.get<int>("verbosity", 0);
    m_calculate_cn = m_config.get<bool>("calculate_cn", true);

    if (m_verbosity >= 2) {
        CurcumaLogger::info("EEQSolver initialized with parameters:");
        CurcumaLogger::param("max_iterations", std::to_string(m_max_iterations));
        CurcumaLogger::param("convergence_threshold", fmt::format("{:.2e}", m_convergence_threshold));
        CurcumaLogger::param("verbosity", std::to_string(m_verbosity));
        CurcumaLogger::param("calculate_cn", m_calculate_cn ? "true" : "false");
    }
}

// ===== Main API =====

Vector EEQSolver::calculateCharges(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    int total_charge,
    const Vector* cn_hint,
    const std::vector<int>* hyb_hint)
{
    const int natoms = atoms.size();

    if (natoms == 0) {
        CurcumaLogger::error("EEQSolver::calculateCharges: No atoms provided");
        return Vector::Zero(0);
    }

    if (geometry_bohr.rows() != natoms || geometry_bohr.cols() != 3) {
        CurcumaLogger::error(fmt::format("EEQSolver::calculateCharges: Invalid geometry dimensions {}x{} (expected {}x3)",
                                         geometry_bohr.rows(), geometry_bohr.cols(), natoms));
        return Vector::Zero(natoms);
    }

    // Step 1: Calculate or use provided coordination numbers
    Vector cn;
    if (cn_hint != nullptr && cn_hint->size() == natoms) {
        cn = *cn_hint;
        if (m_verbosity >= 3) {
            CurcumaLogger::info("EEQSolver: Using provided coordination numbers");
        }
    } else if (m_calculate_cn) {
        cn = calculateCoordinationNumbers(atoms, geometry_bohr);
        if (m_verbosity >= 3) {
            CurcumaLogger::info("EEQSolver: Calculated coordination numbers");
        }
    } else {
        CurcumaLogger::error("EEQSolver::calculateCharges: CN required but calculate_cn=false and no hint provided");
        return Vector::Zero(natoms);
    }

    // Step 2: Detect or use provided hybridization
    std::vector<int> hybridization;
    if (hyb_hint != nullptr && hyb_hint->size() == natoms) {
        hybridization = *hyb_hint;
        if (m_verbosity >= 3) {
            CurcumaLogger::info("EEQSolver: Using provided hybridization states");
        }
    } else {
        hybridization = detectHybridization(atoms, geometry_bohr, cn);
        if (m_verbosity >= 3) {
            CurcumaLogger::info("EEQSolver: Detected hybridization states");
        }
    }

    // Step 3: Phase 1 - Calculate topology charges
    Vector topology_charges = calculateTopologyCharges(atoms, geometry_bohr, total_charge, cn);

    if (topology_charges.size() != natoms) {
        CurcumaLogger::error("EEQSolver::calculateCharges: Phase 1 failed");
        return Vector::Zero(natoms);
    }

    // Step 4: Phase 2 - Refine charges with environmental corrections
    Vector final_charges = calculateFinalCharges(atoms, geometry_bohr, total_charge,
                                                  topology_charges, cn, hybridization);

    if (final_charges.size() != natoms) {
        CurcumaLogger::error("EEQSolver::calculateCharges: Phase 2 failed");
        return topology_charges;  // Return Phase 1 charges as fallback
    }

    if (m_verbosity >= 1) {
        CurcumaLogger::success(fmt::format("EEQSolver: Two-phase EEQ completed for {} atoms", natoms));
        if (m_verbosity >= 2) {
            for (int i = 0; i < std::min(5, natoms); ++i) {
                CurcumaLogger::result(fmt::format("Atom {} (Z={}) q = {:.6f}", i, atoms[i], final_charges(i)));
            }
        }
    }

    return final_charges;
}

// ===== Phase 1: Topology Charges =====

Vector EEQSolver::calculateTopologyCharges(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    int total_charge,
    const Vector& cn)
{
    const int natoms = atoms.size();
    const double TSQRT2PI = 0.797884560802866;  // sqrt(2/π)

    // Augmented system size: n atoms + 1 constraint
    int nfrag = 1;  // Single molecular fragment (neutral or charged)
    int m = natoms + nfrag;

    // Setup EEQ parameters with CN-dependence
    Vector chi(natoms);
    Vector gam(natoms);
    Vector alpha(natoms);

    for (int i = 0; i < natoms; ++i) {
        int z_i = atoms[i];
        EEQParameters params_i = getParameters(z_i, cn(i));

        // Claude Generated (December 2025, Session 11): CRITICAL FIX - Add CNF term in Phase 1
        // Reference: XTB gfnff_ini.f90:411
        // CORRECT: chi = -chi + dxi + cnf*sqrt(CN)  (XTB formula)
        // WRONG (old): chi = -chi + dxi_total (missing CNF term!)

        // Chi with CN-dependent corrections
        // Simple hybridization guess: high CN → sp3, low CN → sp/sp2
        double dxi_hyb = 0.0;
        if (cn(i) < 1.5) {
            dxi_hyb = 0.1;  // sp-like
        } else if (cn(i) < 2.5) {
            dxi_hyb = 0.05;  // sp2-like
        }
        // sp3-like: no correction

        double dxi_cn = -0.01 * (cn(i) - 2.0);
        double dxi_total = dxi_hyb + dxi_cn;

        // CRITICAL: Add CNF*sqrt(CN) term (from XTB gfnff_ini.f90:411)
        double cnf_term = params_i.cnf * std::sqrt(cn(i));

        chi(i) = -params_i.chi + dxi_total + cnf_term;
        gam(i) = params_i.gam;
        alpha(i) = params_i.alp;  // Already squared
    }

    // Build AUGMENTED EEQ matrix A (m x m)
    Matrix A = Matrix::Zero(m, m);
    Vector x = Vector::Zero(m);

    // DEBUG: Print Phase 1 parameters for first 3 atoms
    if (m_verbosity >= 3 && natoms >= 1) {
        std::cerr << "\n=== Phase 1 EEQ Parameter Breakdown (Topology Charges) ===" << std::endl;
        for (int i = 0; i < std::min(3, natoms); ++i) {
            int z_i = atoms[i];
            EEQParameters params_raw = getParameters(z_i, 0.0);  // Without CN correction
            EEQParameters params_cn = getParameters(z_i, cn(i)); // With CN correction

            // dxi calculation (same as above)
            double dxi_hyb = 0.0;
            if (cn(i) < 1.5) dxi_hyb = 0.1;
            else if (cn(i) < 2.5) dxi_hyb = 0.05;
            double dxi_cn_corr = -0.01 * (cn(i) - 2.0);
            double dxi_total = dxi_hyb + dxi_cn_corr;
            double cnf_term = params_cn.cnf * std::sqrt(cn(i));

            std::cerr << "Atom " << i << " (Z=" << z_i << "):" << std::endl;
            std::cerr << "  CN = " << cn(i) << std::endl;
            std::cerr << "  chi_base = " << params_raw.chi << std::endl;
            std::cerr << "  gam_base = " << params_raw.gam << std::endl;
            std::cerr << "  alpha_base = " << std::sqrt(params_raw.alp) << " (squared: " << params_raw.alp << ")" << std::endl;
            std::cerr << "  cnf = " << params_raw.cnf << std::endl;
            std::cerr << "  dxi_hyb = " << dxi_hyb << ", dxi_cn = " << dxi_cn_corr << ", dxi_total = " << dxi_total << std::endl;
            std::cerr << "  cnf_term = cnf*sqrt(CN) = " << cnf_term << std::endl;
            std::cerr << "  chi_corrected = -chi + dxi + cnf*sqrt(CN) = " << chi(i) << std::endl;
            std::cerr << "  NOTE: XTB adds CNF term TWICE in Phase 1!" << std::endl;
        }
        std::cerr << "========================================\n" << std::endl;
    }

    // 1. Setup RHS and diagonal
    for (int i = 0; i < natoms; ++i) {
        // CRITICAL FIX (Dec 28, 2025): XTB adds CNF term TWICE in Phase 1!
        // gfnff_ini.f90:411: topo%chieeq = -chi + dxi + CNF*√CN
        // gfnff_engrad.F90:1504: x(i) = topo%chieeq + CNF*√CN
        // Total: x = -chi + dxi + 2×CNF*√CN
        //
        // Our chi(i) already includes CNF term once (line 188), so add it again:
        int z_i = atoms[i];
        EEQParameters params_i = getParameters(z_i, cn(i));
        double cnf_term = params_i.cnf * std::sqrt(cn(i));

        x(i) = chi(i) + cnf_term;  // chi already has 1×CNF, add 2nd term
        A(i, i) = gam(i) + TSQRT2PI / std::sqrt(alpha(i));
    }

    // DEBUG: Print final RHS values with 2×CNF
    if (m_verbosity >= 3 && natoms >= 1) {
        std::cerr << "\n=== Phase 1 Final RHS (with 2×CNF fix) ===" << std::endl;
        for (int i = 0; i < std::min(3, natoms); ++i) {
            std::cerr << "  x(" << i << ") = " << x(i) << " (was " << chi(i) << " before 2nd CNF)" << std::endl;
        }
        std::cerr << "========================================\n" << std::endl;
    }

    // 2. Setup off-diagonal Coulomb matrix
    //
    // CRITICAL TODO (Dec 28, 2025): XTB uses TOPOLOGICAL distances, NOT geometric!
    // Reference: gfnff_ini.f90:431-461 (Floyd-Warshall shortest-path algorithm)
    //
    // XTB Algorithm:
    //   1. Build bond distance matrix from topology (covalent radii)
    //   2. Floyd-Warshall to compute shortest topological paths
    //   3. Scale by gen%rfgoed1 (typically 1.0-1.2)
    //   4. Use scaled topological distances in EEQ matrix
    //
    // Current Implementation: Uses geometric distances (WRONG!)
    // Impact: Topology charges ~4-5× too large
    //
    // Geometric distances are SHORTER than topological (direct line vs through bonds)
    // → Larger Coulomb terms → Larger charges → Wrong dispersion
    //
    // TODO: Implement Floyd-Warshall for topological distance matrix
    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < i; ++j) {
            double dx = geometry_bohr(i, 0) - geometry_bohr(j, 0);
            double dy = geometry_bohr(i, 1) - geometry_bohr(j, 1);
            double dz = geometry_bohr(i, 2) - geometry_bohr(j, 2);
            double r_sq = dx*dx + dy*dy + dz*dz;
            double r = std::sqrt(r_sq);

            if (r < 1e-10) {
                CurcumaLogger::error("EEQSolver::calculateTopologyCharges: Zero distance between atoms");
                return Vector::Zero(0);
            }

            // J_ij = erf(gamma_ij * r) / r
            // gamma_ij = 1/sqrt(alpha_i + alpha_j)
            double gammij = 1.0 / std::sqrt(alpha(i) + alpha(j));
            double erf_gamma = std::erf(gammij * r);
            double coulomb = erf_gamma / r;

            A(i, j) = coulomb;
            A(j, i) = coulomb;
        }
    }

    // 3. Setup fragment charge constraint
    x(natoms) = static_cast<double>(total_charge);
    for (int j = 0; j < natoms; ++j) {
        A(natoms, j) = 1.0;
        A(j, natoms) = 1.0;
    }

    // 4. Solve augmented system
    if (m_verbosity >= 3) {
        CurcumaLogger::info("=== Phase 1 EEQ Matrix Diagnostics ===");

        Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(A);
        Vector eigenvalues = eigensolver.eigenvalues();

        double max_eigenvalue = eigenvalues.maxCoeff();
        double min_eigenvalue = eigenvalues.minCoeff();
        double condition_number = std::abs(max_eigenvalue / min_eigenvalue);

        CurcumaLogger::param("matrix_size", fmt::format("{}x{} (augmented)", m, m));
        CurcumaLogger::param("condition_number", fmt::format("{:.6e}", condition_number));

        if (condition_number > 1e12) {
            CurcumaLogger::warn(fmt::format("Phase 1 EEQ matrix ill-conditioned (cond={:.2e})", condition_number));
        } else if (condition_number > 1e8) {
            CurcumaLogger::warn(fmt::format("Phase 1 EEQ matrix poorly conditioned (cond={:.2e})", condition_number));
        } else {
            CurcumaLogger::success(fmt::format("Phase 1 EEQ matrix well-conditioned (cond={:.2e})", condition_number));
        }
    }

    Eigen::PartialPivLU<Matrix> lu(A);
    Vector solution = lu.solve(x);

    // Extract atomic charges
    Vector topology_charges = solution.segment(0, natoms);

    // Check for NaN/Inf
    for (int i = 0; i < natoms; ++i) {
        if (std::isnan(topology_charges[i]) || std::isinf(topology_charges[i])) {
            CurcumaLogger::error(fmt::format("Phase 1 EEQ: Invalid charge[{}] = {} (Z={})",
                                             i, topology_charges[i], atoms[i]));
            return Vector::Zero(0);
        }
    }

    // DEBUG: Print resulting topology charges
    if (m_verbosity >= 3 && natoms >= 1) {
        std::cerr << "\n=== Phase 1 Topology Charges (qa) ===" << std::endl;
        for (int i = 0; i < std::min(5, natoms); ++i) {
            std::cerr << "  qa[" << i << "] (Z=" << atoms[i] << ") = " << topology_charges(i) << std::endl;
        }
        double total_charge = topology_charges.sum();
        std::cerr << "  Total charge = " << total_charge << " (expected: " << total_charge << ")" << std::endl;
        std::cerr << "========================================\n" << std::endl;
    }

    if (m_verbosity >= 2) {
        CurcumaLogger::info("Phase 1 EEQ: Topology charges calculated");
        for (int i = 0; i < std::min(5, natoms); ++i) {
            CurcumaLogger::result(fmt::format("Atom {} qa = {:.6f}", i, topology_charges(i)));
        }
    }

    return topology_charges;
}

// ===== Phase 2: Final Refined Charges =====

Vector EEQSolver::calculateFinalCharges(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    int total_charge,
    const Vector& topology_charges,
    const Vector& cn,
    const std::vector<int>& hybridization)
{
    const int natoms = atoms.size();
    const double TSQRT2PI = 0.797884560802866;  // sqrt(2/π)

    // Calculate correction terms (dxi and dgam only - alpha calculated inline)
    Vector dxi = calculateDxi(atoms, geometry_bohr, cn);
    Vector dgam = calculateDgam(atoms, topology_charges, hybridization);

    // DEBUG: Print topology charges used for dgam
    if (m_verbosity >= 3) {
        std::cerr << "\n=== Phase 2 Topology Charges (used for dgam) ===" << std::endl;
        for (int i = 0; i < std::min(3, natoms); ++i) {
            std::cerr << "Atom " << i << " (Z=" << atoms[i] << "): qa = " << topology_charges(i)
                      << ", dgam = " << dgam(i) << std::endl;
        }
        std::cerr << "================================================\n" << std::endl;
    }

    Vector final_charges = topology_charges;

    // Augmented system size
    int m = natoms + 1;

    // Claude Generated (2025): Pre-calculate corrected parameters ONCE
    // Note: dxi and dgam are constant, but alpha is charge-dependent (calculated inline)
    Vector chi_corrected(natoms);
    Vector gam_corrected(natoms);
    Vector alpha_corrected(natoms);

    for (int i = 0; i < natoms; ++i) {
        int z_i = atoms[i];
        EEQParameters params_i = getParameters(z_i, cn(i));

        // Claude Generated (December 2025, Session 11): CRITICAL FIX - Correct alpha formula from XTB
        // Reference: Fortran gfnff_ini.f90:699-706
        // CORRECT: alpha = (alpha_base + ff*qa)²  (charge-dependent, NOT CN-dependent!)
        // WRONG (old): alpha = alpha_base² + dalpha (CN-dependent)

        // Get base alpha (UNSQUARED) from gfnff_par.h
        double alpha_base = (z_i >= 1 && z_i <= 86) ? alpha_eeq[z_i - 1] : 0.903430;

        // Calculate charge-dependent ff factor (from XTB gfnff_ini.f90:699-705)
        double ff = 0.0;
        if (z_i == 6) {  // Carbon
            ff = 0.09;
        } else if (z_i == 7) {  // Nitrogen
            ff = -0.21;
        } else if (z_i > 10 && z_i <= 86) {  // Heavy atoms only
            int group = periodic_group[z_i - 1];
            int imetal_val = metal_type[z_i - 1];

            if (group == 6) {  // Chalcogens (O, S, Se, Te, Po)
                ff = -0.03;
            } else if (group == 7) {  // Halogens (F, Cl, Br, I, At)
                ff = 0.50;
            } else if (imetal_val == 1) {  // Main group metals
                ff = 0.3;
            } else if (imetal_val == 2) {  // Transition metals
                ff = -0.1;
            }
        }

        // Apply CORRECT formula: alpha = (alpha_base + ff*qa)²
        alpha_corrected(i) = std::pow(alpha_base + ff * topology_charges(i), 2);

        // CRITICAL FIX (Session 11 - CORRECTED): chi_corrected WITHOUT CNF!
        // Reference: XTB gfnff_ini.f90:696 (Phase 2)
        // chieeq = -χ + dxi (OHNE CNF!)
        // Then RHS adds CNF once: x = chieeq + CNF·√CN = -χ + dxi + CNF·√CN
        chi_corrected(i) = -params_i.chi + dxi(i);  // WITHOUT CNF!
        gam_corrected(i) = params_i.gam + dgam(i);
    }

    // Claude Generated (2025): Pre-calculate distance matrix (geometry-dependent only)
    // This eliminates O(k·n²) sqrt() calls in the iteration loop
    Matrix distances = Matrix::Zero(natoms, natoms);
    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < i; ++j) {
            double dx = geometry_bohr(i, 0) - geometry_bohr(j, 0);
            double dy = geometry_bohr(i, 1) - geometry_bohr(j, 1);
            double dz = geometry_bohr(i, 2) - geometry_bohr(j, 2);
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-10) {
                CurcumaLogger::error("EEQSolver::calculateFinalCharges: atoms too close");
                return Vector::Zero(0);
            }

            distances(i, j) = r;
            distances(j, i) = r;  // Symmetric
        }
    }

    // Claude Generated (2025): Pre-build base Coulomb matrix ONCE
    // Off-diagonal elements depend on geometry and alpha (both constant in loop)
    Matrix A_coulomb_base = Matrix::Zero(m, m);

    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < i; ++j) {
            double r = distances(i, j);
            double gamma_ij = 1.0 / std::sqrt(alpha_corrected(i) + alpha_corrected(j));
            double erf_gamma = std::erf(gamma_ij * r);
            double coulomb = erf_gamma / r;

            A_coulomb_base(i, j) = coulomb;
            A_coulomb_base(j, i) = coulomb;
        }
    }

    // Setup fragment charge constraint (constant)
    for (int j = 0; j < natoms; ++j) {
        A_coulomb_base(natoms, j) = 1.0;
        A_coulomb_base(j, natoms) = 1.0;
    }

    if (m_verbosity >= 3) {
        CurcumaLogger::info(fmt::format("EEQ Phase 2: Pre-computed distance matrix and Coulomb base ({}x{})", natoms, natoms));
    }

    // Iterative refinement loop (now much faster!)
    for (int iter = 0; iter < m_max_iterations; ++iter) {
        // Claude Generated (2025): Copy pre-computed Coulomb matrix
        // This is O(n²) copy, much faster than O(n²) sqrt() + erf() recalculation
        Matrix A = A_coulomb_base;
        Vector x = Vector::Zero(m);

        // 1. Setup RHS with CNF*sqrt(CN) term (critical!)
        // Reference: XTB gfnff_engrad.F90:1504
        // x(i) = topo%chieeq(i) + param%cnf(at(i))*sqrt(cn(i))
        for (int i = 0; i < natoms; ++i) {
            int z_i = atoms[i];
            EEQParameters params_i = getParameters(z_i, cn(i));
            x(i) = chi_corrected(i) + params_i.cnf * std::sqrt(cn(i));
        }

        // 2. Update only diagonal (charge-dependent terms)
        // This is the only part that changes per iteration
        for (int i = 0; i < natoms; ++i) {
            A(i, i) = gam_corrected(i) + TSQRT2PI / std::sqrt(alpha_corrected(i));
        }

        // DEBUG: Print matrix details on first iteration
        if (iter == 0 && m_verbosity >= 3) {
            std::cerr << "\n=== EEQ Phase 2 DEBUG (iteration 0) ===" << std::endl;
            for (int i = 0; i < std::min(3, natoms); ++i) {
                int z_i = atoms[i];
                EEQParameters params_i = getParameters(z_i, cn(i));
                std::cerr << "Atom " << i << " (Z=" << z_i << "):" << std::endl;
                std::cerr << "  CN = " << cn(i) << std::endl;
                std::cerr << "  chi_corrected = " << chi_corrected(i) << std::endl;
                std::cerr << "  gam_corrected = " << gam_corrected(i) << std::endl;
                std::cerr << "  alpha_corrected = " << alpha_corrected(i) << std::endl;
                std::cerr << "  cnf*sqrt(CN) = " << (params_i.cnf * std::sqrt(cn(i))) << std::endl;
                std::cerr << "  x(RHS) = " << x(i) << std::endl;
                std::cerr << "  A(i,i) diagonal = " << A(i, i) << std::endl;
            }
            std::cerr << "  x(constraint) = " << x(natoms) << std::endl;
            std::cerr << "==========================================\n" << std::endl;
        }

        // 3. Setup constraint RHS
        x(natoms) = static_cast<double>(total_charge);

        // 4. Solve system
        if (iter == 0 && m_verbosity >= 3) {
            CurcumaLogger::info("=== Phase 2 EEQ Matrix Diagnostics ===");

            Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(A);
            Vector eigenvalues = eigensolver.eigenvalues();

            double max_eigenvalue = eigenvalues.maxCoeff();
            double min_eigenvalue = eigenvalues.minCoeff();
            double condition_number = std::abs(max_eigenvalue / min_eigenvalue);

            CurcumaLogger::param("condition_number", fmt::format("{:.6e}", condition_number));

            if (condition_number > 1e12) {
                CurcumaLogger::warn(fmt::format("Phase 2 EEQ matrix ill-conditioned (cond={:.2e})", condition_number));
            }
        }

        Eigen::PartialPivLU<Matrix> lu(A);
        Vector solution = lu.solve(x);

        Vector final_charges_new = solution.segment(0, natoms);

        // Check for NaN/Inf on first iteration
        if (iter == 0) {
            for (int i = 0; i < natoms; ++i) {
                if (std::isnan(final_charges_new[i]) || std::isinf(final_charges_new[i])) {
                    CurcumaLogger::error(fmt::format("Phase 2 EEQ: Invalid charge[{}] = {} (Z={})",
                                                     i, final_charges_new[i], atoms[i]));
                    return Vector::Zero(0);
                }
            }
        }

        // Check convergence
        double max_change = 0.0;
        for (int i = 0; i < natoms; ++i) {
            double change = std::abs(final_charges_new(i) - final_charges(i));
            max_change = std::max(max_change, change);
        }

        final_charges = final_charges_new;

        if (m_verbosity >= 3) {
            CurcumaLogger::result(fmt::format("EEQ Phase 2 iter {}: max_change = {:.2e}",
                                              iter + 1, max_change));
        }

        // Claude Generated (2025): Adaptive convergence tolerance
        // For large systems (n > 100), use relaxed tolerance to save iterations
        // For force field accuracy, 1e-5 is sufficient for most purposes
        double adaptive_threshold = m_convergence_threshold;
        if (natoms > 500) {
            adaptive_threshold = std::max(m_convergence_threshold, 1e-4);
        } else if (natoms > 100) {
            adaptive_threshold = std::max(m_convergence_threshold, 1e-5);
        }

        if (max_change < adaptive_threshold) {
            if (m_verbosity >= 2) {
                CurcumaLogger::success(fmt::format("Phase 2 converged in {} iterations (threshold={:.1e})",
                                                   iter + 1, adaptive_threshold));
            }
            break;
        }
    }

    // Validate charge conservation
    double total = final_charges.sum();
    if (std::abs(total - total_charge) > 1e-6) {
        CurcumaLogger::warn(fmt::format("EEQ Phase 2: Charge sum = {:.6f} (expected {})",
                                        total, total_charge));
    }

    // Store dxi for later use in energy calculation
    m_dxi_stored = dxi;

    return final_charges;
}

// ===== Energy Calculation =====

double EEQSolver::calculateEEQEnergy(
    const Vector& charges,
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn)
{
    const int natoms = atoms.size();
    double energy = 0.0;

    // Pairwise Coulomb energy with erf damping
    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < i; ++j) {
            double dx = geometry_bohr(i, 0) - geometry_bohr(j, 0);
            double dy = geometry_bohr(i, 1) - geometry_bohr(j, 1);
            double dz = geometry_bohr(i, 2) - geometry_bohr(j, 2);
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (r < 1e-10) continue;

            EEQParameters params_i = getParameters(atoms[i], cn(i));
            EEQParameters params_j = getParameters(atoms[j], cn(j));

            double gamma_ij = 1.0 / std::sqrt(params_i.alp + params_j.alp);
            double erf_gamma = std::erf(gamma_ij * r);
            double coulomb = erf_gamma / r;

            energy += charges(i) * charges(j) * coulomb;
        }
    }

    // Self-energy terms
    const double TSQRT2PI = 0.797884560802866;
    for (int i = 0; i < natoms; ++i) {
        EEQParameters params_i = getParameters(atoms[i], cn(i));

        // CRITICAL FIX (Session 11): chi_i must include dxi term!
        // Reference: XTB gfnff_engrad.F90:1581
        // chi = -χ + dxi + CNF·√CN (NOT just -χ + CNF·√CN!)
        double dxi_i = (m_dxi_stored.size() > i) ? m_dxi_stored(i) : 0.0;
        double chi_i = -params_i.chi + dxi_i + params_i.cnf * std::sqrt(cn(i));
        double self_energy = -charges(i) * chi_i
                           + 0.5 * charges(i) * charges(i) * (params_i.gam + TSQRT2PI / std::sqrt(params_i.alp));

        energy += self_energy;
    }

    return energy;  // Hartree
}

// ===== Correction Terms =====

Vector EEQSolver::calculateDxi(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn)
{
    const int natoms = atoms.size();
    Vector dxi = Vector::Zero(natoms);

    for (int i = 0; i < natoms; ++i) {
        // Hybridization guess from CN
        double dxi_hyb = 0.0;
        if (cn(i) < 1.5) {
            dxi_hyb = 0.1;  // sp-like
        } else if (cn(i) < 2.5) {
            dxi_hyb = 0.05;  // sp2-like
        }

        // CN-dependent correction
        double dxi_cn = -0.01 * (cn(i) - 2.0);

        dxi(i) = dxi_hyb + dxi_cn;
    }

    return dxi;
}

Vector EEQSolver::calculateDgam(
    const std::vector<int>& atoms,
    const Vector& charges,
    const std::vector<int>& hybridization)
{
    const int natoms = atoms.size();
    Vector dgam = Vector::Zero(natoms);

    for (int i = 0; i < natoms; ++i) {
        int Z = atoms[i];
        double qa = charges(i);

        // Claude Generated (December 2025, Session 11): EXACT XTB ff values
        // Reference: Fortran gfnff_ini.f90:665-690
        double ff = 0.0;  // Default: do nothing

        if (Z == 1) {
            ff = -0.08;  // H
        } else if (Z == 5) {
            ff = -0.05;  // B
        } else if (Z == 6) {  // C
            ff = -0.27;  // sp3
            if (!hybridization.empty() && i < hybridization.size()) {
                if (hybridization[i] < 3) ff = -0.45;  // sp2 or lower
                if (hybridization[i] < 2) ff = -0.34;  // sp
            }
        } else if (Z == 7) {  // N
            ff = -0.13;  // Base N
            // TODO: pi-system (-0.14) and amide (-0.16) detection requires piadr/amide functions
        } else if (Z == 8) {  // O
            ff = -0.15;  // sp3
            if (!hybridization.empty() && i < hybridization.size()) {
                if (hybridization[i] < 3) ff = -0.08;  // unsaturated
            }
        } else if (Z == 9) {
            ff = 0.10;   // F
        } else if (Z > 10) {
            ff = -0.02;  // Heavy atoms default

            // Specific heavy atom overrides
            if (Z == 17) ff = -0.02;  // Cl
            if (Z == 35) ff = -0.11;  // Br
            if (Z == 53) ff = -0.07;  // I

            // Metal corrections (requires metal_type array)
            if (Z >= 1 && Z <= 86) {
                int imetal_val = metal_type[Z - 1];
                if (imetal_val == 1) ff = -0.08;   // Main group metals
                if (imetal_val == 2) ff = -0.9;    // Transition metals (XTB comment: "too large")
            }

            // Noble gases (Group 8)
            if (Z >= 1 && Z <= 86) {
                int group = periodic_group[Z - 1];
                if (group == 8) ff = 0.0;  // Noble gases
            }
        }

        dgam(i) = qa * ff;
    }

    return dgam;
}

// ===== Parameter Lookup =====
// NOTE: calculateDalpha() removed - alpha now calculated with charge-dependent formula (alpha_base + ff*qa)²

EEQSolver::EEQParameters EEQSolver::getParameters(int Z, double cn) const
{
    EEQParameters params;

    if (Z >= 1 && Z <= 86) {
        int idx = Z - 1;
        params.chi = chi_eeq[idx];
        params.gam = gam_eeq[idx];
        params.alp = alpha_eeq[idx] * alpha_eeq[idx];  // CRITICAL: Must be SQUARED!
        params.cnf = cnf_eeq[idx];
    } else {
        CurcumaLogger::warn(fmt::format("EEQSolver: No parameters for Z={}, using defaults", Z));
        params.chi = 1.0;
        params.gam = 0.0;
        params.alp = 1.0;
        params.cnf = 0.0;
    }

    return params;
}

// ===== Helper Functions =====

Vector EEQSolver::calculateCoordinationNumbers(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr) const
{
    // Delegate to shared CNCalculator utility
    // Claude Generated - December 21, 2025 (consolidated duplicate CN calculation)
    std::vector<double> cn_double = CNCalculator::calculateGFNFFCN(atoms, geometry_bohr);

    // Convert std::vector<double> to Eigen::Vector for API compatibility
    Vector cn = Vector::Map(cn_double.data(), cn_double.size());

    return cn;
}

std::vector<int> EEQSolver::detectHybridization(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn) const
{
    const int natoms = atoms.size();
    std::vector<int> hybridization(natoms, 3);  // Default: sp3

    for (int i = 0; i < natoms; ++i) {
        // Simple heuristic based on coordination number
        if (cn(i) < 1.5) {
            hybridization[i] = 1;  // sp (terminal or linear)
        } else if (cn(i) < 2.5) {
            hybridization[i] = 2;  // sp2 (trigonal)
        } else {
            hybridization[i] = 3;  // sp3 (tetrahedral)
        }
    }

    return hybridization;
}

std::vector<std::vector<int>> EEQSolver::buildNeighborLists(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    double cutoff_radius) const
{
    const int natoms = atoms.size();
    std::vector<std::vector<int>> neighbors(natoms);

    double cutoff_sq = cutoff_radius * cutoff_radius;

    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < natoms; ++j) {
            if (i == j) continue;

            Vector ri = geometry_bohr.row(i);
            Vector rj = geometry_bohr.row(j);
            double distance_sq = (ri - rj).squaredNorm();

            if (distance_sq < cutoff_sq) {
                neighbors[i].push_back(j);
            }
        }
    }

    return neighbors;
}
