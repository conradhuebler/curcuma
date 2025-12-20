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

        chi(i) = -params_i.chi + dxi_total;
        gam(i) = params_i.gam;
        alpha(i) = params_i.alp;  // Already squared
    }

    // Build AUGMENTED EEQ matrix A (m x m)
    Matrix A = Matrix::Zero(m, m);
    Vector x = Vector::Zero(m);

    // 1. Setup RHS and diagonal
    for (int i = 0; i < natoms; ++i) {
        x(i) = chi(i);
        A(i, i) = gam(i) + TSQRT2PI / std::sqrt(alpha(i));
    }

    // 2. Setup off-diagonal Coulomb matrix
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

    // Calculate correction terms
    Vector dxi = calculateDxi(atoms, geometry_bohr, cn);
    Vector dgam = calculateDgam(atoms, topology_charges, hybridization);
    Vector dalpha = calculateDalpha(atoms, geometry_bohr, cn);

    Vector final_charges = topology_charges;

    // Augmented system size
    int m = natoms + 1;

    // Iterative refinement loop
    for (int iter = 0; iter < m_max_iterations; ++iter) {
        // Prepare corrected parameters
        Vector chi_corrected(natoms);
        Vector gam_corrected(natoms);
        Vector alpha_corrected(natoms);

        for (int i = 0; i < natoms; ++i) {
            int z_i = atoms[i];
            EEQParameters params_i = getParameters(z_i, cn(i));

            chi_corrected(i) = -params_i.chi + dxi(i);
            gam_corrected(i) = params_i.gam + dgam(i);
            alpha_corrected(i) = params_i.alp + dalpha(i);
        }

        // Build augmented EEQ matrix
        Matrix A = Matrix::Zero(m, m);
        Vector x = Vector::Zero(m);

        // 1. Setup RHS and diagonal
        for (int i = 0; i < natoms; ++i) {
            x(i) = chi_corrected(i);
            A(i, i) = gam_corrected(i) + TSQRT2PI / std::sqrt(alpha_corrected(i));
        }

        // 2. Setup off-diagonal Coulomb matrix
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

                double gamma_ij = 1.0 / std::sqrt(alpha_corrected(i) + alpha_corrected(j));
                double erf_gamma = std::erf(gamma_ij * r);
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

        if (max_change < m_convergence_threshold) {
            if (m_verbosity >= 2) {
                CurcumaLogger::success(fmt::format("Phase 2 converged in {} iterations", iter + 1));
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

        double chi_i = -params_i.chi + params_i.cnf * std::sqrt(cn(i));
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
        double ff = -0.04;  // Base default

        // Element-specific factors (from Fortran gfnff_ini.f90:677-688)
        if (!hybridization.empty() && i < hybridization.size()) {
            if (hybridization[i] < 3) {
                ff = -0.08;  // Unsaturated
            }
        }

        if (Z == 9) ff = 0.10;   // Fluorine
        if (Z > 10) ff = -0.02;  // Heavy atoms
        if (Z == 17) ff = -0.02; // Chlorine
        if (Z == 35) ff = -0.11; // Bromine
        if (Z == 53) ff = -0.07; // Iodine

        dgam(i) = qa * ff;
    }

    return dgam;
}

Vector EEQSolver::calculateDalpha(
    const std::vector<int>& atoms,
    const Matrix& geometry_bohr,
    const Vector& cn)
{
    const int natoms = atoms.size();
    Vector dalpha = Vector::Zero(natoms);

    for (int i = 0; i < natoms; ++i) {
        // CN-dependent correction
        double dalpha_cn = -0.02 * (cn(i) - 2.0);

        // Simple hybridization guess from CN
        double dalpha_hyb = 0.0;
        if (cn(i) < 1.5) {
            dalpha_hyb = 0.05;  // sp-like
        } else if (cn(i) < 2.5) {
            dalpha_hyb = 0.02;  // sp2-like
        }

        dalpha(i) = dalpha_cn + dalpha_hyb;
    }

    return dalpha;
}

// ===== Parameter Lookup =====

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
    const int natoms = atoms.size();
    const double kn = -7.5;       // Error function steepness
    const double cnmax = 4.4;     // Maximum CN cutoff
    const double threshold = 1600.0;  // Distance threshold (Bohr²)

    // D3 covalent radii (Angstrom) - first 18 elements
    static const std::vector<double> rcov_d3_angstrom = {
        0.32, 0.46,  // H, He
        1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,  // Li-Ne
        1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96   // Na-Ar
    };

    Vector cn = Vector::Zero(natoms);

    // Calculate raw CN using error function
    for (int i = 0; i < natoms; ++i) {
        double cn_i = 0.0;
        for (int j = 0; j < natoms; ++j) {
            if (i == j) continue;

            Vector ri = geometry_bohr.row(i);
            Vector rj = geometry_bohr.row(j);
            double distance_sq = (ri - rj).squaredNorm();

            if (distance_sq > threshold) continue;

            double distance = std::sqrt(distance_sq);

            // Get D3 covalent radii (fallback to Elements::CovalentRadius)
            double rcov_i_angstrom = (atoms[i] >= 1 && atoms[i] <= 18)
                                   ? rcov_d3_angstrom[atoms[i] - 1]
                                   : ((atoms[i] >= 1 && atoms[i] < static_cast<int>(Elements::CovalentRadius.size()))
                                      ? Elements::CovalentRadius[atoms[i]] : 0.7);
            double rcov_j_angstrom = (atoms[j] >= 1 && atoms[j] <= 18)
                                   ? rcov_d3_angstrom[atoms[j] - 1]
                                   : ((atoms[j] >= 1 && atoms[j] < static_cast<int>(Elements::CovalentRadius.size()))
                                      ? Elements::CovalentRadius[atoms[j]] : 0.7);

            // Convert to Bohr (1 Angstrom = 1.8897259886 Bohr)
            double rcov_i_bohr = rcov_i_angstrom * 1.8897259886;
            double rcov_j_bohr = rcov_j_angstrom * 1.8897259886;
            double rcov_ij = rcov_i_bohr + rcov_j_bohr;

            // Error function coordination number
            double dr = (distance - rcov_ij) / rcov_ij;
            double erfCN = 0.5 * (1.0 + std::erf(kn * dr));

            cn_i += erfCN;
        }

        // Apply logarithmic transformation
        cn(i) = std::log(1.0 + std::exp(cnmax)) - std::log(1.0 + std::exp(cnmax - cn_i));
    }

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
