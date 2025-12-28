/*
 * < EEQ (Electronegativity Equalization) Charge Solver >
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
 * This file contains a standalone EEQ (Electronegativity Equalization) solver
 * extracted from GFN-FF implementation for reusability across force field methods.
 *
 * References:
 * - S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665-15673
 *   (GFN-FF with EEQ charge model)
 * - R. T. Sanderson, Chemical Bonds and Bond Energy, Academic Press, 1976
 *   (Original electronegativity equalization concept)
 */

#pragma once

#include "src/core/global.h"
#include "src/core/config_manager.h"
#include "src/core/parameter_macros.h"

#include <Eigen/Dense>
#include <vector>
#include <memory>

/**
 * @brief EEQ (Electronegativity Equalization) charge solver
 *
 * Standalone utility for calculating atomic partial charges using the EEQ method.
 * Implements the two-phase algorithm from GFN-FF:
 * - Phase 1: Topology charges (quick initial charges from augmented linear system)
 * - Phase 2: Final charges (iterative refinement with environmental corrections)
 *
 * Usage:
 *   ConfigManager config("eeq_solver", json{});
 *   EEQSolver solver(config);
 *   Vector charges = solver.calculateCharges(atoms, geometry_bohr, total_charge);
 *
 * Claude Generated - December 2025
 */
class EEQSolver {
public:
    /**
     * @brief Construct EEQ solver from configuration
     * @param config ConfigManager with eeq_solver parameters
     */
    explicit EEQSolver(const ConfigManager& config);

    /**
     * @brief Default destructor
     */
    ~EEQSolver() = default;

    // ===== Main API =====

    /**
     * @brief Calculate atomic partial charges using two-phase EEQ algorithm
     *
     * @param atoms Atomic numbers (Z) for each atom
     * @param geometry_bohr Cartesian coordinates in Bohr (natoms x 3)
     * @param total_charge Total molecular charge (default 0)
     * @param cn_hint Optional pre-calculated coordination numbers (nullptr to calculate internally)
     * @param hyb_hint Optional pre-calculated hybridization states (nullptr to use defaults)
     * @return Vector of atomic partial charges in elementary charge units (e)
     *
     * Automatically uses two-phase algorithm for better accuracy:
     * - Phase 1: Solve augmented EEQ system for topology charges
     * - Phase 2: Iterative refinement with dxi, dgam, dalpha corrections
     */
    Vector calculateCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge = 0,
        const Vector* cn_hint = nullptr,
        const std::vector<int>* hyb_hint = nullptr
    );

    // ===== Advanced API (two-phase workflow) =====

    /**
     * @brief Phase 1: Calculate initial topology charges
     *
     * Solves augmented EEQ linear system:
     * - A_{ii} = gam_i + sqrt(2/π)/sqrt(alpha_i)
     * - A_{ij} = erf(γ_ij * r_ij) / r_ij
     * - Constraint row/column for charge conservation
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param total_charge Total molecular charge
     * @param cn Coordination numbers (for chi corrections)
     * @return Vector of topology charges (qa)
     */
    Vector calculateTopologyCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge,
        const Vector& cn
    );

    /**
     * @brief Phase 2: Refine charges with environmental corrections
     *
     * Iteratively solves EEQ with corrected parameters:
     * - chi' = chi + dxi(environment)
     * - gam' = gam + dgam(charge-dependent)
     * - alpha' = alpha + dalpha(environment)
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param total_charge Total molecular charge
     * @param topology_charges Initial charges from Phase 1
     * @param cn Coordination numbers
     * @param hybridization Hybridization states (1=sp, 2=sp2, 3=sp3)
     * @return Vector of final refined charges
     */
    Vector calculateFinalCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge,
        const Vector& topology_charges,
        const Vector& cn,
        const std::vector<int>& hybridization
    );

    /**
     * @brief Calculate EEQ electrostatic energy from charges
     *
     * @param charges Atomic partial charges
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers (for cnf corrections)
     * @return EEQ energy in Hartree
     */
    double calculateEEQEnergy(
        const Vector& charges,
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn
    );

private:
    // ===== EEQ Parameter Structure =====

    /**
     * @brief Element-specific EEQ parameters
     *
     * From Spicher & Grimme Angew. Chem. Int. Ed. 2020, 59, 15665 (angewChem2020 set)
     */
    struct EEQParameters {
        double chi;    ///< Electronegativity (Hartree)
        double gam;    ///< Chemical hardness (Hartree)
        double alp;    ///< Damping parameter (Bohr^-1, SQUARED!)
        double cnf;    ///< CN correction factor (dimensionless)
    };

    /**
     * @brief Get EEQ parameters for element Z
     *
     * @param Z Atomic number (1-86)
     * @param cn Coordination number (for chi correction)
     * @return EEQParameters structure
     *
     * CRITICAL: alpha value is SQUARED (alpha_eeq[Z-1]²) before storage!
     */
    EEQParameters getParameters(int Z, double cn = 0.0) const;

    // ===== Correction Terms =====

    /**
     * @brief Calculate electronegativity corrections (dxi)
     *
     * Environment-dependent chi corrections based on coordination number
     * and neighbor electronegativity.
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @return Vector of dxi corrections (Hartree)
     */
    Vector calculateDxi(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn
    );

    /**
     * @brief Calculate hardness corrections (dgam)
     *
     * Charge-dependent gamma corrections with element-specific factors.
     * Formula: dgam(i) = qa(i) * ff
     * where ff depends on element and hybridization.
     *
     * @param atoms Atomic numbers
     * @param charges Current charges (from Phase 1)
     * @param hybridization Hybridization states
     * @return Vector of dgam corrections (Hartree)
     */
    Vector calculateDgam(
        const std::vector<int>& atoms,
        const Vector& charges,
        const std::vector<int>& hybridization
    );

    // NOTE: calculateDalpha() removed - alpha now calculated inline with charge-dependent formula

    // ===== Helper Functions =====

    /**
     * @brief Calculate coordination numbers from geometry
     *
     * Uses covalent radii with exponential distance weighting.
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @return Vector of coordination numbers
     */
    Vector calculateCoordinationNumbers(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr
    ) const;

    /**
     * @brief Detect hybridization states from geometry
     *
     * Simple geometric analysis: 1=sp, 2=sp2, 3=sp3
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @return Vector of hybridization states
     */
    std::vector<int> detectHybridization(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn
    ) const;

    /**
     * @brief Build neighbor lists for each atom
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cutoff_radius Cutoff radius in Bohr (default 10.0)
     * @return Vector of neighbor index lists
     */
    std::vector<std::vector<int>> buildNeighborLists(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        double cutoff_radius = 10.0
    ) const;

    // ===== Configuration =====

    ConfigManager m_config;           ///< Configuration manager
    int m_max_iterations;             ///< Maximum iterations for Phase 2 refinement
    double m_convergence_threshold;   ///< Convergence threshold for charge changes (e)
    int m_verbosity;                  ///< Verbosity level (0-3)
    bool m_calculate_cn;              ///< Auto-calculate CN if not provided

    // ===== Cached Data for Energy Calculation =====

    mutable Vector m_dxi_stored;      ///< Stored dxi corrections from last calculateCharges() call
};

// ===== Parameter Definitions =====

BEGIN_PARAMETER_DEFINITION(eeq_solver)
    PARAM(max_iterations, Int, 50,
          "Maximum iterations for EEQ charge refinement (Phase 2)", "Algorithm", {})
    PARAM(convergence_threshold, Double, 1e-6,
          "Convergence threshold for charge changes (e)", "Algorithm", {})
    PARAM(verbosity, Int, 0,
          "Verbosity level (0=silent, 1=basic, 2=detailed, 3=debug)", "Output", {})
    PARAM(calculate_cn, Bool, true,
          "Auto-calculate coordination numbers if not provided", "Algorithm", {})
END_PARAMETER_DEFINITION
