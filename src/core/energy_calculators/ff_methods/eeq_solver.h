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
#include <optional>

/**
 * @brief Distance mode for EEQ matrix construction
 *
 * Controls which type of interatomic distances are used in the EEQ Coulomb matrix:
 * - Topological: Floyd-Warshall bond path distances (for Phase 1 topology charges)
 * - Geometric: Euclidean xyz distances (for Phase 2 energy charges)
 *
 * CRITICAL FIX (January 17, 2026): The Fortran reference uses DIFFERENT distance modes
 * for Phase 1 vs Phase 2:
 * - Phase 1 (goedeckera): Topological distances (`pair(ij)` from Floyd-Warshall)
 * - Phase 2 (goed_gfnff): Geometric distances (`r(ij)` from xyz coordinates)
 *
 * Claude Generated - January 17, 2026
 */
enum class EEQDistanceMode {
    Topological,  ///< Phase 1: Floyd-Warshall bond path distances
    Geometric     ///< Phase 2: Euclidean xyz distances
};

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
     * @brief Topology information for EEQ Phase 1
     *
     * Used to compute topological distances via Floyd-Warshall algorithm
     * instead of geometric (line-of-sight) distances.
     *
     * Topological distances represent shortest paths through the bond graph,
     * where each bond has length = sum of covalent radii. This matches the
     * XTB implementation in gfnff_ini.f90:431-461.
     *
     * Claude Generated December 2025
     */
    struct TopologyInput {
        std::vector<std::vector<int>> neighbor_lists;  // neighbor_lists[i] = atoms bonded to i
        std::vector<double> covalent_radii;            // Covalent radii in ANGSTROM (NOT Bohr!)
        int nfrag = 1;                                 // Number of molecular fragments
        std::vector<int> fraglist;                     // fraglist[i] = fragment ID for atom i (1-indexed)
        std::vector<double> qfrag;                     // qfrag[f] = target charge for fragment f
    };

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
     * @param topology Optional topology information for Floyd-Warshall topological distances
     *                 If provided, uses topological distances (more accurate, matches XTB)
     *                 If omitted, falls back to geometric distances (deprecated)
     * @return Vector of atomic partial charges in elementary charge units (e)
     *
     * Automatically uses two-phase algorithm for better accuracy:
     * - Phase 1: Solve augmented EEQ system for topology charges
     * - Phase 2: Iterative refinement with dxi, dgam, dalpha corrections
     *
     * @note CRITICAL FIX (Dec 28, 2025): Added topology parameter to enable Floyd-Warshall
     *       topological distances in Phase 1, matching XTB reference implementation.
     *       Without topology, geometric distances produce charges 4-5× too large!
     */
    Vector calculateCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge = 0,
        const Vector* cn_hint = nullptr,
        const std::vector<int>* hyb_hint = nullptr,
        const std::optional<TopologyInput>& topology = std::nullopt
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
     * @param topology Optional topology information (neighbor lists, covalent radii)
     *                 If provided, uses Floyd-Warshall topological distances
     *                 If omitted, falls back to geometric distances (deprecated)
     * @return Vector of topology charges (qa)
     *
     * @note Topological distances computed via Floyd-Warshall algorithm match XTB reference
     *       and produce charges ~4-5× smaller than geometric distances.
     *       Reference: XTB gfnff_ini.f90:431-461
     */
    Vector calculateTopologyCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt,
        bool use_corrections = false  // CRITICAL FIX (Jan 4, 2026): default false to match gfnff_final.cpp
    );

    /**
     * @brief Phase 2: Calculate final EEQ charges with environmental corrections
     *
     * Single linear solve (matches XTB gfnff_ini.f90:699-706) with corrected parameters:
     * - chi' = chi + dxi(environment)
     * - gam' = gam + dgam(charge-dependent, uses topology_charges)
     * - alpha' = alpha_base² (BASE VALUE, NOT charge-corrected!)
     *
     * CRITICAL (Jan 7, 2026): Alpha uses BASE values in Phase 2, unlike Phase 1c which
     * applies alpha' = (alpha_base + ff*qa)². This is the key difference between topology
     * and energy charge calculations. The dxi/dgam corrections are applied, but alpha
     * remains constant to maintain numerical stability.
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param total_charge Total molecular charge
     * @param topology_charges Phase-1 charges (qa) used for dgam calculation (NOT for alpha!)
     * @param cn Coordination numbers
     * @param hybridization Hybridization states (1=sp, 2=sp2, 3=sp3)
     * @param topology Optional topology information for integer neighbor counts in CNF calculation
     * @param alpeeq Optional charge-dependent alpha values (squared) from Phase 1
     *               If provided, these pre-computed alpeeq values are used instead of base alpha
     *               Formula: alpeeq(i) = (alpha_base + ff*qa(i))²
     *               Reference: Fortran gfnff_ini.f90:718-725
     * @return Vector of final EEQ charges
     */
    Vector calculateFinalCharges(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        int total_charge,
        const Vector& topology_charges,
        const Vector& cn,
        const std::vector<int>& hybridization,
        const std::optional<TopologyInput>& topology = std::nullopt,
        bool use_corrections = false,  // CRITICAL FIX (Jan 4, 2026): default false to match gfnff_final.cpp
        const std::optional<Vector>& alpeeq = std::nullopt  // Claude Generated (January 2026): Charge-dependent alpha
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
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
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
        const std::vector<int>& hybridization,
        const std::vector<bool>& is_pi_atom = {},
        const std::vector<bool>& is_amide = {}
    );

    // NOTE: calculateDalpha() removed - alpha now calculated inline with charge-dependent formula

    // ===== NEW Single-Solve Helper Functions (Jan 4, 2026) =====

    /**
     * @brief Build EEQ Coulomb matrix with ALL charge-dependent corrections
     *
     * Constructs the augmented EEQ matrix with:
     * - Diagonal: gam_corrected + sqrt(2/π)/sqrt(alpha_corrected)
     * - Off-diagonal: Coulomb terms with corrected alpha
     * - Constraint row/column for charge conservation
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @param current_charges Current charge estimate (used for charge-dependent dgam/alpha)
     * @param dxi Electronegativity corrections
     * @param dgam Hardness corrections (gam - qa*ff)
     * @param hybridization Hybridization states
     * @param topology Optional topology for topological distances
     * @param distance_mode Distance calculation mode (Topological for Phase 1, Geometric for Phase 2)
     *                      CRITICAL (Jan 17, 2026): Fortran uses topological for Phase 1 (goedeckera),
     *                      geometric for Phase 2 (goed_gfnff). Default=Topological for backward compat.
     * @return Augmented EEQ matrix (natoms+1)×(natoms+1)
     */
    Matrix buildCorrectedEEQMatrix(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn,
        const Vector& current_charges,
        const Vector& dxi,
        const Vector& dgam,
        const std::vector<int>& hybridization,
        const std::optional<TopologyInput>& topology,
        EEQDistanceMode distance_mode = EEQDistanceMode::Topological
    );

    /**
     * @brief Build EEQ matrix with intelligent caching for performance optimization
     *
     * Enhanced version of buildCorrectedEEQMatrix that uses intelligent caching
     * to avoid expensive matrix reconstruction when geometry changes are insignificant.
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @param current_charges Current charge estimate (used for charge-dependent dgam/alpha)
     * @param dxi Electronegativity corrections
     * @param dgam Hardness corrections (gam - qa*ff)
     * @param hybridization Hybridization states
     * @param topology Optional topology for topological distances
     * @param distance_mode Distance calculation mode (Topological for Phase 1, Geometric for Phase 2)
     * @return Augmented EEQ matrix (natoms+1)×(natoms+1) with caching
     *
     * Claude Generated - Performance Optimization Implementation
     */
    Matrix buildSmartEEQMatrix(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn,
        const Vector& current_charges,
        const Vector& dxi,
        const Vector& dgam,
        const std::vector<int>& hybridization,
        const std::optional<TopologyInput>& topology,
        EEQDistanceMode distance_mode = EEQDistanceMode::Topological
    );

    /**
     * @brief Solve augmented EEQ linear system with corrected parameters
     *
     * Sets up RHS: x(i) = -chi + dxi [+ CNF*sqrt(nb) if use_cnf_term=true]
     * Solves the augmented system: A*[q; lambda] = [x; total_charge]
     *
     * CRITICAL (Jan 4, 2026): CNF term ONLY in Phase 1 (topology charges)!
     * Phase 2 (final charges with dgam) does NOT use CNF term!
     * Reference: XTB gfnff_ini.f90 lines 563-570 vs 696-707
     *
     * @param A Augmented EEQ matrix from buildCorrectedEEQMatrix()
     * @param atoms Atomic numbers
     * @param cn Coordination numbers
     * @param dxi Electronegativity corrections
     * @param total_charge Total molecular charge
     * @param topology Optional topology for integer neighbor counts in CNF
     * @param use_cnf_term If true, adds CNF*sqrt(nb) to RHS (Phase 1 only!)
     * @return Vector of atomic charges (natoms elements)
     */
    Vector solveEEQ(
        const Matrix& A,
        const std::vector<int>& atoms,
        const Vector& cn,
        const Vector& dxi,
        int total_charge,
        const std::optional<TopologyInput>& topology,
        bool use_cnf_term = true,  // Default true for backward compatibility
        bool use_integer_nb = true // Phase 1 uses integer nb, Phase 2 uses fractional cn
    );

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
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    ) const;

    /**
     * @brief Detect atoms in pi-systems
     *
     * Simplified detection for EEQ corrections: (sp or sp2) AND (C,N,O,F,S)
     * plus picon (sp3 N,O,F bonded to sp/sp2 atom).
     *
     * @param atoms Atomic numbers
     * @param hybridization Hybridization states
     * @param topology Topology information
     * @return Vector of pi-system membership flags
     */
    std::vector<bool> detectPiSystem(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const std::optional<TopologyInput>& topology
    ) const;

    /**
     * @brief Detect amide nitrogens
     *
     * @param atoms Atomic numbers
     * @param hybridization Hybridization states
     * @param is_pi_atom Pi-system membership flags
     * @param topology Topology information
     * @param cn Coordination numbers
     * @return Vector of amide nitrogen flags
     */
    std::vector<bool> detectAmideNitrogens(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const std::vector<bool>& is_pi_atom,
        const std::optional<TopologyInput>& topology,
        const Vector& cn
    ) const;

    /**
     * @brief Detect amide hydrogens
     *
     * @param atoms Atomic numbers
     * @param hybridization Hybridization states
     * @param is_amide Vector of amide nitrogen flags
     * @param topology Topology information
     * @return Vector of amide hydrogen flags
     *
     * Claude Generated January 2026 - Ported from Fortran gfnff_ini2.f90:1566
     */
    std::vector<bool> detectAmideHydrogens(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const std::vector<bool>& is_amide,
        const std::optional<TopologyInput>& topology
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

    /**
     * @brief Compute topological distances via Floyd-Warshall algorithm
     *
     * Computes shortest path distances through the bond graph, where each bond
     * has length = sum of covalent radii. Topological distances are always
     * greater than or equal to geometric distances.
     *
     * Reference: XTB gfnff_ini.f90:431-461
     *
     * @param atoms Atomic numbers
     * @param topology Topology information (neighbor lists, covalent radii)
     * @return Matrix of topological distances in Bohr
     *
     * Claude Generated December 2025
     */
    Matrix computeTopologicalDistances(
        const std::vector<int>& atoms,
        const TopologyInput& topology
    ) const;

    // ===== Configuration =====

    ConfigManager m_config;           ///< Configuration manager
    int m_max_iterations;             ///< Maximum iterations for Phase 2 refinement
    double m_convergence_threshold;   ///< Convergence threshold for charge changes (e)
    int m_verbosity;                  ///< Verbosity level (0-3)
    bool m_calculate_cn;              ///< Auto-calculate CN if not provided

    // ===== Cached Data for Energy Calculation =====

    // EEQ Solver intelligent caching for performance optimization
    class EEQSolverCache {
    private:
        Matrix m_last_geometry;
        Matrix m_last_A_matrix;
        Vector m_last_charges;
        bool m_cache_valid = false;
        double m_change_threshold = 1e-6;

    public:
        bool isGeometryChanged(const Matrix& current_geometry) const {
            if (!m_cache_valid) return true;
            if (m_last_geometry.rows() != current_geometry.rows() ||
                m_last_geometry.cols() != current_geometry.cols()) {
                return true;
            }
            return (m_last_geometry - current_geometry).array().abs().maxCoeff() > m_change_threshold;
        }

        void cacheResults(const Matrix& geometry, const Matrix& A, const Vector& charges) {
            m_last_geometry = geometry;
            m_last_A_matrix = A;
            m_last_charges = charges;
            m_cache_valid = true;
        }

        Matrix getCachedAMatrix() const { return m_last_A_matrix; }
        Vector getCachedCharges() const { return m_last_charges; }
        bool isValid() const { return m_cache_valid; }

        void reset() {
            m_cache_valid = false;
            m_last_geometry = Matrix();
            m_last_A_matrix = Matrix();
            m_last_charges = Vector();
        }
    };

    mutable Vector m_dxi_stored;      ///< Stored dxi corrections from last calculateCharges() call
    mutable Matrix m_cached_topological_distances;  ///< Cached topological distances from Phase 1 for Phase 2 reuse (Jan 2, 2026)

    // Intelligent EEQ matrix caching for performance
    mutable std::unique_ptr<EEQSolverCache> m_eeq_cache;
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
    PARAM(use_iterative_refinement, Bool, false,
          "Use iterative refinement for EEQ Phase 2", "Algorithm", {})
END_PARAMETER_DEFINITION
