/*
 * < EEQ (Electronegativity Equalization) Charge Solver >
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

class CxxThreadPool;  // Forward declaration for pool-based parallelisation

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
 * @brief Linear solve method for EEQ system
 *
 * Controls the algorithm used to solve the augmented EEQ linear system:
 * - LU: PartialPivLU on full augmented matrix (baseline, always works)
 * - SchurCholesky: Exploit SPD structure of NxN Coulomb sub-matrix via Cholesky,
 *   then Schur complement for constraint rows (~2x faster than LU)
 * - PCG: Preconditioned Conjugate Gradient with Schur complement constraint
 *   handling. O(N²·k) where k << N for warm-started MD/optimization (~10-30x for N>500)
 * - Auto: First-call benchmark chooses fastest method (SchurCholesky or PCG),
 *   then caches decision for subsequent calls. Default for optimal performance.
 *
 * Claude Generated - March 2026 (Performance optimization)
 */
enum class EEQSolveMethod {
    LU,             ///< PartialPivLU on full augmented system (baseline)
    SchurCholesky,  ///< Cholesky on NxN SPD block + Schur complement for constraints
    Batched,        ///< Per-fragment Cholesky (GPU only); CPU falls back to SchurCholesky
    PCG,            ///< Preconditioned Conjugate Gradient with warm start
    Auto            ///< Auto-select via first-call benchmark (SchurCholesky vs PCG)
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

    // ===== Convergence Statistics =====

    /**
     * @brief Print summary of PCG convergence statistics since last reset.
     *
     * Prints one line: "EEQ PCG: X/Y converged (worst |r|=Z.ZeZZ)" if any
     * calls did not converge, or nothing if all converged.
     * Resets statistics after printing.
     */
    void printConvergenceSummary();

    /**
     * @brief Reset convergence statistics counters.
     */
    void resetConvergenceStats();

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
        bool use_corrections = false,  // CRITICAL FIX (Jan 4, 2026): default false to match gfnff_final.cpp
        CxxThreadPool* pool = nullptr,  // Claude Generated (WP2, May 2026): pool forwarded to dispatchSolve for Stage-4 batched parallelisation
        int num_threads = 1
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
        const std::optional<Vector>& alpeeq = std::nullopt,  // Claude Generated (January 2026): Charge-dependent alpha
        CxxThreadPool* pool = nullptr,  // Claude Generated (Mar 2026): Pool-based parallelisation
        int num_threads = 1
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

    /**
     * @brief Calculate dgam corrections with full pi-system and amide detection
     *
     * Claude Generated (March 2026): Public interface for computing dgam with
     * complete N-specific corrections (ff=-0.14 for pi-N, ff=-0.16 for amide N).
     * Ensures Coulomb energy parameters are consistent with EEQ solver parameters.
     * Reference: Fortran gfnff_ini.f90:697-724 computes dgam ONCE for both.
     *
     * @param atoms Atomic numbers
     * @param topology_charges Phase-1 charges (qa)
     * @param hybridization Hybridization states
     * @param cn Coordination numbers
     * @param topology Topology information (neighbor lists, etc.)
     * @return Vector of dgam corrections
     */
    Vector calculateDgamFull(
        const std::vector<int>& atoms,
        const Vector& topology_charges,
        const std::vector<int>& hybridization,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    );

    /**
     * @brief Calculate dxi corrections with full environment detection
     *
     * Claude Generated (March 2026): Public interface for computing dxi with
     * complete Fortran-matching corrections (carbene C, free CO, nitro O,
     * polyvalent halogens, pi-system neighbor EN averaging, etc.).
     * Ensures Coulomb energy dxi is consistent with EEQ solver dxi.
     * Reference: Fortran gfnff_ini.f90:358-403
     *
     * @param atoms Atomic numbers
     * @param geometry_bohr Coordinates in Bohr
     * @param cn Coordination numbers
     * @param topology Topology information (neighbor lists, etc.)
     * @return Vector of dxi corrections
     */
    Vector calculateDxiFull(
        const std::vector<int>& atoms,
        const Matrix& geometry_bohr,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    );

    /**
     * @brief Detect amide hydrogens (public wrapper)
     *
     * Claude Generated (March 2026): Public access to amide hydrogen detection
     * for use in Coulomb chi_base correction.
     * Reference: Fortran gfnff_ini.f90:717 - amideH correction to chieeq
     *
     * @param atoms Atomic numbers
     * @param hybridization Hybridization states
     * @param cn Coordination numbers
     * @param topology Topology information
     * @return Vector of amide hydrogen flags
     */
    std::vector<bool> detectAmideHydrogensFull(
        const std::vector<int>& atoms,
        const std::vector<int>& hybridization,
        const Vector& cn,
        const std::optional<TopologyInput>& topology = std::nullopt
    ) const;

    /**
     * @brief Parse solve method string to enum (public so GPU dispatch can reuse).
     * Claude Generated - March 2026; moved to public for WP7-B (May 2026).
     */
    static EEQSolveMethod parseSolveMethod(const std::string& method_str);

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

    /**
     * @brief Generate uniform fallback charges when EEQ solver fails
     * Claude Generated (March 2026): Graceful degradation instead of hard crash
     *
     * @param natoms Number of atoms
     * @param total_charge Total molecular charge
     * @param context Description of failure context (for warning message)
     * @return Vector with q_i = total_charge / natoms
     */
    Vector generateFallbackCharges(int natoms, int total_charge, const std::string& context) const;

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

    /**
     * @brief Compute topological distances via multi-source Dijkstra (O(N·E·logN))
     *
     * Performance replacement for Floyd-Warshall (O(N³)). Uses per-atom Dijkstra
     * on the sparse bond graph with early termination at TDIST_THR cutoff.
     * OpenMP-parallelized over source atoms. Uses float32 arithmetic to match
     * Fortran real(sp) rounding behavior.
     *
     * Claude Generated - March 2026 (Performance optimization)
     *
     * @param atoms Atomic numbers
     * @param topology Topology information (neighbor lists, covalent radii)
     * @return Matrix of topological distances in Bohr (identical to Floyd-Warshall output)
     */
    Matrix computeTopologicalDistancesSparse(
        const std::vector<int>& atoms,
        const TopologyInput& topology
    ) const;

    // ===== Solve Methods =====

    /**
     * @brief Dispatch solve of augmented EEQ system using configured solver method
     * Claude Generated (March 2026): Unified solver dispatch for Phase 1 and Phase 2
     *
     * Routes to SchurCholesky, PCG, LU, or Auto based on m_solve_method.
     * Returns atomic charges (natoms elements), or fallback charges on failure.
     *
     * @param A Full augmented matrix (natoms+nfrag) × (natoms+nfrag)
     * @param x Full RHS vector (natoms+nfrag)
     * @param natoms Number of atoms
     * @param nfrag Number of fragments
     * @param total_charge Total molecular charge (for fallback)
     * @return Vector of atomic charges (natoms elements)
     */
    Vector dispatchSolve(
        const Matrix& A,
        const Vector& x,
        int natoms,
        int nfrag,
        int total_charge,
        CxxThreadPool* pool = nullptr,   // Claude Generated (WP2, May 2026): pool for Stage-4 batched parallelisation
        int num_threads = 1               // Default 1 keeps existing callers (Phase 1, solveEEQ) serial
    );

    /**
     * @brief Solve augmented EEQ system via Schur complement + Cholesky
     *
     * Exploits the SPD structure of the NxN Coulomb sub-matrix A:
     * 1. Cholesky factorize A (O(N³/6) vs O(N³/3) for LU)
     * 2. Solve A·z₁ = b and A·z₂ = Cᵀ via back-substitution
     * 3. Form Schur complement S = C·z₂ (nfrag × nfrag)
     * 4. Solve constraint: λ = S⁻¹·(C·z₁ - d)
     * 5. Final charges: q = z₁ - z₂·λ
     *
     * Falls back to LU if Cholesky fails (matrix not SPD).
     *
     * Claude Generated - March 2026 (Performance optimization)
     *
     * @param A_nn NxN Coulomb+hardness sub-matrix (must be SPD)
     * @param rhs_atoms RHS vector for atoms (N elements)
     * @param C Constraint matrix (nfrag × N)
     * @param rhs_constraints RHS for constraints (nfrag elements)
     * @param natoms Number of atoms
     * @param nfrag Number of fragments
     * @return Charge vector (N elements), or empty vector on failure
     */
    Vector solveWithSchurCholesky(
        const Matrix& A_nn,
        const Vector& rhs_atoms,
        const Matrix& C,
        const Vector& rhs_constraints,
        int natoms,
        int nfrag
    );

    /**
     * @brief Solve A·x = b via Preconditioned Conjugate Gradient
     *
     * Iterative O(N²·k) solver where k is iteration count (typically 10-50 with warm start).
     * Uses Jacobi (diagonal) preconditioner.
     *
     * Claude Generated - March 2026 (Performance optimization)
     *
     * @param A NxN SPD matrix
     * @param b RHS vector
     * @param x0 Initial guess (warm start from previous step)
     * @param max_iter Maximum CG iterations
     * @param tol Convergence tolerance on residual norm
     * @param block_pc Optional block-Jacobi preconditioner (Stage 2, May 2026).
     *                 When non-null, replaces the diagonal Jacobi with per-fragment
     *                 Cholesky inversion — typically reduces k from 30-100 to 2-5.
     * @return Solution vector
     */
    struct BlockJacobiPC;  // Forward decl — defined below
    Vector solveWithPCG(
        const Matrix& A,
        const Vector& b,
        const Vector& x0,
        int max_iter,
        double tol,
        const BlockJacobiPC* block_pc = nullptr
    );

    /**
     * @brief Block-Jacobi preconditioner for multi-fragment EEQ PCG.
     *
     * Stage 2 (May 2026): Replaces the diagonal Jacobi `M_inv = diag(A)^-1` with
     * a block-diagonal `M^-1 = diag(A_ff^-1)`, where A_ff is the sub-block of A_nn
     * for fragment f's atoms. Each block is factorized via dense Cholesky once.
     *
     * Cost:
     *   - Build:  O(Σ atoms_f³)  ≤ O(N · max_frag²)  — negligible for small fragments.
     *   - Apply:  O(Σ atoms_f²)  ≤ O(N · max_frag)   — also negligible.
     *
     * Convergence: PCG sees only the inter-fragment Coulomb coupling as residual,
     * since intra-fragment is exactly inverted. For weakly coupled fragments
     * (e.g. solvated clusters) k typically drops to 2-5 iterations.
     *
     * Falls back to identity (no-op) if any per-fragment block is not SPD.
     */
    struct BlockJacobiPC {
        std::vector<Eigen::LLT<Matrix>> blocks;
        std::vector<std::vector<int>> atom_ids;  // atom_ids[f] = global indices of fragment f's atoms
        bool valid = false;

        /// Apply z = M^-1 r blockwise. Atoms not in any fragment fall back to identity (z = r).
        Vector apply(const Vector& r) const;
    };

    /// Build BlockJacobiPC from A_nn and the fragment constraint matrix C (nfrag × natoms).
    /// Returns a default-constructed (invalid) PC if any fragment block fails Cholesky.
    static BlockJacobiPC buildBlockJacobi(const Matrix& A_nn, const Matrix& C);

    /**
     * @brief Multi-RHS Block-PCG solver — solves A·X = B with B ∈ R^{N × m} simultaneously.
     *
     * Stage 3 (May 2026): Replaces the per-fragment PCG loop in dispatchSolve. The single
     * matvec `A · P` becomes a GEMM (BLAS-3) over an N × m matrix instead of m independent
     * BLAS-2 matvecs, sharing memory bandwidth and benefiting from SIMD/cache reuse.
     *
     * Each column iterates with its own α, β, residual; convergence is checked per column.
     * No deflation — once any column has converged below tolerance, it continues to be
     * updated harmlessly (its residual stays small) until all columns finish or max_iter
     * is reached. Eigen handles the column-wise scalar broadcasting via .array() ops.
     *
     * @param A        N×N SPD matrix (A_nn for EEQ)
     * @param B        N×m matrix of right-hand-sides (column 0 = atomic RHS, cols 1..m-1 = C^T)
     * @param X0       N×m initial guess (zeros or warm-start matrix)
     * @param max_iter Maximum CG iterations
     * @param tol      Convergence tolerance on max-column residual norm
     * @param block_pc Optional block-Jacobi preconditioner (applied column-wise)
     * @return         N×m solution matrix X with A·X ≈ B
     */
    Matrix solveWithPCG_multiRHS(
        const Matrix& A,
        const Matrix& B,
        const Matrix& X0,
        int max_iter,
        double tol,
        const BlockJacobiPC* block_pc = nullptr
    );

    // ===== Configuration =====

    ConfigManager m_config;           ///< Configuration manager
    int m_max_iterations;             ///< Maximum iterations for Phase 2 refinement
    double m_convergence_threshold;   ///< Convergence threshold for charge changes (e)
    int m_verbosity;                  ///< Verbosity level (0-3)
    bool m_calculate_cn;              ///< Auto-calculate CN if not provided
    bool m_allow_unconverged;         ///< Allow continuing with unconverged charges (Claude Generated Apr 2026)
    bool m_skip_phase2;               ///< Skip Phase 2 and use Phase 1 topology charges directly (Claude Generated Apr 2026)
    EEQSolveMethod m_solve_method;    ///< Linear solve algorithm selection

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

    // PCG warm-start cache for iterative EEQ solve
    // Claude Generated - March 2026 (Performance optimization)
    // May 2026: m_pcg_last_Z2 promoted from Vector to Matrix so all nfrag constraint
    // solves warm-start from their previous result, not just fragment 0.
    mutable Vector m_pcg_last_z1;  ///< Previous A⁻¹·b solution (warm start for PCG)
    mutable Matrix m_pcg_last_Z2;  ///< Previous A⁻¹·C^T columns (N × nfrag), one warm start per fragment
    mutable bool m_pcg_cache_valid = false;  ///< Whether PCG warm-start cache is usable

    // Auto-solver benchmark state
    // Claude Generated - March 2026 (Auto-solver selection)
    mutable EEQSolveMethod m_selected_method = EEQSolveMethod::SchurCholesky;  ///< Method chosen by auto-benchmark
    mutable bool m_auto_benchmark_done = false;  ///< Whether auto-benchmark has run

    // WP7-B (May 2026): warn-once flag for CPU "batched" → cholesky fallback
    mutable bool m_batched_cpu_warned = false;

    // Convergence statistics — accumulate across calls, print summary instead of per-call spam
    // Claude Generated (April 2026)
    mutable int m_pcg_total_calls = 0;       ///< Total PCG solve calls
    mutable int m_pcg_nonconv_calls = 0;    ///< PCG calls that did not converge
    mutable int m_pcg_total_iters = 0;      ///< Total PCG iterations across all calls
    mutable double m_pcg_worst_residual = 0.0; ///< Worst |r| among non-converged calls

    // Claude Generated (Apr 2026): Cache for historically implausible Phase 2 results.
    // Once Phase 2 has produced garbage charges for a given system size, skip it forever
    // (until the atom count changes, which indicates a new molecule).
    mutable bool m_phase2_historically_implausible = false;
    mutable int m_phase2_implausible_natoms = 0;

    // ===== Pre-allocated Buffers for calculateFinalCharges (Claude Generated Mar 2026) =====
    // Avoid ~26 MB alloc+free per gradient step for large molecules (N=1280).
    // Buffers are allocated once and reused when atom count stays the same.
    mutable Matrix m_phase2_distances;   ///< N×N distance buffer
    mutable Matrix m_phase2_A;           ///< (N+nfrag)×(N+nfrag) augmented matrix buffer
    mutable Vector m_phase2_rhs;         ///< (N+nfrag) RHS vector buffer
    mutable int m_phase2_buf_natoms = 0; ///< Atom count for current buffer size
    mutable int m_phase2_buf_nfrag = 0;  ///< Fragment count for current buffer size

    /// Last truly successful charges (Phase 2 if available, otherwise Phase 1).
    /// Initialized from topology_charges on first call, then overwritten with the Phase 2
    /// result after a successful solve (line 2919). NOT overwritten on plausibility-check
    /// failure — the cache preserves the last good result. Used as fallback when the solver
    /// fails on subsequent steps. Validated on read: allFinite(), |q|_max < 50, charge sum
    /// within tolerance, to avoid using a corrupted cache.
    mutable Vector m_last_successful_charges;

    /// Validate cached charges are plausible for use as fallback
    bool isValidChargeCache(const Vector& charges, int natoms, int total_charge) const {
        return charges.size() == natoms
            && charges.allFinite()
            && std::abs(charges.sum() - total_charge) < 1.0
            && charges.cwiseAbs().maxCoeff() < 50.0;
    }


    /// Ensure buffers are large enough. Only reallocates if size changed.
    void ensurePhase2Buffers(int natoms, int nfrag) const {
        if (natoms != m_phase2_buf_natoms || nfrag != m_phase2_buf_nfrag) {
            int m = natoms + nfrag;
            m_phase2_distances.resize(natoms, natoms);
            m_phase2_A.resize(m, m);
            m_phase2_rhs.resize(m);
            m_phase2_buf_natoms = natoms;
            m_phase2_buf_nfrag = nfrag;
        }
    }
};

// ===== Parameter Definitions =====

BEGIN_PARAMETER_DEFINITION(eeq_solver)
    PARAM(accuracy, String, "normal", "Accuracy profile: loose|normal|medium|high. Maps to solver tolerances and iteration limits.", "Basic", {})
    PARAM(allow_unconverged_charges, Bool, false, "Allow calculation to continue with unconverged charges (warn instead of abort).", "Advanced", {})
    PARAM(skip_phase2, Bool, false, "Skip Phase 2 EEQ refinement and use Phase 1 topology charges directly. Faster but less accurate.", "Advanced", {})
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
    PARAM(solve_method, String, "cholesky",
          "EEQ linear solve algorithm: cholesky | batched | pcg | auto | lu (legacy). "
          "GPU paths: cholesky → WP5-A/WP7-A (exact, full N×N Cholesky); "
          "batched → WP7-B (per-fragment Cholesky, ignores cross-fragment Coulomb — "
          "use only for well-separated fragments, see eeq_batched_min_distance). "
          "CPU 'batched' logs a warning and falls back to cholesky. "
          "Legacy alias 'schur_cholesky' maps to 'cholesky'.",
          "Algorithm", {})
    PARAM(max_pcg_iterations, Int, 200,
          "Maximum PCG iterations for EEQ solve", "Algorithm", {})
    PARAM(pcg_tolerance, Double, 1e-10,
          "PCG convergence tolerance", "Algorithm", {})
    PARAM(pcg_large_system_iterations, Int, 5000,
          "Max PCG iterations for large systems (N>pcg_large_threshold). Overrides max_pcg_iterations.", "Algorithm", {})
    PARAM(pcg_large_system_scaling, Int, 10,
          "PCG iteration scaling factor for large systems: max_iter = min(scaling*N, pcg_large_system_iterations)", "Algorithm", {})
    PARAM(pcg_large_system_tol_factor, Double, 1e-6,
          "Relative tolerance factor for large systems: tol = max(pcg_tolerance, factor*||rhs||)", "Algorithm", {})
    PARAM(pcg_large_threshold, Int, 500,
          "Atom count threshold above which PCG auto-selection and adaptive scaling activate", "Algorithm", {})
    PARAM(eeq_pcg_nfrag_threshold, Int, 4,
          "Auto solver: prefer SchurCholesky when nfrag >= this threshold. "
          "Cholesky factorizes A_nn once and amortizes the cost over nfrag back-substitutions "
          "via BLAS-3, which beats nfrag independent PCG solves once nfrag is moderate.",
          "Algorithm", {})
    PARAM(eeq_pcg_expected_iters, Int, 30,
          "Auto solver: expected average PCG iteration count. Used in the cost model "
          "k_pcg*(nfrag+1) > N/6 to decide between PCG and SchurCholesky.",
          "Algorithm", {})
    PARAM(eeq_distance_cutoff, Double, 0.0,
          "Distance cutoff in Bohr for Coulomb matrix sparsification (0 = no cutoff, matches Fortran goed_gfnff). Non-zero values violate Hellmann-Feynman vs. the full Coulomb energy and degrade MD energy conservation.", "Advanced", {})
    PARAM(eeq_batched_min_distance, Double, 15.0,
          "Minimum atom-atom distance between different fragments (Bohr) below which "
          "the batched EEQ solver logs a warning. Batched still runs — the warning "
          "alerts the user that cross-fragment Coulomb (which batched ignores) may be "
          "non-negligible. 0 = no warning. Default 15 Bohr ≈ 8 Å.",
          "Advanced", {})
    PARAM(dump_charges, Bool, false,
          "Save Phase 1 and Phase 2 charges to charges_dump_N<size>.json for analysis", "Advanced", {})
END_PARAMETER_DEFINITION
