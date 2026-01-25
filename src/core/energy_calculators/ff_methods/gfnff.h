/*
 * <GFN-FF Implementation for Curcuma>
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
 * GFN-FF (Geometry, Frequency, Noncovalent - Force Field) is a fully
 * automated, quantum chemistry-based force field for the accurate
 * description of structures and dynamics of large molecular systems.
 */

#pragma once

#include "json.hpp"
#include "src/core/config_manager.h"
#include "src/core/energy_calculators/ff_methods/forcefield.h"
#include "src/core/energy_calculators/ff_methods/eeq_solver.h"  // EEQ charge calculation (Dec 2025 - Phase 3)
#include "src/core/energy_calculators/ff_methods/huckel_solver.h"  // Full Hückel calculation (Jan 2026 - Phase 1)
#include "src/core/global.h"
#include "src/core/functional_groups.h"
#include "src/core/periodic_table.h"
#include <utility>
#include <optional>
#include <vector>
#include <memory>

using json = nlohmann::json;

/**
 * @brief GFN-FF Implementation as Standalone Force Field
 *
 * GFN-FF combines quantum chemical accuracy with force field efficiency.
 * It provides:
 * - Automatic parametrization based on extended tight-binding (xTB)
 * - Accurate treatment of non-covalent interactions
 * - Gradients for geometry optimization and dynamics
 * - Coverage of the periodic table up to Z=86
 *
 * References:
 * - Spicher, S.; Grimme, S. "Robust Atomistic Modeling of Materials,
 *   Organometallic, and Biochemical Systems" Angew. Chem. Int. Ed. 59, 15665 (2020)
 *
 * Note: GFNFF is semantically a force field, not a quantum method.
 * It was previously inheriting from QMInterface but has been refactored
 * to standalone class for semantic correctness (Phase 2, November 2025).
 */

/**
 * @brief Triangular indexing function for symmetric matrices
 *
 * Claude Generated (January 10, 2026) - Port from Fortran gfnff_helpers.f90:416-423
 * Converts atom pair (i,j) to linear index for upper-triangular storage
 *
 * @param i First atom index (0-based)
 * @param j Second atom index (0-based)
 * @return Linear index in triangular array
 *
 * Formula (C++ 0-based): lin = min(i,j) + max(i,j)*(max(i,j)+1)/2
 */
inline int lin(int i, int j) {
    int imax = std::max(i, j);
    int imin = std::min(i, j);
    return imin + imax * (imax + 1) / 2;
}

/**
 * @brief Check if an atom is classified as a metal
 *
 * Claude Generated (January 2026) - Phase 4: Metal scaling
 * Reference: gfnff_method.cpp:4382-4386
 *
 * @param atomic_number Atomic number (Z)
 * @return true if atom is a transition metal or lanthanide/actinide
 *
 * Metal classification: Transition metals (Sc-Zn, Y-Cd, La-Hg) + Lanthanides/Actinides
 */
inline bool isMetalAtom(int atomic_number) {
    return (atomic_number >= 21 && atomic_number <= 30) ||   // Sc-Zn
           (atomic_number >= 39 && atomic_number <= 48) ||   // Y-Cd
           (atomic_number >= 57 && atomic_number <= 80) ||   // La-Hg
           (atomic_number >= 89 && atomic_number <= 103);    // Ac-Lr
}

class GFNFF {
public:
    /**
     * @brief Phase 9: Topology information structure (moved here for use in function signatures)
     *
     * Extended for Session 5 (December 2025): Two-phase EEQ system
     * - neighbor_lists: Full neighbor connectivity for functional group detection
     * - functional_groups: Classification of atoms into functional group types
     * - topology_charges: Phase 1 EEQ charges (qa) - used for correction calculations
     * - dxi: Electronegativity corrections per atom
     * - dgam: Hardness corrections per atom (already computed in Phase 1)
     * - dalpha: Polarizability corrections per atom
     */
    struct TopologyInfo {
        Vector coordination_numbers;
        Vector neighbor_counts;  // Simple neighbor counts for XTB-compatible torsion CN correction
        std::vector<int> hybridization;
        std::vector<int> pi_fragments;
        std::vector<int> ring_sizes;
        Vector eeq_charges;
        std::vector<bool> is_metal;
        std::vector<bool> is_aromatic;

        // NEW (Session 5): Topology and correction data
        std::vector<std::vector<int>> neighbor_lists;        // Full neighbor connectivity
        std::vector<FunctionalGroupType> functional_groups;  // Per-atom functional group classification
        Vector topology_charges;                             // Phase 1 EEQ charges (qa) - base topology
        Vector dxi;                                          // Electronegativity corrections
        Vector dgam;                                         // Hardness corrections
        Vector dalpha;                                       // Polarizability corrections

        // Claude Generated (January 2026): Charge-dependent alpha for EEQ
        // Calculated ONCE during topology init using Phase 1 charges
        // Formula: alpeeq(i) = (alpha_base + ff*qa(i))²
        // Reference: Fortran gfnff_ini.f90:718-725, gfnff_data_types.f90:128
        Vector alpeeq;                                       // Squared alpha values, charge-corrected

        // Phase 2.1: Distance caching (Claude Generated - Dec 2025)
        // Computed once per geometry update to eliminate redundant sqrt() calls
        Eigen::MatrixXd distance_matrix;        // N×N distances in Bohr
        Eigen::MatrixXd squared_dist_matrix;    // N×N squared distances (avoids sqrt)

        // Phase 2.2: Adjacency list (Claude Generated - Dec 2025)
        // Per-atom neighbor connectivity built from bond list
        // Used in generateGFNFFAngles() to avoid O(N_bonds) search per atom
        std::vector<std::vector<int>> adjacency_list;

        // Phase 9B: Topological distances (Claude Generated - Dec 24, 2025)
        // N×N matrix of bond counts (shortest path) between atom pairs
        // Used for 1,3 and 1,4 topology factors in non-bonded repulsion
        // topo_distances[i][j] = number of bonds in shortest path between i and j
        //   0 = same atom, 1 = bonded, 2 = separated by 1 bond, 3 = 1,3-pair, 4 = 1,4-pair, etc.
        std::vector<std::vector<int>> topo_distances;

        // Bond type classification (Claude Generated - Jan 2, 2026)
        // Per-bond type following Fortran gfnff_ini.f90:1131-1148
        // btyp = 1: Single bond (default)
        // btyp = 2: Pi bond (sp2-sp2)
        // btyp = 3: Sp bond (linear, no torsion)
        // btyp = 4: Hypervalent
        // btyp = 5: Metal-containing bond
        // btyp = 6: Eta-complex (special metal)
        // btyp = 7: TM metal-metal bond
        std::vector<int> bond_types;

        // Phase 2C: π-bond orders (Claude Generated - January 10, 2026)
        // Simplified approximation based on hybridization (no full Hückel calculation)
        // Stored in triangular format using lin(i,j) indexing
        // pbo[lin(i,j)] = π-bond order between atoms i and j
        std::vector<double> pi_bond_orders;

        // PERFORMANCE OPTIMIZATION (Claude Generated - January 17, 2026)
        // Pre-cached EEQ parameters per atom to avoid repeated lookups in O(N²) loops
        // Individual vectors instead of struct array (struct defined later in class)
        // Populated once during topology initialization, indexed by atom index
        std::vector<double> eeq_chi;  // Electronegativity per atom
        std::vector<double> eeq_gam;  // Chemical hardness per atom
        std::vector<double> eeq_alp;  // Damping parameter (squared) per atom
        std::vector<double> eeq_cnf;  // CN correction factor per atom

        // BF (Bonded ATM/GFN-FF) - Claude Generated (January 17, 2026)
        // GFN-FF bonded ATM (batm) terms - D3-like 3-body dispersion for 1,4-pairs only
        // Reference: external/gfnff/src/gfnff_ini.f90:745-779, gfnff_engrad.F90:562-603
        std::vector<std::vector<int>> bpair;  // N×N topological distance matrix (bonds between atoms)
        std::vector<std::tuple<int,int,int>> b3list;  // Batm triples (i,j,k) for 1,4-pairs
        int nbatm = 0;  // Number of batm triples
    };

    /**
     * @brief Default constructor
     */
    GFNFF();

    /**
     * @brief Constructor with custom parameters
     * @param parameters JSON configuration for GFN-FF
     */
    explicit GFNFF(const json& parameters);

    /**
     * @brief Destructor
     */
    virtual ~GFNFF();

    /**
     * @brief Initialize molecule for GFN-FF calculation from Mol object
     * @param molecule Molecule to initialize
     * @return true if initialization successful
     */
    bool InitialiseMolecule(const Mol& molecule);

    /**
     * @brief Initialize molecule for GFN-FF calculation (parameterless)
     * @return true if initialization successful
     */
    bool InitialiseMolecule();

    /**
     * @brief Update molecular geometry
     * @param geometry New geometry matrix
     * @return true if update successful
     */
    bool UpdateMolecule(const Matrix& geometry);

    /**
     * @brief Update molecular geometry (parameterless)
     * @return true if update successful
     */
    bool UpdateMolecule();

    /**
     * @brief Perform GFN-FF calculation
     * @param gradient Calculate gradients if true
     * @return Total energy in Hartree
     */
    double Calculation(bool gradient = false);

    /**
     * @brief Get analytical gradients
     * @return Gradient matrix (N_atoms x 3) in Hartree/Bohr
     */
    Geometry Gradient() const { return m_gradient; }

    /**
     * @brief Check if gradients are available
     * @return true (GFN-FF always provides gradients)
     */
    bool hasGradient() const { return true; }

    /**
     * @brief Get atomic partial charges (Phase 2 energy charges - nlist%q)
     * @return Vector of atomic charges (final energy charges)
     *
     * Returns Phase 2 EEQ charges used for energy calculation.
     * These are the charges used in gradient calculations and for Coulomb energy.
     */
    Vector Charges() const;

    /**
     * @brief Get topology charges (Phase 1 - topo%qa)
     * @return Vector of topology charges
     *
     * Returns Phase 1 EEQ topology charges used for parameter generation.
     * These charges use integer neighbor count and are used to compute
     * corrections for bonds, angles, and other topological terms.
     *
     * Claude Generated (January 4, 2026)
     */
    Vector getTopologyCharges() const;

    /**
     * @brief Get energy charges (Phase 2 - nlist%q)
     * @return Vector of energy charges
     *
     * Alias for Charges() - returns Phase 2 EEQ charges.
     * Provided for clarity when comparing both charge types.
     *
     * Claude Generated (January 4, 2026)
     */
    Vector getEnergyCharges() const { return Charges(); }

    /**
     * @brief Get bond orders (Wiberg bond orders)
     * @return Vector of bond orders
     */
    Vector BondOrders() const;

    /**
     * @brief Set calculation parameters
     * @param parameters JSON configuration
     */
    void setParameters(const json& parameters);

    /**
     * @brief Get current parameters
     * @return JSON configuration
     */
    json getParameters() const { return m_parameters; }

private:
    /**
     * @brief Initialize GFN-FF force field and generate parameters
     * @return true if successful
     */
    bool initializeForceField();

    /**
     * @brief Generate GFN-FF specific force field parameters
     * @return JSON with GFN-FF parameters
     */
    json generateGFNFFParameters();

    /**
     * @brief Calculate topology and connectivity for GFN-FF
     * @return true if successful
     */
    bool calculateTopology();

    /**
     * @brief Retrieve cached topology information, computing it once if needed
     */
    const TopologyInfo& getCachedTopology() const;

    /**
     * @brief Retrieve cached bond list, computing it once if needed
     */
    const std::vector<std::pair<int,int>>& getCachedBondList() const;

    /**
     * @brief Calculate topological distances (bond counts) between all atom pairs using BFS
     * @param adjacency_list Per-atom neighbor connectivity
     * @return N×N matrix of shortest path lengths (0=same, 1=bonded, 3=1,3-pair, 4=1,4-pair)
     *
     * Claude Generated (Dec 24, 2025): Breadth-First Search for 1,3/1,4 topology factors
     */
    std::vector<std::vector<int>> calculateTopologyDistances(const std::vector<std::vector<int>>& adjacency_list) const;

    /**
     * @brief Classify bond type according to GFN-FF topology rules
     * @param atom_i First atom index
     * @param atom_j Second atom index
     * @param hyb_i Hybridization of atom i (0=sp3, 1=sp, 2=sp2, 3=terminal, 5=hypervalent)
     * @param hyb_j Hybridization of atom j
     * @param is_metal_i True if atom i is a metal
     * @param is_metal_j True if atom j is a metal
     * @return Bond type (btyp): 1=single, 2=pi, 3=sp/linear, 4=hypervalent, 5=metal, 6=eta, 7=TM-TM
     *
     * Claude Generated (Jan 2, 2026): Ported from Fortran gfnff_ini.f90:1131-1148
     * Used for extra torsion filtering (btyp < 5 excludes metal bonds)
     */
    int classifyBondType(int atom_i, int atom_j, int hyb_i, int hyb_j,
                         bool is_metal_i, bool is_metal_j) const;

    /**
     * @brief Validate molecular structure for GFN-FF
     * @return true if molecule is valid
     */
    bool validateMolecule() const;

    /**
     * @brief Convert energy from kcal/mol to Hartree
     * @param energy Energy in kcal/mol
     * @return Energy in Hartree
     */
    double convertToHartree(double energy) const;

    /**
     * @brief Convert gradient from kcal/(mol*Angstrom) to Hartree/Bohr
     * @param gradient Gradient in kcal/(mol*Angstrom)
     * @return Gradient in Hartree/Bohr
     */
    Matrix convertGradientToHartree(const Matrix& gradient) const;

    /**
     * @brief Generate GFN-FF bond parameters from bond detection
     * @return JSON array of bond parameters
     */
    json generateGFNFFBonds() const;

    /**
     * @brief Generate GFN-FF angle parameters from topology
     * @param topo_info Topology information including charges, hybridization, CN, etc.
     * @return JSON array of angle parameters
     */
    json generateGFNFFAngles(const TopologyInfo& topo_info) const;

    /**
     * @brief Generate GFN-FF torsion parameters from topology
     *
     * Claude Generated (2025): Ported from Grimme Lab GFN-FF (Spicher & Grimme 2020)
     * Reference: external/gfnff/src/gfnff_engrad.F90:1041-1122 (egtors subroutine)
     *
     * Implements proper and improper torsion potentials for molecular rotations.
     * See docs/theory/GFNFF_TORSION_THEORY.md for scientific background.
     *
     * @return JSON array of torsion parameters
     */
    json generateGFNFFTorsions() const;

    /**
     * @brief Generate GFN-FF inversion/out-of-plane parameters
     *
     * Claude Generated (2025): Inversion term implementation
     * Reference: external/gfnff/src/gfnff_helpers.f90:427-510 (omega, domegadr)
     *
     * Implements out-of-plane bending potentials for sp² centers (planarity constraints).
     * See docs/theory/GFNFF_INVERSION_THEORY.md for scientific background.
     *
     * @return JSON array of inversion parameters
     */
    json generateGFNFFInversions() const;

    // Phase 4.2: GFN-FF pairwise non-bonded parameter generation (Claude Generated 2025)

    /**
     * @brief Generate EEQ-based Coulomb electrostatics pairwise parameters
     * Formula: E_coul = q_i * q_j * erf(γ_ij * r_ij) / r_ij
     * @return JSON array of Coulomb pair parameters
     */
    json generateGFNFFCoulombPairs() const;

    /**
     * @brief Generate GFN-FF repulsion pairwise parameters
     * Formula: E_rep = repab * exp(-α*r^1.5) / r
     * @return JSON array of repulsion pair parameters
     */
    json generateGFNFFRepulsionPairs() const;

    /**
     * @brief Generate D3/D4 dispersion pairwise parameters with BJ damping
     * Formula: E_disp = -Σ_ij f_damp(r) * (s6*C6/r^6 + s8*C8/r^8)
     * @return JSON array of dispersion pair parameters
     *
     * Claude Generated (December 2025): D3/D4 integration
     * - Calls D4ParameterGenerator if USE_D4 defined (preferred)
     * - Falls back to D3ParameterGenerator if USE_D3 defined
     * - Final fallback to generateFreeAtomDispersion()
     */
    json generateGFNFFDispersionPairs() const;

    /**
     * @brief Generate dispersion parameters using free-atom C6 approximation
     * @return JSON array of dispersion pair parameters
     *
     * Claude Generated (December 2025): Fallback method
     * Legacy implementation extracted from generateGFNFFDispersionPairs().
     * Uses hardcoded free-atom C6 coefficients (geometry-independent).
     * Less accurate than D3/D4 but always available.
     */
    json generateFreeAtomDispersion() const;

    /**
     * @brief Factory method to generate D3 dispersion parameters
     * @return JSON array of D3 dispersion pair parameters
     *
     * Claude Generated (December 2025): Phase 3 - Factory method
     * Encapsulates all D3-specific generation logic with fallback handling.
     * Creates D3ParameterGenerator, runs GenerateParameters(), and converts
     * output to GFN-FF dispersion pair format.
     * Handles exceptions and falls back to free-atom C6 if D3 fails.
     */
    json generateD3Dispersion() const;

    /**
     * @brief Extract D3/D4 configuration from main GFN-FF config
     * @param method Dispersion method: "d3" or "d4"
     * @return ConfigManager for D3/D4 parameter generator
     *
     * Claude Generated (December 2025): Configuration helper
     * Extracts relevant parameters (s6, s8, a1, a2) from main config
     * and creates ConfigManager for D3ParameterGenerator or D4ParameterGenerator.
     */
    ConfigManager extractDispersionConfig(const std::string& method) const;

    /**
     * @brief Get covalent radius for element
     * @param atomic_number Element atomic number
     * @return Covalent radius in Angstrom
     */
    double getCovalentRadius(int atomic_number) const;

    // GFN-FF parameter structures
    struct GFNFFBondParams {
        double force_constant;        // k_b in Fortran (energy scale)
        double equilibrium_distance;  // r₀ reference bond length
        double alpha;                 // α exponential decay parameter (was: anharmonic_factor)
        double rabshift;              // Claude Generated (Dec 2025): vbond(1) = gen%rabshift + shift
        double fqq;                   // Claude Generated (Jan 7, 2026): charge-dependent force constant factor

        // Claude Generated (Jan 18, 2026): Dynamic r0 calculation parameters
        // Reference: Fortran gfnff_rab.f90:147-153 - r0 recalculated at each Calculate()
        // Formula: r0 = (r0_base_i + cnfak_i*cn_i + r0_base_j + cnfak_j*cn_j + rabshift) * ff
        int z_i = 0, z_j = 0;           // Atomic numbers for parameter lookup
        double r0_base_i = 0.0;          // r0_gfnff[z_i-1] (Bohr)
        double r0_base_j = 0.0;          // r0_gfnff[z_j-1] (Bohr)
        double cnfak_i = 0.0;            // cnfak_gfnff[z_i-1]
        double cnfak_j = 0.0;            // cnfak_gfnff[z_j-1]
        double ff = 1.0;                 // EN-correction: 1 - k1*|ΔEN| - k2*ΔEN²
    };

    struct GFNFFAngleParams {
        double force_constant;     // k_ijk in Fortran
        double equilibrium_angle;  // θ₀ reference angle
        // Phase 1.3: Removed c0,c1,c2 Fourier coefficients (were dummy values)
        // GFN-FF uses simple angle bending, not Fourier expansion
    };

    /**
     * @brief GFN-FF torsion parameters
     *
     * Claude Generated (2025): Based on GFN-FF method (Spicher & Grimme 2020)
     *
     * Torsion potential: E = V/2 * [1 - cos(n*(φ - φ₀))] * D(r_ij, r_jk, r_kl)
     *
     * Scientific background:
     * - Describes rotation around central bond j-k in i-j-k-l sequence
     * - Periodicity n determines symmetry (1, 2, or 3)
     * - Barrier height V controls rotation difficulty
     * - Distance damping D couples stretching with rotation
     *
     * Reference: docs/theory/GFNFF_TORSION_THEORY.md
     */
    struct GFNFFTorsionParams {
        double barrier_height;     ///< V_n: Energy barrier in kcal/mol (Spicher & Grimme Eq. 8)
        int periodicity;            ///< n: Rotational symmetry (1, 2, or 3)
        double phase_shift;         ///< φ₀: Reference angle in radians
        bool is_improper;          ///< True for out-of-plane/improper torsions
    };

    /**
     * @brief GFN-FF inversion/out-of-plane parameters
     *
     * Claude Generated (2025): Based on GFN-FF method (Spicher & Grimme 2020)
     *
     * Inversion potential: E = V * [cos(ω) - cos(ω₀)]² * D(r_ij, r_jk, r_jl)
     *
     * Scientific background:
     * - Describes out-of-plane bending for atom i relative to plane j-k-l
     * - Enforces planarity at sp² centers (aromatics, C=C, C=O)
     * - ω = out-of-plane angle ∈ [-π/2, +π/2]
     * - Double-well potential allows ±ω₀ equivalence
     *
     * Reference: docs/theory/GFNFF_INVERSION_THEORY.md
     */
    struct GFNFFInversionParams {
        double barrier_height;     ///< V: Energy barrier in kcal/mol
        double reference_angle;     ///< ω₀: Reference angle in radians (usually 0 for planar)
        int potential_type;        ///< 0: double-well [cos(ω)-cos(ω₀)]², 1: single-well [1-cos(ω)]
    };

    /**
     * @brief EEQ (Electronegativity Equalization) parameters
     *
     * Claude Generated (2025): Phase 3 EEQ charge calculation
     *
     * Parameters from gfnff_param.f90 (angewChem2020 parameter set)
     * Reference: S. Spicher, S. Grimme, Angew. Chem. Int. Ed. 2020, 59, 15665-15673
     *
     * Scientific background:
     * - EEQ method solves linear system A·q = b for atomic charges
     * - chi: atomic electronegativity (controls charge distribution)
     * - gam: chemical hardness (resistance to charge transfer)
     * - alp: damping parameter for Coulomb interaction
     * - cnf: coordination number correction factor
     */
    struct EEQParameters {
        double chi;  ///< Electronegativity (angewChem2020)
        double gam;  ///< Chemical hardness (angewChem2020)
        double alp;  ///< Damping parameter (angewChem2020)
        double cnf;  ///< CN correction factor (angewChem2020)
        double xi_corr;  ///< Environment correction (topology-dependent, optional)
    };

    /**
     * @brief Get GFN-FF bond parameters for element pair
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of second atom
     * @param distance Current bond distance
     * @return GFN-FF bond parameters
     */
    /**
     * @brief Get GFN-FF bond parameters with full topology corrections (Phase 9)
     * @param atom1 First atom index
     * @param atom2 Second atom index
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of second atom
     * @param distance Current bond distance
     * @param topo Topology information (CN, hyb, charges, rings)
     * @return GFN-FF bond parameters with all corrections
     */
    GFNFFBondParams getGFNFFBondParameters(int atom1, int atom2, int z1, int z2,
                                            double distance, const TopologyInfo& topo) const;

    /**
     * @brief Load atomic charges from reference GFN-FF calculation
     * @return true if charges loaded successfully
     */
    bool loadGFNFFCharges();

    /**
     * @brief Get EEQ parameters for an element
     *
     * Claude Generated (2025): Phase 3 EEQ implementation
     * Reference: gfnff_param.f90 (chi/gam/alp/cnf_angewChem2020)
     *
     * Returns EEQ parameters from angewChem2020 parameter set:
     * - chi: atomic electronegativity
     * - gam: chemical hardness
     * - alp: damping parameter for erf(γ*r)/r Coulomb interaction
     * - cnf: coordination number correction factor
     *
     * @param atomic_number Element atomic number (1-86)
     * @return EEQ parameters structure
     */
    EEQParameters getEEQParameters(int atomic_number) const;

    /**
     * @brief Get GFN-FF angle parameters for angle triplet
     * @param atom_i First atom index (needed for hybridization lookup)
     * @param atom_j Center atom index (determines equilibrium angle via hybridization)
     * @param atom_k Third atom index
     * @param current_angle Current angle in radians (used for geometry-dependent overrides)
     * @return GFN-FF angle parameters with topology-based equilibrium angles
     *
     * Claude Generated (Nov 2025): Phase 2 implementation uses topology-aware parameters
     * including charge-dependent corrections (fqq), coordination number scaling (fn),
     * element-specific corrections (f2), and small-angle corrections (fbsmall).
     */
    GFNFFAngleParams getGFNFFAngleParameters(int atom_i, int atom_j, int atom_k,
                                              double current_angle, const TopologyInfo& topo_info) const;

    /**
     * @brief Get GFN-FF torsion parameters for atom quartet
     *
     * Claude Generated (2025): Topology-aware parameter assignment
     * Reference: external/gfnff/src/gfnff_ini2.f90 (torsion setup)
     *
     * Assigns torsion parameters based on:
     * - Hybridization of central atoms j and k
     * - Ring membership (strain corrections)
     * - Conjugation status (planarity preferences)
     *
     * @param z_i Atomic number of first atom
     * @param z_j Atomic number of second atom (central bond)
     * @param z_k Atomic number of third atom (central bond)
     * @param z_l Atomic number of fourth atom
     * @param hyb_j Hybridization of atom j (1=sp, 2=sp2, 3=sp3)
     * @param hyb_k Hybridization of atom k
     * @return GFN-FF torsion parameters
     */
    GFNFFTorsionParams getGFNFFTorsionParameters(int z_i, int z_j, int z_k, int z_l,
                                                  int hyb_j, int hyb_k,
                                                  double qa_j = 0.0, double qa_k = 0.0,
                                                  double cn_i = 2.0, double cn_l = 2.0,
                                                  bool in_ring = false, int ring_size = 0,
                                                  int i_atom_idx = -1, int j_atom_idx = -1,
                                                  int k_atom_idx = -1, int l_atom_idx = -1) const;

    /**
     * @brief Calculate dihedral angle for four atoms
     *
     * Claude Generated (2025): Standard computational chemistry formula
     * Reference: Allen & Tildesley "Computer Simulation of Liquids" (1987)
     *
     * Computes signed dihedral angle φ ∈ [-π, π] between planes i-j-k and j-k-l.
     * Uses atan2 for proper sign handling (crucial for gradient calculation).
     *
     * Physical interpretation:
     * - φ = 0°: cis/eclipsed (atoms i and l on same side)
     * - φ = 180°: trans/anti (atoms i and l on opposite sides)
     *
     * @param i Index of first atom
     * @param j Index of second atom
     * @param k Index of third atom
     * @param l Index of fourth atom
     * @return Dihedral angle in radians [-π, π]
     */
    double calculateDihedralAngle(int i, int j, int k, int l) const;

    /**
     * @brief Calculate derivatives of dihedral angle w.r.t. atomic positions
     *
     * Claude Generated (2025): Analytical gradient for torsions
     * Reference: external/gfnff/src/gfnff_engrad.F90 (dphidr subroutine)
     *
     * Computes ∂φ/∂x_i, ∂φ/∂x_j, ∂φ/∂x_k, ∂φ/∂x_l for chain rule in gradient.
     * Critical for analytical force calculation in molecular dynamics.
     *
     * @param i Index of first atom
     * @param j Index of second atom
     * @param k Index of third atom
     * @param l Index of fourth atom
     * @param phi Current dihedral angle (from calculateDihedralAngle)
     * @param dda Output: ∂φ/∂x_i (3D vector)
     * @param ddb Output: ∂φ/∂x_j (3D vector)
     * @param ddc Output: ∂φ/∂x_k (3D vector)
     * @param ddd Output: ∂φ/∂x_l (3D vector)
     */
    void calculateDihedralGradient(int i, int j, int k, int l, double phi,
                                     Vector& dda, Vector& ddb,
                                     Vector& ddc, Vector& ddd) const;

    /**
     * @brief Calculate damping function for torsion potential
     *
     * Claude Generated (2025): Distance-dependent damping
     * Reference: external/gfnff/src/gfnff_engrad.F90 (gfnffdampt function)
     *
     * Damping function couples bond stretching with torsional motion:
     * D(r) = 1 / [1 + exp(-α*(r/r₀ - 1))]
     *
     * Physical meaning:
     * - Stretched bonds → weaker torsion barrier (easier rotation)
     * - Compressed bonds → stronger torsion barrier
     *
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of second atom
     * @param r_squared Squared distance r² in Bohr²
     * @param damp Output: Damping value D(r)
     * @param damp_deriv Output: Derivative ∂D/∂r
     */
    void calculateTorsionDamping(int z1, int z2, double r_squared,
                                  double& damp, double& damp_deriv) const;

    // =================================================================================
    // INVERSION/OUT-OF-PLANE HELPER FUNCTIONS (Phase 1.2)
    // =================================================================================

    /**
     * @brief Calculate out-of-plane angle (omega) for atom i relative to plane j-k-l
     *
     * Claude Generated (2025): Inversion angle calculation
     * Reference: external/gfnff/src/gfnff_helpers.f90:427-448 (omega function)
     *
     * Computes ω ∈ [-π/2, +π/2] measuring deviation from planarity:
     *   ω = arcsin(n · v̂)
     * where:
     *   n = (r_ij × r_jk) / |r_ij × r_jk|  (normal to plane i-j-k)
     *   v = r_il  (vector from i to l)
     *
     * Physical interpretation:
     * - ω = 0: atom i in plane j-k-l (planar, typical for sp²)
     * - ω = ±π/2: atom i perpendicular to plane (pyramidal)
     *
     * @param i Index of central atom (out-of-plane)
     * @param j Index of first plane atom
     * @param k Index of second plane atom
     * @param l Index of third plane atom
     * @return Out-of-plane angle in radians [-π/2, π/2]
     */
    double calculateOutOfPlaneAngle(int i, int j, int k, int l) const;

    /**
     * @brief Calculate derivatives of out-of-plane angle w.r.t. atomic positions
     *
     * Claude Generated (2025): Analytical gradient for inversions
     * Reference: external/gfnff/src/gfnff_helpers.f90:450-510 (domegadr subroutine)
     *
     * Computes ∂ω/∂x_i, ∂ω/∂x_j, ∂ω/∂x_k, ∂ω/∂x_l for chain rule in gradient.
     * Critical for analytical force calculation in geometry optimization.
     *
     * @param i Index of central atom
     * @param j Index of first plane atom
     * @param k Index of second plane atom
     * @param l Index of third plane atom
     * @param omega Current out-of-plane angle (from calculateOutOfPlaneAngle)
     * @param grad_i Output: ∂ω/∂x_i (3D vector)
     * @param grad_j Output: ∂ω/∂x_j (3D vector)
     * @param grad_k Output: ∂ω/∂x_k (3D vector)
     * @param grad_l Output: ∂ω/∂x_l (3D vector)
     */
    void calculateInversionGradient(int i, int j, int k, int l, double omega,
                                     Vector& grad_i, Vector& grad_j,
                                     Vector& grad_k, Vector& grad_l) const;

    /**
     * @brief Get GFN-FF inversion parameters for atom quartet
     *
     * Claude Generated (2025): Topology-aware inversion parameter assignment
     *
     * Assigns inversion parameters based on:
     * - Hybridization of central atom i (sp² → needs inversion)
     * - Element type (C, N, O, B different barriers)
     * - Pi-system membership (aromatics → higher barriers)
     *
     * @param z_i Atomic number of central atom (out-of-plane)
     * @param z_j Atomic number of plane atom j
     * @param z_k Atomic number of plane atom k
     * @param z_l Atomic number of plane atom l
     * @param hyb_i Hybridization of central atom i (1=sp, 2=sp², 3=sp³)
     * @return GFN-FF inversion parameters
     */
    GFNFFInversionParams getGFNFFInversionParameters(int z_i, int z_j, int z_k, int z_l,
                                                      int hyb_i) const;

    // =================================================================================
    // Advanced GFN-FF Parameter Generation (for future implementation)
    // =================================================================================

    /**
     * @brief Calculate coordination number derivatives for gradients
     * @param cn Coordination numbers
     * @param threshold Coordination number threshold (squared distance in Bohr²)
     * @return 3D tensor of CN derivatives (3 x natoms x natoms)
     */
    std::vector<Matrix> calculateCoordinationNumberDerivatives(const Vector& cn, double threshold = 1600.0) const;  // 40.0² = 1600 (squared)

    /**
     * @brief Determine hybridization states for all atoms
     * @return Vector of hybridization states (1=sp, 2=sp2, 3=sp3, 4=sp3d, 5=sp3d2)
     */
    std::vector<int> determineHybridization() const;

    /**
     * @brief Detect pi-systems and conjugated fragments
     * @param hyb Hybridization states
     * @return Vector mapping atoms to pi-fragment IDs (0 = no pi-system)
     */
    std::vector<int> detectPiSystems(const std::vector<int>& hyb) const;

    /**
     * @brief Find smallest ring size for each atom
     * @return Vector of smallest ring sizes (0 = not in ring)
     */
    std::vector<int> findSmallestRings() const;

    /**
     * @brief Check if two atoms are in the same ring
     * @param i First atom index
     * @param j Second atom index
     * @param ring_size Output: size of the ring they share (0 if not in same ring)
     * @return true if atoms are in the same ring
     */
    bool areAtomsInSameRing(int i, int j, int& ring_size) const;

    /**
     * @brief Calculate topology-dependent electronegativity corrections (dxi)
     *
     * Claude Generated (2025): Element-specific topology corrections for EEQ
     * Reference: gfnff_ini.f90:247-297 (dxicalc subroutine)
     *
     * Implements the missing dxi topology corrections that significantly
     * improve EEQ charge accuracy, especially for heteroatoms.
     *
     * dxi corrections include:
     * - Boron hydrogen-dependent corrections
     * - Carbene carbon corrections
     * - Oxygen nitro group and water corrections
     * - Group 6 (S, Se, etc.) overcoordination corrections
     * - Group 7 (Cl, Br, I) metal-dependent corrections
     *
     * @param atoms Atomic numbers
     * @param coordination_numbers Coordination numbers
     * @param neighbor_list Bond connectivity information
     * @param topo_info Topology information (pi-systems, etc.)
     * @return Vector of dxi corrections for each atom
     */
    Vector calculateDXI(const std::vector<int>& atoms,
                        const Vector& coordination_numbers,
                        const std::vector<std::vector<int>>& neighbor_list,
                        const TopologyInfo& topo_info) const;

    /**
     * @brief Calculate EEQ charges using extended electronegativity equalization
     * @param cn Coordination numbers
     * @param hyb Hybridization states
     * @param rings Ring information
     * @return Vector of EEQ charges
     */
    Vector calculateEEQCharges(const Vector& cn, const std::vector<int>& hyb, const std::vector<int>& rings) const;

    /**
     * @brief Calculate dgam (charge-dependent hardness) corrections
     *
     * Claude Generated (December 2025, Session 6): Extracted from calculateEEQCharges()
     * Reference: external/gfnff/src/gfnff_ini.f90:677-688
     *
     * Calculates charge-dependent gamma corrections that refine the EEQ hardness matrix
     * based on computed atomic charges and element type.
     *
     * @param qa_charges Base EEQ charges (from Phase 3.3 of calculateEEQCharges)
     * @param hybridization Hybridization state per atom (1=sp, 2=sp2, 3=sp3)
     * @param ring_sizes Smallest ring size per atom (0 if not in ring)
     * @return dgam corrections: Delta-gamma values to add to hardness matrix diagonal
     */
    Vector calculateDgam(const Vector& qa_charges,
                        const std::vector<int>& hybridization,
                        const std::vector<int>& ring_sizes) const;

    /**
     * @brief Build per-atom neighbor lists from bond pairs
     *
     * Claude Generated (December 2025, Session 6): Two-phase EEQ support
     * Converts cached bond list into per-atom neighbor connectivity for
     * enhanced topology analysis and future dxi corrections.
     *
     * Creates symmetric neighbor lists: if atom i bonds to j, then both lists updated
     *
     * @return Vector of neighbor lists (one std::vector<int> per atom)
     */
    std::vector<std::vector<int>> buildNeighborLists() const;

    /**
     * @brief Count neighbors within 20 Bohr cutoff (nb20)
     *
     * Claude Generated (January 14, 2026) - Phase 2: Exact nb20 implementation
     * Port from gfnff_ini2.f90 neighbor list generation.
     *
     * Returns the number of atoms within 20 Bohr (≈10.58 Å) of the given atom.
     * This is used for bond fcn correction factors in GFN-FF.
     *
     * @param atom_index Index of atom to count neighbors for
     * @param distance_matrix N×N distance matrix in Bohr
     * @return Number of neighbors within 20 Bohr cutoff
     */
    int countNeighborsWithin20Bohr(int atom_index, const Eigen::MatrixXd& distance_matrix) const;

    /**
     * @brief Calculate simplified π-bond orders for all atom pairs
     *
     * Claude Generated (January 10, 2026) - Phase 2C: π-bond order approximation
     *
     * Claude Generated (January 14, 2026) - Updated for Phase 1: Full Hückel implementation
     *
     * Two modes available (controlled by m_use_full_huckel):
     *
     * **Full Hückel mode (default, m_use_full_huckel=true)**:
     * Uses iterative self-consistent Hückel method from gfnff_ini.f90:928-1062.
     * - P-dependent off-diagonal coupling (prevents over-delocalization)
     * - Charge-dependent diagonal elements
     * - Fermi smearing at 4000K for biradical handling
     * - Exact π-bond orders from density matrix elements
     *
     * **Simplified mode (m_use_full_huckel=false)**:
     * Approximation based on hybridization and bond types:
     * - Single bonds (sp3-sp3): pbo = 0.0
     * - π bonds (sp2-sp2, sp2-sp): pbo = 0.5-1.0
     * - sp bonds (sp-sp): pbo = 1.5
     *
     * @param bond_list Vector of bonded atom pairs
     * @param hybridization Hybridization state per atom
     * @param pi_fragments Pi-system fragment IDs
     * @param charges EEQ atomic charges (needed for full Hückel)
     * @param distances N×N distance matrix in Bohr (needed for full Hückel)
     * @return Vector of π-bond orders in triangular format (access via lin(i,j))
     */
    std::vector<double> calculatePiBondOrders(
        const std::vector<std::pair<int,int>>& bond_list,
        const std::vector<int>& hybridization,
        const std::vector<int>& pi_fragments,
        const std::vector<double>& charges = {},
        const Eigen::MatrixXd& distances = Eigen::MatrixXd()) const;

    /**
     * @brief Calculate EEQ electrostatic energy
     *
     * Claude Generated (2025): Phase 3.2 EEQ energy contribution
     * Reference: gfnff_engrad.F90:1378-1389 (EEQ energy formula)
     *
     * Computes electrostatic energy from EEQ charges:
     * E_EEQ = Σ_ij q_i*q_j*γ_ij(r_ij) + Σ_i q_i*[-χ_i - cnf*√CN_i + 0.5*q_i*(γ_i + √(2π)/√α_i)]
     *
     * Note: Gradients use "frozen charge" approximation (∂q/∂r neglected)
     *
     * @param charges EEQ atomic charges from calculateEEQCharges()
     * @param cn Coordination numbers
     * @return EEQ electrostatic energy in Hartree
     */
    double calculateEEQEnergy(const Vector& charges, const Vector& cn) const;

    /**
     * @brief Generate topology-aware bond parameters
     * @param cn Coordination numbers
     * @param hyb Hybridization states
     * @param charges EEQ charges
     * @param rings Ring information
     * @return JSON with advanced bond parameters
     */
    json generateTopologyAwareBonds(const Vector& cn, const std::vector<int>& hyb,
        const Vector& charges, const std::vector<int>& rings) const;

    /**
     * Claude Generated (Jan 15, 2026): New overload accepting full TopologyInfo with pi_bond_orders
     * @param topo_info Complete topology information including pi_bond_orders
     * @return JSON with advanced bond parameters
     */
    json generateTopologyAwareBonds(const TopologyInfo& topo_info) const;

    /**
     * @brief Generate topology-aware angle parameters
     * @param cn Coordination numbers
     * @param hyb Hybridization states
     * @param charges EEQ charges
     * @param rings Ring information
     * @return JSON with advanced angle parameters
     */
    json generateTopologyAwareAngles(const Vector& cn, const std::vector<int>& hyb,
        const Vector& charges, const std::vector<int>& rings) const;

    /**
     * @brief Detect hydrogen bonds and set up A-H...B interactions
     * @param charges EEQ charges
     * @return JSON with hydrogen bond parameters
     */
    json detectHydrogenBonds(const Vector& charges) const;

    /**
     * @brief Detect halogen bond (XB) interactions
     * @param charges EEQ charges for charge-based criteria
     * @return JSON with halogen bond parameters
     * Claude Generated (2025): Phase 2.2 - XB Detection
     */
    json detectHalogenBonds(const Vector& charges) const;

    // Advanced parameter structures (EEQParameters already defined above at line 298)
    // TopologyInfo now defined at line 51 (public section) for use in function signatures

    /**
     * @brief Calculate full topology information for advanced parametrization
     * @return Complete topology information
     */
    TopologyInfo calculateTopologyInfo() const;

    /**
     * @brief Get EEQ parameters for specific atom with environment corrections
     * @param atom_idx Atom index
     * @param topo_info Topology information
     * @return EEQ parameters for this atom
     */
    EEQParameters getEEQParameters(int atom_idx, const TopologyInfo& topo_info) const;

    // =================================================================================
    // Two-Phase EEQ System Methods (Session 5, December 2025)
    // =================================================================================

    /**
     * @brief Phase 1: Calculate topology charges using base EEQ parameters
     *
     * Solves EEQ using ONLY base parameters without any corrections:
     * - chi = -chi_base (NO dxi corrections)
     * - gamma = gam_base (NO dgam corrections)
     * - alpha = alp_base^2 (NO dalpha corrections)
     *
     * These topology charges (qa) are then used in Phase 2 to calculate
     * the correction terms (dxi, dgam, dalpha).
     *
     * @param cn Coordination numbers
     * @param hyb Hybridization states
     * @param rings Ring information
     * @return Topology charges (qa) from base EEQ
     *
     * Reference: Fortran gfnff_ini.f90:405-421
     */
    Vector calculateTopologyCharges(const Vector& cn, const std::vector<int>& hyb,
                                     const std::vector<int>& rings) const;

    /**
     * @brief Calculate dxi (electronegativity) corrections
     *
     * Applies the full cascade logic from Fortran gfnff_ini.f90:361-403
     * with 30+ lines of element-specific, group-specific, and neighbor-dependent corrections.
     *
     * Corrections include:
     * - Boron: +nh*0.015 (hydrogen neighbors)
     * - Carbon: carbene (-0.15), free CO (+0.15)
     * - Oxygen: nitro (+0.05), water (-0.02), overcoordination (+nn*0.005)
     * - Group 6: overcoordination correction
     * - Group 7: polyvalent halogen corrections
     *
     * @param topo Topology information with neighbor lists and functional groups
     * @param qa_charges Topology charges from Phase 1
     * @return Electronegativity corrections per atom
     *
     * Reference: Fortran gfnff_ini.f90:361-403
     */
    Vector calculateDxi(const TopologyInfo& topo, const Vector& qa_charges) const;

    /**
     * @brief Calculate dalpha (polarizability) corrections
     *
     * Applies charge-dependent polarizability corrections:
     * alpeeq = (alp_base + ff * qa)^2
     *
     * where ff depends on element and group:
     * - C: +0.09, N: -0.21
     * - Group 6: -0.03, Group 7: +0.50
     * - Main-group metals: +0.3, Transition metals: -0.1
     *
     * @param qa_charges Topology charges from Phase 1
     * @return Polarizability corrections per atom (NOT yet squared)
     *
     * Reference: Fortran gfnff_ini.f90:694-707
     */
    Vector calculateDalpha(const Vector& qa_charges) const;

    /**
     * @brief Phase 2: Calculate final charges with all corrections applied
     *
     * Solves second EEQ with corrected parameters:
     * - chi = -chi_base + dxi [+ amide correction]
     * - gamma = gam_base + dgam
     * - alpha = (alp_base + dalpha)^2
     *
     * Also applies special amide hydrogen correction: chi -= 0.02
     *
     * @param topo Topology information with all corrections calculated
     * @return Final EEQ charges (q) to be used for force field calculations
     *
     * Reference: Fortran gfnff_ini.f90:694-707
     */
    Vector calculateFinalCharges(const TopologyInfo& topo) const;

public:
    // =================================================================================
    // Energy Component Access (Claude Generated November 2025)
    // =================================================================================

    /**
     * @brief Get bond energy component
     * @return Bond stretching energy or 0 if not calculated
     */
    double BondEnergy() const;

    /**
     * @brief Get angle energy component
     * @return Angle bending energy or 0 if not calculated
     */
    double AngleEnergy() const;

    /**
     * @brief Get dihedral energy component
     * @return Dihedral torsion energy or 0 if not calculated
     */
    double DihedralEnergy() const;

    /**
     * @brief Get inversion energy component
     * @return Inversion/out-of-plane energy or 0 if not calculated
     */
    double InversionEnergy() const;

    /**
     * @brief Get van der Waals energy component
     * @return Van der Waals interaction energy or 0 if not calculated
     */
    double VdWEnergy() const;

    /**
     * @brief Get repulsion energy component
     * @return Core-core repulsion energy or 0 if not calculated
     */
    double RepulsionEnergy() const;

    /**
     * @brief Get dispersion energy component
     * @return Dispersion correction energy or 0 if not calculated
     */
    double DispersionEnergy() const;

    /**
     * @brief Get Coulomb electrostatic energy component
     * @return Electrostatic energy or 0 if not calculated
     */
    double CoulombEnergy() const;

    /**
     * @brief Get D3 dispersion energy component
     * @return D3 dispersion energy or 0 if not calculated
     * Claude Generated (Jan 2, 2026): D3 dispersion energy accessor
     */
    double D3Energy() const;

    /**
     * @brief Get D4 dispersion energy component
     * @return D4 dispersion energy or 0 if not calculated
     * Claude Generated (Jan 2, 2026): D4 dispersion energy accessor
     */
    double D4Energy() const;

    /**
     * @brief Get batm (bonded ATM) energy component
     * @return Batm energy or 0 if not calculated
     * Claude Generated (Jan 17, 2026): Batm energy accessor for 1,4-pairs
     */
    double BatmEnergy() const;

    // =================================================================================
    // vbond Parameter Access for Verification (Claude Generated November 2025)
    // =================================================================================

    /**
     * @brief Get vbond parameters for verification against reference implementation
     * @param bond_index Index of bond (0-based)
     * @param shift Output: vbond(1,i) - equilibrium distance shift parameter
     * @param alpha Output: vbond(2,i) - exponential decay parameter
     * @param force_constant Output: vbond(3,i) - force constant parameter
     * @return true if parameters were successfully retrieved
     */
    bool getVBondParameters(int bond_index, double& shift, double& alpha, double& force_constant) const;

    /**
     * @brief Get number of bonds in the system
     * @return Number of bonds
     */
    int getBondCount() const;

    // =================================================================================
    // TWO-PHASE EEQ SYSTEM (Claude Generated November 2025, Session 5)
    // =================================================================================

    /**
     * @brief Phase 1: Calculate topology-aware base charges (qa) via EEQ
     *
     * The two-phase EEQ system separates charge calculation into:
     * 1. TOPOLOGY PHASE: Base charges from atomic properties (electronegativity, hardness)
     *    using coordination-dependent parameters
     * 2. CORRECTION PHASE: Apply dxi, dgam, dalpha corrections for refined accuracy
     *
     * Reference: angewChem 2020, GFN-FF parameter set
     *   - chi[z]: Electronegativity for element z
     *   - gam[z]: Chemical hardness for element z
     *   - alp[z]: Damping parameter for erf(γ*r)/r Coulomb
     *   - cnf[z]: Coordination number correction factor
     *
     * Extended Hückel Theory (EHT) approximation:
     *   qa_i = -χ_i - J_ii + Σ_j(1/(2*J_ij) - 1/(2*r_ij))
     *
     * @param topo_info TopologyInfo structure with coordination numbers, hybridization
     * @return true if Phase 1 charges calculated successfully
     *
     * Output written to: topo_info.topology_charges (qa in Hartree)
     */
    bool calculateTopologyCharges(TopologyInfo& topo_info) const;

    /**
     * @brief Calculate dxi (electronegativity) corrections for Phase 2
     *
     * dxi corrects electronegativity based on:
     * - Local environment (neighbor count, hybridization)
     * - Bonding context (pi-systems, heteroatom effects)
     * - Functional group classification
     *
     * Physical meaning: Electronegativity is NOT constant - it depends on
     * chemical context. Atoms in electron-withdrawing groups become more
     * electronegative.
     *
     * @param topo_info TopologyInfo with topology_charges and hybrid classifications
     * @return true if dxi corrections calculated
     *
     * Output written to: topo_info.dxi (corrections to chi in Hartree)
     */
    bool calculateDxi(TopologyInfo& topo_info) const;

    /**
     * @brief Calculate dalpha (polarizability) corrections for Phase 2
     *
     * dalpha corrects the damping parameter (alpha) based on:
     * - Atomic size changes (coordination-dependent)
     * - Electronic environment (hybridization, charge state)
     * - Pi-system participation
     *
     * Physical meaning: Polarizability (and hence Coulomb damping) adapts to
     * local electronic density. More polarizable atoms in electron-rich
     * environments use different damping.
     *
     * @param topo_info TopologyInfo with coordination numbers and charges
     * @return true if dalpha corrections calculated
     *
     * Output written to: topo_info.dalpha (corrections to alpha)
     */
    bool calculateDalpha(TopologyInfo& topo_info) const;

    /**
     * @brief Calculate charge-dependent alpha (alpeeq) for EEQ
     *
     * Claude Generated (January 2026)
     *
     * Computes charge-dependent alpha values used in EEQ matrix construction.
     * This implements the formula from Fortran gfnff_ini.f90:718-725:
     *
     *   alpeeq(i) = (alpha_base + ff*qa(i))²
     *
     * where ff is element-specific:
     *   - Carbon (Z=6): ff = 0.09
     *   - Nitrogen (Z=7): ff = -0.21
     *   - Group 6 (O,S,Se): ff = -0.03
     *   - Group 7 (Halogens): ff = 0.50
     *   - Main group metals: ff = 0.3
     *   - Transition metals: ff = -0.1
     *
     * Physical meaning: Charge state modifies atomic polarizability (Gaussian width).
     * Positive charges increase alpha (softer, more diffuse) for C/main-group metals,
     * negative charges increase alpha for N/halogens/transition metals.
     *
     * CRITICAL: This must be called AFTER topology_charges are computed,
     * and the resulting alpeeq values are used UNCHANGED in all subsequent
     * EEQ calculations (no iteration).
     *
     * Reference: Fortran gfnff_data_types.f90:128 - "atomic alpha for EEQ, squared"
     *
     * @param topo_info TopologyInfo with topology_charges (Phase 1 charges qa)
     * @return true if alpeeq calculated successfully
     *
     * Output written to: topo_info.alpeeq (squared alpha values)
     */
    bool calculateAlpeeq(TopologyInfo& topo_info) const;

    /**
     * @brief Get generated ForceField parameters
     *
     * Returns the full parameter set generated by the ForceField engine,
     * including bonds, angles, dihedrals, electrostatics, dispersion, etc.
     *
     * This is different from getParameters() which only returns the input JSON.
     *
     * @return JSON object containing all force field parameters
     * @return Empty JSON if ForceField is not initialized
     *
     * Claude Generated - December 27, 2025
     */
    json getForceFieldParameters() const {
        if (m_forcefield) {
            return m_forcefield->exportCurrentParameters();
        }
        return json();
    };

    /**
     * @brief Phase 2: Calculate final refined charges by solving corrected EEQ
     *
     * Iteratively solves EEQ with corrections:
     *   qa_i_final = qa_i + dxi_i + (dgam correction)
     *   Then re-solve EEQ with modified parameters
     *
     * The correction application is:
     * 1. Modify electronegativity: χ'_i = χ_i + dxi_i
     * 2. Recalculate Coulomb matrix with corrected polarizabilities: α'_i = α_i + dalpha_i
     * 3. Re-solve the linear EEQ system to consistency
     *
     * Convergence: Typically 1-2 iterations for <0.01 e change per atom
     *
     * @param topo_info TopologyInfo with topology_charges, dxi, dalpha corrections
     * @param max_iterations Maximum iterations for EEQ convergence (default: 10)
     * @param convergence_threshold Threshold for charge change (default: 1e-5 Hartree)
     * @return true if Phase 2 refinement successful
     *
     * Output written to: topo_info.eeq_charges (final qa in Hartree)
     */
    bool calculateFinalCharges(TopologyInfo& topo_info, int max_iterations = 10,
                               double convergence_threshold = 1e-5) const;

    /**
     * @brief Set external charges (for testing/validation)
     * @param charges Atomic partial charges to use instead of EEQ-calculated charges
     *
     * Claude Generated (December 2025): Testing utility for charge-dependent validation
     * Bypasses EEQ calculation and uses provided reference charges directly.
     * Call AFTER InitialiseMolecule() but BEFORE Calculation().
     * This method is intended for validation purposes to isolate energy calculation
     * errors from EEQ charge calculation errors.
     */
    void setCharges(const Vector& charges);

    /**
     * @brief Regenerate GFN-FF parameters using current charges (for testing/validation)
     *
     * Claude Generated (January 2025): Testing utility for charge-dependent validation
     * Regenerates all charge-dependent parameters (bonds, angles, dihedrals, inversions)
     * using the current m_charges (set via setCharges()).
     *
     * Call AFTER setCharges() but BEFORE Calculation() to ensure parameter-charge consistency.
     *
     * This is required because GFN-FF parameters depend on charges:
     * - Bond force constants depend on fqq (charge-dependent correction)
     * - Angle force constants have charge-dependent terms
     * - Electrostatic parameters depend on charge distribution
     *
     * @return true if parameter regeneration successful
     */
    bool regenerateParametersWithCurrentCharges();

    /**
     * @brief Get current bond parameters (for validation)
     *
     * Claude Generated (December 2025): Phase 3 - Parameter validation infrastructure
     * Returns the current bond parameters as JSON for systematic validation.
     * Used by Test 5 to compare generated parameters against XTB reference.
     *
     * @return JSON array of bond parameters with fields: i, j, distance, k, fc, etc.
     */
    json getBondParameters() const;

    /**
     * @brief Get current angle parameters (for validation)
     *
     * Claude Generated (December 2025): Phase 3 - Parameter validation infrastructure
     *
     * @return JSON array of angle parameters
     */
    json getAngleParameters() const;

    /**
     * @brief Set bond parameters for testing (bypasses generation)
     *
     * Claude Generated (December 2025): Phase 3 - Parameter injection for Test 8
     * Allows injection of exact reference parameters to isolate energy calculation
     * errors from parameter generation errors.
     *
     * @param bond_params JSON array of bond parameters in ForceField format
     */
    void setBondParametersForTesting(const json& bond_params);

    /**
     * @brief Set angle parameters for testing (bypasses generation)
     *
     * Claude Generated (December 2025): Phase 3 - Parameter injection for Test 8
     *
     * @param angle_params JSON array of angle parameters in ForceField format
     */
    void setAngleParametersForTesting(const json& angle_params);

private:
    // Molecular structure (formerly from QMInterface base class)
    int m_atomcount = 0; ///< Number of atoms
    Matrix m_geometry; ///< Molecular geometry in Angström
    Matrix m_gradient; ///< Gradient in Hartree/Bohr
    std::vector<int> m_atoms; ///< Atomic numbers (Z values)
    int m_charge = 0; ///< Total molecular charge
    int m_spin = 0; ///< Spin multiplicity

    // GFN-FF specific
    json m_parameters; ///< GFN-FF parameters
    ForceField* m_forcefield; ///< Force field engine using modern structure
    Matrix m_geometry_bohr; ///< Geometry in Bohr (GFN-FF parameters are in Bohr)

    // EEQ charge calculation (Dec 2025 - Phase 3: Extraction and delegation)
    std::unique_ptr<EEQSolver> m_eeq_solver; ///< Standalone EEQ solver (replaces embedded EEQ code)

    // Hückel solver for π-bond orders (Jan 2026 - Phase 1: Full Hückel implementation)
    std::unique_ptr<HuckelSolver> m_huckel_solver; ///< Full iterative Hückel solver
    bool m_use_full_huckel = true; ///< Use full Hückel calculation (default: true, set to false for simplified approximation)

    // ATM three-body dispersion terms (extracted from D3/D4 - Claude Generated Jan 2025)
    mutable json m_atm_triples; ///< Bonded ATM triples from D3/D4 parameter generators

    // Geometry change detection for intelligent caching
    class GeometryChangeDetector {
    private:
        Matrix m_last_geometry;
        double m_change_threshold = 1e-6;

    public:
        bool geometryChanged(const Matrix& new_geometry) const {
            if (m_last_geometry.rows() != new_geometry.rows() ||
                m_last_geometry.cols() != new_geometry.cols()) {
                return true;
            }
            // Only invalidate cache if change exceeds threshold
            return (m_last_geometry - new_geometry).array().abs().maxCoeff() > m_change_threshold;
        }

        void updateGeometry(const Matrix& new_geometry) {
            m_last_geometry = new_geometry;
        }

        void reset() {
            m_last_geometry = Matrix();
        }
    };

    bool m_initialized; ///< Initialization status

    double m_energy_total; ///< Total energy in Hartree
    Vector m_charges; ///< Atomic partial charges
    Vector m_bond_orders; ///< Wiberg bond orders

    // Geometry change detector for intelligent caching
    mutable GeometryChangeDetector m_geometry_tracker;

    // Cached topology and bond detection to avoid redundant calculations
    mutable std::optional<TopologyInfo> m_cached_topology;
    mutable std::optional<std::vector<std::pair<int,int>>> m_cached_bond_list;

    // Conversion factors
    static constexpr double HARTREE_TO_KCAL = 627.5094740631;
    static constexpr double BOHR_TO_ANGSTROM = 0.5291772105638411;
    static constexpr double KCAL_TO_HARTREE = 1.0 / 627.5094740631;
    static constexpr double ANGSTROM_TO_BOHR = 1.0 / 0.5291772105638411;

    /**
     * @brief Element-specific radius scaling factors (fat array from gfnff_ini2.f90:76-97)
     *
     * Claude Generated (January 2026) - Phase 3: Element-specific neighbor detection
     * Applied in bond detection: threshold = 1.3 * (rcov_i + rcov_j) * fat[Z_i] * fat[Z_j]
     *
     * Default: 1.0 for all elements
     * Special adjustments for specific elements to improve bond detection accuracy
     */
    static constexpr double fat[87] = {
        0.0,   // 0: placeholder (atom numbers start at 1)
        1.02,  // 1: H
        1.00,  // 2: He
        1.00,  // 3: Li
        1.03,  // 4: Be
        1.02,  // 5: B
        1.00,  // 6: C
        1.00,  // 7: N
        1.02,  // 8: O
        1.05,  // 9: F
        1.10,  // 10: Ne
        1.01,  // 11: Na
        1.02,  // 12: Mg
        1.00,  // 13: Al
        1.00,  // 14: Si
        0.97,  // 15: P
        1.00,  // 16: S
        1.00,  // 17: Cl
        1.10,  // 18: Ar
        1.02,  // 19: K
        1.02,  // 20: Ca
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  // 21-30: Sc-Zn
        1.00, 1.00, 1.00, 0.99,  // 31-34: Ga-Se (34: Se = 0.99)
        1.00, 1.00, 1.00,        // 35-37: Br-Rb
        1.02,  // 38: Sr
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  // 39-49: Y-In
        1.01,  // 50: Sn
        0.99,  // 51: Sb
        0.95,  // 52: Te
        0.98,  // 53: I
        1.00, 1.00,  // 54-55: Xe-Cs
        1.02,  // 56: Ba
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  // 57-66: La-Dy
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,        // 67-75: Ho-Re
        1.02,  // 76: Os
        1.00, 1.00, 1.00, 1.00, 1.00,  // 77-81: Ir-Tl
        1.06,  // 82: Pb
        0.95,  // 83: Bi
        1.00, 1.00, 1.00   // 84-86: Po-Rn
    };
};

/**
 * @brief GFN-FF Architecture: Two-Phase Implementation
 *
 * PARAMETER GENERATION (this class):
 * - generateTopologyAwareBonds()     → JSON["bonds"]
 * - generateTopologyAwareAngles()    → JSON["angles"]
 * - generateGFNFFTorsions()          → JSON["dihedrals"]
 * - generateGFNFFInversions()        → JSON["inversions"]
 * - generateGFNFFDispersionPairs()   → JSON["gfnff_dispersions"]
 * - generateGFNFFRepulsionPairs()    → JSON["gfnff_repulsions"]
 * - generateGFNFFCoulombPairs()      → JSON["gfnff_coulombs"]
 *
 * TERM CALCULATION (ForceFieldThread):
 * - CalculateGFNFFBondContribution()
 * - CalculateGFNFFAngleContribution()
 * - CalculateGFNFFDihedralContribution()
 * - CalculateGFNFFInversionContribution()
 * - CalculateGFNFFDispersionContribution()
 * - CalculateGFNFFRepulsionContribution()
 * - CalculateGFNFFCoulombContribution()
 *
 * To add new terms: Modify BOTH this class (generation) AND ForceFieldThread (calculation)
 * See: src/core/energy_calculators/ff_methods/CLAUDE.md for detailed checklist
 */
