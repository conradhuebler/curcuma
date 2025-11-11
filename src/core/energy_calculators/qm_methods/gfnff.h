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

#include "interface/abstract_interface.h"
#include "json.hpp"
#include "src/core/energy_calculators/ff_methods/forcefield.h"
#include "src/core/global.h"

using json = nlohmann::json;

/**
 * @brief GFN-FF Implementation as QM Method
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
 */
class GFNFF : public QMInterface {
public:
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
     * @brief Initialize molecule for GFN-FF calculation
     * @return true if initialization successful
     */
    virtual bool InitialiseMolecule() override;

    /**
     * @brief Update molecular geometry
     * @return true if update successful
     */
    virtual bool UpdateMolecule() override;

    /**
     * @brief Perform GFN-FF calculation
     * @param gradient Calculate gradients if true
     * @return Total energy in Hartree
     */
    virtual double Calculation(bool gradient = false) override;

    /**
     * @brief Get analytical gradients
     * @return Gradient matrix (N_atoms x 3) in Hartree/Bohr
     */
    virtual Geometry Gradient() const override { return m_gradient; }

    /**
     * @brief Check if gradients are available
     * @return true (GFN-FF always provides gradients)
     */
    virtual bool hasGradient() const override { return true; }

    /**
     * @brief Get atomic partial charges
     * @return Vector of atomic charges
     */
    virtual Vector Charges() const override;

    /**
     * @brief Get bond orders (Wiberg bond orders)
     * @return Vector of bond orders
     */
    virtual Vector BondOrders() const override;

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
     * @return JSON array of angle parameters
     */
    json generateGFNFFAngles() const;

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
    };

    /**
     * @brief Get GFN-FF bond parameters for element pair
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of second atom
     * @param distance Current bond distance
     * @return GFN-FF bond parameters
     */
    GFNFFBondParams getGFNFFBondParameters(int z1, int z2, double distance) const;

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
     * @brief Get GFN-FF angle parameters for element triplet
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of center atom
     * @param z3 Atomic number of third atom
     * @param current_angle Current angle in radians
     * @return GFN-FF angle parameters
     */
    GFNFFAngleParams getGFNFFAngleParameters(int z1, int z2, int z3, double current_angle) const;

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
                                                  int hyb_j, int hyb_k) const;

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
     * @brief Calculate coordination numbers for all atoms
     * @param threshold Coordination number threshold
     * @return Vector of coordination numbers
     */
    Vector calculateCoordinationNumbers(double threshold = 40.0) const;

    /**
     * @brief Calculate coordination number derivatives for gradients
     * @param cn Coordination numbers
     * @param threshold Coordination number threshold
     * @return 3D tensor of CN derivatives (3 x natoms x natoms)
     */
    std::vector<Matrix> calculateCoordinationNumberDerivatives(const Vector& cn, double threshold = 40.0) const;

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
     * @brief Calculate EEQ charges using extended electronegativity equalization
     * @param cn Coordination numbers
     * @param hyb Hybridization states
     * @param rings Ring information
     * @return Vector of EEQ charges
     */
    Vector calculateEEQCharges(const Vector& cn, const std::vector<int>& hyb, const std::vector<int>& rings) const;

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

    // Advanced parameter structures
    struct EEQParameters {
        double chi; // Electronegativity
        double gam; // Chemical hardness
        double alp; // Polarizability
        double xi_corr; // Environment correction
    };

    struct TopologyInfo {
        Vector coordination_numbers;
        std::vector<int> hybridization;
        std::vector<int> pi_fragments;
        std::vector<int> ring_sizes;
        Vector eeq_charges;
        std::vector<bool> is_metal;
        std::vector<bool> is_aromatic;
    };

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

private:
    json m_parameters; ///< GFN-FF parameters
    ForceField* m_forcefield; ///< Force field engine using modern structure

    bool m_initialized; ///< Initialization status

    double m_energy_total; ///< Total energy in Hartree
    Vector m_charges; ///< Atomic partial charges
    Vector m_bond_orders; ///< Wiberg bond orders

    // Conversion factors
    static constexpr double HARTREE_TO_KCAL = 627.5094740631;
    static constexpr double BOHR_TO_ANGSTROM = 0.5291772105638411;
    static constexpr double KCAL_TO_HARTREE = 1.0 / 627.5094740631;
    static constexpr double ANGSTROM_TO_BOHR = 1.0 / 0.5291772105638411;
};