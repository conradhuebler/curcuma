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
     * @param verbose Enable detailed output
     * @return Total energy in Hartree
     */
    virtual double Calculation(bool gradient = false, bool verbose = false) override;

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
     * @brief Get covalent radius for element
     * @param atomic_number Element atomic number
     * @return Covalent radius in Angstrom
     */
    double getCovalentRadius(int atomic_number) const;

    // GFN-FF parameter structures
    struct GFNFFBondParams {
        double force_constant;
        double equilibrium_distance;
        double anharmonic_factor;
    };

    struct GFNFFAngleParams {
        double force_constant;
        double equilibrium_angle;
        double c0, c1, c2; // Fourier coefficients
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
     * @brief Get GFN-FF angle parameters for element triplet
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of center atom
     * @param z3 Atomic number of third atom
     * @param current_angle Current angle in radians
     * @return GFN-FF angle parameters
     */
    GFNFFAngleParams getGFNFFAngleParameters(int z1, int z2, int z3, double current_angle) const;

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