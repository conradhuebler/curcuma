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

#include "src/core/global.h"
#include "src/core/qm_methods/interface/abstract_interface.h"
#include "src/core/forcefield.h"
#include "json.hpp"

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
    json generateGFNFFBonds();

    /**
     * @brief Generate GFN-FF angle parameters from topology
     * @return JSON array of angle parameters
     */
    json generateGFNFFAngles();

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
     * @brief Get GFN-FF angle parameters for element triplet
     * @param z1 Atomic number of first atom
     * @param z2 Atomic number of center atom
     * @param z3 Atomic number of third atom
     * @param current_angle Current angle in radians
     * @return GFN-FF angle parameters
     */
    GFNFFAngleParams getGFNFFAngleParameters(int z1, int z2, int z3, double current_angle) const;

private:
    json m_parameters;                    ///< GFN-FF parameters
    ForceField* m_forcefield;            ///< Force field engine using modern structure
    
    bool m_initialized;                  ///< Initialization status
    
    double m_energy_total;               ///< Total energy in Hartree
    Vector m_charges;                    ///< Atomic partial charges
    Vector m_bond_orders;                ///< Wiberg bond orders
    
    // Conversion factors
    static constexpr double HARTREE_TO_KCAL = 627.5094740631;
    static constexpr double BOHR_TO_ANGSTROM = 0.5291772105638411;
    static constexpr double KCAL_TO_HARTREE = 1.0 / 627.5094740631;
    static constexpr double ANGSTROM_TO_BOHR = 1.0 / 0.5291772105638411;
};