/*
 * < Unified Computational Method Interface >
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
 */

#pragma once

#include "src/core/global.h"
#include "json.hpp"

using json = nlohmann::json;

/**
 * @brief Unified interface for all computational methods (QM and MM)
 * 
 * This interface provides a common API for quantum mechanical and molecular
 * mechanics methods, enabling polymorphic dispatch and consistent handling
 * in the EnergyCalculator.
 * 
 * Key design principles:
 * - Polymorphic design using virtual functions for flexibility
 * - Thread-safety support for parallel calculations
 * - Unified property access (charges, gradients, bond orders)
 * - Consistent parameter handling via JSON
 * - Error handling and status reporting
 * 
 * Claude Generated: Big-Bang refactoring of EnergyCalculator architecture
 */
class ComputationalMethod {
public:
    virtual ~ComputationalMethod() = default;
    
    // =================================================================================
    // Core Molecular Setup and Calculation
    // =================================================================================
    
    /**
     * @brief Initialize method with molecular structure
     * @param mol Molecule object containing geometry, atoms, charge, spin
     * @return true if initialization successful
     */
    virtual bool setMolecule(const Mol& mol) = 0;
    
    /**
     * @brief Update molecular geometry without full reinitialization
     * @param geometry New atomic coordinates (natoms x 3 matrix)
     * @return true if update successful
     */
    virtual bool updateGeometry(const Matrix& geometry) = 0;

    /**
     * @brief Perform energy and gradient calculation
     * @param gradient Calculate analytical gradients if available and requested
     * @return Total energy in appropriate units (typically Hartree or kcal/mol)
     */
    virtual double calculateEnergy(bool gradient = false) = 0;

    // =================================================================================
    // Property Access
    // =================================================================================
    
    /**
     * @brief Get analytical or numerical gradients
     * @return Gradient matrix (natoms x 3) in appropriate units
     */
    virtual Matrix getGradient() const = 0;
    
    /**
     * @brief Get atomic partial charges
     * @return Vector of atomic charges (length = natoms)
     */
    virtual Vector getCharges() const = 0;
    
    /**
     * @brief Get bond orders between atoms
     * @return Vector or matrix of bond orders (method-dependent format)
     */
    virtual Vector getBondOrders() const = 0;
    
    /**
     * @brief Get molecular dipole moment
     * @return 3D dipole vector in appropriate units
     */
    virtual Position getDipole() const = 0;
    
    /**
     * @brief Check if analytical gradients are available
     * @return true if method can provide analytical gradients
     */
    virtual bool hasGradient() const = 0;
    
    // =================================================================================
    // Method Information and Configuration
    // =================================================================================
    
    /**
     * @brief Get method name for identification
     * @return String identifier (e.g., "eht", "gfn2-xtb", "uff", "cgfnff")
     */
    virtual std::string getMethodName() const = 0;
    
    /**
     * @brief Check if method supports concurrent calculations
     * @return true if method is thread-safe for parallel usage
     */
    virtual bool isThreadSafe() const = 0;
    
    /**
     * @brief Configure number of threads for parallel calculations
     * @param threads Number of threads to use (1 = sequential)
     */
    virtual void setThreadCount(int threads) = 0;
    
    /**
     * @brief Set method-specific parameters
     * @param params JSON configuration with method parameters
     */
    virtual void setParameters(const json& params) = 0;
    
    /**
     * @brief Get current method parameters
     * @return JSON object with current configuration
     */
    virtual json getParameters() const = 0;
    
    // =================================================================================
    // Error Handling and Status
    // =================================================================================
    
    /**
     * @brief Check if method is in error state
     * @return true if last operation failed
     */
    virtual bool hasError() const = 0;
    
    /**
     * @brief Clear error state and reset method
     */
    virtual void clearError() { /* Default: no-op */ }
    
    /**
     * @brief Get human-readable error message
     * @return Error description or empty string if no error
     */
    virtual std::string getErrorMessage() const { return ""; }
    
    // =================================================================================
    // Optional Advanced Features
    // =================================================================================
    
    /**
     * @brief Get orbital energies (QM methods only)
     * @return Vector of orbital energies or empty if not available
     */
    virtual Vector getOrbitalEnergies() const { return Vector{}; }
    
    /**
     * @brief Get orbital occupations (QM methods only)
     * @return Vector of orbital occupations or empty if not available
     */
    virtual Vector getOrbitalOccupations() const { return Vector{}; }
    
    /**
     * @brief Get number of electrons (QM methods only)
     * @return Total number of electrons or 0 if not applicable
     */
    virtual int getNumElectrons() const { return 0; }

    // =================================================================================
    // Force Field Energy Component Access (Claude Generated November 2025)
    // =================================================================================

    /**
     * @brief Get bond energy component (Force field methods only)
     * @return Bond stretching energy or 0 if not applicable
     */
    virtual double getBondEnergy() const { return 0.0; }

    /**
     * @brief Get angle energy component (Force field methods only)
     * @return Angle bending energy or 0 if not applicable
     */
    virtual double getAngleEnergy() const { return 0.0; }

    /**
     * @brief Get dihedral energy component (Force field methods only)
     * @return Dihedral torsion energy or 0 if not applicable
     */
    virtual double getDihedralEnergy() const { return 0.0; }

    /**
     * @brief Get inversion energy component (Force field methods only)
     * @return Inversion/out-of-plane energy or 0 if not applicable
     */
    virtual double getInversionEnergy() const { return 0.0; }

    /**
     * @brief Get van der Waals energy component (Force field methods only)
     * @return Van der Waals interaction energy or 0 if not applicable
     */
    virtual double getVdWEnergy() const { return 0.0; }

    /**
     * @brief Get repulsion energy component (Force field methods only)
     * @return Core-core repulsion energy or 0 if not applicable
     */
    virtual double getRepulsionEnergy() const { return 0.0; }

    /**
     * @brief Get dispersion energy component (Force field methods only)
     * @return Dispersion correction energy or 0 if not applicable
     */
    virtual double getDispersionEnergy() const { return 0.0; }

    /**
     * @brief Get Coulomb electrostatic energy component (Force field methods only)
     * @return Electrostatic energy or 0 if not applicable
     */
    virtual double getCoulombEnergy() const { return 0.0; }

    /**
     * @brief Save method-specific data to file
     * @param filename Output file name
     * @return true if save successful
     */
    virtual bool saveToFile(const std::string& filename) const { return false; }
    
    /**
     * @brief Load method-specific data from file
     * @param filename Input file name
     * @return true if load successful
     */
    virtual bool loadFromFile(const std::string& filename) { return false; }
    
protected:
    // Common data members that derived classes might use
    bool m_initialized = false;           ///< Initialization status
    bool m_has_error = false;            ///< Error state flag
    std::string m_error_message;         ///< Error description
    int m_thread_count = 1;              ///< Number of threads to use
    json m_parameters;                   ///< Method parameters
};

/**
 * @brief Utility functions for ComputationalMethod implementations
 */
namespace ComputationalMethodUtils {
    
    /**
     * @brief Convert energy units (implementation-dependent)
     * @param energy Input energy
     * @param from_unit Source unit ("hartree", "kcal", "ev")
     * @param to_unit Target unit
     * @return Converted energy
     */
    double convertEnergy(double energy, const std::string& from_unit, const std::string& to_unit);
    
    /**
     * @brief Convert gradient units
     * @param gradient Input gradient matrix
     * @param from_unit Source unit ("hartree/bohr", "kcal/angstrom")
     * @param to_unit Target unit
     * @return Converted gradient matrix
     */
    Matrix convertGradient(const Matrix& gradient, const std::string& from_unit, const std::string& to_unit);
    
    /**
     * @brief Validate molecular structure for computational methods
     * @param mol Molecule to validate
     * @return true if molecule is suitable for calculations
     */
    bool validateMolecule(const Mol& mol);
}