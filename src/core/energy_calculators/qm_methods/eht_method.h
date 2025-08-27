/*
 * < EHT Method Wrapper for ComputationalMethod Interface >
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

#include "../computational_method.h"
#include "eht.h"

#include <memory>

/**
 * @brief Extended Hückel Theory wrapper for ComputationalMethod interface
 * 
 * This wrapper adapts the existing EHT implementation to the unified
 * ComputationalMethod interface, providing consistent access to EHT
 * calculations within the new EnergyCalculator architecture.
 * 
 * Features:
 * - Wraps existing EHT semi-empirical quantum method
 * - Provides molecular orbital energies and coefficients
 * - Educational transparency: direct access to quantum chemistry implementation
 * - Thread-safe by design (EHT calculations are independent)
 * 
 * Claude Generated: Big-Bang EnergyCalculator refactoring wrapper
 */
class EHTMethod : public ComputationalMethod {
public:
    /**
     * @brief Constructor with configuration
     * @param config JSON configuration (EHT-specific parameters)
     */
    explicit EHTMethod(const json& config = json{});
    
    /**
     * @brief Destructor
     */
    virtual ~EHTMethod() = default;
    
    // =================================================================================
    // Core ComputationalMethod Interface Implementation
    // =================================================================================
    
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    // Property access
    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override;
    Position getDipole() const override;
    bool hasGradient() const override { return false; } // EHT doesn't provide analytical gradients
    
    // Method information
    std::string getMethodName() const override { return "eht"; }
    bool isThreadSafe() const override { return true; } // EHT calculations are independent
    void setThreadCount(int threads) override; // EHT can be parallelized
    
    // Configuration
    void setParameters(const json& params) override;
    json getParameters() const override;
    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;
    
    // =================================================================================
    // EHT-specific Methods (additional functionality)
    // =================================================================================
    
    /**
     * @brief Get molecular orbital energies
     * @return Vector of orbital energies in eV
     */
    Vector getOrbitalEnergies() const override;
    
    /**
     * @brief Get molecular orbital coefficients matrix
     * @return MO coefficient matrix (nbasis x norbitals)
     */
    Matrix getMolecularOrbitals() const;
    
    /**
     * @brief Get number of electrons
     * @return Total number of valence electrons
     */
    int getNumElectrons() const override;
    
    /**
     * @brief Get HOMO-LUMO gap
     * @return Energy gap between HOMO and LUMO in eV
     */
    double getHOMOLUMOGap() const;
    
    /**
     * @brief Get HOMO energy
     * @return Highest Occupied Molecular Orbital energy in eV
     */
    double getHOMOEnergy() const;
    
    /**
     * @brief Get LUMO energy
     * @return Lowest Unoccupied Molecular Orbital energy in eV
     */
    double getLUMOEnergy() const;
    
    /**
     * @brief Set Wolfsberg-Helmholz constant
     * @param K New K value (typically between 1.5 and 2.0)
     */
    void setWolfsbergHelmholzConstant(double K);
    
    /**
     * @brief Get current Wolfsberg-Helmholz constant
     * @return Current K value
     */
    double getWolfsbergHelmholzConstant() const;
    
    /**
     * @brief Print orbital analysis summary
     * @param num_orbitals_around_gap Number of orbitals to show around HOMO-LUMO gap
     */
    void printOrbitalAnalysis(int num_orbitals_around_gap = 5) const;
    
    /**
     * @brief Save orbital data to file
     * @param filename Output file name
     * @return true if save successful
     */
    bool saveToFile(const std::string& filename) const override;
    
private:
    std::unique_ptr<EHT> m_eht;                    ///< Wrapped EHT implementation
    Mol m_molecule;                                ///< Current molecule
    bool m_calculation_done;                       ///< Flag if calculation was performed
    double m_last_energy;                          ///< Last calculated energy
    
    /**
     * @brief Initialize EHT method with current parameters
     * @return true if initialization successful
     */
    bool initializeEHT();
    
    /**
     * @brief Update EHT parameters from JSON configuration
     */
    void updateEHTParameters();
    
    /**
     * @brief Convert EHT-specific errors to ComputationalMethod format
     */
    void handleEHTError(const std::string& operation);
    
    /**
     * @brief Default EHT configuration
     */
    static json getDefaultConfig();
};

/**
 * @brief Utility functions for EHT method
 */
namespace EHTMethodUtils {
    /**
     * @brief Validate EHT parameters
     * @param config JSON configuration to validate
     * @return true if parameters are valid
     */
    bool validateEHTConfig(const json& config);
    
    /**
     * @brief Get supported elements for EHT
     * @return Vector of atomic numbers supported by EHT
     */
    std::vector<int> getSupportedElements();
    
    /**
     * @brief Check if molecule contains only supported elements
     * @param mol Molecule to check
     * @return true if all elements are supported
     */
    bool isMoleculeSupported(const Mol& mol);
}