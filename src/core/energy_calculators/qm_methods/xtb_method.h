/*
 * < XTB Method Wrapper for ComputationalMethod Interface >
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

#ifdef USE_XTB
#include "xtbinterface.h"
#endif

#include <memory>

/**
 * @brief XTB method wrapper for ComputationalMethod interface
 * 
 * This wrapper adapts the existing XTBInterface to the unified
 * ComputationalMethod interface. XTB provides extended tight-binding
 * methods including GFN-FF, GFN1, and GFN2.
 * 
 * Supported XTB methods:
 * - gfnff: GFN-FF force field (fast, good for large systems)
 * - gfn1/xtb-gfn1: GFN1-xTB tight-binding method
 * - gfn2/xtb-gfn2: GFN2-xTB tight-binding method (most accurate)
 * 
 * Features:
 * - Analytical gradients available for all methods
 * - Atomic charges, bond orders, and dipole moments
 * - Orbital energies and occupations
 * - Temperature and accuracy control
 * - Threading support (configurable)
 * 
 * Claude Generated: Big-Bang EnergyCalculator refactoring wrapper
 */
class XTBMethod : public ComputationalMethod {
public:
    /**
     * @brief Constructor with method name and configuration
     * @param method_name XTB method ("gfnff", "xtb-gfn1", "xtb-gfn2", etc.)
     * @param config JSON configuration (XTB-specific parameters)
     */
    XTBMethod(const std::string& method_name, const json& config = json{});
    
    /**
     * @brief Destructor
     */
    virtual ~XTBMethod() = default;
    
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
    bool hasGradient() const override { return true; } // XTB always provides gradients
    
    // Method information
    std::string getMethodName() const override { return m_method_name; }
    bool isThreadSafe() const override;
    void setThreadCount(int threads) override;
    
    // Configuration
    void setParameters(const json& params) override;
    json getParameters() const override;
    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;

    // Energy decomposition (JSON output - placeholder for native implementation)
    json getEnergyDecomposition() const override;

    // =================================================================================
    // XTB-specific Methods (additional functionality)
    // =================================================================================
    
    /**
     * @brief Get orbital energies
     * @return Vector of orbital energies in Hartree
     */
    Vector getOrbitalEnergies() const override;
    
    /**
     * @brief Get orbital occupations
     * @return Vector of orbital occupations (0-2 for each orbital)
     */
    Vector getOrbitalOccupations() const override;
    
    
    
    
    
    
    
    
    
    /**
     * @brief Get supported XTB methods
     * @return Vector of supported method names
     */
    static std::vector<std::string> getSupportedMethods();
    
    
    /**
     * @brief Check if XTB is available (compilation flag)
     * @return True if XTB was compiled in
     */
    static bool isAvailable();
    
    /**
     * @brief Save XTB results to file
     * @param filename Output file name
     * @return True if save successful
     */
    bool saveToFile(const std::string& filename) const override;
    
private:
#ifdef USE_XTB
    std::unique_ptr<XTBInterface> m_xtb;           ///< Wrapped XTB implementation
#endif
    std::string m_method_name;                     ///< XTB method name
    Mol m_molecule;                                ///< Current molecule
    bool m_calculation_done;                       ///< Flag if calculation was performed
    double m_last_energy;                          ///< Last calculated energy
    
    
    
    
    
    
    
    
};

