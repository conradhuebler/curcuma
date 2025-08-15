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
#include "src/core/qm_methods/xtbinterface.h"
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
    double calculateEnergy(bool gradient = false, bool verbose = false) override;
    
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
     * @brief Set electronic temperature for calculation
     * @param temperature Electronic temperature in Kelvin
     */
    void setTemperature(double temperature);
    
    /**
     * @brief Get current electronic temperature
     * @return Electronic temperature in Kelvin
     */
    double getTemperature() const;
    
    /**
     * @brief Set accuracy level for calculation
     * @param accuracy Accuracy level (1=rough, 2=normal, 3=high)
     */
    void setAccuracy(int accuracy);
    
    /**
     * @brief Get current accuracy level
     * @return Accuracy level (1-3)
     */
    int getAccuracy() const;
    
    /**
     * @brief Set maximum SCF iterations
     * @param maxiter Maximum number of SCF iterations
     */
    void setMaxSCFIterations(int maxiter);
    
    /**
     * @brief Get maximum SCF iterations
     * @return Maximum SCF iterations
     */
    int getMaxSCFIterations() const;
    
    /**
     * @brief Set multiplicity (spin state)
     * @param mult Multiplicity (1=singlet, 2=doublet, 3=triplet, etc.)
     */
    void setMultiplicity(int mult);
    
    /**
     * @brief Get current multiplicity
     * @return Multiplicity
     */
    int getMultiplicity() const;
    
    /**
     * @brief Get supported XTB methods
     * @return Vector of supported method names
     */
    static std::vector<std::string> getSupportedMethods();
    
    /**
     * @brief Check if method is supported by XTB
     * @param method_name Method to check
     * @return True if method is supported
     */
    static bool isMethodSupported(const std::string& method_name);
    
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
    
    /**
     * @brief Initialize XTB interface with current parameters
     * @return True if initialization successful
     */
    bool initializeXTB();
    
    /**
     * @brief Update XTB configuration from JSON parameters
     */
    void updateXTBConfig();
    
    /**
     * @brief Generate XTB settings JSON from parameters
     * @return XTB settings JSON
     */
    json generateXTBSettings() const;
    
    /**
     * @brief Convert method name to XTB-compatible format
     * @param method_name Input method name
     * @return XTB-compatible method name
     */
    std::string convertMethodNameForXTB(const std::string& method_name) const;
    
    /**
     * @brief Handle XTB-specific errors
     * @param operation Description of failed operation
     */
    void handleXTBError(const std::string& operation);
    
    /**
     * @brief Validate method name for XTB compatibility
     * @param method_name Method to validate
     * @return True if method is valid for XTB
     */
    bool validateMethodName(const std::string& method_name) const;
    
    /**
     * @brief Get default configuration for specific XTB method
     * @param method_name XTB method name
     * @return Default JSON configuration
     */
    static json getDefaultConfigForMethod(const std::string& method_name);
    
    /**
     * @brief Check if XTB supports threading (depends on compilation)
     * @return True if XTB was compiled with threading support
     */
    bool supportsThreading() const;
};

/**
 * @brief Utility functions for XTB method
 */
namespace XTBMethodUtils {
    
    /**
     * @brief Validate XTB parameters
     * @param config JSON configuration to validate
     * @return True if parameters are valid
     */
    bool validateXTBConfig(const json& config);
    
    /**
     * @brief Get XTB method information
     * @param method_name XTB method to query
     * @return JSON with method details (accuracy, speed, etc.)
     */
    json getXTBMethodInfo(const std::string& method_name);
    
    /**
     * @brief Estimate XTB calculation time
     * @param mol Molecule for calculation
     * @param method_name XTB method
     * @return Estimated calculation time in seconds
     */
    double estimateCalculationTime(const Mol& mol, const std::string& method_name);
    
    /**
     * @brief Check if molecular system is suitable for XTB
     * @param mol Molecule to validate
     * @return True if molecule is suitable for XTB
     */
    bool isMoleculeSuitableForXTB(const Mol& mol);
    
    /**
     * @brief Get recommended XTB method for molecular system
     * @param mol Molecule to analyze
     * @param priority Priority ("speed", "accuracy", "balance")
     * @return Recommended XTB method name
     */
    std::string getRecommendedMethod(const Mol& mol, const std::string& priority = "balance");
}