/*
 * < ForceField Method Wrapper for ComputationalMethod Interface >
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
#include "forcefield.h"
#include "forcefieldgenerator.h"

#include <memory>

/**
 * @brief ForceField wrapper for ComputationalMethod interface
 * 
 * This wrapper adapts the existing ForceField implementation (UFF, QMDFF, etc.)
 * to the unified ComputationalMethod interface. It maintains the high-performance
 * threading capabilities of the original ForceField using CxxThreadPool.
 * 
 * Key Features:
 * - Maintains existing CxxThreadPool threading system for optimal performance
 * - Supports all ForceField methods: UFF, UFF-D3, QMDFF, cgfnff
 * - Preserves universal parameter caching system (96% speedup)
 * - Thread-safe parameter caching control for concurrent calculations
 * - Analytical gradients always available
 * 
 * Claude Generated: Big-Bang EnergyCalculator refactoring wrapper
 */
class ForceFieldMethod : public ComputationalMethod {
public:
    /**
     * @brief Constructor with method name and configuration
     * @param method_name ForceField method ("uff", "uff-d3", "qmdff", "cgfnff")
     * @param config JSON configuration (ForceField-specific parameters)
     */
    ForceFieldMethod(const std::string& method_name, const json& config = json{});
    
    /**
     * @brief Destructor
     */
    virtual ~ForceFieldMethod() = default;
    
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
    bool hasGradient() const override { return true; } // ForceFields always have gradients
    
    // Method information
    std::string getMethodName() const override { return m_method_name; }
    bool isThreadSafe() const override { return true; } // ForceField is thread-optimized
    void setThreadCount(int threads) override;
    
    // Configuration
    void setParameters(const json& params) override;
    json getParameters() const override;
    bool hasError() const override;
    void clearError() override;
    std::string getErrorMessage() const override;
    
    // =================================================================================
    // ForceField-specific Methods (additional functionality)
    // =================================================================================
    
    /**
     * @brief Enable/disable parameter caching (important for multi-threading)
     * @param enable True to enable caching, false for concurrent calculations
     * 
     * Note: Parameter caching provides 96% speedup but must be disabled
     * when running multiple ForceField instances concurrently
     */
    void setParameterCaching(bool enable);
    
    /**
     * @brief Check if parameter caching is enabled
     * @return True if caching is enabled
     */
    bool isParameterCachingEnabled() const;
    
    /**
     * @brief Load parameters from file
     * @param filename Parameter file path
     * @return True if load successful
     */
    bool loadParametersFromFile(const std::string& filename);
    
    /**
     * @brief Save current parameters to file
     * @param filename Output parameter file path
     * @return True if save successful
     */
    bool saveParametersToFile(const std::string& filename) const;
    
    /**
     * @brief Export current parameters as JSON
     * @return JSON object with current parameters
     */
    json exportCurrentParameters() const;
    
    /**
     * @brief Check if parameters are loaded/available
     * @return True if parameters are available
     */
    bool hasParameters() const;
    
    /**
     * @brief Try to load auto-generated parameters
     * @param method Method name for parameter matching
     * @return True if auto-parameters loaded successfully
     */
    bool tryLoadAutoParameters(const std::string& method = "");
    
    /**
     * @brief Auto-save parameters with intelligent naming
     * @return True if auto-save successful
     */
    bool autoSaveParameters() const;
    
    /**
     * @brief Print parameter summary for analysis
     */
    void printParameterSummary() const;
    
    /**
     * @brief Calculate numerical gradient (for testing/validation)
     * @return Numerical gradient matrix
     */
    Matrix calculateNumericalGradient();
    
    /**
     * @brief Set parameter file path
     * @param filename Path to parameter file
     */
    void setParameterFile(const std::string& filename);
    
    /**
     * @brief Generate ForceField parameters if needed
     * @param mol Molecule for parameter generation
     * @return True if parameters available (generated or cached)
     * @note Claude Generated: UFF parameter generation integration
     */
    bool generateParametersIfNeeded(const Mol& mol);
    
    /**
     * @brief Get supported ForceField methods
     * @return Vector of supported method names
     */
    static std::vector<std::string> getSupportedMethods();

    /**
     * @brief Check if method is supported by ForceField
     * @param method_name Method to check
     * @return True if method is supported
     */
    static bool isMethodSupported(const std::string& method_name);

    // =================================================================================
    // Energy Component Access (Claude Generated November 2025)
    // =================================================================================

    /**
     * @brief Get bond energy component
     * @return Bond stretching energy
     */
    double getBondEnergy() const;

    /**
     * @brief Get angle energy component
     * @return Angle bending energy
     */
    double getAngleEnergy() const;

    /**
     * @brief Get dihedral energy component
     * @return Dihedral torsion energy
     */
    double getDihedralEnergy() const;

    /**
     * @brief Get inversion energy component
     * @return Inversion/out-of-plane energy
     */
    double getInversionEnergy() const;

    /**
     * @brief Get van der Waals energy component
     * @return Van der Waals interaction energy
     */
    double getVdWEnergy() const;

    /**
     * @brief Get repulsion energy component
     * @return Core-core repulsion energy
     */
    double getRepulsionEnergy() const;

    /**
     * @brief Get dispersion energy component
     * @return Dispersion correction energy
     */
    double getDispersionEnergy() const;

    /**
     * @brief Get Coulomb energy component
     * @return Electrostatic energy
     */
    double getCoulombEnergy() const;
    
private:
    std::unique_ptr<ForceField> m_forcefield;      ///< Wrapped ForceField implementation
    std::string m_method_name;                     ///< ForceField method name
    Mol m_molecule;                                ///< Current molecule
    bool m_calculation_done;                       ///< Flag if calculation was performed
    double m_last_energy;                          ///< Last calculated energy
    Matrix m_last_gradient;                        ///< Last calculated gradient
    
    /**
     * @brief Initialize ForceField with current parameters
     * @return True if initialization successful
     */
    bool initializeForceField();
    
    /**
     * @brief Update ForceField configuration from JSON parameters
     */
    void updateForceFieldConfig();
    
    /**
     * @brief Generate ForceField controller JSON from parameters
     * @return Controller JSON for ForceField constructor
     */
    json generateForceFieldController() const;
    
    /**
     * @brief Handle ForceField-specific errors
     * @param operation Description of failed operation
     */
    void handleForceFieldError(const std::string& operation);
    
    /**
     * @brief Validate method name for ForceField compatibility
     * @param method_name Method to validate
     * @return True if method is valid
     */
    bool validateMethodName(const std::string& method_name) const;
    
    /**
     * @brief Get default configuration for specific ForceField method
     * @param method_name ForceField method name
     * @return Default JSON configuration
     */
    static json getDefaultConfigForMethod(const std::string& method_name);
    
    /**
     * @brief Convert method name to ForceField-compatible format
     * @param method_name Input method name
     * @return ForceField-compatible method name
     */
    std::string normalizeMethodName(const std::string& method_name) const;
};

/**
 * @brief Utility functions for ForceField method
 */
namespace ForceFieldMethodUtils {
    
    /**
     * @brief Validate ForceField parameters
     * @param config JSON configuration to validate
     * @return True if parameters are valid
     */
    bool validateForceFieldConfig(const json& config);
    
    /**
     * @brief Get performance-optimized thread count for system
     * @return Recommended number of threads
     */
    int getOptimalThreadCount();
    
    /**
     * @brief Generate parameter file name from geometry file
     * @param geometry_file Geometry file path
     * @return Corresponding parameter file path
     */
    std::string generateParameterFileName(const std::string& geometry_file);
    
    /**
     * @brief Check if molecular system is suitable for ForceField methods
     * @param mol Molecule to validate
     * @return True if molecule is suitable
     */
    bool isMoleculeSuitableForFF(const Mol& mol);
    
    /**
     * @brief Estimate memory usage for ForceField calculation
     * @param mol Molecule for calculation
     * @param threads Number of threads
     * @return Estimated memory usage in MB
     */
    double estimateMemoryUsage(const Mol& mol, int threads = 1);
}