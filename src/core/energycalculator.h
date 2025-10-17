/*
 * < Unified Energy and Gradient Calculator >
 * Copyright (C) 2022 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated: Big-Bang refactoring to unified ComputationalMethod interface
 */

#pragma once

#include "energy_calculators/computational_method.h"
#include "energy_calculators/method_factory.h"
#include "config_manager.h"

#include <memory>
#include <functional>

#include "json.hpp"
using json = nlohmann::json;

/**
 * @brief Unified Energy and Gradient Calculator
 * 
 * This class provides a single, consistent interface for all quantum mechanical
 * and molecular mechanics calculations in Curcuma. It uses the ComputationalMethod
 * interface to provide polymorphic access to all available methods.
 * 
 * Key improvements in refactored version:
 * - Single unified interface instead of dual QMInterface/ForceField system
 * - Priority-based method resolution (gfn2: TBLite > Ulysses > XTB)
 * - Simplified implementation without complex switch statements
 * - Maintained API compatibility with existing code
 * - Enhanced threading support for all methods
 * - Consistent parameter handling and error reporting
 * 
 * Supported methods:
 * - Force Fields: uff, uff-d3, qmdff, cgfnff (native GFN-FF)
 * - Quantum Methods: eht, gfn1, gfn2, ipea1
 * - External Libraries: XTB, TBLite, Ulysses methods
 * - Dispersion: d3, d4 corrections
 * 
 * Claude Generated: Complete architectural overhaul for consistency and extensibility
 */
class EnergyCalculator {
public:
    /**
     * @brief Constructor with method and JSON configuration (backward compatible)
     * @param method Method name (e.g., "gfn2", "uff", "eht", "cgfnff")
     * @param controller JSON configuration
     */
    EnergyCalculator(const std::string& method, const json& controller);

    /**
     * @brief Constructor with basename for parameter caching (backward compatible)
     * @param method Method name
     * @param controller JSON configuration
     * @param basename Base name for parameter file generation
     */
    EnergyCalculator(const std::string& method, const json& controller, const std::string& basename);

    /**
     * @brief Constructor with method and ConfigManager configuration (new, preferred)
     * Claude Generated: Phase 3C - Native ConfigManager support
     * @param method Method name
     * @param config ConfigManager configuration
     */
    EnergyCalculator(const std::string& method, const ConfigManager& config);

    /**
     * @brief Constructor with ConfigManager and basename (new, preferred)
     * Claude Generated: Phase 3C - Native ConfigManager support
     * @param method Method name
     * @param config ConfigManager configuration
     * @param basename Base name for parameter file generation
     */
    EnergyCalculator(const std::string& method, const ConfigManager& config, const std::string& basename);

    /**
     * @brief Destructor
     */
    ~EnergyCalculator();

    // =================================================================================
    // Molecular Setup and Geometry Management
    // =================================================================================
    
    /**
     * @brief Set molecular structure for calculation
     * @param mol Molecule object with geometry, atoms, charge, spin
     */
    void setMolecule(const Mol& mol);
    
    /**
     * @brief Update molecular geometry (raw coordinate array)
     * @param coord Coordinate array [x1,y1,z1,x2,y2,z2,...]
     */
    void updateGeometry(const double* coord);
    
    /**
     * @brief Update molecular geometry (vector format)
     * @param geometry Coordinate vector
     */
    void updateGeometry(const std::vector<double>& geometry);
    
    /**
     * @brief Update molecular geometry (matrix format)
     * @param geometry Coordinate matrix (natoms x 3)
     */
    void updateGeometry(const Matrix& geometry);
    
    /**
     * @brief Update molecular geometry (Eigen vector format)
     * @param geometry Eigen coordinate vector
     */
    void updateGeometry(const Eigen::VectorXd& geometry);

    // =================================================================================
    // Energy and Gradient Calculations
    // =================================================================================

    /**
     * @brief Perform energy and gradient calculation
     * @param gradient Calculate gradients if true
     * @return Total energy in appropriate units
     */
    double CalculateEnergy(bool gradient = false);

    /**
     * @brief Get calculated gradients
     * @return Gradient matrix (natoms x 3)
     */
    Matrix Gradient() const;
    
    /**
     * @brief Calculate numerical gradients (for testing/validation)
     * @return Numerical gradient matrix
     */
    Eigen::MatrixXd NumGrad();

    // =================================================================================
    // Property Access (maintained for API compatibility)
    // =================================================================================
    
    /**
     * @brief Get atomic partial charges
     * @return Vector of atomic charges
     */
    Vector Charges() const;
    
    /**
     * @brief Get molecular dipole moment
     * @return 3D dipole vector
     */
    Position Dipole() const;
    
    /**
     * @brief Get bond orders
     * @return Matrix/vector of bond orders (method-dependent format)
     */
    std::vector<std::vector<double>> BondOrders() const;
    
    /**
     * @brief Get orbital energies (QM methods only)
     * @return Vector of orbital energies
     */
    Vector Energies() const;
    
    /**
     * @brief Get orbital occupations (QM methods only)
     * @return Vector of orbital occupations
     */
    Vector OrbitalOccuptations() const;
    
    /**
     * @brief Get number of electrons (QM methods only)
     * @return Total number of electrons
     */
    int NumElectrons() const;

    // =================================================================================
    // Error Handling and Status
    // =================================================================================
    
    /**
     * @brief Check for NaN values in results
     * @return True if results contain NaN
     */
    bool HasNan() const { return m_containsNaN; }
    
    /**
     * @brief Check if calculator is in error state
     * @return True if error occurred
     */
    bool Error() const { return m_error; }
    
    /**
     * @brief Get error message
     * @return Error description or empty string
     */
    std::string ErrorMessage() const;
    
    /**
     * @brief Clear error state
     */
    void ClearError();

    // =================================================================================
    // Configuration and Parameters
    // =================================================================================
    
    /**
     * @brief Set method-specific parameters
     * @param parameter JSON parameter object
     */
    inline void setParameter(const json& parameter) {
        m_parameter = parameter;
        if (m_method) {
            m_method->setParameters(parameter);
        }
    }
    
    /**
     * @brief Get current parameters
     * @return JSON parameter object
     */
    inline json Parameter() const { return m_parameter; }
    
    /**
     * @brief Set geometry filename for parameter caching
     * @param filename Geometry file path
     */
    void setGeometryFile(const std::string& filename);
    
    /**
     * @brief Set basename for parameter file generation
     * @param basename Base name for parameter files
     */
    void setBasename(const std::string& basename);
    
    /**
     * @brief Set number of threads for calculation
     * @param threads Thread count (method-dependent)
     */
    void setThreadCount(int threads);
    
    /**
     * @brief Get current method name
     * @return Method name string
     */
    std::string getMethodName() const;
    
    /**
     * @brief Check if method supports threading
     * @return True if method is thread-safe
     */
    bool isThreadSafe() const;

    // =================================================================================
    // Legacy Interface Support (for backward compatibility)
    // =================================================================================
    
    /**
     * @brief Get internal interface pointer (legacy support)
     * @return ComputationalMethod pointer (replaces old QMInterface)
     * 
     * @deprecated Use direct EnergyCalculator methods instead
     */
    ComputationalMethod* Interface() const { return m_method.get(); }
    
    /**
     * @brief Print available methods
     */
    static void PrintAvailableMethods();

    // =================================================================================
    // Advanced Features
    // =================================================================================
    
    /**
     * @brief Save calculation results to file
     * @param filename Output file path
     * @return True if save successful
     */
    bool SaveToFile(const std::string& filename) const;
    
    /**
     * @brief Get method information and capabilities
     * @return JSON with method details
     */
    json GetMethodInfo() const;

    // =================================================================================
    // Submodule Verbosity Control (Claude Generated)
    // =================================================================================

    /**
     * @brief Override system verbosity for this EnergyCalculator instance
     *
     * This allows parent modules (like Hessian) to run sub-calculations silently
     * while maintaining their own output level.
     *
     * @param level Verbosity level (0=silent, 1=results, 2=analysis, 3=debug)
     */
    void setVerbosity(int level);

    /**
     * @brief Reset to system-wide verbosity from CurcumaLogger
     */
    void resetVerbosity();

    /**
     * @brief Get current effective verbosity level
     * @return Current verbosity (override or system level)
     */
    int getEffectiveVerbosity() const;

private:
    // =================================================================================
    // Internal State
    // =================================================================================
    
    std::unique_ptr<ComputationalMethod> m_method;    ///< Unified computational method
    json m_controller;                                ///< Configuration parameters
    json m_parameter;                                 ///< Method-specific parameters
    
    std::string m_method_name;                        ///< Method name
    std::string m_geometry_file;                      ///< Geometry file for caching
    std::string m_basename;                           ///< Base name for parameter files
    
    // Cached calculation results
    Mol m_mol;                                        ///< Current molecule
    Matrix m_gradient;                                ///< Last calculated gradient
    double m_energy;                                  ///< Last calculated energy
    
    // Status flags
    bool m_initialized = false;                       ///< Initialization status
    bool m_containsNaN = false;                      ///< NaN detection flag
    bool m_error = false;                            ///< Error state flag
    std::string m_error_message;                     ///< Error description
    
    int m_atoms = 0;                                 ///< Number of atoms
    int m_mult = 1;                                  ///< Multiplicity

    // Verbosity control (Claude Generated)
    int m_verbosity_override = -1; ///< Override verbosity (-1 = use system)

    // =================================================================================
    // Internal Methods
    // =================================================================================
    
    /**
     * @brief Initialize EnergyCalculator with ConfigManager settings (new, preferred)
     * Claude Generated: Phase 3C - Native ConfigManager support
     * @param config ConfigManager configuration
     */
    void initializeCommonFromConfig(const ConfigManager& config);

    /**
     * @brief Initialize EnergyCalculator with JSON settings (backward compatible)
     * @param controller Configuration JSON
     */
    void initializeCommon(const json& controller);
    
    /**
     * @brief Create computational method using factory
     * @param method_name Method to create
     * @param config Configuration for method
     * @return True if method created successfully
     */
    bool createMethod(const std::string& method_name, const json& config);
    
    /**
     * @brief Handle errors from computational method
     * @param operation Description of failed operation
     */
    void handleMethodError(const std::string& operation);
    
    /**
     * @brief Check results for NaN values
     * @param energy Calculated energy
     * @param gradient Calculated gradient
     * @return True if NaN values found
     */
    bool checkForNaN(double energy, const Matrix& gradient = Matrix{});
    
    /**
     * @brief Convert coordinate formats for updateGeometry
     * @param coord Input coordinates
     * @param geometry Output matrix
     */
    void convertCoordinates(const double* coord, Matrix& geometry);
    void convertCoordinates(const std::vector<double>& coord, Matrix& geometry);
    void convertCoordinates(const Eigen::VectorXd& coord, Matrix& geometry);
};