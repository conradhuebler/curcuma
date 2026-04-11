/*
 * <Modern Optimizer Interface - Strategy Pattern Implementation>
 * Copyright (C) 2025 Claude AI - Generated Code
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

#include "json.hpp"
#include "src/core/curcuma_logger.h"
#include <memory>
#include <string>
#include <vector>

using json = nlohmann::json;

// Forward declarations
class Molecule;
class EnergyCalculator;
typedef Eigen::VectorXd Vector;
typedef Eigen::Matrix<double, -1, -1> Matrix;

namespace Optimization {

/**
 * @brief Enum for optimizer type selection - Claude Generated
 * Type-safe method selection replacing magic numbers
 */
enum class OptimizerType {
    LBFGSPP, // "lbfgspp" - External LBFGSpp library
    INTERNAL_LBFGS, // "lbfgs"   - Internal LBFGS implementation
    DIIS, // "diis"    - Direct Inversion of Iterative Subspace
    RFO, // "rfo"     - Rational Function Optimization
    AUTO // "auto"    - Automatic selection based on system
};

/**
 * @brief Parse optimizer type from string - Claude Generated
 */
OptimizerType parseOptimizerType(const std::string& method_name);

/**
 * @brief Convert optimizer type to string - Claude Generated
 */
std::string optimizerTypeToString(OptimizerType type);

/**
 * @brief Optimization result container - Claude Generated
 * Structured result replacing string-based output accumulation
 */
struct OptimizationResult {
    bool success = false;
    std::string error_message;

    // Final state
    Molecule final_molecule;
    double final_energy = 0.0;
    Vector final_gradient;

    // Convergence information
    int iterations_performed = 0;
    double final_energy_change = 0.0; // kJ/mol
    double final_rmsd_change = 0.0; // Å
    double final_gradient_norm = 0.0; // Eh/Bohr

    // Trajectory
    std::vector<Molecule> trajectory;
    std::vector<double> energy_trajectory;

    // Timing
    double optimization_time_seconds = 0.0;

    // Factory methods for common cases
    static OptimizationResult success_result(const Molecule& final_mol, double energy,
        int iterations, double time_s);
    static OptimizationResult failed_result(const std::string& error);
};

/**
 * @brief Abstract optimizer interface - Claude Generated
 * Analog to QMInterface - unified polymorphic interface for all optimization methods
 */
class OptimizerInterface {
public:
    virtual ~OptimizerInterface() = default;

    // Multiple initialization methods (analog to QMInterface)
    virtual bool InitializeOptimization(const Molecule& molecule) = 0;
    virtual bool InitializeOptimization(const Molecule* molecule) = 0;
    virtual bool InitializeOptimization(const double* coordinates, int atom_count) = 0;

    // Geometry updates (analog to QMInterface::UpdateMolecule)
    virtual bool UpdateGeometry(const Molecule& molecule) = 0;
    virtual bool UpdateGeometry(const double* coordinates) = 0;

    // Core optimization method (analog to QMInterface::Calculation)
    virtual OptimizationResult Optimize(bool write_trajectory = false, int verbosity = 1) = 0;

    // Property accessors (analog to QMInterface::Charges, etc.)
    virtual Vector GetCurrentGradient() const = 0;
    virtual double GetCurrentEnergy() const = 0;
    virtual Matrix GetCurrentHessian() const = 0;
    virtual std::vector<Molecule> GetTrajectory() const = 0;

    // Configuration management (analog to QMInterface)
    virtual void LoadConfiguration(const json& config) = 0;
    virtual json GetDefaultConfiguration() const = 0;
    virtual json GetCurrentConfiguration() const = 0;

    // Method identification
    virtual std::string getName() const = 0;
    virtual OptimizerType getType() const = 0;
    virtual bool supportsConstraints() const = 0;
    virtual bool requiresHessian() const = 0;
    virtual std::vector<std::string> getRequiredParameters() const = 0;
};

/**
 * @brief Default configuration for all optimizers - Claude Generated
 * Analog to QMInterfaceJson - consistent defaults across methods
 */
extern const json OptimizerInterfaceJson;

} // namespace Optimization