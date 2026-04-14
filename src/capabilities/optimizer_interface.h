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
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include <memory>
#include <string>
#include <vector>

using json = nlohmann::json;
using curcuma::Molecule;

namespace Optimization {

/**
 * @brief Enum for optimizer type selection - Claude Generated
 * Type-safe method selection replacing magic numbers
 */
enum class OptimizerType {
    LBFGSPP,      // External LBFGSpp library
    NATIVE_LBFGS, // Curcuma's own L-BFGS (Nocedal & Wright)
    NATIVE_DIIS,  // Curcuma's own DIIS (Pulay 1980)
    NATIVE_RFO,   // Curcuma's own RFO (Banerjee et al. 1985)
    ANCOPT,       // Approximate Normal Coordinate Optimizer (Grimme)
    AUTO          // Automatic selection based on system size
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
 * @brief Default configuration for all optimizers - Claude Generated
 */
extern const json OptimizerInterfaceJson;

} // namespace Optimization