/*
 * <Modern Optimizer Interface Implementation>
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

#include "optimizer_interface.h"
#include <algorithm>
#include <chrono>
#include <stdexcept>

namespace Optimization {

// Default configuration analog to QMInterfaceJson
const json OptimizerInterfaceJson{
    // Threading and performance
    { "threads", 1 },
    { "max_iterations", 5000 },

    // Convergence criteria (using unit system from CLAUDE.md)
    { "energy_threshold", 0.1 }, // kJ/mol - relative energy changes
    { "rmsd_threshold", 0.01 }, // Å - geometric changes
    { "gradient_threshold", 5e-4 }, // Eh/Bohr - gradient norm
    { "convergence_count", 11 }, // Bit field for convergence criteria

    // Output control
    { "write_trajectory", true },
    { "verbose", false },
    { "print_output", true },

    // Method selection
    { "method", "lbfgspp" },
    { "optimizer_type", "lbfgspp" },

    // LBFGS-specific parameters (can be overridden by specific strategies)
    { "lbfgs_m", 2000 },
    { "lbfgs_eps_abs", 1e-5 },
    { "lbfgs_eps_rel", 1e-5 },
    { "lbfgs_past", 0 },
    { "lbfgs_delta", 0 },
    { "lbfgs_line_search", 3 },
    { "lbfgs_max_line_search", 20 },
    { "lbfgs_min_step", 1e-20 },
    { "lbfgs_max_step", 1e20 },
    { "lbfgs_ftol", 1e-4 },
    { "lbfgs_wolfe", 0.9 },

    // Advanced options
    { "use_constraints", false },
    { "use_hessian", false },
    { "single_step_mode", 1 },
    { "max_energy_rise", 100 },

    // Internal optimizer parameters (for InternalLBFGS, DIIS, RFO)
    { "diis_history", 5 },
    { "diis_start", 5 },
    { "rfo_lambda", 0.1 },
    { "trust_radius", 0.1 },

    // Molecular properties
    { "charge", 0 },
    { "spin", 0 }
};

// Claude Generated - String to enum conversion
OptimizerType parseOptimizerType(const std::string& method_name)
{
    std::string lower_name = method_name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);

    if (lower_name == "lbfgspp" || lower_name == "lbfgs++" || lower_name == "external_lbfgs") {
        return OptimizerType::LBFGSPP;
    } else if (lower_name == "lbfgs" || lower_name == "internal_lbfgs" || lower_name == "gpt_lbfgs") {
        return OptimizerType::INTERNAL_LBFGS;
    } else if (lower_name == "diis" || lower_name == "direct_inversion") {
        return OptimizerType::DIIS;
    } else if (lower_name == "rfo" || lower_name == "rational_function") {
        return OptimizerType::RFO;
    } else if (lower_name == "ancopt" || lower_name == "anc" || lower_name == "approximate_normal") {
        return OptimizerType::ANCOPT;
    } else if (lower_name == "auto" || lower_name == "automatic") {
        return OptimizerType::AUTO;
    } else {
        throw std::invalid_argument("Unknown optimizer type: " + method_name + ". Valid options: lbfgspp, lbfgs, diis, rfo, ancopt, auto");
    }
}

// Claude Generated - Enum to string conversion
std::string optimizerTypeToString(OptimizerType type)
{
    switch (type) {
    case OptimizerType::LBFGSPP:
        return "lbfgspp";
    case OptimizerType::INTERNAL_LBFGS:
        return "lbfgs";
    case OptimizerType::DIIS:
        return "diis";
    case OptimizerType::RFO:
        return "rfo";
    case OptimizerType::ANCOPT:
        return "ancopt";
    case OptimizerType::AUTO:
        return "auto";
    default:
        return "unknown";
    }
}

// Claude Generated - Factory methods for OptimizationResult
OptimizationResult OptimizationResult::success_result(const Molecule& final_mol,
    double energy,
    int iterations,
    double time_s)
{
    OptimizationResult result;
    result.success = true;
    result.final_molecule = final_mol;
    result.final_energy = energy;
    result.iterations_performed = iterations;
    result.optimization_time_seconds = time_s;
    return result;
}

OptimizationResult OptimizationResult::failed_result(const std::string& error)
{
    OptimizationResult result;
    result.success = false;
    result.error_message = error;
    return result;
}

} // namespace Optimization