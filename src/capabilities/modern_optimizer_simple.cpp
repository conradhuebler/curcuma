/*
 * <Modern Optimizer - Simplified Implementation>
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

#include "modern_optimizer_simple.h"
#include "curcumaopt.h" // For legacy LBFGS functionality
// #include "native_lbfgs_optimizer.h" // For native implementations (disabled for now)
#include <algorithm>
#include <chrono>

using curcuma::Molecule;

namespace ModernOptimization {

// Claude Generated - SimpleOptimizationResult factory methods
SimpleOptimizationResult SimpleOptimizationResult::success_result(const Molecule& mol, double energy,
    const std::string& method, int iterations, double time)
{
    SimpleOptimizationResult result;
    result.success = true;
    result.final_molecule = mol;
    result.final_energy = energy;
    result.method_used = method;
    result.iterations_performed = iterations;
    result.optimization_time_seconds = time;
    return result;
}

SimpleOptimizationResult SimpleOptimizationResult::failed_result(const std::string& error, const std::string& method)
{
    SimpleOptimizationResult result;
    result.success = false;
    result.error_message = error;
    result.method_used = method;
    return result;
}

// Claude Generated - ModernOptimizerDispatcher implementation
ModernOptimizerDispatcher::OptimizerType ModernOptimizerDispatcher::parseOptimizerType(const std::string& method_name)
{
    std::string lower_name = method_name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);

    if (lower_name == "lbfgspp" || lower_name == "lbfgs++" || lower_name == "external") {
        return OptimizerType::LBFGSPP;
    } else if (lower_name == "internal" || lower_name == "gpt") {
        return OptimizerType::INTERNAL;
    } else if (lower_name == "native_lbfgs" || lower_name == "native-lbfgs" || lower_name == "lbfgs") {
        return OptimizerType::NATIVE_LBFGS;
    } else if (lower_name == "native_diis" || lower_name == "native-diis" || lower_name == "diis") {
        return OptimizerType::NATIVE_DIIS;
    } else if (lower_name == "native_rfo" || lower_name == "native-rfo" || lower_name == "rfo") {
        return OptimizerType::NATIVE_RFO;
    } else if (lower_name == "auto" || lower_name == "automatic") {
        return OptimizerType::AUTO;
    } else {
        // Default to native LBFGS for unknown methods
        CurcumaLogger::warn_fmt("Unknown optimizer method '{}', defaulting to native L-BFGS", method_name);
        return OptimizerType::NATIVE_LBFGS;
    }
}

std::map<std::string, std::string> ModernOptimizerDispatcher::getAvailableOptimizers()
{
    return {
        { "lbfgspp", "External LBFGSpp library - robust L-BFGS implementation" },
        { "internal", "Internal LBFGS implementation - custom features (placeholder)" },
        { "lbfgs", "Native L-BFGS implementation (Claude 3.5 generated)" },
        { "diis", "Native DIIS with L-BFGS fallback (Claude 3.5 generated)" },
        { "rfo", "Native Rational Function Optimization (Claude 3.5 generated)" },
        { "auto", "Automatic selection based on system characteristics" }
    };
}

SimpleOptimizationResult ModernOptimizerDispatcher::optimizeStructure(
    Molecule* molecule,
    const std::string& method_name,
    EnergyCalculator* energy_calculator,
    const json& config)
{

    if (!molecule) {
        return SimpleOptimizationResult::failed_result("Null molecule provided", method_name);
    }
    if (!energy_calculator) {
        return SimpleOptimizationResult::failed_result("Null energy calculator provided", method_name);
    }

    try {
        OptimizerType type = parseOptimizerType(method_name);

        // Auto-selection logic
        if (type == OptimizerType::AUTO) {
            if (molecule->AtomCount() > 200) {
                type = OptimizerType::LBFGSPP; // External library for very large systems
                CurcumaLogger::info_fmt("Auto-selected LBFGSpp for large {}-atom system", static_cast<int>(molecule->AtomCount()));
            } else if (molecule->AtomCount() > 50) {
                type = OptimizerType::NATIVE_LBFGS; // Native L-BFGS for medium systems
                CurcumaLogger::info_fmt("Auto-selected native L-BFGS for {}-atom system", static_cast<int>(molecule->AtomCount()));
            } else {
                type = OptimizerType::NATIVE_DIIS; // DIIS for small systems
                CurcumaLogger::info_fmt("Auto-selected native DIIS for small {}-atom system", static_cast<int>(molecule->AtomCount()));
            }
        }

        // Log optimization start
        logOptimizationHeader(method_name, *molecule);

        // Route to specific optimizer
        SimpleOptimizationResult result;
        switch (type) {
        case OptimizerType::LBFGSPP:
            result = optimizeWithLBFGSpp(molecule, energy_calculator, config);
            break;
        case OptimizerType::INTERNAL:
            result = optimizeWithInternal(molecule, energy_calculator, config);
            break;
        case OptimizerType::NATIVE_LBFGS:
            CurcumaLogger::warn("Native L-BFGS temporarily disabled - falling back to LBFGSpp");
            result = optimizeWithLBFGSpp(molecule, energy_calculator, config);
            break;
        case OptimizerType::NATIVE_DIIS:
            CurcumaLogger::warn("Native DIIS temporarily disabled - falling back to LBFGSpp");
            result = optimizeWithLBFGSpp(molecule, energy_calculator, config);
            break;
        case OptimizerType::NATIVE_RFO:
            CurcumaLogger::warn("Native RFO temporarily disabled - falling back to LBFGSpp");
            result = optimizeWithLBFGSpp(molecule, energy_calculator, config);
            break;
        default:
            result = SimpleOptimizationResult::failed_result("Unknown optimizer type", method_name);
            break;
        }

        // Log results
        logOptimizationResult(result);

        return result;

    } catch (const std::exception& e) {
        return SimpleOptimizationResult::failed_result(
            fmt::format("Optimization error: {}", e.what()), method_name);
    }
}

void ModernOptimizerDispatcher::demonstrateModernFeatures(const Molecule& molecule)
{
    CurcumaLogger::header("Modern Optimizer System Features");

    // Show available optimizers
    auto available = getAvailableOptimizers();
    CurcumaLogger::info("Available optimization methods:");
    for (const auto& pair : available) {
        CurcumaLogger::info_fmt("  {} - {}", pair.first, pair.second);
    }

    // Show modern logging features
    CurcumaLogger::info("");
    CurcumaLogger::info("Modern Features Demonstration:");
    CurcumaLogger::success("✓ Strategy Pattern architecture for extensible methods");
    CurcumaLogger::success("✓ CurcumaLogger integration with verbosity control");
    CurcumaLogger::success("✓ Unit-aware output formatting");
    CurcumaLogger::success("✓ Type-safe method selection (no magic numbers)");
    CurcumaLogger::success("✓ Comprehensive error handling and validation");

    // Show unit-aware output examples
    CurcumaLogger::info("");
    CurcumaLogger::info("Unit-aware output examples:");
    CurcumaLogger::energy_abs(-154.123456, "Sample energy");
    CurcumaLogger::energy_rel(0.00123, "Sample energy change");
    CurcumaLogger::length(3.2 / CURCUMA_BOHR_TO_ANGSTROM, "Sample distance");

    // Show parameter examples
    CurcumaLogger::info("");
    json sample_params = {
        { "method", "lbfgspp" },
        { "max_iterations", 5000 },
        { "energy_threshold", 0.1 },
        { "gradient_threshold", 5e-4 }
    };
    CurcumaLogger::param_table(sample_params, "Sample Configuration");
}

// Claude Generated - Implementation methods
SimpleOptimizationResult ModernOptimizerDispatcher::optimizeWithLBFGSpp(Molecule* molecule,
    EnergyCalculator* calc,
    const json& config)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        // Create legacy optimizer to do the actual work
        json legacy_config = config;
        legacy_config["optimethod"] = 0; // Force LBFGSpp mode
        legacy_config["SinglePoint"] = false;

        CurcumaOpt legacy_opt(legacy_config, true); // Silent mode

        // Perform optimization
        std::string output;
        std::vector<Molecule> intermediates;
        Vector charges;

        Molecule result_mol = legacy_opt.LBFGSOptimise(molecule, output, &intermediates, charges, -1, "modern_opt");

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double optimization_time = duration.count() / 1000.0;

        // Parse iteration count from output (simple heuristic)
        int iterations = intermediates.size();
        if (iterations == 0)
            iterations = 1; // At least one iteration

        CurcumaLogger::success("LBFGSpp optimization completed successfully");
        CurcumaLogger::info_fmt("Iterations: {}", iterations);
        CurcumaLogger::info_fmt("Time: {:.3f} seconds", optimization_time);

        return SimpleOptimizationResult::success_result(
            result_mol, result_mol.Energy(), "lbfgspp", iterations, optimization_time);

    } catch (const std::exception& e) {
        return SimpleOptimizationResult::failed_result(
            fmt::format("LBFGSpp optimization failed: {}", e.what()), "lbfgspp");
    }
}

SimpleOptimizationResult ModernOptimizerDispatcher::optimizeWithInternal(Molecule* molecule,
    EnergyCalculator* calc,
    const json& config)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        // Create legacy optimizer to do the actual work
        json legacy_config = config;
        legacy_config["optimethod"] = 1; // Force internal LBFGS mode
        legacy_config["SinglePoint"] = false;

        CurcumaOpt legacy_opt(legacy_config, true); // Silent mode

        // Perform optimization
        std::string output;
        std::vector<Molecule> intermediates;
        Vector charges;

        Molecule result_mol = legacy_opt.GPTLBFGS(molecule, output, &intermediates, charges, -1, "modern_opt");

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double optimization_time = duration.count() / 1000.0;

        // Parse iteration count from output
        int iterations = intermediates.size();
        if (iterations == 0)
            iterations = 1;

        CurcumaLogger::success("Internal LBFGS optimization completed successfully");
        CurcumaLogger::info_fmt("Iterations: {}", iterations);
        CurcumaLogger::info_fmt("Time: {:.3f} seconds", optimization_time);

        return SimpleOptimizationResult::success_result(
            result_mol, result_mol.Energy(), "internal", iterations, optimization_time);

    } catch (const std::exception& e) {
        return SimpleOptimizationResult::failed_result(
            fmt::format("Internal LBFGS optimization failed: {}", e.what()), "internal");
    }
}

// Claude Generated - Helper methods
json ModernOptimizerDispatcher::getDefaultConfig(OptimizerType type)
{
    json config = {
        { "writeXYZ", true },
        { "printOutput", false }, // Suppress legacy output
        { "dE", 0.1 },
        { "dRMSD", 0.01 },
        { "GradNorm", 5e-4 },
        { "MaxIter", 5000 },
        { "threads", 1 }
    };

    switch (type) {
    case OptimizerType::LBFGSPP:
        config["optimethod"] = 0;
        config["LBFGS_m"] = 2000;
        config["LBFGS_eps_abs"] = 1e-5;
        break;
    case OptimizerType::INTERNAL:
        config["optimethod"] = 1;
        config["lambda"] = 0.1;
        config["diis_hist"] = 5;
        break;
    default:
        break;
    }

    return config;
}

void ModernOptimizerDispatcher::logOptimizationHeader(const std::string& method, const Molecule& molecule)
{
    CurcumaLogger::header(fmt::format("Modern Optimization: {}", method));
    CurcumaLogger::param("Atoms", static_cast<int>(molecule.AtomCount()));
    CurcumaLogger::param("Method", method);

    // Show initial energy if available
    if (molecule.Energy() != 0.0) {
        CurcumaLogger::energy_abs(molecule.Energy(), "Initial energy");
    }
}

void ModernOptimizerDispatcher::logOptimizationResult(const SimpleOptimizationResult& result)
{
    if (result.success) {
        CurcumaLogger::success_fmt("Optimization completed with {}", result.method_used);
        CurcumaLogger::energy_abs(result.final_energy, "Final energy");
        CurcumaLogger::param("Iterations", result.iterations_performed);
        CurcumaLogger::param("Time", fmt::format("{:.3f} seconds", result.optimization_time_seconds));
    } else {
        CurcumaLogger::error_fmt("Optimization failed with {}: {}", result.method_used, result.error_message);
    }
}

// Claude Generated - Native optimization methods (temporarily disabled)
SimpleOptimizationResult ModernOptimizerDispatcher::optimizeWithNativeLBFGS(Molecule* molecule,
    EnergyCalculator* calc,
    const json& config)
{
    CurcumaLogger::warn("Native L-BFGS implementation temporarily disabled due to compilation issues");
    CurcumaLogger::info("Falling back to LBFGSpp implementation");
    return optimizeWithLBFGSpp(molecule, calc, config);
}

SimpleOptimizationResult ModernOptimizerDispatcher::optimizeWithNativeDIIS(Molecule* molecule,
    EnergyCalculator* calc,
    const json& config)
{
    CurcumaLogger::warn("Native DIIS implementation temporarily disabled due to compilation issues");
    CurcumaLogger::info("Falling back to LBFGSpp implementation");
    return optimizeWithLBFGSpp(molecule, calc, config);
}

SimpleOptimizationResult ModernOptimizerDispatcher::optimizeWithNativeRFO(Molecule* molecule,
    EnergyCalculator* calc,
    const json& config)
{
    CurcumaLogger::warn("Native RFO implementation temporarily disabled due to compilation issues");
    CurcumaLogger::info("Falling back to LBFGSpp implementation");
    return optimizeWithLBFGSpp(molecule, calc, config);
}

} // namespace ModernOptimization