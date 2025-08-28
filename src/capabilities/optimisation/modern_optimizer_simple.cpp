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
#include "../curcumaopt.h" // For legacy LBFGS functionality
#include "lbfgs.h" // Native LBFGS implementation - Claude Generated
#include <algorithm>
#include <chrono>

using curcuma::Molecule;

namespace ModernOptimization {

// Claude Generated - ModernOptimizerDispatcher CurcumaMethod integration
ModernOptimizerDispatcher::ModernOptimizerDispatcher(const json& controller, bool silent)
    : CurcumaMethod(OptimizerJson, controller, silent)
{
    UpdateController(controller);
    m_controller = MergeJson(m_defaults, controller);
}

void ModernOptimizerDispatcher::start()
{
    // This would be called by CurcumaMethod interface
    // For now, we redirect to static optimizeStructure method
    CurcumaLogger::warn("ModernOptimizerDispatcher::start() called - use static optimizeStructure() instead");
}

void ModernOptimizerDispatcher::printHelp() const 
{
    CurcumaLogger::header("Native Optimization Methods - Modern Computational Chemistry");
    
    CurcumaLogger::info("Available optimization algorithms:");
    auto available = getAvailableOptimizers();
    for (const auto& pair : available) {
        CurcumaLogger::info_fmt("  {} - {}", pair.first, pair.second);
    }
    
    CurcumaLogger::info("");
    CurcumaLogger::success("🧪 Scientific Algorithm Documentation:");
    
    CurcumaLogger::citation("L-BFGS: Nocedal & Wright 'Numerical Optimization' (2006), Chapter 7");
    CurcumaLogger::citation("DIIS: Pulay, P. Chem. Phys. Lett. 73, 393 (1980)");
    CurcumaLogger::citation("RFO: Banerjee et al. J. Phys. Chem. 89, 52 (1985)");
    
    CurcumaLogger::info("");
    CurcumaLogger::success("Configuration Parameters:");
    
    CurcumaLogger::param("General Parameters", "");
    CurcumaLogger::param("  verbosity", "0-3 (0=silent, 1=minimal, 2=standard, 3=verbose)");
    CurcumaLogger::param("  MaxIter", "Maximum iterations (default: 1000)");
    CurcumaLogger::param("  dE", "Energy convergence threshold in Eh (default: 1e-6)");
    CurcumaLogger::param("  GradNorm", "Gradient norm threshold in Eh/Bohr (default: 5e-4)");
    
    CurcumaLogger::param("L-BFGS Specific", "");
    CurcumaLogger::param("  memory_size", "L-BFGS memory size m (default: 10)");
    
    CurcumaLogger::param("DIIS Specific", "");
    CurcumaLogger::param("  diis_hist", "DIIS history size (default: 10)");
    CurcumaLogger::param("  diis_start", "Start DIIS after iteration (default: 5)");
    
    CurcumaLogger::param("RFO Specific", "");
    CurcumaLogger::param("  trust_radius", "Trust radius in Bohr (default: 0.05)");
    CurcumaLogger::param("  energy_threshold", "Energy threshold for line search (default: 1e-6)");
    CurcumaLogger::param("  eigenvalue_shift", "Eigenvalue shift magnitude (default: 1e-4)");
    CurcumaLogger::param("  lambda", "Legacy RFO parameter (default: 0.1)");
    
    CurcumaLogger::info("");
    CurcumaLogger::success("Usage Examples:");
    CurcumaLogger::info("  ./curcuma -opt molecule.xyz -optimizer lbfgs -v 2");
    CurcumaLogger::info("  ./curcuma -opt molecule.xyz -optimizer diis --quiet");
    CurcumaLogger::info("  ./curcuma -opt molecule.xyz -optimizer rfo -verbose");
    CurcumaLogger::info("  ./curcuma -opt molecule.xyz -optimizer auto");
}

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
/**
 * @brief Parse optimizer string to enum type with alias support - Claude Generated
 *
 * STRING TO ALGORITHM MAPPING:
 * - "lbfgspp", "lbfgs++", "external" → LBFGSPP (external library)
 * - "internal", "gpt" → INTERNAL (legacy placeholder)
 * - "native_lbfgs", "native-lbfgs", "lbfgs" → NATIVE_LBFGS (educational L-BFGS)
 * - "native_diis", "native-diis", "diis" → NATIVE_DIIS (acceleration method)
 * - "native_rfo", "native-rfo", "rfo" → NATIVE_RFO (eigenvector following)
 * - "auto", "automatic" → AUTO (system-based selection)
 *
 * DEBUGGING: Add verbosity ≥ 2 to see which algorithm was selected
 */
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

/*
 * ARCHITECTURAL DECISION RECORD: Modern Optimization Dispatcher
 *
 * CONTEXT: Multiple optimization algorithms with different capabilities and performance
 * - Native L-BFGS: Educational implementation with two-loop recursion
 * - Native DIIS: Direct Inversion of Iterative Subspace acceleration
 * - Native RFO: Rational Function Optimization with eigenvector following
 * - LBFGSpp: External high-performance library for large systems
 * - Legacy system: Backward compatibility for existing workflows
 *
 * DECISION: Strategy pattern with automatic algorithm selection
 * - Priority-based selection: user choice > system characteristics > defaults
 * - Educational focus: native algorithms directly accessible and understandable
 * - Unified interface: same SimpleOptimizationResult for all algorithms
 *
 * IMPLEMENTATION CHAIN:
 * 1. main.cpp:692 → ModernOptimizerDispatcher::optimizeStructure()
 * 2. modern_optimizer_simple.cpp:152 → parseOptimizerType() determines algorithm
 * 3. modern_optimizer_simple.cpp:170+ → algorithm-specific optimize*() methods
 * 4. lbfgs.cpp:optimize() OR external LBFGSpp library execution
 *
 * RUNTIME BEHAVIOR:
 * - "native_lbfgs" → optimizeWithNativeLBFGS() → lbfgs.cpp LBFGS class
 * - "native_diis" → optimizeWithNativeDIIS() → lbfgs.cpp DIISStep() method
 * - "native_rfo" → optimizeWithNativeRFO() → lbfgs.cpp RFOStep() method
 * - "lbfgspp" → optimizeWithLBFGSpp() → external LBFGSpp library
 * - "auto" → Heuristic: >200 atoms→LBFGSpp, else native_lbfgs
 *
 * DEBUGGING ENTRY POINTS:
 * - Set config["verbosity"] ≥ 2 to see algorithm selection logic
 * - Each optimizer reports convergence via CurcumaLogger at verbosity ≥ 1
 * - Native algorithms: Check lbfgs.cpp optimize() method for iteration details
 * - External LBFGSpp: Monitor LBFGSpp callback functions for progress
 */
SimpleOptimizationResult ModernOptimizerDispatcher::optimizeStructure(
    Molecule* molecule,
    const std::string& method_name,
    EnergyCalculator* energy_calculator,
    const json& config)
{
    // Parameter validation
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
            // We will have this as default for now
            type = OptimizerType::LBFGSPP; // External library for very large systems
            /*
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
            */
        }
        
       
        // Log optimization start
        // TODO -> Calculate initial energy 

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
            result = optimizeWithNativeLBFGS(molecule, energy_calculator, config);
            break;
        case OptimizerType::NATIVE_DIIS:
            result = optimizeWithNativeDIIS(molecule, energy_calculator, config);
            break;
        case OptimizerType::NATIVE_RFO:
            result = optimizeWithNativeRFO(molecule, energy_calculator, config);
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

// Claude Generated - Native optimization methods - 🧪 REAL COMPUTATIONAL CHEMISTRY
SimpleOptimizationResult ModernOptimizerDispatcher::optimizeWithNativeLBFGS(Molecule* molecule,
    EnergyCalculator* calc,
    const json& config)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Initialize native L-BFGS optimizer - Claude Generated
        // Handle null/empty config gracefully - Claude Generated
        json safe_config = config.is_null() ? json{} : config;
        
        LBFGS optimizer(safe_config.value("memory_size", 10)); // Memory size from config or default
        
        // Setup the optimizer with computational chemistry parameters
        optimizer.setEnergyCalculator(calc);
        optimizer.setOptimizationMethod(LBFGS::Method::LBFGS);
        int verbosity = safe_config.value("verbosity", 1); // Default: minimal output
        optimizer.setVerbosity(verbosity);
        
        // Extract geometry to optimization coordinates
        Vector initial_coords = Vector::Zero(3 * molecule->AtomCount());
        auto geometry = molecule->getGeometry();
        for (size_t i = 0; i < molecule->AtomCount(); ++i) {
            initial_coords[3*i + 0] = geometry(i, 0);
            initial_coords[3*i + 1] = geometry(i, 1);
            initial_coords[3*i + 2] = geometry(i, 2);
        }
        
        // Initialize the scientific algorithm
        optimizer.initialize(molecule->AtomCount(), initial_coords);
        
        // Extract parameters for logging and optimization - Claude Generated
        int max_iterations = safe_config.value("MaxIter", 1000);
        double energy_threshold = safe_config.value("dE", 1e-6);
        double gradient_threshold = safe_config.value("GradNorm", 5e-4);
        
        // Scientific algorithm header with verbosity control - Claude Generated
        if (verbosity >= 1) {
            CurcumaLogger::header("Native L-BFGS Optimization");
        }
        if (verbosity >= 2) {
            CurcumaLogger::citation("Algorithm: Nocedal & Wright 'Numerical Optimization' (2006), Chapter 7");
            CurcumaLogger::param("Atoms", static_cast<int>(molecule->AtomCount()));
            CurcumaLogger::param("Memory size", safe_config.value("memory_size", 10));
            CurcumaLogger::param("Max iterations", max_iterations);
            CurcumaLogger::param("Energy threshold", fmt::format("{:.2e} Eh", energy_threshold));
            CurcumaLogger::param("Gradient threshold", fmt::format("{:.2e} Eh/Bohr", gradient_threshold));
        }
        
        int iteration = 0;
        bool converged = false;
        double previous_energy = 0.0;
        
        while (iteration < max_iterations && !converged && !optimizer.hasError()) {
            Vector new_coords = optimizer.step();
            double current_energy = optimizer.getCurrentEnergy();
            Vector gradient = optimizer.getCurrentGradient();
            
            // Check convergence criteria - Educational transparency
            double energy_change = std::abs(current_energy - previous_energy);
            double gradient_norm = gradient.norm();
            
            if (iteration > 0 && energy_change < energy_threshold && gradient_norm < gradient_threshold) {
                converged = true;
                if (verbosity >= 1) {
                    CurcumaLogger::success("L-BFGS optimization converged!");
                    if (verbosity >= 2) {
                        CurcumaLogger::param("Final energy change", fmt::format("{:.2e} Eh", energy_change));
                        CurcumaLogger::param("Final gradient norm", fmt::format("{:.2e} Eh/Bohr", gradient_norm));
                    }
                }
            }
            
            // Update molecule geometry with optimized coordinates
            Eigen::Matrix<double, Eigen::Dynamic, 3> new_geometry(molecule->AtomCount(), 3);
            for (size_t i = 0; i < molecule->AtomCount(); ++i) {
                new_geometry(i, 0) = new_coords[3*i + 0];
                new_geometry(i, 1) = new_coords[3*i + 1];  
                new_geometry(i, 2) = new_coords[3*i + 2];
            }
            molecule->setGeometry(new_geometry);
            molecule->setEnergy(current_energy);
            
            previous_energy = current_energy;
            iteration++;
            
            // Progress logging with verbosity control - Claude Generated
            if (verbosity >= 2 && iteration % 10 == 0) {
                CurcumaLogger::info_fmt("Step {}: E = {:.6f} Eh, |grad| = {:.2e} Eh/Bohr", 
                                       iteration, current_energy, gradient_norm);
            }
            else if (verbosity >= 3 && iteration % 5 == 0) {
                CurcumaLogger::info_fmt("Step {}: E = {:.6f} Eh, |grad| = {:.2e} Eh/Bohr, ΔE = {:.2e} Eh", 
                                       iteration, current_energy, gradient_norm, energy_change);
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double optimization_time = duration.count() / 1000.0;
        
        if (optimizer.hasError()) {
            return SimpleOptimizationResult::failed_result("Native L-BFGS encountered numerical error", "native_lbfgs");
        }
        
        // Final status reporting with verbosity control - Claude Generated
        if (!converged && verbosity >= 1) {
            CurcumaLogger::warn_fmt("L-BFGS reached maximum iterations ({})", max_iterations);
        }
        
        if (verbosity >= 1) {
            if (converged) {
                CurcumaLogger::success("L-BFGS optimization completed successfully");
            }
            CurcumaLogger::energy_abs(optimizer.getCurrentEnergy(), "Final energy");
            CurcumaLogger::param("Iterations", iteration);
            CurcumaLogger::param("Time", fmt::format("{:.3f} seconds", optimization_time));
        }
        
        return SimpleOptimizationResult::success_result(
            *molecule, optimizer.getCurrentEnergy(), "native_lbfgs", iteration, optimization_time);
            
    } catch (const std::exception& e) {
        return SimpleOptimizationResult::failed_result(
            fmt::format("Native L-BFGS failed: {}", e.what()), "native_lbfgs");
    }
}

SimpleOptimizationResult ModernOptimizerDispatcher::optimizeWithNativeDIIS(Molecule* molecule,
    EnergyCalculator* calc,
    const json& config)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Initialize native DIIS optimizer - Claude Generated
        // Handle null/empty config gracefully - Claude Generated
        json safe_config = config.is_null() ? json{} : config;
        
        LBFGS optimizer(safe_config.value("diis_hist", 10)); // DIIS history size
        
        // Setup the optimizer with DIIS parameters
        optimizer.setEnergyCalculator(calc);
        optimizer.setOptimizationMethod(LBFGS::Method::DIIS);
        optimizer.setDIISParameters(safe_config.value("diis_hist", 10), safe_config.value("diis_start", 5));
        int verbosity = safe_config.value("verbosity", 1); // Default: minimal output
        optimizer.setVerbosity(verbosity);
        
        // Extract geometry to optimization coordinates
        Vector initial_coords = Vector::Zero(3 * molecule->AtomCount());
        auto geometry = molecule->getGeometry();
        for (size_t i = 0; i < molecule->AtomCount(); ++i) {
            initial_coords[3*i + 0] = geometry(i, 0);
            initial_coords[3*i + 1] = geometry(i, 1);
            initial_coords[3*i + 2] = geometry(i, 2);
        }
        
        // Initialize the scientific algorithm
        optimizer.initialize(molecule->AtomCount(), initial_coords);
        
        // Extract parameters for logging and optimization - Claude Generated
        int max_iterations = safe_config.value("MaxIter", 1000);
        double energy_threshold = safe_config.value("dE", 1e-6);
        double gradient_threshold = safe_config.value("GradNorm", 5e-4);
        
        // Scientific algorithm header with verbosity control - Claude Generated
        if (verbosity >= 1) {
            CurcumaLogger::header("Native DIIS Optimization");
        }
        if (verbosity >= 2) {
            CurcumaLogger::citation("Algorithm: Pulay, P. Chem. Phys. Lett. 73, 393 (1980)");
            CurcumaLogger::param("Atoms", static_cast<int>(molecule->AtomCount()));
            CurcumaLogger::param("DIIS history", safe_config.value("diis_hist", 10));
            CurcumaLogger::param("DIIS start", safe_config.value("diis_start", 5));
            CurcumaLogger::param("Max iterations", max_iterations);
            CurcumaLogger::param("Energy threshold", fmt::format("{:.2e} Eh", energy_threshold));
            CurcumaLogger::param("Gradient threshold", fmt::format("{:.2e} Eh/Bohr", gradient_threshold));
        }
        
        int iteration = 0;
        bool converged = false;
        double previous_energy = 0.0;
        
        while (iteration < max_iterations && !converged && !optimizer.hasError()) {
            Vector new_coords = optimizer.step();
            double current_energy = optimizer.getCurrentEnergy();
            Vector gradient = optimizer.getCurrentGradient();
            
            // Check convergence criteria
            double energy_change = std::abs(current_energy - previous_energy);
            double gradient_norm = gradient.norm();
            
            if (iteration > 0 && energy_change < energy_threshold && gradient_norm < gradient_threshold) {
                converged = true;
                if (verbosity >= 1) {
                    CurcumaLogger::success("DIIS optimization converged!");
                    if (verbosity >= 2) {
                        CurcumaLogger::param("Final energy change", fmt::format("{:.2e} Eh", energy_change));
                        CurcumaLogger::param("Final gradient norm", fmt::format("{:.2e} Eh/Bohr", gradient_norm));
                    }
                }
            }
            
            // Update molecule geometry
            Eigen::Matrix<double, Eigen::Dynamic, 3> new_geometry(molecule->AtomCount(), 3);
            for (size_t i = 0; i < molecule->AtomCount(); ++i) {
                new_geometry(i, 0) = new_coords[3*i + 0];
                new_geometry(i, 1) = new_coords[3*i + 1];  
                new_geometry(i, 2) = new_coords[3*i + 2];
            }
            molecule->setGeometry(new_geometry);
            molecule->setEnergy(current_energy);
            
            previous_energy = current_energy;
            iteration++;
            
            // Progress logging with method switching info - Claude Generated
            if (verbosity >= 2 && iteration % 10 == 0) {
                std::string method_info = (iteration >= safe_config.value("diis_start", 5)) ? "DIIS" : "L-BFGS";
                CurcumaLogger::info_fmt("Step {} ({}): E = {:.6f} Eh, |grad| = {:.2e} Eh/Bohr", 
                                       iteration, method_info, current_energy, gradient_norm);
            }
            else if (verbosity >= 3 && iteration % 5 == 0) {
                std::string method_info = (iteration >= safe_config.value("diis_start", 5)) ? "DIIS" : "L-BFGS";
                CurcumaLogger::info_fmt("Step {} ({}): E = {:.6f} Eh, |grad| = {:.2e} Eh/Bohr, ΔE = {:.2e} Eh", 
                                       iteration, method_info, current_energy, gradient_norm, energy_change);
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double optimization_time = duration.count() / 1000.0;
        
        if (optimizer.hasError()) {
            return SimpleOptimizationResult::failed_result("Native DIIS encountered numerical error", "native_diis");
        }
        
        // Final status reporting with verbosity control - Claude Generated
        if (verbosity >= 1) {
            if (converged) {
                CurcumaLogger::success("DIIS optimization completed successfully");
            }
            CurcumaLogger::energy_abs(optimizer.getCurrentEnergy(), "Final energy");
            CurcumaLogger::param("Iterations", iteration);
            CurcumaLogger::param("Time", fmt::format("{:.3f} seconds", optimization_time));
        }
        
        return SimpleOptimizationResult::success_result(
            *molecule, optimizer.getCurrentEnergy(), "native_diis", iteration, optimization_time);
            
    } catch (const std::exception& e) {
        return SimpleOptimizationResult::failed_result(
            fmt::format("Native DIIS failed: {}", e.what()), "native_diis");
    }
}

SimpleOptimizationResult ModernOptimizerDispatcher::optimizeWithNativeRFO(Molecule* molecule,
    EnergyCalculator* calc,
    const json& config)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    
    try {
        // Initialize native RFO optimizer - Claude Generated
        // Handle null/empty config gracefully - Claude Generated
        json safe_config = config.is_null() ? json{} : config;
        
        LBFGS optimizer(10);
        
        // Setup the optimizer with RFO parameters
        optimizer.setEnergyCalculator(calc);
        optimizer.setOptimizationMethod(LBFGS::Method::RFO);
        optimizer.setLambda(safe_config.value("lambda", 0.1)); // RFO lambda parameter (legacy)
        int verbosity = safe_config.value("verbosity", 1); // Default: minimal output
        optimizer.setVerbosity(verbosity);
        
        // Configure RFO-specific parameters from config
        double trust_radius = safe_config.value("trust_radius", 0.05);  // Conservative default
        double energy_threshold = safe_config.value("energy_threshold", 1e-6);  // Stricter threshold
        double eigenvalue_shift = safe_config.value("eigenvalue_shift", 1e-4);  // Smaller shift
        optimizer.setRFOParameters(trust_radius, energy_threshold, eigenvalue_shift);
        
        // Set up masses for mass-weighted coordinates (important for RFO)
        std::vector<double> masses(molecule->AtomCount());
        for (size_t i = 0; i < molecule->AtomCount(); ++i) {
            // Use atomic masses from element types - simplified for now
            int element = molecule->Atoms()[i];
            masses[i] = (element == 1) ? 1.008 : (element == 6) ? 12.011 : (element == 8) ? 15.999 : 12.0; // H, C, O, default
        }
        optimizer.setMasses(masses);
        
        // Extract geometry to optimization coordinates
        Vector initial_coords = Vector::Zero(3 * molecule->AtomCount());
        auto geometry = molecule->getGeometry();
        for (size_t i = 0; i < molecule->AtomCount(); ++i) {
            initial_coords[3*i + 0] = geometry(i, 0);
            initial_coords[3*i + 1] = geometry(i, 1);
            initial_coords[3*i + 2] = geometry(i, 2);
        }
        
        // Initialize the scientific algorithm
        optimizer.initialize(molecule->AtomCount(), initial_coords);
        
        // 🧪 CRITICAL: Initialize Hessian properly for RFO - Claude Generated
        // RFO NEEDS a reasonable Hessian approximation to work properly!
        int coord_size = 3 * molecule->AtomCount();
        Matrix initial_hessian = Matrix::Zero(coord_size, coord_size);
        
        // Use simple diagonal Hessian based on typical bond force constants
        for (int i = 0; i < coord_size; ++i) {
            int atom_index = i / 3;
            int element = molecule->Atoms()[atom_index];
            
            // Typical force constants in atomic units for different elements
            double force_constant = 0.5;  // Default
            if (element == 1)      force_constant = 0.3;  // H
            else if (element == 6) force_constant = 0.5;  // C
            else if (element == 7) force_constant = 0.6;  // N
            else if (element == 8) force_constant = 0.7;  // O
            
            initial_hessian(i, i) = force_constant;
        }
        
        optimizer.setHessian(initial_hessian);
        
        // Extract parameters for logging and optimization - Claude Generated
        int max_iterations = safe_config.value("MaxIter", 500); // RFO typically needs fewer iterations
        double convergence_energy_threshold = safe_config.value("dE", 1e-6);  // For convergence check
        double gradient_threshold = safe_config.value("GradNorm", 5e-4);
        
        // Scientific algorithm header with verbosity control - Claude Generated
        if (verbosity >= 1) {
            CurcumaLogger::header("Native RFO Optimization");
        }
        if (verbosity >= 2) {
            CurcumaLogger::citation("Algorithm: Banerjee et al. J. Phys. Chem. 89, 52 (1985)");
            CurcumaLogger::param("Atoms", static_cast<int>(molecule->AtomCount()));
            CurcumaLogger::param("Trust radius", fmt::format("{:.3f} Bohr", trust_radius));
            CurcumaLogger::param("Energy threshold", fmt::format("{:.2e} Eh", energy_threshold));
            CurcumaLogger::param("Eigenvalue shift", fmt::format("{:.2e}", eigenvalue_shift));
            CurcumaLogger::param("Max iterations", max_iterations);
            CurcumaLogger::param("Convergence energy threshold", fmt::format("{:.2e} Eh", convergence_energy_threshold));
            CurcumaLogger::param("Gradient threshold", fmt::format("{:.2e} Eh/Bohr", gradient_threshold));
        }
        
        int iteration = 0;
        bool converged = false;
        double previous_energy = 0.0;
        
        while (iteration < max_iterations && !converged && !optimizer.hasError()) {
            Vector new_coords = optimizer.step();
            double current_energy = optimizer.getCurrentEnergy();
            Vector gradient = optimizer.getCurrentGradient();
            
            // Check convergence criteria
            double energy_change = std::abs(current_energy - previous_energy);
            double gradient_norm = gradient.norm();
            
            if (iteration > 0 && energy_change < convergence_energy_threshold && gradient_norm < gradient_threshold) {
                converged = true;
                if (verbosity >= 1) {
                    CurcumaLogger::success("RFO optimization converged!");
                    if (verbosity >= 2) {
                        CurcumaLogger::param("Final energy change", fmt::format("{:.2e} Eh", energy_change));
                        CurcumaLogger::param("Final gradient norm", fmt::format("{:.2e} Eh/Bohr", gradient_norm));
                    }
                }
            }
            
            // Update molecule geometry
            Eigen::Matrix<double, Eigen::Dynamic, 3> new_geometry(molecule->AtomCount(), 3);
            for (size_t i = 0; i < molecule->AtomCount(); ++i) {
                new_geometry(i, 0) = new_coords[3*i + 0];
                new_geometry(i, 1) = new_coords[3*i + 1];  
                new_geometry(i, 2) = new_coords[3*i + 2];
            }
            molecule->setGeometry(new_geometry);
            molecule->setEnergy(current_energy);
            
            previous_energy = current_energy;
            iteration++;
            
            // Progress logging with RFO specific info - Claude Generated
            if (verbosity >= 2 && iteration % 10 == 0) {
                CurcumaLogger::info_fmt("Step {} (RFO): E = {:.6f} Eh, |grad| = {:.2e} Eh/Bohr", 
                                       iteration, current_energy, gradient_norm);
            }
            else if (verbosity >= 3 && iteration % 5 == 0) {
                CurcumaLogger::info_fmt("Step {} (RFO): E = {:.6f} Eh, |grad| = {:.2e} Eh/Bohr, ΔE = {:.2e} Eh", 
                                       iteration, current_energy, gradient_norm, energy_change);
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double optimization_time = duration.count() / 1000.0;
        
        if (optimizer.hasError()) {
            return SimpleOptimizationResult::failed_result("Native RFO encountered numerical error", "native_rfo");
        }
        
        // Final status reporting with verbosity control - Claude Generated
        if (verbosity >= 1) {
            if (converged) {
                CurcumaLogger::success("RFO optimization completed successfully");
            }
            CurcumaLogger::energy_abs(optimizer.getCurrentEnergy(), "Final energy");
            CurcumaLogger::param("Iterations", iteration);
            CurcumaLogger::param("Time", fmt::format("{:.3f} seconds", optimization_time));
        }
        
        return SimpleOptimizationResult::success_result(
            *molecule, optimizer.getCurrentEnergy(), "native_rfo", iteration, optimization_time);
            
    } catch (const std::exception& e) {
        return SimpleOptimizationResult::failed_result(
            fmt::format("Native RFO failed: {}", e.what()), "native_rfo");
    }
}

} // namespace ModernOptimization