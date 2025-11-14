/*
 * <Optimizer Factory Implementation>
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

#include "optimizer_factory.h"
#include "lbfgspp_optimizer.h"
#include "ancopt_optimizer.h"
// Future includes:
// #include "internal_lbfgs_optimizer.h"
// #include "diis_optimizer.h"
// #include "rfo_optimizer.h"

namespace Optimization {

// Claude Generated - Factory registry initialization
const std::map<OptimizerType, OptimizerFactory::OptimizerFactoryFunction>
    OptimizerFactory::s_factory_registry = {
        { OptimizerType::LBFGSPP, &OptimizerFactory::createLBFGSpp },
        { OptimizerType::ANCOPT, &OptimizerFactory::createANCOpt },
        // Future implementations:
        // { OptimizerType::INTERNAL_LBFGS, &OptimizerFactory::createInternalLBFGS },
        // { OptimizerType::DIIS, &OptimizerFactory::createDIIS },
        // { OptimizerType::RFO, &OptimizerFactory::createRFO }
    };

const std::map<OptimizerType, std::string> OptimizerFactory::s_optimizer_descriptions = {
    { OptimizerType::LBFGSPP, "External LBFGSpp library - robust L-BFGS implementation" },
    { OptimizerType::ANCOPT, "Approximate Normal Coordinate Optimizer (from XTB by Stefan Grimme)" },
    { OptimizerType::INTERNAL_LBFGS, "Internal LBFGS implementation with custom features" },
    { OptimizerType::DIIS, "Direct Inversion of Iterative Subspace method" },
    { OptimizerType::RFO, "Rational Function Optimization with Hessian updates" },
    { OptimizerType::AUTO, "Automatic selection based on system characteristics" }
};

// Claude Generated - Factory methods
std::unique_ptr<OptimizerInterface> OptimizerFactory::createOptimizer(
    OptimizerType type, EnergyCalculator* energy_calculator)
{

    // Handle AUTO selection
    if (type == OptimizerType::AUTO) {
        // Default to LBFGSpp for now - could be made smarter
        type = OptimizerType::LBFGSPP;
        CurcumaLogger::info("Auto-selected optimizer: LBFGSpp");
    }

    // Find factory function
    auto it = s_factory_registry.find(type);
    if (it == s_factory_registry.end()) {
        throw std::invalid_argument("Optimizer type not implemented: " + optimizerTypeToString(type));
    }

    // Create optimizer
    auto optimizer = it->second();

    // Set energy calculator if provided
    if (energy_calculator && optimizer) {
        optimizer->LoadConfiguration(optimizer->GetDefaultConfiguration());
        // Note: Energy calculator will be set during initialization
    }

    return optimizer;
}

std::unique_ptr<OptimizerInterface> OptimizerFactory::createOptimizer(
    const std::string& method_name, EnergyCalculator* energy_calculator)
{

    OptimizerType type = parseOptimizerType(method_name);
    return createOptimizer(type, energy_calculator);
}

std::unique_ptr<OptimizerInterface> OptimizerFactory::createOptimizer(
    const json& config, EnergyCalculator* energy_calculator)
{

    // Parse optimizer type from config
    OptimizerType type = OptimizerType::LBFGSPP; // Default

    if (config.contains("method")) {
        type = parseOptimizerType(config["method"]);
    } else if (config.contains("optimizer_type")) {
        type = parseOptimizerType(config["optimizer_type"]);
    } else if (config.contains("optimethod")) {
        // Legacy compatibility
        int legacy_method = config["optimethod"];
        type = (legacy_method == 0) ? OptimizerType::LBFGSPP : OptimizerType::INTERNAL_LBFGS;
    }

    auto optimizer = createOptimizer(type, energy_calculator);

    if (optimizer) {
        optimizer->LoadConfiguration(config);
    }

    return optimizer;
}

std::map<OptimizerType, std::string> OptimizerFactory::getAvailableOptimizers()
{
    std::map<OptimizerType, std::string> available;

    for (const auto& pair : s_optimizer_descriptions) {
        if (pair.first == OptimizerType::AUTO || isAvailable(pair.first)) {
            available[pair.first] = pair.second;
        }
    }

    return available;
}

std::string OptimizerFactory::getOptimizerDescription(OptimizerType type)
{
    auto it = s_optimizer_descriptions.find(type);
    return (it != s_optimizer_descriptions.end()) ? it->second : "Unknown optimizer";
}

bool OptimizerFactory::isAvailable(OptimizerType type)
{
    return s_factory_registry.find(type) != s_factory_registry.end();
}

OptimizerType OptimizerFactory::selectOptimalOptimizer(const Molecule& molecule, const json& preferences)
{
    // Simple heuristics for optimizer selection
    int atom_count = molecule.AtomCount();

    // Check user preferences
    if (preferences.contains("prefer_robust")) {
        bool prefer_robust = preferences["prefer_robust"];
        if (prefer_robust) {
            return OptimizerType::LBFGSPP; // Most robust option
        }
    }

    if (preferences.contains("prefer_speed")) {
        bool prefer_speed = preferences["prefer_speed"];
        if (prefer_speed && atom_count < 100) {
            return OptimizerType::INTERNAL_LBFGS; // Faster for small systems
        }
    }

    // Default selection based on system size
    if (atom_count > 500) {
        return OptimizerType::LBFGSPP; // Better for large systems
    } else if (atom_count > 50) {
        return OptimizerType::LBFGSPP; // Good general purpose
    } else {
        return OptimizerType::LBFGSPP; // Default choice
    }
}

// Claude Generated - Individual factory methods
std::unique_ptr<OptimizerInterface> OptimizerFactory::createLBFGSpp()
{
    return std::make_unique<LBFGSppOptimizer>();
}

std::unique_ptr<OptimizerInterface> OptimizerFactory::createANCOpt()
{
    return std::make_unique<ANCOptimizer>();
}

std::unique_ptr<OptimizerInterface> OptimizerFactory::createInternalLBFGS()
{
    throw std::runtime_error("Internal LBFGS optimizer not yet implemented");
    // return std::make_unique<InternalLBFGSOptimizer>();
}

std::unique_ptr<OptimizerInterface> OptimizerFactory::createDIIS()
{
    throw std::runtime_error("DIIS optimizer not yet implemented");
    // return std::make_unique<DIISOptimizer>();
}

std::unique_ptr<OptimizerInterface> OptimizerFactory::createRFO()
{
    throw std::runtime_error("RFO optimizer not yet implemented");
    // return std::make_unique<RFOOptimizer>();
}

// Claude Generated - OptimizationDispatcher implementation
OptimizationResult OptimizationDispatcher::optimizeStructure(
    Molecule* molecule,
    OptimizerType optimizer_type,
    EnergyCalculator* energy_calculator,
    const json& config)
{

    try {
        // Validate inputs
        if (!molecule) {
            return OptimizationResult::failed_result("Null molecule provided");
        }
        if (!energy_calculator) {
            return OptimizationResult::failed_result("Null energy calculator provided");
        }

        // Create optimizer
        auto optimizer = OptimizerFactory::createOptimizer(optimizer_type, energy_calculator);
        if (!optimizer) {
            return OptimizationResult::failed_result("Failed to create optimizer");
        }

        // Configure optimizer
        json merged_config = mergeConfigurations(optimizer->GetDefaultConfiguration(), config);
        optimizer->LoadConfiguration(merged_config);
        // Note: Energy calculator is passed during createOptimizer and set during initialization

        // Initialize optimization
        if (!optimizer->InitializeOptimization(*molecule)) {
            return OptimizationResult::failed_result("Optimizer initialization failed");
        }

        // Extract settings from config
        bool write_trajectory = merged_config.value("write_trajectory", true);
        bool verbose = merged_config.value("verbose", false);

        // Perform optimization
        OptimizationResult result = optimizer->Optimize(write_trajectory, verbose);

        // Update molecule with final structure
        if (result.success) {
            *molecule = result.final_molecule;
        }

        return result;

    } catch (const std::exception& e) {
        return OptimizationResult::failed_result(
            fmt::format("Optimization failed: {}", e.what()));
    }
}

OptimizationResult OptimizationDispatcher::optimizeStructure(
    Molecule* molecule,
    const std::string& method_name,
    EnergyCalculator* energy_calculator,
    const json& config)
{

    try {
        OptimizerType type = parseOptimizerType(method_name);
        return optimizeStructure(molecule, type, energy_calculator, config);
    } catch (const std::exception& e) {
        return OptimizationResult::failed_result(e.what());
    }
}

std::vector<OptimizationResult> OptimizationDispatcher::optimizeBatch(
    const std::vector<Molecule>& molecules,
    OptimizerType optimizer_type,
    EnergyCalculator* energy_calculator,
    const json& config)
{

    std::vector<OptimizationResult> results;
    results.reserve(molecules.size());

    CurcumaLogger::info_fmt("Starting batch optimization of {} structures", molecules.size());

    for (size_t i = 0; i < molecules.size(); ++i) {
        CurcumaLogger::info_fmt("Optimizing structure {} of {}", i + 1, molecules.size());

        Molecule mol_copy = molecules[i];
        OptimizationResult result = optimizeStructure(&mol_copy, optimizer_type,
            energy_calculator, config);
        results.push_back(result);

        if (!result.success) {
            CurcumaLogger::warn_fmt("Structure {} optimization failed: {}",
                i + 1, result.error_message);
        }
    }

    int successful = std::count_if(results.begin(), results.end(),
        [](const OptimizationResult& r) { return r.success; });

    CurcumaLogger::success_fmt("Batch optimization completed: {}/{} successful",
        successful, molecules.size());

    return results;
}

OptimizationResult OptimizationDispatcher::autoOptimize(
    Molecule* molecule,
    EnergyCalculator* energy_calculator,
    const json& preferences)
{

    OptimizerType optimal_type = OptimizerFactory::selectOptimalOptimizer(*molecule, preferences);

    CurcumaLogger::info_fmt("Auto-selected optimizer: {} for {}-atom system",
        optimizerTypeToString(optimal_type), molecule->AtomCount());

    return optimizeStructure(molecule, optimal_type, energy_calculator, preferences);
}

// Claude Generated - Helper methods
json OptimizationDispatcher::mergeConfigurations(const json& base_config, const json& override_config)
{
    json merged = base_config;

    for (auto it = override_config.begin(); it != override_config.end(); ++it) {
        merged[it.key()] = it.value();
    }

    return merged;
}

void OptimizationDispatcher::validateConfiguration(const json& config, OptimizerType type)
{
    // Basic validation - can be extended
    if (config.contains("max_iterations") && config["max_iterations"] <= 0) {
        throw std::invalid_argument("max_iterations must be positive");
    }

    if (config.contains("energy_threshold") && config["energy_threshold"] <= 0) {
        throw std::invalid_argument("energy_threshold must be positive");
    }

    if (config.contains("gradient_threshold") && config["gradient_threshold"] <= 0) {
        throw std::invalid_argument("gradient_threshold must be positive");
    }

    // Type-specific validation could be added here
}

} // namespace Optimization