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
#include "ancopt_optimizer.h"
#include "lbfgspp_optimizer.h"
#include "native_optimizer_adapters.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/intra_parallel_context.h"
#include "src/core/molecule.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

namespace Optimization {

// Claude Generated (Apr 2026) - Factory registry with all optimizer types
const std::map<OptimizerType, OptimizerFactory::OptimizerCreator>
    OptimizerFactory::s_factory_registry = {
        { OptimizerType::LBFGSPP, &OptimizerFactory::createLBFGSpp },
        { OptimizerType::ANCOPT, &OptimizerFactory::createANCOpt },
        { OptimizerType::NATIVE_LBFGS, &OptimizerFactory::createNativeLBFGS },
        { OptimizerType::NATIVE_DIIS, &OptimizerFactory::createNativeDIIS },
        { OptimizerType::NATIVE_RFO, &OptimizerFactory::createNativeRFO },
    };

const std::map<OptimizerType, std::string> OptimizerFactory::s_optimizer_descriptions = {
    { OptimizerType::LBFGSPP, "External LBFGSpp library - robust L-BFGS" },
    { OptimizerType::NATIVE_LBFGS, "Native L-BFGS (Nocedal & Wright)" },
    { OptimizerType::NATIVE_DIIS, "Native DIIS (Pulay 1980)" },
    { OptimizerType::NATIVE_RFO, "Native RFO (Banerjee et al. 1985)" },
    { OptimizerType::ANCOPT, "ANCOpt (Grimme, approximate normal coordinates)" },
    { OptimizerType::AUTO, "Automatic selection based on system size" }
};

// Claude Generated - Factory methods
std::unique_ptr<OptimizerDriver> OptimizerFactory::createOptimizer(
    OptimizerType type, EnergyCalculator* energy_calculator)
{

    // Handle AUTO selection
    if (type == OptimizerType::AUTO) {
        // Default to LBFGSpp for now - could be made smarter
        type = OptimizerType::LBFGSPP;
        CurcumaLogger::info("Auto-selected optimizer: LBFGSpp");
    }

    // Find factory function — fall back to LBFGSpp for types not yet in factory (native_* etc.)
    auto it = s_factory_registry.find(type);
    if (it == s_factory_registry.end()) {
        CurcumaLogger::warn_fmt("Optimizer '{}' not in factory, using LBFGSpp", optimizerTypeToString(type));
        it = s_factory_registry.find(OptimizerType::LBFGSPP);
    }

    // Create optimizer
    auto optimizer = it->second();

    // Set energy calculator if provided (Claude Generated - critical for optimization)
    if (energy_calculator && optimizer) {
        optimizer->setEnergyCalculator(energy_calculator);
        optimizer->LoadConfiguration(optimizer->GetDefaultConfiguration());
    }

    return optimizer;
}

std::unique_ptr<OptimizerDriver> OptimizerFactory::createOptimizer(
    const std::string& method_name, EnergyCalculator* energy_calculator)
{

    OptimizerType type = parseOptimizerType(method_name);
    return createOptimizer(type, energy_calculator);
}

std::unique_ptr<OptimizerDriver> OptimizerFactory::createOptimizer(
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
        type = (legacy_method == 0) ? OptimizerType::LBFGSPP : OptimizerType::NATIVE_LBFGS;
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
            return OptimizerType::NATIVE_LBFGS; // Faster for small systems
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
std::unique_ptr<OptimizerDriver> OptimizerFactory::createLBFGSpp()
{
    return std::make_unique<LBFGSppOptimizer>();
}

std::unique_ptr<OptimizerDriver> OptimizerFactory::createANCOpt()
{
    return std::make_unique<ANCOptimizer>();
}

std::unique_ptr<OptimizerDriver> OptimizerFactory::createNativeLBFGS()
{
    return std::make_unique<NativeLBFGSAdapter>();
}

std::unique_ptr<OptimizerDriver> OptimizerFactory::createNativeDIIS()
{
    return std::make_unique<NativeDIISAdapter>();
}

std::unique_ptr<OptimizerDriver> OptimizerFactory::createNativeRFO()
{
    return std::make_unique<NativeRFOAdapter>();
}

// Claude Generated (May 2026): Build preset config from user-provided preset string
static json buildPresetConfig(const json& config)
{
    if (!config.contains("convergence_preset"))
        return json::object();

    std::string preset = config["convergence_preset"].get<std::string>();
    json preset_config;
    preset_config["convergence_preset"] = preset;
    if (preset == "loose") {
        preset_config["energy_threshold"] = 1.0;
        preset_config["rmsd_threshold"] = 0.05;
        preset_config["gradient_threshold"] = 1e-3;
        preset_config["max_iterations"] = 1000;
    } else if (preset == "normal") {
        preset_config["energy_threshold"] = 0.1;
        preset_config["rmsd_threshold"] = 0.01;
        preset_config["gradient_threshold"] = 5e-4;
        preset_config["max_iterations"] = 5000;
    } else if (preset == "tight") {
        preset_config["energy_threshold"] = 1e-6 * 2625.5; // 1e-6 Eh in kJ/mol
        preset_config["rmsd_threshold"] = 1e-3;
        preset_config["gradient_threshold"] = 1e-5;
        preset_config["max_iterations"] = 10000;
    } else if (preset == "verytight") {
        preset_config["energy_threshold"] = 1e-7 * 2625.5; // 1e-7 Eh in kJ/mol
        preset_config["rmsd_threshold"] = 1e-4;
        preset_config["gradient_threshold"] = 1e-6;
        preset_config["max_iterations"] = 20000;
    }
    return preset_config;
}

// Claude Generated - OptimizationDispatcher implementation
OptimizationResult OptimizationDispatcher::optimizeStructure(
    curcuma::Molecule* molecule,
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
        // Claude Generated (May 2026): Apply convergence preset as intermediate layer:
        // defaults → preset → user-specified individual parameters (only those differing from defaults)
        json preset_config = buildPresetConfig(config);
        json merged_config = mergeConfigurations(optimizer->GetDefaultConfiguration(), preset_config);

        // Override preset only with explicitly user-provided individual parameters.
        // Defaults already in config must not wipe out the preset. We detect explicit
        // user overrides by checking if the value differs from the default.
        json default_config = optimizer->GetDefaultConfiguration();
        for (auto it = config.begin(); it != config.end(); ++it) {
            const std::string& key = it.key();
            if (key == "convergence_preset")
                continue; // Preset itself is already applied above
            if (!default_config.contains(key)) {
                merged_config[key] = it.value(); // New parameter not in defaults
            } else if (default_config[key] != it.value()) {
                merged_config[key] = it.value(); // Explicitly changed from default
            }
        }
        optimizer->LoadConfiguration(merged_config);
        optimizer->setEnergyCalculator(energy_calculator);

        // Initialize optimization
        if (!optimizer->InitializeOptimization(*molecule)) {
            return OptimizationResult::failed_result("Optimizer initialization failed");
        }

        // Extract settings from config
        bool write_trajectory = merged_config.value("write_trajectory", true);
        int verbosity = merged_config.contains("verbosity") ? merged_config["verbosity"].get<int>() : 1;

        // Perform optimization
        OptimizationResult result = optimizer->Optimize(write_trajectory, verbosity);

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
    curcuma::Molecule* molecule,
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

// Claude Generated (Jun 2026): Worker that optimises one frame of a batch.
// Each worker owns its own EnergyCalculator so that the stateful GFN-FF/
// xTB backends are not shared across threads. The intra-molecule parallelism
// is suppressed because the CxxThreadPool already owns molecule-level
// concurrency.
class OptimizationBatchThread : public CxxThread {
public:
    OptimizationBatchThread(
        size_t index,
        const curcuma::Molecule& molecule,
        OptimizerType optimizer_type,
        const std::string& method,
        const json& energy_controller,
        const json& config)
        : m_index(index)
        , m_molecule(molecule)
        , m_optimizer_type(optimizer_type)
        , m_method(method)
        , m_energy_controller(energy_controller)
        , m_config(config)
    {
        // Suppress per-worker output during the parallel phase so that step tables
        // from concurrent optimizations do not interleave on stdout. The main thread
        // prints an ordered summary once all workers are finished. Trajectory files
        // are also disabled because all workers would otherwise write to the same
        // basename simultaneously.
        m_config["verbosity"] = 0;
        m_config["write_trajectory"] = false;
    }

    int execute() override
    {
        curcuma::SuppressIntraParallel intra_guard;

        EnergyCalculator energy_calc(m_method, m_energy_controller);
        energy_calc.setIterativeMode(true);
        if (m_method == "gfn1" || m_method == "gfn2")
            energy_calc.setWarmStart(true);

        m_result = OptimizationDispatcher::optimizeStructure(
            &m_molecule, m_optimizer_type, &energy_calc, m_config);

        // optimizeStructure only updates the molecule on success; keep the
        // final geometry available for the dispatcher in every case.
        if (m_result.final_molecule.AtomCount() == 0 && m_molecule.AtomCount() > 0)
            m_result.final_molecule = m_molecule;

        return 0;
    }

    size_t index() const { return m_index; }
    const OptimizationResult& result() const { return m_result; }

private:
    size_t m_index;
    curcuma::Molecule m_molecule;
    OptimizerType m_optimizer_type;
    std::string m_method;
    json m_energy_controller;
    json m_config;
    OptimizationResult m_result;
};

std::vector<OptimizationResult> OptimizationDispatcher::optimizeBatch(
    const std::vector<curcuma::Molecule>& molecules,
    OptimizerType optimizer_type,
    EnergyCalculator* energy_calculator,
    const json& config,
    int threads,
    const json& energy_controller)
{

    std::vector<OptimizationResult> results;
    results.resize(molecules.size());

    if (threads <= 1 || molecules.size() <= 1) {
        // Sequential path: no thread pool, no progress bar, identical behaviour
        // to the original implementation.
        CurcumaLogger::info_fmt("Starting batch optimization of {} structures", molecules.size());

        for (size_t i = 0; i < molecules.size(); ++i) {
            curcuma::Molecule mol_copy = molecules[i];
            OptimizationResult result = optimizeStructure(
                &mol_copy, optimizer_type, energy_calculator, config);
            results[i] = result;
        }

        int successful = std::count_if(results.begin(), results.end(),
            [](const OptimizationResult& r) { return r.success; });

        CurcumaLogger::success_fmt("Batch optimization completed: {}/{} successful",
            successful, molecules.size());

        CurcumaLogger::info("Batch optimization summary:");
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& r = results[i];
            std::string status = r.success ? "converged" : "not converged";
            CurcumaLogger::result_fmt("  Frame {:3}: {} after {:4} iterations, final energy = {:.8f} Eh",
                i + 1, status, r.iterations_performed, r.final_energy);
            if (!r.success && !r.error_message.empty())
                CurcumaLogger::warn_fmt("    Reason: {}", r.error_message);
        }

        return results;
    }

    // Parallel path: one CxxThreadPool worker per molecule. Each worker builds
    // its own EnergyCalculator so the stateful backends stay thread-local.
    CurcumaLogger::info_fmt("Starting parallel batch optimization of {} structures using {} threads",
        molecules.size(), threads);

    const int worker_count = static_cast<int>(std::min(
        static_cast<size_t>(threads), molecules.size()));

    CxxThreadPool pool;
    pool.setProgressBar(CxxThreadPool::ProgressBarType::Continously);
    pool.setActiveThreadCount(worker_count);

    // Recover the method and full energy controller needed to recreate an
    // EnergyCalculator per worker. The supplied calculator is only used as a
    // template; it is not shared across threads.
    // per_worker_controller must be resolved first because method lives in
    // energy_controller (the full top-level controller), not in config
    // (the optimizer-only JSON that lacks the "method" key).
    json per_worker_controller = energy_controller.empty() ? config : energy_controller;
    std::string method = per_worker_controller.value("method", std::string("gfnff"));

    std::vector<OptimizationBatchThread*> workers;
    workers.reserve(molecules.size());
    for (size_t i = 0; i < molecules.size(); ++i) {
        auto* th = new OptimizationBatchThread(
            i, molecules[i], optimizer_type, method, per_worker_controller, config);
        workers.push_back(th);
        pool.addThread(th);
    }

    // Suppress all CurcumaLogger output during the parallel phase. Because the
    // logger verbosity is a single global static, concurrent workers would
    // otherwise produce interleaved initialization / step messages. The
    // original level is restored immediately after the pool finishes so the
    // ordered final summary can be printed at the user's requested verbosity.
    int saved_global_verbosity = CurcumaLogger::get_verbosity();
    CurcumaLogger::set_verbosity(0);

    pool.StartAndWait();

    CurcumaLogger::set_verbosity(saved_global_verbosity);

    // Collect results in original molecular order, not pool finish order.
    // The CxxThreadPool owns the workers (auto-delete enabled by default),
    // so we must not delete them ourselves.
    int successful = 0;
    for (auto* th : workers) {
        results[th->index()] = th->result();
        if (th->result().success)
            ++successful;
    }

    CurcumaLogger::success_fmt("Batch optimization completed: {}/{} successful",
        successful, molecules.size());

    // Print an input-order summary so the user gets clean, sequential output
    // instead of interleaved per-worker step tables.
    CurcumaLogger::info("Batch optimization summary:");
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        std::string status = r.success ? "converged" : "not converged";
        CurcumaLogger::result_fmt("  Frame {:3}: {} after {:4} iterations, final energy = {:.8f} Eh",
            i + 1, status, r.iterations_performed, r.final_energy);
        if (!r.success && !r.error_message.empty()) {
            CurcumaLogger::warn_fmt("    Reason: {}", r.error_message);
        }
    }

    return results;
}

OptimizationResult OptimizationDispatcher::autoOptimize(
    curcuma::Molecule* molecule,
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