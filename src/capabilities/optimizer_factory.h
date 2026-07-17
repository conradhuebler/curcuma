/*
 * <Optimizer Factory - Strategy Pattern Factory>
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

#include "optimizer_driver.h"
#include <functional>
#include <map>
#include <memory>

namespace Optimization {

/**
 * @brief Factory for creating optimizer strategies - Claude Generated
 * Analog to QM system's method factory pattern
 */
class OptimizerFactory {
public:
    /**
     * @brief Create optimizer by type enum
     */
    static std::unique_ptr<OptimizerDriver> createOptimizer(
        OptimizerType type,
        EnergyCalculator* energy_calculator = nullptr);

    /**
     * @brief Create optimizer by string name
     */
    static std::unique_ptr<OptimizerDriver> createOptimizer(
        const std::string& method_name,
        EnergyCalculator* energy_calculator = nullptr);

    /**
     * @brief Create optimizer from JSON configuration
     */
    static std::unique_ptr<OptimizerDriver> createOptimizer(
        const json& config,
        EnergyCalculator* energy_calculator = nullptr);

    /**
     * @brief Get available optimizer types and descriptions
     */
    static std::map<OptimizerType, std::string> getAvailableOptimizers();

    /**
     * @brief Get optimizer description by type
     */
    static std::string getOptimizerDescription(OptimizerType type);

    /**
     * @brief Check if optimizer type is available
     */
    static bool isAvailable(OptimizerType type);

    /**
     * @brief Auto-select best optimizer for system
     */
    static OptimizerType selectOptimalOptimizer(const curcuma::Molecule& molecule,
        const json& preferences = json{});

private:
    // Factory function type
    using OptimizerCreator = std::function<std::unique_ptr<OptimizerDriver>()>;

    // Factory registry
    static const std::map<OptimizerType, OptimizerCreator> s_factory_registry;
    static const std::map<OptimizerType, std::string> s_optimizer_descriptions;

    // Individual factory methods
    static std::unique_ptr<OptimizerDriver> createLBFGSpp();
    static std::unique_ptr<OptimizerDriver> createANCOpt();
    static std::unique_ptr<OptimizerDriver> createNativeLBFGS();
    static std::unique_ptr<OptimizerDriver> createNativeDIIS();
    static std::unique_ptr<OptimizerDriver> createNativeRFO();
};

/**
 * @brief High-level optimization dispatcher - Claude Generated
 * Analog to EnergyCalculator routing for optimization methods
 */
class OptimizationDispatcher {
public:
    /**
     * @brief Optimize molecule with specified method
     */
    static OptimizationResult optimizeStructure(
        curcuma::Molecule* molecule,
        OptimizerType optimizer_type,
        EnergyCalculator* energy_calculator,
        const json& config = json{});

    /**
     * @brief Optimize molecule with string method name
     */
    static OptimizationResult optimizeStructure(
        curcuma::Molecule* molecule,
        const std::string& method_name,
        EnergyCalculator* energy_calculator,
        const json& config = json{});

    /**
     * @brief Optimize multiple molecules (batch processing)
     * @param threads Number of parallel molecule-level workers. Values <= 1 keep
     *        the original sequential path; values > 1 dispatch frames to a
     *        CxxThreadPool with progress bar.
     * @param energy_controller Full controller used to construct per-worker
     *        EnergyCalculators in parallel mode. When empty, the supplied
     *        energy_calculator is used directly (sequential path only).
     */
    static std::vector<OptimizationResult> optimizeBatch(
        const std::vector<curcuma::Molecule>& molecules,
        OptimizerType optimizer_type,
        EnergyCalculator* energy_calculator,
        const json& config = json{},
        int threads = 1,
        const json& energy_controller = json{});

    /**
     * @brief Auto-optimize with best method selection
     */
    static OptimizationResult autoOptimize(
        curcuma::Molecule* molecule,
        EnergyCalculator* energy_calculator,
        const json& preferences = json{});

private:
    // Helper methods
    static json mergeConfigurations(const json& base_config, const json& override_config);
    static void validateConfiguration(const json& config, OptimizerType type);
};

} // namespace Optimization