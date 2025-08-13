/*
 * <Native LBFGS Optimizer Adapter - Claude Generated>
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
 */

#pragma once

#include "modern_optimizer_simple.h"
#include "optimiser/lbfgs.h"
#include "src/core/curcuma_logger.h"
#include <chrono>

namespace ModernOptimization {

/**
 * @brief Native LBFGS Optimizer Adapter - Claude Generated
 * Adapts the cleaned native LBFGS implementation to our modern optimizer interface
 */
class NativeLBFGSOptimizer {
public:
    explicit NativeLBFGSOptimizer(int memory_size = 10, LBFGS::Method method = LBFGS::Method::LBFGS);
    ~NativeLBFGSOptimizer() = default;

    /**
     * @brief Optimize a molecule using the native LBFGS implementation
     */
    static SimpleOptimizationResult optimizeStructure(
        curcuma::Molecule* molecule,
        EnergyCalculator* energy_calculator,
        const json& config = json{});

    /**
     * @brief Get available native optimization methods
     */
    static std::map<std::string, std::string> getAvailableNativeMethods();

private:
    /**
     * @brief Convert molecule to coordinate vector
     */
    static Vector moleculeToVector(const curcuma::Molecule& molecule);

    /**
     * @brief Convert coordinate vector back to molecule
     */
    static curcuma::Molecule vectorToMolecule(const Vector& coordinates, const curcuma::Molecule& template_mol);

    /**
     * @brief Configure LBFGS optimizer from JSON settings
     */
    static void configureOptimizer(LBFGS& optimizer, const json& config);

    /**
     * @brief Log optimization progress with modern formatting
     */
    static void logOptimizationProgress(int step, double energy, double energy_change,
        double rmsd_change, double gradient_norm);

    int m_memory_size;
    LBFGS::Method m_method;
};

} // namespace ModernOptimization