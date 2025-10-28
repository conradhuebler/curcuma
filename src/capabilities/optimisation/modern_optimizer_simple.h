/*
 * <Modern Optimizer - Simplified Version for Integration>
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
#include "src/core/parameter_macros.h"
#include "src/capabilities/curcumamethod.h"
#include <map>
#include <string>

using json = nlohmann::json;

/* Claude Generated 2025: Modern Optimizer Parameter Registry - replaces static OptimizerJson */
BEGIN_PARAMETER_DEFINITION(modern_optimizer)
    // Convergence criteria
    PARAM(max_iter, Int, 1000, "Maximum optimization iterations", "Convergence", { "MaxIter" })
    PARAM(d_e, Double, 1e-6, "Energy convergence threshold (Hartree)", "Convergence", { "dE" })
    PARAM(grad_norm, Double, 5e-4, "Gradient norm convergence threshold (Hartree/Bohr)", "Convergence", { "GradNorm" })

    // L-BFGS parameters
    PARAM(memory_size, Int, 10, "L-BFGS memory size", "LBFGS", {})

    // DIIS parameters
    PARAM(diis_hist, Int, 10, "DIIS history size", "DIIS", {})
    PARAM(diis_start, Int, 5, "DIIS start iteration", "DIIS", {})

    // RFO (Rational Function Optimizer) parameters
    PARAM(trust_radius, Double, 0.05, "RFO trust radius (Bohr)", "RFO", {})
    PARAM(energy_threshold, Double, 1e-6, "RFO energy convergence threshold (Hartree)", "RFO", {})
    PARAM(eigenvalue_shift, Double, 1e-4, "RFO eigenvalue shift for Hessian regularization", "RFO", {})
    PARAM(lambda, Double, 0.1, "Legacy RFO damping parameter", "RFO", {})

    // Output control
    PARAM(verbosity, Int, 2, "Verbosity level (0-3)", "Output", {})
END_PARAMETER_DEFINITION

// Default JSON configuration for native optimizers - Claude Generated
static const json OptimizerJson = {
    { "verbosity", 2 },           // Standard verbosity level
    { "MaxIter", 1000 },          // Maximum iterations
    { "dE", 1e-6 },              // Energy convergence threshold (Eh)
    { "GradNorm", 5e-4 },        // Gradient norm threshold (Eh/Bohr)
    { "memory_size", 10 },        // L-BFGS memory size
    { "diis_hist", 10 },         // DIIS history size
    { "diis_start", 5 },         // DIIS start iteration
    { "trust_radius", 0.05 },    // RFO trust radius (Bohr)
    { "energy_threshold", 1e-6 }, // RFO energy threshold
    { "eigenvalue_shift", 1e-4 }, // RFO eigenvalue shift
    { "lambda", 0.1 }            // Legacy RFO parameter
};
using curcuma::Molecule;

namespace ModernOptimization {

/**
 * @brief Simple optimization result - Claude Generated
 */
struct SimpleOptimizationResult {
    bool success = false;
    std::string method_used;
    std::string error_message;

    // Final state
    Molecule final_molecule;
    double final_energy = 0.0;
    int iterations_performed = 0;
    double optimization_time_seconds = 0.0;

    // Factory methods
    static SimpleOptimizationResult success_result(const Molecule& mol, double energy,
        const std::string& method, int iterations, double time);
    static SimpleOptimizationResult failed_result(const std::string& error, const std::string& method);
};

/**
 * @brief Modern optimizer dispatcher - Claude Generated
 * Simplified version for integration testing with CurcumaMethod integration
 */
class ModernOptimizerDispatcher : public CurcumaMethod {
public:
    // Constructor for CurcumaMethod integration - Claude Generated
    ModernOptimizerDispatcher(const json& controller, bool silent = false);
    
    // CurcumaMethod interface implementation
    void start() override;
    void printHelp() const override;
    
    // Required pure virtual methods from CurcumaMethod - Claude Generated
    json WriteRestartInformation() override { return json{}; }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return StringList{"modern_optimizer"}; }
    void ReadControlFile() override { }
    void LoadControlJson() override { }
    
    /**
     * @brief Available optimizer types
     */
    enum class OptimizerType {
        LBFGSPP, // External LBFGSpp library (working)
        INTERNAL, // Internal LBFGS (placeholder)
        NATIVE_LBFGS, // Native LBFGS implementation (Claude 3.5)
        NATIVE_DIIS, // Native DIIS implementation (Claude 3.5)
        NATIVE_RFO, // Native RFO implementation (Claude 3.5)
        AUTO // Automatic selection
    };

    /**
     * @brief Parse optimizer type from string
     */
    static OptimizerType parseOptimizerType(const std::string& method_name);

    /**
     * @brief Get available optimizers
     */
    static std::map<std::string, std::string> getAvailableOptimizers();

    /**
     * @brief Optimize molecule with modern system
     */
    static SimpleOptimizationResult optimizeStructure(
        Molecule* molecule,
        const std::string& method_name,
        EnergyCalculator* energy_calculator,
        const json& config = json{});

    /**
     * @brief Demonstrate modern optimizer features
     */
    static void demonstrateModernFeatures(const Molecule& molecule);

private:
    // Implementation methods
    static SimpleOptimizationResult optimizeWithLBFGSpp(Molecule* molecule,
        EnergyCalculator* calc,
        const json& config);
    static SimpleOptimizationResult optimizeWithInternal(Molecule* molecule,
        EnergyCalculator* calc,
        const json& config);

    // Native optimization methods - Claude Generated
    static SimpleOptimizationResult optimizeWithNativeLBFGS(Molecule* molecule,
        EnergyCalculator* calc,
        const json& config);
    static SimpleOptimizationResult optimizeWithNativeDIIS(Molecule* molecule,
        EnergyCalculator* calc,
        const json& config);
    static SimpleOptimizationResult optimizeWithNativeRFO(Molecule* molecule,
        EnergyCalculator* calc,
        const json& config);

    // Helper methods
    static json getDefaultConfig(OptimizerType type);
    static void logOptimizationHeader(const std::string& method, const Molecule& molecule);
    static void logOptimizationResult(const SimpleOptimizationResult& result);
};

} // namespace ModernOptimization