/*
 * <Modern CurcumaOpt - Strategy Pattern Integration>
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

#include "curcumamethod.h"
#include "external/CxxThreadPool/include/CxxThreadPool.hpp"
#include "optimization_validator.h"
#include "optimizer_factory.h"

namespace Optimization {

/**
 * @brief Modern CurcumaOpt with Strategy Pattern - Claude Generated
 * Drop-in replacement for legacy CurcumaOpt with improved architecture
 */
class CurcumaOptModern : public CurcumaMethod {
public:
    CurcumaOptModern(const json& controller, bool silent = false);
    virtual ~CurcumaOptModern() = default;

    // CurcumaMethod interface
    void start() override;

    // Configuration methods
    void setFileName(const std::string& filename);
    void addMolecule(const Molecule& molecule);
    void setSinglePoint(bool sp) { m_singlepoint = sp; }
    void setOptimizerType(const std::string& method_name);
    void setOptimizerType(OptimizerType type);

    // Results access
    const std::vector<Molecule>& getMolecules() const { return m_optimized_molecules; }
    const std::vector<OptimizationResult>& getResults() const { return m_results; }
    std::string getLastOutput() const { return m_last_output; }

    // Legacy compatibility methods
    Molecule optimizeMolecule(Molecule* initial, std::string& output,
        std::vector<Molecule>* intermediate = nullptr,
        Vector* charges = nullptr);
    double calculateSinglePoint(const Molecule* initial, std::string& output,
        Vector* charges = nullptr);

    // File output methods
    std::string getOptFile() const { return Basename() + ".opt.xyz"; }
    std::string getTrajFile() const { return Basename() + ".trj.xyz"; }

private:
    // CurcumaMethod interface implementation
    void LoadControlJson() override;
    void ReadControlFile() override {}
    StringList MethodName() const override { return { "opt", "sp", "optimization" }; }
    json WriteRestartInformation() override { return json{}; }
    bool LoadRestartInformation() override { return true; }

    // Processing methods
    void processInput();
    void processSingleMolecule(const Molecule& molecule, int index = 0);
    void processMultipleMolecules();
    void processSerial();
    void processParallel();

    // Output methods
    void writeResults();
    void writeOptimizedStructures();
    void writeTrajectories();
    void generateSummaryReport();

    // Configuration and state
    OptimizerType m_optimizer_type = OptimizerType::LBFGSPP;
    std::string m_method_name = "lbfgspp";
    json m_optimizer_config;

    // Input molecules and processing flags
    std::vector<Molecule> m_input_molecules;
    std::vector<Molecule> m_optimized_molecules;
    std::vector<OptimizationResult> m_results;

    bool m_file_input = false;
    bool m_molecule_input = false;
    bool m_singlepoint = false;
    bool m_write_xyz = true;
    bool m_write_trajectory = true;

    // Threading
    int m_threads = 1;
    bool m_serial_mode = false;

    // Output
    std::string m_last_output;

    // Energy calculator (inherited from CurcumaMethod)
    EnergyCalculator* getEnergyCalculator();
};

/**
 * @brief Thread-safe optimization worker - Claude Generated
 * Modern replacement for legacy OptThread using new optimizer system
 */
class ModernOptThread : public CxxThread {
public:
    ModernOptThread(CurcumaOptModern* parent, const Molecule& molecule,
        OptimizerType type, const json& config);
    virtual ~ModernOptThread() = default;

    int execute() override;

    // Results access
    OptimizationResult getResult() const { return m_result; }
    Molecule getOptimizedMolecule() const { return m_result.final_molecule; }
    std::string getOutput() const { return m_output; }

private:
    CurcumaOptModern* m_parent;
    Molecule m_input_molecule;
    OptimizerType m_optimizer_type;
    json m_config;

    OptimizationResult m_result;
    std::string m_output;
    std::string m_basename;
};

} // namespace Optimization

// Legacy typedef for backward compatibility
using CurcumaOpt = Optimization::CurcumaOptModern;