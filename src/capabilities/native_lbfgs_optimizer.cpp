/*
 * <Native LBFGS Optimizer Adapter Implementation - Claude Generated>
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

#include "native_lbfgs_optimizer.h"
#include <chrono>
#include <cmath>

namespace ModernOptimization {

// Claude Generated - Constructor
NativeLBFGSOptimizer::NativeLBFGSOptimizer(int memory_size, LBFGS::Method method)
    : m_memory_size(memory_size)
    , m_method(method)
{
}

// Claude Generated - Main optimization interface
SimpleOptimizationResult NativeLBFGSOptimizer::optimizeStructure(
    curcuma::Molecule* molecule,
    EnergyCalculator* energy_calculator,
    const json& config)
{
    if (!molecule) {
        return SimpleOptimizationResult::failed_result("Null molecule provided", "native_lbfgs");
    }
    if (!energy_calculator) {
        return SimpleOptimizationResult::failed_result("Null energy calculator provided", "native_lbfgs");
    }

    try {
        auto start_time = std::chrono::high_resolution_clock::now();

        // Parse configuration
        int memory_size = config.value("lbfgs_memory", 10);
        std::string method_str = config.value("native_method", "lbfgs");
        int max_iterations = config.value("MaxIter", 1000);
        double gradient_threshold = config.value("GradNorm", 1e-4);
        double energy_threshold = config.value("dE", 1e-6);

        // Parse method
        LBFGS::Method method = LBFGS::Method::LBFGS;
        if (method_str == "diis") {
            method = LBFGS::Method::DIIS;
        } else if (method_str == "rfo") {
            method = LBFGS::Method::RFO;
        }

        // Create and configure optimizer
        LBFGS optimizer(memory_size);
        optimizer.setOptimizationMethod(method);
        optimizer.setEnergyCalculator(energy_calculator);
        configureOptimizer(optimizer, config);

        // Initialize with current molecule coordinates
        Vector initial_coords = moleculeToVector(*molecule);
        optimizer.initialize(molecule->AtomCount(), initial_coords);

        // Log start
        CurcumaLogger::header("Native LBFGS Optimization");
        CurcumaLogger::param("Method", method_str);
        CurcumaLogger::param("Memory size", memory_size);
        CurcumaLogger::param("Max iterations", max_iterations);
        CurcumaLogger::energy_abs(optimizer.getCurrentEnergy(), "Initial energy");

        // Optimization loop
        double previous_energy = optimizer.getCurrentEnergy();
        bool converged = false;
        int iteration = 0;

        while (iteration < max_iterations && !converged && !optimizer.hasError()) {
            // Take optimization step
            Vector new_coords = optimizer.step();

            // Check convergence
            double current_energy = optimizer.getCurrentEnergy();
            double energy_change = std::abs(current_energy - previous_energy);
            double gradient_norm = optimizer.getCurrentGradient().norm();

            // Calculate RMSD change
            Vector coord_change = new_coords - initial_coords;
            double rmsd_change = std::sqrt(coord_change.squaredNorm() / (3.0 * molecule->AtomCount()));

            // Log progress
            if (iteration % 10 == 0 || iteration < 10) {
                logOptimizationProgress(iteration + 1, current_energy, energy_change,
                    rmsd_change, gradient_norm);
            }

            // Convergence check
            converged = optimizer.isConverged(energy_threshold, gradient_threshold);
            if (converged) {
                CurcumaLogger::success("Native LBFGS optimization converged!");
                break;
            }

            previous_energy = current_energy;
            iteration++;
        }

        // Get final results
        Vector final_coords = optimizer.getCurrentSolution();
        curcuma::Molecule final_molecule = vectorToMolecule(final_coords, *molecule);
        final_molecule.setEnergy(optimizer.getCurrentEnergy());

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        double optimization_time = duration.count() / 1000.0;

        // Final logging
        CurcumaLogger::energy_abs(optimizer.getCurrentEnergy(), "Final energy");
        CurcumaLogger::param("Iterations", iteration);
        CurcumaLogger::param("Time", fmt::format("{:.3f} seconds", optimization_time));

        if (converged) {
            return SimpleOptimizationResult::success_result(
                final_molecule, optimizer.getCurrentEnergy(),
                "native_" + method_str, iteration, optimization_time);
        } else if (optimizer.hasError()) {
            return SimpleOptimizationResult::failed_result("Optimization error encountered", "native_" + method_str);
        } else {
            return SimpleOptimizationResult::failed_result("Maximum iterations reached", "native_" + method_str);
        }

    } catch (const std::exception& e) {
        return SimpleOptimizationResult::failed_result(
            fmt::format("Native LBFGS error: {}", e.what()), "native_lbfgs");
    }
}

// Claude Generated - Available methods
std::map<std::string, std::string> NativeLBFGSOptimizer::getAvailableNativeMethods()
{
    return {
        { "native_lbfgs", "Native L-BFGS implementation (Claude 3.5)" },
        { "native_diis", "Native DIIS with L-BFGS fallback (Claude 3.5)" },
        { "native_rfo", "Native Rational Function Optimization (Claude 3.5)" }
    };
}

// Claude Generated - Molecule conversion utilities
Vector NativeLBFGSOptimizer::moleculeToVector(const curcuma::Molecule& molecule)
{
    int natoms = molecule.AtomCount();
    Vector coords(3 * natoms);

    for (int i = 0; i < natoms; ++i) {
        auto pos = molecule.Atom(i);
        coords[3 * i + 0] = pos.x();
        coords[3 * i + 1] = pos.y();
        coords[3 * i + 2] = pos.z();
    }

    return coords;
}

curcuma::Molecule NativeLBFGSOptimizer::vectorToMolecule(const Vector& coordinates, const curcuma::Molecule& template_mol)
{
    curcuma::Molecule result = template_mol; // Copy all properties

    int natoms = template_mol.AtomCount();
    for (int i = 0; i < natoms; ++i) {
        Position new_pos(coordinates[3 * i + 0], coordinates[3 * i + 1], coordinates[3 * i + 2]);
        result.setAtom(i, template_mol.Atom(i).first, new_pos);
    }

    return result;
}

// Claude Generated - Configuration helper
void NativeLBFGSOptimizer::configureOptimizer(LBFGS& optimizer, const json& config)
{
    // Set verbosity
    int verbosity = config.value("printOutput", true) ? 2 : 1;
    optimizer.setVerbosity(verbosity);

    // DIIS parameters
    if (config.contains("diis_hist")) {
        int diis_hist = config["diis_hist"];
        int diis_start = config.value("diis_start", diis_hist);
        optimizer.setDIISParameters(diis_hist, diis_start);
    }

    // RFO parameters
    if (config.contains("rfo_lambda")) {
        optimizer.setLambda(config["rfo_lambda"]);
    }

    // Mass weighting (if provided)
    if (config.contains("atomic_masses")) {
        std::vector<double> masses = config["atomic_masses"];
        optimizer.setMasses(masses);
    }
}

// Claude Generated - Progress logging
void NativeLBFGSOptimizer::logOptimizationProgress(int step, double energy, double energy_change,
    double rmsd_change, double gradient_norm)
{
    CurcumaLogger::info_fmt("Step {:3d}: E={:.8f} Eh, ΔE={:.2e} kJ/mol, "
                            "RMSD={:.4f} Å, |g|={:.2e} Eh/Bohr",
        step, energy, energy_change * CURCUMA_EH_TO_KJMOL,
        rmsd_change, gradient_norm);
}

} // namespace ModernOptimization