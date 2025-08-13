/*
 * <Optimizer Driver Implementation - Template Method Pattern>
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

#include "optimizer_driver.h"
#include "src/tools/general.h"
#include <cmath>
#include <fstream>
#include <stdexcept>

namespace Optimization {

// Claude Generated - OptimizationContext implementation
bool OptimizationContext::isValid() const
{
    if (!energy_calculator)
        return false;
    if (energy_threshold <= 0)
        return false;
    if (rmsd_threshold <= 0)
        return false;
    if (gradient_threshold <= 0)
        return false;
    if (max_iterations <= 0)
        return false;
    if (threads <= 0)
        return false;
    return true;
}

std::string OptimizationContext::getValidationErrors() const
{
    std::vector<std::string> errors;

    if (!energy_calculator)
        errors.push_back("Energy calculator not set");
    if (energy_threshold <= 0)
        errors.push_back("Energy threshold must be positive");
    if (rmsd_threshold <= 0)
        errors.push_back("RMSD threshold must be positive");
    if (gradient_threshold <= 0)
        errors.push_back("Gradient threshold must be positive");
    if (max_iterations <= 0)
        errors.push_back("Max iterations must be positive");
    if (threads <= 0)
        errors.push_back("Thread count must be positive");

    std::string result;
    for (size_t i = 0; i < errors.size(); ++i) {
        if (i > 0)
            result += "; ";
        result += errors[i];
    }
    return result;
}

OptimizationContext OptimizationContext::fromJson(const json& config, EnergyCalculator* calc)
{
    OptimizationContext context;
    context.energy_calculator = calc;

    // Load convergence criteria
    if (config.contains("energy_threshold"))
        context.energy_threshold = config["energy_threshold"];
    if (config.contains("rmsd_threshold"))
        context.rmsd_threshold = config["rmsd_threshold"];
    if (config.contains("gradient_threshold"))
        context.gradient_threshold = config["gradient_threshold"];
    if (config.contains("max_iterations"))
        context.max_iterations = config["max_iterations"];
    if (config.contains("convergence_count"))
        context.convergence_count = config["convergence_count"];

    // Load performance settings
    if (config.contains("threads"))
        context.threads = config["threads"];
    if (config.contains("single_step_mode"))
        context.single_step_mode = config["single_step_mode"];
    if (config.contains("max_energy_rise"))
        context.max_energy_rise = config["max_energy_rise"];

    // Load output settings
    if (config.contains("write_trajectory"))
        context.write_trajectory = config["write_trajectory"];
    if (config.contains("verbose"))
        context.verbose = config["verbose"];
    if (config.contains("print_output"))
        context.print_output = config["print_output"];

    // Load molecular properties
    if (config.contains("charge"))
        context.charge = config["charge"];
    if (config.contains("spin"))
        context.spin = config["spin"];

    return context;
}

// Claude Generated - OptimizerDriver implementation
OptimizerDriver::OptimizerDriver()
{
    m_configuration = OptimizerInterfaceJson;
}

bool OptimizerDriver::InitializeOptimization(const Molecule& molecule)
{
    m_molecule = molecule;
    return InitializeOptimization(&m_molecule);
}

bool OptimizerDriver::InitializeOptimization(const Molecule* molecule)
{
    if (!molecule) {
        CurcumaLogger::error("Cannot initialize optimization with null molecule");
        return false;
    }

    m_molecule = *molecule;

    // Validate context
    if (!m_context.isValid()) {
        CurcumaLogger::error_fmt("Invalid optimization context: {}", m_context.getValidationErrors());
        return false;
    }

    // Setup energy calculator
    if (!m_context.energy_calculator) {
        CurcumaLogger::error("Energy calculator not configured");
        return false;
    }

    // Setup RMSD driver for convergence checking
    json rmsd_config = {
        { "reorder", false },
        { "check", false },
        { "heavy", false },
        { "silent", true },
        { "threads", 1 }
    };
    m_context.rmsd_driver = std::make_unique<RMSDDriver>(rmsd_config);

    // Initialize trajectory
    m_trajectory.clear();
    m_energy_trajectory.clear();
    m_current_iteration = 0;
    m_converged = false;

    // Calculate initial energy and gradient
    Vector coordinates = Tools::Mol2Coord(m_molecule);
    if (!evaluateEnergyAndGradient(coordinates, m_current_energy, m_current_gradient)) {
        CurcumaLogger::error("Failed to calculate initial energy and gradient");
        return false;
    }

    m_initial_energy = m_current_energy;
    updateTrajectory(m_molecule, m_current_energy);

    // Method-specific initialization
    if (!InitializeOptimizerInternal()) {
        CurcumaLogger::error("Method-specific initialization failed");
        return false;
    }

    CurcumaLogger::success_fmt("Optimization initialized with {} method", getName());
    CurcumaLogger::energy_abs(m_initial_energy, "Initial energy");
    CurcumaLogger::param("Initial gradient norm", fmt::format("{:.6e} Eh/Bohr", m_current_gradient.norm()));

    return true;
}

bool OptimizerDriver::InitializeOptimization(const double* coordinates, int atom_count)
{
    // Create molecule from coordinate array
    Molecule mol(atom_count, 0);
    for (int i = 0; i < atom_count; ++i) {
        int element = 6; // Default to carbon - should be provided in real implementation
        Vector3 position(coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]);
        mol.addPair({ element, position });
    }

    return InitializeOptimization(mol);
}

bool OptimizerDriver::UpdateGeometry(const Molecule& molecule)
{
    m_molecule = molecule;
    Vector coordinates = Tools::Mol2Coord(m_molecule);
    return evaluateEnergyAndGradient(coordinates, m_current_energy, m_current_gradient);
}

bool OptimizerDriver::UpdateGeometry(const double* coordinates)
{
    // Update molecule coordinates
    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        Vector3 pos(coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]);
        m_molecule.setAtom(i, m_molecule.Atom(i).first, pos);
    }

    return evaluateEnergyAndGradient(Vector::Map(coordinates, 3 * m_molecule.AtomCount()),
        m_current_energy, m_current_gradient);
}

OptimizationResult OptimizerDriver::Optimize(bool write_trajectory, bool verbose)
{
    // Synchronize verbosity settings (analog to recent RMSDTraj fixes)
    CurcumaLogger::set_verbosity(verbose ? 3 : 2);
    m_context.write_trajectory = write_trajectory;
    m_context.verbose = verbose;

    m_start_time = std::chrono::high_resolution_clock::now();

    CurcumaLogger::header(fmt::format("{} Geometry Optimization", getName()));

    // Display optimization parameters (analog to QM system parameter display)
    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::param_table(m_configuration, "Optimization Parameters");
    }

    CurcumaLogger::energy_abs(m_initial_energy, "Initial energy");
    CurcumaLogger::param("Max iterations", m_context.max_iterations);
    CurcumaLogger::param("Convergence criteria",
        fmt::format("ΔE < {:.1e} kJ/mol, ΔRMSD < {:.3f} Å, |∇| < {:.1e} Eh/Bohr",
            m_context.energy_threshold, m_context.rmsd_threshold, m_context.gradient_threshold));

    try {
        // Main optimization loop (Template Method Pattern)
        for (m_current_iteration = 1; m_current_iteration <= m_context.max_iterations; ++m_current_iteration) {

            // Calculate optimization step (method-specific)
            Vector current_coords = Tools::Mol2Coord(m_molecule);
            Vector step = CalculateOptimizationStep(current_coords, m_current_gradient);

            if (step.norm() == 0) {
                CurcumaLogger::warn("Zero optimization step calculated - possible convergence or error");
                break;
            }

            // Apply step and evaluate new energy/gradient
            Vector new_coords = current_coords + step;
            double new_energy;
            Vector new_gradient;

            if (!evaluateEnergyAndGradient(new_coords, new_energy, new_gradient)) {
                CurcumaLogger::error_fmt("Energy evaluation failed at iteration {}", m_current_iteration);
                return OptimizationResult::failed_result("Energy evaluation failed during optimization");
            }

            // Check for energy rise limit
            double energy_change_kjmol = (new_energy - m_current_energy) * CURCUMA_EH_TO_KJMOL;
            if (energy_change_kjmol > m_context.max_energy_rise) {
                CurcumaLogger::warn_fmt("Energy rise ({:.2f} kJ/mol) exceeds limit ({:.1f} kJ/mol)",
                    energy_change_kjmol, m_context.max_energy_rise);
                return OptimizationResult::failed_result("Energy rise exceeded maximum allowed");
            }

            // Update molecule geometry
            Tools::Coord2Mol(new_coords, m_molecule);

            // Calculate RMSD change
            double rmsd_change = 0.0;
            if (m_trajectory.size() > 0) {
                rmsd_change = calculateRMSD(m_molecule, m_trajectory.back());
            }

            // Update state (method-specific)
            UpdateOptimizerState(new_coords, new_gradient, new_energy);

            // Update common state
            m_current_energy = new_energy;
            m_current_gradient = new_gradient;
            updateTrajectory(m_molecule, new_energy);

            // Log progress (analog to recent RMSDTraj improvements)
            if (CurcumaLogger::get_verbosity() >= 3 || m_current_iteration % 10 == 0) {
                logOptimizationStep(m_current_iteration, new_energy, energy_change_kjmol,
                    rmsd_change, new_gradient.norm());
            }

            // Check convergence
            if (checkConvergence(energy_change_kjmol, rmsd_change, new_gradient.norm()) && CheckMethodSpecificConvergence()) {
                m_converged = true;
                m_convergence_reason = "All convergence criteria satisfied";
                break;
            }

            // Single step mode (for debugging/testing)
            if (m_context.single_step_mode) {
                CurcumaLogger::info("Single step mode - stopping after one iteration");
                break;
            }
        }

        // Finalize optimization (method-specific)
        FinalizeOptimizationInternal();

        m_end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(m_end_time - m_start_time);
        double optimization_time = duration.count() / 1000.0;

        // Final reporting
        if (m_converged) {
            CurcumaLogger::success_fmt("Optimization converged after {} iterations", m_current_iteration);
            CurcumaLogger::success_fmt("Convergence reason: {}", m_convergence_reason);
        } else {
            CurcumaLogger::warn_fmt("Optimization did not converge within {} iterations", m_context.max_iterations);
        }

        CurcumaLogger::energy_abs(m_current_energy, "Final energy");
        CurcumaLogger::energy_rel(m_current_energy - m_initial_energy, "Energy change");
        CurcumaLogger::param("Final gradient norm", fmt::format("{:.6e} Eh/Bohr", m_current_gradient.norm()));
        CurcumaLogger::param("Optimization time", fmt::format("{:.3f} seconds", optimization_time));

        if (m_context.write_trajectory && !m_context.trajectory_filename.empty()) {
            // Write trajectory file
            std::ofstream trj_file(m_context.trajectory_filename);
            for (const auto& mol : m_trajectory) {
                trj_file << mol.getXYZString();
            }
            CurcumaLogger::success_fmt("Trajectory written to: {}", m_context.trajectory_filename);
        }

        // Create result
        OptimizationResult result = OptimizationResult::success_result(
            m_molecule, m_current_energy, m_current_iteration, optimization_time);
        result.success = m_converged;
        result.final_gradient = m_current_gradient;
        result.final_energy_change = (m_current_energy - m_initial_energy) * CURCUMA_EH_TO_KJMOL;
        result.final_gradient_norm = m_current_gradient.norm();
        result.trajectory = m_trajectory;
        result.energy_trajectory = m_energy_trajectory;

        if (!m_converged) {
            result.error_message = fmt::format("Did not converge within {} iterations", m_context.max_iterations);
        }

        return result;

    } catch (const std::exception& e) {
        CurcumaLogger::error_fmt("Optimization failed with exception: {}", e.what());
        return OptimizationResult::failed_result(e.what());
    }
}

// Claude Generated - Helper method implementations
bool OptimizerDriver::evaluateEnergyAndGradient(const Vector& coordinates, double& energy, Vector& gradient)
{
    try {
        Tools::Coord2Mol(coordinates, m_molecule);
        m_context.energy_calculator->setMolecule(m_molecule);

        energy = m_context.energy_calculator->CalculateEnergy(true);
        if (std::isnan(energy) || std::isinf(energy)) {
            return false;
        }

        Geometry grad_geom = m_context.energy_calculator->Gradient();
        gradient = Vector::Map(grad_geom.data(), grad_geom.size());

        // Apply constraints if specified
        if (m_context.use_constraints && !m_context.atom_constraints.empty()) {
            for (int i = 0; i < m_molecule.AtomCount(); ++i) {
                if (i < m_context.atom_constraints.size() && m_context.atom_constraints[i] == 0) {
                    gradient[3 * i] = 0.0;
                    gradient[3 * i + 1] = 0.0;
                    gradient[3 * i + 2] = 0.0;
                }
            }
        }

        return true;

    } catch (const std::exception& e) {
        CurcumaLogger::error_fmt("Energy/gradient evaluation failed: {}", e.what());
        return false;
    }
}

bool OptimizerDriver::checkConvergence(double energy_change, double rmsd_change, double gradient_norm) const
{
    int convergence_flags = 0;

    // Energy convergence (bit 0)
    if (std::abs(energy_change) < m_context.energy_threshold) {
        convergence_flags |= 1;
    }

    // RMSD convergence (bit 1)
    if (rmsd_change < m_context.rmsd_threshold) {
        convergence_flags |= 2;
    }

    // Gradient convergence (bit 2)
    if (gradient_norm < m_context.gradient_threshold) {
        convergence_flags |= 4;
    }

    // Check if required convergence criteria are met
    return (convergence_flags & m_context.convergence_count) == m_context.convergence_count;
}

void OptimizerDriver::updateTrajectory(const Molecule& new_structure, double energy)
{
    m_trajectory.push_back(new_structure);
    m_energy_trajectory.push_back(energy);
}

void OptimizerDriver::logOptimizationStep(int iteration, double energy, double energy_change,
    double rmsd_change, double gradient_norm) const
{
    CurcumaLogger::info_fmt("Step {:4d}: E = {:.8f} Eh", iteration, energy);
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::energy_rel(energy_change / CURCUMA_EH_TO_KJMOL, "  ΔE");
        CurcumaLogger::length(rmsd_change / CURCUMA_BOHR_TO_ANGSTROM, "  ΔRMSD");
        CurcumaLogger::param("  |∇|", fmt::format("{:.6e} Eh/Bohr", gradient_norm));
    }
}

double OptimizerDriver::calculateRMSD(const Molecule& mol1, const Molecule& mol2) const
{
    if (!m_context.rmsd_driver)
        return 0.0;

    try {
        m_context.rmsd_driver->setReference(mol1);
        m_context.rmsd_driver->setTarget(mol2);
        m_context.rmsd_driver->start();
        return m_context.rmsd_driver->RMSD();
    } catch (...) {
        return 0.0; // Fallback if RMSD calculation fails
    }
}

void OptimizerDriver::LoadConfiguration(const json& config)
{
    m_configuration = config;
    m_context = OptimizationContext::fromJson(config, m_context.energy_calculator);
}

void OptimizerDriver::setEnergyCalculator(EnergyCalculator* calculator)
{
    m_context.energy_calculator = calculator;
}

void OptimizerDriver::setConstraints(const std::vector<int>& constraints)
{
    m_context.atom_constraints = constraints;
    m_context.use_constraints = !constraints.empty();
}

void OptimizerDriver::setConvergenceCriteria(double energy_thresh, double rmsd_thresh, double grad_thresh)
{
    m_context.energy_threshold = energy_thresh;
    m_context.rmsd_threshold = rmsd_thresh;
    m_context.gradient_threshold = grad_thresh;
}

void OptimizerDriver::setTrajectoryFile(const std::string& filename)
{
    m_context.trajectory_filename = filename;
    m_context.write_trajectory = !filename.empty();
}

void OptimizerDriver::setBasename(const std::string& basename)
{
    m_context.output_basename = basename;
    if (!basename.empty()) {
        m_context.trajectory_filename = basename + ".trj.xyz";
    }
}

} // namespace Optimization