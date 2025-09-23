/*
 * <Casino - Monte Carlo simulation for molecular systems.>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "casino.h"
#include "src/tools/formats.h"
#include "src/core/units.h"
#include "src/core/elements.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

Casino::Casino(const json& controller, bool silent)
    : CurcumaMethod(CasinoJson, controller, silent)
    , m_energy_calculator(nullptr)
    , m_current_step(0)
    , m_total_moves(0)
    , m_accepted_moves(0)
    , m_current_energy(0.0)
    , m_uniform_dist(0.0, 1.0)
    , m_normal_dist(0.0, 1.0)
{
    UpdateController(controller);
    m_config = MergeJson(CasinoJson, controller);

    // Extract simulation parameters
    m_total_steps = Json2KeyWord<int>(m_config, "steps");
    m_temperature = Json2KeyWord<double>(m_config, "temperature");
    m_step_size = Json2KeyWord<double>(m_config, "step_size");
    m_output_frequency = Json2KeyWord<int>(m_config, "output_frequency");
    m_energy_frequency = Json2KeyWord<int>(m_config, "energy_frequency");
    m_output_file = Json2KeyWord<std::string>(m_config, "output_file");
    m_move_type = Json2KeyWord<std::string>(m_config, "move_type");
    m_acceptance_target = Json2KeyWord<double>(m_config, "acceptance_target");
    m_adaptive_step = Json2KeyWord<bool>(m_config, "adaptive_step");
    m_verbose = Json2KeyWord<bool>(m_config, "verbose");

    // Enhanced sampling options
    m_bias_potential = Json2KeyWord<bool>(m_config, "bias_potential");
    m_umbrella_sampling = Json2KeyWord<bool>(m_config, "umbrella_sampling");
    m_metadynamics = Json2KeyWord<bool>(m_config, "metadynamics");

    // Initialize random number generator
    int seed = Json2KeyWord<int>(m_config, "seed");
    m_rng.seed(seed);

    if (!silent) {
        CurcumaLogger::info_fmt("Monte Carlo simulation initialized with {} steps", m_total_steps);
        CurcumaLogger::info_fmt("Temperature: {:.1f} K", m_temperature);
        CurcumaLogger::info_fmt("Initial step size: {:.3f} Å", m_step_size);
    }
}

Casino::~Casino()
{
    delete m_energy_calculator;
}

bool Casino::Initialise()
{
    if (m_filename.empty() && m_molecule.AtomCount() == 0) {
        CurcumaLogger::error("No input molecule provided for Monte Carlo simulation");
        return false;
    }

    // Load molecule if filename provided
    if (!m_filename.empty()) {
        m_molecule = Files::LoadFile(m_filename);
        if (m_molecule.AtomCount() == 0) {
            CurcumaLogger::error_fmt("Could not load molecule from file: {}", m_filename);
            return false;
        }
    }

    // Initialize energy calculator
    std::string method = m_config.value("method", "uff");
    m_energy_calculator = new EnergyCalculator(method, m_config);
    m_energy_calculator->setMolecule(m_molecule.getMolInfo());

    // Calculate initial energy
    m_current_energy = calculateEnergy();

    if (!m_silent) {
        CurcumaLogger::info("=== Monte Carlo Simulation Setup ===");
        CurcumaLogger::param("atoms", fmt::format("{}", m_molecule.AtomCount()));
        CurcumaLogger::param("method", fmt::format("{}", method));
        CurcumaLogger::param("initial_energy", fmt::format("{:.6f} kJ/mol", m_current_energy));

        if (m_molecule.AtomCount() > 0) {
            // Check for CG elements
            bool has_cg = false;
            for (int atom : m_molecule.Atoms()) {
                if (atom == CG_ELEMENT) {
                    has_cg = true;
                    break;
                }
            }
            if (has_cg) {
                CurcumaLogger::param("system_type", "Coarse-grained");
            } else {
                CurcumaLogger::param("system_type", "Atomistic");
            }
        }
    }

    // Reserve space for trajectories
    m_energy_trajectory.reserve(m_total_steps / m_energy_frequency + 100);
    m_acceptance_history.reserve(m_total_steps / 100 + 100);

    return true;
}

void Casino::start()
{
    if (!Initialise()) {
        CurcumaLogger::error("Failed to initialize Casino Monte Carlo simulation");
        return;
    }

    CurcumaLogger::header("Starting Casino Monte Carlo Simulation");
    auto start_time = std::chrono::high_resolution_clock::now();

    // Main MC loop
    for (m_current_step = 1; m_current_step <= m_total_steps; ++m_current_step) {

        // Perform Monte Carlo move
        bool accepted = performMove();

        // Update statistics
        m_total_moves++;
        if (accepted) {
            m_accepted_moves++;
        }

        // Output energy
        if (m_current_step % m_energy_frequency == 0) {
            m_energy_trajectory.push_back(m_current_energy);
            outputStatistics(m_current_step);
        }

        // Write trajectory
        if (m_current_step % m_output_frequency == 0) {
            writeTrajectoryFrame(m_current_step);
        }

        // Adaptive step size adjustment
        if (m_adaptive_step && m_current_step % 1000 == 0) {
            updateStepSize();
        }

        // Progress reporting
        if (!m_silent && m_current_step % (m_total_steps / 10) == 0) {
            double progress = 100.0 * m_current_step / m_total_steps;
            double acceptance = getAcceptanceRatio();
            CurcumaLogger::info_fmt("Progress: {:.1f}% | Acceptance: {:.3f} | Energy: {:.6f} kJ/mol",
                                  progress, acceptance, m_current_energy);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

    // Final statistics
    CurcumaLogger::success("Casino Monte Carlo simulation completed!");
    CurcumaLogger::param("total_steps", fmt::format("{}", m_total_steps));
    CurcumaLogger::param("acceptance_ratio", fmt::format("{:.3f}", getAcceptanceRatio()));
    CurcumaLogger::param("final_energy", fmt::format("{:.6f} kJ/mol", m_current_energy));
    CurcumaLogger::param("simulation_time", fmt::format("{} seconds", duration.count()));
    CurcumaLogger::param("trajectory_file", fmt::format("{}", m_output_file));

    // Write final structure
    std::string final_file = m_output_file;
    final_file.replace(final_file.find(".xyz"), 4, "_final.xyz");
    m_molecule.writeXYZFile(final_file);
    CurcumaLogger::param("final_structure", fmt::format("{}", final_file));
}

bool Casino::performMove()
{
    // Store current state
    m_trial_molecule = m_molecule;
    double energy_old = m_current_energy;

    // Apply random move
    applyRandomMove();

    // Calculate new energy
    m_energy_calculator->setMolecule(m_trial_molecule.getMolInfo());
    double energy_new = calculateEnergy();

    // Add bias potential if enabled
    if (m_bias_potential || m_umbrella_sampling || m_metadynamics) {
        energy_new += calculateBiasPotential();
    }

    // Accept or reject move
    bool accepted = acceptMove(energy_old, energy_new);

    if (accepted) {
        m_molecule = m_trial_molecule;
        m_current_energy = energy_new;

        // Update enhanced sampling
        if (m_metadynamics) {
            updateMetadynamics();
        }
    } else {
        // Restore energy calculator to old state
        m_energy_calculator->setMolecule(m_molecule.getMolInfo());
    }

    return accepted;
}

double Casino::calculateEnergy()
{
    double energy = m_energy_calculator->CalculateEnergy(false);
    return CurcumaUnit::Energy::hartree_to_kjmol(energy); // Convert to kJ/mol
}

void Casino::applyRandomMove()
{
    Geometry geometry = m_trial_molecule.getGeometry();
    int n_atoms = m_trial_molecule.AtomCount();

    if (m_move_type == "translation" || m_move_type == "mixed") {
        // Random translation moves
        for (int i = 0; i < n_atoms; ++i) {
            double dx = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_step_size;
            double dy = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_step_size;
            double dz = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_step_size;

            geometry(i, 0) += dx;
            geometry(i, 1) += dy;
            geometry(i, 2) += dz;
        }
    }

    if (m_move_type == "rotation" || m_move_type == "mixed") {
        // Random rotation around center of mass
        Position com = m_trial_molecule.MassCentroid();

        // Random rotation axis and angle
        double theta = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_step_size * 0.1; // Small rotations
        double phi = m_uniform_dist(m_rng) * 2.0 * M_PI;
        // cos_phi and sin_phi would be used for full 3D rotation
        // Currently implementing simplified z-axis rotation

        // Apply rotation (simplified - around z-axis for now)
        for (int i = 0; i < n_atoms; ++i) {
            double x = geometry(i, 0) - com[0];
            double y = geometry(i, 1) - com[1];

            geometry(i, 0) = com[0] + x * std::cos(theta) - y * std::sin(theta);
            geometry(i, 1) = com[1] + x * std::sin(theta) + y * std::cos(theta);
        }
    }

    m_trial_molecule.setGeometry(geometry);
}

bool Casino::acceptMove(double energy_old, double energy_new)
{
    double delta_E = energy_new - energy_old;

    // Always accept if energy decreases
    if (delta_E <= 0.0) {
        return true;
    }

    // Metropolis acceptance criterion
    double beta = 1.0 / (kB * m_temperature);
    double probability = std::exp(-beta * delta_E);

    return m_uniform_dist(m_rng) < probability;
}

void Casino::updateStepSize()
{
    if (m_total_moves < 100) return; // Need sufficient statistics

    double current_acceptance = static_cast<double>(m_accepted_moves) / m_total_moves;
    double ratio = current_acceptance / m_acceptance_target;

    // Adjust step size to maintain target acceptance ratio
    if (ratio > 1.2) {
        m_step_size *= 1.1; // Increase step size
    } else if (ratio < 0.8) {
        m_step_size *= 0.9; // Decrease step size
    }

    // Reasonable bounds
    m_step_size = std::max(0.001, std::min(m_step_size, 1.0));

    if (m_verbose) {
        CurcumaLogger::info_fmt("Step {} | Acceptance: {:.3f} | Step size: {:.3f} Å",
                              m_current_step, current_acceptance, m_step_size);
    }
}

void Casino::writeTrajectoryFrame(int step)
{
    static bool first_write = true;

    std::ofstream file;
    if (first_write) {
        file.open(m_output_file);
        first_write = false;
    } else {
        file.open(m_output_file, std::ios::app);
    }

    if (file.is_open()) {
        // Write XYZ format with energy in comment
        file << m_molecule.AtomCount() << std::endl;
        file << "Step " << step << " | Energy: " << std::fixed << std::setprecision(6)
             << m_current_energy << " kJ/mol" << std::endl;

        for (int i = 0; i < m_molecule.AtomCount(); ++i) {
            Position pos = m_molecule.getGeometry().row(i);
            file << Elements::ElementAbbr[m_molecule.Atom(i).first] << " "
                 << std::fixed << std::setprecision(6)
                 << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
        }
        file.close();
    }
}

void Casino::outputStatistics(int step)
{
    if (m_verbose && step % (m_energy_frequency * 10) == 0) {
        double acceptance = getAcceptanceRatio();
        CurcumaLogger::param("step_info", fmt::format("Step: {} | Energy: {:.6f} kJ/mol | Acceptance: {:.3f}",
                                step, m_current_energy, acceptance));
    }
}

double Casino::calculateBiasPotential()
{
    // Placeholder for bias potential implementation
    // This would include umbrella sampling, metadynamics, etc.
    return 0.0;
}

void Casino::updateMetadynamics()
{
    // Placeholder for metadynamics implementation
    // Would store current geometry for history-dependent bias
    if (m_current_step % 100 == 0) {
        m_metadynamics_history.push_back(m_molecule.getGeometry());
    }
}

void Casino::printHelp() const
{
    std::cout << "Casino - Monte Carlo Simulation for all molecular systems" << std::endl;
    std::cout << "=========================================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: curcuma -casino file.xyz [options]" << std::endl;
    std::cout << "       curcuma -casino file.vtf [options]" << std::endl;
    std::cout << std::endl;
    std::cout << "Supported Formats: XYZ, VTF, MOL2, SDF, PDB" << std::endl;
    std::cout << "Supported Methods: All Curcuma force fields and QM methods" << std::endl;
    std::cout << std::endl;
    std::cout << "Basic Options:" << std::endl;
    std::cout << "  -steps N                Number of MC steps (default: 10000)" << std::endl;
    std::cout << "  -temperature T          Temperature in Kelvin (default: 300.0)" << std::endl;
    std::cout << "  -step_size S            Maximum displacement in Å (default: 0.1)" << std::endl;
    std::cout << "  -method METHOD          Force field or QM method (default: uff)" << std::endl;
    std::cout << std::endl;
    std::cout << "Move Types:" << std::endl;
    std::cout << "  -move_type translation  Translation moves only" << std::endl;
    std::cout << "  -move_type rotation     Rotation moves only" << std::endl;
    std::cout << "  -move_type mixed        Combined translation and rotation (default)" << std::endl;
    std::cout << std::endl;
    std::cout << "Output Options:" << std::endl;
    std::cout << "  -output_file FILE       Trajectory output file (default: casino_trajectory.xyz)" << std::endl;
    std::cout << "  -output_frequency N     Trajectory output frequency (default: 100)" << std::endl;
    std::cout << "  -energy_frequency N     Energy output frequency (default: 10)" << std::endl;
    std::cout << std::endl;
    std::cout << "Advanced Options:" << std::endl;
    std::cout << "  -adaptive_step true     Adaptive step size adjustment (default)" << std::endl;
    std::cout << "  -acceptance_target R    Target acceptance ratio (default: 0.5)" << std::endl;
    std::cout << "  -seed N                 Random seed for reproducibility (default: 42)" << std::endl;
    std::cout << std::endl;
    std::cout << "Enhanced Sampling (Experimental):" << std::endl;
    std::cout << "  -umbrella_sampling true Enable umbrella sampling" << std::endl;
    std::cout << "  -metadynamics true      Enable metadynamics" << std::endl;
    std::cout << "  -bias_potential true    Enable custom bias potential" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  curcuma -casino molecule.xyz -steps 50000 -temperature 500" << std::endl;
    std::cout << "  curcuma -casino cg_polymer.vtf -method uff -move_type mixed" << std::endl;
    std::cout << "  curcuma -casino protein.pdb -method gfn2 -adaptive_step true" << std::endl;
    std::cout << std::endl;
    std::cout << "Output Files:" << std::endl;
    std::cout << "  casino_trajectory.xyz       Full trajectory with energies" << std::endl;
    std::cout << "  casino_trajectory_final.xyz Final optimized structure" << std::endl;
}