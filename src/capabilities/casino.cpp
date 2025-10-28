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
#include "src/core/parameter_registry.h"

#include <algorithm>
#include <filesystem>
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
    , m_scnp_parser(silent)
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
    m_move_strategy = Json2KeyWord<std::string>(m_config, "move_strategy");
    m_acceptance_target = Json2KeyWord<double>(m_config, "acceptance_target");
    m_adaptive_step = Json2KeyWord<bool>(m_config, "adaptive_step");
    m_verbose = Json2KeyWord<bool>(m_config, "verbose");

    // Enhanced MC parameters - Claude Generated
    m_pivot_moves = Json2KeyWord<bool>(m_config, "pivot_moves");
    m_orientational_moves = Json2KeyWord<bool>(m_config, "orientational_moves");
    m_local_energy_updates = Json2KeyWord<bool>(m_config, "local_energy_updates");

    // Enhanced sampling options
    m_bias_potential = Json2KeyWord<bool>(m_config, "bias_potential");
    m_umbrella_sampling = Json2KeyWord<bool>(m_config, "umbrella_sampling");
    m_metadynamics = Json2KeyWord<bool>(m_config, "metadynamics");

    // Initialize random number generator
    int seed = Json2KeyWord<int>(m_config, "seed");
    m_rng.seed(seed);

    // Initialize move strategy and step sizes - Claude Generated
    m_single_atom_step_size = m_step_size;
    m_orientational_step_size = 0.1; // Radians for orientational moves
    m_last_moved_atom = -1;

    // Convert string parameters to enums - Claude Generated
    if (m_move_strategy == "all_atoms") m_current_strategy = MoveStrategy::ALL_ATOMS;
    else if (m_move_strategy == "single_atom") m_current_strategy = MoveStrategy::SINGLE_ATOM;
    else if (m_move_strategy == "cg_aware") m_current_strategy = MoveStrategy::CG_AWARE;
    else if (m_move_strategy == "chain_segment") m_current_strategy = MoveStrategy::CHAIN_SEGMENT;
    else if (m_move_strategy == "mixed_strategy") m_current_strategy = MoveStrategy::MIXED_STRATEGY;
    else m_current_strategy = MoveStrategy::SINGLE_ATOM; // Default

    if (!silent) {
        CurcumaLogger::info_fmt("Monte Carlo simulation initialized with {} steps", m_total_steps);
        CurcumaLogger::info_fmt("Temperature: {:.1f} K", m_temperature);
        CurcumaLogger::info_fmt("Initial step size: {:.3f} Å", m_step_size);
        CurcumaLogger::info_fmt("Move strategy: {}", m_move_strategy);
        if (m_local_energy_updates) {
            CurcumaLogger::info("Local energy updates enabled for efficiency");
        }
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
            // Analyze system composition - Claude Generated
            analyzeMolecularSystem();

            // Report system type
            if (!m_cg_atom_indices.empty() && !m_atomic_atom_indices.empty()) {
                CurcumaLogger::param("system_type", "Mixed (Atomic + CG)");
                CurcumaLogger::param("cg_atoms", fmt::format("{}", m_cg_atom_indices.size()));
                CurcumaLogger::param("atomic_atoms", fmt::format("{}", m_atomic_atom_indices.size()));
            } else if (!m_cg_atom_indices.empty()) {
                CurcumaLogger::param("system_type", "Coarse-grained");
                CurcumaLogger::param("cg_atoms", fmt::format("{}", m_cg_atom_indices.size()));
            } else {
                CurcumaLogger::param("system_type", "Atomistic");
                CurcumaLogger::param("atomic_atoms", fmt::format("{}", m_atomic_atom_indices.size()));
            }

            // Report detected chains
            if (!m_polymer_chains.empty()) {
                CurcumaLogger::param("polymer_chains", fmt::format("{}", m_polymer_chains.size()));
                CurcumaLogger::param("pivot_moves_available", m_pivot_moves ? "Yes" : "No");
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
    // Determine move strategy dynamically if needed
    if (m_current_strategy == MoveStrategy::MIXED_STRATEGY) {
        m_current_strategy = determineMoveStrategy();
    }

    // Perform move based on strategy
    switch (m_current_strategy) {
        case MoveStrategy::ALL_ATOMS:
            return performLegacyMove(); // Original implementation
        case MoveStrategy::SINGLE_ATOM:
            return performSingleAtomMove();
        case MoveStrategy::CG_AWARE:
            return performCGAwareMove();
        case MoveStrategy::CHAIN_SEGMENT:
            return performChainSegmentMove();
        default:
            return performSingleAtomMove(); // Default fallback
    }
}

bool Casino::performLegacyMove()
{
    // Store current state
    m_trial_molecule = m_molecule;
    double energy_old = m_current_energy;

    // Apply random move (original implementation)
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

bool Casino::performSingleAtomMove()
{
    if (m_molecule.AtomCount() == 0) return false;

    // Select random atom
    int atom_index = selectRandomAtom();

    // Store current state
    m_trial_molecule = m_molecule;
    double energy_old = m_current_energy;

    // Store old position for rollback
    Position old_position = m_molecule.getGeometry().row(atom_index);
    m_last_moved_atom = atom_index;
    m_last_position = old_position;

    // Apply single atom move
    applySingleAtomMove(atom_index);

    // Calculate energy change
    double energy_new;
    if (m_local_energy_updates) {
        Position new_position = m_trial_molecule.getGeometry().row(atom_index);
        double energy_change = calculateLocalEnergyChange(atom_index, old_position, new_position);
        energy_new = energy_old + energy_change;
    } else {
        // Full energy calculation as fallback
        m_energy_calculator->setMolecule(m_trial_molecule.getMolInfo());
        energy_new = calculateEnergy();
    }

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
        // Restore energy calculator to old state if needed
        if (!m_local_energy_updates) {
            m_energy_calculator->setMolecule(m_molecule.getMolInfo());
        }
    }

    return accepted;
}

// ===== SCNP Input File Support - Claude Generated =====

void Casino::setInputFile(const std::string& input_filename)
{
    m_input_filename = input_filename;

    if (!m_silent) {
        CurcumaLogger::info_fmt("Input configuration file set: {}", input_filename);
    }

    // Auto-detect format and load configuration
    if (loadInputConfiguration(input_filename)) {
        if (!m_silent) {
            CurcumaLogger::success("Configuration loaded successfully from input file");
        }
    } else {
        CurcumaLogger::error_fmt("Failed to load configuration from: {}", input_filename);
    }
}

bool Casino::loadInputConfiguration(const std::string& input_filename)
{
    try {
        json input_config;

        // Auto-detect format
        if (ScnpInputParser::isScnpFormat(input_filename)) {
            if (!m_silent) {
                CurcumaLogger::info("Detected SCNP/molsim format - using SCNP parser");
            }
            input_config = m_scnp_parser.parseScnpFile(input_filename);
        } else {
            if (!m_silent) {
                CurcumaLogger::info("Assuming JSON format - using standard parser");
            }
            std::ifstream file(input_filename);
            if (!file.is_open()) {
                CurcumaLogger::error_fmt("Could not open input file: {}", input_filename);
                return false;
            }
            file >> input_config;
            file.close();
        }

        if (input_config.empty()) {
            CurcumaLogger::error("Loaded configuration is empty");
            return false;
        }

        // Merge with existing configuration (input file takes priority)
        json merged_config = MergeJson(m_config, input_config);
        m_config = merged_config;

        // Re-extract parameters from merged configuration
        updateParametersFromConfig();

        if (!m_silent) {
            CurcumaLogger::info("Configuration parameters updated from input file");
            if (ScnpInputParser::isScnpFormat(input_filename)) {
                m_scnp_parser.documentParameters();
            }
        }

        return true;

    } catch (const std::exception& e) {
        CurcumaLogger::error_fmt("Error loading input configuration: {}", e.what());
        return false;
    }
}

void Casino::updateParametersFromConfig()
{
    // Re-extract all parameters from updated configuration
    m_total_steps = Json2KeyWord<int>(m_config, "steps");
    m_temperature = Json2KeyWord<double>(m_config, "temperature");
    m_step_size = Json2KeyWord<double>(m_config, "step_size");
    m_output_frequency = Json2KeyWord<int>(m_config, "output_frequency");
    m_energy_frequency = Json2KeyWord<int>(m_config, "energy_frequency");
    m_output_file = Json2KeyWord<std::string>(m_config, "output_file");
    m_move_type = Json2KeyWord<std::string>(m_config, "move_type");
    m_move_strategy = Json2KeyWord<std::string>(m_config, "move_strategy");
    m_acceptance_target = Json2KeyWord<double>(m_config, "acceptance_target");
    m_adaptive_step = Json2KeyWord<bool>(m_config, "adaptive_step");
    m_verbose = Json2KeyWord<bool>(m_config, "verbose");

    // Enhanced MC parameters
    m_pivot_moves = Json2KeyWord<bool>(m_config, "pivot_moves");
    m_orientational_moves = Json2KeyWord<bool>(m_config, "orientational_moves");
    m_local_energy_updates = Json2KeyWord<bool>(m_config, "local_energy_updates");

    // Enhanced sampling options
    m_bias_potential = Json2KeyWord<bool>(m_config, "bias_potential");
    m_umbrella_sampling = Json2KeyWord<bool>(m_config, "umbrella_sampling");
    m_metadynamics = Json2KeyWord<bool>(m_config, "metadynamics");

    // Update move strategy enum
    if (m_move_strategy == "all_atoms") m_current_strategy = MoveStrategy::ALL_ATOMS;
    else if (m_move_strategy == "single_atom") m_current_strategy = MoveStrategy::SINGLE_ATOM;
    else if (m_move_strategy == "cg_aware") m_current_strategy = MoveStrategy::CG_AWARE;
    else if (m_move_strategy == "chain_segment") m_current_strategy = MoveStrategy::CHAIN_SEGMENT;
    else if (m_move_strategy == "mixed_strategy") m_current_strategy = MoveStrategy::MIXED_STRATEGY;
    else m_current_strategy = MoveStrategy::SINGLE_ATOM;

    // Update step sizes
    m_single_atom_step_size = m_step_size;

    // Handle orientational step size with default value - Claude Generated
    if (m_config.contains("orientational_step_size")) {
        m_orientational_step_size = Json2KeyWord<double>(m_config, "orientational_step_size");
    } else {
        m_orientational_step_size = 0.1; // Default value in radians
    }

    // Handle box length from SCNP files
    if (m_config.contains("box_length")) {
        if (!m_silent) {
            CurcumaLogger::info_fmt("Box length from input file: {}", m_config["box_length"].dump());
            CurcumaLogger::info("Box parameters will be applied to molecule during initialization");
        }
    }
}

// ===== Enhanced Single-Particle Moves Implementation - Claude Generated =====

void Casino::analyzeMolecularSystem()
{
    // Clear previous analysis
    m_cg_atom_indices.clear();
    m_atomic_atom_indices.clear();
    m_cg_shapes.clear();
    m_polymer_chains.clear();
    m_atom_to_chain.clear();

    // Analyze atom types
    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        if (m_molecule.Atom(i).first == CG_ELEMENT) {
            m_cg_atom_indices.push_back(i);
            // Initialize default spherical shape for CG particles
            CGParticleShape shape;
            shape.type = CGParticleShape::SPHERE;
            shape.radii = Eigen::Vector3d(1.0, 1.0, 1.0); // Default 1 Å radius
            shape.orientation = Eigen::Vector3d(0.0, 0.0, 0.0);
            m_cg_shapes.push_back(shape);
        } else {
            m_atomic_atom_indices.push_back(i);
        }
    }

    // Simple chain detection for now (placeholder for more sophisticated topology analysis)
    if (!m_cg_atom_indices.empty() && m_cg_atom_indices.size() > 2) {
        // Assume CG atoms form a linear chain (simplified)
        std::vector<int> chain = m_cg_atom_indices;
        m_polymer_chains.push_back(chain);

        // Map atoms to chain
        for (size_t i = 0; i < chain.size(); ++i) {
            m_atom_to_chain[chain[i]] = 0; // Chain ID 0
        }
    }

    if (!m_silent) {
        CurcumaLogger::info_fmt("System analysis: {} CG atoms, {} atomic atoms, {} chains detected",
                              m_cg_atom_indices.size(), m_atomic_atom_indices.size(), m_polymer_chains.size());
    }
}

MoveStrategy Casino::determineMoveStrategy()
{
    // Adaptive strategy selection based on system composition
    if (!m_cg_atom_indices.empty() && !m_atomic_atom_indices.empty()) {
        // Mixed system: prefer CG-aware moves
        return MoveStrategy::CG_AWARE;
    } else if (!m_cg_atom_indices.empty()) {
        // Pure CG system: check for chains
        if (!m_polymer_chains.empty() && m_pivot_moves) {
            return MoveStrategy::CHAIN_SEGMENT;
        } else {
            return MoveStrategy::CG_AWARE;
        }
    } else {
        // Pure atomic system: use single atom moves
        return MoveStrategy::SINGLE_ATOM;
    }
}

int Casino::selectRandomAtom()
{
    return static_cast<int>(m_uniform_dist(m_rng) * m_molecule.AtomCount());
}

bool Casino::isCGAtom(int atom_index) const
{
    return m_molecule.Atom(atom_index).first == CG_ELEMENT;
}

void Casino::applySingleAtomMove(int atom_index)
{
    Geometry geometry = m_trial_molecule.getGeometry();

    // Apply random displacement
    double dx = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_single_atom_step_size;
    double dy = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_single_atom_step_size;
    double dz = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_single_atom_step_size;

    geometry(atom_index, 0) += dx;
    geometry(atom_index, 1) += dy;
    geometry(atom_index, 2) += dz;

    m_trial_molecule.setGeometry(geometry);
}

double Casino::calculateLocalEnergyChange(int atom_index, const Position& old_pos, const Position& new_pos)
{
    // Placeholder for efficient local energy calculation
    // For now, fall back to full energy calculation
    // TODO: Implement neighbor list and local interaction calculations

    // This would involve:
    // 1. Identify neighbors within cutoff distance
    // 2. Calculate only affected pairwise interactions
    // 3. Return energy difference

    // Fallback to full calculation for safety
    m_energy_calculator->setMolecule(m_trial_molecule.getMolInfo());
    double new_energy = calculateEnergy();

    // Store old geometry temporarily to calculate old energy
    Geometry temp_geometry = m_trial_molecule.getGeometry();
    temp_geometry.row(atom_index) = old_pos;
    m_trial_molecule.setGeometry(temp_geometry);
    m_energy_calculator->setMolecule(m_trial_molecule.getMolInfo());
    double old_energy = calculateEnergy();

    // Restore new geometry
    temp_geometry.row(atom_index) = new_pos;
    m_trial_molecule.setGeometry(temp_geometry);

    return new_energy - old_energy;
}

bool Casino::performCGAwareMove()
{
    if (m_cg_atom_indices.empty()) {
        // Fallback to regular single atom move
        return performSingleAtomMove();
    }

    // Select random CG atom
    int cg_index = static_cast<int>(m_uniform_dist(m_rng) * m_cg_atom_indices.size());
    int atom_index = m_cg_atom_indices[cg_index];

    // Store current state
    m_trial_molecule = m_molecule;
    double energy_old = m_current_energy;

    // Store old position
    Position old_position = m_molecule.getGeometry().row(atom_index);
    m_last_moved_atom = atom_index;
    m_last_position = old_position;

    // Determine move type for CG particle
    bool do_orientational = m_orientational_moves && (m_uniform_dist(m_rng) > 0.5);

    if (do_orientational && !m_cg_shapes[cg_index].isSpherical()) {
        // Apply orientational move
        applyOrientationalMove(atom_index);
    } else {
        // Apply translational move
        applySingleAtomMove(atom_index);
    }

    // Calculate energy change
    double energy_new;
    if (m_local_energy_updates) {
        Position new_position = m_trial_molecule.getGeometry().row(atom_index);
        double energy_change = calculateLocalEnergyChange(atom_index, old_position, new_position);
        energy_new = energy_old + energy_change;
    } else {
        m_energy_calculator->setMolecule(m_trial_molecule.getMolInfo());
        energy_new = calculateEnergy();
    }

    // Accept or reject move
    bool accepted = acceptMove(energy_old, energy_new);

    if (accepted) {
        m_molecule = m_trial_molecule;
        m_current_energy = energy_new;
    } else {
        if (!m_local_energy_updates) {
            m_energy_calculator->setMolecule(m_molecule.getMolInfo());
        }
    }

    return accepted;
}

void Casino::applyOrientationalMove(int cg_atom_index)
{
    // Find CG particle index
    auto it = std::find(m_cg_atom_indices.begin(), m_cg_atom_indices.end(), cg_atom_index);
    if (it == m_cg_atom_indices.end()) return;

    int cg_index = static_cast<int>(it - m_cg_atom_indices.begin());

    // Apply random orientation change (Euler angles)
    double dtheta_x = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_orientational_step_size;
    double dtheta_y = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_orientational_step_size;
    double dtheta_z = (m_uniform_dist(m_rng) - 0.5) * 2.0 * m_orientational_step_size;

    m_cg_shapes[cg_index].orientation[0] += dtheta_x;
    m_cg_shapes[cg_index].orientation[1] += dtheta_y;
    m_cg_shapes[cg_index].orientation[2] += dtheta_z;

    // Note: For now, this only updates the internal orientation state
    // In a full implementation, this would affect the interaction potential
    // based on particle anisotropy
}

bool Casino::performChainSegmentMove()
{
    if (m_polymer_chains.empty()) {
        // Fallback to CG-aware move
        return performCGAwareMove();
    }

    // Placeholder for pivot moves - simplified implementation
    // TODO: Implement real pivot move algorithm

    // For now, just perform a regular single atom move on a chain atom
    int chain_id = static_cast<int>(m_uniform_dist(m_rng) * m_polymer_chains.size());
    const auto& chain = m_polymer_chains[chain_id];

    if (chain.empty()) return false;

    int atom_index_in_chain = static_cast<int>(m_uniform_dist(m_rng) * chain.size());
    int atom_index = chain[atom_index_in_chain];

    // Store current state
    m_trial_molecule = m_molecule;
    double energy_old = m_current_energy;

    // Apply move to selected chain atom
    Position old_position = m_molecule.getGeometry().row(atom_index);
    m_last_moved_atom = atom_index;
    m_last_position = old_position;

    applySingleAtomMove(atom_index);

    // Calculate energy
    double energy_new;
    if (m_local_energy_updates) {
        Position new_position = m_trial_molecule.getGeometry().row(atom_index);
        double energy_change = calculateLocalEnergyChange(atom_index, old_position, new_position);
        energy_new = energy_old + energy_change;
    } else {
        m_energy_calculator->setMolecule(m_trial_molecule.getMolInfo());
        energy_new = calculateEnergy();
    }

    // Accept or reject
    bool accepted = acceptMove(energy_old, energy_new);

    if (accepted) {
        m_molecule = m_trial_molecule;
        m_current_energy = energy_new;
    } else {
        if (!m_local_energy_updates) {
            m_energy_calculator->setMolecule(m_molecule.getMolInfo());
        }
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

    // Adjust step sizes based on current strategy - Claude Generated
    if (ratio > 1.2) {
        // Increase step sizes
        m_step_size *= 1.1;
        m_single_atom_step_size *= 1.1;
        m_orientational_step_size *= 1.1;
    } else if (ratio < 0.8) {
        // Decrease step sizes
        m_step_size *= 0.9;
        m_single_atom_step_size *= 0.9;
        m_orientational_step_size *= 0.9;
    }

    // Apply reasonable bounds
    m_step_size = std::max(0.001, std::min(m_step_size, 1.0));
    m_single_atom_step_size = std::max(0.001, std::min(m_single_atom_step_size, 1.0));
    m_orientational_step_size = std::max(0.001, std::min(m_orientational_step_size, 0.5)); // Radians

    if (m_verbose) {
        CurcumaLogger::info_fmt("Step {} | Acceptance: {:.3f} | Translation: {:.3f} Å | Orientation: {:.3f} rad",
                              m_current_step, current_acceptance, m_single_atom_step_size, m_orientational_step_size);
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
    ParameterRegistry::getInstance().printHelp("casino");

    std::cout << "\nSupported Input Formats:" << std::endl;
    std::cout << "- XYZ, VTF, MOL2, SDF, PDB (molecular structures)" << std::endl;
    std::cout << "- JSON (Curcuma native), SCNP/molsim (auto-converted)" << std::endl;
    std::cout << "\nSCNP Parameter Mapping (common parameters):" << std::endl;
    std::cout << "  txtitle    → simulation_title  (descriptive title)" << std::endl;
    std::cout << "  temp       → temperature        (temperature in K)" << std::endl;
    std::cout << "  nstep1     → steps              (number of MC steps)" << std::endl;
    std::cout << "  nfreq      → output_frequency   (trajectory output frequency)" << std::endl;
    std::cout << "  seed       → seed               (random number seed)" << std::endl;
}

nlohmann::json Casino::WriteRestartInformation()
{
    json restart_data;

    restart_data["metadata"] = {
        {"version", "1.0"},
        {"creation_time", std::time(nullptr)},
        {"simulation_type", "casino_mc"},
        {"total_steps", m_total_steps},
        {"completed_steps", m_current_step}
    };

    restart_data["simulation_state"] = {
        {"current_step", m_current_step},
        {"current_energy", m_current_energy},
        {"total_moves", m_total_moves},
        {"accepted_moves", m_accepted_moves},
        {"step_size", m_step_size},
        {"single_atom_step_size", m_single_atom_step_size},
        {"orientational_step_size", m_orientational_step_size},
        {"last_moved_atom", m_last_moved_atom}
    };

    restart_data["parameters"] = {
        {"temperature", m_temperature},
        {"move_strategy", m_move_strategy},
        {"move_type", m_move_type},
        {"acceptance_target", m_acceptance_target},
        {"adaptive_step", m_adaptive_step},
        {"pivot_moves", m_pivot_moves},
        {"orientational_moves", m_orientational_moves},
        {"local_energy_updates", m_local_energy_updates},
        {"verbose", m_verbose}
    };

    restart_data["random_state"] = {
        {"rng_state", m_rng()},
        {"seed_backup", 42}
    };

    restart_data["statistics"] = {
        {"energy_trajectory", m_energy_trajectory},
        {"acceptance_history", m_acceptance_history}
    };

    if (!m_cg_atom_indices.empty()) {
        restart_data["cg_system"] = exportCGShapes();
        restart_data["cg_system"]["cg_atom_indices"] = m_cg_atom_indices;
        restart_data["cg_system"]["atomic_atom_indices"] = m_atomic_atom_indices;
    }

    if (!m_polymer_chains.empty()) {
        restart_data["polymer_system"] = {
            {"chains", m_polymer_chains},
            {"atom_to_chain_map", json::object()}
        };
        for (const auto& [atom_idx, chain_id] : m_atom_to_chain) {
            restart_data["polymer_system"]["atom_to_chain_map"][std::to_string(atom_idx)] = chain_id;
        }
    }

    restart_data["molecule"] = {
        {"natoms", m_molecule.AtomCount()},
        {"geometry", json::array()},
        {"elements", json::array()}
    };

    for (int i = 0; i < m_molecule.AtomCount(); ++i) {
        json atom_data;
        Position pos = m_molecule.Atom(i).second;
        atom_data["x"] = pos[0];
        atom_data["y"] = pos[1];
        atom_data["z"] = pos[2];
        restart_data["molecule"]["geometry"].push_back(atom_data);
        restart_data["molecule"]["elements"].push_back(m_molecule.Atom(i).first);
    }

    restart_data["file_information"] = {
        {"input_filename", m_input_filename},
        {"output_file", m_output_file},
        {"checkpoint_file", m_checkpoint_file}
    };

    return restart_data;
}

bool Casino::LoadRestartInformation()
{
    if (!validateRestartFile(m_checkpoint_file)) {
        if (m_verbose) {
            std::cout << "[CASINO] No valid restart file found: " << m_checkpoint_file << std::endl;
        }
        return false;
    }

    try {
        std::ifstream restart_file(m_checkpoint_file);
        if (!restart_file.is_open()) {
            if (m_verbose) {
                std::cout << "[CASINO] Cannot open restart file: " << m_checkpoint_file << std::endl;
            }
            return false;
        }

        json restart_data;
        restart_file >> restart_data;
        restart_file.close();

        if (!restart_data.contains("metadata") ||
            restart_data["metadata"]["simulation_type"] != "casino_mc") {
            if (m_verbose) {
                std::cout << "[CASINO] Invalid restart file format" << std::endl;
            }
            return false;
        }

        if (m_verbose) {
            std::cout << "[CASINO] Loading restart from step "
                      << restart_data["simulation_state"]["current_step"] << std::endl;
        }

        m_current_step = restart_data["simulation_state"]["current_step"];
        m_current_energy = restart_data["simulation_state"]["current_energy"];
        m_total_moves = restart_data["simulation_state"]["total_moves"];
        m_accepted_moves = restart_data["simulation_state"]["accepted_moves"];
        m_step_size = restart_data["simulation_state"]["step_size"];

        if (restart_data["simulation_state"].contains("single_atom_step_size")) {
            m_single_atom_step_size = restart_data["simulation_state"]["single_atom_step_size"];
        }
        if (restart_data["simulation_state"].contains("orientational_step_size")) {
            m_orientational_step_size = restart_data["simulation_state"]["orientational_step_size"];
        }
        if (restart_data["simulation_state"].contains("last_moved_atom")) {
            m_last_moved_atom = restart_data["simulation_state"]["last_moved_atom"];
        }

        if (restart_data.contains("statistics")) {
            if (restart_data["statistics"].contains("energy_trajectory")) {
                m_energy_trajectory = restart_data["statistics"]["energy_trajectory"].get<std::vector<double>>();
            }
            if (restart_data["statistics"].contains("acceptance_history")) {
                m_acceptance_history = restart_data["statistics"]["acceptance_history"].get<std::vector<double>>();
            }
        }

        if (restart_data.contains("molecule")) {
            auto mol_data = restart_data["molecule"];
            int natoms = mol_data["natoms"];

            if (natoms != m_molecule.AtomCount()) {
                if (m_verbose) {
                    std::cout << "[CASINO] Warning: Atom count mismatch in restart file" << std::endl;
                }
                return false;
            }

            Geometry new_geometry(natoms, 3);
            for (int i = 0; i < natoms; ++i) {
                auto atom_pos = mol_data["geometry"][i];
                new_geometry(i, 0) = atom_pos["x"];
                new_geometry(i, 1) = atom_pos["y"];
                new_geometry(i, 2) = atom_pos["z"];
            }
            m_molecule.setGeometry(new_geometry);
        }

        if (restart_data.contains("cg_system")) {
            loadCGShapes(restart_data["cg_system"]);
            if (restart_data["cg_system"].contains("cg_atom_indices")) {
                m_cg_atom_indices = restart_data["cg_system"]["cg_atom_indices"].get<std::vector<int>>();
            }
            if (restart_data["cg_system"].contains("atomic_atom_indices")) {
                m_atomic_atom_indices = restart_data["cg_system"]["atomic_atom_indices"].get<std::vector<int>>();
            }
        }

        if (restart_data.contains("polymer_system")) {
            m_polymer_chains = restart_data["polymer_system"]["chains"].get<std::vector<std::vector<int>>>();
            m_atom_to_chain.clear();
            for (const auto& [atom_str, chain_id] : restart_data["polymer_system"]["atom_to_chain_map"].items()) {
                m_atom_to_chain[std::stoi(atom_str)] = chain_id;
            }
        }

        if (restart_data.contains("random_state")) {
            if (restart_data["random_state"].contains("seed_backup")) {
                int seed = restart_data["random_state"]["seed_backup"];
                m_rng.seed(seed + m_current_step);
            }
        }

        if (m_verbose) {
            std::cout << "[CASINO] Successfully loaded restart data:" << std::endl;
            std::cout << "         Current step: " << m_current_step << "/" << m_total_steps << std::endl;
            std::cout << "         Current energy: " << m_current_energy << " Ha" << std::endl;
            std::cout << "         Acceptance ratio: "
                      << (m_total_moves > 0 ? m_accepted_moves / double(m_total_moves) : 0.0) << std::endl;
        }

        return true;
    }
    catch (const std::exception& e) {
        if (m_verbose) {
            std::cout << "[CASINO] Error loading restart file: " << e.what() << std::endl;
        }
        return false;
    }
}

void Casino::performCheckpoint()
{
    try {
        json restart_data = WriteRestartInformation();

        std::string checkpoint_path = m_checkpoint_file;
        if (m_max_checkpoints > 1) {
            checkpoint_path = m_checkpoint_file + ".backup." + std::to_string(m_current_step);
        }

        std::ofstream checkpoint_file(checkpoint_path);
        if (!checkpoint_file.is_open()) {
            if (m_verbose) {
                std::cout << "[CASINO] Warning: Cannot write checkpoint to " << checkpoint_path << std::endl;
            }
            return;
        }

        checkpoint_file << restart_data.dump(4) << std::endl;
        checkpoint_file.close();

        if (m_max_checkpoints > 1) {
            std::filesystem::copy_file(checkpoint_path, m_checkpoint_file,
                                       std::filesystem::copy_options::overwrite_existing);
        }

        if (m_verbose) {
            std::cout << "[CASINO] Checkpoint saved at step " << m_current_step << std::endl;
        }

        if (m_max_checkpoints > 1) {
            cleanupOldCheckpoints();
        }
    }
    catch (const std::exception& e) {
        if (m_verbose) {
            std::cout << "[CASINO] Error during checkpointing: " << e.what() << std::endl;
        }
    }
}

void Casino::cleanupOldCheckpoints()
{
    try {
        std::vector<std::string> checkpoint_files;
        std::string base_name = m_checkpoint_file + ".backup.";

        for (const auto& entry : std::filesystem::directory_iterator(std::filesystem::current_path())) {
            if (entry.path().filename().string().find(base_name) == 0) {
                checkpoint_files.push_back(entry.path().string());
            }
        }

        if (checkpoint_files.size() > static_cast<size_t>(m_max_checkpoints)) {
            std::sort(checkpoint_files.begin(), checkpoint_files.end());

            size_t files_to_remove = checkpoint_files.size() - m_max_checkpoints;
            for (size_t i = 0; i < files_to_remove; ++i) {
                std::filesystem::remove(checkpoint_files[i]);
                if (m_verbose) {
                    std::cout << "[CASINO] Removed old checkpoint: " << checkpoint_files[i] << std::endl;
                }
            }
        }
    }
    catch (const std::exception& e) {
        if (m_verbose) {
            std::cout << "[CASINO] Error cleaning old checkpoints: " << e.what() << std::endl;
        }
    }
}

bool Casino::detectAndLoadRestart()
{
    if (!m_auto_restart) {
        return false;
    }

    std::vector<std::string> potential_restart_files = {
        m_checkpoint_file,
        m_checkpoint_file + ".backup." + std::to_string(m_current_step),
        "casino_checkpoint.json",
        "casino_checkpoint_latest.json"
    };

    for (const auto& file : potential_restart_files) {
        if (validateRestartFile(file)) {
            if (m_verbose) {
                std::cout << "[CASINO] Found restart file: " << file << std::endl;
            }

            std::string backup_checkpoint_file = m_checkpoint_file;
            m_checkpoint_file = file;
            bool success = LoadRestartInformation();

            if (success) {
                m_checkpoint_file = backup_checkpoint_file;
                return true;
            } else {
                m_checkpoint_file = backup_checkpoint_file;
                if (m_verbose) {
                    std::cout << "[CASINO] Failed to load restart from: " << file << std::endl;
                }
            }
        }
    }

    if (m_verbose) {
        std::cout << "[CASINO] No valid restart files found" << std::endl;
    }
    return false;
}

bool Casino::validateRestartFile(const std::string& filename)
{
    try {
        if (!std::filesystem::exists(filename) || std::filesystem::file_size(filename) == 0) {
            return false;
        }

        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }

        json restart_data;
        file >> restart_data;
        file.close();

        if (!restart_data.contains("metadata") ||
            !restart_data["metadata"].contains("simulation_type") ||
            restart_data["metadata"]["simulation_type"] != "casino_mc") {
            return false;
        }

        if (!restart_data.contains("simulation_state") ||
            !restart_data["simulation_state"].contains("current_step") ||
            !restart_data["simulation_state"].contains("current_energy")) {
            return false;
        }

        if (!restart_data.contains("molecule") ||
            !restart_data["molecule"].contains("natoms") ||
            !restart_data["molecule"].contains("geometry")) {
            return false;
        }

        return true;
    }
    catch (const std::exception&) {
        return false;
    }
}

bool Casino::loadFromRestartFile(const std::string& filename)
{
    if (!validateRestartFile(filename)) {
        return false;
    }

    std::string backup_checkpoint_file = m_checkpoint_file;
    m_checkpoint_file = filename;
    bool success = LoadRestartInformation();
    m_checkpoint_file = backup_checkpoint_file;

    return success;
}

void Casino::handleTrajectoryRestart()
{
    if (!m_append_trajectory) {
        return;
    }

    try {
        if (std::filesystem::exists(m_output_file)) {
            if (m_verbose) {
                std::cout << "[CASINO] Existing trajectory found: " << m_output_file << std::endl;
            }

            if (m_current_step == 0) {
                if (m_verbose) {
                    std::cout << "[CASINO] Starting from step 0, overwriting trajectory" << std::endl;
                }
                std::filesystem::remove(m_output_file);
                return;
            }

            std::string backup_name = m_output_file + ".backup." + std::to_string(std::time(nullptr));
            std::filesystem::copy_file(m_output_file, backup_name);

            if (m_verbose) {
                std::cout << "[CASINO] Trajectory backup created: " << backup_name << std::endl;
                std::cout << "[CASINO] Continuing trajectory from step " << m_current_step << std::endl;
            }
        }
    }
    catch (const std::exception& e) {
        if (m_verbose) {
            std::cout << "[CASINO] Error handling trajectory restart: " << e.what() << std::endl;
        }
    }
}

json Casino::exportCGShapes() const
{
    json cg_shapes_data;
    cg_shapes_data["shapes"] = json::array();

    for (size_t i = 0; i < m_cg_shapes.size(); ++i) {
        const auto& shape = m_cg_shapes[i];
        json shape_data;

        shape_data["radii"] = {shape.radii[0], shape.radii[1], shape.radii[2]};
        shape_data["orientation"] = {shape.orientation[0], shape.orientation[1], shape.orientation[2]};

        switch (shape.type) {
            case CGParticleShape::SPHERE:
                shape_data["type"] = "sphere";
                break;
            case CGParticleShape::ELLIPSOID:
                shape_data["type"] = "ellipsoid";
                break;
            case CGParticleShape::CYLINDER:
                shape_data["type"] = "cylinder";
                break;
            default:
                shape_data["type"] = "sphere";
        }

        shape_data["is_spherical"] = shape.isSpherical();
        cg_shapes_data["shapes"].push_back(shape_data);
    }

    cg_shapes_data["num_cg_particles"] = m_cg_shapes.size();
    return cg_shapes_data;
}

void Casino::loadCGShapes(const json& cg_data)
{
    try {
        if (!cg_data.contains("shapes")) {
            return;
        }

        m_cg_shapes.clear();
        auto shapes_array = cg_data["shapes"];

        for (const auto& shape_data : shapes_array) {
            CGParticleShape shape;

            if (shape_data.contains("radii") && shape_data["radii"].is_array()) {
                auto radii = shape_data["radii"];
                shape.radii = Eigen::Vector3d(radii[0], radii[1], radii[2]);
            }

            if (shape_data.contains("orientation") && shape_data["orientation"].is_array()) {
                auto orientation = shape_data["orientation"];
                shape.orientation = Eigen::Vector3d(orientation[0], orientation[1], orientation[2]);
            }

            if (shape_data.contains("type")) {
                std::string type_str = shape_data["type"];
                if (type_str == "sphere") {
                    shape.type = CGParticleShape::SPHERE;
                } else if (type_str == "ellipsoid") {
                    shape.type = CGParticleShape::ELLIPSOID;
                } else if (type_str == "cylinder") {
                    shape.type = CGParticleShape::CYLINDER;
                } else {
                    shape.type = CGParticleShape::SPHERE;
                }
            }

            m_cg_shapes.push_back(shape);
        }

        if (m_verbose && !m_cg_shapes.empty()) {
            std::cout << "[CASINO] Loaded " << m_cg_shapes.size() << " CG particle shapes from restart" << std::endl;
        }
    }
    catch (const std::exception& e) {
        if (m_verbose) {
            std::cout << "[CASINO] Error loading CG shapes: " << e.what() << std::endl;
        }
    }
}