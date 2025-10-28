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

#pragma once

#include <chrono>
#include <random>
#include <vector>

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/parameter_macros.h"
#include "scnp_parser.h"

#include "curcumamethod.h"

/* Claude Generated 2025: Casino Parameter Registry - replaces static CasinoJson */
BEGIN_PARAMETER_DEFINITION(casino)
    // Basic simulation parameters
    PARAM(steps, Int, 10000, "Number of Monte Carlo steps", "Basic", {})
    PARAM(temperature, Double, 300.0, "Temperature in Kelvin", "Basic", {})
    PARAM(step_size, Double, 0.1, "Maximum displacement per move (Angstrom)", "Basic", {})
    PARAM(method, String, "uff", "Computational method for energy/gradient", "Basic", {})
    PARAM(seed, Int, 42, "Random seed for reproducibility", "Basic", {})

    // Output options
    PARAM(output_file, String, "casino_trajectory.xyz", "Trajectory output file", "Output", {})
    PARAM(output_frequency, Int, 100, "Trajectory output frequency", "Output", {})
    PARAM(energy_frequency, Int, 10, "Energy output frequency", "Output", {})

    // Move configuration
    PARAM(move_type, String, "mixed", "Move type: translation|rotation|orientational|pivot|mixed", "Algorithm", {})
    PARAM(move_strategy, String, "single_atom", "Move strategy: all_atoms|single_atom|cg_aware|chain_segment|mixed_strategy", "Algorithm", {})

    // Adaptive sampling
    PARAM(acceptance_target, Double, 0.5, "Target acceptance ratio", "Algorithm", {})
    PARAM(adaptive_step, Bool, true, "Adaptive step size adjustment", "Algorithm", {})

    // Advanced features
    PARAM(local_energy_updates, Bool, true, "Use local energy updates for efficiency", "Advanced", {})
    PARAM(pivot_moves, Bool, false, "Enable pivot moves for polymer chains", "Advanced", {})
    PARAM(orientational_moves, Bool, false, "Enable orientational moves for CG particles", "Advanced", {})
    PARAM(bias_potential, Bool, false, "Enable bias potential", "Advanced", {})
    PARAM(umbrella_sampling, Bool, false, "Enable umbrella sampling", "Advanced", {})
    PARAM(metadynamics, Bool, false, "Enable metadynamics", "Advanced", {})
    PARAM(verbose, Bool, true, "Detailed output", "Output", {})
END_PARAMETER_DEFINITION

// Enhanced move strategies for MC simulation - Claude Generated
enum class MoveStrategy {
    ALL_ATOMS,      // Move all atoms simultaneously (legacy behavior)
    SINGLE_ATOM,    // Move single randomly selected atom
    CG_AWARE,       // CG-particle specific moves with shape awareness
    CHAIN_SEGMENT,  // Polymer chain segment moves
    MIXED_STRATEGY  // Adaptive combination of strategies
};

// Enhanced move types for different particle types - Claude Generated
enum class MoveType {
    TRANSLATION,    // Standard translational displacement
    ROTATION,       // Whole molecule rotation
    ORIENTATIONAL,  // CG particle orientation (Euler angles)
    PIVOT,          // Polymer chain pivot moves
    MIXED           // Combination of move types
};

// CG particle shape information - Claude Generated
struct CGParticleShape {
    Eigen::Vector3d radii;        // x,y,z radii (sphere: all equal)
    Eigen::Vector3d orientation;  // Euler angles for ellipsoids
    enum ShapeType { SPHERE, ELLIPSOID, CYLINDER } type;

    CGParticleShape() : radii(1.0, 1.0, 1.0), orientation(0.0, 0.0, 0.0), type(SPHERE) {}
    bool isSpherical() const { return std::abs(radii[0] - radii[1]) < 1e-6 && std::abs(radii[1] - radii[2]) < 1e-6; }
};

// Default configuration for Casino Monte Carlo simulation - Claude Generated
static const json CasinoJson = {
    { "steps", 10000 },                    // Number of MC steps
    { "temperature", 300.0 },              // Temperature in Kelvin
    { "step_size", 0.1 },                  // Maximum displacement per move (Angstrom)
    { "seed", 42 },                        // Random seed for reproducibility
    { "output_frequency", 100 },           // Frequency for trajectory output
    { "energy_frequency", 10 },            // Frequency for energy output
    { "output_file", "casino_trajectory.xyz" }, // Output trajectory file
    { "move_type", "mixed" },              // Move type: "translation", "rotation", "orientational", "pivot", "mixed"
    { "move_strategy", "single_atom" },    // Move strategy: "all_atoms", "single_atom", "cg_aware", "chain_segment", "mixed_strategy"
    { "acceptance_target", 0.5 },          // Target acceptance ratio
    { "adaptive_step", true },             // Adaptive step size adjustment
    { "pivot_moves", false },              // Enable pivot moves for polymer chains
    { "orientational_moves", false },      // Enable orientational moves for CG particles
    { "local_energy_updates", true },      // Use local energy updates for efficiency
    { "bias_potential", false },           // Enable bias potential
    { "umbrella_sampling", false },        // Enable umbrella sampling
    { "metadynamics", false },             // Enable metadynamics
    { "verbose", true },                   // Detailed output
    // Restart and checkpointing options - Claude Generated
    { "auto_restart", false },             // Automatically detect and load restart files
    { "checkpoint_frequency", 1000 },      // Steps between checkpoints
    { "checkpoint_file", "casino_checkpoint.json" }, // Checkpoint filename
    { "max_checkpoints", 5 },              // Keep last N checkpoints
    { "append_trajectory", true },         // Append to existing trajectory on restart
    { "checkpoint_on_signal", true }       // Checkpoint on SIGTERM/SIGINT
};

/*! \brief Casino - Monte Carlo simulation for molecular systems - Claude Generated
 *
 * Universal Monte Carlo simulation capability that works with:
 * - Atomistic molecules (XYZ, MOL2, SDF, PDB)
 * - Coarse-grained systems (VTF)
 * - Any force field or quantum method supported by Curcuma
 *
 * Features:
 * - Canonical (NVT) ensemble simulation
 * - Multiple move types: translation, rotation, mixed
 * - Adaptive step size for optimal acceptance ratio
 * - Enhanced sampling methods: umbrella sampling, metadynamics
 * - Time-series analysis and trajectory output
 * - Energy statistics and convergence monitoring
 */
class Casino : public CurcumaMethod
{
public:
    Casino(const json& controller, bool silent);
    ~Casino();

    /*! \brief Start Monte Carlo simulation */
    void start() override;

    /*! \brief Set input molecule */
    inline void setMolecule(const Molecule& molecule) { m_molecule = molecule; }

    /*! \brief Set input filename (any supported format) */
    inline void setFileName(const std::string& filename) { m_filename = filename; }

    /*! \brief Set input configuration file (SCNP/JSON auto-detection) */
    void setInputFile(const std::string& input_filename);

    /*! \brief Load configuration from SCNP or JSON input file */
    bool loadInputConfiguration(const std::string& input_filename);

    /*! \brief Print help for Casino Monte Carlo options */
    void printHelp() const override;

    /*! \brief Get final molecule after simulation */
    const Molecule& getFinalMolecule() const { return m_molecule; }

    /*! \brief Get energy trajectory */
    const std::vector<double>& getEnergyTrajectory() const { return m_energy_trajectory; }

    /*! \brief Get acceptance statistics */
    double getAcceptanceRatio() const { return m_accepted_moves / static_cast<double>(m_total_moves); }

    // CurcumaMethod interface implementation
    nlohmann::json WriteRestartInformation() override;
    bool LoadRestartInformation() override;
    StringList MethodName() const override { return {"Casino"}; }
    void ReadControlFile() override {}
    void LoadControlJson() override {}

    // ===== Restart System Methods - Claude Generated =====
    /*! \brief Perform automatic checkpoint during simulation */
    void performCheckpoint();

    /*! \brief Detect and load restart files automatically */
    bool detectAndLoadRestart();

    /*! \brief Validate restart file integrity */
    bool validateRestartFile(const std::string& filename);

    /*! \brief Load restart data from specific file */
    bool loadFromRestartFile(const std::string& filename);

    /*! \brief Handle trajectory file continuation on restart */
    void handleTrajectoryRestart();

    /*! \brief Setup signal handlers for graceful shutdown */
    void setupSignalHandlers();

    /*! \brief Export CG particle shapes to JSON */
    json exportCGShapes() const;

    /*! \brief Load CG particle shapes from JSON */
    void loadCGShapes(const json& cg_data);

    /*! \brief Print restart system help */
    void printRestartHelp() const;

    /*! \brief Clean up old checkpoint files based on max_checkpoints setting */
    void cleanupOldCheckpoints();

private:
    /*! \brief Initialize Monte Carlo simulation */
    bool Initialise() override;

    /*! \brief Perform single Monte Carlo move */
    bool performMove();

    /*! \brief Perform legacy all-atoms move (original implementation) */
    bool performLegacyMove();

    /*! \brief Analyze molecular system for CG atoms and chains */
    void analyzeMolecularSystem();

    /*! \brief Calculate energy of current configuration */
    double calculateEnergy();

    /*! \brief Apply random displacement to molecule */
    void applyRandomMove();

    /*! \brief Accept or reject move based on Metropolis criterion */
    bool acceptMove(double energy_old, double energy_new);

    // ===== Enhanced Single-Particle Moves - Claude Generated =====
    /*! \brief Perform single-atom Monte Carlo move */
    bool performSingleAtomMove();

    /*! \brief Perform CG-aware Monte Carlo move */
    bool performCGAwareMove();

    /*! \brief Perform chain segment move */
    bool performChainSegmentMove();

    /*! \brief Calculate local energy change for single atom move */
    double calculateLocalEnergyChange(int atom_index, const Position& old_pos, const Position& new_pos);

    /*! \brief Apply single atom displacement */
    void applySingleAtomMove(int atom_index);

    /*! \brief Apply orientational move to CG particle */
    void applyOrientationalMove(int cg_atom_index);

    /*! \brief Determine move strategy based on system type */
    MoveStrategy determineMoveStrategy();

    /*! \brief Select random atom for single-particle move */
    int selectRandomAtom();

    /*! \brief Check if atom is part of CG system */
    bool isCGAtom(int atom_index) const;

    /*! \brief Update parameters from loaded configuration */
    void updateParametersFromConfig();

    /*! \brief Update step size for optimal acceptance */
    void updateStepSize();

    /*! \brief Write trajectory frame */
    void writeTrajectoryFrame(int step);

    /*! \brief Output energy and statistics */
    void outputStatistics(int step);

    /*! \brief Calculate bias potential energy */
    double calculateBiasPotential();

    /*! \brief Update metadynamics history */
    void updateMetadynamics();

    // Core simulation data
    Molecule m_molecule;
    Molecule m_trial_molecule;
    EnergyCalculator* m_energy_calculator;
    std::string m_filename;
    std::string m_input_filename;  // SCNP/JSON input file
    json m_config;
    ScnpInputParser m_scnp_parser; // SCNP input parser

    // Simulation parameters
    int m_total_steps;
    double m_temperature;
    double m_step_size;
    int m_output_frequency;
    int m_energy_frequency;
    std::string m_output_file;
    std::string m_move_type;
    std::string m_move_strategy;
    double m_acceptance_target;
    bool m_adaptive_step;
    bool m_verbose;

    // Enhanced MC parameters - Claude Generated
    bool m_pivot_moves;
    bool m_orientational_moves;
    bool m_local_energy_updates;
    MoveStrategy m_current_strategy;
    MoveType m_current_move_type;

    // Restart and checkpointing parameters - Claude Generated
    bool m_auto_restart;
    int m_checkpoint_frequency;
    std::string m_checkpoint_file;
    int m_max_checkpoints;
    bool m_append_trajectory;
    bool m_checkpoint_on_signal;

    // Random number generation
    std::mt19937 m_rng;
    std::uniform_real_distribution<double> m_uniform_dist;
    std::normal_distribution<double> m_normal_dist;

    // Statistics
    int m_current_step;
    int m_total_moves;
    int m_accepted_moves;
    double m_current_energy;
    std::vector<double> m_energy_trajectory;
    std::vector<double> m_acceptance_history;

    // Enhanced sampling
    bool m_bias_potential;
    bool m_umbrella_sampling;
    bool m_metadynamics;
    std::vector<Geometry> m_metadynamics_history;

    // CG and chain analysis - Claude Generated
    std::vector<CGParticleShape> m_cg_shapes;         // Shape information for CG particles
    std::vector<int> m_cg_atom_indices;               // Indices of CG atoms
    std::vector<int> m_atomic_atom_indices;           // Indices of atomic atoms
    std::vector<std::vector<int>> m_polymer_chains;   // Detected polymer chains
    std::map<int, int> m_atom_to_chain;               // Map atom index to chain ID

    // Move selection and efficiency - Claude Generated
    int m_last_moved_atom;                            // Last moved atom for energy updates
    Position m_last_position;                         // Last position for rollback
    double m_single_atom_step_size;                   // Step size for single atom moves
    double m_orientational_step_size;                 // Step size for orientational moves

    // Physical constants
    static constexpr double kB = 8.314462618e-3; // Boltzmann constant in kJ/(mol·K)
};