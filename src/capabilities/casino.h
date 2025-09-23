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

#include "curcumamethod.h"

// Default configuration for Casino Monte Carlo simulation - Claude Generated
static const json CasinoJson = {
    { "steps", 10000 },                    // Number of MC steps
    { "temperature", 300.0 },              // Temperature in Kelvin
    { "step_size", 0.1 },                  // Maximum displacement per move (Angstrom)
    { "seed", 42 },                        // Random seed for reproducibility
    { "output_frequency", 100 },           // Frequency for trajectory output
    { "energy_frequency", 10 },            // Frequency for energy output
    { "output_file", "casino_trajectory.xyz" }, // Output trajectory file
    { "move_type", "translation" },        // Move type: "translation", "rotation", "mixed"
    { "acceptance_target", 0.5 },          // Target acceptance ratio
    { "adaptive_step", true },             // Adaptive step size adjustment
    { "bias_potential", false },           // Enable bias potential
    { "umbrella_sampling", false },        // Enable umbrella sampling
    { "metadynamics", false },             // Enable metadynamics
    { "verbose", true }                    // Detailed output
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

    /*! \brief Print help for Casino Monte Carlo options */
    void printHelp() const override;

    /*! \brief Get final molecule after simulation */
    const Molecule& getFinalMolecule() const { return m_molecule; }

    /*! \brief Get energy trajectory */
    const std::vector<double>& getEnergyTrajectory() const { return m_energy_trajectory; }

    /*! \brief Get acceptance statistics */
    double getAcceptanceRatio() const { return m_accepted_moves / static_cast<double>(m_total_moves); }

    // CurcumaMethod interface implementation
    nlohmann::json WriteRestartInformation() override { return json{}; }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return {"Casino"}; }
    void ReadControlFile() override {}
    void LoadControlJson() override {}

private:
    /*! \brief Initialize Monte Carlo simulation */
    bool Initialise() override;

    /*! \brief Perform single Monte Carlo move */
    bool performMove();

    /*! \brief Calculate energy of current configuration */
    double calculateEnergy();

    /*! \brief Apply random displacement to molecule */
    void applyRandomMove();

    /*! \brief Accept or reject move based on Metropolis criterion */
    bool acceptMove(double energy_old, double energy_new);

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
    json m_config;

    // Simulation parameters
    int m_total_steps;
    double m_temperature;
    double m_step_size;
    int m_output_frequency;
    int m_energy_frequency;
    std::string m_output_file;
    std::string m_move_type;
    double m_acceptance_target;
    bool m_adaptive_step;
    bool m_verbose;

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

    // Physical constants
    static constexpr double kB = 8.314462618e-3; // Boltzmann constant in kJ/(mol·K)
};