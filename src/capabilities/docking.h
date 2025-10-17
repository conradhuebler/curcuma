/*
 * <Docking tool for structures. >
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/capabilities/curcumaopt.h"
#include "src/capabilities/optimiser/LevMarDocking.h"

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include "src/tools/general.h"
#include "src/tools/geometry.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <map>
#include <thread>

#include "json.hpp"
using json = nlohmann::json;

#include "curcumamethod.h"

class DockThread : public CxxThread {
public:
    inline DockThread(const Molecule& host, const Molecule& guest)
        : m_host(host)
        , m_guest(guest)
    {
    }
    ~DockThread() = default;

    inline void setPosition(const Position& position) { m_position = position; }
    inline void setRotation(const Position& rotation) { m_rotation = rotation; }

    inline int execute() override
    {
        std::pair<Position, Position> pair = OptimiseAnchor(&m_host, m_guest, m_position, m_rotation);
        m_last_position = pair.first;
        m_last_rotation = pair.second;
        return 0;
    }

    inline Position InitialPosition() const { return m_position; }
    inline Position InitialRotation() const { return m_rotation; }
    inline Position LastPosition() const { return m_last_position; }
    inline Position LastRotation() const { return m_last_rotation; }

private:
    Position m_position, m_rotation, m_last_position, m_last_rotation;
    Molecule m_host, m_guest;
};

/* Claude Generated 2025: Docking Parameter Registry - replaces static DockingJson */
BEGIN_PARAMETER_DEFINITION(docking)
    PARAM(host, String, "", "Host molecule file.", "Input", {})
    PARAM(guest, String, "", "Guest molecule file.", "Input", {})
    PARAM(complex, String, "", "Pre-assembled complex file.", "Input", {"Complex"})
    PARAM(pos_x, Double, 0.0, "X-position for docking box center.", "Grid", {"Pos_X"})
    PARAM(pos_y, Double, 0.0, "Y-position for docking box center.", "Grid", {"Pos_Y"})
    PARAM(pos_z, Double, 0.0, "Z-position for docking box center.", "Grid", {"Pos_Z"})
    PARAM(step_x, Int, 10, "Number of steps in X direction.", "Grid", {"Step_X"})
    PARAM(step_y, Int, 10, "Number of steps in Y direction.", "Grid", {"Step_Y"})
    PARAM(step_z, Int, 10, "Number of steps in Z direction.", "Grid", {"Step_z"})
    PARAM(auto_pos, Bool, true, "Automatically determine docking box position.", "Grid", {"AutoPos"})
    PARAM(scaling, Double, 1.5, "Scaling factor for docking box.", "Grid", {})
    PARAM(filter, Bool, true, "Filter docking results.", "PostProcessing", {"Filter"})
    PARAM(post_opt, Bool, true, "Post-optimization of results.", "PostProcessing", {"PostOpt"})
    PARAM(no_opt, Bool, false, "Skip optimization step.", "PostProcessing", {"NoOpt"})
    PARAM(centroid_max_distance, Double, 1e5, "Maximum centroid distance.", "Algorithm", {"CentroidMaxDistance"})
    PARAM(centroid_tol_distance, Double, 0.1, "Centroid tolerance distance.", "Algorithm", {"CentroidTolDis"})
    PARAM(centroid_rot_distance, Double, 0.1, "Centroid rotation tolerance.", "Algorithm", {"RotationTolDis"})
    PARAM(energy_threshold, Double, 200.0, "Energy threshold for filtering.", "Filtering", {"EnergyThreshold"})
    PARAM(cycles, Int, 1, "Number of docking cycles.", "Execution", {"Cycles"})
    PARAM(threads, Int, 1, "Number of threads for the main process.", "Performance", {"Threads"})
    PARAM(docking_threads, Int, 1, "Number of threads for docking subprocesses.", "Performance", {"DockingThreads"})
    PARAM(charge, Int, 0, "Total charge of complex.", "Molecule", {"Charge"})
    PARAM(rmsd_method, String, "hybrid", "RMSD calculation method.", "Analysis", {"RMSDMethod"})
    PARAM(rmsd_threads, Int, 1, "Number of RMSD threads.", "Performance", {"RMSDThreads"})
    PARAM(rmsd_element, Int, 7, "Element type for RMSD (atomic number).", "Analysis", {"RMSDElement"})
END_PARAMETER_DEFINITION

class Docking : public CurcumaMethod {
public:
    /**
     * @brief Constructor with JSON configuration (backward compatible)
     * Claude Generated: Phase 4 - ConfigManager Migration
     */
    Docking(const json& controller, bool silent = true);

    /**
     * @brief Constructor with ConfigManager configuration (new, preferred)
     * Claude Generated: Phase 4 - Native ConfigManager support
     */
    Docking(const ConfigManager& config, bool silent = true);

    virtual ~Docking() = default;

    bool Initialise() override;

    void PerformDocking();

    void OptimiseBatch();
    void CollectStructures();
    void FilterStructures();

    void PostOptimise();

    Molecule getMolecule() const { return m_supramol; }

    StringList Files() const { return m_files; }

    void start() override { PerformDocking(); }

private:
    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    bool LoadRestartInformation() override { return true; }

    StringList MethodName() const override { return { std::string("dock") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override {}

    Molecule m_host_structure, m_guest_structure, m_supramol;
    std::vector<Position> m_initial_anchor = { Position{ 0, 0, 0 } };
    int m_step_X = 1, m_step_Y = 1, m_step_Z = 1;
    std::map<double, Vector> m_docking_list;
    std::vector<std::pair<Position, Position>> m_initial_list;
    std::vector<Position> m_anchor_accepted, m_rotation_accepted, m_reuse_anchor;
    std::vector<double> m_fragments_mass;
    bool m_check = false;
    bool m_PostFilter = true, m_PostOptimise = true, m_AutoPos = true, m_NoOpt = false;
    double m_sum_distance = 0;
    double m_scaling = 1.5;
    double m_window_seperator = 0.66666;
    double m_centroid_max_distance = 1e5;
    double m_centroid_tol_distance = 1e-1;
    double m_centroid_rot_distance = 1e-1;
    double m_energy_threshold = 200;
    int m_threads = 1;
    int m_docking_threads = 1;
    int m_charge = 0;
    int m_current_cycle = 0;
    int m_cycles = 1;
    int m_RMSDthreads = 1;
    int m_RMSDElement = 7;
    StringList m_files;
    std::string m_host, m_guest, m_complex, m_RMSDmethod;
    CurcumaOpt *m_optimise, *m_singlepoint;
    std::map<double, Molecule*> m_docking_result, m_optimisation_result, m_result_list, m_final_results, m_temp_results;
};
