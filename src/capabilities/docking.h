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
#include "src/capabilities/optimiser/LBFGSInterface.h"
#include "src/capabilities/optimiser/LevMarDocking.h"

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include "src/tools/general.h"
#include "src/tools/geometry.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

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

    inline int execute()
    {
        std::pair<Position, Position> pair = OptimiseAnchor(&m_host, m_guest, m_position, m_rotation);
        m_last_position = pair.first;
        m_last_rotation = pair.second;
        return 0;
    }

    inline Position InitialPositon() const { return m_position; }
    inline Position InitialRotation() const { return m_rotation; }
    inline Position LastPosition() const { return m_last_position; }
    inline Position LastRotation() const { return m_last_rotation; }

private:
    Position m_position, m_rotation, m_last_position, m_last_rotation;
    Molecule m_host, m_guest;
};

class Thread : public CxxThread {
public:
    Thread() = default;
    ~Thread() = default;

    inline void setMolecule(const Molecule& molecule) { m_molecule = molecule; }
    inline Molecule getMolecule() const { return m_final; }
    inline int execute()
    {
        // OptimiseGeometryThreaded(&m_molecule, &m_result, &m_final, m_controller);
        m_final = CurcumaOpt::LBFGSOptimise(&m_molecule, m_controller);
        return 0;
    }

    inline void setController(const json& controller) { m_controller = controller; }
    std::string Output() const { return m_result; }
private:
    std::string m_result;
    Molecule m_molecule, m_final;
    json m_controller = OptJson;
};

static const json DockingJson = {
    { "Pos_X", 0.0 },
    { "Pos_Y", 0.0 },
    { "Pos_Z", 0.0 },
    { "AutoPos", true },
    { "Filter", true },
    { "PostOpt", true },
    { "Step_X", 10 },
    { "Step_Y", 10 },
    { "Step_z", 10 },
    { "Host", "none" },
    { "Guest", "none" },
    { "Complex", "none" },
    { "scaling", 1.5 },
    { "NoOpt", false },
    { "CentroidMaxDistance", 1e5 },
    { "CentroidTolDis", 1e-1 },
    { "RotationTolDis", 1e-1 },
    { "Threads", 1 },
    { "DockingThreads", 1 },
    { "Charge", 0 }
};

class Docking : public CurcumaMethod {

public:
    Docking(const json& controller = DockingJson, bool silent = true);
    virtual ~Docking() = default;

    bool Initialise() override;

    void PerformDocking();

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

    std::string MethodName() const override { return std::string("dock"); }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override {}

    Molecule m_host_structure, m_guest_structure, m_supramol;
    Position m_initial_anchor = Position{ 0, 0, 0 };
    int m_step_X = 1, m_step_Y = 1, m_step_Z = 1;
    std::map<double, Vector> m_docking_list;
    std::vector<Position> m_initial_list;
    std::map<double, Molecule*> m_result_list;
    std::vector<Position> m_anchor_accepted, m_rotation_accepted;
    std::vector<double> m_fragments_mass;
    bool m_check = false;
    bool m_PostFilter = true, m_PostOptimise = true, m_AutoPos = true, m_NoOpt = false;
    double m_sum_distance = 0;
    double m_scaling = 1.5;
    double m_window_seperator = 0.66666;
    double m_centroid_max_distance = 1e5;
    double m_centroid_tol_distance = 1e-1;
    double m_centroid_rot_distance = 1e-1;
    int m_threads = 1;
    int m_docking_threads = 1;
    int m_charge = 0;
    StringList m_files;
    std::string m_host, m_guest, m_complex;
};
