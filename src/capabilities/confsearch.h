/*
 * <Conformational Search based on Molecular Dynamics>
 * Copyright (C) 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <string>
#include <vector>

#include "src/tools/general.h"

#include "src/core/molecule.h"

#include "json.hpp"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "src/capabilities/curcumamethod.h"
using namespace curcuma;
static const nlohmann::json ConfSearchJson{
    { "method", "uff" },
    { "startT", 600 },
    { "endT", 300 },
    { "deltaT", 50 },
    { "repeat", 10 },
    { "time", 5e4 }, // 10 ps
    { "rmsd", 1.25 },
    { "threads", 1 },
    { "energy_window", 100 },
    { "writeXYZ", true },
    { "printOutput", true },
    { "MaxTime", 5000 },
    { "T", 298.15 },
    { "dt", 1 }, // single step in fs
    { "rm_COM", 100 }, // remove translation and rotation every x fs
    { "charge", 0 },
    { "Spin", 0 },
    { "rmrottrans", 0 },
    { "nocenter", false },
    { "dump", 50 },
    { "print", 1000 },
    { "unique", false },
    { "rmsd", 1.5 },
    { "opt", false },
    { "hmass", 1 },
    { "velo", 1 },
    { "rescue", false },
    { "coupling", 10 },
    { "MaxTopoDiff", 15 },
    { "impuls", 0 },
    { "impuls_scaling", 0.75 },
    { "writeinit", false },
    { "initfile", "none" },
    { "norestart", false },
    { "writerestart", 1000 },
    { "rattle", false },
    { "rattle_tolerance", 1e-2 },
    { "rattle_maxiter", 10 },
    { "thermostat", "csvr" },
    { "respa", 1 },
    { "dipole", false },
    { "seed", -1 },
    { "cleanenergy", false },
    { "wall", "none" }, // can be spheric or rect
    { "wall_type", "logfermi" }, // can be logfermi or harmonic
    { "wall_spheric_radius", 0 },
    { "wall_xl", 0 },
    { "wall_yl", 0 },
    { "wall_zl", 0 },
    { "wall_x_min", 0 },
    { "wall_x_max", 0 },
    { "wall_y_min", 0 },
    { "wall_y_max", 0 },
    { "wall_z_min", 0 },
    { "wall_z_max", 0 },
    { "wall_temp", 298.15 },
    { "wall_beta", 6 }
};

class curcuma::Molecule;

class ConfSearch : public CurcumaMethod {
public:
    ConfSearch(const json& controller, bool silent);
    ~ConfSearch();

    void setFile(const std::string& file);
    virtual bool Initialise() override;

    virtual void start() override;

private:
    std::string PerformMolecularDynamics(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter);

    std::string PerformOptimisation(const std::string& filename, const nlohmann::json& parameter);

    std::string PerformFilter(const std::string& filename, const nlohmann::json& parameter);

    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation() override;

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation() override;

    virtual StringList MethodName() const override
    {
        return { "ConfSearch" };
    }

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile() override;

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson() override;

    StringList m_error_list;
    std::string m_filename, m_method, m_thermostat;
    bool m_silent = true, m_rattle = true;
    double m_dT = 4;
    std::vector<Molecule*> m_in_stack, m_final_stack;
    int m_spin = 0, m_charge = 0, m_repeat = 5, m_threads = 1;
    double m_time = 1e4, m_startT = 500, m_endT = 300, m_deltaT = 50, m_currentT = 0, m_rmsd = 1.25, m_energy_window = 100;
    Matrix m_topo_matrix;
};
