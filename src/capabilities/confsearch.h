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

static const nlohmann::json ConfSearchJson{
    { "gfn", 66 },
    { "charge", 0 },
    { "Spin", 0 },
    { "startT", 500 },
    { "endT", 300 },
    { "deltaT", 50 },
    { "repeat", 5 },
    { "time", 1e4 } // 10 ps
};

class Molecule;

class ConfSearch : public CurcumaMethod {
public:
    ConfSearch(const json& controller, bool silent);
    ~ConfSearch();

    void setFile(const std::string& file);
    virtual bool Initialise();

    virtual void start();

private:
    std::vector<Molecule*> PerformMolecularDynamics(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter);

    std::vector<Molecule*> PerformOptimisation(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter);

    void PerformFilter(const std::vector<Molecule*>& molecules, const nlohmann::json& parameter);

    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation();

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation();

    virtual StringList MethodName() const
    {
        return { "ConfSearch" };
    }

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile();

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson();

    StringList m_error_list;
    std::string m_filename, m_basename;
    bool m_silent = true;

    std::vector<Molecule*> m_in_stack, m_final_stack;
    int m_gfn = 66, m_spin = 0, m_charge = 0, m_repeat = 5;
    double m_time = 1e4, m_startT = 500, m_endT = 300, m_deltaT = 50;
};
