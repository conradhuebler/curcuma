/*
 * <Statistics of conformers >
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

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "curcumamethod.h"

static json ConfStatJson{
    "Cutoff", 10,
    "Temp", 298.5,
    "Threshold", 0.5
};

class ConfStat : public CurcumaMethod {
public:
    ConfStat(const json& controller = ConfStatJson, bool silent = true);
    void setFileName(const std::string& filename)
    {
        m_filename = filename;
    }

    void start() override; // TODO make pure virtual and move all main action here

private:
    /* Lets have this for all modules */
    inline nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    inline bool LoadRestartInformation() override { return true; }

    inline std::string MethodName() const override { return std::string("ConfStat"); }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override{};

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    std::string m_filename;
    std::vector<double> m_energies;
    double m_temp = 298.5;
    double m_cutoff = 10;
    double m_print_threshold = 0.5;
};
