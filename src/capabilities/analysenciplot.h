/*
 * <Post-process NCI RDG vs sign(lambda2)rho plots>
 * Copyright (C) 2021 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "curcumamethod.h"

static json AnalyseNCIPlotJson{
    { "bins", 1000 },
    { "scale_d1", 1 },
    { "scale_d2", 1 },
    { "local_distance", false }
};

struct NCIplot {
    std::map<double, double> m_plot;
    double m_RDG_max = 0, m_RDG_min = 0;
};

class AnalyseNCIPlot : public CurcumaMethod {
public:
    AnalyseNCIPlot(const json& controller);

    void start() override; // TODO make pure virtual and move all main action here

    inline void setFiles(const std::string& f1, const std::string& f2)
    {
        m_file1 = f1;
        m_file2 = f2;
    }

private:
    double getDistance(const std::map<double, double>& map, const std::map<double, double>::iterator& begin, const std::map<double, double>::iterator& end);

    /* Lets have this for all modules */
    inline nlohmann::json WriteRestartInformation() override { return json(); }

    /* Lets have this for all modules */
    inline bool LoadRestartInformation() override { return true; }

    inline StringList MethodName() const override { return { std::string("nci") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override{};

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    std::map<double, double> LoadFile(const std::string& file);

    std::string m_file1, m_file2;
    std::map<double, double> m_NCI1, m_NCI2;

    double m_bins = 1000, m_scale_d1 = 1, m_scale_d2 = 1;
    bool m_local_distance = false;
};
