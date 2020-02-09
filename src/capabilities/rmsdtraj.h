/*
 * <Trajectory RMSD Analyse. >
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

class RMSDTraj {
public:
    RMSDTraj();
    void setFile(const std::string& filename) { m_filename = filename; }
    void AnalyseTrajectory();

    void setReferenceStructure(const std::string& reference) { m_reference = reference; }

    /*! \brief Set the index of the fragment that is used for rmsd calculation/atom reordering */
    inline void setFragment(int fragment) { m_fragment = fragment; }

private:
    std::string m_filename, m_reference;
    std::ofstream m_rmsd_file;
    std::vector<Molecule> m_stored_structures;
    int m_fragment = -1;
};
