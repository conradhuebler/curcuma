/*
 * <Scan and judge conformers from different input. >
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

#include "src/core/molecule.h"

class ConfScan {
public:
    ConfScan();
    ~ConfScan();

    void setFileName(const std::string& filename) { m_filename = filename; }
    void scan();

    inline void setMaxRank(int rank) { m_maxrank = rank; }
    inline void setWriteXYZ(bool writeXYZ) { m_writeXYZ = writeXYZ; }

private:
    bool openFile();

    std::string m_filename;
    std::map<double, int> m_ordered_list;
    std::vector<std::pair<std::string, Molecule*>> m_molecules;
    double m_energy_threshold = 1.0, m_rmsd_threshold = 1.0;
    std::vector<Molecule*> m_result;
    int m_maxrank = 10000;
    bool m_writeXYZ = false;
};
