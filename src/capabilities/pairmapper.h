/*
 * <Collect HBond Pairs in Trajectories. >
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

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#pragma once

#include "src/core/molecule.h"

class PairMapper {
public:
    PairMapper();
    void setFile(const std::string& filename) { m_filename = filename; }
    void FindPairs();

private:
    void InitialisePairs(const Molecule* molecule);
    void ScanPairs(const Molecule* molecule);
    void addPair(std::pair<int, int> pair);
    std::string m_filename;
    std::vector<std::pair<int, int>> m_pairs;
    bool m_intramolecular = false, m_intermolecule = true;
    double m_cutoff = 2.5;
};
