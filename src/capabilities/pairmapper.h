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

#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>


#include "src/core/molecule.h"

class PairMapper {
public:
    PairMapper();
    void setFile(const std::string& filename) { m_filename = filename; }
    void FindPairs();
    inline void addPair(std::pair<int, int> pair) { addPair(pair, m_user_pairs); }
    inline void addElementPair(std::pair<int, int> pair) { addPair(pair, m_element_pairs); }

private:
    void InitialisePairs(const Molecule* molecule);
    void ScanPairs(const Molecule* molecule);
    void BlackListProtons(const Molecule* molecule);
    void addPair(std::pair<int, int> pair, std::vector<std::pair<int, int>>& pairs);
    std::string m_filename;
    std::vector<std::pair<int, int>> m_intra_pairs, m_inter_pairs, m_user_pairs, m_element_pairs;
    std::vector<int> m_proton_blacklist;
    bool m_intramolecular = false, m_intermolecule = true;
    double m_cutoff = 2.5, m_scaling = 1.3;
    std::ofstream m_intermol_file, m_intramol_file, m_centroid_file, m_user_file, m_pair_file;
};
