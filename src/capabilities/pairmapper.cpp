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
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "src/core/elements.h"
#include "src/core/molecule.h"

#include "src/tools/geometry.h"

#include "pairmapper.h"

PairMapper::PairMapper()
{
}

void PairMapper::FindPairs()
{
    string outfile = m_filename;
    for (int i = 0; i < 4; ++i)
        outfile.pop_back();

    m_intermol_file.open(outfile + "_intermol.dat");
    m_user_file.open(outfile + "_user.dat");
    m_intramol_file.open(outfile + "_intramol.dat");
    m_centroid_file.open(outfile + "_centroid.dat");
    m_pair_file.open(outfile + "_pairs.dat");

    std::ifstream input(m_filename);
    std::vector<std::string> lines;
    int atoms = 0;
    int index = 0;
    int i = 0;
    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;
    Molecule mol(atoms, 0);
    for (std::string line; getline(input, line);) {
        if (index == 0 && xyzfile) {
            atoms = stoi(line);
            mol = Molecule(atoms, 0);
        }
        if (xyzfile) {
            if (i > 1) {
                mol.setXYZ(line, i - 2);
            }
            if (i - 1 == atoms) {
                InitialisePairs(&mol);
                i = -1;
                mol = Molecule(atoms, 0);
            }
            ++i;
        } else {
            mol.setAtom(line, i);
        }
        index++;
    }

    m_intermol_file << "# ";
    for (const std::pair<int, int> p : m_inter_pairs) {
        m_intermol_file << "(" << std::setprecision(6) << p.first + 1 << "-" << p.second + 1 << ")   ";
    }
    m_intermol_file << std::endl;

    for (const std::pair<int, int> p : m_inter_pairs) {
        m_pair_file << "" << std::setprecision(6) << p.first + 1 << " " << p.second + 1 << std::endl;
    }

    m_intramol_file << "# ";
    for (const std::pair<int, int> p : m_intra_pairs) {
        m_intramol_file << "(" << std::setprecision(6) << p.first + 1 << "-" << p.second + 1 << ")   ";
    }
    m_intramol_file << std::endl;

    m_user_file << "# ";
    for (const std::pair<int, int> p : m_user_pairs) {
        m_user_file << "(" << std::setprecision(6) << p.first + 1 << "-" << p.second + 1 << ")   ";
    }
    m_user_file << std::endl;

    std::ifstream input2(m_filename);

    for (std::string line; getline(input2, line);) {
        if (index == 0 && xyzfile) {
            atoms = stoi(line);
            mol = Molecule(atoms, 0);
        }
        if (xyzfile) {
            if (i > 1) {
                mol.setXYZ(line, i - 2);
            }
            if (i - 1 == atoms) {
                ScanPairs(&mol);
                i = -1;
                mol = Molecule(atoms, 0);
            }
            ++i;
        } else {
            mol.setAtom(line, i);
        }
        index++;
    }

    m_intermol_file.close();
    m_intramol_file.close();
    m_centroid_file.close();
}

void PairMapper::InitialisePairs(const Molecule* molecule)
{

    for (std::size_t i = 0; i < molecule->AtomCount(); ++i) {
        const int element_A = molecule->Atom(i).first;
        for (std::size_t j = i + 1; j < molecule->AtomCount(); ++j) {
            const int element_B = molecule->Atom(j).first;
            for (const std::pair<int, int>& pair : m_element_pairs) {
                if ((element_A == pair.first && element_B == pair.second) || (element_B == pair.first && element_A == pair.second))
                    addPair(std::pair<int, int>(i, j));
            }
        }
    }

    if (m_proton_blacklist.size() == 0) {
        BlackListProtons(molecule);
    }

    std::vector<std::vector<int>> f = molecule->GetFragments();

    for (std::size_t i = 0; i < f.size(); ++i) {
        for (std::size_t k = 0; k < f[i].size(); ++k) {
            for (std::size_t l = k + 1; l < f[i].size(); ++l) {

                int a = f[i][k];
                int b = f[i][l];
                if ((std::find(m_proton_blacklist.begin(), m_proton_blacklist.end(), a) != m_proton_blacklist.end()) || (std::find(m_proton_blacklist.begin(), m_proton_blacklist.end(), b) != m_proton_blacklist.end()))
                    continue;
                double distance = molecule->Distance(a, b);
                if (distance < m_cutoff && distance > (Elements::CovalentRadius[molecule->Atom(a).first] + Elements::CovalentRadius[molecule->Atom(b).first]) * m_scaling) {
                    if ((molecule->Atom(a).first == 1 && (molecule->Atom(b).first == 7 || molecule->Atom(b).first == 8)) || ((molecule->Atom(a).first == 7 || molecule->Atom(a).first == 8) && molecule->Atom(b).first == 1)) {
                        addPair(std::pair<double, double>(a, b), m_intra_pairs);
                    }
                }
            }
        }
    }

    for (std::size_t i = 0; i < molecule->GetFragments().size(); ++i) {
        for (std::size_t j = i + 1; j < molecule->GetFragments().size(); ++j) {

            for (int a : f[i]) {

                for (int b : f[j]) {
                    if ((std::find(m_proton_blacklist.begin(), m_proton_blacklist.end(), a) != m_proton_blacklist.end()) || (std::find(m_proton_blacklist.begin(), m_proton_blacklist.end(), b) != m_proton_blacklist.end()))
                        continue;
                    double distance = molecule->Distance(a, b);
                    if (distance < m_cutoff) {
                        if ((molecule->Atom(a).first == 1 && (molecule->Atom(b).first == 7 || molecule->Atom(b).first == 8)) || ((molecule->Atom(b).first == 7 || molecule->Atom(b).first == 8) && molecule->Atom(b).first == 1)) {
                            addPair(std::pair<double, double>(a, b), m_inter_pairs);
                        }
                    }
                }
            }
        }
    }
}

void PairMapper::ScanPairs(const Molecule* molecule)
{
    for (const std::pair<int, int> p : m_intra_pairs) {
        m_intramol_file << molecule->Distance(p.first, p.second) << "    ";
    }
    m_intramol_file << std::endl;

    for (const std::pair<int, int> p : m_inter_pairs) {
        m_intermol_file << molecule->Distance(p.first, p.second) << "    ";
    }
    m_intermol_file << std::endl;

    for (const std::pair<int, int> p : m_user_pairs) {
        m_user_file << molecule->Distance(p.first, p.second) << "    ";
    }
    m_user_file << std::endl;

    std::vector<std::vector<int>> f = molecule->GetFragments();
    for (std::size_t i = 0; i < molecule->GetFragments().size(); ++i) {
        for (std::size_t j = i + 1; j < molecule->GetFragments().size(); ++j) {
            Geometry geom1 = molecule->getGeometryByFragment(i);
            Geometry geom2 = molecule->getGeometryByFragment(j);
            double distance = GeometryTools::Distance(GeometryTools::Centroid(geom1), GeometryTools::Centroid(geom2));
            m_centroid_file << distance << "    ";
        }
    }
    m_centroid_file << std::endl;
}

void PairMapper::addPair(std::pair<int, int> pair, std::vector<std::pair<int, int>>& pairs)
{
    bool exist = false;
    for (const std::pair<int, int> p : pairs)
        if (pair == p) {
            exist = true;
            break;
        }

    if (!exist)
        pairs.push_back(pair);
}

void PairMapper::BlackListProtons(const Molecule* molecule)
{
    for (std::size_t i = 0; i < molecule->AtomCount(); ++i) {
        if (molecule->Atom(i).first != 1)
            continue;
        double distance = 1e8;
        int element = 0;
        for (std::size_t j = 0; j < molecule->AtomCount(); ++j) {
            if (i == j)
                continue;
            double d = molecule->Distance(i, j);
            if (d < distance)
                element = molecule->Atom(j).first;
            distance = min(d, distance);
        }
        if (element == 6)
            m_proton_blacklist.push_back(i);
    }
}
