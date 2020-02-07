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

#include "src/core/molecule.h"

#include "pairmapper.h"

PairMapper::PairMapper()
{
}

void PairMapper::FindPairs()
{
    std::ifstream input(m_filename);
    std::vector<std::string> lines;
    int atoms = 0;
    int index = 0;
    int i = 0;
    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos;
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

    for (const std::pair<int, int> p : m_pairs) {
        std::cout << std::setprecision(6) << p.first << " - " << p.second << "   ";
    }
    std::cout << std::endl;

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
}

void PairMapper::InitialisePairs(const Molecule* molecule)
{
    if (m_intramolecular) // can be better, I know ....
    {
        for (std::size_t i = 0; i < molecule->GetFragments().size(); ++i) {
            for (int a : molecule->GetFragments()[i]) {
                for (int b : molecule->GetFragments()[i]) {
                    double distance = molecule->Distance(a, b);
                    if (distance < m_cutoff) {
                        if ((molecule->Atom(a).first == 1 && molecule->Atom(b).first != 1) || (molecule->Atom(a).first != 1 && molecule->Atom(b).first == 1)) {
                            addPair(std::pair<double, double>(a, b));
                        }
                        //std::cout <<  std::setprecision(6) << distance << " " << a << "(" << Atom(a).first << ") - " << b << "(" << Atom(b).first << ")" << std::endl;
                    }
                }
            }
        }
    }

    if (m_intermolecule) {
        std::vector<std::vector<int>> f = molecule->GetFragments();
        for (std::size_t i = 0; i < molecule->GetFragments().size(); ++i) {
            for (std::size_t j = i + 1; j < molecule->GetFragments().size(); ++j) {

                for (int a : f[i]) {

                    for (int b : f[j]) {
                        double distance = molecule->Distance(a, b);
                        if (distance < m_cutoff) {
                            if ((molecule->Atom(a).first == 1 && molecule->Atom(b).first != 1) || (molecule->Atom(a).first != 1 && molecule->Atom(b).first == 1)) {
                                addPair(std::pair<double, double>(a, b));
                            }
                            //std::cout <<  std::setprecision(6) << distance << " " << a << "(" << Atom(a).first << ") - " << b << "(" << Atom(b).first << ")" << std::endl;
                        }
                    }
                }
            }
        }
    }
}

void PairMapper::ScanPairs(const Molecule* molecule)
{
    for (const std::pair<int, int> p : m_pairs) {
        std::cout << molecule->Distance(p.first, p.second) << "    ";
    }
    std::cout << std::endl;
}

void PairMapper::addPair(std::pair<int, int> pair)
{
    bool exist = false;
    for (const std::pair<int, int> p : m_pairs)
        if (pair == p) {
            exist = true;
            break;
        }

    if (!exist)
        m_pairs.push_back(pair);
}
