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

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "src/capabilities/rmsd.h"

#include "src/tools/general.h"

#include "confscan.h"

ConfScan::ConfScan()
{
}

ConfScan::~ConfScan()
{
    for (auto i : m_molecules) {
        delete i.second;
    }
}

bool ConfScan::openFile()
{

    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos;

    if (xyzfile == false)
        throw 1;

    std::vector<std::string> lines;
    std::ifstream input(m_filename);

    int atoms = 0;
    int index = 0;
    int i = 0;
    int molecule = 0;
    Molecule* mol = new Molecule(atoms, 0);
    for (std::string line; getline(input, line);) {
        if (i == 0 && xyzfile) {
            atoms = stoi(line);
            if (atoms < 1)
                throw 2;
            std::pair<std::string, Molecule*> pair(mol->Name(), mol);
            m_molecules.push_back(pair);
            m_ordered_list.insert(std::pair<double, int>(mol->Energy(), molecule));
            molecule++;
            mol = new Molecule(atoms, 0);
        }
        if (i == 1) {
            StringList list = Tools::SplitString(line);
            mol->setName(list[0]);
            mol->setEnergy(std::stod((list[3])));
            // std::cout << mol->Name() << " " << mol->Energy() << std::endl;
        }
        if (i > 1) {
            mol->setXYZ(line, i - 2);
        }
        index++;
        ++i;
        if (i - 2 == atoms)
            i = 0;
    }

    return true;
}

void ConfScan::scan()
{
    openFile();
    Molecule* mol1;
    bool ok = true;

    std::vector<std::string> filtered;
    std::size_t start = 0;
    std::size_t ende = m_ordered_list.size();
    for (const auto& i : m_ordered_list) {

        int index = i.second;
        mol1 = m_molecules.at(index).second;
        if (mol1->AtomCount() == 0)
            continue;
        mol1->CalculateRotationalConstants();
        std::cout << start / double(ende) * 100;
        for (const auto* mol2 : m_result) {
            double difference = abs(mol1->Energy() - mol2->Energy()) * 2625.5;
            double rmsd = 0;
            double Ia = abs(mol1->Ia() - mol2->Ia()) / mol2->Ia();
            double Ib = abs(mol1->Ib() - mol2->Ib()) / mol2->Ib();
            double Ic = abs(mol1->Ic() - mol2->Ic()) / mol2->Ic();
            double diff_rot = Ia + Ib + Ic;
            RMSDDriver* driver = new RMSDDriver(mol1, mol2);
            // driver->setForceReorder(reorder);
            //driver->setFragment(fragment);
            //driver->setCheckConnections(check_connect);
            driver->AutoPilot();
            rmsd = driver->RMSD();
            // std::cout << mol2->Ia() << " " << mol2->Ib() << " " << mol2->Ic() << std::endl;

            //std::cout << Ia/mol2->Ia() << " " << Ib/mol2->Ib() << " " << Ic/mol2->Ic() << std::endl;
            delete driver;
            if ((difference < m_energy_threshold && rmsd < m_rmsd_threshold && diff_rot < 0.3) || m_result.size() >= m_maxrank) {
                ok = false;
                filtered.push_back(mol1->Name());
            }
            std::cout << ".";
        }
        if (ok)
            m_result.push_back(mol1);
        ok = true;
        start++;
        std::cout << " " << start / double(ende) * 100 << "% done!" << std::endl;
    }

    for (const auto molecule : m_result)
        molecule->print_geom(false);

    if (m_writeXYZ) {
        for (const auto molecule : m_result)
            molecule->writeXYZFile();
    }

    std::sort(filtered.begin(), filtered.end());
    std::vector<std::string>::iterator iterator;
    std::cerr << m_result.size() << " structures were kept - of " << m_molecules.size() << " total!" << std::endl;
    iterator = std::unique(filtered.begin(), filtered.end());
    filtered.resize(std::distance(filtered.begin(), iterator));
    std::cerr << "List of filtered names ... " << std::endl;
    for (const auto& element : filtered)
        std::cout << element << std::endl;
    std::cout << " done :-) " << std::endl;
}
