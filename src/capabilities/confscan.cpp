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
            if (atoms < 1) {
                ++i;
                continue;
            }
            //if (mol->AtomCount())
            //{
            std::pair<std::string, Molecule*> pair(mol->Name(), mol);
            m_molecules.push_back(pair);
            m_ordered_list.insert(std::pair<double, int>(mol->Energy(), molecule));
            molecule++;
            //}

            mol = new Molecule(atoms, 0);
        }
        if (i == 1) {
            StringList list = Tools::SplitString(line);
            if (list.size() == 4) {
                mol->setName(list[0]);
                //mol->setName("Molecule " + std::to_string(molecule));
                if (list[3] == "")
                    mol->setEnergy(0);
                else {
                    try {
                        mol->setEnergy(std::stod((list[3])));
                    } catch (const std::string& what_arg) {
                        mol->setEnergy(0);
                    }
                }
            } else if (list.size() == 1) {
                try {
                    mol->setEnergy(std::stod((list[0])));
                } catch (const std::string& what_arg) {
                }
                mol->setName("Molecule " + std::to_string(molecule));
            }
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
    std::size_t fail = 0;
    std::size_t start = 0;
    std::size_t ende = m_ordered_list.size();
    for (auto& i : m_ordered_list) {

        int index = i.second;
        mol1 = m_molecules.at(index).second;
        if (mol1->AtomCount() == 0) {
            fail++;
            continue;
        }

        if (mol1->Energy() == 0) {
            filtered.push_back(mol1->Name());
            ok = false;
        } else {
            mol1->CalculateRotationalConstants();
            //std::cout << start / double(ende) * 100;
            for (auto* mol2 : m_result) {
                std::vector< std::string >::iterator find = std::find(filtered.begin(), filtered.end(), mol2->Name());
                if(find != filtered.end())
                {
                    std::cout << "We already had this structure on blacklist" << std::endl;
                    continue;
                }
                std::cout << std::endl
                          << std::endl
                          << std::endl
                          << "Reference Molecule:" << mol1->Name() << "        Target Molecule " << mol2->Name() << std::endl;
                double difference = abs(mol1->Energy() - mol2->Energy()) * 2625.5;

                double rmsd = 0;
                double Ia = abs(mol1->Ia() - mol2->Ia()) / mol2->Ia();
                double Ib = abs(mol1->Ib() - mol2->Ib()) / mol2->Ib();
                double Ic = abs(mol1->Ic() - mol2->Ic()) / mol2->Ic();

                double diff_rot = (Ia + Ib + Ic) * 0.33333;
                std::cout << "Energy Difference: " << difference << "        Average Difference in rot constant " << diff_rot << std::endl;

                RMSDDriver* driver = new RMSDDriver(mol1, mol2);
                driver->setSilent(true);
                driver->setProtons(!m_heavy);
                driver->setForceReorder(ForceReorder());
                driver->setCheckConnections(CheckConnections());
                //driver->setFragment(fragment);
                driver->AutoPilot();
                rmsd = driver->RMSD();
                if (rmsd > 0.5 && difference < 1 && diff_rot < 0.1 && diff_rot > 0.01) {
                    std::cout << std::endl
                              << std::endl
                              << "*** Reordering forced as energies and rotational constants are too close and rmsd (" << rmsd << ") is too different! ***" << std::endl
                              << std::endl;
                    driver->setForceReorder(true);
                    driver->AutoPilot();
                    double rmsd_tmp = driver->RMSD();
                    std::cout << "New rmsd is " << rmsd_tmp << ". Old was " << rmsd << std::endl;
                    rmsd = rmsd_tmp;
                } else {
                    std::cout << "RMSD is " << rmsd << std::endl;
                }
                // std::cout << mol2->Ia() << " " << mol2->Ib() << " " << mol2->Ic() << std::endl;

                //std::cout << Ia/mol2->Ia() << " " << Ib/mol2->Ib() << " " << Ic/mol2->Ic() << std::endl;
                delete driver;

                //if(mol1->Energy() > mol2->Energy() && rmsd < m_rmsd_threshold)
                //    std::swap(mol1, mol2);

                if ((difference < m_energy_threshold && rmsd < m_rmsd_threshold && diff_rot < 0.3) || m_result.size() >= m_maxrank) {
                    ok = false;
                    filtered.push_back(mol1->Name());
                    std::cout << "  ** Rejecting structure **" << std::endl;
                }
                if (diff_rot < 0.01) {
                    ok = false;
                    filtered.push_back(mol1->Name());
                    std::cout << "  ** Rejecting structure **" << std::endl;
                }
                //std::cout << ".";
            }
        }
        if (ok)
            m_result.push_back(mol1);
        ok = true;
        start++;
        std::cout << std::endl
                  << std::endl
                  << "             ###   " << start / double(ende) * 100 << "% done!   ###" << std::endl;
    }
    if (m_result.size() == 0)
        return;

    for (const auto molecule : m_result)
        molecule->print_geom(false);

    if (m_writeXYZ) {
        for (const auto molecule : m_result)
            molecule->writeXYZFile();
    }

    std::sort(filtered.begin(), filtered.end());
    std::vector<std::string>::iterator iterator;
    std::cout << m_result.size() << " structures were kept - of " << m_molecules.size() - fail << " total!" << std::endl;

    std::cout << "Best structure is " << m_result[0]->Name() << " Energy = " << m_result[0]->Energy() << std::endl;
    iterator = std::unique(filtered.begin(), filtered.end());
    filtered.resize(std::distance(filtered.begin(), iterator));
    std::cout << "List of filtered names ... " << std::endl;
    for (const auto& element : filtered)
        std::cout << element << std::endl;
    std::cout << " done :-) " << std::endl;
}
