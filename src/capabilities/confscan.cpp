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
#include <iomanip>
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

    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;

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
            if (list.size() == 7) {
                try {
                    mol->setEnergy(std::stod((list[4])));
                } catch (const std::string& what_arg) {
                    mol->setEnergy(0);
                }
            } else if (list.size() == 4) {
                if (list[0].compare("SCF") == 0 && list[1].compare("done") == 0) {
                    mol->setName("Molecule " + std::to_string(molecule));
                    try {
                        mol->setEnergy(std::stod((list[2])));
                    } catch (const std::string& what_arg) {
                        mol->setEnergy(0);
                    }
                } else {
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
                }
            } else if (list.size() == 1) {
                try {
                    mol->setEnergy(std::stod((list[0])));
                } catch (const std::string& what_arg) {
                }
                mol->setName("Molecule " + std::to_string(molecule));
            } else {
                for (const string& s : list) {
                    double energy = 0;
                    if (Tools::isDouble(s)) {
                        energy = std::stod(s);
                        mol->setEnergy(energy);
                        break;
                    }
                }
                mol->setName(NamePattern(molecule));
            }
            if (m_noname)
                mol->setName(NamePattern(molecule));
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

    std::map<std::string, std::vector<std::string>> filtered;
    std::size_t fail = 0;
    std::size_t start = 0;
    std::size_t ende = m_ordered_list.size();

    std::string result_name = m_filename;
    result_name.erase(result_name.end() - 4, result_name.end());
    std::string nearly_missed = result_name + "_missed.xyz";
    std::string failed = result_name + "_failed.xyz";

    result_name += +"_filter.xyz";

    std::ofstream result_file;
    result_file.open(result_name);
    result_file.close();

    std::ofstream missed_file;
    missed_file.open(nearly_missed);
    missed_file.close();

    std::ofstream failed_file;
    missed_file.open(failed);
    missed_file.close();

    std::vector<std::vector<int>> rules_list;
    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << "" << std::endl;

    if (m_heavy)
        std::cout << "    RMSD Calculation will be performed only on heavy atoms! " << std::endl;
    else
        std::cout << "    RMSD Calculation will be performed on all atoms! " << std::endl;

    std::cout << "    RMSD Threshold set to: " << m_rmsd_threshold << " Angstrom" << std::endl;
    std::cout << "    Energy Threshold set to: " << m_energy_threshold << " kJ/mol" << std::endl;
    std::cout << "    Average Difference in rot constants: " << std::endl;
    std::cout << "    Loose Threshold: " << m_diff_rot_loose << std::endl;
    std::cout << "    Tight Threshold: " << m_diff_rot_tight << std::endl;
    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << std::endl;

    for (auto& i : m_ordered_list) {

        int index = i.second;
        mol1 = m_molecules.at(index).second;
        if (mol1->AtomCount() == 0) {
            fail++;
            continue;
        }
        if (m_result.size() >= m_maxrank) {
            ok = false;
            continue;
        }
        /* Lets postpone the reordering */
        std::vector<Molecule*> temp_list;
        if (mol1->Energy() == 0) {
            filtered[mol1->Name()].push_back("Empty");
            ok = false;
        } else {
            mol1->CalculateRotationalConstants();
            for (auto* mol2 : m_result) {
                if (filtered.count(mol1->Name())) {
                    std::cout << mol1->Name() << " already rejected. Skipping check against " << mol2->Name() << std::endl;
                    continue;
                }
                RMSDDriver* driver = new RMSDDriver;
                driver->setSilent(true);
                driver->setProtons(!m_heavy);
                driver->setForceReorder(ForceReorder());
                driver->setCheckConnections(CheckConnections());

                std::cout << std::endl
                          << std::setprecision(10)
                          << std::endl
                          << std::endl
                          << "Reference Molecule:" << mol1->Name() << " (" << mol1->Energy() << " Eh)        Target Molecule " << mol2->Name() << " (" << mol2->Energy() << " Eh)" << std::endl;
                double difference = abs(mol1->Energy() - mol2->Energy()) * 2625.5;
                if (m_energy_cutoff > 0)
                    if (difference < m_energy_cutoff) {
                        ok = false;
                        continue;
                    }
                double rmsd = 0;
                double Ia = abs(mol1->Ia() - mol2->Ia()) / mol2->Ia();
                double Ib = abs(mol1->Ib() - mol2->Ib()) / mol2->Ib();
                double Ic = abs(mol1->Ic() - mol2->Ic()) / mol2->Ic();

                double diff_rot = (Ia + Ib + Ic) * 0.33333;

                driver->setReference(mol1);
                driver->setTarget(mol2);

                driver->AutoPilot();
                rmsd = driver->RMSD();

                std::cout << "Energy Difference: " << std::setprecision(2) << difference << " kJ/mol" << std::endl;
                std::cout << "Average Difference in rot constant " << std::setprecision(4) << diff_rot << " some unit" << std::endl;
                std::cout << "RMSD = " << std::setprecision(5) << rmsd << " A" << std::endl;

                if (rmsd > m_rmsd_threshold && difference < 1 && diff_rot < 0.1 && diff_rot > 0.01) {
                    if (!m_prevent_reorder) {
                        temp_list.push_back(mol2);
                        std::cout << "~~ Reordering forced as energies and rotational constants are too close and rmsd is too different! ~~" << std::endl;
                        std::cout << "*** Adding " << mol2->Name() << " to the list as RMSD is " << rmsd << "! ***" << std::endl;
                    } else {
                        std::cout << "~~ Reordering prevented, so no more will be done regarding these structures. ~~" << std::endl;
                    }
                } else {
                    if ((difference < m_energy_threshold && rmsd < m_rmsd_threshold && diff_rot < m_diff_rot_loose)) {
                        ok = false;
                        filtered[mol1->Name()].push_back(mol2->Name());
                        std::cout << "  ** Rejecting structure **" << std::endl;
                        if (rmsd <= m_rmsd_threshold * m_nearly_missed) {
                            std::cout << " Nearly missed for " << mol1->Name() << std::endl;
                            m_nearly.push_back(mol1);
                        }
                        continue;
                    }
                    if (diff_rot < m_diff_rot_tight && difference < m_energy_threshold) {
                        ok = false;
                        filtered[mol1->Name()].push_back(mol2->Name());
                        std::cout << "  ** Rejecting structure **" << std::endl;
                        continue;
                    }
                }
                delete driver;
            }
        }

        if (temp_list.size()) {
            std::cout << std::endl
                      << "             ###   " << std::setprecision(4) << start / double(ende) * 100 << "% done!   ###" << std::endl;
            std::cout << "*** Start Reordering Block **" << std::endl;
        }
        if (!filtered.count(mol1->Name())) {
            for (const auto mol2 : temp_list) {
                if (filtered.count(mol1->Name())) {
                    std::cout << mol1->Name() << " already rejected. Skipping check against " << mol2->Name() << std::endl;
                    continue;
                }

                std::cout << "*** Reordering " << mol2->Name() << " now! ***" << std::endl;
                double difference = abs(mol1->Energy() - mol2->Energy()) * 2625.5;

                double Ia = abs(mol1->Ia() - mol2->Ia()) / mol2->Ia();
                double Ib = abs(mol1->Ib() - mol2->Ib()) / mol2->Ib();
                double Ic = abs(mol1->Ic() - mol2->Ic()) / mol2->Ic();

                double diff_rot = (Ia + Ib + Ic) * 0.33333;

                RMSDDriver* driver = new RMSDDriver;
                driver->setSilent(true);
                driver->setProtons(!m_heavy);
                driver->setForceReorder(ForceReorder());
                driver->setCheckConnections(CheckConnections());

                driver->setReference(mol1);
                driver->setTarget(mol2);

                driver->AutoPilot();
                double rmsd = driver->RMSD();

                std::cout << std::endl
                          << std::endl
                          << "*** Reordering forced as energies and rotational constants are too close and rmsd (" << rmsd << ") is too different! ***" << std::endl
                          << std::endl;

                bool digDeeper = true;
                std::cout << "*** Checking old reordering solutions first. ***" << std::endl;
                for (const auto rule : rules_list) {
                    double tmp_rmsd = driver->Rules2RMSD(rule);
                    std::cout << tmp_rmsd << " A ";
                    if (tmp_rmsd < m_rmsd_threshold) {
                        digDeeper = false;
                        std::cout << std::endl
                                  << std::endl
                                  << "*** Old reordering solution worked here! ***" << std::endl
                                  << std::endl;
                        ok = false;
                        filtered[mol1->Name()].push_back(mol2->Name());
                        break;
                    }
                }
                std::cout << std::endl;
                if (digDeeper) {
                    std::cout << "*** Nothing fitting found, reorder now! ***" << std::endl;
                    driver->setForceReorder(true);
                    driver->AutoPilot();
                    driver->setForceReorder(ForceReorder());

                    double rmsd_tmp = driver->RMSD();
                    rules_list.push_back(driver->ReorderRules());

                    std::cout << "New rmsd is " << rmsd_tmp << ". Old was " << rmsd << std::endl;
                    if (rmsd_tmp > rmsd) {
                        std::cout << "Reordering failed, a better fit was not found! This may happen." << std::endl;
                        m_failed.push_back(mol1);
                        m_failed.push_back(mol2);
                    }
                    rmsd = rmsd_tmp;
                }
                if ((difference < m_energy_threshold && rmsd < m_rmsd_threshold && diff_rot < m_diff_rot_loose)) {
                    ok = false;
                    filtered[mol1->Name()].push_back(mol2->Name());
                    std::cout << "  ** Rejecting structure **" << std::endl;
                    if (rmsd <= m_rmsd_threshold * m_nearly_missed) {
                        std::cout << " Nearly missed for " << mol1->Name() << std::endl;
                        m_nearly.push_back(mol1);
                    }
                    continue;
                }
                if (diff_rot < m_diff_rot_tight && difference < m_energy_threshold) {
                    ok = false;
                    filtered[mol1->Name()].push_back(mol2->Name());
                    std::cout << "  ** Rejecting structure **" << std::endl;
                    continue;
                }
                delete driver;
            }
        } else {
            std::cout << "Skipping check, as structure already rejected." << std::endl;
        }
        if (ok /* && m_result.size() < m_maxrank*/) {
            m_result.push_back(mol1);
            mol1->appendXYZFile(result_name);

            std::cout << std::endl
                      << "Having now " << m_result.size() << " structures accepted." << std::endl
                      << "               ** Accepting " << mol1->Name() << " **" << std::endl;
            /* Write each accepted structure immediately */
        }
        ok = true;
        start++;
        std::cout << std::endl
                  << std::endl
                  << "             ###   " << std::setprecision(4) << start / double(ende) * 100 << "% done!   ###" << std::endl;
    }

    if (m_result.size() == 0)
        return;

    if (m_writeXYZ) {
        for (const auto molecule : m_result)
            molecule->writeXYZFile();
    }

    for (const auto molecule : m_nearly) {
        molecule->appendXYZFile(nearly_missed);
    }

    for (const auto molecule : m_failed) {
        molecule->appendXYZFile(failed);
    }

    std::cout << m_result.size() << " structures were kept - of " << m_molecules.size() - fail << " total!" << std::endl;

    std::cout << "Best structure is " << m_result[0]->Name() << " Energy = " << std::setprecision(8) << m_result[0]->Energy() << std::endl;
    std::cout << "List of filtered names ... " << std::endl;
    for (const auto& element : filtered) {
        std::cout << element.first << " rejected due to: ";
        for (const std::string& str : element.second)
            std::cout << str << " ";
        std::cout << std::endl;
    }
    std::cout << " done :-) " << std::endl;
}
