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

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"

#include "src/tools/general.h"

#include "json.hpp"

// for convenience
using json = nlohmann::json;

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

    int molecule = 0;

    FileIterator file(m_filename);
    while (!file.AtEnd()) {
        Molecule* mol = new Molecule(file.Next());
        m_ordered_list.insert(std::pair<double, int>(mol->Energy(), molecule));
        molecule++;
        if (m_noname)
            mol->setName(NamePattern(molecule));
        std::pair<std::string, Molecule*> pair(mol->Name(), mol);
        m_molecules.push_back(pair);
    }

    return true;
}

void ConfScan::setMolecules(const std::map<double, Molecule*>& molecules)
{
    std::cout << "Loading molecules" << std::endl;
    int molecule = 0;
    for (const auto m : molecules) {
        Molecule* mol = new Molecule(m.second);
        if (m_noname)
            mol->setName(NamePattern(molecule));
        std::pair<std::string, Molecule*> pair(mol->Name(), mol);
        m_molecules.push_back(pair);
        m_ordered_list.insert(std::pair<double, int>(mol->Energy(), molecule));
        molecule++;
    }
    m_writeFiles = false;
    std::cout << "Loading molecules done!" << std::endl;
    m_filename = "blocklist.xyz";
}

bool ConfScan::LoadRestartInformation()
{
    if (!Restart())
        return false;
    StringList files = RestartFiles();
    if (!m_prevent_reorder)
        m_prevent_reorder = files.size() > 1;
    int error = 0;
    for (const auto& f : files) {
        std::vector<std::vector<int>> reorder_cached;

        std::cout << "Reading file " << f << std::endl;
        std::ifstream file(f);
        json restart;
        try {
            file >> restart;
        } catch (json::type_error& e) {
            error++;
            continue;

        } catch (json::parse_error& e) {
            error++;
            continue;
        }

        json confscan;
        try {
            confscan = restart[MethodName()];
        } catch (json::type_error& e) {
            error++;
            continue;
        }

        try {
            reorder_cached = Tools::String2VectorVector(confscan["ReorderRules"]);
        } catch (json::type_error& e) {
        }
        try {
            m_reference_last_energy = confscan["ReferenceLastEnergy"];
        } catch (json::type_error& e) {
        }
        try {
            m_target_last_energy = confscan["TargetLastEnergy"];
        } catch (json::type_error& e) {
        }
        for (const auto& vector : reorder_cached)
            if (std::find(m_reorder_rules.begin(), m_reorder_rules.end(), vector) == m_reorder_rules.end())
                m_reorder_rules.push_back(vector);
    }
    if (files.size() != 1) {
        m_reference_last_energy = 0;
        m_target_last_energy = 0;
    }

    m_useRestart = files.size() == 1 && error != int(files.size());

    std::cout << "Starting with " << m_reorder_rules.size() << " initial reorder rules." << std::endl;
    return true;
}

nlohmann::json ConfScan::WriteRestartInformation()
{
    json block;
    block["ReorderRules"] = Tools::VectorVector2String(m_reorder_rules);
    block["ReferenceLastEnergy"] = m_reference_last_energy;
    block["TargetLastEnergy"] = m_target_last_energy;
    return block;
}

void ConfScan::scan()
{
    LoadRestartInformation();

    m_silent = m_ordered_list.size() > 500;
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
    if (m_writeFiles) {
        result_file.open(result_name);
        result_file.close();
    }

    std::ofstream missed_file;
    if (m_writeFiles) {
        missed_file.open(nearly_missed);
        missed_file.close();
    }

    std::ofstream failed_file;
    if (m_writeFiles) {
        missed_file.open(failed);
        missed_file.close();
    }

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

    /*
     * rejected, accepted are clear
     * reordered - for every reorder calculation increase this
     * reordered_worked every time one reordered did help
     * reordered_failed_complete got a worse rmsd after reordering
     * reordered_reused every time an old solution did the trick
     */
    int rejected = 0, accepted = 0, reordered = 0, reordered_worked = 0, reordered_failed_completely = 0, reordered_reused = 0;
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
                    if (!m_silent)
                        std::cout << mol1->Name() << " already rejected. Skipping check against " << mol2->Name() << std::endl;
                    continue;
                }
                RMSDDriver* driver = new RMSDDriver;
                driver->setSilent(true);
                driver->setProtons(!m_heavy);
                driver->setForceReorder(ForceReorder());
                driver->setCheckConnections(CheckConnections());

                if (!m_silent) {
                    std::cout << std::endl
                              << std::setprecision(10)
                              << std::endl
                              << std::endl
                              << "Reference Molecule:" << mol1->Name() << " (" << mol1->Energy() << " Eh)        Target Molecule " << mol2->Name() << " (" << mol2->Energy() << " Eh)" << std::endl;
                }
                if (!m_useRestart) {
                    m_reference_last_energy = mol1->Energy();
                    m_target_last_energy = mol2->Energy();
                }

                double difference = abs(mol1->Energy() - mol2->Energy()) * 2625.5;
                if (m_energy_cutoff > 0)
                    if (difference > m_energy_cutoff) {
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
                if (!m_silent) {
                    std::cout << "Energy Difference: " << std::setprecision(2) << difference << " kJ/mol" << std::endl;
                    std::cout << "Average Difference in rot constant " << std::setprecision(4) << diff_rot << " some unit" << std::endl;
                    std::cout << "RMSD = " << std::setprecision(5) << rmsd << " A" << std::endl;
                }
                if (rmsd > m_rmsd_threshold && difference < 1 && diff_rot < 0.1 && diff_rot > 0.01) {
                    temp_list.push_back(mol2);
                    if (!m_silent) {
                        std::cout << "~~ Reordering forced as energies and rotational constants are too close and rmsd is too different! ~~" << std::endl;
                        std::cout << "*** Adding " << mol2->Name() << " to the list as RMSD is " << rmsd << "! ***" << std::endl;
                    }
                } else {
                    if ((difference < m_energy_threshold && rmsd < m_rmsd_threshold && diff_rot < m_diff_rot_loose)) {
                        ok = false;
                        filtered[mol1->Name()].push_back(mol2->Name());
                        if (!m_silent) {
                            std::cout << "  ** Rejecting structure **" << std::endl;
                        }
                        if (rmsd <= m_rmsd_threshold * m_nearly_missed) {
                            if (!m_silent) {
                                std::cout << " Nearly missed for " << mol1->Name() << std::endl;
                            }
                            m_nearly.push_back(mol1);
                        }
                        continue;
                    }
                    if (diff_rot < m_diff_rot_tight && difference < m_energy_threshold) {
                        ok = false;
                        filtered[mol1->Name()].push_back(mol2->Name());
                        if (!m_silent) {
                            std::cout << "  ** Rejecting structure **" << std::endl;
                        }
                        continue;
                    }
                }
                delete driver;
            }
        }
        TriggerWriteRestart();

        if (temp_list.size()) {
            if (!m_silent) {
                std::cout << std::endl
                          << "             ###   " << std::setprecision(4) << start / double(ende) * 100 << "% done!   ###" << std::endl;
                std::cout << "*** Start Reordering Block **" << std::endl;
            }
        }
        if (!filtered.count(mol1->Name())) {
            for (const auto mol2 : temp_list) {
                if (filtered.count(mol1->Name())) {
                    if (!m_silent) {
                        std::cout << mol1->Name() << " already rejected. Skipping check against " << mol2->Name() << std::endl;
                    }
                    continue;
                }
                if (!m_silent) {
                    std::cout << "*** Reordering " << mol2->Name() << " now! ***" << std::endl;
                }
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
                if (!m_silent) {
                    std::cout << std::endl
                              << std::endl
                              << "*** Reordering forced as energies and rotational constants are too close and rmsd (" << rmsd << ") is too different! ***" << std::endl
                              << std::endl;
                }
                bool digDeeper = true;
                if (!m_silent) {
                    std::cout << "*** Checking old reordering solutions first. ***" << std::endl;
                }
                for (const auto& rule : m_reorder_rules) {
                    double tmp_rmsd = driver->Rules2RMSD(rule);
                    if (!m_silent) {
                        std::cout << tmp_rmsd << " A ";
                    }
                    if (tmp_rmsd < m_rmsd_threshold) {
                        digDeeper = false;
                        if (!m_silent) {
                            std::cout << std::endl
                                      << std::endl
                                      << "*** Old reordering solution worked here! ***" << std::endl
                                      << std::endl;
                        }
                        ok = false;
                        filtered[mol1->Name()].push_back(mol2->Name());
                        reordered_reused++;
                        break;
                    }
                }
                if (!m_silent) {
                    std::cout << std::endl;
                }
                if (m_reference_last_energy - mol1->Energy() > 1e-5 && m_useRestart) {
                    {
                        if (!m_silent) {
                            std::cout << "*** Skip reordering, as we are in Restart Modus *** " << std::endl;
                        }
                        digDeeper = false;
                    }
                } else {
                    m_useRestart = false;
                }
                if (m_prevent_reorder) {
                    digDeeper = false;
                    if (!m_silent) {
                        std::cout << "~~ Reordering prevented, so no more will be done regarding these structures. ~~" << std::endl;
                    }
                }
                if (digDeeper) {

                    {
                        m_reference_last_energy = mol1->Energy();
                        m_target_last_energy = mol2->Energy();
                    }
                    if (!m_silent) {
                        std::cout << "*** Nothing fitting found, reorder now! ***" << std::endl;
                    }
                    driver->setForceReorder(true);
                    driver->AutoPilot();
                    driver->setForceReorder(ForceReorder());
                    reordered++;
                    double rmsd_tmp = driver->RMSD();
                    if (std::find(m_reorder_rules.begin(), m_reorder_rules.end(), driver->ReorderRules()) == m_reorder_rules.end()) {
                        if (!m_silent) {
                            std::cout << "*** Newly obtained reorder solution added to heap of information ***" << std::endl;
                        }
                        m_reorder_rules.push_back(driver->ReorderRules());
                    }
                    TriggerWriteRestart();
                    if (!m_silent) {
                        std::cout << "New rmsd is " << rmsd_tmp << ". Old was " << rmsd << std::endl;
                    }
                    if (rmsd_tmp > rmsd) {
                        if (!m_silent) {
                            std::cout << "Reordering failed, a better fit was not found! This may happen." << std::endl;
                        }
                        reordered_failed_completely++;
                        m_failed.push_back(mol1);
                        m_failed.push_back(mol2);
                    }
                    rmsd = rmsd_tmp;
                }
                if ((difference < m_energy_threshold && rmsd < m_rmsd_threshold && diff_rot < m_diff_rot_loose)) {
                    reordered_worked++;
                    ok = false;
                    filtered[mol1->Name()].push_back(mol2->Name());
                    if (!m_silent) {
                        std::cout << "  ** Rejecting structure **" << std::endl;
                    }
                    if (rmsd <= m_rmsd_threshold * m_nearly_missed) {
                        if (!m_silent) {
                            std::cout << " Nearly missed for " << mol1->Name() << std::endl;
                        }
                        m_nearly.push_back(mol1);
                    }
                    continue;
                }
                if (diff_rot < m_diff_rot_tight && difference < m_energy_threshold) {
                    ok = false;
                    filtered[mol1->Name()].push_back(mol2->Name());
                    if (!m_silent) {
                        std::cout << "  ** Rejecting structure **" << std::endl;
                    }
                    continue;
                }
                delete driver;
            }
        } else {
            if (!m_silent) {
                std::cout << "Skipping check, as structure already rejected." << std::endl;
            }
        }
        if (ok && accepted < m_maxrank && m_maxrank > 0) {
            m_result.push_back(mol1);
            if (m_writeFiles)
                mol1->appendXYZFile(result_name);

            if (!m_silent) {
                std::cout << std::endl
                          << "Having now " << m_result.size() << " structures accepted." << std::endl
                          << "               ** Accepting " << mol1->Name() << " **" << std::endl;
            }
            accepted++;
            /* Write each accepted structure immediately */
        } else
            rejected++;
        ok = true;
        start++;
        std::cout << std::endl
                  << "             ###   " << std::setprecision(4) << start / double(ende) * 100 << "% done!   ###" << std::endl;
        std::cout << "# Accepted : " << accepted << "     ";
        std::cout << "# Rejected : " << rejected << "     ";
        std::cout << "# Reordered : " << reordered << "     ";
        std::cout << "# Successfully : " << reordered_worked << "    ";
        std::cout << "# Not at all : " << reordered_failed_completely << "     ";
        std::cout << "# Reused Results : " << reordered_reused << std::endl;
    }

    if (m_result.size() == 0)
        return;
    if (m_writeFiles) {
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
    }
    {
        m_reference_last_energy = 0;
        m_target_last_energy = 0;
    }
    TriggerWriteRestart();

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
