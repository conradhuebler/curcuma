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

#include <fmt/color.h>
#include <fmt/core.h>

#include "src/capabilities/confstat.h"
#include "src/capabilities/rmsd.h"

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

#include "src/tools/general.h"

#include "json.hpp"
using json = nlohmann::json;

#include "confscan.h"

ConfScan::ConfScan(const json& controller, bool silent)
    : CurcumaMethod(ConfScanJson, controller, silent)
{
    UpdateController(controller);
}

ConfScan::~ConfScan()
{
    for (auto i : m_molecules) {
        delete i.second;
    }
}

void ConfScan::LoadControlJson()
{
    m_noname = Json2KeyWord<bool>(m_defaults, "noname");
    m_restart = Json2KeyWord<bool>(m_defaults, "restart");

    m_heavy = Json2KeyWord<bool>(m_defaults, "heavy");

    m_rmsd_threshold = Json2KeyWord<double>(m_defaults, "rmsd");

    if (m_heavy && m_rmsd_threshold == -1)
        m_rmsd_threshold = 0.75;
    else if (!m_heavy && m_rmsd_threshold == -1)
        m_rmsd_threshold = 1;

    m_maxrank = Json2KeyWord<double>(m_defaults, "rank");
    m_writeXYZ = Json2KeyWord<bool>(m_defaults, "writeXYZ");
    m_force_reorder = Json2KeyWord<bool>(m_defaults, "forceReorder");
    m_check_connections = Json2KeyWord<bool>(m_defaults, "check");
    m_energy_threshold = Json2KeyWord<double>(m_defaults, "energy");
    m_energy_cutoff = Json2KeyWord<double>(m_defaults, "maxenergy");
    m_prevent_reorder = Json2KeyWord<bool>(m_defaults, "preventReorder");
    m_force_silent = Json2KeyWord<bool>(m_defaults, "silent");
    m_scale_loose = Json2KeyWord<double>(m_defaults, "scaleLoose");
    m_scale_tight = Json2KeyWord<double>(m_defaults, "scaleTight");
    m_skip = Json2KeyWord<int>(m_defaults, "skip");
    m_allxyz = Json2KeyWord<bool>(m_defaults, "allxyz");
    m_update = Json2KeyWord<bool>(m_defaults, "update");
    m_maxParam = Json2KeyWord<int>(m_defaults, "MaxParam");
    m_useorders = Json2KeyWord<int>(m_defaults, "UseOrders");
    m_MaxHTopoDiff = Json2KeyWord<int>(m_defaults, "MaxHTopoDiff");
    m_RMSDthreads = Json2KeyWord<int>(m_defaults, "RMSDThreads");

    std::string method = Json2KeyWord<std::string>(m_defaults, "RMSDMethod");

    if (method.compare("template") == 0) {
        m_RMSDmethod = 2;
        if (m_useorders == -1)
            m_useorders = 0;
    } else {
        m_RMSDmethod = 1;
        if (m_useorders == -1)
            m_useorders = 10;
    }
    m_gfn = Json2KeyWord<int>(m_defaults, "GFN");
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
        double energy = mol->Energy();
        if (std::abs(energy) < 1e-5 || m_gfn != -1) {
            XTBInterface interface; // As long as xtb leaks, we have to put it heare
                // I might not leak really, but was unable to clear everything
            if (m_gfn == -1)
                m_gfn = 2;
            interface.InitialiseMolecule(mol);
            energy = interface.GFNCalculation(m_gfn, 0);
            //interface.clear();
        }
        m_ordered_list.insert(std::pair<double, int>(energy, molecule));
        molecule++;
        if (m_noname)
            mol->setName(NamePattern(molecule));
        mol->CalculateRotationalConstants();
        std::pair<std::string, Molecule*> pair(mol->Name(), mol);
        m_molecules.push_back(pair);
    }

    return true;
}

void ConfScan::ReadControlFile()
{
    json control;
    try {
        control = LoadControl();
    } catch (int error) {
        if (error == 404) // No control given
            return;
    }
    json confscan;
    try {
        confscan = control[MethodName()];
    } catch (json::type_error& e) {
        return; // File does not contain control information for ConfScan
    }

    try {
        m_maxrank = confscan["MaxRank"];
    } catch (json::type_error& e) {
    }

    try {
        m_rmsd_threshold = confscan["RMSDThreshold"];
    } catch (json::type_error& e) {
    }
    try {
        m_diff_rot_abs_tight = confscan["AbsRotationalTight"];
        m_parameter_loaded = true;
        m_internal_parametrised = true;
    } catch (json::type_error& e) {
    }
    try {
        m_diff_rot_abs_loose = confscan["AbsRotationalLoose"];
        m_parameter_loaded = true;
        m_internal_parametrised = true;
    } catch (json::type_error& e) {
    }
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
    return true;
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
            m_diff_rot_abs_tight = confscan["AbsRotationalTight"];
            m_parameter_loaded = true;
            m_internal_parametrised = true;
        } catch (json::type_error& e) {
        }
        try {
            m_diff_rot_abs_loose = confscan["AbsRotationalLoose"];
            m_parameter_loaded = true;
            m_internal_parametrised = true;
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
    block["AbsRotationalTight"] = m_diff_rot_abs_tight;
    block["AbsRotationalLoose"] = m_diff_rot_abs_loose;

    return block;
}

void ConfScan::ParametriseRotationalCutoffs()
{
    if (m_prevent_reorder)
        return;
    json rmsd = RMSDJson;
    rmsd["silent"] = true;
    rmsd["reorder"] = false;
    rmsd["noreorder"] = true;
    rmsd["check"] = false;
    rmsd["heavy"] = m_heavy;
    std::cout << "Parametrise cutoff of rotational constants for reordering decison tree!" << std::endl;

    m_internal_parametrised = true;
    std::vector<double> accepted, rejected;
    double accepted_single = 0;
    int counter = 0;
    if (m_maxParam == -1)
        m_maxParam = m_ordered_list.size() * (m_ordered_list.size() - 1);
    if (m_maxParam > 1e4)
        m_maxParam = 1e4;
    std::cout << m_maxParam << " RMSD calculation will be performed." << std::endl;

    for (int i = 0; i < m_ordered_list.size() && counter < m_maxParam; ++i) {
        Molecule* mol1 = m_molecules.at(i).second;
        if (mol1->Check() == 1)
            continue;

        for (int j = i + 1; j < m_ordered_list.size() && counter < m_maxParam; ++j) {
            RMSDDriver* driver = new RMSDDriver(rmsd);
            Molecule* mol2 = m_molecules.at(j).second;
            if (mol2->Check() == 1)
                continue;

            double Ia = abs(mol1->Ia() - mol2->Ia());
            double Ib = abs(mol1->Ib() - mol2->Ib());
            double Ic = abs(mol1->Ic() - mol2->Ic());

            double diff_rot = (Ia + Ib + Ic) * third;
            driver->setReference(mol1);
            driver->setTarget(mol2);

            driver->start();
            double rmsd = driver->RMSD();
            delete driver;
            if (rmsd < m_rmsd_threshold) {
                accepted.push_back(diff_rot);
                accepted_single = std::max(accepted_single, diff_rot);
            } else
                rejected.push_back(diff_rot);
            counter++;
        }
    }
    if (accepted.size() < 2 || rejected.size() < 2) {
        std::cout << " To few samples. Parametrisation rejected!" << std::endl;
        if (!m_force_reorder)
            m_prevent_reorder = true;
        return;
    }
    m_diff_rot_abs_tight = Tools::median(accepted);
    m_diff_rot_abs_loose = Tools::median(rejected);
}

int ConfScan::AcceptRotationalConstant(double constant)
{
    return 0;
    /* This has to be compressed to logic gatters to avoid branching */
    if (m_internal_parametrised) {
        if (constant < m_diff_rot_abs_tight * m_scale_tight)
            return -1; // Assume identical without reordering
        else if (m_diff_rot_abs_tight * m_scale_tight < constant && constant < m_diff_rot_abs_loose * m_scale_loose)
            return 0; // Lets reorder it
        else
            return 1 && m_prevent_reorder == false; // Do not reorder it
    } else {
        if (constant < m_diff_rot_rel_tight)
            return -1; // Assume identical without reordering
        else if (m_diff_rot_rel_tight < constant && constant < m_diff_rot_rel_loose)
            return 0; // Lets reorder it
        else
            return 1 && m_prevent_reorder == false; // Do not reorder it
    }
}

void ConfScan::SetUp()
{
    ReadControlFile();
    LoadRestartInformation();

    m_silent = m_ordered_list.size() > 500 || m_force_silent;

    m_fail = 0;
    m_start = 0;
    m_end = m_ordered_list.size();

    m_result_basename = m_filename;
    m_result_basename.erase(m_result_basename.end() - 4, m_result_basename.end());

    m_accepted_filename = m_result_basename + ".accepted.xyz";
    m_rejected_filename = m_result_basename + ".rejected.xyz";

    std::ofstream result_file;
    if (m_writeFiles) {
        result_file.open(m_accepted_filename);
        result_file.close();
    }

    std::ofstream failed_file;
    if (m_writeFiles) {
        failed_file.open(m_rejected_filename);
        failed_file.close();
    }
    /*
    if (!m_parameter_loaded)
        ParametriseRotationalCutoffs();
*/

    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << "" << std::endl;

    if (m_heavy)
        std::cout << "    RMSD Calculation will be performed only on heavy atoms! " << std::endl;
    else
        std::cout << "    RMSD Calculation will be performed on all atoms! " << std::endl;

    std::cout << "    RMSD Threshold set to: " << m_rmsd_threshold << " Angstrom" << std::endl;
    std::cout << "    Energy Threshold set to: " << m_energy_threshold << " kJ/mol" << std::endl;
    std::cout << "    Thresholds in rotational constants (averaged over Ia, Ib and Ic): " << std::endl;
    std::cout << "    Loose Threshold: " << m_diff_rot_abs_loose << " MHz" << std::endl;
    std::cout << "    Tight Threshold: " << m_diff_rot_abs_tight << " MHz" << std::endl;
    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << std::endl;
}

void ConfScan::start()
{
    PrintController(m_controller);
    SetUp();

    fmt::print("\n\n1st Pass\nPerforming RMSD calculation without reordering now!\n\n");
    RunTimer timer(false);
    CheckRMSD();
    fmt::print("\n1st Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);
    if (m_prevent_reorder == false && CheckStop() == false) {
        timer.Reset();
        fmt::print("\n\n2nd Pass\nPerforming RMSD calculation with reordering now!\n\n");
        ReorderCheck(false, false);
        fmt::print("\n2nd Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);

        timer.Reset();
        fmt::print("\n\n3rd Pass\nPerforming RMSD calculation with reordering, but only reuse previouse reordering rules.\n\n");
        ReorderCheck(true, true);
        fmt::print("\n3rd Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);
    }

    Finalise();
}

void ConfScan::CheckRMSD()
{
    m_maxmol = m_ordered_list.size();

    json rmsd = RMSDJson;
    rmsd["silent"] = true;
    rmsd["check"] = CheckConnections();
    rmsd["heavy"] = m_heavy;
    rmsd["noreorder"] = true;

    for (auto& i : m_ordered_list) {
        if (m_skip) {
            m_skip--;
            continue;
        }
        if (m_prevent_reorder && m_maxrank <= m_accepted)
            continue;

        int index = i.second;
        Molecule* mol1 = m_molecules.at(index).second;
        if (mol1->Check() == 1) {
            m_rejected++;
            m_start++;
            PrintStatus();
            continue;
        }
        if (m_result.size() == 0) {
            m_result.push_back(mol1);
            continue;
        }
        if (m_lowest_energy > 0)
            m_lowest_energy = mol1->Energy();

        bool keep_molecule = true;
        RMSDDriver* driver = new RMSDDriver(rmsd);
        for (const auto& mol2 : m_result) {
            if (CheckStop()) {
                fmt::print("\n\n** Found stop file, will end now! **\n\n");
                delete driver;
                return;
            }
            keep_molecule = SingleCheckRMSD(mol1, mol2, driver);
            if (keep_molecule == false)
                break;
        }
        if (keep_molecule) {
            m_result.push_back(mol1);
            m_accepted++;
        } else {
            m_rejected_structures.push_back(mol1);
            m_rejected++;
        }
        PrintStatus();
        delete driver;
    }
}

bool ConfScan::SingleCheckRMSD(const Molecule* mol1, const Molecule* mol2, RMSDDriver* driver)
{
    bool keep_molecule = true;
    if (!m_silent) {
        std::cout << std::endl
                  << std::setprecision(10)
                  << std::endl
                  << std::endl
                  << "Reference Molecule:" << mol1->Name() << " (" << mol1->Energy() << " Eh)        Target Molecule " << mol2->Name() << " (" << mol2->Energy() << " Eh)" << std::endl;
    }
    //if (!m_useRestart) {
    m_reference_last_energy = mol1->Energy();
    m_target_last_energy = mol2->Energy();
    //}

    double rmsd = 0;
    double Ia = abs(mol1->Ia() - mol2->Ia());
    double Ib = abs(mol1->Ib() - mol2->Ib());
    double Ic = abs(mol1->Ic() - mol2->Ic());

    double diff_rot = (Ia + Ib + Ic) * third;

    driver->setReference(mol1);
    driver->setTarget(mol2);

    driver->start();
    rmsd = driver->RMSD();
    double difference = abs(mol1->Energy() - m_lowest_energy) * 2625.5;
    if (difference > m_energy_cutoff && m_energy_cutoff != -1) {
        keep_molecule = false;
    }

    if (!m_silent) {
        std::cout << "Energy Difference: " << std::setprecision(2) << difference << " kJ/mol" << std::endl;
        std::cout << "Difference in Ia " << std::setprecision(2) << Ia << " MHz" << std::endl;
        std::cout << "Difference in Ib " << std::setprecision(2) << Ib << " MHz" << std::endl;
        std::cout << "Difference in Ic " << std::setprecision(2) << Ic << " MHz" << std::endl;

        std::cout << "RMSD = " << std::setprecision(5) << rmsd << " A" << std::endl;
    }

    if (rmsd <= m_rmsd_threshold && (m_MaxHTopoDiff == -1 || driver->HBondTopoDifference() <= m_MaxHTopoDiff) /*|| accepted_rotational == -1*/) {
        keep_molecule = false;
        std::string reject_reason = mol2->Name() + " [I] RMSD = " + std::to_string(rmsd) + "; dE = " + std::to_string(difference) + "; dIx = " + std::to_string(diff_rot);
        m_filtered[mol1->Name()].push_back(reject_reason);
        m_accept_rmsd.push_back(diff_rot);
        if (!m_silent) {
            std::cout << "  ** Rejecting structure **" << std::endl;
        }
    } else if (rmsd > m_rmsd_threshold) {
        m_reject_rmsd.push_back(diff_rot);
    }
    return keep_molecule;
}

void ConfScan::ReorderCheck(bool reuse_only, bool limit)
{
    m_maxmol = m_result.size();

    m_diff_rot_abs_tight = Tools::median(m_accept_rmsd);
    m_diff_rot_abs_loose = Tools::median(m_reject_rmsd);

    // To be finalised and tested

    /*
    fmt::print(
        "'{0:'^{1}}'\n"
        "'{2: ^{1}}'\n"
        "'{0: ^{1}}'\n"
        "*{3: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{12: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{4: ^{1}}*\n"
        "*{5: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{6: ^{1}}*\n"
        "*{7: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{8: ^{1}}*\n"
        "*{9: ^{1}}*\n"
        "*{10: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{11: ^{1}}*\n"
        "*{0: ^{1}}*\n"
        "*{0:*^{1}}*\n",
        "", 60,
        "Thresholds in rotational constants (averaged over Ia, Ib and Ic):");

    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << "    Thresholds in rotational constants (averaged over Ia, Ib and Ic): " << std::endl;
    std::cout << "    Loose Threshold: " << m_diff_rot_abs_loose << " MHz" << std::endl;
    std::cout << "    Tight Threshold: " << m_diff_rot_abs_tight << " MHz" << std::endl;
    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << std::endl;
*/
    m_rejected = 0, m_accepted = 0, m_reordered = 0, m_reordered_worked = 0, m_reordered_failed_completely = 0, m_reordered_reused = 0;

    json rmsd = RMSDJson;
    rmsd["silent"] = true;
    rmsd["reorder"] = true;
    rmsd["check"] = CheckConnections();
    rmsd["heavy"] = m_heavy;
    rmsd["threads"] = m_RMSDthreads;
    if (m_RMSDmethod == 2)
        rmsd["method"] = "template";
    else
        rmsd["method"] = "incr";

    std::vector<Molecule*> cached = m_result;
    m_result.clear();
    for (Molecule* mol1 : cached) {
        if (m_result.size() == 0) {
            m_result.push_back(mol1);
            continue;
        }

        bool keep_molecule = true;
        RMSDDriver* driver = new RMSDDriver(rmsd);
        for (const auto& mol2 : m_result) {
            if (CheckStop()) {
                delete driver;
                fmt::print("\n\n** Found stop file, will end now! **\n\n");
                return;
            }
            if (keep_molecule == false)
                break;
            keep_molecule = SingleReorderRMSD(mol1, mol2, driver, reuse_only);
            if (keep_molecule == false)
                break;
        }
        if (keep_molecule) {
            m_result.push_back(mol1);
            m_accepted++;
        } else {
            m_rejected_structures.push_back(mol1);
            m_rejected++;
        }

        PrintStatus();
        TriggerWriteRestart();
        delete driver;
        if ((m_result.size() >= m_maxrank && limit) || (m_result.size() >= 2 * m_maxrank && !limit))
            break;
        double difference = abs(mol1->Energy() - m_lowest_energy) * 2625.5;
        if (difference > m_energy_cutoff && m_energy_cutoff != -1) {
            break;
        }
    }
}

bool ConfScan::SingleReorderRMSD(const Molecule* mol1, const Molecule* mol2, RMSDDriver* driver, bool reuse_only)
{
    bool keep_molecule = true;
    if (!m_silent) {
        std::cout << std::endl
                  << std::setprecision(10)
                  << std::endl
                  << std::endl
                  << "Reference Molecule:" << mol1->Name() << " (" << mol1->Energy() << " Eh)        Target Molecule " << mol2->Name() << " (" << mol2->Energy() << " Eh)" << std::endl;
    }

    m_reference_last_energy = mol1->Energy();
    m_target_last_energy = mol2->Energy();

    bool allow_reorder = true;

    double rmsd = 0;
    driver->setReference(mol1);
    driver->setTarget(mol2);

    for (const auto& rule : m_reorder_rules) {
        double tmp_rmsd = driver->Rules2RMSD(rule);
        if (!m_silent) {
            std::cout << tmp_rmsd << " A ";
        }
        if (tmp_rmsd < m_rmsd_threshold && (m_MaxHTopoDiff == -1 || driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
            if (!m_silent) {
                std::cout << std::endl
                          << std::endl
                          << "*** Old reordering solution worked here! ***" << std::endl
                          << std::endl;
            }
            keep_molecule = false;
            allow_reorder = false;

            /*
            if (reuse_only) {
                mol1->writeXYZFile(m_result_basename + ".M1_X" + std::to_string(m_reordered_reused) + ".xyz");
                mol2->writeXYZFile(m_result_basename + ".M2_X" + std::to_string(m_reordered_reused) + ".xyz");
            }
            */
            rmsd = tmp_rmsd;
            m_reordered_reused++;
            break;
        }
    }

    double Ia = abs(mol1->Ia() - mol2->Ia());
    double Ib = abs(mol1->Ib() - mol2->Ib());
    double Ic = abs(mol1->Ic() - mol2->Ic());

    double diff_rot = (Ia + Ib + Ic) * third;

    if (m_reference_last_energy - mol1->Energy() > 1e-5 && m_useRestart) {
        {
            if (!m_silent) {
                std::cout << "*** Skip reordering, as we are in Restart Modus *** " << std::endl;
            }
            allow_reorder = false;
        }
    } else {
        m_useRestart = false;
    }

    if (diff_rot < m_diff_rot_abs_loose && allow_reorder && reuse_only == false) {
        driver->setReference(mol1);
        driver->setTarget(mol2);
        driver->start();
        rmsd = driver->RMSD();
        m_reordered++;
        if (!m_silent) {
            std::cout << "Difference in Ia " << std::setprecision(2) << Ia << " MHz" << std::endl;
            std::cout << "Difference in Ib " << std::setprecision(2) << Ib << " MHz" << std::endl;
            std::cout << "Difference in Ic " << std::setprecision(2) << Ic << " MHz" << std::endl;

            std::cout << "RMSD = " << std::setprecision(5) << rmsd << " A" << std::endl;
        }

        if (rmsd <= m_rmsd_threshold && (m_MaxHTopoDiff == -1 || driver->HBondTopoDifference() <= m_MaxHTopoDiff) /*|| accepted_rotational == -1*/) {
            keep_molecule = false;
            m_reordered_worked++;
            std::string reject_reason = mol2->Name() + " [I] RMSD = " + std::to_string(rmsd) + "; dIx = " + std::to_string(diff_rot);
            m_filtered[mol1->Name()].push_back(reject_reason);
            m_accept_rmsd.push_back(diff_rot);
            if (!m_silent) {
                std::cout << "  ** Rejecting structure **" << std::endl;
            }
            AddRules(driver->ReorderRules());
        } else if (rmsd > m_rmsd_threshold) {
            m_reject_rmsd.push_back(diff_rot);
        }
    }
    return keep_molecule;
}

void ConfScan::Finalise()
{
    int i = 0;
    for (const auto molecule : m_result) {
        double difference = abs(molecule->Energy() - m_lowest_energy) * 2625.5;
        if (i >= m_maxrank && m_maxrank != -1) {
            molecule->appendXYZFile(m_rejected_filename);
            continue;
        }

        if (difference > m_energy_cutoff && m_energy_cutoff != -1) {
            molecule->appendXYZFile(m_rejected_filename);
            continue;
        }
        molecule->appendXYZFile(m_accepted_filename);
        i++;
    }
    for (const auto molecule : m_rejected_structures) {
        molecule->appendXYZFile(m_rejected_filename);
    }
    std::cout << m_result.size() << " structures were kept - of " << m_molecules.size() - m_fail << " total!" << std::endl;
    /*
    //  std::cout << "Best structure is " << m_result[0]->Name() << " Energy = " << std::setprecision(8) << m_result[0]->Energy() << std::endl;
    std::cout << "List of filtered names ... " << std::endl;
    for (const auto& element : m_filtered) {
        std::cout << element.first << " rejected due to: ";
        for (const std::string& str : element.second)
            std::cout << str << " ";
        std::cout << std::endl;
    }
    std::cout << " done :-) " << std::endl;
    */
}

void ConfScan::AddRules(const std::vector<int>& rules)
{
    if (std::find(m_reorder_rules.begin(), m_reorder_rules.end(), rules) == m_reorder_rules.end()) {
        if (!m_silent) {
            std::cout << "*** Newly obtained reorder solution added to heap of information ***" << std::endl;
        }
        m_reorder_rules.push_back(rules);
    }
}


void ConfScan::PrintStatus()
{
    std::cout << std::endl
              << "             ###   " << std::setprecision(4) << (m_result.size() + m_rejected) / double(m_maxmol) * 100 << "% done!   ###" << std::endl;
    std::cout << "# Accepted : " << m_result.size() << "     ";
    std::cout << "# Rejected : " << m_rejected << "     ";
    std::cout << "# Reordered : " << m_reordered << "     ";
    std::cout << "# Successfully : " << m_reordered_worked << "    ";
    std::cout << "# Not at all : " << m_reordered_failed_completely << "     ";
    std::cout << "# Reused Results : " << m_reordered_reused << "     ";
    std::cout << "# Current Energy [kJ/mol] : " << (m_target_last_energy - m_lowest_energy) * 2625.5 << std::endl;
}
