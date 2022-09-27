/*
 * <Scan and judge conformers from different input. >
 * Copyright (C) 2020 - 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "src/capabilities/persistentdiagram.h"
#include "src/capabilities/rmsd.h"

#include "src/core/fileiterator.h"

#include "src/core/tbliteinterface.h"

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
    m_scale_loose = Json2KeyWord<double>(m_defaults, "scaleLoose");
    m_scale_tight = Json2KeyWord<double>(m_defaults, "scaleTight");
    m_last_dE = Json2KeyWord<double>(m_defaults, "lastdE");

    m_skip = Json2KeyWord<int>(m_defaults, "skip");
    m_allxyz = Json2KeyWord<bool>(m_defaults, "allxyz");
    m_update = Json2KeyWord<bool>(m_defaults, "update");
    m_maxParam = Json2KeyWord<int>(m_defaults, "MaxParam");
    m_useorders = Json2KeyWord<int>(m_defaults, "UseOrders");
    m_MaxHTopoDiff = Json2KeyWord<int>(m_defaults, "MaxHTopoDiff");
    m_threads = Json2KeyWord<int>(m_defaults, "threads");
    m_RMSDElement = Json2KeyWord<int>(m_defaults, "RMSDElement");
    m_prev_accepted = Json2KeyWord<std::string>(m_defaults, "accepted");

    m_RMSDmethod = Json2KeyWord<std::string>(m_defaults, "RMSDMethod");
    if (m_useorders == -1)
        m_useorders = 10;

    m_gfn = Json2KeyWord<int>(m_defaults, "GFN");
}

bool ConfScan::openFile()
{
    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;

    if (xyzfile == false)
        throw 1;

    int molecule = 0;
    PersistentDiagram diagram;
    FileIterator file(m_filename);
    while (!file.AtEnd()) {
        Molecule* mol = new Molecule(file.Next());
        double energy = mol->Energy();
        if (std::abs(energy) < 1e-5 || m_gfn != -1) {
            // XTBInterface interface; // As long as xtb leaks, we have to put it heare
            TBLiteInterface interface;
            // I might not leak really, but was unable to clear everything
            if (m_gfn == -1)
                m_gfn = 2;
            interface.InitialiseMolecule(mol);
            energy = interface.GFNCalculation(m_gfn, 0);
        }
        m_ordered_list.insert(std::pair<double, int>(energy, molecule));
        molecule++;
        if (m_noname)
            mol->setName(NamePattern(molecule));
        mol->CalculateRotationalConstants();

        diagram.setDimension(2);
        diagram.setDistanceMatrix(mol->LowerDistanceVector());
        mol->setPersisentImage(diagram.generateImage(diagram.generatePairs()));
        std::pair<std::string, Molecule*> pair(mol->Name(), mol);
        m_molecules.push_back(pair);
    }

    if (m_prev_accepted != "") {
        double min_energy = 0;
        bool xyzfile = std::string(m_prev_accepted).find(".xyz") != std::string::npos || std::string(m_prev_accepted).find(".trj") != std::string::npos;

        if (xyzfile == false)
            throw 1;

        int molecule = 0;

        FileIterator file(m_prev_accepted);
        while (!file.AtEnd()) {
            Molecule* mol = new Molecule(file.Next());
            double energy = mol->Energy();
            if (std::abs(energy) < 1e-5 || m_gfn != -1) {
                TBLiteInterface interface;
                // I might not leak really, but was unable to clear everything
                if (m_gfn == -1)
                    m_gfn = 2;
                interface.InitialiseMolecule(mol);
                energy = interface.GFNCalculation(m_gfn, 0);
            }
            min_energy = std::min(min_energy, energy);
            mol->CalculateRotationalConstants();

            diagram.setDimension(2);
            diagram.setDistanceMatrix(mol->LowerDistanceVector());
            mol->setPersisentImage(diagram.generateImage(diagram.generatePairs()));
            m_previously_accepted.push_back(mol);
        }
        m_lowest_energy = min_energy;
        m_result = m_previously_accepted;
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
        confscan = control[MethodName()[0]];
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
            confscan = restart[MethodName()[0]];
        } catch (json::type_error& e) {
            error++;
            continue;
        }

        try {
            reorder_cached = Tools::String2VectorVector(confscan["ReorderRules"]);
        } catch (json::type_error& e) {
        }
        try {
            m_reference_restored_energy = confscan["ReferenceLastEnergy"];
        } catch (json::type_error& e) {
        }
        try {
            m_target_restored_energy = confscan["TargetLastEnergy"];
        } catch (json::type_error& e) {
        }
        if (m_last_dE < 0) {
            try {
                m_last_dE = confscan["deltaE"];
            } catch (json::type_error& e) {
            }
        }
        for (const auto& vector : reorder_cached)
            if (std::find(m_reorder_rules.begin(), m_reorder_rules.end(), vector) == m_reorder_rules.end())
                m_reorder_rules.push_back(vector);
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
    block["deltaE"] = m_dE;
    return block;
}

void ConfScan::SetUp()
{
    ReadControlFile();
    LoadRestartInformation();

    m_fail = 0;
    m_start = 0;
    m_end = m_ordered_list.size();

    m_result_basename = m_filename;
    m_result_basename.erase(m_result_basename.end() - 4, m_result_basename.end());

    m_accepted_filename = m_result_basename + ".accepted.xyz";
    m_1st_filename = m_result_basename + ".1st.xyz";
    m_2nd_filename = m_result_basename + ".2nd.xyz";

    m_rejected_filename = m_result_basename + ".rejected.xyz";
    m_statistic_filename = m_result_basename + ".statistic.log";
    m_joined_filename = m_result_basename + ".joined.xyz";
    m_threshold_filename = m_result_basename + ".thresh.xyz";
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

    std::ofstream statistic_file;
    if (m_writeFiles) {
        statistic_file.open(m_statistic_filename);
        statistic_file.close();
    }

    std::ofstream thresh_file;
    if (m_writeFiles) {
        thresh_file.open(m_threshold_filename);
        thresh_file.close();
    }

    if (m_previously_accepted.size()) {
        std::ofstream joined_file;
        joined_file.open(m_joined_filename);
        joined_file.close();
    }

    std::ofstream st_file;
    if (m_writeFiles) {
        st_file.open(m_1st_filename);
        st_file.close();
    }

    std::ofstream nd_file;
    if (m_writeFiles) {
        nd_file.open(m_2nd_filename);
        nd_file.close();
    }

    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << "" << std::endl;

    if (m_heavy)
        std::cout << "    RMSD Calculation will be performed only on heavy atoms! " << std::endl;
    else
        std::cout << "    RMSD Calculation will be performed on all atoms! " << std::endl;

    std::cout << "    RMSD Threshold set to: " << m_rmsd_threshold << " Angstrom" << std::endl;
    std::cout << "    Highest energy conformer allowed: " << m_energy_cutoff << " kJ/mol " << std::endl;
    //  std::cout << "    Energy Threshold set to: " << m_energy_threshold << " kJ/mol" << std::endl;
    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << std::endl;
}

void ConfScan::AcceptMolecule(Molecule* molecule)
{
    m_result.push_back(molecule);
    m_stored_structures.push_back(molecule);
    m_accepted++;
}
void ConfScan::RejectMolecule(Molecule* molecule)
{
    m_rejected_structures.push_back(molecule);
    m_rejected++;
}

void ConfScan::start()
{
    PrintController(m_controller);
    SetUp();

    fmt::print("\n\n1st Pass\nPerforming RMSD calculation without reordering now!\n\n");
    RunTimer timer(false);
    std::ofstream result_file;
    result_file.open(m_statistic_filename, std::ios_base::app);
    result_file << "Results of 1st Pass" << std::endl;
    result_file.close();
    CheckRMSD();

    for (const auto molecule : m_stored_structures) {
        molecule->appendXYZFile(m_1st_filename);
    }

    fmt::print("\n1st Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);
    if (m_prevent_reorder == false) {
        if (!CheckStop()) {
            timer.Reset();
            fmt::print("\n\n2nd Pass\nPerforming RMSD calculation with reordering now!\n\n");
            result_file.open(m_statistic_filename, std::ios_base::app);
            result_file << "Results of 2nd Pass" << std::endl;
            result_file.close();
            ReorderCheck(false, false);
            for (const auto molecule : m_stored_structures) {
                molecule->appendXYZFile(m_2nd_filename);
            }
            fmt::print("\n2nd Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);
            timer.Reset();
        }
        if (!CheckStop()) {
            result_file.open(m_statistic_filename, std::ios_base::app);
            result_file << "Results of 3rd Pass" << std::endl;
            result_file.close();
            fmt::print("\n\n3rd Pass\nPerforming RMSD calculation with reordering, but only reuse previouse reordering rules.\n\n");
            ReorderCheck(true, true);
            fmt::print("\n3rd Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);
        }
    }
    if (!CheckStop())
        m_dE = -1;

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
        if (/*m_prevent_reorder && */ m_maxrank <= m_accepted && m_maxrank > -1)
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
            AcceptMolecule(mol1);
            m_lowest_energy = mol1->Energy();
            continue;
        }
        m_current_energy = mol1->Energy();
        m_dE = (m_current_energy - m_lowest_energy) * 2625.5;
        bool keep_molecule = true;
        RMSDDriver* driver = new RMSDDriver(rmsd);
        for (const auto& mol2 : m_result) {
            if (CheckStop()) {
                fmt::print("\n\n** Found stop file, will end now! **\n\n");
                delete driver;
                return;
            }
            keep_molecule = SingleCheckRMSD(mol1, mol2, driver);
            if (keep_molecule == false) {
                writeStatisticFile(driver->ReferenceAlignedReference(), driver->TargetAlignedReference(), driver->RMSD());
                break;
            }
        }
        if (keep_molecule) {
            AcceptMolecule(mol1);
        } else {
            RejectMolecule(mol1);
        }
        PrintStatus();
        delete driver;
    }
}

bool ConfScan::SingleCheckRMSD(const Molecule* mol1, const Molecule* mol2, RMSDDriver* driver)
{
    bool keep_molecule = true;

    double rmsd = 0;
    double Ia = abs(mol1->Ia() - mol2->Ia());
    double Ib = abs(mol1->Ib() - mol2->Ib());
    double Ic = abs(mol1->Ic() - mol2->Ic());

    double diff_rot = (Ia + Ib + Ic) * third;
    if (m_dE > m_energy_cutoff && m_energy_cutoff != -1) {
        keep_molecule = false;

        m_reference_last_energy = mol1->Energy();
        m_target_last_energy = mol2->Energy();

        return false;
    }

    driver->setReference(mol1);
    driver->setTarget(mol2);

    driver->start();
    rmsd = driver->RMSD();
    double diff = (mol1->getPersisentImage() - mol2->getPersisentImage()).cwiseAbs().sum();
    if (rmsd <= m_scale_tight * m_rmsd_threshold) {
        m_diff_rot_threshold_tight = std::max(m_diff_rot_threshold_tight, diff_rot);
        m_diff_ripser_threshold_tight = std::max(m_diff_ripser_threshold_tight, diff);
    } else if (rmsd <= m_scale_loose * m_rmsd_threshold && rmsd > m_scale_tight * m_rmsd_threshold) {
        m_diff_rot_threshold_loose = std::max(m_diff_rot_threshold_loose, diff_rot);
        m_diff_ripser_threshold_loose = std::max(m_diff_ripser_threshold_loose, diff);
    }
    /*
        std::ofstream result_file;
        result_file.open("ripser.dat", std::ios_base::app);
        result_file << rmsd << "\t" << diff << "\t" << diff_rot << std::endl;
        result_file.close();
    */

    if (rmsd <= m_rmsd_threshold && (m_MaxHTopoDiff == -1 || driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
        keep_molecule = false;
    }

    m_reference_last_energy = mol1->Energy();
    m_target_last_energy = mol2->Energy();

    return keep_molecule;
}

void ConfScan::ReorderCheck(bool reuse_only, bool limit)
{
    m_maxmol = m_stored_structures.size();

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
            "Thresholds in rotational constants (averaged over Ia, Ib and Ic) and Ripser Image:");
    */
    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << "    Thresholds in rotational constants (averaged over Ia, Ib and Ic): " << std::endl;
    std::cout << "    Loose Threshold: " << m_diff_rot_threshold_loose << " MHz" << std::endl;
    std::cout << "    Tight Threshold: " << m_diff_rot_threshold_tight << " MHz" << std::endl;

    std::cout << "    Thresholds in difference of ripser images: " << std::endl;
    std::cout << "    Loose Threshold: " << m_diff_ripser_threshold_loose << " " << std::endl;
    std::cout << "    Tight Threshold: " << m_diff_ripser_threshold_tight << " " << std::endl;

    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << std::endl;

    m_rejected = 0, m_accepted = 0, m_reordered = 0, m_reordered_worked = 0, m_reordered_reused = 0;

    json rmsd = RMSDJson;
    rmsd["silent"] = true;
    rmsd["reorder"] = true;
    rmsd["check"] = CheckConnections();
    rmsd["heavy"] = m_heavy;
    rmsd["threads"] = m_threads;
    rmsd["method"] = m_RMSDmethod;
    rmsd["element"] = m_RMSDElement;

    std::vector<Molecule*> cached = m_stored_structures;
    m_result = m_previously_accepted;
    m_stored_structures.clear();
    for (Molecule* mol1 : cached) {
        if (m_result.size() == 0) {
            AcceptMolecule(mol1);
            continue;
        }
        m_current_energy = mol1->Energy();
        m_dE = (m_current_energy - m_lowest_energy) * 2625.5;

        bool keep_molecule = true;
        RMSDDriver* driver = new RMSDDriver(rmsd);
        for (const auto& mol2 : m_result) {
            if (CheckStop()) {
                delete driver;
                fmt::print("\n\n** Found stop file, will end now! **\n\n");
                // TriggerWriteRestart();
                return;
            }

            double Ia = abs(mol1->Ia() - mol2->Ia());
            double Ib = abs(mol1->Ib() - mol2->Ib());
            double Ic = abs(mol1->Ic() - mol2->Ic());

            double diff_rot = (Ia + Ib + Ic) * third;
            double diff = (mol1->getPersisentImage() - mol2->getPersisentImage()).cwiseAbs().sum();

            if (diff_rot > m_diff_rot_threshold_loose && diff > m_diff_ripser_threshold_loose) {
                keep_molecule = true;
                break;
            }
            if (diff_rot < m_diff_rot_threshold_tight && diff < m_diff_ripser_threshold_tight) {
                std::cout << "Differences " << diff_rot << " MHz and " << diff << " below tight threshold, reject molecule directly!" << std::endl;
                m_last_diff = diff_rot;
                m_last_ripser = diff;
                keep_molecule = false;
                writeStatisticFile(mol1, mol2, driver->RMSD(), false);
                m_threshold.push_back(mol2);
                break;
            }

            keep_molecule = SingleReorderRMSD(mol1, mol2, driver, reuse_only);
            if (keep_molecule == false) {
                writeStatisticFile(driver->ReferenceAlignedReference(), driver->TargetAlignedReference(), driver->RMSD(), true);
                /*
                std::ofstream result_file;
                result_file.open("ripser.dat", std::ios_base::app);
                result_file << driver->RMSD() << "\t" << diff << "\t" << diff_rot << std::endl;
                result_file.close();
                */
                break;
            }
        }
        if (keep_molecule) {
            AcceptMolecule(mol1);
        } else {
            RejectMolecule(mol1);
        }

        PrintStatus();
        delete driver;
        if ((m_result.size() >= m_maxrank && limit) || (m_result.size() >= 2 * m_maxrank && !limit))
            break;
        if (m_dE > m_energy_cutoff && m_energy_cutoff != -1) {
            break;
        }
    }
}

bool ConfScan::SingleReorderRMSD(const Molecule* mol1, const Molecule* mol2, RMSDDriver* driver, bool reuse_only)
{
    bool keep_molecule = true;
    bool allow_reorder = true;

    double rmsd = 0;
    driver->setReference(mol1);
    driver->setTarget(mol2);
    for (const auto& rule : m_reorder_rules) {
        if (rule.size() != mol1->AtomCount())
            continue;

        double tmp_rmsd = driver->Rules2RMSD(rule);
        if (tmp_rmsd < m_rmsd_threshold && (m_MaxHTopoDiff == -1 || driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
            keep_molecule = false;
            allow_reorder = false;
            rmsd = tmp_rmsd;
            m_reordered_reused++;
            writeStatisticFile(driver->ReferenceAlignedReference(), driver->TargetAlignedReference(), driver->RMSD(), true);
            break;
        }
    }
    if (m_useRestart && m_dE < m_last_dE) {
        allow_reorder = false;
    } else {
        m_useRestart = false;
    }

    if (allow_reorder && reuse_only == false) {
        driver->setReference(mol1);
        driver->setTarget(mol2);
        driver->start();
        rmsd = driver->RMSD();

        m_reordered++;

        if (rmsd <= m_rmsd_threshold && (m_MaxHTopoDiff == -1 || driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
            keep_molecule = false;
            AddRules(driver->ReorderRules());
            m_reordered_worked++;
        }
    }

    m_reference_last_energy = mol1->Energy();
    m_target_last_energy = mol2->Energy();

    return keep_molecule;
}

void ConfScan::Finalise()
{
    TriggerWriteRestart();

    int i = 0;
    for (const auto molecule : m_stored_structures) {
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
        if (m_previously_accepted.size()) {
            molecule->appendXYZFile(m_joined_filename);
        }
        i++;
    }

    for (const auto molecule : m_previously_accepted) {
        molecule->appendXYZFile(m_joined_filename);
    }

    for (const auto molecule : m_rejected_structures) {
        molecule->appendXYZFile(m_rejected_filename);
    }

    for (const auto molecule : m_threshold)
        molecule->appendXYZFile(m_threshold_filename);

    std::cout << m_stored_structures.size() << " structures were kept - of " << m_molecules.size() - m_fail << " total!" << std::endl;
}

void ConfScan::AddRules(const std::vector<int>& rules)
{
    if (std::find(m_reorder_rules.begin(), m_reorder_rules.end(), rules) == m_reorder_rules.end()) {
        m_reorder_rules.push_back(rules);
    }
}


void ConfScan::PrintStatus()
{
    std::cout << std::endl
              << "             ###   " << std::setprecision(4) << (m_stored_structures.size() + m_rejected) / double(m_maxmol) * 100 << "% done!   ###" << std::endl;
    std::cout << "# Accepted : " << m_stored_structures.size() << "     ";
    std::cout << "# Rejected : " << m_rejected << "     ";
    std::cout << "# Reordered : " << m_reordered << "     ";
    std::cout << "# Successfully : " << m_reordered_worked << "    ";
    std::cout << "# Reused Results : " << m_reordered_reused << "     ";
    std::cout << "# Current Energy [kJ/mol] : " << m_dE << std::endl;
}

void ConfScan::writeStatisticFile(const Molecule* mol1, const Molecule* mol2, double rmsd, bool reason)
{
    std::ofstream result_file;
    result_file.open(m_statistic_filename, std::ios_base::app);
    if (reason)
        result_file << "Molecule got rejected due to small rmsd " << rmsd << " with and energy difference of " << m_dE << std::endl;
    else
        result_file << "Molecule got rejected as differences " << m_last_diff << " MHz and " << m_last_ripser << " are below the estimated thresholds;  with and energy difference of " << std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 << std::endl;

    result_file << mol1->XYZString();
    result_file << mol2->XYZString();
    result_file << std::endl;
    result_file.close();
}
