/*
 * <Scan and judge conformers from different input. >
 * Copyright (C) 2020 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include <fmt/color.h>
#include <fmt/core.h>

#include "src/capabilities/confstat.h"
#include "src/capabilities/persistentdiagram.h"
#include "src/capabilities/rmsd.h"

#include "src/core/fileiterator.h"

#include "src/core/energycalculator.h"

#include "src/tools/general.h"

#include "json.hpp"
using json = nlohmann::json;

#include "confscan.h"

int ConfScanThread::execute()
{
    m_driver->setThreads(m_threads);
    m_driver->setReference(m_reference);
    m_driver->setTarget(m_target);
    m_keep_molecule = true;
    m_break_pool = false;
    m_reorder_worked = false;
    m_reused_worked = false;
    m_reorder_rule.clear();

    double Ia = abs(m_reference.Ia() - m_target.Ia());
    double Ib = abs(m_reference.Ib() - m_target.Ib());
    double Ic = abs(m_reference.Ic() - m_target.Ic());

#ifdef WriteMoreInfo
    m_input.dIa = Ia;
    m_input.dIb = Ib;
    m_input.dIc = Ic;
    const Matrix m = (m_reference.getPersisentImage() - m_target.getPersisentImage());
    double diff_ripser = m.cwiseAbs().sum();
    m_input.dH = diff_ripser;
    m_input.dHM = m;
    m_input.dE = std::abs(m_reference.Energy() - m_target.Energy()) * 2625.5;
#endif

    for (int i = 0; i < m_reorder_rules.size(); ++i) {
        if (m_reorder_rules[i].size() != m_reference.AtomCount() || m_reorder_rules[i].size() == 0)
            continue;

        double tmp_rmsd = m_driver->Rules2RMSD(m_reorder_rules[i]);
        if (tmp_rmsd < m_rmsd_threshold && (m_MaxHTopoDiff == -1 || m_driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
            m_keep_molecule = false;
            m_break_pool = true;
            m_reused_worked = true;
            m_rmsd = tmp_rmsd;

#ifdef WriteMoreInfo
            m_input.rmsd = m_rmsd;
#endif
            m_reorder_rule = m_reorder_rules[i];
            return 0;
        }
    }

    if (m_reuse_only) {
        return 0;
    }

    m_driver->start();
    m_rmsd = m_driver->RMSD();

#ifdef WriteMoreInfo
    m_input.rmsd = m_rmsd;
#endif

    if (m_rmsd <= m_rmsd_threshold && (m_MaxHTopoDiff == -1 || m_driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
        m_keep_molecule = false;
        m_break_pool = true;
        m_reorder_worked = true;

        m_reorder_rule = m_driver->ReorderRules();
    }
    m_driver->clear();
    return 0;
}

int ConfScanThreadNoReorder::execute()
{
    m_driver->setReference(m_reference);
    m_driver->setTarget(m_target);
    m_keep_molecule = true;

    m_driver->start();
    m_rmsd = m_driver->RMSD();
#ifdef WriteMoreInfo
    m_input.rmsd = m_rmsd;
#endif

    double Ia = abs(m_reference.Ia() - m_target.Ia());
    double Ib = abs(m_reference.Ib() - m_target.Ib());
    double Ic = abs(m_reference.Ic() - m_target.Ic());
#ifdef WriteMoreInfo
    m_input.dIa = Ia;
    m_input.dIb = Ib;
    m_input.dIc = Ic;
#endif
    m_DI = (Ia + Ib + Ic) * third;
    const Matrix m = (m_reference.getPersisentImage() - m_target.getPersisentImage());
    m_DH = m.cwiseAbs().sum();
#ifdef WriteMoreInfo
    m_input.dH = m_DH;
    m_input.dHM = m;
    m_input.dE = std::abs(m_reference.Energy() - m_target.Energy()) * 2625.5;
#endif
    if (m_rmsd <= m_rmsd_threshold && (m_MaxHTopoDiff == -1 || m_driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
        m_keep_molecule = false;
        m_break_pool = true;
    }

    m_driver->clear();
    return 0;
}

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
        m_rmsd_threshold = 0.9;

    m_maxrank = Json2KeyWord<double>(m_defaults, "rank");
    m_writeXYZ = Json2KeyWord<bool>(m_defaults, "writeXYZ");
    m_force_reorder = Json2KeyWord<bool>(m_defaults, "forceReorder");
    m_check_connections = Json2KeyWord<bool>(m_defaults, "check");
    m_energy_cutoff = Json2KeyWord<double>(m_defaults, "maxenergy");
    m_prevent_reorder = Json2KeyWord<bool>(m_defaults, "preventReorder");

    if (m_defaults["sLE"].is_number())
        m_sLE = { Json2KeyWord<double>(m_defaults, "sLE") };
    else if (m_defaults["sLE"].is_string())
        m_sLE = Tools::String2DoubleVec(Json2KeyWord<std::string>(m_defaults, "sLE"), ";");

    if (m_defaults["sLI"].is_number())
        m_sLI = { Json2KeyWord<double>(m_defaults, "sLI") };
    else if (m_defaults["sLI"].is_string())
        m_sLI = Tools::String2DoubleVec(Json2KeyWord<std::string>(m_defaults, "sLI"), ";");

    if (m_defaults["sLH"].is_number())
        m_sLH = { Json2KeyWord<double>(m_defaults, "sLH") };
    else if (m_defaults["sLH"].is_string())
        m_sLH = Tools::String2DoubleVec(Json2KeyWord<std::string>(m_defaults, "sLH"), ";");

    if (m_sLE.size() != m_sLI.size() || m_sLE.size() != m_sLH.size()) {
        std::cout << "Inconsistent length of steps requested, will abort now." << std::endl;
        exit(1);
    }

    m_noreorderpass = Json2KeyWord<bool>(m_defaults, "noreorderpass");
    n_doreusepass = Json2KeyWord<bool>(m_defaults, "doreuse");

    m_sTE = Json2KeyWord<double>(m_defaults, "sTE");
    m_sTI = Json2KeyWord<double>(m_defaults, "sTI");
    m_sTH = Json2KeyWord<double>(m_defaults, "sTH");

    m_lastdE = Json2KeyWord<double>(m_defaults, "lastdE");
    m_domolalign = Json2KeyWord<double>(m_defaults, "domolalign");

    m_skip = Json2KeyWord<int>(m_defaults, "skip");
    m_allxyz = Json2KeyWord<bool>(m_defaults, "allxyz");
    m_reduced_file = Json2KeyWord<bool>(m_defaults, "fewerFile");

    m_update = Json2KeyWord<bool>(m_defaults, "update");
    m_maxParam = Json2KeyWord<int>(m_defaults, "MaxParam");
    m_useorders = Json2KeyWord<int>(m_defaults, "UseOrders");
    m_MaxHTopoDiff = Json2KeyWord<int>(m_defaults, "MaxHTopoDiff");
    m_threads = m_defaults["threads"].get<int>();
    m_skipfirst = Json2KeyWord<bool>(m_defaults, "skipfirst");
    m_RMSDmethod = Json2KeyWord<std::string>(m_defaults, "RMSDMethod");
    if (m_RMSDmethod == "molalign") {
        fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nPlease cite the follow research report!\nJ. Chem. Inf. Model. 2023, 63, 4, 1157–1165 - DOI: 10.1021/acs.jcim.2c01187\n\n");
        m_domolalign = -1;
    }
    m_ignoreRotation = Json2KeyWord<bool>(m_defaults, "ignoreRotation");
    m_ignoreBarCode = Json2KeyWord<bool>(m_defaults, "ignoreBarCode");
    m_update_rotation = Json2KeyWord<bool>(m_defaults, "update-rotation");
    m_split = Json2KeyWord<bool>(m_defaults, "split");
    m_write = Json2KeyWord<bool>(m_defaults, "writefiles");

    m_nomunkres = Json2KeyWord<bool>(m_defaults, "nomunkres");
    m_molalign = Json2KeyWord<std::string>(m_defaults, "molalignbin");
    m_molaligntol = Json2KeyWord<int>(m_defaults, "molaligntol");

    m_openLoop = true;
    m_closeLoop = false;

    m_looseThresh = m_defaults["looseThresh"];
    m_tightThresh = m_defaults["tightThresh"];

#pragma message("these hacks to overcome the json stuff are not nice, TODO!")
    try {
        m_rmsd_element_templates = m_defaults["RMSDElement"].get<std::string>();
        StringList elements = Tools::SplitString(m_rmsd_element_templates, ";");
        for (const std::string& str : elements)
            m_rmsd_element_templates.push_back(std::stod(str));

        if (m_element_templates.size())
            m_RMSDElement = m_element_templates[0];

    } catch (const nlohmann::detail::type_error& error) {
        m_RMSDElement = Json2KeyWord<int>(m_defaults, "RMSDElement");
        m_rmsd_element_templates.push_back(m_RMSDElement);
    }
    if (m_RMSDmethod.compare("hybrid") == 0 /* && m_rmsd_element_templates.compare("") == 0 */) {
        std::cout << "Reordering method hybrid has to be combined with element types. I will chose for you nitrogen and oxygen!" << std::endl;
        std::cout << "This is equivalent to adding:\' -rmsdelement 7;8 \' to your argument list!" << std::endl;
        m_rmsd_element_templates = "7;8";
    }
    m_prev_accepted = Json2KeyWord<std::string>(m_defaults, "accepted");

    if (m_useorders == -1)
        m_useorders = 10;
}

bool ConfScan::openFile()
{
    bool xyzfile = std::string(m_filename).find(".xyz") != std::string::npos || std::string(m_filename).find(".trj") != std::string::npos;

    if (xyzfile == false)
        throw 1;

    int molecule = 0;
    PersistentDiagram diagram(m_defaults);
    FileIterator file(m_filename);
    while (!file.AtEnd()) {
        Molecule* mol = new Molecule(file.Next());
        double energy = mol->Energy();
        if (std::abs(energy) < 1e-5 || m_method.compare("") != 0) {
            // XTBInterface interface; // As long as xtb leaks, we have to put it heare
            if (m_method == "")
                m_method = "gfn2";
            EnergyCalculator interface(m_method, m_controller);
            // I might not leak really, but was unable to clear everything

            interface.setMolecule(mol);
            energy = interface.CalculateEnergy(false);
        }
        m_ordered_list.insert(std::pair<double, int>(energy, molecule));
        molecule++;
        if (m_noname)
            mol->setName(NamePattern(molecule));

        mol->CalculateRotationalConstants();
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
            if (std::abs(energy) < 1e-5 || m_method.compare("") != 0) {
                // XTBInterface interface; // As long as xtb leaks, we have to put it heare
                if (m_method == "")
                    m_method = "gfn2";
                EnergyCalculator interface(m_method, m_controller);
                // I might not leak really, but was unable to clear everything

                interface.setMolecule(mol);
                energy = interface.CalculateEnergy(false);
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
        if (m_lastdE < 0) {
            try {
                m_lastdE = confscan["deltaE"];
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
    block["dLI"] = m_dLI;
    block["dLH"] = m_dLH;
    block["dLE"] = m_dLE;
    block["dTI"] = m_dTI;
    block["dTH"] = m_dTH;
    block["dTE"] = m_dTE;

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
    m_1st_filename = m_result_basename + ".initial.xyz";
    m_2nd_filename = m_result_basename + ".reorder";
    m_3rd_filename = m_result_basename + ".reuse.xyz";

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
    if (m_writeFiles && !m_reduced_file) {
        failed_file.open(m_rejected_filename);
        failed_file.close();
    }

    std::ofstream statistic_file;
    if (m_writeFiles && !m_reduced_file) {
        statistic_file.open(m_statistic_filename);
        statistic_file.close();
    }

    std::ofstream thresh_file;
    if (m_writeFiles && !m_reduced_file) {
        thresh_file.open(m_threshold_filename);
        thresh_file.close();
    }

    if (m_previously_accepted.size()) {
        std::ofstream joined_file;
        joined_file.open(m_joined_filename);
        joined_file.close();
    }

    std::ofstream st_file;
    if (m_writeFiles && !m_reduced_file) {
        st_file.open(m_1st_filename);
        st_file.close();
    }

    std::ofstream nd_file;
    if (m_writeFiles && !m_reduced_file) {
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
    std::cout << "    Threshold multipliers are loose / tight " << std::endl;

    std::cout << "    Ripser Persitance Diagrams definition for loose ";
    for (const auto& d : m_sLH)
        std::cout << d << " ";
    std::cout << " and tight thresholds " << m_sTH << std::endl;
    std::cout << "    Rotational Constants definition for loose ";
    for (const auto& d : m_sLI)
        std::cout << d << " ";
    std::cout << " and tight thresholds " << m_sTI << std::endl;

    std::cout << "    Energy definition for loose ";
    for (const auto& d : m_sLE)
        std::cout << d << " ";
    std::cout << " and tight thresholds " << m_sTE << std::endl;

    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << std::endl;
}

void ConfScan::AcceptMolecule(Molecule* molecule)
{
    m_result.push_back(molecule);
    m_stored_structures.push_back(molecule);
    m_accepted++;
    if (m_writeFiles && !m_reduced_file && m_current_filename.length()) {
        molecule->appendXYZFile(m_current_filename);
    }
}

void ConfScan::RejectMolecule(Molecule* molecule)
{
    m_rejected_structures.push_back(molecule);
    m_rejected++;
}

void ConfScan::WriteDotFile(const std::string& filename, const std::string& content)
{
    std::ofstream result_file;
    result_file.open(filename);
    result_file << "digraph graphname \n {\n";
    result_file << content << std::endl;
    result_file << "}";
}

void ConfScan::start()
{
    PrintController(m_controller);
    SetUp();

    fmt::print("\n\n1st Pass\nPerforming RMSD calculation without reordering now!\n\n");
    RunTimer timer(false);
    m_current_filename = m_1st_filename;
    std::ofstream result_file;
    if (m_writeFiles && !m_reduced_file) {
        result_file.open(m_statistic_filename, std::ios_base::app);
        result_file << "Results of 1st Pass" << std::endl;
        result_file.close();
    }
    if (!m_skipfirst) {
        CheckOnly(m_sLE[0], m_sLI[0], m_sLH[0]);
        PrintStatus("Result 1st pass:");
        WriteDotFile(m_result_basename + ".initial.dot", m_first_content);
        m_sLE[0] = 1;
        m_sLI[0] = 1;
        m_sLH[0] = 1;
        m_collective_content = "edge [color=green];\n";
        m_collective_content += m_first_content + "\n";
    } else {
        for (const auto& i : m_ordered_list)
            m_stored_structures.push_back(m_molecules.at(i.second).second);
    }
#ifdef WriteMoreInfo
    std::cout << m_dnn_data.size() << std::endl;
#endif
    fmt::print("\nInitial Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);

    if (!m_noreorderpass) {
        for (int run = 0; run < m_sLE.size(); ++run) {
            double dLI = m_dLI;
            double dLH = m_dLH;
            double dLE = m_dLE;

            dLI = m_dLI * m_sLI[run];
            dLH = m_dLH * m_sLH[run];
            dLE = m_dLE * m_sLE[run];

            if (!CheckStop()) {
                timer.Reset();
                m_current_filename = m_2nd_filename + "." + std::to_string(run + 1) + ".xyz";
                fmt::print("\n\nReorder Pass\nPerforming RMSD calculation with reordering now!\n\n");
                if (m_writeFiles && !m_reduced_file) {
                    result_file.open(m_statistic_filename, std::ios_base::app);
                    result_file << "Results of Reorder Pass #" << run + 1 << std::endl;
                    result_file.close();
                }
                Reorder(dLE, dLI, dLH, false);
                PrintStatus("Result Reorder pass:");
                WriteDotFile(m_result_basename + ".reorder." + std::to_string(run + 1) + ".dot", m_second_content);
                fmt::print("\nReorder Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);
                timer.Reset();
                m_collective_content += "edge [color=red];\n";
                m_collective_content += m_second_content + "\n";
                m_second_content.clear();
            }
        }
    } else
        fmt::print("\nReorder Pass skipped!\n");

    if (n_doreusepass) {
        if (!CheckStop()) {
            timer.Reset();
            m_current_filename = m_3rd_filename;
            fmt::print("\n\nReuse Pass\nPerforming RMSD calculation with stored reorder rules using all structures.\n\n");
            if (m_writeFiles && !m_reduced_file) {
                result_file.open(m_statistic_filename, std::ios_base::app);
                PrintStatus("Result Reuse pass:");
                result_file.close();
            }
            Reorder(-1, -1, -1, true);
            PrintStatus("Result reuse pass:");
            WriteDotFile(m_result_basename + ".reuse.dot", m_second_content);
            fmt::print("\nReuse Pass finished after {} seconds!\n", timer.Elapsed() / 1000.0);
            timer.Reset();
            m_collective_content += "edge [color=blue];\n";
            m_collective_content += m_second_content + "\n";
        }
    } else
        fmt::print("\Reuse Pass skipped!\n");

#ifdef WriteMoreInfo
    int index = 1;
    for (const auto& i : m_dnn_data) {
        json data;
        data["xcount"] = 5;
        data["ycount"] = 1;
        data["Xcount"] = 1;
        data["y1"] = i.rmsd;
        data["x1"] = i.dE;
        data["x2"] = i.dIa;
        data["x3"] = i.dIb;
        data["x4"] = i.dIc;
        data["x5"] = i.dH;
        data["X1"] = Tools::Matrix2String(i.dHM);
        std::ofstream restart_file("confscan_saved_" + std::to_string(index) + ".json");
        restart_file << data;
        index++;
    }
#endif
    std::ofstream dotfile;
    dotfile.open(m_result_basename + ".dot");
    dotfile << "digraph graphname \n {\n";
    dotfile << m_collective_content;
    dotfile << "}";
    Finalise();
}

void ConfScan::CheckOnly(double sLE, double sLI, double sLH)
{
    std::string content_dot;
    std::string laststring;
    m_maxmol = m_ordered_list.size();

    json rmsd = RMSDJson;
    rmsd["silent"] = true;
    rmsd["check"] = CheckConnections();
    rmsd["heavy"] = m_heavy;
    rmsd["noreorder"] = true;

    std::vector<ConfScanThreadNoReorder*> threads;
    std::vector<std::vector<int>> rules;
    CxxThreadPool* p = new CxxThreadPool;
    p->setActiveThreadCount(m_threads);

    for (auto& i : m_ordered_list) {
        if (m_skip) {
            m_skip--;
            continue;
        }
        if (m_maxrank <= m_accepted && m_maxrank > -1)
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
            ConfScanThreadNoReorder* thread = addThreadNoreorder(mol1, rmsd);
            threads.push_back(thread);
            p->addThread(thread);
            m_all_structures.push_back(mol1);

            m_lowest_energy = mol1->Energy();
            continue;
        }
        p->Reset();
        m_current_energy = mol1->Energy();
        m_dE = (m_current_energy - m_lowest_energy) * 2625.5;

        bool keep_molecule = true;
        for (int i = 0; i < threads.size(); ++i) {
            threads[i]->setTarget(mol1);
        }
        p->StaticPool();
        p->StartAndWait();

        for (auto* t : threads) {

            m_dTI = std::max(m_dTI, t->DI() * (t->RMSD() <= (m_sTI * m_rmsd_threshold)));
            m_dTH = std::max(m_dTH, t->DH() * (t->RMSD() <= (m_sTH * m_rmsd_threshold)));
            m_dTE = std::max(m_dTE, (std::abs(t->Reference()->Energy() - mol1->Energy()) * 2625.5) * (t->RMSD() <= (m_sTE * m_rmsd_threshold)));

            m_dLI = std::max(m_dLI, t->DI() * (t->RMSD() <= (sLI * m_rmsd_threshold)));
            m_dLH = std::max(m_dLH, t->DH() * (t->RMSD() <= (sLH * m_rmsd_threshold)));
            m_dLE = std::max(m_dLE, (std::abs(t->Reference()->Energy() - mol1->Energy()) * 2625.5) * (t->RMSD() <= (sLE * m_rmsd_threshold)));

            if (t->KeepMolecule() == false) {
                keep_molecule = false;
                writeStatisticFile(t->Reference(), mol1, t->RMSD());

                if (laststring.compare("") != 0 && laststring.compare(t->Reference()->Name()) != 0)
                    m_first_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\"[style=dotted];\n";

                m_first_content += "\"" + t->Reference()->Name() + "\" [shape=box, label=\"" + t->Reference()->Name() + "\"];\n";
                m_first_content += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\"];\n";
                m_first_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(t->RMSD()) + "];\n";
                laststring = t->Reference()->Name();

#ifdef WriteMoreInfo
                m_dnn_data.push_back(t->getDNNInput());
#endif
                break;
            }
        }

        if (keep_molecule) {
            ConfScanThreadNoReorder* thread = addThreadNoreorder(mol1, rmsd);
            threads.push_back(thread);
            p->addThread(thread);
            AcceptMolecule(mol1);
        } else {
            RejectMolecule(mol1);
        }
        PrintStatus();
        m_all_structures.push_back(mol1);
    }
    p->clear();
    delete p;
}

void ConfScan::PrintSetUp(double dLE, double dLI, double dLH)
{
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
    std::cout << "    Loose Threshold: " << dLI << " MHz" << std::endl;
    std::cout << "    Tight Threshold: " << m_dTI << " MHz" << std::endl;
    std::cout << "    Thresholds in difference of ripser images: " << std::endl;
    std::cout << "    Loose Threshold: " << dLH << " " << std::endl;
    std::cout << "    Tight Threshold: " << m_dTH << " " << std::endl;
    std::cout << "    Thresholds for energy differences in kJ/mol: " << std::endl;
    std::cout << "    Loose Threshold: " << dLE << " " << std::endl;
    std::cout << "    Tight Threshold: " << m_dTE << " " << std::endl;

    std::cout << "" << std::endl
              << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl
              << std::endl;
}

void ConfScan::Reorder(double dLE, double dLI, double dLH, bool reuse_only)
{
    std::string content_dot;
    std::string laststring;
    PrintSetUp(dLE, dLI, dLH);
    /* if ignore rotation or ignore barcode are true, the tight cutoffs are set to -1, that no difference is smaller
     * than the tight threshold and the loose threshold is set to a big number, that every structure has to be reordered*/
    m_skiped = 0;
    m_rejected_directly = 0;
    if (m_ignoreRotation == true) {
        m_dLI = 1e10;
        m_dTI = -1;
    }
    if (m_ignoreBarCode == true) {
        m_dLH = 1e10;
        m_dTH = -1;
    }

    TriggerWriteRestart();

    m_rejected = 0, m_accepted = 0, m_reordered = 0, m_reordered_worked = 0, m_reordered_reused = 0;

    json rmsd = RMSDJson;
    rmsd["silent"] = true;
    rmsd["reorder"] = true;
    rmsd["check"] = CheckConnections();
    rmsd["heavy"] = m_heavy;
    rmsd["method"] = m_RMSDmethod;
    rmsd["element"] = m_rmsd_element_templates;
    rmsd["update-rotation"] = m_update_rotation;
    rmsd["damping"] = m_damping;
    rmsd["nomunkres"] = m_nomunkres;
    rmsd["molalignbin"] = m_molalign;
    rmsd["molaligntol"] = m_molaligntol;

    std::vector<Molecule*> cached;
    if (dLE < 0 && dLI < 0 && dLH < 0)
        cached = m_all_structures;
    else
        cached = m_stored_structures;
    m_maxmol = cached.size();

    m_result = m_previously_accepted;
    m_stored_structures.clear();
    std::vector<ConfScanThread*> threads;
    std::vector<std::vector<int>> rules;
    CxxThreadPool* p = new CxxThreadPool;
    p->setActiveThreadCount(m_threads);

    for (Molecule* mol1 : cached) {
        if (m_result.size() == 0) {
            AcceptMolecule(mol1);
            ConfScanThread* thread = addThread(mol1, rmsd, reuse_only);
            threads.push_back(thread);
            p->addThread(thread);
            m_lowest_energy = mol1->Energy();
            continue;
        }
        p->Reset();
        m_current_energy = mol1->Energy();
        m_dE = (m_current_energy - m_lowest_energy) * 2625.5;

        bool keep_molecule = true;
        bool reorder = false;
        for (int t = 0; t < threads.size(); ++t) {
            if (CheckStop()) {
                fmt::print("\n\n** Found stop file, will end now! **\n\n");
                // TriggerWriteRestart();
                return;
            }
            const Molecule* mol2 = threads[t]->Reference();
            double Ia = abs(mol1->Ia() - mol2->Ia());
            double Ib = abs(mol1->Ib() - mol2->Ib());
            double Ic = abs(mol1->Ic() - mol2->Ic());

            double dI = (Ia + Ib + Ic) * third;
            double dH = (mol1->getPersisentImage() - mol2->getPersisentImage()).cwiseAbs().sum();

            /* rotational = 1
             * ripser     = 2
             * energy     = 4 */
            int looseThresh = 1 * (dI < dLI) + 2 * (dH < dLH) + 4 * (std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 < dLE);
            if ((looseThresh & m_looseThresh) == m_looseThresh || (dLI <= 1e-8 && dLH <= 1e-8 && dLE <= 1e-8)) {
                reorder = true;
                threads[t]->setEnabled(m_openLoop);
                int tightThresh = 1 * (dI < m_dTI) + 2 * (dH < m_dTH) + 4 * ((std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 < m_dTE));

                if ((tightThresh & m_tightThresh) == m_tightThresh) {
                    std::cout << "Differences " << dI << " MHz and " << dH << " below tight threshold, reject molecule directly!" << std::endl;
                    m_lastDI = dI;
                    m_lastDH = dH;
                    writeStatisticFile(mol1, mol2, -1, false);
                    m_threshold.push_back(mol2);
                    m_rejected_directly++;
                    reorder = false;
                    keep_molecule = false;
                    break;
                }
                m_exclude_list.push_back(std::pair<std::string, std::string>(mol1->Name(), mol2->Name()));

            } else {
                threads[t]->setEnabled(m_closeLoop);
            }
        }

        if (reorder && keep_molecule) {
            int free_threads = m_threads;
            if (threads.size())
                free_threads /= threads.size();

            if (free_threads < 1)
                free_threads = 1;
            for (int i = 0; i < threads.size(); ++i) {
                threads[i]->setTarget(mol1);
                threads[i]->setReorderRules(m_reorder_rules);
                threads[i]->setThreads(free_threads);
                for (int j = 0; j < rules.size(); ++j)
                    threads[i]->addReorderRule(rules[j]);
            }

            if (m_RMSDmethod.compare("molalign") != 0 || m_threads == 1) {
                p->StaticPool();
                p->StartAndWait();
            } else {
                for (auto* t : threads) {
                    t->execute();
                    if (t->KeepMolecule() == false)
                        break;
                }
            }

            for (auto* t : threads) {
                if (!t->isEnabled()) {
                    t->setEnabled(true);
                    m_skiped++;
                    continue;
                }
#ifdef WriteMoreInfo
                m_dnn_data.push_back(t->getDNNInput());
#endif
                m_reordered++;
                if (t->KeepMolecule() == false) {
                    keep_molecule = false;
                    m_reordered_worked += t->ReorderWorked();
                    m_reordered_reused += t->ReusedWorked();
                    if (AddRules(t->ReorderRule()))
                        rules.push_back(t->ReorderRule());

                    if (laststring.compare("") != 0 && laststring.compare(t->Reference()->Name()) != 0)
                        m_second_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\"[style=dotted];\n";

                    m_second_content += "\"" + t->Reference()->Name() + "\" [shape=box, label=\"" + t->Reference()->Name() + "\"];\n";
                    m_second_content += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\"];\n";
                    m_second_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(t->RMSD()) + "];\n";
                    laststring = t->Reference()->Name();

                    writeStatisticFile(t->Reference(), mol1, t->RMSD(), true, t->ReorderRule());
                    mol1->ApplyReorderRule(t->ReorderRule());
                    break;
                } else {
                    if ((m_domolalign > 1) && t->RMSD() < m_domolalign * m_rmsd_threshold) {
                        fmt::print(fg(fmt::color::yellow) | fmt::emphasis::bold, "Starting molalign for more precise reordering ...\n");
                        json molalign = rmsd;
                        molalign["method"] = "molalign";
                        m_molalign_count++;
                        RMSDDriver driver(molalign);
                        driver.setReference(t->Reference());
                        driver.setTarget(mol1);
                        driver.start();
                        mol1->LoadMolecule(driver.TargetReorderd());
                        if (driver.RMSD() < m_rmsd_threshold) {
                            m_molalign_success++;
                            keep_molecule = false;
                            m_reordered_worked += 1;
                            m_reordered_reused += 1;
                            writeStatisticFile(t->Reference(), mol1, driver.RMSD(), true, { 0, 0 });

                            if (laststring.compare("") != 0 && laststring.compare(t->Reference()->Name()) != 0)
                                m_second_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\"[style=dotted];\n";

                            m_second_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\";\n";
                            m_second_content += "\"" + mol1->Name() + "\" [shape=box, style=filled,color=\".7 .3 1.0\"];\n";
                            m_second_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(driver.RMSD()) + "];\n";
                            laststring = t->Reference()->Name();

                            fmt::print(fg(fmt::color::yellow) | fmt::emphasis::bold, "... success!\n");
                        } else
                            fmt::print(fg(fmt::color::yellow) | fmt::emphasis::bold, "... without effect!\n");
                    }
                }
            }
        } else
            m_skiped += threads.size();

        if (keep_molecule) {
            AcceptMolecule(mol1);
            ConfScanThread* thread = addThread(mol1, rmsd, reuse_only);
            p->addThread(thread);
            threads.push_back(thread);
        } else {
            RejectMolecule(mol1);
        }
        PrintStatus();
        if (m_result.size() >= m_maxrank)
            break;
        if (m_dE > m_energy_cutoff && m_energy_cutoff != -1) {
            break;
        }
    }
    p->clear();
    delete p;
}

void ConfScan::ReuseOnly()
{
    bool reuse_only = false;
    m_molalign_count = 0;
    m_molalign_success = 0;

    m_maxmol = m_stored_structures.size();
    std::string content_dot;
    std::string laststring;

    /* if ignore rotation or ignore barcode are true, the tight cutoffs are set to -1, that no difference is smaller
     * than the tight threshold and the loose threshold is set to a big number, that every structure has to be reordered*/
    m_skiped = 0;
    m_rejected_directly = 0;
    if (m_ignoreRotation == true) {
        m_dLI = 1e10;
        m_dTI = -1;
    }
    if (m_ignoreBarCode == true) {
        m_dLH = 1e10;
        m_dTH = -1;
    }
    m_dLI *= 2;
    m_dLH *= 2;
    m_dLE *= 2;

    PrintSetUp(m_dLI, m_dLH, m_dLE);

    TriggerWriteRestart();

    m_rejected = 0, m_accepted = 0, m_reordered = 0, m_reordered_worked = 0, m_reordered_reused = 0;

    json rmsd = RMSDJson;
    rmsd["silent"] = true;
    rmsd["reorder"] = true;
    rmsd["check"] = CheckConnections();
    rmsd["heavy"] = m_heavy;
    rmsd["method"] = m_RMSDmethod;
    rmsd["element"] = m_rmsd_element_templates;
    rmsd["update-rotation"] = m_update_rotation;
    rmsd["damping"] = m_damping;
    rmsd["nomunkres"] = m_nomunkres;
    rmsd["molalignbin"] = m_molalign;
    std::vector<Molecule*> cached = m_stored_structures;
    m_result = m_previously_accepted;
    m_stored_structures.clear();
    std::vector<ConfScanThread*> threads;
    std::vector<std::vector<int>> rules;
    CxxThreadPool* p = new CxxThreadPool;
    p->setActiveThreadCount(m_threads);

    for (Molecule* mol1 : cached) {
        if (m_result.size() == 0) {
            AcceptMolecule(mol1);
            ConfScanThread* thread = addThread(mol1, rmsd, reuse_only);
            threads.push_back(thread);
            p->addThread(thread);
            m_lowest_energy = mol1->Energy();
            continue;
        }
        p->Reset();
        m_current_energy = mol1->Energy();
        m_dE = (m_current_energy - m_lowest_energy) * 2625.5;

        bool keep_molecule = true;
        bool reorder = false;
        for (int t = 0; t < threads.size(); ++t) {
            if (CheckStop()) {
                fmt::print("\n\n** Found stop file, will end now! **\n\n");
                // TriggerWriteRestart();
                return;
            }
            const Molecule* mol2 = threads[t]->Reference();
            std::pair<std::string, std::string> names(mol1->Name(), mol2->Name());
            if (std::find(m_exclude_list.begin(), m_exclude_list.end(), names) != m_exclude_list.end()) {
                m_duplicated++;
                continue;
            }
            double Ia = abs(mol1->Ia() - mol2->Ia());
            double Ib = abs(mol1->Ib() - mol2->Ib());
            double Ic = abs(mol1->Ic() - mol2->Ic());

            double dI = (Ia + Ib + Ic) * third;

            double dH = (mol1->getPersisentImage() - mol2->getPersisentImage()).cwiseAbs().sum();

            /* rotational = 1
             * ripser     = 2
             * energy     = 4 */
            int looseThresh = 1 * (dI < m_dLI) + 2 * (dH < m_dLH) + 4 * (std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 < m_dLE);
            if ((looseThresh & m_looseThresh) == m_looseThresh /* || (m_dLI == 0 && m_dLH == 0 && m_dLE == 0)*/) {
                reorder = true;
                threads[t]->setEnabled(m_openLoop);
                int tightThresh = 1 * (dI < m_dTI) + 2 * (dH < m_dTH) + 4 * ((std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 < m_dTE));

                if ((tightThresh & m_tightThresh) == m_tightThresh) {
                    std::cout << "Differences " << dI << " MHz and " << dH << " below tight threshold, reject molecule directly!" << std::endl;
                    m_lastDI = dI;
                    m_lastDH = dH;
                    writeStatisticFile(mol1, mol2, -1, false);
                    m_threshold.push_back(mol2);
                    m_rejected_directly++;
                    reorder = false;
                    keep_molecule = false;
                    break;
                }
            } else {
                threads[t]->setEnabled(m_closeLoop);
            }
        }

        if (reorder && keep_molecule) {
            int free_threads = m_threads;
            if (threads.size())
                free_threads /= threads.size();

            if (free_threads < 1)
                free_threads = 1;
            for (int i = 0; i < threads.size(); ++i) {
                threads[i]->setTarget(mol1);
                threads[i]->setReorderRules(m_reorder_rules);
                threads[i]->setThreads(free_threads);
                for (int j = 0; j < rules.size(); ++j)
                    threads[i]->addReorderRule(rules[j]);
            }

            if (m_RMSDmethod.compare("molalign") != 0 || m_threads == 1) {
                p->StaticPool();
                p->StartAndWait();
            } else {
                for (auto* t : threads) {
                    t->execute();
                    if (t->KeepMolecule() == false)
                        break;
                }
            }

            for (auto* t : threads) {
                if (!t->isEnabled()) {
                    t->setEnabled(true);
                    m_skiped++;
                    continue;
                }

                m_reordered++;
                if (t->KeepMolecule() == false) {
                    keep_molecule = false;
                    m_reordered_worked += t->ReorderWorked();
                    m_reordered_reused += t->ReusedWorked();
                    if (AddRules(t->ReorderRule()))
                        rules.push_back(t->ReorderRule());

                    writeStatisticFile(t->Reference(), mol1, t->RMSD(), true, t->ReorderRule());
                    mol1->ApplyReorderRule(t->ReorderRule());

                    if (laststring.compare("") != 0 && laststring.compare(t->Reference()->Name()) != 0)
                        m_third_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\"[style=dotted];\n";

                    m_third_content += "\"" + t->Reference()->Name() + "\" [shape=box, label=\"" + t->Reference()->Name() + "\"];\n";
                    m_third_content += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\"];\n";
                    m_third_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(t->RMSD()) + "];\n";

                    laststring = t->Reference()->Name();
                    break;
                } else {
                    if ((m_domolalign > 1) && t->RMSD() < m_domolalign * m_rmsd_threshold) {
                        fmt::print(fg(fmt::color::yellow) | fmt::emphasis::bold, "Starting molalign for more precise reordering ...\n");
                        json molalign = rmsd;
                        molalign["method"] = "molalign";
                        m_molalign_count++;
                        RMSDDriver driver(molalign);
                        driver.setReference(t->Reference());
                        driver.setTarget(mol1);
                        driver.start();
                        mol1->LoadMolecule(driver.TargetReorderd());
                        if (driver.RMSD() < m_rmsd_threshold) {
                            m_molalign_success++;
                            keep_molecule = false;

                            if (laststring.compare("") != 0 && laststring.compare(t->Reference()->Name()) != 0)
                                m_third_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\"[style=dotted];\n";

                            m_third_content += "\"" + mol1->Name() + "\" [shape=box, style=filled,color=\".7 .3 1.0\"];\n";
                            m_third_content += "\"" + t->Reference()->Name() + "\" [shape=box];\n";
                            m_third_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(driver.RMSD()) + "];\n";

                            laststring = t->Reference()->Name();
                            m_reordered_worked += 1;
                            m_reordered_reused += 1;
                            writeStatisticFile(t->Reference(), mol1, driver.RMSD(), true, { 0, 0 });
                            fmt::print(fg(fmt::color::yellow) | fmt::emphasis::bold, "... success!\n");
                        } else
                            fmt::print(fg(fmt::color::yellow) | fmt::emphasis::bold, "... without effect!\n");
                    }
                }
            }
        } else
            m_skiped += threads.size();

        if (keep_molecule) {
            AcceptMolecule(mol1);

            ConfScanThread* thread = addThread(mol1, rmsd, reuse_only);
            p->addThread(thread);
            threads.push_back(thread);
        } else {
            RejectMolecule(mol1);
        }

        PrintStatus();
        if (m_result.size() >= 2 * m_maxrank)
            break;
        if (m_dE > m_energy_cutoff && m_energy_cutoff != -1) {
            break;
        }
    }
    p->clear();
    delete p;
}
/*
void ConfScan::ReuseOnly()
{
    m_stored_structures.clear();
    m_maxmol = m_final_batch.size();
    m_rejected = 0;
    m_reordered_reused = 0;
    m_reordered_worked = 0;
    m_skiped = 0;
    std::string content_dot;
    std::string laststring;
    m_maxmol = m_ordered_list.size();

    json rmsd = RMSDJson;
    rmsd["silent"] = true;
    rmsd["check"] = CheckConnections();
    rmsd["heavy"] = m_heavy;
    rmsd["noreorder"] = true;

    std::vector<ConfScanThreadNoReorder*> threads;
    std::vector<std::vector<int>> rules;
    CxxThreadPool* p = new CxxThreadPool;
    p->setActiveThreadCount(m_threads);

    for (auto& i : m_final_batch) {
        if (m_skip) {
            m_skip--;
            continue;
        }
        if (m_maxrank <= m_accepted && m_maxrank > -1)
            continue;

        Molecule* mol1 = i.second;
        if (mol1->Check() == 1) {
            m_rejected++;
            m_start++;
            PrintStatus();
            continue;
        }
        if (m_result.size() == 0) {
            AcceptMolecule(mol1);
            ConfScanThreadNoReorder* thread = addThreadNoreorder(mol1, rmsd);
            threads.push_back(thread);
            p->addThread(thread);

            m_lowest_energy = mol1->Energy();
            continue;
        }
        p->Reset();
        m_current_energy = mol1->Energy();
        m_dE = (m_current_energy - m_lowest_energy) * 2625.5;

        bool keep_molecule = true;
        for (int i = 0; i < threads.size(); ++i) {
            threads[i]->setTarget(mol1);
        }
        p->StaticPool();
        p->StartAndWait();

        for (auto* t : threads) {

            if (t->KeepMolecule() == false) {
                keep_molecule = false;
                writeStatisticFile(t->Reference(), mol1, t->RMSD());

                if (laststring.compare("") != 0 && laststring.compare(t->Reference()->Name()) != 0)
                    m_4th_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\"[style=dotted];\n";

                m_4th_content += "\"" + t->Reference()->Name() + "\" [shape=box, label=\"" + t->Reference()->Name() + "\"];\n";
                m_4th_content += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\"];\n";
                m_4th_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(t->RMSD()) + "];\n";
                laststring = t->Reference()->Name();

                break;
            }
        }

        if (keep_molecule) {
            ConfScanThreadNoReorder* thread = addThreadNoreorder(mol1, rmsd);
            threads.push_back(thread);
            p->addThread(thread);
            AcceptMolecule(mol1);
        } else {
            RejectMolecule(mol1);
        }
        PrintStatus();
    }
    p->clear();
    delete p;
}
*/
ConfScanThreadNoReorder* ConfScan::addThreadNoreorder(const Molecule* reference, const json& config)
{
    ConfScanThreadNoReorder* thread = new ConfScanThreadNoReorder(m_rmsd_threshold, m_MaxHTopoDiff, config);
    thread->setReference(*reference);
    return thread;
}

ConfScanThread* ConfScan::addThread(const Molecule* reference, const json& config, bool reuse_only)
{
    ConfScanThread* thread = new ConfScanThread(m_reorder_rules, m_rmsd_threshold, m_MaxHTopoDiff, reuse_only, config);
    thread->setReference(*reference);
    return thread;
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
    if (m_writeFiles && !m_reduced_file) {
        for (const auto molecule : m_rejected_structures) {
            molecule->appendXYZFile(m_rejected_filename);
        }

        for (const auto molecule : m_threshold)
            molecule->appendXYZFile(m_threshold_filename);
    }
    std::cout << m_stored_structures.size() << " structures were kept - of " << m_molecules.size() - m_fail << " total!" << std::endl;
}

bool ConfScan::AddRules(const std::vector<int>& rules)
{
    if (rules.size() == 0)
        return false;

    if (std::find(m_reorder_rules.begin(), m_reorder_rules.end(), rules) == m_reorder_rules.end()) {
        m_reorder_rules.push_back(rules);
    }
    return true;
}

void ConfScan::PrintStatus(const std::string& info)
{
    std::cout << std::endl;
    std::cout << "             ###   " << std::setprecision(4) << (m_stored_structures.size() + m_rejected) / double(m_maxmol) * 100 << "% done!   ###" << std::endl;
    if (info.compare("") != 0)
        std::cout << info;
    std::cout << "# Accepted : " << m_stored_structures.size() << "     ";
    std::cout << "# Rejected : " << m_rejected << "     ";
    std::cout << "# Reordered : " << m_reordered << " (+ " << m_molalign_count << ")"
              << "     ";
    std::cout << "# Successfully : " << m_reordered_worked << " (+ " << m_molalign_success << ")"
              << "    ";
    std::cout << "# Reused Results : " << m_reordered_reused << "     ";
    std::cout << "# Reordering Skipped : " << m_skiped << " (+ " << m_duplicated << ")";
    std::cout << "# Rejected Directly : " << m_rejected_directly << "     ";

    std::cout << "# Current Energy [kJ/mol] : " << m_dE << std::endl;
}

void ConfScan::writeStatisticFile(const Molecule* mol1, const Molecule* mol2, double rmsd, bool reason, const std::vector<int>& rule)
{
    if (!(m_writeFiles && !m_reduced_file))
        return;
    std::ofstream result_file;
    result_file.open(m_statistic_filename, std::ios_base::app);
    if (reason)
        result_file << "Molecule got rejected due to small rmsd " << rmsd << " with and energy difference of " << std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 << " kJ/mol." << std::endl;
    else
        result_file << "Molecule got rejected as differences " << m_lastDI << " MHz and " << m_lastDH << " are below the estimated thresholds;  with and energy difference of " << std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 << " kJ/mol." << std::endl;
    if (rule.size())
        for (auto i : rule)
            result_file << i << "|";
    result_file << std::endl;
    result_file << mol1->XYZString();
    result_file << mol2->XYZString();
    result_file << std::endl;
    result_file.close();

    if (m_write && rule.size()) {
        mol1->writeXYZFile("A" + std::to_string(m_rejected) + ".xyz");
        mol2->writeXYZFile("B" + std::to_string(m_rejected) + ".xyz");
    }
}
