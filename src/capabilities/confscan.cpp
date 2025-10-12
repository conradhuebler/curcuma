/*
 * <Scan and judge conformers from different input. >
 * Copyright (C) 2020 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include <fmt/color.h>
#include <fmt/core.h>

#include "src/core/global.h" // For CurcumaLogger - Claude Generated
#include "src/global_config.h" // For new logging system - Claude Generated
#include "src/core/parameter_registry.h" // For ParameterRegistry - Claude Generated 2025

#include "src/capabilities/confstat.h"
#include "src/capabilities/persistentdiagram.h"
#include "src/capabilities/rmsd.h"

#include "src/core/fileiterator.h"

#include "src/core/energycalculator.h"

#include "src/tools/general.h"

#include "json.hpp"
using json = nlohmann::json;

#include "confscan.h"

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 2)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

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

    m_input.dIa = Ia;
    m_input.dIb = Ib;
    m_input.dIc = Ic;
    const Matrix m = (m_reference.getPersisentImage() - m_target.getPersisentImage());
    double diff_ripser = m.cwiseAbs().sum();
    m_input.dH = diff_ripser;
    m_input.dHM = m;
    m_input.dE = std::abs(m_reference.Energy() - m_target.Energy()) * 2625.5;

    m_old_rmsd = m_driver->BestFitRMSD();
    if (m_old_rmsd < m_rmsd_threshold) {
        m_rmsd = m_old_rmsd;
        m_keep_molecule = false;
        m_break_pool = true;
        return 0;
    }

    for (int i = 0; i < m_reorder_rules.size(); ++i) {
        if (m_reorder_rules[i].size() != m_reference.AtomCount() || m_reorder_rules[i].size() == 0)
            continue;

        double tmp_rmsd = m_driver->Rules2RMSD(m_reorder_rules[i]);
        if (tmp_rmsd < m_rmsd_threshold && (m_MaxHTopoDiff == -1 || m_driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
            m_keep_molecule = false;
            /* early breaks affect the number of finally accepted structures */
            m_break_pool = (m_earlybreak & 1) == 0;
            m_reused_worked = true;
            m_rmsd = tmp_rmsd;

            m_input.rmsd = m_rmsd;
            m_reorder_rule = m_reorder_rules[i];
            if (m_verbosity >= 3) {
                if (m_earlybreak)
                    CurcumaLogger::info_fmt("Reuse: {} {} {:.6f} Early break", m_reference.Name(), m_target.Name(), m_rmsd);
                else
                    CurcumaLogger::info_fmt("Reuse: {} {} {:.6f}", m_reference.Name(), m_target.Name(), m_rmsd);
            }
            m_driver->clear();
            return 0;
        }
    }

    if (m_reuse_only) {
        m_driver->clear();
        return 0;
    }
    /* In case multiple threads are working, the target molecule gets reset / wrongly overwritten or something
        similar. This is a workaround to prevent this.
        So using a different number of threads effects the number of finally accepted structures.
    */

    m_driver->start();
    m_rmsd = m_driver->RMSD();

    m_input.rmsd = m_rmsd;

    if (m_rmsd <= m_rmsd_threshold && (m_MaxHTopoDiff == -1 || m_driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
        m_keep_molecule = false;
        /* early breaks affect the number of finally accepted structures */

        m_break_pool = (m_earlybreak & 2) == 0;
        m_reorder_worked = true;

        m_reorder_rule = m_driver->ReorderRules();
    }
    if (m_verbosity >= 3) {
        if (m_break_pool)
            CurcumaLogger::info_fmt("Permutation: {} {} {:.6f} Early break", m_reference.Name(), m_target.Name(), m_rmsd);
        else
            CurcumaLogger::info_fmt("Permutation: {} {} {:.6f}", m_reference.Name(), m_target.Name(), m_rmsd);
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
    m_input.rmsd = m_rmsd;

    double Ia = abs(m_reference.Ia() - m_target.Ia());
    double Ib = abs(m_reference.Ib() - m_target.Ib());
    double Ic = abs(m_reference.Ic() - m_target.Ic());
    m_input.dIa = Ia;
    m_input.dIb = Ib;
    m_input.dIc = Ic;
    m_DI = (Ia + Ib + Ic) * third;
    const Matrix m = (m_reference.getPersisentImage() - m_target.getPersisentImage());
    m_DH = m.cwiseAbs().sum();
    m_input.dH = m_DH;
    m_input.dHM = m;
    m_input.dE = std::abs(m_reference.Energy() - m_target.Energy()) * 2625.5;
    if (m_rmsd <= m_rmsd_threshold && (m_MaxHTopoDiff == -1 || m_driver->HBondTopoDifference() <= m_MaxHTopoDiff)) {
        m_keep_molecule = false;
        m_break_pool = true;
    }

    m_driver->clear();
    return 0;
}

ConfScan::ConfScan(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("confscan"), controller, silent),
      m_config(std::vector<std::string>{"confscan", "rmsd"}, controller)  // Claude Generated 2025: Multi-Module ConfigManager
{
    LoadControlJson();
}

ConfScan::~ConfScan()
{
    for (auto i : m_molecules) {
        delete i.second;
    }
}

void ConfScan::LoadControlJson()
{
    // Claude Generated 2025: Migrated to ConfigManager - Multi-Module architecture

    // ConfScan-specific parameters
    m_noname = m_config.get<bool>("noname");
    m_restart = m_config.get<bool>("restart");

    // RMSD parameters from rmsd module (dot notation)
    // Note: RMSD has "protons" (include protons), ConfScan expects "heavy" (exclude protons) - inverse logic!
    bool protons_included = m_config.get<bool>("rmsd.protons");
    m_heavy = !protons_included;  // Invert: protons=false means heavy-only

    // Filtering & thresholds
    m_rmsd_threshold = m_config.get<double>("rmsd");
    if (m_config.get<bool>("get_rmsd")) {
        m_rmsd_threshold = -1;
        m_rmsd_set = false;
        CurcumaLogger::warn("RMSD value is not set, will obtain it from ensemble");
    }
    if (m_rmsd_threshold < 0) {
        m_rmsd_set = false;
        m_rmsd_threshold = 1e5;
        CurcumaLogger::warn("RMSD value is not set, will obtain it from ensemble");
    }
    m_maxrank = m_config.get<double>("rank");
    m_writeXYZ = m_config.get<bool>("write_xyz");
    m_force_reorder = m_config.get<bool>("force_reorder");
    m_check_connections = m_config.get<bool>("check_connections");
    m_energy_cutoff = m_config.get<double>("max_energy");

    // Multi-type threshold parameters (simplified with ConfigManager) - Claude Generated 2025
    // sLX can be "default", "1.5", or "1.0,2.0,3.0"
    std::string slx = m_config.get<std::string>("slx");
    if (slx == "default") {
        if (m_verbosity >= 2)
            CurcumaLogger::info("Using default values for the steps");
        m_sLE = { 1.0, 2.0 };
        m_sLI = { 1.0, 2.0 };
        m_sLH = { 1.0, 2.0 };
    } else if (slx.find(',') != std::string::npos) {
        // CSV: "1.0,2.0,3.0"
        m_sLE = Tools::String2DoubleVec(slx, ",");
        m_sLI = m_sLE;
        m_sLH = m_sLE;
    } else {
        // Single number: "1.5"
        double val = std::stod(slx);
        m_sLE = { val };
        m_sLI = { val };
        m_sLH = { val };
    }

    // Individual overrides (sLE, sLI, sLH can override sLX)
    std::string sle = m_config.get<std::string>("sle");
    if (sle != "default") {
        if (sle.find(',') != std::string::npos) {
            m_sLE = Tools::String2DoubleVec(sle, ",");
        } else {
            m_sLE = { std::stod(sle) };
        }
    }

    std::string sli = m_config.get<std::string>("sli");
    if (sli != "default") {
        if (sli.find(',') != std::string::npos) {
            m_sLI = Tools::String2DoubleVec(sli, ",");
        } else {
            m_sLI = { std::stod(sli) };
        }
    }

    std::string slh = m_config.get<std::string>("slh");
    if (slh != "default") {
        if (slh.find(',') != std::string::npos) {
            m_sLH = Tools::String2DoubleVec(slh, ",");
        } else {
            m_sLH = { std::stod(slh) };
        }
    }
    if (m_sLE.size() != m_sLI.size() || m_sLE.size() != m_sLH.size()) {
        CurcumaLogger::error("Inconsistent length of steps requested, will abort now");
        exit(1);
    }

    if (m_verbosity >= 2) {
        CurcumaLogger::info("Using the following steps for the thresholds:");
        for (std::size_t i = 0; i < m_sLE.size(); ++i) {
            CurcumaLogger::info_fmt("sLE: {}, sLI: {}, sLH: {}",
                to_string_with_precision(m_sLE[i]),
                to_string_with_precision(m_sLI[i]),
                to_string_with_precision(m_sLH[i]));
        }
    }

    // Tight thresholds
    m_sTE = m_config.get<double>("ste");
    m_sTI = m_config.get<double>("sti");
    m_sTH = m_config.get<double>("sth");

    // Workflow control
    m_reset = m_config.get<bool>("reset");
    m_analyse = m_config.get<bool>("analyse");
    m_skipinit = m_config.get<bool>("skip_init");
    m_skipreorder = m_config.get<bool>("skip_reorder");
    m_skipreuse = m_config.get<bool>("skip_reuse");
    m_mapped = m_config.get<bool>("mapped");
    m_skip_orders = m_config.get<bool>("skip_orders");

    // Algorithm parameters
    m_lastdE = m_config.get<double>("last_de");
    m_getrmsd_scale = m_config.get<double>("getrmsd_scale");
    m_getrmsd_thresh = m_config.get<double>("getrmsd_thresh");
    m_skip = m_config.get<int>("skip");
    m_cycles = m_config.get<int>("cycles");
    m_maxParam = m_config.get<int>("max_param");
    m_useorders = m_config.get<int>("use_orders");

    // Output control
    m_allxyz = m_config.get<bool>("all_xyz");
    m_reduced_file = m_config.get<bool>("fewer_file");
    m_update = m_config.get<bool>("update");
    m_split = m_config.get<bool>("split");
    m_write = m_config.get<bool>("write_files");

    // Performance
    m_threads = m_config.get<int>("threads");  // ConfScan ensemble threads

    // RMSD parameters from rmsd module (dot notation) - corrected names
    m_RMSDmethod = m_config.get<std::string>("rmsd.method");
    m_update_rotation = m_config.get<bool>("rmsd.update_rotation");
    m_nomunkres = m_config.get<bool>("rmsd.nomunkres", false);  // May not exist in RMSD
    m_molalign = m_config.get<std::string>("rmsd.molalign_bin");
    m_molaligntol = m_config.get<int>("rmsd.molalign_tolerance");
    m_domolalign = m_config.get<double>("rmsd.do_molalign", -1.0);  // May not exist
    m_ignoreRotation = m_config.get<bool>("rmsd.ignore_rotation", false);  // May not exist

    // Ripser/Topology parameters
    m_MaxHTopoDiff = m_config.get<int>("ripser.max_topo_diff", -1);  // Optional, may not be in rmsd
    m_ignoreBarCode = m_config.get<bool>("ripser.ignore_barcode", false);  // Optional

    // Threshold bit flags
    m_looseThresh = m_config.get<int>("loose_thresh");
    m_tightThresh = m_config.get<int>("tight_thresh");
    m_earlybreak = m_config.get<int>("early_break");

    // Early break flags logging
    if ((m_earlybreak & 1) == 0)
        if (m_verbosity >= 2)
            CurcumaLogger::info("Early break in reuse part is enabled");
    if ((m_earlybreak & 2) == 0)
        if (m_verbosity >= 2)
            CurcumaLogger::info("Early break in reorder part is enabled");

    // RMSD Element handling (simplified with ConfigManager - Claude Generated 2025)
    // No more try-catch needed! #pragma message removed!
    m_rmsd_element_templates = m_config.get<std::string>("rmsd.element", "7");  // Default to nitrogen
    if (!m_rmsd_element_templates.empty() && m_rmsd_element_templates != "0") {
        if (m_rmsd_element_templates.find(',') != std::string::npos) {
            // CSV: "7,8"
            StringList elements = Tools::SplitString(m_rmsd_element_templates, ",");
            for (const std::string& str : elements)
                m_element_templates.push_back(std::stod(str));

            if (m_element_templates.size())
                m_RMSDElement = m_element_templates[0];
        } else {
            // Single number: "7"
            m_RMSDElement = std::stoi(m_rmsd_element_templates);
            m_element_templates.push_back(m_RMSDElement);
        }
    }

    // Hybrid method validation
    if (m_RMSDmethod == "hybrid" && m_rmsd_element_templates.empty()) {
        CurcumaLogger::warn("Reordering method hybrid has to be combined with element types. I will choose for you nitrogen and oxygen!");
        CurcumaLogger::info("This is equivalent to adding: '-rmsd.element 7,8' to your argument list");
        m_rmsd_element_templates = "7,8";
        m_element_templates = {7.0, 8.0};
        m_RMSDElement = 7;
    }

    // Molalign citation
    if (m_RMSDmethod == "molalign") {
        if (m_verbosity >= 1) {
            CurcumaLogger::citation("J. Chem. Inf. Model. 2023, 63, 4, 1157–1165 - DOI: 10.1021/acs.jcim.2c01187");
        }
        m_domolalign = -1;
    }

    // Previous accepted structures
    m_prev_accepted = m_config.get<std::string>("accepted");

    // RMSD method logging
    if (m_verbosity >= 2) {
        CurcumaLogger::info_fmt("Permutation of atomic indices performed according to {}", m_RMSDmethod);
    }

    if (m_useorders == -1)
        m_useorders = 10;

    if (m_verbosity >= 2) {
        CurcumaLogger::info("Current Configuration:");
        CurcumaLogger::param("Threads", m_threads);
        CurcumaLogger::param("Molalign Tolerance", m_molaligntol);
        CurcumaLogger::param("Force Reorder", m_force_reorder);
        CurcumaLogger::param("Verbosity", m_verbosity);
        CurcumaLogger::param("Write", m_write);
        CurcumaLogger::param("Update Rotation", m_update_rotation);
        CurcumaLogger::param("Split", m_split);
        CurcumaLogger::param("Method", m_method);
    }
}

bool ConfScan::openFile()
{
    bool xyzfile = std::string(Filename()).find(".xyz") != std::string::npos || std::string(Filename()).find(".trj") != std::string::npos;
    if (xyzfile == false)
        throw 1;

    int molecule = 0;
    // Claude Generated 2025: Export flattened config for legacy PersistentDiagram API
    PersistentDiagram diagram(m_config.exportConfig());
    FileIterator file(Filename());
    int calcH = 0;
    int calcI = 0;
    // std::cout << m_looseThresh <<" "<<int((m_looseThresh & 1) == 1) << " " << int((m_looseThresh & 2) == 2) << std::endl;

    if (m_verbosity >= 2) {
        CurcumaLogger::info("Calculation of descriptors:");
        if ((m_looseThresh & 1) == 1)
            CurcumaLogger::info("  - rotational constants");
        if ((m_looseThresh & 2) == 2)
            CurcumaLogger::info("  - ripser barcodes");
        CurcumaLogger::info("required");
    }
    while (!file.AtEnd()) {
        Molecule* mol = new Molecule(file.Next());
        double energy = mol->Energy();
        if (std::abs(energy) < 1e-5 || m_method.compare("") != 0) {
            // XTBInterface interface; // As long as xtb leaks, we have to put it heare
            if (m_method == "")
                m_method = "gfn2";
            // Claude Generated: Use new constructor with basename for parameter caching
            EnergyCalculator interface(m_method, m_controller, Basename());
            // I might not leak really, but was unable to clear everything

            interface.setMolecule(mol->getMolInfo());
            energy = interface.CalculateEnergy(false);
        }
        m_ordered_list.insert(std::pair<double, int>(energy, molecule));
        molecule++;
        if (m_noname)
            mol->setName(NamePattern(molecule));

        auto rot = std::chrono::system_clock::now();
        if ((m_looseThresh & 1) == 1)
            mol->CalculateRotationalConstants();
        auto ripser = std::chrono::system_clock::now();
        if ((m_looseThresh & 2) == 2) {
            diagram.setDistanceMatrix(mol->LowerDistanceVector());
            mol->setPersisentImage(diagram.generateImage(diagram.generatePairs()));
        }
        calcH += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - ripser).count();
        calcI += std::chrono::duration_cast<std::chrono::milliseconds>(ripser - rot).count();

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
                // Claude Generated: Use new constructor with basename for parameter caching
                EnergyCalculator interface(m_method, m_controller, Basename());
                // I might not leak really, but was unable to clear everything

                interface.setMolecule(mol->getMolInfo());
                energy = interface.CalculateEnergy(false);
            }
            min_energy = std::min(min_energy, energy);
            auto rot = std::chrono::system_clock::now();
            if ((m_looseThresh & 1) == 1)
                mol->CalculateRotationalConstants();
            auto ripser = std::chrono::system_clock::now();

            // diagram.setDimension(2);
            if ((m_looseThresh & 2) == 2) {
                diagram.setDistanceMatrix(mol->LowerDistanceVector());
                mol->setPersisentImage(diagram.generateImage(diagram.generatePairs()));
            }
            calcH += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - ripser).count();
            calcI += std::chrono::duration_cast<std::chrono::milliseconds>(ripser - rot).count();

            m_previously_accepted.push_back(mol);
        }
        m_lowest_energy = min_energy;
        m_result = m_previously_accepted;
    }
    m_timing_rot = calcI;
    m_timing_ripser = calcH;

    if (m_verbosity >= 2) {
        CurcumaLogger::info("Time for calculating descriptors:");
        CurcumaLogger::info_fmt("Rotational constants: {:.3f} seconds", m_timing_rot / 1000.0);
        CurcumaLogger::info_fmt("Ripser bar code: {:.3f} seconds", m_timing_ripser / 1000.0);
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

    int error = 0;
    for (const auto& f : files) {
        std::vector<std::vector<int>> reorder_cached;

#ifdef CURCUMA_DEBUG
        if (m_verbosity >= 3)
            CurcumaLogger::debug(1, fmt::format("Reading file {}", f));
#endif
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
        if (m_restart) {
            for (const auto& vector : reorder_cached)
                if (std::find(m_reorder_rules.begin(), m_reorder_rules.end(), vector) == m_reorder_rules.end())
                    m_reorder_rules.push_back(vector);
        }
    }
    m_useRestart = files.size() == 1 && error != int(files.size());
    if (m_verbosity >= 2)
        CurcumaLogger::info_fmt("Starting with {} initial reorder rules", m_reorder_rules.size());
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

    m_result_basename = Filename();
    m_result_basename.erase(m_result_basename.end() - 4, m_result_basename.end());

    m_accepted_filename = m_result_basename + ".accepted.xyz";
    m_1st_filename = m_result_basename + ".initial.xyz";
    m_2nd_filename = m_result_basename + ".reorder";
    m_3rd_filename = m_result_basename + ".reuse.xyz";
    m_rejected_filename = m_result_basename + ".rejected.xyz";
    m_statistic_filename = m_result_basename + ".statistic.log";
    m_joined_filename = m_result_basename + ".joined.xyz";
    m_threshold_filename = m_result_basename + ".thresh.xyz";
    m_param_file = m_result_basename + ".param.dat";
    m_skip_file = m_result_basename + ".param.skip.dat";
    m_perform_file = m_result_basename + ".param.perf.dat";
    m_success_file = m_result_basename + ".param.success.dat";
    m_limit_file = m_result_basename + ".param.limit.dat";

    auto createFile = [](const std::string& filename) {
        std::ofstream file(filename);
        file.close();
    };

    if (m_writeFiles) {
        createFile(m_accepted_filename);
        if (!m_reduced_file) {
            createFile(m_rejected_filename);
            createFile(m_statistic_filename);
            createFile(m_threshold_filename);
            createFile(m_1st_filename);
        }
    }

    if (!m_previously_accepted.empty()) {
        createFile(m_joined_filename);
    }

    if (m_analyse) {
        createFile(m_success_file);
        std::ofstream parameters_file(m_success_file);
        parameters_file << "# RMSD(new)\tRMSD(old)\tDelta E\tDelta H\tDelta I" << std::endl;
        parameters_file.close();

        createFile(m_skip_file);
        createFile(m_perform_file);
        createFile(m_limit_file);
    }

    if (m_verbosity >= 1) {
        CurcumaLogger::info(std::string(70, '='));
        CurcumaLogger::info_fmt("RMSD Calculation will be performed on {}", m_heavy ? "heavy atoms" : "all atoms");
        CurcumaLogger::info_fmt("RMSD Threshold set to: {:.6f} Angstrom", m_rmsd_threshold);
        CurcumaLogger::info_fmt("Highest energy conformer allowed: {:.2f} kJ/mol", m_energy_cutoff);
        CurcumaLogger::info("Threshold multipliers are loose / tight");
    }

    auto printThresholds = [](const std::string& label, const std::vector<double>& loose, double tight) {
        // Claude Generated: Convert to new logging system
        std::ostringstream loose_values;
        for (const auto& d : loose)
            loose_values << d << " ";
        CurcumaLogger::info_fmt("    {} definition for loose {} and tight thresholds {}",
            label, loose_values.str(), tight);
    };

    printThresholds("Ripser Persistence Diagrams", m_sLH, m_sTH);
    printThresholds("Rotational Constants", m_sLI, m_sTI);
    printThresholds("Energy", m_sLE, m_sTE);

    // Claude Generated: Convert end separator to new logging system
    if (m_verbosity >= 1) {
        CurcumaLogger::info("");
        CurcumaLogger::info("''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''");
        CurcumaLogger::info("");
    }
}

void ConfScan::AcceptMolecule(Molecule* molecule)
{
    m_result.push_back(molecule);
    m_stored_structures.push_back(molecule);
    m_accepted++;
    if (m_writeFiles && !m_reduced_file && m_current_filename.length()) {
        molecule->appendXYZFile(m_current_filename);
    }
    // Claude Generated: Always show accept info for each structure at default verbosity
    if (m_verbosity >= 1) { // Show at SMALL_PRINT level and above
        // Ensure logger verbosity matches local verbosity for accept messages
        int old_verbosity = CurcumaLogger::get_verbosity();
        CurcumaLogger::set_verbosity(m_verbosity);
        CurcumaLogger::success_fmt("Accept {}", molecule->Name());
        CurcumaLogger::set_verbosity(old_verbosity);
    }
}

void ConfScan::RejectMolecule(Molecule* molecule)
{
    m_rejected_structures.push_back(molecule);
    m_rejected++;
    // Claude Generated: Always show reject info for each structure at default verbosity
    if (m_verbosity >= 1) { // Show at SMALL_PRINT level and above
        // Ensure logger verbosity matches local verbosity for reject messages
        int old_verbosity = CurcumaLogger::get_verbosity();
        CurcumaLogger::set_verbosity(m_verbosity);
        CurcumaLogger::warn_fmt("Reject {}", molecule->Name());
        CurcumaLogger::set_verbosity(old_verbosity);
    }
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
    // Level 1: Display input parameters nicely formatted - Claude Generated
    if (m_verbosity >= 1) {
        CurcumaLogger::header("Conformational Scanning");

        // Display parameter comparison table showing defaults vs current settings at level 1
        // Claude Generated 2025: Use ConfigManager export for comparison table
        CurcumaLogger::param_comparison_table(m_config.exportConfig(), m_controller, "ConfScan Configuration");
    }

    SetUp();
    RunTimer timer(false);
    std::ofstream result_file;

    if (!m_skipinit) {
        if (m_verbosity >= 1) {
            CurcumaLogger::success("Initial Pass: Performing RMSD calculation without reordering");
        }
        m_current_filename = m_1st_filename;
        if (m_writeFiles && !m_reduced_file) {
            result_file.open(m_statistic_filename, std::ios_base::app);
            result_file << "Results of 1st Pass" << std::endl;
            result_file.close();
        }
        CheckOnly(m_sLE[0], m_sLI[0], m_sLH[0]);
        PrintStatus("Result initial pass:");
        if (m_analyse) {
            WriteDotFile(m_result_basename + ".initial.dot", m_first_content);
        }
        // m_sLE[0] = 1;
        // m_sLI[0] = 1;
        // m_sLH[0] = 1;
        if (m_analyse) {
            m_collective_content = "edge [color=green];\n";
            for (auto node : m_nodes)
                m_collective_content += node.second;
            m_nodes.clear();
            m_collective_content += m_first_content + "\n";

            std::ofstream parameters_file;
            parameters_file.open(m_param_file);
            parameters_file << "# RMSD(old)\tDelta E\tDelta H\tDelta I" << std::endl;
            bool breakline = false;
            for (const auto& i : m_listThresh) {
                if (i.first > m_rmsd_threshold && !breakline) {
                    parameters_file << std::endl;
                    breakline = true;
                }
                parameters_file << i.first << " " << i.second[0] << " " << i.second[1] << " " << i.second[2] << std::endl;
            }
        }
        if (m_verbosity >= 1) {
            CurcumaLogger::success_fmt("Initial Pass finished after {:.3f} seconds", timer.Elapsed() / 1000.0);
        }
    } else {
        if (m_verbosity >= 1) {
            CurcumaLogger::info("Skipping initial pass - Setting thresholds to high value");
        }

        for (const auto& i : m_ordered_list)
            m_stored_structures.push_back(m_molecules.at(i.second).second);
        if (m_rmsd_set) {
            m_dLI = 1e23;
            m_dLH = 1e23;
            m_dLE = 1e23;
        } else {
            m_dLI = 0;
            m_dLH = 0;
            m_dLE = 0;
        }
        m_looseThresh = 0;
        m_skipreuse = true;
    }

    if (!m_skipreorder) {
        std::ofstream parameters_skip;
        std::ofstream parameters_performed;

        if (m_analyse) {
            parameters_skip.open(m_skip_file, std::ios_base::app);
            parameters_skip << "# RMSD(old)\tDelta E\tDelta H\tDelta I" << std::endl;
            parameters_performed.open(m_perform_file, std::ios_base::app);
            parameters_performed << "# RMSD(old)\tDelta E\tDelta H\tDelta I" << std::endl;
        }

        for (int run = 0; run < m_sLE.size(); ++run) {
            m_current_filename = m_2nd_filename + "." + std::to_string(run + 1) + ".xyz";

            std::ofstream nd_file;
            if (m_writeFiles && !m_reduced_file) {
                nd_file.open(m_current_filename);
                nd_file.close();
            }
            double dLI = m_dLI;
            double dLH = m_dLH;
            double dLE = m_dLE;
            if (!m_mapped) {
                dLI = m_dLI * m_sLI[run];
                dLH = m_dLH * m_sLH[run];
                dLE = m_dLE * m_sLE[run];
                m_print_rmsd = m_rmsd_threshold;
            } else {
                m_print_rmsd = m_sLI[run] * m_rmsd_threshold;
                for (const auto& i : m_listThresh) {
                    if (i.first <= m_sLI[run] * m_rmsd_threshold)
                        dLI = std::max(dLI, i.second[2]);

                    if (i.first <= m_sLH[run] * m_rmsd_threshold)
                        dLH = std::max(dLH, i.second[1]);

                    if (i.first <= m_sLE[run] * m_rmsd_threshold)
                        dLE = std::max(dLE, i.second[0]);
                }
            }
            if (!CheckStop()) {
                timer.Reset();
                if (m_verbosity >= 1) {
                    CurcumaLogger::success("Reorder Pass: Performing RMSD calculation with reordering");
                }
                if (m_writeFiles && !m_reduced_file) {
                    result_file.open(m_statistic_filename, std::ios_base::app);
                    result_file << "Results of Reorder Pass #" << run + 1 << std::endl;
                    result_file.close();
                }
                if (m_analyse) {
                    std::ofstream parameters_success;
                    parameters_success.open(m_success_file, std::ios_base::app);
                    parameters_success << std::endl
                                       << "# " << run << " run" << std::endl;
                }
                Reorder(dLE, dLI, dLH, false);
                PrintStatus("Result Reorder pass:");

                if (m_verbosity >= 1) {
                    CurcumaLogger::success_fmt("Reorder Pass finished after {:.3f} seconds", timer.Elapsed() / 1000.0);
                }
                timer.Reset();
                if (m_analyse) {
                    WriteDotFile(m_result_basename + ".reorder." + std::to_string(run + 1) + ".dot", m_second_content);

                    m_collective_content += "edge [color=red];\n";
                    for (auto node : m_nodes)
                        m_collective_content += node.second;
                    m_nodes.clear();
                    m_collective_content += m_second_content + "\n";
                    m_second_content.clear();
                    if (m_analyse) {
                        parameters_performed << "# " << run << " run" << std::endl
                                             << std::endl;
                        parameters_skip << "# " << run << " run" << std::endl
                                        << std::endl;

                        for (const auto& i : m_listThresh) {
                            std::vector<double> element = { i.second[0], i.second[1], i.second[2] };
                            auto position = std::find(m_list_skipped.begin(), m_list_skipped.end(), element);
                            if (position != m_list_skipped.end()) {
                                parameters_skip << i.first << " " << i.second[0] << " " << i.second[1] << " " << i.second[2] << std::endl;
                                m_list_skipped.erase(position);
                                continue;
                            }

                            position = std::find(m_list_performed.begin(), m_list_performed.end(), element);
                            if (position != m_list_performed.end()) {
                                parameters_performed << i.first << " " << i.second[0] << " " << i.second[1] << " " << i.second[2] << std::endl;
                                m_list_performed.erase(position);
                            }
                        }
                        parameters_performed << std::endl;
                        parameters_skip << std::endl;
                    }
                }
            }
        }
        if (m_analyse) {
            parameters_skip.close();
            parameters_performed.close();
        }
    } else if (m_verbosity >= 1)
        CurcumaLogger::success("Reorder Pass skipped");
    if (!m_skipreuse) {
        if (!CheckStop()) {
            timer.Reset();
            m_current_filename = m_3rd_filename;
            if (m_verbosity >= 1) {
                if (m_reset)
                    CurcumaLogger::success("Reuse Pass: Performing RMSD calculation with stored reorder rules using all structures");
                else
                    CurcumaLogger::success("Reuse Pass: Performing RMSD calculation with stored reorder rules using previously accepted structures");
            }

            if (m_writeFiles && !m_reduced_file) {
                result_file.open(m_statistic_filename, std::ios_base::app);
                PrintStatus("Result Reuse pass:");
                result_file.close();
            }
            m_exclude_list.clear();
            if (m_analyse) {
                std::ofstream parameters_success;
                parameters_success.open(m_success_file, std::ios_base::app);
                parameters_success << std::endl
                                   << "# reuse run" << std::endl;
            }
            Reorder(-1, -1, -1, true, m_reset);
            PrintStatus("Result reuse pass:");

            if (m_verbosity >= 1) {
                CurcumaLogger::success_fmt("Reuse Pass finished after {:.3f} seconds", timer.Elapsed() / 1000.0);
            }
            timer.Reset();

            if (m_analyse) {
                WriteDotFile(m_result_basename + ".reuse.dot", m_second_content);

                m_collective_content += "edge [color=blue];\n";
                for (auto node : m_nodes)
                    m_collective_content += node.second;
                m_nodes.clear();
                m_collective_content += m_second_content + "\n";
            }
        }
    }
    if (m_analyse) {
        std::ofstream energy;
        energy.open(m_result_basename + ".energy.gnuplot");
        energy << "scale = 4063.0/800.0" << std::endl;
        energy << "set terminal pngcairo  transparent size 600*scale,400*scale transparent font \"Noto Sans\" fontscale scale linewidth scale pointscale scale" << std::endl;
        energy << "set encoding utf8" << std::endl;
        energy << "set output '" << m_result_basename << ".energy.png'" << std::endl;
        energy << "set xlabel \"RMSD [Å]\"" << std::endl;
        energy << "set ylabel \"Energy ΔE [kJ/mol]\"" << std::endl;
        energy << "set key left horizontal  font \"Helvetica, 10\" maxrows 1 outside" << std::endl;
        energy << "plot '" << m_result_basename << ".param.dat' using 1:2 pt 10 ps 0.5 lt rgb \"grey\" title \"Parametrisation\", '" << m_result_basename << ".param.skip.dat' using 1:2 pt 10 ps 0.1 lt rgb \"blue\" title \"Reorder skipped\", '" << m_result_basename << ".param.perf.dat' using 1:2 pt 10 ps 0.1 lt rgb \"yellow\" title \"Reorder performed\", '" << m_result_basename << ".param.success.dat' using 2:3 pt 10 ps 0.5 lt rgb \"red\" title \"Reorder successful\", '" << m_result_basename << ".param.limit.dat' using 1:2 with linespoints linestyle 1 notitle" << std::endl;
        energy.close();

        std::ofstream ripser;
        ripser.open(m_result_basename + ".ripser.gnuplot");
        ripser << "scale = 4063.0/800.0" << std::endl;
        ripser << "set terminal pngcairo  transparent size 600*scale,400*scale transparent font \"Noto Sans\" fontscale scale linewidth scale pointscale scale" << std::endl;
        ripser << "set encoding utf8" << std::endl;
        ripser << "set xlabel \"RMSD [Å]\"" << std::endl;
        ripser << "set ylabel \"ΔH\"" << std::endl;
        ripser << "set output '" << m_result_basename << ".ripser.png'" << std::endl;
        ripser << "set key left horizontal  font \"Helvetica, 10\" maxrows 1 outside" << std::endl;
        ripser << "plot '" << m_result_basename << ".param.dat' using 1:3 pt 10 ps 0.5 lt rgb \"grey\" title \"Parametrisation\", '" << m_result_basename << ".param.skip.dat' using 1:3 pt 10 ps 0.1 lt rgb \"blue\" title \"Reorder skipped\", '" << m_result_basename << ".param.perf.dat' using 1:3 pt 10 ps 0.1 lt rgb \"yellow\" title \"Reorder performed\",  '" << m_result_basename << ".param.success.dat' using 2:4 pt 10 ps 0.5 lt rgb \"red\" title \"Reorder successful\", '" << m_result_basename << ".param.limit.dat' using 1:3 with linespoints linestyle 1 notitle" << std::endl;
        ripser.close();

        std::ofstream rotational;
        rotational.open(m_result_basename + ".rotational.gnuplot");
        rotational << "scale = 4063.0/800.0" << std::endl;
        rotational << "set terminal pngcairo  transparent size 600*scale,400*scale transparent font \"Noto Sans\" fontscale scale linewidth scale pointscale scale" << std::endl;
        rotational << "set encoding utf8" << std::endl;
        rotational << "set xlabel \"RMSD [Å]\"" << std::endl;
        rotational << "set ylabel \"ΔI [MHz]\"" << std::endl;
        rotational << "set output '" << m_result_basename << ".rotational.png'" << std::endl;
        rotational << "set key left horizontal  font \"Helvetica, 10\" maxrows 1 outside" << std::endl;
        rotational << "plot '" << m_result_basename << ".param.dat' using 1:4 pt 10 ps 0.5 lt rgb \"grey\" title \"Parametrisation\", '" << m_result_basename << ".param.skip.dat' using 1:4 pt 10 ps 0.1 lt rgb \"blue\" title \"Reorder skipped\", '" << m_result_basename << ".param.perf.dat' using 1:4 pt 10 ps 0.1 lt rgb \"yellow\" title \"Reorder performed\",  '" << m_result_basename << ".param.success.dat' using 2:5 pt 10 ps 0.5 lt rgb \"red\" title \"Reorder successful\", '" << m_result_basename << ".param.limit.dat' using 1:4 with linespoints linestyle 1 notitle" << std::endl;
        rotational.close();
    }

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
    Finalise();
    if (m_analyse) {
        std::ofstream dotfile;
        dotfile.open(m_result_basename + ".dot");
        dotfile << "digraph graphname \n {\n";
        dotfile << m_collective_content;
        dotfile << "}";
    }
}

void ConfScan::CheckOnly(double sLE, double sLI, double sLH)
{
    std::string content_dot;
    std::string laststring;
    m_maxmol = m_ordered_list.size();

    // Export RMSD config from ConfigManager - Claude Generated 2025
    json rmsd = m_config.exportModule("rmsd");
    rmsd["verbosity"] = 0;  // Override for thread silence
    rmsd["silent"] = true;
    rmsd["check"] = CheckConnections();
    rmsd["noreorder"] = true;

    std::vector<ConfScanThreadNoReorder*> threads;
    std::vector<std::vector<int>> rules;
    m_energies.clear();

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
            m_first_node = mol1->Name();
            ConfScanThreadNoReorder* thread = addThreadNoreorder(mol1, rmsd);
            threads.push_back(thread);
            p->addThread(thread);
            m_all_structures.push_back(mol1);

            m_lowest_energy = mol1->Energy();
            if (m_analyse) {
                laststring = mol1->Name();

                std::string node = "\"" + mol1->Name() + "\" [shape=box, label=\"" + mol1->Name() + "\\n0.00 kj/mol\",id=\"" + mol1->Name() + "\"];\n";
                node += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\"];\n";
                m_nodes.insert(std::pair<double, std::string>(mol1->Energy(), node));
            }
            continue;
        }
        if (m_analyse) {
            std::string node = "\"" + mol1->Name() + "\" [shape=box, label=\"" + mol1->Name() + "\\n" + to_string_with_precision((mol1->Energy() - m_lowest_energy) * 2625.5) + " kJ/mol\",id=\"" + mol1->Name() + "\"];\n";
            node += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\\n" + to_string_with_precision((mol1->Energy() - m_lowest_energy) * 2625.5) + " kJ/mol\"];\n";
            m_nodes.insert(std::pair<double, std::string>(mol1->Energy(), node));
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
        double min_rmsd = 1e4;
        for (auto* t : threads) {
            m_listThresh.insert(std::pair<double, std::vector<double>>(t->RMSD(), { (std::abs(t->Reference()->Energy() - mol1->Energy()) * 2625.5), t->DH(), t->DI() }));

            if (!m_rmsd_set) {
                min_rmsd = std::min(min_rmsd, t->RMSD());
                keep_molecule = true;

            } else {
                if (!m_mapped) {
                    //    m_dLI = std::max(m_dLI, t->DI() * (t->RMSD() <= (sLI * m_rmsd_threshold)));
                    //    m_dLH = std::max(m_dLH, t->DH() * (t->RMSD() <= (sLH * m_rmsd_threshold)));
                    //    m_dLE = std::max(m_dLE, (std::abs(t->Reference()->Energy() - mol1->Energy()) * 2625.5) * (t->RMSD() <= (sLE * m_rmsd_threshold)));
                }

                m_dTI = std::max(m_dTI, t->DI() * (t->RMSD() <= (m_sTI * m_rmsd_threshold)));
                m_dTH = std::max(m_dTH, t->DH() * (t->RMSD() <= (m_sTH * m_rmsd_threshold)));
                m_dTE = std::max(m_dTE, (std::abs(t->Reference()->Energy() - mol1->Energy()) * 2625.5) * (t->RMSD() <= (m_sTE * m_rmsd_threshold)));

                if (t->KeepMolecule() == false) {
                    keep_molecule = false;
                    writeStatisticFile(t->Reference(), mol1, t->RMSD());
                    if (m_analyse) {
                        if (laststring.compare("") != 0 && laststring.compare(t->Reference()->Name()) != 0)
                            m_first_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\" [style=dotted,arrowhead=onormal];\n";

                        std::string node = "\"" + t->Reference()->Name() + "\" [shape=box, label=\"" + t->Reference()->Name() + "\\n" + to_string_with_precision(m_dE) + " kJ/mol\",id=\"" + t->Reference()->Name() + "\"];\n";
                        node += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\\n" + to_string_with_precision((mol1->Energy() - m_lowest_energy) * 2625.5) + " kJ/mol\"];\n";
                        m_nodes.insert(std::pair<double, std::string>(t->Reference()->Energy(), node));
                        m_first_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(t->RMSD()) + "];\n";
                        //  m_nodes_list.push_back(t->Reference()->Name());
                    }
                    laststring = t->Reference()->Name();

#ifdef WriteMoreInfo
                m_dnn_data.push_back(t->getDNNInput());
#endif
                break;
                }
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
        if (!m_rmsd_set) {
            m_rmsd_threshold = std::min(min_rmsd, m_rmsd_threshold);
            // std::cout << "RMSD threshold set to " << m_rmsd_threshold << " Å" << "obtained (" <<  min_rmsd << ")" << std::endl;
        }
        PrintStatus();
        m_all_structures.push_back(mol1);
    }
    p->clear();
    delete p;
    if (!m_rmsd_set) {
        if (m_verbosity >= 2)
            CurcumaLogger::info_fmt("RMSD threshold set to {:.6f} Å", m_rmsd_threshold);
        for (const auto& i : m_listThresh) {
            if (i.first > m_getrmsd_thresh)
                break;
            m_dLI = std::max(m_dLI, i.second[2]);
            m_dLH = std::max(m_dLH, i.second[1]);
            m_dLE = std::max(m_dLE, i.second[0]);
            // std::cout << i.first << " " << i.second[0] << " " << i.second[1] << " " << i.second[2] << std::endl;
        }
    }

    for (const auto& i : m_listThresh) {
        if (i.first > m_getrmsd_thresh)
            break;
        m_dLI = std::max(m_dLI, i.second[2]);
        m_dLH = std::max(m_dLH, i.second[1]);
        m_dLE = std::max(m_dLE, i.second[0]);
        // std::cout << i.first << " " << i.second[0] << " " << i.second[1] << " " << i.second[2] << std::endl;
    }

    m_rmsd_set = true;

    // ConfStat will be called once at the end in Finalise() to avoid interference with status updates
}

void ConfScan::PrintSetUp(double dLE, double dLI, double dLH)
{
    // Claude Generated: Convert threshold display to new logging system
    if (m_verbosity >= 2) {
        CurcumaLogger::info("```");
        CurcumaLogger::info("* Thresholds in Delta I (averaged over Ia, Ib and Ic):");
        CurcumaLogger::info_fmt("  Loose Threshold: {:.2f} MHz \t Tight Threshold: {:.2f} MHz", dLI, m_dTI);
        CurcumaLogger::info("* Thresholds Delta H:");
        CurcumaLogger::info_fmt("  Loose Threshold: {:.2f} \t Tight Threshold: {:.2f}", dLH, m_dTH);
        CurcumaLogger::info("* Thresholds Delta E:");
        CurcumaLogger::info_fmt("  Loose Threshold: {:.2f} kJ/mol \t Tight Threshold: {:.2f} kJ/mol", dLE, m_dTE);
        CurcumaLogger::info("```");
    }

    if (dLE > 0 || dLH > 0 || dLI > 0) {
        std::ofstream parameters_limit;
        if (m_analyse) {
            parameters_limit.open(m_limit_file, std::ios_base::app);
            parameters_limit << "0\t" << dLE << "\t" << dLH << "\t" << dLI << std::endl;
            parameters_limit << std::prev(m_listThresh.end())->first << "\t" << dLE << "\t" << dLH << "\t" << dLI << std::endl;
            parameters_limit << std::endl;
            parameters_limit << m_print_rmsd << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
            parameters_limit << m_print_rmsd << "\t" << dLE << "\t" << dLH << "\t" << dLI << std::endl;
            parameters_limit << std::endl;
            parameters_limit.close();
        }
    }
}

void ConfScan::Reorder(double dLE, double dLI, double dLH, bool reuse_only, bool reset)
{
    std::string content_dot;
    std::string laststring;
    PrintSetUp(dLE, dLI, dLH);
    /* if ignore rotation or ignore barcode are true, the tight cutoffs are set to -1, that no difference is smaller
     * than the tight threshold and the loose threshold is set to a big number, that every structure has to be reordered*/

    m_rejected_directly = 0;
    m_duplicated = 0;
    if (m_ignoreRotation == true) {
        m_dLI = 1e10;
        m_dTI = -1;
    }
    if (m_ignoreBarCode == true) {
        m_dLH = 1e10;
        m_dTH = -1;
    }

    TriggerWriteRestart();
    m_reorder_count += m_reordered;
    m_reorder_successfull_count += m_reordered_worked;
    m_skipped_count += m_skiped;
    m_rejected = 0, m_accepted = 0, m_reordered = 0, m_reordered_worked = 0, m_reordered_reused = 0, m_skiped = 0;
    // Export RMSD config from ConfigManager - Claude Generated 2025
    json rmsd = m_config.exportModule("rmsd");
    rmsd["verbosity"] = 0;  // Override for thread silence
    rmsd["silent"] = true;
    rmsd["reorder"] = true;
    rmsd["threads"] = 1;  // Override: 1 thread per RMSD calculation
    std::vector<Molecule*> cached;
    if (reset)
        cached = m_all_structures;
    else
        cached = m_stored_structures;
    m_maxmol = cached.size();

    m_result = m_previously_accepted;
    m_stored_structures.clear();
    m_energies.clear();
    std::vector<ConfScanThread*> threads;
    std::vector<std::vector<int>> rules;
    CxxThreadPool* p = new CxxThreadPool;
    p->setActiveThreadCount(m_threads);

    std::ofstream parameters_success;
    if (m_analyse) {
        parameters_success.open(m_success_file, std::ios_base::app);
    }

    for (Molecule* mol1 : cached) {
        if (m_result.size() == 0) {
            AcceptMolecule(mol1);
            ConfScanThread* thread = addThread(mol1, rmsd, reuse_only);
            threads.push_back(thread);
            p->addThread(thread);
            m_lowest_energy = mol1->Energy();
            if (m_analyse) {
                std::string node = "\"" + mol1->Name() + "\" [shape=box, label=\"" + mol1->Name() + "\\n" + to_string_with_precision((mol1->Energy() - m_lowest_energy) * 2625.5) + " kJ/mol\",id=\"" + mol1->Name() + "\"];\n";
                node += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\\n" + to_string_with_precision((mol1->Energy() - m_lowest_energy) * 2625.5) + " kJ/mol\"];\n";
                m_nodes.insert(std::pair<double, std::string>(mol1->Energy(), node));
            }
            continue;
        }
        if (m_analyse) {
            std::string node = "\"" + mol1->Name() + "\" [shape=box, label=\"" + mol1->Name() + "\\n" + to_string_with_precision((mol1->Energy() - m_lowest_energy) * 2625.5) + " kJ/mol\",id=\"" + mol1->Name() + "\"];\n";
            node += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\\n" + to_string_with_precision((mol1->Energy() - m_lowest_energy) * 2625.5) + " kJ/mol\"];\n";
            m_nodes.insert(std::pair<double, std::string>(mol1->Energy(), node));
        }
        p->Reset();
        m_current_energy = mol1->Energy();
        m_dE = (m_current_energy - m_lowest_energy) * 2625.5;
        if (m_dE > m_energy_cutoff && m_energy_cutoff != -1) {
            if (m_verbosity >= 1)
                CurcumaLogger::warn_fmt("Energy of {} is {:.2f} kJ/mol, above cutoff of {:.2f} kJ/mol, skipping", mol1->Name(), m_current_energy, m_energy_cutoff);
            break;
        }
        bool keep_molecule = true;
        bool reorder = false;
        for (int t = 0; t < threads.size(); ++t) {
            if (CheckStop()) {
                CurcumaLogger::warn("Found stop file, will end now!");
                TriggerWriteRestart();
                return;
            }
            const Molecule* mol2 = threads[t]->Reference();
            std::pair<std::string, std::string> names(mol1->Name(), mol2->Name());

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
                if (std::find(m_exclude_list.begin(), m_exclude_list.end(), names) != m_exclude_list.end()) {
                    m_duplicated++;
                    m_list_performed.push_back({ std::abs(mol1->Energy() - mol2->Energy()) * 2625.5, dH, dI });
                    continue;
                }
                reorder = true;
                threads[t]->setEnabled(true);
                int tightThresh = 1 * (dI < m_dTI) + 2 * (dH < m_dTH) + 4 * ((std::abs(mol1->Energy() - mol2->Energy()) * 2625.5 < m_dTE));

                if ((tightThresh & m_tightThresh) == m_tightThresh) {
                    if (m_verbosity >= 1)
                        CurcumaLogger::info_fmt("Differences {:.3f} MHz and {:.3f} below tight threshold, reject molecule directly!", dI, dH);
                    m_lastDI = dI;
                    m_lastDH = dH;
                    writeStatisticFile(mol1, mol2, -1, false);
                    m_threshold.push_back(mol2);
                    m_rejected_directly++;
                    reorder = false;
                    keep_molecule = false;
                    break;
                }
                m_list_performed.push_back({ std::abs(mol1->Energy() - mol2->Energy()) * 2625.5, dH, dI });
                m_exclude_list.push_back(std::pair<std::string, std::string>(mol1->Name(), mol2->Name()));
            } else {
                threads[t]->setEnabled(false);
                m_list_skipped.push_back({ std::abs(mol1->Energy() - mol2->Energy()) * 2625.5, dH, dI });
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

            if (m_RMSDmethod.compare("molalign") != 0 || m_threads != 1) {
                // if (m_threads > 2)
                //     p->StaticPool();
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
                    m_reordered_worked += t->ReorderWorked();
                    m_reordered_reused += t->ReusedWorked();
                    if (AddRules(t->ReorderRule()))
                        rules.push_back(t->ReorderRule());
                    if (keep_molecule) {
                        if (m_analyse) {
                            if (laststring.compare("") != 0 && laststring.compare(t->Reference()->Name()) != 0)
                                m_second_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\" [style=dotted,arrowhead=onormal];\n";

                            std::string node = "\"" + t->Reference()->Name() + "\" [shape=box, label=\"" + t->Reference()->Name() + "\\n" + to_string_with_precision(m_dE) + " kJ/mol\",id=\"" + t->Reference()->Name() + "\"];\n";
                            node += "\"" + mol1->Name() + "\" [label=\"" + mol1->Name() + "\\n" + to_string_with_precision((mol1->Energy() - m_lowest_energy) * 2625.5) + " kJ/mol\",id = \"" + mol1->Name() + "\"];\n";
                            m_nodes.insert(std::pair<double, std::string>(t->Reference()->Energy(), node));

                            m_second_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(t->RMSD()) + "];\n";
                            laststring = t->Reference()->Name();
                            auto i = t->getDNNInput();

                            parameters_success << t->RMSD() << " " << t->OldRMSD() << " " << i.dE << " " << i.dH << " " << (i.dIa + i.dIb + i.dIc) * third << std::endl;
                            // m_nodes_list.push_back(t->Reference()->Name());
                        }
                        writeStatisticFile(t->Reference(), mol1, t->RMSD(), true, t->ReorderRule());
                        mol1->ApplyReorderRule(t->ReorderRule());
                    }
                    keep_molecule = false;

                    // break;
                } else { /* Only if additional molaign is invoked */
                    if ((m_domolalign > 1) && t->RMSD() < m_domolalign * m_rmsd_threshold) {
                        if (m_verbosity >= 1)
                            CurcumaLogger::info("Starting molalign for more precise reordering...");
                        json molalign = rmsd;
                        molalign["method"] = "molalign";
                        molalign["verbosity"] = 0;
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
                                m_second_content += "\"" + laststring + "\" -> \"" + t->Reference()->Name() + "\" [style=dotted,arrowhead=onormal];\n";

                            m_second_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\";\n";
                            m_second_content += "\"" + mol1->Name() + "\" [shape=box, style=filled,color=\".7 .3 1.0\"];\n";
                            m_second_content += "\"" + t->Reference()->Name() + "\" -> \"" + mol1->Name() + "\" [style=bold,label=" + std::to_string(driver.RMSD()) + "];\n";
                            laststring = t->Reference()->Name();

                            if (m_verbosity >= 2)
                                CurcumaLogger::success("... success!");
                        } else if (m_verbosity >= 2)
                            CurcumaLogger::warn("... without effect!");
                    }
                }
                m_kuhn_munkres_iterations += t->getKuhnMunkresIterations();
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
    }
    if (m_analyse) {
        parameters_success.close();
    }
    p->clear();
    delete p;

    // ConfStat will be called once at the end in Finalise() to avoid interference with status updates
}

ConfScanThreadNoReorder* ConfScan::addThreadNoreorder(const Molecule* reference, const json& config)
{
    ConfScanThreadNoReorder* thread = new ConfScanThreadNoReorder(m_rmsd_threshold, m_MaxHTopoDiff, config, m_verbosity);
    thread->setReference(*reference);
    m_energies.push_back(reference->Energy());
    return thread;
}

ConfScanThread* ConfScan::addThread(const Molecule* reference, const json& config, bool reuse_only)
{
    ConfScanThread* thread = new ConfScanThread(m_reorder_rules, m_rmsd_threshold, m_MaxHTopoDiff, reuse_only, config, m_verbosity);
    thread->setReference(*reference);
    thread->setEarlyBreak(m_earlybreak);
    thread->setVerbose(m_analyse);
    m_energies.push_back(reference->Energy());

    // thread->setVerbose(false);
    return thread;
}

void ConfScan::Finalise()
{
    TriggerWriteRestart();

    if (m_verbosity >= 1) {
        // Ensure logger verbosity matches local verbosity for these messages
        int old_verbosity = CurcumaLogger::get_verbosity();
        CurcumaLogger::set_verbosity(m_verbosity);

        CurcumaLogger::info("Final Statistics:");
        CurcumaLogger::info_fmt("Rotational constants: {} ms", m_timing_rot);
        CurcumaLogger::info_fmt("Ripser bar code difference: {} ms", m_timing_ripser);

        // Calculate success rate with bounds checking
        double success_rate = (m_reorder_count > 0) ? (m_reorder_successfull_count / double(m_reorder_count) * 100.0) : 0.0;
        CurcumaLogger::info_fmt("Success rate: {:.2f}%", success_rate);

        // Calculate efficiency with bounds checking
        double efficiency = (m_reorder_successfull_count > 0) ? (m_skipped_count / double(m_reorder_successfull_count)) : 0.0;
        CurcumaLogger::info_fmt("Efficiency: {:.3f}", efficiency);

        CurcumaLogger::info_fmt("Kuhn-Munkres iterations: {}", m_kuhn_munkres_iterations);

        // Restore original logger verbosity
        CurcumaLogger::set_verbosity(old_verbosity);
    }

    int i = 0;
    m_collective_content += "subgraph cluster_bevor {\nrank = same;\nstyle= invis;\n";
    std::string content_after;
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
        if (m_analyse) {
            std::string content = "\"" + molecule->Name() + "\" [shape=box, label=\"" + molecule->Name() + "\\n" + to_string_with_precision(difference) + " kJ/mol\", fontcolor=\"orange\", fontname=\"times-bold\",id=\"" + molecule->Name() + "\"];\n";
            content_after += content;
            if (std::find(m_nodes_list.begin(), m_nodes_list.end(), molecule->Name()) != m_nodes_list.end()) {
                //    std::cout << molecule->Name() << " inside ";
            } else {
                std::string content = "\"" + molecule->Name() + "\";\n";
                m_collective_content += content;
                content_after += "\"" + molecule->Name() + "\" -> \"" + m_first_node + "\" [style=invis];\n";
            }
        }
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
    m_collective_content += "}\n";
    m_collective_content += "\"" + m_first_node + "\";\n";
    m_collective_content += content_after;

    // Always show final result regardless of verbosity - this is critical information
    CurcumaLogger::success_fmt("{} structures were kept - of {} total!",
        m_stored_structures.size(), m_molecules.size() - m_fail);

    // Show conformer statistics at the end - this is the final analysis
    if (!m_energies.empty()) {
        ConfStat stat;
        stat.setEnergies(m_energies);
        stat.start();
    }
}

bool ConfScan::AddRules(const std::vector<int>& rules)
{
    if (rules.size() == 0 || m_skip_orders)
        return false;

    if (std::find(m_reorder_rules.begin(), m_reorder_rules.end(), rules) == m_reorder_rules.end()) {
        m_reorder_rules.push_back(rules);
    }
    return true;
}

void ConfScan::PrintStatus(const std::string& info)
{
    // Always show critical status updates for user feedback during calculations
    if (m_verbosity >= 1) {
        if (!info.empty()) {
            fmt::print("        {}\n", info);
        }

        // Show progress percentage
        fmt::print("        Progress: {:.1f}% done\n",
            (m_stored_structures.size() + m_rejected) / double(m_maxmol) * 100);

        // Show detailed status line - this is the key information users expect to see
        std::string status = fmt::format("Accepted: {}, Rejected: {}, Reordered: {} (+{})",
            m_stored_structures.size(), m_rejected, m_reordered, m_molalign_count);
        status += fmt::format(" | Successfully: {} (+{}), Reused Results: {}",
            m_reordered_worked, m_molalign_success, m_reordered_reused);
        status += fmt::format(" | Reordering Skipped: {} (+{}), Rejected Directly: {}",
            m_skiped, m_duplicated, m_rejected_directly);
        status += fmt::format(" | Current Energy: {:.2f} kJ/mol", m_dE);
        fmt::print("        {}\n", status);

        std::cout << std::endl; // Force immediate flush for real-time feedback
    }
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
    m_nodes_list.push_back(mol1->Name());
}
