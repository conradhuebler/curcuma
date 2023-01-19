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

#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "src/capabilities/rmsd.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#include "src/core/molecule.h"

#include "curcumamethod.h"

constexpr double third = 1 / 3.0;

static const json ConfScanJson = {
    { "noname", true },
    { "restart", true },
    { "heavy", false },
    { "rmsd", -1 },
    { "rank", -1 },
    { "writeXYZ", false },
    { "forceReorder", false },
    { "check", false },
    { "energy", 1.0 },
    { "maxenergy", -1.0 },
    { "preventreorder", false },
    { "scaleLoose", 1.5 },
    { "scaleTight", 0.1 },
    { "skip", 0 },
    { "allxyz", false },
    { "update", false },
    { "MaxParam", -1 },
    { "UseOrders", -1 },
    { "RMSDMethod", "incr" },
    { "MaxHTopoDiff", -1 },
    { "threads", 1 },
    { "RMSDElement", 7 },
    { "accepted", "" },
    { "method", "" },
    { "lastdE", -1 },
    { "fewerFile", false },
    { "dothird", false },
    { "skipfirst", false }
};

class ConfScanThread : public CxxThread {
public:
    ConfScanThread(const std::vector<std::vector<int>>& reorder_rules, double rmsd_threshold, int MaxHTopoDiff, bool reuse_only, const json& config)
    {
        m_driver = new RMSDDriver(config, true);
        m_config = config;
        m_reuse_only = reuse_only;
        m_reorder_rules = reorder_rules;
        m_rmsd_threshold = rmsd_threshold;
        m_MaxHTopoDiff = MaxHTopoDiff;
        setAutoDelete(false);
    }

    ~ConfScanThread()
    {
        delete m_driver;
    }

    virtual int execute() override;
    virtual bool BreakThreadPool() const override { return false; }

    bool KeepMolecule() const { return m_keep_molecule; }
    bool ReorderWorked() const { return m_reorder_worked; }
    bool ReusedWorked() const { return m_reused_worked; }

    void setReference(const Molecule& molecule)
    {
        m_reference = molecule;
        m_target = molecule;
    }
    void setTarget(const Molecule* molecule)
    {
        // m_target.setGeometry(molecule->getGeometry());
        m_target = *molecule;
    }
    std::vector<int> ReorderRule() const { return m_reorder_rule; }
    void setReorderRules(const std::vector<std::vector<int>>& reorder_rules)
    {
        m_reorder_rules = reorder_rules;
    }
    void addReorderRule(const std::vector<int>& rule)
    {
        m_reorder_rules.push_back(rule);
    }
    void setThreads(int threads) { m_threads = threads; }

private:
    bool m_keep_molecule = true, m_break_pool = false, m_reorder_worked = false, m_reuse_only = false, m_reused_worked = false;
    Molecule m_reference, m_target;
    double m_rmsd = 0, m_rmsd_threshold = 1;
    int m_MaxHTopoDiff;
    int m_threads = 1;
    std::vector<int> m_reorder_rule;
    std::vector<std::vector<int>> m_reorder_rules;
    RMSDDriver* m_driver;
    json m_config;
};

class ConfScan : public CurcumaMethod {
public:
    ConfScan(const json& controller = ConfScanJson, bool silent = true);
    virtual ~ConfScan();

    void setFileName(const std::string& filename)
    {
        m_filename = filename;
        openFile();
    }

    // void setMolecules(const std::map<double, Molecule*>& molecules);

    /*! \brief Force Connectivitiy Check */
    inline bool CheckConnections() const { return m_check_connections; }

    /*! \brief Check, if Reordering is forced */
    inline bool ForceReorder() const { return m_force_reorder; }

    /*! \brief Check, if Reordering is forced */
    inline bool PreventReorder() const { return m_prevent_reorder; }

    inline std::string NamePattern(int index) const { return "input_" + std::to_string(index); }

    std::vector<Molecule*> Result() const { return m_result; }
    // std::vector<Molecule*> Failed() const { return m_failed; }

    void ParametriseRotationalCutoffs();

    void start() override; // TODO make pure virtual and move all main action here
    ConfScanThread* addThread(const Molecule* reference, const json& config, bool reuse_only);

private:
    void SetUp();

    void CheckRMSD();
    bool SingleCheckRMSD(const Molecule* mol1, const Molecule* mol2, RMSDDriver* driver);

    void ReorderCheck(bool reuse_only = false, bool limit = false);
    bool SingleReorderRMSD(const Molecule* mol1, const Molecule* mol2, RMSDDriver* driver, bool reuse_only);

    void writeStatisticFile(const Molecule* mol1, const Molecule* mol2, double rmsd, bool reason = true);

    void Finalise();

    void AcceptMolecule(Molecule* molecule);
    void RejectMolecule(Molecule* molecule);

    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override;

    /* Lets have this for all modules */
    bool LoadRestartInformation() override;

    StringList MethodName() const override { return { std::string("ConfScan") }; }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override;

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    bool AddRules(const std::vector<int>& rules);

    bool openFile();

    std::vector<std::vector<int>> m_reorder_rules;

    void PrintStatus();

    std::map<std::string, std::vector<std::string>> m_filtered;
    bool m_ok;
    std::size_t m_fail = 0, m_start = 0, m_end;
    std::vector<Molecule*> m_global_temp_list;
    int m_rejected = 0, m_accepted = 0, m_reordered = 0, m_reordered_worked = 0, m_reordered_failed_completely = 0, m_reordered_reused = 0, m_skip = 0;

    std::string m_filename, m_accepted_filename, m_1st_filename, m_2nd_filename, m_rejected_filename, m_result_basename, m_statistic_filename, m_prev_accepted, m_joined_filename, m_threshold_filename, m_current_filename;
    std::map<double, int> m_ordered_list;
    std::vector<std::pair<std::string, Molecule*>> m_molecules;
    double m_energy_threshold = 1.0, m_rmsd_threshold = 1.0, m_diff_rot_rel_loose = 0.3, m_diff_rot_rel_tight = 0.01, m_nearly_missed = 0.8, m_energy_cutoff = -1, m_reference_last_energy = 0, m_target_last_energy = 0, m_lowest_energy = 1, m_current_energy = 0;
    double m_diff_rot_abs_tight = 0, m_diff_rot_abs_loose = 0, m_scale_tight = 0.5, m_scale_loose = 2;
    double m_reference_restored_energy = -1e10, m_target_restored_energy = -1e10;
    double m_diff_rot_threshold_loose = 0.0, m_diff_ripser_threshold_loose = 0.0, m_diff_rot_threshold_tight = 0.0, m_diff_ripser_threshold_tight = 0.0;
    std::vector<Molecule*> m_result, m_rejected_structures, m_stored_structures, m_previously_accepted, m_threshold;
    std::vector<int> m_element_templates;
    std::string m_rmsd_element_templates;
    std::string m_method = "";

    double m_last_diff = 0.0, m_last_ripser = 0.0, m_last_dE = -1, m_dE = -1;
    int m_maxmol = 0;
    int m_maxrank = 10000;
    int m_maxParam = -1;
    int m_useorders = 10;
    std::string m_RMSDmethod = "incr";
    int m_MaxHTopoDiff = -1;
    int m_threads = 1;
    int m_RMSDElement = 7;
    bool m_writeXYZ = false;
    bool m_check_connections = false;
    bool m_force_reorder = false, m_prevent_reorder = false;
    bool m_heavy = false;
    bool m_noname = false;
    bool m_writeFiles = true;
    bool m_useRestart = false;
    bool m_internal_parametrised = false;
    bool m_parameter_loaded = false;
    bool m_force_silent = false;
    bool m_allxyz = false;
    bool m_update = false;
    bool m_reduced_file = false;
    bool m_do_third = false;
    bool m_skipfirst = false;
};
