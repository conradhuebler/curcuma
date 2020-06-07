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

#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "src/core/molecule.h"

#include "curcumamethod.h"

constexpr double third = 1 / 3.0;

static const json ConfScanJson = {
    { "noname", true },
    { "restart", true },
    { "heavy", false },
    { "rmsd", 1 },
    { "rank", -1 },
    { "writeXYZ", false },
    { "forceReorder", false },
    { "check", false },
    { "energy", 1.0 },
    { "maxenergy", -1.0 },
    { "preventreorder", false },
    { "silent", false },
    { "scaleLoose", 2 },
    { "scaleTight", 0.5 }
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

    void setMolecules(const std::map<double, Molecule*>& molecules);

    /*! \brief Force Connectivitiy Check */
    inline bool CheckConnections() const { return m_check_connections; }

    /*! \brief Check, if Reordering is forced */
    inline bool ForceReorder() const { return m_force_reorder; }

    /*! \brief Check, if Reordering is forced */
    inline bool PreventReorder() const { return m_prevent_reorder; }

    inline std::string NamePattern(int index) const { return "input_" + std::to_string(index); }

    std::vector<Molecule*> Result() const { return m_result; }
    std::vector<Molecule*> Failed() const { return m_failed; }

    void ParametriseRotationalCutoffs();

    int AcceptRotationalConstant(double constant);

    void start() override; // TODO make pure virtual and move all main action here

private:
    /* Lets have this for all modules */
    nlohmann::json WriteRestartInformation() override;

    /* Lets have this for all modules */
    bool LoadRestartInformation() override;

    std::string MethodName() const override { return std::string("ConfScan"); }

    /* Lets have all methods read the input/control file */
    void ReadControlFile() override;

    /* Read Controller has to be implemented for all */
    void LoadControlJson() override;

    bool openFile();

    std::vector<std::vector<int>> m_reorder_rules;

    int PreCheckAgainstAccepted(int index);
    int CheckTempList(int index);
    void PrintStatus();

    std::map<std::string, std::vector<std::string>> m_filtered;
    bool m_ok;
    std::size_t m_fail = 0, m_start = 0, m_end;
    std::vector<Molecule*> m_global_temp_list;
    int m_rejected = 0, m_accepted = 0, m_reordered = 0, m_reordered_worked = 0, m_reordered_failed_completely = 0, m_reordered_reused = 0;

    std::string m_filename;
    std::map<double, int> m_ordered_list;
    std::vector<std::pair<std::string, Molecule*>> m_molecules;
    double m_energy_threshold = 1.0, m_rmsd_threshold = 1.0, m_diff_rot_rel_loose = 0.3, m_diff_rot_rel_tight = 0.01, m_nearly_missed = 0.8, m_energy_cutoff = -1, m_reference_last_energy = 0, m_target_last_energy = 0;
    double m_diff_rot_abs_tight = 0, m_diff_rot_abs_loose = 0, m_scale_tight = 0.5, m_scale_loose = 2;
    std::vector<Molecule*> m_result, m_nearly, m_failed;
    int m_maxrank = 10000;
    bool m_writeXYZ = false;
    bool m_check_connections = false;
    bool m_force_reorder = false, m_prevent_reorder = false;
    bool m_heavy = false;
    bool m_noname = false;
    bool m_writeFiles = true;
    bool m_useRestart = false;
    bool m_silent = false;
    bool m_internal_parametrised = false;
    bool m_parameter_loaded = false;
    bool m_force_silent = false;
};
