/*
 * <Abstract Curcuma Method, please try to subclass from that!>
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

#pragma once

#include "src/tools/general.h"

#include <string>

#include "json.hpp"

// Claude Generated 2025: Restart file validation result
struct RestartValidationResult {
    bool valid;
    std::string error_message;
};

class CurcumaMethod {
public:
    CurcumaMethod(const json& defaults, const json& controller, bool silent); // Legacy constructor
    CurcumaMethod(const json& defaults, const json& controller, int verbosity); // New verbosity constructor - Claude Generated
    CurcumaMethod()
    {
        m_help = true;
        m_verbosity = 1; // Default: Normal Print
    }
    ~CurcumaMethod();

    inline void setRestart(bool restart) { m_restart = restart; }

    inline bool Restart() const { return m_restart; }

    inline void setController(const json& controller) { m_controller = controller; }

    virtual bool Initialise() { return true; } // TODO make pure virtual

    void printError()
    {
        for (const auto& error : m_error_list)
            std::cerr << error << std::endl;
        m_error_list.clear();
    }

    void UpdateController(const json& controller);

    virtual void start() = 0; // TODO make pure virtual and move all main action here
    virtual void printHelp() const { std::cout << "No help available for this method." << std::endl; };

    bool CheckStop() const;

    std::string Basename() const { return m_basename; }
    std::string Filename() const { return m_filename; }
    void getBasename(const std::string& filename);
    void overrideBasename(const std::string& basename) { m_basename = basename; }
    virtual void setFile(const std::string& filename);

protected:
    void checkHelp();
    void TriggerWriteRestart();

    StringList RestartFiles() const;

    nlohmann::json LoadControl() const;

    // Claude Generated 2025: Universal restart file validation
    size_t computeRestartChecksum(const json& state, const std::vector<std::string>& fields) const;
    RestartValidationResult validateRestartData(const json& state,
                                                 const std::vector<std::string>& required_fields,
                                                 const std::vector<std::string>& checksum_fields) const;
    bool isValidDoubleString(const std::string& str) const;

    json m_defaults, m_controller;
    bool m_restart = true;

    void AppendError(const std::string& error) { m_error_list.push_back(error); }

    // Logging system integration - Claude Generated
    // Claude Generated (Updated 2025-11-03): Sync local verbosity with global logger
    void setVerbosity(int level) {
        m_verbosity = level;
        CurcumaLogger::set_verbosity(level);  // Always sync to logger
    }
    int getVerbosity() const { return m_verbosity; }

    //std::filebuf m_curcuma_progress;
    bool m_silent = true; // Legacy - kept for backwards compatibility
    bool m_verbose = false; // Legacy - kept for backwards compatibility
    bool m_help = false;
    int m_verbosity = 1; // New verbosity system: 0=Silent, 1=Small, 2=Normal, 3=Informative
    int m_threads = 1; // Number of threads for parallel processing - Claude Generated

private:
    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation() = 0;

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation() = 0;

    virtual StringList MethodName() const = 0;

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile() = 0;

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson() = 0;

    StringList m_error_list;

    std::string m_basename;
    std::string m_filename;
};
