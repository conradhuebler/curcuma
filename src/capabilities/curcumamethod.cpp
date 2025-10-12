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

#include "src/core/curcuma_logger.h"
#include "src/core/global.h"
#include "src/global_config.h"

#include "src/tools/general.h"

#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdio.h>

#ifdef C17
#ifndef _WIN32
#include <filesystem>
namespace fs = std::filesystem;
#endif
#endif
#include "curcumamethod.h"

CurcumaMethod::CurcumaMethod(const json& defaults, const json& controller, bool silent)
    : m_defaults(defaults)
    , m_controller(controller)
    , m_silent(silent)
{
    // Legacy constructor - convert boolean silent to verbosity levels
    if (controller.count("verbose") > 0) {
        m_silent = false;
        m_verbose = true;
        m_verbosity = 3; // Verbose = Informative Print
        CurcumaLogger::set_verbosity(m_verbosity); // Claude Generated - Sync to logger
    } else {
        m_verbosity = silent ? 0 : 1; // Silent or Normal Print
        CurcumaLogger::set_verbosity(m_verbosity); // Claude Generated - Sync to logger
    }

    // Check for explicit verbosity level in CLI arguments - Claude Generated
    if (controller.count("verbosity") > 0) {
        try {
            m_verbosity = controller["verbosity"].get<int>();
            // Clamp to valid range 0-3
            m_verbosity = std::max(0, std::min(3, m_verbosity));
            // Update legacy flags for backwards compatibility
            m_silent = (m_verbosity == 0);
            m_verbose = (m_verbosity >= 3);
            CurcumaLogger::set_verbosity(m_verbosity); // Claude Generated - Sync to logger
        } catch (const std::exception& e) {
            // Invalid verbosity value, keep current setting
        }
    }

    // Check for threads setting in CLI arguments - Claude Generated
    if (controller.count("threads") > 0) {
        try {
            m_threads = controller["threads"].get<int>();
            m_threads = std::max(1, m_threads); // At least 1 thread
        } catch (const std::exception& e) {
            // Invalid threads value, keep default
        }
    }

    controller.count("help") > 0 ? m_help = true : m_help = false;

    // Initialize CurcumaLogger with this method's verbosity - Claude Generated
    CurcumaLogger::set_verbosity(m_verbosity);

    // m_curcuma_progress.open("curcuma_progress", std::ios::out);
}

// New verbosity constructor - Claude Generated
CurcumaMethod::CurcumaMethod(const json& defaults, const json& controller, int verbosity)
    : m_defaults(defaults)
    , m_controller(controller)
    , m_verbosity(verbosity)
{
    // Set legacy flags for backwards compatibility
    m_silent = (verbosity == 0);
    m_verbose = (verbosity >= 3);
    CurcumaLogger::set_verbosity(m_verbosity); // Claude Generated - Sync to logger

    // Check for verbosity override in controller
    if (controller.count("verbosity") > 0) {
        try {
            m_verbosity = controller["verbosity"].get<int>();
            // Clamp to valid range 0-3
            m_verbosity = std::max(0, std::min(3, m_verbosity));
            m_silent = (m_verbosity == 0);
            m_verbose = (m_verbosity >= 3);
            CurcumaLogger::set_verbosity(m_verbosity); // Claude Generated - Sync to logger
        } catch (const std::exception& e) {
            // Invalid verbosity value, keep parameter setting
        }
    }

    // Check for threads setting in CLI arguments - Claude Generated
    if (controller.count("threads") > 0) {
        try {
            m_threads = controller["threads"].get<int>();
            m_threads = std::max(1, m_threads); // At least 1 thread
        } catch (const std::exception& e) {
            // Invalid threads value, keep default
        }
    }

    controller.count("help") > 0 ? m_help = true : m_help = false;

    // Initialize CurcumaLogger with this method's verbosity - Claude Generated
    CurcumaLogger::set_verbosity(m_verbosity);

    //m_curcuma_progress.open("curcuma_progress", std::ios::out);
}

CurcumaMethod::~CurcumaMethod()
{
}

// Claude Generated 2025: Enhanced restart writing with automatic checksum and version
void CurcumaMethod::TriggerWriteRestart()
{
    std::ofstream restart_file("curcuma_restart.json");
    nlohmann::json restart;
    try {
        nlohmann::json state = WriteRestartInformation();

        // Add metadata
        state["format_version"] = "1.0";
        state["timestamp"] = std::time(nullptr);

        // Compute and add checksum for critical string fields
        std::vector<std::string> checksum_fields;
        for (const auto& item : state.items()) {
            if (item.value().is_string()) {
                const std::string& str = item.value().get<std::string>();
                // Only checksum fields that look like data (contain '|' separator)
                if (str.find('|') != std::string::npos) {
                    checksum_fields.push_back(item.key());
                }
            }
        }

        if (!checksum_fields.empty()) {
            size_t checksum = computeRestartChecksum(state, checksum_fields);
            state["checksum"] = checksum;
            state["checksum_fields"] = checksum_fields;
        }

        restart[MethodName()[0]] = state;
    } catch (nlohmann::json::type_error& e) {
    }
    try {
        restart_file << restart << std::endl;
    } catch (nlohmann::json::type_error& e) {
    }
}

StringList CurcumaMethod::RestartFiles() const
{
    StringList file_list;

#ifdef C17
#ifndef _WIN32
    for (auto& p : fs::directory_iterator(".")) {
        std::string file(p.path());
        if (file.find("curcuma_restart") != std::string::npos)
            file_list.push_back(file);
    }
#endif
#else
    std::ifstream test_file("curcuma_restart.json");
    if (test_file.is_open())
        file_list.push_back("curcuma_restart.json");
#endif

    return file_list;
}

nlohmann::json CurcumaMethod::LoadControl() const
{
    nlohmann::json control;
    std::ifstream restart_file("curcuma_control.json");
    try {
        restart_file >> control;
    } catch (nlohmann::json::type_error& e) {
        throw 404;
    } catch (nlohmann::json::parse_error& e) {
        throw 404;
    }
    std::cout << control << std::endl;
    return control;
}

void CurcumaMethod::UpdateController(const json& controller)
{
    json method;
    for (const auto& s : MethodName()) {
        try {
            method = Json2KeyWord<json>(controller, s);
        } catch (int i) {
            method = controller;
        }
        m_defaults = MergeJson(m_defaults, method);
    }
    if (!m_silent)
        PrintController(m_defaults);
    LoadControlJson();
}

void CurcumaMethod::checkHelp()
{
    if (m_help) {
        printHelp();
        exit(0);
    }
}

bool CurcumaMethod::CheckStop() const
{
#ifdef C17
#ifndef _WIN32
    return std::filesystem::exists("stop");
#endif
    std::ifstream test_file("stop");
    bool result = test_file.is_open();
    test_file.close();
    return result;
#else
    std::ifstream test_file("stop");
    bool result = test_file.is_open();
    test_file.close();
    return result;
#endif
}

void CurcumaMethod::getBasename(const std::string& filename)
{
    std::string name = std::string(filename);
    for (int i = 0; i < 4; ++i)
        name.pop_back();
    m_basename = name;
}

void CurcumaMethod::setFile(const std::string& filename)
{
    m_filename = filename;
    getBasename(filename);
}

// Claude Generated 2025: Compute checksum over specified fields for restart validation
size_t CurcumaMethod::computeRestartChecksum(const json& state, const std::vector<std::string>& fields) const
{
    std::string data_to_hash;
    for (const auto& field : fields) {
        if (state.contains(field) && state[field].is_string()) {
            data_to_hash += state[field].get<std::string>();
        }
    }
    return std::hash<std::string>{}(data_to_hash);
}

// Claude Generated 2025: Validate pipe-separated double strings (e.g., "1.5|2.3|3.1")
bool CurcumaMethod::isValidDoubleString(const std::string& str) const
{
    if (str.empty()) return true;

    std::istringstream ss(str);
    std::string token;
    while (std::getline(ss, token, '|')) {
        // Trim whitespace
        token.erase(0, token.find_first_not_of(" \t\n\r"));
        token.erase(token.find_last_not_of(" \t\n\r") + 1);

        if (token.empty()) continue;

        // Check for truncated NaN/Inf patterns
        if (token.size() < 3) {
            if (token.find("na") != std::string::npos ||
                token.find("in") != std::string::npos ||
                token == "-n" || token == "-i") {
                return false; // Likely truncated "-nan" or "-inf"
            }
        }

        // Try to parse as double
        try {
            std::stod(token);
        } catch (const std::invalid_argument&) {
            return false;
        } catch (const std::out_of_range&) {
            return false;
        }
    }
    return true;
}

// Claude Generated 2025: Universal restart data validation for all capabilities
RestartValidationResult CurcumaMethod::validateRestartData(
    const json& state,
    const std::vector<std::string>& required_fields,
    const std::vector<std::string>& checksum_fields) const
{
    // Check format version (warning only, not fatal)
    if (state.contains("format_version")) {
        std::string version = state["format_version"];
        if (version != "1.0") {
            std::cerr << "\033[1;33m[WARNING]\033[0m Restart file format version " << version
                      << " may be incompatible with current version 1.0" << std::endl;
        }
    }

    // Check required fields
    for (const auto& field : required_fields) {
        if (!state.contains(field)) {
            return {false, "Missing required field: " + field};
        }
    }

    // Validate checksum if present
    if (state.contains("checksum") && state.contains("checksum_fields")) {
        size_t stored_checksum = state["checksum"].get<size_t>();
        std::vector<std::string> stored_fields = state["checksum_fields"].get<std::vector<std::string>>();

        // Recompute checksum
        size_t computed_checksum = computeRestartChecksum(state, stored_fields);

        if (computed_checksum != stored_checksum) {
            return {false, "Checksum mismatch - restart file may be corrupted (disk error or manual edit)"};
        }
    }

    // Validate string fields containing doubles
    for (const auto& field : checksum_fields) {
        if (state.contains(field) && state[field].is_string()) {
            std::string value = state[field];
            if (!isValidDoubleString(value)) {
                return {false, "Field '" + field + "' contains malformed numerical data (e.g., truncated NaN)"};
            }
        }
    }

    return {true, ""};
}
