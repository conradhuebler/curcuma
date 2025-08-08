/*
 * <Abstract Curcuma Method, please try to subclass from that!>
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

#include "src/core/global.h"
#include "src/global_config.h"

#include "src/tools/general.h"

#include <fstream>
#include <iostream>
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
    } else {
        m_verbosity = silent ? 0 : 1; // Silent or Normal Print
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

    // Check for verbosity override in controller
    if (controller.count("verbosity") > 0) {
        try {
            m_verbosity = controller["verbosity"].get<int>();
            // Clamp to valid range 0-3
            m_verbosity = std::max(0, std::min(3, m_verbosity));
            m_silent = (m_verbosity == 0);
            m_verbose = (m_verbosity >= 3);
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

void CurcumaMethod::TriggerWriteRestart()
{
    std::ofstream restart_file("curcuma_restart.json");
    nlohmann::json restart;
    try {
        restart[MethodName()[0]] = WriteRestartInformation();
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
