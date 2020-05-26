/*
 * <Abstract Curcuma Method, please try to subclass from that!>
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

#include "src/global_config.h"

#include "src/tools/general.h"

#include <fstream>
#include <iostream>

#ifdef C17
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include "curcumamethod.h"

CurcumaMethod::CurcumaMethod(const json controller)
    : m_controller(controller)
{
}

void CurcumaMethod::TriggerWriteRestart()
{
    std::ofstream restart_file("curcuma_restart.json");
    nlohmann::json restart;
    try {
        restart[MethodName()] = WriteRestartInformation();

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
    for (auto& p : fs::directory_iterator(".")) {
        std::string file(p.path());
        if (file.find("curcuma_restart") != string::npos)
            file_list.push_back(file);
    }
#else
    ifstream test_file("curcuma_restart.json");
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
