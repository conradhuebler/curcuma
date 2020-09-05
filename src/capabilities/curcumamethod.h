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

#pragma once

#include "src/tools/general.h"

#include <string>

#include "json.hpp"

class CurcumaMethod {
public:
    CurcumaMethod(const json& defaults, const json& controller, bool silent);

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

protected:
    void TriggerWriteRestart();

    StringList RestartFiles() const;

    nlohmann::json LoadControl() const;

    json m_defaults, m_controller;
    bool m_restart = true;

    void AppendError(const std::string& error) { m_error_list.push_back(error); }
    std::filebuf m_curcuma_progress;

private:
    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation() = 0;

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation() = 0;

    virtual std::string MethodName() const = 0;

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile() = 0;

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson() = 0;

    StringList m_error_list;
    bool m_silent = true;
};
