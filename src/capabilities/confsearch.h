/*
 * <Conformational Search based on Molecular Dynamics>
 * Copyright (C) 2022 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/capabilities/curcumamethod.h"

static const json ConfSearchJson = {

};

class ConfSearch : public CurcumaMethod {
public:
    ConfSearch(const json& controller, bool silent);
    ~ConfSearch();

    virtual bool Initialise();

    virtual void start();

private:
    /* Lets have this for all modules */
    virtual nlohmann::json WriteRestartInformation();

    /* Lets have this for all modules */
    virtual bool LoadRestartInformation();

    virtual StringList MethodName() const
    {
        return { "ConfSearch" };
    }

    /* Lets have all methods read the input/control file */
    virtual void ReadControlFile();

    /* Read Controller has to be implemented for all */
    virtual void LoadControlJson();

    StringList m_error_list;
    bool m_silent = true;
};
