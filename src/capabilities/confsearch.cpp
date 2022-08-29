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

#include "src/global_config.h"

#include "src/tools/general.h"

#include <fstream>
#include <iostream>
#include <stdio.h>

#include "confsearch.h"

ConfSearch::ConfSearch(const json& controller, bool silent)
    : CurcumaMethod(ConfSearchJson, controller, silent)
{
    UpdateController(controller);
}

ConfSearch::~ConfSearch()
{
}

bool ConfSearch::Initialise()
{
    return true;
}

void ConfSearch::start()
{
}

nlohmann::json ConfSearch::WriteRestartInformation()
{
    nlohmann::json restart;
    return restart;
}

bool ConfSearch::LoadRestartInformation()
{
    return true;
}

void ConfSearch::ReadControlFile()
{
}

void ConfSearch::LoadControlJson()
{
}
