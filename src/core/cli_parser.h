/*
 * <Command-line parsing: argv -> controller JSON>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

// Claude Generated (Nov 2026): moved unchanged from main.cpp so the routing can be tested without
// the whole program (docs/debt_2026-11/PARAMS.md, migration step 0).
#pragma once

#include "json.hpp"

#include <string>

using json = nlohmann::json;

/*! \brief Translate argv into the controller JSON: controller[module][param] plus globals. */
json CLI2Json(int argc, char** argv);

/*! \brief Set target["a"]["b"] for dotKey "a.b" (flat key when there is no dot). */
void setNestedJsonValue(json& target, const std::string& dotKey, const json& value);

/*! \brief Recursive deep merge; values already present in dst win (CLI over file). */
void mergeJsonRecursive(json& dst, const json& src);
