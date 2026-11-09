/*
 * <Normalisation of parameter keys>
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

#pragma once

#include <algorithm>
#include <cctype>
#include <string>

namespace curcuma {

/// Claude Generated (Nov 2026): the one rule by which parameter keys are compared in the
/// registry, in ConfigManager and in Json2KeyWord: lower case, hyphens removed. "Max-Iter" and
/// "maxiter" are the same key; "max_iter" is a different one (underscores are kept).
inline std::string normalizeKey(std::string key)
{
    std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    key.erase(std::remove(key.begin(), key.end(), '-'), key.end());
    return key;
}

} // namespace curcuma
