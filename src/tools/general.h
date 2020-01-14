/*
 * <Geometry tools for chemical structures.>
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

#include <Eigen/Dense>

#include "src/core/global.h"
#include "src/core/molecule.h"

namespace Tools {

int VectorDifference(const std::vector<int>& tmp_a, const std::vector<int>& tmp_b)
{
    int difference = 0;
    std::vector<int> a, b;
    if (tmp_a.size() < tmp_b.size()) {
        a = tmp_a;
        b = tmp_b;
    } else {
        b = tmp_a;
        a = tmp_b;
    }

    difference += b.size() - a.size();
    for (int i : a) {
        std::vector<int>::iterator it = std::find(b.begin(), b.end(), i);
        difference += it == b.end();
    }
    return difference;
}

}
