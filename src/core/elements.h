/*
 * <Some globale definition for chemical structures.>
 * Copyright (C) 2019 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <array>
#include <set>
#include <string>
#include <vector>

namespace Elements {

static const std::vector<std::string> ElementAbbr = { "xXx", "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar" };

static const std::vector<std::string> ElementAbbr_Low = { "h", "he",
    "li", "be", "b", "c", "n", "o", "f", "ne",
    "na", "mg", "al", "si", "p", "s", "cl", "ar" };

static const std::vector<double> AtomicMass = { 1.008, 4.003 };

/* for now, taken from here
 * https://de.wikipedia.org/wiki/Van-der-Waals-Radius */
static const std::vector<double> VanDerWaalsRadius = { 1.10, 1.40,
    1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
    2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88 };

static const std::vector<double> CovalentRadius = { 0.32, 0.32,
    1.34, 0.9, 0.82, 0.77, 0.71, 0.73, 0.71, 0.69,
    1.54, 1.30, 1.18, 1.11, 1.06, 1.02, 0.99, 0.97 };

static int String2Element(const std::string& string)
{
    int element = 0;

    for (int i = 0; i < ElementAbbr.size(); ++i) {
        if (string.find(ElementAbbr[i]) != std::string::npos) {
            element = i;
            return element;
        }
    }
    return element;
}

}
