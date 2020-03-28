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
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba" };

static const std::vector<std::string> ElementAbbr_Low = { "h", "he",
    "li", "be", "b", "c", "n", "o", "f", "ne",
    "na", "mg", "al", "si", "p", "s", "cl", "ar",
    "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr",
    "rb", "sr", "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", "sb", "te", "i", "xe",
    "cs", "ba" };

static const std::vector<double> AtomicMass = { 1.0079, 4.00326,
    6.941, 9.0122, 10.811, 12.011, 14.007, 15.999, 18.998, 20.180 };

/* for now, taken from here
 * https://de.wikipedia.org/wiki/Van-der-Waals-Radius */
static const std::vector<double> VanDerWaalsRadius = { 1.10, 1.40,
    1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
    2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88,
    2.75, 2.31, 2, 2, 2, 2, 2, 2, 2, 1.62, 1.40, 1.39, 1.87, 2.11, 1.85, 1.90, 1.85, 2.02,
    3.03, 2.49, 2, 2, 2, 2, 2, 2, 2, 1.63, 1.72, 1.58, 1.93, 2.17, 2.06, 2.06, 1.98, 2.16 };

static const std::vector<double> CovalentRadius = { 0.32, 0.32,
    1.34, 0.9, 0.82, 0.77, 0.71, 0.73, 0.71, 0.69,
    1.54, 1.30, 1.18, 1.11, 1.06, 1.02, 0.99, 0.97,
    1.96, 1.74, 1.44, 1.36, 1.25, 1.27, 1.39, 1.25, 1.26, 1.21, 1.38, 1.31, 1.26, 1.22, 1.21, 1.16, 1.14, 1.10,
    2.11, 1.92, 1.62, 1.48, 1.37, 1.45, 1.31, 1.26, 1.35, 1.31, 1.53, 1.48, 1.44, 1.41, 1.38, 1.35, 1.33, 1.30 };

static int String2Element(const std::string& string)
{
    int element = 0;

    for (int i = 0; i < ElementAbbr.size(); ++i) {
        if (string.compare(ElementAbbr[i]) == 0) {
            element = i;
            return element;
        }
    }
    return element;
}

}
