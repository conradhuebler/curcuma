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

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <vector>

namespace Elements {

static const std::vector<std::string> ElementAbbr = {
    "xXx", // leading index
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn"
    // La - Block not written down here
};

static const std::vector<std::string> ElementAbbr_Low = {
    "xXx", // leading index
    "h",
    "he",
    "li",
    "be",
    "b",
    "c",
    "n",
    "o",
    "f",
    "ne",
    "na",
    "mg",
    "al",
    "si",
    "p",
    "s",
    "cl",
    "ar",
    "k",
    "ca",
    "sc",
    "ti",
    "v",
    "cr",
    "mn",
    "fe",
    "co",
    "ni",
    "cu",
    "zn",
    "ga",
    "ge",
    "as",
    "se",
    "br",
    "kr",
    "rb",
    "sr",
    "y",
    "zr",
    "nb",
    "mo",
    "tc",
    "ru",
    "rh",
    "pd",
    "ag",
    "cd",
    "in",
    "sn",
    "sb",
    "te",
    "i",
    "xe",
    "cs",
    "ba",
    "la",
    "ce",
    "pr",
    "nd",
    "pm",
    "sm",
    "eu",
    "gd",
    "tb",
    "dy",
    "ho",
    "er",
    "tm",
    "yb",
    "lu",
    "hf",
    "ta",
    "w",
    "re",
    "os",
    "ir",
    "pt",
    "au",
    "hg",
    "tl",
    "pb",
    "bi",
    "po",
    "at",
    "rn"
    // La - Block not written down here
};

static const std::vector<double> AtomicMass = {
    -1, // leading index
    1.007940, // H
    4.002602, // He
    6.941000, // Li
    9.012182, // Be
    10.81100, // B
    12.01100, // C
    14.00674, // N
    15.99940, // O
    18.99840, // F
    20.17970, // Ne
    22.98977, // Na
    24.30500, // Mg
    26.98114, // Al
    28.08550, // S
    30.97376, // P
    32.06600, // S
    35.45270, // Cl
    39.94800, // Ar
    39.09830, // K
    40.07800, // Ca
    44.95591, // Sc
    47.88000, // Ti
    50.94150, // V
    51.99610, // Cr
    54.93805, // Mn
    55.84700, // Fe
    58.93320, // Co
    58.69340, // Ni
    63.54600, // Cu
    65.39000, // Zn
    69.72300, // Ga
    72.61000, // Ge
    74.92159, // As
    78.96000, // Se
    79.90400, // Br
    83.80000, // Kr
    85.46780, // Rb
    87.62000, // Sr
    88.90585, // Y
    91.22400, // Zr
    92.90638, // Nb
    95.94000, // Mo
    97.90721, // Tc
    101.0700, // Ru
    102.9055, // Rh
    106.4200, // Pd
    107.8682, // Ag
    112.4110, // Cd
    114.8280, // In
    118.7100, // Sn
    121.7570, // Sb
    127.6000, // Te
    126.9045, // I
    131.2900, // Xe
    132.9054, // Cs
    137.3270, // Ba
    138.9055, // La
    140.116, // Ce
    140.908, // Pr
    144.242, // Nd
    145.000, // Pm
    150.36, // Sm
    151.964, // Eu
    157.25, // Gd
    158.925, // Tb
    162.5, // Dy
    164.93, // Ho
    167.259, // Er
    168.934, // Tm
    173.054, // Yb
    174.967, // Lu
    178.49, // Hf
    180.948, // Ta
    183.84, // W
    186.207, // Re
    190.23, // Os
    192.217, // Ir
    195.084, // Pt
    196.967, // Au
    200.59, // Hg
    204.383, // Tl
    207.2, // Pb
    208.98, // Bi
    209, // Po
    210 // At
};

/*
 * Taken from Cramer and Truhlar J. Phys. Chem. A 2009, 113, 19, 5806-5812
 * 10.1021/jp8111556
 * d Block elements are not available in that work, there it is set to two ....
 */

static const std::vector<double> VanDerWaalsRadius = {
    -1, // leading index
    1.10, // H
    1.40, // He
    1.81, // Li
    1.53, // Be
    1.92, // B
    1.70, // C
    1.55, // N
    1.52, // O
    1.47, // F
    1.54, // Ne
    2.27, // Na
    1.73, // Mg
    1.84, // Al
    2.10, // S
    1.80, // P
    1.80, // S
    1.75, // Cl
    1.88, // Ar
    2.75, // K
    2.31, // Ca
    2, // Sc
    2, // Ti
    2, // V
    2, // Cr
    2, // Mn
    2, // Fe
    2, // Co
    2, // Ni
    2, // Cu
    2, // Zn
    1.87, // Ga
    2.11, // Ge
    1.85, // As
    1.90, // Se
    1.83, // Br
    2.02, // Kr
    3.03, // Rb
    2.49, // Sr
    2, // Y
    2, // Zr
    2, // Nb
    2, // Mo
    2, // Tc
    2, // Ru
    2, // Rh
    2, // Pd
    2, // Ag
    2, // Cd
    1.93, // In
    2.17, // Sn
    2.06, // Sb
    2.06, // Te
    1.98, // I
    2.16, // Xe
    3.43, // Cs
    2.68, // Ba
    2.5, // La
    2.48, // Ce
    2.47, // Pr
    2.45, // Nd
    2.43, // Pm
    2.42, // Sm
    2.40, // Eu
    2.38, // Gd
    2.37, // Tb
    2.35, // Dy
    2.33, // Ho
    2.32, // Er
    2.30, // Tm
    2.28, // Yb
    2.27, // Lu
    2.25, // Hf
    2.20, // Ta
    2.10, // W
    2.05, // Re
    2.00, // Os
    2.00, // Ir
    2.05, // Pt
    2.10, // Au
    2.05, // Hg
    2.20, // Tl
    2.30, // Pb
    2.30, // Bi
    2.00, // Po
    2.00 // At
};

/*
 * Taken from: Pyykkö and Asumi Chem. Eur. J. 2009, 15, 186 – 197
 * 10.1002/chem.200800987
 */
static const std::vector<double> CovalentRadius = {
    -1, // leading index
    0.32, // H
    0.46, // He
    1.33, // Li
    1.02, // Be
    0.85, // B
    0.75, // C
    0.71, // N
    0.63, // O
    0.64, // F
    0.67, // Ne
    1.55, // Na
    1.39, // Mg
    1.26, // Al
    1.16, // S
    1.11, // P
    1.03, // S
    0.99, // Cl
    0.96, // Ar
    1.96, // K
    1.71, // Ca
    1.48, // Sc
    1.36, // Ti
    1.34, // V
    1.22, // Cr
    1.19, // Mn
    1.16, // Fe
    1.11, // Co
    1.10, // Ni
    1.12, // Cu
    1.18, // Zn
    1.24, // Ga
    1.21, // Ge
    1.21, // As
    1.16, // Se
    1.14, // Br
    1.17, // Kr
    2.10, // Rb
    1.85, // Sr
    1.63, // Y
    1.54, // Zr
    1.47, // Nb
    1.38, // Mo
    1.28, // Tc
    1.25, // Ru
    1.25, // Rh
    1.20, // Pd
    1.28, // Ag
    1.36, // Cd
    1.42, // In
    1.40, // Sn
    1.40, // Sb
    1.36, // Te
    1.33, // I
    1.31, // Xe
    2.32, // Cs
    1.96, // Ba
    1.69, // La
    1.69, // Ce - taken from La
    1.69, // Pr - taken from La
    1.69, // Nd - taken from La
    1.69, // Pm - taken from La
    1.69, // Sm - taken from La
    1.69, // Eu - taken from La
    1.69, // Gd - taken from La
    1.69, // Tb - taken from La
    1.69, // Dy - taken from La
    1.69, // Ho - taken from La
    1.69, // Er - taken from La
    1.69, // Tm - taken from La
    1.69, // Yb - taken from La
    1.60, // Lu
    1.50, // Hf
    1.38, // Ta
    1.46, // W
    1.59, // Re
    1.28, // Os
    1.37, // Ir
    1.28, // Pt
    1.44, // Au
    1.49, // Hg
    1.48, // Tl
    1.47, // Pb
    1.46, // Bi
    1.46, // Po - taken from Bi
    1.46 // At - taken from Bi
};

// The old ones
/* for now, taken from here
 * https://de.wikipedia.org/wiki/Van-der-Waals-Radius */

/*
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
*/
static int String2Element(std::string string)
{
    transform(string.begin(), string.end(), string.begin(), ::tolower);
    int element = 0;

    for (int i = 0; i < ElementAbbr_Low.size(); ++i) {
        if (string.compare(ElementAbbr_Low[i]) == 0) {
            element = i;
            return element;
        }
    }
    return element;
}

}
