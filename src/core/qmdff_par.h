/*
 * <QMDFF implementation for Cucuma. >
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

/*

 */

#pragma once

#include "json.hpp"
#include <Eigen/Dense>
#include <vector>
using json = nlohmann::json;

/* Parameter taken from publication */
const double kEN = -0.164;
const double ka2 = 0.221;
const double ka13 = 2.81;
const double kb13 = 0.53;
const double k13r = 0.7;
const double kdmp = 0.11;
const double kovlp = 0.5;
const double a1 = 0.45;
const double a2 = 4.0;
const double s8 = 2.7;
const double epsilon_ES = 0.85;
const double epsilon_disp_res = 0.5;
const double beta_rep = 16.5;
const double khbnd_N = 0.8;
const double khbnd_O = 0.3;
const double khbnd_F = 0.1;
const double khbnd_X = 2.0;
const double kXCl = 0.3;
const double kXBr = 0.6;
const double kXI = 0.8;
const double kXAt = 1;
const double kq1_hbnd = 10;
const double kq1_xbnd = -6.5;
const double kq2_hbnd = 5;
const double kq2_xbnd = 1;
const double kq = 1.15;

struct QMDFFBond {
    int a, b, distance;
    double reAB, kAB, exponA;
};

struct QMDFFAngle {
    int a, b, c;
    double thetae, kabc, reAB, reAC;
};

inline double ka(int element)
{
    if (element == 1 || element == 2) // H, He
        return 1.755;
    else if (element == 5) // B
        return 2.287;
    else if (element == 6) // C
        return 2.463;
    else if (element == 7) // N
        return 2.559;
    else if (element == 8) // O
        return 2.579;
    else if (element == 9 || element == 10) // F, Ne
        return 2.2;
    else if (element == 13) // Al
        return 2.559 + ka2;
    else if (element == 14) // Si
        return 2.463 + ka2;
    else if (element == 15) // P
        return 2.559 + ka2;
    else if (element == 16) // S
        return 2.579 + ka2;
    else if (element == 17) // Cl
        return 2.2 + ka2;
    else
        return 1; // todo
}

inline double AllenEN(int element)
{
    if (element == 1 || element == 2) // H, He
        return 2.30;
    else if (element == 5) // B
        return 2.051;
    else if (element == 6) // C
        return 2.544;
    else if (element == 7) // N
        return 3.066;
    else if (element == 8) // O
        return 3.610;
    else if (element == 9) // F, Ne
        return 4.193;
    else if (element == 10)
        return 4.787;
    else if (element == 13) // Al
        return 1.613;
    else if (element == 14) // Si
        return 1.916;
    else if (element == 15) // P
        return 2.253;
    else if (element == 16) // S
        return 2.589;
    else if (element == 17) // Cl
        return 2.869;
    else
        return 1; // todo
}

inline double kZ(int element)
{
    if (element <= 2) // row 1
        return 2.35;
    else if (element <= 10) // row 2
        return 0.95;
    else if (element <= 18) // row 3
        return 0.75;
    else if (element <= 36) // row 4
        return 0.65;
    else // row 5
        return 0.6;
}

const json QMDFFParameterJson{
    { "bond_scaling", 1 },
    { "angle_scaling", 1 },
    { "dihedral_scaling", 1 },
    { "inversion_scaling", 1 },
    { "vdw_scaling", 1 },
    { "rep_scaling", 1 },
    { "coulomb_scaling", 1 },
    { "bond_force", 664.12 },
    { "angle_force", 664.12 },
    { "differential", 1e-7 },
    { "h4_scaling", 0 },
    { "hh_scaling", 0 },
    { "h4_oh_o", 2.32 },
    { "h4_oh_n", 3.10 },
    { "h4_nh_o", 1.07 },
    { "h4_nh_n", 2.01 },
    { "h4_wh_o", 0.42 },
    { "h4_nh4", 3.61 },
    { "h4_coo", 1.41 },
    { "hh_rep_k", 0.42 },
    { "hh_rep_e", 12.7 },
    { "hh_rep_r0", 2.3 },
    { "d4", 0 },
    { "d3", 0 },
    { "d_s6", 1.00 },
    { "d_s8", 1.20065498 },
    { "d_s10", 0 },
    { "d_s9", 1 },
    { "d_a1", 0.40085597 },
    { "d_a2", 5.02928789 },
    { "d_alp", 16 },
    { "d_func", "pbe0" },
    { "d_atm", true },
    { "param_file", "none" },
    { "ff_file", "none" },
    { "writeparam", "none" },
    { "writeff", "none" },
    { "verbose", false },
    { "rings", false },
    { "threads", 1 },
    { "gradient", 0 },
    { "variable", true },
    { "const", true }
};
