/*
 * <COSMO van-der-Waals radii for the native CPCM solvation model>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (June 2026): COSMO-optimised van-der-Waals radii used to
 * build the molecular cavity for the native CPCM (ddCOSMO) solvation model.
 * Transcribed verbatim from tblite
 * (external/tblite/src/tblite/solvation/data.f90, `vdw_rad_cosmo(1:94)`).
 */

#pragma once

#include <cstddef>

namespace Curcuma {
namespace Solvation {

/// Default radius for elements without an optimised COSMO radius (Angstrom).
inline constexpr double kCosmoStub = 2.223;

/// CODATA-2018 Angstrom -> Bohr (matches Curcuma's ANGSTROM_TO_BOHR and the
/// AA_TO_AU used throughout the native xTB path, so the cavity is built in the
/// same length unit as the rest of the SCF).
inline constexpr double kCosmoAngstromToBohr = 1.0 / 0.529177210903;

/**
 * @brief COSMO-optimised van-der-Waals radii in Angstrom for Z = 1..94.
 *
 * Index i holds the radius of element Z = i+1. Elements without an optimised
 * value carry @ref kCosmoStub. Source: tblite data.f90:106 (`vdw_rad_cosmo`),
 * which itself stores them in Bohr via `aatoau * [...]`; here they are kept in
 * Angstrom and converted on access by @ref getVdwRadCosmo.
 */
inline constexpr double kVdwRadCosmoAngstrom[94] = {
    1.3000, 1.6380, 1.5700, 1.0530,        // H  He Li Be
    2.0480, 2.0000, 1.8300, 1.7200,        // B  C  N  O
    1.7200, 1.8018, 1.8000, 1.6380,        // F  Ne Na Mg
    2.1530, 2.2000, 2.1060, 2.1600,        // Al Si P  S
    2.0500, 2.2000, 2.2230, kCosmoStub,    // Cl Ar K  Ca
    kCosmoStub, 2.2930, kCosmoStub, kCosmoStub, // Sc Ti V  Cr
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Mn Fe Co Ni
    kCosmoStub, 1.6260, kCosmoStub, 2.7000,     // Cu Zn Ga Ge
    2.3500, 2.2000, 2.1600, 2.3630,        // As Se Br Kr
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Rb Sr Y  Zr
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Nb Mo Tc Ru
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Rh Pd Ag Cd
    2.2580, 2.5500, 2.4100, 2.4100,        // In Sn Sb Te
    2.3200, 2.5270, kCosmoStub, kCosmoStub,     // I  Xe Cs Ba
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // La Ce Pr Nd
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Pm Sm Eu Gd
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Tb Dy Ho Er
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Tm Yb Lu Hf
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Ta W  Re Os
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Ir Pt Au Hg
    kCosmoStub, 2.3600, 2.4220, 2.3050,    // Tl Pb Bi Po
    2.3630, 2.5740, kCosmoStub, kCosmoStub,     // At Rn Fr Ra
    kCosmoStub, kCosmoStub, kCosmoStub, kCosmoStub, // Ac Th Pa U
    kCosmoStub, kCosmoStub                       // Np Pu
};

/**
 * @brief COSMO van-der-Waals radius in Bohr for atomic number Z.
 * @param z Atomic number (1-based). Out-of-range Z returns the stub radius.
 */
inline double getVdwRadCosmo(int z)
{
    const double r_ang = (z >= 1 && z <= 94) ? kVdwRadCosmoAngstrom[z - 1] : kCosmoStub;
    return r_ang * kCosmoAngstromToBohr;
}

} // namespace Solvation
} // namespace Curcuma
