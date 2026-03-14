/*
 * <Periodic Table Utility for Curcuma>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Periodic Table Implementation
 *
 * Claude Generated (Session 5, December 2025)
 * Translated from Fortran gfnff_param.f90:408-424
 */

#include "periodic_table.h"

namespace PeriodicTable {

/**
 * Metal classification (0/1/2)
 * Index: atomic_number - 1 (C++ 0-indexed)
 * Value: 0 = non-metal, 1 = main group metal, 2 = transition metal
 *
 * Fortran gfnff_param.f90:408-415
 */
static constexpr int METAL_TABLE[86] = {
    0,  0,                                                   // H(1), He(2)
    1,  1,  0,  0,  0,  0,  0,  0,                           // Li-Ne (3-10)
    1,  1,  1,  0,  0,  0,  0,  0,                           // Na-Ar (11-18)
    1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  0,  0,  0,  0,  0,  // K-Kr (19-36)
    1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  1,  0,  0,  0,  0,  // Rb-Xe (37-54)
    1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,      // Cs-Lu (55-71)
    2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  1,  1,  1,  0,  0              // Hf-Rn (72-86)
};

/**
 * Periodic table group classification
 * Index: atomic_number - 1 (C++ 0-indexed)
 * Value: 1-8 for main groups, negative (-3 to -12) for transition metals
 *
 * Fortran gfnff_param.f90:417-424
 */
static constexpr int GROUP_TABLE[86] = {
    1,  8,                                                   // H(1), He(2)
    1,  2,  3,  4,  5,  6,  7,  8,                           // Li-Ne (3-10)
    1,  2,  3,  4,  5,  6,  7,  8,                           // Na-Ar (11-18)
    1,  2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,  3,  4,  5,  6,  7,  8,  // K-Kr (19-36)
    1,  2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,  3,  4,  5,  6,  7,  8,  // Rb-Xe (37-54)
    1,  2, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,      // Cs-Lu (55-71)
    -4, -5, -6, -7, -8, -9,-10,-11,-12,  3,  4,  5,  6,  7,  8              // Hf-Rn (72-86)
};

// ============================================================================
// PUBLIC INTERFACE
// ============================================================================

int getGroup(int atomic_number) {
    if (atomic_number < 1 || atomic_number > 86) {
        return 0;  // Invalid
    }
    return GROUP_TABLE[atomic_number - 1];
}

int getMetalType(int atomic_number) {
    if (atomic_number < 1 || atomic_number > 86) {
        return 0;  // Invalid - treat as non-metal
    }
    return METAL_TABLE[atomic_number - 1];
}

bool isHalogen(int atomic_number) {
    int group = getGroup(atomic_number);
    // Group 7 for main group: F(9), Cl(17), Br(35), I(53), At(85)
    return (group == 7);
}

bool isChalcogen(int atomic_number) {
    int group = getGroup(atomic_number);
    // Group 6 for main group: O(8), S(16), Se(34), Te(52), Po(84)
    return (group == 6);
}

bool isTransitionMetal(int atomic_number) {
    int metal = getMetalType(atomic_number);
    return (metal == 2);
}

bool isMainGroupMetal(int atomic_number) {
    int metal = getMetalType(atomic_number);
    return (metal == 1);
}

} // namespace PeriodicTable
