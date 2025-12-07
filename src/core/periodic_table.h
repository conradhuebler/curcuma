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
 * Periodic Table Utilities
 * ========================
 * Centralized classification system for chemical elements.
 * Provides group and metal type information for all 86 elements.
 *
 * Claude Generated (Session 5, December 2025)
 * Translated from Fortran gfnff_param.f90:408-424
 */

#pragma once

namespace PeriodicTable {

/**
 * @brief Get the periodic table group for an element
 * @param atomic_number Atomic number (1-86)
 * @return Group number:
 *         - 1-8: Main group elements
 *         - 0: Invalid/unknown
 *         - Negative (-3 to -12): Transition metals and lanthanides
 *
 * Reference: Fortran gfnff_param.f90:417-424
 */
int getGroup(int atomic_number);

/**
 * @brief Get the metal classification for an element
 * @param atomic_number Atomic number (1-86)
 * @return Metal type:
 *         - 0: Non-metal (H, C, N, O, halogens, etc.)
 *         - 1: Main group metal (Li, Be, Al, Ga, In, Sn, Pb, Bi)
 *         - 2: Transition metal or lanthanide
 *
 * Reference: Fortran gfnff_param.f90:408-415
 */
int getMetalType(int atomic_number);

/**
 * @brief Check if an element is a halogen (Group 7)
 * @param atomic_number Atomic number (1-86)
 * @return true if element is F, Cl, Br, I, or At
 */
bool isHalogen(int atomic_number);

/**
 * @brief Check if an element is a chalcogen (Group 6)
 * @param atomic_number Atomic number (1-86)
 * @return true if element is O, S, Se, Te, or Po
 */
bool isChalcogen(int atomic_number);

/**
 * @brief Check if an element is a transition metal
 * @param atomic_number Atomic number (1-86)
 * @return true if element is in d-block (Sc-Zn, Y-Cd, La-Hg, Ac-Cn)
 */
bool isTransitionMetal(int atomic_number);

/**
 * @brief Check if an element is a main group metal
 * @param atomic_number Atomic number (1-86)
 * @return true if metal type is 1 (main group metal)
 */
bool isMainGroupMetal(int atomic_number);

} // namespace PeriodicTable
