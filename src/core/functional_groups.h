/*
 * <Functional Group Detection for Curcuma>
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
 * Functional Group Detection Framework
 * ====================================
 * Extensible pattern matching system for chemical groups.
 * Identifies functional groups based on topology, hybridization, and pi-system membership.
 *
 * Claude Generated (Session 5, December 2025)
 * Translated from Fortran gfnff_ini2.f90:1461-1511
 */

#pragma once

#include <vector>
#include <memory>

/**
 * @brief Enumeration for functional group types
 *
 * Classification of recognized functional groups.
 * Used for applying group-specific dxi corrections in EEQ calculation.
 */
enum class FunctionalGroupType {
    None = 0,      // No functional group
    Amide,         // Amide group (N-C=O): for dxi corrections
    AmideH,        // Hydrogen bonded to amide nitrogen: chi -= 0.02
    Nitro,         // Nitro group (N-O): oxygen in pi-bonded nitro
    Water,         // Water molecule (H2O): dxi -= 0.02
    Carbonyl,      // Carbonyl group (C=O): carbon with pi-bonded oxygen
    FreeCarbon     // Free CO molecule: C and O both with CN=1
};

/**
 * @brief Functional Group Detector Class
 *
 * Extensible pattern matching system for identifying chemical functional groups.
 * Uses topology information (neighbors, hybridization, pi-system membership) to
 * classify atoms into functional group categories.
 *
 * Reference: Fortran gfnff_ini2.f90:1461-1511
 *
 * Claude Generated (Session 5, December 2025)
 */
class FunctionalGroupDetector {
public:
    /**
     * @brief Constructor
     *
     * @param atom_count Number of atoms in molecule
     * @param atomic_numbers Element atomic numbers [atom_count]
     * @param neighbor_lists Full neighbor connectivity [atom_count]
     * @param hybridization Hybridization states [atom_count] (1=sp, 2=sp², 3=sp³, 5=sp³d²)
     * @param pi_fragments Pi-system membership [atom_count] (0=not in pi, >0=pi fragment ID)
     */
    FunctionalGroupDetector(
        int atom_count,
        const std::vector<int>& atomic_numbers,
        const std::vector<std::vector<int>>& neighbor_lists,
        const std::vector<int>& hybridization,
        const std::vector<int>& pi_fragments
    );

    /**
     * @brief Detect all functional groups in the molecule
     *
     * @return Vector of functional group classifications [atom_count]
     *         Each element corresponds to the functional group type for that atom
     *
     * Pattern detection order:
     * 1. Amide nitrogen and hydrogen
     * 2. Nitro group oxygen
     * 3. Water molecule
     * 4. Free carbon monoxide
     * 5. Carbonyl oxygen (C=O)
     */
    std::vector<FunctionalGroupType> detectAllGroups() const;

    /**
     * @brief Check if atom is an amide nitrogen
     *
     * Pattern: N(sp³, in pi-system) bonded to C(pi) which has C=O
     *
     * Criteria:
     * - Element: Nitrogen (Z=7)
     * - Hybridization: sp³
     * - Pi-system: Must be in pi-system
     * - Neighbors: One pi-bonded carbon
     * - That carbon: Has C=O (carbonyl oxygen with pi, CN=1)
     *
     * Reference: Fortran gfnff_ini2.f90:1476-1503
     */
    bool isAmideNitrogen(int atom) const;

    /**
     * @brief Check if atom is hydrogen bonded to amide nitrogen
     *
     * Pattern: H bonded to amide nitrogen
     *
     * Criteria:
     * - Element: Hydrogen (Z=1)
     * - Coordination: CN=1 (single bond)
     * - Neighbor: Must be amide nitrogen
     * - Additional: Nitrogen has exactly 1 sp³ carbon neighbor
     *
     * Reference: Fortran gfnff_ini2.f90:1461-1475
     */
    bool isAmideHydrogen(int atom) const;

    /**
     * @brief Check if atom is oxygen in nitro group
     *
     * Pattern: O in N-O (pi-bonded)
     *
     * Criteria:
     * - Element: Oxygen (Z=8)
     * - Coordination: CN=1
     * - Pi-system: In pi-system
     * - Neighbor: Pi-bonded nitrogen (also in pi-system)
     */
    bool isNitroOxygen(int atom) const;

    /**
     * @brief Check if atom is oxygen in water molecule
     *
     * Pattern: H-O-H
     *
     * Criteria:
     * - Element: Oxygen (Z=8)
     * - Coordination: CN=2
     * - Both neighbors: Hydrogen atoms
     */
    bool isWaterOxygen(int atom) const;

    /**
     * @brief Check if atom is oxygen in carbonyl group
     *
     * Pattern: C=O (pi-bonded)
     *
     * Criteria:
     * - Element: Oxygen (Z=8)
     * - Coordination: CN=1
     * - Pi-system: In pi-system
     * - Neighbor: Pi-bonded carbon
     */
    bool isCarbonylOxygen(int atom) const;

    /**
     * @brief Check if atom is carbon in free CO molecule
     *
     * Pattern: C with CN=1 bonded to O with CN=1
     *
     * Criteria:
     * - Element: Carbon (Z=6)
     * - Coordination: CN=1
     * - Neighbor: Oxygen with CN=1
     */
    bool isFreeCarbon(int atom) const;

    /**
     * @brief Check if carbon atom has pi-bonded carbonyl oxygen
     *
     * Helper method for amide detection.
     * Checks if this carbon has exactly one neighbor that is:
     * - Oxygen (Z=8)
     * - Pi-bonded (in pi-system)
     * - CN=1 (terminal)
     *
     * @param carbon_atom Index of carbon atom to check
     * @return true if carbon has C=O
     */
    bool hasCarbonylOxygen(int carbon_atom) const;

private:
    int m_natoms;
    const std::vector<int>& m_atomic_numbers;
    const std::vector<std::vector<int>>& m_neighbor_lists;
    const std::vector<int>& m_hybridization;
    const std::vector<int>& m_pi_fragments;

    // Helper to get coordination number
    int getCoordinationNumber(int atom) const;
};
