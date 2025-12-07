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
 * Functional Group Detection Implementation
 *
 * Claude Generated (Session 5, December 2025)
 * Translated from Fortran gfnff_ini2.f90:1461-1511
 */

#include "functional_groups.h"
#include <algorithm>

// ============================================================================
// Constructor
// ============================================================================

FunctionalGroupDetector::FunctionalGroupDetector(
    int atom_count,
    const std::vector<int>& atomic_numbers,
    const std::vector<std::vector<int>>& neighbor_lists,
    const std::vector<int>& hybridization,
    const std::vector<int>& pi_fragments)
    : m_natoms(atom_count),
      m_atomic_numbers(atomic_numbers),
      m_neighbor_lists(neighbor_lists),
      m_hybridization(hybridization),
      m_pi_fragments(pi_fragments)
{
}

// ============================================================================
// Main Detection Method
// ============================================================================

std::vector<FunctionalGroupType> FunctionalGroupDetector::detectAllGroups() const
{
    std::vector<FunctionalGroupType> groups(m_natoms, FunctionalGroupType::None);

    for (int i = 0; i < m_natoms; ++i) {
        // Priority order of detection (first match wins)

        // 1. Amide hydrogen (must check before amide nitrogen)
        if (isAmideHydrogen(i)) {
            groups[i] = FunctionalGroupType::AmideH;
            continue;
        }

        // 2. Amide nitrogen
        if (isAmideNitrogen(i)) {
            groups[i] = FunctionalGroupType::Amide;
            continue;
        }

        // 3. Nitro group oxygen
        if (isNitroOxygen(i)) {
            groups[i] = FunctionalGroupType::Nitro;
            continue;
        }

        // 4. Water molecule
        if (isWaterOxygen(i)) {
            groups[i] = FunctionalGroupType::Water;
            continue;
        }

        // 5. Free carbon monoxide
        if (isFreeCarbon(i)) {
            groups[i] = FunctionalGroupType::FreeCarbon;
            continue;
        }

        // 6. Carbonyl oxygen
        if (isCarbonylOxygen(i)) {
            groups[i] = FunctionalGroupType::Carbonyl;
            continue;
        }
    }

    return groups;
}

// ============================================================================
// Individual Detection Methods
// ============================================================================

int FunctionalGroupDetector::getCoordinationNumber(int atom) const
{
    if (atom < 0 || atom >= m_natoms) {
        return 0;
    }
    return m_neighbor_lists[atom].size();
}

// ============================================================================
// Amide Detection (Fortran gfnff_ini2.f90:1461-1511)
// ============================================================================

bool FunctionalGroupDetector::isAmideNitrogen(int atom) const
{
    // Pattern: N(sp³, in pi-system) bonded to C(pi) which has C=O

    // Check element and hybridization
    if (m_atomic_numbers[atom] != 7) return false;     // Not nitrogen
    if (m_hybridization[atom] != 3) return false;      // Not sp³
    if (m_pi_fragments[atom] == 0) return false;       // Not in pi-system

    int nn = getCoordinationNumber(atom);
    if (nn == 0) return false;

    // Count pi-bonded carbons and check for carbonyl
    int pi_carbons = 0;
    int carbonyl_oxygens = 0;

    for (int neighbor : m_neighbor_lists[atom]) {
        if (m_atomic_numbers[neighbor] == 6 && m_pi_fragments[neighbor] != 0) {
            pi_carbons++;
            // Check if this carbon has C=O
            if (hasCarbonylOxygen(neighbor)) {
                carbonyl_oxygens++;
            }
        }
    }

    // Amide nitrogen: exactly 1 pi-bonded carbon with C=O
    return (pi_carbons == 1 && carbonyl_oxygens == 1);
}

bool FunctionalGroupDetector::isAmideHydrogen(int atom) const
{
    // Pattern: H bonded to amide nitrogen

    // Check element and coordination
    if (m_atomic_numbers[atom] != 1) return false;     // Not hydrogen
    if (getCoordinationNumber(atom) != 1) return false; // Must have exactly 1 bond

    // Get the nitrogen neighbor
    int nitrogen = m_neighbor_lists[atom][0];

    // Check if neighbor is amide nitrogen
    if (!isAmideNitrogen(nitrogen)) {
        return false;
    }

    // Additional check: nitrogen must have exactly 1 sp³ carbon neighbor
    int sp3_carbons = 0;
    for (int neighbor : m_neighbor_lists[nitrogen]) {
        if (m_atomic_numbers[neighbor] == 6 && m_hybridization[neighbor] == 3) {
            sp3_carbons++;
        }
    }

    return (sp3_carbons == 1);
}

bool FunctionalGroupDetector::hasCarbonylOxygen(int carbon_atom) const
{
    // Check if carbon has exactly one pi-bonded oxygen with CN=1

    int carbonyl_count = 0;

    for (int neighbor : m_neighbor_lists[carbon_atom]) {
        // Check if neighbor is oxygen
        if (m_atomic_numbers[neighbor] != 8) continue;

        // Check if oxygen is pi-bonded
        if (m_pi_fragments[neighbor] == 0) continue;

        // Check if oxygen has CN=1 (terminal)
        if (getCoordinationNumber(neighbor) != 1) continue;

        carbonyl_count++;
    }

    return (carbonyl_count >= 1);  // At least one C=O
}

// ============================================================================
// Nitro Group Detection
// ============================================================================

bool FunctionalGroupDetector::isNitroOxygen(int atom) const
{
    // Pattern: O in N-O (pi-bonded)

    // Check element, coordination, and pi-system
    if (m_atomic_numbers[atom] != 8) return false;      // Not oxygen
    if (getCoordinationNumber(atom) != 1) return false;  // CN=1 only
    if (m_pi_fragments[atom] == 0) return false;        // Must be in pi-system

    // Check neighbor: must be pi-bonded nitrogen also in pi-system
    int neighbor = m_neighbor_lists[atom][0];
    if (m_atomic_numbers[neighbor] != 7) return false;           // Not nitrogen
    if (m_pi_fragments[neighbor] == 0) return false;             // Not in pi-system

    return true;
}

// ============================================================================
// Water Detection
// ============================================================================

bool FunctionalGroupDetector::isWaterOxygen(int atom) const
{
    // Pattern: H-O-H

    // Check element
    if (m_atomic_numbers[atom] != 8) return false;

    // Check coordination: exactly 2
    if (getCoordinationNumber(atom) != 2) return false;

    // Both neighbors must be hydrogen
    for (int neighbor : m_neighbor_lists[atom]) {
        if (m_atomic_numbers[neighbor] != 1) return false;  // Not hydrogen
    }

    return true;
}

// ============================================================================
// Carbonyl Detection
// ============================================================================

bool FunctionalGroupDetector::isCarbonylOxygen(int atom) const
{
    // Pattern: C=O (pi-bonded)

    // Check element, coordination, and pi-system
    if (m_atomic_numbers[atom] != 8) return false;      // Not oxygen
    if (getCoordinationNumber(atom) != 1) return false;  // CN=1 only
    if (m_pi_fragments[atom] == 0) return false;        // Must be in pi-system

    // Check neighbor: must be pi-bonded carbon
    int neighbor = m_neighbor_lists[atom][0];
    if (m_atomic_numbers[neighbor] != 6) return false;           // Not carbon
    if (m_pi_fragments[neighbor] == 0) return false;             // Not in pi-system

    return true;
}

// ============================================================================
// Free Carbon Monoxide Detection
// ============================================================================

bool FunctionalGroupDetector::isFreeCarbon(int atom) const
{
    // Pattern: C with CN=1 bonded to O with CN=1

    // Check element and coordination
    if (m_atomic_numbers[atom] != 6) return false;      // Not carbon
    if (getCoordinationNumber(atom) != 1) return false;  // CN=1 only

    // Check neighbor: must be oxygen with CN=1
    int neighbor = m_neighbor_lists[atom][0];
    if (m_atomic_numbers[neighbor] != 8) return false;           // Not oxygen
    if (getCoordinationNumber(neighbor) != 1) return false;       // Oxygen CN=1

    return true;
}
