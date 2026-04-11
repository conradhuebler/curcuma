/*
 * <Advanced GFN-FF Parameter Generation for Curcuma>
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
 */

#pragma once

#include "json.hpp"
#include "src/core/global.h"

using json = nlohmann::json;

namespace GFNFFAdvanced {

/**
 * @brief Calculate EEQ charges using iterative method
 * @param natoms Number of atoms
 * @param atoms Atomic numbers
 * @param geometry Geometry matrix
 * @param chi Electronegativity parameters
 * @param gam Hardness parameters
 * @param alp Polarizability parameters
 * @param fragment_charges Fragment charge constraints
 * @return EEQ charges
 */
Vector calculateEEQCharges(int natoms, const std::vector<int>& atoms,
    const Matrix& geometry, const Vector& chi,
    const Vector& gam, const Vector& alp,
    const std::vector<double>& fragment_charges = {});

/**
 * @brief Detect molecular fragments using connectivity
 * @param natoms Number of atoms
 * @param atoms Atomic numbers
 * @param geometry Geometry matrix
 * @param bond_threshold Bond detection threshold
 * @return Fragment assignment for each atom
 */
std::vector<int> detectFragments(int natoms, const std::vector<int>& atoms,
    const Matrix& geometry, double bond_threshold = 1.3);

/**
 * @brief Calculate topological distances using Floyd-Warshall algorithm
 * @param natoms Number of atoms
 * @param atoms Atomic numbers
 * @param neighbors Neighbor list
 * @param covalent_radii Covalent radii
 * @return Topological distance matrix
 */
Matrix calculateTopologicalDistances(int natoms, const std::vector<int>& atoms,
    const std::vector<std::vector<int>>& neighbors,
    const Vector& covalent_radii);

/**
 * @brief Find rings in molecular graph using depth-first search
 * @param natoms Number of atoms
 * @param neighbors Neighbor list
 * @param max_ring_size Maximum ring size to detect
 * @return Ring information for each atom
 */
std::vector<std::vector<int>> findRings(int natoms,
    const std::vector<std::vector<int>>& neighbors,
    int max_ring_size = 10);

/**
 * @brief Determine hybridization using geometry analysis
 * @param natoms Number of atoms
 * @param atoms Atomic numbers
 * @param geometry Geometry matrix
 * @param neighbors Neighbor list
 * @return Hybridization states
 */
std::vector<int> determineHybridization(int natoms, const std::vector<int>& atoms,
    const Matrix& geometry,
    const std::vector<std::vector<int>>& neighbors);

/**
 * @brief Detect pi-conjugated systems
 * @param natoms Number of atoms
 * @param atoms Atomic numbers
 * @param hybridization Hybridization states
 * @param neighbors Neighbor list
 * @return Pi-system assignment
 */
std::vector<int> detectPiSystems(int natoms, const std::vector<int>& atoms,
    const std::vector<int>& hybridization,
    const std::vector<std::vector<int>>& neighbors);

/**
 * @brief Build neighbor list from geometry
 * @param natoms Number of atoms
 * @param atoms Atomic numbers
 * @param geometry Geometry matrix
 * @param covalent_radii Covalent radii
 * @param scale_factor Scaling factor for bond detection
 * @return Neighbor list
 */
std::vector<std::vector<int>> buildNeighborList(int natoms, const std::vector<int>& atoms,
    const Matrix& geometry,
    const Vector& covalent_radii,
    double scale_factor = 1.3);

/**
 * @brief Generate environment-dependent force constants
 * @param atom1 First atom type
 * @param atom2 Second atom type
 * @param distance Bond distance
 * @param cn1 Coordination number of atom 1
 * @param cn2 Coordination number of atom 2
 * @param hyb1 Hybridization of atom 1
 * @param hyb2 Hybridization of atom 2
 * @param ring_size Ring size (0 = no ring)
 * @return Force constant
 */
double calculateEnvironmentBondForceConstant(int atom1, int atom2, double distance,
    double cn1, double cn2, int hyb1, int hyb2,
    int ring_size = 0);

/**
 * @brief Generate environment-dependent angle force constants
 * @param atom1 First atom type
 * @param atom2 Center atom type
 * @param atom3 Third atom type
 * @param angle Current angle
 * @param cn2 Coordination number of center atom
 * @param hyb2 Hybridization of center atom
 * @param ring_size Ring size (0 = no ring)
 * @return Force constant
 */
double calculateEnvironmentAngleForceConstant(int atom1, int atom2, int atom3,
    double angle, double cn2, int hyb2,
    int ring_size = 0);

/**
 * @brief TODO: Implementation roadmap for advanced features
 *
 * Priority 1 (High):
 * - EEQ charge calculation with matrix solver
 * - Coordination number calculation with derivatives
 * - Basic hybridization detection
 *
 * Priority 2 (Medium):
 * - Ring detection algorithm
 * - Environment-dependent parameters
 * - Fragment detection and constraints
 *
 * Priority 3 (Low):
 * - Pi-system detection
 * - Hydrogen bond detection
 * - Metallic character corrections
 * - Torsion and out-of-plane parameters
 */

} // namespace GFNFFAdvanced