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

#include "gfnff.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <set>

namespace GFNFFAdvanced {

Vector calculateEEQCharges(int natoms, const std::vector<int>& atoms,
    const Matrix& geometry, const Vector& chi,
    const Vector& gam, const Vector& alp,
    const std::vector<double>& fragment_charges)
{
    // TODO: Implement full EEQ charge calculation
    // This requires:
    // 1. Build EEQ matrix A with hardness and coulomb terms
    // 2. Build constraint matrix for fragments
    // 3. Solve linear system A*q = b
    // 4. Apply damping and convergence criteria

    Vector charges = Vector::Zero(natoms);

    // Placeholder: Simple electronegativity-based charges
    for (int i = 0; i < natoms; ++i) {
        int z = atoms[i];
        if (z == 1)
            charges[i] = 0.1; // H
        else if (z == 6)
            charges[i] = -0.1; // C
        else if (z == 7)
            charges[i] = -0.2; // N
        else if (z == 8)
            charges[i] = -0.3; // O
        else if (z == 9)
            charges[i] = -0.4; // F
    }

    return charges;
}

std::vector<int> detectFragments(int natoms, const std::vector<int>& atoms,
    const Matrix& geometry, double bond_threshold)
{
    std::vector<int> fragments(natoms, -1);
    std::vector<bool> visited(natoms, false);
    int fragment_id = 0;

    // Build adjacency list
    std::vector<std::vector<int>> neighbors(natoms);
    for (int i = 0; i < natoms; ++i) {
        for (int j = i + 1; j < natoms; ++j) {
            Vector ri = geometry.row(i);
            Vector rj = geometry.row(j);
            double distance = (ri - rj).norm();

            // Simple bond detection (should use covalent radii)
            if (distance < bond_threshold) {
                neighbors[i].push_back(j);
                neighbors[j].push_back(i);
            }
        }
    }

    // DFS to find connected components
    for (int i = 0; i < natoms; ++i) {
        if (!visited[i]) {
            std::queue<int> queue;
            queue.push(i);
            visited[i] = true;

            while (!queue.empty()) {
                int atom = queue.front();
                queue.pop();
                fragments[atom] = fragment_id;

                for (int neighbor : neighbors[atom]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        queue.push(neighbor);
                    }
                }
            }
            fragment_id++;
        }
    }

    return fragments;
}

Matrix calculateTopologicalDistances(int natoms, const std::vector<int>& atoms,
    const std::vector<std::vector<int>>& neighbors,
    const Vector& covalent_radii)
{
    Matrix distances = Matrix::Constant(natoms, natoms, 1e6);

    // Initialize direct neighbors
    for (int i = 0; i < natoms; ++i) {
        distances(i, i) = 0.0;
        for (int j : neighbors[i]) {
            double rcov_i = covalent_radii[i];
            double rcov_j = covalent_radii[j];
            distances(i, j) = rcov_i + rcov_j;
        }
    }

    // Floyd-Warshall algorithm
    for (int k = 0; k < natoms; ++k) {
        for (int i = 0; i < natoms; ++i) {
            for (int j = 0; j < natoms; ++j) {
                if (distances(i, k) + distances(k, j) < distances(i, j)) {
                    distances(i, j) = distances(i, k) + distances(k, j);
                }
            }
        }
    }

    return distances;
}

std::vector<std::vector<int>> findRings(int natoms,
    const std::vector<std::vector<int>>& neighbors,
    int max_ring_size)
{
    std::vector<std::vector<int>> rings;

    // TODO: Implement ring detection algorithm
    // This is complex and should use sophisticated algorithms like
    // the smallest set of smallest rings (SSSR) algorithm

    return rings;
}

std::vector<int> determineHybridization(int natoms, const std::vector<int>& atoms,
    const Matrix& geometry,
    const std::vector<std::vector<int>>& neighbors)
{
    std::vector<int> hybridization(natoms, 3); // Default sp3

    for (int i = 0; i < natoms; ++i) {
        int z = atoms[i];
        int neighbor_count = neighbors[i].size();

        // Basic hybridization assignment
        if (neighbor_count <= 1) {
            hybridization[i] = 1; // sp
        } else if (neighbor_count == 2) {
            // Check if linear (sp) or bent (sp2)
            if (neighbors[i].size() == 2) {
                Vector ri = geometry.row(i);
                Vector rj = geometry.row(neighbors[i][0]);
                Vector rk = geometry.row(neighbors[i][1]);

                Vector v1 = rj - ri;
                Vector v2 = rk - ri;
                double angle = acos(v1.dot(v2) / (v1.norm() * v2.norm()));

                if (angle > 2.8) { // Close to 180°
                    hybridization[i] = 1; // sp
                } else {
                    hybridization[i] = 2; // sp2
                }
            } else {
                hybridization[i] = 2; // sp2
            }
        } else if (neighbor_count == 3) {
            hybridization[i] = 2; // sp2
        } else {
            hybridization[i] = 3; // sp3
        }
    }

    return hybridization;
}

std::vector<int> detectPiSystems(int natoms, const std::vector<int>& atoms,
    const std::vector<int>& hybridization,
    const std::vector<std::vector<int>>& neighbors)
{
    std::vector<int> pi_systems(natoms, 0);

    // TODO: Implement pi-system detection
    // This should identify connected sp/sp2 systems

    return pi_systems;
}

std::vector<std::vector<int>> buildNeighborList(int natoms, const std::vector<int>& atoms,
    const Matrix& geometry,
    const Vector& covalent_radii,
    double scale_factor)
{
    std::vector<std::vector<int>> neighbors(natoms);

    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < natoms; ++j) {
            if (i == j)
                continue;

            Vector ri = geometry.row(i);
            Vector rj = geometry.row(j);
            double distance = (ri - rj).norm();

            double rcov_i = covalent_radii[i];
            double rcov_j = covalent_radii[j];

            if (distance < scale_factor * (rcov_i + rcov_j)) {
                neighbors[i].push_back(j);
            }
        }
    }

    return neighbors;
}

double calculateEnvironmentBondForceConstant(int atom1, int atom2, double distance,
    double cn1, double cn2, int hyb1, int hyb2,
    int ring_size)
{
    // TODO: Implement environment-dependent force constants
    // This should include:
    // - Base force constant from bond parameters
    // - CN corrections
    // - Hybridization corrections
    // - Ring strain corrections

    return 1.0; // Placeholder
}

double calculateEnvironmentAngleForceConstant(int atom1, int atom2, int atom3,
    double angle, double cn2, int hyb2,
    int ring_size)
{
    // TODO: Implement environment-dependent angle force constants
    // Similar to bond force constants but for angles

    return 1.0; // Placeholder
}

} // namespace GFNFFAdvanced