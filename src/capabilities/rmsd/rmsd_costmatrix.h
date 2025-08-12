/*
 * <Cost Matrix Calculator for RMSD calculations>
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
 * Claude Generated - Refactored from rmsd.cpp for better modularity
 */

#pragma once

#include "src/core/molecule.h"
#include "src/tools/geometry.h"
#include <utility>
#include <vector>

/**
 * @brief Configuration for cost matrix calculation
 * Claude Generated - Extracted from RMSDDriver parameters
 */
struct CostMatrixConfig {
    int cost_function = 1; // Cost function type (1-6)
    double scaling = 1.5; // Scaling factor
    bool check_connectivity = false; // Perform connectivity checks

    CostMatrixConfig() = default;
    CostMatrixConfig(int cost_func, double scale = 1.5, bool check = false)
        : cost_function(cost_func)
        , scaling(scale)
        , check_connectivity(check)
    {
    }
};

/**
 * @brief Calculates cost matrices for RMSD atom assignment problems
 * Claude Generated - Extracted from RMSDDriver::MakeCostMatrix methods
 */
class CostMatrixCalculator {
public:
    CostMatrixCalculator() = default;
    ~CostMatrixCalculator() = default;

    /**
     * @brief Calculate cost matrix from geometries with specific atoms
     * Claude Generated - Refactored from RMSDDriver::MakeCostMatrix
     */
    static std::pair<double, Matrix> calculate(const Geometry& reference,
        const Geometry& target,
        const std::vector<int>& reference_atoms,
        const std::vector<int>& target_atoms,
        const CostMatrixConfig& config = CostMatrixConfig());

    /**
     * @brief Calculate cost matrix from full geometries
     * Claude Generated - Refactored from RMSDDriver::MakeCostMatrix
     */
    static std::pair<double, Matrix> calculate(const Geometry& reference,
        const Geometry& target,
        const CostMatrixConfig& config = CostMatrixConfig());

    /**
     * @brief Calculate cost matrix from atom index mappings
     * Claude Generated - Refactored from RMSDDriver::MakeCostMatrix
     */
    static std::pair<double, Matrix> calculate(const std::vector<int>& reference_indices,
        const std::vector<int>& target_indices,
        const Molecule& reference_mol,
        const Molecule& target_mol,
        const CostMatrixConfig& config = CostMatrixConfig());

    /**
     * @brief Calculate cost matrix from permutation vector
     * Claude Generated - Refactored from RMSDDriver::MakeCostMatrix
     */
    static std::pair<double, Matrix> calculate(const std::vector<int>& permutation,
        const Molecule& reference_mol,
        const Molecule& target_mol,
        const CostMatrixConfig& config = CostMatrixConfig());

private:
    /**
     * @brief Cost function calculation - supports multiple cost types
     * Claude Generated - Extracted from RMSDDriver::Cost static method
     */
    static inline double calculateCost(double distance, double norm, int cost_function)
    {
        switch (cost_function) {
        case 1:
            return distance * distance;
        case 2:
            return distance;
        case 3:
            return distance + norm;
        case 4:
            return distance * distance + norm * norm;
        case 5:
            return distance * norm;
        case 6:
            return distance * distance * norm * norm;
        default:
            return distance * distance;
        }
    }
};