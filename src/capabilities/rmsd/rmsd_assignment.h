/*
 * <Assignment Problem Solver for RMSD calculations>
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
#include <vector>

/**
 * @brief Configuration for assignment solver
 * Claude Generated - Extracted from RMSDDriver parameters
 */
struct AssignmentConfig {
    int max_iterations = 100; // Maximum iterations for convergence
    double convergence_threshold = 1e-3; // Convergence threshold
    bool use_munkres = true; // Use Munkres algorithm (vs other methods)

    AssignmentConfig() = default;
    AssignmentConfig(int max_iter, double threshold = 1e-3, bool munkres = true)
        : max_iterations(max_iter)
        , convergence_threshold(threshold)
        , use_munkres(munkres)
    {
    }
};

/**
 * @brief Abstract base class for assignment problem solvers
 * Claude Generated - Strategy pattern for different assignment algorithms
 */
class AssignmentSolver {
public:
    virtual ~AssignmentSolver() = default;

    /**
     * @brief Solve assignment problem given cost matrix
     * @param cost_matrix The cost matrix to solve
     * @param config Configuration for the solver
     * @return Vector of assignments (target indices for each reference atom)
     */
    virtual std::vector<int> solve(Matrix& cost_matrix, const AssignmentConfig& config = AssignmentConfig()) = 0;

    /**
     * @brief Get number of iterations used in last solve
     */
    virtual int getIterations() const = 0;
};

/**
 * @brief Munkres (Hungarian) algorithm implementation for assignment problems
 * Claude Generated - Refactored from RMSDDriver::SolveCostMatrix
 */
class MunkresAssignmentSolver : public AssignmentSolver {
public:
    MunkresAssignmentSolver() = default;
    ~MunkresAssignmentSolver() = default;

    /**
     * @brief Solve assignment using Munkres algorithm
     * Claude Generated - Refactored from RMSDDriver::SolveCostMatrix
     */
    std::vector<int> solve(Matrix& cost_matrix, const AssignmentConfig& config = AssignmentConfig()) override;

    /**
     * @brief Get number of iterations used in last solve
     */
    int getIterations() const override { return m_last_iterations; }

private:
    int m_last_iterations = 0;
};