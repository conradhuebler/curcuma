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

#include "rmsd_assignment.h"
#include "munkres.h" // For MunkressAssign function
#include "src/core/global.h" // For logging
#include "src/global_config.h"

extern "C" {
#include "src/capabilities/c_code/interface.h" // For C assignment function
}

#include <fmt/core.h>

// Claude Generated - Refactored from RMSDDriver::SolveCostMatrix
std::vector<int> MunkresAssignmentSolver::solve(Matrix& cost_matrix, const AssignmentConfig& config)
{
    std::vector<int> assignment;
    assignment.resize(cost_matrix.rows());

    double difference = cost_matrix.cwiseAbs().sum();
    int iter = 0;

    for (iter = 0; iter < config.max_iterations && difference > config.convergence_threshold; ++iter) {
        const int dim = cost_matrix.rows();

        if (config.use_munkres) {
            // Use C++ Munkres implementation
            auto result = MunkressAssign(cost_matrix);

            // Extract assignment from result matrix
            for (int i = 0; i < result.cols(); ++i) {
                for (int j = 0; j < result.rows(); ++j) {
                    if (result(i, j) == 1) {
                        assignment[i] = j;
                        break;
                    }
                }
            }
        } else {
            // Use C implementation as fallback
            double* table = new double[dim * dim];
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    table[i * dim + j] = cost_matrix(i, j);
                }
            }

            int* order = new int[dim];
            assign(dim, table, order);

            for (int i = 0; i < dim; ++i) {
                assignment[i] = order[i];
            }

            delete[] table;
            delete[] order;
        }

        // For iterative solving, the caller should update the cost matrix
        // This simple implementation returns after one iteration
        break;
    }

    m_last_iterations = iter;

    // Debug: Munkres assignment completed

    return assignment;
}