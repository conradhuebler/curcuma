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

#include "rmsd_costmatrix.h"
#include <algorithm>
#include <cmath>

// Claude Generated - Refactored from RMSDDriver::MakeCostMatrix(Geometry, Geometry, vector, vector, int)
std::pair<double, Matrix> CostMatrixCalculator::calculate(const Geometry& reference,
    const Geometry& target,
    const std::vector<int>& reference_atoms,
    const std::vector<int>& target_atoms,
    const CostMatrixConfig& config)
{
    const double penalty = 1e23;

    if (reference.rows() != target.rows()) {
        return std::pair<double, Matrix>(penalty, Matrix());
    }

    if (reference_atoms.size() != target_atoms.size()) {
        return std::pair<double, Matrix>(penalty, Matrix());
    }

    const size_t atom_count = reference_atoms.size();
    Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(atom_count, atom_count);
    double sum = 0;

    for (size_t i = 0; i < atom_count; ++i) {
        double min = penalty;

        for (size_t j = 0; j < atom_count; ++j) {
            double d = (target.row(j) - reference.row(i)).norm();
            double norm = (target.row(j).norm() - reference.row(i).norm());

            distance(i, j) = calculateCost(d, norm, config.cost_function);

            // Add penalty for different atom types
            distance(i, j) += penalty * (reference_atoms[i] != target_atoms[j]);

            min = std::min(min, distance(i, j));
        }
        sum += min;
    }

    return std::pair<double, Matrix>(sum, distance);
}

// Claude Generated - Refactored from RMSDDriver::MakeCostMatrix(Geometry, Geometry)
std::pair<double, Matrix> CostMatrixCalculator::calculate(const Geometry& reference,
    const Geometry& target,
    const CostMatrixConfig& config)
{
    const double penalty = 1e23;

    if (reference.rows() != target.rows()) {
        return std::pair<double, Matrix>(penalty, Matrix());
    }

    const int atom_count = reference.rows();
    Eigen::MatrixXd distance = Eigen::MatrixXd::Zero(atom_count, atom_count);
    double sum = 0;

    for (int i = 0; i < atom_count; ++i) {
        double min = penalty;

        for (int j = 0; j < atom_count; ++j) {
            double d = (target.row(j) - reference.row(i)).norm();
            double norm = (target.row(j).norm() - reference.row(i).norm());

            distance(i, j) = calculateCost(d, norm, config.cost_function);

            min = std::min(min, distance(i, j));
        }
        sum += min;
    }

    return std::pair<double, Matrix>(sum, distance);
}

// Claude Generated - Refactored from RMSDDriver::MakeCostMatrix(vector, vector)
std::pair<double, Matrix> CostMatrixCalculator::calculate(const std::vector<int>& reference_indices,
    const std::vector<int>& target_indices,
    const Molecule& reference_mol,
    const Molecule& target_mol,
    const CostMatrixConfig& config)
{
    // Extract geometries for the specified atoms
    Geometry ref_geom(reference_indices.size(), 3);
    Geometry tar_geom(target_indices.size(), 3);
    std::vector<int> ref_atoms, tar_atoms;

    for (size_t i = 0; i < reference_indices.size(); ++i) {
        ref_geom.row(i) = reference_mol.Atom(reference_indices[i]).second;
        ref_atoms.push_back(reference_mol.Atom(reference_indices[i]).first);
    }

    for (size_t i = 0; i < target_indices.size(); ++i) {
        tar_geom.row(i) = target_mol.Atom(target_indices[i]).second;
        tar_atoms.push_back(target_mol.Atom(target_indices[i]).first);
    }

    return calculate(ref_geom, tar_geom, ref_atoms, tar_atoms, config);
}

// Claude Generated - Refactored from RMSDDriver::MakeCostMatrix(vector)
std::pair<double, Matrix> CostMatrixCalculator::calculate(const std::vector<int>& permutation,
    const Molecule& reference_mol,
    const Molecule& target_mol,
    const CostMatrixConfig& config)
{
    std::vector<int> reference_indices(permutation.size());
    for (size_t i = 0; i < permutation.size(); ++i) {
        reference_indices[i] = i;
    }

    return calculate(reference_indices, permutation, reference_mol, target_mol, config);
}