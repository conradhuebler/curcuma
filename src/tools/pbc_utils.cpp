/*
 * <Periodic Boundary Conditions Utilities - orthorhombic cell helpers>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED A1): the out-of-line part of the orthorhombic
 * helpers declared and documented in pbc_utils.h (too large to inline; -Winline).
 */

#include "src/tools/pbc_utils.h"

#include <fmt/format.h>

namespace PBCUtils {

OrthoCell makeOrthoCell(const Eigen::Matrix3d& cell, std::string& error)
{
    OrthoCell oc;
    error.clear();
    if (!isOrthorhombic(cell)) {
        error = fmt::format(
            "periodic cell is not orthorhombic (cell matrix rows: [{:.6f} {:.6f} {:.6f}] [{:.6f} {:.6f} {:.6f}] "
            "[{:.6f} {:.6f} {:.6f}]); only rectangular cells are supported. Triclinic cells also hit a known "
            "convention mismatch in tools/pbc_utils.h (buildLatticeVectors() stores lattice vectors as rows, "
            "applyMinimumImage() uses cell_inv * r as if they were columns).",
            cell(0, 0), cell(0, 1), cell(0, 2), cell(1, 0), cell(1, 1), cell(1, 2), cell(2, 0), cell(2, 1), cell(2, 2));
        return oc;
    }
    const Eigen::Matrix3d inv = cell.inverse();  // same inverse as the GPU upload
    for (int k = 0; k < 3; ++k) {
        oc.L[k] = cell(k, k);
        oc.inv_L[k] = inv(k, k);
    }
    oc.active = true;
    return oc;
}

std::string checkCutoffsAgainstCell(const OrthoCell& cell,
    const std::vector<std::pair<std::string, double>>& cutoffs, const std::string& unit)
{
    if (!cell.active)
        return std::string();
    const double half = 0.5 * cell.minLength();
    std::string msg;
    for (const auto& c : cutoffs) {
        if (!(c.second < half)) {
            msg += fmt::format("\n  {} = {:.4f} {} is not below L_min/2 = {:.4f} {}",
                c.first, c.second, unit, half, unit);
        }
    }
    if (msg.empty())
        return msg;
    return fmt::format("periodic cell too small for the cutoffs in use (L = {:.4f} x {:.4f} x {:.4f} {}):{}",
        cell.L[0], cell.L[1], cell.L[2], unit, msg);
}

int wrapFragmentsIntoCell(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& geometry,
    const std::vector<double>& masses, const std::vector<std::vector<int>>& fragments,
    const Eigen::Matrix3d& cell, const Eigen::Matrix3d& cell_inv)
{
    int moved = 0;
    for (const std::vector<int>& fragment : fragments) {
        double mass_sum = 0.0;
        Eigen::Vector3d com = Eigen::Vector3d::Zero();
        for (int a : fragment) {
            mass_sum += masses[a];
            com += masses[a] * Eigen::Vector3d(geometry(a, 0), geometry(a, 1), geometry(a, 2));
        }
        if (!(mass_sum > 0.0))
            continue;
        com /= mass_sum;
        const Eigen::Vector3d frac = cell_inv * com;
        const Eigen::Vector3d n(-std::floor(frac(0)), -std::floor(frac(1)), -std::floor(frac(2)));
        if (n(0) == 0.0 && n(1) == 0.0 && n(2) == 0.0)
            continue;
        const Eigen::Vector3d shift = cell * n;
        for (int a : fragment) {
            geometry(a, 0) += shift(0);
            geometry(a, 1) += shift(1);
            geometry(a, 2) += shift(2);
        }
        ++moved;
    }
    return moved;
}

} // namespace PBCUtils
