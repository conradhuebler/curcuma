/*
 * <Spatial Cell List for Neighbor Queries>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (April 2026): Minimal spatial cell list for O(N) neighbor queries.
 * Replaces O(N^2) double loops in non-bonded detection (HB, XB, dispersion pairs).
 * Cell size equals cutoff so only 27 neighbour cells need to be checked.
 */

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <functional>

/*!
 * @brief Minimal spatial cell-list for O(N) neighbor queries.
 *
 * Header-only, no PBC (molecular systems).  Cell size == cutoff,
 * therefore only 27 neighbor cells have to be inspected.
 * For N < ~800 the build overhead usually does not pay off;
 * callers should keep a sequential fallback path.
 */
class SpatialCellList {
public:
    SpatialCellList() = default;

    /*! @brief Build the cell list from Cartesian coordinates.
     *  @param coords  Nx3 matrix, rows = atom positions (any unit).
     *  @param cutoff  Pair cutoff (same unit as coords).  Cell size = cutoff.
     */
    void build(const Eigen::MatrixXd& coords, double cutoff)
    {
        m_coords = coords;
        m_natoms = static_cast<int>(coords.rows());
        m_cutoff = cutoff;
        m_cells.clear();

        if (m_natoms == 0 || cutoff <= 0.0)
            return;

        // Bounding box
        m_min = coords.row(0);
        m_max = coords.row(0);
        for (int i = 1; i < m_natoms; ++i) {
            for (int d = 0; d < 3; ++d) {
                m_min[d] = std::min(m_min[d], coords(i, d));
                m_max[d] = std::max(m_max[d], coords(i, d));
            }
        }

        // Number of cells in each dimension (at least 1)
        m_inv_cell_size = 1.0 / cutoff;
        m_ncells_x = std::max(1, static_cast<int>(std::ceil((m_max[0] - m_min[0]) * m_inv_cell_size)));
        m_ncells_y = std::max(1, static_cast<int>(std::ceil((m_max[1] - m_min[1]) * m_inv_cell_size)));
        m_ncells_z = std::max(1, static_cast<int>(std::ceil((m_max[2] - m_min[2]) * m_inv_cell_size)));

        m_cells.resize(static_cast<size_t>(m_ncells_x) * m_ncells_y * m_ncells_z);

        // Bin atoms
        for (int i = 0; i < m_natoms; ++i) {
            int ci = atomCellIndex(i);
            m_cells[ci].push_back(i);
        }
    }

    /*! @brief True if no atoms were inserted. */
    bool empty() const { return m_natoms == 0; }

    /*! @brief Number of atoms stored. */
    int natoms() const { return m_natoms; }

    /*! @brief Cutoff used during build. */
    double cutoff() const { return m_cutoff; }

    /*! @brief Iterate over all unique pairs (i < j) with distance <= cutoff.
     *  @param callback  Invoked as callback(i, j, r2) for every pair.
     *                   i < j is guaranteed; r2 is the squared distance.
     */
    template <typename Func>
    void forEachPair(Func&& callback) const
    {
        if (m_natoms == 0)
            return;

        const double cutoff_sq = m_cutoff * m_cutoff;

        for (int i = 0; i < m_natoms; ++i) {
            int cix, ciy, ciz;
            atomCellCoords(i, cix, ciy, ciz);

            for (int dz = -1; dz <= 1; ++dz) {
                int nz = ciz + dz;
                if (nz < 0 || nz >= m_ncells_z)
                    continue;
                for (int dy = -1; dy <= 1; ++dy) {
                    int ny = ciy + dy;
                    if (ny < 0 || ny >= m_ncells_y)
                        continue;
                    for (int dx = -1; dx <= 1; ++dx) {
                        int nx = cix + dx;
                        if (nx < 0 || nx >= m_ncells_x)
                            continue;
                        int nci = cellIndex(nx, ny, nz);

                        // Eliminate duplicate pairs: only visit cells with
                        // flat index >= current cell index.  When the cell
                        // is the same, enforce j > i.
                        if (nci < cellIndex(cix, ciy, ciz))
                            continue;

                        const auto& cell = m_cells[nci];
                        for (int j : cell) {
                            if (nci == cellIndex(cix, ciy, ciz) && j <= i)
                                continue;
                            double r2 = (m_coords.row(i) - m_coords.row(j)).squaredNorm();
                            if (r2 <= cutoff_sq) {
                                callback(i, j, r2);
                            }
                        }
                    }
                }
            }
        }
    }

    /*! @brief Iterate over neighbors of atom @p i within a (possibly smaller) cutoff.
     *  @param i          Central atom.
     *  @param cutoff_sq  Squared cutoff distance (may be <= m_cutoff^2).
     *  @param callback   Invoked as callback(j, r2) for every neighbor j != i.
     */
    template <typename Func>
    void forEachNeighbor(int i, double cutoff_sq, Func&& callback) const
    {
        if (m_natoms == 0 || i < 0 || i >= m_natoms)
            return;

        int cix, ciy, ciz;
        atomCellCoords(i, cix, ciy, ciz);

        for (int dz = -1; dz <= 1; ++dz) {
            int nz = ciz + dz;
            if (nz < 0 || nz >= m_ncells_z)
                continue;
            for (int dy = -1; dy <= 1; ++dy) {
                int ny = ciy + dy;
                if (ny < 0 || ny >= m_ncells_y)
                    continue;
                for (int dx = -1; dx <= 1; ++dx) {
                    int nx = cix + dx;
                    if (nx < 0 || nx >= m_ncells_x)
                        continue;
                    int nci = cellIndex(nx, ny, nz);
                    for (int j : m_cells[nci]) {
                        if (j == i)
                            continue;
                        double r2 = (m_coords.row(i) - m_coords.row(j)).squaredNorm();
                        if (r2 <= cutoff_sq) {
                            callback(j, r2);
                        }
                    }
                }
            }
        }
    }

private:
    int m_natoms = 0;
    double m_cutoff = 0.0;
    double m_inv_cell_size = 0.0;
    int m_ncells_x = 0, m_ncells_y = 0, m_ncells_z = 0;
    Eigen::Vector3d m_min = Eigen::Vector3d::Zero();
    Eigen::Vector3d m_max = Eigen::Vector3d::Zero();
    std::vector<std::vector<int>> m_cells;
    Eigen::MatrixXd m_coords;

    inline int cellIndex(int cx, int cy, int cz) const
    {
        return cx + m_ncells_x * (cy + m_ncells_y * cz);
    }

    inline int atomCellIndex(int atom) const
    {
        int cx = static_cast<int>((m_coords(atom, 0) - m_min[0]) * m_inv_cell_size);
        int cy = static_cast<int>((m_coords(atom, 1) - m_min[1]) * m_inv_cell_size);
        int cz = static_cast<int>((m_coords(atom, 2) - m_min[2]) * m_inv_cell_size);
        cx = std::clamp(cx, 0, m_ncells_x - 1);
        cy = std::clamp(cy, 0, m_ncells_y - 1);
        cz = std::clamp(cz, 0, m_ncells_z - 1);
        return cellIndex(cx, cy, cz);
    }

    inline void atomCellCoords(int atom, int& cx, int& cy, int& cz) const
    {
        cx = static_cast<int>((m_coords(atom, 0) - m_min[0]) * m_inv_cell_size);
        cy = static_cast<int>((m_coords(atom, 1) - m_min[1]) * m_inv_cell_size);
        cz = static_cast<int>((m_coords(atom, 2) - m_min[2]) * m_inv_cell_size);
        cx = std::clamp(cx, 0, m_ncells_x - 1);
        cy = std::clamp(cy, 0, m_ncells_y - 1);
        cz = std::clamp(cz, 0, m_ncells_z - 1);
    }
};
