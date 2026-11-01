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

#include "src/tools/pbc_utils.h"

/*!
 * @brief Minimal spatial cell-list for O(N) neighbor queries.
 *
 * Header-only.  Cell size == cutoff, therefore only 27 neighbor cells have to be inspected.
 * For N < ~800 the build overhead usually does not pay off;
 * callers should keep a sequential fallback path.
 *
 * Periodic mode (Claude Generated, Nov 2026, WP-PERIODIC-NONBONDED A2): build() with an active
 * PBCUtils::OrthoCell bins the atoms by their wrapped fractional coordinate into a grid that
 * spans the cell (cell size >= cutoff), the neighbour cells wrap around, and every squared
 * distance handed to a callback is that of the minimum-image vector
 * (PBCUtils::applyMinimumImageOrtho of x_i - x_j). Requires cutoff < L_min/2, so at most one image
 * of an atom is within the cutoff; the caller checks that (PBCUtils::checkCutoffsAgainstCell).
 * The periodic queries visit pairs in a different order than the non-periodic ones; without a
 * cell (or with an inactive one) the class behaves exactly as before.
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
        m_periodic = false;  // Nov 2026: a rebuild without a cell is non-periodic again
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

        // Number of cells in each dimension (at least 1).
        //
        // Claude Generated (Sep 2026): the extent is data, and on a diverging trajectory it
        // grows without bound. The cell count is its cube, so the resize below asks for an
        // impossible allocation and the process dies - observed as a crash inside
        // SpatialCellList::build(), reached from GFNFF::detectHydrogenBondsNative(), when a
        // 7320-atom GFN-FF MD ran at too large a step. Reproduced standalone: with the
        // previous code an extent of 2e4 A in all three directions already aborts (6.4e10
        // cells requested), 2e3 A still works. A crash is the worst outcome, because the run
        // is lost together with its last snapshot, where the ordinary instability handling
        // would have reported and stopped cleanly.
        //
        // (A non-finite extent was the first suspicion - static_cast<int> of NaN is undefined
        // - but it turned out to be harmless by accident: the cast yields a negative value and
        // std::max(1, ...) collapses the grid to a single cell. It is guarded anyway, since
        // relying on that is relying on undefined behaviour.)
        //
        // Enlarging the cells is always SAFE for correctness - forEachNeighbor scans the
        // 3x3x3 block and filters on the true squared distance, so bigger cells only cost
        // time, while smaller ones would miss pairs. So when the requested grid does not fit
        // the budget, coarsen until it does; in the worst case that is one cell, i.e. the
        // O(N^2) fallback.
        const double extent_x = m_max[0] - m_min[0];
        const double extent_y = m_max[1] - m_min[1];
        const double extent_z = m_max[2] - m_min[2];
        constexpr double kMaxCells = 2.0e6;

        if (!std::isfinite(extent_x) || !std::isfinite(extent_y) || !std::isfinite(extent_z)) {
            m_inv_cell_size = 0.0;
            m_ncells_x = m_ncells_y = m_ncells_z = 1;
            m_cells.assign(1, {});
            m_cells[0].reserve(static_cast<size_t>(m_natoms));
            for (int i = 0; i < m_natoms; ++i)
                m_cells[0].push_back(i);
            return;
        }

        double cell_size = cutoff;
        double nx = 1.0, ny = 1.0, nz = 1.0;
        for (;;) {
            nx = std::max(1.0, std::ceil(extent_x / cell_size));
            ny = std::max(1.0, std::ceil(extent_y / cell_size));
            nz = std::max(1.0, std::ceil(extent_z / cell_size));
            if (nx * ny * nz <= kMaxCells)
                break;
            cell_size *= 2.0;
        }
        m_inv_cell_size = 1.0 / cell_size;
        m_ncells_x = static_cast<int>(nx);
        m_ncells_y = static_cast<int>(ny);
        m_ncells_z = static_cast<int>(nz);

        m_cells.resize(static_cast<size_t>(m_ncells_x) * m_ncells_y * m_ncells_z);

        // Bin atoms
        for (int i = 0; i < m_natoms; ++i) {
            int ci = atomCellIndex(i);
            m_cells[ci].push_back(i);
        }
    }

    /*! @brief Build a periodic cell list - Claude Generated (Nov 2026, A2).
     *  @param coords  Nx3 matrix (same unit as the cell), atoms may lie outside the cell.
     *  @param cutoff  Pair cutoff, must be < cell.minLength()/2 (not checked here).
     *  @param cell    Orthorhombic cell; inactive -> identical to build(coords, cutoff).
     */
    void build(const Eigen::MatrixXd& coords, double cutoff, const PBCUtils::OrthoCell& cell)
    {
        if (!cell.active) {
            m_periodic = false;
            build(coords, cutoff);
            return;
        }
        m_periodic = true;
        m_pbc = cell;
        m_coords = coords;
        m_natoms = static_cast<int>(coords.rows());
        m_cutoff = cutoff;
        m_cells.clear();
        if (m_natoms == 0 || cutoff <= 0.0)
            return;

        // n cells per axis with L/n >= cutoff, so a pair within the cutoff is in the same or an
        // adjacent (wrapped) cell. Same allocation budget as the non-periodic build: halving n
        // only enlarges the cells, which is always correct.
        constexpr double kMaxCells = 2.0e6;
        double n[3];
        for (int k = 0; k < 3; ++k)
            n[k] = std::max(1.0, std::floor(cell.L[k] / cutoff));
        while (n[0] * n[1] * n[2] > kMaxCells)
            for (int k = 0; k < 3; ++k)
                n[k] = std::max(1.0, std::floor(n[k] / 2.0));
        m_ncells_x = static_cast<int>(n[0]);
        m_ncells_y = static_cast<int>(n[1]);
        m_ncells_z = static_cast<int>(n[2]);
        m_cells.resize(static_cast<size_t>(m_ncells_x) * m_ncells_y * m_ncells_z);
        for (int i = 0; i < m_natoms; ++i) {
            int cx, cy, cz;
            atomCellCoords(i, cx, cy, cz);
            m_cells[cellIndex(cx, cy, cz)].push_back(i);
        }
    }

    /*! @brief True when the list was built with an active periodic cell. */
    bool periodic() const { return m_periodic; }

    /*! @brief Pair displacement x_i - x_j, minimum image in periodic mode - Claude Generated (Nov 2026).
     *  The same vector whose squared norm the queries hand to their callbacks. */
    void displacement(int i, int j, double& dx, double& dy, double& dz) const
    {
        dx = m_coords(i, 0) - m_coords(j, 0);
        dy = m_coords(i, 1) - m_coords(j, 1);
        dz = m_coords(i, 2) - m_coords(j, 2);
        if (m_periodic)
            PBCUtils::applyMinimumImageOrtho(m_pbc, dx, dy, dz);
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
        for (int i = 0; i < m_natoms; ++i)
            forEachPairOfAtom(i, callback);
    }

    /*! @brief The pairs forEachPair() visits for one first atom @p i, in the same order.
     *
     *  Claude Generated (Oct 2026): forEachPair() is the loop over i of this function, so a
     *  caller can split the atoms into ranges, work on the ranges in parallel and concatenate
     *  their results in range order to get exactly the serial sequence of pairs.
     *  @param callback  Invoked as callback(i, j, r2); every pair has i < j.
     */
    template <typename Func>
    void forEachPairOfAtom(int i, Func&& callback) const
    {
        if (m_natoms == 0 || i < 0 || i >= m_natoms)
            return;

        const double cutoff_sq = m_cutoff * m_cutoff;
        int cix, ciy, ciz;
        atomCellCoords(i, cix, ciy, ciz);

        if (m_periodic) {
            // Every distinct wrapped neighbour cell once; j > i gives each pair once.
            forEachPeriodicNeighbourCell(cix, ciy, ciz, [&](int nci) {
                for (int j : m_cells[nci]) {
                    if (j <= i)
                        continue;
                    double ddx, ddy, ddz;
                    displacement(i, j, ddx, ddy, ddz);
                    const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
                    if (r2 <= cutoff_sq)
                        callback(i, j, r2);
                }
            });
            return;
        }

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

    /*! @brief Copy of this list that keeps only the atoms with keep[atom] != 0.
     *
     *  Claude Generated (Oct 2026): same grid (same bounding box, cell size and cell
     *  numbering), every cell holds the kept atoms in their original order. Every query on the
     *  copy therefore visits the kept atoms in the order the full list visits them and skips
     *  the others, without computing their distances.
     */
    SpatialCellList subset(const std::vector<char>& keep) const
    {
        SpatialCellList out = *this;
        for (auto& cell : out.m_cells) {
            size_t w = 0;
            for (int a : cell)
                if (a >= 0 && a < static_cast<int>(keep.size()) && keep[a]) cell[w++] = a;
            cell.resize(w);
        }
        return out;
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

        if (m_periodic) {
            forEachPeriodicNeighbourCell(cix, ciy, ciz, [&](int nci) {
                for (int j : m_cells[nci]) {
                    if (j == i)
                        continue;
                    double ddx, ddy, ddz;
                    displacement(i, j, ddx, ddy, ddz);
                    const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
                    if (r2 <= cutoff_sq)
                        callback(j, r2);
                }
            });
            return;
        }

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
    bool m_periodic = false;      ///< built with an active cell (Nov 2026, A2)
    PBCUtils::OrthoCell m_pbc;    ///< that cell
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

    /*! @brief Claude Generated (Sep 2026): coordinate -> cell index along one axis, without
     *  undefined behaviour. The cast is only reached for a value that is already known to be
     *  inside [0, n-1]; the `!(t >= 0.0)` form catches NaN as well as negatives. */
    inline int axisCell(double x, double lo, int n) const
    {
        const double t = (x - lo) * m_inv_cell_size;
        if (!(t >= 0.0))
            return 0;
        const double f = std::floor(t);
        if (f >= static_cast<double>(n - 1))
            return n - 1;
        return static_cast<int>(f);
    }

    inline int atomCellIndex(int atom) const
    {
        return cellIndex(axisCell(m_coords(atom, 0), m_min[0], m_ncells_x),
                         axisCell(m_coords(atom, 1), m_min[1], m_ncells_y),
                         axisCell(m_coords(atom, 2), m_min[2], m_ncells_z));
    }

    /*! @brief Claude Generated (Nov 2026, A2): cell along one axis of the periodic grid from the
     *  wrapped fractional coordinate s - floor(s) in [0,1); NaN goes to cell 0. */
    static inline int axisCellPeriodic(double x, double inv_L, int n)
    {
        double s = x * inv_L;
        s -= std::floor(s);
        if (s >= 1.0)  // s = -1e-17 rounds to 1.0 above; that atom belongs to cell 0
            s = 0.0;
        const double t = s * static_cast<double>(n);
        if (!(t >= 0.0))
            return 0;
        const double f = std::floor(t);
        if (f >= static_cast<double>(n - 1))
            return n - 1;
        return static_cast<int>(f);
    }

    /*! @brief Claude Generated (Nov 2026, A2): the distinct cells c-1, c, c+1 (mod n) of one
     *  axis; fewer than three when n < 3, so no cell is visited twice. Returns their count. */
    static inline int periodicAxisNeighbours(int c, int n, int out[3])
    {
        if (n == 1) {
            out[0] = 0;
            return 1;
        }
        if (n == 2) {
            out[0] = 0;
            out[1] = 1;
            return 2;
        }
        out[0] = (c + n - 1) % n;
        out[1] = c;
        out[2] = (c + 1) % n;
        return 3;
    }

    template <typename Func>
    void forEachPeriodicNeighbourCell(int cix, int ciy, int ciz, Func&& f) const
    {
        int xs[3], ys[3], zs[3];
        const int nx = periodicAxisNeighbours(cix, m_ncells_x, xs);
        const int ny = periodicAxisNeighbours(ciy, m_ncells_y, ys);
        const int nz = periodicAxisNeighbours(ciz, m_ncells_z, zs);
        for (int a = 0; a < nz; ++a)
            for (int b = 0; b < ny; ++b)
                for (int c = 0; c < nx; ++c)
                    f(cellIndex(xs[c], ys[b], zs[a]));
    }

    inline void atomCellCoords(int atom, int& cx, int& cy, int& cz) const
    {
        if (m_periodic) {
            cx = axisCellPeriodic(m_coords(atom, 0), m_pbc.inv_L[0], m_ncells_x);
            cy = axisCellPeriodic(m_coords(atom, 1), m_pbc.inv_L[1], m_ncells_y);
            cz = axisCellPeriodic(m_coords(atom, 2), m_pbc.inv_L[2], m_ncells_z);
            return;
        }
        cx = axisCell(m_coords(atom, 0), m_min[0], m_ncells_x);
        cy = axisCell(m_coords(atom, 1), m_min[1], m_ncells_y);
        cz = axisCell(m_coords(atom, 2), m_min[2], m_ncells_z);
    }
};
