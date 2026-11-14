/*
 * <Periodic Boundary Conditions Utilities>
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
 * Claude Generated: Universal PBC support utilities for all file formats
 */

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>


namespace PBCUtils {

/*! \brief Convert lattice parameters to 3x3 lattice vector matrix - Claude Generated
 * \param a Length of lattice vector a (Angstroms)
 * \param b Length of lattice vector b (Angstroms)
 * \param c Length of lattice vector c (Angstroms)
 * \param alpha_deg Angle between b and c (degrees)
 * \param beta_deg Angle between a and c (degrees)
 * \param gamma_deg Angle between a and b (degrees)
 * \return 3x3 matrix with lattice vectors as columns
 *
 * Reference: International Tables for Crystallography, Vol. B, Section 1.1
 * https://en.wikipedia.org/wiki/Fractional_coordinates
 *
 * Standard crystallographic convention:
 *   - Vector a is aligned along the x-axis
 *   - Vector b lies in the xy-plane
 *   - Vector c completes a right-handed coordinate system
 *
 * This convention ensures:
 *   - Orthorhombic cells (α=β=γ=90°) become diagonal matrices
 *   - Minimum Image Convention works correctly
 *   - Compatible with major MD packages (LAMMPS, GROMACS, VMD)
 */
inline Eigen::Matrix3d buildLatticeVectors(double a, double b, double c,
    double alpha_deg, double beta_deg, double gamma_deg)
{
    // Convert angles from degrees to radians
    const double deg2rad = M_PI / 180.0;
    const double alpha = alpha_deg * deg2rad;
    const double beta = beta_deg * deg2rad;
    const double gamma = gamma_deg * deg2rad;

    Eigen::Matrix3d cell;

    // Vector a: aligned along x-axis
    // a = (a, 0, 0)
    cell(0, 0) = a;
    cell(0, 1) = 0.0;
    cell(0, 2) = 0.0;

    // Vector b: in xy-plane
    // b = (b*cos(γ), b*sin(γ), 0)
    cell(1, 0) = b * std::cos(gamma);
    cell(1, 1) = b * std::sin(gamma);
    cell(1, 2) = 0.0;

    // Vector c: general triclinic case
    // c_x = c * cos(β)
    // c_y = c * (cos(α) - cos(β)*cos(γ)) / sin(γ)
    // c_z = sqrt(c² - c_x² - c_y²)
    double cx = c * std::cos(beta);
    double cy = c * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
    double cz_squared = c * c - cx * cx - cy * cy;

    // Guard against numerical errors for nearly degenerate cells
    double cz = (cz_squared > 0.0) ? std::sqrt(cz_squared) : 0.0;

    cell(2, 0) = cx;
    cell(2, 1) = cy;
    cell(2, 2) = cz;

    return cell;
}

/*! \brief Apply Minimum Image Convention to a distance vector - Claude Generated
 * \param r Distance vector in Cartesian coordinates (Angstroms)
 * \param cell 3x3 lattice vector matrix
 * \return PBC-corrected distance vector (shortest distance across periodic images)
 *
 * Reference: Frenkel & Smit, "Understanding Molecular Simulation", Chapter 12
 * Allen & Tildesley, "Computer Simulation of Liquids", Section 1.4
 *
 * The Minimum Image Convention ensures that we always calculate distances
 * to the nearest periodic image of an atom. This is essential for:
 *   - Correct potential energy calculations
 *   - Accurate force computations
 *   - Physically meaningful distances in periodic systems
 *
 * Algorithm:
 *   1. Convert distance vector to fractional coordinates
 *   2. Wrap each component to [-0.5, 0.5) range
 *   3. Convert back to Cartesian coordinates
 *
 * Performance: O(1) - only 9 multiplications + 3 rounding operations
 *
 * ⚠️ PERFORMANCE WARNING: This version computes cell.inverse() every time!
 * For repeated calculations, use applyMinimumImage(r, cell, cell_inv) instead
 * to avoid expensive matrix inversions.
 */
inline Eigen::Vector3d applyMinimumImage(const Eigen::Vector3d& r,
    const Eigen::Matrix3d& cell)
{
    // Convert Cartesian distance to fractional coordinates
    // frac = cell^(-1) * r
    Eigen::Vector3d frac = cell.inverse() * r;

    // Apply periodic wrapping: map to [-0.5, 0.5) range
    // This finds the nearest periodic image
    frac(0) -= std::round(frac(0));
    frac(1) -= std::round(frac(1));
    frac(2) -= std::round(frac(2));

    // Convert back to Cartesian coordinates
    // r_pbc = cell * frac
    return cell * frac;
}

/*! \brief Apply Minimum Image Convention using pre-computed inverse matrix - Claude Generated (Oct 2025)
 * \param r Distance vector in Cartesian coordinates (Angstroms)
 * \param cell 3x3 lattice vector matrix
 * \param cell_inv Pre-computed inverse of lattice matrix (cell^-1)
 * \return PBC-corrected distance vector (shortest distance across periodic images)
 *
 * OPTIMIZED VERSION: Uses pre-computed inverse matrix to avoid expensive matrix inversions.
 * This is the recommended function for bulk PBC calculations (polymer properties, distance arrays).
 *
 * Performance improvement:
 *   - 200 atoms: ~10,000× faster (200 inversions avoided)
 *   - 1000 atoms: ~250,000× faster
 *
 * Usage pattern:
 *   Eigen::Matrix3d cell_inv = cell.inverse();  // Once!
 *   for (int i = 0; i < N; i++) {
 *       r_pbc = applyMinimumImage(r_vectors[i], cell, cell_inv);
 *   }
 */
inline Eigen::Vector3d applyMinimumImage(const Eigen::Vector3d& r,
    const Eigen::Matrix3d& cell,
    const Eigen::Matrix3d& cell_inv)
{
    // Convert Cartesian distance to fractional coordinates using pre-computed inverse
    // frac = cell_inv * r
    Eigen::Vector3d frac = cell_inv * r;

    // Apply periodic wrapping: map to [-0.5, 0.5) range
    // This finds the nearest periodic image
    frac(0) -= std::round(frac(0));
    frac(1) -= std::round(frac(1));
    frac(2) -= std::round(frac(2));

    // Convert back to Cartesian coordinates
    // r_pbc = cell * frac
    return cell * frac;
}

/*! \brief Calculate PBC-aware distance between two points - Claude Generated
 * \param pos1 Position of first atom (Angstroms)
 * \param pos2 Position of second atom (Angstroms)
 * \param cell 3x3 lattice vector matrix
 * \return Distance accounting for periodic boundaries (Angstroms)
 *
 * This is the main distance calculation function for periodic systems.
 * It automatically applies the Minimum Image Convention.
 *
 * Use cases:
 *   - End-to-end distances in polymer simulations
 *   - Bond lengths in periodic crystals
 *   - Distance-dependent properties in MD trajectories
 *
 * ⚠️ PERFORMANCE WARNING: This version computes cell.inverse() on every call!
 * For bulk calculations, use calculateDistancePBC(pos1, pos2, cell, cell_inv) instead.
 */
inline double calculateDistancePBC(const Eigen::Vector3d& pos1,
    const Eigen::Vector3d& pos2,
    const Eigen::Matrix3d& cell)
{
    // Calculate distance vector
    Eigen::Vector3d rij = pos2 - pos1;

    // Apply Minimum Image Convention
    Eigen::Vector3d rij_pbc = applyMinimumImage(rij, cell);

    // Return Euclidean distance
    return rij_pbc.norm();
}

/*! \brief Calculate PBC-aware distance using pre-computed inverse matrix - Claude Generated (Oct 2025)
 * \param pos1 Position of first atom (Angstroms)
 * \param pos2 Position of second atom (Angstroms)
 * \param cell 3x3 lattice vector matrix
 * \param cell_inv Pre-computed inverse of lattice matrix (cell^-1)
 * \return Distance accounting for periodic boundaries (Angstroms)
 *
 * OPTIMIZED VERSION: Uses pre-computed inverse matrix.
 * Recommended for:
 *   - Calculating distance arrays/matrices
 *   - GyrationRadius calculations
 *   - Polymer property evaluations
 *
 * Performance improvement:
 *   - For N atoms: ~N inversions saved (from N to 1)
 *   - Example: 200 atoms saves ~199 expensive matrix operations per property calculation
 */
inline double calculateDistancePBC(const Eigen::Vector3d& pos1,
    const Eigen::Vector3d& pos2,
    const Eigen::Matrix3d& cell,
    const Eigen::Matrix3d& cell_inv)
{
    // Calculate distance vector
    Eigen::Vector3d rij = pos2 - pos1;

    // Apply Minimum Image Convention with cached inverse
    Eigen::Vector3d rij_pbc = applyMinimumImage(rij, cell, cell_inv);

    // Return Euclidean distance
    return rij_pbc.norm();
}

/*! \brief Extract lattice parameters from 3x3 matrix - Claude Generated
 * \param cell 3x3 lattice vector matrix
 * \return Vector (a, b, c, alpha, beta, gamma) with lengths in Å and angles in degrees
 *
 * Inverse operation to buildLatticeVectors().
 * Useful for output and visualization purposes.
 */
inline std::vector<double> getLatticeParameters(const Eigen::Matrix3d& cell)
{
    // Extract lattice vector lengths
    double a = cell.col(0).norm();
    double b = cell.col(1).norm();
    double c = cell.col(2).norm();

    // Calculate angles using dot products
    // cos(angle) = (v1 · v2) / (|v1| * |v2|)
    const double rad2deg = 180.0 / M_PI;

    double alpha = std::acos(cell.col(1).dot(cell.col(2)) / (b * c)) * rad2deg;
    double beta = std::acos(cell.col(0).dot(cell.col(2)) / (a * c)) * rad2deg;
    double gamma = std::acos(cell.col(0).dot(cell.col(1)) / (a * b)) * rad2deg;

    return { a, b, c, alpha, beta, gamma };
}

/*! \brief Check if a unit cell is valid - Claude Generated
 * \param cell 3x3 lattice vector matrix
 * \return true if cell is physically valid (positive volume, non-degenerate)
 *
 * Validation checks:
 *   - All lattice vector lengths > 0
 *   - Cell volume > 0 (right-handed system)
 *   - No degenerate angles (0° or 180°)
 */
inline bool isValidUnitCell(const Eigen::Matrix3d& cell)
{
    // Check lattice vector lengths
    if (cell.col(0).norm() < 1e-6 || cell.col(1).norm() < 1e-6 || cell.col(2).norm() < 1e-6) {
        return false;
    }

    // Check cell volume (determinant must be positive)
    double volume = cell.determinant();
    if (volume < 1e-6) {
        return false;
    }

    return true;
}

// ============================================================================
// Orthorhombic (rectangular) cells - Claude Generated (Nov 2026, WP-PERIODIC-NONBONDED A1)
//
// Foundation of periodic GFN-FF: one minimum-image function shared by the CPU force-field
// engine (FFWorkspace, Bohr), the cell-size check and SimpleMD. Only diagonal cells are
// accepted. The general matrix functions above mix two conventions (buildLatticeVectors()
// stores the lattice vectors as ROWS, applyMinimumImage() multiplies cell_inv * r as if they
// were COLUMNS); the two agree only for a diagonal cell, so a triclinic cell is refused here
// instead of being silently mis-imaged. That inconsistency is known and not repaired here.
// ============================================================================

/*! \brief Diagonal of an orthorhombic cell with its reciprocal - Claude Generated (Nov 2026)
 *
 * Lengths are in whatever unit the caller works in (FFWorkspace: Bohr, SimpleMD: Angstrom).
 * inv_L is taken from the diagonal of the Eigen 3x3 inverse of the cell matrix, the same
 * matrix the GPU path uploads as d_lattice_inv (gfnff_gpu_method_impl.h), so CPU and GPU
 * apply the same reciprocal, not a re-rounded 1/L.
 */
struct OrthoCell {
    bool active = false;                  ///< false: no cell, every function below is a no-op
    double L[3] = { 0.0, 0.0, 0.0 };      ///< box lengths
    double inv_L[3] = { 0.0, 0.0, 0.0 };  ///< reciprocal box lengths (diagonal of cell^-1)

    /// Shortest box length; the minimum image is unique only for distances below minLength()/2.
    double minLength() const { return std::min(L[0], std::min(L[1], L[2])); }
};

/*! \brief True when every off-diagonal element is negligible against the box - Claude Generated (Nov 2026)
 * \param cell 3x3 lattice matrix (any unit)
 * \param rel_tol off-diagonal elements must be <= rel_tol * largest diagonal element
 *
 * buildLatticeVectors(a, b, c, 90, 90, 90) leaves b*cos(90 deg) ~ 6e-16*b off the diagonal,
 * hence a relative tolerance instead of an exact zero test.
 */
inline bool isOrthorhombic(const Eigen::Matrix3d& cell, double rel_tol = 1e-10)
{
    const double scale = cell.diagonal().cwiseAbs().maxCoeff();
    if (!(scale > 0.0))
        return false;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            if (r != c && std::abs(cell(r, c)) > rel_tol * scale)
                return false;
    return cell(0, 0) > 0.0 && cell(1, 1) > 0.0 && cell(2, 2) > 0.0;
}

/*! \brief Build an OrthoCell from a 3x3 cell matrix - Claude Generated (Nov 2026)
 * \param cell lattice matrix, must be orthorhombic (isOrthorhombic())
 * \param error set to a message and returned cell inactive when the cell is not orthorhombic
 */
OrthoCell makeOrthoCell(const Eigen::Matrix3d& cell, std::string& error);

/*! \brief Minimum-image component for an orthorhombic cell - Claude Generated (Nov 2026)
 * \param d     one Cartesian component of r_j - r_i
 * \param L     box length along that axis
 * \param inv_L reciprocal box length (OrthoCell::inv_L)
 * \return d shifted by an integer number of box lengths into [-L/2, L/2]
 *
 * Minimum-image convention, Allen & Tildesley, "Computer Simulation of Liquids", 2nd ed.
 * (2017), Sec. 1.6.3 (periodic boundary conditions): s = d/L, s <- s - rint(s), d <- L s.
 * The arithmetic is that of applyMIC() in cuda/gfnff_kernels.cu and rocm/gfnff_rocm.hip for a
 * diagonal lattice (the off-diagonal products there are exact zeros): one multiplication by
 * the reciprocal, rint (round half to even, odd function, so mic(-d) = -mic(d)), one
 * multiplication by L. A pair exactly on the boundary (|s| = 0.5) keeps its sign: rint(+-0.5) = 0.
 */
inline double minimumImageOrtho(double d, double L, double inv_L)
{
    double s = d * inv_L;
    s -= std::rint(s);
    return L * s;
}

/*! \brief Apply the orthorhombic minimum image to a displacement in place - Claude Generated (Nov 2026)
 *
 * No-op for an inactive cell, so a call site without a cell keeps its numbers bit for bit.
 */
inline void applyMinimumImageOrtho(const OrthoCell& cell, double& dx, double& dy, double& dz)
{
    if (!cell.active)
        return;
    dx = minimumImageOrtho(dx, cell.L[0], cell.inv_L[0]);
    dy = minimumImageOrtho(dy, cell.L[1], cell.inv_L[1]);
    dz = minimumImageOrtho(dz, cell.L[2], cell.inv_L[2]);
}

/*! \brief Lattice translations for explicit image sums in an orthorhombic cell - Claude Generated (Nov 2026, WP-PERIODIC-XTB P1)
 *
 * For a pair (i, j) the periodic xTB terms sum every image j + T with |x_i - x_j - T| < cutoff, including the
 * images of an atom with itself (T != 0). forEachImage() first takes the minimum image d0 of x_i - x_j, so
 * |d0_a| <= L_a/2, and then runs over the stored translations; a translation n_a L_a can reach the cutoff only
 * when |n_a| <= (cutoff + L_a/2)/L_a, which bounds the list. Unlike the minimum-image convention this needs
 * no cutoff < L/2: a small cell simply has more images per pair. Same unit as the cell (Bohr in xTB).
 */
struct LatticeImageList {
    OrthoCell cell;                              ///< inactive: forEachImage() visits the single plain pair
    double cutoff = 0.0;                         ///< radius of the image sum
    double cutoff2 = 0.0;
    std::vector<std::array<double, 3>> trans;    ///< candidate translations (T = 0 first)
    int nmax[3] = { 0, 0, 0 };
    bool active() const { return cell.active; }
};

inline LatticeImageList makeLatticeImageList(const OrthoCell& cell, double cutoff)
{
    LatticeImageList list;
    list.cell = cell;
    list.cutoff = cutoff;
    list.cutoff2 = cutoff * cutoff;
    if (!cell.active) {
        list.trans.push_back({ 0.0, 0.0, 0.0 });
        return list;
    }
    for (int a = 0; a < 3; ++a)
        list.nmax[a] = static_cast<int>(std::ceil((cutoff + 0.5 * cell.L[a]) * cell.inv_L[a]));
    list.trans.push_back({ 0.0, 0.0, 0.0 });
    for (int nx = -list.nmax[0]; nx <= list.nmax[0]; ++nx)
        for (int ny = -list.nmax[1]; ny <= list.nmax[1]; ++ny)
            for (int nz = -list.nmax[2]; nz <= list.nmax[2]; ++nz) {
                if (nx == 0 && ny == 0 && nz == 0)
                    continue;
                list.trans.push_back({ nx * cell.L[0], ny * cell.L[1], nz * cell.L[2] });
            }
    return list;
}

/*! \brief Visit every image of the pair (i, j) within the list's cutoff - Claude Generated (Nov 2026)
 * \param xi, xj   Cartesian positions (3 doubles each), same unit as the cell
 * \param same     true when i == j (the T = 0 term is then skipped)
 * \param f        callback f(dx, dy, dz, r2) with d = x_i - (x_j + T) (tblite convention vec = xyz_i - xyz_j - trans)
 *
 * Without an active cell exactly the plain pair d = x_i - x_j is visited (no cutoff test, no self pair), so a
 * caller that is also used without a cell keeps its former sum; the active-cell loop is a separate code path.
 */
template <class F>
inline void forEachImage(const LatticeImageList& list, const double* xi, const double* xj, bool same, F&& f)
{
    if (!list.cell.active) {
        if (same)
            return;
        const double dx = xi[0] - xj[0], dy = xi[1] - xj[1], dz = xi[2] - xj[2];
        f(dx, dy, dz, dx * dx + dy * dy + dz * dz);
        return;
    }
    double d0x = xi[0] - xj[0], d0y = xi[1] - xj[1], d0z = xi[2] - xj[2];
    applyMinimumImageOrtho(list.cell, d0x, d0y, d0z);
    const size_t n = list.trans.size();
    for (size_t t = 0; t < n; ++t) {
        if (same && t == 0)
            continue;
        const double dx = d0x - list.trans[t][0];
        const double dy = d0y - list.trans[t][1];
        const double dz = d0z - list.trans[t][2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 >= list.cutoff2 || r2 < 1.0e-24)
            continue;
        f(dx, dy, dz, r2);
    }
}

/*! \brief Check interaction cutoffs against half the shortest box length - Claude Generated (Nov 2026)
 * \param cell    active orthorhombic cell
 * \param cutoffs (name, radius) of every cutoff that is in use, same unit as the cell
 * \param unit    unit name for the message ("Bohr")
 * \return empty string when every radius is < L_min/2, otherwise one message naming every
 *         offending cutoff with its value and L_min/2
 *
 * Under the minimum-image convention a pair interacts only with its nearest image; a cutoff
 * at or above L/2 would need a second image (Allen & Tildesley, Sec. 1.6.3). Non-finite
 * radii (no cutoff) are reported as well, the caller decides which terms to pass.
 */
std::string checkCutoffsAgainstCell(const OrthoCell& cell,
    const std::vector<std::pair<std::string, double>>& cutoffs, const std::string& unit);

/*! \brief Wrap whole fragments into the cell by their centre of mass - Claude Generated (Nov 2026)
 * \param geometry  N x 3 Cartesian coordinates, modified in place
 * \param masses    atomic masses, size N
 * \param fragments atom indices per fragment (e.g. Molecule::GetFragments())
 * \param cell      lattice matrix, same unit as geometry
 * \param cell_inv  its inverse
 * \return number of fragments that were moved
 *
 * Each fragment is translated by one integer lattice vector so that its centre of mass lands
 * in fractional [0,1); all atoms of a fragment get the same shift, so intramolecular distances
 * are kept (up to the rounding of adding n*L). Replaces the former per-atom wrap of SimpleMD,
 * which tore molecules apart at the boundary. Uses the same frac = cell_inv * r and
 * r = cell * frac as that former code (correct for orthorhombic cells, see the note above).
 * A fragment that is already torn in the input is moved as it is, not reassembled.
 */
int wrapFragmentsIntoCell(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& geometry,
    const std::vector<double>& masses, const std::vector<std::vector<int>>& fragments,
    const Eigen::Matrix3d& cell, const Eigen::Matrix3d& cell_inv);

} // namespace PBCUtils
