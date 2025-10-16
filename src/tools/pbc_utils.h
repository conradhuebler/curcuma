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
#include <cmath>

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
 * Performance note: For repeated distance calculations, consider caching
 * the cell.inverse() matrix to avoid redundant matrix inversions.
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

} // namespace PBCUtils
