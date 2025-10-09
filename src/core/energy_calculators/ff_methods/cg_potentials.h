/*
 * < CG Potentials for Curcuma Coarse-Graining >
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
 */

#pragma once

#include "forcefieldthread.h"
#include <Eigen/Dense>

using Vector3d = Eigen::Vector3d;
using Matrix3d = Eigen::Matrix3d;

namespace CGPotentials {

/*! \brief Check if shape represents a sphere
 * \param shape Vector3d with x,y,z radii
 * \param tolerance Numerical tolerance for equality check
 * \return true if all radii are equal (sphere)
 */
bool isSpherical(const Vector3d& shape, double tolerance = 1e-6);

/*! \brief Check if shape represents an ellipsoid (non-spherical)
 * \param shape Vector3d with x,y,z radii
 * \param tolerance Numerical tolerance for equality check
 * \return true if radii are different (ellipsoid)
 */
bool isEllipsoidal(const Vector3d& shape, double tolerance = 1e-6);

/*! \brief Lennard-Jones (6,12) potential for standard applications
 * \param r Distance between particles
 * \param sigma LJ size parameter
 * \param epsilon LJ energy parameter
 * \return LJ energy
 */
double calculateLJ_612(double r, double sigma, double epsilon);

/*! \brief Lennard-Jones (1,6,12) potential for SCNP simulations
 * \param r Distance between particles
 * \param sigma LJ size parameter
 * \param epsilon LJ energy parameter
 * \return Modified LJ energy with additional repulsion
 */
double calculateLJ_1612(double r, double sigma, double epsilon);

/*! \brief Convert Euler angles to rotation matrix (Z-Y-X convention)
 * \param euler_angles Vector3d with (alpha, beta, gamma) Euler angles
 * \return 3x3 rotation matrix
 * \note Prepared for ellipsoid rotations, unused for spheres
 */
Matrix3d eulerToRotationMatrix(const Vector3d& euler_angles);

/*! \brief Calculate effective radius of ellipsoid in given direction
 * \param axes Vector3d with (a,b,c) semi-axes of ellipsoid
 * \param direction Unit vector for direction of interest
 * \return Effective radius in given direction
 * \note For ellipsoid contact calculations
 */
double calculateEllipsoidRadius(const Vector3d& axes, const Vector3d& direction);

/*! \brief Calculate effective contact distance between two ellipsoids
 * \param shape_i Semi-axes of particle i
 * \param shape_j Semi-axes of particle j
 * \param orient_i Euler angles of particle i
 * \param orient_j Euler angles of particle j
 * \param contact_vector Vector from i to j
 * \return Effective contact distance for angle-dependent interactions
 * \note Prepared for ellipsoids, returns sphere equivalent for spheres
 */
double calculateEffectiveDistance(const Vector3d& shape_i, const Vector3d& shape_j,
    const Vector3d& orient_i, const Vector3d& orient_j,
    const Vector3d& contact_vector);

/*! \brief Main CG pair energy calculation (pragmatic sphere implementation)
 * \param pair vdW struct containing all CG parameters
 * \param pos_i Position of particle i
 * \param pos_j Position of particle j
 * \return Interaction energy between CG particles
 *
 * \note Current implementation:
 * - Spheres: Direct LJ calculation using pair.sigma
 * - Ellipsoids: Prepared structure with placeholder for future implementation
 */
double calculateCGPairEnergy(const vdW& pair, const Vector3d& pos_i, const Vector3d& pos_j);

} // namespace CGPotentials