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

#include "cg_potentials.h"
#include <cmath>

namespace CGPotentials {

bool isSpherical(const Vector3d& shape, double tolerance)
{
    return std::abs(shape[0] - shape[1]) < tolerance && std::abs(shape[1] - shape[2]) < tolerance;
}

bool isEllipsoidal(const Vector3d& shape, double tolerance)
{
    return !isSpherical(shape, tolerance);
}

double calculateLJ_612(double r, double sigma, double epsilon)
{
    // Standard Lennard-Jones (6,12) potential
    double sigma_r = sigma / r;
    double r6 = std::pow(sigma_r, 6);
    double r12 = r6 * r6;
    return 4.0 * epsilon * (r12 - r6);
}

double calculateLJ_1612(double r, double sigma, double epsilon)
{
    // Modified LJ (1,6,12) potential for SCNP simulations
    // Additional linear repulsion term for cross-linking prevention
    double sigma_r = sigma / r;
    double r6 = std::pow(sigma_r, 6);
    double r12 = r6 * r6;
    return epsilon * (r12 - r6 + 1.0); // +1 linear term
}

Matrix3d eulerToRotationMatrix(const Vector3d& euler)
{
    // Convert Euler angles (Z-Y-X convention) to rotation matrix
    // Prepared for ellipsoid rotations, unused for current sphere implementation
    double alpha = euler[0], beta = euler[1], gamma = euler[2];

    Matrix3d Rz, Ry, Rx;

    // Z rotation (yaw)
    Rz << std::cos(alpha), -std::sin(alpha), 0,
        std::sin(alpha), std::cos(alpha), 0,
        0, 0, 1;

    // Y rotation (pitch)
    Ry << std::cos(beta), 0, std::sin(beta),
        0, 1, 0,
        -std::sin(beta), 0, std::cos(beta);

    // X rotation (roll)
    Rx << 1, 0, 0,
        0, std::cos(gamma), -std::sin(gamma),
        0, std::sin(gamma), std::cos(gamma);

    return Rz * Ry * Rx;
}

double calculateEllipsoidRadius(const Vector3d& axes, const Vector3d& direction)
{
    // Calculate radius of ellipsoid in given direction
    // Formula: r(θ,φ) = abc / sqrt((bc·x)² + (ac·y)² + (ab·z)²)
    double a = axes[0], b = axes[1], c = axes[2];
    double x = direction[0], y = direction[1], z = direction[2];

    double denominator = std::pow(b * c * x, 2) + std::pow(a * c * y, 2) + std::pow(a * b * z, 2);

    // Avoid division by zero
    if (denominator < 1e-12) {
        return std::max({ a, b, c }); // Return largest axis
    }

    return (a * b * c) / std::sqrt(denominator);
}

double calculateEffectiveDistance(const Vector3d& shape_i, const Vector3d& shape_j,
    const Vector3d& orient_i, const Vector3d& orient_j,
    const Vector3d& contact_vector)
{
    // Calculate effective contact distance between two particles
    // Handles both spheres and ellipsoids

    Vector3d normalized_contact = contact_vector.normalized();

    // For spherical particles, use average radius
    if (isSpherical(shape_i) && isSpherical(shape_j)) {
        double radius_i = shape_i[0]; // All components equal for spheres
        double radius_j = shape_j[0];
        return radius_i + radius_j;
    }

    // For ellipsoids: calculate orientation-dependent contact distance
    // Rotation matrices for each particle
    Matrix3d R_i = eulerToRotationMatrix(orient_i);
    Matrix3d R_j = eulerToRotationMatrix(orient_j);

    // Effective radii in contact direction (accounting for orientation)
    double r_eff_i = calculateEllipsoidRadius(shape_i, R_i * normalized_contact);
    double r_eff_j = calculateEllipsoidRadius(shape_j, R_j * (-normalized_contact));

    return r_eff_i + r_eff_j;
}

double calculateCGPairEnergy(const vdW& pair, const Vector3d& pos_i, const Vector3d& pos_j)
{
    // Only process CG interactions
    if (pair.type != 3)
        return 0.0;

    Vector3d contact_vector = pos_j - pos_i;
    double distance = contact_vector.norm();

    // Avoid singularities at very short distances
    if (distance < 0.1) {
        return 1e6; // Large repulsive energy
    }

    double sigma_effective;

    // === PRAGMATIC IMPLEMENTATION: Sphere-sphere interaction ===
    if (isSpherical(pair.shape_i) && isSpherical(pair.shape_j)) {
        // Simple sphere interaction - direct use of sigma parameter
        sigma_effective = pair.sigma;
    }
    // === FUTURE IMPLEMENTATION: Ellipsoid interactions ===
    else {
        /*
        // TODO: Implement angle-dependent logic for ellipsoids here
        // This structure is prepared but not yet implemented (pragmatic approach)

        double contact_distance = calculateEffectiveDistance(
            pair.shape_i, pair.shape_j,
            pair.orient_i, pair.orient_j,
            contact_vector
        );
        sigma_effective = pair.sigma * (contact_distance / 4.0);  // Scaling factor
        */

        // Temporary fallback: treat ellipsoids as equivalent spheres
        double avg_radius_i = (pair.shape_i[0] + pair.shape_i[1] + pair.shape_i[2]) / 3.0;
        double avg_radius_j = (pair.shape_j[0] + pair.shape_j[1] + pair.shape_j[2]) / 3.0;
        sigma_effective = pair.sigma * (avg_radius_i + avg_radius_j) / 4.0;
    }

    // Calculate energy based on potential type
    switch (pair.cg_potential_type) {
    case 1: // LJ (1,6,12) - for SCNP cross-linking prevention
        return calculateLJ_1612(distance, sigma_effective, pair.epsilon);

    case 2: // Standard LJ (6,12) - for general CG simulations
        return calculateLJ_612(distance, sigma_effective, pair.epsilon);

    default:
        return 0.0; // Unknown potential type
    }
}

} // namespace CGPotentials