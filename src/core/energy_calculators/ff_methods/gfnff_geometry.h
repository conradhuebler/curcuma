/*
 * <GFN-FF Geometry Functions for Force Field Threads>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * =============================================================================
 * ACKNOWLEDGMENT OF ORIGINAL WORK
 * =============================================================================
 *
 * This implementation is based on the GFN-FF force field method developed by:
 *   - Prof. Stefan Grimme (University of Bonn)
 *   - Dr. Sebastian Spicher (University of Bonn)
 *
 * Original Publication:
 *   S. Spicher, S. Grimme
 *   "Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems"
 *   Angew. Chem. Int. Ed. 59, 15665-15676 (2020)
 *   DOI: 10.1002/anie.202004239
 *
 * Original Fortran Implementation:
 *   external/gfnff/src/gfnff_helpers.f90 (valijklff, omega, dphidr functions)
 *   Developed by: Grimme Group (Sebastian Ehlert, Sebastian Spicher, Stefan Grimme)
 *
 * Claude Generated (2025): Standalone geometry functions for multi-threaded force field
 */

#pragma once

#include "src/core/global.h"
#include <Eigen/Dense>
#include <cmath>

namespace GFNFF_Geometry {

/**
 * @brief Calculate dihedral angle for four atoms i-j-k-l (GFN-FF method)
 *
 * Reference: external/gfnff/src/gfnff_helpers.f90:319-358 (valijklff function)
 *           src/core/energy_calculators/qm_methods/gfnff_torsions.cpp:102-167
 *
 * Physical Interpretation:
 *   φ = 0°:   cis/eclipsed   (i and l on same side)
 *   φ = 180°: trans/anti     (i and l on opposite sides)
 *   φ = ±60°: gauche         (staggered conformations)
 *
 * Mathematical Formula:
 *   φ = atan2(|v₂|·(n1·v3), n1·n2)
 *   where n1 = v1 × v2 (normal to plane i-j-k)
 *         n2 = v2 × v3 (normal to plane j-k-l)
 *         v1 = r_j - r_i, v2 = r_k - r_j, v3 = r_l - r_k
 *
 * @param r_i Position of atom i (Eigen::Vector3d or row from Matrix)
 * @param r_j Position of atom j
 * @param r_k Position of atom k
 * @param r_l Position of atom l
 * @param gradient Output: Matrix(4,3) with ∂φ/∂r_i, ∂φ/∂r_j, ∂φ/∂r_k, ∂φ/∂r_l
 * @param calculate_gradient If true, compute analytical gradients
 * @return Dihedral angle in radians [-π, π]
 *
 * Claude Generated (2025): Extracted from GFNFF class for thread-safe usage
 */
inline double calculateDihedralAngle(
    const Eigen::Vector3d& r_i,
    const Eigen::Vector3d& r_j,
    const Eigen::Vector3d& r_k,
    const Eigen::Vector3d& r_l,
    Matrix& gradient,
    bool calculate_gradient = false)
{
    // Calculate bond vectors
    Eigen::Vector3d v1 = r_j - r_i;  // i→j bond vector
    Eigen::Vector3d v2 = r_k - r_j;  // j→k bond vector (central bond)
    Eigen::Vector3d v3 = r_l - r_k;  // k→l bond vector

    // Calculate cross products (plane normals)
    Eigen::Vector3d n1 = v1.cross(v2);  // Normal to plane i-j-k
    Eigen::Vector3d n2 = v2.cross(v3);  // Normal to plane j-k-l

    // Calculate magnitudes
    double n1_norm = n1.norm();
    double n2_norm = n2.norm();
    double v2_norm = v2.norm();

    // Check for degenerate geometries (colinear atoms)
    const double epsilon = 1.0e-10;
    if (n1_norm < epsilon || n2_norm < epsilon || v2_norm < epsilon) {
        // Colinear atoms → dihedral angle undefined, return 0
        if (calculate_gradient) {
            gradient = Matrix::Zero(4, 3);
        }
        return 0.0;
    }

    // Normalize the normal vectors
    Eigen::Vector3d n1_normalized = n1 / n1_norm;
    Eigen::Vector3d n2_normalized = n2 / n2_norm;

    // Calculate dihedral angle using atan2 for proper sign
    // Reference: standard atan2 definition for dihedrals
    // CRITICAL FIX (Jan 25, 2026): Correct sin_phi normalization.
    // Standard sin_phi = (v1 x v2) . (v2 x v3) x v2 / (|v1 x v2| * |v2 x v3| * |v2|)
    // In terms of normals: sin_phi = (n1 x n2) . v2_unit / (|n1| * |n2|)
    // Formula below: v2_norm * (n1 . v3) / (|n1| * |n2|) is equivalent.
    double cos_phi = n1_normalized.dot(n2_normalized);
    double sin_phi = (v2_norm * n1_normalized.dot(v3)) / n2_norm;

    // CRITICAL FIX (Feb 2026): Clamp sin_phi to [-1, 1] to avoid numerical edge cases
    // Prevents NaN in gradient calculations if sin_phi slightly exceeds bounds
    sin_phi = std::max(-1.0, std::min(1.0, sin_phi));

    // Use atan2 for correct quadrant (returns angle in [-π, π])
    double phi = atan2(sin_phi, cos_phi);

    if (!calculate_gradient) {
        return phi;
    }

    // =========================================================================
    // ANALYTICAL GRADIENT CALCULATION
    // =========================================================================
    // Reference: external/gfnff/src/gfnff_helpers.f90:514-583 (dphidr subroutine)
    //
    // Gradient formula (from Fortran):
    //   ∂φ/∂r = (1/(nan*nbn*sinφ)) * [cross product expressions]
    //
    // where nan = |n1|, nbn = |n2|

    gradient = Matrix::Zero(4, 3);

    double sin_phi_val = sin(phi);
    // Guard matching Fortran dphidr (gfnff_helpers.f90:546): eps=1.d-14
    // nenner = nan*nbn*sinphi — only guard against true numerical zero.
    // Previous threshold of 1e-4 was too aggressive: zeroed gradients for
    // many near-planar dihedrals (aromatic rings), losing physical contributions.
    double nenner = n1_norm * n2_norm * sin_phi_val;
    if (std::abs(nenner) < 1e-14) {
        gradient.setZero();
        return phi;
    }

    // Normalization factor
    double onenner = 1.0 / nenner;

    // Calculate various cross products needed for gradients
    Eigen::Vector3d rab = n1.cross(v2);  // n1 × v2
    Eigen::Vector3d rbb = n2.cross(v2);  // n2 × v2
    Eigen::Vector3d rac = n1.cross(v3);  // n1 × v3
    Eigen::Vector3d rbc = n2.cross(v3);  // n2 × v3
    Eigen::Vector3d rba = n2.cross(v1);  // n2 × v1
    Eigen::Vector3d raa = n1.cross(v1);  // n1 × v1

    Eigen::Vector3d rapb = v1 + v2;  // r_i-r_j + r_j-r_k
    Eigen::Vector3d rbpc = v2 + v3;  // r_j-r_k + r_k-r_l

    Eigen::Vector3d rapba = rapb.cross(n1);  // (v1+v2) × n1
    Eigen::Vector3d rapbb = rapb.cross(n2);  // (v1+v2) × n2
    Eigen::Vector3d rbpca = rbpc.cross(n1);  // (v2+v3) × n1
    Eigen::Vector3d rbpcb = rbpc.cross(n2);  // (v2+v3) × n2

    // Gradient components (Fortran dphidr formulas)
    // ∂φ/∂r_i
    gradient.row(0) = onenner * (cos_phi * n2_norm / n1_norm * rab - rbb);

    // ∂φ/∂r_j
    gradient.row(1) = onenner * (cos_phi * (n2_norm / n1_norm * rapba + n1_norm / n2_norm * rbc)
                                 - (rac + rapbb));

    // ∂φ/∂r_k
    gradient.row(2) = onenner * (cos_phi * (n2_norm / n1_norm * raa + n1_norm / n2_norm * rbpcb)
                                 - (rba + rbpca));

    // ∂φ/∂r_l
    gradient.row(3) = onenner * (cos_phi * n1_norm / n2_norm * rbb - rab);

    return phi;
}

/**
 * @brief Calculate out-of-plane angle (omega) for inversion term
 *
 * Reference: external/gfnff/src/gfnff_helpers.f90:427-448 (omega function)
 *           src/core/energy_calculators/qm_methods/gfnff_inversions.cpp
 *
 * Physical Interpretation:
 *   ω = 0:     atom i in plane j-k-l (planar, typical for sp²)
 *   ω = ±π/2:  atom i perpendicular to plane (pyramidal)
 *
 * Mathematical Formula:
 *   ω = arcsin(n · v̂)
 *   where n = (r_ij × r_ik) / |r_ij × r_ik|  (normal to plane i-j-k)
 *         v = r_il  (vector from i to l)
 *
 * @param r_i Position of central atom (out-of-plane)
 * @param r_j Position of first plane atom
 * @param r_k Position of second plane atom
 * @param r_l Position of third plane atom
 * @param gradient Output: Matrix(4,3) with ∂ω/∂r_i, ∂ω/∂r_j, ∂ω/∂r_k, ∂ω/∂r_l
 * @param calculate_gradient If true, compute analytical gradients
 * @return Out-of-plane angle in radians [-π/2, π/2]
 *
 * Claude Generated (2025): Extracted from GFNFF class for thread-safe usage
 */
inline double calculateOutOfPlaneAngle(
    const Eigen::Vector3d& r_i,
    const Eigen::Vector3d& r_j,
    const Eigen::Vector3d& r_k,
    const Eigen::Vector3d& r_l,
    Matrix& gradient,
    bool calculate_gradient = false)
{
    // Fortran vectors (gfnff_helpers.f90:436-440, 465-471)
    // re = center - nb1, rd = nb2 - nb1, rv = nb3 - center
    Eigen::Vector3d re = r_i - r_j;
    Eigen::Vector3d rd = r_k - r_j;
    Eigen::Vector3d rv = r_l - r_i;

    // Normal vector to plane (gfnff_helpers.f90:441)
    Eigen::Vector3d rn = re.cross(rd);
    double rnn = rn.norm();
    double rvn = rv.norm();

    const double epsilon = 1.0e-10;
    if (rnn < epsilon || rvn < epsilon) {
        if (calculate_gradient) {
            gradient = Matrix::Zero(4, 3);
        }
        return 0.0;
    }

    // omega = asin(rn_hat · rv_hat) — Fortran convention (gfnff_helpers.f90:442-446)
    // Sign-flipped vs previous C++ convention; energy unaffected (cos is even)
    Eigen::Vector3d rn_hat = rn / rnn;
    Eigen::Vector3d rv_hat = rv / rvn;
    double sin_omega = rn_hat.dot(rv_hat);
    sin_omega = std::max(-1.0, std::min(1.0, sin_omega));
    double omega = asin(sin_omega);

    if (!calculate_gradient) {
        return omega;
    }

    // =========================================================================
    // ANALYTICAL GRADIENT (gfnff_helpers.f90:450-510, domegadr subroutine)
    // =========================================================================
    // Exact port of Fortran domegadr — proven translation-invariant:
    //   Σ(dω/dr_i + dω/dr_j + dω/dr_k + dω/dr_l) = 0 (algebraically verified)

    gradient = Matrix::Zero(4, 3);

    double cos_omega_val = cos(omega);
    double nenner = rnn * rvn * cos_omega_val;
    if (std::abs(nenner) < 1e-14) {
        gradient.setZero();
        return omega;
    }

    double onenner = 1.0 / nenner;
    double sin_omega_val = sin(omega);

    Eigen::Vector3d rdme = rd - re;  // = r_k - r_i

    // Cross products (gfnff_helpers.f90:477-482)
    Eigen::Vector3d rve   = rv.cross(re);
    Eigen::Vector3d rne   = rn.cross(re);
    Eigen::Vector3d rdv   = rd.cross(rv);
    Eigen::Vector3d rdn   = rd.cross(rn);
    Eigen::Vector3d rvdme = rv.cross(rdme);
    Eigen::Vector3d rndme = rn.cross(rdme);

    // Gradient formulas (gfnff_helpers.f90:489-499)
    // row(0) = dω/dr_i (center), row(1) = dω/dr_j (nb1),
    // row(2) = dω/dr_k (nb2),    row(3) = dω/dr_l (nb3)
    gradient.row(0) = onenner * (rdv - rn - sin_omega_val * (rvn / rnn * rdn - rnn / rvn * rv));
    gradient.row(1) = onenner * (rvdme - sin_omega_val * rvn / rnn * rndme);
    gradient.row(2) = onenner * (rve   - sin_omega_val * rvn / rnn * rne);
    gradient.row(3) = onenner * (rn    - sin_omega_val * rnn / rvn * rv);

    return omega;
}

/**
 * @brief Calculate torsion damping function (GFN-FF specific)
 *
 * Reference: external/gfnff/src/gfnff_engrad.F90 (gfnffdampt function)
 *           src/core/energy_calculators/qm_methods/gfnff_torsions.cpp:193-278
 *
 * Physical meaning:
 *   D(r) = 1 / [1 + exp(-α*(r/r₀ - 1))]
 *
 *   - Stretched bonds (r > r₀) → D → 1 (weak damping, easier rotation)
 *   - Compressed bonds (r < r₀) → D → 0 (strong damping, harder rotation)
 *   - Equilibrium (r = r₀) → D = 0.5
 *
 * This couples bond stretching with torsional motion, representing the
 * physical reality that stretched bonds have lower torsion barriers.
 *
 * @param r_squared Squared distance r² in Angstrom²
 * @param r0 Equilibrium bond length in Angstrom
 * @param damp Output: Damping value D(r) ∈ [0, 1]
 * @param damp_deriv Output: Derivative dD/dr
 *
 * Claude Generated (2025): Standalone version for ForceFieldThread
 */
inline void calculateTorsionDamping(
    double r_squared,
    double r0,
    double& damp,
    double& damp_deriv)
{
    // Current distance
    double r = sqrt(r_squared);

    // Steepness parameter α (from Fortran gfnffdampt)
    // Typical value: 16-20 (controls how rapidly damping changes with distance)
    const double alpha = 16.0;

    // Calculate normalized distance deviation
    // x = r/r₀ - 1
    double x = (r / r0) - 1.0;

    // Exponential argument: -α·x
    // Capped at ±20 to avoid numerical overflow/underflow
    double exp_arg = -alpha * x;

    if (exp_arg > 20.0) {
        // Very negative x → strongly compressed bond
        // exp(-α·x) → ∞, so D → 0 (maximum damping)
        damp = 0.0;
        damp_deriv = 0.0;
        return;
    }

    if (exp_arg < -20.0) {
        // Very positive x → strongly stretched bond
        // exp(-α·x) → 0, so D → 1 (minimum damping)
        damp = 1.0;
        damp_deriv = 0.0;
        return;
    }

    // Calculate damping function
    // D(r) = 1 / [1 + exp(-α·x)]
    double exp_val = exp(exp_arg);
    damp = 1.0 / (1.0 + exp_val);

    // Calculate derivative dD/dr
    // dD/dr = α·D·(1 - D) / r₀
    damp_deriv = alpha * damp * (1.0 - damp) / r0;

    // Numerical stability check
    if (damp < 1.0e-10 || damp > (1.0 - 1.0e-10)) {
        damp_deriv = 0.0;
    }
}

} // namespace GFNFF_Geometry
