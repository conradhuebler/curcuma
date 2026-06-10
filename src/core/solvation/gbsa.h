/*
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This file is part of Curcuma - Native Solvation Module
 *
 * Extracted and adapted from Ulysses (Copyright (C) 2023- Filipe Menezes et al.)
 * Implemented for Curcuma by Claude (Anthropic AI Assistant)
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */

#ifndef GBSA_H
#define GBSA_H

#include <vector>
#include <string>
#include <cmath>

#include "solvation_parameters.h"
#include "lebedev_grid.h"

namespace Curcuma {
namespace Solvation {

/**
 * @brief Native GBSA (Generalized Born + Surface Area) solvation model
 *
 * Implements implicit solvation using the GBSA approach:
 *   E_solv = E_GB + E_SA
 *
 * Where:
 *   - E_GB: Generalized Born electrostatic solvation energy
 *   - E_SA: Solvent-accessible surface area (hydrophobic) term
 *
 * References:
 *   - GBOBC-II Model: Onufriev, Bashford, Case, Proteins 55, 383 (2004)
 *   - GBSA: Hawkins et al., J. Phys. Chem. 100, 19824 (1996)
 *   - Implementation: Ulysses QC.hpp calcBornRadii(), SolvationEnergy()
 *
 * Algorithm Flow:
 *   1. Calculate Born radii for each atom (GBOBC-II)
 *   2. Calculate solvent-accessible surface areas (SASA)
 *   3. Compute Born matrix (pairwise effective Coulomb interactions)
 *   4. Calculate solvation energy: E_GB + γ·SASA
 *
 * @note Currently implements energy-only (no gradients)
 * @todo Add analytical gradients for geometry optimization
 * @todo Extract full GFN2 solvation parameters from Ulysses
 */
class GBSA {
public:
    /**
     * @brief Construct GBSA solvation model
     *
     * @param solvent Solvent name (e.g., "water", "dmso", "acetone")
     * @param method QM method name (e.g., "GFN2", "GFN1") for parameter selection
     */
    GBSA(const std::string& solvent = "water", const std::string& method = "GFN2")
        : m_solvent(solvent)
        , m_method(method)
        , m_epsilon(SolventDielectricConstant(solvent, method))
        , m_surface_tension(SurfaceTension(solvent))
        , m_short_range_cutoff(5.0)   // Angstrom
        , m_long_range_cutoff(35.0)   // Angstrom
    {
    }

    /**
     * @brief Calculate GBSA solvation energy
     *
     * Main entry point for solvation energy calculation.
     *
     * @param atomic_numbers Atomic numbers (Z) for all atoms
     * @param positions Cartesian coordinates in Angstrom [Natoms × 3]
     * @param charges Atomic partial charges (from QM calculation)
     * @return Solvation energy in Hartree
     *
     * @note Positions must be in Angstrom, returns energy in Hartree
     */
    double calculateEnergy(
        const std::vector<int>& atomic_numbers,
        const std::vector<std::array<double, 3>>& positions,
        const std::vector<double>& charges
    );

    /**
     * @brief Calculate GBSA solvation energy and gradients
     *
     * Computes both energy and analytical gradients for geometry optimization.
     *
     * @param atomic_numbers Atomic numbers (Z) for all atoms
     * @param positions Cartesian coordinates in Angstrom [Natoms × 3]
     * @param charges Atomic partial charges (from QM calculation)
     * @param gradients Output gradients in Hartree/Angstrom [Natoms × 3]
     * @return Solvation energy in Hartree
     *
     * @note Gradients are added to existing values (not cleared)
     * @note Claude Generated: Analytical GBSA gradients (November 2025)
     */
    double calculateEnergyAndGradients(
        const std::vector<int>& atomic_numbers,
        const std::vector<std::array<double, 3>>& positions,
        const std::vector<double>& charges,
        std::vector<std::array<double, 3>>& gradients
    );

    /**
     * @brief Get Born radii for each atom
     *
     * Returns effective Born radii after last energy calculation.
     * Useful for analysis and debugging.
     *
     * @return Vector of Born radii in Angstrom
     */
    const std::vector<double>& getBornRadii() const { return m_born_radii; }

    /**
     * @brief Get solvent-accessible surface areas
     *
     * Returns SASA for each atom after last energy calculation.
     *
     * @return Vector of SASA values in Angstrom²
     */
    const std::vector<double>& getSASA() const { return m_sasa; }

    /**
     * @brief Get GB energy contribution
     *
     * @return Generalized Born energy in Hartree
     */
    double getGBEnergy() const { return m_energy_gb; }

    /**
     * @brief Get surface area energy contribution
     *
     * @return Surface area energy (γ·SASA) in Hartree
     */
    double getSAEnergy() const { return m_energy_sa; }

private:
    // Solvent properties
    std::string m_solvent;
    std::string m_method;
    double m_epsilon;           // Dielectric constant
    double m_surface_tension;   // γ for SASA term (kcal/mol/Ų)

    // Cutoff distances
    double m_short_range_cutoff;  // For neighbor lists (Angstrom)
    double m_long_range_cutoff;   // For Born radii calculation (Angstrom)

    // Cached results
    std::vector<double> m_born_radii;  // Effective Born radii (Angstrom)
    std::vector<double> m_sasa;        // Solvent-accessible surface areas (Ų)
    double m_energy_gb;                // GB energy contribution (Hartree)
    double m_energy_sa;                // SA energy contribution (Hartree)

    // Cached gradient data (for efficiency)
    std::vector<std::array<double, 3>> m_born_radii_gradients;  // ∂R_i/∂r [Natoms × 3]
    std::vector<std::array<double, 3>> m_sasa_gradients;        // ∂A_i/∂r [Natoms × 3]

    /**
     * @brief Calculate Born radii using GBOBC-II algorithm
     *
     * Implements the GBOBC-II (Generalized Born with Optimal Coulomb approximation, version II)
     * algorithm from Onufriev et al. 2004.
     *
     * Algorithm:
     *   1. Calculate pairwise descreening integrals I_ij
     *   2. Sum integrals to get Ψ_i = Σ I_ij for each atom
     *   3. Apply non-linear scaling: R_i = 1/(ρ̃_i⁻¹ - ρ_i⁻¹·tanh(α·Ψ_i + β·Ψ_i² - γ·Ψ_i³))
     *
     * Parameters (GBOBC-II):
     *   - α = 1.0
     *   - β = 0.8
     *   - γ = 4.85
     *
     * @param atomic_numbers Atomic numbers (Z) for all atoms
     * @param positions Cartesian coordinates in Angstrom [Natoms × 3]
     */
    void calculateBornRadii(
        const std::vector<int>& atomic_numbers,
        const std::vector<std::array<double, 3>>& positions
    );

    /**
     * @brief Calculate solvent-accessible surface areas using Lebedev quadrature
     *
     * Integrates over spherical grid around each atom to determine exposed surface.
     * Uses smooth overlap function to handle atom-atom overlaps.
     *
     * Algorithm:
     *   1. For each atom i with vdW radius r_i:
     *   2.   Generate Lebedev grid points on sphere of radius r_i + r_probe
     *   3.   For each grid point:
     *   4.     Check overlap with neighboring atoms
     *   5.     Apply smooth switching function
     *   6.   Integrate: SASA_i = 4π(r_i + r_probe)² · Σ w_k · S_k
     *
     * @param atomic_numbers Atomic numbers (Z) for all atoms
     * @param positions Cartesian coordinates in Angstrom [Natoms × 3]
     */
    void calculateSASA(
        const std::vector<int>& atomic_numbers,
        const std::vector<std::array<double, 3>>& positions
    );

    /**
     * @brief Calculate Born energy from Born matrix
     *
     * Computes the generalized Born electrostatic energy using effective
     * Coulomb operator with Born radii.
     *
     * E_GB = -½ · (1/ε_in - 1/ε_out) · Σ_i Σ_j q_i q_j / f_GB(r_ij, R_i, R_j)
     *
     * Where f_GB is the effective Coulomb operator:
     *   f_GB = √(r_ij² + R_i·R_j·exp(-r_ij²/(4·R_i·R_j)))  (Still et al. formulation)
     *
     * @param charges Atomic partial charges from QM calculation
     * @param positions Cartesian coordinates in Angstrom [Natoms × 3]
     * @return GB energy in Hartree
     */
    double calculateBornEnergy(
        const std::vector<double>& charges,
        const std::vector<std::array<double, 3>>& positions
    );

    /**
     * @brief Calculate overlap integral for Born radii (non-overlapping case)
     *
     * Helper function for GBOBC-II Born radii calculation.
     * Computes integral when atom spheres don't overlap.
     *
     * @param r_AB Distance between atoms A and B (Angstrom)
     * @param rho Descreening radius (Angstrom)
     * @return Integral contribution
     */
    double integralNonOverlap(double r_AB, double rho) const;

    /**
     * @brief Calculate overlap integral for Born radii (overlapping case)
     *
     * Helper function for GBOBC-II Born radii calculation.
     * Computes integral when atom spheres overlap.
     *
     * @param r_AB Distance between atoms A and B (Angstrom)
     * @param rho Descreening radius (Angstrom)
     * @param vdW_radius van der Waals radius (Angstrom)
     * @return Integral contribution
     */
    double integralOverlap(double r_AB, double rho, double vdW_radius) const;

    /**
     * @brief Calculate Born radii with gradients
     *
     * Extended version that computes both Born radii and their derivatives.
     *
     * @param atomic_numbers Atomic numbers (Z) for all atoms
     * @param positions Cartesian coordinates in Angstrom [Natoms × 3]
     *
     * @note Updates m_born_radii and m_born_radii_gradients
     * @note Claude Generated: Gradient implementation (November 2025)
     */
    void calculateBornRadiiWithGradients(
        const std::vector<int>& atomic_numbers,
        const std::vector<std::array<double, 3>>& positions
    );

    /**
     * @brief Calculate SASA with gradients
     *
     * Extended version that computes both SASA and their derivatives.
     *
     * @param atomic_numbers Atomic numbers (Z) for all atoms
     * @param positions Cartesian coordinates in Angstrom [Natoms × 3]
     *
     * @note Updates m_sasa and m_sasa_gradients
     * @note Claude Generated: Gradient implementation (November 2025)
     */
    void calculateSASAWithGradients(
        const std::vector<int>& atomic_numbers,
        const std::vector<std::array<double, 3>>& positions
    );

    /**
     * @brief Calculate Born energy with gradients
     *
     * Extended version that computes both energy and gradients.
     *
     * @param charges Atomic partial charges from QM calculation
     * @param positions Cartesian coordinates in Angstrom [Natoms × 3]
     * @param gradients Output gradients (added to existing values)
     * @return GB energy in Hartree
     *
     * @note Requires m_born_radii and m_born_radii_gradients to be calculated first
     * @note Claude Generated: Gradient implementation (November 2025)
     */
    double calculateBornEnergyWithGradients(
        const std::vector<double>& charges,
        const std::vector<std::array<double, 3>>& positions,
        std::vector<std::array<double, 3>>& gradients
    );
};

} // namespace Solvation
} // namespace Curcuma

#endif // GBSA_H
