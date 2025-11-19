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

#include "gbsa.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace Curcuma {
namespace Solvation {

// Physical constants
namespace {
    const double ANGSTROM_TO_BOHR = 1.8897259886;     // Conversion factor
    const double KCAL_TO_HARTREE = 0.00159360144;      // kcal/mol to Hartree
    const double SQRT_PI = 1.7724538509055159;         // √π
}

double GBSA::calculateEnergy(
    const std::vector<int>& atomic_numbers,
    const std::vector<std::array<double, 3>>& positions,
    const std::vector<double>& charges
) {
    const size_t natoms = atomic_numbers.size();

    if (positions.size() != natoms || charges.size() != natoms) {
        throw std::invalid_argument("GBSA::calculateEnergy: Inconsistent array sizes");
    }

    // Initialize result arrays
    m_born_radii.resize(natoms);
    m_sasa.resize(natoms);

    // Step 1: Calculate Born radii (GBOBC-II)
    calculateBornRadii(atomic_numbers, positions);

    // Step 2: Calculate solvent-accessible surface areas
    calculateSASA(atomic_numbers, positions);

    // Step 3: Calculate GB electrostatic energy
    m_energy_gb = calculateBornEnergy(charges, positions);

    // Step 4: Calculate surface area (hydrophobic) energy
    // E_SA = γ · Σ SASA_i
    // γ in kcal/(mol·Ų), convert to Hartree
    double total_sasa = 0.0;
    for (double sasa : m_sasa) {
        total_sasa += sasa;
    }
    m_energy_sa = m_surface_tension * total_sasa * KCAL_TO_HARTREE;

    // Total solvation energy
    return m_energy_gb + m_energy_sa;
}

void GBSA::calculateBornRadii(
    const std::vector<int>& atomic_numbers,
    const std::vector<std::array<double, 3>>& positions
) {
    const size_t natoms = atomic_numbers.size();

    // GBOBC-II parameters (Onufriev et al. 2004)
    const double alpha_LR = 1.0;
    const double beta_LR = 0.8;
    const double gamma_LR = 4.85;
    const double born_offset = BornOffset();  // 0.09 Angstrom

    // Initialize descreening integrals Ψ
    std::vector<double> psi(natoms, 0.0);

    // Step 1: Calculate pairwise descreening integrals
    for (size_t i = 0; i < natoms; ++i) {
        const int atom_i = atomic_numbers[i];
        const double vdW_i = vdWaalsRadius(atom_i);
        const double rho_i = vdW_i * DeScreening(m_solvent, atom_i);

        for (size_t j = 0; j < i; ++j) {
            const int atom_j = atomic_numbers[j];

            // Calculate distance r_ij
            double r_ij_sq = 0.0;
            for (int k = 0; k < 3; ++k) {
                const double delta = positions[i][k] - positions[j][k];
                r_ij_sq += delta * delta;
            }

            // Apply long-range cutoff
            if (r_ij_sq > m_long_range_cutoff * m_long_range_cutoff) {
                continue;
            }

            const double r_ij = std::sqrt(r_ij_sq);
            const double vdW_j = vdWaalsRadius(atom_j);
            const double rho_j = vdW_j * DeScreening(m_solvent, atom_j);

            // Check for overlaps
            const bool AB_overlap = (r_ij < vdW_i + rho_j);
            const bool BA_overlap = (r_ij < vdW_j + rho_i);

            double integral_i = 0.0;  // I_ij contribution to atom i
            double integral_j = 0.0;  // I_ji contribution to atom j

            if (!AB_overlap && !BA_overlap) {
                // Neither sphere overlaps
                integral_i = integralNonOverlap(r_ij, rho_j);
                integral_j = integralNonOverlap(r_ij, rho_i);
            } else if (!AB_overlap && BA_overlap) {
                // j overlaps i, but i doesn't overlap j
                integral_i = integralNonOverlap(r_ij, rho_j);
                if (r_ij + rho_i > vdW_j) {
                    integral_j = integralOverlap(r_ij, rho_i, vdW_j);
                }
            } else if (AB_overlap && !BA_overlap) {
                // i overlaps j, but j doesn't overlap i
                if (r_ij + rho_j > vdW_i) {
                    integral_i = integralOverlap(r_ij, rho_j, vdW_i);
                }
                integral_j = integralNonOverlap(r_ij, rho_i);
            } else {
                // Both spheres overlap
                if (r_ij + rho_j > vdW_i) {
                    integral_i = integralOverlap(r_ij, rho_j, vdW_i);
                }
                if (r_ij + rho_i > vdW_j) {
                    integral_j = integralOverlap(r_ij, rho_i, vdW_j);
                }
            }

            // Accumulate integrals
            psi[i] += integral_i;
            psi[j] += integral_j;
        }
    }

    // Step 2: Convert integrals to Born radii
    for (size_t i = 0; i < natoms; ++i) {
        const int atom_i = atomic_numbers[i];
        const double vdW_i = vdWaalsRadius(atom_i);
        const double rho_tilde_i = vdW_i - born_offset;  // Modified vdW radius

        // Scale and invert
        const double inv_vdW = 1.0 / vdW_i;
        const double inv_rho_tilde = 1.0 / rho_tilde_i;

        // Ψ_i scaled by ρ̃_i
        const double psi_scaled = 0.5 * psi[i] * rho_tilde_i;

        // Non-linear scaling function: tanh(α·Ψ + β·Ψ² - γ·Ψ³)
        const double psi_arg = alpha_LR * psi_scaled +
                               beta_LR * psi_scaled * psi_scaled -
                               gamma_LR * psi_scaled * psi_scaled * psi_scaled;
        const double tanh_psi = std::tanh(psi_arg);

        // Born radius: R_i = 1/(1/ρ̃_i - (1/ρ_i)·tanh(...))
        // Apply solvent-specific scaling
        const double born_radius_scale = BornRadiusScale(m_solvent);
        m_born_radii[i] = born_radius_scale / (inv_rho_tilde - inv_vdW * tanh_psi);
    }
}

void GBSA::calculateSASA(
    const std::vector<int>& atomic_numbers,
    const std::vector<std::array<double, 3>>& positions
) {
    const size_t natoms = atomic_numbers.size();
    const double probe_radius = 1.4;  // Water probe radius (Angstrom)

    // Get Lebedev grid for integration
    auto grid = LebedevGrid::generate110();

    // Calculate SASA for each atom
    for (size_t i = 0; i < natoms; ++i) {
        const int atom_i = atomic_numbers[i];
        const double vdW_i = vdWaalsRadius(atom_i);
        const double r_probe_i = vdW_i + probe_radius;  // Effective radius
        const double r_probe_i_sq = r_probe_i * r_probe_i;

        double sasa_sum = 0.0;

        // Integrate over grid points
        for (const auto& grid_point : grid) {
            // Grid point in real space: center_i + r_probe_i * grid_direction
            std::array<double, 3> probe_pos;
            for (int k = 0; k < 3; ++k) {
                probe_pos[k] = positions[i][k] + r_probe_i *
                               (k == 0 ? grid_point.x : (k == 1 ? grid_point.y : grid_point.z));
            }

            // Check overlap with other atoms using smooth switching function
            double overlap_factor = 1.0;

            for (size_t j = 0; j < natoms; ++j) {
                if (i == j) continue;

                const int atom_j = atomic_numbers[j];
                const double vdW_j = vdWaalsRadius(atom_j);
                const double r_probe_j = vdW_j + probe_radius;

                // Distance from probe position to atom j center
                double dist_sq = 0.0;
                for (int k = 0; k < 3; ++k) {
                    const double delta = probe_pos[k] - positions[j][k];
                    dist_sq += delta * delta;
                }
                const double dist = std::sqrt(dist_sq);

                // Smooth switching function to avoid discontinuities
                // S(d) = 1 if d > r_j, 0 if d < r_j, smooth transition in between
                const double switch_width = 0.5;  // Angstrom
                if (dist < r_probe_j - switch_width) {
                    overlap_factor = 0.0;
                    break;  // Fully buried
                } else if (dist < r_probe_j + switch_width) {
                    // Smooth interpolation
                    const double x = (dist - r_probe_j) / switch_width;
                    overlap_factor *= 0.5 * (1.0 + x);  // Linear for simplicity
                }
            }

            // Add weighted contribution
            sasa_sum += grid_point.weight * overlap_factor;
        }

        // SASA = (radius)² · 4π · (grid_sum / 4π) = radius² · grid_sum
        m_sasa[i] = r_probe_i_sq * sasa_sum;
    }
}

double GBSA::calculateBornEnergy(
    const std::vector<double>& charges,
    const std::vector<std::array<double, 3>>& positions
) {
    const size_t natoms = charges.size();

    // Prefactor: -½ · (1/ε_in - 1/ε_out)
    // ε_in = 1 (vacuum inside molecule), ε_out = ε_solvent
    const double epsilon_factor = -0.5 * (1.0 - 1.0 / m_epsilon);

    // Still et al. parameter for effective Coulomb operator
    const double cc = 4.0;  // Empirical parameter

    double energy = 0.0;

    // Self-energy terms (i = j)
    for (size_t i = 0; i < natoms; ++i) {
        const double q_i = charges[i];
        const double R_i = m_born_radii[i];

        // Self-energy: q_i² / R_i
        energy += q_i * q_i / R_i;
    }

    // Pairwise interaction terms (i ≠ j)
    for (size_t i = 0; i < natoms; ++i) {
        const double q_i = charges[i];
        const double R_i = m_born_radii[i];

        for (size_t j = 0; j < i; ++j) {
            const double q_j = charges[j];
            const double R_j = m_born_radii[j];

            // Distance r_ij
            double r_ij_sq = 0.0;
            for (int k = 0; k < 3; ++k) {
                const double delta = positions[i][k] - positions[j][k];
                r_ij_sq += delta * delta;
            }
            const double r_ij = std::sqrt(r_ij_sq);

            // Effective Coulomb operator (Still et al. formulation)
            // f_GB = √(r² + R_i·R_j·exp(-r²/(cc·R_i·R_j)))
            const double R_prod = R_i * R_j;
            const double exp_arg = -r_ij_sq / (cc * R_prod);
            const double f_GB = std::sqrt(r_ij_sq + R_prod * std::exp(exp_arg));

            // Interaction energy: 2 · q_i·q_j / f_GB  (factor 2 for i<j symmetry)
            energy += 2.0 * q_i * q_j / f_GB;
        }
    }

    // Apply prefactor and convert to Hartree
    // Energy is in e²/Angstrom, convert to Hartree
    const double COULOMB_ANGSTROM_TO_HARTREE = 1.0 / ANGSTROM_TO_BOHR;  // 1/a₀
    return epsilon_factor * energy * COULOMB_ANGSTROM_TO_HARTREE;
}

double GBSA::integralNonOverlap(double r_AB, double rho) const {
    // Non-overlapping case: simple 1/r⁴ integral
    // I = (1/rho³) · integral formula

    const double r_inv = 1.0 / r_AB;
    const double rho_inv = 1.0 / rho;

    // Analytical formula from Onufriev et al. 2004
    // I = (1/4) · (1/rho³) · [(r - rho)/(r + rho)]
    const double numerator = r_AB - rho;
    const double denominator = r_AB + rho;

    if (denominator < 1.0e-8) {
        return 0.0;
    }

    const double ratio = numerator / denominator;
    const double log_term = std::log(ratio);

    // Full formula includes additional terms for proper integration
    // Simplified version for quick-win implementation
    return 0.25 * rho_inv * rho_inv * rho_inv * ratio;
}

double GBSA::integralOverlap(double r_AB, double rho, double vdW_radius) const {
    // Overlapping case: more complex integral with cutoff at vdW surface

    const double r_inv = 1.0 / r_AB;
    const double vdW_inv = 1.0 / vdW_radius;

    // Analytical formula for overlap region
    // Includes geometric corrections for sphere intersection
    const double upper_limit = vdW_radius;
    const double lower_limit = std::abs(r_AB - rho);

    if (upper_limit <= lower_limit) {
        return 0.0;
    }

    // Simplified overlap integral
    // Full implementation would include complete Onufriev et al. 2004 formulas
    const double delta = upper_limit - lower_limit;
    return 0.5 * delta / (rho * rho * rho);
}

// =================================================================================
// Gradient Calculations (Claude Generated: November 2025)
// =================================================================================

double GBSA::calculateEnergyAndGradients(
    const std::vector<int>& atomic_numbers,
    const std::vector<std::array<double, 3>>& positions,
    const std::vector<double>& charges,
    std::vector<std::array<double, 3>>& gradients
) {
    const size_t natoms = atomic_numbers.size();

    if (positions.size() != natoms || charges.size() != natoms) {
        throw std::invalid_argument("GBSA::calculateEnergyAndGradients: Inconsistent array sizes");
    }

    // Ensure gradients array is properly sized
    if (gradients.size() != natoms) {
        gradients.resize(natoms);
    }

    // Initialize result arrays
    m_born_radii.resize(natoms);
    m_sasa.resize(natoms);
    m_born_radii_gradients.resize(natoms);
    m_sasa_gradients.resize(natoms);

    // Clear gradient caches
    for (size_t i = 0; i < natoms; ++i) {
        m_born_radii_gradients[i] = {0.0, 0.0, 0.0};
        m_sasa_gradients[i] = {0.0, 0.0, 0.0};
    }

    // Step 1: Calculate Born radii with gradients
    calculateBornRadiiWithGradients(atomic_numbers, positions);

    // Step 2: Calculate SASA with gradients
    calculateSASAWithGradients(atomic_numbers, positions);

    // Step 3: Calculate GB energy and gradients
    m_energy_gb = calculateBornEnergyWithGradients(charges, positions, gradients);

    // Step 4: Calculate SA energy and gradients
    // ∂E_SA/∂r = γ · ∂SASA/∂r
    double total_sasa = 0.0;
    for (size_t i = 0; i < natoms; ++i) {
        total_sasa += m_sasa[i];

        // Add SASA gradient contribution
        for (int k = 0; k < 3; ++k) {
            gradients[i][k] += m_surface_tension * m_sasa_gradients[i][k] * KCAL_TO_HARTREE;
        }
    }
    m_energy_sa = m_surface_tension * total_sasa * KCAL_TO_HARTREE;

    // Total solvation energy
    return m_energy_gb + m_energy_sa;
}

void GBSA::calculateBornRadiiWithGradients(
    const std::vector<int>& atomic_numbers,
    const std::vector<std::array<double, 3>>& positions
) {
    // Simplified gradient calculation for Born radii
    // Full implementation would require derivatives of all integrals
    // Quick-win: Finite difference approximation for now

    const size_t natoms = atomic_numbers.size();
    const double h = 1.0e-5;  // Finite difference step (Angstrom)

    // First, calculate Born radii at current positions
    calculateBornRadii(atomic_numbers, positions);
    std::vector<double> R0 = m_born_radii;

    // Calculate numerical derivatives using finite differences
    // ∂R_i/∂x_j ≈ (R_i(x_j + h) - R_i(x_j - h)) / (2h)
    for (size_t j = 0; j < natoms; ++j) {
        for (int coord = 0; coord < 3; ++coord) {
            // Perturb coordinate
            std::vector<std::array<double, 3>> pos_plus = positions;
            std::vector<std::array<double, 3>> pos_minus = positions;
            pos_plus[j][coord] += h;
            pos_minus[j][coord] -= h;

            // Calculate Born radii at perturbed positions
            calculateBornRadii(atomic_numbers, pos_plus);
            std::vector<double> R_plus = m_born_radii;

            calculateBornRadii(atomic_numbers, pos_minus);
            std::vector<double> R_minus = m_born_radii;

            // Numerical derivative for all atoms
            for (size_t i = 0; i < natoms; ++i) {
                double dR_dr = (R_plus[i] - R_minus[i]) / (2.0 * h);
                m_born_radii_gradients[i][coord] += dR_dr;
            }
        }
    }

    // Restore original Born radii
    m_born_radii = R0;
}

void GBSA::calculateSASAWithGradients(
    const std::vector<int>& atomic_numbers,
    const std::vector<std::array<double, 3>>& positions
) {
    // Simplified SASA gradient calculation
    // Quick-win: Finite difference approximation

    const size_t natoms = atomic_numbers.size();
    const double h = 1.0e-5;  // Finite difference step (Angstrom)

    // Calculate SASA at current positions
    calculateSASA(atomic_numbers, positions);
    std::vector<double> A0 = m_sasa;

    // Calculate numerical derivatives
    for (size_t j = 0; j < natoms; ++j) {
        for (int coord = 0; coord < 3; ++coord) {
            // Perturb coordinate
            std::vector<std::array<double, 3>> pos_plus = positions;
            std::vector<std::array<double, 3>> pos_minus = positions;
            pos_plus[j][coord] += h;
            pos_minus[j][coord] -= h;

            // Calculate SASA at perturbed positions
            calculateSASA(atomic_numbers, pos_plus);
            std::vector<double> A_plus = m_sasa;

            calculateSASA(atomic_numbers, pos_minus);
            std::vector<double> A_minus = m_sasa;

            // Numerical derivative for all atoms
            for (size_t i = 0; i < natoms; ++i) {
                double dA_dr = (A_plus[i] - A_minus[i]) / (2.0 * h);
                m_sasa_gradients[i][coord] += dA_dr;
            }
        }
    }

    // Restore original SASA
    m_sasa = A0;
}

double GBSA::calculateBornEnergyWithGradients(
    const std::vector<double>& charges,
    const std::vector<std::array<double, 3>>& positions,
    std::vector<std::array<double, 3>>& gradients
) {
    const size_t natoms = charges.size();

    // Prefactor: -½ · (1/ε_in - 1/ε_out)
    const double epsilon_factor = -0.5 * (1.0 - 1.0 / m_epsilon);
    const double cc = 4.0;  // Still et al. parameter

    double energy = 0.0;

    // Self-energy terms and gradients
    for (size_t i = 0; i < natoms; ++i) {
        const double q_i = charges[i];
        const double R_i = m_born_radii[i];

        // Self-energy: q_i² / R_i
        energy += q_i * q_i / R_i;

        // Self-energy gradient: -q_i² / R_i² · ∂R_i/∂r
        const double dE_dR = -q_i * q_i / (R_i * R_i);
        for (int k = 0; k < 3; ++k) {
            gradients[i][k] += epsilon_factor * dE_dR * m_born_radii_gradients[i][k]
                             * (1.0 / ANGSTROM_TO_BOHR);
        }
    }

    // Pairwise interaction terms and gradients
    for (size_t i = 0; i < natoms; ++i) {
        const double q_i = charges[i];
        const double R_i = m_born_radii[i];

        for (size_t j = 0; j < i; ++j) {
            const double q_j = charges[j];
            const double R_j = m_born_radii[j];

            // Distance and direction
            std::array<double, 3> r_ij;
            double r_ij_sq = 0.0;
            for (int k = 0; k < 3; ++k) {
                r_ij[k] = positions[i][k] - positions[j][k];
                r_ij_sq += r_ij[k] * r_ij[k];
            }
            const double r = std::sqrt(r_ij_sq);

            // Effective Coulomb operator
            const double R_prod = R_i * R_j;
            const double exp_arg = -r_ij_sq / (cc * R_prod);
            const double exp_term = std::exp(exp_arg);
            const double f_GB = std::sqrt(r_ij_sq + R_prod * exp_term);

            // Interaction energy
            const double E_pair = 2.0 * q_i * q_j / f_GB;
            energy += E_pair;

            // Gradient w.r.t. distance: ∂E/∂r_ij
            const double df_dr = (r - R_prod * exp_term * r / (cc * R_prod)) / f_GB;
            const double dE_dr = -E_pair * df_dr / f_GB;

            // Add distance-dependent gradient contribution
            for (int k = 0; k < 3; ++k) {
                const double grad_component = dE_dr * r_ij[k] / r;
                gradients[i][k] += epsilon_factor * grad_component * (1.0 / ANGSTROM_TO_BOHR);
                gradients[j][k] -= epsilon_factor * grad_component * (1.0 / ANGSTROM_TO_BOHR);
            }

            // Gradient w.r.t. Born radii: ∂E/∂R_i and ∂E/∂R_j
            const double dE_dRi = -E_pair * R_j * exp_term * (1.0 - r_ij_sq / (cc * R_prod)) / (2.0 * f_GB * f_GB);
            const double dE_dRj = -E_pair * R_i * exp_term * (1.0 - r_ij_sq / (cc * R_prod)) / (2.0 * f_GB * f_GB);

            for (int k = 0; k < 3; ++k) {
                gradients[i][k] += epsilon_factor * dE_dRi * m_born_radii_gradients[i][k]
                                 * (1.0 / ANGSTROM_TO_BOHR);
                gradients[j][k] += epsilon_factor * dE_dRj * m_born_radii_gradients[j][k]
                                 * (1.0 / ANGSTROM_TO_BOHR);
            }
        }
    }

    // Apply prefactor and convert to Hartree
    return epsilon_factor * energy * (1.0 / ANGSTROM_TO_BOHR);
}

} // namespace Solvation
} // namespace Curcuma
