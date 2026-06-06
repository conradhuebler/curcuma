/*
 * D4 electronegativity-weighted covalent coordination number
 *
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License v3 (or later).
 *
 * Claude Generated 2026 — exact port of dftd4's get_coordination_number
 * (electronegativity-weighted covalent CN) used for the D4 C6 interpolation.
 *
 * This is the coordination number that tblite's GFN2-xTB feeds into its D4
 * dispersion model (tblite/disp/d4.f90: get_coordination_number(mol, …, rcov,
 * en, …)). It differs from the GFN-FF erf-CN (CNCalculator::calculateGFNFFCN)
 * in two ways that matter for polar molecules:
 *   1. a per-pair electronegativity-difference weight
 *        den = k4·exp(−(|EN_i−EN_j| + k5)² / k6)
 *   2. NO log-compression cap (the GFN-FF CN is capped at CNMAX=4.4; the dftd4
 *      D4 dispersion CN is uncapped — the cn_max=8 cut applies only to the
 *      separate EEQ CN, not this one).
 * Both effects vanish for homonuclear H₂ and grow for O–H/N–H/C–H, which is
 * exactly the residual pattern this module fixes.
 *
 * Reference: external/cpp-d4/src/dftd_ncoord.cpp (ncoord_d4 / dncoord_d4),
 *            E. Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019).
 */

#pragma once

#include "src/core/global.h"

#include <vector>

namespace curcuma::dispersion {

/**
 * @brief EN-weighted covalent coordination number (dftd4 GFN2 D4 reference).
 *
 * @param atoms      Atomic numbers (1-based Z), size N
 * @param geom_bohr  Coordinates in Bohr, N×3 (row per atom)
 * @param cutoff_bohr Radial cutoff for the pair sum (dftd4 default ≈ 30 Bohr)
 * @return Per-atom coordination number (size N), uncapped.
 */
std::vector<double> computeD4CovalentCN(const std::vector<int>& atoms,
                                        const Matrix& geom_bohr,
                                        double cutoff_bohr = 30.0);

/**
 * @brief Fold dE/dCN into a Cartesian gradient using this CN's own derivative.
 *
 * grad_out(a) += Σ_pairs (dEdcn_i + dEdcn_j) · d(countf_ij)/dr · (R_i − R_j)/r
 * (Hartree/Bohr). grad_out is N×3 and is *added to*, not overwritten. Use this
 * instead of routing dE_D4/dCN through the GFN2-Hamiltonian CN loop so the D4
 * dispersion gradient stays self-consistent with the D4 CN definition.
 *
 * @param atoms      Atomic numbers (1-based Z), size N
 * @param geom_bohr  Coordinates in Bohr, N×3 (must match computeD4CovalentCN)
 * @param dEdcn      Per-atom dE_D4/dCN (size N)
 * @param grad_out   N×3 Cartesian gradient accumulator (Hartree/Bohr)
 * @param cutoff_bohr Radial cutoff (must match computeD4CovalentCN)
 */
void addD4CovalentCNGradient(const std::vector<int>& atoms,
                             const Matrix& geom_bohr,
                             const Vector& dEdcn,
                             Matrix& grad_out,
                             double cutoff_bohr = 30.0);

}  // namespace curcuma::dispersion
