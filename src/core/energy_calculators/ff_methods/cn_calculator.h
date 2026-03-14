/*
 * Coordination Number Calculator for D3 Dispersion
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
 * Claude Generated 2025 - Consolidated CN calculation utility
 */

#pragma once

#include <vector>
#include <Eigen/Dense>

/**
 * @brief Coordination Number Calculator
 *
 * Calculates D3-style coordination numbers using exponential counting function.
 * Shared utility for D3, D4, and GFN-FF dispersion methods.
 *
 * Formula: CN_i = sum_{j≠i} 1 / (1 + exp(-k1 * (k2 * (R_cov_i + R_cov_j) / R_ij - 1)))
 * Reference: Grimme et al., J. Chem. Phys. 132, 154104 (2010)
 *
 * Claude Generated - December 20, 2025
 */
class CNCalculator {
public:
    /**
     * Calculate D3 coordination numbers
     *
     * @param atoms Atomic numbers (1-based: H=1, C=6, O=8, etc.)
     * @param geometry Molecular geometry in Ångström
     * @param k1 Steepness parameter (default: 16.0)
     * @param k2 Scaling factor (default: 4.0/3.0)
     * @return Vector of coordination numbers (one per atom)
     */
    static std::vector<double> calculateD3CN(
        const std::vector<int>& atoms,
        const Eigen::MatrixXd& geometry,
        double k1 = 16.0,
        double k2 = 4.0 / 3.0
    );

    /**
     * Calculate GFN-FF coordination numbers
     *
     * Purpose: EEQ parameter corrections and topology classification
     * Formula: Error function with logarithmic transformation
     *   CN_erf = 0.5 * (1 + erf(kn * (r_cov/r_ij - 1)))
     *   CN_log = log(1 + e^cnmax) - log(1 + e^(cnmax - CN_erf))
     *
     * @param atoms Atomic numbers (1-based: H=1, C=6, O=8, etc.)
     * @param geometry_bohr Molecular geometry in Bohr
     * @param threshold Distance cutoff in Bohr² (default: 1600.0)
     * @param kn Steepness parameter (default: -7.5)
     * @param cnmax Max CN for log compression (default: 4.4)
     * @return Vector of coordination numbers (one per atom, max 4.4)
     *
     * Reference: gfnff_cn.f90:66-126, Spicher & Grimme, J. Chem. Theory Comput. 2020
     * Claude Generated - December 21, 2025
     */
    static std::vector<double> calculateGFNFFCN(
        const std::vector<int>& atoms,
        const Eigen::MatrixXd& geometry_bohr,
        double threshold = 1600.0,
        double kn = -7.5,
        double cnmax = 4.4
    );

    /**
     * Get covalent radius for an element
     *
     * @param atomic_number 1-based atomic number (H=1, C=6, etc.)
     * @return Covalent radius in Ångström
     */
    static double getCovalentRadius(int atomic_number);

private:
    // Covalent radii from simple-dftd3 (Angstrom)
    // Data source: s-dftd3 atomic radii - 18 elements (H through Ar)
    static const std::vector<double> COVALENT_RADII;
};
