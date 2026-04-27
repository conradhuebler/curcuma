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

class ConfigManager;  // Forward declaration — avoids circular include

/**
 * @brief Coordination Number Calculator
 *
 * Calculates coordination numbers using different methods:
 * - GFN-FF erf-based CN with configurable neighbor list cutoff (P2b)
 * - D3 exponential CN
 *
 * Claude Generated 2025 - Consolidated CN calculation utility
 * P2b (Apr 2026) - Added neighbor list mode and ConfigManager integration
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
     * Calculate GFN-FF coordination numbers (original threshold-based API)
     *
     * @param atoms Atomic numbers (1-based)
     * @param geometry_bohr Molecular geometry in Bohr
     * @param threshold Distance cutoff in Bohr² (default: 1600.0 = 40²)
     * @param kn Steepness parameter (default: -7.5)
     * @param cnmax Max CN for log compression (default: 4.4)
     * @return Vector of coordination numbers (one per atom, max 4.4)
     */
    static std::vector<double> calculateGFNFFCN(
        const std::vector<int>& atoms,
        const Eigen::MatrixXd& geometry_bohr,
        double threshold = 1600.0,
        double kn = -7.5,
        double cnmax = 4.4
    );

    /**
     * Calculate GFN-FF coordination numbers with configurable cutoff mode (P2b)
     *
     * Three modes selected by parameters:
     *   cn_cutoff_bohr > 0: Neighbor-list mode — O(N*k) where k is avg neighbors
     *   cn_cutoff_bohr = 0, cn_accuracy > 0: Fortran accuracy-based threshold
     *     cnthr = 100 - log10(acc)*50, threshold = cnthr Bohr²
     *   cn_cutoff_bohr = 0, cn_accuracy = 0: Full O(N²) reference mode (threshold = inf)
     *
     * @param atoms Atomic numbers (1-based)
     * @param geometry_bohr Molecular geometry in Bohr
     * @param cn_cutoff_bohr Neighbor list cutoff in Bohr (>0: neighbor list, 0: threshold mode)
     * @param cn_accuracy Accuracy for threshold mode (only used when cn_cutoff_bohr = 0)
     * @param kn Steepness parameter (default: -7.5)
     * @param cnmax Max CN for log compression (default: 4.4)
     * @return Vector of coordination numbers (one per atom)
     */
    static std::vector<double> calculateGFNFFCN(
        const std::vector<int>& atoms,
        const Eigen::MatrixXd& geometry_bohr,
        double cn_cutoff_bohr,
        double cn_accuracy,
        double kn,
        double cnmax
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
    // Data source: s-dftd3 atomic radii - 86 elements (H through Rn)
    static const std::vector<double> COVALENT_RADII;
};
