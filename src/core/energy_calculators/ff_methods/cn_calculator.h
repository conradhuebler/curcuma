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

#include "src/core/global.h"
#include "src/core/math_compat.h"   // curcuma_erf / curcuma_exp (portable-math aware)

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
     * @param cn_cutoff_bohr Real-space cutoff for the CN sum in Bohr
     *        (default: 25.0, matching tblite's d3 dispersion container,
     *        tblite/disp/d3.f90: realspace_cutoff(cn=25.0)). <=0 disables.
     *        NOTE: the s-dftd3 *standalone* default is 40.0; tblite overrides it.
     * @return Vector of coordination numbers (one per atom)
     */
    static std::vector<double> calculateD3CN(
        const std::vector<int>& atoms,
        const Eigen::MatrixXd& geometry,
        double k1 = 16.0,
        double k2 = 4.0 / 3.0,
        double cn_cutoff_bohr = 25.0
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
        double cnmax = 4.4,
        std::vector<double>* out_cn_raw = nullptr  ///< WP-D (May 2026): optional raw CN output for reuse in dcn calculation
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

    /// WP-D Stage C (May 2026): CN + cn_raw + symmetric neighbor list in one pass.
    /// Returns all three so dcn can skip its own N²-erf loop AND its O(N²) pair scan.
    struct CNResult {
        std::vector<double> cn_values;           ///< post-log CN (size N)
        std::vector<double> cn_raw;              ///< pre-log raw erf-sum (size N)
        std::vector<std::vector<int>> neighbors; ///< symmetric: neighbors[i] = all j within cutoff
        double cutoff_sq = 0.0;                  ///< cutoff² in Bohr²
    };

    /// WP-D Stage C: compute CN values, raw CN, and neighbor list in a single O(N²) pass.
    static CNResult calculateGFNFFCNWithNeighbors(
        const std::vector<int>& atoms,
        const Eigen::MatrixXd& geometry_bohr,
        double cn_cutoff_bohr,
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

    /// rev-gfnff stage 3a(i) (Claude Generated, Sep 2026): the radii and the pair term the CN
    /// itself is built from, exposed so that a caller which needs "what would this pair
    /// contribute to my coordination number" gets the SAME numbers as the CN, not a second
    /// approximation. calculateGFNFFCN()/calculateGFNFFCNWithNeighbors() below use both.
    ///
    /// GFN-FF CN (Fortran gfnff_param.f90:551 / gfnff_rab.f90): the raw per-pair count is
    ///     c_ij(r) = 0.5 * (1 + erf(kn * (r - R)/R)),   R = rcov_i + rcov_j,  kn = -7.5
    /// with r, R in Bohr and rcov = 4/3 * covalent_radii(Z) converted to Bohr, and the atom's
    /// CN is the log-compressed sum of c_ij over its neighbours.

    /// rcov of one atom in Bohr, exactly as the CN build computes it. An element outside the
    /// table keeps the 0.0 sentinel the CN build uses to exclude the atom.
    static double gfnffCNRadiusBohr(int atomic_number)
    {
        constexpr double k_scaled = 4.0 / 3.0;      // gfnff_param.f90 covalentRadD3 scaling
        constexpr double ANG2BOHR = 1.8897259886;   // as in calculateGFNFFCN()
        const int idx = atomic_number - 1;
        if (idx < 0 || idx >= static_cast<int>(COVALENT_RADII.size()))
            return 0.0;
        return k_scaled * COVALENT_RADII[idx] * ANG2BOHR;
    }

    /// this pair's own contribution to the raw CN of either atom (symmetric in i, j)
    static double pairCNContribution(double r_bohr, double rcov_sum_bohr, double kn = -7.5)
    {
        const double dr = (r_bohr - rcov_sum_bohr) / rcov_sum_bohr;
        return 0.5 * (1.0 + curcuma_erf(kn * dr));
    }

    /// d c_ij / dr of the same expression (Bohr^-1):  erf'(kn*dr) * kn / R
    static double pairCNContributionDerivative(double r_bohr, double rcov_sum_bohr, double kn = -7.5)
    {
        constexpr double INV_SQRTPI = 0.56418958354775628695;  // 1/sqrt(pi)
        const double dr = (r_bohr - rcov_sum_bohr) / rcov_sum_bohr;
        return (kn * INV_SQRTPI) * curcuma_exp(-(kn * dr) * (kn * dr)) / rcov_sum_bohr;
    }

    /**
     * Add D3 CN chain-rule gradient to Cartesian gradient matrix.
     * Claude Generated (May 2026): Analytical gradient for D3 dispersion CN dependence.
     *
     * Uses the D3 exponential counting function:
     *   CN_i = sum_j 1/(1+exp(-k1*(k2*(R_cov_i+R_cov_j)/R_ij - 1)))
     *
     * @param atoms Atomic numbers (1-based)
     * @param geometry Molecular geometry in Ångström
     * @param dEdcn Per-atom dE/dCN values [N]
     * @param gradient_out [N,3] Cartesian gradient in Eh/Angstrom (accumulated)
     * @param k1 Steepness parameter (default: 16.0)
     * @param k2 Scaling factor (default: 4.0/3.0)
     * @param cn_cutoff_bohr Real-space cutoff for the CN sum in Bohr
     *        (default: 25.0, matching tblite's d3 dispersion container). <=0
     *        disables. Must match calculateD3CN so energy and gradient agree.
     */
    static void addD3CNGradient(
        const std::vector<int>& atoms,
        const Matrix& geometry,
        const Vector& dEdcn,
        Matrix& gradient_out,
        double k1 = 16.0,
        double k2 = 4.0 / 3.0,
        double distance_unit_to_bohr = 1.0,  // multiply output by this (e.g. au=1.8897)
        double cn_cutoff_bohr = 25.0
    );

private:
    // Covalent radii from simple-dftd3 (Angstrom)
    // Data source: s-dftd3 atomic radii - 86 elements (H through Rn)
    static const std::vector<double> COVALENT_RADII;
};
