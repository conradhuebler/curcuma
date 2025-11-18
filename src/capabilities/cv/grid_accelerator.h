/*
 * <Grid Acceleration for Multi-CV Metadynamics>
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
 * ---
 * Claude Generated (November 2025)
 *
 * GRID ACCELERATION FOR METADYNAMICS
 *
 * This class implements grid-based bias calculation for multi-dimensional
 * metadynamics, providing O(1) bias evaluation instead of O(N_gaussians).
 *
 * ALGORITHM:
 * ----------
 * 1. Create regular grid spanning CV space
 * 2. Precompute bias potential at all grid points
 * 3. Use multilinear interpolation for bias/gradient evaluation
 * 4. Update grid incrementally when new Gaussians deposited
 *
 * PERFORMANCE GAINS:
 * -----------------
 * - Without grid: O(N_gaussians) per MD step
 * - With grid: O(1) per MD step (after O(N_gaussians) grid setup)
 * - Break-even point: ~100-500 Gaussians (depends on dimensionality)
 *
 * MEMORY USAGE:
 * ------------
 * - 1D: N_bins points (negligible, ~1 KB for 100 points)
 * - 2D: N_bins^2 points (moderate, ~80 KB for 100×100 grid)
 * - 3D: N_bins^3 points (significant, ~8 MB for 100×100×100 grid)
 *
 * ACCURACY:
 * --------
 * - Grid spacing should be σ/4 to σ/2 for good interpolation accuracy
 * - Adaptive grid refinement in high-gradient regions (future feature)
 *
 * REFERENCES:
 * ----------
 * [1] Branduardi, D. et al. (2012). Metadynamics with adaptive Gaussians.
 *     J. Chem. Theory Comput. 8, 2247-2254.
 *     DOI: 10.1021/ct3002464
 *
 * [2] Tiwary, P. & Parrinello, M. (2013). From metadynamics to dynamics.
 *     Phys. Rev. Lett. 111, 230602.
 *     DOI: 10.1103/PhysRevLett.111.230602
 */

#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace CV {

/**
 * @class GridAccelerator
 * @brief Fast grid-based bias evaluation for metadynamics
 *
 * Uses regular grid with multilinear interpolation to achieve O(1) bias
 * evaluation instead of O(N_gaussians) direct summation.
 *
 * USAGE:
 * ------
 * ```cpp
 * GridAccelerator grid(num_cvs);
 * grid.setGridBounds(0, -3.0, 3.0, 100);  // CV0: -3 to 3 Å, 100 bins
 * grid.setGridBounds(1, 0.0, 180.0, 90);   // CV1: 0 to 180°, 90 bins
 * grid.initialize();
 *
 * // After depositing Gaussians, update grid
 * grid.addGaussianToGrid(centers, heights, sigmas);
 *
 * // Fast bias evaluation
 * double bias = grid.interpolateBias(cv_values);
 * Eigen::VectorXd grad = grid.interpolateGradient(cv_values);
 * ```
 */
class GridAccelerator {
public:
    /**
     * @brief Constructor
     * @param num_cvs Number of collective variables (1-3 supported)
     */
    GridAccelerator(int num_cvs)
        : m_num_cvs(num_cvs)
        , m_initialized(false)
        , m_grid_size(0)
    {
        if (num_cvs < 1 || num_cvs > 3) {
            throw std::runtime_error("GridAccelerator: Only 1-3 CVs supported");
        }

        m_grid_min.resize(num_cvs);
        m_grid_max.resize(num_cvs);
        m_grid_bins.resize(num_cvs);
        m_grid_spacing.resize(num_cvs);
    }

    /**
     * @brief Set grid bounds for a CV dimension
     * @param cv_index CV index (0, 1, or 2)
     * @param min_val Minimum CV value
     * @param max_val Maximum CV value
     * @param num_bins Number of grid bins
     */
    void setGridBounds(int cv_index, double min_val, double max_val, int num_bins) {
        if (cv_index < 0 || cv_index >= m_num_cvs) {
            throw std::out_of_range("GridAccelerator: Invalid CV index");
        }

        m_grid_min[cv_index] = min_val;
        m_grid_max[cv_index] = max_val;
        m_grid_bins[cv_index] = num_bins;
        m_grid_spacing[cv_index] = (max_val - min_val) / (num_bins - 1);
    }

    /**
     * @brief Initialize grid (allocate memory and set to zero)
     *
     * Must be called after setting all grid bounds and before using grid.
     */
    void initialize() {
        // Compute total grid size
        m_grid_size = 1;
        for (int d = 0; d < m_num_cvs; ++d) {
            m_grid_size *= m_grid_bins[d];
        }

        // Allocate grid storage (bias potential)
        m_grid_data.resize(m_grid_size, 0.0);

        // Allocate gradient grids (one per CV dimension)
        m_grid_gradients.resize(m_num_cvs);
        for (int d = 0; d < m_num_cvs; ++d) {
            m_grid_gradients[d].resize(m_grid_size, 0.0);
        }

        m_initialized = true;

        // Memory usage warning for 3D
        if (m_num_cvs == 3 && m_grid_size > 1000000) {
            double memory_mb = (m_grid_size * (1 + m_num_cvs) * sizeof(double)) / (1024.0 * 1024.0);
            // Note: Would use CurcumaLogger here if included
            // For now, just set a flag or use std::cerr
        }
    }

    /**
     * @brief Add a Gaussian to the grid
     *
     * Incrementally updates grid with new Gaussian contribution.
     *
     * @param centers CV values at Gaussian center (size = num_cvs)
     * @param height Gaussian height
     * @param sigmas Gaussian widths (size = num_cvs)
     *
     * OPTIMIZATION: Only updates grid points within 6σ of center (>99.9% of Gaussian)
     */
    void addGaussianToGrid(const std::vector<double>& centers,
                          double height,
                          const std::vector<double>& sigmas)
    {
        if (!m_initialized) {
            throw std::runtime_error("GridAccelerator: Must call initialize() before addGaussianToGrid()");
        }

        if (centers.size() != static_cast<size_t>(m_num_cvs) ||
            sigmas.size() != static_cast<size_t>(m_num_cvs)) {
            throw std::invalid_argument("GridAccelerator: centers/sigmas size mismatch with num_cvs");
        }

        // Determine grid region to update (6σ cutoff)
        std::vector<int> min_idx(m_num_cvs), max_idx(m_num_cvs);
        for (int d = 0; d < m_num_cvs; ++d) {
            double cutoff = 6.0 * sigmas[d];
            min_idx[d] = std::max(0, gridIndex(d, centers[d] - cutoff));
            max_idx[d] = std::min(m_grid_bins[d] - 1, gridIndex(d, centers[d] + cutoff));
        }

        // Loop over relevant grid points
        if (m_num_cvs == 1) {
            for (int i = min_idx[0]; i <= max_idx[0]; ++i) {
                double s = gridValue(0, i);
                double ds = s - centers[0];
                double arg = ds * ds / (2.0 * sigmas[0] * sigmas[0]);

                if (arg < 20.0) {  // Avoid exp underflow
                    double gaussian = height * std::exp(-arg);
                    int flat_idx = i;
                    m_grid_data[flat_idx] += gaussian;

                    // Gradient: ∂V/∂s = -V * (s - s0) / σ²
                    m_grid_gradients[0][flat_idx] += -gaussian * ds / (sigmas[0] * sigmas[0]);
                }
            }
        } else if (m_num_cvs == 2) {
            for (int i = min_idx[0]; i <= max_idx[0]; ++i) {
                double s0 = gridValue(0, i);
                double ds0 = s0 - centers[0];

                for (int j = min_idx[1]; j <= max_idx[1]; ++j) {
                    double s1 = gridValue(1, j);
                    double ds1 = s1 - centers[1];

                    double arg = (ds0 * ds0) / (2.0 * sigmas[0] * sigmas[0]) +
                                (ds1 * ds1) / (2.0 * sigmas[1] * sigmas[1]);

                    if (arg < 20.0) {
                        double gaussian = height * std::exp(-arg);
                        int flat_idx = flatIndex2D(i, j);
                        m_grid_data[flat_idx] += gaussian;

                        m_grid_gradients[0][flat_idx] += -gaussian * ds0 / (sigmas[0] * sigmas[0]);
                        m_grid_gradients[1][flat_idx] += -gaussian * ds1 / (sigmas[1] * sigmas[1]);
                    }
                }
            }
        } else if (m_num_cvs == 3) {
            for (int i = min_idx[0]; i <= max_idx[0]; ++i) {
                double s0 = gridValue(0, i);
                double ds0 = s0 - centers[0];

                for (int j = min_idx[1]; j <= max_idx[1]; ++j) {
                    double s1 = gridValue(1, j);
                    double ds1 = s1 - centers[1];

                    for (int k = min_idx[2]; k <= max_idx[2]; ++k) {
                        double s2 = gridValue(2, k);
                        double ds2 = s2 - centers[2];

                        double arg = (ds0 * ds0) / (2.0 * sigmas[0] * sigmas[0]) +
                                    (ds1 * ds1) / (2.0 * sigmas[1] * sigmas[1]) +
                                    (ds2 * ds2) / (2.0 * sigmas[2] * sigmas[2]);

                        if (arg < 20.0) {
                            double gaussian = height * std::exp(-arg);
                            int flat_idx = flatIndex3D(i, j, k);
                            m_grid_data[flat_idx] += gaussian;

                            m_grid_gradients[0][flat_idx] += -gaussian * ds0 / (sigmas[0] * sigmas[0]);
                            m_grid_gradients[1][flat_idx] += -gaussian * ds1 / (sigmas[1] * sigmas[1]);
                            m_grid_gradients[2][flat_idx] += -gaussian * ds2 / (sigmas[2] * sigmas[2]);
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief Interpolate bias potential at CV point
     *
     * Uses multilinear interpolation (bilinear for 2D, trilinear for 3D).
     *
     * @param cv_values Current CV values (size = num_cvs)
     * @return Interpolated bias potential
     */
    double interpolateBias(const std::vector<double>& cv_values) const {
        if (!m_initialized) {
            throw std::runtime_error("GridAccelerator: Must call initialize() before interpolateBias()");
        }

        if (cv_values.size() != static_cast<size_t>(m_num_cvs)) {
            throw std::invalid_argument("GridAccelerator: cv_values size mismatch with num_cvs");
        }

        // Check bounds (extrapolate to zero outside grid)
        for (int d = 0; d < m_num_cvs; ++d) {
            if (cv_values[d] < m_grid_min[d] || cv_values[d] > m_grid_max[d]) {
                return 0.0;  // Outside grid bounds
            }
        }

        if (m_num_cvs == 1) {
            return interpolate1D(cv_values[0]);
        } else if (m_num_cvs == 2) {
            return interpolate2D(cv_values[0], cv_values[1]);
        } else {  // m_num_cvs == 3
            return interpolate3D(cv_values[0], cv_values[1], cv_values[2]);
        }
    }

    /**
     * @brief Interpolate bias gradient at CV point
     *
     * @param cv_values Current CV values (size = num_cvs)
     * @return Interpolated gradient vector (size = num_cvs)
     */
    Eigen::VectorXd interpolateGradient(const std::vector<double>& cv_values) const {
        if (!m_initialized) {
            throw std::runtime_error("GridAccelerator: Must call initialize() before interpolateGradient()");
        }

        Eigen::VectorXd gradient = Eigen::VectorXd::Zero(m_num_cvs);

        // Check bounds
        for (int d = 0; d < m_num_cvs; ++d) {
            if (cv_values[d] < m_grid_min[d] || cv_values[d] > m_grid_max[d]) {
                return gradient;  // Zero gradient outside grid
            }
        }

        // Interpolate each gradient component
        for (int d = 0; d < m_num_cvs; ++d) {
            if (m_num_cvs == 1) {
                gradient[d] = interpolate1D_fromGrid(m_grid_gradients[d], cv_values[0]);
            } else if (m_num_cvs == 2) {
                gradient[d] = interpolate2D_fromGrid(m_grid_gradients[d], cv_values[0], cv_values[1]);
            } else {
                gradient[d] = interpolate3D_fromGrid(m_grid_gradients[d], cv_values[0], cv_values[1], cv_values[2]);
            }
        }

        return gradient;
    }

    /**
     * @brief Get grid memory usage in bytes
     */
    size_t getMemoryUsage() const {
        return m_grid_size * (1 + m_num_cvs) * sizeof(double);
    }

    /**
     * @brief Check if grid is initialized
     */
    bool isInitialized() const { return m_initialized; }

    /**
     * @brief Get grid size (total number of points)
     */
    int getGridSize() const { return m_grid_size; }

private:
    int m_num_cvs;                           ///< Number of CVs
    bool m_initialized;                      ///< Initialization flag
    int m_grid_size;                         ///< Total grid points
    std::vector<double> m_grid_min;          ///< Grid minimum bounds
    std::vector<double> m_grid_max;          ///< Grid maximum bounds
    std::vector<int> m_grid_bins;            ///< Number of bins per dimension
    std::vector<double> m_grid_spacing;      ///< Grid spacing per dimension
    std::vector<double> m_grid_data;         ///< Grid values (bias potential)
    std::vector<std::vector<double>> m_grid_gradients;  ///< Grid gradients (one per CV)

    /**
     * @brief Convert CV value to grid index
     */
    int gridIndex(int dim, double value) const {
        int idx = static_cast<int>((value - m_grid_min[dim]) / m_grid_spacing[dim]);
        return std::max(0, std::min(m_grid_bins[dim] - 1, idx));
    }

    /**
     * @brief Convert grid index to CV value
     */
    double gridValue(int dim, int idx) const {
        return m_grid_min[dim] + idx * m_grid_spacing[dim];
    }

    /**
     * @brief Flatten 2D index to 1D
     */
    int flatIndex2D(int i, int j) const {
        return i * m_grid_bins[1] + j;
    }

    /**
     * @brief Flatten 3D index to 1D
     */
    int flatIndex3D(int i, int j, int k) const {
        return (i * m_grid_bins[1] + j) * m_grid_bins[2] + k;
    }

    /**
     * @brief 1D linear interpolation
     */
    double interpolate1D(double s) const {
        double frac = (s - m_grid_min[0]) / m_grid_spacing[0];
        int i0 = static_cast<int>(frac);
        int i1 = i0 + 1;

        if (i1 >= m_grid_bins[0]) {
            i1 = m_grid_bins[0] - 1;
            i0 = i1 - 1;
        }

        double t = frac - i0;
        return (1.0 - t) * m_grid_data[i0] + t * m_grid_data[i1];
    }

    /**
     * @brief 1D interpolation from specific grid
     */
    double interpolate1D_fromGrid(const std::vector<double>& grid, double s) const {
        double frac = (s - m_grid_min[0]) / m_grid_spacing[0];
        int i0 = static_cast<int>(frac);
        int i1 = i0 + 1;

        if (i1 >= m_grid_bins[0]) {
            i1 = m_grid_bins[0] - 1;
            i0 = i1 - 1;
        }

        double t = frac - i0;
        return (1.0 - t) * grid[i0] + t * grid[i1];
    }

    /**
     * @brief 2D bilinear interpolation
     */
    double interpolate2D(double s0, double s1) const {
        // Find grid cell
        double frac0 = (s0 - m_grid_min[0]) / m_grid_spacing[0];
        double frac1 = (s1 - m_grid_min[1]) / m_grid_spacing[1];

        int i0 = static_cast<int>(frac0);
        int j0 = static_cast<int>(frac1);
        int i1 = std::min(i0 + 1, m_grid_bins[0] - 1);
        int j1 = std::min(j0 + 1, m_grid_bins[1] - 1);

        double t = frac0 - i0;
        double u = frac1 - j0;

        // Bilinear interpolation
        double v00 = m_grid_data[flatIndex2D(i0, j0)];
        double v10 = m_grid_data[flatIndex2D(i1, j0)];
        double v01 = m_grid_data[flatIndex2D(i0, j1)];
        double v11 = m_grid_data[flatIndex2D(i1, j1)];

        return (1.0 - t) * (1.0 - u) * v00 +
               t * (1.0 - u) * v10 +
               (1.0 - t) * u * v01 +
               t * u * v11;
    }

    /**
     * @brief 2D bilinear interpolation from specific grid
     */
    double interpolate2D_fromGrid(const std::vector<double>& grid, double s0, double s1) const {
        double frac0 = (s0 - m_grid_min[0]) / m_grid_spacing[0];
        double frac1 = (s1 - m_grid_min[1]) / m_grid_spacing[1];

        int i0 = static_cast<int>(frac0);
        int j0 = static_cast<int>(frac1);
        int i1 = std::min(i0 + 1, m_grid_bins[0] - 1);
        int j1 = std::min(j0 + 1, m_grid_bins[1] - 1);

        double t = frac0 - i0;
        double u = frac1 - j0;

        double v00 = grid[flatIndex2D(i0, j0)];
        double v10 = grid[flatIndex2D(i1, j0)];
        double v01 = grid[flatIndex2D(i0, j1)];
        double v11 = grid[flatIndex2D(i1, j1)];

        return (1.0 - t) * (1.0 - u) * v00 +
               t * (1.0 - u) * v10 +
               (1.0 - t) * u * v01 +
               t * u * v11;
    }

    /**
     * @brief 3D trilinear interpolation
     */
    double interpolate3D(double s0, double s1, double s2) const {
        double frac0 = (s0 - m_grid_min[0]) / m_grid_spacing[0];
        double frac1 = (s1 - m_grid_min[1]) / m_grid_spacing[1];
        double frac2 = (s2 - m_grid_min[2]) / m_grid_spacing[2];

        int i0 = static_cast<int>(frac0);
        int j0 = static_cast<int>(frac1);
        int k0 = static_cast<int>(frac2);
        int i1 = std::min(i0 + 1, m_grid_bins[0] - 1);
        int j1 = std::min(j0 + 1, m_grid_bins[1] - 1);
        int k1 = std::min(k0 + 1, m_grid_bins[2] - 1);

        double t = frac0 - i0;
        double u = frac1 - j0;
        double v = frac2 - k0;

        // Trilinear interpolation (8 corners of cube)
        double v000 = m_grid_data[flatIndex3D(i0, j0, k0)];
        double v100 = m_grid_data[flatIndex3D(i1, j0, k0)];
        double v010 = m_grid_data[flatIndex3D(i0, j1, k0)];
        double v110 = m_grid_data[flatIndex3D(i1, j1, k0)];
        double v001 = m_grid_data[flatIndex3D(i0, j0, k1)];
        double v101 = m_grid_data[flatIndex3D(i1, j0, k1)];
        double v011 = m_grid_data[flatIndex3D(i0, j1, k1)];
        double v111 = m_grid_data[flatIndex3D(i1, j1, k1)];

        return (1.0 - t) * (1.0 - u) * (1.0 - v) * v000 +
               t * (1.0 - u) * (1.0 - v) * v100 +
               (1.0 - t) * u * (1.0 - v) * v010 +
               t * u * (1.0 - v) * v110 +
               (1.0 - t) * (1.0 - u) * v * v001 +
               t * (1.0 - u) * v * v101 +
               (1.0 - t) * u * v * v011 +
               t * u * v * v111;
    }

    /**
     * @brief 3D trilinear interpolation from specific grid
     */
    double interpolate3D_fromGrid(const std::vector<double>& grid, double s0, double s1, double s2) const {
        double frac0 = (s0 - m_grid_min[0]) / m_grid_spacing[0];
        double frac1 = (s1 - m_grid_min[1]) / m_grid_spacing[1];
        double frac2 = (s2 - m_grid_min[2]) / m_grid_spacing[2];

        int i0 = static_cast<int>(frac0);
        int j0 = static_cast<int>(frac1);
        int k0 = static_cast<int>(frac2);
        int i1 = std::min(i0 + 1, m_grid_bins[0] - 1);
        int j1 = std::min(j0 + 1, m_grid_bins[1] - 1);
        int k1 = std::min(k0 + 1, m_grid_bins[2] - 1);

        double t = frac0 - i0;
        double u = frac1 - j0;
        double v = frac2 - k0;

        double v000 = grid[flatIndex3D(i0, j0, k0)];
        double v100 = grid[flatIndex3D(i1, j0, k0)];
        double v010 = grid[flatIndex3D(i0, j1, k0)];
        double v110 = grid[flatIndex3D(i1, j1, k0)];
        double v001 = grid[flatIndex3D(i0, j0, k1)];
        double v101 = grid[flatIndex3D(i1, j0, k1)];
        double v011 = grid[flatIndex3D(i0, j1, k1)];
        double v111 = grid[flatIndex3D(i1, j1, k1)];

        return (1.0 - t) * (1.0 - u) * (1.0 - v) * v000 +
               t * (1.0 - u) * (1.0 - v) * v100 +
               (1.0 - t) * u * (1.0 - v) * v010 +
               t * u * (1.0 - v) * v110 +
               (1.0 - t) * (1.0 - u) * v * v001 +
               t * (1.0 - u) * v * v101 +
               (1.0 - t) * u * v * v011 +
               t * u * v * v111;
    }
};

}  // namespace CV
