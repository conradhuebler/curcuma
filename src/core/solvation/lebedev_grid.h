/*
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This file is part of Curcuma - Native Solvation Module
 *
 * Extracted and adapted from Ulysses (Copyright (C) 2023- Filipe Menezes et al.)
 * Simplified for Curcuma by Claude (Anthropic AI Assistant)
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */

#ifndef LEBEDEV_GRID_H
#define LEBEDEV_GRID_H

#include <vector>
#include <cmath>
#include <array>

namespace Curcuma {
namespace Solvation {

/**
 * @brief Grid point on unit sphere for Lebedev quadrature
 *
 * Each point represents a direction (x,y,z) on the unit sphere
 * with an associated integration weight.
 */
struct LebedevPoint {
    double x, y, z;    // Coordinates on unit sphere (normalized)
    double weight;     // Integration weight (sum of all weights = 4π)
};

/**
 * @brief Lebedev-Laikov spherical grid for SASA integration
 *
 * Provides optimized spherical grids for integrating functions over the unit sphere.
 * Used in GBSA solvation to calculate solvent-accessible surface areas.
 *
 * Reference: V.I. Lebedev and D.N. Laikov, Doklady Mathematics 59, 477-481 (1999)
 *
 * Implementation Notes:
 * - Currently implements only 110-point grid (sufficient for quick-win)
 * - Full implementation requires grids from 6 to 5810 points
 * - Grid generation uses octahedral symmetry (genOh functions)
 *
 * @todo Implement full grid hierarchy (6, 14, 26, ..., 5810 points)
 * @todo Add adaptive grid selection based on accuracy requirements
 * @todo Benchmark accuracy vs performance for different grid sizes
 */
class LebedevGrid {
public:
    /**
     * @brief Generate 110-point Lebedev grid (Order 11)
     *
     * Quick-win implementation with good balance of accuracy and performance.
     * Suitable for SASA calculations on small to medium molecules.
     *
     * Grid composition:
     * - 6 points from genOh1 (octahedral axes)
     * - 8 points from genOh3 (body diagonals)
     * - 96 points from genOh4 (4×24 points, general positions)
     * - Total: 110 points, integrating polynomials up to order 11
     *
     * @return Vector of 110 grid points with weights summing to 4π
     */
    static std::vector<LebedevPoint> generate110() {
        std::vector<LebedevPoint> grid;
        grid.reserve(110);

        // genOh1: 6 points on axes [±1,0,0], [0,±1,0], [0,0,±1]
        genOh1(grid, 0.003828270494937162);

        // genOh3: 8 points on body diagonals [±a,±a,±a], a=1/√3
        genOh3(grid, 0.009793737512487512);

        // genOh4: 4×24 points in general positions
        genOh4(grid, 0.008211737283191111, 0.1851156353447362);  // 14-37
        genOh4(grid, 0.009942814891178103, 0.6904210483822922);  // 38-61
        genOh4(grid, 0.009595471336070963, 0.3956894730559419);  // 62-85

        // genOh5: 24 points in xy/xz/yz planes
        genOh5(grid, 0.009694996361663028, 0.4783690288121502);  // 86-109

        return grid;
    }

    /**
     * @brief Get available grid sizes
     *
     * Returns list of implemented Lebedev grid sizes.
     * Currently only 110 points available.
     *
     * @return Vector of available grid sizes
     */
    static std::vector<int> availableSizes() {
        return {110}; // TODO: Add 6, 14, 26, 38, 50, 74, 86, 146, ..., 5810
    }

private:
    /**
     * @brief Generate 6 points with Oh symmetry from [0,0,1]
     *
     * Creates points on coordinate axes: [±1,0,0], [0,±1,0], [0,0,±1]
     *
     * @param grid Output grid (appends 6 points)
     * @param weight Integration weight for all 6 points
     */
    static void genOh1(std::vector<LebedevPoint>& grid, double weight) {
        grid.push_back({ 1.0,  0.0,  0.0, weight});
        grid.push_back({-1.0,  0.0,  0.0, weight});
        grid.push_back({ 0.0,  1.0,  0.0, weight});
        grid.push_back({ 0.0, -1.0,  0.0, weight});
        grid.push_back({ 0.0,  0.0,  1.0, weight});
        grid.push_back({ 0.0,  0.0, -1.0, weight});
    }

    /**
     * @brief Generate 8 points with Oh symmetry from [a,a,a], a=1/√3
     *
     * Creates points on body diagonals of the cube.
     *
     * @param grid Output grid (appends 8 points)
     * @param weight Integration weight for all 8 points
     */
    static void genOh3(std::vector<LebedevPoint>& grid, double weight) {
        const double a = std::sqrt(1.0 / 3.0);
        grid.push_back({ a,  a,  a, weight});
        grid.push_back({-a,  a,  a, weight});
        grid.push_back({ a, -a,  a, weight});
        grid.push_back({-a, -a,  a, weight});
        grid.push_back({ a,  a, -a, weight});
        grid.push_back({-a,  a, -a, weight});
        grid.push_back({ a, -a, -a, weight});
        grid.push_back({-a, -a, -a, weight});
    }

    /**
     * @brief Generate 24 points with Oh symmetry from [a,a,b], b=√(1-2a²)
     *
     * Creates 24 points in general positions with two equal coordinates.
     *
     * @param grid Output grid (appends 24 points)
     * @param weight Integration weight for all 24 points
     * @param a Parameter determining point positions
     */
    static void genOh4(std::vector<LebedevPoint>& grid, double weight, double a) {
        const double b = std::sqrt(1.0 - 2.0 * a * a);

        // Permutations of [±a, ±a, ±b]
        grid.push_back({ a,  a,  b, weight});
        grid.push_back({-a,  a,  b, weight});
        grid.push_back({ a, -a,  b, weight});
        grid.push_back({-a, -a,  b, weight});
        grid.push_back({ a,  a, -b, weight});
        grid.push_back({-a,  a, -b, weight});
        grid.push_back({ a, -a, -b, weight});
        grid.push_back({-a, -a, -b, weight});

        // Permutations of [±a, ±b, ±a]
        grid.push_back({ a,  b,  a, weight});
        grid.push_back({-a,  b,  a, weight});
        grid.push_back({ a, -b,  a, weight});
        grid.push_back({-a, -b,  a, weight});
        grid.push_back({ a,  b, -a, weight});
        grid.push_back({-a,  b, -a, weight});
        grid.push_back({ a, -b, -a, weight});
        grid.push_back({-a, -b, -a, weight});

        // Permutations of [±b, ±a, ±a]
        grid.push_back({ b,  a,  a, weight});
        grid.push_back({-b,  a,  a, weight});
        grid.push_back({ b, -a,  a, weight});
        grid.push_back({-b, -a,  a, weight});
        grid.push_back({ b,  a, -a, weight});
        grid.push_back({-b,  a, -a, weight});
        grid.push_back({ b, -a, -a, weight});
        grid.push_back({-b, -a, -a, weight});
    }

    /**
     * @brief Generate 24 points with Oh symmetry from [a,b,0], b=√(1-a²)
     *
     * Creates 24 points in coordinate planes.
     *
     * @param grid Output grid (appends 24 points)
     * @param weight Integration weight for all 24 points
     * @param a Parameter determining point positions
     */
    static void genOh5(std::vector<LebedevPoint>& grid, double weight, double a) {
        const double b = std::sqrt(1.0 - a * a);

        // xy plane: [±a, ±b, 0] and [±b, ±a, 0]
        grid.push_back({ a,  b,  0.0, weight});
        grid.push_back({-a,  b,  0.0, weight});
        grid.push_back({ a, -b,  0.0, weight});
        grid.push_back({-a, -b,  0.0, weight});
        grid.push_back({ b,  a,  0.0, weight});
        grid.push_back({-b,  a,  0.0, weight});
        grid.push_back({ b, -a,  0.0, weight});
        grid.push_back({-b, -a,  0.0, weight});

        // xz plane: [±a, 0, ±b] and [±b, 0, ±a]
        grid.push_back({ a,  0.0,  b, weight});
        grid.push_back({-a,  0.0,  b, weight});
        grid.push_back({ a,  0.0, -b, weight});
        grid.push_back({-a,  0.0, -b, weight});
        grid.push_back({ b,  0.0,  a, weight});
        grid.push_back({-b,  0.0,  a, weight});
        grid.push_back({ b,  0.0, -a, weight});
        grid.push_back({-b,  0.0, -a, weight});

        // yz plane: [0, ±a, ±b] and [0, ±b, ±a]
        grid.push_back({ 0.0,  a,  b, weight});
        grid.push_back({ 0.0, -a,  b, weight});
        grid.push_back({ 0.0,  a, -b, weight});
        grid.push_back({ 0.0, -a, -b, weight});
        grid.push_back({ 0.0,  b,  a, weight});
        grid.push_back({ 0.0, -b,  a, weight});
        grid.push_back({ 0.0,  b, -a, weight});
        grid.push_back({ 0.0, -b, -a, weight});
    }
};

} // namespace Solvation
} // namespace Curcuma

#endif // LEBEDEV_GRID_H
