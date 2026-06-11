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

// ============================================================================
// tblite-faithful Lebedev-Laikov grid API (Claude Generated, June 2026)
//
// Ported one-to-one from external/tblite/src/tblite/mesh/lebedev.f90 for the
// native CPCM (ddCOSMO) solvation model, which builds a per-atom angular grid.
// Kept separate from the LebedevGrid class above (used by ALPB/GBSA SASA) so
// those callers are unaffected. Weights are the raw Lebedev weights (sum = 1),
// exactly as tblite hands them to new_domain_decomposition.
// ============================================================================

/// Available Lebedev grid point counts (tblite grid_size, 32 entries).
inline constexpr int kLebedevGridSize[32] = {
    6, 14, 26, 38, 50, 74, 86, 110,
    146, 170, 194, 230, 266, 302, 350, 434,
    590, 770, 974, 1202, 1454, 1730, 2030, 2354,
    2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
};

namespace lebedev_detail {

// Octahedral point generators (tblite gen_oh1..gen_oh6, exact ordering).
inline void genOh1(std::vector<LebedevPoint>& g, double v)
{
    const double a = 1.0;
    g.push_back({ a, 0, 0, v}); g.push_back({-a, 0, 0, v});
    g.push_back({ 0, a, 0, v}); g.push_back({ 0, -a, 0, v});
    g.push_back({ 0, 0, a, v}); g.push_back({ 0, 0, -a, v});
}
inline void genOh2(std::vector<LebedevPoint>& g, double v)
{
    const double a = std::sqrt(0.5);
    g.push_back({0, a, a, v}); g.push_back({0, -a, a, v}); g.push_back({0, a, -a, v}); g.push_back({0, -a, -a, v});
    g.push_back({a, 0, a, v}); g.push_back({-a, 0, a, v}); g.push_back({a, 0, -a, v}); g.push_back({-a, 0, -a, v});
    g.push_back({a, a, 0, v}); g.push_back({-a, a, 0, v}); g.push_back({a, -a, 0, v}); g.push_back({-a, -a, 0, v});
}
inline void genOh3(std::vector<LebedevPoint>& g, double v)
{
    const double a = std::sqrt(1.0 / 3.0);
    g.push_back({a, a, a, v}); g.push_back({-a, a, a, v}); g.push_back({a, -a, a, v}); g.push_back({-a, -a, a, v});
    g.push_back({a, a, -a, v}); g.push_back({-a, a, -a, v}); g.push_back({a, -a, -a, v}); g.push_back({-a, -a, -a, v});
}
inline void genOh4(std::vector<LebedevPoint>& g, double a, double v)
{
    const double b = std::sqrt(1.0 - 2.0 * a * a);
    g.push_back({a, a, b, v}); g.push_back({-a, a, b, v}); g.push_back({a, -a, b, v}); g.push_back({-a, -a, b, v});
    g.push_back({a, a, -b, v}); g.push_back({-a, a, -b, v}); g.push_back({a, -a, -b, v}); g.push_back({-a, -a, -b, v});
    g.push_back({a, b, a, v}); g.push_back({-a, b, a, v}); g.push_back({a, -b, a, v}); g.push_back({-a, -b, a, v});
    g.push_back({a, b, -a, v}); g.push_back({-a, b, -a, v}); g.push_back({a, -b, -a, v}); g.push_back({-a, -b, -a, v});
    g.push_back({b, a, a, v}); g.push_back({-b, a, a, v}); g.push_back({b, -a, a, v}); g.push_back({-b, -a, a, v});
    g.push_back({b, a, -a, v}); g.push_back({-b, a, -a, v}); g.push_back({b, -a, -a, v}); g.push_back({-b, -a, -a, v});
}
inline void genOh5(std::vector<LebedevPoint>& g, double a, double v)
{
    const double b = std::sqrt(1.0 - a * a);
    g.push_back({a, b, 0, v}); g.push_back({-a, b, 0, v}); g.push_back({a, -b, 0, v}); g.push_back({-a, -b, 0, v});
    g.push_back({b, a, 0, v}); g.push_back({-b, a, 0, v}); g.push_back({b, -a, 0, v}); g.push_back({-b, -a, 0, v});
    g.push_back({a, 0, b, v}); g.push_back({-a, 0, b, v}); g.push_back({a, 0, -b, v}); g.push_back({-a, 0, -b, v});
    g.push_back({b, 0, a, v}); g.push_back({-b, 0, a, v}); g.push_back({b, 0, -a, v}); g.push_back({-b, 0, -a, v});
    g.push_back({0, a, b, v}); g.push_back({0, -a, b, v}); g.push_back({0, a, -b, v}); g.push_back({0, -a, -b, v});
    g.push_back({0, b, a, v}); g.push_back({0, -b, a, v}); g.push_back({0, b, -a, v}); g.push_back({0, -b, -a, v});
}

} // namespace lebedev_detail

/**
 * @brief Locate the grid index pos with grid_size[pos] <= val < grid_size[pos+1].
 *        1-based to match tblite list_bisection. Clamps to [1, 32].
 */
inline int lebedevListBisection(int val)
{
    const int n = 32;
    if (val <= kLebedevGridSize[0]) return 1;
    if (val >= kLebedevGridSize[n - 1]) return n;
    int lower = 0, current = n + 1;
    while ((current - lower) > 1) {
        const int upper = (current + lower) / 2;       // 1-based midpoint
        if (val >= kLebedevGridSize[upper - 1])
            lower = upper;
        else
            current = upper;
    }
    return lower;
}

/**
 * @brief Generate the Lebedev grid for the 1-based size index nang_index.
 *        Implements indices 1..8 (6..110 points), which covers the CPCM default
 *        (grid_size(6) = 74). Returns an empty grid for unsupported indices.
 * @param nang_index 1-based index into kLebedevGridSize.
 * @param ok         [out, optional] set to false if the index is unsupported.
 */
inline std::vector<LebedevPoint> lebedevAngularGrid(int nang_index, bool* ok = nullptr)
{
    using namespace lebedev_detail;
    std::vector<LebedevPoint> g;
    if (ok) *ok = true;
    switch (nang_index) {
    case 1: // ld0006
        genOh1(g, 0.1666666666666667e+0);
        break;
    case 2: // ld0014
        genOh1(g, 0.6666666666666667e-1);
        genOh3(g, 0.7500000000000000e-1);
        break;
    case 3: // ld0026
        genOh1(g, 0.4761904761904762e-1);
        genOh2(g, 0.3809523809523810e-1);
        genOh3(g, 0.3214285714285714e-1);
        break;
    case 4: // ld0038
        genOh1(g, 0.9523809523809524e-2);
        genOh3(g, 0.3214285714285714e-1);
        genOh5(g, 0.4597008433809831e+0, 0.2857142857142857e-1);
        break;
    case 5: // ld0050
        genOh1(g, 0.1269841269841270e-1);
        genOh2(g, 0.2257495590828924e-1);
        genOh3(g, 0.2109375000000000e-1);
        genOh4(g, 0.3015113445777636e+0, 0.2017333553791887e-1);
        break;
    case 6: // ld0074
        genOh1(g, 0.5130671797338464e-3);
        genOh2(g, 0.1660406956574204e-1);
        genOh3(g, -0.2958603896103896e-1);
        genOh4(g, 0.4803844614152614e+0, 0.2657620708215946e-1);
        genOh5(g, 0.3207726489807764e+0, 0.1652217099371571e-1);
        break;
    case 7: // ld0086
        genOh1(g, 0.1154401154401154e-1);
        genOh3(g, 0.1194390908585628e-1);
        genOh4(g, 0.3696028464541502e+0, 0.1111055571060340e-1);
        genOh4(g, 0.6943540066026664e+0, 0.1187650129453714e-1);
        genOh5(g, 0.3742430390903412e+0, 0.1181230374690448e-1);
        break;
    case 8: // ld0110
        genOh1(g, 0.3828270494937162e-2);
        genOh3(g, 0.9793737512487512e-2);
        genOh4(g, 0.1851156353447362e+0, 0.8211737283191111e-2);
        genOh4(g, 0.6904210483822922e+0, 0.9942814891178103e-2);
        genOh4(g, 0.3956894730559419e+0, 0.9595471336070963e-2);
        genOh5(g, 0.4783690288121502e+0, 0.9694996361663028e-2);
        break;
    default:
        if (ok) *ok = false;
        break;
    }
    return g;
}

} // namespace Solvation
} // namespace Curcuma

#endif // LEBEDEV_GRID_H
