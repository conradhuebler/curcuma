/*
 * <D4 pair C6 and dC6/dCN from the per-atom half-contraction>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#pragma once

/**
 * @file d4_pair_c6.h
 * @brief The two per-pair D4 formulas of GFN-FF, shared by the D4 generator and the
 *        FFWorkspace energy loop (Claude Generated, Oct 2026).
 *
 * D4 C6 with CN-only Gaussian weights w (gfnff_gdisp0.f90:405 weight_cn):
 *   C6(i,j)       = sum_a sum_b w_i[a] w_j[b] C6ref(Z_i,a; Z_j,b)
 * With the per-atom half-contraction g_i[Z][b] = sum_a w_i[a] C6ref(Z_i,a; Z,b) and
 * dg_i[Z][b] = sum_a dw_i[a]/dCN_i C6ref(Z_i,a; Z,b) (D4ParameterGenerator::
 * precomputeC6HalfContraction) every pair needs only a sum over the <= 7 references of j:
 *   C6(i,j)       = sum_b g_i[Z_j][b]  w_j[b]
 *   dC6/dCN_i     = sum_b dg_i[Z_j][b] w_j[b]
 *   dC6/dCN_j     = sum_b g_i[Z_j][b]  dw_j[b]
 * Reference for the derivative use: gfnff_gdisp0.f90:262-305 (dc6dcn), 382-395 (gradient).
 *
 * Both the generator (getChargeWeightedC6, computeDC6DCN) and the fused energy loop
 * (FFWorkspace::calcD4Dispersion) call these functions, so the operation order - and with
 * it every bit of the result - is the same in both places.
 */
namespace D4PairC6 {

/// C6(i,j) = sum_{b < nb} g_i[b] * w_j[b]   (g_i: row of atom i for the element of j)
inline double c6(const double* g_i, const double* w_j, int nb)
{
    double c6 = 0.0;
    for (int b = 0; b < nb; ++b) c6 += g_i[b] * w_j[b];
    return c6;
}

/// dC6(i,j)/dCN_i and dC6(i,j)/dCN_j, same loop as D4ParameterGenerator::computeDC6DCN
inline void dc6dcn(const double* g_i, const double* dg_i, const double* w_j, const double* dw_j,
                   int nb, double& dc6_di, double& dc6_dj)
{
    dc6_di = 0.0;
    dc6_dj = 0.0;
    for (int b = 0; b < nb; ++b) {
        dc6_di += dg_i[b] * w_j[b];
        dc6_dj += g_i[b] * dw_j[b];
    }
}

} // namespace D4PairC6
