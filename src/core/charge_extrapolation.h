/*
 * <Multi-step charge/density extrapolation weights — shared ASPC / least-squares>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <string>

// Multi-step extrapolation weights w_j for predicting the next value of a slowly
// varying vector from its history: x_pred = Σ_j w_j * x_hist[j], with x_hist[0] the
// most recent. The weights are charge/trace-conserving (Σ w_j = 1): if every history
// entry shares the same sum (e.g. total charge), so does the prediction. Shared by the
// native GFN SCF (SCC vector) and the GFN-FF EEQ PCG warm-start (z1/Z2 solution).
// Claude Generated.
//
//  mode = "aspc"  : Always Stable Predictor-Corrector (Kolafa, J. Comput. Chem. 25,
//                   335 (2004)). m = min(order+2, n_hist) points, fixed binomial
//                   coefficients w_j = (-1)^j C(m, j+1). Time-reversible for a fixed
//                   step (MD); order 0 -> the 2-point linear predictor [2, -1].
//  mode = "gauss" : least-squares polynomial of degree d = min(order, m-1) fit to the
//                   history at abscissae x_j = -(j+1) and extrapolated to x = 0
//                   (w = V·a with (VᵀV)a = e_0). Robust for the irregular steps of a
//                   geometry optimisation.
//
// Returns an empty vector when there is too little history (caller keeps its 1-step
// warm-start) or on a numerical failure.
namespace curcuma::extrapolation {

inline Eigen::VectorXd weights(const std::string& mode, int order, int n_hist)
{
    if (n_hist < 2)
        return {};

    if (mode == "aspc") {
        const int m = std::min(order + 2, n_hist);
        if (m < 2) return {};
        // Binomial C(m, i) via the multiplicative recurrence (exact for small m).
        Eigen::VectorXd w(m);
        long long c = 1; // C(m,0)
        for (int j = 0; j < m; ++j) {
            c = c * (m - j) / (j + 1); // advance C(m,j) -> C(m,j+1)
            w(j) = (j % 2 == 0 ? 1.0 : -1.0) * static_cast<double>(c);
        }
        return w;
    }

    if (mode == "gauss") {
        const int m = n_hist;
        const int d = std::min(order, m - 1); // polynomial degree actually fit
        if (d < 1) return {};
        Eigen::MatrixXd V(m, d + 1);
        for (int j = 0; j < m; ++j) {
            const double x = -static_cast<double>(j + 1);
            double xp = 1.0;
            for (int p = 0; p <= d; ++p) { V(j, p) = xp; xp *= x; }
        }
        Eigen::VectorXd e0 = Eigen::VectorXd::Zero(d + 1);
        e0(0) = 1.0;
        Eigen::MatrixXd VtV = V.transpose() * V;
        Eigen::VectorXd a = VtV.ldlt().solve(e0);
        Eigen::VectorXd w = V * a;
        if (!w.allFinite())
            return {};
        return w;
    }

    return {}; // unknown mode -> caller falls back
}

} // namespace curcuma::extrapolation
