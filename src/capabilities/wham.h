/*
 * <Self-consistent WHAM solver for 1D umbrella sampling with harmonic windows.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested).
 * Header-only so that NebMD and the standalone test share exactly one implementation.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace curcuma {

/*! \brief WHAM (Kumar et al., J. Comput. Chem. 13, 1011 (1992)) for windows V_i = k_i/2 (s-s0_i)^2.
 *
 *   P(s)        = sum_i H_i(s) / sum_i N_i exp(-(V_i(s) - f_i)/kT)
 *   exp(-f_i/kT) = sum_s P(s) exp(-V_i(s)/kT)
 *
 * \param hist  hist[i][b]: counts of window i in bin b, bins equidistant over s in [0,1]
 * \param n     total samples per window (normally the row sum of hist)
 * \param s0    window centres
 * \param k     stiffness PER WINDOW, in the units of kT / s^2 (the k the dynamics used)
 * \param kT    thermal energy, same energy unit as k
 * \return      unnormalised P(s) per bin; G(s) = -kT ln P(s). f_i returned via f_out if given.
 */
inline std::vector<double> whamSolve(const std::vector<std::vector<int>>& hist,
    const std::vector<double>& n, const std::vector<double>& s0,
    const std::vector<double>& k, double kT, std::vector<double>* f_out = nullptr)
{
    const int W = static_cast<int>(hist.size());
    const int B = W > 0 ? static_cast<int>(hist[0].size()) : 0;
    std::vector<double> f(W, 0.0), P(B, 0.0);
    auto bias = [&](int i, int b) {
        const double d = (b + 0.5) / B - s0[i];
        return 0.5 * k[i] * d * d;
    };
    for (int iter = 0; iter < 2000; ++iter) {
        for (int b = 0; b < B; ++b) {
            double num = 0.0, den = 0.0;
            for (int i = 0; i < W; ++i) {
                num += hist[i][b];
                den += n[i] * std::exp(-(bias(i, b) - f[i]) / kT);
            }
            P[b] = (den > 1e-300) ? num / den : 0.0;
        }
        double maxdev = 0.0;
        for (int i = 0; i < W; ++i) {
            double sum = 0.0;
            for (int b = 0; b < B; ++b)
                sum += P[b] * std::exp(-bias(i, b) / kT);
            const double fnew = (sum > 1e-300) ? -kT * std::log(sum) : f[i];
            maxdev = std::max(maxdev, std::abs(fnew - f[i]));
            f[i] = fnew;
        }
        if (maxdev < 1e-10)
            break;
    }
    if (f_out)
        *f_out = f;
    return P;
}

} // namespace curcuma
