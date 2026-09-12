/*
 * < rev-gfnff: continuous bond order and over-coordination helpers >
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated (Sep 2026, rev-gfnff stage 1)
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

#include <cmath>

/**
 * @file rev_bond_order.h
 * @brief The two smooth functions of the reactive GFN-FF stage 1.
 *
 * Continuous bond order (the "weight" every bonded term is multiplied with):
 *
 *     b(r) = 1/2 (1 + erf(k (r - R) / R)),   k < 0
 *
 * the erf counting function of the GFN-FF coordination number, with the switching
 * radius R = f_b (rcov_i + rcov_j) fat_i fat_j - the same threshold expression the
 * react topology scan uses - and a steepness k. For k = -7.5 and f_b = 2.0 the weight is
 * 1 to 1e-6 at the bond length, 1/2 at r = R and 0.05 at about 1.15 R, so the equilibrium
 * region of every term is untouched and a term switches off smoothly where the bond list
 * would otherwise flip it off in one step (docs/REV_GFNFF_ROADMAP.md, WP3).
 *
 * Over-coordination energy (replaces the hard valence cap of the react scan):
 *
 *     E_over,i = p_Z  sp(sum_j b_ij BO_ij - Val_Z)^2,   sp(x) = ln(1 + e^{k x}) / k
 *
 * a smooth penalty that rises once the bond-order sum of atom i exceeds its nominal
 * valence; sp is C-infinity and exponentially small for negative arguments.
 */
namespace RevGFNFF {

/// Continuous bond order and its radial derivative. R and r in the same unit.
inline double bondOrder(double r, double R, double k, double* dbdr = nullptr)
{
    const double x = k * (r - R) / R;
    if (dbdr)
        *dbdr = k / (R * 1.7724538509055159) * std::exp(-x * x); // k/(R sqrt(pi)) e^{-x^2}
    return 0.5 * (1.0 + std::erf(x));
}

/// softplus with steepness k, and its derivative (the logistic function)
inline double softplus(double x, double k, double* dsp = nullptr)
{
    const double kx = k * x;
    double sp;
    if (kx > 30.0) {
        sp = x;
        if (dsp) *dsp = 1.0;
    } else if (kx < -30.0) {
        sp = std::exp(kx) / k;
        if (dsp) *dsp = std::exp(kx);
    } else {
        sp = std::log1p(std::exp(kx)) / k;
        if (dsp) *dsp = 1.0 / (1.0 + std::exp(-kx));
    }
    return sp;
}

} // namespace RevGFNFF
