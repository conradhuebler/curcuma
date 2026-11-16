/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2022 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#pragma once

#include "qmdff_par.h"
#include "src/core/global.h"
#include "uff_par.h"

#include <cmath>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"

namespace UFF {
inline double BondStretching(const Vector& i, const Vector& j, Matrix& derivate, bool gradient)
{
    Vector ij = i - j;
    double distance = (ij).norm();
    if (!gradient)
        return distance;
    derivate = Matrix::Zero(2, 3);
    derivate.row(0) = ij / distance;
    derivate.row(1) = -ij / distance;
    return distance;
}

inline double AngleBending(const Vector& i, const Vector& j, const Vector& k, Matrix& derivate, bool gradient)
{
    Eigen::Vector3d rij = i - j;
    Eigen::Vector3d rkj = k - j;
    double rij_norm = rij.norm();
    double rkj_norm = rkj.norm();

    // Claude Generated (Jul 2026, F1): guard the singularities that produced NaN
    // gradients (1/sin(theta), 1/|rij|, 1/|rkj|, acos out-of-range). Mirrors the
    // guards already present in calculateDihedralAngle (gfnff_geometry.h:89-96,140-143).
    // Without this, a near-collinear angle (current geometry collinear, theta0 not)
    // underflowed sin(theta) to 0 -> 1/0 = Inf -> 0*Inf = NaN in derivate, which
    // propagated to every downstream per-term gradient and forced the energy to 0.0.
    if (rij_norm < 1e-10 || rkj_norm < 1e-10) {
        if (gradient) derivate = Matrix::Zero(3, 3);
        return 1.0;  // degenerate: treat as 0 angle (cos=1)
    }
    Eigen::Vector3d nij = rij / rij_norm;
    Eigen::Vector3d nkj = rkj / rkj_norm;
    double costheta = rij.dot(rkj) / (rij_norm * rkj_norm);
    // clamp before acos (round-off can push |costheta| past 1)
    costheta = std::max(-1.0, std::min(1.0, costheta));

    if (!gradient)
        return costheta;
    derivate = Matrix::Zero(3, 3);

    double sintheta = std::sqrt(std::max(0.0, 1.0 - costheta * costheta));
    if (sintheta < 1e-14) {
        // collinear: the 1/sin(theta) factor is singular. The caller multiplies
        // derivate by sin(theta) (dedtheta = 2k sin(theta)(cos0-cos)), so the
        // product is 0 * Inf = NaN. Return a zero derivate; the angle-term energy
        // is still finite via the clamped costheta, only its gradient is undefined.
        return costheta;
    }
    double dThetadCosTheta = 1.0 / sintheta;
    derivate.row(0) = -dThetadCosTheta * (nkj - nij * costheta) / rij_norm;
    derivate.row(2) = -dThetadCosTheta * (nij - nkj * costheta) / rkj_norm;
    derivate.row(1) = -derivate.row(0) - derivate.row(2);

    return costheta;
}

inline Eigen::Vector3d NormalVector(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k)
{
    return (j - i).cross(j - k);
}

}
