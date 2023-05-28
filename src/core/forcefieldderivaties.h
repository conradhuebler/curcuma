/*
 * < General derivates for force field components. >
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/global.h"
#include "src/tools/general.h"

#include <Eigen/Dense>

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
    auto nij = rij / rij.norm();
    Eigen::Vector3d rkj = k - j;
    auto nkj = rkj / rkj.norm();
    double costheta = (rij.dot(rkj) / (sqrt(rij.dot(rij) * rkj.dot(rkj))));

    if (!gradient)
        return costheta;
    derivate = Matrix::Zero(3, 3);

    double sintheta = sin(acos(costheta));
    double dThetadCosTheta = 1 / sintheta;
    derivate = Matrix::Zero(3, 3);
    derivate.row(0) = -dThetadCosTheta * (nkj - nij * costheta) / (rij.norm());
    derivate.row(2) = -dThetadCosTheta * (nij - nkj * costheta) / (rkj.norm());
    derivate.row(1) = -derivate.row(0) - derivate.row(2);

    return costheta;
}
