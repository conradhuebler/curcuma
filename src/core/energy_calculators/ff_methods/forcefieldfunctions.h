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

inline Eigen::Vector3d NormalVector(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k)
{
    return (j - i).cross(j - k);
}

inline double Dihedral(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, const Eigen::Vector3d& l, double V, double n, double phi0)
{
    Eigen::Vector3d nijk = NormalVector(i, j, k);
    Eigen::Vector3d njkl = NormalVector(j, k, l);
    double n_ijk = (nijk).norm();
    double n_jkl = (njkl).norm();
    double dotpr = nijk.dot(njkl);
    Eigen::Vector3d ji = j - i;
    double sign = (-1 * ji).dot(njkl) < 0 ? -1 : 1;
    double phi = pi + sign * acos(dotpr / (n_ijk * n_jkl));
    double energy = (1 / 2.0 * V * (1 - cos(n * phi0) * cos(n * phi)));
    if (std::isnan(energy))
        return 0;
    else
        return energy;
}

inline double Torsion(const Vector& i, const Vector& j, const Vector& k, const Vector& l, Matrix& derivate, bool gradient)
{
    Eigen::Vector3d rij = i - j;
    Eigen::Vector3d rjk = j - k;
    Eigen::Vector3d rkl = k - l;

    // Calculate normal vectors
    Eigen::Vector3d n1 = rij.cross(rjk);
    Eigen::Vector3d n2 = rjk.cross(rkl);

    double n1_norm = n1.norm();
    double n2_norm = n2.norm();

    if (n1_norm < 1e-8 || n2_norm < 1e-8) {
        if (gradient) {
            derivate = Matrix::Zero(4, 3);
        }
        return 0.0; // Linear arrangement
    }

    n1 /= n1_norm;
    n2 /= n2_norm;

    // Calculate torsion angle
    double cos_phi = n1.dot(n2);
    cos_phi = std::max(-1.0, std::min(1.0, cos_phi)); // Clamp to [-1,1]

    // Determine sign using triple scalar product
    double sign = (rij.cross(rjk)).dot(rkl) >= 0 ? 1.0 : -1.0;
    double phi = sign * acos(cos_phi);

    if (!gradient) {
        return phi;
    }

    // Calculate analytical derivatives
    derivate = Matrix::Zero(4, 3);

    double sin_phi = sin(phi);
    if (abs(sin_phi) < 1e-8) {
        // At extrema, derivatives are zero
        return phi;
    }

    double rjk_norm = rjk.norm();
    double rjk_norm2 = rjk_norm * rjk_norm;

    // Derivatives for each atom
    Eigen::Vector3d dphi_di = -rjk_norm / (n1_norm * n1_norm) * n1;
    Eigen::Vector3d dphi_dl = rjk_norm / (n2_norm * n2_norm) * n2;

    Eigen::Vector3d dphi_dj = -(1.0 - rij.dot(rjk) / rjk_norm2) * dphi_di
        - (rkl.dot(rjk) / rjk_norm2) * dphi_dl;

    Eigen::Vector3d dphi_dk = -(1.0 - rkl.dot(rjk) / rjk_norm2) * dphi_dl
        - (rij.dot(rjk) / rjk_norm2) * dphi_di;

    derivate.row(0) = dphi_di.transpose();
    derivate.row(1) = dphi_dj.transpose();
    derivate.row(2) = dphi_dk.transpose();
    derivate.row(3) = dphi_dl.transpose();

    return phi;
}

inline double OutOfPlane(const Vector& i, const Vector& j, const Vector& k, const Vector& l, Matrix& derivate, bool gradient)
{
    // Out-of-plane angle: angle between atom i and plane formed by j,k,l
    // l is the central atom, i is out-of-plane

    Eigen::Vector3d rlj = l - j;
    Eigen::Vector3d rlk = l - k;
    Eigen::Vector3d rli = l - i;

    // Normal vector to plane j-l-k
    Eigen::Vector3d normal = rlj.cross(rlk);
    double normal_norm = normal.norm();

    if (normal_norm < 1e-8) {
        if (gradient) {
            derivate = Matrix::Zero(4, 3);
        }
        return 0.0; // Degenerate plane
    }

    normal /= normal_norm;

    // Out-of-plane angle
    double rli_norm = rli.norm();
    if (rli_norm < 1e-8) {
        if (gradient) {
            derivate = Matrix::Zero(4, 3);
        }
        return 0.0;
    }

    double sin_theta = normal.dot(rli) / rli_norm;
    sin_theta = std::max(-1.0, std::min(1.0, sin_theta));
    double theta = asin(sin_theta);

    if (!gradient) {
        return theta;
    }

    // Calculate analytical derivatives
    derivate = Matrix::Zero(4, 3);

    double cos_theta = cos(theta);
    if (abs(cos_theta) < 1e-8) {
        // At 90 degrees, derivatives become complex
        return theta;
    }

    // Simplified derivatives (TODO: Complete implementation)
    // For now, using numerical approximation structure
    Eigen::Vector3d dtheta_di = cos_theta * normal / rli_norm;

    derivate.row(0) = dtheta_di.transpose(); // atom i
    // derivate.row(1) = ...; // atom j (TODO)
    // derivate.row(2) = ...; // atom k (TODO)
    // derivate.row(3) = ...; // atom l (TODO)

    return theta;
}
}

namespace QMDFF {

inline double LJStretchEnergy(double r_ij, double r0_ij, double fc, double exponent)
{
    const double ratio = r0_ij / r_ij;
    double energy = fc * (1 + pow(ratio, exponent) - 2 * pow(ratio, exponent * 0.5));
    if (std::isnan(energy))
        return 0;
    else
        return energy;
}

inline double AngleDamping(double r_ij, double r_ik, double r0_ij, double r0_ik)
{
    double kdamp = 1;
    double ratioij = r_ij / r0_ij;
    double ratioik = r_ik / r0_ik;

    double fijinv = 1 + kdamp * pow(ratioij, 4);
    double fikinv = 1 + kdamp * pow(ratioik, 4);
    double finv = fijinv * fikinv;
    return 1 / finv;
}

inline double AngleBend(const Eigen::Vector3d& j, const Eigen::Vector3d& i, const Eigen::Vector3d& k, double theta0_ijk, double fc, double r0_ij, double r0_ik, Matrix& derivate, bool gradient)
{
    Eigen::Vector3d ij = i - j;
    Eigen::Vector3d ik = i - k;
    double damp = 1; // AngleDamping(ij.norm(), ik.norm(), r0_ij, r0_ik);
    double costheta = (ij.dot(ik) / (sqrt(ij.dot(ij) * ik.dot(ik))));
    double costheta0_ijk = cos(theta0_ijk * pi / 180.0);
    double energy = (fc * damp * (costheta0_ijk - costheta) * (costheta0_ijk - costheta));

    if (std::isnan(energy))
        return 0;
    else
        return energy;
}

inline double LinearAngleBend(const Eigen::Vector3d& j, const Eigen::Vector3d& i, const Eigen::Vector3d& k, double theta0_ijk, double fc, double r0_ij, double r0_ik, Matrix& derivate, bool gradient)
{
#pragma message("check")

    Eigen::Vector3d ij = i - j;
    Eigen::Vector3d ik = i - k;
    double damp = AngleDamping(ij.norm(), ik.norm(), r0_ij, r0_ik);
    double theta = acos((ij.dot(ik) / (sqrt(ij.dot(ij) * ik.dot(ik)))));
    double energy = (fc * damp * (theta0_ijk - theta) * (theta0_ijk - theta));

    if (std::isnan(energy))
        return 0;
    else
        return energy;
}
}
