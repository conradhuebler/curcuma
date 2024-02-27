/*
 * <QMDFF-Terms for force field calculation>
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

#include <Eigen/Dense>
/*
namespace QMDFFTerms {

double LJStretchEnergy(double r_ij, double r0_ij, double fc, double exponent)
{
    const double ratio = r0_ij / r_ij;
    double energy = fc * (1 + pow(ratio, exponent) - 2 * pow(ratio, exponent * 0.5));
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double AngleDamping(double r_ij, double r_ik, double r0_ij, double r0_ik)
{
    double kdamp = 1;
    double ratio_ij = r_ij / r0_ij;
    double ratio_ik = r_ik / r0_ik;

    double f_ij_inv = 1 + kdamp * pow(ratio_ij, 4);
    double f_ik_inv = 1 + kdamp * pow(ratio_ik, 4);
    double finv = f_ij_inv * f_ik_inv;
    return 1 / finv;
}

double AngleBend(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, double thetae, double fc, double r0_ij, double r0_ik)
{
    Eigen::Vector3d ij = i - j;
    Eigen::Vector3d ik = i - k;
    double damp = AngleDamping((ij).norm(), (ik).norm()), r0_ij, r0_ik);
    double costheta = (ij.dot(ik) / (sqrt(ij.dot(ij) * ik.dot(ik))));
    double costhetae = cos(thetae);
    double energy = (fc * damp * (costhetae - costheta) * (costhetae - costheta));
    if (isnan(energy))
        return 0;
    else
        return energy;
}

double LinearAngleBend(const Eigen::Vector3d& i, const Eigen::Vector3d& j, const Eigen::Vector3d& k, double thetae, double fc, double r0_ij, double r0_ik)
{
    if (fc < 0)
        return 0;
    Eigen::Vector3d ij = i - j;
    Eigen::Vector3d ik = i - k;
    double damp = AngleDamping((ij).norm(), (ik).norm()), r0_ij, r0_ik);
    double theta = acos(ij.dot(ik) / (sqrt(ij.dot(ij) * ik.dot(ik))));
    double energy = (fc * damp * (thetae - theta) * (thetae - theta));
    if (isnan(energy))
        return 0;
    else
        return energy;
}

}
*/
