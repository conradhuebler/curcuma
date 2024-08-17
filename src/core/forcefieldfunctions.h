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

#include "src/core/global.h"
#include "src/core/qmdff_par.h"
#include "src/core/uff_par.h"

#include <cmath>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
/* thanks to ChatGPT Maths module */
// Funktion zur Berechnung des Kosinus des Winkels zwischen zwei Vektoren
inline double cosTheta(const Eigen::Vector3d& AB, const Eigen::Vector3d& AC)
{
    double dotProduct = AB.dot(AC);
    double norms = AB.norm() * AC.norm();
    return dotProduct / norms;
}

// Funktion zur Berechnung des Gradienten von cos(theta) bezüglich eines Punktes
inline Eigen::Vector3d gradientCosTheta(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C)
{
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AC = C - A;

    double cosThetaValue = cosTheta(AB, AC);
    double normAB = AB.norm();
    double normAC = AC.norm();

    // Berechnung der Gradientenkomponenten
    Eigen::Vector3d gradA = (AC / (normAB * normAC)) - (cosThetaValue / normAB) * (AB / normAB)
        - (cosThetaValue / normAC) * (AC / normAC);

    return gradA;
}

// Funktion zur Berechnung des Gradienten von E
inline double computeGradient(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C,
    double fc, double cosTheta0,
    Eigen::Vector3d& gradA, Eigen::Vector3d& gradB, Eigen::Vector3d& gradC)
{
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AC = C - A;

    // Berechne cos(theta)
    double cosThetaValue = cosTheta(AB, AC);
    // std::cout << acos(cosThetaValue)*180/pi << " " << acos(cosTheta0)*180/pi<< std::endl;
    //  Berechne den Faktor (cos(theta) - cos(theta_0))
    double factor = 2 * fc * (cosThetaValue - cosTheta0);

    // Berechne die Gradienten in Bezug auf A, B und C
    gradA = factor * gradientCosTheta(A, B, C);
    gradB = -factor * gradientCosTheta(B, A, C); // wegen AB wird zu BA und umgekehrt
    gradC = -factor * gradientCosTheta(C, B, A); // wegen AC wird zu CA und umgekehrt
    return fc * (cosThetaValue - cosTheta0) * (cosThetaValue - cosTheta0);
}

inline double harmonicPotential(double theta, double theta0, double k)
{
    return 0.5 * k * (theta - theta0) * (theta - theta0);
}

inline double computeTheta(const Eigen::Vector3d& AB, const Eigen::Vector3d& AC)
{
    Eigen::Vector3d crossProduct = AB.cross(AC);
    double sinTheta = crossProduct.norm() / (AB.norm() * AC.norm());
    return std::asin(sinTheta); // kleine Winkelannahme
}

inline double computeHarmonicGradient(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C,
    double k, double theta0,
    Eigen::Vector3d& gradA, Eigen::Vector3d& gradB, Eigen::Vector3d& gradC)
{
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AC = C - A;

    double theta = computeTheta(AB, AC);
    double factor = -k * (theta - theta0);

    // Annahme: grad_theta/grad_r proportional zu senkrechtem Abstand
    Eigen::Vector3d crossProduct = AB.cross(AC);
    gradA = factor * crossProduct.normalized();
    gradB = -gradA; // Wegen der Definition des Abstands
    gradC = -gradA;
    return harmonicPotential(theta, theta0, k);
}

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
    // std::cout << acos(costheta)*180/pi << " " << theta0_ijk << std::endl;
    // std::cout << damp << " " << costheta << " " << costheta0_ijk << " " << energy << std::endl;

    if (!gradient) {
        if (std::isnan(energy))
            return 0;
        else
            return energy;
    }
    Eigen::Vector3d rij = i - j;
    auto nij = rij / rij.norm();
    Eigen::Vector3d rkj = k - j;
    auto nkj = rkj / rkj.norm();
    derivate = Matrix::Zero(3, 3);

    double sintheta = sin(acos(costheta));
    double dThetadCosTheta = 1 / sintheta;
    double dE = fc * (sintheta * (2 * costheta + costheta0_ijk));
    derivate = Matrix::Zero(3, 3);
    derivate.row(0) = -dThetadCosTheta * (nkj - nij * costheta) / (rij.norm()) * dE;
    derivate.row(2) = -dThetadCosTheta * (nij - nkj * costheta) / (rkj.norm()) * dE;
    derivate.row(1) = -derivate.row(0) - derivate.row(2) * dE;

    if (std::isnan(energy))
        return 0;
    else
        return energy;
}

inline double LinearAngleBend(const Eigen::Vector3d& j, const Eigen::Vector3d& i, const Eigen::Vector3d& k, double theta0_ijk, double fc, double r0_ij, double r0_ik, Matrix& derivate, bool gradient)
{
#pragma message("check")
    /*
    if (kabc < 0)
        return 0;
*/
    Eigen::Vector3d ij = i - j;
    Eigen::Vector3d ik = i - k;
    double damp = AngleDamping(ij.norm(), ik.norm(), r0_ij, r0_ik);
    double theta = acos((ij.dot(ik) / (sqrt(ij.dot(ij) * ik.dot(ik)))));
    double energy = (fc * damp * (theta0_ijk - theta) * (theta0_ijk - theta));

    if (!gradient) {
        if (std::isnan(energy))
            return 0;
        else
            return energy;
    }
    Eigen::Vector3d rij = i - j;
    auto nij = rij / rij.norm();
    Eigen::Vector3d rkj = k - j;
    auto nkj = rkj / rkj.norm();
    derivate = Matrix::Zero(3, 3);
    double costheta = (rij.dot(rkj) / (sqrt(rij.dot(rij) * rkj.dot(rkj))));

    double sintheta = sin(acos(costheta));
    double dThetadCosTheta = 1 / sintheta;
    double dE = 2 * fc * (theta0_ijk - theta);
    derivate = Matrix::Zero(3, 3);
    derivate.row(0) = -dThetadCosTheta * (nkj - nij * costheta) / (rij.norm()) * dE;
    derivate.row(2) = -dThetadCosTheta * (nij - nkj * costheta) / (rkj.norm()) * dE;
    derivate.row(1) = -derivate.row(0) - derivate.row(2) * dE;

    if (std::isnan(energy))
        return 0;
    else
        return energy;
}
}
