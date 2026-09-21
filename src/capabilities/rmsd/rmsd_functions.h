/*
 * <Collection of functions to calculate rmsd. >
 * Copyright (C) 2020 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "src/core/molecule.h"

#include "src/tools/geometry.h"

#include <Eigen/Dense>

// Claude Generated 2026 - GCC-only optimization attribute. MSVC has no
// __attribute__ at all (misparses the declaration -> C3861/C4430), and Clang has
// no 'optimize' attribute either (warns "unknown attribute ignored"). Guard it so
// only real GCC emits the attribute; everyone else gets the plain declaration.
#if defined(__GNUC__) && !defined(__clang__)
#define CURCUMA_NO_TREE_VECTORIZE __attribute__((optimize("no-tree-vectorize")))
#else
#define CURCUMA_NO_TREE_VECTORIZE
#endif

namespace RMSDFunctions {

/*! \brief Calculate the best fit rotation of two sets of coordinates, both have to be centered already */

CURCUMA_NO_TREE_VECTORIZE inline Eigen::Matrix3d BestFitRotation(const Geometry& reference, const Geometry& target, int factor = 1)
{
    /* The rmsd kabsch algorithmn was adopted from here:
     * https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
     * The specific git commit was
     * https://github.com/oleg-alexandrov/projects/blob/e7b1eb7a4d83d41af563c24859072e4ddd9b730b/eigen/Kabsch.cpp
     */
    // Direkte Berechnung der Kovarianzmatrix ohne temporäre Matrizen
    Eigen::Matrix3d Cov;
    Cov.setZero();

    const int n = reference.rows();
    for (int i = 0; i < n; ++i) {
        // Verwende Eigen's Vektorisierung für diese Operation
        Cov += reference.row(i).transpose() * target.row(i);
    }

    // Sichere Determinantenberechnung für numerische Stabilität
    double detCov = Cov.determinant();

    // Spezieller Fall: Wenn Cov nahe an einer Nullmatrix ist
    if (std::abs(detCov) < 1e-10) {
        // Rückgabe einer Matrix, die nah an der Identität ist
        return factor * Eigen::Matrix3d::Identity();
    }

    // Optimierte SVD für 3x3 Matrix
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(Cov, Eigen::ComputeFullU | Eigen::ComputeFullV);

    // Direkte Konstruktion der Rotationsmatrix
    Eigen::Matrix3d V = svd.matrixV();
    Eigen::Matrix3d U = svd.matrixU();

    // Handhabung von Reflexionen für eine echte Rotationsmatrix
    double d = (V * U.transpose()).determinant();
    if (d < 0) {
        // Nur die letzte Spalte von V modifizieren
        V.col(2) *= -1.0;
    }

    return factor * V * U.transpose();
}

inline Eigen::Matrix3d BestFitRotation(const Molecule& reference, const Molecule& target, int factor = 1)
{
    return BestFitRotation(reference.getGeometry(), target.getGeometry(), factor);
}

inline Geometry applyRotation(const Geometry& geometry, const Eigen::Matrix3d& rotation)
{
    return geometry * rotation;
}

inline Geometry getAligned(const Geometry& reference, const Geometry& target, int factor)
{
    Eigen::Matrix3d rotation = BestFitRotation(reference, target, factor);
    return applyRotation(target, rotation);
}

inline Molecule getAligned(const Molecule& reference, const Molecule& target, int factor)
{
    Molecule result = target;
    result.setGeometry(getAligned(reference.getGeometry(), target.getGeometry(), factor));
    return result;
}

inline double getRMSD(const Geometry& reference, const Geometry& target)
{
    double rmsd = 0.0;
    for (int i = 0; i < target.rows(); ++i) {
        rmsd += (target(i, 0) - reference(i, 0)) * (target(i, 0) - reference(i, 0))
            + (target(i, 1) - reference(i, 1)) * (target(i, 1) - reference(i, 1))
            + (target(i, 2) - reference(i, 2)) * (target(i, 2) - reference(i, 2));
    }
    rmsd = sqrt(rmsd / double(target.rows()));
    return rmsd;
}

/* Claude Generated (Jun 2026): weighted Kabsch variants for flexibility-weighted RMSD.
 * With uniform weights (all equal) these reduce EXACTLY to the unweighted functions above,
 * so the default ConfSearch/MTD path stays bit-identical. Per-atom weight w_i down-weights
 * floppy atoms (large positional fluctuation) so the RMSD reflects rigid-core conformer
 * identity. The rotation derivative is neglected (same approximation as the analytic
 * RMSDDriver::Gradient()). */

inline Eigen::Vector3d WeightedCentroid(const Geometry& g, const std::vector<double>& w)
{
    Eigen::Vector3d c = Eigen::Vector3d::Zero();
    double wsum = 0.0;
    for (int i = 0; i < g.rows(); ++i) {
        c += w[i] * g.row(i).transpose();
        wsum += w[i];
    }
    if (wsum > 0)
        c /= wsum;
    return c;
}

CURCUMA_NO_TREE_VECTORIZE inline Eigen::Matrix3d BestFitRotationW(const Geometry& reference, const Geometry& target, const std::vector<double>& w)
{
    Eigen::Matrix3d Cov;
    Cov.setZero();
    const int n = reference.rows();
    for (int i = 0; i < n; ++i)
        Cov += w[i] * reference.row(i).transpose() * target.row(i);

    if (std::abs(Cov.determinant()) < 1e-10)
        return Eigen::Matrix3d::Identity();

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(Cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d V = svd.matrixV();
    Eigen::Matrix3d U = svd.matrixU();
    if ((V * U.transpose()).determinant() < 0)
        V.col(2) *= -1.0;
    return V * U.transpose();
}

inline double getRMSDW(const Geometry& reference, const Geometry& target, const std::vector<double>& w)
{
    double rmsd = 0.0, wsum = 0.0;
    for (int i = 0; i < target.rows(); ++i) {
        double d2 = (target(i, 0) - reference(i, 0)) * (target(i, 0) - reference(i, 0))
            + (target(i, 1) - reference(i, 1)) * (target(i, 1) - reference(i, 1))
            + (target(i, 2) - reference(i, 2)) * (target(i, 2) - reference(i, 2));
        rmsd += w[i] * d2;
        wsum += w[i];
    }
    return sqrt(rmsd / (wsum > 0 ? wsum : 1.0));
}
};
