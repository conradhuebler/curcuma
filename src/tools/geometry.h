/*
 * <Geometry tools for chemical structures.>
 * Copyright (C) 2019 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/global.h"
#include "src/core/molecule.h"
using namespace curcuma;
namespace GeometryTools {

inline double degreesToRadians(double angleDegrees) { return angleDegrees * pi / 180.0; }

inline double radiansToDegrees(double angleRadians) { return angleRadians * 180.0 / pi; }

inline double Distance(const Position& a, const Position& b)
{
    double distance = 0;
    distance = sqrt((a(0) - b(0)) * (a(0) - b(0)) + (a(1) - b(1)) * (a(1) - b(1)) + (a(2) - b(2)) * (a(2) - b(2)));
    return distance;
}

/*! \brief Angle a-b-c at the middle atom b, in degrees (0 ... 180).
 *
 * Claude Generated (Sep 2026). The cosine is clamped before acos: rounding can push
 * the dot product a hair outside [-1, 1] for a straight or fully folded arrangement,
 * which would return NaN.
 */
inline double Angle(const Position& a, const Position& b, const Position& c)
{
    const Position u = a - b;
    const Position v = c - b;
    const double nu = u.norm(), nv = v.norm();
    if (nu < 1e-12 || nv < 1e-12)
        return 0.0;
    double cosine = u.dot(v) / (nu * nv);
    cosine = std::max(-1.0, std::min(1.0, cosine));
    return radiansToDegrees(std::acos(cosine));
}

/*! \brief Dihedral a-b-c-d in degrees, signed (-180 ... 180, IUPAC convention).
 *
 * Claude Generated (Sep 2026). Standard atan2 form over the two plane normals, so
 * the sign distinguishes the two enantiomeric arrangements instead of folding them
 * onto one another.
 */
inline double Dihedral(const Position& a, const Position& b, const Position& c, const Position& d)
{
    const Position b1 = b - a, b2 = c - b, b3 = d - c;
    const double n2 = b2.norm();
    if (n2 < 1e-12)
        return 0.0;
    const Position nn1 = b1.cross(b2);
    const Position nn2 = b2.cross(b3);
    const Position m = nn1.cross(b2 / n2);
    return radiansToDegrees(std::atan2(m.dot(nn2), nn1.dot(nn2)));
}

inline Position Centroid(const Geometry& geom)
{
    Position position{ 0, 0, 0 };

    for (int i = 0; i < geom.rows(); ++i) {
        position += geom.row(i);
    }

    position /= double(geom.rows());

    return position;
}

/*! \brief Radius of gyration of a bare geometry, in the geometry's length unit.
 *
 * Claude Generated (Sep 2026). Unweighted: every position counts once, so this is a
 * pure shape measure of a coordinate set and needs no element information. For the
 * mass-weighted radius (and the PBC-aware and per-fragment variants) use
 * Molecule::GyrationRadius, which has the masses.
 */
inline double GyrationRadius(const Geometry& geom)
{
    if (geom.rows() == 0)
        return 0.0;
    const Position centre = Centroid(geom);
    double sum = 0.0;
    for (int i = 0; i < geom.rows(); ++i) {
        const Position d = Position(geom.row(i)) - centre;
        sum += d.squaredNorm();
    }
    return std::sqrt(sum / double(geom.rows()));
}

inline Geometry RotationX(double alpha)
{
    const double radian = degreesToRadians(alpha);
    Geometry rotation = Geometry::Zero(3, 3);

    rotation(0, 0) = 1;
    rotation(1, 1) = std::cos(radian);
    rotation(1, 2) = -1 * std::sin(radian);
    rotation(2, 1) = std::sin(radian);
    rotation(2, 2) = std::cos(radian);

    return rotation;
}

inline Geometry RotationY(double alpha)
{
    const double radian = degreesToRadians(alpha);
    Geometry rotation = Geometry::Zero(3, 3);

    rotation(0, 0) = std::cos(radian);
    rotation(0, 2) = std::sin(radian);
    rotation(1, 1) = 1;
    rotation(2, 0) = -1 * std::sin(radian);
    rotation(2, 2) = std::cos(radian);

    return rotation;
}

inline Geometry RotationZ(double alpha)
{
    const double radian = degreesToRadians(alpha);
    Geometry rotation = Geometry::Zero(3, 3);

    rotation(0, 0) = std::cos(radian);
    rotation(0, 1) = -1 * std::sin(radian);
    rotation(1, 0) = std::sin(radian);
    rotation(1, 1) = std::cos(radian);
    rotation(2, 2) = 1;

    return rotation;
}

inline Position rotateX(const Position& position, double alpha)
{
    const double radian = degreesToRadians(alpha);

    const double X = position(0);
    const double Y = position(1) * std::cos(radian) - position(2) * std::sin(radian);
    const double Z = position(1) * std::sin(radian) + position(2) * std::cos(radian);

    return Position{ X, Y, Z };
}

inline Position rotateY(const Position& position, double alpha)
{
    const double radian = degreesToRadians(alpha);

    const double X = position(0) * std::cos(radian) + position(2) * std::sin(radian);
    const double Y = position(1);
    const double Z = -1 * position(0) * std::sin(radian) + position(2) * std::cos(radian);

    return Position{ X, Y, Z };
}

inline Position rotateZ(const Position& position, double alpha)
{
    const double radian = degreesToRadians(alpha);

    const double X = position(0) * std::cos(radian) - position(1) * std::sin(radian);
    const double Y = position(0) * std::sin(radian) + position(1) * std::cos(radian);
    const double Z = position(2);

    return Position{ X, Y, Z };
}

inline Geometry TranslateMolecule(const Molecule& molecule, const Position& start, const Position& destination)
{
    Geometry geom = molecule.getGeometry();
    Position direction = destination - start;
    for(int i = 0; i < geom.rows(); ++i)
        geom.row(i) += direction;

    return geom;
}

inline Geometry TranslateMolecule(const Molecule& molecule, const Position& translate)
{
    Geometry geom = molecule.getGeometry();
    for (int i = 0; i < geom.rows(); ++i)
        geom.row(i) += translate;

    return geom;
}

inline Geometry TranslateGeometry(const Geometry& geom, const Position& start, const Position& destination)
{
    Position direction = destination - start;
    Geometry temp = geom;
    for (int i = 0; i < geom.rows(); ++i)
        temp.row(i) += direction;
    return temp;
}

inline Geometry TranslateGeometry(const Geometry& geom, const Position& translate)
{
    Geometry temp = geom;
    for (int i = 0; i < geom.rows(); ++i)
        temp.row(i) += translate;

    return temp;
}

inline Geometry TranslateAndRotate(const Geometry& geom, const Position& start, const Position& destination, const Position& rotation)
{
    Geometry temp = GeometryTools::TranslateGeometry(geom, start, Position{ 0, 0, 0 });
    Geometry rot = GeometryTools::RotationX(rotation(0));
    rot *= GeometryTools::RotationY(rotation(1));
    rot *= GeometryTools::RotationZ(rotation(2));
    temp = temp * rot;
    return GeometryTools::TranslateGeometry(temp, Position{ 0, 0, 0 }, destination);
}
}
