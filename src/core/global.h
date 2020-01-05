/*
 * <Some globale definition for chemical structures.>
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

const double pi = 3.14159265359;

typedef Eigen::MatrixXd Geometry;
typedef Eigen::Vector3d Position;
typedef Eigen::VectorXd Vector;

inline Vector PositionPair2Vector(const std::pair<Position, Position>& pair)
{
    Vector vector = Vector::Zero(6);
    vector(0) = pair.first(0);
    vector(1) = pair.first(1);
    vector(2) = pair.first(2);
    vector(3) = pair.second(0);
    vector(4) = pair.second(1);
    vector(5) = pair.second(2);
    return vector;
}

inline Vector PositionPair2Vector(const Position& first, const Position& second)
{
    Vector vector = Vector::Zero(6);
    vector(0) = first(0);
    vector(1) = first(1);
    vector(2) = first(2);
    vector(3) = second(0);
    vector(4) = second(1);
    vector(5) = second(2);
    return vector;
}

inline std::pair<Position, Position> Vector2PositionPair(const Vector& vector)
{
    return std::pair<Position, Position>(Position{ vector(0), vector(1), vector(2) }, Position{ vector(3), vector(4), vector(5) });
}
