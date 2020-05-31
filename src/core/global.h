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

#include "src/global_config.h"
#include "src/version.h"

#include "json.hpp"

// for convenience
using json = nlohmann::json;

const double pi = 3.14159265359;
const double au = 0.52917721092;
const double amu2au = 1822.8884850;
const double kb = 3.166811e-6;
const double fs2amu = 41.34137314;

typedef Eigen::MatrixXd Geometry;
typedef Eigen::Vector3d Position;
typedef Eigen::Vector4d Vector4d;

typedef Eigen::VectorXd Vector;
typedef std::pair<int, int> IntPair;
typedef std::vector<std::string> StringList;

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

inline json CLI2Json(int argc, char** argv)
{
    json controller;
    json key;
    if (argc < 2)
        return controller;
    std::string keyword = argv[1];
    keyword.erase(0, 1);
    for (int i = 2; i < argc; ++i) {
        std::string current = argv[i];
        std::string sub = current.substr(0, 1);
        if ((i + 1) >= argc) {
            current.erase(0, 1);
            if (sub.compare("-") == 0)
                key[current] = true;
            else
                key[current] = false;
        } else {
            if (sub.compare("-") == 0 && ((i + 1) < argc)) {
                std::string next = argv[i + 1];
                std::string next_sub = next.substr(0, 1);
                if (next_sub.compare("-") == 0) {
                    current.erase(0, 1);
                    key[current] = true;
                    continue;
                } else if (next_sub.compare("+") == 0) {
                    current.erase(0, 1);
                    key[current] = false;
                    continue;
                } else {
                    current.erase(0, 1);
                    try {
                        key[current] = std::stoi(argv[i + 1]);
                    } catch (const std::invalid_argument& error) {
                        try {
                            key[current] = std::stod(argv[i + 1]);
                        } catch (const std::invalid_argument& error) {
                            key[current] = argv[i + 1];
                        }
                    }

                    ++i;
                }
            }
        }
    }
    controller[keyword] = key;
    return controller;
}

template <class T>
inline T Json2KeyWord(const json& controller, std::string name)
{
    T temp;
    bool found = false;
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    for (auto& el : controller.items()) {
        std::string key = el.key();
        transform(key.begin(), key.end(), key.begin(), ::tolower);
        if (key.compare(name) == 0) {
            temp = el.value();
            found = true;
        }
    }
    if (found)
        return temp;
    else
        throw -1;
}
