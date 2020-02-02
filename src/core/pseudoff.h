/*
 * <Pseudo ForceField Implementation.>
 * Copyright (C) 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include "src/tools/geometry.h"

class PseudoFF {
public:
    PseudoFF();

    /*! \brief Calculate the LennardJones Potential between two molecules */
    static double LennardJones(const Molecule& a, const Molecule& b);

    inline static double LJPotential(double distance, double vdW_A, double vdW_B)
    {
        const double potenz = pow((vdW_A + vdW_B) / (distance), 6);
        return potenz * (potenz - 2);
    }

    inline static double LJPotential(const std::pair<int, Position>& a, const std::pair<int, Position>& b)
    {
        const double distance = GeometryTools::Distance(a.second, b.second);
        const double vdW_A = Elements::VanDerWaalsRadius[a.first];
        const double vdW_B = Elements::VanDerWaalsRadius[b.first];

        return LJPotential(distance, vdW_A, vdW_B);
    }

    inline static double DistancePenalty(const std::pair<int, Position>& a, const std::pair<int, Position>& b, double scaling = 1.5)
    {
        double error = 0;

        const double distance = GeometryTools::Distance(a.second, b.second);
        const double cov_A = Elements::CovalentRadius[a.first];
        const double cov_B = Elements::CovalentRadius[b.first];
        double mind = (cov_A + cov_B) * scaling;

        error = (distance < mind) * (mind - distance);
        return error;
    }
    /*
    inline static Vector4d DistancePenaltyDiff(const std::pair<int, Position>& a, const std::pair<int, Position>& b, double scaling = 1.5)
    {
        double distance = 0;
        double dx = 0, dy = 0, dz = 0;
        dx = (a.second(0) - b.second(0)) * (a.second(0) - b.second(0));
        dy = (a.second(1) - b.second(1)) * (a.second(1) - b.second(1));
        dz = (a.second(2) - b.second(2)) * (a.second(2) - b.second(2));
        distance = sqrt(dx + dy + dz);
        double ddx = 0, ddy = 0, ddz = 0;
        double error = 0;
        const double cov_A = Elements::CovalentRadius[a.first];
        const double cov_B = Elements::CovalentRadius[b.first];
        double mind = (cov_A + cov_B) * scaling;

        error = (distance < mind) * (mind - distance);
        ddx = (2 * dx / distance < mind) * (mind - distance);
        ddy = (2 * dy / distance < mind) * (mind - distance);
        ddz = (2 * dz / distance < mind) * (mind - distance);

        return Vector4d{ error, ddx, ddy, ddz };
    }
    */
    inline static Vector4d DistancePenaltyDiff(const std::pair<int, Position>& a, const std::pair<int, Position>& b, double scaling = 1.5)
    {
        double distance = 0;
        double dx = 0, dy = 0, dz = 0;
        const double vdW_A = Elements::VanDerWaalsRadius[a.first];
        const double vdW_B = Elements::VanDerWaalsRadius[b.first];

        dx = (a.second(0) - b.second(0)) * (a.second(0) - b.second(0));
        dy = (a.second(1) - b.second(1)) * (a.second(1) - b.second(1));
        dz = (a.second(2) - b.second(2)) * (a.second(2) - b.second(2));

        distance = sqrt(dx + dy + dz);
        double ddx = 0, ddy = 0, ddz = 0;
        double error = 2 * pow((vdW_A + vdW_B) / (distance), 2);

        ddx = dx * 48 * pow((vdW_A + vdW_B) / (sqrt(dx + dy + dz)), 3) / (vdW_A + vdW_B);
        ddy = dy * 48 * pow((vdW_A + vdW_B) / (sqrt(dy + dx + dz)), 3) / (vdW_A + vdW_B);
        ddz = dz * 48 * pow((vdW_A + vdW_B) / (sqrt(dz + dx + dy)), 3) / (vdW_A + vdW_B);

        return Vector4d{ error, ddx, ddy, ddz };
    }
};
