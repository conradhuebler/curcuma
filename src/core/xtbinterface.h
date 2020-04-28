/*
 * <Docking tool for structures. >
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

#include "external/xtb/include/xtb.h"

#include "src/core/molecule.h"

class XTBInterface {
public:
    XTBInterface();
    double GFN2Energy(const Molecule& molecule);
    double GFN1Energy(const Molecule& molecule);
    double GFN0Energy(const Molecule& molecule);

    double GFN2Energy(const int* attyp, const double* coord, const int natoms, const double charge);
    double GFN1Energy(const int* attyp, const double* coord, const int natoms, const double charge);
    double GFN0Energy(const int* attyp, const double* coord, const int natoms, const double charge);

    double GFN2Gradient(const int* attyp, const double* coord, const int natoms, const double charge, double* grad = 0);
    double GFN1Gradient(const int* attyp, const double* coord, const int natoms, const double charge, double* grad = 0);
    double GFN0Gradient(const int* attyp, const double* coord, const int natoms, const double charge, double* grad = 0);

private:
    double m_thr = 1.0e-10;
};
