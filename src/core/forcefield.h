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

#include "hbonds.h"

#include "external/CxxThreadPool/include/CxxThreadPool.h"

#ifdef USE_D3
#include "src/core/dftd3interface.h"
#endif

#ifdef USE_D4
#include "src/core/dftd4interface.h"
#endif

#include "src/core/qmdff_par.h"
#include "src/core/uff_par.h"

#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

struct Bond {
    int i = 0, j = 0, k = 0, distance = 0;
    double fc = 0, exponent = 0, r0_ij = 0, r0_ik = 0;
};

struct Angle {
    int i = 0, j = 0, k = 0;
    double fc = 0, r0_ij = 0, r0_ik = 0, theta0_ijk = 0;
};

class ForceFieldThread : public CxxThread {

public:
    ForceFieldThread();

private:
    double CalculateBondContribution();

    double HarmonicBondStretching();

    double LJBondStretching();

    double m_final_factor = 1;
    double m_bond_scaling = 1;
    double m_au = 1;
    double m_d = 1e-3;
    int m_calc_gradient = 1;
    bool m_calculate_gradient = true;

    Matrix m_geometry, m_gradient;

    std::vector<Bond> m_bonds;
};

class ForceField {

public:
    ForceField(const json& controller);
    ~ForceField();

    void UpdateGeometry(const Matrix& geometry);
    inline void UpdateGeometry(const double* coord);
    inline void UpdateGeometry(const std::vector<std::array<double, 3>>& geometry);

    double Calculate(bool gradient = true, bool verbose = false);

    Matrix Gradient() const { return m_gradient; }

    void setParameter(const json& parameter);

private:
    Matrix m_geometry, m_gradient;
    int m_natoms = 0;
};
