/*
 * < Generic force field class for curcuma . >
 * Copyright (C) 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "forcefieldthread.h"

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

#include <functional>
#include <set>
#include <vector>

#include <Eigen/Dense>

#include "json.hpp"
using json = nlohmann::json;

class ForceField {

public:
    ForceField(const json& controller);
    ~ForceField();

    inline void setAtomCount(int atom) { m_natoms = atom; }
    void UpdateGeometry(const Matrix& geometry);
    inline void UpdateGeometry(const double* coord);
    inline void UpdateGeometry(const std::vector<std::array<double, 3>>& geometry);

    double Calculate(bool gradient = true, bool verbose = false);

    Matrix Gradient() const { return m_gradient; }

    void setParameter(const json& parameter);
    void setParameterFile(const std::string& file);
    Eigen::MatrixXd NumGrad();

private:
    void AutoRanges();
    void setBonds(const json& bonds);
    void setAngles(const json& angles);
    void setDihedrals(const json& dihedrals);
    void setInversions(const json& inversions);

    std::vector<ForceFieldThread*> m_stored_threads;
    CxxThreadPool* m_threadpool;
    void setvdWs(const json& vdws);

    Matrix m_geometry, m_gradient;
    int m_natoms = 0;
    int m_threads = 1;
    std::vector<Bond> m_bonds;
    std::vector<Angle> m_angles;
    std::vector<Dihedral> m_dihedrals;
    std::vector<Inversion> m_inversions;
    std::vector<vdW> m_vdWs;
    std::vector<EQ> m_EQs;
};
