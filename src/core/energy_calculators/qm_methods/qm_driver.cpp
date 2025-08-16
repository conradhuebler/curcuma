/*
 * <Extendend Hückel Theory Implementation in Cucuma. >
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "interface/abstract_interface.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"

#include "ParallelEigenSolver.hpp"

#include <set>
#include <vector>

#include <Eigen/Dense>

#include "GTOIntegrals.hpp"
#include "STOIntegrals.hpp"
#include "basissetparser.hpp"

#include "json.hpp"

#include "eht.h"

#include <iomanip>

QMDriver::QMDriver()
{
}

bool QMDriver::InitialiseMolecule()
{
    m_num_electrons = 0;
    m_mo = Matrix::Zero(m_atoms.size(), m_atoms.size());
    m_energies = Vector::Zero(m_atoms.size());
    return true;
}
