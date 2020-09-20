/*
 * <Geometrie optimisation using external LBFGS and xtb. >
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

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <iomanip>
#include <iostream>

#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

#include "src/tools/geometry.h"

#include "json.hpp"
#include <LBFGSB.h>
using json = nlohmann::json;

const json OptJson2{
    { "writeXYZ", false },
    { "printOutput", true },
    { "dE", 0.75 },
    { "dRMSD", 0.01 },
    { "GFN", 2 },
    { "InnerLoop", 20 },
    { "OuterLoop", 100 },
    { "LBFGS_eps", 1e-5 }
};

using Eigen::VectorXd;
using namespace LBFGSpp;
