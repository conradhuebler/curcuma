/*
 * <Geometrie optimisation using external LBFGS and xtb. >
 * Copyright (C) 2020 - 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <iomanip>
#include <iostream>

#include "src/core/energycalculator.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

#include "json.hpp"
#include <LBFGS.h>
#include <LBFGSB.h>

using json = nlohmann::json;

using Eigen::VectorXd;
using namespace LBFGSpp;

class LBFGSInterface {
public:
    LBFGSInterface(int n_)
    //: n(n_)
    {
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        double charge = 0;
        m_interface->updateGeometry(x);
        fx = m_interface->CalculateEnergy(true);
        auto gradient = m_interface->getGradient();
        m_error = std::isnan(fx);

        for (int i = 0; i < m_atoms; ++i) {
            grad[3 * i + 0] = gradient[i][0] * (m_constrains[i]);
            grad[3 * i + 1] = gradient[i][1] * (m_constrains[i]);
            grad[3 * i + 2] = gradient[i][2] * (m_constrains[i]);
        }
        m_energy = fx;
        m_parameter = x;
        return fx;
    }

    double LastEnergy() const { return m_energy; }

    double m_energy = 0, m_last_change = 0, m_last_rmsd = 0;
    Vector Parameter() const { return m_parameter; }
    void setMolecule(const Molecule* molecule)
    {
        m_molecule = molecule;
        m_atoms = m_molecule->AtomCount();
        for (int i = 0; i < m_atoms; ++i)
            m_constrains.push_back(1);
    }
    void setConstrains(const std::vector<int> constrains) { m_constrains = constrains; }
    void setInterface(EnergyCalculator* interface) { m_interface = interface; }
    void setMethod(int method) { m_method = method; }
    bool isError() const { return m_error; }

private:
    //  int m_iter = 0;
    int m_atoms = 0;
    //  int n;
    int m_method = 2;
    std::vector<int> m_constrains;
    EnergyCalculator* m_interface;
    Vector m_parameter;
    const Molecule* m_molecule;
    bool m_error = false;
};
