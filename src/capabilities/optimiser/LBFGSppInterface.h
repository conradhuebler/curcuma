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

#include <iomanip>
#include <iostream>

#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/tbliteinterface.h"
#include "src/core/xtbinterface.h"

#include "json.hpp"
#include <LBFGS.h>
#include <LBFGSB.h>

using json = nlohmann::json;

using Eigen::VectorXd;
using namespace LBFGSpp;

class LBFGSInterface {
public:
    LBFGSInterface(int n_)
        : n(n_)
    {
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        int attyp[m_atoms];
        double charge = 0;
        std::vector<int> atoms = m_molecule->Atoms();
        double coord[3 * m_atoms];
        double gradient[3 * m_atoms];
        double dist_gradient[3 * m_atoms];

        for (int i = 0; i < m_atoms; ++i) {
            attyp[i] = m_molecule->Atoms()[i];
            coord[3 * i + 0] = x(3 * i + 0) / au;
            coord[3 * i + 1] = x(3 * i + 1) / au;
            coord[3 * i + 2] = x(3 * i + 2) / au;
        }
        m_interface->UpdateMolecule(coord);
        fx = m_interface->GFNCalculation(m_method, gradient);

        for (int i = 0; i < m_atoms; ++i) {
            grad[3 * i + 0] = gradient[3 * i + 0];
            grad[3 * i + 1] = gradient[3 * i + 1];
            grad[3 * i + 2] = gradient[3 * i + 2];
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
    }
    void setInterface(TBLiteInterface* interface) { m_interface = interface; }
    void setMethod(int method) { m_method = method; }

private:
    int m_iter = 0;
    int m_atoms = 0;
    int n;
    int m_method = 2;
    // XTBInterface* m_interface;
    TBLiteInterface* m_interface;
    Vector m_parameter;
    const Molecule* m_molecule;
};
