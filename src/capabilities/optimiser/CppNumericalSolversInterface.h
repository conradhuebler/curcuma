/*
 * <Interface to LBFGS Implementaton from Patrick Wieschollek. >
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

// see
// https://github.com/PatWie/CppNumericalSolvers for more information

#pragma once

#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

#include <functional>
#include <iostream>
#include <list>

#include "function.h"
#include "solver/lbfgs.h"

class CppNumSolvInterface : public cppoptlib::function::Function<double> {
public:
    using vector_t = typename cppoptlib::function::Function<double>::vector_t;

    double operator()(const vector_t& x) const override
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
        Molecule host(m_molecule);
        Geometry geometry = host.getGeometry();
        for (int i = 0; i < host.AtomCount(); ++i) {
            geometry(i, 0) = x(3 * i);
            geometry(i, 1) = x(3 * i + 1);
            geometry(i, 2) = x(3 * i + 2);
        }
        host.setGeometry(geometry);
        host.setEnergy(fx);
        host.appendXYZFile("move_host_eval.xyz");

        m_energy = fx;
        m_parameter = x;
        return fx;
    }

    void Gradient(const vector_t& x, vector_t* grad) const override
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

        Molecule host(m_molecule);
        Geometry geometry = host.getGeometry();
        for (int i = 0; i < host.AtomCount(); ++i) {
            geometry(i, 0) = x(3 * i);
            geometry(i, 1) = x(3 * i + 1);
            geometry(i, 2) = x(3 * i + 2);
        }
        host.setGeometry(geometry);
        host.setEnergy(fx);
        host.appendXYZFile("move_host_grad.xyz");

        for (int i = 0; i < m_atoms; ++i) {
            (*grad)[3 * i + 0] = gradient[3 * i + 0];
            (*grad)[3 * i + 1] = gradient[3 * i + 1];
            (*grad)[3 * i + 2] = gradient[3 * i + 2];
        }
        m_energy = fx;
        m_parameter = x;
    }

    double LastEnergy() const { return m_energy; }

    mutable double m_energy = 0, m_last_change = 0, m_last_rmsd = 0;
    Vector Parameter() const { return m_parameter; }
    void setMolecule(const Molecule* molecule)
    {
        m_molecule = molecule;
        m_atoms = m_molecule->AtomCount();
    }
    void setInterface(XTBInterface* interface) { m_interface = interface; }
    void setMethod(int method) { m_method = method; }

private:
    int m_iter = 0;
    int m_atoms = 0;
    int n;
    int m_method = 2;
    XTBInterface* m_interface;
    mutable Vector m_parameter;
    const Molecule* m_molecule;
};
