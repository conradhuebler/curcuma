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
#include <unsupported/Eigen/NonLinearOptimization>

#include <iomanip>
#include <iostream>

#include "src/capabilities/rmsd.h"

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/xtbinterface.h"

#include "src/tools/geometry.h"

#include <external/LBFGSpp/include/LBFGS.h>

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
        fx = interface->GFN2Gradient(attyp, coord, m_atoms, charge, gradient);
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

private:
    int m_iter = 0;
    int m_atoms = 0;
    int n;
    XTBInterface* interface;
    Vector m_parameter;
    const Molecule* m_molecule;
};

inline Molecule OptimiseGeometry(const Molecule* host, bool writeXYZ = true, bool printOutput = true, double dE = 0.75, double dRMSD = 0.01)
{
    Geometry geometry = host->getGeometry();
    Molecule tmp(host);
    Molecule h(host);
    Vector parameter(3 * host->AtomCount());

    for (int i = 0; i < host->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now(), end;
    std::cout << "Step\tCurrent Energy [Eh]\tEnergy Change\tRMSD Change\tt [s]" << std::endl;
    int maxouter = 100;

    XTBInterface interface;
    double final_energy = interface.GFN2Energy(*host);

    LBFGSParam<double> param;
    param.epsilon = 1e-5;
    param.max_iterations = 20;

    LBFGSSolver<double> solver(param);
    LBFGSInterface fun(3 * host->AtomCount());
    fun.setMolecule(host);
    double fx;

    RMSDDriver* driver = new RMSDDriver;

    driver->setSilent(true);
    driver->setProtons(true);
    driver->setForceReorder(false);
    driver->setCheckConnections(false);

    int atoms_count = host->AtomCount();
    for (int outer = 0; outer < maxouter; ++outer) {
        int niter = solver.minimize(fun, parameter, fx);
        parameter = fun.Parameter();

        for (int i = 0; i < atoms_count; ++i) {
            geometry(i, 0) = parameter(3 * i);
            geometry(i, 1) = parameter(3 * i + 1);
            geometry(i, 2) = parameter(3 * i + 2);
        }
        h.setGeometry(geometry);

        driver->setReference(tmp);
        driver->setTarget(h);
        driver->AutoPilot();
        if (printOutput) {
            end = std::chrono::system_clock::now();
            std::cout << outer << "\t" << std::setprecision(9) << fun.m_energy << "\t\t" << std::setprecision(5) << (fun.m_energy - final_energy) * 2625.5 << "\t\t" << std::setprecision(6) << driver->RMSD() << "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0 << std::endl;
            start = std::chrono::system_clock::now();
        }
        final_energy = fun.m_energy;
        tmp = h;
        if (writeXYZ) {
            h.setEnergy(final_energy);
            h.appendXYZFile("curcuma_optim.xyz");
        }
        if (((fun.m_energy - final_energy) * 2625.5 < dE && driver->RMSD() < dRMSD) || niter < param.max_iterations)
            break;
    }


    for (int i = 0; i < host->AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }
    h.setEnergy(final_energy);
    h.setGeometry(geometry);
    return h;
}

inline void OptimiseGeometryThreaded(const Molecule* host, std::string* result_string, Molecule* result_molecule, double dE = 0.75, double dRMSD = 0.01)
{
    stringstream ss;

    Geometry geometry = host->getGeometry();
    Molecule tmp(host);
    Molecule h(host);
    Vector parameter(3 * host->AtomCount());

    for (int i = 0; i < host->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }
    std::chrono::time_point<std::chrono::system_clock> init = std::chrono::system_clock::now(), start = std::chrono::system_clock::now(), end;

    ss << "Step\tCurrent Energy [Eh]\tEnergy Change\tRMSD Change\tt [s]" << std::endl;
    int maxouter = 100;
    XTBInterface interface;
    double final_energy = interface.GFN2Energy(*host);

    LBFGSParam<double> param;
    param.epsilon = 1e-5;
    param.max_iterations = 20;

    LBFGSSolver<double> solver(param);
    LBFGSInterface fun(3 * host->AtomCount());
    fun.setMolecule(host);
    double fx;

    RMSDDriver* driver = new RMSDDriver;

    driver->setSilent(true);
    driver->setProtons(true);
    driver->setForceReorder(false);
    driver->setCheckConnections(false);

    int atoms_count = host->AtomCount();
    for (int outer = 0; outer < maxouter; ++outer) {
        int niter = solver.minimize(fun, parameter, fx);
        parameter = fun.Parameter();

        for (int i = 0; i < atoms_count; ++i) {
            geometry(i, 0) = parameter(3 * i);
            geometry(i, 1) = parameter(3 * i + 1);
            geometry(i, 2) = parameter(3 * i + 2);
        }
        h.setGeometry(geometry);

        driver->setReference(tmp);
        driver->setTarget(h);
        driver->AutoPilot();
        end = std::chrono::system_clock::now();
        ss << outer << "\t" << std::setprecision(9) << fun.m_energy << "\t\t" << std::setprecision(5) << (fun.m_energy - final_energy) * 2625.5 << "\t\t" << std::setprecision(6) << driver->RMSD() << "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0 << std::endl;
        start = std::chrono::system_clock::now();

        final_energy = fun.m_energy;
        tmp = h;
        if (((fun.m_energy - final_energy) * 2625.5 < dE && driver->RMSD() < dRMSD) || niter < param.max_iterations)
            break;
    }

    for (int i = 0; i < host->AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }
    h.setEnergy(final_energy);
    h.setGeometry(geometry);
    ss << "Final"
       << "\t" << std::setprecision(9) << final_energy << "\t\t" << std::setprecision(5) << (fun.m_energy - final_energy) * 2625.5 << "\t\t" << std::setprecision(6) << driver->RMSD() << "\t" << std::chrono::duration_cast<std::chrono::milliseconds>(end - init).count() / 1000.0 << std::endl;

    result_string->append(ss.str());
    result_molecule->LoadMolecule(h);
}
