/*
 * <LevenbergMarquardt Anchor Optimisation for Docking. >
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
#include "src/core/pseudoff.h"
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
        VectorXd error(n + 1);
        for (int i = 0; i < n; ++i)
            error[i] = -1;
        error[n] = 0;
        m_diis_last = error;

        m_driver = new RMSDDriver;
        m_driver->setSilent(true);
        m_driver->setProtons(true);
        m_driver->setForceReorder(false);
        m_driver->setCheckConnections(false);
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        VectorXd v = x;
        parameter = x;
        /*
        if(m_iter >= 2)
        {
            VectorXd error(n+1);
            for(int i = 0; i < n; ++i)
                error[i] = abs(x[i] - m_previous[i]);
            error[n] = -1;
            m_error.push_back(error);
            if(m_error.size() == 8)
                m_error.erase(m_error.begin());
            Eigen::MatrixXd matrix(m_error.size() + 1, n + 1);
            Eigen::VectorXd vector(m_error.size() + 1);
            for(int i = 0; i < m_error.size(); ++i)
            {
              matrix.row(i) = m_error[i];
              vector(i) = 0;
            }
            matrix.row(m_error.size()) = m_diis_last;
            vector(m_error.size()) = -1;
            std::cout << matrix << std::endl;
            std::cout << vector << std::endl;
            Eigen::VectorXd solved = matrix.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(vector);
            std::cout << solved.transpose() << std::endl;
            VectorXd zero = VectorXd::Zero(x.size());
            for(int i = 0; i < m_coords.size(); ++i)
                zero += m_coords[i] * solved(i);
            std::cout << zero.transpose() << std::endl << std::endl << v.transpose() << std::endl;
            v = zero;
        }
        m_coords.push_back(v);
        if(m_coords.size() == 8)
            m_coords.erase(m_coords.begin());
        */
        double fx = 0.0;

        Molecule host = m_host;

        Geometry geometry = host.getGeometry();
        int natoms = host.AtomCount();
        int attyp[host.AtomCount()];
        double charge = 0;
        std::vector<int> atoms = host.Atoms();
        double coord[3 * natoms];
        double gradient[3 * natoms];
        double dist_gradient[3 * natoms];

        for (int i = 0; i < m_host->AtomCount(); ++i) {
            //geometry(i, 0) = v(3 * i);
            //geometry(i, 1) = v(3 * i + 1);
            //geometry(i, 2) = v(3 * i + 2);
            attyp[i] = host.Atoms()[i];
            coord[3 * i + 0] = v(3 * i + 0) / au;
            coord[3 * i + 1] = v(3 * i + 1) / au;
            coord[3 * i + 2] = v(3 * i + 2) / au;
        }

        //host.setGeometry(geometry);
        /*
        m_driver->setReference(m_last);
        m_driver->setTarget(host);
        m_driver->AutoPilot();*/
        double Energy = interface->GFN2Energy(attyp, coord, natoms, charge, gradient);
        fx = Energy;
        for (int i = 0; i < host.AtomCount(); ++i) {
            grad[3 * i + 0] = gradient[3 * i + 0];
            grad[3 * i + 1] = gradient[3 * i + 1];
            grad[3 * i + 2] = gradient[3 * i + 2];
        }
        //host.setEnergy(fx);
        //host.appendXYZFile("move_host.xyz");
        //m_iter++;
        //m_last_change =  (Energy - m_energy)*2625.5;

        //std::cout << m_iter << "\t"<< std::setprecision(9) << fx << "\t\t" << std::setprecision(5) << (Energy - m_energy)*2625.5 << "\t\t" << std::setprecision(5)<< m_driver->RMSD() <<std::endl;
        m_energy = Energy;
        //m_last_rmsd = m_driver->RMSD();
        //m_previous = v;
        m_last = host;
        return fx;
    }

    const Molecule* m_host;
    double LastEnergy() const { return m_energy; }
    std::vector<VectorXd> m_error;
    std::vector<VectorXd> m_coords;

    VectorXd m_previous, m_diis_last;
    Molecule m_last;
    double m_energy = 0, m_last_change = 0, m_last_rmsd = 0;
    Vector parameter;

private:
    int m_iter = 0;
    int n;
    XTBInterface* interface;
    RMSDDriver* m_driver;
};

Molecule OptimiseGeometry(const Molecule* host, bool writeXYZ = true, bool printOutput = true)
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

    std::cout << "Step\tCurrent Energy [Eh]\tEnergy Change\tRMSD Change" << std::endl;
    int maxouter = 100;
    double final_energy = 0;

    LBFGSParam<double> param;
    param.epsilon = 1e-5;
    param.max_iterations = 20;

    LBFGSSolver<double> solver(param);
    LBFGSInterface fun(3 * host->AtomCount());
    fun.m_host = host;
    fun.m_last = host;
    double fx;

    RMSDDriver* driver = new RMSDDriver;

    driver->setSilent(true);
    driver->setProtons(true);
    driver->setForceReorder(false);
    driver->setCheckConnections(false);

    int atoms_count = host->AtomCount();
    for (int outer = 0; outer < maxouter; ++outer) {
        int niter = solver.minimize(fun, parameter, fx);
        parameter = fun.parameter;

        for (int i = 0; i < atoms_count; ++i) {
            geometry(i, 0) = parameter(3 * i);
            geometry(i, 1) = parameter(3 * i + 1);
            geometry(i, 2) = parameter(3 * i + 2);
        }
        h.setGeometry(geometry);

        driver->setReference(tmp);
        driver->setTarget(h);
        driver->AutoPilot();

        if (printOutput)
            std::cout << outer << "\t" << std::setprecision(9) << fun.m_energy << "\t\t" << std::setprecision(5) << (fun.m_energy - final_energy) * 2625.5 << "\t\t" << std::setprecision(5) << driver->RMSD() << std::endl;

        final_energy = fun.m_energy;
        tmp = h;
        if (writeXYZ) {
            h.setEnergy(final_energy);
            h.appendXYZFile("curcuma_optim.xyz");
        }
        if (((fun.m_energy - final_energy) * 2625.5 < 0.75 && driver->RMSD() < 0.01) || niter < param.max_iterations)
            break;
    }

    /*
    LBFGSParam<double> param;
    param.epsilon = 1e-5;
    param.max_iterations = 2000;

    LBFGSSolver<double> solver(param);
    LBFGSInterface fun(3 * host->AtomCount());
    fun.m_host = host;
    fun.m_last = host;
    double fx;
    int niter = solver.minimize(fun, parameter, fx);
    */
    for (int i = 0; i < host->AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }
    h.setEnergy(final_energy);
    h.setGeometry(geometry);
    return h;
}
