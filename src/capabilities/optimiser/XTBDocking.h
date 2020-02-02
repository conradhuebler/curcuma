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

#include <iostream>

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/pseudoff.h"
#include "src/core/xtbinterface.h"

#include "src/tools/geometry.h"

/*
template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>

struct LevMarXTBDockingFunctor {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    inline LevMarXTBDockingFunctor(int inputs, int values)
        : m_inputs(inputs)
        , m_values(values)
    {
    }

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

struct LevMarXTBDockingFunction : LevMarXTBDockingFunctor<double> {
    inline LevMarXTBDockingFunction(int inputs, int values)
        : LevMarXTBDockingFunctor(inputs, values)
        , no_parameter(inputs)
        , no_points(values)
    {
        interface = new XTBInterface;
    }
    inline ~LevMarXTBDockingFunction() {delete interface;}
    inline int operator()(const Eigen::VectorXd& position, Eigen::VectorXd& fvec) const
    {
        Molecule host = m_host;
        Molecule guest = m_guest;

        Geometry geometry = host.getGeometry();
        int natoms = host.AtomCount();
        int attyp[host.AtomCount()];
        double charge = 0;
        std::vector<int> atoms = host.Atoms();
        double coord[3*natoms];

        for(int i = 0; i < m_host->AtomCount(); ++i)
        {
            geometry(i, 0) = position(3*i);
            geometry(i, 1) = position(3*i + 1);
            geometry(i, 2) = position(3*i + 2);
            attyp[i] = host.Atoms()[i];
            coord[3*i+0] = position(3*i+0)/au;
            coord[3*i+1] = position(3*i+1)/au;
            coord[3*i+2] = position(3*i+2)/au;
        }

        host.setGeometry(geometry);
        // host.appendXYZFile("move_host.xyz");

        double Energy = interface->GFN2Energy(attyp, coord, natoms, charge);

        double distance = 0;
        for (int i = 0; i < m_host->AtomCount(); ++i) {
            for (int j = 0; j < guest.AtomCount(); ++j) {
                distance += PseudoFF::LJPotential(m_host->Atom(i), guest.Atom(j)) + PseudoFF::DistancePenalty(m_host->Atom(i), guest.Atom(j));
            }
        }
        std::cout << Energy << " " << distance << std::endl;
        fvec(0) = Energy + distance;
        return 0;
    }
    int no_parameter;
    int no_points;
    const Molecule* m_host;
    const Molecule* m_guest;;
    XTBInterface *interface;
    int inputs() const { return no_parameter; }
    int values() const { return no_points; }
};

struct LevMarXTBDockingFunctionNumericalDiff : Eigen::NumericalDiff<LevMarXTBDockingFunction> {
};
*/

#include <external/LBFGSpp/include/LBFGS.h>

using Eigen::VectorXd;
using namespace LBFGSpp;

class XTBPotential {

public:
    XTBPotential(int n_)
        : n(n_)
    {
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {

        double fx = 0.0;

        Molecule host = m_host;
        Molecule guest = m_guest;

        Geometry geometry = host.getGeometry();
        int natoms = host.AtomCount();
        int attyp[host.AtomCount()];
        double charge = 0;
        std::vector<int> atoms = host.Atoms();
        double coord[3 * natoms];
        double gradient[3 * natoms];
        double dist_gradient[3 * natoms];

        for (int i = 0; i < m_host->AtomCount(); ++i) {
            geometry(i, 0) = x(3 * i);
            geometry(i, 1) = x(3 * i + 1);
            geometry(i, 2) = x(3 * i + 2);
            attyp[i] = host.Atoms()[i];
            coord[3 * i + 0] = x(3 * i + 0) / au;
            coord[3 * i + 1] = x(3 * i + 1) / au;
            coord[3 * i + 2] = x(3 * i + 2) / au;
        }

        host.setGeometry(geometry);

        double scal = 0.9999;
        double Energy = interface->GFN2Energy(attyp, coord, natoms, charge, gradient);

        double distance = 0;
        double max_x = 0, max_y = 0, max_z = 0;
        for (int i = 0; i < m_host->AtomCount(); ++i) {
            Position g{ 0, 0, 0 };
            for (int j = 0; j < guest.AtomCount(); ++j) {
                Vector4d indiv_grad = PseudoFF::DistancePenaltyDiff(host.Atom(i), guest.Atom(j));
                distance += indiv_grad(0);
                max_x = std::max(indiv_grad(0), max_x);
                max_y = std::max(indiv_grad(1), max_y);
                max_z = std::max(indiv_grad(2), max_z);

                g(0) += indiv_grad(1);
                g(1) += indiv_grad(2);
                g(2) += indiv_grad(3);
            }
            dist_gradient[3 * i + 0] = g(0);
            dist_gradient[3 * i + 1] = g(1);
            dist_gradient[3 * i + 2] = g(2);
        }
        std::cout << Energy << " " << distance << std::endl;
        fx = scal * Energy + (1 - scal) * distance;
        for (int i = 0; i < host.AtomCount(); ++i) {
            grad[3 * i + 0] = scal * gradient[3 * i + 0] + (1 - scal) * dist_gradient[3 * i + 0];
            grad[3 * i + 1] = scal * gradient[3 * i + 1] + (1 - scal) * dist_gradient[3 * i + 1];
            grad[3 * i + 2] = scal * gradient[3 * i + 2] + (1 - scal) * dist_gradient[3 * i + 2];
        }
        host.setEnergy(fx);
        host.appendXYZFile("move_host.xyz");

        Molecule result = host;
        for (int i = 0; i < guest.AtomCount(); ++i)
            result.addPair(guest.Atom(i));
        result.appendXYZFile("test_xxx.xyz");
        return fx;
    }

    const Molecule* m_host;
    const Molecule* m_guest;

private:
    int n;

    XTBInterface* interface;
};

Geometry PrepareHost(const Molecule* host, const Molecule* guest)
{
    Geometry geometry = host->getGeometry();
    Molecule h(host);
    Vector parameter(3 * host->AtomCount());

    for (int i = 0; i < host->AtomCount(); ++i) {
        parameter(3 * i) = geometry(i, 0);
        parameter(3 * i + 1) = geometry(i, 1);
        parameter(3 * i + 2) = geometry(i, 2);
    }

    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    LBFGSSolver<double> solver(param);
    XTBPotential fun(3 * host->AtomCount());
    fun.m_host = host;
    fun.m_guest = guest;
    double fx;
    int niter = solver.minimize(fun, parameter, fx);
    for (int i = 0; i < host->AtomCount(); ++i) {
        geometry(i, 0) = parameter(3 * i);
        geometry(i, 1) = parameter(3 * i + 1);
        geometry(i, 2) = parameter(3 * i + 2);
    }

    return geometry;
}
