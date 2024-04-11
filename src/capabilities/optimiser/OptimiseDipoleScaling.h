/*
 * <LevenbergMarquardt Optimisation for Dipole Calculation from Partial Charges. >
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
 *               2024 Gerd Gehrisch
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

#include "src/tools/geometry.h"

#include "json.hpp"
#include <LBFGS.h>
#include <LBFGSB.h>

using json = nlohmann::json;

using Eigen::VectorXd;
using namespace LBFGSpp;

class LBFGSDipoleInterface {
public:
    explicit LBFGSDipoleInterface(int n_)
    {
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx;
        double dx = m_dx;
        std::vector<double> scaling(x.size());
        for (int i = 0; i < scaling.size(); ++i) {
            scaling[i] = x(i);
        }

        auto dipole_vector = m_molecule->CalculateDipoleMoment(scaling);
        auto dipole = dipole_vector.norm();
        fx = std::abs((m_dipole.norm() - dipole_vector.norm()));
        // std::cout << dipole << " " << m_dipole << " "<< fx<< std::endl;
        for (int i = 0; i < scaling.size(); ++i) {
            // if(std::abs(fx) < std::abs(m_smaller))
                 std::cout << scaling[i] << " ";
            scaling[i] += dx;
            dipole_vector = m_molecule->CalculateDipoleMoment(scaling);
            double dipolep = dipole_vector.norm();
            scaling[i] -= 2 * dx;
            dipole_vector = m_molecule->CalculateDipoleMoment(scaling);
            double dipolem = dipole_vector.norm();

            grad[i] = (std::abs(dipolep - dipole) - std::abs(dipolem - dipole)) / (2.0 * dx);
            scaling[i] += dx;
        }
        // if(std::abs(fx) < std::abs(m_smaller))
             std::cout << std::endl;
        // m_smaller = std::abs(fx);
        m_parameter = x;
        return fx;
    }

    Vector Parameter() const { return m_parameter; }
    void setMolecule(const Molecule* molecule)
    {
        m_molecule = molecule;
    }

    Position m_dipole;
    double m_dx = 5e-1;

private:
    int m_atoms = 0;
    Vector m_parameter;
    const Molecule* m_molecule;
};

inline std::pair<Position, std::vector<double>> OptimiseScaling(const Molecule* molecule, Position dipole, double initial, double threshold, int maxiter, double dx)
{
    Vector parameter(molecule->AtomCount()); //dec parameter as array and init size of parameter to size of atoms
    for (int i = 0; i < molecule->AtomCount(); ++i)
        parameter(i) = initial; //init every parameter as initial

    //Vector old_param = parameter;
    // std::cout << parameter.transpose() << std::endl;

    LBFGSParam<double> param;
    param.epsilon = 1e-5;
    LBFGSSolver<double> solver(param);
    LBFGSDipoleInterface fun(parameter.size());
    fun.m_dipole = dipole;
    fun.m_dx = dx;
    fun.setMolecule(molecule);
    int iteration = 0;
    // std::cout << parameter.transpose() << std::endl;
    double fx;
    int converged = solver.InitializeSingleSteps(fun, parameter, fx);

    for (iteration = 1; iteration <= maxiter && fx > threshold; ++iteration) {
        try {
            solver.SingleStep(fun, parameter, fx);
            std::cout << "Iteration: " << iteration << " Difference: " << fx << std::endl;
        }
        catch (const std::logic_error& e) {
            std::cerr << e.what();
        }

    }
    std::vector<double> scaling(parameter.size());
    for (int i = 0; i < scaling.size(); ++i)
        scaling[i] = parameter(i);
    auto dipole_vector = molecule->CalculateDipoleMoment(scaling);
    dipole = dipole_vector;

    return std::pair<Position, std::vector<double>>(dipole, scaling);
}
