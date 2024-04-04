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
    LBFGSDipoleInterface(int n_)
    {
    }
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        double dx = 5e-1;
        std::vector<double> scaling(x.size());
        for (int i = 0; i < scaling.size(); ++i) {
            scaling[i] = x(i);
        }

        auto dipole_vector = m_molecule->CalculateDipoleMoment(scaling);
        double dipole = sqrt(dipole_vector[0] * dipole_vector[0] + dipole_vector[1] * dipole_vector[1] + dipole_vector[2] * dipole_vector[2]) * 2.5418;
        fx = std::abs(m_dipole - dipole);
        // std::cout << dipole << " " << m_dipole << " "<< fx<< std::endl;
        for (int i = 0; i < scaling.size(); ++i) {
            // if(std::abs(fx) < std::abs(m_smaller))
            //     std::cout << scaling[i] << " ";
            scaling[i] += dx;
            dipole_vector = m_molecule->CalculateDipoleMoment(scaling);
            double dipolep = sqrt(dipole_vector[0] * dipole_vector[0] + dipole_vector[1] * dipole_vector[1] + dipole_vector[2] * dipole_vector[2]) * 2.5418;
            scaling[i] -= 2 * dx;
            dipole_vector = m_molecule->CalculateDipoleMoment(scaling);
            double dipolem = sqrt(dipole_vector[0] * dipole_vector[0] + dipole_vector[1] * dipole_vector[1] + dipole_vector[2] * dipole_vector[2]) * 2.5418;

            grad[i] = (std::abs(dipolep - dipole) - std::abs(dipolem - dipole)) / (2.0 * dx);
            scaling[i] += dx;
        }
        // if(std::abs(fx) < std::abs(m_smaller))
        //     std::cout << std::endl;
        m_smaller = std::abs(fx);
        m_parameter = x;
        return fx;
    }

    Vector Parameter() const { return m_parameter; }
    void setMolecule(const Molecule* molecule)
    {
        m_molecule = molecule;
    }

    double m_dipole;
    double m_smaller = 1;

private:
    int m_atoms = 0;
    Vector m_parameter;
    const Molecule* m_molecule;
};

inline std::pair<double, std::vector<double>> OptimiseScaling(const Molecule* molecule, double dipole, double initial, double threshold, int maxiter)
{
    Vector parameter(molecule->AtomCount());
    for (int i = 0; i < molecule->AtomCount(); ++i)
        parameter(i) = initial;

    Vector old_param = parameter;
    // std::cout << parameter.transpose() << std::endl;

    LBFGSParam<double> param;
    param.epsilon = 1e-5;
    LBFGSSolver<double> solver(param);
    LBFGSDipoleInterface fun(parameter.size());
    fun.m_dipole = dipole * 2.5418;
    fun.setMolecule(molecule);
    int iteration = 0;
    // std::cout << parameter.transpose() << std::endl;
    double fx;
    int converged = solver.InitializeSingleSteps(fun, parameter, fx);

    for (iteration = 1; iteration <= maxiter && fx > threshold; ++iteration) {
        solver.SingleStep(fun, parameter, fx);
        std::cout << "Iteration: " << iteration << " Difference: " << fx << std::endl;
    }
    std::vector<double> scaling(parameter.size());
    for (int i = 0; i < scaling.size(); ++i)
        scaling[i] = parameter(i);
    auto dipole_vector = molecule->CalculateDipoleMoment(scaling);
    dipole = sqrt(dipole_vector[0] * dipole_vector[0] + dipole_vector[1] * dipole_vector[1] + dipole_vector[2] * dipole_vector[2]);

    return std::pair<double, std::vector<double>>(dipole, scaling);
}
