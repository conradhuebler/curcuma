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
#include <unsupported/Eigen/NonLinearOptimization>

#include <iostream>

#include "src/core/global.h"
#include "src/core/molecule.h"


template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct TFunctor {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    TFunctor(int inputs, int values)
        : m_inputs(inputs)
        , m_values(values)
    {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

struct OptDipoleFunctor : TFunctor<double> {
    OptDipoleFunctor(int inputs, int values)
        : TFunctor(inputs, values)
        , no_parameter(inputs)
        , no_points(values) {
    }
    ~OptDipoleFunctor() = default;
    int operator()(const Vector& scaling, Eigen::VectorXd& fvec) const {
        // calculation of residuals
        for (int i = 0; i < m_conformers.size(); ++i) {
            auto conf = m_conformers.at(i);
            fvec(i) = (conf.getDipole() - conf.CalculateDipoleMoment(scaling, m_bond)).norm();
        }
        return 0;
    }
    int no_parameter;
    int no_points;
    std::vector<Molecule> m_conformers;
    bool m_bond;

    int inputs() const { return no_parameter; }
    int values() const { return no_points; }
};

struct OptDipoleFunctorNumericalDiff : Eigen::NumericalDiff<OptDipoleFunctor> {};

inline Vector OptimiseDipoleScaling(const std::vector<Molecule>& conformers, Vector scaling, const bool bond = false) {

    OptDipoleFunctor functor(2, conformers.size());
    functor.m_conformers = conformers;
    functor.m_bond = bond;
    Eigen::NumericalDiff numDiff(functor);
    Eigen::LevenbergMarquardt lm(numDiff);

    /*
    lm.parameters.factor = config["LevMar_Factor"].toInt(); //step bound for the diagonal shift, is this related to damping parameter, lambda?
    lm.parameters.maxfev = config["LevMar_MaxFEv"].toDouble(); //max number of function evaluations
    lm.parameters.xtol = config["LevMar_Xtol"].toDouble(); //tolerance for the norm of the solution vector
    lm.parameters.ftol = config["LevMar_Ftol"].toDouble(); //tolerance for the norm of the vector function
    lm.parameters.gtol = config["LevMar_Gtol"].toDouble(); // tolerance for the norm of the gradient of the error vector
    lm.parameters.epsfcn = config["LevMar_epsfcn"].toDouble(); //error precision
    */

    Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(scaling);

    constexpr int MaxIter = 3000;
    Vector old_param = scaling;

    for (int iter = 0; iter < MaxIter; ++iter) {
        status = lm.minimizeOneStep(scaling);

        if ((old_param - scaling).norm() < 1e-6)
            break;

        old_param = scaling;
    }

    return scaling;
}

inline Matrix DipoleScalingCalculation(const std::vector<Molecule>& conformers)
{
    std::vector<Position> y; //Dipole moments
    std::vector<Geometry> F; // Geometry multiplied with partial Charge
    const auto para_size = conformers[0].AtomCount();
    const auto conformer_size = conformers.size();
    Matrix FTF = Matrix::Zero(para_size, para_size);
    Matrix FTy = Matrix::Zero(para_size, 1);
    for (int i = 0; i < conformer_size; ++i) {
        y.push_back(conformers[i].getDipole()); // in eA
        F.push_back(conformers[i].ChargeDistribution());
    }
    for (auto k : F) {
        for (int i = 0; i < para_size; ++i) {
            for (int j = 0; j < para_size; ++j) {
                FTF(i, j) += k(i, 0) * k(j, 0) + k(i, 1) * k(j, 1) + k(i, 2) * k(j, 2);
            }
        }
    }

    for (int i = 0; i < conformer_size; ++i) {
        for (int j = 0; j < para_size; ++j) {
            FTy(j) += y[i](0) * F[i](j, 0) + y[i](1) * F[i](j, 1) + y[i](2) * F[i](j, 2);
        }
    }

    Matrix Theta = FTF.colPivHouseholderQr().solve(FTy);

    return Theta;
}