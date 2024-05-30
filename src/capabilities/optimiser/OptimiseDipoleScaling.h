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

#include "src/capabilities/optimiser/LevMarDocking.h"
#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/pseudoff.h"

#include "src/tools/geometry.h"

template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>

// Implement argmin x: F(x) = sum_i^N (y_i - f_i(x))**2
// f_i(x) = sum_j (x*q_j*r_j), N... number of confomere
// r_j=[ [xyz]_1, [xyz]_2, ..., [xyz]_m] m... number of atoms/parameter
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

    inline TFunctor(int inputs, int values)
        : m_inputs(inputs)
        , m_values(values)
    {
    }

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

struct OptDipoleFunctor : TFunctor<double> {
    inline OptDipoleFunctor(int inputs, int values)
        : TFunctor(inputs, values)
        , no_parameter(inputs)
        , no_points(values)
    {
    }
    inline ~OptDipoleFunctor() = default;
    inline int operator()(const Vector& scaling, Eigen::VectorXd& fvec) const
    {
        for (int i = 0; i < m_conformers.size(); ++i ){
            auto conf = m_conformers.at(i);
            fvec(i) = (conf.getDipole() - conf.CalculateDipoleMoment(scaling)).norm() ;

        }

        return 0;
    }
    int no_parameter;
    int no_points;
    std::vector<Molecule> m_conformers;

    int inputs() const { return no_parameter; }
    int values() const { return no_points; }
};

struct OptDipoleFunctorNumericalDiff : Eigen::NumericalDiff<OptDipoleFunctor> {
};

inline Vector OptimiseDipoleScaling(const std::vector<Molecule>& conformers, Vector scaling)
{
    OptDipoleFunctor functor(6,conformers.size());
    functor.m_conformers = conformers;
    Eigen::NumericalDiff<OptDipoleFunctor> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<OptDipoleFunctor>> lm(numDiff);


    /*
    lm.parameters.factor = config["LevMar_Factor"].toInt(); //step bound for the diagonal shift, is this related to damping parameter, lambda?
    lm.parameters.maxfev = config["LevMar_MaxFEv"].toDouble(); //max number of function evaluations
    lm.parameters.xtol = config["LevMar_Xtol"].toDouble(); //tolerance for the norm of the solution vector
    lm.parameters.ftol = config["LevMar_Ftol"].toDouble(); //tolerance for the norm of the vector function
    lm.parameters.gtol = config["LevMar_Gtol"].toDouble(); // tolerance for the norm of the gradient of the error vector
    lm.parameters.epsfcn = config["LevMar_epsfcn"].toDouble(); //error precision
    */

    Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(scaling);

    int MaxIter = 3000;
    Vector old_param = scaling;

    for (int iter = 0; iter < MaxIter; ++iter) {
        status = lm.minimizeOneStep(scaling);

        if ((old_param - scaling).norm() < 1e-5)
            break;

        old_param = scaling;
    }

    return scaling;
}


inline Matrix DipoleScalingCalculation(const std::vector<Molecule>& conformers){
    std::vector<Position> y;
    std::vector<Geometry> F;
    auto para_size = conformers[0].AtomCount();
    auto conformer_size = conformers.size();
    Matrix FTF(para_size, para_size);
    Vector FTy(para_size);
    for (const auto & conformer : conformers) {
        y.push_back(conformer.getDipole());//TODO Einheit überprüfen
        F.push_back(conformer.ChargeDistribution());
    }

    for (int i = 0; i < para_size; ++i){
        for (int j = 0; j < para_size; ++j){
            for (auto & k : F)
                FTF(i, j) += k(i, 0) * k(j, 0) + k(i, 1) * k(j, 1) + k(i, 2) * k(j, 2);
        }
    }

    for (int j = 0; j < para_size; ++j) {
        for (int i = 0; i < conformer_size; ++i){
            FTy(j) += y[i](0) * F[i](j, 0) + y[i](1) * F[i](j, 1) + y[i](2) * F[i](j, 2);
        }
    }


    Matrix Theta(para_size,1);

    Theta = FTF.inverse() * FTy;
    

    //inv(F.t@F)@(F.t@y);
    return Theta;
}