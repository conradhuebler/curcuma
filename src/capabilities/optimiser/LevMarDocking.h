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

#include "src/tools/geometry.h"

template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>

struct Functor {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    inline Functor(int inputs, int values)
        : m_inputs(inputs)
        , m_values(values)
    {
    }

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

struct MyFunctor : Functor<double> {
    inline MyFunctor(int inputs, int values)
        : Functor(inputs, values)
        , no_parameter(inputs)
        , no_points(values)
    {
    }
    inline ~MyFunctor() {}
    inline int operator()(const Eigen::VectorXd& position, Eigen::VectorXd& fvec) const
    {
        Molecule guest = m_guest;
        guest.setGeometry(
            GeometryTools::TranslateAndRotate(
                guest.getGeometry(),
                guest.Centroid(),
                Position{ position(0), position(1), position(2) },
                Position{ position(3), position(4), position(5) }));

        for (int i = 0; i < m_host->AtomCount(); ++i) {
            fvec(i) = 0;
            for (int j = 0; j < guest.AtomCount(); ++j) {
                fvec(i) += PseudoFF::LJPotential(m_host->Atom(i), guest.Atom(j)); // + PseudoFF::DistancePenalty(m_host->Atom(i), guest.Atom(j));
            }
        }

        return 0;
    }
    int no_parameter;
    int no_points;
    const Molecule* m_host;

    Molecule m_guest;

    int inputs() const { return no_parameter; }
    int values() const { return no_points; }
};

struct MyFunctorNumericalDiff : Eigen::NumericalDiff<MyFunctor> {
};

inline std::pair<Position, Position> OptimiseAnchor(const Molecule* host, const Molecule& guest, Position anchor, Position rotation)
{
    Vector parameter = PositionPair2Vector(anchor, rotation);

    MyFunctor functor(6, host->AtomCount());
    functor.m_host = host;
    functor.m_guest = guest;
    Eigen::NumericalDiff<MyFunctor> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MyFunctor>> lm(numDiff);
    int iter = 0;

    /*
    lm.parameters.factor = config["LevMar_Factor"].toInt(); //step bound for the diagonal shift, is this related to damping parameter, lambda?
    lm.parameters.maxfev = config["LevMar_MaxFEv"].toDouble(); //max number of function evaluations
    lm.parameters.xtol = config["LevMar_Xtol"].toDouble(); //tolerance for the norm of the solution vector
    lm.parameters.ftol = config["LevMar_Ftol"].toDouble(); //tolerance for the norm of the vector function
    lm.parameters.gtol = config["LevMar_Gtol"].toDouble(); // tolerance for the norm of the gradient of the error vector
    lm.parameters.epsfcn = config["LevMar_epsfcn"].toDouble(); //error precision
    */

    Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(parameter);

    int MaxIter = 3000;
    Vector old_param = parameter;
    // std::cout << parameter.transpose() << std::endl;
    for (; iter < MaxIter /*&& ((qAbs(error_0 - error_2) > ErrorConvergence) || norm > DeltaParameter)*/; ++iter) {
        status = lm.minimizeOneStep(parameter);

        if ((old_param - parameter).norm() < 1e-5)
            break;

        old_param = parameter;
    }
    // std::cout << parameter.transpose() << std::endl;
    // std::cout << "took " << iter << " optimisation steps." << std::endl;

    return Vector2PositionPair(parameter);
}
