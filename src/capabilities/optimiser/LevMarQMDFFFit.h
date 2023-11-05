/*
 * <LevenbergMarquardt QMDFF FC Fit. >
 * Copyright (C) 2023 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "src/capabilities/hessian.h"

#include "src/core/elements.h"
#include "src/core/global.h"
#include "src/core/molecule.h"
#include "src/core/pseudoff.h"

#include "src/tools/geometry.h"

template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>

struct FCFunctor {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    inline FCFunctor(int inputs, int values)
        : m_inputs(inputs)
        , m_values(values)
    {
    }

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

struct MyFCFunctor : FCFunctor<double> {
    inline MyFCFunctor(int inputs, int values)
        : FCFunctor(inputs, values)
        , no_parameter(inputs)
        , no_points(values)
    {
    }
    inline ~MyFCFunctor() {}
    inline int operator()(const Eigen::VectorXd& fc, Eigen::VectorXd& fvec) const
    {
        // std::cout << fc.transpose() << std::endl;
        json hc = HessianJson;
        hc["method"] = "qmdff";
        hc["threads"] = m_threads;
        Hessian he2(hc, true);

        json bonds = m_parameter["bonds"];
        json angles = m_parameter["angles"];
        int index = 0;
        for (int i = 0; i < bonds.size(); ++i) {
            bonds[i]["kAB"] = fc(index);
            index++;
        }
        for (int i = 0; i < angles.size(); ++i) {
            angles[i]["kabc"] = fc(index);
            index++;
        }
        json parameter;
        parameter["bonds"] = bonds;
        parameter["angles"] = angles;
        he2.setMolecule(m_molecule);
        he2.setParameter(parameter);
        he2.start();
        Matrix hessian = he2.getHessian();

        index = 0;
        for (int i = 0; i < m_hessian.rows(); ++i) {
            for (int j = 0; j < m_hessian.cols(); ++j) {
                fvec(index++) = (m_hessian(i, j) - hessian(i, j)) * (m_hessian(i, j) - hessian(i, j)); // + PseudoFF::DistancePenalty(m_host->Atom(i), guest.Atom(j));
            }
        }

        return 0;
    }
    int no_parameter;
    int no_points;
    int m_threads = 1;

    int inputs() const { return no_parameter; }
    int values() const { return no_points; }
    Matrix m_hessian;
    json m_parameter;
    Molecule m_molecule;
};

struct MyFunctorNumericalDiff : Eigen::NumericalDiff<MyFCFunctor> {
};

inline Vector OptimiseFC(const Molecule& molecule, const Matrix& hessian, const Vector& fc, const json& parameters, int threads)
{

    // Vector parameter = PositionPair2Vector(anchor, rotation);
    Vector parameter = fc;
    MyFCFunctor functor(fc.size(), hessian.cols() * hessian.rows());
    functor.m_hessian = hessian;
    functor.m_molecule = molecule;
    functor.m_parameter = parameters;
    functor.m_threads = threads;

    Eigen::NumericalDiff<MyFCFunctor> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MyFCFunctor>> lm(numDiff);
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
    // std::cout << parameter.transpose() << std::endl;
    int MaxIter = 300;
    Vector old_param = parameter;
    // std::cout << parameter.transpose() << std::endl;
    for (; iter < MaxIter /*&& ((qAbs(error_0 - error_2) > ErrorConvergence) || norm > DeltaParameter)*/; ++iter) {
        std::cout << "Step " << iter << " of maximal " << MaxIter << std::endl;
        status = lm.minimizeOneStep(parameter);

        if ((old_param - parameter).norm() < 1e-8)
            break;
        std::cout << parameter.transpose() << std::endl;
        old_param = parameter;
    }
    // std::cout << parameter.transpose() << std::endl;
    // std::cout << "took " << iter << " optimisation steps." << std::endl;

    return old_param;
}
