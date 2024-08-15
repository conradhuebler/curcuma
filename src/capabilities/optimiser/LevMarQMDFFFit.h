/*
 * <LevenbergMarquardt QMDFF FC Fit. >
 * Copyright (C) 2023 - 2024 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "json.hpp"
using json = nlohmann::json;

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
        Hessian he2(m_controller, true);

        json bonds = m_parameter["bonds"];
        json angles = m_parameter["angles"];
        int index = 0;
        std::vector<std::pair<int, int>> pairs(bonds.size());
        for (int i = 0; i < bonds.size(); ++i) {
            bonds[i]["fc"] = fc(index);
            index++;
            pairs[i] = std::pair<int, int>(bonds[i]["i"], bonds[i]["j"]);
        }

        for (int i = 0; i < angles.size(); ++i) {
            angles[i]["fc"] = fc(index);
            index++;
        }
        json parameter = m_parameter;
        parameter["bonds"] = bonds;
        parameter["angles"] = angles;
        he2.setMolecule(m_molecule);
        he2.setParameter(parameter);
        he2.start();
        Matrix hessian = he2.getHessian();

        index = 0;
        double diff = 0;

                for (int i = 0; i < m_hessian.rows(); ++i) {
                    for (int j = i; j < m_hessian.cols(); ++j) {
                        if (index >= fvec.size())
                            break;
                        fvec(index) = (m_hessian(j, i) - (hessian(j, i) + m_const_hessian(j, i))) + (m_hessian(i, j) - (hessian(i, j) + m_const_hessian(i, j)));
                        diff += (m_hessian(i, j) - (hessian(i, j) + m_const_hessian(i, j)));
                        index++;
                    }
                }
                /*
                for (auto pair : pairs) {
                    int i = pair.first;
                    int j = pair.second;
                    if (index >= fvec.size()) {
                        std::cout << "mist" << index << " " << i << " " << j << std::endl;
                        break;
                    }
                    for (int c = 0; c < 3; ++c) {
                        for (int d = 0; d < 3; ++d) {
                            fvec(index) = m_hessian(3 * i + c, 3 * j + d) - (hessian(3 * i + c, 3 * j + d) + m_const_hessian(3 * i + c, 3 * j + d));
                            diff += m_hessian(3 * i + c, 3 * j + d) - (hessian(3 * i + c, 3 * j + d) + m_const_hessian(3 * i + c, 3 * j + d));
                            index++;
                        }
                    }
                }
                */
                // std::cout << index << " " << fvec.size()<< " " << diff << " ";
                return 0;
    }

    void Controller(const json& controller)
    {
        m_controller = MergeJson(HessianJson, m_controller);
        m_controller["method"] = "uff";
    }

    int no_parameter;
    int no_points;

    int inputs() const { return no_parameter; }
    int values() const { return no_points; }
    Matrix m_hessian, m_const_hessian;
    json m_parameter, m_controller;
    Molecule m_molecule;
};

struct MyFunctorNumericalDiff : Eigen::NumericalDiff<MyFCFunctor> {
};

inline Vector OptimiseFC(const Molecule& molecule, const Matrix& hessian, const Matrix& const_hessian, const Vector& fc, const json& parameters, const json& controller)
{
    Vector parameter = fc;
    MyFCFunctor functor(fc.size(), hessian.cols() * hessian.rows() / 2 + hessian.cols() / 2 /* 6*parameters["bonds"].size() */);
    // MyFCFunctor functor(fc.size(), 9 * parameters["bonds"].size());

    functor.m_hessian = hessian;
    functor.m_const_hessian = const_hessian;
    functor.m_molecule = molecule;
    functor.m_parameter = parameters;
    functor.m_controller = controller;

    std::cout << controller << std::endl;

    Eigen::NumericalDiff<MyFCFunctor> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MyFCFunctor>> lm(numDiff);
    int iter = 0;
    /*
        std::cout << "Target hessian " << std::endl
                  << hessian << std::endl;
        std::cout << "Const part of qmdff hessian " << std::endl
                  << const_hessian << std::endl;
    */
    /*
    lm.parameters.factor = config["LevMar_Factor"].toInt(); //step bound for the diagonal shift, is this related to damping parameter, lambda?
    lm.parameters.maxfev = config["LevMar_MaxFEv"].toDouble(); //max number of function evaluations
    lm.parameters.xtol = config["LevMar_Xtol"].toDouble(); //tolerance for the norm of the solution vector
    lm.parameters.ftol = config["LevMar_Ftol"].toDouble(); //tolerance for the norm of the vector function
    lm.parameters.gtol = config["LevMar_Gtol"].toDouble(); // tolerance for the norm of the gradient of the error vector
    lm.parameters.epsfcn = config["LevMar_epsfcn"].toDouble(); //error precision
    */

    Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(parameter);
    int MaxIter = 30;
    Vector old_param = parameter;
    for (; iter < MaxIter /*&& ((qAbs(error_0 - error_2) > ErrorConvergence) || norm > DeltaParameter)*/; ++iter) {
        std::cout << "Step " << iter << " of maximal " << MaxIter << std::endl;
        status = lm.minimizeOneStep(parameter);
        std::cout << "Norm " << (old_param - parameter).norm() << std::endl;
        if ((old_param - parameter).norm() < 1e-12)
            break;
        std::cout << parameter.transpose() << std::endl;
        old_param = parameter;
    }

    return old_param;
}
