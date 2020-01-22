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

struct LevMarNEBPseudoFFBaseFunctor {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    inline LevMarNEBPseudoFFBaseFunctor(int inputs, int values)
        : m_inputs(inputs)
        , m_values(values)
    {
    }

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

struct LevMarNEBPseudoFF : LevMarNEBPseudoFFBaseFunctor<double> {
    inline LevMarNEBPseudoFF(int inputs, int values, int atoms)
        : LevMarNEBPseudoFFBaseFunctor(inputs, values)
        , no_parameter(inputs)
        , no_points(values)
    {
        m_first = Molecule(0, atoms);
        m_second = Molecule(0, atoms);
    }
    inline ~LevMarNEBPseudoFF() {}
    inline int operator()(const Eigen::VectorXd& position, Eigen::VectorXd& fvec) const
    {
        Molecule first = m_first;
        first.GetFragments();
        Molecule second = m_second;
        second.GetFragments();

        std::cout << position << std::endl;

        std::cout << std::endl;

        for (std::size_t i = 0; i < m_fragments_a; ++i) {
            Geometry geom = first.getGeometryByFragment(i + 1);
            first.setGeometryByFragment(GeometryTools::TranslateGeometry(
                                            geom,
                                            GeometryTools::Centroid(geom),
                                            Position{ position(3 * i), position(3 * i + 1), position(3 * i + 2) }),
                i + 1);
            std::cout << Position{ position(3 * i), position(3 * i + 1), position(3 * i + 2) }.transpose() << std::endl;
        }

        for (std::size_t i = 0; i < m_fragments_b; ++i) {
            Geometry geom = second.getGeometryByFragment(i + 1);
            second.setGeometryByFragment(GeometryTools::TranslateGeometry(
                                             geom,
                                             GeometryTools::Centroid(geom),
                                             Position{ position(3 * (i + m_fragments_a)), position(3 * (i + m_fragments_a) + 1), position(3 * (i + m_fragments_a) + 2) }),
                i + 1);
            std::cout << Position{ position(3 * (i + m_fragments_a)), position(3 * (i + m_fragments_a) + 1), position(3 * (i + m_fragments_a) + 2) }.transpose() << std::endl;
        }

        if (m_protons) {
            for (int i = 0; i < first.AtomCount(); ++i)
                fvec(i) = GeometryTools::Distance(first.Atom(i).second, second.Atom(i).second) * GeometryTools::Distance(first.Atom(i).second, second.Atom(i).second);
        } else {
            for (int i = 0; i < first.AtomCount(); ++i)
                if (m_first.Atom(i).first != 1)
                    fvec(i) = GeometryTools::Distance(first.Atom(i).second, second.Atom(i).second) * GeometryTools::Distance(first.Atom(i).second, second.Atom(i).second);
        }
        second.setName("neb_ende");
        second.appendXYZFile();

        return 0;
    }
    int no_parameter;
    int no_points;

    Molecule m_first;
    Molecule m_second;
    int m_fragments_a = 0, m_fragments_b = 0, m_fragments = 0;

    bool m_protons;
    int inputs() const { return no_parameter; }
    int values() const { return no_points; }
};

struct MyFunctorNumericalDiff : Eigen::NumericalDiff<LevMarNEBPseudoFF> {
};

std::pair<Molecule, Molecule> OptimiseAtoms(const Molecule& first, const Molecule& second)
{
    Molecule f = first;
    f.GetFragments();
    f.setName("neb_start");
    Molecule s = second;
    s.GetFragments();
    s.setName("neb_end");

    bool protons = false;
    int fragments_a = first.GetFragments().size() - 1;
    int fragments_b = second.GetFragments().size() - 1;
    int fragments = fragments_a + fragments_b;
    Vector parameter(fragments * 3);
    for (int i = 0; i < fragments * 3; ++i)
        parameter(i, 0) = 0.1;

    LevMarNEBPseudoFF functor(fragments * 3, first.AtomCount() * 3, first.AtomCount());
    functor.m_first.LoadMolecule(first);
    functor.m_second.LoadMolecule(second);
    functor.m_protons = protons;
    functor.m_fragments_a = fragments_a;
    functor.m_fragments_b = fragments_b;
    functor.m_fragments = fragments;

    Eigen::NumericalDiff<LevMarNEBPseudoFF> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LevMarNEBPseudoFF>> lm(numDiff);
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

    for (std::size_t i = 0; i < fragments_a; ++i) {
        Geometry geom = f.getGeometryByFragment(i + 1);
        f.setGeometryByFragment(GeometryTools::TranslateGeometry(
                                    geom,
                                    GeometryTools::Centroid(geom),
                                    Position{ parameter(3 * i), parameter(3 * i + 1), parameter(3 * i + 2) }),
            i + 1);
        std::cout << Position{ parameter(3 * i), parameter(3 * i + 1), parameter(3 * i + 2) }.transpose() << std::endl;
    }

    for (std::size_t i = 0; i < fragments_b; ++i) {
        Geometry geom = s.getGeometryByFragment(i + 1);
        s.setGeometryByFragment(GeometryTools::TranslateGeometry(
                                    geom,
                                    GeometryTools::Centroid(geom),
                                    Position{ parameter(3 * (i + fragments_a)), parameter(3 * (i + fragments_a) + 1), parameter(3 * (i + fragments_a) + 2) }),
            i + 1);
        std::cout << Position{ parameter(3 * (i + fragments_a)), parameter(3 * (i + fragments_a) + 1), parameter(3 * (i + fragments_a) + 2) }.transpose() << std::endl;
    }
    return std::pair<Molecule, Molecule>(f, s);
}
