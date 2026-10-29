/*
 * <Backtracking line search of LBFGSpp with an optional cap on the displacement of one atom>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Derived from LBFGSpp/LineSearchBacktracking.h, Copyright (C) 2016-2023 Yixuan Qiu, MIT license.
 */
#pragma once

#include <Eigen/Core>
#include <LBFGSpp/Param.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Optimization {

/**
 * @brief Cap on the Cartesian displacement of any single atom in one L-BFGS step (Angstrom, 0 = off).
 *
 * Claude Generated (Oct 2026). thread_local: the driver sets it before every step in the thread that
 * runs the line search, so optimisations running in parallel threads do not share it.
 * Read by LineSearchBacktrackingCapped.
 */
struct LBFGSAtomStepCap {
    static double& value()
    {
        thread_local double v = 0.0;
        return v;
    }
};

/**
 * @brief LBFGSpp's LineSearchBacktracking with one addition: no trial point lies further than
 *        LBFGSAtomStepCap from the start point for any atom.
 *
 * Claude Generated (Oct 2026). Why: LBFGSpp starts from H0 = (s'y/y'y) I, learned from the soft modes
 * (water), and its line search (Strong Wolfe) looks at the total energy only. In a solvated system the
 * first steps overshoot stiff bonds (C=O, 1.21 -> 1.99 A) while the soft part gains much more energy;
 * the bounded GFN-FF bond term has no restoring force from about 1.9 A (docs/LBFGS_ATOM_STEP_CAP.md).
 * With the cap, the trial step length is limited to cap / max_i |d_i| (d_i the 3-vector of atom i),
 * including the expansion of the Wolfe branch. Cap = 0: the arithmetic is that of the original, step
 * for step. The loop is a copy of LineSearchBacktracking::LineSearch (state of external/LBFGSpp, Oct 2026).
 */
template <typename Scalar>
class LineSearchBacktrackingCapped {
private:
    using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

public:
    template <typename Foo>
    static void LineSearch(Foo& f, const LBFGSpp::LBFGSParam<Scalar>& param,
        const Vector& xp, const Vector& drt, const Scalar& /*step_max*/,
        Scalar& step, Scalar& fx, Vector& grad, Scalar& dg, Vector& x)
    {
        const Scalar dec = 0.5;
        const Scalar inc = 2.1;

        if (step <= Scalar(0))
            throw std::invalid_argument("'step' must be positive");

        // Largest step length that keeps every atom within the cap (infinite = no cap).
        Scalar cap_step = std::numeric_limits<Scalar>::infinity();
        const double cap = LBFGSAtomStepCap::value();
        if (cap > 0.0 && drt.size() % 3 == 0) {
            Scalar max_d2 = 0;
            for (Eigen::Index i = 0; i < drt.size(); i += 3)
                max_d2 = std::max(max_d2, drt.segment(i, 3).squaredNorm());
            if (max_d2 > 0)
                cap_step = Scalar(cap) / std::sqrt(max_d2);
        }

        const Scalar fx_init = fx;
        const Scalar dg_init = grad.dot(drt);
        if (dg_init > 0)
            throw std::logic_error("the moving direction increases the objective function value");

        const Scalar test_decr = param.ftol * dg_init;
        Scalar width;

        int iter;
        for (iter = 0; iter < param.max_linesearch; iter++) {
            // Capped: never evaluate beyond the cap.
            const bool at_cap = step >= cap_step;
            if (at_cap)
                step = cap_step;

            x.noalias() = xp + step * drt;
            fx = f(x, grad);

            if (fx > fx_init + step * test_decr || (fx != fx)) {
                width = dec;
            } else {
                dg = grad.dot(drt);

                if (param.linesearch == LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
                    break;

                if (dg < param.wolfe * dg_init) {
                    // The step wants to grow, but the cap is reached: accept the capped step.
                    if (at_cap)
                        break;
                    width = inc;
                } else {
                    if (param.linesearch == LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_WOLFE)
                        break;

                    if (dg > -param.wolfe * dg_init) {
                        width = dec;
                    } else {
                        break;
                    }
                }
            }

            if (step < param.min_step)
                throw std::runtime_error("the line search step became smaller than the minimum value allowed");

            if (step > param.max_step)
                throw std::runtime_error("the line search step became larger than the maximum value allowed");

            step *= width;
        }
    }
};

} // namespace Optimization
