/*
 * <Path collective variables (Branduardi s,z) for reaction-coordinate sampling.>
 * Copyright (C) 2019 - 2026 Conrad Huebler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated 2026 (AI-generated, machine-tested pending - not TESTED/APPROVED).
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
/*
 * Path collective variables - Branduardi, Gervasio & Parrinello,
 * J. Chem. Phys. 126, 054103 (2007).
 *
 * Given a reference path R_1 ... R_N (e.g. an optimised NEB band), two scalar
 * functions of the CURRENT geometry R are defined via the RMSD to each reference
 * frame:
 *
 *     s(R) = 1/(N-1) * sum_i (i-1) exp(-lambda d_i^2) / sum_i exp(-lambda d_i^2)
 *     z(R) = -1/lambda * ln sum_i exp(-lambda d_i^2)
 *
 * with d_i = RMSD(R, R_i). `s` runs from 0 (reactant) to 1 (product) and measures
 * progress ALONG the path; `z` measures the distance FROM the path.
 *
 * Why this matters here: unlike the arclength of a live NEB-MD band, `s` is an
 * explicit function of R alone - the reference frames are frozen, so the umbrella
 * window does not migrate while sampling. That removes the systematic error the
 * spring-free control experiment exposed (see docs/NEB_MD.md), and it makes
 * standard umbrella sampling + WHAM applicable.
 *
 * The construction is system-agnostic: it needs a path, not a hand-picked internal
 * coordinate, so it works wherever a NEB/IDPP path can be produced.
 *
 * lambda should be chosen so that the exponentials of NEIGHBOURING frames overlap;
 * the usual rule of thumb is lambda ~ 1 / <d_i,i+1>^2 (the mean squared spacing of
 * the reference frames). computeLambdaFromPath() implements that.
 */

#pragma once

#include <vector>

#include "src/capabilities/rmsd.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

class PathCV {
public:
    /*! \brief Build the CV from a frozen reference path (>= 2 frames). */
    PathCV(const std::vector<Geometry>& path, const Molecule& templ, double lambda = -1.0);

    /*! \brief Progress along the path, s in [0,1]. */
    double s(const Geometry& R) const;
    /*! \brief Distance from the path (same units as RMSD^2 -> Angstrom^2). */
    double z(const Geometry& R) const;

    /*! \brief s, z and their cartesian gradients (natoms,3).
     *
     *  ds/dR = sum_i [ i/(N-1) - s ] * w_i * (-2 lambda d_i) * dd_i/dR
     *  with w_i = exp(-lambda d_i^2)/sum_j exp(-lambda d_j^2), i.e. the analytic
     *  derivative of the softmax-weighted index. dd_i/dR comes from
     *  RMSDDriver::Gradient() (best-fit RMSD derivative).
     *
     *  dz/dR = sum_i w_i * 2 d_i * dd_i/dR   with w_i the same softmax weights, since
     *  z = -1/lambda ln sum_i exp(-lambda d_i^2) and
     *  dz/dR = sum_i [exp(-lambda d_i^2)/sum_j ...] * 2 d_i dd_i/dR.
     *  Pass dz_dR = nullptr when the perpendicular gradient is not needed. */
    void evaluate(const Geometry& R, double* s_out, double* z_out, Geometry* ds_dR,
                  Geometry* dz_dR = nullptr) const;

    /*! \brief Rule-of-thumb lambda from the mean squared spacing of the frames. */
    static double computeLambdaFromPath(const std::vector<Geometry>& path);

    double lambda() const { return m_lambda; }
    int frames() const { return static_cast<int>(m_path.size()); }

private:
    /*! \brief Best-fit RMSD of R to reference frame i, optionally with its gradient. */
    double rmsdTo(int i, const Geometry& R, Geometry* grad) const;

    std::vector<Geometry> m_path;
    Molecule m_template;   ///< supplies the atom labels for the RMSD driver
    double m_lambda = 1.0;
};
