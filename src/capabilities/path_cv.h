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
 *
 * Atom subset (Claude Generated 2026, -nebmd.cv_atoms): optionally d_i is the RMSD over a
 * SUBSET A of the atoms only (e.g. the heavy atoms), with centroid and Kabsch rotation
 * taken over A as well:
 *
 *     d_i^2 = 1/|A| sum_{a in A} |Rot (p_a - <p>_A) - (q_ia - <q_i>_A)|^2
 *
 * Atoms outside A do not enter s or z, so their gradient is exactly zero; they stay
 * dynamic. Why: with all atoms counted equally, thermal motion of light peripheral atoms
 * (methyl H) adds to z, and a large displacement of one heavy atom is diluted by the
 * average over all atoms. An empty selection means all atoms (the original code path).
 */

#pragma once

#include <string>
#include <vector>

#include "src/capabilities/rmsd.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

class PathCV {
public:
    /*! \brief Build the CV from a frozen reference path (>= 2 frames).
     *  \param atoms 0-based indices of the atoms the RMSD runs over (empty = all atoms).
     *  lambda <= 0: computeLambdaFromPath(path, atoms). */
    PathCV(const std::vector<Geometry>& path, const Molecule& templ, double lambda = -1.0,
           const std::vector<int>& atoms = {});

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

    /*! \brief Flat-bottom wall on z ("tube potential"), Claude Generated 2026.
     *
     *      W(z) = 0                               z <= z_wall
     *      W(z) = k/2 (z - z_wall)^2              z >  z_wall
     *      dW/dR = k (z - z_wall) dz/dR           (0 inside the tube)
     *
     *  z in Angstrom^2, k in Eh/Angstrom^4, W in Eh, dW/dR in Eh/Angstrom (laboratory frame,
     *  like dz/dR). `grad` is always resized and set (zero inside the tube). Returns W. */
    double zWallEnergy(const Geometry& R, double z_wall, double k, Geometry* grad,
                       double* z_out = nullptr) const;

    /*! \brief Default wall stiffness, Claude Generated 2026: W = n_kT * kT at a penetration
     *  delta = 0.2 z_wall, i.e. k = 2 n_kT kT / (0.2 z_wall)^2 (n_kT = 3). */
    static double zWallDefaultK(double kT, double z_wall);

    /*! \brief Rule-of-thumb lambda from the mean squared spacing of the frames: plain (not
     *  best-fit) differences, mean square over `atoms` (empty = all atoms). */
    static double computeLambdaFromPath(const std::vector<Geometry>& path,
                                        const std::vector<int>& atoms = {});

    /*! \brief Parse a cv_atoms specification, Claude Generated 2026.
     *
     *  "all" -> empty vector (all atoms); "heavy" -> atoms with Z > 1; otherwise a list of
     *  1-based indices and ranges, e.g. "1-6,9". Result: sorted 0-based indices. Returns
     *  false with a message in `error` on an index out of range, a duplicate, a malformed
     *  token or fewer than 3 selected atoms (a CV on 2 atoms cannot see rotations about
     *  their axis, the Kabsch fit is then not unique). */
    static bool parseAtomSelection(const std::string& spec, const Molecule& mol,
                                   std::vector<int>& atoms, std::string& error);

    /*! \brief Per-atom weights of a selection, the form the string_path helpers take (Claude
     *  Generated 2026): 1 on the selected atoms, `other` on the rest (empty for an empty
     *  selection). other = 0: the selected-atom metric. other = 1e-6 (alignment only): the
     *  selected atoms decide the fit, the others only fix a rotation the selection leaves free
     *  (a collinear selection such as Cl-C-Cl cannot fix the rotation about its axis). */
    static std::vector<double> selectionWeights(const std::vector<int>& atoms, int natoms, double other = 0.0);

    double lambda() const { return m_lambda; }
    int frames() const { return static_cast<int>(m_path.size()); }
    const std::vector<int>& atoms() const { return m_atoms; }

private:
    /*! \brief Best-fit RMSD of R to reference frame i, optionally with its gradient. */
    double rmsdTo(int i, const Geometry& R, Geometry* grad) const;
    /*! \brief Rows of g that belong to the selection (g itself for an empty selection). */
    Geometry selectRows(const Geometry& g) const;

    std::vector<Geometry> m_path;
    std::vector<Geometry> m_path_sel; ///< reference frames restricted to the selected atoms
    std::vector<int> m_atoms;         ///< selected atoms, 0-based (empty = all)
    Molecule m_template;   ///< supplies the atom labels for the RMSD driver
    double m_lambda = 1.0;
};
