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
 *
 * Distance-space metric (Claude Generated 2026, -nebmd.cv_metric distances): instead of the
 * best-fit RMSD, d_i is the root-mean-square difference of a set of P atom-pair coordinates
 *
 *     x_p(R) = |R_a - R_b|          (transform 'd')   or   x_p(R) = 1/|R_a - R_b|  ('inv')
 *     d_i^2  = 1/P sum_p (x_p(R) - x_p(C_i))^2        (C_i = CV frame i, Cartesian-stored)
 *
 * No Kabsch fit is needed: pair distances are invariant under rotation and translation, and
 * atoms that enter no pair (e.g. methyl H) are invisible to s and z. s, z and their
 * gradients follow from d_i with the same formulas as above; z is in A^2 ('d') or A^-2 ('inv').
 * Limit: distances are also invariant under reflection, so a path between mirror images
 * (enantiomerisation) is not resolved in distance space.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "src/capabilities/rmsd.h"
#include "src/core/global.h"
#include "src/core/molecule.h"

class PathCV {
public:
    /*! \brief Metric of d_i (Claude Generated 2026): best-fit RMSD or atom-pair distances. */
    enum class Metric { Rmsd, Distances };
    /*! \brief Pair coordinate of the distance metric: x = d or x = 1/d (Claude Generated 2026). */
    enum class PairTransform { Distance, Inverse };
    /*! \brief Atom pair, 0-based, first < second. */
    using AtomPair = std::pair<int, int>;

    /*! \brief Build the CV from a frozen reference path (>= 2 frames).
     *  \param atoms 0-based indices of the atoms the RMSD runs over (empty = all atoms;
     *         ignored by the distance metric).
     *  \param metric Rmsd (default, the original code path) or Distances.
     *  \param pairs atom pairs of the distance metric (>= 1 required for Distances).
     *  \param transform x = d or x = 1/d for the distance metric.
     *  lambda <= 0: computeLambdaFromPath(path, atoms) resp. computeLambdaFromPairs(). */
    PathCV(const std::vector<Geometry>& path, const Molecule& templ, double lambda = -1.0,
           const std::vector<int>& atoms = {}, Metric metric = Metric::Rmsd,
           const std::vector<AtomPair>& pairs = {}, PairTransform transform = PairTransform::Distance);

    /*! \brief Progress along the path, s in [0,1]. */
    double s(const Geometry& R) const;
    /*! \brief Distance from the path (units of d_i^2: Angstrom^2; A^-2 for the 'inv' pair metric). */
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

    /* --- Distance-space metric helpers (Claude Generated 2026) --- */

    /*! \brief Pair coordinates x_p(R) = |R_a - R_b| (Distance) or 1/|R_a - R_b| (Inverse). */
    static Vector pairValues(const Geometry& R, const std::vector<AtomPair>& pairs, PairTransform t);
    /*! \brief Squared distance-space spacing 1/P sum_p (x_p(A) - x_p(B))^2. */
    static double pairSpacing2(const Geometry& A, const Geometry& B, const std::vector<AtomPair>& pairs,
                               PairTransform t);
    /*! \brief lambda = 1 / <d_{k,k+1}^2> over consecutive frames in the distance metric. */
    static double computeLambdaFromPairs(const std::vector<Geometry>& path, const std::vector<AtomPair>& pairs,
                                         PairTransform t);
    /*! \brief Parse an explicit pair list "1-4,1-9,4-9" (1-based; "i-j" is the PAIR (i,j), not a
     *  range). Errors: malformed token, i == j, index outside 1..natoms, duplicate pair (also
     *  reversed), fewer than 2 pairs. Result normalised to first < second, in input order. */
    static bool parsePairList(const std::string& spec, int natoms, std::vector<AtomPair>& pairs,
                              std::string& error);
    /*! \brief 'auto' pair rule: all pairs of `candidates` (0-based atoms) whose distance spans
     *  (max - min over `frames`) more than `threshold` Angstrom. spans (optional) receives the
     *  span of every selected pair. No minimum count is enforced here. */
    static std::vector<AtomPair> selectPairsBySpan(const std::vector<Geometry>& frames,
                                                   const std::vector<int>& candidates, double threshold,
                                                   std::vector<double>* spans = nullptr);
    /*! \brief Re-sample a Cartesian polyline onto nout frames of equal DISTANCE-SPACE arclength
     *  (endpoints kept, linear Cartesian interpolation inside a segment). The arclength of a
     *  segment is integrated over `sub` linear sub-steps because x_p is not linear along a
     *  straight Cartesian segment. */
    static std::vector<Geometry> resampleByPairs(const std::vector<Geometry>& path, int nout,
                                                 const std::vector<AtomPair>& pairs, PairTransform t,
                                                 int sub = 16);

    /*! \brief Index of the nearest CV frame (argmin_i d_i) and optionally z, without gradients. */
    int nearestFrame(const Geometry& R, double* z_out = nullptr) const;

    double lambda() const { return m_lambda; }
    int frames() const { return static_cast<int>(m_path.size()); }
    const std::vector<int>& atoms() const { return m_atoms; }
    Metric metric() const { return m_metric; }
    const std::vector<AtomPair>& pairs() const { return m_pairs; }
    PairTransform transform() const { return m_transform; }
    /*! \brief Unit of z: "A^2" (rmsd, transform d) or "A^-2" (transform inv). */
    std::string zUnit() const
    {
        return (m_metric == Metric::Distances && m_transform == PairTransform::Inverse) ? "A^-2" : "A^2";
    }

private:
    /*! \brief Best-fit RMSD of R to reference frame i, optionally with its gradient. */
    double rmsdTo(int i, const Geometry& R, Geometry* grad) const;
    /*! \brief Distance-space d_i of R to reference frame i; grad (same sign convention as
     *  rmsdTo: MINUS d_i dd_i/dR) = -(1/P) sum_p (x_p - x_p,i) dx_p/dR. Claude Generated 2026. */
    double distTo(int i, const Geometry& R, const Vector& x, Geometry* grad) const;
    /*! \brief d_i of every frame (and optionally its gradient), dispatching on the metric. */
    void allDistances(const Geometry& R, std::vector<double>& d, std::vector<Geometry>* gd) const;
    /*! \brief Rows of g that belong to the selection (g itself for an empty selection). */
    Geometry selectRows(const Geometry& g) const;

    std::vector<Geometry> m_path;
    std::vector<Geometry> m_path_sel; ///< reference frames restricted to the selected atoms
    std::vector<int> m_atoms;         ///< selected atoms, 0-based (empty = all)
    Molecule m_template;   ///< supplies the atom labels for the RMSD driver
    double m_lambda = 1.0;
    Metric m_metric = Metric::Rmsd;                       ///< Claude Generated 2026
    std::vector<AtomPair> m_pairs;                        ///< pairs of the distance metric
    PairTransform m_transform = PairTransform::Distance;
    std::vector<Vector> m_path_x;                         ///< x_p of every CV frame (distance metric)
};
