/*
 * <Torsion space: rotatable bonds, rotamer states, dihedral manipulation>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (Jul 2026)
 */

#pragma once

#include <string>
#include <vector>

#include "src/core/global.h"
#include "src/core/molecule.h"

/**
 * @brief The conformational search space of a molecule, reduced to its soft coordinates.
 *
 * WHY THIS EXISTS (the dimensionality argument):
 * A molecule has 3N-6 internal degrees of freedom, but conformers differ only in the SOFT ones --
 * the torsions about rotatable single bonds. Bond lengths, bond angles and ring puckers are stiff:
 * distort them and the geometry optimisation puts them back. So the space that actually distinguishes
 * conformers is
 *
 *     3N-6  (100-250 for a drug-sized molecule)   ->   n_rot  (typically 5-15)
 *
 * and because the values each torsion takes in a real ensemble cluster into a few groups (sp3-sp3:
 * gauche+/gauche-/anti, amide: cis/trans, ...), the continuous torus reduces once more to a DISCRETE
 * state vector s = (s_1 ... s_n_rot). That is the representation every conformer generator uses
 * (RDKit ETKDG, OpenBabel confab, CREST's Z-matrix crossing) and it is the same construction as a
 * side-chain rotamer library in protein modelling.
 *
 * This module provides exactly that reduction and its inverse:
 *   - detectTorsions()   : topology -> list of rotatable torsions (+ which atoms a rotation moves)
 *   - dihedral()         : geometry -> angle of one torsion
 *   - clusterStates()    : observed angles of one torsion -> its rotamer state centres
 *   - stateVector()      : geometry -> discrete state vector
 *   - setDihedral()      : state vector -> geometry (the inverse map, used to build proposals)
 *
 * Everything here works on the BOND TOPOLOGY only (Molecule::DistanceMatrix()), no force field and no
 * energy evaluation, so it is cheap and independent of the method.
 */
namespace curcuma::TorsionSpace {

/**
 * @brief One rotatable torsion i-j-k-l; the rotatable bond is the central pair j-k.
 *
 * `moving` holds the atoms that a rotation about j-k displaces: the connected component of k after
 * the bond j-k is cut, minus k itself (k sits on the rotation axis). The torsion is normalised so
 * that this is the SMALLER of the two sides -- the dihedral value is invariant under reversal
 * (phi(i,j,k,l) == phi(l,k,j,i)), so this only decides which half of the molecule is moved.
 */
struct Torsion {
    int i = -1, j = -1, k = -1, l = -1;
    std::vector<int> moving;    ///< atoms rotated when this torsion changes (k-side, without k)
    int heavy_moving = 0;       ///< non-hydrogen atoms in `moving` (size of the moved fragment)

    /// Human-readable label, e.g. "C3-C4" (the rotatable bond), for reports and CSV headers.
    std::string label(const curcuma::Molecule& mol) const;
};

/**
 * @brief Find all rotatable torsions from the bond topology.
 *
 * A bond j-k is rotatable when
 *   (1) neither atom is hydrogen,
 *   (2) both atoms have at least two heavy neighbours -- otherwise the rotation only spins a
 *       terminal group (methyl, -OH, -NH2) and produces no new conformer,
 *   (3) the bond is not part of a ring: cutting it must disconnect the molecule. This is tested
 *       exactly by a breadth-first search (no ring perception data needed), which also yields the
 *       `moving` set for free.
 *
 * Multiple bonds are NOT filtered out here (no bond orders in the plain topology): a C=C or amide
 * bond is reported as rotatable, and the ensemble statistics will simply show it in one or two
 * states. That is deliberate -- an amide cis/trans flip IS a conformational degree of freedom.
 *
 * @param mol molecule with a valid geometry (only its connectivity is used)
 * @return torsions ordered by the index pair of their rotatable bond (deterministic)
 */
std::vector<Torsion> detectTorsions(const curcuma::Molecule& mol);

/// Dihedral angle of `t` in degrees, in (-180, 180]. Geometry in Angstrom.
double dihedral(const Geometry& geometry, const Torsion& t);

/// All torsion angles of one structure, in degrees (same order as `torsions`).
std::vector<double> dihedrals(const Geometry& geometry, const std::vector<Torsion>& torsions);

/**
 * @brief Rotate `t.moving` about the j-k axis so that the dihedral becomes `target_deg`.
 * @return the modified geometry (Angstrom); the input is left untouched.
 *
 * The sign convention of the rotation relative to the dihedral definition is VERIFIED numerically
 * (one extra dihedral evaluation) rather than derived, so a convention mismatch cannot silently
 * produce mirror-image structures.
 */
Geometry setDihedral(const Geometry& geometry, const Torsion& t, double target_deg);

/// Shortest signed difference a-b of two angles in degrees, wrapped to (-180, 180].
double angleDifference(double a_deg, double b_deg);

/**
 * @brief Cluster the observed values of ONE torsion into rotamer state centres.
 *
 * Periodic single-linkage clustering: values whose circular gap is below `tolerance_deg` join the
 * same state; the centre is the circular mean of the group. Handles the wrap-around at +-180 (an
 * anti state scattered over 175 ... -175 degrees is ONE state, not two).
 *
 * @return state centres in degrees, sorted ascending; empty if `angles_deg` is empty
 */
std::vector<double> clusterStates(const std::vector<double>& angles_deg, double tolerance_deg = 40.0);

/// Index of the closest state centre (circular distance); -1 if `centres_deg` is empty.
int assignState(double angle_deg, const std::vector<double>& centres_deg);

/**
 * @brief Discrete state vector of one structure.
 * @param centres per-torsion state centres (outer index = torsion, as returned by clusterStates)
 */
std::vector<int> stateVector(const Geometry& geometry, const std::vector<Torsion>& torsions,
    const std::vector<std::vector<double>>& centres);

/// Number of positions in which two state vectors differ (Hamming distance); -1 on size mismatch.
int hammingDistance(const std::vector<int>& a, const std::vector<int>& b);

} // namespace curcuma::TorsionSpace
