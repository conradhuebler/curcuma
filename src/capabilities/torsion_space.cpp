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

#include "torsion_space.h"

#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/energy_calculators/ff_methods/gfnff_geometry.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <fmt/core.h>
#include <numeric>

namespace curcuma::TorsionSpace {

namespace {

    constexpr double rad2deg = 180.0 / M_PI;
    constexpr double deg2rad = M_PI / 180.0;

    /// Wrap an angle in degrees into (-180, 180].
    double wrap180(double a)
    {
        while (a > 180.0)
            a -= 360.0;
        while (a <= -180.0)
            a += 360.0;
        return a;
    }

    /// Neighbour lists from the 0/1 bond matrix of the molecule.
    std::vector<std::vector<int>> neighbourList(const curcuma::Molecule& mol)
    {
        const int n = mol.AtomCount();
        const Matrix bonds = mol.DistanceMatrix().second;
        std::vector<std::vector<int>> nb(n);
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                if (a != b && bonds(a, b) > 0.5)
                    nb[a].push_back(b);
        return nb;
    }

    /**
     * Breadth-first search from `start` over the neighbour list, refusing to traverse the single
     * bond (block_a, block_b). Returns the visited set (including `start`).
     */
    std::vector<int> componentWithoutBond(const std::vector<std::vector<int>>& nb, int start,
        int block_a, int block_b)
    {
        std::vector<char> seen(nb.size(), 0);
        std::vector<int> out;
        std::deque<int> queue{ start };
        seen[start] = 1;
        while (!queue.empty()) {
            const int cur = queue.front();
            queue.pop_front();
            out.push_back(cur);
            for (int next : nb[cur]) {
                if ((cur == block_a && next == block_b) || (cur == block_b && next == block_a))
                    continue; // the cut bond
                if (!seen[next]) {
                    seen[next] = 1;
                    queue.push_back(next);
                }
            }
        }
        return out;
    }

    /// Reference neighbour for the torsion definition: prefer a heavy atom, then the lowest index.
    int referenceNeighbour(const std::vector<int>& neighbours, int exclude, const curcuma::Molecule& mol)
    {
        int best = -1;
        bool best_heavy = false;
        for (int cand : neighbours) {
            if (cand == exclude)
                continue;
            const bool heavy = mol.Atom(cand).first != 1;
            if (best < 0 || (heavy && !best_heavy) || (heavy == best_heavy && cand < best)) {
                best = cand;
                best_heavy = heavy;
            }
        }
        return best;
    }

} // namespace

std::string Torsion::label(const curcuma::Molecule& mol) const
{
    return fmt::format("{}{}-{}{}",
        Elements::ElementAbbr[mol.Atom(j).first], j,
        Elements::ElementAbbr[mol.Atom(k).first], k);
}

std::vector<Torsion> detectTorsions(const curcuma::Molecule& mol)
{
    const int n = mol.AtomCount();
    std::vector<Torsion> torsions;
    if (n < 4)
        return torsions;

    const std::vector<std::vector<int>> nb = neighbourList(mol);

    // Heavy-neighbour count decides whether a rotation produces a new conformer at all.
    std::vector<int> heavy_degree(n, 0);
    for (int a = 0; a < n; ++a)
        for (int b : nb[a])
            if (mol.Atom(b).first != 1)
                heavy_degree[a]++;

    for (int a = 0; a < n; ++a) {
        if (mol.Atom(a).first == 1)
            continue;
        for (int b : nb[a]) {
            if (b <= a || mol.Atom(b).first == 1)
                continue;
            // (2) terminal group -> rotation only spins hydrogens / a single terminal atom
            if (heavy_degree[a] < 2 || heavy_degree[b] < 2)
                continue;
            // (3) ring bond -> cutting it leaves the two atoms connected
            std::vector<int> side_b = componentWithoutBond(nb, b, a, b);
            if (std::find(side_b.begin(), side_b.end(), a) != side_b.end())
                continue;

            Torsion t;
            // Normalise: the SMALLER side is the one that gets rotated. phi(i,j,k,l) is invariant
            // under reversal, so this is purely a choice of which half of the molecule moves.
            std::vector<int> side_a = componentWithoutBond(nb, a, a, b);
            if (side_b.size() <= side_a.size()) {
                t.j = a;
                t.k = b;
                t.moving = std::move(side_b);
            } else {
                t.j = b;
                t.k = a;
                t.moving = std::move(side_a);
            }
            t.i = referenceNeighbour(nb[t.j], t.k, mol);
            t.l = referenceNeighbour(nb[t.k], t.j, mol);
            if (t.i < 0 || t.l < 0)
                continue; // no reference atom -> dihedral undefined (should not happen after (2))

            // k lies on the rotation axis and must not be moved.
            t.moving.erase(std::remove(t.moving.begin(), t.moving.end(), t.k), t.moving.end());
            if (t.moving.empty())
                continue;
            for (int idx : t.moving)
                if (mol.Atom(idx).first != 1)
                    t.heavy_moving++;

            torsions.push_back(std::move(t));
        }
    }
    return torsions;
}

double dihedral(const Geometry& geometry, const Torsion& t)
{
    Matrix dummy;
    const double phi = GFNFF_Geometry::calculateDihedralAngle(
        Eigen::Vector3d(geometry.row(t.i)), Eigen::Vector3d(geometry.row(t.j)),
        Eigen::Vector3d(geometry.row(t.k)), Eigen::Vector3d(geometry.row(t.l)), dummy, false);
    return wrap180(phi * rad2deg);
}

std::vector<double> dihedrals(const Geometry& geometry, const std::vector<Torsion>& torsions)
{
    std::vector<double> out;
    out.reserve(torsions.size());
    for (const auto& t : torsions)
        out.push_back(dihedral(geometry, t));
    return out;
}

Geometry setDihedral(const Geometry& geometry, const Torsion& t, double target_deg)
{
    const double current = dihedral(geometry, t);
    const double delta = wrap180(target_deg - current);
    if (std::abs(delta) < 1e-9)
        return geometry;

    const Eigen::Vector3d origin = geometry.row(t.j);
    Eigen::Vector3d axis = Eigen::Vector3d(geometry.row(t.k)) - origin;
    const double axis_norm = axis.norm();
    if (axis_norm < 1e-10)
        return geometry; // degenerate bond, nothing sensible to do
    axis /= axis_norm;

    auto rotate = [&](double angle_deg) {
        Geometry out = geometry;
        const Eigen::AngleAxisd rot(angle_deg * deg2rad, axis);
        for (int idx : t.moving)
            out.row(idx) = (origin + rot * (Eigen::Vector3d(geometry.row(idx)) - origin)).transpose();
        return out;
    };

    // The sign relating a right-handed rotation about j->k to the dihedral definition is verified,
    // not derived: rotate, measure, and flip if the angle moved the wrong way. Costs one extra
    // dihedral evaluation and is immune to a convention mismatch in either routine.
    Geometry candidate = rotate(delta);
    if (std::abs(angleDifference(dihedral(candidate, t), target_deg)) > 1e-4)
        candidate = rotate(-delta);
    return candidate;
}

double angleDifference(double a_deg, double b_deg)
{
    return wrap180(a_deg - b_deg);
}

std::vector<double> clusterStates(const std::vector<double>& angles_deg, double tolerance_deg)
{
    if (angles_deg.empty())
        return {};

    std::vector<double> sorted = angles_deg;
    for (double& a : sorted)
        a = wrap180(a);
    std::sort(sorted.begin(), sorted.end());

    // Claude Generated (Aug 2026) -- TWO defects fixed here, both of which turned a well-sampled
    // torsion into "one state, carries no information":
    //
    // (a) WRAP MERGE. The old code split at every gap and then unconditionally merged the group at
    //     -180 with the one at +180 whenever they touched across the seam. For a torsion whose
    //     values cover the whole circle with a single interior gap, that split it there and glued it
    //     back together -- one "state" spanning 360 degrees, whose arithmetic mean is meaningless.
    //     Measured on a 142-structure peptide ensemble: torsion C1-C3, values from -179.6 to +179.8
    //     with a standard deviation of 108 degrees, was reported as a single state at -65.3 degrees.
    //     19 of 29 torsions were affected, i.e. exactly the BEST sampled ones were declared rigid.
    //     Fix: rotate the sorted sequence so it STARTS after the largest circular gap. The seam then
    //     lies where the data is sparsest and no wrap merge is needed at all.
    //
    // (b) CHAINING. Single linkage joins A-B and B-C into one group however far A and C are apart.
    //     A densely covered torsion therefore chains into one state even without (a). Fix: a state
    //     may not span more than twice the tolerance; wider groups are split recursively at their
    //     largest internal gap.
    const std::size_t n = sorted.size();
    std::size_t start = 0;
    double largest_gap = (sorted.front() + 360.0) - sorted.back(); // the seam gap
    for (std::size_t idx = 1; idx < n; ++idx) {
        const double gap = sorted[idx] - sorted[idx - 1];
        if (gap > largest_gap) {
            largest_gap = gap;
            start = idx;
        }
    }
    // Unrolled sequence beginning after the largest gap; values may exceed +180 and are wrapped back
    // when the centre is formed.
    std::vector<double> circle;
    circle.reserve(n);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t idx = (start + k) % n;
        double value = sorted[idx];
        if (idx < start)
            value += 360.0;
        circle.push_back(value);
    }

    std::vector<std::vector<double>> groups;
    groups.push_back({ circle.front() });
    for (std::size_t idx = 1; idx < circle.size(); ++idx) {
        if (circle[idx] - circle[idx - 1] <= tolerance_deg)
            groups.back().push_back(circle[idx]);
        else
            groups.push_back({ circle[idx] });
    }

    // Split over-wide groups at their largest internal gap until every state is narrow enough.
    const double max_span = 2.0 * tolerance_deg;
    for (std::size_t g = 0; g < groups.size(); ++g) {
        if (groups[g].size() < 2 || groups[g].back() - groups[g].front() <= max_span)
            continue;
        std::size_t cut = 1;
        double widest = -1.0;
        for (std::size_t idx = 1; idx < groups[g].size(); ++idx) {
            const double gap = groups[g][idx] - groups[g][idx - 1];
            if (gap > widest) {
                widest = gap;
                cut = idx;
            }
        }
        std::vector<double> tail(groups[g].begin() + cut, groups[g].end());
        groups[g].resize(cut);
        groups.insert(groups.begin() + g + 1, std::move(tail));
        --g; // re-examine the left part, it may still be too wide
    }

    std::vector<double> centres;
    centres.reserve(groups.size());
    for (const auto& g : groups)
        centres.push_back(wrap180(std::accumulate(g.begin(), g.end(), 0.0) / static_cast<double>(g.size())));
    std::sort(centres.begin(), centres.end());
    return centres;
}

int assignState(double angle_deg, const std::vector<double>& centres_deg)
{
    int best = -1;
    double best_dist = std::numeric_limits<double>::infinity();
    for (int s = 0; s < static_cast<int>(centres_deg.size()); ++s) {
        const double d = std::abs(angleDifference(angle_deg, centres_deg[s]));
        if (d < best_dist) {
            best_dist = d;
            best = s;
        }
    }
    return best;
}

std::vector<int> stateVector(const Geometry& geometry, const std::vector<Torsion>& torsions,
    const std::vector<std::vector<double>>& centres)
{
    std::vector<int> states(torsions.size(), -1);
    for (std::size_t idx = 0; idx < torsions.size(); ++idx) {
        if (idx >= centres.size())
            break;
        states[idx] = assignState(dihedral(geometry, torsions[idx]), centres[idx]);
    }
    return states;
}

int hammingDistance(const std::vector<int>& a, const std::vector<int>& b)
{
    if (a.size() != b.size())
        return -1;
    int d = 0;
    for (std::size_t idx = 0; idx < a.size(); ++idx)
        if (a[idx] != b[idx])
            d++;
    return d;
}

} // namespace curcuma::TorsionSpace
