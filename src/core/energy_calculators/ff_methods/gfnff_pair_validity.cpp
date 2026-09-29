/*
 * < GFN-FF: rev-gfnff pair-validity gate >
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
 * Claude Generated (Sep 29, 2026) - AI-generated, machine-tested only; human production
 * testing pending. Design: test_cases/revgfnff/_log/FABLE_BOND_STATE_2.md section 2.1-rev
 * (the formal VALID(i,j,b) rule) and BOND_VALIDITY_GATE_SWEEP_STATUS.md (the independently
 * verified acceptance numbers this C++ port must reproduce). Implementation status/acceptance
 * table: test_cases/revgfnff/_log/PAIR_VALIDITY_IMPL_STATUS.md.
 *
 * -gfnff.rev_pair_validity
 * ------------------------
 * A perceived bonded pair (i, j) is not always a bond: GFN-FF's discrete getnb() criterion
 * sometimes lists a closed-shell repulsion instead (six F...F contacts alongside the four B-F
 * bonds of a compressed BF4-; a geminal H...H "bond" between two already-saturated hydrogens of
 * the same carbon in a hot reactive trajectory). Both cases share one signature: neither end has
 * a free valence slot, no lone pair can bridge them, and no acceptor/electron-count exception
 * applies - so the pair is energetically a bond in the force field (it gets a well, it
 * re-parametrises hybridisation/angles/torsions around it) but chemically nothing. VALID(i,j,b)
 * below is that admission test, evaluated once per topology "corner" b from integers and the
 * corner's own per-atom budget cap (no geometry, no new global constant):
 *
 *     N(i)          listed partners of i in corner b;  deg(i) = |N(i)|
 *     metal(i)      GFN-FF metal_type(Z_i) > 0                       (topo.is_metal[i])
 *     cap_i         the conserving-share budget cap of this corner (FFWorkspace::shareCapForAtom,
 *                   the SAME formula prepareConservingShare() uses); Val_i = Val_Z(i) + cap_i
 *     deficient(i)  cap_i >= 0.5  or  metal(i)
 *     purebridge(k;i,j)  k in N(i) ∩ N(j)  and every partner of k other than i, j is H
 *     bridge(i,j)   deficient(i) and deficient(j) and #{k : purebridge(k;i,j)} >= 2
 *     n_other(i;j)  #{m in N(i) \ {j} : not bridge(i,m)}
 *     free(i;j)     metal(i)  or  n_other(i;j) < Val_i
 *     lp(i;j)       Z_i not in {H, C}  and  ve(Z_i) - n_other(i;j) >= 2   (main-group valence
 *                   electrons, GFNFFParameters::periodic_group; undefined for d/f-block -> false)
 *     acc(k)        metal(k)  or  deg(k) < Val_k
 *     qloc(i,j)     sum of Phase-1 topology_charges over {i,j} u N(i) u N(j) u
 *                   {H partners of any of those}
 *
 *     VALID(i,j,b)  iff  free(i;j) or free(j;i)
 *                    or  (Z_i = H and lp(j;i)) or (Z_j = H and lp(i;j))
 *                    or  bridge(i,j)
 *                    or  exists k in N(i) ∩ N(j) with acc(k)
 *                    or  qloc(i,j) >= +0.5
 *
 * An invalid pair is not merely zeroed in the bond energy: the corner is REGENERATED from the
 * bond list with that pair removed (generateGFNFFParameterSet() below), so hybridisation, rings,
 * pi-systems, angles, torsions and the Phase-1 EEQ are all freshly derived exactly as if the pair
 * had never been perceived - "aliased to the corner with the invalid pair removed"
 * (FABLE_BOND_STATE_2.md section 2). This is what lets the compressed-BF4- probe reproduce its
 * own 4-bond topology's energy bit-for-bit instead of merely cancelling one term.
 *
 * DEFAULT OFF. When on and every pair is already valid, findInvalidPairValidityPairs() returns
 * empty and generateGFNFFParameterSet() takes exactly the one call it always took - no forcing,
 * no second pass, bit-identical to the flag being off.
 */

#include "gfnff.h"
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"
#include "src/core/curcuma_logger.h"

#include <algorithm>
#include <unordered_set>

#include <fmt/format.h>

namespace {

/// Main-group valence electron count from GFNFFParameters::periodic_group (1-8; H=1, C=4, N=5,
/// O=6, halogens=7). Returns -1 (undefined) for a d/f-block element or an out-of-range Z - the
/// gate's lp() clause then reads as false there, the same "no rescue, but nothing decided by it
/// either" reading BOND_VALIDITY_GATE_SWEEP_STATUS.md section 4 measured for the reference sets.
inline int mainGroupValenceElectrons(int Z)
{
    if (Z < 1 || Z > 86)
        return -1;
    const int g = GFNFFParameters::periodic_group[Z - 1];
    return (g >= 1 && g <= 8) ? g : -1;
}

} // namespace

std::vector<std::pair<int, int>> GFNFF::findInvalidPairValidityPairs(const TopologyInfo& topo) const
{
    std::vector<std::pair<int, int>> invalid;
    const int N = m_atomcount;
    if (N == 0 || static_cast<int>(topo.neighbor_lists.size()) != N
        || static_cast<int>(topo.is_metal.size()) != N || topo.topology_charges.size() != N
        || static_cast<int>(m_atoms.size()) != N)
        return invalid; // inconsistent topology snapshot - do not gate rather than gate wrongly

    const auto& nb = topo.neighbor_lists;
    const auto& is_metal = topo.is_metal;

    // --- per-atom: qgroup, donor, cap, Val, deficient -------------------------------------
    // Verbatim the same three quantities prepareConservingShare() (ff_workspace_gfnff.cpp)
    // builds for the "conserving" valence share, evaluated here directly from the corner's own
    // topology (no FFWorkspace instance exists yet for a corner under test).
    std::vector<double> qgroup(N, 0.0);
    for (int i = 0; i < N; ++i) {
        qgroup[i] = topo.topology_charges(i);
        for (int j : nb[i])
            if (j >= 0 && j < N && m_atoms[j] == 1)
                qgroup[i] += topo.topology_charges(j);
    }
    std::vector<double> valz(N, 0.0);
    for (int i = 0; i < N; ++i)
        valz[i] = revValence(m_atoms[i]);
    std::vector<char> donor(N, 0);
    if (m_rev_share_donor_rule) {
        auto acceptor = [&](int j) {
            const int Zj = m_atoms[j];
            if (Zj >= 1 && Zj <= 86 && GFNFFParameters::periodic_group[Zj - 1] == 3)
                return true; // group 13 (curcuma main-group numbering): an empty p orbital
            return static_cast<double>(nb[j].size()) < valz[j] - 0.5; // a free coordination site
        };
        for (int i = 0; i < N; ++i)
            for (int j : nb[i])
                if (j >= 0 && j < N && acceptor(j)) {
                    donor[i] = 1;
                    break;
                }
    }
    std::vector<double> cap(N, 0.0), Val(N, 0.0);
    std::vector<char> deficient(N, 0);
    for (int i = 0; i < N; ++i) {
        bool delivered_growth = false;
        cap[i] = FFWorkspace::shareCapForAtom(m_atoms[i], qgroup[i], donor[i] != 0, valz[i],
                                               m_rev_budget_fix_h, delivered_growth);
        Val[i] = valz[i] + cap[i]; // meaningful only when not metal(i); every VALID() clause that
                                   // reads Val_i is OR'd with metal(i), so a delivered-growth
                                   // metal's Val is never actually dereferenced below
        deficient[i] = (cap[i] >= 0.5) || is_metal[i];
    }

    // --- bridge(a,b) for every listed edge (a,b), canonical a < b --------------------------
    auto canon = [](int a, int b) { return a < b ? std::make_pair(a, b) : std::make_pair(b, a); };
    struct PairHash {
        size_t operator()(const std::pair<int, int>& p) const noexcept
        {
            return (static_cast<size_t>(p.first) << 20) ^ static_cast<size_t>(p.second);
        }
    };
    std::unordered_set<std::pair<int, int>, PairHash> edges, bridgeEdges;
    for (int i = 0; i < N; ++i)
        for (int j : nb[i])
            if (j > i)
                edges.insert({ i, j });

    for (const auto& e : edges) {
        const int a = e.first, b = e.second;
        if (!deficient[a] || !deficient[b])
            continue; // bridge() needs both ends deficient; skip the purebridge scan otherwise
        int pure_count = 0;
        // shared neighbours of a and b (excluding a, b themselves)
        std::unordered_set<int> nb_b(nb[b].begin(), nb[b].end());
        for (int k : nb[a]) {
            if (k == a || k == b || nb_b.find(k) == nb_b.end())
                continue;
            bool pure = true;
            for (int p : nb[k])
                if (p != a && p != b && (p < 0 || p >= N || m_atoms[p] != 1)) {
                    pure = false;
                    break;
                }
            if (pure)
                ++pure_count;
        }
        if (pure_count >= 2)
            bridgeEdges.insert(e);
    }
    auto isBridge = [&](int a, int b) { return bridgeEdges.find(canon(a, b)) != bridgeEdges.end(); };

    // --- n_other(i;j) for every ordered (i,j) with j in N(i) -------------------------------
    // Stored as n_other[i][slot], slot matching nb[i]'s own order (avoids a map of pairs).
    std::vector<std::vector<int>> n_other(N);
    for (int i = 0; i < N; ++i) {
        n_other[i].resize(nb[i].size(), 0);
        for (size_t s = 0; s < nb[i].size(); ++s) {
            const int j = nb[i][s];
            int count = 0;
            for (int m : nb[i])
                if (m != j && !isBridge(i, m))
                    ++count;
            n_other[i][s] = count;
        }
    }
    auto nOtherOf = [&](int i, int j) -> int {
        for (size_t s = 0; s < nb[i].size(); ++s)
            if (nb[i][s] == j)
                return n_other[i][s];
        return static_cast<int>(nb[i].size()); // j not actually listed - conservative (no rescue)
    };

    // --- acc(k) = metal(k) or deg(k) < Val_k -----------------------------------------------
    std::vector<char> acc(N, 0);
    for (int k = 0; k < N; ++k)
        acc[k] = is_metal[k] || (static_cast<double>(nb[k].size()) < Val[k]);

    // --- VALID(i,j) for every listed edge ---------------------------------------------------
    for (const auto& e : edges) {
        const int i = e.first, j = e.second;
        const int ni = nOtherOf(i, j), nj = nOtherOf(j, i);
        const bool free_ij = is_metal[i] || (static_cast<double>(ni) < Val[i]);
        const bool free_ji = is_metal[j] || (static_cast<double>(nj) < Val[j]);
        bool valid = free_ij || free_ji;
        if (!valid) {
            const int ve_j = mainGroupValenceElectrons(m_atoms[j]);
            const int ve_i = mainGroupValenceElectrons(m_atoms[i]);
            const bool lp_j = (m_atoms[j] != 1 && m_atoms[j] != 6) && ve_j >= 0 && (ve_j - nj >= 2);
            const bool lp_i = (m_atoms[i] != 1 && m_atoms[i] != 6) && ve_i >= 0 && (ve_i - ni >= 2);
            valid = (m_atoms[i] == 1 && lp_j) || (m_atoms[j] == 1 && lp_i);
        }
        if (!valid)
            valid = isBridge(i, j);
        if (!valid) {
            // exists k in N(i) ∩ N(j) with acc(k)
            std::unordered_set<int> nb_j(nb[j].begin(), nb[j].end());
            for (int k : nb[i])
                if (k != i && k != j && nb_j.find(k) != nb_j.end() && k >= 0 && k < N && acc[k]) {
                    valid = true;
                    break;
                }
        }
        if (!valid) {
            // qloc(i,j) = sum of topology_charges over {i,j} u N(i) u N(j) u {H partners of those}
            std::vector<char> in_shell(N, 0);
            std::vector<int> shell;
            auto addAtom = [&](int a) {
                if (a >= 0 && a < N && !in_shell[a]) {
                    in_shell[a] = 1;
                    shell.push_back(a);
                }
            };
            addAtom(i);
            addAtom(j);
            for (int m : nb[i])
                addAtom(m);
            for (int m : nb[j])
                addAtom(m);
            const size_t base_count = shell.size();
            for (size_t s = 0; s < base_count; ++s)
                for (int h : nb[shell[s]])
                    if (h >= 0 && h < N && m_atoms[h] == 1)
                        addAtom(h);
            double qloc = 0.0;
            for (int a : shell)
                qloc += topo.topology_charges(a);
            valid = qloc >= 0.5;
        }
        if (!valid)
            invalid.emplace_back(i, j);
    }
    return invalid;
}

// Claude Generated (Sep 2026): the gate wrapper around generateGFNFFParameterSetImpl()
// (gfnff_method.cpp) - see gfnff.h for both declarations. Off by default (-gfnff.rev_pair_validity
// false), in which case this is exactly the single call generateGFNFFParameterSet() always made.
GFNFFParameterSet GFNFF::generateGFNFFParameterSet()
{
    if (!m_rev_pair_validity)
        return generateGFNFFParameterSetImpl();

    std::vector<std::pair<int, int>> invalid, pruned;
    int n_edges_before = 0;
    {
        const TopologyInfo& topo = getCachedTopology();
        invalid = findInvalidPairValidityPairs(topo);
        if (invalid.empty())
            return generateGFNFFParameterSetImpl(); // bit-identical: no forcing, no second pass
        const int N = static_cast<int>(topo.neighbor_lists.size());
        for (int i = 0; i < N; ++i)
            for (int j : topo.neighbor_lists[i])
                if (j > i) {
                    ++n_edges_before;
                    const auto pr = std::make_pair(i, j);
                    if (std::find(invalid.begin(), invalid.end(), pr) == invalid.end())
                        pruned.push_back(pr);
                }
        // `topo` is a reference into m_cached_topology, which the forced regeneration below will
        // reassign - nothing below this block may read it again.
    }

    const std::vector<std::pair<int, int>> saved_forced = m_forced_bonds;
    const bool saved_owns = m_react_owns_bonds;
    // Claude Generated (Sep 2026): whether m_forced_bonds/m_react_owns_bonds were ALREADY under
    // an explicit caller's control before this call. Found by testing on a live react-mode
    // trajectory (test_cases/revgfnff/_log/PAIR_VALIDITY_IMPL_STATUS.md): prepareTransitionCorners()
    // / rebuildReactiveTopology() set m_forced_bonds to a corner's own bond list BEFORE calling
    // generateGFNFFParameterSet(), then immediately read m_cached_topology again afterward
    // (captureCornerEEQ()) EXPECTING it to describe the corner that was just generated - if the
    // gate restores m_forced_bonds and invalidates the cache in between, that read silently
    // recomputes from the RESTORED (un-gated) bond list instead, so the corner's captured EEQ/
    // hybridisation disagree with the bonds actually installed into it (measured: the geminal
    // c2h6 H...H frame's blended energy came out IDENTICAL gate on/off, because the "new" corner's
    // EEQ snapshot silently reverted to the 8-bond state even though its bond list was the gated
    // 7-bond one). The static/default path is the opposite case: NOTHING re-manages
    // m_forced_bonds afterward, so leaving it at the pruned list would freeze every later
    // geometry (a subsequent -opt/-md step, or the next frame of a -batch run) to this one gated
    // decision forever - there restoring is required for correctness.
    // Resolution: restore+reset ONLY when the caller had NOT already taken explicit ownership of
    // the bond source (the static/default virgin state); otherwise leave the gated state in
    // place; the caller (prepareTransitionCorners()/rebuildReactiveTopology()) already re-sets
    // m_forced_bonds itself before its own next need, exactly as it does between every corner.
    const bool caller_owns_bonds = !saved_forced.empty() || saved_owns;
    auto forceAndReset = [&](std::vector<std::pair<int, int>> bonds, bool owns) {
        // Claude Generated (Sep 2026): the exact reset sequence prepareTransitionCorners() /
        // rebuildReactiveTopology() already use to force an explicit bond list through
        // getCachedBondList() - see GFNFF::getCachedBondList() ("Use externally provided bond
        // list exclusively"). owns=true so a corner that legitimately loses every pair (all of
        // its listed bonds invalid) is still authoritative rather than falling back to geometric
        // re-detection (same reasoning as react mode's m_react_owns_bonds).
        m_forced_bonds = std::move(bonds);
        m_react_owns_bonds = owns;
        m_cached_bond_list.reset();
        m_geometry_tracker.reset();
        m_static_topology_valid = false;
    };
    forceAndReset(pruned, true);

    GFNFFParameterSet gated;
    bool ok = true;
    try {
        gated = generateGFNFFParameterSetImpl();
    } catch (const std::exception& e) {
        ok = false;
        CurcumaLogger::warn(fmt::format(
            "rev-gfnff pair-validity gate: regeneration on the pruned bond list failed ({}), "
            "keeping the ungated corner", e.what()));
    }

    if (!caller_owns_bonds || !ok)
        forceAndReset(saved_forced, saved_owns);

    if (!ok)
        return generateGFNFFParameterSetImpl(); // fall back to the ORIGINAL (ungated) corner

    if (CurcumaLogger::get_verbosity() >= 2)
        CurcumaLogger::info(fmt::format(
            "rev-gfnff pair-validity gate: {} pair(s) invalid, corner topology regenerated "
            "({} -> {} bonds)",
            invalid.size(), n_edges_before, pruned.size()));
    return gated;
}
