/*
 * < GFN-FF: chemistry-aware, continuous fragment-charge placement >
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
 * Claude Generated (Sep 24, 2026) - AI-generated, machine-tested only; human production
 * testing pending. Design and validation: test_cases/revgfnff/_log/FRAG_CHARGE_STATUS.md.
 *
 * -gfnff.frag_charge_model ensemble
 * ---------------------------------
 * The reference GFN-FF (pprcht/xtb) puts the whole net charge of a molecule that it perceives as
 * several fragments on fragment 0 - the fragment of atom 1 - and pins each fragment's Phase-1 and
 * Phase-2 EEQ charge sum to that integer. That rule is (1) index dependent (H3O+(H2O)6 with a
 * water written first gets a charged WATER) and (2) discrete: at the pass-1 bond threshold the
 * charge of a symmetric X2- pair jumps from (-1/2,-1/2) to (-1,0), a 100-200 kcal/mol step.
 *
 * The model here keeps the reference's physics INSIDE each charge state and changes only which
 * states are used and how they are weighted:
 *
 *   variant v = (grouping of the base fragments into EEQ constraint groups, integer charge per
 *               group); E_v = the complete reference GFN-FF energy for that assignment, with all
 *               qa-dependent parameters consistent with it (a GFNFF instance with an override).
 *
 *   contact edge (f,g): s_ij = r_ij / r_thr_ij, r_thr = the pass-1 (qa = 0) getnb threshold, so
 *               s = 1 is exactly where pass 1 splits the pair. L(s) = smootherstep on [1, s_max],
 *               lambda_fg = prod_ij L(s_ij).  lambda = 0: in contact (merged), 1: separated.
 *
 *   corner c   = subset of the window edges that are merged (2^k, as in the rev-gfnff stage-1
 *               topology corners): W_c = prod_{e in c}(1 - lambda_e) prod_{e not in c} lambda_e.
 *
 *   placement p = integer charges on the corner's groups; E_c = sum_p omega_p E_p with
 *               omega = softmax(-E_p / tau) (smooth stand-in for the PPLB lowest-state limit).
 *
 *   E = sum_c W_c E_c,
 *   dE/dx = sum_c [ W_c sum_p omega_p (1 - (E_p - E_c)/tau) dE_p/dx ] + sum_e (sum_c E_c dW_c/dlambda_e) dlambda_e/dx
 *
 * Limits: all separated -> lowest-energy integer placement (reference physics, carrier chosen by
 * energy, not by index); symmetric -> exact average of mirror states (correct dissociation limit
 * with consistent parameters); at the pass-1 split lambda -> 0 -> only the merged corner, which is
 * the one-fragment charge model of the other side of the threshold (continuous).
 */

#include "gfnff.h"
#include "src/core/energy_calculators/ff_methods/gfnff_par.h"
#include "src/core/curcuma_logger.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>

#include <fmt/format.h>

namespace {

/// smootherstep L(s) on [1, s_max] and its derivative dL/ds (C2 at both ends)
inline void fragWindowL(double s, double s_max, double& L, double& dL)
{
    if (s <= 1.0) { L = 0.0; dL = 0.0; return; }
    if (s >= s_max) { L = 1.0; dL = 0.0; return; }
    const double w = s_max - 1.0;
    const double t = (s - 1.0) / w;
    L = t * t * t * (10.0 - 15.0 * t + 6.0 * t * t);
    dL = 30.0 * t * t * (1.0 - t) * (1.0 - t) / w;
}

/// Experimental atomic electron affinities in eV (NIST / Andersen, Haugen, Hotop, J. Phys. Chem.
/// Ref. Data 28 (1999) 1511; rounded), for -gfnff.frag_charge_atomic_ea (opt-in). 0 = anion not
/// bound (He, Be, N, Ne, Mg, Ar, Zn, Kr, Cd, Xe, Hg, Rn); NaN = not tabulated here (transition
/// metals, lanthanides, ...), which makes the rule fall back to the free-charge weighting.
/// Claude Generated (Sep 2026), test_cases/revgfnff/_log/I2_CLF_STATUS.md.
double atomicElectronAffinityEV(int z)
{
    switch (z) {
    case 1: return 0.754195;  case 2: return 0.0;       case 3: return 0.618049;  case 4: return 0.0;
    case 5: return 0.279723;  case 6: return 1.262119;  case 7: return 0.0;       case 8: return 1.461105;
    case 9: return 3.401190;  case 10: return 0.0;      case 11: return 0.547926; case 12: return 0.0;
    case 13: return 0.43283;  case 14: return 1.389521; case 15: return 0.746607; case 16: return 2.077103;
    case 17: return 3.612724; case 18: return 0.0;      case 19: return 0.501459; case 30: return 0.0;
    case 32: return 1.232712; case 33: return 0.8048;   case 34: return 2.020605; case 35: return 3.363588;
    case 36: return 0.0;      case 37: return 0.485916; case 48: return 0.0;      case 50: return 1.112066;
    case 51: return 1.047401; case 52: return 1.970875; case 53: return 3.059047; case 54: return 0.0;
    case 55: return 0.471626; case 80: return 0.0;      case 82: return 0.356743; case 83: return 0.942363;
    case 86: return 0.0;
    default: return std::numeric_limits<double>::quiet_NaN();
    }
}

/// scoped per-thread verbosity (the variants run silently)
struct FragQuiet {
    int saved;
    explicit FragQuiet(int level) : saved(CurcumaLogger::thread_verbosity()) { CurcumaLogger::set_thread_verbosity(level); }
    ~FragQuiet() { CurcumaLogger::set_thread_verbosity(saved); }
};

/// all compositions of `units` identical units over `m` slots
void compositions(int units, int m, std::vector<int>& cur, std::vector<std::vector<int>>& out)
{
    if (m == 0) {
        if (units == 0) out.push_back(cur);
        return;
    }
    if (m == 1) {
        cur.push_back(units);
        out.push_back(cur);
        cur.pop_back();
        return;
    }
    for (int k = units; k >= 0; --k) {
        cur.push_back(k);
        compositions(units - k, m - 1, cur, out);
        cur.pop_back();
    }
}

void addScaled(FFEnergyComponents& a, const FFEnergyComponents& b, double w)
{
    a.bond += w * b.bond; a.angle += w * b.angle; a.dihedral += w * b.dihedral; a.inversion += w * b.inversion;
    a.dispersion += w * b.dispersion; a.vdw += w * b.vdw; a.rep += w * b.rep;
    a.bonded_rep += w * b.bonded_rep; a.nonbonded_rep += w * b.nonbonded_rep;
    a.coulomb += w * b.coulomb; a.hbond += w * b.hbond; a.xbond += w * b.xbond;
    a.atm += w * b.atm; a.batm += w * b.batm; a.stors += w * b.stors; a.over_coord += w * b.over_coord;
    a.sqe_hardness += w * b.sqe_hardness;
    a.hbond_case1 += w * b.hbond_case1; a.hbond_case2 += w * b.hbond_case2;
    a.hbond_case3 += w * b.hbond_case3; a.hbond_case4 += w * b.hbond_case4;
}

long binom(int n, int k)
{
    if (k < 0 || k > n) return 0;
    long r = 1;
    for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

} // namespace

void GFNFF::setFragmentOverride(const std::vector<int>& fraglist, int nfrag, const std::vector<double>& qfrag)
{
    m_frag_override.active = true;
    m_frag_override.fraglist = fraglist;
    m_frag_override.nfrag = nfrag;
    m_frag_override.qfrag = qfrag;
    m_frag_is_variant = true;
}

double GFNFF::fragPass1Threshold(int i, int j) const
{
    // Identical to perceiveGeometricBonds() at qa = 0 (pass 1 of the q-loop, gfnff_ini2.f90:111-126
    // + getnb icase 1): threshold = fm_i fm_j rthr rab(Z_i,Z_j,normcn) fat_i fat_j.
    using GFNFFParameters::metal_type;
    using GFNFFParameters::normcn;
    constexpr double rthr = 1.25;
    constexpr double rthr2 = 1.00;
    const int zi = m_atoms[i], zj = m_atoms[j];
    const double ncn_i = (zi >= 1 && zi <= 86) ? static_cast<double>(normcn[zi - 1]) : 4.0;
    const double ncn_j = (zj >= 1 && zj <= 86) ? static_cast<double>(normcn[zj - 1]) : 4.0;
    auto fm = [&](int z) {
        const int mt = (z >= 1 && z <= 86) ? metal_type[z - 1] : 0;
        return mt == 2 ? rthr2 : (mt == 1 ? rthr2 + 0.025 : 1.0);
    };
    double rco = GFNFFParameters::computeRabEstimate(zi, zj, ncn_i, ncn_j);
    // Fragments split by the charge-shrunk pass 2 (cations: rtmp -= qa*fq, gfnff_ini2.f90:122):
    // anchor the window at THAT threshold, with the pass-1 charges it was built from.
    if (m_frag_split_pass == 2 && static_cast<int>(m_bond_qa.size()) == m_atomcount) {
        constexpr double rqshrink = 0.23;
        auto fq = [&](int z) { const int mt = (z >= 1 && z <= 86) ? metal_type[z - 1] : 0; return rqshrink * (mt > 0 ? 2.0 : 1.0); };
        rco -= m_bond_qa[i] * fq(zi) + m_bond_qa[j] * fq(zj);
    }
    rco *= fat[zi] * fat[zj];
    double thr = fm(zi) * fm(zj) * rthr * rco;
    // rev_excess_bond_extend (Sep 27, 2026): a 2c-3e candidate pair splits at the EXTENDED
    // threshold, so its window starts there. Isolation is tested with the ordinary pass-1
    // criterion against every other atom (the same "no other bond" condition as the rule).
    if (m_rev_bond_extend > 1.0 && revX2PairExtendable(i, j))
        return m_rev_bond_extend * getnbThresholdPass1(i, j);   // charge-independent, as the rule
    return thr;
}

std::vector<GFNFF::FragEdge> GFNFF::fragWindowEdges(const std::vector<int>& fraglist, int nfrag) const
{
    std::map<std::pair<int, int>, FragEdge> edges;
    for (int i = 0; i < m_atomcount; ++i) {
        const int fi = fraglist[i] - 1;
        for (int j = i + 1; j < m_atomcount; ++j) {
            const int fj = fraglist[j] - 1;
            if (fi == fj || fi < 0 || fj < 0 || fi >= nfrag || fj >= nfrag)
                continue;
            const double r = (m_geometry_bohr.row(i) - m_geometry_bohr.row(j)).norm();
            const double thr = fragPass1Threshold(i, j);
            if (thr <= 0.0 || r >= m_frag_s_max * thr)
                continue;
            double L, dL;
            fragWindowL(r / thr, m_frag_s_max, L, dL);
            const auto key = std::make_pair(std::min(fi, fj), std::max(fi, fj));
            FragEdge& e = edges[key];
            e.f = key.first;
            e.g = key.second;
            e.contacts.push_back({ i, j, thr, L, dL });
        }
    }
    std::vector<FragEdge> out;
    for (auto& kv : edges) {
        FragEdge e = std::move(kv.second);
        e.lambda = 1.0;
        for (const auto& c : e.contacts)
            e.lambda *= c.L;
        out.push_back(std::move(e));
    }
    return out;
}

GFNFF* GFNFF::fragVariant(const std::string& key, const std::vector<int>& group_of_atom, int ngroups,
                          const std::vector<double>& qgroup)
{
    for (auto& v : m_frag_variants)
        if (v.key == key)
            return v.ff.get();

    json p = m_parameters;
    p.erase("gfnff");          // its keys are already promoted; it would re-enable the ensemble
    p.erase("geometry_file");  // no on-disk topology cache shared with the master
    p["cache_topology"] = false;
    p["frag_charge_model"] = "reference";
    p["print_timing"] = false;
    auto sub = std::make_unique<GFNFF>(p);
    sub->setThreadCount(m_threads);
    sub->setFragmentOverride(group_of_atom, ngroups, qgroup);

    Mol mol;
    mol.m_number_atoms = m_atomcount;
    mol.m_atoms = m_atoms;
    mol.m_geometry = m_geometry_bohr * BOHR_TO_ANGSTROM;
    mol.m_charge = m_charge;
    mol.m_spin = m_spin;
    mol.m_bonds = m_forced_bonds;
    mol.m_has_pbc = m_has_pbc;
    mol.m_unit_cell = m_unit_cell;
    bool ok = false;
    {
        FragQuiet quiet(0);
        ok = sub->InitialiseMolecule(mol) && !sub->eeqSolveFailed();
    }
    if (!ok) {
        CurcumaLogger::warn(fmt::format("GFNFF frag_charge_model ensemble: variant {} failed to initialise - skipped", key));
        return nullptr;
    }
    FragVariant v;
    v.key = key;
    v.ff = std::move(sub);
    m_frag_variants.push_back(std::move(v));
    return m_frag_variants.back().ff.get();
}

double GFNFF::fragEnsembleBlend(bool gradient, double e_master)
{
    const TopologyInfo& topo = getCachedTopology();
    if (static_cast<int>(topo.fraglist.size()) != m_atomcount || static_cast<int>(topo.nb_full.size()) != m_atomcount) {
        m_frag_variants.clear();
        return e_master;
    }
    if (m_topology_version != m_frag_variants_topo_version) {
        m_frag_variants.clear();
        m_frag_free_q = Vector();
        m_frag_variants_topo_version = m_topology_version;
    }

    // ---- base fragments at THIS geometry -------------------------------------------------
    // Components of the kept bond graph restricted to the pairs a fresh pass-1 perception would
    // still bond (s <= 1). With a fresh topology this is the reference fragmentation; with a
    // topology carried over from another geometry (MD, optimisation, batch reuse) a bond that has
    // meanwhile been stretched past the split is cut here, so the energy does not depend on where
    // the topology was built (package 25 section 4, topology-history dependence).
    std::vector<int> fraglist(m_atomcount, 0);
    int F = 0;
    {
        std::vector<int> parent(m_atomcount);
        std::iota(parent.begin(), parent.end(), 0);
        std::function<int(int)> root = [&](int a) { return parent[a] == a ? a : (parent[a] = root(parent[a])); };
        for (int i = 0; i < m_atomcount; ++i)
            for (int j : topo.nb_full[i]) {
                if (j <= i || j >= m_atomcount) continue;
                const double r = (m_geometry_bohr.row(i) - m_geometry_bohr.row(j)).norm();
                if (r <= fragPass1Threshold(i, j))
                    parent[root(i)] = root(j);
            }
        std::map<int, int> lab;
        for (int i = 0; i < m_atomcount; ++i) {
            auto it = lab.find(root(i));
            if (it == lab.end())
                it = lab.emplace(root(i), static_cast<int>(lab.size()) + 1).first;
            fraglist[i] = it->second;
        }
        F = static_cast<int>(lab.size());
    }
    // exactly the master's own charge model -> nothing to blend (bit-identical to the reference)
    if (F < 2 && topo.nfrag < 2)
        return e_master;

    // ---- contact edges inside the window -------------------------------------------------
    std::vector<FragEdge> edges = fragWindowEdges(fraglist, F);
    if (static_cast<int>(edges.size()) > m_frag_max_edges) {
        std::stable_sort(edges.begin(), edges.end(), [](const FragEdge& a, const FragEdge& b) { return a.lambda < b.lambda; });
        if (!m_frag_warned_edges) {
            CurcumaLogger::warn(fmt::format("GFNFF frag_charge_model ensemble: {} fragment contacts inside the window, blending the {} closest "
                                            "(frag_charge_max_edges); the rest are treated as separated",
                edges.size(), m_frag_max_edges));
            m_frag_warned_edges = true;
        }
        edges.resize(m_frag_max_edges);
    }
    const int k = static_cast<int>(edges.size());

    // ---- free single-constraint Phase-1 charges (placement pre-selection only) ------------
    auto freeCharges = [&]() -> const Vector& {
        if (m_frag_free_q.size() == m_atomcount)
            return m_frag_free_q;
        EEQSolver::TopologyInput ti;
        ti.neighbor_lists = topo.neighbor_lists;
        ti.nfrag = 1;
        ti.fraglist.assign(m_atomcount, 1);
        ti.qfrag = { static_cast<double>(m_charge) };
        ti.itag = topo.itag;
        ti.covalent_radii.resize(m_atomcount);
        for (int i = 0; i < m_atomcount; ++i) {
            const int z = m_atoms[i];
            ti.covalent_radii[i] = (z >= 1 && z <= static_cast<int>(GFNFFParameters::covalent_radii.size()))
                ? GFNFFParameters::covalent_radii[z - 1] : 1.0;
        }
        FragQuiet quiet(0);
        m_eeq_solver->invalidateCholeskyCache();
        m_frag_free_q = m_eeq_solver->calculateTopologyCharges(m_atoms, m_geometry_bohr, m_charge,
            topo.coordination_numbers, ti, true, threadPool(), m_threads);
        m_eeq_solver->invalidateCholeskyCache();
        m_eeq_solver->invalidateMatrixCache();
        if (m_frag_free_q.size() != m_atomcount)
            m_frag_free_q = Vector::Zero(m_atomcount);
        return m_frag_free_q;
    };

    const int nq = std::abs(m_charge);
    const double unit = (m_charge > 0) ? 1.0 : -1.0;

    // ---- enumerate corners and placements --------------------------------------------------
    // Placement rule (FRAG_CHARGE_STATUS.md section 2): GFN-FF's absolute energies of differently
    // charged fragments are NOT comparable (its "electron affinities" are -554 .. -648 kcal/mol and
    // ordered C6H6 > CH4 > Cl > H2O > F), so the carrier is chosen by chemistry, not by energy:
    //   1. parity: no group stripped of all its electrons if avoidable, then the fewest
    //      odd-electron (radical) groups - a closed-shell ion next to closed-shell neutrals always
    //      wins (H3O+ in water, formate...HF, Cl- + CH4);
    //   2. among those, classes of chemically identical carriers (same element multiset and charge)
    //      are weighted by the free single-constraint Phase-1 EEQ charge they carry (topology
    //      constant, index-free): Omega = softmax(S / sigma);
    //   3. within a class (e.g. the two ends of a symmetric X...X- pair) the members are weighted by
    //      their energies, softmax(-E/tau): there the fragment biases cancel and the energy only
    //      measures the environment.
    struct Place { int variant; int cls; double Omega; };
    std::vector<double> Wc;                                   // weight per active corner
    std::vector<unsigned> corner_mask;
    std::vector<std::vector<Place>> corner_uses;              // per active corner
    const Vector& qfree = freeCharges();
    // A one-group corner of a one-fragment master IS the master's own evaluation: reuse it
    // (exact, and saves one full GFN-FF evaluation per step).
    int master_idx = -1;
    if (topo.nfrag == 1) {
        for (int v = 0; v < static_cast<int>(m_frag_variants.size()); ++v)
            if (m_frag_variants[v].key == "master") master_idx = v;
        if (master_idx < 0) {
            FragVariant mv;
            mv.key = "master";
            m_frag_variants.push_back(std::move(mv));
            master_idx = static_cast<int>(m_frag_variants.size()) - 1;
        }
    }
    for (unsigned mask = 0; mask < (1u << k); ++mask) {
        double W = 1.0;
        for (int e = 0; e < k; ++e)
            W *= ((mask >> e) & 1u) ? (1.0 - edges[e].lambda) : edges[e].lambda;
        if (W <= 0.0)
            continue;   // exact: see FRAG_CHARGE_STATUS.md section 1.3 (its dW/dx vanishes too)
        // groups = connected components of the merged edges over the base fragments
        std::vector<int> parent(F);
        std::iota(parent.begin(), parent.end(), 0);
        std::function<int(int)> root = [&](int a) { return parent[a] == a ? a : (parent[a] = root(parent[a])); };
        for (int e = 0; e < k; ++e)
            if ((mask >> e) & 1u)
                parent[root(edges[e].f)] = root(edges[e].g);
        std::vector<int> gid(F, -1);
        std::map<int, int> rootg;
        for (int f = 0; f < F; ++f) {
            const int r = root(f);
            auto it = rootg.find(r);
            if (it == rootg.end())
                it = rootg.emplace(r, static_cast<int>(rootg.size())).first;
            gid[f] = it->second;
        }
        const int G = static_cast<int>(rootg.size());
        std::vector<int> group_of_atom(m_atomcount);
        std::vector<long> nelec(G, 0);                        // electrons of the NEUTRAL group
        std::vector<double> qfree_g(G, 0.0);
        std::vector<std::vector<int>> zlist(G);
        for (int i = 0; i < m_atomcount; ++i) {
            const int g = gid[fraglist[i] - 1];
            group_of_atom[i] = g + 1;
            nelec[g] += m_atoms[i];
            qfree_g[g] += qfree(i);
            zlist[g].push_back(m_atoms[i]);
        }
        for (auto& z : zlist) std::sort(z.begin(), z.end());

        // all compositions of the net charge over the groups, then the parity filter
        std::vector<std::vector<int>> comps;
        std::vector<int> cur;
        compositions(nq, G, cur, comps);
        // rank = bare_nuclei * (G + 1) + radicals: a placement that strips a group of ALL its
        // electrons (a bare proton next to CO in a pass-2-split CH2O2+) is never preferred to one
        // that does not, whatever the radical count; then the fewest odd-electron groups.
        struct Cand { std::vector<double> qg; int radicals; double S; std::string cls; double ea; };
        std::vector<Cand> cands;
        int rmin = std::numeric_limits<int>::max();
        for (const auto& cmp : comps) {
            Cand c;
            c.qg.assign(G, 0.0);
            c.radicals = 0;
            c.S = 0.0;
            c.ea = std::numeric_limits<double>::quiet_NaN();   // atomic EA of a single-atom -1 carrier
            std::vector<std::string> parts;
            for (int g = 0; g < G; ++g) {
                c.qg[g] = unit * cmp[g];
                const long ne = nelec[g] - static_cast<long>(c.qg[g]);
                if (ne & 1L) ++c.radicals;
                if (ne <= 0) c.radicals += G + 1;
                c.S += c.qg[g] * qfree_g[g] / static_cast<double>(nq);
                if (m_charge == -1 && cmp[g] == 1 && zlist[g].size() == 1)
                    c.ea = atomicElectronAffinityEV(zlist[g][0]);
                if (cmp[g] != 0) {
                    std::string zs = fmt::format("{:+d}:", static_cast<int>(c.qg[g]));
                    for (int z : zlist[g]) zs += fmt::format("{}.", z);
                    parts.push_back(zs);
                }
            }
            std::sort(parts.begin(), parts.end());
            for (const auto& ps : parts) c.cls += ps + "/";
            rmin = std::min(rmin, c.radicals);
            cands.push_back(std::move(c));
        }
        std::vector<Cand> kept;
        for (auto& c : cands)
            if (c.radicals == rmin) kept.push_back(std::move(c));
        if (static_cast<int>(kept.size()) > m_frag_max_placements) {
            std::stable_sort(kept.begin(), kept.end(), [](const Cand& x, const Cand& y) { return x.S > y.S; });
            kept.resize(m_frag_max_placements);
        }
        // class weights Omega = softmax(S_class / sigma), S_class = max over its members
        // opt-in frag_charge_atomic_ea (I2_CLF_STATUS.md): when EVERY kept candidate puts the -1 on
        // one single atom with a tabulated EA, rank by that EA instead (the asymptote gap of
        // A- + B vs A + B- is exactly EA(B) - EA(A)); otherwise the free-charge rule below.
        bool use_ea = m_frag_atomic_ea && !kept.empty();
        for (const auto& c : kept)
            if (!std::isfinite(c.ea)) use_ea = false;
        const double sig = use_ea ? m_frag_ea_sigma : m_frag_sigma;
        std::map<std::string, double> cls_S;
        for (const auto& c : kept) {
            const double sc = use_ea ? c.ea : c.S;
            auto it = cls_S.find(c.cls);
            if (it == cls_S.end()) cls_S[c.cls] = sc;
            else it->second = std::max(it->second, sc);
        }
        double smax = -std::numeric_limits<double>::infinity();
        for (const auto& kv : cls_S) smax = std::max(smax, kv.second);
        std::map<std::string, double> cls_W;
        std::map<std::string, int> cls_id;
        double zc = 0.0;
        for (const auto& kv : cls_S) { cls_W[kv.first] = std::exp((kv.second - smax) / sig); zc += cls_W[kv.first]; }
        for (auto& kv : cls_W) { kv.second /= zc; cls_id.emplace(kv.first, static_cast<int>(cls_id.size())); }

        std::string gkey = "a";
        for (int i = 0; i < m_atomcount; ++i)
            gkey += fmt::format(",{}", group_of_atom[i]);
        std::vector<Place> uses;
        if (G == 1 && master_idx >= 0) {
            uses.push_back({ master_idx, 0, 1.0 });
            kept.clear();
        }
        for (const auto& c : kept) {
            if (cls_W[c.cls] < 1e-14)
                continue;   // negligible class: not worth a GFN-FF evaluation (weight is a topology constant)
            std::string key = gkey + "|q";
            for (int g = 0; g < G; ++g)
                key += fmt::format(",{}", static_cast<int>(c.qg[g]));
            GFNFF* ff = fragVariant(key, group_of_atom, G, c.qg);
            if (!ff)
                continue;
            for (int v = 0; v < static_cast<int>(m_frag_variants.size()); ++v)
                if (m_frag_variants[v].key == key) { uses.push_back({ v, cls_id[c.cls], cls_W[c.cls] }); break; }
        }
        if (uses.empty())
            continue;
        Wc.push_back(W);
        corner_mask.push_back(mask);
        corner_uses.push_back(std::move(uses));
    }
    if (Wc.empty()) {
        CurcumaLogger::warn("GFNFF frag_charge_model ensemble: no usable charge variant - keeping the reference energy");
        return e_master;
    }
    // renormalise only if a corner had to be dropped (a variant failed); otherwise sum W = 1
    const double Wsum = std::accumulate(Wc.begin(), Wc.end(), 0.0);

    // ---- evaluate the needed variants ------------------------------------------------------
    std::vector<char> needed(m_frag_variants.size(), 0);
    for (const auto& u : corner_uses)
        for (const auto& pl : u)
            needed[pl.variant] = 1;
    const Matrix geom_ang = m_geometry_bohr * BOHR_TO_ANGSTROM;
    for (size_t v = 0; v < m_frag_variants.size(); ++v) {
        if (!needed[v])
            continue;
        FragVariant& fv = m_frag_variants[v];
        if (!fv.ff) {   // the master itself (see master_idx)
            fv.energy = e_master;
            fv.ok = std::isfinite(e_master) && !m_eeq_solve_failed;
            if (gradient)
                fv.gradient = m_gradient;
            fv.charges = m_charges;
            if (m_workspace)
                fv.comp = m_workspace->energyComponents();
            continue;
        }
        FragQuiet quiet(0);
        fv.ff->UpdateMolecule(geom_ang);
        fv.energy = fv.ff->Calculation(gradient);
        fv.ok = std::isfinite(fv.energy) && !fv.ff->eeqSolveFailed();
        if (gradient)
            fv.gradient = fv.ff->Gradient();
        fv.charges = fv.ff->Charges();
        if (fv.ff->m_workspace)
            fv.comp = fv.ff->m_workspace->energyComponents();
    }

    // ---- blend ---------------------------------------------------------------------------
    // E_c = sum_class Omega_class sum_{p in class} omega_p E_p, omega = softmax(-E/tau) in the class;
    // dE_c/dx = sum Omega omega_p (1 - (E_p - E_class)/tau) dE_p/dx (Omega is a topology constant).
    const double tau = m_frag_tau_eh;
    const int nv = static_cast<int>(m_frag_variants.size());
    std::vector<double> a(nv, 0.0);       // energy weight per variant (sum a = 1)
    std::vector<double> b(nv, 0.0);       // gradient weight per variant
    std::vector<double> Ec(Wc.size(), 0.0);
    double E = 0.0;
    std::vector<char> corner_ok(Wc.size(), 0);
    double Wok = 0.0;
    for (size_t c = 0; c < Wc.size(); ++c) {
        std::map<int, std::vector<size_t>> by_cls;
        for (size_t p = 0; p < corner_uses[c].size(); ++p)
            if (m_frag_variants[corner_uses[c][p].variant].ok)
                by_cls[corner_uses[c][p].cls].push_back(p);
        if (by_cls.empty())
            continue;
        double Om_sum = 0.0;
        for (const auto& kv : by_cls) Om_sum += corner_uses[c][kv.second.front()].Omega;
        double e_c = 0.0;
        std::vector<std::pair<int, double>> wa, wb;   // (variant, a-weight), (variant, b-weight) inside the corner
        for (const auto& kv : by_cls) {
            const double Om = corner_uses[c][kv.second.front()].Omega / Om_sum;
            double emin = std::numeric_limits<double>::infinity();
            for (size_t p : kv.second) emin = std::min(emin, m_frag_variants[corner_uses[c][p].variant].energy);
            double Z = 0.0;
            std::vector<double> om;
            for (size_t p : kv.second) { om.push_back(std::exp(-(m_frag_variants[corner_uses[c][p].variant].energy - emin) / tau)); Z += om.back(); }
            double e_cls = 0.0;
            for (size_t n = 0; n < om.size(); ++n) { om[n] /= Z; e_cls += om[n] * m_frag_variants[corner_uses[c][kv.second[n]].variant].energy; }
            e_c += Om * e_cls;
            for (size_t n = 0; n < om.size(); ++n) {
                const int v = corner_uses[c][kv.second[n]].variant;
                wa.push_back({ v, Om * om[n] });
                wb.push_back({ v, Om * om[n] * (1.0 - (m_frag_variants[v].energy - e_cls) / tau) });
            }
        }
        Ec[c] = e_c;
        corner_ok[c] = 1;
        Wok += Wc[c];
        for (const auto& x : wa) a[x.first] += Wc[c] * x.second;
        for (const auto& x : wb) b[x.first] += Wc[c] * x.second;
        E += Wc[c] * e_c;
    }
    if (Wok <= 0.0) {
        CurcumaLogger::warn("GFNFF frag_charge_model ensemble: every charge variant failed - keeping the reference energy");
        return e_master;
    }
    // a corner whose variants all failed is dropped and the rest renormalised (warned once above)
    E /= Wok;
    for (int v = 0; v < nv; ++v) { a[v] /= Wok; b[v] /= Wok; }
    (void)Wsum;

    // components and charges: the same weights as the energy
    FFEnergyComponents comp;
    Vector q = Vector::Zero(m_atomcount);
    for (int v = 0; v < nv; ++v) {
        if (a[v] == 0.0) continue;
        addScaled(comp, m_frag_variants[v].comp, a[v]);
        if (m_frag_variants[v].charges.size() == m_atomcount)
            q += a[v] * m_frag_variants[v].charges;
    }

    if (gradient) {
        Matrix g = Matrix::Zero(m_atomcount, 3);
        for (int v = 0; v < nv; ++v)
            if (b[v] != 0.0 && m_frag_variants[v].gradient.rows() == m_atomcount)
                g += b[v] * m_frag_variants[v].gradient;
        // dW_c/dlambda_e contribution
        for (int e = 0; e < k; ++e) {
            double coef = 0.0;   // dE/dlambda_e
            for (size_t c = 0; c < Wc.size(); ++c) {
                if (!corner_ok[c]) continue;
                double d = ((corner_mask[c] >> e) & 1u) ? -1.0 : 1.0;
                for (int e2 = 0; e2 < k; ++e2) {
                    if (e2 == e) continue;
                    d *= ((corner_mask[c] >> e2) & 1u) ? (1.0 - edges[e2].lambda) : edges[e2].lambda;
                }
                coef += d * (Ec[c] - E) / Wok;   // d(sum W E / sum W)/dlambda
            }
            if (coef == 0.0) continue;
            const auto& cs = edges[e].contacts;
            for (size_t m = 0; m < cs.size(); ++m) {
                if (cs[m].dL == 0.0) continue;
                double dlam = cs[m].dL;
                for (size_t m2 = 0; m2 < cs.size(); ++m2)
                    if (m2 != m) dlam *= cs[m2].L;
                if (dlam == 0.0) continue;
                const Eigen::RowVector3d rij = m_geometry_bohr.row(cs[m].i) - m_geometry_bohr.row(cs[m].j);
                const double r = rij.norm();
                const Eigen::RowVector3d ds = rij / (r * cs[m].thr);
                g.row(cs[m].i) += coef * dlam * ds;
                g.row(cs[m].j) -= coef * dlam * ds;
            }
        }
        m_gradient = g;
    }

    m_energy_total = E;
    m_charges = q;
    m_frag_blend_comp = comp;
    m_frag_blend_valid = true;
    m_frag_last_nvariants = 0;
    for (int v = 0; v < nv; ++v)
        if (a[v] > 0.0) ++m_frag_last_nvariants;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format("GFN-FF frag_charge_model ensemble: {} fragments, {} window contact(s), {} corner(s), {} variant(s)",
            F, k, Wc.size(), m_frag_last_nvariants));
        for (int e = 0; e < k; ++e)
            CurcumaLogger::info(fmt::format("  contact {}-{}: lambda = {:.6f} ({} atom pair(s))", edges[e].f + 1, edges[e].g + 1,
                edges[e].lambda, edges[e].contacts.size()));
        for (int v = 0; v < nv; ++v)
            if (a[v] > 0.0)
                CurcumaLogger::info(fmt::format("  variant {:<40s} weight {:.6f}  E - E(reference rule) = {:+.4f} kcal/mol",
                    m_frag_variants[v].key, a[v], (m_frag_variants[v].energy - e_master) * 627.5094740631));
        CurcumaLogger::info(fmt::format("  blended E = {:.10f} Eh, E - E(reference rule) = {:+.4f} kcal/mol", E, (E - e_master) * 627.5094740631));
    }
    return E;
}
