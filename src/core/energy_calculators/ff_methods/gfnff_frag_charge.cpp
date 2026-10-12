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
 * testing pending. Documentation: docs/FRAG_CHARGE_MODEL.md; design and validation record on the
 * reactff2-llm branch (test_cases/revgfnff/_log/FRAG_CHARGE_STATUS.md there).
 *
 * Port to master (Sep 27, 2026, Known Issue #31): taken verbatim from reactff2-llm commit
 * 837266d8 except (a) FragQuiet scopes the process-wide CurcumaLogger verbosity (master has no
 * per-thread verbosity) and (b) addScaled() omits the rev-gfnff-only energy components.
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
 *   corner c   = subset of the window edges that are merged (2^k, the same pattern as the reactive
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
#include <cfloat>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>

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

/// scoped verbosity (the variants run silently). master's CurcumaLogger has only the process-wide
/// level, so this is the same save/restore pattern as the other nested silent helpers (Known Issue #3).
struct FragQuiet {
    int saved;
    explicit FragQuiet(int level) : saved(CurcumaLogger::get_verbosity()) { CurcumaLogger::set_verbosity(level); }
    ~FragQuiet() { CurcumaLogger::set_verbosity(saved); }
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

/// rank contribution of one group carrying c units: [odd electron count] + (G+1)*[no electron left].
/// Same expressions as the pre-Oct-2026 enumeration (qg = unit * c as double, ne as long).
inline int fragRankTerm(long nelec, double unit, int c, int G)
{
    const double qg = unit * c;
    const long ne = nelec - static_cast<long>(qg);
    int r = 0;
    if (ne & 1L) ++r;
    if (ne <= 0) r += G + 1;
    return r;
}

/// score contribution (free-charge pre-selection) of one group carrying c units, same float
/// expression as the enumeration: (unit*c) * qfree / nq, added to S in group order.
inline double fragScoreTerm(double unit, int c, double qfree, int nq)
{
    const double qg = unit * c;
    return qg * qfree / static_cast<double>(nq);
}

/// sparse placement prefix: (group, units) for the groups with c > 0, groups ascending
using FragNZ = std::vector<std::pair<int, int>>;

/// dense(a) > dense(b) lexicographically, i.e. a comes first in the old enumeration order
/// (compositions() runs k = units .. 0 at the first slot, recursively)
bool fragLexGreater(const FragNZ& a, const FragNZ& b)
{
    size_t i = 0;
    for (; i < a.size() && i < b.size(); ++i) {
        if (a[i].first != b[i].first)
            return a[i].first < b[i].first;   // a is non-zero at an earlier group where b is zero
        if (a[i].second != b[i].second)
            return a[i].second > b[i].second;
    }
    return i < a.size() && i >= b.size();
}

void addScaled(FFEnergyComponents& a, const FFEnergyComponents& b, double w)
{
    a.bond += w * b.bond; a.angle += w * b.angle; a.dihedral += w * b.dihedral; a.inversion += w * b.inversion;
    a.dispersion += w * b.dispersion; a.vdw += w * b.vdw; a.rep += w * b.rep;
    a.bonded_rep += w * b.bonded_rep; a.nonbonded_rep += w * b.nonbonded_rep;
    a.coulomb += w * b.coulomb; a.hbond += w * b.hbond; a.xbond += w * b.xbond;
    a.atm += w * b.atm; a.batm += w * b.batm; a.stors += w * b.stors;
    a.hbond_case1 += w * b.hbond_case1; a.hbond_case2 += w * b.hbond_case2;
    a.hbond_case3 += w * b.hbond_case3; a.hbond_case4 += w * b.hbond_case4;
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
    return fm(zi) * fm(zj) * rthr * rco;
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

std::vector<std::vector<int>> GFNFF::fragEnumeratePlacements(int nq, double unit, const std::vector<long>& nelec,
    const std::vector<double>& qfree_g, int max_placements)
{
    // The pre-Oct-2026 code of fragEnsembleBlend(), unchanged in its arithmetic: all
    // compositions, minimum rank, stable sort by S (only if more than max_placements), truncate.
    const int G = static_cast<int>(nelec.size());
    std::vector<std::vector<int>> comps;
    std::vector<int> cur;
    compositions(nq, G, cur, comps);
    struct Cand { size_t idx; int radicals; double S; };
    std::vector<Cand> cands;
    int rmin = std::numeric_limits<int>::max();
    for (size_t n = 0; n < comps.size(); ++n) {
        Cand c{ n, 0, 0.0 };
        for (int g = 0; g < G; ++g) {
            c.radicals += fragRankTerm(nelec[g], unit, comps[n][g], G);
            c.S += fragScoreTerm(unit, comps[n][g], qfree_g[g], nq);
        }
        rmin = std::min(rmin, c.radicals);
        cands.push_back(c);
    }
    std::vector<Cand> kept;
    for (const auto& c : cands)
        if (c.radicals == rmin) kept.push_back(c);
    if (static_cast<int>(kept.size()) > max_placements) {
        std::stable_sort(kept.begin(), kept.end(), [](const Cand& x, const Cand& y) { return x.S > y.S; });
        kept.resize(max_placements);
    }
    std::vector<std::vector<int>> out;
    for (const auto& c : kept) out.push_back(comps[c.idx]);
    return out;
}

std::vector<std::vector<int>> GFNFF::fragSelectPlacements(int nq, double unit, const std::vector<long>& nelec,
    const std::vector<double>& qfree_g, int max_placements, bool* exact_tiebreak)
{
    // Claude Generated (Oct 2026). Replaces the full enumeration, which held C(nq+G-1, nq)
    // vectors of length G (497 groups, nq = 3: 2e7 vectors, ~41 GB) before any filter acted.
    //
    // Both the rank and the score S are sums of per-group terms, so the placements are paths
    // through the states (g, u) = (groups decided, units used). Rank: exact minimum by DP
    // (suffix minima suf[g][r], prefix minima pre[u]). Only states on a minimum-rank path are
    // kept; every prefix there has the same rank, so prefixes in one state differ only in S and
    // in their enumeration order (lexicographic, decided inside the prefix).
    // Selection: prefix B in a state can be dropped if at least K = max_placements other
    // prefixes A beat it for EVERY common suffix. A beats B for every suffix if
    //   (S_A >= S_B and A before B in enumeration order)   - float addition is monotone, a tie
    //                                                          is then broken by the order, or
    //   (S_A - S_B > delta)                                  - the gap survives all remaining
    //                                                          roundings (delta bound below).
    // Then the true top K are never dropped, and a final sort with the old comparator over the
    // survivors gives exactly the old list. If there are at most K placements nothing is
    // dropped and the list is returned in enumeration order (the old code did not sort then).
    if (exact_tiebreak) *exact_tiebreak = true;
    const int G = static_cast<int>(nelec.size());
    const int K = std::max(1, max_placements);
    if (G == 0 || nq < 0)
        return {};
    const int Q = nq + 1;
    constexpr int INF = std::numeric_limits<int>::max() / 4;

    // suffix minima: suf[g*Q + r] = minimum rank of r units over groups g .. G-1
    std::vector<int> suf(static_cast<size_t>(G + 1) * Q, INF);
    suf[static_cast<size_t>(G) * Q + 0] = 0;
    for (int g = G - 1; g >= 0; --g)
        for (int r = 0; r <= nq; ++r) {
            int best = INF;
            for (int c = 0; c <= r; ++c) {
                const int tail = suf[static_cast<size_t>(g + 1) * Q + (r - c)];
                if (tail < INF) best = std::min(best, fragRankTerm(nelec[g], unit, c, G) + tail);
            }
            suf[static_cast<size_t>(g) * Q + r] = best;
        }
    const int rmin = suf[nq];

    // delta: every partial sum and term is bounded by M = sum |qfree_g| (c <= nq); each of the
    // <= G remaining additions shifts a difference by at most 2 u M (u = DBL_EPSILON/2). Doubled.
    double M = 0.0;
    for (double q : qfree_g) M += std::abs(q);
    const double delta = 4.0 * (G + 1) * DBL_EPSILON * M * 1.01;
    const size_t cap = static_cast<size_t>(64) * K + 64;   // per state; only near-tie clusters reach it

    struct Entry { double S; FragNZ nz; };
    auto beats = [](const Entry& a, const Entry& b) {   // the old comparator: S descending, then enumeration order
        if (a.S > b.S) return true;
        if (b.S > a.S) return false;
        return fragLexGreater(a.nz, b.nz);
    };
    std::vector<int> pre(Q, INF), npre(Q);
    pre[0] = 0;
    std::vector<unsigned long long> cnt(Q, 0), ncnt(Q);
    cnt[0] = 1;
    constexpr unsigned long long CNT_SAT = 1ull << 62;
    std::vector<std::vector<Entry>> cur(Q), nxt(Q);
    cur[0].push_back({ 0.0, {} });
    for (int g = 0; g < G; ++g) {
        const size_t sg = static_cast<size_t>(g) * Q, sg1 = static_cast<size_t>(g + 1) * Q;
        std::fill(npre.begin(), npre.end(), INF);
        for (int u = 0; u <= nq; ++u) {
            if (pre[u] >= INF || pre[u] + suf[sg + (nq - u)] != rmin) continue;   // not on a minimum path
            for (int c = 0; u + c <= nq; ++c)
                npre[u + c] = std::min(npre[u + c], pre[u] + fragRankTerm(nelec[g], unit, c, G));
        }
        for (int t = 0; t <= nq; ++t) { nxt[t].clear(); ncnt[t] = 0; }
        for (int u = 0; u <= nq; ++u) {
            if (pre[u] >= INF || pre[u] + suf[sg + (nq - u)] != rmin) continue;
            for (int c = 0; u + c <= nq; ++c) {
                const int t = u + c;
                if (npre[t] >= INF || npre[t] + suf[sg1 + (nq - t)] != rmin) continue;
                if (pre[u] + fragRankTerm(nelec[g], unit, c, G) != npre[t]) continue;
                ncnt[t] = std::min(CNT_SAT, ncnt[t] + cnt[u]);
                const double term = fragScoreTerm(unit, c, qfree_g[g], nq);
                for (const Entry& e : cur[u]) {
                    Entry x{ e.S + term, e.nz };
                    if (c > 0) x.nz.emplace_back(g, c);
                    nxt[t].push_back(std::move(x));
                }
            }
        }
        // prune every state: drop an entry that at least K others beat for every suffix
        for (int t = 0; t <= nq; ++t) {
            std::vector<Entry>& L = nxt[t];
            if (L.size() <= static_cast<size_t>(K)) continue;
            std::vector<char> keep(L.size(), 1);
            for (size_t b = 0; b < L.size(); ++b) {
                int dom = 0;
                for (size_t a = 0; a < L.size() && dom < K; ++a)
                    if (a != b && ((L[a].S >= L[b].S && fragLexGreater(L[a].nz, L[b].nz)) || (L[a].S - L[b].S > delta)))
                        ++dom;
                keep[b] = dom < K;
            }
            std::vector<Entry> kept;
            for (size_t b = 0; b < L.size(); ++b)
                if (keep[b]) kept.push_back(std::move(L[b]));
            if (kept.size() > cap) {   // near-tie cluster wider than the cap: keep the best by the old comparator
                std::sort(kept.begin(), kept.end(), beats);
                kept.resize(cap);
                if (exact_tiebreak) *exact_tiebreak = false;
            }
            L = std::move(kept);
        }
        std::swap(cur, nxt);
        pre = npre;
        cnt = ncnt;
    }
    std::vector<Entry>& fin = cur[nq];
    if (cnt[nq] <= static_cast<unsigned long long>(K))
        std::sort(fin.begin(), fin.end(), [](const Entry& a, const Entry& b) { return fragLexGreater(a.nz, b.nz); });
    else {
        std::sort(fin.begin(), fin.end(), beats);
        if (fin.size() > static_cast<size_t>(K)) fin.resize(K);
    }
    std::vector<std::vector<int>> out;
    out.reserve(fin.size());
    for (const Entry& e : fin) {
        std::vector<int> v(G, 0);
        for (const auto& p : e.nz) v[p.first] = p.second;
        out.push_back(std::move(v));
    }
    return out;
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
    // the topology was built (topology-history dependence).
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
    // Placement rule (docs/FRAG_CHARGE_MODEL.md): GFN-FF's absolute energies of differently
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
    // ---- active corners and the variant memory budget (Claude Generated, Oct 2026) ------------
    // Each variant is a complete GFNFF instance: measured 64-74 bytes * N^2 (VmRSS per variant,
    // N = 2822: ~0.5 GB, N = 6200: ~2.6 GB), estimated as 72 bytes * N^2. The budget
    // frag_charge_variant_memory_mb gives the number of variants held at once (>= 1); placements
    // per corner and, if needed, the lightest corners are cut to fit (warning once).
    std::vector<unsigned> act_mask;
    std::vector<double> act_W;
    for (unsigned mask = 0; mask < (1u << k); ++mask) {
        double W = 1.0;
        for (int e = 0; e < k; ++e)
            W *= ((mask >> e) & 1u) ? (1.0 - edges[e].lambda) : edges[e].lambda;
        if (W <= 0.0)
            continue;   // exact: a corner with W = 0 has dW/dx = 0 too (smootherstep is C2 at both ends)
        act_mask.push_back(mask);
        act_W.push_back(W);
    }
    const double variant_bytes = std::max(1.0, 72.0 * static_cast<double>(m_atomcount) * static_cast<double>(m_atomcount));
    const double budget_bytes = static_cast<double>(m_frag_variant_memory_mb) * 1024.0 * 1024.0;
    const int max_held = static_cast<int>(std::min(1.0e6, std::max(1.0, std::floor(budget_bytes / variant_bytes))));
    std::vector<char> act_use(act_mask.size(), 1);
    bool mem_limited = false;
    if (static_cast<int>(act_mask.size()) > max_held) {
        // keep the max_held heaviest corners (ties: lower mask); the rest is renormalised away (Wok)
        std::vector<size_t> order(act_mask.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](size_t x, size_t y) { return act_W[x] > act_W[y]; });
        std::fill(act_use.begin(), act_use.end(), 0);
        for (int i = 0; i < max_held; ++i) act_use[order[i]] = 1;
        mem_limited = true;
    }
    const int n_corners = std::min(static_cast<int>(act_mask.size()), max_held);
    const int max_place = std::max(1, std::min(m_frag_max_placements, n_corners > 0 ? max_held / n_corners : max_held));

    // ---- plan: the placements every corner needs, before any variant is built ----------------
    struct Cand { std::vector<double> qg; int radicals; double S; std::string cls; };
    struct CornerPlan {
        unsigned mask = 0;
        double W = 0.0;
        int G = 0;
        bool master = false;                 // the one-group corner of a one-fragment master
        std::vector<int> group_of_atom;
        std::vector<Cand> kept;
        std::map<std::string, double> cls_W;
        std::map<std::string, int> cls_id;
        std::vector<std::string> keys;       // per kept candidate, "" = negligible class (skipped)
    };
    std::vector<CornerPlan> plans;
    for (size_t ci = 0; ci < act_mask.size(); ++ci) {
        if (!act_use[ci])
            continue;
        const unsigned mask = act_mask[ci];
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

        // placements: minimum rank, then at most frag_charge_max_placements by the free-charge
        // score S. rank = bare_nuclei * (G + 1) + radicals: a placement that strips a group of ALL
        // its electrons (a bare proton next to CO in a pass-2-split CH2O2+) is never preferred to
        // one that does not, whatever the radical count; then the fewest odd-electron groups.
        // fragSelectPlacements() gives the list of the former full enumeration (same members,
        // same order) without enumerating all C(nq+G-1, nq) compositions (Oct 2026).
        std::vector<std::vector<int>> comps;
        bool exact_tiebreak = true;
        bool enumerate = m_frag_enumerate;
        if (enumerate) {
            double ncomp = 1.0;
            for (int t = 1; t <= nq; ++t) ncomp = ncomp * (G - 1 + t) / t;
            if (ncomp > 2.0e6) {
                if (!m_frag_warned_select) {
                    CurcumaLogger::warn(fmt::format("GFNFF frag_charge_enumerate: {:.3g} charge distributions over {} groups exceed the "
                                                    "enumeration limit (2e6) - using the bounded selection", ncomp, G));
                    m_frag_warned_select = true;
                }
                enumerate = false;
            }
        }
        if (enumerate)
            comps = fragEnumeratePlacements(nq, unit, nelec, qfree_g, m_frag_max_placements);
        else
            comps = fragSelectPlacements(nq, unit, nelec, qfree_g, m_frag_max_placements, &exact_tiebreak);
        if (!exact_tiebreak && !m_frag_warned_select) {
            CurcumaLogger::warn("GFNFF frag_charge_model ensemble: many placements with free-charge scores equal within rounding; "
                                "their tie order may differ from the full enumeration");
            m_frag_warned_select = true;
        }
        std::vector<Cand> kept;
        for (const auto& cmp : comps) {
            Cand c;
            c.qg.assign(G, 0.0);
            c.radicals = 0;
            c.S = 0.0;
            std::vector<std::string> parts;
            for (int g = 0; g < G; ++g) {
                c.qg[g] = unit * cmp[g];
                c.radicals += fragRankTerm(nelec[g], unit, cmp[g], G);
                c.S += fragScoreTerm(unit, cmp[g], qfree_g[g], nq);
                if (cmp[g] != 0) {
                    std::string zs = fmt::format("{:+d}:", static_cast<int>(c.qg[g]));
                    for (int z : zlist[g]) zs += fmt::format("{}.", z);
                    parts.push_back(zs);
                }
            }
            std::sort(parts.begin(), parts.end());
            for (const auto& ps : parts) c.cls += ps + "/";
            kept.push_back(std::move(c));
        }
        // memory budget: fewer placements per corner, the lowest score S dropped first - exactly
        // what frag_charge_max_placements = max_place would select
        if (static_cast<int>(kept.size()) > max_place) {
            std::stable_sort(kept.begin(), kept.end(), [](const Cand& x, const Cand& y) { return x.S > y.S; });
            kept.resize(max_place);
            mem_limited = true;
        }
        // class weights Omega = softmax(S_class / sigma), S_class = max over its members
        std::map<std::string, double> cls_S;
        for (const auto& c : kept) {
            auto it = cls_S.find(c.cls);
            if (it == cls_S.end()) cls_S[c.cls] = c.S;
            else it->second = std::max(it->second, c.S);
        }
        double smax = -std::numeric_limits<double>::infinity();
        for (const auto& kv : cls_S) smax = std::max(smax, kv.second);
        std::map<std::string, double> cls_W;
        std::map<std::string, int> cls_id;
        double zc = 0.0;
        for (const auto& kv : cls_S) { cls_W[kv.first] = std::exp((kv.second - smax) / m_frag_sigma); zc += cls_W[kv.first]; }
        for (auto& kv : cls_W) { kv.second /= zc; cls_id.emplace(kv.first, static_cast<int>(cls_id.size())); }

        CornerPlan plan;
        plan.mask = mask;
        plan.W = act_W[ci];
        plan.G = G;
        if (G == 1 && master_idx >= 0) {
            plan.master = true;
            kept.clear();
        }
        std::string gkey = "a";
        for (int i = 0; i < m_atomcount; ++i)
            gkey += fmt::format(",{}", group_of_atom[i]);
        for (const auto& c : kept) {
            std::string key;
            if (cls_W[c.cls] >= 1e-14) {   // negligible class: not worth a GFN-FF evaluation (weight is a topology constant)
                key = gkey + "|q";
                for (int g = 0; g < G; ++g)
                    key += fmt::format(",{}", static_cast<int>(c.qg[g]));
            }
            plan.keys.push_back(std::move(key));
        }
        plan.group_of_atom = std::move(group_of_atom);
        plan.kept = std::move(kept);
        plan.cls_W = std::move(cls_W);
        plan.cls_id = std::move(cls_id);
        plans.push_back(std::move(plan));
    }
    if (mem_limited && !m_frag_warned_memory) {
        CurcumaLogger::warn(fmt::format("GFNFF frag_charge_model ensemble: the charge variants of {} atoms (~{:.0f} MB each) exceed "
                                        "frag_charge_variant_memory_mb = {}; evaluating {} corner(s) with at most {} placement(s) each "
                                        "instead of {} corner(s) x {}",
            m_atomcount, variant_bytes / (1024.0 * 1024.0), m_frag_variant_memory_mb, n_corners, max_place, act_mask.size(),
            m_frag_max_placements));
        m_frag_warned_memory = true;
    }

    // ---- release variants this step does not need, if keeping them would exceed the budget -----
    {
        std::set<std::string> needed_keys;
        for (const auto& pl : plans)
            for (const auto& key : pl.keys)
                if (!key.empty()) needed_keys.insert(key);
        int held = 0, held_unneeded = 0;
        std::set<std::string> present;
        for (const auto& v : m_frag_variants)
            if (v.ff) {
                ++held;
                present.insert(v.key);
                if (!needed_keys.count(v.key)) ++held_unneeded;
            }
        int to_build = 0;
        for (const auto& key : needed_keys)
            if (!present.count(key)) ++to_build;
        if (held_unneeded > 0 && held + to_build > max_held) {
            m_frag_variants.erase(std::remove_if(m_frag_variants.begin(), m_frag_variants.end(),
                                      [&](const FragVariant& v) { return v.ff && !needed_keys.count(v.key); }),
                m_frag_variants.end());
            master_idx = -1;
            for (int v = 0; v < static_cast<int>(m_frag_variants.size()); ++v)
                if (m_frag_variants[v].key == "master") master_idx = v;
        }
    }

    // ---- build the variants and the corner list --------------------------------------------
    for (const auto& pl : plans) {
        std::vector<Place> uses;
        if (pl.master)
            uses.push_back({ master_idx, 0, 1.0 });
        for (size_t n = 0; n < pl.kept.size(); ++n) {
            const Cand& c = pl.kept[n];
            const std::string& key = pl.keys[n];
            if (key.empty())
                continue;
            GFNFF* ff = fragVariant(key, pl.group_of_atom, pl.G, c.qg);
            if (!ff)
                continue;
            for (int v = 0; v < static_cast<int>(m_frag_variants.size()); ++v)
                if (m_frag_variants[v].key == key) { uses.push_back({ v, pl.cls_id.at(c.cls), pl.cls_W.at(c.cls) }); break; }
        }
        if (uses.empty())
            continue;
        Wc.push_back(pl.W);
        corner_mask.push_back(pl.mask);
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
