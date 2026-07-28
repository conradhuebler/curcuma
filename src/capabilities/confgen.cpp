/*
 * <Conformer proposals from an energy-decomposed ensemble>
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

#include "confgen.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/fileiterator.h"
#include "src/capabilities/optimisation/geometry_restraints.h"
#include "src/capabilities/optimizer_factory.h"
#include "src/capabilities/rmsd/rmsd_functions.h"
#include "src/core/elements.h"

#include <algorithm>
#include <cmath>
#include <fmt/core.h>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <numeric>
#include <set>

namespace {
/// Hartree -> kJ/mol (same literal as the rest of the conformer code, see confsearch.cpp).
constexpr double kEh2kJ = 2625.5;
/// Boltzmann constant in Hartree/Kelvin.
constexpr double kBoltzmannEh = 3.166811563e-6;
}

ConfGen::ConfGen(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("confgen"), controller, silent)
    , m_config("confgen", controller)
{
    UpdateController(controller);
}

void ConfGen::setFile(const std::string& filename)
{
    CurcumaMethod::setFile(filename);
}

void ConfGen::LoadControlJson()
{
    m_method = m_config.get<std::string>("method");
    m_charge = m_config.get<int>("charge");
    m_spin = m_config.get<int>("spin");
    m_threads = m_config.get<int>("threads");
    m_state_tolerance = m_config.get<double>("state_tolerance");
    m_temperature = m_config.get<double>("temperature");
    m_min_pairs = m_config.get<int>("min_pairs");
    m_report_threshold = m_config.get<double>("report_threshold");
    m_cv_folds = m_config.get<int>("cv_folds");
    m_couplings = m_config.get<bool>("couplings");
    m_generate = m_config.get<bool>("generate");
    m_max_proposals = m_config.get<int>("max_proposals");
    m_proposal_templates = m_config.get<int>("proposal_templates");
    m_proposal_depth = m_config.get<int>("proposal_depth");
    m_clash_factor = m_config.get<double>("clash_factor");
    m_new_rmsd = m_config.get<double>("new_rmsd");
    m_topology_factor = m_config.get<double>("topology_factor");
    // Claude Generated (Jul 2026): restrained build (P0)
    m_restrained_build = m_config.get<bool>("restrained_build");
    m_restraint_force = m_config.get<double>("restraint_force");
    m_restraint_max_iterations = m_config.get<int>("restraint_max_iterations");
}

bool ConfGen::analyseEnsemble()
{
    // ---- 1. read the ensemble -------------------------------------------------------------------
    std::vector<Molecule> input;
    {
        FileIterator file(Filename());
        while (!file.AtEnd()) {
            Molecule mol = file.Next();
            if (mol.AtomCount() > 0) {
                mol.setCharge(m_charge);
                mol.setSpin(m_spin);
                input.push_back(mol);
            }
        }
    }
    if (input.size() < 2) {
        CurcumaLogger::error_fmt("ConfGen: need at least 2 structures in {} (found {})",
            Filename(), static_cast<int>(input.size()));
        return false;
    }
    CurcumaLogger::result_fmt("ConfGen: {} structures read from {}", static_cast<int>(input.size()), Filename());

    // ---- 2. torsion space from the FIRST structure (the reference topology) --------------------
    m_torsions = TorsionSpace::detectTorsions(input.front());
    if (m_torsions.empty()) {
        CurcumaLogger::warn("ConfGen: no rotatable torsion found -- nothing to analyse. A rigid "
                            "molecule has no torsion states to recombine.");
        return false;
    }
    CurcumaLogger::result_fmt("ConfGen: {} rotatable torsion(s) detected", static_cast<int>(m_torsions.size()));

    // ---- 3. one single point per member, ONE shared calculator ---------------------------------
    // Sharing the calculator is essential, not an optimisation: it sets up topology, parameters and
    // (for GFN-FF) the charge model once from the reference structure and only updates coordinates.
    // Every frame therefore reports the SAME terms, which is what makes their differences comparable.
    // A frame whose bond topology differs from the reference is dropped -- its terms would describe a
    // different molecule.
    json calc_config;
    calc_config["method"] = m_method;
    calc_config["threads"] = m_threads;
    calc_config["charge"] = m_charge;
    calc_config["spin"] = m_spin;
    m_calculator = std::make_unique<EnergyCalculator>(m_method, calc_config);
    EnergyCalculator& calculator = *m_calculator;
    calculator.setMolecule(input.front().getMolInfo());

    // Claude Generated (Jul 2026): use the EXPLICIT topology factor here too, not
    // Molecule::DistanceMatrix()'s 1.5. ConfGen already knew 1.5 is unusable for this decision (it
    // is why -topology_factor exists and why the proposal checks use it), but its own ensemble
    // filter still used the default -- and 1.5 counts a compressed 1-3 contact as a bond, so it
    // flips on a two-degree bending fluctuation. Measured on a 142-structure peptide ensemble:
    // 28 structures dropped as "topology differs" at 1.5, ZERO at 1.3 -- including the reference
    // structure this analysis was run to compare against.
    const std::vector<std::pair<int, int>> reference_topology
        = topologyFingerprint(input.front(), m_topology_factor);
    int rejected_topology = 0;

    m_frames.clear();
    m_frames.reserve(input.size());
    for (const Molecule& mol : input) {
        Frame frame;
        frame.molecule = mol;

        if (mol.AtomCount() != input.front().AtomCount()) {
            rejected_topology++;
            continue;
        }
        if (topologyFingerprint(mol, m_topology_factor) != reference_topology) {
            rejected_topology++;
            continue;
        }

        calculator.updateGeometry(mol.getGeometry());
        frame.energy = calculator.CalculateEnergy(false);
        if (calculator.Error()) {
            CurcumaLogger::error("ConfGen: " + calculator.ErrorMessage());
            return false;
        }
        const json decomposition = calculator.getEnergyDecomposition();
        for (auto& [key, value] : decomposition.items())
            if (value.is_number())
                frame.terms[key] = value.get<double>();

        frame.angles = TorsionSpace::dihedrals(mol.getGeometry(), m_torsions);
        frame.valid = true;
        m_frames.push_back(std::move(frame));
    }
    if (rejected_topology > 0)
        CurcumaLogger::warn_fmt("ConfGen: {} structure(s) dropped -- their bond topology differs from "
                                "the reference structure, so their energy terms are not comparable",
            rejected_topology);
    if (m_frames.size() < 2) {
        CurcumaLogger::error("ConfGen: fewer than 2 usable structures left after the topology check");
        return false;
    }

    // Stable term order for all tables: the keys of the first frame, sorted, non-zero ones first.
    std::vector<std::string> keys;
    for (const auto& [key, value] : m_frames.front().terms)
        keys.push_back(key);
    std::sort(keys.begin(), keys.end());
    m_term_names = keys;

    // Consistency check: for a force field the components must add up to the total energy. If they
    // do not, a term is missing from the decomposition and every attribution below is incomplete --
    // so this is reported rather than silently ignored.
    double sum_terms = 0.0;
    for (const auto& name : m_term_names)
        sum_terms += m_frames.front().terms.at(name);
    const double mismatch = std::abs(sum_terms - m_frames.front().energy) * kEh2kJ;
    if (mismatch > 1.0)
        CurcumaLogger::warn_fmt("ConfGen: the energy decomposition does not sum to the total energy "
                                "({:.2f} vs {:.2f} kJ/mol, difference {:.2f}). Some contribution is not "
                                "exposed as a term; per-term attributions below are incomplete.",
            sum_terms * kEh2kJ, m_frames.front().energy * kEh2kJ, mismatch);
    else
        CurcumaLogger::info_fmt("ConfGen: energy decomposition sums to the total energy (difference {:.4f} kJ/mol)",
            mismatch);

    // ---- 4. rotamer states per torsion --------------------------------------------------------
    m_state_centres.assign(m_torsions.size(), {});
    for (std::size_t t = 0; t < m_torsions.size(); ++t) {
        std::vector<double> observed;
        observed.reserve(m_frames.size());
        for (const Frame& f : m_frames)
            observed.push_back(f.angles[t]);
        m_state_centres[t] = TorsionSpace::clusterStates(observed, m_state_tolerance);
    }
    for (Frame& f : m_frames) {
        f.states.assign(m_torsions.size(), -1);
        for (std::size_t t = 0; t < m_torsions.size(); ++t)
            f.states[t] = TorsionSpace::assignState(f.angles[t], m_state_centres[t]);
    }

    m_reference_energy = std::min_element(m_frames.begin(), m_frames.end(),
        [](const Frame& a, const Frame& b) { return a.energy < b.energy; })
                             ->energy;
    return true;
}

std::vector<ConfGen::Transition> ConfGen::matchedPairs() const
{
    // Key: (torsion, low state, high state). Canonical direction low -> high, so dE is unambiguous.
    std::map<std::tuple<int, int, int>, std::vector<std::pair<double, std::map<std::string, double>>>> collected;
    // Distinct structures per side: with only one structure in a state, N "pairs" are N repetitions
    // of the SAME measurement, not N independent ones. Tracked so the report can say so.
    std::map<std::tuple<int, int, int>, std::pair<std::set<int>, std::set<int>>> members;

    for (std::size_t a = 0; a < m_frames.size(); ++a) {
        for (std::size_t b = a + 1; b < m_frames.size(); ++b) {
            if (TorsionSpace::hammingDistance(m_frames[a].states, m_frames[b].states) != 1)
                continue;
            // locate the single differing torsion
            int t = -1;
            for (std::size_t idx = 0; idx < m_frames[a].states.size(); ++idx)
                if (m_frames[a].states[idx] != m_frames[b].states[idx]) {
                    t = static_cast<int>(idx);
                    break;
                }
            if (t < 0)
                continue;
            const int s_a = m_frames[a].states[t], s_b = m_frames[b].states[t];
            if (s_a < 0 || s_b < 0)
                continue;
            // orient the difference from the LOWER to the HIGHER state index
            const std::size_t low = (s_a < s_b) ? a : b;
            const std::size_t high = (s_a < s_b) ? b : a;
            std::map<std::string, double> d_terms;
            for (const auto& name : m_term_names)
                d_terms[name] = m_frames[high].terms.at(name) - m_frames[low].terms.at(name);
            const auto key = std::make_tuple(t, std::min(s_a, s_b), std::max(s_a, s_b));
            collected[key].emplace_back(m_frames[high].energy - m_frames[low].energy, std::move(d_terms));
            members[key].first.insert(static_cast<int>(low));
            members[key].second.insert(static_cast<int>(high));
        }
    }

    std::vector<Transition> transitions;
    for (const auto& [key, samples] : collected) {
        Transition tr;
        std::tie(tr.torsion, tr.state_from, tr.state_to) = key;
        tr.pairs = static_cast<int>(samples.size());
        if (tr.pairs < m_min_pairs)
            continue;
        tr.d_total_min = samples.front().first;
        tr.d_total_max = samples.front().first;
        for (const auto& [d_total, d_terms] : samples) {
            tr.d_total_mean += d_total;
            tr.d_total_min = std::min(tr.d_total_min, d_total);
            tr.d_total_max = std::max(tr.d_total_max, d_total);
            for (const auto& [name, value] : d_terms)
                tr.d_terms_mean[name] += value;
        }
        tr.d_total_mean /= static_cast<double>(tr.pairs);
        for (auto& [name, value] : tr.d_terms_mean)
            value /= static_cast<double>(tr.pairs);
        if (tr.pairs > 1) {
            double sq = 0.0;
            for (const auto& [d_total, d_terms] : samples)
                sq += (d_total - tr.d_total_mean) * (d_total - tr.d_total_mean);
            tr.d_total_sd = std::sqrt(sq / static_cast<double>(tr.pairs - 1));
        }
        tr.distinct_from = static_cast<int>(members.at(key).first.size());
        tr.distinct_to = static_cast<int>(members.at(key).second.size());
        transitions.push_back(std::move(tr));
    }
    // strongest effect first
    std::sort(transitions.begin(), transitions.end(), [](const Transition& x, const Transition& y) {
        return std::abs(x.d_total_mean) > std::abs(y.d_total_mean);
    });
    return transitions;
}

std::vector<int> ConfGen::informativeTorsions() const
{
    std::vector<int> out;
    for (std::size_t t = 0; t < m_torsions.size(); ++t) {
        std::set<int> populated;
        for (const Frame& f : m_frames)
            if (f.states[t] >= 0)
                populated.insert(f.states[t]);
        if (populated.size() >= 2)
            out.push_back(static_cast<int>(t));
    }
    return out;
}

std::vector<ConfGen::Coupling> ConfGen::doubleMutantCycles() const
{
    // A cycle needs four members that agree in every torsion except two, and cover all four
    // combinations of those two state changes. Index by "everything else" so the rectangles can be
    // found without an O(N^4) scan: for each torsion pair, group the frames by the remaining states.
    const std::vector<int> informative = informativeTorsions();
    std::map<std::tuple<int, int, int, int, int, int>, std::vector<double>> j_total;
    std::map<std::tuple<int, int, int, int, int, int>, std::vector<std::map<std::string, double>>> j_terms;

    for (std::size_t ia = 0; ia < informative.size(); ++ia) {
        for (std::size_t ib = ia + 1; ib < informative.size(); ++ib) {
            const int ta = informative[ia], tb = informative[ib];
            // group frames by the state of all OTHER torsions
            std::map<std::vector<int>, std::map<std::pair<int, int>, int>> groups; // context -> (sa,sb) -> frame
            for (std::size_t f = 0; f < m_frames.size(); ++f) {
                std::vector<int> context;
                context.reserve(m_frames[f].states.size());
                for (std::size_t t = 0; t < m_frames[f].states.size(); ++t)
                    if (static_cast<int>(t) != ta && static_cast<int>(t) != tb)
                        context.push_back(m_frames[f].states[t]);
                const std::pair<int, int> key{ m_frames[f].states[ta], m_frames[f].states[tb] };
                // keep the lowest-energy representative if a combination occurs twice
                auto& slot = groups[context];
                auto it = slot.find(key);
                if (it == slot.end() || m_frames[f].energy < m_frames[it->second].energy)
                    slot[key] = static_cast<int>(f);
            }

            for (const auto& [context, combos] : groups) {
                if (combos.size() < 4)
                    continue;
                // every rectangle (a_from,a_to) x (b_from,b_to) present in this context
                std::set<int> states_a, states_b;
                for (const auto& [key, frame] : combos) {
                    states_a.insert(key.first);
                    states_b.insert(key.second);
                }
                for (auto a_it = states_a.begin(); a_it != states_a.end(); ++a_it)
                    for (auto a2 = std::next(a_it); a2 != states_a.end(); ++a2)
                        for (auto b_it = states_b.begin(); b_it != states_b.end(); ++b_it)
                            for (auto b2 = std::next(b_it); b2 != states_b.end(); ++b2) {
                                const std::pair<int, int> ac{ *a_it, *b_it }, bc{ *a2, *b_it },
                                    ad{ *a_it, *b2 }, bd{ *a2, *b2 };
                                if (!combos.count(ac) || !combos.count(bc) || !combos.count(ad) || !combos.count(bd))
                                    continue;
                                const Frame& f_ac = m_frames[combos.at(ac)];
                                const Frame& f_bc = m_frames[combos.at(bc)];
                                const Frame& f_ad = m_frames[combos.at(ad)];
                                const Frame& f_bd = m_frames[combos.at(bd)];
                                const double j = (f_bd.energy - f_ad.energy) - (f_bc.energy - f_ac.energy);
                                const auto key = std::make_tuple(ta, tb, *a_it, *a2, *b_it, *b2);
                                j_total[key].push_back(j);
                                std::map<std::string, double> terms;
                                for (const auto& name : m_term_names)
                                    terms[name] = (f_bd.terms.at(name) - f_ad.terms.at(name))
                                        - (f_bc.terms.at(name) - f_ac.terms.at(name));
                                j_terms[key].push_back(std::move(terms));
                            }
            }
        }
    }

    std::vector<Coupling> couplings;
    for (const auto& [key, values] : j_total) {
        Coupling c;
        std::tie(c.torsion_a, c.torsion_b, c.a_from, c.a_to, c.b_from, c.b_to) = key;
        c.cycles = static_cast<int>(values.size());
        c.j_mean = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(c.cycles);
        if (c.cycles > 1) {
            double sq = 0.0;
            for (double v : values)
                sq += (v - c.j_mean) * (v - c.j_mean);
            c.j_sd = std::sqrt(sq / static_cast<double>(c.cycles - 1));
        }
        for (const auto& terms : j_terms.at(key))
            for (const auto& [name, value] : terms)
                c.j_terms_mean[name] += value / static_cast<double>(c.cycles);
        couplings.push_back(std::move(c));
    }
    std::sort(couplings.begin(), couplings.end(), [](const Coupling& x, const Coupling& y) {
        return std::abs(x.j_mean) > std::abs(y.j_mean);
    });
    return couplings;
}

ConfGen::ModelFit ConfGen::fitModel(int level) const
{
    ModelFit fit;
    fit.name = (level == 0) ? "constant" : (level == 1) ? "additive (marginals)" : "additive + pair couplings";

    const std::vector<int> informative = informativeTorsions();
    // Column layout: [intercept][per torsion: states 1..k-1][per torsion pair: (a,b) with a,b >= 1]
    // State 0 is the reference level, so the columns stay linearly independent by construction.
    struct Column {
        int torsion_a = -1, state_a = -1, torsion_b = -1, state_b = -1;
    };
    std::vector<Column> columns;
    if (level >= 1)
        for (int t : informative)
            for (std::size_t s = 1; s < m_state_centres[t].size(); ++s)
                columns.push_back({ t, static_cast<int>(s), -1, -1 });
    if (level >= 2)
        for (std::size_t ia = 0; ia < informative.size(); ++ia)
            for (std::size_t ib = ia + 1; ib < informative.size(); ++ib)
                for (std::size_t sa = 1; sa < m_state_centres[informative[ia]].size(); ++sa)
                    for (std::size_t sb = 1; sb < m_state_centres[informative[ib]].size(); ++sb)
                        columns.push_back({ informative[ia], static_cast<int>(sa),
                            informative[ib], static_cast<int>(sb) });

    const int n = static_cast<int>(m_frames.size());
    const int p = static_cast<int>(columns.size()) + 1; // + intercept
    fit.columns = p;

    auto buildRow = [&](const Frame& f, Eigen::VectorXd& row) {
        row.setZero(p);
        row(0) = 1.0;
        for (int c = 0; c < static_cast<int>(columns.size()); ++c) {
            const Column& col = columns[c];
            const bool a_on = (f.states[col.torsion_a] == col.state_a);
            const bool on = (col.torsion_b < 0) ? a_on : (a_on && f.states[col.torsion_b] == col.state_b);
            row(c + 1) = on ? 1.0 : 0.0;
        }
    };

    Eigen::MatrixXd X(n, p);
    Eigen::VectorXd y(n);
    Eigen::VectorXd row(p);
    for (int i = 0; i < n; ++i) {
        buildRow(m_frames[i], row);
        X.row(i) = row.transpose();
        y(i) = (m_frames[i].energy - m_reference_energy) * kEh2kJ;
    }

    // Complete orthogonal decomposition: the design matrix IS rank deficient whenever a state
    // combination never occurs in the ensemble (an empty column) -- COD gives the minimum-norm
    // solution instead of exploding, and its rank tells us how much the ensemble can actually resolve.
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
    cod.setThreshold(1e-10);
    fit.rank = static_cast<int>(cod.rank());
    const Eigen::VectorXd beta_in = cod.solve(y);
    fit.rmse_in = std::sqrt((X * beta_in - y).squaredNorm() / static_cast<double>(n));

    // k-fold cross-validation, deterministic split (fold = index % k). The ensemble is usually
    // energy-sorted, so modulo assignment stratifies over the energy range automatically.
    const int folds = std::max(2, m_cv_folds);
    std::vector<double> residuals;
    residuals.reserve(n);
    for (int k = 0; k < folds; ++k) {
        std::vector<int> train, test;
        for (int i = 0; i < n; ++i)
            ((i % folds) == k ? test : train).push_back(i);
        if (train.size() < static_cast<std::size_t>(p) || test.empty())
            continue;
        Eigen::MatrixXd Xtr(train.size(), p);
        Eigen::VectorXd ytr(train.size());
        for (std::size_t r = 0; r < train.size(); ++r) {
            Xtr.row(r) = X.row(train[r]);
            ytr(r) = y(train[r]);
        }
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod_tr(Xtr);
        cod_tr.setThreshold(1e-10);
        const Eigen::VectorXd beta = cod_tr.solve(ytr);
        for (int i : test)
            residuals.push_back(X.row(i) * beta - y(i));
    }
    if (!residuals.empty()) {
        double sq = 0.0;
        for (double r : residuals)
            sq += r * r;
        fit.rmse_cv = std::sqrt(sq / static_cast<double>(residuals.size()));
        std::vector<double> abs_res;
        abs_res.reserve(residuals.size());
        for (double r : residuals)
            abs_res.push_back(std::abs(r));
        std::sort(abs_res.begin(), abs_res.end());
        fit.mae_cv = abs_res[abs_res.size() / 2];
    }
    return fit;
}

void ConfGen::writeCouplingTable(const std::string& path, const std::vector<Coupling>& couplings) const
{
    std::ofstream out(path);
    out << "# torsion_a,bond_a,a_from,a_to,torsion_b,bond_b,b_from,b_to,cycles,J_mean_kJ,J_sd_kJ";
    for (const auto& name : m_term_names)
        out << ",J_" << name << "_kJ";
    out << "\n";
    for (const auto& c : couplings) {
        out << c.torsion_a << "," << m_torsions[c.torsion_a].label(m_frames.front().molecule) << ","
            << c.a_from << "," << c.a_to << "," << c.torsion_b << ","
            << m_torsions[c.torsion_b].label(m_frames.front().molecule) << "," << c.b_from << ","
            << c.b_to << "," << c.cycles << "," << fmt::format("{:.3f}", c.j_mean * kEh2kJ) << ","
            << fmt::format("{:.3f}", c.j_sd * kEh2kJ);
        for (const auto& name : m_term_names)
            out << "," << fmt::format("{:.3f}", c.j_terms_mean.at(name) * kEh2kJ);
        out << "\n";
    }
}

void ConfGen::reportCouplings(const std::vector<Coupling>& couplings) const
{
    CurcumaLogger::header("=== ConfGen: torsion-torsion couplings (double-mutant cycles) ===");
    if (couplings.empty()) {
        CurcumaLogger::warn("ConfGen: no double-mutant cycle in the ensemble -- no four structures form a "
                            "rectangle in state space (two torsions changed, everything else identical). "
                            "Couplings can then only be estimated by the regression below, not measured.");
        return;
    }
    int shown = 0;
    for (const auto& c : couplings) {
        const double j_kj = c.j_mean * kEh2kJ;
        if (std::abs(j_kj) < m_report_threshold || shown >= 10)
            continue;
        std::vector<std::pair<std::string, double>> parts;
        for (const auto& [name, value] : c.j_terms_mean)
            if (std::abs(value * kEh2kJ) >= 0.5)
                parts.emplace_back(name, value * kEh2kJ);
        std::sort(parts.begin(), parts.end(), [](const auto& x, const auto& y) {
            return std::abs(x.second) > std::abs(y.second);
        });
        std::string breakdown;
        for (const auto& [name, value] : parts)
            breakdown += fmt::format("{}{:+.1f} {}", breakdown.empty() ? "" : ", ", value, name);

        CurcumaLogger::result_fmt("  {:<10} ({}->{}) x {:<10} ({}->{}): J = {:+.2f} kJ/mol"
                                  "{}  [{}]  from {} cycle(s)",
            m_torsions[c.torsion_a].label(m_frames.front().molecule), c.a_from, c.a_to,
            m_torsions[c.torsion_b].label(m_frames.front().molecule), c.b_from, c.b_to,
            j_kj, c.cycles > 1 ? fmt::format(" +- {:.2f}", c.j_sd * kEh2kJ)
                               : std::string(" (single cycle, no error estimate)"),
            breakdown.empty() ? "no single term above 0.5 kJ/mol" : breakdown, c.cycles);
        shown++;
    }
    CurcumaLogger::info("J is the non-additivity: J = 0 means the two state changes are independent and "
                        "their effects simply add. |J| of the size of the effects themselves means the "
                        "torsions must be chosen together, not one at a time.");
}

void ConfGen::reportModelComparison(const std::vector<ModelFit>& fits) const
{
    CurcumaLogger::header("=== ConfGen: can the ensemble energies be modelled at all? ===");
    CurcumaLogger::result_fmt("  {:<28} {:>7} {:>7} {:>12} {:>12} {:>10} {:>12}",
        "model", "params", "rank", "RMSE_cv", "medAE_cv", "R2_cv", "RMSE_in");
    const double null_rmse = fits.empty() ? 0.0 : fits.front().rmse_cv;
    for (const auto& f : fits) {
        // R2_cv against the constant model: the fraction of the energy variation that the model
        // predicts on data it has not seen. This is the number that decides whether the torsion-state
        // picture describes this ensemble at all.
        const double r2 = (null_rmse > 0.0) ? 1.0 - (f.rmse_cv * f.rmse_cv) / (null_rmse * null_rmse) : 0.0;
        CurcumaLogger::result_fmt("  {:<28} {:>7} {:>7} {:>9.2f} kJ {:>9.2f} kJ {:>9.0f} % {:>9.2f} kJ",
            f.name, f.columns, f.rank, f.rmse_cv, f.mae_cv, r2 * 100.0, f.rmse_in);
    }
    CurcumaLogger::info("RMSE_cv/medAE_cv are out-of-sample (k-fold). RMSE_in always improves with more "
                        "parameters and is shown for reference only.");

    if (fits.size() < 3)
        return;
    const double null_err = fits[0].rmse_cv, add_err = fits[1].rmse_cv, pair_err = fits[2].rmse_cv;
    const double add_gain = (null_err > 0.0) ? (1.0 - add_err / null_err) * 100.0 : 0.0;
    const double pair_gain = (add_err > 0.0) ? (1.0 - pair_err / add_err) * 100.0 : 0.0;

    const double r2_add = (null_err > 0.0) ? 1.0 - (add_err * add_err) / (null_err * null_err) : 0.0;
    const double r2_pair = (null_err > 0.0) ? 1.0 - (pair_err * pair_err) / (null_err * null_err) : 0.0;

    // Verdict. Two guards against reading noise as signal:
    //   - a few percent of RMSE is within the fold-to-fold scatter of a k-fold CV, so a minimum
    //     improvement is required before anything is called informative;
    //   - RMSE reacts to outliers. If the MEDIAN error moves the other way, the "improvement" is a
    //     handful of structures, not a better description of the ensemble.
    constexpr double kMinRelGain = 0.05; // 5 % of the RMSE
    const bool add_useful = (add_err < null_err * (1.0 - kMinRelGain)) && (fits[1].mae_cv <= fits[0].mae_cv);
    const bool pair_useful = (pair_err < add_err * (1.0 - kMinRelGain)) && (fits[2].mae_cv <= fits[1].mae_cv);

    CurcumaLogger::result_fmt("ConfGen: torsion states explain {:.0f}% of the energy variation out of "
                              "sample, adding pair couplings takes it to {:.0f}%",
        r2_add * 100.0, r2_pair * 100.0);

    if (!add_useful)
        CurcumaLogger::warn_fmt("ConfGen: the additive model is NOT a useful description of this ensemble "
                                "(RMSE {:.2f} -> {:.2f} kJ/mol = {:.0f}%, median error {:.2f} -> {:.2f}). "
                                "Most of the energy differences come from degrees of freedom the torsion "
                                "states do not capture, so proposals ranked by these marginals would be "
                                "little better than guessing.",
            null_err, add_err, add_gain, fits[0].mae_cv, fits[1].mae_cv);
    else
        CurcumaLogger::success_fmt("ConfGen: the additive model removes {:.0f}% of the energy scatter "
                                   "({:.2f} -> {:.2f} kJ/mol out of sample)",
            add_gain, null_err, add_err);

    if (!pair_useful)
        CurcumaLogger::warn_fmt("ConfGen: pair couplings give NO reliable improvement (RMSE {:.2f} -> "
                                "{:.2f} kJ/mol = {:.0f}%, median error {:.2f} -> {:.2f}; {} of {} columns "
                                "resolvable). Either the ensemble is too small to fit them, or the "
                                "non-additivity is not of pairwise form.",
            add_err, pair_err, pair_gain, fits[1].mae_cv, fits[2].mae_cv, fits[2].rank, fits[2].columns);
    else
        CurcumaLogger::success_fmt("ConfGen: pair couplings improve the out-of-sample prediction by a "
                                   "further {:.0f}% ({:.2f} -> {:.2f} kJ/mol) -- the state vector must be "
                                   "optimised as a whole, not torsion by torsion.",
            pair_gain, add_err, pair_err);

    if (!add_useful && !pair_useful)
        CurcumaLogger::warn("ConfGen: on this ensemble the torsion-state model does not support "
                            "generating proposals. More conformers, or torsions that the search actually "
                            "varied, are needed before recombination can be expected to work.");
}

namespace {
/// Best-fit (Kabsch) RMSD in Angstrom between two geometries in the same atom order.
/// Twin of the lambda in ConfSearch::PermRMSD (confsearch.cpp) -- kept local and permutation-free
/// because ConfGen compares structures of one molecule in a fixed order.
double bestFitRMSD(const Geometry& reference, const Geometry& target)
{
    const int n = reference.rows();
    if (n == 0 || target.rows() != n)
        return std::numeric_limits<double>::infinity();
    Geometry r = reference, t = target;
    Eigen::Vector3d cr = Eigen::Vector3d::Zero(), ct = Eigen::Vector3d::Zero();
    for (int i = 0; i < n; ++i) {
        cr += r.row(i).transpose();
        ct += t.row(i).transpose();
    }
    cr /= n;
    ct /= n;
    for (int i = 0; i < n; ++i) {
        r.row(i) -= cr.transpose();
        t.row(i) -= ct.transpose();
    }
    const Eigen::Matrix3d rot = RMSDFunctions::BestFitRotation(r, t, 1);
    return RMSDFunctions::getRMSD(r, RMSDFunctions::applyRotation(t, rot));
}
}

std::vector<std::vector<double>> ConfGen::additiveCoefficients() const
{
    // Same design as fitModel(level=1) but returning the coefficients instead of a score. Used ONLY
    // to order the proposals -- the model explains ~15 % of the energy variation, which is enough to
    // decide what to try first and far too little to decide what is good.
    const std::vector<int> informative = informativeTorsions();
    std::vector<std::pair<int, int>> columns; // (torsion, state)
    for (int t : informative)
        for (std::size_t st = 1; st < m_state_centres[t].size(); ++st)
            columns.emplace_back(t, static_cast<int>(st));

    std::vector<std::vector<double>> coefficients(m_torsions.size());
    for (std::size_t t = 0; t < m_torsions.size(); ++t)
        coefficients[t].assign(m_state_centres[t].size(), 0.0);
    if (columns.empty())
        return coefficients;

    const int n = static_cast<int>(m_frames.size());
    const int p = static_cast<int>(columns.size()) + 1;
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n, p);
    Eigen::VectorXd y(n);
    for (int i = 0; i < n; ++i) {
        X(i, 0) = 1.0;
        for (int c = 0; c < static_cast<int>(columns.size()); ++c)
            X(i, c + 1) = (m_frames[i].states[columns[c].first] == columns[c].second) ? 1.0 : 0.0;
        y(i) = (m_frames[i].energy - m_reference_energy) * kEh2kJ;
    }
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
    cod.setThreshold(1e-10);
    const Eigen::VectorXd beta = cod.solve(y);
    for (int c = 0; c < static_cast<int>(columns.size()); ++c)
        coefficients[columns[c].first][columns[c].second] = beta(c + 1);
    return coefficients;
}

std::vector<std::pair<int, int>> ConfGen::topologyFingerprint(const Molecule& mol, double factor) const
{
    std::vector<std::pair<int, int>> bonds;
    const Geometry geom = mol.getGeometry();
    for (int i = 0; i < mol.AtomCount(); ++i)
        for (int j = i + 1; j < mol.AtomCount(); ++j) {
            const double cutoff = factor
                * (Elements::CovalentRadius[mol.Atom(i).first] + Elements::CovalentRadius[mol.Atom(j).first]);
            if ((Eigen::Vector3d(geom.row(i)) - Eigen::Vector3d(geom.row(j))).norm() <= cutoff)
                bonds.emplace_back(i, j);
        }
    return bonds; // already sorted by construction
}

bool ConfGen::hasClash(const Molecule& mol, double factor) const
{
    // Claude Generated (Jul 2026): the reference bond list decides which pairs are ALLOWED to be
    // close. Built with the explicit topology factor, like every other topology decision here --
    // Molecule::DistanceMatrix()'s 1.5 also exempts compressed 1-3 contacts, which are not bonds.
    const std::vector<std::pair<int, int>> ref = topologyFingerprint(m_frames.front().molecule, m_topology_factor);
    std::set<std::pair<int, int>> bonded(ref.begin(), ref.end());
    const Geometry geom = mol.getGeometry();
    for (int i = 0; i < mol.AtomCount(); ++i) {
        for (int j = i + 1; j < mol.AtomCount(); ++j) {
            if (bonded.count({ i, j }))
                continue; // bonded pairs are supposed to be close
            const double limit = factor
                * (Elements::CovalentRadius[mol.Atom(i).first] + Elements::CovalentRadius[mol.Atom(j).first]);
            if ((Eigen::Vector3d(geom.row(i)) - Eigen::Vector3d(geom.row(j))).norm() < limit)
                return true;
        }
    }
    return false;
}

// Claude Generated (Jul 2026, roadmap item P0): restrained build -- see the header for the why.
bool ConfGen::restrainedBuild(const Proposal& p, Molecule& driven) const
{
    if (!m_calculator)
        return false;

    // One restraint per torsion whose state differs from the template's.
    std::vector<Optimization::DihedralRestraint> restraints;
    for (std::size_t t = 0; t < m_torsions.size(); ++t) {
        if (p.states[t] == m_frames[p.template_frame].states[t] || p.states[t] < 0)
            continue;
        if (p.states[t] >= static_cast<int>(m_state_centres[t].size()))
            continue;
        Optimization::DihedralRestraint r;
        r.i = m_torsions[t].i;
        r.j = m_torsions[t].j;
        r.k = m_torsions[t].k;
        r.l = m_torsions[t].l;
        r.target = m_state_centres[t][p.states[t]] * M_PI / 180.0; // state centres are in degrees
        r.force = m_restraint_force;
        restraints.push_back(r);
    }
    if (restraints.empty())
        return false;

    json opt_config;
    opt_config["method"] = m_method;
    opt_config["threads"] = m_threads;
    opt_config["charge"] = m_charge;
    opt_config["spin"] = m_spin;
    opt_config["verbosity"] = 0;
    opt_config["write_trajectory"] = false;
    opt_config["max_iterations"] = m_restraint_max_iterations;
    opt_config["dihedral_restraints"] = Optimization::GeometryRestraints::toJson(restraints);

    // Start from the TEMPLATE geometry -- clash-free and with the correct topology. The restrained
    // optimisation itself performs the rotation.
    Molecule mol = m_frames[p.template_frame].molecule;
    auto result = Optimization::OptimizationDispatcher::optimizeStructure(
        &mol, Optimization::OptimizerType::LBFGSPP, m_calculator.get(), opt_config);
    if (result.final_molecule.AtomCount() == 0)
        return false; // hard failure; a non-converged but finite geometry is still usable below

    // Did the torsions actually arrive? A restraint competes with the force field, so a target that
    // is sterically impossible ends up somewhere else -- and then this is not the proposed state
    // vector at all and must not be reported as one. 30 degrees is the tolerance: it is well inside
    // the default state_tolerance of 40 degrees used to assign states.
    Optimization::GeometryRestraints check;
    for (const auto& r : restraints)
        check.add(r);
    const double worst = check.maxDeviationDegrees(result.final_molecule.getGeometry());
    if (worst > 30.0)
        return false;

    driven = result.final_molecule;
    return true;
}

std::vector<ConfGen::Proposal> ConfGen::generateProposals() const
{
    std::vector<Proposal> proposals;
    const std::vector<int> informative = informativeTorsions();
    if (informative.empty()) {
        CurcumaLogger::warn("ConfGen: no torsion has more than one populated state -- there is nothing "
                            "to recombine. The search never varied these torsions.");
        return proposals;
    }

    // Everything the ensemble already contains; proposals must be outside this set.
    std::set<std::vector<int>> known;
    for (const Frame& f : m_frames)
        known.insert(f.states);

    // Templates: the lowest-energy members (their geometry is the starting point).
    std::vector<int> order(m_frames.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [this](int a, int b) { return m_frames[a].energy < m_frames[b].energy; });
    const int n_templates = std::min<int>(std::max(1, m_proposal_templates), static_cast<int>(order.size()));

    const std::vector<std::vector<double>> coefficients = additiveCoefficients();
    std::set<std::vector<int>> proposed;

    // Enumerate mutations up to the requested Hamming depth around each template.
    std::function<void(int, std::vector<int>&, int, int)> mutate =
        [&](int template_frame, std::vector<int>& states, int start, int depth) {
            if (depth > 0 && !known.count(states) && !proposed.count(states)) {
                Proposal p;
                p.states = states;
                p.template_frame = template_frame;
                p.distance = depth;
                for (std::size_t t = 0; t < states.size(); ++t)
                    if (states[t] >= 0 && states[t] < static_cast<int>(coefficients[t].size()))
                        p.predicted += coefficients[t][states[t]];
                proposed.insert(states);
                proposals.push_back(std::move(p));
            }
            if (depth >= m_proposal_depth)
                return;
            for (std::size_t idx = start; idx < informative.size(); ++idx) {
                const int t = informative[idx];
                const int original = states[t];
                for (int st = 0; st < static_cast<int>(m_state_centres[t].size()); ++st) {
                    if (st == original)
                        continue;
                    states[t] = st;
                    mutate(template_frame, states, idx + 1, depth + 1);
                }
                states[t] = original;
            }
        };

    for (int i = 0; i < n_templates; ++i) {
        std::vector<int> states = m_frames[order[i]].states;
        mutate(order[i], states, 0, 0);
    }

    // Order by the model estimate and keep the requested number.
    std::sort(proposals.begin(), proposals.end(),
        [](const Proposal& a, const Proposal& b) { return a.predicted < b.predicted; });
    if (static_cast<int>(proposals.size()) > m_max_proposals)
        proposals.resize(m_max_proposals);
    return proposals;
}

void ConfGen::optimiseProposals(std::vector<Proposal>& proposals) const
{
    json calc_config;
    calc_config["method"] = m_method;
    calc_config["threads"] = m_threads;
    calc_config["charge"] = m_charge;
    calc_config["spin"] = m_spin;

    json opt_config = calc_config;
    opt_config["verbosity"] = 0;
    opt_config["write_trajectory"] = false;

    if (!m_calculator)
        return; // analyseEnsemble must have run first
    EnergyCalculator& calculator = *m_calculator;
    const std::vector<std::pair<int, int>> reference_bonds
        = topologyFingerprint(m_frames.front().molecule, m_topology_factor);

    for (Proposal& p : proposals) {
        if (p.geometry.AtomCount() == 0)
            continue; // rejected by the clash filter, never built
        Molecule mol = p.geometry;
        auto result = Optimization::OptimizationDispatcher::optimizeStructure(
            &mol, Optimization::OptimizerType::LBFGSPP, &calculator, opt_config);
        if (!result.success)
            continue;
        p.optimised = true;
        p.geometry = result.final_molecule;
        p.energy = result.final_energy;
        p.states_after = TorsionSpace::stateVector(p.geometry.getGeometry(), m_torsions, m_state_centres);

        // MANDATORY topology check. A built structure that brings two atoms close makes the force
        // field derive a NEW BOND; the optimisation then relaxes into a different molecule whose
        // energy is far below any conformer (measured: 2649 kJ/mol "below the minimum", with 3-6
        // changed bonds). Without this check the generator reports reactions as brilliant conformers.
        p.topology_ok = (topologyFingerprint(p.geometry, m_topology_factor) == reference_bonds);

        // Novelty: closest input structure by best-fit RMSD. The force field, not the model, decides.
        p.min_rmsd_to_ensemble = std::numeric_limits<double>::infinity();
        for (const Frame& f : m_frames)
            p.min_rmsd_to_ensemble = std::min(p.min_rmsd_to_ensemble,
                bestFitRMSD(f.molecule.getGeometry(), p.geometry.getGeometry()));
        p.is_new = p.topology_ok && p.min_rmsd_to_ensemble > m_new_rmsd;
    }
}

double ConfGen::referenceEnergyOptimised(double& worst_gain_kJ) const
{
    json calc_config;
    calc_config["method"] = m_method;
    calc_config["threads"] = m_threads;
    calc_config["charge"] = m_charge;
    calc_config["spin"] = m_spin;
    json opt_config = calc_config;
    opt_config["verbosity"] = 0;
    opt_config["write_trajectory"] = false;

    if (!m_calculator)
        return std::numeric_limits<double>::infinity();
    EnergyCalculator& calculator = *m_calculator;

    std::vector<int> order(m_frames.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [this](int a, int b) { return m_frames[a].energy < m_frames[b].energy; });
    const int n = std::min<int>(std::max(1, m_proposal_templates), static_cast<int>(order.size()));

    double best = std::numeric_limits<double>::infinity();
    worst_gain_kJ = 0.0;
    for (int i = 0; i < n; ++i) {
        Molecule mol = m_frames[order[i]].molecule;
        auto result = Optimization::OptimizationDispatcher::optimizeStructure(
            &mol, Optimization::OptimizerType::LBFGSPP, &calculator, opt_config);
        if (!result.success)
            continue;
        best = std::min(best, result.final_energy);
        worst_gain_kJ = std::max(worst_gain_kJ, (m_frames[order[i]].energy - result.final_energy) * kEh2kJ);
    }
    return best;
}

void ConfGen::reportProposals(const std::vector<Proposal>& proposals, const std::string& base) const
{
    CurcumaLogger::header("=== ConfGen: proposals (recombined state vectors) ===");
    if (proposals.empty()) {
        CurcumaLogger::warn("ConfGen: no proposal could be generated -- every reachable state "
                            "combination is already in the ensemble, or no torsion has a second state.");
        return;
    }

    // Like-for-like energy reference: the templates optimised with exactly the settings the proposals
    // were optimised with.
    double template_gain_kJ = 0.0;
    const double reference_opt = referenceEnergyOptimised(template_gain_kJ);

    int built = 0, optimised = 0, kept_states = 0, collapsed = 0, novel = 0, reacted = 0;
    double best_new = std::numeric_limits<double>::infinity();
    for (const Proposal& p : proposals) {
        if (p.geometry.AtomCount() > 0)
            built++;
        if (!p.optimised)
            continue;
        optimised++;
        if (!p.topology_ok) {
            reacted++;
            continue;
        }
        if (p.states_after == p.states)
            kept_states++;
        if (p.states_after == m_frames[p.template_frame].states)
            collapsed++;
        if (p.is_new) {
            novel++;
            best_new = std::min(best_new, p.energy);
        }
    }

    CurcumaLogger::result_fmt("ConfGen: {} state vector(s) proposed, {} built (rest rejected by the "
                              "clash filter), {} optimised successfully",
        static_cast<int>(proposals.size()), built, optimised);
    CurcumaLogger::result_fmt("ConfGen: {} optimised proposal(s) changed their bond topology and are NOT "
                              "conformers (a bond formed or broke -- rejected)",
        reacted);
    CurcumaLogger::result_fmt("ConfGen: of the {} chemically valid ones, {} kept their target states, {} "
                              "fell back into their template's state vector",
        optimised - reacted, kept_states, collapsed);
    CurcumaLogger::result_fmt("ConfGen: {} of {} chemically valid proposals are NEW conformers "
                              "(best-fit RMSD > {:.2f} A to every input structure)",
        novel, optimised - reacted, m_new_rmsd);

    if (template_gain_kJ > 5.0)
        CurcumaLogger::warn_fmt("ConfGen: re-optimising the input structures lowers them by up to {:.1f} "
                                "kJ/mol -- the ensemble is NOT at the minimum of '{}' (different method or "
                                "convergence). Energies below are therefore compared against re-optimised "
                                "templates, not against the file's own values.",
            template_gain_kJ, m_method);

    if (novel > 0) {
        const double delta = (best_new - reference_opt) * kEh2kJ;
        if (delta < 0.0)
            CurcumaLogger::success_fmt("ConfGen: the best new conformer is {:.2f} kJ/mol BELOW the "
                                       "re-optimised ensemble minimum ({:.6f} Eh) -- the search had "
                                       "missed it",
                -delta, best_new);
        else
            CurcumaLogger::result_fmt("ConfGen: the best new conformer sits {:.2f} kJ/mol above the "
                                      "re-optimised ensemble minimum",
                delta);
        CurcumaLogger::result_fmt("ConfGen: new structures written to {}.proposals.new.xyz", base);
    } else if (optimised > 0) {
        CurcumaLogger::warn("ConfGen: every optimised proposal fell back onto a structure that is "
                            "already in the ensemble. Recombination found nothing new here -- either "
                            "the ensemble is already complete in this torsion space, or the built "
                            "geometries relax straight back (a restrained pre-optimisation would be "
                            "the next lever).");
    }

    // How good was the ordering heuristic? Correlation between predicted and optimised energy.
    std::vector<std::pair<double, double>> pairs;
    for (const Proposal& p : proposals)
        if (p.optimised && p.topology_ok)
            pairs.emplace_back(p.predicted, (p.energy - m_reference_energy) * kEh2kJ);
    if (pairs.size() >= 3) {
        double mx = 0, my = 0;
        for (const auto& [x, y] : pairs) { mx += x; my += y; }
        mx /= pairs.size(); my /= pairs.size();
        double sxy = 0, sxx = 0, syy = 0;
        for (const auto& [x, y] : pairs) {
            sxy += (x - mx) * (y - my);
            sxx += (x - mx) * (x - mx);
            syy += (y - my) * (y - my);
        }
        const double r = (sxx > 0 && syy > 0) ? sxy / std::sqrt(sxx * syy) : 0.0;
        CurcumaLogger::result_fmt("ConfGen: ordering heuristic vs. optimised energy: correlation r = {:.2f} "
                                  "over {} proposals (this is the only claim the weak model makes)",
            r, static_cast<int>(pairs.size()));
    }
}

void ConfGen::writeFrameTable(const std::string& path) const
{
    std::ofstream out(path);
    out << "# structure";
    for (std::size_t t = 0; t < m_torsions.size(); ++t)
        out << ",angle_" << m_torsions[t].label(m_frames.front().molecule)
            << ",state_" << m_torsions[t].label(m_frames.front().molecule);
    out << ",energy_Eh,rel_energy_kJ";
    for (const auto& name : m_term_names)
        out << "," << name << "_Eh";
    out << "\n";
    for (std::size_t f = 0; f < m_frames.size(); ++f) {
        out << f + 1;
        for (std::size_t t = 0; t < m_torsions.size(); ++t)
            out << "," << fmt::format("{:.2f}", m_frames[f].angles[t]) << "," << m_frames[f].states[t];
        out << "," << fmt::format("{:.8f}", m_frames[f].energy)
            << "," << fmt::format("{:.3f}", (m_frames[f].energy - m_reference_energy) * kEh2kJ);
        for (const auto& name : m_term_names)
            out << "," << fmt::format("{:.8f}", m_frames[f].terms.at(name));
        out << "\n";
    }
}

void ConfGen::writeStateStatistics(const std::string& path) const
{
    std::ofstream out(path);
    out << "# torsion,bond,state,centre_deg,population,rel_energy_min_kJ,rel_energy_boltzmann_kJ\n";
    const double kT = kBoltzmannEh * m_temperature;
    for (std::size_t t = 0; t < m_torsions.size(); ++t) {
        for (std::size_t s = 0; s < m_state_centres[t].size(); ++s) {
            int population = 0;
            double e_min = std::numeric_limits<double>::infinity();
            double weight_sum = 0.0, weighted_energy = 0.0;
            for (const Frame& f : m_frames) {
                if (f.states[t] != static_cast<int>(s))
                    continue;
                population++;
                const double rel = f.energy - m_reference_energy;
                e_min = std::min(e_min, rel);
                const double w = std::exp(-rel / kT);
                weight_sum += w;
                weighted_energy += w * rel;
            }
            const double boltz = (weight_sum > 0.0) ? weighted_energy / weight_sum : 0.0;
            out << t << "," << m_torsions[t].label(m_frames.front().molecule) << "," << s << ","
                << fmt::format("{:.2f}", m_state_centres[t][s]) << "," << population << ","
                << (population > 0 ? fmt::format("{:.3f}", e_min * kEh2kJ) : "NA") << ","
                << (population > 0 ? fmt::format("{:.3f}", boltz * kEh2kJ) : "NA") << "\n";
        }
    }
}

void ConfGen::writeTransitionTable(const std::string& path, const std::vector<Transition>& transitions) const
{
    std::ofstream out(path);
    out << "# torsion,bond,state_from,centre_from_deg,state_to,centre_to_deg,pairs,"
           "structures_from,structures_to,dE_total_mean_kJ,dE_total_sd_kJ,dE_total_min_kJ,"
           "dE_total_max_kJ,range_kJ";
    for (const auto& name : m_term_names)
        out << ",dE_" << name << "_kJ";
    out << "\n";
    for (const auto& tr : transitions) {
        out << tr.torsion << "," << m_torsions[tr.torsion].label(m_frames.front().molecule) << ","
            << tr.state_from << "," << fmt::format("{:.2f}", m_state_centres[tr.torsion][tr.state_from]) << ","
            << tr.state_to << "," << fmt::format("{:.2f}", m_state_centres[tr.torsion][tr.state_to]) << ","
            << tr.pairs << "," << tr.distinct_from << "," << tr.distinct_to << ","
            << fmt::format("{:.3f}", tr.d_total_mean * kEh2kJ) << ","
            << fmt::format("{:.3f}", tr.d_total_sd * kEh2kJ) << ","
            << fmt::format("{:.3f}", tr.d_total_min * kEh2kJ) << ","
            << fmt::format("{:.3f}", tr.d_total_max * kEh2kJ) << ","
            << fmt::format("{:.3f}", tr.range() * kEh2kJ);
        for (const auto& name : m_term_names)
            out << "," << fmt::format("{:.3f}", tr.d_terms_mean.at(name) * kEh2kJ);
        out << "\n";
    }
}

void ConfGen::reportTransitions(const std::vector<Transition>& transitions) const
{
    CurcumaLogger::header("=== ConfGen: torsion states in the ensemble ===");
    for (std::size_t t = 0; t < m_torsions.size(); ++t) {
        std::string line = fmt::format("  {:>3} {:<10}", t, m_torsions[t].label(m_frames.front().molecule));
        for (std::size_t s = 0; s < m_state_centres[t].size(); ++s) {
            int population = 0;
            for (const Frame& f : m_frames)
                if (f.states[t] == static_cast<int>(s))
                    population++;
            line += fmt::format("  [{}] {:>7.1f} deg n={:<4}", s, m_state_centres[t][s], population);
        }
        if (m_state_centres[t].size() < 2)
            line += "   (single state -- no contrast, this torsion carries no information)";
        CurcumaLogger::result(line);
    }

    CurcumaLogger::header("=== ConfGen: measured contribution per state change (matched pairs) ===");
    if (transitions.empty()) {
        CurcumaLogger::warn("ConfGen: no matched pairs found. The ensemble contains no two structures "
                            "that differ in exactly one torsion state, so no single-torsion contribution "
                            "can be measured. More conformers, or a coarser -state_tolerance, would help.");
        return;
    }
    int shown = 0;
    for (const auto& tr : transitions) {
        const double d_total_kj = tr.d_total_mean * kEh2kJ;
        if (std::abs(d_total_kj) < m_report_threshold)
            continue;
        // Term breakdown, largest contributions first, only those above 0.5 kJ/mol.
        std::vector<std::pair<std::string, double>> parts;
        for (const auto& [name, value] : tr.d_terms_mean)
            if (std::abs(value * kEh2kJ) >= 0.5)
                parts.emplace_back(name, value * kEh2kJ);
        std::sort(parts.begin(), parts.end(), [](const auto& x, const auto& y) {
            return std::abs(x.second) > std::abs(y.second);
        });
        std::string breakdown;
        for (const auto& [name, value] : parts)
            breakdown += fmt::format("{}{:+.1f} {}", breakdown.empty() ? "" : ", ", value, name);
        if (breakdown.empty())
            breakdown = "no single term above 0.5 kJ/mol";

        // A state represented by a single structure gives one measurement, however many pairs it
        // forms -- flagged so the numbers are not over-read.
        const bool thin = (tr.distinct_from < 2 || tr.distinct_to < 2);
        CurcumaLogger::result_fmt("  {:<10} state {} ({:.0f} deg) -> state {} ({:.0f} deg): "
                                  "dE = {:+.2f} +- {:.2f} kJ/mol  [{}]  {} pair(s) over {}/{} structures{}",
            m_torsions[tr.torsion].label(m_frames.front().molecule),
            tr.state_from, m_state_centres[tr.torsion][tr.state_from],
            tr.state_to, m_state_centres[tr.torsion][tr.state_to],
            d_total_kj, tr.d_total_sd * kEh2kJ, breakdown, tr.pairs,
            tr.distinct_from, tr.distinct_to,
            thin ? "  <- one side has a single structure, not independent" : "");
        shown++;
    }
    if (shown == 0)
        CurcumaLogger::result_fmt("ConfGen: no transition exceeds the report threshold of {:.1f} kJ/mol "
                                  "(all {} of them are in the CSV)",
            m_report_threshold, static_cast<int>(transitions.size()));
    // The scientifically decisive number: how does the scatter compare to the effect itself?
    double sum_abs_mean = 0.0, sum_sd = 0.0;
    int counted = 0;
    for (const auto& tr : transitions) {
        if (tr.pairs < 2 || tr.distinct_from < 2 || tr.distinct_to < 2)
            continue;
        sum_abs_mean += std::abs(tr.d_total_mean);
        sum_sd += tr.d_total_sd;
        counted++;
    }
    if (counted > 0) {
        const double ratio = (sum_abs_mean > 0.0) ? sum_sd / sum_abs_mean : 0.0;
        CurcumaLogger::result_fmt("ConfGen: additivity check over {} well-sampled transition(s): mean |dE| "
                                  "= {:.2f} kJ/mol, mean scatter = {:.2f} kJ/mol (ratio {:.1f})",
            counted, sum_abs_mean / counted * kEh2kJ, sum_sd / counted * kEh2kJ, ratio);
        if (ratio > 1.0)
            CurcumaLogger::warn("ConfGen: the scatter EXCEEDS the effect -- on this molecule a single "
                                "torsion state has no environment-independent contribution. Recombining "
                                "from these marginals alone will be unreliable; the pairwise coupling "
                                "stage is required, not optional.");
        else
            CurcumaLogger::result("ConfGen: scatter below the effect -- the state contributions are "
                                 "approximately additive, so marginal-based recombination is meaningful.");
    }
}

void ConfGen::start()
{
    const std::string base = Basename();
    CurcumaLogger::header("=== ConfGen: ensemble analysis (matched pairs) ===");
    CurcumaLogger::result_fmt("ConfGen: method={}, state_tolerance={} deg, temperature={} K",
        m_method, m_state_tolerance, m_temperature);

    if (!analyseEnsemble())
        return;

    const std::vector<Transition> transitions = matchedPairs();

    writeFrameTable(outputPath(base + ".torsions.csv"));
    writeStateStatistics(outputPath(base + ".torsion_states.csv"));
    writeTransitionTable(outputPath(base + ".matched_pairs.csv"), transitions);

    reportTransitions(transitions);

    if (m_couplings) {
        const std::vector<Coupling> couplings = doubleMutantCycles();
        writeCouplingTable(outputPath(base + ".couplings.csv"), couplings);
        reportCouplings(couplings);

        if (m_cv_folds >= 2) {
            std::vector<ModelFit> fits;
            for (int level = 0; level <= 2; ++level)
                fits.push_back(fitModel(level));
            reportModelComparison(fits);
        }
    }

    if (m_generate) {
        std::vector<Proposal> proposals = generateProposals();

        // Build the geometries: start from the template and set every torsion whose state differs.
        const std::vector<std::pair<int, int>> reference_bonds
            = topologyFingerprint(m_frames.front().molecule, m_topology_factor);
        int clashes = 0, restrained = 0, restrained_ok = 0;
        for (Proposal& p : proposals) {
            Geometry geom = m_frames[p.template_frame].molecule.getGeometry();
            for (std::size_t t = 0; t < m_torsions.size(); ++t)
                if (p.states[t] != m_frames[p.template_frame].states[t] && p.states[t] >= 0)
                    geom = TorsionSpace::setDihedral(geom, m_torsions[t], m_state_centres[t][p.states[t]]);
            Molecule mol = m_frames[p.template_frame].molecule;
            mol.setGeometry(geom);
            mol.setName(fmt::format("proposal_from_{}_d{}", p.template_frame + 1, p.distance));
            // Two gates before a built structure is handed to the force field:
            //   (a) no non-bonded contact inside clash_factor x covalent radii, and
            //   (b) the bond topology must ALREADY match the reference.
            // (b) is the strict one: a rigidly set torsion that creates a bond can never become a
            // conformer anyway, so optimising it is wasted work -- and it is exactly the kind of
            // pathological input that has been observed to crash the GFN-FF parameter generation
            // (wild pointer in getGFNFFBondParameters, intermittent). Screening here keeps such
            // geometries out of the force field entirely instead of relying on it to cope.
            if (hasClash(mol, m_clash_factor)
                || topologyFingerprint(mol, m_topology_factor) != reference_bonds) {
                // Claude Generated (Jul 2026, P0): the rigid build failed -- try the restrained one.
                // Rigidly rotating a torsion on a fixed template moves a whole fragment into whatever
                // happens to be there; on a compact molecule that lost 72 % of all proposals. The
                // restrained build never creates that geometry in the first place: it starts from the
                // clash-free TEMPLATE and drives the torsions to their targets with a harmonic
                // restraint, so the rest of the molecule relaxes out of the way while they turn.
                if (!m_restrained_build) {
                    clashes++;
                    continue; // leave p.geometry empty -> counted as "not built"
                }
                restrained++;
                Molecule driven;
                if (!restrainedBuild(p, driven)) {
                    clashes++;
                    continue;
                }
                // The driven structure has to pass exactly the same gates -- the restraint is a way
                // to REACH the geometry, never a licence to skip the chemistry check.
                if (hasClash(driven, m_clash_factor)
                    || topologyFingerprint(driven, m_topology_factor) != reference_bonds) {
                    clashes++;
                    continue;
                }
                driven.setName(fmt::format("proposal_from_{}_d{}_restrained", p.template_frame + 1, p.distance));
                p.geometry = driven;
                p.restrained_build = true;
                restrained_ok++;
                continue;
            }
            p.geometry = mol;
        }
        if (restrained > 0)
            CurcumaLogger::result_fmt("ConfGen: rigid build failed for {} proposal(s) -- restrained build "
                                      "recovered {} of them (k = {} Eh/rad^2)",
                restrained, restrained_ok, m_restraint_force);
        if (clashes > 0)
            CurcumaLogger::result_fmt("ConfGen: {} built structure(s) rejected before optimisation "
                                      "(clash or changed bond topology)", clashes);

        // Write what was built, then optimise and write the survivors.
        bool first = true;
        for (const Proposal& p : proposals) {
            if (p.geometry.AtomCount() == 0)
                continue;
            if (first) { p.geometry.writeXYZFile(outputPath(base + ".proposals.xyz")); first = false; }
            else          p.geometry.appendXYZFile(outputPath(base + ".proposals.xyz"));
        }

        optimiseProposals(proposals);

        bool first_opt = true, first_new = true;
        for (const Proposal& p : proposals) {
            if (!p.optimised)
                continue;
            Molecule out = p.geometry;
            out.setEnergy(p.energy);
            if (first_opt) { out.writeXYZFile(outputPath(base + ".proposals.opt.xyz")); first_opt = false; }
            else             out.appendXYZFile(outputPath(base + ".proposals.opt.xyz"));
            if (p.is_new) {
                if (first_new) { out.writeXYZFile(outputPath(base + ".proposals.new.xyz")); first_new = false; }
                else             out.appendXYZFile(outputPath(base + ".proposals.new.xyz"));
            }
        }
        reportProposals(proposals, base);
    }

    CurcumaLogger::success_fmt("ConfGen: wrote {}.torsions.csv (per structure), {}.torsion_states.csv "
                               "(states + populations), {}.matched_pairs.csv ({} transition(s))",
        base, base, base, static_cast<int>(transitions.size()));
}
