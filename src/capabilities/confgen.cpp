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

#include <algorithm>
#include <cmath>
#include <fmt/core.h>
#include <fstream>
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
    EnergyCalculator calculator(m_method, calc_config);
    calculator.setMolecule(input.front().getMolInfo());

    const Matrix reference_topology = input.front().DistanceMatrix().second;
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
        if ((reference_topology - mol.DistanceMatrix().second).cwiseAbs().sum() > 1e-4) {
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

    CurcumaLogger::success_fmt("ConfGen: wrote {}.torsions.csv (per structure), {}.torsion_states.csv "
                               "(states + populations), {}.matched_pairs.csv ({} transition(s))",
        base, base, base, static_cast<int>(transitions.size()));
}
