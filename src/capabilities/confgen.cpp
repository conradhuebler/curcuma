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
#include <random>
#include <set>
#include <sstream>

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
    m_analysis_file = m_config.get<std::string>("analysis_file");
    m_proposal_ranking = m_config.get<std::string>("proposal_ranking");
    m_concerted_max = m_config.get<int>("concerted_max");
    m_isomerise_max = m_config.get<int>("isomerise_max");
    m_isomerise_scan_steps = m_config.get<int>("isomerise_scan_steps");
    m_isomerise_min_separation = m_config.get<double>("isomerise_min_separation");
    m_isomerise_max_rise = m_config.get<double>("isomerise_max_rise");
    m_proposal_novelty_weight = m_config.get<double>("proposal_novelty_weight");
    m_proposal_depth = m_config.get<int>("proposal_depth");
    m_proposal_candidate_cap = m_config.get<int>("proposal_candidate_cap");
    m_proposal_seed = m_config.get<int>("proposal_seed");
    m_crossover_max = m_config.get<int>("crossover_max");
    m_crossover_window = m_config.get<int>("crossover_window");
    m_mode_max = m_config.get<int>("mode_max");
    m_mode_count = m_config.get<int>("mode_count");
    m_mode_amplitude = m_config.get<double>("mode_amplitude");
    m_path_max = m_config.get<int>("path_max");
    m_path_images = m_config.get<int>("path_images");
    m_path_min_distance = m_config.get<int>("path_min_distance");
    m_nci_staged = m_config.get<bool>("nci_staged");
    m_eval_method = m_config.get<std::string>("eval_method");
    /* Claude Generated (Aug 2026): "none" is the caller-side sentinel for "no split" (that is the
     * default of ConfSearch's -confgen_eval_method). Normalise it here as well, so a caller that
     * forwards its sentinel verbatim cannot make this class build a calculator for a method
     * literally named "none" -- which returns a finite garbage energy with a NaN gradient and
     * silently kills every proposal. Fixed on the ConfSearch side too; this is the second net. */
    if (m_eval_method == "none" || m_eval_method == "None" || m_eval_method == "NONE")
        m_eval_method.clear();
    m_clash_factor = m_config.get<double>("clash_factor");
    m_new_rmsd = m_config.get<double>("new_rmsd");
    m_topology_factor = m_config.get<double>("topology_factor");
    // Claude Generated (Jul 2026): restrained build (P0)
    m_restrained_build = m_config.get<bool>("restrained_build");
    m_restraint_force = m_config.get<double>("restraint_force");
    m_restraint_max_iterations = m_config.get<int>("restraint_max_iterations");
    // Claude Generated (Jul 2026): NCI pattern
    m_nci_analysis = m_config.get<bool>("nci_analysis");
    m_nci_hbond_distance = m_config.get<double>("nci_hbond_distance");
    m_nci_hbond_angle = m_config.get<double>("nci_hbond_angle");
    m_nci_xbond_distance = m_config.get<double>("nci_xbond_distance");
    m_nci_xbond_angle = m_config.get<double>("nci_xbond_angle");
    m_nci_contact_distance = m_config.get<double>("nci_contact_distance");
    m_nci_charge_product = m_config.get<double>("nci_charge_product");
    m_nci_min_population = m_config.get<int>("nci_min_population");
    m_nci_generate = m_config.get<bool>("nci_generate");
    m_nci_max_proposals = m_config.get<int>("nci_max_proposals");
    m_nci_depth = m_config.get<int>("nci_depth");
    m_nci_form_distance = m_config.get<double>("nci_form_distance");
    m_nci_break_distance = m_config.get<double>("nci_break_distance");
    m_nci_restraint_force = m_config.get<double>("nci_restraint_force");
    m_consensus_build = m_config.get<bool>("consensus_build");
    m_consensus_max = m_config.get<int>("consensus_max");
    m_proposal_memory_file = m_config.get<std::string>("proposal_memory_file");
}

std::string ConfGen::NCIContact::label() const
{
    const char* tag = (kind == HBond) ? "HB" : (kind == XBond) ? "XB" : (kind == Ionic) ? "EL" : "CT";
    if (kind == HBond)
        return fmt::format("{} {}-{}...{}", tag, first + 1, hydrogen + 1, second + 1);
    return fmt::format("{} {}...{}", tag, first + 1, second + 1);
}

/**
 * @brief Detect the non-covalent interactions of one structure.
 *
 * Four kinds, each with the criterion that the corresponding GFN-FF term uses in spirit:
 *
 *  - HBond   D-H...A with D,A electronegative, r(H...A) below nci_hbond_distance and the
 *            D-H...A angle above nci_hbond_angle. Directional, because the term is.
 *  - XBond   C-X...B with X a heavy halogen: the sigma hole sits opposite the covalent bond, so
 *            the C-X...B angle must be near linear (nci_xbond_angle).
 *  - Ionic   a close heavy-atom pair whose partial charges attract each other strongly
 *            (q_i q_j < -nci_charge_product). This is where the Coulomb term lives.
 *  - Contact any remaining close heavy-atom pair -- dispersion and repulsion.
 *
 * Pairs closer than four bonds are excluded throughout: those are 1-2, 1-3 and 1-4 relations, whose
 * energy the bonded terms already carry. The bond graph comes from the same explicit covalent-radius
 * criterion as everywhere else in this class.
 */
std::vector<ConfGen::NCIContact> ConfGen::detectNCI(const Molecule& mol, const std::vector<double>& charges) const
{
    const int n = mol.AtomCount();
    const Geometry geom = mol.getGeometry();
    auto element = [&](int i) { return mol.Atom(i).first; };
    auto position = [&](int i) { return Eigen::Vector3d(geom.row(i)); };
    auto distance = [&](int i, int j) { return (position(i) - position(j)).norm(); };
    auto angle = [&](int i, int j, int k) { // angle at j, degrees
        const Eigen::Vector3d a = position(i) - position(j), b = position(k) - position(j);
        const double denominator = a.norm() * b.norm();
        if (denominator < 1e-12)
            return 0.0;
        return std::acos(std::max(-1.0, std::min(1.0, a.dot(b) / denominator))) * 180.0 / M_PI;
    };

    // Bond graph and pairwise bond-count distance (BFS, capped at 4 -- more is not needed).
    std::vector<std::vector<int>> neighbours(n);
    for (const auto& [i, j] : topologyFingerprint(mol, m_topology_factor)) {
        neighbours[i].push_back(j);
        neighbours[j].push_back(i);
    }
    std::vector<std::vector<int>> bond_distance(n, std::vector<int>(n, 99));
    for (int start = 0; start < n; ++start) {
        bond_distance[start][start] = 0;
        std::vector<int> frontier { start };
        for (int depth = 1; depth <= 4 && !frontier.empty(); ++depth) {
            std::vector<int> next;
            for (int atom : frontier)
                for (int nb : neighbours[atom])
                    if (bond_distance[start][nb] > depth) {
                        bond_distance[start][nb] = depth;
                        next.push_back(nb);
                    }
            frontier.swap(next);
        }
    }

    auto is_donor_acceptor = [](int z) { return z == 7 || z == 8 || z == 9 || z == 16; };
    auto is_halogen = [](int z) { return z == 17 || z == 35 || z == 53; };
    auto is_heavy = [](int z) { return z > 1; };

    std::vector<NCIContact> contacts;
    std::set<std::pair<int, int>> claimed; // heavy pairs already described by a specific interaction

    // --- hydrogen bonds ---------------------------------------------------------------------------
    for (int h = 0; h < n; ++h) {
        if (element(h) != 1 || neighbours[h].empty())
            continue;
        const int donor = neighbours[h].front();
        if (!is_donor_acceptor(element(donor)))
            continue;
        for (int a = 0; a < n; ++a) {
            if (a == donor || !is_donor_acceptor(element(a)) || bond_distance[donor][a] <= 3)
                continue;
            if (distance(h, a) > m_nci_hbond_distance || angle(donor, h, a) < m_nci_hbond_angle)
                continue;
            NCIContact c;
            c.kind = NCIContact::HBond;
            c.first = donor;
            c.second = a;
            c.hydrogen = h;
            contacts.push_back(c);
            claimed.insert({ std::min(donor, a), std::max(donor, a) });
        }
    }

    // --- halogen bonds ----------------------------------------------------------------------------
    for (int x = 0; x < n; ++x) {
        if (!is_halogen(element(x)) || neighbours[x].empty())
            continue;
        const int carbon = neighbours[x].front();
        for (int b = 0; b < n; ++b) {
            if (b == x || !(is_donor_acceptor(element(b)) || is_halogen(element(b))))
                continue;
            if (bond_distance[x][b] <= 3)
                continue;
            if (distance(x, b) > m_nci_xbond_distance || angle(carbon, x, b) < m_nci_xbond_angle)
                continue;
            NCIContact c;
            c.kind = NCIContact::XBond;
            c.first = x;
            c.second = b;
            contacts.push_back(c);
            claimed.insert({ std::min(x, b), std::max(x, b) });
        }
    }

    // --- electrostatic and plain contacts ---------------------------------------------------------
    const bool have_charges = static_cast<int>(charges.size()) == n;
    for (int i = 0; i < n; ++i) {
        if (!is_heavy(element(i)))
            continue;
        for (int j = i + 1; j < n; ++j) {
            if (!is_heavy(element(j)) || bond_distance[i][j] <= 3)
                continue;
            if (distance(i, j) > m_nci_contact_distance)
                continue;
            if (claimed.count({ i, j }))
                continue;
            NCIContact c;
            c.first = i;
            c.second = j;
            c.kind = (have_charges && charges[i] * charges[j] < -m_nci_charge_product)
                ? NCIContact::Ionic
                : NCIContact::Contact;
            contacts.push_back(c);
        }
    }
    return contacts;
}

/**
 * @brief Turn the per-structure contact lists into one binary presence vector per structure.
 *
 * The construction mirrors the torsion states exactly: a global, sorted list of everything observed
 * anywhere (m_nci_pairs) plus a 0/1 vector per frame. Contacts seen in fewer than nci_min_population
 * structures are dropped -- like a torsion with only one state they carry no contrast, and with a few
 * hundred structures there are many of them.
 */
void ConfGen::buildNCISpace()
{
    std::map<std::tuple<int, int, int, int>, int> population; // (kind, first, second, hydrogen) -> count
    std::vector<std::vector<NCIContact>> per_frame(m_frames.size());
    for (std::size_t f = 0; f < m_frames.size(); ++f) {
        per_frame[f] = detectNCI(m_frames[f].molecule, m_frames[f].charges);
        for (const NCIContact& c : per_frame[f])
            population[{ static_cast<int>(c.kind), c.first, c.second, -1 }]++;
    }

    m_nci_pairs.clear();
    std::map<std::tuple<int, int, int, int>, int> index;
    for (const auto& [key, count] : population) {
        if (count < m_nci_min_population || count == static_cast<int>(m_frames.size()))
            continue; // never seen enough, or always present -- neither distinguishes conformers
        NCIContact c;
        c.kind = static_cast<NCIContact::Kind>(std::get<0>(key));
        c.first = std::get<1>(key);
        c.second = std::get<2>(key);
        index[key] = static_cast<int>(m_nci_pairs.size());
        m_nci_pairs.push_back(c);
    }

    for (std::size_t f = 0; f < m_frames.size(); ++f) {
        m_frames[f].nci.assign(m_nci_pairs.size(), 0);
        for (const NCIContact& c : per_frame[f]) {
            auto it = index.find({ static_cast<int>(c.kind), c.first, c.second, -1 });
            if (it != index.end()) {
                m_frames[f].nci[it->second] = 1;
                if (c.kind == NCIContact::HBond)
                    m_nci_pairs[it->second].hydrogen = c.hydrogen;
            }
        }
    }
}

bool ConfGen::analyseEnsemble()
{
    // ---- 1. read the ensemble -------------------------------------------------------------------
    std::vector<Molecule> input;
    std::vector<bool> is_template;      // parallel to `input`: may this frame serve as a template?
    {
        FileIterator file(Filename());
        while (!file.AtEnd()) {
            Molecule mol = file.Next();
            if (mol.AtomCount() > 0) {
                mol.setCharge(m_charge);
                mol.setSpin(m_spin);
                input.push_back(mol);
                is_template.push_back(true);
            }
        }
    }
    // Claude Generated (Aug 2026): additional structures for the DESCRIPTION only. The state and
    // contact statistics are estimated from a handful of structures when only one cycle is available
    // -- 6 in one measured case, against 29 torsions and over 100 contacts -- while the cumulative
    // pool of the run holds one to two orders of magnitude more. Templates stay with the file the
    // run was called with; the novelty check profits as well, because a proposal that matches an
    // older structure is then correctly recognised as known.
    if (!m_analysis_file.empty() && m_analysis_file != Filename()) {
        std::ifstream check(m_analysis_file);
        if (check.good()) {
            int added = 0;
            FileIterator file(m_analysis_file);
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                if (mol.AtomCount() > 0) {
                    mol.setCharge(m_charge);
                    mol.setSpin(m_spin);
                    input.push_back(mol);
                    is_template.push_back(false);
                    added++;
                }
            }
            if (added > 0 && m_verbosity >= 1)
                CurcumaLogger::result_fmt("ConfGen: description built from {} structure(s) of this call plus "
                                          "{} from {} -- templates remain the ones of this call",
                    static_cast<int>(is_template.size()) - added, added, m_analysis_file);
        } else {
            CurcumaLogger::warn("ConfGen: -analysis_file not readable, using the input structures only: "
                + m_analysis_file);
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
    // Claude Generated (Aug 2026): freeze the GFN-FF topology, as every other child computation of
    // the search does (ConfSearch::ChildConfig). Without it the "auto" mode re-derives the topology
    // whenever an atom moved far enough, and each re-derivation shifts the ENERGY SCALE -- with one
    // calculator shared across all proposals (reuse_calculator) that shift then survives into every
    // later structure. Measured on a 107-atom peptide: a proposal was written with -18.968697 Eh
    // while the same geometry recomputes to -18.577015 Eh, a difference of 1028 kJ/mol. It passed
    // the topology check, became the cycle's "new best", and its energy window then excluded every
    // real conformer -- the cycle ended with 0 structures in the ensemble.
    {
        const json full = m_config.exportConfig();
        if (full.contains("gfnff") && full["gfnff"].is_object())
            calc_config["gfnff"] = full["gfnff"];
        if (!(calc_config.contains("gfnff") && calc_config["gfnff"].contains("topology_mode")))
            calc_config["gfnff"]["topology_mode"] = "constant";
    }
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
    for (std::size_t frame_index = 0; frame_index < input.size(); ++frame_index) {
        const Molecule& mol = input[frame_index];
        Frame frame;
        frame.molecule = mol;
        frame.from_input = frame_index < is_template.size() ? is_template[frame_index] : true;

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
        // Claude Generated (Jul 2026): keep the partial charges of THIS single point -- the NCI
        // detection needs them to tell an electrostatic contact from a plain close one, and taking
        // them here avoids a second single point per structure later.
        const Vector charges = calculator.Charges();
        frame.charges.assign(charges.data(), charges.data() + charges.size());

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

/**
 * @brief Which energy term drives the spread of the ensemble.
 *
 * Because the total is the sum of its terms, the variance decomposes exactly:
 *     Var(E) = sum_i Cov(E_i, E),
 * so Cov(E_i, E)/Var(E) is the share of the ensemble's energy spread attributable to term i. The
 * shares add to one and may be negative (a term that anticorrelates with the total damps the spread).
 * This is a statement about the ENSEMBLE, not about a single structure pair, and it is the cheapest
 * possible answer to "what makes these conformers differ in energy".
 */
void ConfGen::reportTermVariance() const
{
    const int n = static_cast<int>(m_frames.size());
    if (n < 3 || m_term_names.empty())
        return;
    double mean_total = 0.0;
    for (const Frame& f : m_frames)
        mean_total += f.energy;
    mean_total /= n;
    double variance = 0.0;
    for (const Frame& f : m_frames)
        variance += (f.energy - mean_total) * (f.energy - mean_total);
    variance /= (n - 1);
    if (variance <= 0.0)
        return;

    std::vector<std::pair<std::string, double>> shares; // term -> Cov(E_i, E)/Var(E)
    std::vector<std::pair<std::string, double>> spreads; // term -> own standard deviation, kJ/mol
    for (const std::string& name : m_term_names) {
        double mean_term = 0.0;
        for (const Frame& f : m_frames)
            mean_term += f.terms.at(name);
        mean_term /= n;
        double covariance = 0.0, own = 0.0;
        for (const Frame& f : m_frames) {
            const double dt = f.terms.at(name) - mean_term;
            covariance += dt * (f.energy - mean_total);
            own += dt * dt;
        }
        covariance /= (n - 1);
        own /= (n - 1);
        shares.emplace_back(name, covariance / variance);
        spreads.emplace_back(name, std::sqrt(own) * kEh2kJ);
    }
    std::sort(shares.begin(), shares.end(),
        [](const auto& a, const auto& b) { return std::abs(a.second) > std::abs(b.second); });

    CurcumaLogger::result_fmt("ConfGen: energy spread of the ensemble: {:.2f} kJ/mol (standard deviation)",
        std::sqrt(variance) * kEh2kJ);
    CurcumaLogger::result("ConfGen: which term carries that spread (Cov(term,total)/Var(total)):");
    for (const auto& [name, share] : shares) {
        if (std::abs(share) < 0.005)
            continue;
        const double own = std::find_if(spreads.begin(), spreads.end(),
            [&](const auto& s) { return s.first == name; })
                               ->second;
        CurcumaLogger::result_fmt("            {:<14} {:+6.1f} %   (own spread {:6.2f} kJ/mol)",
            name, 100.0 * share, own);
    }
}

/**
 * @brief Report the NCI pattern as a conformer description and compare it with the torsion states.
 *
 * The two descriptions are printed side by side on purpose. A description is only useful if it
 * SEPARATES the ensemble: if many structures share a description, it cannot be the basis for telling
 * conformers apart, let alone for proposing new ones.
 */
void ConfGen::reportNCISpace() const
{
    if (m_nci_pairs.empty()) {
        CurcumaLogger::warn("ConfGen: no variable non-covalent contact found -- the NCI pattern is "
                            "identical in every structure and carries no information");
        return;
    }
    int hb = 0, xb = 0, el = 0, ct = 0;
    for (const NCIContact& c : m_nci_pairs)
        switch (c.kind) {
        case NCIContact::HBond: hb++; break;
        case NCIContact::XBond: xb++; break;
        case NCIContact::Ionic: el++; break;
        default: ct++; break;
        }
    CurcumaLogger::result_fmt("ConfGen: NCI pattern: {} variable contact(s) -- {} H-bond, {} halogen "
                              "bond, {} electrostatic, {} close contact",
        static_cast<int>(m_nci_pairs.size()), hb, xb, el, ct);

    auto hamming = [](const std::vector<int>& a, const std::vector<int>& b) {
        int d = 0;
        for (std::size_t i = 0; i < a.size() && i < b.size(); ++i)
            d += (a[i] != b[i]) ? 1 : 0;
        return d;
    };
    std::set<std::vector<int>> nci_patterns, torsion_patterns;
    for (const Frame& f : m_frames) {
        nci_patterns.insert(f.nci);
        torsion_patterns.insert(f.states);
    }
    const int n = static_cast<int>(m_frames.size());
    int nci_pairs_one = 0, torsion_pairs_one = 0, identical_torsions = 0, identical_nci = 0;
    long long sum_nci = 0, sum_torsion = 0, count = 0;
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b) {
            const int dn = hamming(m_frames[a].nci, m_frames[b].nci);
            const int dt = hamming(m_frames[a].states, m_frames[b].states);
            nci_pairs_one += (dn == 1);
            torsion_pairs_one += (dt == 1);
            identical_nci += (dn == 0);
            identical_torsions += (dt == 0);
            sum_nci += dn;
            sum_torsion += dt;
            count++;
        }
    CurcumaLogger::result_fmt("ConfGen: distinct descriptions of {} structures: {} torsion vector(s), "
                              "{} NCI pattern(s)",
        n, static_cast<int>(torsion_patterns.size()), static_cast<int>(nci_patterns.size()));
    CurcumaLogger::result_fmt("ConfGen: mean Hamming distance: {:.1f} (torsions) vs {:.1f} (NCI); "
                              "pairs that are identical: {} vs {}",
        static_cast<double>(sum_torsion) / std::max<long long>(1, count),
        static_cast<double>(sum_nci) / std::max<long long>(1, count),
        identical_torsions, identical_nci);
    CurcumaLogger::result_fmt("ConfGen: pairs at Hamming distance 1 (the ones a matched-pair analysis "
                              "needs): {} in torsion space, {} in NCI space",
        torsion_pairs_one, nci_pairs_one);

    // The full contact map is high-dimensional and separates almost anything, which makes "it
    // separates the ensemble" a weak statement. The hydrogen-bond subset is the chemically specific
    // and much smaller one -- if it already separates structures that the torsion vector cannot tell
    // apart, the conclusion does not rest on the dimensionality of the descriptor.
    std::vector<int> hb_index;
    for (std::size_t k = 0; k < m_nci_pairs.size(); ++k)
        if (m_nci_pairs[k].kind == NCIContact::HBond)
            hb_index.push_back(static_cast<int>(k));
    if (hb_index.empty())
        return;
    auto hb_pattern = [&](const Frame& f) {
        std::vector<int> pattern;
        pattern.reserve(hb_index.size());
        for (int k : hb_index)
            pattern.push_back(k < static_cast<int>(f.nci.size()) ? f.nci[k] : 0);
        return pattern;
    };
    std::set<std::vector<int>> hb_patterns;
    for (const Frame& f : m_frames)
        hb_patterns.insert(hb_pattern(f));
    int hb_identical = 0, hb_identical_and_torsion_identical = 0, torsion_identical = 0;
    for (int a = 0; a < n; ++a)
        for (int b = a + 1; b < n; ++b) {
            const bool same_torsion = (m_frames[a].states == m_frames[b].states);
            const bool same_hb = (hb_pattern(m_frames[a]) == hb_pattern(m_frames[b]));
            hb_identical += same_hb;
            torsion_identical += same_torsion;
            hb_identical_and_torsion_identical += (same_hb && same_torsion);
        }
    CurcumaLogger::result_fmt("ConfGen: hydrogen-bond subset alone ({} bonds): {} distinct pattern(s); "
                              "{} identical pair(s) versus {} in torsion space",
        static_cast<int>(hb_index.size()), static_cast<int>(hb_patterns.size()), hb_identical, torsion_identical);
    if (torsion_identical > 0)
        CurcumaLogger::result_fmt("ConfGen: of the {} pair(s) with an identical torsion vector, {} also "
                                  "share their hydrogen-bond pattern -- the remaining {} are structures "
                                  "the torsion description cannot tell apart but the H-bond pattern can",
            torsion_identical, hb_identical_and_torsion_identical,
            torsion_identical - hb_identical_and_torsion_identical);
}

void ConfGen::writeNCITable(const std::string& path) const
{
    if (m_nci_pairs.empty())
        return;
    std::ofstream out(path);
    out << "# structure,energy_Eh";
    for (const NCIContact& c : m_nci_pairs)
        out << "," << c.label();
    out << "\n";
    for (std::size_t f = 0; f < m_frames.size(); ++f) {
        out << (f + 1) << "," << std::fixed << std::setprecision(8) << m_frames[f].energy;
        for (std::size_t k = 0; k < m_nci_pairs.size(); ++k)
            out << "," << (k < m_frames[f].nci.size() ? m_frames[f].nci[k] : 0);
        out << "\n";
    }
    CurcumaLogger::info_fmt("ConfGen: NCI pattern written to {}", path);
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
    // Reference for the term columns: the deepest structure of the ensemble. Any fixed member does,
    // the choice only shifts the intercept -- but a member (rather than the mean) keeps every column
    // a difference between two REAL structures, which is what the matched-pair analysis reports too.
    if (m_reference_terms.empty() && !m_frames.empty()) {
        const Frame* deepest = &m_frames.front();
        for (const Frame& f : m_frames)
            if (f.energy < deepest->energy)
                deepest = &f;
        m_reference_terms = deepest->terms;
    }
    switch (level) {
    case 0: fit.name = "constant"; break;
    case 1: fit.name = "torsions (marginals)"; break;
    case 2: fit.name = "torsions + pair couplings"; break;
    case 3: fit.name = "NCI pattern"; break;
    // Claude Generated (Aug 2026): these two are CONSISTENCY CHECKS, not predictions --
    // E = sum of its own terms is an identity, so a perfect score is the expected result and
    // says only that the decomposition is complete and the differences are formed correctly.
    case 5: fit.name = "[Pruefung] Energieterme (Identitaet)"; break;
    case 6: fit.name = "[Pruefung] Torsionen + NCI + Terme"; break;
    default: fit.name = "torsions + NCI pattern"; break;
    }

    const std::vector<int> informative = informativeTorsions();
    // Column layout: [intercept][per torsion: states 1..k-1][per torsion pair: (a,b) with a,b >= 1]
    // [per NCI contact: present/absent]. State 0 is the reference level, so the torsion columns stay
    // linearly independent by construction; an NCI column is an indicator and needs no reference.
    struct Column {
        int torsion_a = -1, state_a = -1, torsion_b = -1, state_b = -1;
        int nci = -1; ///< index into m_nci_pairs; >= 0 marks an NCI column
        /** Claude Generated (Aug 2026): name of an energy term of the decomposition. Such a column
         *  holds the term's VALUE, not an indicator -- the description then contains the physical
         *  quantities themselves instead of only the discrete pattern that produces them. Ten
         *  columns for the whole force field, so it costs nothing in degrees of freedom, and the
         *  terms are what the variance attribution identified as the carriers of the spread. */
        std::string term;
    };
    std::vector<Column> columns;
    const bool use_torsions = (level == 1 || level == 2 || level == 4 || level == 6);
    const bool use_nci = (level == 3 || level == 4 || level == 6);
    const bool use_terms = (level == 5 || level == 6);
    if (use_torsions)
        for (int t : informative)
            for (std::size_t s = 1; s < m_state_centres[t].size(); ++s)
                columns.push_back({ t, static_cast<int>(s), -1, -1, -1 });
    if (level == 2)
        for (std::size_t ia = 0; ia < informative.size(); ++ia)
            for (std::size_t ib = ia + 1; ib < informative.size(); ++ib)
                for (std::size_t sa = 1; sa < m_state_centres[informative[ia]].size(); ++sa)
                    for (std::size_t sb = 1; sb < m_state_centres[informative[ib]].size(); ++sb)
                        columns.push_back({ informative[ia], static_cast<int>(sa),
                            informative[ib], static_cast<int>(sb), -1 });
    if (use_nci) {
        // Column budget. There are typically far more variable contacts than structures, and an
        // indicator model with p > n interpolates every training set exactly -- its cross-validation
        // is then either impossible or meaningless. The contacts are therefore ranked by CONTRAST,
        // min(population, n - population): an interaction present in half the ensemble separates it,
        // one present in two structures does not. Only the most contrasting ones enter the model, at
        // most a third of the structure count, which is the same budget the torsion model gets.
        std::vector<std::pair<int, int>> ranked; // (contrast, index)
        for (std::size_t k = 0; k < m_nci_pairs.size(); ++k) {
            int population = 0;
            for (const Frame& f : m_frames)
                population += (k < f.nci.size() && f.nci[k]) ? 1 : 0;
            const int contrast = std::min(population, static_cast<int>(m_frames.size()) - population);
            if (contrast > 0)
                ranked.emplace_back(contrast, static_cast<int>(k));
        }
        std::sort(ranked.begin(), ranked.end(),
            [](const auto& a, const auto& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
        const std::size_t budget = std::max<std::size_t>(1, m_frames.size() / 3);
        if (ranked.size() > budget)
            ranked.resize(budget);
        for (const auto& [contrast, k] : ranked)
            columns.push_back({ -1, -1, -1, -1, k, {} });
    }
    if (use_terms)
        for (const auto& [name, value] : m_frames.front().terms)
            columns.push_back({ -1, -1, -1, -1, -1, name });

    const int n = static_cast<int>(m_frames.size());
    const int p = static_cast<int>(columns.size()) + 1; // + intercept
    fit.columns = p;

    auto buildRow = [&](const Frame& f, Eigen::VectorXd& row) {
        row.setZero(p);
        row(0) = 1.0;
        for (int c = 0; c < static_cast<int>(columns.size()); ++c) {
            const Column& col = columns[c];
            if (!col.term.empty()) {
                // Claude Generated (Aug 2026): the CHANGE of the term, not its value. The raw terms
                // are dominated by their offset -- the bond term sits at -50360 kJ/mol and varies by
                // 13 -- so as absolute columns they are numerically indistinguishable from the
                // intercept and the decomposition drops them (measured: rank 1 for twelve columns).
                // Referenced to the deepest structure of the ensemble, the column holds exactly what
                // distinguishes the conformers, in the same units as the target.
                const auto it = f.terms.find(col.term);
                const auto ir = m_reference_terms.find(col.term);
                row(c + 1) = ((it != f.terms.end()) ? it->second : 0.0)
                    - ((ir != m_reference_terms.end()) ? ir->second : 0.0);
                continue;   // ohne dies liefe der Torsionszweig darunter weiter -- und der liest
                            // f.states[-1] und ueberschreibt den Wert (Ursache des Rangs 1)
            } else if (col.nci >= 0) {
                row(c + 1) = (col.nci < static_cast<int>(f.nci.size()) && f.nci[col.nci]) ? 1.0 : 0.0;
                continue;
            }
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
    if (residuals.empty()) {
        // Every fold was skipped because the training part had fewer rows than the model has
        // parameters. Reporting rmse_cv = 0 here would look like a perfect model; it means the
        // opposite -- the model could not be tested at all.
        fit.evaluated = false;
        return fit;
    }
    fit.evaluated = true;
    {
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
        if (!f.evaluated) {
            CurcumaLogger::result_fmt("  {:<28} {:>7} {:>7} {:>34} {:>9.2f} kJ",
                f.name, f.columns, f.rank, "not testable (more parameters than data)", f.rmse_in);
            continue;
        }
        // R2_cv against the constant model: the fraction of the energy variation that the model
        // predicts on data it has not seen. This is the number that decides whether the torsion-state
        // picture describes this ensemble at all.
        const double r2 = (null_rmse > 0.0) ? 1.0 - (f.rmse_cv * f.rmse_cv) / (null_rmse * null_rmse) : 0.0;
        CurcumaLogger::result_fmt("  {:<28} {:>7} {:>7} {:>9.2f} kJ {:>9.2f} kJ {:>9.0f} % {:>9.2f} kJ",
            f.name, f.columns, f.rank, f.rmse_cv, f.mae_cv, r2 * 100.0, f.rmse_in);
    }
    CurcumaLogger::info("RMSE_cv/medAE_cv are out-of-sample (k-fold). RMSE_in always improves with more "
                        "parameters and is shown for reference only.");
    CurcumaLogger::info("Rows marked [Pruefung] are not models: the energy is the SUM of its own terms, so "
                        "reproducing it from them is an identity. A perfect score there confirms that the "
                        "decomposition is complete -- it predicts nothing.");

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
    const bool add_useful = fits[1].evaluated
        && (add_err < null_err * (1.0 - kMinRelGain)) && (fits[1].mae_cv <= fits[0].mae_cv);
    const bool pair_useful = fits[2].evaluated
        && (pair_err < add_err * (1.0 - kMinRelGain)) && (fits[2].mae_cv <= fits[1].mae_cv);

    // Claude Generated (Aug 2026): a model that could not be cross-validated (more parameters than
    // training rows in every fold) has rmse_cv = 0, which would print as "explains 100 %" -- the
    // exact opposite of what it means. Such a model gets no verdict at all.
    // Report each level only if it was actually tested. A model whose folds were all skipped has
    // rmse_cv = 0 and would otherwise print as "explains 100 %" -- the opposite of what it means.
    auto not_testable = [&](const ModelFit& f) {
        CurcumaLogger::warn_fmt("ConfGen: the '{}' model could not be tested on this ensemble -- {} "
                                "parameters for {} structures, so every cross-validation fold had fewer "
                                "training rows than parameters. More structures, or a coarser "
                                "-state_tolerance, are needed before it can be judged at all.",
            f.name, f.columns, static_cast<int>(m_frames.size()));
    };
    if (!fits[1].evaluated) {
        not_testable(fits[1]);
        if (!fits[2].evaluated)
            not_testable(fits[2]);
        return;
    }
    if (fits[2].evaluated)
        CurcumaLogger::result_fmt("ConfGen: torsion states explain {:.0f}% of the energy variation out of "
                                  "sample, adding pair couplings takes it to {:.0f}%",
            r2_add * 100.0, r2_pair * 100.0);
    else {
        CurcumaLogger::result_fmt("ConfGen: torsion states explain {:.0f}% of the energy variation out of "
                                  "sample", r2_add * 100.0);
        not_testable(fits[2]);
    }

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
bool ConfGen::restrainedBuild(const Proposal& p, Molecule& driven, const Molecule* start) const
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
    /* Claude Generated (Aug 2026): explicit target angles for the isomerisation move. They are not
     * state centres -- the whole point is that the ensemble never showed this value, so it cannot
     * be addressed by a state index. */
    for (const auto& at : p.angle_targets) {
        if (at.first < 0 || at.first >= static_cast<int>(m_torsions.size()))
            continue;
        Optimization::DihedralRestraint r;
        r.i = m_torsions[at.first].i;
        r.j = m_torsions[at.first].j;
        r.k = m_torsions[at.first].k;
        r.l = m_torsions[at.first].l;
        r.target = at.second * M_PI / 180.0;
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
    // Claude Generated (Aug 2026): ONE shared calculator for the whole run (see m_calculator) --
    // the optimiser must not re-derive the force field per structure. Saves a full GFN-FF setup
    // per proposal and keeps the run out of the intermittent parameter-generation crash.
    opt_config["reuse_calculator"] = true;
    opt_config["write_trajectory"] = false;
    opt_config["max_iterations"] = m_restraint_max_iterations;
    opt_config["dihedral_restraints"] = Optimization::GeometryRestraints::toJson(restraints);

    // Default: start from the TEMPLATE -- clash-free and with the correct topology, so the restrained
    // optimisation performs the rotation itself. With an explicit start geometry (the rigidly built
    // one) the same restraints act as a clash repair instead; see the parameter documentation.
    Molecule mol = start ? *start : m_frames[p.template_frame].molecule;
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


void ConfGen::loadProposalMemory()
{
    m_proposed_before.clear();
    if (m_proposal_memory_file.empty())
        return;
    std::ifstream in(m_proposal_memory_file);
    if (!in.good())
        return;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')
            continue;
        // Claude Generated (Aug 2026): one line holds BOTH descriptions, separated by '|' --
        // the torsion state vector and the contact pattern of the structure that was built. Both
        // move sets read both, so a torsion proposal knows which patterns already exist and an NCI
        // move knows which state vectors do. Files written before this change carry only the state
        // vector and are still read correctly.
        const std::size_t bar = line.find('|');
        auto parse = [](const std::string& s) {
            std::vector<int> v;
            std::istringstream is(s);
            int x;
            while (is >> x) v.push_back(x);
            return v;
        };
        std::vector<int> states = parse(bar == std::string::npos ? line : line.substr(0, bar));
        if (!states.empty())
            m_proposed_before.insert(states);
        if (bar != std::string::npos) {
            std::vector<int> pattern = parse(line.substr(bar + 1));
            if (!pattern.empty())
                m_patterns_before.insert(pattern);
        }
    }
    if (!m_proposed_before.empty())
        CurcumaLogger::result_fmt("ConfGen: {} state vector(s) were already proposed earlier and are skipped",
            static_cast<int>(m_proposed_before.size()));
}

void ConfGen::appendProposalMemory(
    const std::vector<std::pair<std::vector<int>, std::vector<int>>>& entries) const
{
    if (m_proposal_memory_file.empty() || entries.empty())
        return;
    std::ofstream out(m_proposal_memory_file, std::ios::app);
    for (const auto& [states, pattern] : entries) {
        for (std::size_t i = 0; i < states.size(); ++i)
            out << states[i] << (i + 1 < states.size() ? " " : "");
        if (!pattern.empty()) {
            out << " |";
            for (int v : pattern)
                out << " " << v;
        }
        out << "\n";
    }
}

std::vector<std::vector<double>> ConfGen::stateEnergies() const
{
    const double kT = kBoltzmannEh * m_temperature;
    std::vector<std::vector<double>> energies(m_torsions.size());
    for (std::size_t t = 0; t < m_torsions.size(); ++t) {
        energies[t].assign(m_state_centres[t].size(), std::numeric_limits<double>::quiet_NaN());
        for (std::size_t s = 0; s < m_state_centres[t].size(); ++s) {
            double weight_sum = 0.0, weighted = 0.0;
            for (const Frame& f : m_frames) {
                if (f.states[t] != static_cast<int>(s))
                    continue;
                const double rel = f.energy - m_reference_energy;
                const double w = std::exp(-rel / kT);
                weight_sum += w;
                weighted += w * rel;
            }
            if (weight_sum > 0.0)
                energies[t][s] = weighted / weight_sum * kEh2kJ;
        }
    }
    return energies;
}

std::vector<ConfGen::Proposal> ConfGen::generateConsensusProposals() const
{
    std::vector<Proposal> proposals;
    const std::vector<std::vector<double>> energies = stateEnergies();

    // The consensus vector: per torsion the state with the lowest Boltzmann-weighted mean energy.
    // Torsions with a single populated state contribute that state -- they carry no choice.
    std::vector<int> consensus(m_torsions.size(), 0);
    std::vector<std::pair<double, int>> margins; // (energy gap best->second best, torsion)
    for (std::size_t t = 0; t < m_torsions.size(); ++t) {
        int best = -1, second = -1;
        for (std::size_t s = 0; s < energies[t].size(); ++s) {
            if (std::isnan(energies[t][s]))
                continue;
            if (best < 0 || energies[t][s] < energies[t][best]) {
                second = best;
                best = static_cast<int>(s);
            } else if (second < 0 || energies[t][s] < energies[t][second]) {
                second = static_cast<int>(s);
            }
        }
        if (best < 0)
            return proposals; // no populated state at all -- nothing to assemble
        consensus[t] = best;
        if (second >= 0)
            margins.emplace_back(energies[t][second] - energies[t][best], static_cast<int>(t));
    }

    std::set<std::vector<int>> known;
    for (const Frame& f : m_frames)
        known.insert(f.states);

    // Template = the ensemble member CLOSEST to the target vector, so the build has to drive as few
    // torsions as possible. Deliberately not the lowest-energy structure: here the geometry is a
    // starting point, not a candidate.
    auto make = [&](const std::vector<int>& states, const std::string& tag) {
        if (known.count(states))
            return;
        int best_frame = 0, best_distance = std::numeric_limits<int>::max();
        for (std::size_t i = 0; i < m_frames.size(); ++i) {
            const int d = TorsionSpace::hammingDistance(m_frames[i].states, states);
            if (d < best_distance) {
                best_distance = d;
                best_frame = static_cast<int>(i);
            }
        }
        Proposal p;
        p.states = states;
        p.template_frame = best_frame;
        p.distance = best_distance;
        p.nci_label = tag; // reused as a label; nci_targets stays empty, so this is a torsion build
        proposals.push_back(std::move(p));
    };

    make(consensus, "consensus");

    // Claude Generated (Aug 2026): walk AWAY from the ensemble, not around the consensus.
    //
    // Flipping the least certain torsions of the consensus one at a time produced assemblies that
    // still sat one torsion away from an existing structure -- the consensus itself is usually a
    // sampled structure, so its neighbours are too, and the whole point of a de-novo build is lost.
    // What is needed are vectors far from EVERYTHING, and they are reachable: the missing reference
    // structure of the test system differs from every found conformer in at least 15 of 29 torsions,
    // yet all 29 of its states occur individually in the ensemble. So the states are there; only
    // their combination is not.
    //
    // Each step therefore takes the flip that gains the most distance to the nearest ensemble member,
    // breaking ties by the smallest energy penalty, and emits the result. The assemblies form a
    // trajectory of increasing novelty instead of a cloud around the average structure.
    auto min_distance = [&](const std::vector<int>& v) {
        int best = std::numeric_limits<int>::max();
        for (const Frame& f : m_frames)
            best = std::min(best, TorsionSpace::hammingDistance(f.states, v));
        return best;
    };

    std::vector<int> current = consensus;
    while (static_cast<int>(proposals.size()) < m_consensus_max) {
        int best_torsion = -1, best_state = -1, best_gain = -1;
        double best_penalty = std::numeric_limits<double>::infinity();
        const int base_distance = min_distance(current);
        for (std::size_t t = 0; t < m_torsions.size(); ++t) {
            for (std::size_t s = 0; s < energies[t].size(); ++s) {
                if (std::isnan(energies[t][s]) || static_cast<int>(s) == current[t])
                    continue;
                std::vector<int> trial = current;
                trial[t] = static_cast<int>(s);
                const int gain = min_distance(trial) - base_distance;
                const double penalty = energies[t][s] - energies[t][consensus[t]];
                if (gain > best_gain || (gain == best_gain && penalty < best_penalty)) {
                    best_gain = gain;
                    best_penalty = penalty;
                    best_torsion = static_cast<int>(t);
                    best_state = static_cast<int>(s);
                }
            }
        }
        if (best_torsion < 0)
            break;
        current[best_torsion] = best_state;
        make(current, fmt::format("de-novo, {} torsion(s) from the consensus, {} from the nearest structure",
                 TorsionSpace::hammingDistance(current, consensus), min_distance(current)));
    }
    return proposals;
}

/* Claude Generated (Aug 2026): CROSSOVER -- transfer a CONNECTED WINDOW of torsions from one
 * conformer into another, instead of changing d torsions independently.
 *
 * Why this move exists. The reference conformer of the measured 107-atom peptide differs from the
 * nearest of ~4600 sampled structures in 6-7 of 29 torsions by more than 60 degrees -- further than
 * any mutation-style move reaches in one step. Raising -proposal_depth to 5 or 7 does reach that far
 * and is useless in practice: measured, NOT ONE of the built structures kept its target state vector
 * through the free optimisation, because a random combination of seven torsions is almost never near
 * a minimum. "Six to seven torsions at once" is not seven independent events, it is the signature of
 * one collective rearrangement.
 *
 * A window transfer respects that. Every transferred block comes from a geometry that was already
 * viable, so the local fold it carries is a real one and only its context changes. This is fragment
 * insertion as used in protein folding, with the ensemble itself as the fragment library -- crossover
 * rather than mutation.
 *
 * Connectivity is defined ON THE TORSIONS: two are neighbours when their central bonds share an
 * atom, and a window is grown breadth-first from a seed. A window in storage order would be an
 * arbitrary set of unrelated dihedrals.
 */
/* Claude Generated (Aug 2026): see the doc comment in the header. */
std::vector<ConfGen::Proposal> ConfGen::generateIsomerisationProposals() const
{
    std::vector<Proposal> out;
    if (m_isomerise_max <= 0 || m_frames.empty() || !m_calculator)
        return out;

    /* Which torsions did the ensemble only ever show in ONE state? That is the whole detection,
     * and it is read off the rotamer analysis -- no bond orders, no element list, nothing about
     * amides or guanidinium groups. A single-state torsion is a degree of freedom the sampling
     * never opened; whether that is a rotation barrier or plain bad luck, the generator cannot
     * recombine what it has not seen. */
    std::vector<int> frozen;
    for (std::size_t t = 0; t < m_torsions.size(); ++t)
        if (m_state_centres[t].size() == 1)
            frozen.push_back(static_cast<int>(t));
    if (frozen.empty()) {
        CurcumaLogger::result("ConfGen: no isomerisation candidate -- every torsion was sampled in "
                              "more than one state");
        return out;
    }

    /* WHERE to flip to is measured, not assumed. Earlier this used "the opposite planar value",
     * which is chemical knowledge smuggled in after seeing one system. Instead the frozen torsion
     * is scanned rigidly (the moving side is rotated, everything else held) and the energy profile
     * is read: every local minimum that is far enough from the state the ensemble knows is a
     * candidate target. That finds the cis form of a conjugated bond without being told such a
     * thing exists, and equally a non-planar second minimum, which the planarity rule would miss. */
    const Molecule& probe = m_frames.front().molecule;
    m_calculator->setMolecule(probe.getMolInfo());
    const int steps = std::max(6, m_isomerise_scan_steps);
    const double dphi = 360.0 / steps;
    /* Claude Generated (Aug 2026): the scan is the only part of this move whose cost grows with the
     * system -- frozen torsions x steps single points, on the DESCRIPTION surface (the force field),
     * never on the ranking one. Reported, because "how does this scale" must be answerable from a
     * log rather than from the source. */
    const auto scan_start = std::chrono::steady_clock::now();

    struct Target { int torsion; double angle; double barrier_kJ; double depth_kJ; };
    std::vector<Target> targets;
    for (int t : frozen) {
        const double known = m_state_centres[t][0];
        std::vector<double> ang(steps), en(steps);
        Molecule scan = probe;
        for (int s = 0; s < steps; ++s) {
            ang[s] = known + s * dphi;
            while (ang[s] > 180.0) ang[s] -= 360.0;
            Geometry g = TorsionSpace::setDihedral(probe.getGeometry(), m_torsions[t], ang[s]);
            scan.setGeometry(g);
            m_calculator->setMolecule(scan.getMolInfo());
            en[s] = m_calculator->CalculateEnergy(false);
        }
        const double e0 = en[0];
        // local minima of the cyclic profile, excluding the neighbourhood of the known state
        for (int s = 0; s < steps; ++s) {
            const double prev = en[(s - 1 + steps) % steps], next = en[(s + 1) % steps];
            if (!(en[s] < prev && en[s] < next))
                continue;
            double sep = std::fabs(ang[s] - known);
            if (sep > 180.0) sep = 360.0 - sep;
            if (sep < m_isomerise_min_separation)
                continue;
            /* The rigid profile's maximum is NOT a barrier -- nothing relaxes along the way, so
             * it is dominated by atoms driven into each other. Kept only as a coarse flag, and
             * named for what it is. */
            double crest = 0.0;
            for (int q = 1; q < steps; ++q)
                crest = std::max(crest, en[q] - e0);
            const double rise = (en[s] - e0) * 2625.4996;
            if (rise > m_isomerise_max_rise)
                continue;   // a collision, not a conformer -- see isomerise_max_rise
            targets.push_back({ t, ang[s], crest * 2625.4996, rise });
        }
    }
    const double scan_seconds
        = std::chrono::duration<double>(std::chrono::steady_clock::now() - scan_start).count();
    CurcumaLogger::result_fmt("ConfGen: torsion scan: {} frozen torsion(s) x {} points = {} single "
                              "point(s) at {} in {:.1f} s ({:.0f} ms each)",
        static_cast<int>(frozen.size()), steps, static_cast<int>(frozen.size()) * steps, m_method,
        scan_seconds, frozen.empty() ? 0.0 : 1000.0 * scan_seconds / (frozen.size() * steps));
    if (targets.empty()) {
        CurcumaLogger::result_fmt("ConfGen: {} torsion(s) sampled in a single state, but the scan found "
                                  "no usable second minimum ({:.0f} deg away, at most {:.0f} kJ/mol up) "
                                  "-- nothing to isomerise",
            static_cast<int>(frozen.size()), m_isomerise_min_separation, m_isomerise_max_rise);
        return out;
    }
    // The interesting ones are those whose second minimum is not much worse: a target 200 kJ/mol up
    // is a real minimum and still useless.
    std::sort(targets.begin(), targets.end(),
        [](const Target& a, const Target& b) { return a.depth_kJ < b.depth_kJ; });

    CurcumaLogger::result_fmt("ConfGen: {} torsion(s) held in ONE state by the whole ensemble; the scan "
                              "finds {} second minimum/minima -- states no recombination can invent",
        static_cast<int>(frozen.size()), static_cast<int>(targets.size()));
    for (const Target& g : targets)
        CurcumaLogger::result_fmt("ConfGen:   {} known at {:+.0f} deg -> {:+.0f} deg "
                                  "({:+.1f} kJ/mol above it; rigid crest {:.0f} kJ/mol, not a barrier)",
            m_torsions[g.torsion].label(probe), m_state_centres[g.torsion][0], g.angle,
            g.depth_kJ, g.barrier_kJ);

    /* Templates: the lowest-energy structures. A barrier crossing is worth trying from a good
     * basin -- the flip changes one degree of freedom, everything else stays the template's. */
    std::vector<int> order(m_frames.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
        [this](int a, int b) { return m_frames[a].energy < m_frames[b].energy; });

    for (std::size_t g = 0; g < targets.size() && static_cast<int>(out.size()) < m_isomerise_max; ++g) {
        Proposal p;
        p.states = m_frames[order.front()].states;  // unchanged -- the target angle has no state
        p.template_frame = order.front();
        p.distance = 0;
        p.novelty = 0;
        p.angle_targets.emplace_back(targets[g].torsion, targets[g].angle);
        out.push_back(std::move(p));
    }
    return out;
}

std::vector<ConfGen::Proposal> ConfGen::generateCrossoverProposals() const
{
    std::vector<Proposal> proposals;
    const std::vector<int> informative = informativeTorsions();
    if (m_crossover_max <= 0 || informative.size() < 2 || m_frames.size() < 2)
        return proposals;

    std::set<std::vector<int>> known;
    for (const Frame& f : m_frames)
        known.insert(f.states);
    known.insert(m_proposed_before.begin(), m_proposed_before.end());

    std::vector<int> order;
    for (std::size_t i = 0; i < m_frames.size(); ++i)
        if (m_frames[i].from_input)
            order.push_back(static_cast<int>(i));
    if (order.empty())
        for (std::size_t i = 0; i < m_frames.size(); ++i)
            order.push_back(static_cast<int>(i));
    std::sort(order.begin(), order.end(), [this](int a, int b) { return m_frames[a].energy < m_frames[b].energy; });
    const int n_templates = std::min<int>(std::max(1, m_proposal_templates), static_cast<int>(order.size()));

    // Torsion adjacency over the central bonds.
    std::vector<std::vector<int>> neighbours(informative.size());
    for (std::size_t a = 0; a < informative.size(); ++a)
        for (std::size_t b = a + 1; b < informative.size(); ++b) {
            const auto& ta = m_torsions[informative[a]];
            const auto& tb = m_torsions[informative[b]];
            if (ta.j == tb.j || ta.j == tb.k || ta.k == tb.j || ta.k == tb.k) {
                neighbours[a].push_back(static_cast<int>(b));
                neighbours[b].push_back(static_cast<int>(a));
            }
        }

    // Donors: every structure of this call, thinned evenly so a large cumulative ensemble does not
    // blow the candidate count up (the thinning keeps the energy range, unlike taking the N lowest).
    std::vector<int> donors = order;
    const int donor_cap = 50;
    if (static_cast<int>(donors.size()) > donor_cap) {
        std::vector<int> thinned;
        const double step = static_cast<double>(donors.size()) / donor_cap;
        for (int n = 0; n < donor_cap; ++n)
            thinned.push_back(donors[static_cast<std::size_t>(n * step)]);
        donors = std::move(thinned);
    }

    const std::vector<std::vector<double>> coefficients = additiveCoefficients();
    const int window = std::max(1, m_crossover_window);
    std::set<std::vector<int>> seen_here;
    for (int t = 0; t < n_templates; ++t) {
        const std::vector<int>& base = m_frames[order[t]].states;
        for (int donor : donors) {
            if (donor == order[t] || m_frames[donor].states.size() != base.size())
                continue;
            const std::vector<int>& from = m_frames[donor].states;
            for (std::size_t seed = 0; seed < informative.size(); ++seed) {
                std::vector<char> in_window(informative.size(), 0);
                std::vector<std::size_t> queue = { seed };
                in_window[seed] = 1;
                for (std::size_t q = 0; q < queue.size() && static_cast<int>(queue.size()) < window; ++q)
                    for (int nb : neighbours[queue[q]]) {
                        if (in_window[nb])
                            continue;
                        in_window[nb] = 1;
                        queue.push_back(static_cast<std::size_t>(nb));
                        if (static_cast<int>(queue.size()) >= window)
                            break;
                    }
                std::vector<int> states = base;
                int changed = 0;
                for (std::size_t idx : queue) {
                    const int tor = informative[idx];
                    if (tor >= static_cast<int>(states.size()) || from[tor] == states[tor])
                        continue;
                    states[tor] = from[tor];
                    changed++;
                }
                if (changed == 0 || known.count(states) || seen_here.count(states))
                    continue;
                seen_here.insert(states);
                Proposal p;
                p.states = states;
                p.template_frame = order[t];
                p.distance = changed;
                p.nci_label = fmt::format("crossover, {} torsion(s) from structure {}", changed, donor + 1);
                for (std::size_t k = 0; k < states.size(); ++k)
                    if (states[k] >= 0 && states[k] < static_cast<int>(coefficients[k].size()))
                        p.predicted += coefficients[k][states[k]];
                proposals.push_back(std::move(p));
            }
        }
    }
    if (proposals.empty())
        return proposals;

    // Same ordering rule as the torsion proposals: energy estimate and coverage, mixed.
    std::vector<const std::vector<int>*> seen;
    for (const Frame& f : m_frames)
        seen.push_back(&f.states);
    for (const std::vector<int>& v : m_proposed_before)
        seen.push_back(&v);
    for (Proposal& p : proposals) {
        int best = std::numeric_limits<int>::max();
        for (const std::vector<int>* v : seen) {
            if (v->size() != p.states.size())
                continue;
            int d = 0;
            for (std::size_t i = 0; i < v->size() && d < best; ++i)
                d += ((*v)[i] != p.states[i]);
            best = std::min(best, d);
        }
        p.novelty = (best == std::numeric_limits<int>::max()) ? 0 : best;
    }
    const double w = (m_proposal_ranking == "energy") ? 0.0
        : (m_proposal_ranking == "coverage") ? 1.0
        : std::max(0.0, std::min(1.0, m_proposal_novelty_weight));
    std::stable_sort(proposals.begin(), proposals.end(), [w](const Proposal& a, const Proposal& b) {
        return (1.0 - w) * a.predicted - w * a.novelty < (1.0 - w) * b.predicted - w * b.novelty;
    });
    if (static_cast<int>(proposals.size()) > m_crossover_max)
        proposals.resize(m_crossover_max);
    return proposals;
}

/* Claude Generated (Aug 2026): COLLECTIVE MODES -- displace along the principal components of the
 * ensemble, computed IN TORSION SPACE.
 *
 * The crossover move transfers folds the ensemble already contains. This one is the complement: it
 * is not restricted to observed combinations. The measured problem is that the reference conformer
 * sits 6-7 torsions away and that six or seven torsions do not move independently -- a backbone
 * rearranges collectively. The leading eigenvectors of the ensemble's own covariance ARE those
 * collective directions, and they cost nothing because the structures are already there.
 *
 * Why torsion space and not Cartesian space: measured, a displacement of 2 sigma along a Cartesian
 * principal component had ALL 20 proposals rejected before optimisation -- a linear displacement of
 * coordinates does not preserve bond lengths, so it produces stretched bonds and clashes rather than
 * a conformer. In torsion space bonds and angles are preserved by construction.
 *
 * Periodicity is handled by working on (cos, sin) of each dihedral instead of the angle itself; a
 * PCA over raw angles would be dominated by the 359 -> 1 degree wrap. The displaced angles are then
 * snapped to the nearest populated state, so the result is a legal state vector that the ordinary
 * build path can realise -- the mode supplies the DIRECTION, the state table the vocabulary.
 */
std::vector<ConfGen::Proposal> ConfGen::generateModeProposals() const
{
    std::vector<Proposal> proposals;
    const std::vector<int> informative = informativeTorsions();
    if (m_mode_max <= 0 || m_frames.size() < 4 || informative.empty())
        return proposals;

    std::vector<int> order;
    for (std::size_t i = 0; i < m_frames.size(); ++i)
        if (m_frames[i].from_input)
            order.push_back(static_cast<int>(i));
    if (order.empty())
        for (std::size_t i = 0; i < m_frames.size(); ++i)
            order.push_back(static_cast<int>(i));
    std::sort(order.begin(), order.end(), [this](int a, int b) { return m_frames[a].energy < m_frames[b].energy; });

    const int T = static_cast<int>(m_torsions.size());
    const int dim = 2 * T;   // (cos, sin) per torsion
    std::vector<Eigen::VectorXd> rows;
    rows.reserve(m_frames.size());
    for (const Frame& f : m_frames) {
        if (static_cast<int>(f.angles.size()) != T)
            continue;
        Eigen::VectorXd v(dim);
        for (int t = 0; t < T; ++t) {
            const double a = f.angles[t] * M_PI / 180.0;
            v(2 * t) = std::cos(a);
            v(2 * t + 1) = std::sin(a);
        }
        rows.push_back(std::move(v));
    }
    if (static_cast<int>(rows.size()) < 4)
        return proposals;

    Eigen::VectorXd mean = Eigen::VectorXd::Zero(dim);
    for (const auto& v : rows)
        mean += v;
    mean /= static_cast<double>(rows.size());
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(dim, dim);
    for (const auto& v : rows) {
        const Eigen::VectorXd d = v - mean;
        cov.noalias() += d * d.transpose();
    }
    cov /= static_cast<double>(rows.size() - 1);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(cov);
    if (solver.info() != Eigen::Success)
        return proposals;
    const Eigen::VectorXd& lambda = solver.eigenvalues();    // ascending
    const Eigen::MatrixXd& vectors = solver.eigenvectors();
    const double total = lambda.sum();
    if (total <= 0)
        return proposals;

    std::set<std::vector<int>> known;
    for (const Frame& f : m_frames)
        known.insert(f.states);
    known.insert(m_proposed_before.begin(), m_proposed_before.end());
    std::set<std::vector<int>> seen_here;

    const int modes = std::max(1, std::min(m_mode_count, dim));
    const int n_templates = std::min<int>(std::max(1, m_proposal_templates), static_cast<int>(order.size()));
    for (int t = 0; t < n_templates; ++t) {
        const Frame& tmpl = m_frames[order[t]];
        if (static_cast<int>(tmpl.angles.size()) != T)
            continue;
        Eigen::VectorXd base(dim);
        for (int k = 0; k < T; ++k) {
            const double a = tmpl.angles[k] * M_PI / 180.0;
            base(2 * k) = std::cos(a);
            base(2 * k + 1) = std::sin(a);
        }
        for (int m = 0; m < modes; ++m) {
            const int col = dim - 1 - m;
            const double sigma = std::sqrt(std::max(0.0, lambda(col)));
            if (sigma < 1e-8)
                continue;
            for (int sign = -1; sign <= 1; sign += 2) {
                Eigen::VectorXd moved = base + sign * m_mode_amplitude * sigma * vectors.col(col);
                std::vector<int> states = tmpl.states;
                int changed = 0;
                for (int k = 0; k < T; ++k) {
                    if (m_state_centres[k].size() < 2)
                        continue;
                    double angle = std::atan2(moved(2 * k + 1), moved(2 * k)) * 180.0 / M_PI;
                    const int state = TorsionSpace::assignState(angle, m_state_centres[k]);
                    if (state >= 0 && state != states[k]) {
                        states[k] = state;
                        changed++;
                    }
                }
                if (changed == 0 || known.count(states) || seen_here.count(states))
                    continue;
                seen_here.insert(states);
                Proposal p;
                p.states = states;
                p.template_frame = order[t];
                p.distance = changed;
                p.nci_label = fmt::format("collective mode {} ({:.0f} % of the torsion variance), {}{:.1f} sigma",
                    m + 1, 100.0 * lambda(col) / total, sign > 0 ? "+" : "-", m_mode_amplitude);
                proposals.push_back(std::move(p));
            }
        }
    }
    if (static_cast<int>(proposals.size()) > m_mode_max)
        proposals.resize(m_mode_max);
    return proposals;
}

/* Claude Generated (Aug 2026): PATH IMAGES -- propose the structures that lie BETWEEN two known
 * conformers instead of around one.
 *
 * Every move set so far starts at one structure and steps away from it. The measured gap is
 * different in kind: the reference conformer is 6-7 torsions from the nearest of ~4600 sampled
 * structures, and a step of that size does not survive optimisation (measured: not one deep-mutation
 * proposal kept its target state vector). But two ensemble members that are themselves 6-7 torsions
 * APART have that same region between them -- and a point halfway along is only 3-4 torsions from
 * either end, which is a distance the machinery handles well.
 *
 * The path here is the straight line in torsion space (circular interpolation per dihedral, snapped
 * to the populated states), not a minimum-energy path. That is deliberate: it costs nothing and
 * tests the hypothesis -- do the missing structures lie between the known ones? -- before the NEB
 * machinery is brought in to replace the straight line with the real path.
 */
std::vector<ConfGen::Proposal> ConfGen::generatePathProposals() const
{
    std::vector<Proposal> proposals;
    if (m_path_max <= 0 || m_frames.size() < 2)
        return proposals;

    std::vector<int> order;
    for (std::size_t i = 0; i < m_frames.size(); ++i)
        if (m_frames[i].from_input)
            order.push_back(static_cast<int>(i));
    if (order.empty())
        for (std::size_t i = 0; i < m_frames.size(); ++i)
            order.push_back(static_cast<int>(i));
    std::sort(order.begin(), order.end(), [this](int a, int b) { return m_frames[a].energy < m_frames[b].energy; });

    std::set<std::vector<int>> known;
    for (const Frame& f : m_frames)
        known.insert(f.states);
    known.insert(m_proposed_before.begin(), m_proposed_before.end());
    std::set<std::vector<int>> seen_here;

    const int T = static_cast<int>(m_torsions.size());
    const int images = std::max(1, m_path_images);
    const int n_templates = std::min<int>(std::max(1, m_proposal_templates), static_cast<int>(order.size()));
    int pairs = 0;
    for (int a = 0; a < n_templates; ++a) {
        const Frame& A = m_frames[order[a]];
        if (static_cast<int>(A.angles.size()) != T)
            continue;
        for (int idx : order) {
            const Frame& B = m_frames[idx];
            if (idx == order[a] || static_cast<int>(B.angles.size()) != T)
                continue;
            const int gap = TorsionSpace::hammingDistance(A.states, B.states);
            if (gap < m_path_min_distance)
                continue;
            pairs++;
            for (int im = 1; im <= images; ++im) {
                const double f = static_cast<double>(im) / (images + 1);
                std::vector<int> states = (f < 0.5) ? A.states : B.states;
                const int template_frame = (f < 0.5) ? order[a] : idx;
                int changed = 0;
                for (int t = 0; t < T; ++t) {
                    if (m_state_centres[t].size() < 2)
                        continue;
                    // circular interpolation along the shorter arc
                    double d = B.angles[t] - A.angles[t];
                    while (d > 180.0) d -= 360.0;
                    while (d < -180.0) d += 360.0;
                    const double angle = A.angles[t] + f * d;
                    const int state = TorsionSpace::assignState(angle, m_state_centres[t]);
                    if (state >= 0 && state != states[t]) {
                        states[t] = state;
                        changed++;
                    }
                }
                if (changed == 0 || known.count(states) || seen_here.count(states))
                    continue;
                seen_here.insert(states);
                Proposal p;
                p.states = states;
                p.template_frame = template_frame;
                p.distance = changed;
                p.nci_label = fmt::format("path image {}/{} between structures {} and {} ({} torsions apart)",
                    im, images, order[a] + 1, idx + 1, gap);
                proposals.push_back(std::move(p));
            }
        }
    }
    if (proposals.empty())
        return proposals;

    std::vector<const std::vector<int>*> seen;
    for (const Frame& f : m_frames)
        seen.push_back(&f.states);
    for (const std::vector<int>& v : m_proposed_before)
        seen.push_back(&v);
    for (Proposal& p : proposals) {
        int best = std::numeric_limits<int>::max();
        for (const std::vector<int>* v : seen) {
            if (v->size() != p.states.size())
                continue;
            int d = 0;
            for (std::size_t i = 0; i < v->size() && d < best; ++i)
                d += ((*v)[i] != p.states[i]);
            best = std::min(best, d);
        }
        p.novelty = (best == std::numeric_limits<int>::max()) ? 0 : best;
    }
    // Coverage first here: the point of a path image is to sit where nothing has been.
    std::stable_sort(proposals.begin(), proposals.end(),
        [](const Proposal& x, const Proposal& y) { return x.novelty > y.novelty; });
    if (static_cast<int>(proposals.size()) > m_path_max)
        proposals.resize(m_path_max);
    return proposals;
}


/* Claude Generated (Aug 2026): the CONCERTED move -- one torsion and one hydrogen bond changed in
 * the SAME restrained optimisation.
 *
 * Until now this existed only as a parameter and a promise: -concerted_max was read and never used,
 * and the two move sets ran strictly separately. That separation is the problem they were supposed
 * to solve. A torsion move ignores the term that carries most of the energy spread; an NCI move
 * re-ties a bridge but leaves the backbone where it was, and measured, it then holds its target
 * pattern in 7 of 13 cases while landing in an already known minimum in 6 of those 7 -- because
 * changing one bond does not by itself change the basin. The backbone has to follow.
 *
 * The coupling is geometric, not statistical: a torsion is chosen whose rotation moves EXACTLY ONE
 * of the two bridge partners (its `moving` set contains one of them). Rotating it therefore changes
 * their relative position, which is precisely what a re-tied bridge needs. A torsion that moves both
 * partners or neither would leave the bridge geometry untouched and reduce the move to the two
 * separate ones again.
 *
 * The build needs no new code: restrainedBuildNCI already applies the dihedral restraints of a
 * changed state vector together with the distance restraints of the bridge -- it just never received
 * a proposal that carried both.
 */
std::vector<ConfGen::Proposal> ConfGen::generateConcertedProposals() const
{
    std::vector<Proposal> proposals;
    if (m_concerted_max <= 0 || m_nci_pairs.empty())
        return proposals;
    const std::vector<Proposal> nci = generateNCIProposals();
    if (nci.empty())
        return proposals;

    std::set<std::vector<int>> known;
    for (const Frame& f : m_frames)
        known.insert(f.states);
    known.insert(m_proposed_before.begin(), m_proposed_before.end());

    for (const Proposal& base : nci) {
        if (static_cast<int>(proposals.size()) >= m_concerted_max)
            break;
        if (base.template_frame < 0 || base.template_frame >= static_cast<int>(m_frames.size()))
            continue;
        // The heavy partners of every bridge this move touches.
        std::vector<int> partners;
        for (const auto& [k, on] : base.nci_targets) {
            if (k < 0 || k >= static_cast<int>(m_nci_pairs.size()))
                continue;
            partners.push_back(m_nci_pairs[k].first);
            partners.push_back(m_nci_pairs[k].second);
        }
        if (partners.size() < 2)
            continue;

        // A torsion is coupled to the bridge when rotating it moves one partner and not the other.
        int chosen = -1, chosen_state = -1, best_population = -1;
        for (std::size_t t = 0; t < m_torsions.size(); ++t) {
            if (t >= m_state_centres.size() || m_state_centres[t].size() < 2)
                continue;
            const std::vector<int>& moving = m_torsions[t].moving;
            bool couples = false;
            for (std::size_t k = 0; k + 1 < partners.size(); k += 2) {
                const bool a = std::find(moving.begin(), moving.end(), partners[k]) != moving.end();
                const bool b = std::find(moving.begin(), moving.end(), partners[k + 1]) != moving.end();
                if (a != b) { couples = true; break; }
            }
            if (!couples)
                continue;
            // Target state: the most populated one that is NOT the template's -- a state the
            // ensemble actually visits, so the restraint pulls towards a real basin.
            const int current = (t < base.states.size()) ? base.states[t] : -1;
            std::vector<int> population(m_state_centres[t].size(), 0);
            for (const Frame& f : m_frames)
                if (t < f.states.size() && f.states[t] >= 0 && f.states[t] < static_cast<int>(population.size()))
                    population[f.states[t]]++;
            for (std::size_t st = 0; st < population.size(); ++st) {
                if (static_cast<int>(st) == current)
                    continue;
                if (population[st] > best_population) {
                    best_population = population[st];
                    chosen = static_cast<int>(t);
                    chosen_state = static_cast<int>(st);
                }
            }
        }
        if (chosen < 0 || best_population <= 0)
            continue;

        Proposal p = base;
        p.states[chosen] = chosen_state;
        if (known.count(p.states))
            continue;
        p.distance = 1;   // one torsion away from its template, plus the bridge change
        p.nci_label = fmt::format("{} + torsion {} -> {:.0f} deg (concerted)", base.nci_label,
            m_torsions[chosen].label(m_frames[base.template_frame].molecule),
            m_state_centres[chosen][chosen_state]);
        proposals.push_back(std::move(p));
    }
    return proposals;
}

std::vector<ConfGen::Proposal> ConfGen::generateNCIProposals() const
{
    std::vector<Proposal> proposals;
    // Only hydrogen bonds are moved. They are the directional, chemically specific interactions and
    // the ones the variance attribution points at; a "close contact" has no well-defined target
    // geometry to restrain towards.
    std::vector<int> movable;
    for (std::size_t k = 0; k < m_nci_pairs.size(); ++k)
        if (m_nci_pairs[k].kind == NCIContact::HBond && m_nci_pairs[k].hydrogen >= 0)
            movable.push_back(static_cast<int>(k));
    if (movable.empty())
        return proposals;

    std::vector<int> population(m_nci_pairs.size(), 0);
    for (const Frame& f : m_frames)
        for (std::size_t k = 0; k < m_nci_pairs.size() && k < f.nci.size(); ++k)
            population[k] += f.nci[k];

    // Templates: the same lowest-energy members the torsion stage uses, so both move sets start from
    // the same structures and their yields are comparable.
    std::vector<int> order;
    for (std::size_t i = 0; i < m_frames.size(); ++i)
        if (m_frames[i].from_input)
            order.push_back(static_cast<int>(i));
    if (order.empty())
        for (std::size_t i = 0; i < m_frames.size(); ++i)
            order.push_back(static_cast<int>(i));
    std::sort(order.begin(), order.end(), [this](int a, int b) { return m_frames[a].energy < m_frames[b].energy; });
    const int n_templates = std::min<int>(std::max(1, m_proposal_templates), static_cast<int>(order.size()));

    // Every pattern the ensemble already has -- a move that lands on one of them proposes nothing new.
    std::set<std::vector<int>> known;
    for (const Frame& f : m_frames)
        known.insert(f.nci);
    std::set<std::vector<int>> proposed;

    for (int t = 0; t < n_templates; ++t) {
        const int frame = order[t];
        const std::vector<int>& have = m_frames[frame].nci;
        // Depth 1: every single flip. Depth 2: one break together with one form -- the concerted move
        // that a hydrogen bond actually performs, since an opened donor usually finds a new acceptor.
        std::vector<int> present, absent;
        for (int k : movable)
            (k < static_cast<int>(have.size()) && have[k] ? present : absent).push_back(k);

        auto emit = [&](const std::vector<std::pair<int, int>>& targets) {
            std::vector<int> pattern = have;
            for (const auto& [k, on] : targets)
                if (k < static_cast<int>(pattern.size()))
                    pattern[k] = on;
            if (known.count(pattern) || proposed.count(pattern))
                return;
            proposed.insert(pattern);
            Proposal p;
            p.states = m_frames[frame].states; // unchanged: this move does not touch the torsions
            p.template_frame = frame;
            p.distance = static_cast<int>(targets.size());
            p.nci_targets = targets;
            for (const auto& [k, on] : targets)
                p.nci_label += fmt::format("{}{} {}", p.nci_label.empty() ? "" : ", ",
                    on ? "form" : "break", m_nci_pairs[k].label());
            proposals.push_back(std::move(p));
        };

        for (int k : present)
            emit({ { k, 0 } });
        for (int k : absent)
            emit({ { k, 1 } });
        if (m_nci_depth >= 2)
            for (int b : present)
                for (int f : absent)
                    emit({ { b, 0 }, { f, 1 } });
    }

    // Ranking without an energy model: prefer moves whose TARGET bond is well populated in the
    // ensemble (a bond many structures realise is a plausible one to ask for) and, at equal
    // population, the smaller move. Deliberately not the additive model -- it was measured to
    // explain ~15 % of the energy variation and has no say over the NCI pattern at all.
    // Claude Generated (Aug 2026): same two-sided ordering as the torsion proposals. The population
    // term is the "is this pattern plausible" half -- a bond many structures form is easier to form
    // again -- and it is exactly the half that was measured to drive TOWARDS over-bridged structures
    // and AWAY from the reference (9 bridges built, the reference has 6). The coverage term is the
    // counterweight: distance in contact space to every pattern the run has already produced.
    {
        std::vector<double> pop(proposals.size()), nov(proposals.size());
        for (std::size_t i = 0; i < proposals.size(); ++i) {
            int s = 0;
            for (const auto& [k, on] : proposals[i].nci_targets)
                s += on ? population[k] : -population[k];
            pop[i] = -static_cast<double>(s);       // klein ist gut, daher Vorzeichen umdrehen
            // Muster dieses Vorschlags, und sein Abstand zu allem bereits Beobachteten
            std::vector<int> pattern = m_frames[proposals[i].template_frame].nci;
            for (const auto& [k, on] : proposals[i].nci_targets)
                if (k < static_cast<int>(pattern.size()))
                    pattern[k] = on;
            int best = std::numeric_limits<int>::max();
            for (const Frame& f : m_frames) {
                if (f.nci.size() != pattern.size())
                    continue;
                int d = 0;
                for (std::size_t j = 0; j < pattern.size() && d < best; ++j)
                    d += (f.nci[j] != pattern[j]);
                best = std::min(best, d);
            }
            proposals[i].novelty = (best == std::numeric_limits<int>::max()) ? 0 : best;
            nov[i] = proposals[i].novelty;
        }
        auto zscore = [](std::vector<double> v) {
            const double mean = v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
            double var = 0.0;
            for (double x : v) var += (x - mean) * (x - mean);
            const double sd = v.size() > 1 ? std::sqrt(var / (v.size() - 1)) : 0.0;
            for (double& x : v) x = (sd > 1e-12) ? (x - mean) / sd : 0.0;
            return v;
        };
        pop = zscore(pop); nov = zscore(nov);
        const double w = (m_proposal_ranking == "energy") ? 0.0
            : (m_proposal_ranking == "coverage") ? 1.0
            : std::max(0.0, std::min(1.0, m_proposal_novelty_weight));
        std::vector<std::size_t> order(proposals.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return (1.0 - w) * pop[a] - w * nov[a] < (1.0 - w) * pop[b] - w * nov[b];
        });
        std::vector<Proposal> sorted;
        sorted.reserve(proposals.size());
        for (std::size_t i : order)
            sorted.push_back(std::move(proposals[i]));
        proposals = std::move(sorted);
    }
    if (static_cast<int>(proposals.size()) > m_nci_max_proposals)
        proposals.resize(m_nci_max_proposals);
    return proposals;
}

bool ConfGen::restrainedBuildNCI(const Proposal& p, Molecule& driven) const
{
    if (!m_calculator || p.nci_targets.empty())
        return false;

    std::vector<Optimization::DistanceRestraint> restraints;
    for (const auto& [k, on] : p.nci_targets) {
        if (k < 0 || k >= static_cast<int>(m_nci_pairs.size()))
            continue;
        const NCIContact& c = m_nci_pairs[k];
        if (c.hydrogen < 0)
            continue;
        Optimization::DistanceRestraint r;
        r.i = c.hydrogen; // the H...acceptor distance is what the detection criterion uses
        r.j = c.second;
        r.target = on ? m_nci_form_distance : m_nci_break_distance;
        r.force = m_nci_restraint_force;
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
    // Claude Generated (Aug 2026): ONE shared calculator for the whole run (see m_calculator) --
    // the optimiser must not re-derive the force field per structure. Saves a full GFN-FF setup
    // per proposal and keeps the run out of the intermittent parameter-generation crash.
    opt_config["reuse_calculator"] = true;
    opt_config["write_trajectory"] = false;
    opt_config["max_iterations"] = m_restraint_max_iterations;
    opt_config["distance_restraints"] = Optimization::GeometryRestraints::toJson(restraints);

    // Claude Generated (Aug 2026): CONCERTED move. When the proposal also changes torsion states
    // relative to its template, their dihedral restraints act in the SAME optimisation as the
    // distance restraints. The two move sets otherwise pull one after the other and each undoes part
    // of the other's work: a new hydrogen bond usually needs the backbone to turn, and a turned
    // backbone usually breaks the old bond. Only both together express what actually happens.
    {
        std::vector<Optimization::DihedralRestraint> dihedrals;
        const std::vector<int>& from = m_frames[p.template_frame].states;
        for (std::size_t k = 0; k < m_torsions.size() && k < p.states.size(); ++k) {
            if (k < from.size() && p.states[k] == from[k])
                continue;
            if (p.states[k] < 0 || p.states[k] >= static_cast<int>(m_state_centres[k].size()))
                continue;
            Optimization::DihedralRestraint d;
            d.i = m_torsions[k].i; d.j = m_torsions[k].j;
            d.k = m_torsions[k].k; d.l = m_torsions[k].l;
            d.target = m_state_centres[k][p.states[k]] * M_PI / 180.0;  // Zustandsmitten in Grad
            d.force = m_restraint_force;
            dihedrals.push_back(d);
        }
        if (!dihedrals.empty())
            opt_config["dihedral_restraints"] = Optimization::GeometryRestraints::toJson(dihedrals);
    }

    Molecule mol = m_frames[p.template_frame].molecule;

    /* Claude Generated (Aug 2026): STAGED realisation of a multi-bridge move.
     *
     * Pulling several distance restraints at once is what the feasibility test did that failed:
     * taking the ensemble structure closest to a reference conformer (2.61 A) and closing all six of
     * its target hydrogen bonds simultaneously left the RMSD at 2.55-2.64 A and the restrained
     * optimisation stalled. Six distances are six conditions in 315 degrees of freedom, and the
     * optimiser satisfies them with a local compromise instead of refolding.
     *
     * Staged means: close ONE bridge, let the molecule relax around it, keep it restrained, add the
     * next. The order is greedy -- always the bridge whose current distance is nearest its target,
     * so each stage is the smallest possible step and the backbone reorganises gradually instead of
     * being torn in six directions at once. The final stage is the same simultaneous optimisation as
     * before, so the acceptance test below is unchanged. */
    if (m_nci_staged && restraints.size() > 1) {
        std::vector<Optimization::DistanceRestraint> pending = restraints, active;
        while (!pending.empty()) {
            const Geometry current = mol.getGeometry();
            std::size_t best = 0;
            double best_gap = std::numeric_limits<double>::max();
            for (std::size_t k = 0; k < pending.size(); ++k) {
                const double d = (Eigen::Vector3d(current.row(pending[k].i))
                    - Eigen::Vector3d(current.row(pending[k].j))).norm();
                const double gap = std::abs(d - pending[k].target);
                if (gap < best_gap) { best_gap = gap; best = k; }
            }
            active.push_back(pending[best]);
            pending.erase(pending.begin() + best);
            json stage = opt_config;
            stage["distance_restraints"] = Optimization::GeometryRestraints::toJson(active);
            auto step = Optimization::OptimizationDispatcher::optimizeStructure(
                &mol, Optimization::OptimizerType::LBFGSPP, m_calculator.get(), stage);
            if (step.final_molecule.AtomCount() == 0)
                return false;
            mol = step.final_molecule;
        }
    }

    auto result = Optimization::OptimizationDispatcher::optimizeStructure(
        &mol, Optimization::OptimizerType::LBFGSPP, m_calculator.get(), opt_config);
    if (result.final_molecule.AtomCount() == 0)
        return false;

    // Did the move arrive? A restraint competes with the force field, so a target that the molecule
    // cannot reach ends up somewhere else -- and then this is not the proposed pattern. Half an
    // Angstrom is the tolerance; the form and break targets are 1.6 A apart, so it cannot confuse
    // the two.
    const Geometry geom = result.final_molecule.getGeometry();
    for (const auto& r : restraints) {
        const double d = (Eigen::Vector3d(geom.row(r.i)) - Eigen::Vector3d(geom.row(r.j))).norm();
        if (std::abs(d - r.target) > 0.5)
            return false;
    }
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
    // Claude Generated (Aug 2026): combinations an earlier call already built count as known -- a
    // proposal that was tried and rejected is not in the ensemble, so without this memory every
    // repetition rebuilds and re-optimises it.
    known.insert(m_proposed_before.begin(), m_proposed_before.end());

    // Templates: the lowest-energy members OF THIS CALL (their geometry is the starting point).
    // Frames contributed by -analysis_file describe, they do not seed -- see Frame::from_input.
    std::vector<int> order;
    for (std::size_t i = 0; i < m_frames.size(); ++i)
        if (m_frames[i].from_input)
            order.push_back(static_cast<int>(i));
    if (order.empty())
        for (std::size_t i = 0; i < m_frames.size(); ++i)
            order.push_back(static_cast<int>(i));
    std::sort(order.begin(), order.end(), [this](int a, int b) { return m_frames[a].energy < m_frames[b].energy; });
    const int n_templates = std::min<int>(std::max(1, m_proposal_templates), static_cast<int>(order.size()));

    const std::vector<std::vector<double>> coefficients = additiveCoefficients();
    std::set<std::vector<int>> proposed;

    // Claude Generated (Aug 2026): how large is the neighbourhood we are about to walk? The number of
    // candidates at exactly Hamming distance d is the elementary symmetric polynomial e_d of the
    // per-torsion alternative counts (how many OTHER states each torsion offers) -- computable in
    // O(n * depth) before a single candidate is built. Measured on a 107-atom peptide with 29
    // torsions: ~1e2 at depth 1, 4e3 at depth 2, 1e5 at depth 3, 3e6 at depth 4 and 5e7 at depth 5.
    // The last one needs about 10 GB for the state vectors alone and aborted with std::bad_alloc,
    // which is why -proposal_depth was effectively capped at 3.
    std::vector<double> alternatives;
    alternatives.reserve(informative.size());
    for (int t : informative)
        alternatives.push_back(std::max<int>(0, static_cast<int>(m_state_centres[t].size()) - 1));
    const int max_depth = std::max(1, m_proposal_depth);
    std::vector<double> shell(max_depth + 1, 0.0);
    shell[0] = 1.0;
    for (double a : alternatives)
        for (int d = max_depth; d >= 1; --d)
            shell[d] += shell[d - 1] * a;
    double ball = 0.0;
    for (int d = 1; d <= max_depth; ++d)
        ball += shell[d];
    const double cap = std::max(1, m_proposal_candidate_cap);
    const bool sampled = ball > cap;

    // One candidate, whichever way it was produced: enumerated or sampled.
    auto record = [&](int template_frame, const std::vector<int>& states, int depth) {
        if (depth <= 0 || known.count(states) || proposed.count(states))
            return false;
        Proposal p;
        p.states = states;
        p.template_frame = template_frame;
        p.distance = depth;
        for (std::size_t t = 0; t < states.size(); ++t)
            if (states[t] >= 0 && states[t] < static_cast<int>(coefficients[t].size()))
                p.predicted += coefficients[t][states[t]];
        proposed.insert(states);
        proposals.push_back(std::move(p));
        return true;
    };

    // Enumerate mutations up to the requested Hamming depth around each template.
    std::function<void(int, std::vector<int>&, int, int)> mutate =
        [&](int template_frame, std::vector<int>& states, int start, int depth) {
            record(template_frame, states, depth);
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

    // Claude Generated (Aug 2026): the same neighbourhood, drawn instead of walked. Depth is drawn
    // proportional to the shell sizes and the torsions proportional to how many alternatives they
    // offer, so the sample follows the ball it replaces rather than favouring the near shells. The
    // draws are independent, hence duplicates -- they cost one set lookup and are dropped.
    const unsigned int seed = (m_proposal_seed != 0)
        ? static_cast<unsigned int>(m_proposal_seed)
        : static_cast<unsigned int>(1234567u + 977u * m_frames.size() + 131u * m_proposed_before.size());
    std::mt19937 rng(seed);
    std::discrete_distribution<int> depth_dist(shell.begin() + 1, shell.end());
    auto draw = [&](int template_frame, const std::vector<int>& base, int budget) {
        std::vector<double> weight;
        for (int attempt = 0; attempt < budget; ++attempt) {
            const int d = depth_dist(rng) + 1;
            weight.assign(alternatives.begin(), alternatives.end());
            std::vector<int> states = base;
            int changed = 0;
            for (int k = 0; k < d; ++k) {
                const double total = std::accumulate(weight.begin(), weight.end(), 0.0);
                if (total <= 0.0)
                    break;
                double x = std::uniform_real_distribution<double>(0.0, total)(rng);
                std::size_t idx = 0;
                for (; idx + 1 < weight.size(); ++idx) {
                    x -= weight[idx];
                    if (x <= 0.0)
                        break;
                }
                weight[idx] = 0.0;                      // each torsion at most once per candidate
                const int t = informative[idx];
                const int n_states = static_cast<int>(m_state_centres[t].size());
                if (n_states < 2 || t >= static_cast<int>(states.size()))
                    continue;
                const int original = states[t];
                int st = 0;
                if (original >= 0 && original < n_states) {
                    st = std::uniform_int_distribution<int>(0, n_states - 2)(rng);
                    if (st >= original)
                        st += 1;                         // skip the state the template already has
                } else {
                    st = std::uniform_int_distribution<int>(0, n_states - 1)(rng);
                }
                states[t] = st;
                changed++;
            }
            record(template_frame, states, changed);
        }
    };

    // Claude Generated (Aug 2026): a template whose whole neighbourhood has already been built
    // contributes nothing but consumes one of the few template slots. Since the deepest structures of
    // a cycle are largely the same ones from repetition to repetition, the enumeration otherwise
    // drills the same Hamming ball again and again while the memory quietly rejects every candidate.
    // Walk the energy ranking and count only the templates that actually produced something.
    const int per_template = std::max(1, static_cast<int>(cap) / std::max(1, n_templates));
    int used = 0, skipped = 0;
    for (std::size_t i = 0; i < order.size() && used < n_templates; ++i) {
        const std::size_t before = proposals.size();
        std::vector<int> states = m_frames[order[i]].states;
        if (sampled)
            draw(order[i], states, per_template);
        else
            mutate(order[i], states, 0, 0);
        if (proposals.size() > before)
            used++;
        else
            skipped++;
    }
    if (m_verbosity >= 1) {
        if (sampled)
            CurcumaLogger::result_fmt("ConfGen: the neighbourhood up to depth {} holds {:.3g} combinations -- "
                                      "too many to enumerate; {} distinct unknown ones were drawn from it "
                                      "(seed {}), the bound is -proposal_candidate_cap",
                max_depth, ball, proposals.size(), seed);
        else if (m_verbosity >= 2)
            CurcumaLogger::info_fmt("ConfGen: neighbourhood up to depth {} enumerated exactly ({:.0f} combinations, "
                                    "{} of them not yet known)",
                max_depth, ball, proposals.size());
    }
    if (skipped > 0 && m_verbosity >= 1)
        CurcumaLogger::result_fmt("ConfGen: {} template(s) skipped -- every combination in their "
                                  "neighbourhood had already been built; {} productive template(s) used",
            skipped, used);

    // Claude Generated (Aug 2026): order by energy AND by coverage. The two properties of the
    // description are not equally good: predicting the energy from it works poorly (cross-validated
    // ~30 % of the spread; a delta model between the two surfaces even makes the ranking worse),
    // while separating structures works very well (141 of 142 carry a distinct contact pattern).
    // Ordering by the model alone therefore uses the weak property and ignores the strong one --
    // but the terms DO relate to the final energy, so dropping them would be equally one-sided.
    // Both contributions are standardised over the candidate set and mixed; the old behaviour is
    // -proposal_ranking energy.
    {
        // "already seen" = every state vector of the description ensemble plus everything the run
        // has tried before (the memory), so a later repetition does not re-explore the same shell.
        std::vector<const std::vector<int>*> seen;
        seen.reserve(m_frames.size() + m_proposed_before.size());
        for (const Frame& f : m_frames)
            seen.push_back(&f.states);
        for (const std::vector<int>& v : m_proposed_before)
            seen.push_back(&v);
        for (Proposal& p : proposals) {
            int best = std::numeric_limits<int>::max();
            for (const std::vector<int>* v : seen) {
                if (v->size() != p.states.size())
                    continue;
                int d = 0;
                for (std::size_t i = 0; i < v->size() && d < best; ++i)
                    d += ((*v)[i] != p.states[i]);
                best = std::min(best, d);
            }
            p.novelty = (best == std::numeric_limits<int>::max()) ? 0 : best;
        }
        auto zscore = [](std::vector<double> v) {
            const double mean = v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
            double var = 0.0;
            for (double x : v) var += (x - mean) * (x - mean);
            const double sd = v.size() > 1 ? std::sqrt(var / (v.size() - 1)) : 0.0;
            for (double& x : v) x = (sd > 1e-12) ? (x - mean) / sd : 0.0;
            return v;
        };
        std::vector<double> e, d;
        e.reserve(proposals.size()); d.reserve(proposals.size());
        for (const Proposal& p : proposals) { e.push_back(p.predicted); d.push_back(p.novelty); }
        e = zscore(e); d = zscore(d);
        const double w = (m_proposal_ranking == "energy") ? 0.0
            : (m_proposal_ranking == "coverage") ? 1.0
            : std::max(0.0, std::min(1.0, m_proposal_novelty_weight));
        std::vector<std::size_t> order(proposals.size());
        std::iota(order.begin(), order.end(), 0);
        std::vector<double> score(proposals.size());
        for (std::size_t i = 0; i < proposals.size(); ++i)
            score[i] = (1.0 - w) * e[i] - w * d[i];   // klein ist gut: tiefe Energie, grosser Abstand
        std::stable_sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) { return score[a] < score[b]; });
        std::vector<Proposal> sorted;
        sorted.reserve(proposals.size());
        for (std::size_t i : order)
            sorted.push_back(std::move(proposals[i]));
        proposals = std::move(sorted);
        if (m_verbosity >= 1 && !proposals.empty())
            CurcumaLogger::result_fmt("ConfGen: candidate ordering '{}' (novelty weight {:.2f}); the {} kept "
                                      "proposals have Hamming distance {} to {} from everything seen so far",
                m_proposal_ranking, w, std::min<int>(m_max_proposals, static_cast<int>(proposals.size())),
                std::min_element(proposals.begin(), proposals.begin() + std::min<std::size_t>(m_max_proposals, proposals.size()),
                    [](const Proposal& a, const Proposal& b) { return a.novelty < b.novelty; })->novelty,
                std::max_element(proposals.begin(), proposals.begin() + std::min<std::size_t>(m_max_proposals, proposals.size()),
                    [](const Proposal& a, const Proposal& b) { return a.novelty < b.novelty; })->novelty);
    }
    if (static_cast<int>(proposals.size()) > m_max_proposals)
        proposals.resize(m_max_proposals);
    return proposals;
}


/* Claude Generated (Aug 2026): the calculator that JUDGES a proposal, which need not be the one that
 * DESCRIBES the ensemble.
 *
 * ConfGen used one method for everything: the per-term decomposition (which only a force field
 * provides), the geometry building, and the optimisation + energy that decide whether a proposal is
 * any good. That mixes two surfaces which were measured NOT to agree -- within one cycle of a
 * 107-atom peptide GFN-FF and GFN2 rank the same structures at r = -0.32 ... -0.46, and every move
 * set built here was therefore selected on a surface that points elsewhere. The clearest case: the
 * collective-mode move produced the deepest GFN-FF structure this generator ever made (20.4 kJ/mol
 * below the ensemble minimum) and the same structures sit 87 to 134 kJ/mol above the reference on
 * GFN2.
 *
 * -eval_method splits the two. The description stays on the force field, so the term decomposition,
 * the matched pairs and the additive model keep working; the proposals are optimised and compared on
 * the accurate surface, which is the one that decides. Empty = same method for both (old behaviour).
 */
bool ConfGen::useSplitEvaluation() const
{
    if (m_eval_method.empty() || m_eval_method == m_method)
        return false;
    evaluationCalculator(); // creates and PROBES the calculator once; may set m_eval_unusable
    return !m_eval_unusable;
}

EnergyCalculator* ConfGen::evaluationCalculator() const
{
    if (m_eval_method.empty() || m_eval_method == m_method || m_eval_unusable)
        return m_calculator.get();
    if (!m_eval_calculator) {
        json cfg;
        cfg["method"] = m_eval_method;
        cfg["threads"] = m_threads;
        cfg["charge"] = m_charge;
        cfg["spin"] = m_spin;
        cfg["verbosity"] = 0;
        try {
            m_eval_calculator = std::make_unique<EnergyCalculator>(m_eval_method, cfg);
            if (!m_frames.empty()) {
                Molecule ref = m_frames.front().molecule;
                ref.setCharge(m_charge);
                ref.setSpin(m_spin);
                m_eval_calculator->setMolecule(ref.getMolInfo());
                /* Claude Generated (Aug 2026): PROVE the calculator works before every proposal is
                 * staked on it. An unusable method name does not necessarily throw -- it can return
                 * a finite garbage energy with a NaN gradient, and then the optimiser refuses to
                 * start for EVERY proposal while the phase reports "no new conformer". One energy +
                 * gradient on the reference structure costs nothing next to one optimisation per
                 * proposal and converts a silent total loss into a warning plus the old behaviour. */
                const double e_probe = m_eval_calculator->CalculateEnergy(true);
                const Matrix g_probe = m_eval_calculator->Gradient();
                if (!std::isfinite(e_probe) || !g_probe.allFinite() || m_eval_calculator->Error()) {
                    CurcumaLogger::error(fmt::format(
                        "ConfGen: the evaluation method '{}' does not produce a usable energy and "
                        "gradient (E = {}, gradient finite: {}) -- falling back to '{}'. Every "
                        "proposal would otherwise fail its optimisation and the phase would report "
                        "'no new conformer'.",
                        m_eval_method, e_probe, g_probe.allFinite() ? "yes" : "no", m_method));
                    m_eval_calculator.reset();
                    m_eval_unusable = true;
                    return m_calculator.get();
                }
            }
            CurcumaLogger::result_fmt("ConfGen: proposals are optimised and judged with '{}' while the "
                                      "description stays on '{}' -- the two surfaces disagree, and the "
                                      "one that decides should be the accurate one",
                m_eval_method, m_method);
        } catch (...) {
            CurcumaLogger::warn_fmt("ConfGen: could not create the evaluation calculator '{}' -- falling "
                                    "back to '{}' for the proposals as well", m_eval_method, m_method);
            m_eval_calculator.reset();
            m_eval_unusable = true;
            return m_calculator.get();
        }
    }
    return m_eval_calculator ? m_eval_calculator.get() : m_calculator.get();
}

void ConfGen::optimiseProposals(std::vector<Proposal>& proposals) const
{
    const bool split = useSplitEvaluation();
    json calc_config;
    calc_config["method"] = split ? m_eval_method : m_method;
    calc_config["threads"] = m_threads;
    calc_config["charge"] = m_charge;
    calc_config["spin"] = m_spin;

    json opt_config = calc_config;
    opt_config["verbosity"] = 0;
    // Claude Generated (Aug 2026): ONE shared calculator for the whole run (see m_calculator) --
    // the optimiser must not re-derive the force field per structure. Saves a full GFN-FF setup
    // per proposal and keeps the run out of the intermittent parameter-generation crash.
    opt_config["reuse_calculator"] = true;
    opt_config["write_trajectory"] = false;

    if (!m_calculator)
        return; // analyseEnsemble must have run first
    EnergyCalculator* judge = evaluationCalculator();
    if (!judge)
        return;
    EnergyCalculator& calculator = *judge;
    const std::vector<std::pair<int, int>> reference_bonds
        = topologyFingerprint(m_frames.front().molecule, m_topology_factor);

    /* Claude Generated (Aug 2026): the optimisations run in PARALLEL, the checks and the logging
     * afterwards do not.
     *
     * They used to be strictly serial, because ConfGen shares ONE calculator across the whole run
     * (reuse_calculator, introduced so the force field is not re-derived per structure and to dodge
     * an intermittent GFN-FF parametrisation crash) -- and one calculator cannot serve several
     * threads. Measured cost of that: one repetition of a 107-atom peptide spent 28 minutes on 30
     * proposals while ten threads idled, and over 35 repetitions the production run put ~16 of its
     * 105 hours into this one phase. It also inflated the +29 % that the ConfGen A/B measured.
     *
     * The fix keeps the property that mattered: each WORKER gets its own calculator and reuses it
     * across its own share of the proposals, so nothing is re-derived per structure -- there are
     * simply W independent chains instead of one. Each optimisation then runs single-threaded,
     * which is the same trade ConfSearch already makes for its optimisation batches ("single
     * threaded per optimisation when the caller parallelises externally").
     *
     * Only the optimisation is parallel. Topology check, novelty and the log lines stay in the
     * original serial loop below: they touch shared state (m_frames, the logger) and they are
     * cheap, so there is nothing to win and a race to lose. */
    std::vector<int> todo;
    for (int i = 0; i < static_cast<int>(proposals.size()); ++i)
        if (proposals[i].geometry.AtomCount() > 0)
            todo.push_back(i);
    const int workers = std::max(1, std::min<int>(m_threads, static_cast<int>(todo.size())));
    if (workers > 1) {
        json worker_config = opt_config;
        worker_config["threads"] = 1;   // W structures at once, one thread each
        const auto t_start = std::chrono::steady_clock::now();
        /* Claude Generated (Aug 2026): the child optimisations run at verbosity 0, and the logger
         * level is a shared static (see the known issue in the ConfSearch docs) -- they leave it at
         * 0 and every message after the join is swallowed. Measured: the summary below simply never
         * appeared, which is what made this parallelisation look like it was not running at all.
         * Save and restore around the batch, exactly as ConfSearch::PerformOptimisation does. */
        const int verbosity_before = CurcumaLogger::get_verbosity();
        std::vector<std::thread> pool;
        for (int w = 0; w < workers; ++w) {
            pool.emplace_back([&, w]() {
                // One calculator per worker, reused across its whole share -- the point of
                // reuse_calculator, preserved.
                json cfg = calc_config;
                cfg["threads"] = 1;
                EnergyCalculator local(split ? m_eval_method : m_method, cfg);
                if (!m_frames.empty())
                    local.setMolecule(m_frames.front().molecule.getMolInfo());
                for (std::size_t idx = w; idx < todo.size(); idx += workers) {
                    Proposal& q = proposals[todo[idx]];
                    Molecule mol = q.geometry;
                    auto r = Optimization::OptimizationDispatcher::optimizeStructure(
                        &mol, Optimization::OptimizerType::LBFGSPP, &local, worker_config);
                    if (!r.success)
                        continue;
                    q.optimised = true;
                    q.geometry = r.final_molecule;
                    q.energy = r.final_energy;
                }
            });
        }
        for (auto& t : pool)
            t.join();
        CurcumaLogger::set_verbosity(verbosity_before);
        CurcumaLogger::result_fmt("ConfGen: {} proposal(s) optimised on {} worker(s) in {:.1f} s "
                                  "(one calculator each, reused across its share)",
            static_cast<int>(todo.size()), workers,
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count());
    }

    const auto serial_start = std::chrono::steady_clock::now();
    int serial_done = 0;
    /* Claude Generated (Aug 2026): a failed proposal optimisation used to be a bare `continue` --
     * counted nowhere, reported nowhere. The proposal then simply never appears in
     * <base>.proposals.opt.xyz, which looks exactly like "the generator had nothing to offer".
     * Measured: a production run built 89 proposals in three cycles and optimised none of them,
     * and the only visible message was "no new conformer this cycle". The first failure now names
     * its reason, and the total is reported below. */
    int opt_failed = 0;
    std::string first_failure;
    for (Proposal& p : proposals) {
        if (p.geometry.AtomCount() == 0)
            continue; // rejected by the clash filter, never built
        if (workers > 1) {
            if (!p.optimised) {
                ++opt_failed;
                continue; // already done in parallel above (or failed there)
            }
        } else {
            Molecule mol = p.geometry;
            auto result = Optimization::OptimizationDispatcher::optimizeStructure(
                &mol, Optimization::OptimizerType::LBFGSPP, &calculator, opt_config);
            if (!result.success) {
                ++opt_failed;
                if (first_failure.empty())
                    first_failure = result.error_message.empty() ? "no reason reported"
                                                                 : result.error_message;
                continue;
            }
            p.optimised = true;
            p.geometry = result.final_molecule;
            p.energy = result.final_energy;
            ++serial_done;
        }
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
        /* Claude Generated (Aug 2026): novelty is geometric everywhere else, and for an
         * isomerisation that is the wrong test. Flipping a torsion the ensemble never opened
         * produces a structure in a rotamer state that NO ensemble member occupies -- new by
         * construction -- while moving few enough atoms to stay inside the RMSD threshold.
         * Measured: the guanidinium cis form (+21.8 deg against +174 in all 99 members) was
         * discarded as a duplicate by the RMSD gate. So an isomerisation counts as new when its
         * target torsion actually ARRIVED, i.e. sits at least isomerise_min_separation away from
         * the single state the ensemble knows. The topology gate still applies unchanged. */
        p.is_new = p.topology_ok && p.min_rmsd_to_ensemble > m_new_rmsd;
        if (p.topology_ok && !p.is_new && p.isomerisation() && p.geometry.AtomCount() > 0) {
            const int t = p.angle_targets.front().first;
            if (t >= 0 && t < static_cast<int>(m_torsions.size()) && m_state_centres[t].size() == 1) {
                const double now = TorsionSpace::dihedral(p.geometry.getGeometry(), m_torsions[t]);
                double sep = std::fabs(now - m_state_centres[t][0]);
                if (sep > 180.0) sep = 360.0 - sep;
                if (sep >= m_isomerise_min_separation) {
                    p.is_new = true;
                    CurcumaLogger::result_fmt("ConfGen: isomerisation kept despite RMSD {:.2f} A -- {} now "
                                              "at {:+.0f} deg against {:+.0f} in every input structure, "
                                              "a state the ensemble does not contain",
                        p.min_rmsd_to_ensemble, m_torsions[t].label(m_frames.front().molecule), now,
                        m_state_centres[t][0]);
                }
            }
        }
    }
    /* Claude Generated (Aug 2026): the serial path is timed too, so the two are comparable from a
     * log without rebuilding. It only fires with one worker -- and then each optimisation still has
     * all the threads internally, which is what the old code did. */
    if (workers <= 1 && serial_done > 0)
        CurcumaLogger::result_fmt("ConfGen: {} proposal(s) optimised serially in {:.1f} s "
                                  "(one shared calculator, {} thread(s) inside each optimisation)",
            serial_done, std::chrono::duration<double>(std::chrono::steady_clock::now() - serial_start).count(),
            m_threads);
    if (opt_failed > 0) {
        const int built_total = static_cast<int>(std::count_if(proposals.begin(), proposals.end(),
            [](const Proposal& p) { return p.geometry.AtomCount() > 0; }));
        // ALL of them failing is a different event from a few of them failing: nothing reaches the
        // ensemble, and every downstream message ("no new conformer") is then misleading.
        if (opt_failed == built_total)
            CurcumaLogger::error(fmt::format("ConfGen: ALL {} built proposal(s) failed their optimisation "
                                             "-- this phase produced nothing{}",
                built_total, first_failure.empty() ? "" : ". First reason: " + first_failure));
        else
            CurcumaLogger::warn_fmt("ConfGen: {} of {} built proposal(s) failed their optimisation and "
                                    "were dropped{}",
                opt_failed, built_total,
                first_failure.empty() ? "" : ". First reason: " + first_failure);
    }
}

double ConfGen::referenceEnergyOptimised(double& worst_gain_kJ) const
{
    // The comparison base has to sit on the SAME surface as the proposals, otherwise "x kJ/mol below
    // the ensemble minimum" compares two different energy scales. Claude Generated (Aug 2026).
    const bool split = useSplitEvaluation();
    json calc_config;
    calc_config["method"] = split ? m_eval_method : m_method;
    calc_config["threads"] = m_threads;
    calc_config["charge"] = m_charge;
    calc_config["spin"] = m_spin;
    json opt_config = calc_config;
    opt_config["verbosity"] = 0;
    // Claude Generated (Aug 2026): ONE shared calculator for the whole run (see m_calculator) --
    // the optimiser must not re-derive the force field per structure. Saves a full GFN-FF setup
    // per proposal and keeps the run out of the intermittent parameter-generation crash.
    opt_config["reuse_calculator"] = true;
    opt_config["write_trajectory"] = false;

    if (!m_calculator)
        return std::numeric_limits<double>::infinity();
    EnergyCalculator* judge = evaluationCalculator();
    if (!judge)
        return std::numeric_limits<double>::infinity();
    EnergyCalculator& calculator = *judge;

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
        /* Claude Generated (Aug 2026): the "how far from its own minimum is the input" check needs a
         * baseline on the SAME surface. With -eval_method the frame energy comes from the
         * description method and the optimised energy from the evaluation method -- subtracting them
         * produced a nonsensical 375203 kJ/mol and named the wrong method in the warning. One single
         * point of the judging method on the unoptimised geometry is the correct baseline. */
        double before = m_frames[order[i]].energy;
        if (split) {
            Molecule original = m_frames[order[i]].molecule;
            original.setCharge(m_charge);
            original.setSpin(m_spin);
            calculator.setMolecule(original.getMolInfo());
            const double e = calculator.CalculateEnergy(false);
            before = (calculator.HasNan() || calculator.Error() || !std::isfinite(e)) ? result.final_energy : e;
        }
        worst_gain_kJ = std::max(worst_gain_kJ, (before - result.final_energy) * kEh2kJ);
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

    // Claude Generated (Aug 2026): the two move sets are reported separately. They address different
    // descriptions -- torsions vs the non-covalent pattern -- so a joint yield would hide which one
    // is actually producing, and that is the question the NCI stage exists to answer.
    int nci_total = 0, nci_opt = 0, nci_new = 0;
    double nci_best = std::numeric_limits<double>::infinity();
    std::string nci_best_label;
    for (const Proposal& p : proposals) {
        if (!p.nci_move())
            continue;
        nci_total++;
        if (p.optimised && p.topology_ok)
            nci_opt++;
        if (p.is_new) {
            nci_new++;
            if (p.energy < nci_best) {
                nci_best = p.energy;
                nci_best_label = p.nci_label;
            }
        }
    }
    if (nci_total > 0) {
        CurcumaLogger::result_fmt("ConfGen:   of those, {} came from NCI moves ({} built and valid of {} "
                                  "proposed), {} from torsion recombination",
            nci_new, nci_opt, nci_total, novel - nci_new);
        if (nci_new > 0)
            CurcumaLogger::result_fmt("ConfGen:   best NCI move: {} -> {:.6f} Eh", nci_best_label, nci_best);

        // Did the moved bond SURVIVE the release? This separates the two ways an NCI move can fail:
        // the restrained geometry was never a minimum (the bond snaps back), or it was a minimum but
        // an already known one. Only the first is a defect of the move; the second means the ensemble
        // already covers that pattern. The H-bond criterion is purely geometric, so no charges are
        // needed here.
        int kept = 0, snapped = 0;
        for (const Proposal& p : proposals) {
            if (!p.nci_move() || !p.optimised || !p.topology_ok)
                continue;
            const std::vector<NCIContact> after = detectNCI(p.geometry, {});
            bool all_ok = true;
            for (const auto& [k, on] : p.nci_targets) {
                if (k < 0 || k >= static_cast<int>(m_nci_pairs.size()))
                    continue;
                const NCIContact& target = m_nci_pairs[k];
                const bool present = std::any_of(after.begin(), after.end(), [&](const NCIContact& c) {
                    return c.kind == NCIContact::HBond && c.first == target.first && c.second == target.second;
                });
                if (present != (on != 0))
                    all_ok = false;
            }
            all_ok ? kept++ : snapped++;
        }
        CurcumaLogger::result_fmt("ConfGen:   after releasing the restraint, {} of {} NCI move(s) kept "
                                  "their target bond pattern, {} snapped back",
            kept, kept + snapped, snapped);
    }

    if (template_gain_kJ > 5.0)
        CurcumaLogger::warn_fmt("ConfGen: re-optimising the input structures lowers them by up to {:.1f} "
                                "kJ/mol -- the ensemble is NOT at the minimum of '{}' (different method or "
                                "convergence). Energies below are therefore compared against re-optimised "
                                "templates, not against the file's own values.",
            template_gain_kJ, useSplitEvaluation() ? m_eval_method : m_method);

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
    // Claude Generated (Aug 2026): the four atom indices of every torsion, written next to the state
    // table. Without them the CSV cannot be used to reproduce a geometry -- a torsion is identified
    // by its central bond in the label, but the two outer atoms are a choice of the detection, and
    // guessing them wrong silently produces different dihedrals (measured: 0 of 29 reconstructed
    // correctly from the label alone).
    {
        std::ofstream def(path.substr(0, path.rfind(".torsion_states.csv")) + ".torsion_definitions.csv");
        def << "# torsion,bond,atom_i,atom_j,atom_k,atom_l  (0-based indices into the XYZ)\n";
        for (std::size_t t = 0; t < m_torsions.size(); ++t)
            def << t << "," << m_torsions[t].label(m_frames.front().molecule) << ","
                << m_torsions[t].i << "," << m_torsions[t].j << ","
                << m_torsions[t].k << "," << m_torsions[t].l << "\n";
    }
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

    loadProposalMemory();
    if (!analyseEnsemble())
        return;

    // Claude Generated (Jul 2026): the non-covalent description. Built before the torsion analysis is
    // reported so that both descriptions can be compared in one place.
    if (m_nci_analysis) {
        buildNCISpace();
        writeNCITable(outputPath(base + ".nci.csv"));
    }

    const std::vector<Transition> transitions = matchedPairs();

    writeFrameTable(outputPath(base + ".torsions.csv"));
    writeStateStatistics(outputPath(base + ".torsion_states.csv"));
    writeTransitionTable(outputPath(base + ".matched_pairs.csv"), transitions);

    reportTermVariance();
    if (m_nci_analysis)
        reportNCISpace();
    reportTransitions(transitions);

    if (m_couplings) {
        const std::vector<Coupling> couplings = doubleMutantCycles();
        writeCouplingTable(outputPath(base + ".couplings.csv"), couplings);
        reportCouplings(couplings);

        if (m_cv_folds >= 2) {
            std::vector<ModelFit> fits;
            for (int level = 0; level <= 2; ++level)
                fits.push_back(fitModel(level));
            // Levels 3 and 4 use the NCI pattern; on the same folds, so the comparison with level 1
            // is like-for-like and answers which description carries the energy.
            if (m_nci_analysis && !m_nci_pairs.empty()) {
                fits.push_back(fitModel(3));
                fits.push_back(fitModel(4));
            }
            // Claude Generated (Aug 2026): the energy terms themselves as columns. Ten continuous
            // numbers against hundreds of indicators -- and unlike the indicators they carry the
            // physical quantity rather than the pattern that produces it. Whether that is enough is
            // decided on the same folds as everything else.
            if (!m_frames.empty() && !m_frames.front().terms.empty()) {
                fits.push_back(fitModel(5));
                fits.push_back(fitModel(6));
            }
            reportModelComparison(fits);
        }
    }

    if (m_generate) {
        std::vector<Proposal> proposals = generateProposals();
        // Record what we are about to try, so a later call does not repeat it.
        {
            // Claude Generated (Aug 2026): remember BOTH descriptions. The state vector is known
            // before the build, the contact pattern only after it -- and exactly that pattern is
            // what the other move set needs in order not to propose the same thing again.
            std::vector<std::pair<std::vector<int>, std::vector<int>>> tried;
            for (const Proposal& pr : proposals) {
                std::vector<int> pattern;
                if (pr.geometry.AtomCount() > 0 && !m_nci_pairs.empty()) {
                    const std::vector<NCIContact> found = detectNCI(pr.geometry, {});
                    pattern.assign(m_nci_pairs.size(), 0);
                    std::map<std::string, int> index;
                    for (std::size_t k = 0; k < m_nci_pairs.size(); ++k)
                        index[m_nci_pairs[k].label()] = static_cast<int>(k);
                    for (const NCIContact& c : found) {
                        const auto it = index.find(c.label());
                        if (it != index.end() && it->second < static_cast<int>(pattern.size()))
                            pattern[it->second] = 1;
                    }
                }
                tried.emplace_back(pr.states, std::move(pattern));
            }
            appendProposalMemory(tried);
        }
        // Claude Generated (Aug 2026): the second move set. Appended to the same list, so both kinds
        // pass through identical build gates, optimisation, topology and novelty checks -- an NCI
        // proposal is never privileged over a torsion one.
        const int n_torsion_proposals = static_cast<int>(proposals.size());
        // Claude Generated (Aug 2026): the de-novo assembly. Not a mutation of an existing structure
        // but the composition of the individually best elements -- it reaches vectors no short walk
        // reaches. Runs before the NCI moves so the report reads in order of increasing locality.
        int n_consensus = 0;
        if (m_consensus_build) {
            const std::vector<Proposal> consensus = generateConsensusProposals();
            n_consensus = static_cast<int>(consensus.size());
            for (const Proposal& c : consensus)
                CurcumaLogger::result_fmt("ConfGen: de-novo assembly ({}) differs from its nearest ensemble "
                                          "member in {} torsion(s)", c.nci_label, c.distance);
            proposals.insert(proposals.end(), consensus.begin(), consensus.end());
        }
        // Claude Generated (Aug 2026): the crossover move set, budgeted separately like the others.
        if (m_crossover_max > 0) {
            const std::vector<Proposal> cross = generateCrossoverProposals();
            if (!cross.empty()) {
                int min_c = cross.front().distance, max_c = cross.front().distance;
                for (const Proposal& c : cross) {
                    min_c = std::min(min_c, c.distance);
                    max_c = std::max(max_c, c.distance);
                }
                CurcumaLogger::result_fmt("ConfGen: {} crossover proposal(s), each transferring a connected "
                                          "window of {} to {} torsion(s) from another conformer -- more "
                                          "torsions at once than a mutation, but every transferred block "
                                          "comes from a structure that exists",
                    static_cast<int>(cross.size()), min_c, max_c);
                proposals.insert(proposals.end(), cross.begin(), cross.end());
            } else if (m_verbosity >= 1) {
                CurcumaLogger::result_fmt("ConfGen: crossover produced nothing -- either every window "
                                          "transfer reproduces a known state vector, or fewer than two "
                                          "structures carry different states");
            }
        }
        if (m_mode_max > 0) {
            const std::vector<Proposal> modes = generateModeProposals();
            if (!modes.empty()) {
                CurcumaLogger::result_fmt("ConfGen: {} collective-mode proposal(s) -- displacements along the "
                                          "principal components of the ensemble itself, which move five to ten "
                                          "torsions together and are not restricted to observed combinations",
                    static_cast<int>(modes.size()));
                for (const Proposal& m : modes)
                    if (m_verbosity >= 2)
                        CurcumaLogger::info_fmt("ConfGen:   {} -> {} torsion(s) away from its template",
                            m.nci_label, m.distance);
                proposals.insert(proposals.end(), modes.begin(), modes.end());
            }
        }
        if (m_path_max > 0) {
            const std::vector<Proposal> path = generatePathProposals();
            if (!path.empty()) {
                CurcumaLogger::result_fmt("ConfGen: {} path image(s) -- structures BETWEEN two known "
                                          "conformers that are at least {} torsions apart, where a "
                                          "halfway point is only half that distance from either end",
                    static_cast<int>(path.size()), m_path_min_distance);
                proposals.insert(proposals.end(), path.begin(), path.end());
            }
        }
        /* Claude Generated (Aug 2026): the isomerisation move. Placed before the NCI moves
         * because it is the only one that can ADD a state to the space -- every other move set
         * recombines what the ensemble already showed. */
        if (m_isomerise_max > 0) {
            const std::vector<Proposal> iso = generateIsomerisationProposals();
            if (!iso.empty()) {
                CurcumaLogger::result_fmt("ConfGen: {} isomerisation proposal(s) -- flipping a torsion "
                                          "the whole ensemble holds in one planar state",
                    static_cast<int>(iso.size()));
                proposals.insert(proposals.end(), iso.begin(), iso.end());
            }
        }
        if (m_concerted_max > 0 && m_nci_generate && !m_nci_pairs.empty()) {
            const std::vector<Proposal> concerted = generateConcertedProposals();
            if (!concerted.empty()) {
                CurcumaLogger::result_fmt("ConfGen: {} concerted proposal(s) -- one torsion AND one "
                                          "hydrogen bond changed in the same restrained optimisation, "
                                          "the torsion chosen so that it moves one bridge partner and "
                                          "not the other",
                    static_cast<int>(concerted.size()));
                for (const Proposal& c : concerted)
                    if (m_verbosity >= 2)
                        CurcumaLogger::info_fmt("ConfGen:   {}", c.nci_label);
                proposals.insert(proposals.end(), concerted.begin(), concerted.end());
            }
        }
        if (m_nci_generate && !m_nci_pairs.empty()) {
            const std::vector<Proposal> nci = generateNCIProposals();
            CurcumaLogger::result_fmt("ConfGen: {} torsion proposal(s) + {} de-novo assembly/assemblies + {} "
                                      "NCI move(s) to build",
                n_torsion_proposals, n_consensus, static_cast<int>(nci.size()));
            proposals.insert(proposals.end(), nci.begin(), nci.end());
        }

        // Build the geometries: start from the template and set every torsion whose state differs.
        const std::vector<std::pair<int, int>> reference_bonds
            = topologyFingerprint(m_frames.front().molecule, m_topology_factor);
        int clashes = 0, restrained = 0, restrained_ok = 0, nci_built = 0, nci_failed = 0, nci_unreached = 0;
        for (Proposal& p : proposals) {
            // An NCI move has no rigid-build variant: there is no single coordinate to set. It goes
            // straight to the restrained build, which is the whole point -- the concerted relaxation
            // around a re-tied hydrogen bond is what a torsion move cannot express.
            if (p.nci_move()) {
                Molecule driven;
                if (!restrainedBuildNCI(p, driven)) {
                    nci_unreached++;
                    continue;
                }
                if (hasClash(driven, m_clash_factor)
                    || topologyFingerprint(driven, m_topology_factor) != reference_bonds) {
                    nci_failed++;
                    continue;
                }
                driven.setName(fmt::format("nci_from_{}_{}", p.template_frame + 1,
                    p.nci_targets.size() == 1 ? "d1" : "d2"));
                p.geometry = driven;
                p.restrained_build = true;
                nci_built++;
                continue;
            }
            // Claude Generated (Aug 2026): a collective-mode displacement arrives with its geometry
            // already set -- there is no torsion target to drive towards. It still has to pass the
            // clash and topology gates below, so it is not privileged, only pre-built.
            if (p.prebuilt) {
                if (p.geometry.AtomCount() == 0
                    || hasClash(p.geometry, m_clash_factor)
                    || topologyFingerprint(p.geometry, m_topology_factor) != reference_bonds) {
                    p.geometry = Molecule();
                    clashes++;
                }
                continue;
            }
            Geometry geom = m_frames[p.template_frame].molecule.getGeometry();
            for (std::size_t t = 0; t < m_torsions.size(); ++t)
                if (p.states[t] != m_frames[p.template_frame].states[t] && p.states[t] >= 0)
                    geom = TorsionSpace::setDihedral(geom, m_torsions[t], m_state_centres[t][p.states[t]]);
            /* Claude Generated (Aug 2026): the rigid build knows only STATE INDICES, so an
             * isomerisation proposal -- whose target angle deliberately has no state -- would be
             * built as an unchanged copy of its template and then pass every gate as a "new"
             * structure that is not new at all. Measured before this line existed: the proposals
             * came out at -157 deg, i.e. the template's own value, instead of the requested +24.
             * The explicit angles have to be applied here too, exactly as restrainedBuild does. */
            for (const auto& at : p.angle_targets)
                if (at.first >= 0 && at.first < static_cast<int>(m_torsions.size()))
                    geom = TorsionSpace::setDihedral(geom, m_torsions[at.first], at.second);
            Molecule mol = m_frames[p.template_frame].molecule;
            mol.setGeometry(geom);
            /* Claude Generated (Aug 2026): the rigid path needs the same honest name as the
             * restrained one -- an isomerisation has Hamming distance 0 by construction, so
             * "proposal_from_N_d0" hides the only move that adds a state to the space. */
            if (p.isomerisation())
                mol.setName(fmt::format("isomerise_from_{}_{}_{:+.0f}deg", p.template_frame + 1,
                    m_torsions[p.angle_targets.front().first].label(m_frames.front().molecule),
                    p.angle_targets.front().second));
            else
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
                // Claude Generated (Aug 2026): TWO ways to use the restraints, and which one works
                // depends on how far the proposal is from its template.
                //
                //   REPAIR (mol as start): the rigid build already produced the target FOLD -- only
                //   some atoms overlap. The restraints hold the torsions where they are while the
                //   optimiser relieves the clash. This is the only route that scales: measured on a
                //   107-atom peptide, rigidly setting all 29 dihedrals to a target vector reproduces
                //   that structure to 0.41 A, and after the free optimisation to 0.75 A.
                //
                //   DRIVE (template as start): the restraints must perform the rotation themselves.
                //   Fine for one or two torsions; for 29 at once it stalls at 74 degrees worst
                //   deviation and never arrives.
                //
                // So repair first, drive second -- the old order was drive-only, which is why the
                // de-novo assemblies could not be built at all.
                Molecule driven;
                bool ok = restrainedBuild(p, driven, &mol)
                    && !hasClash(driven, m_clash_factor)
                    && topologyFingerprint(driven, m_topology_factor) == reference_bonds;
                std::string how = "repaired";
                if (!ok) {
                    ok = restrainedBuild(p, driven)
                        && !hasClash(driven, m_clash_factor)
                        && topologyFingerprint(driven, m_topology_factor) == reference_bonds;
                    how = "driven";
                }
                if (!ok) {
                    clashes++;
                    continue;
                }
                /* Claude Generated (Aug 2026): name the isomerisation for what it is -- its Hamming
                 * distance is 0 (the state vector is the template's), so the generic name would
                 * read "proposal_from_N_d0" and hide the only move that changes the state space. */
                if (p.isomerisation())
                    driven.setName(fmt::format("isomerise_from_{}_{}_{:+.0f}deg_{}", p.template_frame + 1,
                        m_torsions[p.angle_targets.front().first].label(m_frames.front().molecule),
                        p.angle_targets.front().second, how));
                else
                    driven.setName(fmt::format("proposal_from_{}_d{}_{}", p.template_frame + 1, p.distance, how));
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
        if (nci_built + nci_failed + nci_unreached > 0)
            CurcumaLogger::result_fmt("ConfGen: NCI moves: {} built, {} did not reach their restraint "
                                      "(sterically impossible), {} rejected by the clash/topology gate",
                nci_built, nci_unreached, nci_failed);

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
