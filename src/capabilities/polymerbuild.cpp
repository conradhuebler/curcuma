#include "polymerbuild.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energycalculator.h"
#include "src/core/pseudoff.h"
#include "src/tools/geometry.h"
#include "curcumaopt.h"
#include "simplemd.h"
#include "src/tools/general.h"
#include "src/core/elements.h"

#include <regex>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <Eigen/Geometry>
#include <unsupported/Eigen/NonLinearOptimization>
#include <fmt/format.h>
#include "../../external/json.hpp"

/**
 * @brief Check if a Geometry matrix contains any NaN values
 * @param geom The geometry matrix to check (atoms x 3 coordinates)
 * @return true if any coordinate is NaN, false otherwise
 */
bool containsNaN(const Geometry& geom) {
    return geom.array().isNaN().any();
}

// ============================================================================
// Claude Generated: PolymerPlacementFunctor for LM optimization
// ============================================================================

/**
 * @brief Functor for Levenberg-Marquardt optimization of fragment placement.
 *
 * Optimizes 6 DOF (3 translation + 3 rotation) to minimize steric energy
 * between polymer and fragment atoms, with a constraint on the interface bond length.
 *
 * Residues:
 * - 0..N-1: LJ energy contribution from each polymer atom (sum over fragment atoms)
 * - N: Bond-length constraint penalty
 *
 * Based on MyFunctor from LevMarDocking.h - follows same interface pattern.
 */
struct PolymerPlacementFunctor {
    typedef double Scalar;
    enum {
        InputsAtCompileTime = Eigen::Dynamic,
        ValuesAtCompileTime = Eigen::Dynamic
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    const Molecule* m_polymer;          ///< Host (polymer, Xx→I so chain ends are visible/fixed)
    Molecule m_fragment;                 ///< Guest (fragment, Xx→cap_intermediate)
    Position m_anchor_polymer;           ///< Position of active_out in polymer
    int m_anchor_polymer_idx;            ///< Index of active_out atom in polymer
    int m_anchor_fragment_idx;           ///< Index of active_in atom in fragment
    int m_polymer_xx_idx;                ///< Index of polymer interface Xx (→I); excluded from LJ vs fragment active/xx
    int m_frag_xx_idx;                   ///< Index of fragment interface Xx (→cap); excluded from LJ vs polymer Xx
    double m_bond_length;                ///< Target bond length
    double m_k_bond;                     ///< Constraint weight
    Position m_fragment_centroid;        ///< Fragment centroid (for rotation pivot)

    inline PolymerPlacementFunctor(int inputs, int values)
        : m_inputs(inputs)
        , m_values(values)
        , m_anchor_polymer_idx(-1)
        , m_anchor_fragment_idx(-1)
        , m_polymer_xx_idx(-1)
        , m_frag_xx_idx(-1)
        , m_k_bond(100000.0)
    {
    }

    inline ~PolymerPlacementFunctor() {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

    /**
     * @brief Compute residuals for LM optimization.
     *
     * Parameters:
     *   p(0), p(1), p(2) = translation vector (Å)
     *   p(3), p(4), p(5) = rotation angles around centroid (degrees)
     *
     * Residuals:
     *   fvec(i) for i < N = sqrt(LJ energy) contribution from polymer atom i
     *   fvec(N) = sqrt(k_bond) * (distance - bond_length)
     */
    inline int operator()(const Eigen::VectorXd& p, Eigen::VectorXd& fvec) const
    {
        // Extract parameters
        Position translation{p(0), p(1), p(2)};
        Position rotation{p(3), p(4), p(5)};  // degrees

        // Transform fragment: translate to origin, rotate, translate to target
        Geometry frag_geom = m_fragment.getGeometry();
        for (int i = 0; i < frag_geom.rows(); ++i) {
            Position atom_pos = frag_geom.row(i).transpose();
            // Translate to origin (centroid)
            atom_pos -= m_fragment_centroid;
            // Rotate around origin
            atom_pos = rotatePoint(atom_pos, rotation);
            // Translate to final position
            atom_pos += translation;
            frag_geom.row(i) = atom_pos.transpose();
        }

        // Compute residuals
        int n_polymer = m_polymer->AtomCount();

        // Residues 0..N-1: LJ repulsion from each polymer atom over all fragment atoms.
        // The interface atom pair (m_anchor_polymer_idx, m_anchor_fragment_idx) is excluded
        // because those atoms will be bonded — including them would fight the bond constraint.
        for (int i = 0; i < n_polymer; ++i) {
            fvec(i) = 0;
            int elem_i = m_polymer->Atom(i).first;
            if (elem_i <= 0) continue;

            Position pos_i = m_polymer->Atom(i).second;
            double r_vdw_i = Elements::VanDerWaalsRadius[elem_i];

            for (int j = 0; j < m_fragment.AtomCount(); ++j) {
                // Skip topologically connected interface pairs to avoid fighting the bond constraint:
                //   polymer_active  ↔ frag_active   (the new bond)
                // When interface Xx indices are -1, they've been removed before LM, so no need to skip
                if (i == m_anchor_polymer_idx && j == m_anchor_fragment_idx) continue;
                if (m_polymer_xx_idx >= 0 && i == m_polymer_xx_idx && j == m_anchor_fragment_idx) continue;
                if (m_polymer_xx_idx >= 0 && m_frag_xx_idx >= 0 && i == m_polymer_xx_idx && j == m_frag_xx_idx) continue;
                if (m_frag_xx_idx >= 0 && i == m_anchor_polymer_idx && j == m_frag_xx_idx) continue;

                int elem_j = m_fragment.Atom(j).first;
                if (elem_j <= 0) continue;

                Position pos_j = frag_geom.row(j).transpose();
                double dist = (pos_i - pos_j).norm();
                if (dist < 1e-6) dist = 1e-6;

                double r_vdw = r_vdw_i + Elements::VanDerWaalsRadius[elem_j];

                if (dist < r_vdw * 1.5) {
                    double ratio = r_vdw / dist;
                    double lj = std::pow(ratio, 12) - 2.0 * std::pow(ratio, 6);
                    fvec(i) += lj;
                }
            }
            if (fvec(i) > 0) fvec(i) = std::sqrt(fvec(i));
        }

        // Residue N: Bond-length constraint
        // Position of anchor_in in transformed fragment
        Position anchor_frag_pos = frag_geom.row(m_anchor_fragment_idx).transpose();
        double dist = (anchor_frag_pos - m_anchor_polymer).norm();
        fvec(n_polymer) = std::sqrt(m_k_bond) * std::abs(dist - m_bond_length);

        return 0;
    }

    /**
     * @brief Rotate a point around all three axes.
     */
    inline Position rotatePoint(const Position& p, const Position& rotation_deg) const
    {
        Position result = p;

        // Rotation around X
        double rx = rotation_deg(0) * M_PI / 180.0;
        double y1 = result(1) * std::cos(rx) - result(2) * std::sin(rx);
        double z1 = result(1) * std::sin(rx) + result(2) * std::cos(rx);
        result(1) = y1;
        result(2) = z1;

        // Rotation around Y
        double ry = rotation_deg(1) * M_PI / 180.0;
        double x2 = result(0) * std::cos(ry) + result(2) * std::sin(ry);
        double z2 = -result(0) * std::sin(ry) + result(2) * std::cos(ry);
        result(0) = x2;
        result(2) = z2;

        // Rotation around Z
        double rz = rotation_deg(2) * M_PI / 180.0;
        double x3 = result(0) * std::cos(rz) - result(1) * std::sin(rz);
        double y3 = result(0) * std::sin(rz) + result(1) * std::cos(rz);
        result(0) = x3;
        result(1) = y3;

        return result;
    }
};

// ============================================================================
// PolymerBuild Implementation
// ============================================================================

PolymerBuild::PolymerBuild(const json& controller, bool silent)
    : m_config("polymerbuild", controller)
    , m_silent(silent)
{
    // Set up verbosity for CurcumaLogger
    int verbosity = silent ? 0 : 1;
    if (controller.contains("verbosity"))
        verbosity = std::max(0, std::min(3, controller["verbosity"].get<int>()));
    else if (m_config.get<bool>("verbose", true))
        verbosity = 2;
    CurcumaLogger::set_verbosity(verbosity);

    // Parse fragments map
    if (controller.contains("fragments") && controller["fragments"].is_object()) {
        for (auto& [name, path] : controller["fragments"].items()) {
            m_fragments[name] = path.get<std::string>();
        }
    } else {
        std::string fragments_str = m_config.get<std::string>("fragments", "{}");
        try {
            json frags = json::parse(fragments_str);
            for (auto& [name, path] : frags.items()) {
                m_fragments[name] = path.get<std::string>();
            }
        } catch (...) {
             if (fragments_str != "{}")
                CurcumaLogger::error("Failed to parse fragments JSON map: " + fragments_str);
        }
    }
}

PolymerBuild::~PolymerBuild() {}

void PolymerBuild::start()
{
    std::string sequence_str = m_config.get<std::string>("sequence", "");
    if (sequence_str.empty()) {
        CurcumaLogger::error("No polymer sequence provided. Use -sequence \"(A)10-B\"");
        return;
    }

    CurcumaLogger::info("Starting PolymerBuild");
    CurcumaLogger::param("sequence", sequence_str);

    std::vector<SequenceEntry> sequence = parseSequence(sequence_str);
    if (sequence.empty()) {
        CurcumaLogger::error("Parsed sequence is empty or invalid");
        return;
    }

    // Log parsed sequence with Xx selection info — Claude Generated
    std::string seq_str;
    for (const auto& e : sequence)
        seq_str += e.fragment_name + std::string(e.xx_selection, '\'') + " ";
    CurcumaLogger::info(fmt::format("Parsed sequence ({} monomers): {}", sequence.size(), seq_str));

    assemblePolymer(sequence);
}

void PolymerBuild::printHelp() const
{
    // Handled by ParameterRegistry
}

/// Claude Generated: Parse inner token string into SequenceEntry list
/// e.g., "AA'" → [{A,0},{A,1}], "B''" → [{B,2}], "pdmaema" → [{pdmaema,0}]
/// Convention: Each uppercase letter starts a new token. Lowercase/digits/underscores
/// continue the current token. Use '-' to force a new token boundary for names
/// starting with lowercase (rare). This allows (AA')3 to mean alternating Xx selection.
std::vector<SequenceEntry> PolymerBuild::parseTokens(const std::string& content)
{
    std::vector<SequenceEntry> tokens;
    size_t pos = 0;
    while (pos < content.size()) {
        // Skip optional separator
        if (content[pos] == '-') { ++pos; continue; }

        // Skip unexpected characters
        if (!std::isalnum(content[pos]) && content[pos] != '_' && content[pos] != '\'') {
            ++pos; continue;
        }

        // Start a new token: first char must be alphanumeric or underscore
        if (content[pos] == '\'') { ++pos; continue; }  // stray prime, skip

        size_t name_start = pos;
        ++pos;  // consume first character

        // Continue token with lowercase, digits, underscores (NOT uppercase — that starts a new token)
        while (pos < content.size() && (std::islower(content[pos]) || std::isdigit(content[pos]) || content[pos] == '_'))
            ++pos;

        std::string name = content.substr(name_start, pos - name_start);

        // Count trailing prime characters for xx_selection
        int primes = 0;
        while (pos < content.size() && content[pos] == '\'') {
            ++primes;
            ++pos;
        }

        tokens.push_back({name, primes});
    }
    return tokens;
}

/// Claude Generated: Two-level parser supporting prime notation for Xx selection
/// e.g., "(AA')3-B''" → [{A,0},{A,1},{A,0},{A,1},{A,0},{A,1},{B,2}]
std::vector<SequenceEntry> PolymerBuild::parseSequence(const std::string& sequence)
{
    std::vector<SequenceEntry> result;

    // Outer regex: matches either "(content)N" repetition groups or bare tokens with primes
    std::regex re(R"((\(([^)]+)\)(\d+))|([A-Za-z0-9_]+'*))");
    auto it = std::sregex_iterator(sequence.begin(), sequence.end(), re);
    auto end = std::sregex_iterator();

    for (; it != end; ++it) {
        std::smatch match = *it;
        if (match[1].matched) {
            // Repetition group: (content)N
            std::string inner = match[2].str();
            int count = std::stoi(match[3].str());
            std::vector<SequenceEntry> inner_tokens = parseTokens(inner);
            for (int i = 0; i < count; ++i) {
                for (const auto& tok : inner_tokens)
                    result.push_back(tok);
            }
        } else if (match[4].matched) {
            // Bare token (possibly with primes): e.g., "A", "B'", "C''"
            std::vector<SequenceEntry> tokens = parseTokens(match[4].str());
            for (const auto& tok : tokens)
                result.push_back(tok);
        }
    }

    return result;
}

int PolymerBuild::findBondedAtom(const Molecule& mol, int xx_idx) const
{
    Position xx_pos = mol.Atom(xx_idx).second;

    // Step 1: Find all atoms within Xx-bond cutoff (use 1.8 Å as max Xx-bond)
    std::vector<std::pair<double, int>> candidates;
    for (int i = 0; i < mol.AtomCount(); ++i) {
        if (i == xx_idx) continue;
        if (mol.Atom(i).first == 0) continue; // skip other Xx
        double d = (mol.Atom(i).second - xx_pos).norm();
        candidates.push_back({d, i});
    }
    std::sort(candidates.begin(), candidates.end());

    // Step 2: Walk the bond path from Xx → nearest atom → ... → first heavy atom
    for (const auto& [d, idx] : candidates) {
        if (d > 2.5) break; // too far for any bond
        int elem = mol.Atom(idx).first;

        if (elem > 1) {
            // Direct bond to heavy atom
            return idx;
        }

        if (elem == 1) {
            // H atom — find the heavy atom this H is bonded to
            Position h_pos = mol.Atom(idx).second;
            double best_hd = 1e9;
            int best_heavy = -1;
            for (int j = 0; j < mol.AtomCount(); ++j) {
                if (j == idx || j == xx_idx) continue;
                int ej = mol.Atom(j).first;
                if (ej <= 1) continue; // only heavy atoms
                double hd = (mol.Atom(j).second - h_pos).norm();
                double cutoff = (Elements::CovalentRadius[1] + Elements::CovalentRadius[ej]) * 1.15;
                if (hd < cutoff && hd < best_hd) {
                    best_hd = hd;
                    best_heavy = j;
                }
            }
            if (best_heavy >= 0) return best_heavy;
        }
    }

    // Fallback: closest heavy atom
    double min_dist = std::numeric_limits<double>::max();
    int closest = -1;
    for (int i = 0; i < mol.AtomCount(); ++i) {
        if (i == xx_idx) continue;
        if (mol.Atom(i).first <= 1) continue;
        double d = (mol.Atom(i).second - xx_pos).norm();
        if (d < min_dist) { min_dist = d; closest = i; }
    }
    return closest;
}

/**
 * @brief LM optimization for fragment placement.
 *
 * Claude Generated: Uses Levenberg-Marquardt to minimize steric energy + bond constraint.
 */
std::pair<Position, Position> PolymerBuild::optimizeFragmentPlacement(
    const Molecule& polymer,
    const Molecule& fragment,
    const Position& anchor_polymer,
    int anchor_polymer_idx,
    int anchor_fragment_idx,
    double bond_length,
    const Position& initial_translation,
    const Position& initial_rotation,
    int step_number,
    int polymer_xx_idx,
    int frag_xx_idx)
{
    int n_polymer = polymer.AtomCount();

    // Setup functor
    PolymerPlacementFunctor functor(6, n_polymer + 1);
    functor.m_polymer = &polymer;
    functor.m_fragment = fragment;
    functor.m_anchor_polymer = anchor_polymer;
    functor.m_anchor_polymer_idx = anchor_polymer_idx;
    functor.m_anchor_fragment_idx = anchor_fragment_idx;
    functor.m_polymer_xx_idx = polymer_xx_idx;
    functor.m_frag_xx_idx = frag_xx_idx;
    functor.m_bond_length = bond_length;

    // Compute fragment centroid for rotation pivot
    functor.m_fragment_centroid = fragment.Centroid();

    // Initial parameters: translation(0-2), rotation(3-5)
    Vector parameter(6);
    parameter << initial_translation(0), initial_translation(1), initial_translation(2),
                 initial_rotation(0), initial_rotation(1), initial_rotation(2);

    // LM optimization
    Eigen::NumericalDiff<PolymerPlacementFunctor> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<PolymerPlacementFunctor>> lm(numDiff);

    int max_iter = m_config.get<int>("lm_max_iter", 500);
    double tolerance = m_config.get<double>("lm_tolerance", 1e-6);

    // Always write trajectory — one file per LM call, one frame per iteration.
    // Claude Generated: helper lambda to render a single XYZ frame
    std::string traj_filename = fmt::format("{}_lm_{:03d}.xyz", polymer.Name(), step_number);
    std::ofstream traj_file(traj_filename);

    Geometry frag_geom_orig = fragment.getGeometry();

    auto writeFrame = [&](const Vector& p, const std::string& label) {
        Position trans{p(0), p(1), p(2)};
        Position rot{p(3), p(4), p(5)};
        traj_file << polymer.AtomCount() + fragment.AtomCount() << "\n"
                  << label << "\n";
        for (int i = 0; i < polymer.AtomCount(); ++i) {
            auto [elem, pos] = polymer.Atom(i);
            traj_file << Elements::ElementAbbr[elem] << " "
                      << std::fixed << std::setprecision(6)
                      << pos(0) << " " << pos(1) << " " << pos(2) << "\n";
        }
        for (int i = 0; i < fragment.AtomCount(); ++i) {
            Position fp = frag_geom_orig.row(i).transpose();
            fp -= functor.m_fragment_centroid;
            fp = functor.rotatePoint(fp, rot);
            fp += trans;
            traj_file << Elements::ElementAbbr[fragment.Atom(i).first] << " "
                      << std::fixed << std::setprecision(6)
                      << fp(0) << " " << fp(1) << " " << fp(2) << "\n";
        }
    };

    writeFrame(parameter, "LM step 0 (initial)");

    lm.minimizeInit(parameter);

    int iter = 0;
    Vector old_param = parameter;
    for (; iter < max_iter; ++iter) {
        lm.minimizeOneStep(parameter);
        writeFrame(parameter, fmt::format("LM step {}", iter + 1));
        if ((old_param - parameter).norm() < tolerance)
            break;
        old_param = parameter;
    }

    traj_file.close();
    CurcumaLogger::info(fmt::format(
        "LM step {:02d}: {} iterations, traj → {}", step_number, iter + 1, traj_filename));

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "LM optimization converged in {} iterations (norm = {:.2e})",
            iter, (old_param - parameter).norm()));
    }

    Position translation{parameter(0), parameter(1), parameter(2)};
    Position rotation{parameter(3), parameter(4), parameter(5)};

    return {translation, rotation};
}

/// Claude Generated: Thin wrapper — loads fragment from file, delegates to connectMolecule
ConnectionResult PolymerBuild::connectFragment(
    const Molecule& current_polymer,
    const std::string& fragment_file,
    const std::vector<std::pair<int, int>>& prev_tracked_xx,
    const std::vector<std::pair<int, int>>& prev_interface_bonds,
    int step_number,
    int xx_selection)
{
    Molecule next;
    next.LoadMolecule(fragment_file);
    return connectMolecule(current_polymer, next, prev_tracked_xx, prev_interface_bonds, step_number, {}, xx_selection);
}

/// Claude Generated: Core connection logic accepting a Molecule directly
ConnectionResult PolymerBuild::connectMolecule(
    const Molecule& current_polymer,
    const Molecule& next,
    const std::vector<std::pair<int, int>>& prev_tracked_xx,
    const std::vector<std::pair<int, int>>& prev_interface_bonds,
    int step_number,
    const std::vector<std::pair<int, int>>& fragment_internal_bonds,
    int xx_selection)
{
    ConnectionResult result;

    // Find Xx atoms in polymer (from prev_tracked_xx)
    if (prev_tracked_xx.empty()) {
        CurcumaLogger::error("connectMolecule: No tracked Xx atoms available for connection");
        result.polymer = current_polymer;
        return result;
    }

    // Find Xx atoms in next fragment
    std::vector<int> next_xx;
    for (int i = 0; i < next.AtomCount(); ++i) {
        if (next.Atom(i).first == 0) next_xx.push_back(i);
    }

    if (next_xx.empty()) {
        CurcumaLogger::error("connectFragment: Next fragment has no Xx connection points");
        result.polymer = current_polymer;
        return result;
    }

    // Report fragment topology
    if (CurcumaLogger::get_verbosity() >= 2) {
        std::map<int, std::vector<int>> active_to_xx;
        for (int xx : next_xx) {
            int active = findBondedAtom(next, xx);
            if (active >= 0) {
                active_to_xx[active].push_back(xx);
                CurcumaLogger::info(fmt::format(
                    "Fragment topology: Xx{} → {}{} (dist {:.3f} Å)",
                    xx, Elements::ElementAbbr[next.Atom(active).first], active,
                    (next.Atom(xx).second - next.Atom(active).second).norm()));
            }
        }
        for (const auto& [atom, xxs] : active_to_xx) {
            if (xxs.size() > 1) {
                CurcumaLogger::info(fmt::format(
                    "Fragment topology: {} Xx atoms connect to {}{} (branching backbone)",
                    xxs.size(), Elements::ElementAbbr[next.Atom(atom).first], atom));
            }
        }
    }

    // Select chain-end Xx from tracked list
    // Use the last tracked Xx for linear chain extension
    int polymer_xx_idx = prev_tracked_xx.back().first;
    int polymer_active_idx = prev_tracked_xx.back().second;

    CurcumaLogger::info(fmt::format(
        "connectFragment: polymer Xx idx={}, active={} ({})",
        polymer_xx_idx, polymer_active_idx,
        Elements::ElementAbbr[current_polymer.Atom(polymer_active_idx).first]));

    // Compute polymer chain direction for initial fragment placement
    Position polymer_xx_pos = current_polymer.Atom(polymer_xx_idx).second;
    Position polymer_active_pos = current_polymer.Atom(polymer_active_idx).second;
    Position v_out = (polymer_xx_pos - polymer_active_pos).normalized();

    // Claude Generated: Select Xx in next fragment by index (xx_selection from prime notation)
    // next_xx is built by scanning atoms in index order, so next_xx[0]=first Xx, next_xx[1]=second, etc.
    int next_xx_idx;
    if (xx_selection >= 0 && xx_selection < (int)next_xx.size()) {
        next_xx_idx = next_xx[xx_selection];
        CurcumaLogger::info(fmt::format(
            "connectMolecule: Using Xx#{} (atom idx {}) from prime notation",
            xx_selection, next_xx_idx));
    } else {
        CurcumaLogger::error(fmt::format(
            "Requested Xx#{} but fragment has only {} Xx — using Xx#0",
            xx_selection, next_xx.size()));
        next_xx_idx = next_xx[0];
    }
    int next_active_idx = findBondedAtom(next, next_xx_idx);
    if (next_active_idx < 0) {
        CurcumaLogger::error("connectFragment: Could not find bonded atom for Xx in next fragment");
        result.polymer = current_polymer;
        return result;
    }

    CurcumaLogger::info(fmt::format(
        "connectFragment: fragment Xx idx={}, active={} ({})",
        next_xx_idx, next_active_idx,
        Elements::ElementAbbr[next.Atom(next_active_idx).first]));

    // Compute target bond length
    double bond_length = (Elements::CovalentRadius[current_polymer.Atom(polymer_active_idx).first] +
                          Elements::CovalentRadius[next.Atom(next_active_idx).first])
                       * m_config.get<double>("bond_distance_scaling", 1.0);

    // Collect all Xx indices for polymer and fragment
    std::vector<int> polymer_xx_indices;
    for (int i = 0; i < current_polymer.AtomCount(); ++i)
        if (current_polymer.Atom(i).first == 0)
            polymer_xx_indices.push_back(i);

    std::vector<int> fragment_xx_indices;
    for (int i = 0; i < next.AtomCount(); ++i)
        if (next.Atom(i).first == 0)
            fragment_xx_indices.push_back(i);

    // Determine cap_intermediate element for Xx replacement (used by both polymer and fragment)
    std::string ci_lm_str = m_config.get<std::string>("cap_intermediate", "H");
    int ci_lm_elem = Elements::String2Element(ci_lm_str);
    if (ci_lm_elem <= 0) ci_lm_elem = 1;

    // Prepare polymer for LM:
    // - Interface Xx (polymer_xx_idx) is REMOVED before LM (only used for alignment)
    // - Non-interface Xx are replaced with cap_intermediate
    Molecule polymer_lm = current_polymer;
    {
        Mol poly_mol = polymer_lm.getMolInfo();
        for (int xx : polymer_xx_indices) {
            if (xx == polymer_xx_idx) continue;  // Skip interface Xx (will be removed)
            poly_mol.m_atoms[xx] = ci_lm_elem;  // Non-interface Xx → cap_intermediate
        }
        polymer_lm.LoadMolecule(poly_mol);
    }
    // Remove interface Xx from polymer
    polymer_lm = polymer_lm.AtomsRemoved({polymer_xx_idx});
    // Adjust polymer_active_idx after removal
    int polymer_active_idx_lm = polymer_active_idx;
    if (polymer_active_idx_lm > polymer_xx_idx) polymer_active_idx_lm--;

    Position anchor_polymer_pos = polymer_lm.Atom(polymer_active_idx_lm).second;

    // Prepare fragment for LM:
    // - Interface Xx (next_xx_idx) is REMOVED before LM (only used for alignment)
    // - Non-interface Xx are replaced with cap_intermediate
    // Adjust fragment_xx_idx after interface removal (computed below)
    int next_active_idx_lm = next_active_idx;

    Molecule fragment_lm = next;
    {
        Mol frag_mol = fragment_lm.getMolInfo();
        for (int xx : fragment_xx_indices) {
            if (xx == next_xx_idx) continue;  // Skip interface Xx (will be removed)
            frag_mol.m_atoms[xx] = ci_lm_elem;  // Non-interface Xx → cap_intermediate
        }
        fragment_lm.LoadMolecule(frag_mol);
    }
    // Remove interface Xx from fragment
    fragment_lm = fragment_lm.AtomsRemoved({next_xx_idx});
    // Adjust next_active_idx_lm after removal
    if (next_active_idx_lm > next_xx_idx) next_active_idx_lm--;

    // Adjust fragment_xx_indices for the new molecule (after removal)
    std::vector<int> fragment_xx_indices_lm;
    for (int xx : fragment_xx_indices) {
        if (xx == next_xx_idx) continue;  // Skip interface Xx (removed)
        fragment_xx_indices_lm.push_back(xx > next_xx_idx ? xx - 1 : xx);
    }

    // Compute initial placement — using original molecule for alignment vectors
    Position fragment_centroid = fragment_lm.Centroid();
    Position anchor_atom_orig = fragment_lm.Atom(next_active_idx_lm).second;
    Position target_anchor_pos = anchor_polymer_pos + v_out * bond_length;

    // Compute initial rotation: align v_in_fragment with -v_out
    // Use ORIGINAL molecule (next) for alignment — interface Xx is still there
    Position next_xx_pos_orig = next.Atom(next_xx_idx).second;
    Position next_active_pos_in_fragment = next.Atom(next_active_idx).second;
    Position v_in_fragment = (next_xx_pos_orig - next_active_pos_in_fragment).normalized();

    Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(v_in_fragment, -v_out);
    Eigen::Matrix3d R = q.toRotationMatrix();

    // Extract Euler angles (Rz*Ry*Rx convention used in rotatePoint)
    double ry = -std::asin(std::clamp(R(2, 0), -1.0, 1.0));
    double rx, rz;
    if (std::abs(std::cos(ry)) > 1e-6) {
        rx = std::atan2(R(2, 1), R(2, 2));
        rz = std::atan2(R(1, 0), R(0, 0));
    } else {
        rx = std::atan2(-R(1, 2), R(1, 1));
        rz = 0.0;
    }
    Position initial_rotation{rx * 180.0 / M_PI, ry * 180.0 / M_PI, rz * 180.0 / M_PI};

    // Recompute translation so rotated anchor lands exactly on target
    Position anchor_relative = anchor_atom_orig - fragment_centroid;
    Position anchor_rotated = R * anchor_relative;
    Position initial_translation = target_anchor_pos - anchor_rotated;

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "LM initial: anchor_polymer=({:.3f}, {:.3f}, {:.3f}), v_out=({:.3f}, {:.3f}, {:.3f})",
            anchor_polymer_pos(0), anchor_polymer_pos(1), anchor_polymer_pos(2),
            v_out(0), v_out(1), v_out(2)));
        CurcumaLogger::info(fmt::format(
            "LM initial: fragment_centroid=({:.3f}, {:.3f}, {:.3f}), anchor_atom_orig=({:.3f}, {:.3f}, {:.3f})",
            fragment_centroid(0), fragment_centroid(1), fragment_centroid(2),
            anchor_atom_orig(0), anchor_atom_orig(1), anchor_atom_orig(2)));
        CurcumaLogger::info(fmt::format(
            "LM initial: target_anchor=({:.3f}, {:.3f}, {:.3f}), initial_translation=({:.3f}, {:.3f}, {:.3f})",
            target_anchor_pos(0), target_anchor_pos(1), target_anchor_pos(2),
            initial_translation(0), initial_translation(1), initial_translation(2)));
        CurcumaLogger::info(fmt::format(
            "LM initial: bond_length={:.3f}, v_in=({:.3f}, {:.3f}, {:.3f})",
            bond_length, v_in_fragment(0), v_in_fragment(1), v_in_fragment(2)));
    }

    // Run LM optimization (interface Xx already removed from both molecules)
    // Pass -1 for xx indices since no interface atoms exist in LM molecules
    auto [opt_translation, opt_rotation] = optimizeFragmentPlacement(
        polymer_lm, fragment_lm,
        anchor_polymer_pos, polymer_active_idx_lm, next_active_idx_lm,
        bond_length, initial_translation, initial_rotation,
        step_number, -1, -1);  // No interface Xx in LM molecules

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "LM result: translation=({:.3f}, {:.3f}, {:.3f}), rotation=({:.3f}, {:.3f}, {:.3f})",
            opt_translation(0), opt_translation(1), opt_translation(2),
            opt_rotation(0), opt_rotation(1), opt_rotation(2)));
    }

    // Apply optimized transform to ALL fragment atoms (incl. I-substituted Xx positions)
    PolymerPlacementFunctor func(6, 1);
    Geometry frag_geom = fragment_lm.getGeometry();
    for (int i = 0; i < frag_geom.rows(); ++i) {
        Position p = frag_geom.row(i).transpose();
        p -= fragment_centroid;
        p = func.rotatePoint(p, opt_rotation);
        p += opt_translation;
        frag_geom.row(i) = p.transpose();
    }

    // RC4 fix: Dihedral scan around interface bond axis — Claude Generated
    // The LM optimizer has no torsional DOF around the new bond. Scan 12×30° rotations
    // of the fragment around the polymer_active → fragment_active axis, pick lowest clash.
    {
        Position anchor_poly = polymer_lm.Atom(polymer_active_idx_lm).second;
        Position anchor_frag = frag_geom.row(next_active_idx_lm).transpose();
        Position bond_axis = (anchor_frag - anchor_poly);
        if (bond_axis.norm() > 1e-6) {
            bond_axis.normalize();
            const int n_steps = 12;
            const double angle_step = 2.0 * M_PI / n_steps;
            double best_clash_energy = std::numeric_limits<double>::max();
            int best_step = 0;
            Geometry best_geom = frag_geom;

            for (int s = 0; s < n_steps; ++s) {
                double angle = s * angle_step;
                // Rodrigues rotation around bond_axis through anchor_frag
                Geometry rotated = frag_geom;
                double cos_a = std::cos(angle);
                double sin_a = std::sin(angle);
                for (int i = 0; i < rotated.rows(); ++i) {
                    if (i == next_active_idx_lm) continue;  // Pivot atom stays
                    Position p = rotated.row(i).transpose() - anchor_frag;
                    Position p_rot = p * cos_a + bond_axis.cross(p) * sin_a
                                   + bond_axis * (bond_axis.dot(p)) * (1.0 - cos_a);
                    rotated.row(i) = (p_rot + anchor_frag).transpose();
                }

                // Compute total LJ clash energy between polymer and rotated fragment
                double clash_energy = 0.0;
                for (int i = 0; i < polymer_lm.AtomCount(); ++i) {
                    int ei = polymer_lm.Atom(i).first;
                    if (ei <= 0) continue;
                    Position pi = polymer_lm.Atom(i).second;
                    double rvdw_i = Elements::VanDerWaalsRadius[ei];
                    for (int j = 0; j < fragment_lm.AtomCount(); ++j) {
                        if (i == polymer_active_idx_lm && j == next_active_idx_lm) continue;
                        int ej = fragment_lm.Atom(j).first;
                        if (ej <= 0) continue;
                        Position pj = rotated.row(j).transpose();
                        double dist = (pi - pj).norm();
                        if (dist < 1e-6) dist = 1e-6;
                        double rvdw = rvdw_i + Elements::VanDerWaalsRadius[ej];
                        if (dist < rvdw) {
                            double ratio = rvdw / dist;
                            clash_energy += ratio * ratio * ratio * ratio;  // (r_vdw/r)^4 repulsion
                        }
                    }
                }

                if (clash_energy < best_clash_energy) {
                    best_clash_energy = clash_energy;
                    best_step = s;
                    best_geom = rotated;
                }
            }

            if (best_step != 0) {
                frag_geom = best_geom;
                if (CurcumaLogger::get_verbosity() >= 2) {
                    CurcumaLogger::info(fmt::format(
                        "Dihedral scan: best orientation at {}° (clash energy {:.2f} → {:.2f})",
                        best_step * 30, best_clash_energy, best_clash_energy));
                }
            }
        }
    }

    // Build combined: polymer_lm (interface Xx already removed, other Xx as cap_intermediate)
    // + fragment_lm (interface Xx already removed, other Xx as cap_intermediate)
    Molecule combined = polymer_lm;
    int offset = combined.AtomCount();

    for (int i = 0; i < fragment_lm.AtomCount(); ++i)
        combined.addAtom(fragment_lm.Atom(i));

    Geometry combined_geom = combined.getGeometry();
    for (int i = 0; i < fragment_lm.AtomCount(); ++i)
        combined_geom.row(offset + i) = frag_geom.row(i);
    combined.setGeometry(combined_geom);

    // No need to remove interface Xx — already removed before LM
    // No need to convert I back to Xx — fragment Xx are already cap_intermediate

    // Check bond length after LM optimization
    Position anchor_polymer_final = combined.Atom(polymer_active_idx_lm).second;
    Position anchor_fragment_final = combined.Atom(offset + next_active_idx_lm).second;
    double actual_bond_length = (anchor_fragment_final - anchor_polymer_final).norm();
    double bond_tolerance = 0.1;  // Å - warn if bond exceeds target by more than this
    if (actual_bond_length > bond_length + bond_tolerance) {
        CurcumaLogger::warn(fmt::format(
            "LM bond length {:.3f} Å exceeds target {:.3f} Å by {:.3f} Å (tolerance: {:.3f} Å)",
            actual_bond_length, bond_length, actual_bond_length - bond_length, bond_tolerance));
    } else if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "LM bond length: {:.3f} Å (target: {:.3f} Å)",
            actual_bond_length, bond_length));
    }

    // Set up interface bond between polymer_active and fragment_active
    // Adjust indices: polymer_active_idx in original → polymer_active_idx_lm in polymer_lm
    //                 next_active_idx in original → next_active_idx_lm in fragment_lm → +offset in combined
    // Adjust previous interface bond indices for removed polymer interface Xx — Claude Generated
    result.interface_bonds.reserve(prev_interface_bonds.size() + 1);
    for (const auto& [a, b] : prev_interface_bonds) {
        int adj_a = (a > polymer_xx_idx) ? a - 1 : a;
        int adj_b = (b > polymer_xx_idx) ? b - 1 : b;
        result.interface_bonds.push_back({adj_a, adj_b});
    }
    result.interface_bonds.push_back({polymer_active_idx_lm, offset + next_active_idx_lm});

    // Carry forward fragment's internal interface bonds (for sub-chain merging) — Claude Generated
    // These are bonds WITHIN the fragment sub-chain, offset to combined indices.
    // The fragment's interface Xx was removed → adjust indices in fragment_internal_bonds.
    for (const auto& [a, b] : fragment_internal_bonds) {
        int adj_a = (a >= next_xx_idx) ? a - 1 : a;
        int adj_b = (b >= next_xx_idx) ? b - 1 : b;
        result.interface_bonds.push_back({offset + adj_a, offset + adj_b});
    }

    // Build tracked_xx list for next iteration
    // Old polymer Xx (not used for connection): adjust indices for removed interface Xx
    for (const auto& [pxx, pactive] : prev_tracked_xx) {
        if (pxx == polymer_xx_idx) continue;  // Skip interface Xx (removed)
        int adj_pxx = (pxx > polymer_xx_idx) ? pxx - 1 : pxx;  // Adjust for removed polymer interface Xx
        int adj_pactive = pactive;
        if (adj_pactive > polymer_xx_idx) adj_pactive--;
        result.tracked_xx.push_back({adj_pxx, adj_pactive});
    }

    // New fragment Xx (not used for connection): map from original indices to combined indices
    // fragment_xx_indices_lm are indices in fragment_lm (already adjusted for removed interface Xx)
    for (int i = 0; i < (int)fragment_xx_indices_lm.size(); ++i) {
        int xx_lm = fragment_xx_indices_lm[i];
        // Find active atom in fragment_lm
        int active_lm = -1;
        // Need to map from xx_lm back to original to find active, then back to lm
        for (int j = 0; j < (int)fragment_xx_indices.size(); ++j) {
            if (fragment_xx_indices[j] == next_xx_idx) continue;  // Skip interface
            int orig_xx = fragment_xx_indices[j];
            int lm_xx = (orig_xx > next_xx_idx) ? orig_xx - 1 : orig_xx;
            if (lm_xx == xx_lm) {
                // Found matching Xx, now find its active atom in original
                for (int orig_active : {findBondedAtom(next, orig_xx)}) {
                    if (orig_active < 0) continue;
                    active_lm = (orig_active > next_xx_idx) ? orig_active - 1 : orig_active;
                    break;
                }
                break;
            }
        }
        result.tracked_xx.push_back({offset + xx_lm, (active_lm >= 0) ? offset + active_lm : -1});
    }

    // Restore Xx element for tracked connection points — Claude Generated
    // Non-interface Xx were replaced with cap_intermediate during LM preparation;
    // restore them so sub-chain merging can find Xx for the next connection.
    {
        Mol combined_info = combined.getMolInfo();
        for (const auto& [xx_idx, _active] : result.tracked_xx) {
            if (xx_idx >= 0 && xx_idx < (int)combined_info.m_atoms.size())
                combined_info.m_atoms[xx_idx] = 0;  // Xx = element 0
        }
        combined.LoadMolecule(combined_info);
    }

    result.polymer = combined;
    result.fragment_offset = offset;
    result.removed_polymer_xx_idx = polymer_xx_idx;
    result.removed_fragment_xx_idx = next_xx_idx;
    return result;
}

int PolymerBuild::checkDistances(const Molecule& mol,
                                   const std::string& tag,
                                   const std::vector<std::pair<int,int>>& interface_bonds) const
{
    int problems = 0;
    int n = mol.AtomCount();

    // Prefer persistent topology (set by assemblePolymer after each rebuild);
    // fall back to the explicit interface_bonds list when topology is not yet available
    // (e.g. the first post-connect check runs before the topology rebuild block).
    bool has_topo = mol.hasPersistentTopology();
    Matrix topo;
    if (has_topo)
        topo = mol.getTopologyMatrix();

    for (int i = 0; i < n; ++i) {
        int ei = mol.Atom(i).first;
        if (ei <= 0) continue;  // skip Xx dummy atoms
        for (int j = i + 1; j < n; ++j) {
            int ej = mol.Atom(j).first;
            if (ej <= 0) continue;

            // Skip bonded pairs: topology takes priority, interface_bonds as fallback
            bool is_bonded = false;
            if (has_topo) {
                is_bonded = (topo(i, j) > 0.5);
            } else {
                for (const auto& b : interface_bonds)
                    if ((b.first == i && b.second == j) || (b.first == j && b.second == i)) {
                        is_bonded = true; break;
                    }
            }
            if (is_bonded) continue;

            double d = (mol.Atom(i).second - mol.Atom(j).second).norm();
            double r_cov = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
            double r_vdw = Elements::VanDerWaalsRadius[ei] + Elements::VanDerWaalsRadius[ej];

            if (d < r_cov * 1.1) {
                ++problems;
                if (CurcumaLogger::get_verbosity() >= 3)
                    CurcumaLogger::warn(fmt::format(
                        "{}: FUSION {}{}–{}{}: {:.3f} Å (r_cov {:.3f}, not in topology)",
                        tag, Elements::ElementAbbr[ei], i, Elements::ElementAbbr[ej], j, d, r_cov));
            } else if (d < r_vdw * 0.65) {
                ++problems;
                if (CurcumaLogger::get_verbosity() >= 3)
                    CurcumaLogger::warn(fmt::format(
                        "{}: clash {}{}–{}{}: {:.3f} Å (0.65·r_vdw {:.3f})",
                        tag, Elements::ElementAbbr[ei], i, Elements::ElementAbbr[ej], j, d, r_vdw * 0.65));
            }
        }
    }
    if (CurcumaLogger::get_verbosity() >= 2) {
        if (problems == 0)
            CurcumaLogger::info(fmt::format("{}: all distances OK", tag));
        else
            CurcumaLogger::warn(fmt::format("{}: {} problematic contact(s)", tag, problems));
    }
    return problems;
}

/// Claude Generated 2026 — Local LJ-LM refinement for post-optimization fragment reorientation
int PolymerBuild::localLJRefinement(Molecule& polymer,
                                      const std::vector<int>& atom_monomer_id,
                                      const std::vector<std::pair<int,int>>& interface_bonds,
                                      int current_monomer,
                                      int step_number)
{
    // 1. Separate polymer into "old polymer" and "new fragment"
    std::vector<int> fragment_indices;
    std::vector<int> polymer_indices;
    for (int i = 0; i < polymer.AtomCount(); ++i) {
        if (atom_monomer_id[i] == current_monomer)
            fragment_indices.push_back(i);
        else
            polymer_indices.push_back(i);
    }

    if (fragment_indices.empty() || polymer_indices.empty()) {
        CurcumaLogger::warn("localLJRefinement: empty fragment or polymer, skipping");
        return checkDistances(polymer, "localLJRefinement-empty", interface_bonds);
    }

    CurcumaLogger::info(fmt::format(
        "localLJRefinement step {:02d}: fragment has {} atoms, polymer has {} atoms",
        step_number, fragment_indices.size(), polymer_indices.size()));

    // 2. Build sub-molecules using AtomsRemoved
    Molecule polymer_sub = polymer.AtomsRemoved(fragment_indices);
    Molecule fragment_sub = polymer.AtomsRemoved(polymer_indices);

    // Replace Xx atoms with cap_intermediate for LJ calculation
    std::string ci_str = m_config.get<std::string>("cap_intermediate", "H");
    int ci_elem = Elements::String2Element(ci_str);
    if (ci_elem <= 0) ci_elem = 1;

    {
        Mol m_poly = polymer_sub.getMolInfo();
        for (int k = 0; k < polymer_sub.AtomCount(); ++k)
            if (m_poly.m_atoms[k] == 0) m_poly.m_atoms[k] = ci_elem;
        polymer_sub.LoadMolecule(m_poly);
    }
    {
        Mol m_frag = fragment_sub.getMolInfo();
        for (int k = 0; k < fragment_sub.AtomCount(); ++k)
            if (m_frag.m_atoms[k] == 0) m_frag.m_atoms[k] = ci_elem;
        fragment_sub.LoadMolecule(m_frag);
    }

    // 3. Build index maps: combined-molecule idx -> sub-molecule idx
    std::map<int, int> polymer_index_map;
    std::map<int, int> fragment_index_map;
    {
        int sub_idx = 0;
        for (int i = 0; i < polymer.AtomCount(); ++i) {
            if (atom_monomer_id[i] != current_monomer)
                polymer_index_map[i] = sub_idx++;
        }
    }
    {
        int sub_idx = 0;
        for (int i = 0; i < polymer.AtomCount(); ++i) {
            if (atom_monomer_id[i] == current_monomer)
                fragment_index_map[i] = sub_idx++;
        }
    }

    // 4. Find the single interface bond connecting fragment to polymer
    //    Only the LAST interface bond (the one just created) is used as constraint.
    //    All earlier interface bonds are within the fixed old polymer.
    int anchor_poly_combined = -1, anchor_frag_combined = -1;
    for (int bi = static_cast<int>(interface_bonds.size()) - 1; bi >= 0; --bi) {
        const auto& [a, b] = interface_bonds[bi];
        bool a_in_frag = fragment_index_map.count(a) > 0;
        bool b_in_frag = fragment_index_map.count(b) > 0;
        if ((a_in_frag && !b_in_frag) || (!a_in_frag && b_in_frag)) {
            anchor_frag_combined = a_in_frag ? a : b;
            anchor_poly_combined = a_in_frag ? b : a;
            break;
        }
    }

    if (anchor_poly_combined < 0 || anchor_frag_combined < 0) {
        CurcumaLogger::warn("localLJRefinement: no interface bond found for current monomer");
        return checkDistances(polymer, "localLJRefinement-no-bond", interface_bonds);
    }

    int anchor_poly_sub = polymer_index_map[anchor_poly_combined];
    int anchor_frag_sub = fragment_index_map[anchor_frag_combined];
    Position anchor_poly_pos = polymer_sub.Atom(anchor_poly_sub).second;

    double bond_length = (Elements::CovalentRadius[polymer.Atom(anchor_poly_combined).first] +
                            Elements::CovalentRadius[polymer.Atom(anchor_frag_combined).first]) *
                           m_config.get<double>("bond_distance_scaling", 1.0);

    // 5. Initial placement: fragment is already in place, so translation = centroid, rotation = 0
    //    The functor uses: new_pos = R * (old_pos - centroid) + translation
    //    With translation = centroid and rotation = 0, the fragment stays at its current position.
    Position frag_centroid = fragment_sub.Centroid();
    Position initial_translation = frag_centroid;
    Position initial_rotation = Position::Zero();

    // 6. Run optimizeFragmentPlacement (reuse PolymerPlacementFunctor + LM)
    auto [opt_translation, opt_rotation] = optimizeFragmentPlacement(
        polymer_sub, fragment_sub,
        anchor_poly_pos, anchor_poly_sub, anchor_frag_sub,
        bond_length, initial_translation, initial_rotation,
        step_number, -1, -1);

    CurcumaLogger::info(fmt::format(
        "localLJRefinement step {:02d}: LM result, "
        "translation=({:.3f},{:.3f},{:.3f}), rotation=({:.3f},{:.3f},{:.3f})",
        step_number,
        opt_translation(0), opt_translation(1), opt_translation(2),
        opt_rotation(0), opt_rotation(1), opt_rotation(2)));

    // 7. Apply optimized transform to fragment atoms in the combined polymer
    //    Same transform as in connectMolecule: p = R * (old - centroid) + translation
    Geometry geom = polymer.getGeometry();
    PolymerPlacementFunctor func(6, 1);  // just for rotatePoint
    for (int idx : fragment_indices) {
        Position p = geom.row(idx).transpose();
        p -= frag_centroid;
        p = func.rotatePoint(p, opt_rotation);
        p += opt_translation;
        geom.row(idx) = p.transpose();
    }
    polymer.setGeometry(geom);

    // Write debug trajectory if enabled
    if (m_config.get<bool>("write_intermediates", false)) {
        Molecule debug_mol = polymer;
        Mol debug_info = debug_mol.getMolInfo();
        for (int k = 0; k < debug_mol.AtomCount(); ++k)
            if (debug_info.m_atoms[k] == 0) debug_info.m_atoms[k] = ci_elem;
        debug_mol.LoadMolecule(debug_info);
        debug_mol.writeXYZFile(fmt::format("{}_ljlm_{:02d}.xyz", polymer.Name(), step_number));
    }

    // 8. Check distances and report
    return checkDistances(polymer,
                           fmt::format("step {:02d} post-LJ-LM", step_number),
                           interface_bonds);
}

/// Claude Generated — Validate topology against expected structure
int PolymerBuild::validateTopology(const Molecule& mol,
                                    const std::vector<int>& atom_monomer_id,
                                    const std::vector<std::pair<int,int>>& interface_bonds,
                                    const std::string& tag) const
{
    if (!mol.hasPersistentTopology()) {
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::warn(fmt::format("{}: no persistent topology to validate", tag));
        return -1;
    }

    const Matrix topo = mol.getTopologyMatrix();
    int n = mol.AtomCount();
    int spurious = 0;
    int intra_monomer = 0, cross_interface = 0;

    // Build set of allowed cross-monomer bonds (interface bonds)
    std::set<std::pair<int,int>> allowed_cross;
    for (const auto& b : interface_bonds) {
        auto canonical = b.first < b.second ? b : std::make_pair(b.second, b.first);
        allowed_cross.insert(canonical);
    }

    // Per-atom monomer lookup — Claude Generated
    auto get_monomer = [&](int atom_idx) -> int {
        if (atom_idx >= 0 && atom_idx < (int)atom_monomer_id.size())
            return atom_monomer_id[atom_idx];
        return 0;
    };

    for (int ii = 0; ii < n; ++ii) {
        for (int jj = ii + 1; jj < n; ++jj) {
            if (topo(ii, jj) < 0.5) continue;

            if (get_monomer(ii) == get_monomer(jj)) {
                ++intra_monomer;
            } else {
                // Cross-monomer bond — must be an interface bond
                auto canonical = std::make_pair(ii, jj);
                if (allowed_cross.count(canonical)) {
                    ++cross_interface;
                } else {
                    ++spurious;
                    CurcumaLogger::error(fmt::format(
                        "{}: SPURIOUS bond {}{}-{}{} (mono {}-{}, dist {:.3f} Å)",
                        tag,
                        Elements::ElementAbbr[mol.Atom(ii).first], ii,
                        Elements::ElementAbbr[mol.Atom(jj).first], jj,
                        get_monomer(ii), get_monomer(jj),
                        (mol.Atom(ii).second - mol.Atom(jj).second).norm()));
                }
            }
        }
    }

    // Detect unbound atoms: heavy atoms (or H) with zero bonds in topology — Claude Generated
    // An atom floating in space without any bond is always an error.
    int unbound = 0;
    for (int ii = 0; ii < n; ++ii) {
        int ei = mol.Atom(ii).first;
        if (ei == 0) continue;  // Xx atoms are expected to be unbound
        int bond_count = 0;
        for (int jj = 0; jj < n; ++jj)
            if (topo(ii, jj) > 0.5) ++bond_count;
        if (bond_count == 0) {
            ++unbound;
            CurcumaLogger::error(fmt::format(
                "{}: UNBOUND atom {}{}(monomer {}) at ({:.2f}, {:.2f}, {:.2f}) — no bonds in topology",
                tag, Elements::ElementAbbr[ei], ii, get_monomer(ii),
                mol.Atom(ii).second(0), mol.Atom(ii).second(1), mol.Atom(ii).second(2)));
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2 || spurious > 0 || unbound > 0) {
        CurcumaLogger::info(fmt::format(
            "{}: {} intra-monomer, {} interface, {} spurious, {} unbound",
            tag, intra_monomer, cross_interface, spurious, unbound));
    }

    if (spurious > 0)
        CurcumaLogger::error(fmt::format("{}: {} SPURIOUS cross-monomer bond(s)!", tag, spurious));
    if (unbound > 0)
        CurcumaLogger::error(fmt::format("{}: {} UNBOUND atom(s) with no bonds!", tag, unbound));

    return spurious + unbound;
}

/// Claude Generated: Find expected bond partner for unbound atoms, enforce bond, move atom
int PolymerBuild::repairUnboundAtoms(Molecule& mol,
                                      const std::vector<int>& atom_monomer_id,
                                      const std::vector<std::pair<int,int>>& interface_bonds,
                                      const std::string& tag)
{
    if (!mol.hasPersistentTopology()) return 0;

    Matrix topo = mol.getTopologyMatrix();
    Geometry geom = mol.getGeometry();
    int n = mol.AtomCount();
    int repaired = 0;

    for (int ii = 0; ii < n; ++ii) {
        int ei = mol.Atom(ii).first;
        if (ei == 0) continue;  // Skip Xx

        // Count bonds for this atom
        int bond_count = 0;
        for (int jj = 0; jj < n; ++jj)
            if (topo(ii, jj) > 0.5) ++bond_count;
        if (bond_count > 0) continue;  // Atom is properly bonded

        // Unbound atom found — search for best bond partner
        // Prefer same-monomer atoms; use smallest dist/cov_sum ratio
        int best_partner = -1;
        double best_ratio = 1e9;
        int monomer_ii = (ii < (int)atom_monomer_id.size()) ? atom_monomer_id[ii] : -1;

        for (int jj = 0; jj < n; ++jj) {
            if (jj == ii) continue;
            int ej = mol.Atom(jj).first;
            if (ej == 0) continue;  // Skip Xx

            double cov_sum = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
            double dist = (geom.row(ii) - geom.row(jj)).norm();
            if (dist < 1e-6) continue;
            double ratio = dist / cov_sum;

            // Same-monomer partners get priority (ratio halved for ranking)
            int monomer_jj = (jj < (int)atom_monomer_id.size()) ? atom_monomer_id[jj] : -1;
            double effective_ratio = (monomer_ii == monomer_jj) ? ratio : ratio + 10.0;

            if (effective_ratio < best_ratio) {
                best_ratio = effective_ratio;
                best_partner = jj;
            }
        }

        if (best_partner < 0) continue;

        int ej = mol.Atom(best_partner).first;
        double cov_sum = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
        double current_dist = (geom.row(ii) - geom.row(best_partner)).norm();

        // Set bond in topology
        topo(ii, best_partner) = topo(best_partner, ii) = 1.0;

        // Move unbound atom to correct bonding distance along current direction
        Position dir = (geom.row(ii).transpose() - geom.row(best_partner).transpose());
        if (dir.norm() < 1e-6)
            dir = Position(1.0, 0.0, 0.0);  // Fallback direction
        dir.normalize();
        Position new_pos = geom.row(best_partner).transpose() + dir * cov_sum;
        geom.row(ii) = new_pos.transpose();

        ++repaired;
        CurcumaLogger::warn(fmt::format(
            "{}: repaired UNBOUND {}{} (monomer {}) → bonded to {}{} (monomer {}), "
            "moved {:.2f} → {:.2f} Å (cov {:.2f} Å)",
            tag, Elements::ElementAbbr[ei], ii,
            (ii < (int)atom_monomer_id.size()) ? atom_monomer_id[ii] : -1,
            Elements::ElementAbbr[ej], best_partner,
            (best_partner < (int)atom_monomer_id.size()) ? atom_monomer_id[best_partner] : -1,
            current_dist, cov_sum, cov_sum));
    }

    if (repaired > 0) {
        mol.setTopologyMatrix(topo);
        mol.setGeometry(geom);
        CurcumaLogger::info(fmt::format("{}: repaired {} unbound atom(s)", tag, repaired));
    }

    return repaired;
}

/// Claude Generated: Extract fragment topology template from loaded molecule
FragmentTopologyTemplate PolymerBuild::extractFragmentTemplate(const Molecule& mol, const std::string& fragment_name) const
{
    FragmentTopologyTemplate tmpl;
    tmpl.fragment_name = fragment_name;

    // Get atoms and identify Xx atoms
    int n = mol.AtomCount();
    std::vector<int> heavy_indices;  // maps fragment-local index to original index (skipping Xx)
    std::vector<int> idx_to_local(n, -1);  // maps original index to fragment-local index (skipping Xx)

    for (int i = 0; i < n; ++i) {
        int elem = mol.Atom(i).first;
        if (elem == 0) continue;  // skip Xx
        idx_to_local[i] = static_cast<int>(heavy_indices.size());
        heavy_indices.push_back(i);
    }
    tmpl.n_atoms = static_cast<int>(heavy_indices.size());

    // Compute distance matrix for bond detection
    auto [dist, _topo] = mol.DistanceMatrix();

    // Extract bonds between heavy atoms (skip Xx-Xx and heavy-Xx bonds)
    for (int i = 0; i < n; ++i) {
        if (idx_to_local[i] < 0) continue;  // skip Xx
        int elem_i = mol.Atom(i).first;

        for (int j = i + 1; j < n; ++j) {
            if (idx_to_local[j] < 0) continue;  // skip Xx
            int elem_j = mol.Atom(j).first;

            double cov_sum = Elements::CovalentRadius[elem_i] + Elements::CovalentRadius[elem_j];
            if (dist(i, j) > 1e-6 && dist(i, j) <= cov_sum * 1.15) {
                // Bond found - store with local indices
                tmpl.internal_bonds.push_back({idx_to_local[i], idx_to_local[j]});
            }
        }
    }

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info(fmt::format(
            "Fragment template '{}': {} heavy atoms, {} internal bonds",
            fragment_name, tmpl.n_atoms, tmpl.internal_bonds.size()));
    }

    return tmpl;
}

/// Claude Generated: Rebuild topology from stored templates
/// Claude Generated: Build expected topology, diff against actual, and fix discrepancies
int PolymerBuild::verifyAndFixTopology(
    Molecule& mol,
    const std::vector<int>& atom_monomer_id,
    const std::vector<std::string>& monomer_fragment_type,
    const std::vector<std::pair<int, int>>& interface_bonds,
    const std::vector<int>& monomer_start_atoms,
    const std::string& context_label) const
{
    int n = mol.AtomCount();
    Matrix expected = Matrix::Zero(n, n);

    // Build monomer atom ranges from monomer_start_atoms
    std::vector<std::pair<int, int>> monomer_ranges;  // (start, end_exclusive) for each monomer
    for (size_t m = 0; m < monomer_start_atoms.size(); ++m) {
        int start = monomer_start_atoms[m];
        int end = (m + 1 < monomer_start_atoms.size()) ? static_cast<int>(monomer_start_atoms[m + 1]) : n;
        monomer_ranges.push_back({start, end});
    }

    // For each monomer, stamp template bonds (adjusted for atom indices)
    int monomer_count = static_cast<int>(monomer_fragment_type.size());
    bool all_templates_available = true;
    for (int m = 0; m < monomer_count; ++m) {
        const std::string& frag_name = monomer_fragment_type[m];
        auto it = m_fragment_templates.find(frag_name);
        if (it == m_fragment_templates.end()) {
            all_templates_available = false;
            CurcumaLogger::warn(fmt::format(
                "verifyAndFixTopology: no template for fragment '{}', using distance-based fallback",
                frag_name));

            // Distance-based fallback for this monomer
            int start = monomer_ranges[m].first;
            int end = monomer_ranges[m].second;
            auto [dist2, _u] = mol.DistanceMatrix();
            for (int ii = start; ii < end; ++ii) {
                int ei = mol.Atom(ii).first;
                if (ei == 0) continue;
                for (int jj = ii + 1; jj < end; ++jj) {
                    int ej = mol.Atom(jj).first;
                    if (ej == 0) continue;
                    double co = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                    if (dist2(ii, jj) > 1e-3 && dist2(ii, jj) <= co * 1.15)
                        expected(ii, jj) = expected(jj, ii) = 1.0;
                }
            }
            continue;
        }

        const FragmentTopologyTemplate& tmpl = it->second;
        int start_idx = monomer_ranges[m].first;

        for (const auto& bond : tmpl.internal_bonds) {
            int i = start_idx + bond.first;
            int j = start_idx + bond.second;
            if (i < n && j < n) {
                expected(i, j) = expected(j, i) = 1.0;
            }
        }
    }

    // Add interface bonds
    for (const auto& bond : interface_bonds) {
        if (bond.first < n && bond.second < n) {
            expected(bond.first, bond.second) = expected(bond.second, bond.first) = 1.0;
        }
    }

    // Get actual topology
    Matrix actual = mol.hasPersistentTopology() ? mol.getTopologyMatrix() : Matrix::Zero(n, n);
    if (!mol.hasPersistentTopology()) {
        auto [_, topo] = mol.DistanceMatrix();
        actual = topo;
    }

    // Compute diff: missing bonds and extra bonds
    int missing = 0;
    int extra = 0;
    std::vector<std::pair<int,int>> missing_bonds;
    std::vector<std::pair<int,int>> extra_bonds;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            bool in_expected = (expected(i, j) > 0.5);
            bool in_actual = (actual(i, j) > 0.5);
            if (in_expected && !in_actual) {
                ++missing;
                missing_bonds.push_back({i, j});
            } else if (!in_expected && in_actual) {
                ++extra;
                extra_bonds.push_back({i, j});
            }
        }
    }

    // Log diagnostics
    if (missing > 0 || extra > 0) {
        CurcumaLogger::warn(fmt::format(
            "verifyAndFixTopology{}: {} missing bonds, {} extra bonds",
            context_label.empty() ? "" : " [" + context_label + "]",
            missing, extra));

        if (CurcumaLogger::get_verbosity() >= 2) {
            for (const auto& [i, j] : missing_bonds) {
                CurcumaLogger::info(fmt::format(
                    "  MISSING: {}({}) - {}({})",
                    i, Elements::ElementAbbr[mol.Atom(i).first],
                    j, Elements::ElementAbbr[mol.Atom(j).first]));
            }
            for (const auto& [i, j] : extra_bonds) {
                CurcumaLogger::info(fmt::format(
                    "  EXTRA:  {}({}) - {}({}) [cross-monomer: {}]",
                    i, Elements::ElementAbbr[mol.Atom(i).first],
                    j, Elements::ElementAbbr[mol.Atom(j).first],
                    (i < (int)atom_monomer_id.size() && j < (int)atom_monomer_id.size() &&
                     atom_monomer_id[i] != atom_monomer_id[j]) ? "yes" : "no"));
            }
        }
    } else if (CurcumaLogger::get_verbosity() >= 2) {
        int bond_count = 0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (expected(i, j) > 0.5) ++bond_count;
        CurcumaLogger::info(fmt::format(
            "verifyAndFixTopology{}: topology correct ({} bonds, {} monomers, {} interface{})",
            context_label.empty() ? "" : " [" + context_label + "]",
            bond_count, monomer_count, interface_bonds.size(),
            all_templates_available ? ", all templates" : ", partial distance-based"));
    }

    // Set expected topology as the new persistent topology
    mol.setTopologyMatrix(expected);

    return missing + extra;
}

/// Claude Generated: Validate topology consistency across identical monomers
int PolymerBuild::validateTopologyConsistency(
    const Molecule& mol,
    const std::vector<int>& atom_monomer_id,
    const std::vector<std::string>& monomer_fragment_type,
    const std::vector<std::pair<int, int>>& interface_bonds,
    const std::string& tag) const
{
    if (!mol.hasPersistentTopology()) {
        if (CurcumaLogger::get_verbosity() >= 2)
            CurcumaLogger::warn(fmt::format("{}: no persistent topology for consistency check", tag));
        return -1;
    }

    const Matrix topo = mol.getTopologyMatrix();
    int n = mol.AtomCount();

    // Build set of interface bonds for exclusion
    std::set<std::pair<int, int>> interface_set;
    for (const auto& b : interface_bonds) {
        auto canonical = (b.first < b.second) ? b : std::make_pair(b.second, b.first);
        interface_set.insert(canonical);
    }

    // Build monomer atom lists
    std::map<int, std::vector<int>> monomer_atoms;  // monomer_id -> list of atom indices
    for (int i = 0; i < n; ++i) {
        if (i < (int)atom_monomer_id.size()) {
            monomer_atoms[atom_monomer_id[i]].push_back(i);
        }
    }

    // Count actual intra-monomer bonds for each monomer
    std::map<std::string, std::vector<int>> bonds_by_type;  // fragment_type -> list of bond counts per instance
    int inconsistencies = 0;

    for (const auto& [monomer_id, atoms] : monomer_atoms) {
        if (monomer_id < 0 || monomer_id >= (int)monomer_fragment_type.size()) continue;

        const std::string& frag_type = monomer_fragment_type[monomer_id];
        int intra_bonds = 0;

        // Count bonds within this monomer (excluding interface bonds)
        for (size_t i = 0; i < atoms.size(); ++i) {
            for (size_t j = i + 1; j < atoms.size(); ++j) {
                auto canonical = (atoms[i] < atoms[j]) ? std::make_pair(atoms[i], atoms[j])
                                                        : std::make_pair(atoms[j], atoms[i]);
                if (interface_set.count(canonical)) continue;  // exclude interface bonds
                if (topo(atoms[i], atoms[j]) > 0.5) ++intra_bonds;
            }
        }

        bonds_by_type[frag_type].push_back(intra_bonds);
    }

    // Compare against expected bond counts from templates
    for (const auto& [frag_type, bond_counts] : bonds_by_type) {
        auto it = m_fragment_templates.find(frag_type);
        if (it == m_fragment_templates.end()) continue;

        int expected = static_cast<int>(it->second.internal_bonds.size());

        for (int count : bond_counts) {
            if (count != expected) {
                ++inconsistencies;
                CurcumaLogger::warn(fmt::format(
                    "{}: TOPOLOGY INCONSISTENCY: monomer type '{}' has {} bonds, expected {} (from template)",
                    tag, frag_type, count, expected));
            }
        }
    }

    if (inconsistencies == 0 && CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::success(fmt::format("{}: all monomer topologies consistent with templates", tag));
    }

    return inconsistencies;
}

/// RC5 fix: Raised threshold from 0.85×cov to 1.15×cov, exclude bonded pairs — Claude Generated
void PolymerBuild::resolveOverlaps(Molecule& mol,
                                    const std::vector<std::pair<int,int>>& bonded_pairs,
                                    int max_steps, double max_displacement)
{
    int n = mol.AtomCount();
    int clashes_initial = 0;

    // Build bonded pair lookup for O(1) exclusion
    std::set<std::pair<int,int>> bonded_set;
    for (const auto& [a, b] : bonded_pairs) {
        bonded_set.insert({std::min(a, b), std::max(a, b)});
    }
    // Also include topology bonds if available
    if (mol.hasPersistentTopology()) {
        const Matrix& topo = mol.getTopologyMatrix();
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (topo(i, j) > 0.5)
                    bonded_set.insert({i, j});
    }

    for (int step = 0; step < max_steps; ++step) {
        Geometry coords = mol.Coords();
        Geometry grad = Geometry::Zero(n, 3);
        int clashes = 0;

        for (int i = 0; i < n; ++i) {
            int ei = mol.Atom(i).first;
            if (ei == 0) continue; // skip Xx
            for (int j = i + 1; j < n; ++j) {
                int ej = mol.Atom(j).first;
                if (ej == 0) continue;

                // Skip bonded pairs — they are supposed to be close
                if (bonded_set.count({i, j})) continue;

                Position diff = coords.row(i).transpose() - coords.row(j).transpose();
                double dist = diff.norm();
                if (dist < 1e-6) dist = 1e-6;

                double r_cov = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                double r_cut = r_cov * 1.15;  // RC5: raised from 0.85 to match bond detection cutoff
                if (dist >= r_cut) continue;

                clashes++;
                double ratio = r_cut / dist;
                double r4 = ratio * ratio * ratio * ratio;
                double gscale = 4.0 * r4 / (dist * dist);
                Position g = gscale * diff;
                grad.row(i) += g.transpose();
                grad.row(j) -= g.transpose();
            }
        }

        if (step == 0) clashes_initial = clashes;
        if (clashes == 0) {
            if (step > 0)
                CurcumaLogger::info(fmt::format("resolveOverlaps: resolved {} clashes in {} steps",
                                                  clashes_initial, step));
            break;
        }

        double gmax = 0.0;
        for (int i = 0; i < n; ++i)
            gmax = std::max(gmax, grad.row(i).norm());
        if (gmax < 1e-9) break;
        double scale = std::min(max_displacement / gmax, max_displacement);

        for (int i = 0; i < n; ++i) {
            if (mol.Atom(i).first == 0) continue;
            coords.row(i) += scale * grad.row(i);
        }
        mol.setGeometry(coords);
    }

    Geometry coords = mol.Coords();
    double min_d = std::numeric_limits<double>::max();
    for (int i = 0; i < n; ++i) {
        if (mol.Atom(i).first <= 1) continue;
        for (int j = i+1; j < n; ++j) {
            if (mol.Atom(j).first <= 1) continue;
            double d = (coords.row(i) - coords.row(j)).norm();
            if (d < min_d) min_d = d;
        }
    }
    CurcumaLogger::info(fmt::format(
        "resolveOverlaps: residual heavy-atom min-distance = {:.3f} Å", min_d));
}

/// Claude Generated: Build a sub-chain from a list of fragment names
SubchainResult PolymerBuild::buildSubchain(
    const std::vector<SequenceEntry>& fragment_entries,
    int monomer_id_offset,
    int subchain_index)
{
    SubchainResult sc_result;
    if (fragment_entries.empty()) return sc_result;

    // Load first fragment
    std::string first_name = fragment_entries[0].fragment_name;
    if (m_fragments.find(first_name) == m_fragments.end()) {
        CurcumaLogger::error("Fragment not found in map: " + first_name);
        return sc_result;
    }

    Molecule polymer;
    polymer.LoadMolecule(m_fragments[first_name]);
    polymer.setName(fmt::format("subchain_{:02d}", subchain_index));
    polymer.setCharge(0);

    // Debug: Show raw first fragment with Xx — appended to consolidated _debug.xyz
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("DEBUG STEP 1: Raw first fragment '{}' ({} atoms, with Xx)",
                                        first_name, polymer.AtomCount()));
        // Claude Generated: one debug file per step, all frames within a step have the same atom count
        std::string debug_file = fmt::format("{}_debug_01.xyz", polymer.Name());
        polymer.appendXYZFile(debug_file, fmt::format("step 01: raw_Xx ({} atoms)", polymer.AtomCount()));
        CurcumaLogger::info(fmt::format("DEBUG: Appended step1 raw_Xx frame to {}", debug_file));

        Molecule capped = polymer;
        Mol capped_info = capped.getMolInfo();
        for (int k = 0; k < capped.AtomCount(); ++k) {
            if (capped_info.m_atoms[k] == 0) capped_info.m_atoms[k] = 1;
        }
        if (capped_info.m_charge > 1000 || capped_info.m_charge < -1000) {
            capped_info.m_charge = 0;
        }
        capped.LoadMolecule(capped_info);
        capped.appendXYZFile(debug_file, fmt::format("step 01: capped_H ({} atoms)", capped.AtomCount()));
        CurcumaLogger::info(fmt::format("DEBUG: Appended step1 capped_H frame to {}", debug_file));
    }

    // Initialize Xx tracking: find all Xx atoms and their bonded atoms
    std::vector<std::pair<int, int>> tracked_xx;
    for (int i = 0; i < polymer.AtomCount(); ++i) {
        if (polymer.Atom(i).first == 0) {
            int active = findBondedAtom(polymer, i);
            tracked_xx.push_back({i, active});
        }
    }

    CurcumaLogger::info(fmt::format("First fragment: {} Xx atoms tracked", tracked_xx.size()));
    for (const auto& [xx, active] : tracked_xx) {
        CurcumaLogger::info(fmt::format("  Xx{} → active atom {} ({})",
            xx, active, Elements::ElementAbbr[polymer.Atom(active).first]));
    }

    std::vector<std::pair<int, int>> interface_bonds;
    // Per-atom monomer ID: atom_monomer_id[i] = which monomer atom i belongs to — Claude Generated
    std::vector<int> atom_monomer_id(polymer.AtomCount(), monomer_id_offset);
    // Claude Generated: Track monomer types for topology validation
    std::vector<std::string> monomer_fragment_type;
    std::vector<int> monomer_start_atoms;  // starting atom index for each monomer

    // Extract and store template for first fragment — Claude Generated
    if (m_fragment_templates.find(first_name) == m_fragment_templates.end()) {
        m_fragment_templates[first_name] = extractFragmentTemplate(polymer, first_name);
    }
    monomer_fragment_type.push_back(first_name);
    monomer_start_atoms.push_back(0);  // first monomer starts at atom 0

    int chunk_size = m_config.get<int>("chunk_size", 1);
    if (chunk_size < 1) chunk_size = 1;
    bool optimize = m_config.get<bool>("optimize", true);

    // Connect remaining fragments
    for (size_t i = 1; i < fragment_entries.size(); ++i) {
        const auto& entry = fragment_entries[i];
        std::string next_name = entry.fragment_name;
        if (m_fragments.find(next_name) == m_fragments.end()) {
            CurcumaLogger::error("Fragment not found in map: " + next_name);
            continue;
        }

        CurcumaLogger::info(fmt::format("[sc{:02d}] Adding fragment {}/{}: {}{}",
            subchain_index, i, fragment_entries.size() - 1, next_name,
            std::string(entry.xx_selection, '\'')));

        // Use LM-based connectFragment with explicit Xx selection — Claude Generated
        ConnectionResult conn_result = connectFragment(
            polymer, m_fragments[next_name], tracked_xx, interface_bonds,
            static_cast<int>(i), entry.xx_selection);

        polymer = conn_result.polymer;
        tracked_xx = conn_result.tracked_xx;
        interface_bonds = conn_result.interface_bonds;

        // Claude Generated: Extract template for new fragment type
        if (m_fragment_templates.find(next_name) == m_fragment_templates.end()) {
            Molecule frag_tmp;
            frag_tmp.LoadMolecule(m_fragments[next_name]);
            m_fragment_templates[next_name] = extractFragmentTemplate(frag_tmp, next_name);
        }

        // Claude Generated: Track monomer type and start atom
        monomer_fragment_type.push_back(next_name);
        monomer_start_atoms.push_back(conn_result.fragment_offset);

        // Update atom_monomer_id: remove entry for the deleted interface Xx, append new fragment — Claude Generated
        if (conn_result.removed_polymer_xx_idx >= 0
            && conn_result.removed_polymer_xx_idx < (int)atom_monomer_id.size())
            atom_monomer_id.erase(atom_monomer_id.begin() + conn_result.removed_polymer_xx_idx);
        int new_frag_atoms = polymer.AtomCount() - conn_result.fragment_offset;
        for (int k = 0; k < new_frag_atoms; ++k)
            atom_monomer_id.push_back(monomer_id_offset + static_cast<int>(i));

        checkDistances(polymer,
                       fmt::format("step {:02d} post-connect", i),
                       interface_bonds);

        // Debug output — appended to consolidated _debug.xyz
        if (CurcumaLogger::get_verbosity() >= 3) {
            int xx_count = 0;
            for (int k = 0; k < polymer.AtomCount(); ++k)
                if (polymer.Atom(k).first == 0) xx_count++;
            CurcumaLogger::info(fmt::format("DEBUG: After connectFragment - {} atoms, {} Xx remaining",
                                            polymer.AtomCount(), xx_count));
            // Claude Generated: one debug file per step — all frames have the same atom count
            std::string debug_file = fmt::format("{}_debug_{:02d}.xyz", polymer.Name(), i);
            polymer.appendXYZFile(debug_file, fmt::format("step {:02d}: connected ({} atoms)", i, polymer.AtomCount()));
            CurcumaLogger::info(fmt::format("DEBUG: Appended step {:02d} connected frame to {}", i, debug_file));
        }

        // Claude Generated: Verify and fix topology via expected-vs-actual diff
        verifyAndFixTopology(polymer, atom_monomer_id, monomer_fragment_type,
                              interface_bonds, monomer_start_atoms,
                              fmt::format("step {:02d} post-connect", i));

        int topo_issues = validateTopology(polymer, atom_monomer_id, interface_bonds,
                         fmt::format("step {:02d} post-connect", i));
        int repaired = repairUnboundAtoms(polymer, atom_monomer_id, interface_bonds,
                           fmt::format("step {:02d} post-connect", i));

        // Claude Generated: Validate topology consistency across identical monomers
        int consistency_issues = validateTopologyConsistency(polymer, atom_monomer_id,
            monomer_fragment_type, interface_bonds,
            fmt::format("step {:02d} consistency", i));

        // Claude Generated: per-step connectivity status visible at verbose ≥1 for easy scanning
        if (topo_issues == 0 && repaired == 0)
            CurcumaLogger::success(fmt::format("step {:02d}: connectivity OK", i));
        else if (topo_issues > 0)
            CurcumaLogger::error(fmt::format("step {:02d}: {} connectivity issue(s) — check above", i, topo_issues));
        else
            CurcumaLogger::warn(fmt::format("step {:02d}: {} atom(s) repaired", i, repaired));

        // Manual bond count verification — Claude Generated
        {
            int n = polymer.AtomCount();
            auto [dist_v, _uv] = polymer.DistanceMatrix();

            // Count bonds in topology matrix
            int topo_bonds = 0;
            const Matrix& topo_check = polymer.getTopologyMatrix();
            for (int ii = 0; ii < n; ++ii)
                for (int jj = ii + 1; jj < n; ++jj)
                    if (topo_check(ii, jj) > 0.5) ++topo_bonds;

            // Count ALL geometric bonds (ignoring monomer restriction)
            int geo_bonds_all = 0;
            for (int ii = 0; ii < n; ++ii) {
                int ei = polymer.Atom(ii).first;
                if (ei == 0) continue;
                for (int jj = ii + 1; jj < n; ++jj) {
                    int ej = polymer.Atom(jj).first;
                    if (ej == 0) continue;
                    double cutoff = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                    if (dist_v(ii, jj) > 1e-3 && dist_v(ii, jj) <= cutoff * 1.15)
                        ++geo_bonds_all;
                }
            }

            // Count per-monomer geometric bonds
            std::map<int, int> geo_per_monomer;
            for (int ii = 0; ii < n; ++ii) {
                int ei = polymer.Atom(ii).first;
                if (ei == 0) continue;
                for (int jj = ii + 1; jj < n; ++jj) {
                    int ej = polymer.Atom(jj).first;
                    if (ej == 0) continue;
                    if (atom_monomer_id[ii] != atom_monomer_id[jj]) continue;
                    double cutoff = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                    if (dist_v(ii, jj) > 1e-3 && dist_v(ii, jj) <= cutoff * 1.15)
                        ++geo_per_monomer[atom_monomer_id[ii]];
                }
            }

            // Count atoms per monomer
            std::map<int, int> atoms_per_monomer;
            for (int ii = 0; ii < n; ++ii)
                ++atoms_per_monomer[atom_monomer_id[ii]];

            // Count cross-monomer close contacts (would-be spurious) — Claude Generated
            // Use 1.3*cov cutoff (same as GFN-FF geometric detection) to catch all problematic contacts
            int cross_geo = 0;
            for (int ii = 0; ii < n; ++ii) {
                int ei = polymer.Atom(ii).first;
                if (ei == 0) continue;
                for (int jj = ii + 1; jj < n; ++jj) {
                    int ej = polymer.Atom(jj).first;
                    if (ej == 0) continue;
                    if (atom_monomer_id[ii] == atom_monomer_id[jj]) continue;
                    double co = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                    if (dist_v(ii, jj) > 1e-3 && dist_v(ii, jj) <= co * 1.3) {
                        ++cross_geo;
                        bool is_iface = false;
                        for (const auto& b : interface_bonds) {
                            if ((b.first == ii && b.second == jj) || (b.first == jj && b.second == ii)) {
                                is_iface = true; break;
                            }
                        }
                        if (!is_iface && CurcumaLogger::get_verbosity() >= 2) {
                            CurcumaLogger::warn(fmt::format(
                                "step {:02d} BOND-CHECK: cross-monomer overlap {}{}(m{})-{}{}(m{}): {:.3f} Å ({:.2f}×cov)",
                                i, Elements::ElementAbbr[ei], ii, atom_monomer_id[ii],
                                Elements::ElementAbbr[ej], jj, atom_monomer_id[jj],
                                dist_v(ii, jj), dist_v(ii, jj) / co));
                        }
                    }
                }
            }

            int spurious_cross = cross_geo - (int)interface_bonds.size();
            if (spurious_cross > 0) {
                // Pre-optimization overlap: optimizer handles this, only show at verbose ≥2
                if (CurcumaLogger::get_verbosity() >= 2)
                    CurcumaLogger::warn(fmt::format(
                        "step {:02d} pre-opt: {} cross-monomer overlaps (will optimize), topo={}, geo_all={}",
                        i, spurious_cross, topo_bonds, geo_bonds_all));
            } else {
                CurcumaLogger::info(fmt::format(
                    "step {:02d} BOND-CHECK: OK — topo={}, geo_all={}, atoms={}, monomers={}",
                    i, topo_bonds, geo_bonds_all, n, (int)atoms_per_monomer.size()));
            }

            if (CurcumaLogger::get_verbosity() >= 2) {
                for (const auto& [mono, cnt] : atoms_per_monomer) {
                    CurcumaLogger::info(fmt::format(
                        "  monomer {}: {} atoms, {} geo bonds",
                        mono, cnt, geo_per_monomer[mono]));
                }
            }
        }

        // Optimization at chunk boundaries (including last fragment before capping)
        bool should_optimize = (static_cast<int>(i) % chunk_size == 0);
        bool ljlm_enabled = m_config.get<bool>("ljlm_refine", true);
        bool ljlm_disable_md = m_config.get<bool>("ljlm_disable_md", true);

        if (optimize && should_optimize) {
            const int max_opt_retries = m_config.get<int>("overlap_retries", 3);
            std::string opt_method = m_config.get<std::string>("opt_method", "gfnff");
            std::string ci_str = m_config.get<std::string>("cap_intermediate", "H");
            int cap_intermediate_elem = Elements::String2Element(ci_str);
            if (cap_intermediate_elem <= 0) cap_intermediate_elem = 1;

            for (int opt_attempt = 0; opt_attempt < max_opt_retries; ++opt_attempt) {

                if (CurcumaLogger::get_verbosity() >= 3) {
                    CurcumaLogger::info(fmt::format("DEBUG: Before resolveOverlaps (fragment {}, attempt {}), {} atoms",
                                                    i, opt_attempt, polymer.AtomCount()));
                }

                resolveOverlaps(polymer, interface_bonds);
                checkDistances(polymer,
                               fmt::format("step {:02d} post-resolveOverlaps", i),
                               interface_bonds);

                if (CurcumaLogger::get_verbosity() >= 3) {
                    // Claude Generated: same per-step debug file — postLJ has same atom count as connected
                    std::string debug_file = fmt::format("{}_debug_{:02d}.xyz", polymer.Name(), i);
                    polymer.appendXYZFile(debug_file, fmt::format("step {:02d}: postLJ ({} atoms)", i, polymer.AtomCount()));
                    CurcumaLogger::info(fmt::format("DEBUG: Appended step {:02d} postLJ frame to {}", i, debug_file));
                }

                int opt_max_iter = m_config.get<int>("opt_max_iter", 0);
                json opt_ctrl;
                opt_ctrl["opt"]["method"] = opt_method;
                opt_ctrl["opt"]["printOutput"] = false;
                opt_ctrl["opt"]["writeXYZ"] = true;
		opt_ctrl["opt"]["silent"] = true;
                if (opt_max_iter > 0)
                    opt_ctrl["opt"]["MaxIter"] = opt_max_iter;
                if (opt_attempt == 0)
                    CurcumaLogger::info(fmt::format("Running {} optimization (after fragment {})...", opt_method, i));
                else
                    CurcumaLogger::info(fmt::format("Re-running {} optimization (fragment {}, attempt {})...", opt_method, i, opt_attempt + 1));
                if (CurcumaLogger::get_verbosity() >= 2)
                    CurcumaLogger::info(fmt::format("opt_ctrl: {}", opt_ctrl.dump()));

                // Replace Xx with cap_intermediate for optimization; enforce charge = 0 (GFN-FF EEQ requires it)
                Molecule polymer_opt = polymer;
                polymer_opt.setCharge(0);
                Mol m_opt = polymer_opt.getMolInfo();
                for (int k = 0; k < polymer_opt.AtomCount(); ++k) {
                    if (m_opt.m_atoms[k] == 0)
                        m_opt.m_atoms[k] = cap_intermediate_elem;
                }
                // Pass complete bond list from topology to GFN-FF — Claude Generated
                if (polymer.hasPersistentTopology()) {
                    const Matrix topo = polymer.getTopologyMatrix();
                    m_opt.m_bonds.clear();
                    int n = topo.rows();
                    int topo_bond_count = 0;
                    for (int ii = 0; ii < n; ++ii)
                        for (int jj = ii + 1; jj < n; ++jj)
                            if (topo(ii, jj) > 0.5) {
                                m_opt.m_bonds.push_back({ii, jj});
                                ++topo_bond_count;
                            }

                    // Claude Generated: Debug - verify bond list matches topology
                    if (CurcumaLogger::get_verbosity() >= 2) {
                        // Count bonds per atom in topology
                        std::vector<int> topo_bonds_per_atom(n, 0);
                        for (int ii = 0; ii < n; ++ii)
                            for (int jj = 0; jj < n; ++jj)
                                if (ii != jj && topo(ii, jj) > 0.5)
                                    ++topo_bonds_per_atom[ii];

                        // Count bonds per atom in passed list
                        std::vector<int> passed_bonds_per_atom(n, 0);
                        for (const auto& b : m_opt.m_bonds) {
                            ++passed_bonds_per_atom[b.first];
                            ++passed_bonds_per_atom[b.second];
                        }

                        // Compare and report mismatches
                        int mismatches = 0;
                        for (int ii = 0; ii < n; ++ii) {
                            if (topo_bonds_per_atom[ii] != passed_bonds_per_atom[ii]) {
                                ++mismatches;
                                if (mismatches <= 5) {
                                    CurcumaLogger::warn(fmt::format(
                                        "TOPOLOGY MISMATCH: Atom {} has {} bonds in topology but {} in passed list",
                                        ii, topo_bonds_per_atom[ii], passed_bonds_per_atom[ii]));
                                }
                            }
                        }

                        if (mismatches > 0)
                            CurcumaLogger::error(fmt::format(
                                "TOPOLOGY MISMATCH: {} atoms have bond count mismatch!", mismatches));
                        else
                            CurcumaLogger::info(fmt::format(
                                "Passing {} bonds from topology to GFN-FF (verified, no mismatches)", topo_bond_count));
                    }
                }
                polymer_opt.LoadMolecule(m_opt);
                polymer_opt.setCharge(0);

                if (CurcumaLogger::get_verbosity() >= 3) {
                    // Claude Generated: same per-step debug file — forFF (Xx→H) has same atom count
                    std::string debug_file = fmt::format("{}_debug_{:02d}.xyz", polymer.Name(), i);
                    polymer_opt.appendXYZFile(debug_file, fmt::format("step {:02d}: forFF ({} atoms)", i, polymer_opt.AtomCount()));
                    CurcumaLogger::info(fmt::format("DEBUG: Appended step {:02d} forFF frame to {}", i, debug_file));
                }

                // Save current good structure before optimization
                Molecule polymer_before_opt = polymer;
                Geometry geom_before_opt = polymer.getGeometry();

                // CurcumaOpt(silent=true) resets global CurcumaLogger verbosity → save/restore
                int saved_verbosity = CurcumaLogger::get_verbosity();
                CurcumaOpt opt(opt_ctrl, true);
                opt.overrideBasename(fmt::format("{}_opt_{:02d}", polymer.Name(), i));
                opt.addMolecule(polymer_opt);
                opt.start();
                CurcumaLogger::set_verbosity(saved_verbosity);
                const std::vector<Molecule>* finals = opt.Molecules();
                if (finals && !finals->empty()) {
                    // Check if result contains NaN before using it
                    Geometry result_geom = finals->back().Coords();
                    if (!containsNaN(result_geom)) {
                        polymer.setGeometry(result_geom);

                        if (CurcumaLogger::get_verbosity() >= 3) {
                            // Claude Generated: same per-step debug file — postFF has same atom count
                            std::string debug_file = fmt::format("{}_debug_{:02d}.xyz", polymer.Name(), i);
                            polymer.appendXYZFile(debug_file, fmt::format("step {:02d}: postFF ({} atoms)", i, polymer.AtomCount()));
                            CurcumaLogger::info(fmt::format("DEBUG: Appended step {:02d} postFF frame to {}", i, debug_file));
                        }
                    } else {
                        // Optimization produced invalid structure - restore and abort fragment
                        CurcumaLogger::warn("Optimization returned invalid structure (NaN detected)");
                        CurcumaLogger::info("Using last good structure and terminating polymer assembly early");
                        polymer = polymer_before_opt;  // Restore good structure
                        break;  // Exit fragment assembly loop
                    }
                }

                // Verify and fix topology after optimization — Claude Generated
                verifyAndFixTopology(polymer, atom_monomer_id, monomer_fragment_type,
                                      interface_bonds, monomer_start_atoms,
                                      fmt::format("step {:02d} post-FF (attempt {})", i, opt_attempt + 1));

                validateTopology(polymer, atom_monomer_id, interface_bonds,
                                 fmt::format("step {:02d} post-FF (attempt {})", i, opt_attempt + 1));
                repairUnboundAtoms(polymer, atom_monomer_id, interface_bonds,
                                   fmt::format("step {:02d} post-FF (attempt {})", i, opt_attempt + 1));

                checkDistances(polymer,
                               fmt::format("step {:02d} post-FF", i),
                               interface_bonds);

                // Net-repulsion displacement: quick first pass to resolve small overlaps — Claude Generated 2026
                {
                    const int max_disp_rounds = 10;
                    int overlaps_found = 0;
                    for (int disp_round = 0; disp_round < max_disp_rounds; ++disp_round) {
                        int n_ov = polymer.AtomCount();
                        auto [dist_fresh, _uf] = polymer.DistanceMatrix();
                        Geometry geom = polymer.getGeometry();
                        const Matrix& topo_post = polymer.getTopologyMatrix();
                        std::vector<Position> displacement(n_ov, Position::Zero());
                        overlaps_found = 0;
                        for (int ii = 0; ii < n_ov; ++ii) {
                            int ei = polymer.Atom(ii).first;
                            if (ei == 0) continue;
                            for (int jj = ii + 1; jj < n_ov; ++jj) {
                                int ej = polymer.Atom(jj).first;
                                if (ej == 0) continue;
                                if (atom_monomer_id[ii] == atom_monomer_id[jj]) continue;
                                if (topo_post(ii, jj) > 0.5) continue;
                                bool is_iface = false;
                                for (const auto& b : interface_bonds)
                                    if ((b.first == ii && b.second == jj) || (b.first == jj && b.second == ii))
                                        { is_iface = true; break; }
                                if (is_iface) continue;
                                double d = dist_fresh(ii, jj);
                                double co = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                                if (d < 1e-6 || d >= co * 1.3) continue;
                                ++overlaps_found;
                                double overlap_ratio = (co * 1.15 - d) / (co * 1.15);
                                Position dir = (geom.row(ii).transpose() - geom.row(jj).transpose()).normalized();
                                displacement[ii] += dir * overlap_ratio * 0.3;
                                displacement[jj] -= dir * overlap_ratio * 0.3;
                            }
                        }
                        if (overlaps_found == 0) break;
                        for (int k = 0; k < n_ov; ++k)
                            if (displacement[k].norm() > 1e-6)
                                geom.row(k) += displacement[k].transpose();
                        polymer.setGeometry(geom);
                    }
                }

                // LJ-LM iterative refinement: reorient fragment to minimize steric clashes — Claude Generated 2026
                if (ljlm_enabled) {
                    // Iterative LJ-LM → GeoOpt until self-consistency or max iterations
                    int ljlm_max = m_config.get<int>("ljlm_max_iter", 3);
                    int current_monomer_id = monomer_id_offset + static_cast<int>(i);

                    for (int ljlm_iter = 0; ljlm_iter < ljlm_max; ++ljlm_iter) {
                        int overlaps = localLJRefinement(polymer, atom_monomer_id,
                                                          interface_bonds, current_monomer_id,
                                                          static_cast<int>(i));

                        verifyAndFixTopology(polymer, atom_monomer_id, monomer_fragment_type,
                                              interface_bonds, monomer_start_atoms,
                                              fmt::format("step {:02d} LJ-LM iter {}", i, ljlm_iter));

                        if (overlaps == 0) {
                            CurcumaLogger::success(fmt::format(
                                "step {:02d}: LJ-LM refinement converged at iteration {}", i, ljlm_iter));
                            break;
                        }

                        // Re-optimize with GFN-FF (only if not converged and not last iteration)
                        if (ljlm_iter < ljlm_max - 1) {
                            CurcumaLogger::info(fmt::format(
                                "step {:02d}: {} overlaps after LJ-LM iter {}, re-optimizing...",
                                i, overlaps, ljlm_iter));

                            Molecule polymer_opt = polymer;
                            polymer_opt.setCharge(0);
                            Mol m_opt = polymer_opt.getMolInfo();
                            for (int k = 0; k < polymer_opt.AtomCount(); ++k) {
                                if (m_opt.m_atoms[k] == 0)
                                    m_opt.m_atoms[k] = cap_intermediate_elem;
                            }
                            if (polymer.hasPersistentTopology()) {
                                const Matrix& topo = polymer.getTopologyMatrix();
                                m_opt.m_bonds.clear();
                                int n = topo.rows();
                                for (int ii2 = 0; ii2 < n; ++ii2)
                                    for (int jj2 = ii2 + 1; jj2 < n; ++jj2)
                                        if (topo(ii2, jj2) > 0.5)
                                            m_opt.m_bonds.push_back({ii2, jj2});
                            }
                            m_opt.m_charge = 0;
                            polymer_opt.LoadMolecule(m_opt);

                            json opt_ctrl;
                            opt_ctrl["opt"]["method"] = opt_method;
                            opt_ctrl["opt"]["printOutput"] = false;
                            opt_ctrl["opt"]["silent"] = true;
                            opt_ctrl["opt"]["writeXYZ"] = m_config.get<bool>("write_intermediates", false);

                            CurcumaOpt optimizer(opt_ctrl, true);
                            optimizer.overrideBasename(fmt::format("{}_ljlm_opt_{:02d}_{}", polymer.Name(), i, ljlm_iter));
                            optimizer.addMolecule(polymer_opt);
                            optimizer.start();

                            const std::vector<Molecule>* ljlm_finals = optimizer.Molecules();
                            Geometry result_geom;
                            if (ljlm_finals && !ljlm_finals->empty())
                                result_geom = ljlm_finals->back().Coords();
                            if (!containsNaN(result_geom)) {
                                polymer.setGeometry(result_geom);
                            }

                            verifyAndFixTopology(polymer, atom_monomer_id, monomer_fragment_type,
                                                   interface_bonds, monomer_start_atoms,
                                                   fmt::format("step {:02d} post-LJ-LM-opt iter {}", i, ljlm_iter));
                        }
                    }
                } else {
                    // Fallback: net-repulsion displacement (original behavior)
                    const int max_disp_rounds = 10;
                    int overlaps_found = 0;

                    for (int disp_round = 0; disp_round < max_disp_rounds; ++disp_round) {
                        int n_ov = polymer.AtomCount();
                        auto [dist_fresh, _uf] = polymer.DistanceMatrix();
                        Geometry geom = polymer.getGeometry();
                        const Matrix& topo_post = polymer.getTopologyMatrix();

                        std::vector<Position> displacement(n_ov, Position::Zero());
                        overlaps_found = 0;

                        for (int ii = 0; ii < n_ov; ++ii) {
                            int ei = polymer.Atom(ii).first;
                            if (ei == 0) continue;
                            for (int jj = ii + 1; jj < n_ov; ++jj) {
                                int ej = polymer.Atom(jj).first;
                                if (ej == 0) continue;
                                if (atom_monomer_id[ii] == atom_monomer_id[jj]) continue;
                                if (topo_post(ii, jj) > 0.5) continue;
                                bool is_iface = false;
                                for (const auto& b : interface_bonds)
                                    if ((b.first == ii && b.second == jj) || (b.first == jj && b.second == ii))
                                        { is_iface = true; break; }
                                if (is_iface) continue;

                                double d = dist_fresh(ii, jj);
                                double co = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                                if (d < 1e-6 || d >= co * 1.3) continue;

                                ++overlaps_found;
                                double overlap_ratio = (co * 1.15 - d) / (co * 1.15);
                                Position dir = (geom.row(ii).transpose() - geom.row(jj).transpose()).normalized();
                                displacement[ii] += dir * overlap_ratio * 0.3;
                                displacement[jj] -= dir * overlap_ratio * 0.3;
                            }
                        }

                        if (overlaps_found == 0) {
                            if (opt_attempt > 0 || disp_round > 0)
                                CurcumaLogger::success(fmt::format("step {:02d}: overlaps resolved (attempt {}, round {})",
                                                                    i, opt_attempt + 1, disp_round + 1));
                            break;
                        }

                        for (int k = 0; k < n_ov; ++k)
                            if (displacement[k].norm() > 1e-6)
                                geom.row(k) += displacement[k].transpose();

                        polymer.setGeometry(geom);
                        if (disp_round < max_disp_rounds - 1)
                            CurcumaLogger::info(fmt::format(
                                "step {:02d}: {} overlaps, net-repulsion round {}/{}",
                                i, overlaps_found, disp_round + 1, max_disp_rounds));
                    }

                    if (overlaps_found > 0 && opt_attempt == max_opt_retries - 1 && !ljlm_disable_md) {
                        // Cold MD fallback only if LJ-LM is disabled
                        CurcumaLogger::warn(fmt::format(
                            "step {:02d}: {} overlaps remain after all retries, running cold MD...",
                            i, overlaps_found));
                        Molecule md_mol = polymer;
                        Mol md_info = md_mol.getMolInfo();
                        for (int k = 0; k < md_mol.AtomCount(); ++k)
                            if (md_info.m_atoms[k] == 0)
                                md_info.m_atoms[k] = cap_intermediate_elem;
                        if (polymer.hasPersistentTopology()) {
                            const Matrix& t = polymer.getTopologyMatrix();
                            md_info.m_bonds.clear();
                            int n_md = polymer.AtomCount();
                            for (int ii2 = 0; ii2 < n_md; ++ii2)
                                for (int jj2 = ii2 + 1; jj2 < n_md; ++jj2)
                                    if (t(ii2, jj2) > 0.5)
                                        md_info.m_bonds.push_back({ii2, jj2});
                        }
                        md_info.m_charge = 0;
                        md_mol.LoadMolecule(md_info);
                        md_mol.setCharge(0);

                        double cold_T  = m_config.get<double>("cold_md_temperature", 10.0);
                        double cold_t  = m_config.get<double>("cold_md_time", 200.0);
                        double cold_dt = m_config.get<double>("cold_md_time_step", 0.5);

                        json cold_md_ctrl;
                        cold_md_ctrl["simplemd"]["method"] = opt_method;
                        cold_md_ctrl["simplemd"]["max_time"] = cold_t;
                        cold_md_ctrl["simplemd"]["temperature"] = cold_T;
                        cold_md_ctrl["simplemd"]["time_step"] = cold_dt;
                        cold_md_ctrl["simplemd"]["thermostat"] = "csvr";
                        cold_md_ctrl["simplemd"]["coupling"] = 1.0;
                        cold_md_ctrl["simplemd"]["verbosity"] = 0;
                        cold_md_ctrl["simplemd"]["write_xyz"] = false;
                        cold_md_ctrl["simplemd"]["no_center"] = true;
                        cold_md_ctrl["simplemd"]["dump"] = 1;
                        cold_md_ctrl["simplemd"]["restart"] = false;
                        cold_md_ctrl["simplemd"]["print_frequency"] = 1000;

                        int sv = CurcumaLogger::get_verbosity();
                        Molecule polymer_before_md = polymer;

                        SimpleMD cold_md(cold_md_ctrl, true);
                        cold_md.setMolecule(md_mol);
                        cold_md.Initialise();
                        cold_md.start();
                        CurcumaLogger::set_verbosity(sv);

                        bool md_was_stable = cold_md.wasStable();
                        Geometry result_geom = cold_md.CurrentMolecule().Coords();
                        bool result_has_nan = containsNaN(result_geom);

                        if (!md_was_stable || result_has_nan) {
                            CurcumaLogger::warn(fmt::format(
                                "step {:02d}: Cold MD {} - restoring last good structure",
                                i, !md_was_stable ? "unstable" : "NaN"));
                            polymer = polymer_before_md;
                        } else {
                            polymer.setGeometry(result_geom);
                        }
                    }
                }
            } // end retry loop
        }

        if (m_config.get<bool>("dynamics", false) && should_optimize && !(ljlm_enabled && ljlm_disable_md)) {
            CurcumaLogger::info("Running intermediate dynamics...");
            json md_ctrl;
            md_ctrl["simplemd"]["max_time"] = m_config.get<int>("md_steps", 1000);
            md_ctrl["simplemd"]["method"] = m_config.get<std::string>("md_method", "gfnff");
            md_ctrl["simplemd"]["temperature"] = m_config.get<double>("md_temperature", 300.0);
            md_ctrl["simplemd"]["time_step"] = m_config.get<double>("md_time_step", 1.0);
            md_ctrl["simplemd"]["thermostat"] = "csvr";
            md_ctrl["simplemd"]["coupling"] = 10.0;
            md_ctrl["simplemd"]["verbosity"] = 0;
            md_ctrl["simplemd"]["write_xyz"] = false;
            md_ctrl["simplemd"]["restart"] = false;
            md_ctrl["simplemd"]["no_center"] = true;

            // Prepare molecule for MD: Xx → cap_intermediate, topology as bonds — Claude Generated
            std::string ci_md_str = m_config.get<std::string>("cap_intermediate", "H");
            int ci_md_elem = Elements::String2Element(ci_md_str);
            if (ci_md_elem <= 0) ci_md_elem = 1;
            Molecule md_mol = polymer;
            Mol md_info = md_mol.getMolInfo();
            for (int k = 0; k < md_mol.AtomCount(); ++k)
                if (md_info.m_atoms[k] == 0)
                    md_info.m_atoms[k] = ci_md_elem;
            if (polymer.hasPersistentTopology()) {
                const Matrix& t = polymer.getTopologyMatrix();
                md_info.m_bonds.clear();
                int nn = t.rows();
                for (int ii2 = 0; ii2 < nn; ++ii2)
                    for (int jj2 = ii2 + 1; jj2 < nn; ++jj2)
                        if (t(ii2, jj2) > 0.5)
                            md_info.m_bonds.push_back({ii2, jj2});
            }
            md_info.m_charge = 0;
            md_mol.LoadMolecule(md_info);
            md_mol.setCharge(0);

            // Save current good structure before regular MD
            Molecule polymer_before_md = polymer;

            int sv = CurcumaLogger::get_verbosity();
            SimpleMD md(md_ctrl, true);
            md.setMolecule(md_mol);
            md.Initialise();
            md.start();
            CurcumaLogger::set_verbosity(sv);

            // Check if MD completed without instability — Claude Generated
            bool md_was_stable = md.wasStable();
            Geometry md_result_geom = md.CurrentMolecule().Coords();
            bool md_has_nan = containsNaN(md_result_geom);

            if (!md_was_stable || md_has_nan) {
                CurcumaLogger::warn(fmt::format(
                    "step {:02d}: Regular MD {} - using last good structure",
                    i, !md_was_stable ? "became unstable (NaN/Inf velocities)" : "returned NaN coordinates"));
                polymer = polymer_before_md;
                break;
            }

            // MD was stable — apply geometry back to polymer (keeps Xx elements)
            polymer.setGeometry(md_result_geom);
        }

        if (m_config.get<bool>("write_intermediates", false)) {
            Molecule debug_mol = polymer;
            Mol debug_info = debug_mol.getMolInfo();
            for (int k = 0; k < debug_mol.AtomCount(); ++k) {
                if (debug_info.m_atoms[k] == 0) {
                    debug_info.m_atoms[k] = 1;
                }
            }
            debug_mol.LoadMolecule(debug_info);
            std::string step_file = fmt::format("{}_step_{:03d}.xyz", polymer.Name(), i);
            debug_mol.writeXYZFile(step_file);
            CurcumaLogger::info(fmt::format("Wrote intermediate (Xx→H): {}", step_file));
        }
    }

    // Return sub-chain result for use by assemblePolymer — Claude Generated
    sc_result.molecule = polymer;
    sc_result.tracked_xx = tracked_xx;
    sc_result.interface_bonds = interface_bonds;
    sc_result.atom_monomer_id = atom_monomer_id;
    sc_result.monomer_fragment_type = monomer_fragment_type;
    sc_result.monomer_start_atoms = monomer_start_atoms;
    return sc_result;
}

/// Claude Generated: Assemble polymer — dispatches between sequential and divide-and-conquer
void PolymerBuild::assemblePolymer(const std::vector<SequenceEntry>& sequence)
{
    int scs = m_config.get<int>("subchain_size", 0);

    if (scs >= 2 && scs < (int)sequence.size()) {
        // ====== Divide-and-conquer assembly ======
        CurcumaLogger::info(fmt::format(
            "Divide-and-conquer assembly: {} monomers, sub-chain size {}",
            sequence.size(), scs));

        // Split sequence into sub-sequences
        std::vector<std::vector<SequenceEntry>> sub_sequences;
        for (size_t i = 0; i < sequence.size(); i += scs) {
            size_t end = std::min(i + (size_t)scs, sequence.size());
            sub_sequences.push_back(
                std::vector<SequenceEntry>(sequence.begin() + i, sequence.begin() + end));
        }
        CurcumaLogger::info(fmt::format("Split into {} sub-chains", sub_sequences.size()));

        // Phase 1: Build unique sub-chains, clone duplicates — Claude Generated
        // Identical fragment sequences produce identical geometry, so we build once and
        // clone with remapped monomer IDs.  For a homopolymer (A)30 with scs=5 this
        // reduces 6 independent builds to 1 build + 5 clones.
        std::vector<SubchainResult> subchains;
        std::map<std::string, std::pair<SubchainResult, int>> subchain_cache;  // key → (result, original_offset)

        for (size_t sc = 0; sc < sub_sequences.size(); ++sc) {
            int monomer_offset = static_cast<int>(sc * scs);

            // Build a cache key from the fragment name list including Xx selection — Claude Generated
            std::string cache_key;
            for (const auto& entry : sub_sequences[sc])
                cache_key += entry.fragment_name + std::string(entry.xx_selection, '\'') + ",";

            auto it = subchain_cache.find(cache_key);
            if (it != subchain_cache.end()) {
                // Clone cached sub-chain with remapped monomer IDs
                SubchainResult clone = it->second.first;
                int id_delta = monomer_offset - it->second.second;
                for (auto& id : clone.atom_monomer_id)
                    id += id_delta;
                subchains.push_back(std::move(clone));
                CurcumaLogger::info(fmt::format(
                    "=== Sub-chain {}/{}: cloned from cache (identical sequence) ===",
                    sc + 1, sub_sequences.size()));
            } else {
                // Build new sub-chain and cache it
                CurcumaLogger::info(fmt::format(
                    "=== Building sub-chain {}/{} ({} monomers) ===",
                    sc + 1, sub_sequences.size(), sub_sequences[sc].size()));
                SubchainResult result = buildSubchain(sub_sequences[sc], monomer_offset, static_cast<int>(sc));
                subchain_cache[cache_key] = {result, monomer_offset};
                subchains.push_back(std::move(result));
            }
        }

        // Phase 2: Connect sub-chains sequentially
        CurcumaLogger::info("=== Connecting sub-chains ===");
        Molecule polymer = subchains[0].molecule;
        std::vector<std::pair<int, int>> tracked_xx = subchains[0].tracked_xx;
        std::vector<std::pair<int, int>> interface_bonds = subchains[0].interface_bonds;
        std::vector<int> atom_monomer_id = subchains[0].atom_monomer_id;
        // Claude Generated: Track monomer types for topology validation
        std::vector<std::string> monomer_fragment_type = subchains[0].monomer_fragment_type;
        std::vector<int> monomer_start_atoms = subchains[0].monomer_start_atoms;

        std::string opt_method = m_config.get<std::string>("opt_method", "gfnff");
        std::string ci_str = m_config.get<std::string>("cap_intermediate", "H");
        int cap_intermediate_elem = Elements::String2Element(ci_str);
        if (cap_intermediate_elem <= 0) cap_intermediate_elem = 1;

        for (size_t sc = 1; sc < subchains.size(); ++sc) {
            CurcumaLogger::info(fmt::format(
                "Connecting sub-chain {}/{} ({} atoms) to polymer ({} atoms)...",
                sc + 1, subchains.size(),
                subchains[sc].molecule.AtomCount(), polymer.AtomCount()));

            // Connect using connectMolecule, passing sub-chain's internal interface bonds
            ConnectionResult conn = connectMolecule(
                polymer,
                subchains[sc].molecule,
                tracked_xx,
                interface_bonds,
                static_cast<int>(1000 + sc),  // step_number offset to distinguish from sub-chain steps
                subchains[sc].interface_bonds);

            polymer = conn.polymer;
            tracked_xx = conn.tracked_xx;
            interface_bonds = conn.interface_bonds;

            // Update atom_monomer_id: remove deleted Xx, append sub-chain B's IDs
            if (conn.removed_polymer_xx_idx >= 0
                && conn.removed_polymer_xx_idx < (int)atom_monomer_id.size())
                atom_monomer_id.erase(atom_monomer_id.begin() + conn.removed_polymer_xx_idx);

            // Append sub-chain B's atom_monomer_id (already globally unique from monomer_offset)
            // One Xx was removed from sub-chain B during connection — skip that entry
            const auto& sc_ids = subchains[sc].atom_monomer_id;
            int removed_xx = conn.removed_fragment_xx_idx;
            for (int k = 0; k < (int)sc_ids.size(); ++k) {
                if (k == removed_xx) continue;  // This Xx was removed during connection
                atom_monomer_id.push_back(sc_ids[k]);
            }

            // Claude Generated: Merge monomer type tracking
            for (const auto& ft : subchains[sc].monomer_fragment_type)
                monomer_fragment_type.push_back(ft);
            // RC1 fix: Adjust sub-chain B's monomer_start_atoms by fragment_offset
            // and account for the removed interface Xx atom — Claude Generated
            for (int start : subchains[sc].monomer_start_atoms) {
                int adjusted = conn.fragment_offset + start;
                // If the removed fragment Xx was before this monomer's start, shift down by 1
                if (conn.removed_fragment_xx_idx >= 0 && start > conn.removed_fragment_xx_idx)
                    adjusted--;
                monomer_start_atoms.push_back(adjusted);
            }

            // Sanity check: atom_monomer_id size must match polymer atom count — Claude Generated
            if ((int)atom_monomer_id.size() != polymer.AtomCount()) {
                CurcumaLogger::error(fmt::format(
                    "subchain-join {:02d}: atom_monomer_id size mismatch: {} vs {} atoms — topology may be incorrect",
                    sc, atom_monomer_id.size(), polymer.AtomCount()));
            }

            // Verify and fix topology after sub-chain connection — Claude Generated
            verifyAndFixTopology(polymer, atom_monomer_id, monomer_fragment_type,
                                  interface_bonds, monomer_start_atoms,
                                  fmt::format("subchain-join {:02d}", sc));

            validateTopology(polymer, atom_monomer_id, interface_bonds,
                             fmt::format("subchain-join {:02d}", sc));
            repairUnboundAtoms(polymer, atom_monomer_id, interface_bonds,
                               fmt::format("subchain-join {:02d}", sc));
            // Claude Generated: Validate topology consistency across identical monomers
            validateTopologyConsistency(polymer, atom_monomer_id, monomer_fragment_type,
                interface_bonds, fmt::format("subchain-join {:02d}", sc));

            // Light FF optimization after sub-chain connection
            resolveOverlaps(polymer, interface_bonds);

            Molecule polymer_opt = polymer;
            polymer_opt.setCharge(0);
            Mol m_opt = polymer_opt.getMolInfo();
            for (int k = 0; k < polymer_opt.AtomCount(); ++k) {
                if (m_opt.m_atoms[k] == 0)
                    m_opt.m_atoms[k] = cap_intermediate_elem;
            }
            if (polymer.hasPersistentTopology()) {
                const Matrix topo = polymer.getTopologyMatrix();
                m_opt.m_bonds.clear();
                int n = topo.rows();
                for (int ii = 0; ii < n; ++ii)
                    for (int jj = ii + 1; jj < n; ++jj)
                        if (topo(ii, jj) > 0.5)
                            m_opt.m_bonds.push_back({ii, jj});
            }
            polymer_opt.LoadMolecule(m_opt);
            polymer_opt.setCharge(0);

            CurcumaLogger::info(fmt::format("Running {} optimization after sub-chain join...", opt_method));
            int saved_verbosity = CurcumaLogger::get_verbosity();
            int opt_max_iter = m_config.get<int>("opt_max_iter", 0);
            json opt_ctrl;
            opt_ctrl["opt"]["method"] = opt_method;
            opt_ctrl["opt"]["printOutput"] = false;
            opt_ctrl["opt"]["writeXYZ"] = true;
            if (opt_max_iter > 0)
                opt_ctrl["opt"]["MaxIter"] = opt_max_iter;
            CurcumaOpt opt(opt_ctrl, true);
            opt.overrideBasename(fmt::format("polymer_scjoin_{:02d}", sc));
            opt.addMolecule(polymer_opt);
            opt.start();
            CurcumaLogger::set_verbosity(saved_verbosity);
            const std::vector<Molecule>* finals = opt.Molecules();
            if (finals && !finals->empty()) {
                Geometry result_geom = finals->back().Coords();
                if (!containsNaN(result_geom))
                    polymer.setGeometry(result_geom);
                else
                    CurcumaLogger::warn("Sub-chain join optimization returned NaN — keeping pre-opt geometry");
            }

            checkDistances(polymer,
                           fmt::format("subchain-join {:02d} post-opt", sc),
                           interface_bonds);
        }

        // Phase 3: Capping + final output
        polymer.setName("polymer_A");
        applyCapping(polymer);
        polymer.writeXYZFile(polymer.Name() + "_final.xyz");
        CurcumaLogger::success(fmt::format(
            "D&C polymer assembly completed: {}_final.xyz ({} atoms)",
            polymer.Name(), polymer.AtomCount()));

    } else {
        // ====== Sequential assembly (original behavior) ======
        SubchainResult result = buildSubchain(sequence, 0, 0);
        Molecule polymer = result.molecule;
        polymer.setName("polymer_A");
        applyCapping(polymer);
        polymer.writeXYZFile(polymer.Name() + "_final.xyz");
        CurcumaLogger::success("Polymer assembly completed: " + polymer.Name() + "_final.xyz");
    }
}

void PolymerBuild::applyCapping(Molecule& mol)
{
    // cap sets both ends uniformly; cap_start / cap_end override per end.
    // Empty cap_start or cap_end falls back to cap.
    std::string cap_default = m_config.get<std::string>("cap", "H");
    std::string cap_start_str = m_config.get<std::string>("cap_start", "");
    std::string cap_end_str   = m_config.get<std::string>("cap_end",   "");
    if (cap_start_str.empty()) cap_start_str = cap_default;
    if (cap_end_str.empty())   cap_end_str   = cap_default;

    std::vector<int> xx_indices;
    for (int i = 0; i < mol.AtomCount(); ++i)
        if (mol.Atom(i).first == 0) xx_indices.push_back(i);

    if (xx_indices.empty()) return;

    CurcumaLogger::info(fmt::format(
        "Capping {} chain end(s): start={}, end={}",
        xx_indices.size(), cap_start_str, cap_end_str));

    Mol m = mol.getMolInfo();
    for (size_t i = 0; i < xx_indices.size(); ++i) {
        // First Xx → cap_start, last Xx → cap_end, any middle ones → cap_default
        std::string cap_str = (i == 0) ? cap_start_str
                            : (i == xx_indices.size() - 1) ? cap_end_str
                            : cap_default;
        if (cap_str == "none") continue;

        int element = Elements::String2Element(cap_str);
        if (element <= 0) element = 1;  // fallback to H for unrecognised strings
        m.m_atoms[xx_indices[i]] = element;
    }
    mol.LoadMolecule(m);
}
