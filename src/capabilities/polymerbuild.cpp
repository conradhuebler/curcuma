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
        verbosity = 3;
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

    std::vector<std::string> sequence = parseSequence(sequence_str);
    if (sequence.empty()) {
        CurcumaLogger::error("Parsed sequence is empty or invalid");
        return;
    }

    assemblePolymer(sequence);
}

void PolymerBuild::printHelp() const
{
    // Handled by ParameterRegistry
}

std::vector<std::string> PolymerBuild::parseSequence(const std::string& sequence)
{
    std::vector<std::string> result;
    // Simple parser for (A)n and A-B-C
    std::regex re("(\\(([A-Za-z0-9_]+)\\)([0-9]+))|([A-Za-z0-9_]+)");
    auto it = std::sregex_iterator(sequence.begin(), sequence.end(), re);
    auto end = std::sregex_iterator();

    for (; it != end; ++it) {
        std::smatch match = *it;
        if (match[1].matched) { // (NAME)NUMBER
            std::string name = match[2].str();
            int count = std::stoi(match[3].str());
            for (int i = 0; i < count; ++i) {
                result.push_back(name);
            }
        } else if (match[4].matched) { // NAME
            result.push_back(match[4].str());
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
    std::string traj_filename = fmt::format("lm_trajectory_step{:02d}.xyz", step_number);
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

ConnectionResult PolymerBuild::connectFragment(
    const Molecule& current_polymer,
    const std::string& fragment_file,
    const std::vector<std::pair<int, int>>& prev_tracked_xx,
    const std::vector<std::pair<int, int>>& prev_interface_bonds,
    int step_number)
{
    ConnectionResult result;

    // Load next fragment
    Molecule next;
    next.LoadMolecule(fragment_file);

    // Find Xx atoms in polymer (from prev_tracked_xx)
    if (prev_tracked_xx.empty()) {
        CurcumaLogger::error("connectFragment: No tracked Xx atoms available for connection");
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

    // Find best-matching Xx in next fragment
    // Use direction-based matching: pick Xx whose bond direction best aligns
    // with -v_out (opposite to polymer's Xx→active direction)
    Position polymer_xx_pos = current_polymer.Atom(polymer_xx_idx).second;
    Position polymer_active_pos = current_polymer.Atom(polymer_active_idx).second;
    Position v_out = (polymer_xx_pos - polymer_active_pos).normalized();

    int best_next_xx = next_xx.front();
    double best_dot = -2.0;
    for (int xx : next_xx) {
        int active = findBondedAtom(next, xx);
        if (active < 0) continue;
        Position v_in = (next.Atom(xx).second - next.Atom(active).second).normalized();
        double dot = v_in.dot(-v_out);
        if (dot > best_dot) {
            best_dot = dot;
            best_next_xx = xx;
        }
    }

    int next_xx_idx = best_next_xx;
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
    result.interface_bonds = prev_interface_bonds;
    result.interface_bonds.push_back({polymer_active_idx_lm, offset + next_active_idx_lm});

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

    result.polymer = combined;
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
                // Non-bonded pair at bonding distance → topology inconsistency / nuclear fusion
                CurcumaLogger::error(fmt::format(
                    "{}: FUSION {}{}–{}{}: {:.3f} Å (r_cov {:.3f}, not in topology)",
                    tag, Elements::ElementAbbr[ei], i, Elements::ElementAbbr[ej], j, d, r_cov));
                ++problems;
            } else if (d < r_vdw * 0.65) {
                // Non-bonded pair closer than 65 % of vdW sum → significant steric clash
                // (threshold chosen to avoid false positives for 1,3-pairs: e.g. H-C-H ~1.8 Å
                //  vs. 0.65 * vdW(H+H) = 0.65 * 2.40 = 1.56 Å — safely below typical 1,3 distances)
                CurcumaLogger::warn(fmt::format(
                    "{}: clash {}{}–{}{}: {:.3f} Å (0.65·r_vdw {:.3f})",
                    tag, Elements::ElementAbbr[ei], i, Elements::ElementAbbr[ej], j, d, r_vdw * 0.65));
                ++problems;
            }
        }
    }
    if (problems == 0)
        CurcumaLogger::info(fmt::format("{}: all distances OK", tag));
    else
        CurcumaLogger::warn(fmt::format("{}: {} problematic contact(s)", tag, problems));
    return problems;
}

void PolymerBuild::resolveOverlaps(Molecule& mol, int max_steps, double max_displacement)
{
    int n = mol.AtomCount();
    int clashes_initial = 0;

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

                Position diff = coords.row(i).transpose() - coords.row(j).transpose();
                double dist = diff.norm();
                if (dist < 1e-6) dist = 1e-6;

                double r_cov = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                double r_cut = r_cov * 0.85;
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

void PolymerBuild::assemblePolymer(const std::vector<std::string>& sequence)
{
    if (sequence.empty()) return;

    // Load first fragment
    std::string first_name = sequence[0];
    if (m_fragments.find(first_name) == m_fragments.end()) {
        CurcumaLogger::error("Fragment not found in map: " + first_name);
        return;
    }

    Molecule polymer;
    polymer.LoadMolecule(m_fragments[first_name]);
    polymer.setName("polymer_" + first_name);
    polymer.setCharge(0);  // Molecule() = default → m_charge uninitialized; GFN-FF needs charge = 0

    // Debug: Show raw first fragment with Xx
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info(fmt::format("DEBUG STEP 1: Raw first fragment '{}' ({} atoms, with Xx)",
                                        first_name, polymer.AtomCount()));
        std::string debug_file = fmt::format("{}_step1_00_raw_Xx.xyz", polymer.Name());
        polymer.writeXYZFile(debug_file);
        CurcumaLogger::info(fmt::format("DEBUG: Wrote {}", debug_file));

        Molecule capped = polymer;
        Mol capped_info = capped.getMolInfo();
        for (int k = 0; k < capped.AtomCount(); ++k) {
            if (capped_info.m_atoms[k] == 0) capped_info.m_atoms[k] = 1;
        }
        if (capped_info.m_charge > 1000 || capped_info.m_charge < -1000) {
            capped_info.m_charge = 0;
        }
        capped.LoadMolecule(capped_info);
        std::string capped_file = fmt::format("{}_step1_01_capped_H.xyz", polymer.Name());
        capped.writeXYZFile(capped_file);
        CurcumaLogger::info(fmt::format("DEBUG: Wrote (Xx→H capped) {}", capped_file));
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
    int chunk_size = m_config.get<int>("chunk_size", 1);
    if (chunk_size < 1) chunk_size = 1;
    bool optimize = m_config.get<bool>("optimize", true);

    // Connect remaining fragments
    for (size_t i = 1; i < sequence.size(); ++i) {
        std::string next_name = sequence[i];
        if (m_fragments.find(next_name) == m_fragments.end()) {
            CurcumaLogger::error("Fragment not found in map: " + next_name);
            continue;
        }

        CurcumaLogger::info(fmt::format("Adding fragment {}/{}: {}", i, sequence.size() - 1, next_name));

        // Use LM-based connectFragment (pass step_number for debug output)
        ConnectionResult conn_result = connectFragment(
            polymer, m_fragments[next_name], tracked_xx, interface_bonds, static_cast<int>(i));

        polymer = conn_result.polymer;
        tracked_xx = conn_result.tracked_xx;
        interface_bonds = conn_result.interface_bonds;

        checkDistances(polymer,
                       fmt::format("step {:02d} post-connect", i),
                       interface_bonds);

        // Debug output
        if (CurcumaLogger::get_verbosity() >= 3) {
            int xx_count = 0;
            for (int k = 0; k < polymer.AtomCount(); ++k)
                if (polymer.Atom(k).first == 0) xx_count++;
            CurcumaLogger::info(fmt::format("DEBUG: After connectFragment - {} atoms, {} Xx remaining",
                                            polymer.AtomCount(), xx_count));
            std::string debug_file = fmt::format("{}_debug_{:02d}_connected.xyz", polymer.Name(), i);
            polymer.writeXYZFile(debug_file);
            CurcumaLogger::info(fmt::format("DEBUG: Wrote {}", debug_file));
        }

        // Rebuild topology
        {
            int n = polymer.AtomCount();
            Matrix topo = Matrix::Zero(n, n);
            auto [dist, _unused] = polymer.DistanceMatrix();
            for (int ii = 0; ii < n; ++ii) {
                for (int jj = ii + 1; jj < n; ++jj) {
                    int ei = polymer.Atom(ii).first;
                    int ej = polymer.Atom(jj).first;
                    if (ei == 0 || ej == 0) continue;
                    double cutoff = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                    if (dist(ii, jj) > 1e-3 && dist(ii, jj) <= cutoff * 1.15) {
                        topo(ii, jj) = topo(jj, ii) = 1.0;
                    }
                }
            }
            for (const auto& bond : interface_bonds) {
                if (bond.first < n && bond.second < n)
                    topo(bond.first, bond.second) = topo(bond.second, bond.first) = 1.0;
            }
            polymer.setTopologyMatrix(topo);
        }

        // Optimization at chunk boundaries (but NOT for last fragment - post-capping handles that)
        bool should_optimize = (static_cast<int>(i) % chunk_size == 0)
                               && (i != sequence.size() - 1);

        if (optimize && should_optimize) {
            if (CurcumaLogger::get_verbosity() >= 3) {
                CurcumaLogger::info(fmt::format("DEBUG: Before resolveOverlaps (fragment {}), {} atoms",
                                                i, polymer.AtomCount()));
            }

            resolveOverlaps(polymer);
            checkDistances(polymer,
                           fmt::format("step {:02d} post-resolveOverlaps", i),
                           interface_bonds);

            if (CurcumaLogger::get_verbosity() >= 3) {
                std::string debug_file = fmt::format("{}_debug_{:02d}_postLJ_preH.xyz", polymer.Name(), i);
                polymer.writeXYZFile(debug_file);
                CurcumaLogger::info(fmt::format("DEBUG: Wrote {}", debug_file));
            }

            CurcumaLogger::info(fmt::format("Running UFF optimization (after fragment {})...", i));
            json opt_ctrl;
            opt_ctrl["opt"]["method"] = m_config.get<std::string>("opt_method", "uff");
            opt_ctrl["opt"]["printOutput"] = false;
            opt_ctrl["opt"]["writeXYZ"] = false;

            // Replace Xx with cap_intermediate for optimization; enforce charge = 0 (GFN-FF EEQ requires it)
            std::string ci_str = m_config.get<std::string>("cap_intermediate", "H");
            int cap_intermediate_elem = Elements::String2Element(ci_str);
            if (cap_intermediate_elem <= 0) cap_intermediate_elem = 1;

            Molecule polymer_opt = polymer;
            polymer_opt.setCharge(0);
            std::vector<int> xx_indices;
            Mol m_opt = polymer_opt.getMolInfo();
            for (int k = 0; k < polymer_opt.AtomCount(); ++k) {
                if (m_opt.m_atoms[k] == 0) {
                    xx_indices.push_back(k);
                    m_opt.m_atoms[k] = cap_intermediate_elem;
                }
            }
            polymer_opt.LoadMolecule(m_opt);
            polymer_opt.setCharge(0);  // LoadMolecule may overwrite from Mol, reset again

            if (CurcumaLogger::get_verbosity() >= 3) {
                std::string debug_file = fmt::format("{}_debug_{:02d}_forFF.xyz", polymer.Name(), i);
                polymer_opt.writeXYZFile(debug_file);
                CurcumaLogger::info(fmt::format("DEBUG: Wrote {}", debug_file));
            }

            // CurcumaOpt(silent=true) resets global CurcumaLogger verbosity → save/restore
            int saved_verbosity = CurcumaLogger::get_verbosity();
            CurcumaOpt opt(opt_ctrl, true);
            opt.addMolecule(polymer_opt);
            opt.start();
            CurcumaLogger::set_verbosity(saved_verbosity);
            const std::vector<Molecule>* finals = opt.Molecules();
            if (finals && !finals->empty()) {
                const Molecule& opt_mol = finals->back();
                // Copy ALL optimized coordinates back, including Xx placeholder positions.
                // Xx were replaced by H in polymer_opt; the optimized H position is a
                // better estimate for v_out in the next LM step than the stale pre-FF position.
                // Atom element types in polymer are NOT changed (Xx remains Xx).
                polymer.setGeometry(opt_mol.Coords());

                if (CurcumaLogger::get_verbosity() >= 3) {
                    std::string debug_file = fmt::format("{}_debug_{:02d}_postFF.xyz", polymer.Name(), i);
                    polymer.writeXYZFile(debug_file);
                    CurcumaLogger::info(fmt::format("DEBUG: Wrote {}", debug_file));
                }

                // Restore topology
                int n = polymer.AtomCount();
                Matrix topo = Matrix::Zero(n, n);
                auto [dist2, _u2] = polymer.DistanceMatrix();
                for (int ii = 0; ii < n; ++ii)
                    for (int jj = ii+1; jj < n; ++jj) {
                        int ei = polymer.Atom(ii).first, ej = polymer.Atom(jj).first;
                        if (ei == 0 || ej == 0) continue;
                        double co = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                        if (dist2(ii,jj) > 1e-3 && dist2(ii,jj) <= co * 1.15)
                            topo(ii,jj) = topo(jj,ii) = 1.0;
                    }
                for (const auto& b : interface_bonds)
                    if (b.first < n && b.second < n)
                        topo(b.first, b.second) = topo(b.second, b.first) = 1.0;
                polymer.setTopologyMatrix(topo);

                checkDistances(polymer,
                               fmt::format("step {:02d} post-FF", i),
                               interface_bonds);
            }
        }

        if (m_config.get<bool>("dynamics", false) && should_optimize) {
            CurcumaLogger::info("Running intermediate dynamics...");
            json md_ctrl;
            md_ctrl["md"]["max_time"] = m_config.get<int>("md_steps", 1000);
            md_ctrl["md"]["method"] = m_config.get<std::string>("md_method", "uff");
            md_ctrl["md"]["verbosity"] = 0;

            SimpleMD md(md_ctrl, true);
            md.setMolecule(polymer);
            md.Initialise();
            md.start();
            auto unique = md.UniqueMolecules();
            if (!unique.empty()) {
                polymer = *unique.back();
            }
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

    applyCapping(polymer);

    // Post-capping optimization
    if (m_config.get<bool>("optimize", true)) {
        resolveOverlaps(polymer);

        CurcumaLogger::info("Running post-capping optimization...");
        json opt_ctrl;
        opt_ctrl["opt"]["method"] = m_config.get<std::string>("opt_method", "uff");
        opt_ctrl["opt"]["printOutput"] = false;
        opt_ctrl["opt"]["writeXYZ"] = false;

        Molecule polymer_opt = polymer;
        polymer_opt.setCharge(0);
        // cap=="none" may leave residual Xx — substitute them for the FF run
        std::string ci_str2 = m_config.get<std::string>("cap_intermediate", "H");
        int cap_intermediate_elem2 = Elements::String2Element(ci_str2);
        if (cap_intermediate_elem2 <= 0) cap_intermediate_elem2 = 1;
        std::vector<int> xx_indices;
        Mol m_opt = polymer_opt.getMolInfo();
        for (int i = 0; i < polymer_opt.AtomCount(); ++i) {
            if (m_opt.m_atoms[i] == 0) {
                xx_indices.push_back(i);
                m_opt.m_atoms[i] = cap_intermediate_elem2;
            }
        }
        polymer_opt.LoadMolecule(m_opt);

        int saved_verbosity_cap = CurcumaLogger::get_verbosity();
        CurcumaOpt opt(opt_ctrl, true);
        opt.addMolecule(polymer_opt);
        opt.start();
        CurcumaLogger::set_verbosity(saved_verbosity_cap);
        const std::vector<Molecule>* finals = opt.Molecules();
        if (finals && !finals->empty()) {
            const Molecule& opt_mol = finals->back();
            Geometry opt_coords = opt_mol.Coords();
            Geometry polymer_coords = polymer.Coords();
            for (int i = 0; i < polymer.AtomCount() && i < opt_mol.AtomCount(); ++i) {
                if (polymer.Atom(i).first != 0) {
                    polymer_coords.row(i) = opt_coords.row(i);
                }
            }
            polymer.setGeometry(polymer_coords);
        }
    }

    // Final topology rebuild
    {
        int n = polymer.AtomCount();
        Matrix topo = Matrix::Zero(n, n);
        auto [dist, _unused] = polymer.DistanceMatrix();
        for (int ii = 0; ii < n; ++ii) {
            for (int jj = ii + 1; jj < n; ++jj) {
                int ei = polymer.Atom(ii).first;
                int ej = polymer.Atom(jj).first;
                if (ei == 0 || ej == 0) continue;
                double cutoff = Elements::CovalentRadius[ei] + Elements::CovalentRadius[ej];
                if (dist(ii, jj) > 1e-3 && dist(ii, jj) <= cutoff * 1.15) {
                    topo(ii, jj) = topo(jj, ii) = 1.0;
                }
            }
        }
        for (const auto& bond : interface_bonds) {
            if (bond.first < n && bond.second < n)
                topo(bond.first, bond.second) = topo(bond.second, bond.first) = 1.0;
        }
        polymer.setTopologyMatrix(topo);
    }

    polymer.writeXYZFile(polymer.Name() + "_final.xyz");
    CurcumaLogger::success("Polymer assembly completed: " + polymer.Name() + "_final.xyz");
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