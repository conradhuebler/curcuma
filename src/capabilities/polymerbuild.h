/*
 * <Capability for assembling polymers from fragments.>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated: Implementation of iterative polymer assembly with topological matrix support
 */

#pragma once

#include "curcumamethod.h"
#include "src/core/config_manager.h"
#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"

#include <string>
#include <vector>
#include <map>

/**
 * @brief Result of connecting a fragment to a polymer.
 *
 * Claude Generated: Tracks all Xx connection points and their active atom assignments.
 */
struct ConnectionResult {
    Molecule polymer;                                      ///< The polymer with the new fragment added
    std::vector<std::pair<int, int>> tracked_xx;           ///< All (Xx index, active atom index) pairs for chain extension
    std::vector<std::pair<int, int>> interface_bonds;      ///< Bonds created at interface (for topology)
};

/**
 * @brief PolymerBuild Capability
 *
 * Assembles complex polymers from fragment structures using a sequence string.
 * Connection points are defined by dummy atoms (Xx) in the fragments.
 *
 * Claude Generated: Uses Levenberg-Marquardt optimization for fragment placement
 * with bond-length constraints, replacing the previous grid-scan approach.
 */
class PolymerBuild : public CurcumaMethod {
public:
    PolymerBuild(const json& controller, bool silent);
    ~PolymerBuild();

    /**
     * @brief Start the polymer assembly process
     */
    void start() override;

    /**
     * @brief Print help for polymer build options
     */
    void printHelp() const override;

    // CurcumaMethod interface implementation
    nlohmann::json WriteRestartInformation() override { return json{}; }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return { "PolymerBuild" }; }
    void ReadControlFile() override {}
    void LoadControlJson() override {}

private:
    /**
     * @brief Parse the polymer sequence string
     * @param sequence The sequence string (e.g., (A)10-B-(C)5)
     * @return List of fragment names in order
     */
    std::vector<std::string> parseSequence(const std::string& sequence);

    /**
     * @brief Assemble the polymer from the sequence of fragments
     * @param sequence List of fragment names
     */
    void assemblePolymer(const std::vector<std::string>& sequence);

    /**
     * @brief Connect a new fragment to the current polymer using LM optimization.
     *
     * Uses Levenberg-Marquardt optimization with 6 DOF (3 translation + 3 rotation)
     * and a bond-length constraint to find optimal fragment placement.
     * The steric energy is minimized using LJ potential between polymer and fragment atoms.
     *
     * @param current_polymer The existing polymer (will NOT be modified)
     * @param fragment_file Path to the fragment XYZ file
     * @param prev_tracked_xx Previously tracked Xx atoms (from last fragment)
     * @param prev_interface_bonds Previous interface bonds (for topology continuity)
     * @param step_number Fragment connection step (for debug output filenames)
     * @return ConnectionResult with new polymer, tracked Xx, and interface bonds
     *
     * Claude Generated
     */
    ConnectionResult connectFragment(
        const Molecule& current_polymer,
        const std::string& fragment_file,
        const std::vector<std::pair<int, int>>& prev_tracked_xx,
        const std::vector<std::pair<int, int>>& prev_interface_bonds,
        int step_number = 0);

    /**
     * @brief Optimize fragment placement using Levenberg-Marquardt.
     *
     * Minimizes LJ steric energy between polymer and fragment atoms, with a
     * constraint on the interface bond length (active_out to active_in distance).
     *
     * @param polymer Heavy atoms of polymer (no Xx)
     * @param fragment Fragment molecule (no Xx)
     * @param anchor_polymer Position of active_out atom in polymer
     * @param anchor_fragment_idx Index of active_in atom in fragment
     * @param bond_length Target bond length (sum of covalent radii)
     * @param initial_translation Starting translation (fragment centroid displacement)
     * @param initial_rotation Starting rotation angles (degrees)
     * @param step_number Fragment connection step (for debug output filenames)
     * @return Optimized (translation, rotation) pair
     *
     * Claude Generated
     */
    std::pair<Position, Position> optimizeFragmentPlacement(
        const Molecule& polymer,
        const Molecule& fragment,
        const Position& anchor_polymer,
        int anchor_polymer_idx,
        int anchor_fragment_idx,
        double bond_length,
        const Position& initial_translation,
        const Position& initial_rotation,
        int step_number = 0,
        int polymer_xx_idx = -1,
        int frag_xx_idx = -1);

    /**
     * @brief Apply capping to the terminal connection points
     * @param mol The molecule to cap
     */
    void applyCapping(Molecule& mol);

    /**
     * @brief Pre-relax clashing atoms with a pure pairwise repulsion SD.
     *
     * Resolves short contacts that make UFF/LBFGS diverge on step 1.
     * Uses only (r_cut/r)^4 repulsion — no bonded terms.
     * Claude Generated
     *
     * @param mol               Molecule to relax in-place
     * @param max_steps         Maximum steepest-descent iterations (default 300)
     * @param max_displacement  Max atom displacement per step in Å (default 0.05)
     */
    void resolveOverlaps(Molecule& mol, int max_steps = 300, double max_displacement = 0.05);

    /**
     * @brief Check all pairwise distances and report close contacts.
     *
     * Categories:
     *   - "FUSION"  : dist < 0.5 * r_cov_sum  → FF will certainly crash
     *   - "CLASH"   : dist < 0.75 * r_cov_sum → likely problematic
     *
     * Bonded pairs (within interface_bonds) are excluded from non-bonded checks.
     *
     * @param mol           Molecule to check
     * @param tag           Label printed in log messages
     * @param interface_bonds  Known bonded pairs to skip
     * @return Number of problematic contacts found
     *
     * Claude Generated
     */
    int checkDistances(const Molecule& mol,
                       const std::string& tag,
                       const std::vector<std::pair<int,int>>& interface_bonds = {}) const;

    /**
     * @brief Find the heavy atom bonded to an Xx via topology analysis.
     *
     * Uses distance-based topology to find the atom that Xx connects to.
     * Falls back to closest heavy atom heuristic if topology gives no result.
     *
     * @param mol The molecule
     * @param xx_idx Index of Xx atom
     * @return Index of bonded heavy atom, or -1 if not found
     *
     * Claude Generated
     */
    int findBondedAtom(const Molecule& mol, int xx_idx) const;

    ConfigManager m_config;
    bool m_silent;
    std::map<std::string, std::string> m_fragments;

    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    BEGIN_PARAMETER_DEFINITION(polymerbuild)

    // Sequence options
    PARAM(sequence, String, "", "Polymer sequence string (e.g., (A)10-B-(C)5)", "Assembly", { "seq" })
    PARAM(fragments, String, "{}", "JSON map of fragment names to file paths", "Assembly", { "frag" })
    PARAM(cap, String, "H", "Capping element for all chain ends (H, F, Cl, none) — overridden per end by cap_start/cap_end", "Assembly", { "c" })
    PARAM(cap_start, String, "", "Override cap for the chain start (defaults to cap)", "Assembly", {})
    PARAM(cap_end, String, "", "Override cap for the chain end (defaults to cap)", "Assembly", {})
    PARAM(cap_intermediate, String, "H", "Element substituted for Xx during intermediate FF optimization (default H; use I or F for stronger LJ repulsion at chain ends)", "Refinement", { "ci" })

    // Optimization options
    PARAM(optimize, Bool, true, "Enable intermediate optimization steps", "Refinement", { "opt" })
    PARAM(opt_method, String, "gfnff", "Method for optimization (gfnff, uff, gfn2, etc.)", "Refinement", {})
    PARAM(lm_max_iter, Int, 500, "Maximum LM iterations for fragment placement", "Refinement", {})
    PARAM(lm_tolerance, Double, 1e-6, "Convergence tolerance for LM optimization", "Refinement", {})

    // Dynamics options
    PARAM(dynamics, Bool, false, "Enable intermediate molecular dynamics", "Refinement", { "md" })
    PARAM(md_steps, Int, 1000, "Number of MD steps for intermediate dynamics", "Refinement", {})
    PARAM(md_method, String, "gfnff", "Method for dynamics", "Refinement", {})

    // Geometric options
    PARAM(bond_distance_scaling, Double, 1.0, "Scaling factor for interface bond lengths", "Assembly", { "scaling" })
    PARAM(chunk_size, Int, 1, "Number of fragments to place before running optimization (1=after each)", "Refinement", {})
    PARAM(write_intermediates, Bool, false, "Write XYZ file after each fragment addition for debugging", "Output", { "intermediates" })
    PARAM(verbose, Bool, true, "Detailed output", "Output", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};