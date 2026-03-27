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

/// Claude Generated: Sequence entry with explicit Xx connection point selection via prime notation
struct SequenceEntry {
    std::string fragment_name;  ///< Fragment key (without primes) — key into m_fragments
    int xx_selection;           ///< 0=first Xx, 1=second Xx, 2=third, ... (from prime count)
};

/**
 * @brief Result of connecting a fragment to a polymer.
 *
 * Claude Generated: Tracks all Xx connection points and their active atom assignments.
 */
struct ConnectionResult {
    Molecule polymer;                                      ///< The polymer with the new fragment added
    std::vector<std::pair<int, int>> tracked_xx;           ///< All (Xx index, active atom index) pairs for chain extension
    std::vector<std::pair<int, int>> interface_bonds;      ///< Bonds created at interface (for topology)
    int fragment_offset = 0;                               ///< Start index of new fragment atoms in combined polymer
    int removed_polymer_xx_idx = -1;                       ///< Index of removed interface Xx in original polymer (for atom_monomer_id update)
    int removed_fragment_xx_idx = -1;                      ///< Index of removed interface Xx in original fragment (for D&C monomer ID mapping)
};

/// Claude Generated: Result of building a sub-chain for divide-and-conquer assembly
struct SubchainResult {
    Molecule molecule;
    std::vector<std::pair<int, int>> tracked_xx;
    std::vector<std::pair<int, int>> interface_bonds;
    std::vector<int> atom_monomer_id;
    std::vector<std::string> monomer_fragment_type;  ///< Claude Generated: fragment name for each monomer index
    std::vector<int> monomer_start_atoms;             ///< Claude Generated: starting atom index for each monomer
};

/**
 * @brief Template for fragment topology (bonds within fragment, excluding Xx).
 *
 * Claude Generated: Stores expected intra-monomer bonds for topology consistency validation.
 * Used to verify that identical monomers have identical internal topology.
 */
struct FragmentTopologyTemplate {
    std::string fragment_name;                          ///< Fragment identifier
    int n_atoms;                                         ///< Number of heavy atoms (excluding Xx)
    std::vector<std::pair<int, int>> internal_bonds;    ///< Bonds between heavy atoms (indices relative to fragment, Xx removed)
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

    friend class TestablePolymerBuild;  // Claude Generated: test access to parsing methods

private:
    /**
     * @brief Parse the polymer sequence string
     * @param sequence The sequence string (e.g., (A)10-B-(C)5)
     * @return List of sequence entries with fragment names and Xx selection
     */
    std::vector<SequenceEntry> parseSequence(const std::string& sequence);

    /// Claude Generated: Parse inner token string into SequenceEntry list (e.g., "AA'" → [{A,0},{A,1}])
    std::vector<SequenceEntry> parseTokens(const std::string& content);

    /**
     * @brief Assemble the polymer from the sequence of fragments
     * @param sequence List of fragment names
     */
    void assemblePolymer(const std::vector<SequenceEntry>& sequence);

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
        int step_number = 0,
        int xx_selection = 0);

    /// Claude Generated: Core connection logic accepting a Molecule directly (used by connectFragment and D&C merge)
    ConnectionResult connectMolecule(
        const Molecule& current_polymer,
        const Molecule& next_fragment,
        const std::vector<std::pair<int, int>>& prev_tracked_xx,
        const std::vector<std::pair<int, int>>& prev_interface_bonds,
        int step_number = 0,
        const std::vector<std::pair<int, int>>& fragment_internal_bonds = {},
        int xx_selection = 0);

    /// Claude Generated: Build a sub-chain from a list of sequence entries (for divide-and-conquer)
    SubchainResult buildSubchain(
        const std::vector<SequenceEntry>& fragment_entries,
        int monomer_id_offset,
        int subchain_index);

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
     * @brief Validate topology: check that no spurious cross-block bonds exist.
     *
     * Verifies that the topology matrix contains only intra-block bonds and
     * explicit interface bonds. Reports any spurious cross-monomer bonds.
     *
     * @return Number of spurious bonds found (0 = topology correct)
     *
     * Claude Generated
     */
    int validateTopology(const Molecule& mol,
                         const std::vector<int>& atom_monomer_id,
                         const std::vector<std::pair<int,int>>& interface_bonds,
                         const std::string& tag) const;

    /**
     * @brief Repair unbound atoms: find expected bond partner, enforce bond, move atom.
     *
     * For each non-Xx atom with zero bonds in the topology matrix:
     * 1. Find the best bond partner (smallest dist/cov_sum ratio within same monomer)
     * 2. Set the bond in the topology matrix
     * 3. Move the atom to the correct covalent bond distance from the partner
     *
     * @param mol               Molecule to repair in-place
     * @param atom_monomer_id   Per-atom monomer ID assignment
     * @param interface_bonds   Known interface bonds (for partner search across monomers)
     * @param tag               Label for log messages
     * @return Number of repaired atoms
     *
     * Claude Generated
     */
    int repairUnboundAtoms(Molecule& mol,
                           const std::vector<int>& atom_monomer_id,
                           const std::vector<std::pair<int,int>>& interface_bonds,
                           const std::string& tag);

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

    /**
     * @brief Extract fragment topology template from a loaded molecule.
     *
     * Extracts internal bonds (excluding Xx) for use as template during topology rebuild.
     * Indices are relative to the fragment (Xx atoms removed).
     *
     * @param mol Fragment molecule with Xx atoms still present
     * @param fragment_name Name of the fragment for storage
     * @return FragmentTopologyTemplate with bonds and metadata
     *
     * Claude Generated
     */
    FragmentTopologyTemplate extractFragmentTemplate(const Molecule& mol, const std::string& fragment_name) const;

    /**
     * @brief Rebuild topology from stored templates instead of distance-based detection.
     *
     * For each monomer, copies bonds from the stored template (adjusted for atom indices).
     * Interface bonds are added explicitly from the interface_bonds list.
     *
     * @param mol Polymer molecule
     * @param atom_monomer_id Per-atom monomer ID assignment
     * @param monomer_fragment_type Fragment name for each monomer index
     * @param interface_bonds Known interface bonds
     * @param monomer_start_atoms Starting atom index for each monomer
     *
     * Claude Generated
     */
    void rebuildTopologyFromTemplates(
        Molecule& mol,
        const std::vector<int>& atom_monomer_id,
        const std::vector<std::string>& monomer_fragment_type,
        const std::vector<std::pair<int, int>>& interface_bonds,
        const std::vector<int>& monomer_start_atoms) const;

    /**
     * @brief Validate topology consistency across identical monomers.
     *
     * Compares the actual intra-monomer bond counts against the expected counts
     * from stored templates. Reports discrepancies.
     *
     * @param mol Polymer molecule
     * @param atom_monomer_id Per-atom monomer ID assignment
     * @param monomer_fragment_type Fragment name for each monomer index
     * @param interface_bonds Known interface bonds (to exclude from count)
     * @param tag Label for log messages
     * @return Number of inconsistencies found (0 = all consistent)
     *
     * Claude Generated
     */
    int validateTopologyConsistency(
        const Molecule& mol,
        const std::vector<int>& atom_monomer_id,
        const std::vector<std::string>& monomer_fragment_type,
        const std::vector<std::pair<int, int>>& interface_bonds,
        const std::string& tag) const;

    ConfigManager m_config;
    bool m_silent;
    std::map<std::string, std::string> m_fragments;
    std::map<std::string, FragmentTopologyTemplate> m_fragment_templates;  ///< Claude Generated: Cached topology templates by fragment name

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
    PARAM(opt_max_iter, Int, 0, "Max iterations for FF optimization (0=until convergence)", "Refinement", {})
    PARAM(lm_max_iter, Int, 500, "Maximum LM iterations for fragment placement", "Refinement", {})
    PARAM(lm_tolerance, Double, 1e-6, "Convergence tolerance for LM optimization", "Refinement", {})
    PARAM(overlap_retries, Int, 3, "Max optimization retries to resolve cross-monomer overlaps", "Refinement", { "retries" })

    // Dynamics options
    PARAM(dynamics, Bool, false, "Enable intermediate molecular dynamics after each fragment", "Refinement", { "md" })
    PARAM(md_steps, Int, 1000, "Max simulation time for intermediate MD [fs]", "Refinement", {})
    PARAM(md_temperature, Double, 300.0, "Temperature for intermediate MD [K]", "Refinement", {})
    PARAM(md_time_step, Double, 1.0, "Integration time step for intermediate MD [fs]", "Refinement", {})
    PARAM(md_method, String, "gfnff", "Method for dynamics", "Refinement", {})

    // Cold MD fallback (overlap resolution)
    PARAM(cold_md_temperature, Double, 10.0, "Temperature for cold MD overlap resolution [K]", "Refinement", {})
    PARAM(cold_md_time, Double, 200.0, "Max simulation time for cold MD [fs]", "Refinement", {})
    PARAM(cold_md_time_step, Double, 0.5, "Integration time step for cold MD [fs]", "Refinement", {})

    // Geometric options
    PARAM(bond_distance_scaling, Double, 1.0, "Scaling factor for interface bond lengths", "Assembly", { "scaling" })
    PARAM(chunk_size, Int, 1, "Number of fragments to place before running optimization (1=after each)", "Refinement", {})
    PARAM(subchain_size, Int, 0, "Sub-chain size for divide-and-conquer assembly (0=sequential)", "Assembly", { "scs" })
    PARAM(write_intermediates, Bool, false, "Write XYZ file after each fragment addition for debugging", "Output", { "intermediates" })
    PARAM(verbose, Bool, true, "Detailed output", "Output", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};
