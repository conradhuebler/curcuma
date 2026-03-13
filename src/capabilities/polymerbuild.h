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
 * @brief PolymerBuild Capability
 *
 * Assembles complex polymers from fragment structures using a sequence string.
 * Connection points are defined by dummy atoms (Xx) in the fragments.
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
     * @brief Connect a new fragment to the current molecule
     * @param current The existing molecule (will be updated)
     * @param fragment_name Name of the fragment to add
     * @param fragment_file Path to the fragment XYZ file
     * @param interface_bonds List of bonds between atom indices (will be updated)
     */
    void connectFragment(Molecule& current, const std::string& fragment_name,
                         const std::string& fragment_file,
                         std::vector<std::pair<int, int>>& interface_bonds);

    /**
     * @brief Apply capping to the terminal connection points
     * @param mol The molecule to cap
     */
    void applyCapping(Molecule& mol);

    ConfigManager m_config;
    bool m_silent;
    std::map<std::string, std::string> m_fragments;

    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    BEGIN_PARAMETER_DEFINITION(polymerbuild)

    // Sequence options
    PARAM(sequence, String, "", "Polymer sequence string (e.g., (A)10-B-(C)5)", "Assembly", { "seq" })
    PARAM(fragments, String, "{}", "JSON map of fragment names to file paths", "Assembly", { "frag" })
    PARAM(cap_start, String, "H", "Functional group for the start of the chain (H, CH3, or none)", "Assembly", {})
    PARAM(cap_end, String, "H", "Functional group for the end of the chain (H, CH3, or none)", "Assembly", {})

    // Optimization options
    PARAM(optimize, Bool, true, "Enable intermediate optimization steps", "Refinement", { "opt" })
    PARAM(opt_method, String, "uff", "Method for optimization (uff, gfn2, etc.)", "Refinement", {})

    // Dynamics options
    PARAM(dynamics, Bool, false, "Enable intermediate molecular dynamics", "Refinement", { "md" })
    PARAM(md_steps, Int, 1000, "Number of MD steps for intermediate dynamics", "Refinement", {})
    PARAM(md_method, String, "uff", "Method for dynamics", "Refinement", {})

    // Geometric options
    PARAM(bond_distance_scaling, Double, 1.0, "Scaling factor for interface bond lengths", "Assembly", { "scaling" })
    PARAM(verbose, Bool, true, "Detailed output", "Output", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};
