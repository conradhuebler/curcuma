/*
 * <Unified molecular analysis for all structure types and formats.>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"
#include "src/tools/formats.h"
#include "src/tools/trajectory_writer.h"
#include "tda_engine.h"

#include "curcumamethod.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

/*! \brief Unified molecular analysis for all structure types and formats - Claude Generated
 *
 * Provides comprehensive molecular property analysis that works transparently with:
 * - Atomistic structures (XYZ, MOL2, SDF, PDB)
 * - Coarse-grained structures (VTF)
 * - Single structures or full trajectories
 *
 * Available properties:
 * - Basic: Mass, atoms count, fragments, formula
 * - Geometric: Gyration radius, center of mass, inertia constants
 * - Polymer/CG: End-to-end distance, Rout (COM to outermost)
 * - Topological: Persistent homology analysis (ripser)
 */
class UnifiedAnalysis : public CurcumaMethod
{
public:
    UnifiedAnalysis(const json& controller, bool silent);
    ~UnifiedAnalysis();

    /*! \brief Start unified molecular analysis */
    void start() override;

    /*! \brief Set input filename (any supported format) */
    inline void setFileName(const std::string& filename) { m_filename = filename; }

    /*! \brief Print help for analysis options */
    void printHelp() const override;

    /*! \brief Print comprehensive help for enhanced topological analysis (dMatrix alternative)
     *
     * Displays detailed documentation for the enhanced topological analysis features
     * that replace the legacy -dMatrix command line interface. This includes:
     * - Distance matrix export options (.dMat files)
     * - Persistence diagram generation (.PD files)
     * - Persistence image creation (.PI files)
     * - Advanced image formats and post-processing
     * - Batch trajectory analysis capabilities
     */
    void printEnhancedTDAHelp() const;

    // CurcumaMethod interface implementation
    nlohmann::json WriteRestartInformation() override { return json{}; }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return {"Analysis"}; }
    void ReadControlFile() override {}
    void LoadControlJson() override {}

private:
    /*! \brief Analyze single molecule structure */
    json analyzeMolecule(const Molecule& mol, int timestep = -1);

    /*! \brief Calculate basic molecular properties */
    json calculateBasicProperties(const Molecule& mol);

    /*! \brief Calculate geometric properties */
    json calculateGeometricProperties(const Molecule& mol);

    /*! \brief Calculate polymer/chain properties */
    json calculatePolymerProperties(const Molecule& mol);

    /*! \brief Calculate topological properties (ripser) */
    json calculateTopologicalProperties(const Molecule& mol);

    /*! \brief Detect if molecule has chain/polymer structure */
    bool detectChainStructure(const Molecule& mol);

    /*! \brief Calculate scattering properties P(q) and S(q) - Claude Generated 2025 */
    json calculateScatteringProperties(const Molecule& mol);

    /*! \brief Calculate radial distribution function g(r) - Claude Generated 2025 */
    json calculatePairDistribution(const Molecule& mol);

    /*! \brief Output results in requested format */
    void outputResults(const json& results);

    /*! \brief Output results to file */
    void outputToFile(const json& results, const std::string& filename);

    std::string m_filename;
    ConfigManager m_config;       // Claude Generated 2025: Modern configuration manager
    json m_config_legacy;         // Claude Generated 2025: Legacy JSON for TDAEngine compatibility
    bool m_silent;

    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    BEGIN_PARAMETER_DEFINITION(analysis)

    // Analysis options - Claude Generated
    PARAM(properties, String, "all", "Properties to calculate: all|basic|geometric|cg|topology", "Analysis", {})
    PARAM(output_format, String, "human", "Output format: human|json|csv", "Output", { "format" })
    PARAM(output_file, String, "", "Optional output file path", "Output", { "out" })
    PARAM(fragments, Bool, true, "Include per-fragment analysis", "Analysis", {})
    PARAM(ripser, Bool, false, "Include basic topological analysis", "Topology", {})
    PARAM(verbose, Bool, true, "Detailed output", "Output", {})

    // Trajectory statistics options - Claude Generated
    PARAM(metrics, String, "gyration,rout,end2end", "Comma-separated trajectory metrics: gyration|rout|end2end|com|inertia|mass|all", "Trajectory", {})
    PARAM(statistics, String, "none", "Statistics mode: none|cumulative|moving|all", "Trajectory", { "stats" })
    PARAM(window, Int, 10, "Moving average window size", "Trajectory", {})

    // Frame selection and stride options - Claude Generated 2026
    PARAM(frames, String, "", "Frame selection (e.g., '1:5,8,10:12', 'last' or '1:-1'=all, 'N:N'=single)", "Trajectory", {})
    PARAM(stride, Int, 1, "Analyze every N-th frame from selection", "Trajectory", {})
    PARAM(frame_range_start, Int, 0, "Alternative: Start frame (used if 'frames' empty)", "Trajectory", {"start_frame"})
    PARAM(frame_range_end, Int, -1, "Alternative: End frame, -1=last (used if 'frames' empty)", "Trajectory", {"end_frame"})

    // Enhanced topological analysis options - Claude Generated
    PARAM(topological_save_distance_matrix, Bool, false, "Save distance matrix (.dMat files)", "Topology", {})
    PARAM(topological_save_persistence_pairs, Bool, false, "Save persistence pairs (.pairs files)", "Topology", {})
    PARAM(topological_save_persistence_diagram, Bool, false, "Save persistence diagram (.PD files)", "Topology", {})
    PARAM(topological_save_persistence_image, Bool, false, "Save persistence images (.PI files)", "Topology", {})
    PARAM(topological_exclude_bonds, Bool, false, "Exclude bonds from distance matrix", "Topology", {})
    PARAM(topological_exclude_hydrogen, Bool, false, "Exclude hydrogen atoms from analysis", "Topology", {})
    PARAM(topological_print_elements, Bool, false, "Include element symbols in output", "Topology", {})
    PARAM(topological_print_energy, Bool, false, "Include molecular energies in output", "Topology", {})
    PARAM(topological_image_format, String, "png", "Image output format: png|jpg|bmp|tga", "Topology", {})
    PARAM(topological_colormap, String, "hot", "Image colormap: grayscale|jet|hot|viridis|coolwarm", "Topology", {})
    PARAM(topological_resolution, String, "800x800", "Image resolution (e.g., 800x800, 1024x1024)", "Topology", {})
    PARAM(topological_post_processing, String, "none", "Post-processing: none|adaptive|ring_focused", "Topology", {})
    PARAM(topological_temperature, Double, 2.0, "Enhancement temperature parameter", "Topology Advanced", {})
    PARAM(topological_damping, Double, 1.5, "Damping strength for processing", "Topology Advanced", {})
    PARAM(topological_preserve_structure, Bool, true, "Maintain original structure during processing", "Topology Advanced", {})
    PARAM(topological_atom_selection, String, "", "Atom indices for selective analysis (e.g., 1;5-10;15)", "Topology", {})

    // Scattering analysis options - Claude Generated 2025
    PARAM(scattering_enable, Bool, false, "Enable P(q)/S(q) scattering analysis", "Scattering", {})
    PARAM(scattering_q_min, Double, 0.01, "Minimum q value (Å⁻¹)", "Scattering", {"qmin"})
    PARAM(scattering_q_max, Double, 2.0, "Maximum q value (Å⁻¹)", "Scattering", {"qmax"})
    PARAM(scattering_q_steps, Int, 100, "Number of q points", "Scattering", {"qsteps"})
    PARAM(scattering_form_factor, String, "auto", "Form factor type: auto|cromer_mann|cg_sphere", "Scattering", {"ff"})
    PARAM(scattering_cg_radius, Double, 3.0, "CG bead radius for sphere form factor (Å)", "Scattering Advanced", {})
    PARAM(scattering_angular_samples, Int, 50, "Angular samples for S(q) spherical averaging", "Scattering Advanced", {})

    // Per-frame scattering file output options - Claude Generated 2026
    PARAM(scattering_per_frame_files, Bool, false, "Generate separate CSV files for each frame", "Scattering", {})
    PARAM(scattering_q_values, String, "", "Comma-separated q-values to export (e.g., 0.0,0.1,0.5,1.0,2.0)", "Scattering", {})
    PARAM(scattering_output_directory, String, ".", "Directory for per-frame files", "Scattering", {})
    PARAM(scattering_file_prefix, String, "scattering_frame", "Prefix for generated files", "Scattering", {})

    // Radial distribution function options - Claude Generated 2025
    PARAM(rdf_enable, Bool, false, "Enable g(r) radial distribution function calculation", "RDF", {})
    PARAM(rdf_r_max, Double, 15.0, "Maximum r for g(r) (Å)", "RDF", {"rmax"})
    PARAM(rdf_bin_width, Double, 0.05, "Histogram bin width (Å)", "RDF", {"dr"})
    PARAM(rdf_coordination_shells, Bool, false, "Calculate running coordination number n(r)", "RDF", {})

    // Shape descriptor options - Claude Generated 2025
    PARAM(shape_asphericity, Bool, false, "Calculate asphericity descriptor (b = λ₃ - 0.5(λ₁+λ₂))", "Shape", {})
    PARAM(shape_acylindricity, Bool, false, "Calculate acylindricity descriptor (c = λ₂ - λ₁)", "Shape", {})
    PARAM(shape_anisotropy, Bool, false, "Calculate relative shape anisotropy (κ²)", "Shape", {})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
};