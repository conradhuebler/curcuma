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
#include "tda_engine.h"

#include "curcumamethod.h"

// Default configuration for unified analysis - Claude Generated
static const json UnifiedAnalysisJson = {
    { "properties", "all" },           // Which properties to calculate: "all", "basic", "cg", "topology"
    { "output_format", "human" },      // Output format: "human", "json", "csv"
    { "output_file", "" },             // Optional output file
    { "trajectory", false },           // Analyze full trajectory vs single structure
    { "fragments", true },             // Include per-fragment analysis
    { "ripser", false },               // Include basic topological analysis
    { "verbose", true },               // Detailed output
    // Enhanced topological analysis options (dMatrix integration)
    { "topological", {
        { "save_distance_matrix", false },
        { "save_persistence_pairs", false },
        { "save_persistence_diagram", false },
        { "save_persistence_image", false },
        { "exclude_bonds", false },
        { "exclude_hydrogen", false },
        { "print_elements", false },
        { "print_energy", false },
        { "image_format", "png" },
        { "colormap", "hot" },
        { "resolution", "800x800" },
        { "post_processing", "none" },
        { "temperature", 2.0 },
        { "damping", 1.5 },
        { "preserve_structure", true },
        { "atom_selection", "" }     // New: Atom indices for selective analysis (e.g., "1,5:10,15")
    }}
};

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

    /*! \brief Output results in requested format */
    void outputResults(const json& results);

    /*! \brief Output results to file */
    void outputToFile(const json& results, const std::string& filename);

    std::string m_filename;
    json m_config;
    bool m_silent;
};