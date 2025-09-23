/*
 * <Topological Data Analysis Engine - Integration of dMatrix functionality>
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

#include <string>
#include <vector>
#include <Eigen/Dense>

#include "src/core/molecule.h"
#include "persistentdiagram.h"
#include "json.hpp"
#include "src/core/imagewriter.hpp"

using json = nlohmann::json;

// Configuration for TDA Engine - Claude Generated
static const json TDAEngineJson = {
    {"save_distance_matrix", false},
    {"save_persistence_pairs", false},
    {"save_persistence_diagram", false},
    {"save_persistence_image", false},
    {"exclude_bonds", false},
    {"exclude_hydrogen", false},
    {"print_elements", false},
    {"print_energy", false},
    {"image_format", "png"},
    {"colormap", "hot"},
    {"resolution", "800x800"},
    {"post_processing", "none"},
    {"temperature", 2.0},
    {"damping", 1.5},
    {"preserve_structure", true},
    {"output_prefix", ""}
};

/*! \brief Comprehensive Topological Data Analysis Engine - Claude Generated
 *
 * **MODERN REPLACEMENT FOR -dMatrix FUNCTIONALITY**
 *
 * This engine provides research-grade topological data analysis functionality
 * equivalent to the legacy -dMatrix command but with modern integration into
 * the unified analysis system. All original dMatrix features are preserved
 * and enhanced with additional capabilities.
 *
 * **MIGRATION FROM -dMatrix:**
 * Old: curcuma -dMatrix molecule.xyz
 * New: curcuma -analysis molecule.xyz -topological.save_persistence_image true
 *
 * **CORE CAPABILITIES:**
 * - Distance matrix calculation and export (.dMat files)
 * - Persistent homology analysis with persistence diagrams (.PD files)
 * - Persistence image generation (.PI files)
 * - Advanced image processing with multiple colormaps and post-processing
 * - Batch trajectory analysis with statistical processing
 * - Multi-format image export (PNG, JPG, BMP, TGA)
 * - Configurable filtering (hydrogen exclusion, bond exclusion)
 *
 * **RESEARCH APPLICATIONS:**
 * - Molecular shape analysis and comparison
 * - Protein folding and conformational studies
 * - Materials science and crystal structure analysis
 * - Drug discovery and molecular recognition
 * - Machine learning feature generation
 *
 * **INTEGRATION ARCHITECTURE:**
 * This engine encapsulates all dMatrix functionality and integrates seamlessly
 * with the unified analysis system while maintaining research-grade quality
 * suitable for publication in computational chemistry and materials science.
 *
 * **CITATION:**
 * When using TDA features, please cite:
 * Townsend, J., Micucci, C.P., Hymel, J.H. et al. Nat Commun 11, 3230 (2020)
 * https://doi.org/10.1038/s41467-020-17035-5
 */
class TDAEngine
{
public:
    /*! \brief Constructor with configuration */
    explicit TDAEngine(const json& config = TDAEngineJson);

    /*! \brief Destructor */
    ~TDAEngine() = default;

    /*! \brief Perform complete TDA analysis on a single molecule */
    json analyzeMolecule(const Molecule& mol, int index = 0);

    /*! \brief Process multiple molecules/trajectory with statistical analysis */
    json analyzeTrajectory(const std::vector<Molecule>& molecules);

    /*! \brief Set output file prefix for exports */
    void setOutputPrefix(const std::string& prefix) { m_output_prefix = prefix; }

    /*! \brief Get current configuration */
    const json& getConfig() const { return m_config; }

    /*! \brief Update configuration */
    void updateConfig(const json& new_config);

private:
    /*! \brief Calculate and optionally save distance matrix */
    json processDistanceMatrix(const Molecule& mol, int index);

    /*! \brief Convert colormap string to enum */
    EigenImageWriter::ColorMap parseColormap(const std::string& cmap);

    /*! \brief Convert postprocessing string to enum */
    EigenImageWriter::PostProcessing parsePostProcessing(const std::string& pp);

    /*! \brief Generate persistence diagram and handle exports */
    json processPersistenceDiagram(const Molecule& mol, int index);

    /*! \brief Generate persistence image and handle exports */
    json processPersistenceImage(const Molecule& mol, int index);

    /*! \brief Save distance matrix to file */
    bool saveDistanceMatrix(const Molecule& mol, int index);

    /*! \brief Save persistence pairs to file */
    bool savePersistencePairs(const dpairs& pairs, int index);

    /*! \brief Save persistence diagram as text */
    bool savePersistenceDiagramText(const Eigen::MatrixXd& pd_matrix, int index);

    /*! \brief Save persistence diagram as image */
    bool savePersistenceDiagramImage(const Eigen::MatrixXd& pd_matrix, int index);

    /*! \brief Save persistence image as text */
    bool savePersistenceImageText(const Eigen::MatrixXd& pi_matrix, int index);

    /*! \brief Save persistence image as image */
    bool savePersistenceImageImage(const Eigen::MatrixXd& pi_matrix, int index);

    /*! \brief Calculate statistical analysis for trajectory */
    json calculateTrajectoryStatistics(const std::vector<Eigen::MatrixXd>& pd_matrices,
                                      const std::vector<Eigen::MatrixXd>& pi_matrices);

    /*! \brief Save trajectory statistics (average, stddev images) */
    void saveTrajectoryStatistics(const Eigen::MatrixXd& avg_pd, const Eigen::MatrixXd& std_pd,
                                 const Eigen::MatrixXd& avg_pi, const Eigen::MatrixXd& std_pi);

    /*! \brief Parse resolution string (e.g., "800x600") */
    std::pair<int, int> parseResolution(const std::string& resolution);

    /*! \brief Get full output filename for given type and index */
    std::string getOutputFilename(const std::string& type, int index = -1, const std::string& suffix = "");

    // Configuration and state
    json m_config;
    std::string m_output_prefix;

    // Image export parameters (parsed from config)
    std::string m_image_format;
    std::string m_colormap;
    int m_width, m_height;
    std::string m_post_processing;
    double m_temperature;
    double m_damping;
    bool m_preserve_structure;
};