/*
 * <Topological Data Analysis Engine Implementation>
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

#include "tda_engine.h"
#include "src/core/imagewriter.hpp"
#include "src/core/curcuma_logger.h"

#include <fstream>
#include <iomanip>
#include <sstream>

TDAEngine::TDAEngine(const json& config)
    : m_config(config)
    , m_output_prefix("")
{
    // Parse image parameters from config
    m_image_format = m_config.value("image_format", "png");
    m_colormap = m_config.value("colormap", "hot");
    m_post_processing = m_config.value("post_processing", "none");
    m_temperature = m_config.value("temperature", 2.0);
    m_damping = m_config.value("damping", 1.5);
    m_preserve_structure = m_config.value("preserve_structure", true);

    // Parse resolution
    std::string resolution = m_config.value("resolution", "800x800");
    auto [width, height] = parseResolution(resolution);
    m_width = width;
    m_height = height;

    if (m_config.contains("output_prefix")) {
        m_output_prefix = m_config["output_prefix"].get<std::string>();
    }
}

void TDAEngine::updateConfig(const json& new_config)
{
    // Merge new config with existing
    for (auto it = new_config.begin(); it != new_config.end(); ++it) {
        m_config[it.key()] = it.value();
    }

    // Re-parse parameters
    m_image_format = m_config.value("image_format", "png");
    m_colormap = m_config.value("colormap", "hot");
    auto [width, height] = parseResolution(m_config.value("resolution", "800x800"));
    m_width = width;
    m_height = height;
}

json TDAEngine::analyzeMolecule(const Molecule& mol, int index)
{
    json result;

    try {
        // Claude Generated: Input validation to prevent crashes
        if (mol.AtomCount() == 0) {
            result["error"] = "Empty molecule - no atoms to analyze";
            return result;
        }

        if (mol.AtomCount() == 1) {
            result["warning"] = "Single atom molecule - limited topological analysis possible";
        }

        // Basic molecule info
        result["molecule_index"] = index;
        result["atom_count"] = mol.AtomCount();
        result["energy"] = mol.Energy();

        // Distance matrix analysis
        if (m_config.value("save_distance_matrix", false) ||
            m_config.value("save_persistence_pairs", false) ||
            m_config.value("save_persistence_diagram", false) ||
            m_config.value("save_persistence_image", false)) {

            result["distance_matrix"] = processDistanceMatrix(mol, index);
        }

        // Persistent homology analysis - Claude Generated: Only run if specifically requested
        if (m_config.value("save_persistence_pairs", false) ||
            m_config.value("save_persistence_diagram", false)) {
            result["persistence_diagram"] = processPersistenceDiagram(mol, index);
        }

        if (m_config.value("save_persistence_image", false)) {
            result["persistence_image"] = processPersistenceImage(mol, index);
        }

        // Basic topological invariants
        result["topology_summary"] = {
            {"connectivity_components", 1}, // Simplified - would need complex calculation
            {"cycles", 0}                   // Simplified - would need complex calculation
        };

    } catch (const std::exception& e) {
        result["error"] = e.what();
        CurcumaLogger::error("TDA analysis failed: " + std::string(e.what()));
    }

    return result;
}

json TDAEngine::analyzeTrajectory(const std::vector<Molecule>& molecules)
{
    json result;
    std::vector<Eigen::MatrixXd> pd_matrices;
    std::vector<Eigen::MatrixXd> pi_matrices;

    result["trajectory_length"] = molecules.size();
    result["frames"] = json::array();

    // Process each frame
    for (size_t i = 0; i < molecules.size(); ++i) {
        json frame_result = analyzeMolecule(molecules[i], static_cast<int>(i));
        result["frames"].push_back(frame_result);

        // Collect matrices for statistical analysis
        if (frame_result.contains("persistence_diagram") &&
            frame_result["persistence_diagram"].contains("matrix_norm")) {
            // Would need to store actual matrices for statistics
            // This is a simplified version
        }
    }

    // Calculate trajectory-wide statistics
    if (molecules.size() > 1) {
        result["trajectory_statistics"] = calculateTrajectoryStatistics(pd_matrices, pi_matrices);
    }

    return result;
}

json TDAEngine::processDistanceMatrix(const Molecule& mol, int index)
{
    json result;

    // Claude Generated: Validate distance matrix data
    auto distance_vector = mol.LowerDistanceVector(m_config.value("exclude_hydrogen", false));
    if (distance_vector.empty()) {
        result["error"] = "Empty distance matrix - cannot compute topological properties";
        return result;
    }

    bool exclude_bonds = m_config.value("exclude_bonds", false);
    bool print_elements = m_config.value("print_elements", false);

    // Save distance matrix if requested - Claude Generated: Include legacy functionality
    if (m_config.value("save_distance_matrix", false)) {
        if (saveDistanceMatrix(mol, index)) {
            result["distance_matrix_file"] = getOutputFilename("dMat", index);
        }
    }

    // Calculate basic distance matrix statistics using already validated distance_vector
    result["distance_count"] = distance_vector.size();

    if (!distance_vector.empty()) {
        double min_dist = *std::min_element(distance_vector.begin(), distance_vector.end());
        double max_dist = *std::max_element(distance_vector.begin(), distance_vector.end());
        double avg_dist = std::accumulate(distance_vector.begin(), distance_vector.end(), 0.0) / distance_vector.size();

        result["distance_statistics"] = {
            {"min_distance", min_dist},
            {"max_distance", max_dist},
            {"average_distance", avg_dist}
        };
    }

    return result;
}

json TDAEngine::processPersistenceDiagram(const Molecule& mol, int index)
{
    json result;

    try {
        // Claude Generated: Validate input data for persistence diagram
        auto distance_vector = mol.LowerDistanceVector(m_config.value("exclude_hydrogen", false));
        if (distance_vector.empty()) {
            result["error"] = "Empty distance matrix - cannot compute persistence diagram";
            return result;
        }

        auto en_scaling = mol.DeltaEN();
        if (en_scaling.empty()) {
            result["error"] = "Empty EN scaling data - cannot compute persistence pairs";
            return result;
        }

        // Create persistence diagram with current config
        PersistentDiagram diagram(m_config);
        diagram.setDistanceMatrix(distance_vector);
        diagram.setENScaling(en_scaling);

        // Generate persistence pairs
        auto pairs = diagram.generatePairs();
        result["pairs_count"] = pairs.size();

        // Claude Generated: Check if pairs generation succeeded
        if (pairs.empty()) {
            result["warning"] = "No persistence pairs generated - topology may be trivial";
        }

        // Save pairs if requested
        if (m_config.value("save_persistence_pairs", false)) {
            if (savePersistencePairs(pairs, index)) {
                result["pairs_file"] = getOutputFilename("pairs", index);
            }
        }

        // Generate persistence diagram image
        Eigen::MatrixXd pd_matrix = diagram.generateImage(pairs);

        // Claude Generated: Validate generated matrix
        if (pd_matrix.size() == 0) {
            result["warning"] = "Empty persistence diagram matrix generated";
            result["matrix_norm"] = 0.0;
            result["matrix_size"] = {0, 0};
        } else {
            result["matrix_norm"] = pd_matrix.norm();
            result["matrix_size"] = {pd_matrix.rows(), pd_matrix.cols()};
        }

        // Save persistence diagram as text
        if (m_config.value("save_persistence_diagram", false)) {
            if (savePersistenceDiagramText(pd_matrix, index)) {
                result["diagram_text_file"] = getOutputFilename("PD", index);
            }
        }

        // Save persistence diagram as image
        if (m_config.value("save_persistence_image", false)) {
            if (savePersistenceDiagramImage(pd_matrix, index)) {
                result["diagram_image_file"] = getOutputFilename("PD", index, m_image_format);
            }
        }

        // Calculate basic persistence statistics
        if (!pairs.empty()) {
            std::vector<double> lifetimes;
            for (const auto& pair : pairs) {
                lifetimes.push_back(pair.second - pair.first);
            }

            double max_lifetime = *std::max_element(lifetimes.begin(), lifetimes.end());
            double avg_lifetime = std::accumulate(lifetimes.begin(), lifetimes.end(), 0.0) / lifetimes.size();

            result["persistence_statistics"] = {
                {"max_lifetime", max_lifetime},
                {"average_lifetime", avg_lifetime},
                {"persistent_features", pairs.size()}
            };
        }

    } catch (const std::exception& e) {
        result["error"] = e.what();
    }

    return result;
}

json TDAEngine::processPersistenceImage(const Molecule& mol, int index)
{
    json result;

    try {
        // Claude Generated: Validate input data for persistence image
        auto distance_vector = mol.LowerDistanceVector(m_config.value("exclude_hydrogen", false));
        if (distance_vector.empty()) {
            result["error"] = "Empty distance matrix - cannot compute persistence image";
            return result;
        }

        auto en_scaling = mol.DeltaEN();
        if (en_scaling.empty()) {
            result["error"] = "Empty EN scaling data - cannot compute persistence image";
            return result;
        }

        // Create fresh persistence diagram for PI generation - Claude Generated
        PersistentDiagram diagram(m_config);
        // Must reset distance matrix as it may have been moved by previous operations
        auto fresh_distance_vector = mol.LowerDistanceVector(m_config.value("exclude_hydrogen", false));
        diagram.setDistanceMatrix(fresh_distance_vector);
        diagram.setENScaling(en_scaling);  // Claude Generated: Add missing EN scaling for generateTriples()

        // Generate persistence image (EN scaled bond topology)
        auto triples = diagram.generateTriples();
        result["triples_count"] = triples.size();

        // Claude Generated: Check if triples generation succeeded
        if (triples.empty()) {
            result["warning"] = "No persistence triples generated - topology may be trivial";
            result["matrix_norm"] = 0.0;
            result["matrix_size"] = {0, 0};
        } else {
            Eigen::MatrixXd pi_matrix = diagram.generateImage(triples);

            // Claude Generated: Validate generated matrix
            if (pi_matrix.size() == 0) {
                result["warning"] = "Empty persistence image matrix generated";
                result["matrix_norm"] = 0.0;
                result["matrix_size"] = {0, 0};
            } else {
                result["matrix_norm"] = pi_matrix.norm();
                result["matrix_size"] = {pi_matrix.rows(), pi_matrix.cols()};

                // Save persistence image files only if matrix is valid
                if (m_config.value("save_persistence_image", false)) {
                    if (savePersistenceImageText(pi_matrix, index)) {
                        result["image_text_file"] = getOutputFilename("PI", index);
                    }

                    if (savePersistenceImageImage(pi_matrix, index)) {
                        result["image_file"] = getOutputFilename("PI", index, m_image_format);
                    }
                }
            }
        }

    } catch (const std::exception& e) {
        result["error"] = e.what();
    }

    return result;
}

bool TDAEngine::saveDistanceMatrix(const Molecule& mol, int index)
{
    try {
        std::string filename = getOutputFilename("dMat", index);
        std::ofstream file(filename);

        if (m_config.value("print_energy", false)) {
            file << std::setprecision(10) << mol.Energy() << std::endl;
        }

        bool exclude_bonds = m_config.value("exclude_bonds", false);
        bool print_elements = m_config.value("print_elements", false);
        file << mol.DistanceMatrixString(exclude_bonds, print_elements);

        file.close();
        return true;
    } catch (const std::exception& e) {
        CurcumaLogger::error("Failed to save distance matrix: " + std::string(e.what()));
        return false;
    }
}

bool TDAEngine::savePersistencePairs(const dpairs& pairs, int index)
{
    try {
        std::string filename = getOutputFilename("pairs", index);
        std::ofstream file(filename);

        for (const auto& pair : pairs) {
            file << pair.first << " " << pair.second << std::endl;
        }

        file.close();
        return true;
    } catch (const std::exception& e) {
        CurcumaLogger::error("Failed to save persistence pairs: " + std::string(e.what()));
        return false;
    }
}

bool TDAEngine::savePersistenceDiagramText(const Eigen::MatrixXd& pd_matrix, int index)
{
    try {
        std::string filename = getOutputFilename("PD", index);
        std::ofstream file(filename);
        file << pd_matrix;
        file.close();
        return true;
    } catch (const std::exception& e) {
        CurcumaLogger::error("Failed to save persistence diagram text: " + std::string(e.what()));
        return false;
    }
}

bool TDAEngine::savePersistenceDiagramImage(const Eigen::MatrixXd& pd_matrix, int index)
{
    try {
        std::string filename = getOutputFilename("PD", index, m_image_format);
        EigenImageWriter::saveMatrix(pd_matrix, filename, parseColormap(m_colormap), 90, false,
                                    m_width, m_height, parsePostProcessing( m_post_processing ),
                                    m_temperature, m_damping, m_preserve_structure);
        return true;
    } catch (const std::exception& e) {
        CurcumaLogger::error("Failed to save persistence diagram image: " + std::string(e.what()));
        return false;
    }
}

bool TDAEngine::savePersistenceImageText(const Eigen::MatrixXd& pi_matrix, int index)
{
    try {
        std::string filename = getOutputFilename("PI", index);
        std::ofstream file(filename);
        file << pi_matrix;
        file.close();
        return true;
    } catch (const std::exception& e) {
        CurcumaLogger::error("Failed to save persistence image text: " + std::string(e.what()));
        return false;
    }
}

bool TDAEngine::savePersistenceImageImage(const Eigen::MatrixXd& pi_matrix, int index)
{
    try {
        std::string filename = getOutputFilename("PI", index, m_image_format);
        EigenImageWriter::saveMatrix(pi_matrix, filename, parseColormap(m_colormap), 90, false,
                                    m_width, m_height, parsePostProcessing(m_post_processing),
                                    m_temperature, m_damping, m_preserve_structure);
        return true;
    } catch (const std::exception& e) {
        CurcumaLogger::error("Failed to save persistence image: " + std::string(e.what()));
        return false;
    }
}

json TDAEngine::calculateTrajectoryStatistics(const std::vector<Eigen::MatrixXd>& pd_matrices,
                                             const std::vector<Eigen::MatrixXd>& pi_matrices)
{
    json result;

    // This would implement the statistical analysis from dMatrix
    // For now, return basic info
    result["frames_processed"] = pd_matrices.size();
    result["note"] = "Trajectory statistics calculation - implementation pending";

    return result;
}

std::pair<int, int> TDAEngine::parseResolution(const std::string& resolution)
{
    size_t pos = resolution.find('x');
    if (pos == std::string::npos) {
        return {800, 800}; // Default
    }

    try {
        int width = std::stoi(resolution.substr(0, pos));
        int height = std::stoi(resolution.substr(pos + 1));
        return {width, height};
    } catch (...) {
        return {800, 800}; // Default on error
    }
}

std::string TDAEngine::getOutputFilename(const std::string& type, int index, const std::string& suffix)
{
    std::stringstream ss;
    ss << m_output_prefix;

    if (index >= 0) {
        ss << "_" << index;
    }

    ss << "." << type;

    if (!suffix.empty()) {
        ss << "." << suffix;
    }

    return ss.str();
}

EigenImageWriter::ColorMap TDAEngine::parseColormap(const std::string& cmap)
{
    if (cmap == "jet") return EigenImageWriter::JET;
    if (cmap == "hot") return EigenImageWriter::HOT;
    if (cmap == "viridis") return EigenImageWriter::VIRIDIS;
    if (cmap == "coolwarm") return EigenImageWriter::COOLWARM;
    return EigenImageWriter::GRAYSCALE; // Default
}

EigenImageWriter::PostProcessing TDAEngine::parsePostProcessing(const std::string& pp)
{    
    //enum PostProcessing { NONE,
    //    ADAPTIVE,
    //    RING_FOCUSED };

    if (pp == "adaptive") return EigenImageWriter::ADAPTIVE;
    if (pp == "ring_focused") return EigenImageWriter::RING_FOCUSED;   
    return EigenImageWriter::NONE; // Default
}