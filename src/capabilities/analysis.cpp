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

#include "analysis.h"
#include "persistentdiagram.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/units.h"
#include "src/core/parameter_registry.h"
#include "src/tools/pbc_utils.h"
#include "trajectory_statistics.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

UnifiedAnalysis::UnifiedAnalysis(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("analysis"), controller, silent)
    , m_config("analysis", controller)  // Claude Generated 2025: ConfigManager accepts full controller
    , m_silent(silent)
{
    // Claude Generated 2025: Analysis needs verbosity >=2 for output
    // If user didn't explicitly set verbosity, default to 2 (informative)
    if (!controller.contains("verbosity") && !controller.contains("verbose") && !silent) {
        m_verbosity = 2;
        CurcumaLogger::set_verbosity(m_verbosity);
    }

    UpdateController(controller);

    // Claude Generated 2025: Legacy JSON export for TDAEngine compatibility
    m_config_legacy = m_config.exportConfig();
}

UnifiedAnalysis::~UnifiedAnalysis()
{
}

// Helper function: Parse comma-separated metrics list - Claude Generated 2025
std::vector<std::string> parseMetricsList(const std::string& metrics_str)
{
    std::vector<std::string> metrics;

    if (metrics_str == "all") {
        return { "gyration", "rout", "end2end", "com", "mass" };
    }

    std::stringstream ss(metrics_str);
    std::string metric;
    while (std::getline(ss, metric, ',')) {
        // Trim whitespace
        metric.erase(0, metric.find_first_not_of(" \t"));
        metric.erase(metric.find_last_not_of(" \t") + 1);
        if (!metric.empty()) {
            metrics.push_back(metric);
        }
    }
    return metrics;
}

// Helper function: Extract metric value from JSON - Claude Generated 2025
double extractMetricValue(const json& timestep_result, const std::string& metric, bool mass_weighted = false)
{
    if (metric == "gyration") {
        if (timestep_result.contains("geometric") && timestep_result["geometric"].contains("gyration_radius")) {
            const auto& gyr = timestep_result["geometric"]["gyration_radius"];
            return mass_weighted ? gyr["mass_weighted"].get<double>() : gyr["unweighted"].get<double>();
        }
    } else if (metric == "rout") {
        if (timestep_result.contains("polymer") && timestep_result["polymer"].contains("rout")) {
            return timestep_result["polymer"]["rout"].get<double>();
        }
    } else if (metric == "end2end") {
        if (timestep_result.contains("polymer") && timestep_result["polymer"].contains("end_to_end_distance")) {
            return timestep_result["polymer"]["end_to_end_distance"].get<double>();
        }
    } else if (metric == "mass") {
        if (timestep_result.contains("basic") && timestep_result["basic"].contains("mass")) {
            return timestep_result["basic"]["mass"].get<double>();
        }
    } else if (metric == "com") {
        // Return COM magnitude (distance from origin)
        if (timestep_result.contains("geometric") && timestep_result["geometric"].contains("center_of_mass")) {
            const auto& com = timestep_result["geometric"]["center_of_mass"];
            double x = com[0].get<double>();
            double y = com[1].get<double>();
            double z = com[2].get<double>();
            return std::sqrt(x * x + y * y + z * z);
        }
    }
    return 0.0;
}

void UnifiedAnalysis::start()
{
    if (m_filename.empty()) {
        CurcumaLogger::error("No input file specified for analysis");
        return;
    }

    if (!m_silent) {
        CurcumaLogger::info_fmt("Starting unified molecular analysis of: {}", m_filename);
    }

    // Claude Generated 2025: Unified approach - ALWAYS use FileIterator
    // FileIterator automatically stops after 1 structure for single-structure files
    // No need for separate trajectory vs single-structure code paths
    FileIterator file_iter(m_filename, m_silent);

    json results;
    results["filename"] = m_filename;
    // Note: MaxMolecules() may return 0 for some formats (XYZ), but that's OK
    // We'll count structures as we iterate instead
    results["total_timesteps"] = 0; // Will be updated after loop
    results["timesteps"] = json::array();

    // Parse statistics configuration - Claude Generated 2025
    std::string metrics_str = m_config.get<std::string>("metrics");
    std::string statistics_mode = m_config.get<std::string>("statistics");
    int window_size = m_config.get<int>("window");

    std::vector<std::string> enabled_metrics = parseMetricsList(metrics_str);
    bool enable_cumulative = (statistics_mode == "cumulative" || statistics_mode == "all");
    bool enable_moving = (statistics_mode == "moving" || statistics_mode == "all");

    // Create statistics tracker if needed
    std::unique_ptr<TrajectoryStatistics> stats;
    if (enable_cumulative || enable_moving) {
        stats = std::make_unique<TrajectoryStatistics>(window_size);
    }

    // ALWAYS store configuration in results for output formatting
    results["statistics_config"] = {
        { "metrics", enabled_metrics },
        { "mode", statistics_mode },
        { "window", window_size },
        { "enable_cumulative", enable_cumulative },
        { "enable_moving", enable_moving }
    };

    // Process all structures (1 for single-structure, N for trajectory)
    int timestep = 0;

    while (!file_iter.AtEnd()) {
        Molecule mol = file_iter.Next();

        // Claude Generated 2025: Log PBC detection at verbosity level 2+ (first structure only)
        if (timestep == 0 && mol.hasPBC() && CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Detected periodic boundary conditions");
            auto params = PBCUtils::getLatticeParameters(mol.getUnitCell());
            CurcumaLogger::param("lattice_a", fmt::format("{:.4f} Å", params[0]));
            CurcumaLogger::param("lattice_b", fmt::format("{:.4f} Å", params[1]));
            CurcumaLogger::param("lattice_c", fmt::format("{:.4f} Å", params[2]));
            CurcumaLogger::param("lattice_alpha", fmt::format("{:.2f}°", params[3]));
            CurcumaLogger::param("lattice_beta", fmt::format("{:.2f}°", params[4]));
            CurcumaLogger::param("lattice_gamma", fmt::format("{:.2f}°", params[5]));
        }

        json timestep_result = analyzeMolecule(mol, timestep);

        // Collect statistics for enabled metrics - Claude Generated 2025
        if (stats) {
            for (const auto& metric : enabled_metrics) {
                if (metric == "gyration") {
                    // Track both unweighted and mass-weighted
                    double gyr_u = extractMetricValue(timestep_result, "gyration", false);
                    double gyr_m = extractMetricValue(timestep_result, "gyration", true);
                    stats->addValue("gyration_unweighted", gyr_u);
                    stats->addValue("gyration_mass", gyr_m);
                } else {
                    double value = extractMetricValue(timestep_result, metric);
                    stats->addValue(metric, value);
                }
            }

            // Add statistics to timestep result
            timestep_result["statistics"] = json::object();
            for (const auto& metric : enabled_metrics) {
                if (metric == "gyration") {
                    // Both variants
                    for (const auto& variant : { "gyration_unweighted", "gyration_mass" }) {
                        if (enable_cumulative) {
                            timestep_result["statistics"][variant]["mean"] = stats->getMean(variant);
                            timestep_result["statistics"][variant]["std"] = stats->getStdDev(variant);
                        }
                        if (enable_moving) {
                            timestep_result["statistics"][variant]["moving_avg"] = stats->getMovingAverage(variant);
                        }
                    }
                } else {
                    if (enable_cumulative) {
                        timestep_result["statistics"][metric]["mean"] = stats->getMean(metric);
                        timestep_result["statistics"][metric]["std"] = stats->getStdDev(metric);
                    }
                    if (enable_moving) {
                        timestep_result["statistics"][metric]["moving_avg"] = stats->getMovingAverage(metric);
                    }
                }
            }
        }

        results["timesteps"].push_back(timestep_result);

        if (!m_silent && timestep % 10 == 0 && timestep > 0) {
            // Only show progress for multi-structure files
            int total = file_iter.MaxMolecules();
            if (total > 1) {
                CurcumaLogger::info_fmt("Analyzed timestep {}/{}", timestep, total);
            }
        }
        timestep++;
    }

    // Update total count after iterating - Claude Generated 2025
    results["total_timesteps"] = timestep;

    if (!m_silent && timestep > 1) {
        CurcumaLogger::info_fmt("Processed {} structures with {} metrics selected",
            timestep, enabled_metrics.size());
    }

    outputResults(results);
}

json UnifiedAnalysis::analyzeMolecule(const Molecule& mol, int timestep)
{
    json result;

    if (timestep >= 0) {
        result["timestep"] = timestep;
    }

    // Claude Generated 2025: Use metrics parameter for selective computation
    // Check if metrics parameter exists, otherwise fall back to properties for backward compatibility
    std::string metrics_str = m_config.get<std::string>("metrics");
    std::string properties = m_config.get<std::string>("properties");

    // Determine what to compute based on metrics (preferred) or properties (fallback)
    // Claude Generated 2025: ALWAYS compute basic properties (mass, atom count, etc.)
    bool compute_basic = true; // Basic properties are always needed
    bool compute_geometric = false;
    bool compute_polymer = false;
    bool compute_topology = false;

    if (!metrics_str.empty() && metrics_str != "none") {
        // Metrics-based computation (new way) - Claude Generated 2025
        std::vector<std::string> metrics = parseMetricsList(metrics_str);

        for (const auto& metric : metrics) {
            if (metric == "gyration" || metric == "com" || metric == "inertia") {
                compute_geometric = true;
            } else if (metric == "rout" || metric == "end2end") {
                compute_polymer = true;
            }
            // Note: "mass" is always in basic, so no need to check for it
        }

        // If "all" was specified, compute everything
        if (metrics_str == "all") {
            compute_geometric = compute_polymer = true;
        }
    } else {
        // Properties-based computation (old way for backward compatibility)
        compute_geometric = (properties == "all" || properties == "geometric");
        compute_polymer = (properties == "all" || properties == "cg" || properties == "polymer");
    }

    // Topology computation is independent of metrics/properties - Claude Generated 2025
    // Check both properties setting AND if any topological options are enabled
    compute_topology = (properties == "all" || properties == "topology");

    // Calculate only the needed properties - Claude Generated 2025
    if (compute_basic) {
        result["basic"] = calculateBasicProperties(mol);
    }

    if (compute_geometric) {
        result["geometric"] = calculateGeometricProperties(mol);
    }

    if (compute_polymer) {
        result["polymer"] = calculatePolymerProperties(mol);
    }

    // Calculate topological properties if requested
    if (compute_topology) {
        bool enable_ripser = m_config.get<bool>("ripser");

        // Also enable if any enhanced topological options are set - Claude Generated 2025
        bool enhanced_tda = (m_config.get<bool>("topological.save_distance_matrix") ||
                            m_config.get<bool>("topological.save_persistence_pairs") ||
                            m_config.get<bool>("topological.save_persistence_diagram") ||
                            m_config.get<bool>("topological.save_persistence_image"));


        if (enable_ripser || enhanced_tda) {
            result["topology"] = calculateTopologicalProperties(mol);
        }
    }

    return result;
}

json UnifiedAnalysis::calculateBasicProperties(const Molecule& mol)
{
    json basic;

    basic["atom_count"] = mol.AtomCount();
    basic["mass"] = mol.Mass();
    basic["charge"] = mol.Charge();
    basic["fragment_count"] = mol.GetFragments().size();

    if (!mol.Formula().empty()) {
        basic["formula"] = mol.Formula();
    }

    // Check if this might be a CG structure
    bool has_cg_elements = false;
    for (int atom : mol.Atoms()) {
        if (atom == CG_ELEMENT) {
            has_cg_elements = true;
            break;
        }
    }
    basic["coarse_grained"] = has_cg_elements;

    return basic;
}

json UnifiedAnalysis::calculateGeometricProperties(const Molecule& mol)
{
    json geometric;

    // Gyration radius (both mass-weighted and unweighted)
    // Claude Generated 2025: Use PBC-aware calculation when periodic boundaries are present
    auto gyr = mol.hasPBC()
        ? const_cast<Molecule&>(mol).GyrationRadiusPBC()
        : const_cast<Molecule&>(mol).GyrationRadius();
    geometric["gyration_radius"] = {
        {"unweighted", gyr.first},
        {"mass_weighted", gyr.second}
    };

    // Center of mass
    Position com = mol.MassCentroid();
    geometric["center_of_mass"] = {com[0], com[1], com[2]};

    // Geometric center
    Position centroid = mol.Centroid();
    geometric["centroid"] = {centroid[0], centroid[1], centroid[2]};

    // Inertia constants (if molecule has been analyzed)
    try {
        // Note: This requires CalculateRotationalConstants() to have been called
        geometric["inertia_constants"] = {
            {"Ia", mol.Ia()},
            {"Ib", mol.Ib()},
            {"Ic", mol.Ic()}
        };
    } catch (...) {
        // Inertia constants not available
    }

    return geometric;
}

json UnifiedAnalysis::calculatePolymerProperties(const Molecule& mol)
{
    json polymer;

    // End-to-end distance (for chain structures)
    // Claude Generated 2025: Use PBC-aware distance calculation when periodic boundaries are present
    // This ensures correct end-to-end distances for chains that cross periodic boundaries
    double end_to_end = mol.EndToEndDistancePBC();
    if (end_to_end > 0.0) {
        polymer["end_to_end_distance"] = end_to_end;
    }

    // Rout - maximum distance from COM to outermost atom
    // Claude Generated 2025: Use PBC-aware distance calculation when periodic boundaries are present
    double rout = mol.hasPBC() ? mol.RoutPBC() : mol.Rout();
    polymer["rout"] = rout;

    // Per-fragment analysis if multiple fragments
    if (mol.GetFragments().size() > 1) {
        json fragments = json::array();
        for (int i = 0; i < mol.GetFragments().size(); ++i) {
            json frag;
            frag["fragment_id"] = i;

            // Claude Generated 2025: Use PBC-aware calculation for fragments
            auto frag_gyr = mol.hasPBC()
                ? const_cast<Molecule&>(mol).GyrationRadiusPBC(1.0, true, i)
                : const_cast<Molecule&>(mol).GyrationRadius(1.0, true, i);
            frag["gyration_radius"] = frag_gyr.first;

            // Claude Generated 2025: Use PBC-aware distance for fragment end-to-end
            double frag_end_to_end = mol.hasPBC()
                ? mol.EndToEndDistancePBC(i)
                : mol.EndToEndDistance(i);
            if (frag_end_to_end > 0.0) {
                frag["end_to_end_distance"] = frag_end_to_end;
            }

            // Claude Generated 2025: Use PBC-aware Rout for fragments
            double frag_rout = mol.hasPBC() ? mol.RoutPBC(i) : mol.Rout(i);
            frag["rout"] = frag_rout;

            fragments.push_back(frag);
        }
        polymer["fragments"] = fragments;
    }

    return polymer;
}

json UnifiedAnalysis::calculateTopologicalProperties(const Molecule& mol)
{
    json topology;

    try {
        // Check if enhanced topological analysis is requested - Claude Generated 2025
        bool enhanced_tda = (m_config.get<bool>("topological.save_distance_matrix") ||
                            m_config.get<bool>("topological.save_persistence_pairs") ||
                            m_config.get<bool>("topological.save_persistence_diagram") ||
                            m_config.get<bool>("topological.save_persistence_image"));

        // Parse atom selection if provided - Claude Generated
        std::vector<int> selected_indices;
        std::string atom_selection = m_config.get<std::string>("topological.atom_selection");
        if (!atom_selection.empty()) {
            try {
                selected_indices = Tools::ParseStringToVector(atom_selection);

                // Validate indices against molecule size - Claude Generated
                if (!selected_indices.empty()) {
                    std::vector<int> valid_indices;
                    std::vector<int> invalid_indices;

                    for (int idx : selected_indices) {
                        if (idx >= 0 && idx < mol.AtomCount()) {
                            valid_indices.push_back(idx);
                        } else {
                            invalid_indices.push_back(idx);
                        }
                    }

                    if (!invalid_indices.empty() && !m_silent) {
                        std::string invalid_str = "";
                        for (size_t i = 0; i < invalid_indices.size(); ++i) {
                            if (i > 0) invalid_str += ",";
                            invalid_str += std::to_string(invalid_indices[i]);
                        }
                        CurcumaLogger::warn_fmt("Ignoring invalid atom indices (molecule has {} atoms): {}",
                                                mol.AtomCount(), invalid_str);
                    }

                    selected_indices = valid_indices;

                    if (!m_silent) {
                        if (!selected_indices.empty()) {
                            // Detailed selection feedback - Claude Generated
                            std::string indices_str = "[";
                            std::map<int, int> element_count;
                            double min_dist = std::numeric_limits<double>::max();
                            double max_dist = 0.0;

                            // Build indices string and collect element statistics
                            for (size_t i = 0; i < selected_indices.size(); ++i) {
                                if (i > 0) indices_str += ",";
                                int idx = selected_indices[i];
                                indices_str += std::to_string(idx);

                                // Count elements
                                int element = mol.Atom(idx).first;
                                element_count[element]++;

                                // Calculate spatial extent
                                for (size_t j = i + 1; j < selected_indices.size(); ++j) {
                                    double dist = mol.CalculateDistance(selected_indices[i], selected_indices[j]);
                                    min_dist = std::min(min_dist, dist);
                                    max_dist = std::max(max_dist, dist);
                                }
                            }
                            indices_str += "]";

                            // Build element composition string
                            std::string elements_str;
                            for (const auto& elem : element_count) {
                                if (!elements_str.empty()) elements_str += ", ";
                                std::string symbol = (elem.first > 0 && elem.first < Elements::ElementAbbr.size()) ? Elements::ElementAbbr[elem.first] : "X";
                                elements_str += symbol + "(" + std::to_string(elem.second) + ")";
                            }

                            // Main feedback message
                            double percentage = (100.0 * selected_indices.size()) / mol.AtomCount();
                            CurcumaLogger::success_fmt("✓ Atom selection: {}/{} atoms ({:.1f}%) selected: {}",
                                                     selected_indices.size(), mol.AtomCount(), percentage, indices_str);

                            // Additional details
                            if (selected_indices.size() > 1 && max_dist > 0.0) {
                                CurcumaLogger::info_fmt("  Elements: {} | Spatial extent: {:.2f} Å",
                                                      elements_str, max_dist);
                            } else if (selected_indices.size() == 1) {
                                CurcumaLogger::info_fmt("  Single atom: {} at index {}",
                                                      (mol.Atom(selected_indices[0]).first > 0 && mol.Atom(selected_indices[0]).first < Elements::ElementAbbr.size()) ? Elements::ElementAbbr[mol.Atom(selected_indices[0]).first] : "X",
                                                      selected_indices[0]);
                            } else {
                                CurcumaLogger::info_fmt("  Elements: {}", elements_str);
                            }

                            // Helpful hints
                            if (percentage < 10.0) {
                                CurcumaLogger::warn("  Small selection (<10%) - consider including more atoms for meaningful topology");
                            } else if (percentage > 90.0) {
                                CurcumaLogger::info("  Large selection (>90%) - consider full molecule analysis");
                            }

                            // Hydrogen exclusion info - Claude Generated 2025
                            bool has_hydrogen = false;
                            for (int idx : selected_indices) {
                                if (mol.Atom(idx).first == 1) {
                                    has_hydrogen = true;
                                    break;
                                }
                            }
                            if (has_hydrogen && m_config.get<bool>("topological.exclude_hydrogen")) {
                                CurcumaLogger::info("  Note: Hydrogen atoms will be excluded from distance calculations");
                            }

                        } else {
                            CurcumaLogger::warn("No valid atom indices found - using full molecule analysis");
                        }
                    }
                }
            } catch (const std::exception& e) {
                if (!m_silent) {
                    CurcumaLogger::error_fmt("Failed to parse atom selection '{}': {} - using full molecule analysis",
                                            atom_selection, e.what());
                    CurcumaLogger::info("  Syntax help: Use format \"1;5-10;15\" (semicolon-separated, dash for ranges)");
                    CurcumaLogger::info("  Examples: \"0;2;4\" or \"1-5;10-15\" or \"0;3-7;12\"");
                }
                selected_indices.clear();
            }
        }

        if (enhanced_tda) {
            // Use comprehensive TDA Engine for full dMatrix functionality - Claude Generated 2025
            // Create TDA config using ConfigManager with elegant dot notation
            json tda_config = {
                {"save_distance_matrix", m_config.get<bool>("topological.save_distance_matrix")},
                {"save_persistence_pairs", m_config.get<bool>("topological.save_persistence_pairs")},
                {"save_persistence_diagram", m_config.get<bool>("topological.save_persistence_diagram")},
                {"save_persistence_image", m_config.get<bool>("topological.save_persistence_image")},
                {"exclude_bonds", m_config.get<bool>("topological.exclude_bonds")},
                {"exclude_hydrogen", m_config.get<bool>("topological.exclude_hydrogen")},
                {"print_elements", m_config.get<bool>("topological.print_elements")},
                {"print_energy", m_config.get<bool>("topological.print_energy")},
                {"image_format", m_config.get<std::string>("topological.image_format")},
                {"colormap", m_config.get<std::string>("topological.colormap")},
                {"resolution", m_config.get<std::string>("topological.resolution")},
                {"post_processing", m_config.get<std::string>("topological.post_processing")},
                {"temperature", m_config.get<double>("topological.temperature")},
                {"damping", m_config.get<double>("topological.damping")},
                {"preserve_structure", m_config.get<bool>("topological.preserve_structure")}
            };

            // Set output prefix from filename if available
            if (!m_filename.empty()) {
                std::string prefix = m_filename;
                // Remove extension
                size_t lastdot = prefix.find_last_of(".");
                if (lastdot != std::string::npos) {
                    prefix = prefix.substr(0, lastdot);
                }
                tda_config["output_prefix"] = prefix;
            }

            // Auto-detect multi-structure files for TDA (legacy dMatrix compatibility) - Claude Generated
            FileIterator structure_detector(m_filename, true); // silent mode
            int total_structures = structure_detector.MaxMolecules();

            TDAEngine tda_engine(tda_config);

            if (total_structures > 1) {
                // Multi-structure file: Process all structures like legacy dMatrix
                if (!m_silent) {
                    CurcumaLogger::info_fmt("Detected multi-structure file with {} structures - processing all (legacy dMatrix compatibility)", total_structures);
                }

                FileIterator file_iter(m_filename, m_silent);
                json all_structures = json::array();

                int structure_index = 0;
                while (!file_iter.AtEnd()) {
                    Molecule current_mol = file_iter.Next();
                    json structure_result = selected_indices.empty() ? tda_engine.analyzeMolecule(current_mol, structure_index) : tda_engine.analyzeMolecule(current_mol, selected_indices, structure_index);
                    structure_result["structure_index"] = structure_index;
                    all_structures.push_back(structure_result);

                    if (!m_silent && structure_index % 10 == 0) {
                        CurcumaLogger::info_fmt("Processed structure {}/{}", structure_index + 1, total_structures);
                    }
                    structure_index++;
                }

                topology["multi_structure_analysis"] = all_structures;
                topology["total_structures"] = total_structures;
                topology["analysis_type"] = "multi_structure_tda";
            } else {
                // Single structure: Use existing logic
                topology = selected_indices.empty() ? tda_engine.analyzeMolecule(mol, 0) : tda_engine.analyzeMolecule(mol, selected_indices, 0);
                topology["analysis_type"] = "single_structure_tda";
            }

            // Add citation info for enhanced analysis
            topology["citation"] = "Townsend, J., Micucci, C.P., Hymel, J.H. et al. "
                                 "Representation of molecular structures with persistent homology "
                                 "for machine learning applications in chemistry. "
                                 "Nat Commun 11, 3230 (2020). https://doi.org/10.1038/s41467-020-17035-5";

            if (!m_silent) {
                CurcumaLogger::success("Enhanced topological analysis completed with dMatrix functionality");
            }

        } else {
            // Basic topological analysis (original implementation)
            PersistentDiagram pd(m_config_legacy);  // Legacy uses JSON

            // Use atom selection if provided, otherwise use all atoms - Claude Generated 2025
            if (!selected_indices.empty()) {
                bool exclude_hydrogen = m_config.get<bool>("topological.exclude_hydrogen");
                pd.setDistanceMatrix(mol.LowerDistanceVector(exclude_hydrogen, selected_indices));
                if (!m_silent) {
                    CurcumaLogger::info("Basic topological analysis with atom selection");
                }
            } else {
                pd.setDistanceMatrix(mol.LowerDistanceVector());
            }

            auto pairs = pd.generatePairs();
            topology["persistent_pairs_count"] = pairs.size();

            // Generate persistence image
            Eigen::MatrixXd image = pd.generateImage(pairs);
            topology["persistence_image_norm"] = image.norm();

            // Basic topological invariants (simplified)
            topology["betti_numbers"] = {
                {"b0", 1},  // Always at least one connected component
                {"b1", 0}   // Simplified - would need complex calculation from pairs
            };
        }

    } catch (const std::exception& e) {
        topology["error"] = e.what();
        CurcumaLogger::error("Topological analysis failed: " + std::string(e.what()));
    }

    return topology;
}

bool UnifiedAnalysis::detectChainStructure(const Molecule& mol)
{
    // Simple heuristic: if molecule is linear (bonds form a chain)
    // or if it's a CG structure, consider it a polymer/chain

    // Check if contains CG elements
    for (int atom : mol.Atoms()) {
        if (atom == CG_ELEMENT) {
            return true;  // CG structures are likely polymers
        }
    }

    // Check if structure is elongated (simple chain detection)
    auto gyr = const_cast<Molecule&>(mol).GyrationRadius();
    double rout = mol.Rout();

    // If Rout is significantly larger than gyration radius, likely chain-like
    return (rout > 2.0 * gyr.first);
}

void UnifiedAnalysis::outputResults(const json& results)
{
    std::string output_format = m_config.get<std::string>("output_format");
    std::string output_file = m_config.get<std::string>("output_file");

    if (!output_file.empty()) {
        outputToFile(results, output_file);
        return;
    }

    if (output_format == "json") {
        std::cout << std::setw(2) << results << std::endl;
    } else if (output_format == "csv") {
        // Use TrajectoryWriter for CSV output
        json writer_config;
        writer_config["default_format"] = "CSV";
        writer_config["precision"] = 3;
        TrajectoryWriter writer(writer_config);
        writer.writeCSV(std::cout, results);
        return;
    }

    if (output_format == "human") {
        // Human-readable output
        CurcumaLogger::info("=== Molecular Analysis Results ===");

        // Create modifiable copy for single-structure detailed view
        json display_results = results;

        if (display_results.contains("filename")) {
            CurcumaLogger::info_fmt("File: {}", display_results["filename"].get<std::string>());
        }

        // For single-structure files, extract timestep[0] for detailed view
        // For multi-structure files (trajectories), show table instead
        if (display_results.contains("timesteps") && display_results["timesteps"].size() == 1) {
            const auto& timestep = display_results["timesteps"][0];
            // Copy timestep properties to top level for detailed display below
            if (timestep.contains("basic"))
                display_results["basic"] = timestep["basic"];
            if (timestep.contains("geometric"))
                display_results["geometric"] = timestep["geometric"];
            if (timestep.contains("polymer"))
                display_results["polymer"] = timestep["polymer"];
            if (timestep.contains("topology"))
                display_results["topology"] = timestep["topology"];
        }

        if (display_results.contains("basic")) {
            const auto& basic = display_results["basic"];
            CurcumaLogger::info("");
            CurcumaLogger::info("Basic Properties:");
            CurcumaLogger::param("atoms", fmt::format("{}", basic["atom_count"].get<int>()));
            CurcumaLogger::param("mass", fmt::format("{:.3f} amu", basic["mass"].get<double>()));
            CurcumaLogger::param("charge", fmt::format("{}", basic["charge"].get<int>()));
            CurcumaLogger::param("fragments", fmt::format("{}", basic["fragment_count"].get<int>()));

            if (basic.contains("formula")) {
                CurcumaLogger::param("formula", fmt::format("{}", basic["formula"].get<std::string>()));
            }

            if (basic["coarse_grained"].get<bool>()) {
                CurcumaLogger::param("type", "Coarse-grained structure");
            }
        }

        if (display_results.contains("geometric")) {
            const auto& geom = display_results["geometric"];
            CurcumaLogger::info("");
            CurcumaLogger::info("Geometric Properties:");

            const auto& gyr = geom["gyration_radius"];
            CurcumaLogger::info_fmt("  Gyration radius: {:.3f} Å (unweighted)",
                gyr["unweighted"].get<double>());
            CurcumaLogger::info_fmt("                   {:.3f} Å (mass-weighted)",
                gyr["mass_weighted"].get<double>());

            const auto& com = geom["center_of_mass"];
            CurcumaLogger::info_fmt("  Center of mass: ({:.3f}, {:.3f}, {:.3f}) Å",
                com[0].get<double>(), com[1].get<double>(), com[2].get<double>());
        }

        if (display_results.contains("polymer")) {
            const auto& polymer = display_results["polymer"];
            CurcumaLogger::info("");
            CurcumaLogger::info("Polymer/Chain Properties:");

            if (polymer.contains("end_to_end_distance")) {
                CurcumaLogger::info_fmt("  End-to-end distance: {:.3f} Å",
                    polymer["end_to_end_distance"].get<double>());
            }

            CurcumaLogger::info_fmt("  Rout (COM to outermost): {:.3f} Å",
                polymer["rout"].get<double>());
        }

        if (display_results.contains("topology")) {
            const auto& topo = display_results["topology"];
            CurcumaLogger::info("");
            CurcumaLogger::info("Topological Properties:");

            if (topo.contains("persistent_pairs_count")) {
                CurcumaLogger::param("  Persistent pairs: {}",
                                       topo["persistent_pairs_count"].get<int>());
            }
        }

        // Trajectory tabular output - Use TrajectoryWriter for consistency
        // Only show table for multi-structure files (trajectories)
        if (results.contains("timesteps") && results["timesteps"].size() > 1) {
            // Create TrajectoryWriter with appropriate config for human table
            json writer_config;
            writer_config["default_format"] = "HumanTable";
            writer_config["precision"] = 3;

            // Configure statistics for trajectory data
            if (results.contains("statistics_config")) {
                writer_config["enable_cumulative"] = results["statistics_config"]["enable_cumulative"];
                writer_config["enable_moving"] = results["statistics_config"]["enable_moving"];
            }

            TrajectoryWriter writer(writer_config);
            writer.writeHumanTable(std::cout, results);

            CurcumaLogger::info("");
            CurcumaLogger::info_fmt("Trajectory analysis completed: {} timesteps",
                                  results["total_timesteps"].get<int>());
        }
    }
}

void UnifiedAnalysis::outputToFile(const json& results, const std::string& filename)
{
    std::string output_format = m_config.get<std::string>("output_format");

    // Create TrajectoryWriter with appropriate config
    json writer_config;
    writer_config["default_format"] = (output_format == "csv" ? "CSV" : output_format == "json" ? "JSON" : "HumanTable");
    writer_config["precision"] = 3;

    // Configure statistics for trajectory data
    if (results.contains("statistics_config")) {
        writer_config["enable_cumulative"] = results["statistics_config"]["enable_cumulative"];
        writer_config["enable_moving"] = results["statistics_config"]["enable_moving"];
    }

    TrajectoryWriter writer(writer_config);

    // Map format strings to TrajectoryWriter::Format
    TrajectoryWriter::Format format = TrajectoryWriter::Format::HumanTable; // default
    if (output_format == "csv") {
        format = TrajectoryWriter::Format::CSV;
    } else if (output_format == "json") {
        format = TrajectoryWriter::Format::JSON;
    } else if (output_format == "dat") {
        format = TrajectoryWriter::Format::DAT;
    }

    // Output using TrajectoryWriter
    writer.writeToFile(filename, format, results);
    CurcumaLogger::success_fmt("Analysis results saved to: {} (format: {})", filename, output_format);
}

void UnifiedAnalysis::printHelp() const
{
    ParameterRegistry::getInstance().printHelp("analysis");
}

void UnifiedAnalysis::printEnhancedTDAHelp() const
{
    std::cout << "Enhanced Topological Data Analysis (dMatrix Alternative)" << std::endl;
    std::cout << "=========================================================" << std::endl;
    std::cout << std::endl;
    std::cout << "OVERVIEW:" << std::endl;
    std::cout << "The enhanced TDA functionality provides comprehensive topological data analysis" << std::endl;
    std::cout << "for molecular structures, equivalent to the legacy -dMatrix command but with" << std::endl;
    std::cout << "modern integration into the unified analysis system." << std::endl;
    std::cout << std::endl;
    std::cout << "MIGRATION FROM -dMatrix:" << std::endl;
    std::cout << "Old usage: curcuma -dMatrix molecule.xyz" << std::endl;
    std::cout << "New usage: curcuma -analysis molecule.xyz -topological.save_persistence_image true" << std::endl;
    std::cout << std::endl;
    std::cout << "CORE FEATURES:" << std::endl;
    std::cout << std::endl;
    std::cout << "1. Distance Matrix Analysis" << std::endl;
    std::cout << "   -topological.save_distance_matrix true" << std::endl;
    std::cout << "   • Exports complete distance matrix (.dMat files)" << std::endl;
    std::cout << "   • Includes atomic elements and energies if requested" << std::endl;
    std::cout << "   • Supports bond exclusion for topological analysis" << std::endl;
    std::cout << std::endl;
    std::cout << "2. Persistent Homology Analysis" << std::endl;
    std::cout << "   -topological.save_persistence_pairs true" << std::endl;
    std::cout << "   • Generates persistence pairs (.pairs files)" << std::endl;
    std::cout << "   • Calculates birth-death persistence data" << std::endl;
    std::cout << "   • Provides topological invariants and statistics" << std::endl;
    std::cout << std::endl;
    std::cout << "3. Persistence Diagrams" << std::endl;
    std::cout << "   -topological.save_persistence_diagram true" << std::endl;
    std::cout << "   • Creates persistence diagram matrices (.PD files)" << std::endl;
    std::cout << "   • Exports both text and image formats" << std::endl;
    std::cout << "   • Supports multiple colormaps and post-processing" << std::endl;
    std::cout << std::endl;
    std::cout << "4. Persistence Images" << std::endl;
    std::cout << "   -topological.save_persistence_image true" << std::endl;
    std::cout << "   • Generates persistence images (.PI files)" << std::endl;
    std::cout << "   • EN-scaled bond topology analysis" << std::endl;
    std::cout << "   • Advanced image processing and enhancement" << std::endl;
    std::cout << std::endl;
    std::cout << "ADVANCED CONFIGURATION:" << std::endl;
    std::cout << std::endl;
    std::cout << "Molecular Filtering:" << std::endl;
    std::cout << "  -topological.exclude_hydrogen true     # Remove H atoms from analysis" << std::endl;
    std::cout << "  -topological.exclude_bonds true        # Exclude bonded atom pairs" << std::endl;
    std::cout << "  -topological.print_elements true       # Include element symbols" << std::endl;
    std::cout << "  -topological.print_energy true         # Include molecular energies" << std::endl;
    std::cout << std::endl;
    std::cout << "Image Generation:" << std::endl;
    std::cout << "  -topological.image_format png|jpg|bmp|tga" << std::endl;
    std::cout << "  -topological.colormap grayscale|jet|hot|viridis|coolwarm" << std::endl;
    std::cout << "  -topological.resolution 800x600|1024x1024|2048x2048" << std::endl;
    std::cout << std::endl;
    std::cout << "Post-Processing:" << std::endl;
    std::cout << "  -topological.post_processing none|adaptive|ring_focused" << std::endl;
    std::cout << "  -topological.temperature 2.0           # Enhancement temperature parameter" << std::endl;
    std::cout << "  -topological.damping 1.5               # Damping strength for processing" << std::endl;
    std::cout << "  -topological.preserve_structure true   # Maintain original structure" << std::endl;
    std::cout << std::endl;
    std::cout << "MULTI-STRUCTURE & TRAJECTORY ANALYSIS:" << std::endl;
    std::cout << "The enhanced TDA supports batch processing of molecular trajectories and" << std::endl;
    std::cout << "multi-structure files with automatic detection:" << std::endl;
    std::cout << std::endl;
    std::cout << "Supported Multi-Structure Formats:" << std::endl;
    std::cout << "• XYZ trajectories (.trj.xyz, multi-frame .xyz)" << std::endl;
    std::cout << "• VTF trajectories (.vtf with multiple timesteps)" << std::endl;
    std::cout << "• MOL2, SDF multi-structure files" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage Examples:" << std::endl;
    std::cout << "  curcuma -analysis trajectory.trj.xyz -topological.save_distance_matrix true" << std::endl;
    std::cout << "  curcuma -analysis dynamics.vtf -topological.save_persistence_image true" << std::endl;
    std::cout << "  curcuma -analysis multi_conf.xyz -topological.save_persistence_pairs true" << std::endl;
    std::cout << std::endl;
    std::cout << "Auto-Detection Features:" << std::endl;
    std::cout << "• Automatic multi-structure detection (legacy dMatrix compatibility)" << std::endl;
    std::cout << "• Per-frame analysis with automatic indexing" << std::endl;
    std::cout << "• VTF timestep parsing with structure template preservation" << std::endl;
    std::cout << "• Progress reporting for large trajectory files" << std::endl;
    std::cout << "• Statistical analysis across trajectory (average, standard deviation)" << std::endl;
    std::cout << "• Correlation analysis between frames" << std::endl;
    std::cout << "• Batch export of all persistence data" << std::endl;
    std::cout << std::endl;
    std::cout << "OUTPUT FILES:" << std::endl;
    std::cout << "The system generates indexed files for each analysis type:" << std::endl;
    std::cout << std::endl;
    std::cout << "• prefix_N.dMat     - Distance matrix (text format)" << std::endl;
    std::cout << "• prefix_N.pairs    - Persistence pairs (birth-death values)" << std::endl;
    std::cout << "• prefix_N.PD       - Persistence diagram (matrix format)" << std::endl;
    std::cout << "• prefix_N.PD.png   - Persistence diagram (image)" << std::endl;
    std::cout << "• prefix_N.PI       - Persistence image (matrix format)" << std::endl;
    std::cout << "• prefix_N.PI.png   - Persistence image (visualization)" << std::endl;
    std::cout << std::endl;
    std::cout << "Where N is the frame index (0 for single structures)." << std::endl;
    std::cout << std::endl;
    std::cout << "RESEARCH APPLICATIONS:" << std::endl;
    std::cout << "This implementation is suitable for:" << std::endl;
    std::cout << "• Molecular shape analysis and comparison" << std::endl;
    std::cout << "• Protein folding and conformational studies" << std::endl;
    std::cout << "• Materials science and crystal structure analysis" << std::endl;
    std::cout << "• Drug discovery and molecular recognition" << std::endl;
    std::cout << "• Machine learning feature generation" << std::endl;
    std::cout << std::endl;
    std::cout << "CITATION:" << std::endl;
    std::cout << "When using enhanced TDA features, please cite:" << std::endl;
    std::cout << "Townsend, J., Micucci, C.P., Hymel, J.H. et al." << std::endl;
    std::cout << "Representation of molecular structures with persistent homology" << std::endl;
    std::cout << "for machine learning applications in chemistry." << std::endl;
    std::cout << "Nat Commun 11, 3230 (2020). https://doi.org/10.1038/s41467-020-17035-5" << std::endl;
    std::cout << std::endl;
    std::cout << "EXAMPLES:" << std::endl;
    std::cout << std::endl;
    std::cout << "Basic persistence analysis:" << std::endl;
    std::cout << "  curcuma -analysis molecule.xyz -topological.save_persistence_pairs true" << std::endl;
    std::cout << std::endl;
    std::cout << "Complete TDA workflow:" << std::endl;
    std::cout << "  curcuma -analysis structure.xyz \\" << std::endl;
    std::cout << "          -topological.save_distance_matrix true \\" << std::endl;
    std::cout << "          -topological.save_persistence_diagram true \\" << std::endl;
    std::cout << "          -topological.save_persistence_image true" << std::endl;
    std::cout << std::endl;
    std::cout << "High-quality images without hydrogen:" << std::endl;
    std::cout << "  curcuma -analysis protein.pdb \\" << std::endl;
    std::cout << "          -topological.save_persistence_image true \\" << std::endl;
    std::cout << "          -topological.exclude_hydrogen true \\" << std::endl;
    std::cout << "          -topological.colormap viridis \\" << std::endl;
    std::cout << "          -topological.resolution 2048x2048 \\" << std::endl;
    std::cout << "          -topological.image_format png" << std::endl;
    std::cout << std::endl;
    std::cout << "Trajectory analysis with post-processing:" << std::endl;
    std::cout << "  curcuma -analysis dynamics.trj.xyz \\" << std::endl;
    std::cout << "          -trajectory true \\" << std::endl;
    std::cout << "          -topological.save_persistence_image true \\" << std::endl;
    std::cout << "          -topological.post_processing adaptive \\" << std::endl;
    std::cout << "          -topological.temperature 3.0" << std::endl;
    std::cout << std::endl;
    std::cout << "VTF trajectory analysis (coarse-grained systems):" << std::endl;
    std::cout << "  curcuma -analysis polymer_dynamics.vtf \\" << std::endl;
    std::cout << "          -topological.save_distance_matrix true \\" << std::endl;
    std::cout << "          -topological.save_persistence_image true \\" << std::endl;
    std::cout << "          -topological.exclude_bonds true" << std::endl;
    std::cout << std::endl;
    std::cout << "Multi-structure file with automatic detection:" << std::endl;
    std::cout << "  curcuma -analysis conformers.xyz \\" << std::endl;
    std::cout << "          -topological.save_persistence_pairs true" << std::endl;
    std::cout << "  # Automatically processes all structures in the file" << std::endl;
    std::cout << std::endl;
}