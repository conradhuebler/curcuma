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
#include "analysis_output.h"  // Claude Generated 2026: Output dispatcher
#include "persistentdiagram.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"
#include "src/core/units.h"
#include "src/core/parameter_registry.h"
#include "src/core/form_factors.h"
#include "src/tools/pbc_utils.h"
#include "src/tools/general.h"  // For Tools::CreateList - Claude Generated 2026
#include "trajectory_statistics.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

UnifiedAnalysis::UnifiedAnalysis(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("analysis"), controller, silent)
    , m_config("analysis", controller)  // Claude Generated 2025: ConfigManager accepts full controller
    , m_silent(silent)
    , m_analysis_config(m_config)  // Claude Generated 2026: Initialize configuration state
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

// AnalysisConfig implementation - Claude Generated 2026
UnifiedAnalysis::AnalysisConfig::AnalysisConfig(const ConfigManager& config)
{
    // Frame selection parameters
    frames_str = config.get<std::string>("frames", "");
    stride = config.get<int>("stride", 1);
    start_frame = config.get<int>("frame_range_start", 0);
    end_frame = config.get<int>("frame_range_end", -1);

    // Scattering parameters
    scattering_enabled = config.get<bool>("scattering_enable", false);
    scattering_per_frame_files = config.get<bool>("scattering_per_frame_files", false);
    scattering_q_values = config.get<std::string>("scattering_q_values", "");
    scattering_output_directory = config.get<std::string>("scattering_output_directory", ".");
    scattering_file_prefix = config.get<std::string>("scattering_file_prefix", "scattering_frame");
    scattering_stats_include_median = config.get<bool>("scattering_stats_include_median", true);

    // RDF parameters
    rdf_enabled = config.get<bool>("rdf_enable", false);
    rdf_r_max = config.get<double>("rdf_r_max", 15.0);
    rdf_bin_width = config.get<double>("rdf_bin_width", 0.05);
    rdf_coordination_shells = config.get<bool>("rdf_coordination_shells", false);

    // Output parameters
    output_format = config.get<std::string>("output_format", "human");
    output_file = config.get<std::string>("output_file", "");
    metrics = config.get<std::string>("metrics", "gyration,rout,end2end");
    statistics_mode = config.get<std::string>("statistics", "none");
    window_size = config.get<int>("window", 10);

    // Validate configuration
    validateConfiguration();
}

void UnifiedAnalysis::AnalysisConfig::validateConfiguration()
{
    // Positive stride
    if (stride <= 0) {
        CurcumaLogger::warn("stride must be positive, using 1");
        stride = 1;
    }

    // Dependency check: per_frame_files requires scattering
    if (scattering_per_frame_files && !scattering_enabled) {
        CurcumaLogger::warn("scattering_per_frame_files requires scattering_enable=true, disabling per_frame_files");
        scattering_per_frame_files = false;
    }

    // Window size for moving statistics
    if (window_size < 2 && statistics_mode.find("moving") != std::string::npos) {
        CurcumaLogger::warn("window size must be >= 2 for moving statistics, using 10");
        window_size = 10;
    }
}

bool UnifiedAnalysis::AnalysisConfig::hasScatteringOutput() const
{
    return scattering_enabled;
}

bool UnifiedAnalysis::AnalysisConfig::hasPerFrameFiles() const
{
    return scattering_enabled && scattering_per_frame_files;
}

bool UnifiedAnalysis::AnalysisConfig::hasRDFOutput() const
{
    return rdf_enabled;
}

json UnifiedAnalysis::AnalysisConfig::toJSON() const
{
    json config_json;

    // Frame selection
    config_json["frames"] = frames_str;
    config_json["stride"] = stride;
    config_json["frame_range_start"] = start_frame;
    config_json["frame_range_end"] = end_frame;

    // Scattering
    if (scattering_enabled) {
        config_json["scattering_enable"] = true;
        config_json["scattering_per_frame_files"] = scattering_per_frame_files;
        config_json["scattering_q_values"] = scattering_q_values;
        config_json["scattering_output_directory"] = scattering_output_directory;
        config_json["scattering_file_prefix"] = scattering_file_prefix;
        config_json["scattering_stats_include_median"] = scattering_stats_include_median;
    }

    // RDF
    if (rdf_enabled) {
        config_json["rdf_enable"] = true;
        config_json["rdf_r_max"] = rdf_r_max;
        config_json["rdf_bin_width"] = rdf_bin_width;
        config_json["rdf_coordination_shells"] = rdf_coordination_shells;
    }

    // Output
    config_json["output_format"] = output_format;
    config_json["output_file"] = output_file;
    config_json["metrics"] = metrics;
    config_json["statistics"] = statistics_mode;
    config_json["window"] = window_size;

    return config_json;
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

    // Parse statistics configuration - Claude Generated 2026 (Refactored to use AnalysisConfig)
    std::vector<std::string> enabled_metrics = parseMetricsList(m_analysis_config.metrics);
    bool enable_cumulative = (m_analysis_config.statistics_mode == "cumulative" || m_analysis_config.statistics_mode == "all");
    bool enable_moving = (m_analysis_config.statistics_mode == "moving" || m_analysis_config.statistics_mode == "all");

    // Create statistics tracker if needed
    std::unique_ptr<TrajectoryStatistics> stats;
    if (enable_cumulative || enable_moving) {
        stats = std::make_unique<TrajectoryStatistics>(m_analysis_config.window_size);
    }

    // ALWAYS store configuration in results for output formatting
    results["statistics_config"] = {
        { "metrics", enabled_metrics },
        { "mode", m_analysis_config.statistics_mode },
        { "window", m_analysis_config.window_size },
        { "enable_cumulative", enable_cumulative },
        { "enable_moving", enable_moving }
    };

    // Store scattering configuration for per-frame file output - Claude Generated 2026 (Refactored)
    if (m_analysis_config.hasPerFrameFiles()) {
        results["config"] = {
            { "scattering_per_frame_files", m_analysis_config.scattering_per_frame_files },
            { "scattering_q_values", m_analysis_config.scattering_q_values }
        };
    }

    std::vector<int> selected_frames;
    bool use_frame_list = false;

    // Two-pass approach: First count total frames to resolve -1 - Claude Generated 2026
    int total_frames = 0;
    {
        FileIterator counter(m_filename, true); // silent mode
        while (!counter.AtEnd()) {
            counter.Next();
            total_frames++;
        }
    }

    // Resolve end_frame for range-based selection - Claude Generated 2026
    int end_frame = m_analysis_config.end_frame;
    if (end_frame == -1) {
        end_frame = total_frames;
    }

    // Resolve frame selection - Claude Generated 2026 (Refactored to use AnalysisConfig)
    if (!m_analysis_config.frames_str.empty()) {
        use_frame_list = true;

        // Special cases: "-1" or "last" alone means only the last frame
        if (m_analysis_config.frames_str == "-1" || m_analysis_config.frames_str == "last") {
            selected_frames.push_back(total_frames - 1); // 0-based
        } else {
            // Replace -1 with actual last frame number (1-based for user) in ranges
            std::string resolved_frames = m_analysis_config.frames_str;
            if (resolved_frames.find("-1") != std::string::npos) {
                size_t pos = 0;
                std::string from = "-1";
                std::string to = std::to_string(total_frames);
                while ((pos = resolved_frames.find(from, pos)) != std::string::npos) {
                    resolved_frames.replace(pos, from.length(), to);
                    pos += to.length();
                }
            }

            // Parse frame selection using Tools::CreateList (1-based input)
            selected_frames = Tools::CreateList(resolved_frames);

            // Convert 1-based indexing to 0-based for internal use
            for (auto& frame : selected_frames) {
                frame -= 1;
            }
        }
    }

    // Main analysis loop with frame filtering - Claude Generated 2026
    int current_frame = 0;
    int analyzed_frames = 0;
    bool first_analyzed = true; // Track first analyzed frame for PBC logging
    
    std::cout << "Selected frames for analysis: ";
    for (const auto i : selected_frames) {
        std::cout << i << " ";
        CurcumaLogger::info_fmt("Selected frame for analysis: {}", i);
        if (i < 0 || i >= total_frames) {
            CurcumaLogger::error_fmt("Selected frame {} is out of bounds (0 to {})", i, total_frames - 1);
            return;
        }
    }

    while (!file_iter.AtEnd()) {
        Molecule mol = file_iter.Next();

        bool should_analyze = false;

        if (use_frame_list) {
            // Check if current frame is in selected list
            auto it = std::find(selected_frames.begin(), selected_frames.end(), current_frame);
            if (it != selected_frames.end()) {
                // Apply stride filter to the list index
                int list_index = std::distance(selected_frames.begin(), it);
                should_analyze = (list_index % m_analysis_config.stride == 0);
            }
        } else {
            // Simple range check with stride (using local end_frame from line 316)
            if (current_frame >= m_analysis_config.start_frame && current_frame < end_frame) {
                should_analyze = ((current_frame - m_analysis_config.start_frame) % m_analysis_config.stride == 0);
            }
        }
        if(!should_analyze){
            current_frame++;
            continue;
        }
        if (should_analyze) {
            // Claude Generated 2025: Log PBC detection at verbosity level 2+ (first analyzed frame only)
            if (first_analyzed && mol.hasPBC() && CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::info("Detected periodic boundary conditions");
                auto params = PBCUtils::getLatticeParameters(mol.getUnitCell());
                CurcumaLogger::param("lattice_a", fmt::format("{:.4f} Å", params[0]));
                CurcumaLogger::param("lattice_b", fmt::format("{:.4f} Å", params[1]));
                CurcumaLogger::param("lattice_c", fmt::format("{:.4f} Å", params[2]));
                CurcumaLogger::param("lattice_alpha", fmt::format("{:.2f}°", params[3]));
                CurcumaLogger::param("lattice_beta", fmt::format("{:.2f}°", params[4]));
                CurcumaLogger::param("lattice_gamma", fmt::format("{:.2f}°", params[5]));
                first_analyzed = false;
            }

            json timestep_result = analyzeMolecule(mol, current_frame);

            // Store actual frame index in result for traceability
            timestep_result["frame_index"] = current_frame;

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
            analyzed_frames++;

            // Progress reporting (every 10 analyzed frames) - Claude Generated 2026
            if (!m_silent && analyzed_frames % 10 == 0 && analyzed_frames > 0) {
                if (use_frame_list) {
                    CurcumaLogger::info_fmt("Analyzed {} frames (current: {}/{})",
                                           analyzed_frames, current_frame + 1, total_frames);
                } else if (total_frames > 1) {
                    CurcumaLogger::info_fmt("Analyzed frame {}/{}",
                                           current_frame + 1, total_frames);
                }
            }
        }

        current_frame++;
    }

    // Update total counts - Claude Generated 2026
    results["total_timesteps"] = analyzed_frames;
    results["total_frames_in_file"] = total_frames;

    if (!m_silent && analyzed_frames > 0) {
        if (use_frame_list) {
            CurcumaLogger::info_fmt("Processed {} frames from {} total frames",
                analyzed_frames, total_frames);
        } else if (analyzed_frames > 1) {
            CurcumaLogger::info_fmt("Processed {} structures with {} metrics selected",
                analyzed_frames, enabled_metrics.size());
        }
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

    // Calculate scattering properties if enabled - Claude Generated 2025
    // NOTE: Only skip for trajectories in non-JSON formats (Human/CSV causes crashes)
    bool scattering_enabled = false;
    try {
        scattering_enabled = m_config.get<bool>("scattering_enable");
    } catch (...) {
        scattering_enabled = false;
    }
    if (scattering_enabled) {
        bool is_trajectory = (timestep >= 0);
        std::string output_format = "human";
        try {
            output_format = m_config.get<std::string>("output_format");
        } catch (...) {
            output_format = "human";
        }
        bool is_json_output = (output_format == "json");

        // NOTE: With enhanced trajectory writer, scattering analysis now works for all formats
        result["scattering"] = calculateScatteringProperties(mol);
    }

    // Calculate radial distribution function if enabled - Claude Generated 2025
    // NOTE: Only skip for trajectories in non-JSON formats (Human/CSV causes crashes)
    bool rdf_enabled = false;
    try {
        rdf_enabled = m_config.get<bool>("rdf_enable");
    } catch (...) {
        rdf_enabled = false;
    }
    if (rdf_enabled) {
        bool is_trajectory = (timestep >= 0);
        std::string output_format_rdf = "human";
        try {
            output_format_rdf = m_config.get<std::string>("output_format");
        } catch (...) {
            output_format_rdf = "human";
        }
        bool is_json_output_rdf = (output_format_rdf == "json");

        // Only skip for trajectories in non-JSON formats (Human/CSV crash)
        if (is_trajectory && !is_json_output_rdf) {
            static bool warned_rdf = false;
            if (!warned_rdf) {
                CurcumaLogger::warn("RDF analysis on trajectories requires -output_format json");
                warned_rdf = true;
            }
        } else {
            result["rdf"] = calculatePairDistribution(mol);
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

    // Extended shape descriptors - Claude Generated 2025
    // Check if any shape parameters are enabled
    bool calc_asphericity = false;
    bool calc_acylindricity = false;
    bool calc_anisotropy = false;

    try {
        calc_asphericity = m_config.get<bool>("shape_asphericity");
        calc_acylindricity = m_config.get<bool>("shape_acylindricity");
        calc_anisotropy = m_config.get<bool>("shape_anisotropy");
    } catch (...) {
        // Use defaults if parameters not available
    }

    if (calc_asphericity || calc_acylindricity || calc_anisotropy) {
        // Compute gyration tensor: Gₐᵦ = (1/N) Σᵢ (rᵢₐ - COM_a)(rᵢᵦ - COM_b)
        Position com = mol.MassCentroid();
        const int N = mol.AtomCount();

        Eigen::Matrix3d gyration_tensor = Eigen::Matrix3d::Zero();

        for (int i = 0; i < N; ++i) {
            Position pos = mol.Atom(i).second;

            // Relative position from COM
            double dx = pos.x() - com[0];
            double dy = pos.y() - com[1];
            double dz = pos.z() - com[2];

            // Build symmetric tensor
            gyration_tensor(0, 0) += dx * dx;
            gyration_tensor(0, 1) += dx * dy;
            gyration_tensor(0, 2) += dx * dz;
            gyration_tensor(1, 0) += dy * dx;
            gyration_tensor(1, 1) += dy * dy;
            gyration_tensor(1, 2) += dy * dz;
            gyration_tensor(2, 0) += dz * dx;
            gyration_tensor(2, 1) += dz * dy;
            gyration_tensor(2, 2) += dz * dz;
        }

        gyration_tensor /= N;

        // Eigenvalue decomposition
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(gyration_tensor);
        Eigen::Vector3d eigenvalues = solver.eigenvalues();

        // Sort eigenvalues: λ₁ ≥ λ₂ ≥ λ₃
        std::sort(eigenvalues.data(), eigenvalues.data() + 3, std::greater<double>());
        double lambda1 = eigenvalues[0];
        double lambda2 = eigenvalues[1];
        double lambda3 = eigenvalues[2];

        json shape_descriptors;
        shape_descriptors["eigenvalues"] = {lambda1, lambda2, lambda3};

        // Asphericity: b = λ₃ - 0.5(λ₁ + λ₂)
        if (calc_asphericity) {
            double asphericity = lambda3 - 0.5 * (lambda1 + lambda2);
            shape_descriptors["asphericity"] = asphericity;

            // Interpretation
            if (asphericity > 0.1) {
                shape_descriptors["shape_type"] = "prolate";  // Cigar-like
            } else if (asphericity < -0.1) {
                shape_descriptors["shape_type"] = "oblate";   // Disk-like
            } else {
                shape_descriptors["shape_type"] = "spherical";
            }
        }

        // Acylindricity: c = λ₂ - λ₁
        if (calc_acylindricity) {
            double acylindricity = lambda2 - lambda1;
            shape_descriptors["acylindricity"] = acylindricity;
        }

        // Relative shape anisotropy: κ² = 1 - 3⟨I₁I₂I₃⟩/⟨I₁+I₂+I₃⟩²
        if (calc_anisotropy) {
            double lambda_sum = lambda1 + lambda2 + lambda3;
            double lambda_product = lambda1 * lambda2 * lambda3;

            double anisotropy = 0.0;
            if (lambda_sum > 0.0) {
                anisotropy = 1.0 - 3.0 * lambda_product / (lambda_sum * lambda_sum);
            }

            shape_descriptors["relative_anisotropy"] = anisotropy;
            // κ² ∈ [0, 1]: 0 = sphere, 1 = rod or disk
        }

        geometric["shape_descriptors"] = shape_descriptors;
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

    // Phase 4: ALWAYS write files (derive basename from output_file OR input filename)
    std::string basename;
    if (!output_file.empty()) {
        // Use explicit output_file as basename
        basename = output_file;
    } else {
        // Derive basename from input filename: "input.xyz" → "input"
        basename = m_filename;
        size_t last_dot = basename.find_last_of(".");
        if (last_dot != std::string::npos) {
            basename = basename.substr(0, last_dot);
        }
    }

    // ALWAYS dispatch to file handlers (write specialized files)
    outputToFile(results, basename);

    // Continue with console output for immediate feedback
    // (no more early return - console output ADDITIONAL to files)

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

            if (polymer.contains("end_to_end_distance") && !polymer["end_to_end_distance"].is_null()) {
                CurcumaLogger::info_fmt("  End-to-end distance: {:.3f} Å",
                    polymer["end_to_end_distance"].get<double>());
            }

            if (polymer.contains("rout") && !polymer["rout"].is_null()) {
                CurcumaLogger::info_fmt("  Rout (COM to outermost): {:.3f} Å",
                    polymer["rout"].get<double>());
            }
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

            // NOTE: TrajectoryWriter now supports scattering data (2026 enhancement)
            TrajectoryWriter writer(writer_config);
            writer.writeHumanTable(std::cout, results);

            CurcumaLogger::info("");
            CurcumaLogger::info_fmt("Trajectory analysis completed: {} timesteps",
                                  results["total_timesteps"].get<int>());
        }
    }

    // Phase 4: File output already handled by outputToFile() call above
    // No need for separate dispatcher call - handlers write directly to files
}

void UnifiedAnalysis::outputToFile(const json& results, const std::string& basename)
{
    // Phase 4: Only dispatch to handlers - no main file write
    // Each handler creates its own specialized CSV file (general, scattering, rdf, etc.)

    // Create TrajectoryWriter with CSV config for handler use
    json writer_config;
    writer_config["default_format"] = "CSV";
    writer_config["precision"] = 3;
    TrajectoryWriter writer(writer_config);

    // Dispatch to all registered handlers
    // Each handler writes its own file: basename.type.csv or basename.NNN.type.csv
    AnalysisOutputDispatcher dispatcher(m_analysis_config, writer);
    dispatcher.dispatchToFile(results, basename);

    // Handlers log their own success messages
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

// =====================================================================
// Scattering Analysis Implementation - Claude Generated 2025
// =====================================================================

json UnifiedAnalysis::calculateScatteringProperties(const Molecule& mol)
{
    json scattering;

    // Get parameters with defaults - Claude Generated 2025
    double q_min = 0.01;
    double q_max = 2.0;
    int q_steps = 100;
    std::string ff_type = "auto";
    double cg_radius = 3.0;
    int angular_samples = 50;

    try {
        q_min = m_config.get<double>("scattering_q_min");
        q_max = m_config.get<double>("scattering_q_max");
        q_steps = m_config.get<int>("scattering_q_steps");
        ff_type = m_config.get<std::string>("scattering_form_factor");
        cg_radius = m_config.get<double>("scattering_cg_radius");
        angular_samples = m_config.get<int>("scattering_angular_samples");
    } catch (...) {
        // Use defaults if parameters not available
    }


    // Auto-detect system type
    bool is_cg = mol.isCGSystem();
    if (ff_type == "auto") {
        ff_type = is_cg ? "cg_sphere" : "cromer_mann";
    }

    scattering["system_type"] = is_cg ? "cg" : "atomic";
    scattering["form_factor_type"] = ff_type;
    if (is_cg) {
        scattering["cg_radius"] = cg_radius;
    }

    // Claude Generated 2026: Build q-vector with log or linear spacing
    std::vector<double> q_values;
    std::string q_spacing = m_config.get<std::string>("scattering_q_spacing");

    if (q_spacing == "log" || q_spacing == "logarithmic") {
        // Logarithmic spacing: q[i] = q_min * (q_max/q_min)^(i/(N-1))
        // Provides better resolution at low q (Guinier region)

        if (q_min <= 0.0 || q_max <= 0.0) {
            CurcumaLogger::error("Logarithmic q-spacing requires q_min > 0 and q_max > 0");
            CurcumaLogger::warn("Falling back to linear spacing");
            q_spacing = "linear";
        } else if (q_max <= q_min) {
            CurcumaLogger::error_fmt("Invalid q-range: q_max ({}) must be > q_min ({})", q_max, q_min);
            return scattering;  // Empty result
        } else {
            double log_ratio = std::log(q_max / q_min);
            for (int i = 0; i < q_steps; ++i) {
                double exponent = static_cast<double>(i) / std::max(1, q_steps - 1);
                q_values.push_back(q_min * std::exp(exponent * log_ratio));
            }
            CurcumaLogger::info_fmt("Using logarithmic q-spacing ({} points from {} to {} Å⁻¹)",
                                    q_steps, q_min, q_max);
        }
    }

    if (q_spacing == "linear" || q_values.empty()) {
        // Linear spacing: q[i] = q_min + i * (q_max - q_min) / (N-1)
        double q_step = (q_max - q_min) / std::max(1, q_steps - 1);
        q_values.clear();
        for (int i = 0; i < q_steps; ++i) {
            q_values.push_back(q_min + i * q_step);
        }
        if (q_spacing == "linear") {
            CurcumaLogger::info_fmt("Using linear q-spacing ({} points from {} to {} Å⁻¹)",
                                    q_steps, q_min, q_max);
        }
    }

    scattering["q_min"] = q_min;
    scattering["q_max"] = q_max;
    scattering["q_steps"] = q_steps;
    scattering["q_values"] = q_values;

    // Pre-calculate distance matrix (reuse if available)
    const int N = mol.AtomCount();
    std::vector<std::vector<double>> distances(N, std::vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dist = mol.hasPBC()
                ? mol.CalculateDistancePBC(i, j)
                : mol.CalculateDistance(i, j);
            distances[i][j] = dist;
            distances[j][i] = dist;
        }
    }

    // ===== FORM FACTOR P(q) - Debye Formula =====
    // P(q) = (1/N²) Σᵢⱼ fᵢ(q)fⱼ(q) sin(qrᵢⱼ)/(qrᵢⱼ)

    std::vector<double> Pq_values;
    Pq_values.reserve(q_values.size());

    for (double q : q_values) {
        double Pq_sum = 0.0;

        // Pre-compute form factors for all atoms at this q
        std::vector<double> form_factors(N);
        for (int i = 0; i < N; ++i) {
            int element = mol.Atom(i).first;
            if (ff_type == "cg_sphere") {
                form_factors[i] = FormFactors::getCGSphereFormFactor(cg_radius, q);
            } else {
                form_factors[i] = FormFactors::getAtomicFormFactor(element, q);
            }
        }

        // Debye double sum
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double fi = form_factors[i];
                double fj = form_factors[j];

                if (i == j) {
                    // Diagonal: sin(0)/0 = 1
                    Pq_sum += fi * fj;
                } else {
                    double rij = distances[i][j];
                    double qr = q * rij;

                    // Sinc function: sin(qr)/(qr)
                    double sinc_qr = (qr < 1e-6) ? 1.0 : std::sin(qr) / qr;
                    Pq_sum += fi * fj * sinc_qr;
                }
            }
        }

        // Normalize by N²
        double Pq = Pq_sum / (N * N);
        Pq_values.push_back(Pq);
    }

    scattering["P_q"] = Pq_values;

    // Guinier analysis: Extract Rg from P(q) at small q
    // ln(P(q)) ≈ -(q²Rg²)/3 for qRg < 1
    if (q_values.size() >= 5) {
        // Use first 5 points for Guinier fit
        double sum_q2 = 0.0, sum_lnP = 0.0, sum_q2_lnP = 0.0, sum_q4 = 0.0;
        int guinier_points = 0;

        for (int i = 0; i < std::min(5, (int)q_values.size()); ++i) {
            if (Pq_values[i] > 0.0) {
                double q2 = q_values[i] * q_values[i];
                double lnP = std::log(Pq_values[i]);

                sum_q2 += q2;
                sum_lnP += lnP;
                sum_q2_lnP += q2 * lnP;
                sum_q4 += q2 * q2;
                guinier_points++;
            }
        }

        if (guinier_points >= 3) {
            // Linear fit: ln(P) = a + b*q²
            // b = -Rg²/3  =>  Rg = √(-3b)
            double b = (guinier_points * sum_q2_lnP - sum_q2 * sum_lnP) /
                      (guinier_points * sum_q4 - sum_q2 * sum_q2);

            if (b < 0.0) {
                double Rg_guinier = std::sqrt(-3.0 * b);
                scattering["guinier_rg"] = Rg_guinier;
            }
        }
    }

    // ===== STRUCTURE FACTOR S(q) - Spherical Averaging =====
    // S(q) = (1/N) |Σᵢ fᵢ(q) exp(i q·rᵢ)|²
    // Average over many q-directions at same |q|

    std::vector<double> Sq_values;
    Sq_values.reserve(q_values.size());

    // Generate Fibonacci sphere sampling for uniform angular coverage
    std::vector<std::array<double, 3>> q_directions;
    q_directions.reserve(angular_samples);

    const double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
    for (int k = 0; k < angular_samples; ++k) {
        double theta = 2.0 * M_PI * k / golden_ratio;
        double phi = std::acos(1.0 - 2.0 * (k + 0.5) / angular_samples);

        q_directions.push_back({
            std::sin(phi) * std::cos(theta),
            std::sin(phi) * std::sin(theta),
            std::cos(phi)
        });
    }

    scattering["angular_samples"] = angular_samples;

    for (double q : q_values) {
        double Sq_avg = 0.0;

        // Average over all angular samples
        for (const auto& q_dir : q_directions) {
            // Complex amplitude: Σᵢ fᵢ(q) exp(i q·rᵢ)
            double real_sum = 0.0;
            double imag_sum = 0.0;

            for (int i = 0; i < N; ++i) {
                int element = mol.Atom(i).first;
                double fi = (ff_type == "cg_sphere")
                    ? FormFactors::getCGSphereFormFactor(cg_radius, q)
                    : FormFactors::getAtomicFormFactor(element, q);

                // Get atom position
                Position pos = mol.Atom(i).second;

                // q·r dot product
                double q_dot_r = q * (q_dir[0] * pos.x() +
                                     q_dir[1] * pos.y() +
                                     q_dir[2] * pos.z());

                // exp(i q·r) = cos(q·r) + i sin(q·r)
                real_sum += fi * std::cos(q_dot_r);
                imag_sum += fi * std::sin(q_dot_r);
            }

            // |Σ exp(i q·r)|² = real² + imag²
            double amplitude_sq = real_sum * real_sum + imag_sum * imag_sum;
            Sq_avg += amplitude_sq / N;  // Normalize by N
        }

        // Average over angular samples
        double Sq = Sq_avg / angular_samples;
        Sq_values.push_back(Sq);
    }

    scattering["S_q"] = Sq_values;

    return scattering;
}

// =====================================================================
// Radial Distribution Function g(r) Implementation - Claude Generated 2025
// =====================================================================

json UnifiedAnalysis::calculatePairDistribution(const Molecule& mol)
{
    json rdf;

    // Get parameters with defaults - Claude Generated 2025
    double r_max = 15.0;
    double bin_width = 0.05;
    bool calc_coordination = false;

    try {
        r_max = m_config.get<double>("rdf_r_max");
        bin_width = m_config.get<double>("rdf_bin_width");
        calc_coordination = m_config.get<bool>("rdf_coordination_shells");
    } catch (...) {
        // Use defaults if parameters not available
    }

    // Setup histogram bins
    int num_bins = static_cast<int>(std::ceil(r_max / bin_width));
    std::vector<double> histogram(num_bins, 0.0);
    std::vector<double> r_values(num_bins);

    // Fill r_values (bin centers)
    for (int i = 0; i < num_bins; ++i) {
        r_values[i] = (i + 0.5) * bin_width;  // Bin center
    }

    rdf["r_max"] = r_max;
    rdf["bin_width"] = bin_width;
    rdf["num_bins"] = num_bins;
    rdf["r_values"] = r_values;

    // Count atom pairs and fill histogram
    const int N = mol.AtomCount();
    int pair_count = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            // Use PBC-aware distance if periodic boundary conditions present
            double dist = mol.hasPBC()
                ? mol.CalculateDistancePBC(i, j)
                : mol.CalculateDistance(i, j);

            // Only count pairs within r_max
            if (dist < r_max) {
                int bin_idx = static_cast<int>(dist / bin_width);
                if (bin_idx >= 0 && bin_idx < num_bins) {
                    histogram[bin_idx] += 2.0;  // Count both (i,j) and (j,i)
                    pair_count += 2;
                }
            }
        }
    }

    // Calculate system volume
    double volume;
    bool has_pbc = mol.hasPBC();

    if (has_pbc) {
        // Volume from cell matrix determinant
        Eigen::Matrix3d cell = mol.getUnitCell();
        volume = std::abs(cell.determinant());
        rdf["volume_source"] = "periodic_cell";
    } else {
        // Estimate from Rout (spherical approximation)
        double Rout = mol.Rout();
        volume = (4.0 / 3.0) * M_PI * Rout * Rout * Rout;
        rdf["volume_source"] = "spherical_estimate";
        rdf["rout"] = Rout;
    }

    rdf["volume"] = volume;
    rdf["has_pbc"] = has_pbc;

    // Normalize histogram to g(r)
    // g(r) = (V / (4πr²Δr N(N-1))) × histogram(r)
    std::vector<double> gr_values(num_bins);
    std::vector<double> coordination_number;

    if (calc_coordination) {
        coordination_number.resize(num_bins, 0.0);
    }

    double number_density = N / volume;
    double running_coord = 0.0;

    for (int i = 0; i < num_bins; ++i) {
        double r = r_values[i];
        double shell_volume = 4.0 * M_PI * r * r * bin_width;

        // Ideal gas pair count for this shell
        double ideal_pairs = number_density * shell_volume * (N - 1);

        // g(r) = actual_pairs / ideal_pairs
        if (ideal_pairs > 0.0) {
            gr_values[i] = histogram[i] / ideal_pairs;
        } else {
            gr_values[i] = 0.0;
        }

        // Running coordination number: n(r) = 4πρ ∫₀ʳ r²g(r)dr
        // Discrete: n(r) += 4πρ r² g(r) Δr
        if (calc_coordination) {
            running_coord += 4.0 * M_PI * number_density * r * r * gr_values[i] * bin_width;
            coordination_number[i] = running_coord;
        }
    }

    rdf["g_r"] = gr_values;

    if (calc_coordination) {
        rdf["coordination_number"] = coordination_number;
    }

    // Find first peak (important structural feature)
    double max_gr = 0.0;
    double peak_position = 0.0;

    for (int i = 0; i < num_bins; ++i) {
        if (gr_values[i] > max_gr) {
            max_gr = gr_values[i];
            peak_position = r_values[i];
        }
    }

    rdf["first_peak_position"] = peak_position;
    rdf["first_peak_height"] = max_gr;

    return rdf;
}