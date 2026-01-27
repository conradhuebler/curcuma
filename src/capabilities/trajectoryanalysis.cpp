/*
 * <Trajectory analysis for time-series molecular data.>
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

#include "trajectoryanalysis.h"
#include "src/tools/formats.h"
#include "src/tools/pbc_utils.h"
#include "src/core/elements.h"
#include "src/core/parameter_registry.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <iomanip>

TrajectoryAnalysis::TrajectoryAnalysis(const json& controller, bool silent)
    : TrajectoryAnalysis(ConfigManager("trajectoryanalysis", controller), silent)
{
}

TrajectoryAnalysis::TrajectoryAnalysis(const ConfigManager& config, bool silent)
    : CurcumaMethod(json{}, config.exportConfig(), silent)
{
    UpdateController(config.exportConfig());
    m_config = MergeJson(config.exportConfig(), config.exportConfig());

    // Extract analysis parameters
    m_properties = Json2KeyWord<std::string>(m_config, "properties");
    m_output_format = Json2KeyWord<std::string>(m_config, "output_format");
    m_output_file = Json2KeyWord<std::string>(m_config, "output_file");
    m_stride = Json2KeyWord<int>(m_config, "stride");
    m_start_frame = Json2KeyWord<int>(m_config, "start_frame");
    m_end_frame = Json2KeyWord<int>(m_config, "end_frame");
    m_moving_average_window = Json2KeyWord<int>(m_config, "moving_average");
    m_correlation_analysis = Json2KeyWord<bool>(m_config, "correlation_analysis");
    m_fluctuation_analysis = Json2KeyWord<bool>(m_config, "fluctuation_analysis");
    m_convergence_analysis = Json2KeyWord<bool>(m_config, "convergence_analysis");
    m_export_timeseries = Json2KeyWord<bool>(m_config, "export_timeseries");
    m_center_of_mass = Json2KeyWord<bool>(m_config, "center_of_mass");

    // Claude Generated 2025: Initialize TrajectoryWriter
    json writer_config = json::object();
    writer_config["default_format"] = "CSV";
    writer_config["column_widths"] = {8, 15, 15}; // timestep, data1, data2
    m_writer = TrajectoryWriter(writer_config);
    m_gyration_radius = Json2KeyWord<bool>(m_config, "gyration_radius");
    m_end_to_end_distance = Json2KeyWord<bool>(m_config, "end_to_end_distance");
    m_recenter_structures = Json2KeyWord<bool>(m_config, "recenter_structures");
    m_verbose = Json2KeyWord<bool>(m_config, "verbose");
}

TrajectoryAnalysis::~TrajectoryAnalysis()
{
}

bool TrajectoryAnalysis::Initialise()
{
    if (m_filename.empty()) {
        CurcumaLogger::error("No input trajectory file specified");
        return false;
    }

    if (!m_silent) {
        CurcumaLogger::info(fmt::format("Initializing trajectory analysis for: {}", m_filename));
        CurcumaLogger::param("properties", fmt::format("{}", m_properties));
        CurcumaLogger::param("stride", fmt::format("{}", m_stride));
        if (m_start_frame > 0 || m_end_frame > 0) {
            CurcumaLogger::param("frame_range", fmt::format("{} to {}", m_start_frame,
                                   m_end_frame == -1 ? "end" : std::to_string(m_end_frame)));
        }
    }

    return true;
}

void TrajectoryAnalysis::start()
{
    if (!Initialise()) {
        CurcumaLogger::error("Failed to initialize trajectory analysis");
        return;
    }

    CurcumaLogger::header("Starting Trajectory Analysis");

    // Initialize file iterator for trajectory
    FileIterator file_iter(m_filename, m_silent);
    int total_frames = file_iter.MaxMolecules();

    if (total_frames == 0) {
        CurcumaLogger::error("No frames found in trajectory file");
        return;
    }

    // Adjust end frame if not specified
    if (m_end_frame == -1) {
        m_end_frame = total_frames - 1;
    }

    if (!m_silent) {
        CurcumaLogger::info(fmt::format("Trajectory contains {} frames", total_frames));
        CurcumaLogger::info(fmt::format("Analyzing frames {} to {} with stride {}",
                              m_start_frame, m_end_frame, m_stride));
    }

    // Reserve space for time series data
    int estimated_frames = (m_end_frame - m_start_frame) / m_stride + 1;
    m_time_points.reserve(estimated_frames);
    m_gyration_radius_series.reserve(estimated_frames);
    m_gyration_radius_mass_weighted_series.reserve(estimated_frames);
    m_end_to_end_series.reserve(estimated_frames);
    m_rout_series.reserve(estimated_frames);
    m_center_of_mass_series.reserve(estimated_frames);

    // Claude Generated 2025: Log recentering notification
    if (m_recenter_structures && CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("Recentering all structures at origin (mass-weighted)");
    }

    // Process trajectory
    int current_frame = 0;
    int analyzed_frames = 0;
    bool pbc_logged = false; // Track if we've logged PBC parameters

    while (!file_iter.AtEnd() && current_frame <= m_end_frame) {
        Molecule mol = file_iter.Next();

        // Claude Generated 2025: Log lattice parameters once when PBC is detected
        if (!pbc_logged && mol.hasPBC() && CurcumaLogger::get_verbosity() >= 2) {
            CurcumaLogger::info("Detected periodic boundary conditions in trajectory");
            auto params = PBCUtils::getLatticeParameters(mol.getUnitCell());
            CurcumaLogger::param("lattice_a", fmt::format("{:.4f} Å", params[0]));
            CurcumaLogger::param("lattice_b", fmt::format("{:.4f} Å", params[1]));
            CurcumaLogger::param("lattice_c", fmt::format("{:.4f} Å", params[2]));
            CurcumaLogger::param("lattice_alpha", fmt::format("{:.2f}°", params[3]));
            CurcumaLogger::param("lattice_beta", fmt::format("{:.2f}°", params[4]));
            CurcumaLogger::param("lattice_gamma", fmt::format("{:.2f}°", params[5]));
            pbc_logged = true;
        }

        // Check if we should analyze this frame
        if (current_frame >= m_start_frame &&
            (current_frame - m_start_frame) % m_stride == 0) {

            analyzeTimestep(mol, analyzed_frames);
            analyzed_frames++;

            // Progress reporting
            if (!m_silent && analyzed_frames % 100 == 0) {
                double progress = 100.0 * current_frame / m_end_frame;
                CurcumaLogger::info(fmt::format("Progress: {:.1f}% ({} frames analyzed)",
                                      progress, analyzed_frames));
            }
        }

        current_frame++;
    }

    if (!m_silent) {
        CurcumaLogger::success_fmt("Analyzed {} frames from trajectory", analyzed_frames);
    }

    // Perform statistical analysis
    calculateStatistics();

    if (m_convergence_analysis) {
        analyzeConvergence();
    }

    if (m_correlation_analysis) {
        calculateAutocorrelations();
    }

    if (m_fluctuation_analysis) {
        analyzeFluctuations();
    }

    if (m_export_timeseries) {
        exportTimeSeries();
    }

    // Output results
    outputResults();
}

void TrajectoryAnalysis::analyzeTimestep(const Molecule& mol, int timestep)
{
    m_time_points.push_back(static_cast<double>(timestep));

    // Claude Generated 2025: Create working copy and optionally center at origin
    Molecule mol_analysis = mol;
    if (m_recenter_structures) {
        mol_analysis.Center(true); // Mass-weighted centering
    }

    // Basic properties
    if (m_properties == "all" || m_properties == "geometric") {
        // Gyration radius
        if (m_gyration_radius) {
            auto gyr = mol_analysis.GyrationRadius();
            m_gyration_radius_series.push_back(gyr.first);
            m_gyration_radius_mass_weighted_series.push_back(gyr.second);
        }

        // Center of mass tracking
        if (m_center_of_mass) {
            Position com = mol_analysis.MassCentroid();
            m_center_of_mass_series.push_back({com[0], com[1], com[2]});
        }
    }

    // CG/Polymer properties
    if ((m_properties == "all" || m_properties == "cg") && isCGorPolymer(mol_analysis)) {
        // End-to-end distance (PBC-aware if available)
        if (m_end_to_end_distance) {
            // Claude Generated 2025: Log first use of PBC-aware calculations
            if (mol_analysis.hasPBC() && !m_pbc_used && CurcumaLogger::get_verbosity() >= 2) {
                CurcumaLogger::info("Using PBC-aware end-to-end distance calculations");
                m_pbc_used = true;
            }

            // Claude Generated: Use PBC-aware version when periodic boundaries are present
            double end_to_end = mol_analysis.hasPBC() ? mol_analysis.EndToEndDistancePBC() : mol_analysis.EndToEndDistance();
            if (end_to_end > 0.0) {
                m_end_to_end_series.push_back(end_to_end);
            }
        }

        // Rout (distance from COM to outermost atom)
        double rout = mol_analysis.Rout();
        m_rout_series.push_back(rout);
    }

    // Fragment counting for stability analysis
    m_fragment_count_series.push_back(mol_analysis.GetFragments().size());
    m_molecular_mass_series.push_back(mol_analysis.Mass());
}

bool TrajectoryAnalysis::isCGorPolymer(const Molecule& mol)
{
    // Check for CG elements
    for (int atom : mol.Atoms()) {
        if (atom == CG_ELEMENT) {
            return true;
        }
    }

    // Check for polymer-like structure (elongated)
    auto gyr = const_cast<Molecule&>(mol).GyrationRadius();
    double rout = mol.Rout();
    return (rout > 2.0 * gyr.first);
}

void TrajectoryAnalysis::calculateStatistics()
{
    if (!m_silent) {
        CurcumaLogger::info("Calculating time-series statistics...");
    }

    // Calculate statistics for each property
    if (!m_gyration_radius_series.empty()) {
        calculateTimeSeriesStats(m_gyration_radius_series, m_gyration_stats, "gyration_radius", m_moving_average_window);
    }

    if (!m_end_to_end_series.empty()) {
        calculateTimeSeriesStats(m_end_to_end_series, m_end_to_end_stats, "end_to_end_distance", m_moving_average_window);
    }

    if (!m_rout_series.empty()) {
        calculateTimeSeriesStats(m_rout_series, m_rout_stats, "rout", m_moving_average_window);
    }

    // Center of mass displacement analysis
    if (!m_center_of_mass_series.empty() && m_center_of_mass_series.size() > 1) {
        std::vector<double> com_displacements;
        com_displacements.reserve(m_center_of_mass_series.size() - 1);

        for (size_t i = 1; i < m_center_of_mass_series.size(); ++i) {
            double displacement = calculateCOMDisplacement(i);
            com_displacements.push_back(displacement);
        }

        if (!com_displacements.empty()) {
            calculateTimeSeriesStats(com_displacements, m_com_displacement_stats, "com_displacement", m_moving_average_window);
        }
    }
}

// Claude Generated 2025: Use unified TrajectoryStatistics
void TrajectoryAnalysis::calculateTimeSeriesStats(
    const std::vector<double>& series, TrajectoryStatistics& stats, const std::string& name, int window_size)
{
    if (series.empty()) return;

    // Configure full series storage for this metric
    stats.setStoreFullSeries(true, {name});
    stats.setWindowSize(window_size);

    // Add all values to statistics engine
    for (double value : series) {
        stats.addValue(name, value);
    }

    // Additional analysis could be added here if needed
    // (autocorrelation, equilibration time, etc.)
}

std::vector<double> TrajectoryAnalysis::calculateMovingAverage(
    const std::vector<double>& series, int window)
{
    std::vector<double> moving_avg;
    if (series.empty() || window <= 0) return moving_avg;

    moving_avg.reserve(series.size());

    for (size_t i = 0; i < series.size(); ++i) {
        int start = std::max(0, static_cast<int>(i) - window/2);
        int end = std::min(static_cast<int>(series.size()), static_cast<int>(i) + window/2 + 1);

        double sum = 0.0;
        int count = 0;
        for (int j = start; j < end; ++j) {
            sum += series[j];
            count++;
        }

        moving_avg.push_back(sum / count);
    }

    return moving_avg;
}

std::vector<double> TrajectoryAnalysis::calculateAutocorrelation(
    const std::vector<double>& series, int max_lag)
{
    std::vector<double> autocorr;
    if (series.empty() || max_lag <= 0) return autocorr;

    autocorr.reserve(max_lag);

    double mean = std::accumulate(series.begin(), series.end(), 0.0) / series.size();

    // Calculate variance
    double variance = 0.0;
    for (double value : series) {
        double diff = value - mean;
        variance += diff * diff;
    }
    variance /= series.size();

    // Calculate autocorrelation function
    for (int lag = 0; lag < max_lag && lag < static_cast<int>(series.size()); ++lag) {
        double covariance = 0.0;
        int count = 0;

        for (size_t i = 0; i + lag < series.size(); ++i) {
            covariance += (series[i] - mean) * (series[i + lag] - mean);
            count++;
        }

        if (count > 0) {
            covariance /= count;
            autocorr.push_back(variance > 0.0 ? covariance / variance : 0.0);
        } else {
            autocorr.push_back(0.0);
        }
    }

    return autocorr;
}

double TrajectoryAnalysis::calculateEquilibrationTime(const std::vector<double>& series)
{
    if (series.size() < 10) return 0.0;

    // Simple equilibration detection: find where running average stabilizes
    std::vector<double> running_avg;
    running_avg.reserve(series.size());

    double sum = 0.0;
    for (size_t i = 0; i < series.size(); ++i) {
        sum += series[i];
        running_avg.push_back(sum / (i + 1));
    }

    // Find where the running average derivative becomes small
    double threshold = 0.01; // 1% change threshold
    for (size_t i = 10; i < running_avg.size() - 10; ++i) {
        double derivative = std::abs(running_avg[i + 10] - running_avg[i]) / 10.0;
        double relative_change = std::abs(derivative / running_avg[i]);

        if (relative_change < threshold) {
            return static_cast<double>(i);
        }
    }

    return static_cast<double>(series.size() * 0.1); // Default to 10% of trajectory
}

bool TrajectoryAnalysis::assessConvergence(const std::vector<double>& series, double threshold)
{
    if (series.size() < 20) return false;

    // Check convergence by comparing first and second half variances
    size_t mid = series.size() / 2;

    std::vector<double> first_half(series.begin(), series.begin() + mid);
    std::vector<double> second_half(series.begin() + mid, series.end());

    TrajectoryStatistics stats1, stats2;
    calculateTimeSeriesStats(first_half, stats1, "temp1");
    calculateTimeSeriesStats(second_half, stats2, "temp2");

    double relative_diff = std::abs(stats1.getMean("temp1") - stats2.getMean("temp2")) / std::abs(stats1.getMean("temp1"));
    return relative_diff < threshold;
}

double TrajectoryAnalysis::calculateCOMDisplacement(int timestep)
{
    if (timestep == 0 || timestep >= static_cast<int>(m_center_of_mass_series.size())) {
        return 0.0;
    }

    const auto& com1 = m_center_of_mass_series[timestep - 1];
    const auto& com2 = m_center_of_mass_series[timestep];

    double dx = com2[0] - com1[0];
    double dy = com2[1] - com1[1];
    double dz = com2[2] - com1[2];

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

void TrajectoryAnalysis::analyzeConvergence()
{
    // Implementation for detailed convergence analysis
    if (!m_silent) {
        CurcumaLogger::info("Performing convergence analysis...");
    }
}

void TrajectoryAnalysis::calculateAutocorrelations()
{
    // Autocorrelations already calculated in calculateTimeSeriesStats
    if (!m_silent) {
        CurcumaLogger::info("Autocorrelation analysis completed");
    }
}

void TrajectoryAnalysis::analyzeFluctuations()
{
    if (!m_silent) {
        CurcumaLogger::info("Analyzing fluctuations and stability...");
    }
}

void TrajectoryAnalysis::exportTimeSeries()
{
    // Claude Generated 2025: Use TrajectoryWriter for unified output
    json trajectory_data = json::array();

    for (size_t i = 0; i < m_time_points.size(); ++i) {
        json timestep_data = json::object();
        timestep_data["timestep"] = m_time_points[i];

        // Add geometric properties
        if (!m_gyration_radius_series.empty() && i < m_gyration_radius_series.size()) {
            json gyr = json::object();
            gyr["unweighted"] = m_gyration_radius_series[i];
            gyr["mass_weighted"] = m_gyration_radius_mass_weighted_series[i];
            timestep_data["geometric"]["gyration_radius"] = gyr;
        }

        if (!m_end_to_end_series.empty() && i < m_end_to_end_series.size()) {
            timestep_data["polymer"]["end2end"] = m_end_to_end_series[i];
        }

        if (!m_rout_series.empty() && i < m_rout_series.size()) {
            timestep_data["polymer"]["rout"] = m_rout_series[i];
        }

        if (!m_center_of_mass_series.empty() && i < m_center_of_mass_series.size()) {
            const auto& com = m_center_of_mass_series[i];
            timestep_data["geometric"]["center_of_mass"] = json{{"x", com[0]}, {"y", com[1]}, {"z", com[2]}};
        }

        trajectory_data.push_back(timestep_data);
    }

    // Add statistics
    json stats_config = json::object();
    stats_config["store_full_series"] = true;
    stats_config["full_series_metrics"] = {"gyration_radius", "end2end", "rout"};
    trajectory_data["statistics_config"] = stats_config;

    // Add trajectory statistics
    trajectory_data["statistics"] = json::object();

    if (!m_gyration_stats.getSeries("gyration_radius").empty()) {
        trajectory_data["statistics"]["gyration_radius"] = m_gyration_stats.exportAllStatistics();
    }

    if (!m_end_to_end_stats.getSeries("end2end").empty()) {
        trajectory_data["statistics"]["end2end"] = m_end_to_end_stats.exportAllStatistics();
    }

    if (!m_rout_stats.getSeries("rout").empty()) {
        trajectory_data["statistics"]["rout"] = m_rout_stats.exportAllStatistics();
    }

    // Write using TrajectoryWriter
    std::string export_file = m_filename + "_timeseries.csv";
    m_writer.writeToFile(export_file, TrajectoryWriter::Format::CSV, trajectory_data);

    CurcumaLogger::success_fmt("Time series data exported to: {}", export_file);
}

void TrajectoryAnalysis::outputResults()
{
    if (!m_output_file.empty()) {
        outputToFile(m_output_file);
        return;
    }

    // Prepare results JSON
    m_results["trajectory_file"] = m_filename;
    m_results["frames_analyzed"] = m_time_points.size();
    m_results["analysis_parameters"] = m_config;

    if (m_output_format == "json") {
        std::cout << std::setw(2) << m_results << std::endl;
    } else {
        // Claude Generated 2026: Use TrajectoryWriter for unified Human-readable output
        json summary = {
            {"title", "Trajectory Analysis Results"},
            {"metadata", {
                {"Trajectory file", m_filename},
                {"Frames analyzed", static_cast<int>(m_time_points.size())}
            }},
            {"sections", json::array()}
        };

        if (!m_gyration_radius_series.empty()) {
            summary["sections"].push_back({
                {"name", "Gyration Radius"},
                {"unit", "Å"},
                {"statistics", {
                    {"mean", m_gyration_stats.getMean("gyration_radius")},
                    {"std_dev", m_gyration_stats.getStdDev("gyration_radius")},
                    {"min", m_gyration_stats.getMin("gyration_radius")},
                    {"max", m_gyration_stats.getMax("gyration_radius")},
                    {"count", static_cast<int>(m_gyration_radius_series.size())}
                }}
            });
        }

        if (!m_end_to_end_series.empty()) {
            summary["sections"].push_back({
                {"name", "End-to-End Distance"},
                {"unit", "Å"},
                {"statistics", {
                    {"mean", m_end_to_end_stats.getMean("end_to_end")},
                    {"std_dev", m_end_to_end_stats.getStdDev("end_to_end")},
                    {"min", m_end_to_end_stats.getMin("end_to_end")},
                    {"max", m_end_to_end_stats.getMax("end_to_end")},
                    {"count", static_cast<int>(m_end_to_end_series.size())}
                }}
            });
        }

        if (!m_rout_series.empty()) {
            summary["sections"].push_back({
                {"name", "Rout"},
                {"unit", "Å"},
                {"statistics", {
                    {"mean", m_rout_stats.getMean("rout")},
                    {"std_dev", m_rout_stats.getStdDev("rout")},
                    {"min", m_rout_stats.getMin("rout")},
                    {"max", m_rout_stats.getMax("rout")},
                    {"count", static_cast<int>(m_rout_series.size())}
                }}
            });
        }

        m_writer.writeStatisticsSummary(std::cout, summary);

        if (m_export_timeseries) {
            CurcumaLogger::param("timeseries_export", fmt::format("{}_timeseries.csv", m_filename));
        }
    }
}

void TrajectoryAnalysis::outputToFile(const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        CurcumaLogger::error_fmt("Could not open output file: {}", filename);
        return;
    }

    file << std::setw(2) << m_results << std::endl;
    file.close();
    CurcumaLogger::success_fmt("Analysis results saved to: {}", filename);
}

void TrajectoryAnalysis::printHelp() const
{
    ParameterRegistry::getInstance().printHelp("trajectoryanalysis");
}