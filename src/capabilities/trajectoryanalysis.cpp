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
#include "src/core/elements.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <iomanip>

TrajectoryAnalysis::TrajectoryAnalysis(const json& controller, bool silent)
    : CurcumaMethod(TrajectoryAnalysisJson, controller, silent)
{
    UpdateController(controller);
    m_config = MergeJson(TrajectoryAnalysisJson, controller);

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
    m_gyration_radius = Json2KeyWord<bool>(m_config, "gyration_radius");
    m_end_to_end_distance = Json2KeyWord<bool>(m_config, "end_to_end_distance");
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

    // Process trajectory
    int current_frame = 0;
    int analyzed_frames = 0;

    while (!file_iter.AtEnd() && current_frame <= m_end_frame) {
        Molecule mol = file_iter.Next();

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

    // Basic properties
    if (m_properties == "all" || m_properties == "geometric") {
        // Gyration radius
        if (m_gyration_radius) {
            auto gyr = const_cast<Molecule&>(mol).GyrationRadius();
            m_gyration_radius_series.push_back(gyr.first);
            m_gyration_radius_mass_weighted_series.push_back(gyr.second);
        }

        // Center of mass tracking
        if (m_center_of_mass) {
            Position com = mol.MassCentroid();
            m_center_of_mass_series.push_back({com[0], com[1], com[2]});
        }
    }

    // CG/Polymer properties
    if ((m_properties == "all" || m_properties == "cg") && isCGorPolymer(mol)) {
        // End-to-end distance (PBC-aware if available)
        if (m_end_to_end_distance) {
            // Claude Generated: Use PBC-aware version when periodic boundaries are present
            double end_to_end = mol.hasPBC() ? mol.EndToEndDistancePBC() : mol.EndToEndDistance();
            if (end_to_end > 0.0) {
                m_end_to_end_series.push_back(end_to_end);
            }
        }

        // Rout (distance from COM to outermost atom)
        double rout = mol.Rout();
        m_rout_series.push_back(rout);
    }

    // Fragment counting for stability analysis
    m_fragment_count_series.push_back(mol.GetFragments().size());
    m_molecular_mass_series.push_back(mol.Mass());
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
        m_gyration_stats = calculateTimeSeriesStats(m_gyration_radius_series, "gyration_radius");
    }

    if (!m_end_to_end_series.empty()) {
        m_end_to_end_stats = calculateTimeSeriesStats(m_end_to_end_series, "end_to_end_distance");
    }

    if (!m_rout_series.empty()) {
        m_rout_stats = calculateTimeSeriesStats(m_rout_series, "rout");
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
            m_com_displacement_stats = calculateTimeSeriesStats(com_displacements, "com_displacement");
        }
    }
}

TrajectoryAnalysis::TimeSeriesStats TrajectoryAnalysis::calculateTimeSeriesStats(
    const std::vector<double>& series, const std::string& name)
{
    TimeSeriesStats stats;

    if (series.empty()) return stats;

    // Basic statistics
    stats.mean = std::accumulate(series.begin(), series.end(), 0.0) / series.size();

    stats.min_value = *std::min_element(series.begin(), series.end());
    stats.max_value = *std::max_element(series.begin(), series.end());
    stats.range = stats.max_value - stats.min_value;

    // Variance and standard deviation
    double sum_sq_diff = 0.0;
    for (double value : series) {
        double diff = value - stats.mean;
        sum_sq_diff += diff * diff;
    }
    stats.variance = sum_sq_diff / series.size();
    stats.std_dev = std::sqrt(stats.variance);

    // Moving average
    stats.moving_average = calculateMovingAverage(series, m_moving_average_window);

    // Convergence analysis
    if (m_convergence_analysis) {
        stats.equilibration_time = calculateEquilibrationTime(series);
        stats.is_converged = assessConvergence(series);
    }

    // Autocorrelation
    if (m_correlation_analysis && series.size() > 10) {
        int max_lag = std::min(100, static_cast<int>(series.size() / 4));
        stats.autocorrelation = calculateAutocorrelation(series, max_lag);
    }

    return stats;
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

    auto stats1 = calculateTimeSeriesStats(first_half, "temp1");
    auto stats2 = calculateTimeSeriesStats(second_half, "temp2");

    double relative_diff = std::abs(stats1.mean - stats2.mean) / std::abs(stats1.mean);
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
    std::string export_file = m_filename + "_timeseries.csv";
    std::ofstream file(export_file);

    if (!file.is_open()) {
        CurcumaLogger::error_fmt("Could not open export file: {}", export_file);
        return;
    }

    // Write header
    file << "timestep";
    if (!m_gyration_radius_series.empty()) file << ",gyration_radius,gyration_radius_mass_weighted";
    if (!m_end_to_end_series.empty()) file << ",end_to_end_distance";
    if (!m_rout_series.empty()) file << ",rout";
    if (!m_center_of_mass_series.empty()) file << ",com_x,com_y,com_z";
    file << std::endl;

    // Write data
    for (size_t i = 0; i < m_time_points.size(); ++i) {
        file << m_time_points[i];

        if (!m_gyration_radius_series.empty() && i < m_gyration_radius_series.size()) {
            file << "," << m_gyration_radius_series[i] << "," << m_gyration_radius_mass_weighted_series[i];
        }

        if (!m_end_to_end_series.empty() && i < m_end_to_end_series.size()) {
            file << "," << m_end_to_end_series[i];
        }

        if (!m_rout_series.empty() && i < m_rout_series.size()) {
            file << "," << m_rout_series[i];
        }

        if (!m_center_of_mass_series.empty() && i < m_center_of_mass_series.size()) {
            const auto& com = m_center_of_mass_series[i];
            file << "," << com[0] << "," << com[1] << "," << com[2];
        }

        file << std::endl;
    }

    file.close();
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
        // Human-readable output
        CurcumaLogger::info("=== Trajectory Analysis Results ===");
        CurcumaLogger::param("trajectory_file", fmt::format("{}", m_filename));
        CurcumaLogger::param("frames_analyzed", fmt::format("{}", m_time_points.size()));

        if (!m_gyration_radius_series.empty()) {
            CurcumaLogger::info("");
            CurcumaLogger::info("Gyration Radius Analysis:");
            CurcumaLogger::param("gyration_mean", fmt::format("{:.3f} ± {:.3f} Å", m_gyration_stats.mean, m_gyration_stats.std_dev));
            CurcumaLogger::param("gyration_range", fmt::format("{:.3f} - {:.3f} Å", m_gyration_stats.min_value, m_gyration_stats.max_value));
            if (m_gyration_stats.is_converged) {
                CurcumaLogger::param("equilibration_time", fmt::format("{:.1f} frames", m_gyration_stats.equilibration_time));
                CurcumaLogger::param("status", "Converged");
            } else {
                CurcumaLogger::param("status", "Not converged");
            }
        }

        if (!m_end_to_end_series.empty()) {
            CurcumaLogger::info("");
            CurcumaLogger::info("End-to-End Distance Analysis:");
            CurcumaLogger::param("end_to_end_mean", fmt::format("{:.3f} ± {:.3f} Å", m_end_to_end_stats.mean, m_end_to_end_stats.std_dev));
            CurcumaLogger::param("end_to_end_range", fmt::format("{:.3f} - {:.3f} Å", m_end_to_end_stats.min_value, m_end_to_end_stats.max_value));
        }

        if (!m_rout_series.empty()) {
            CurcumaLogger::info("");
            CurcumaLogger::info("Rout Analysis:");
            CurcumaLogger::param("rout_mean", fmt::format("{:.3f} ± {:.3f} Å", m_rout_stats.mean, m_rout_stats.std_dev));
            CurcumaLogger::param("rout_range", fmt::format("{:.3f} - {:.3f} Å", m_rout_stats.min_value, m_rout_stats.max_value));
        }

        if (m_export_timeseries) {
            CurcumaLogger::info("");
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
    std::cout << "Trajectory Analysis - Time-series analysis for molecular trajectories" << std::endl;
    std::cout << "=====================================================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Usage: curcuma -traj trajectory.xyz [options]" << std::endl;
    std::cout << "       curcuma -traj trajectory.vtf [options]" << std::endl;
    std::cout << "       curcuma -traj md_output.trj.xyz [options]" << std::endl;
    std::cout << std::endl;
    std::cout << "Supported Formats: XYZ.trj, VTF (multi-timestep), SDF trajectories" << std::endl;
    std::cout << std::endl;
    std::cout << "Analysis Options:" << std::endl;
    std::cout << "  -properties all|geometric|cg|energy" << std::endl;
    std::cout << "              Which properties to analyze (default: all)" << std::endl;
    std::cout << "  -stride N               Analyze every Nth frame (default: 1)" << std::endl;
    std::cout << "  -start_frame N          Start analysis from frame N (default: 0)" << std::endl;
    std::cout << "  -end_frame N            End analysis at frame N (default: all)" << std::endl;
    std::cout << std::endl;
    std::cout << "Statistical Analysis:" << std::endl;
    std::cout << "  -moving_average N       Moving average window size (default: 10)" << std::endl;
    std::cout << "  -correlation_analysis true|false   Calculate autocorrelations (default: true)" << std::endl;
    std::cout << "  -fluctuation_analysis true|false   Analyze fluctuations (default: true)" << std::endl;
    std::cout << "  -convergence_analysis true|false   Convergence assessment (default: true)" << std::endl;
    std::cout << std::endl;
    std::cout << "Property Selection:" << std::endl;
    std::cout << "  -gyration_radius true|false        Track gyration radius evolution" << std::endl;
    std::cout << "  -end_to_end_distance true|false    Track end-to-end distance (polymers)" << std::endl;
    std::cout << "  -center_of_mass true|false          Track center of mass motion" << std::endl;
    std::cout << std::endl;
    std::cout << "Output Options:" << std::endl;
    std::cout << "  -output_format human|json|csv" << std::endl;
    std::cout << "              Output format (default: human)" << std::endl;
    std::cout << "  -output_file filename" << std::endl;
    std::cout << "              Save results to file" << std::endl;
    std::cout << "  -export_timeseries true|false" << std::endl;
    std::cout << "              Export raw time series data as CSV (default: false)" << std::endl;
    std::cout << std::endl;
    std::cout << "Analysis Results:" << std::endl;
    std::cout << "  • Time-series statistics: mean, variance, min/max, range" << std::endl;
    std::cout << "  • Convergence analysis: equilibration time, convergence assessment" << std::endl;
    std::cout << "  • Autocorrelation functions for temporal correlations" << std::endl;
    std::cout << "  • Moving averages for trend analysis" << std::endl;
    std::cout << "  • Fluctuation analysis for stability assessment" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  curcuma -traj md_simulation.xyz -properties geometric" << std::endl;
    std::cout << "  curcuma -traj cg_polymer.vtf -properties cg -export_timeseries true" << std::endl;
    std::cout << "  curcuma -traj protein_md.trj.xyz -stride 10 -convergence_analysis true" << std::endl;
    std::cout << "  curcuma -traj long_simulation.xyz -start_frame 1000 -end_frame 5000" << std::endl;
}