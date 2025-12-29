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

#pragma once

#include <vector>
#include <string>
#include <fstream>

#include "src/core/fileiterator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/parameter_macros.h"
#include "src/core/config_manager.h"

#include "curcumamethod.h"

/* Claude Generated 2025: TrajectoryAnalysis Parameter Registry - replaces static TrajectoryAnalysisJson */
BEGIN_PARAMETER_DEFINITION(trajectoryanalysis)
    PARAM(properties, String, "all", "Properties to calculate: all|basic|geometric|cg.", "Analysis", {})
    PARAM(output_format, String, "human", "Output format: human|json|csv.", "Output", {})
    PARAM(output_file, String, "", "File to save results.", "Output", {})
    PARAM(stride, Int, 1, "Analyze every N-th frame.", "Input", {})
    PARAM(start_frame, Int, 0, "Frame to start analysis from.", "Input", {})
    PARAM(end_frame, Int, -1, "Frame to end analysis at (-1 for end).", "Input", {})
    PARAM(moving_average, Int, 10, "Window size for moving average.", "Analysis", {})
    PARAM(correlation_analysis, Bool, true, "Calculate autocorrelations.", "Analysis", {})
    PARAM(fluctuation_analysis, Bool, true, "Calculate fluctuations and variances.", "Analysis", {})
    PARAM(convergence_analysis, Bool, true, "Analyze convergence of properties.", "Analysis", {})
    PARAM(export_timeseries, Bool, false, "Export raw time series data.", "Output", {})
    PARAM(center_of_mass, Bool, true, "Track center of mass evolution.", "Tracking", {})
    PARAM(gyration_radius, Bool, true, "Track gyration radius evolution.", "Tracking", {})
    PARAM(end_to_end_distance, Bool, true, "Track end-to-end distance (polymers).", "Tracking", {})
    PARAM(recenter_structures, Bool, false, "Center structures at origin before analysis.", "Processing", {})
    PARAM(verbose, Bool, true, "Detailed output.", "Output", {})
    PARAM(verbosity, Int, 2, "Verbosity level (0-3)", "Output", {})
END_PARAMETER_DEFINITION

/*! \brief Trajectory analysis for time-series molecular data - Claude Generated
 *
 * Comprehensive trajectory analysis that works transparently with:
 * - All trajectory formats: XYZ.trj, VTF multi-timestep, SDF trajectories
 * - Atomistic and coarse-grained systems
 * - Any molecular property time series
 *
 * Analysis capabilities:
 * - Time-series statistics: mean, variance, autocorrelation
 * - Convergence analysis and equilibration detection
 * - Moving averages and trend analysis
 * - Property fluctuations and stability assessment
 * - Correlation analysis between different properties
 * - CG-specific analysis: gyration radius, end-to-end distance, Rout
 * - Export to various formats for external analysis
 */
class TrajectoryAnalysis : public CurcumaMethod
{
public:
    /**
     * @brief Constructor with JSON configuration (backward compatible)
     * Claude Generated: Phase 4 - ConfigManager Migration
     */
    TrajectoryAnalysis(const json& controller, bool silent);

    /**
     * @brief Constructor with ConfigManager configuration (new, preferred)
     * Claude Generated: Phase 4 - Native ConfigManager support
     */
    TrajectoryAnalysis(const ConfigManager& config, bool silent);

    ~TrajectoryAnalysis();

    /*! \brief Start trajectory analysis */
    void start() override;

    /*! \brief Set input filename (any supported trajectory format) */
    inline void setFileName(const std::string& filename) { m_filename = filename; }

    /*! \brief Print help for trajectory analysis options */
    void printHelp() const override;

    /*! \brief Get analysis results as JSON */
    const json& getResults() const { return m_results; }

    // CurcumaMethod interface implementation
    nlohmann::json WriteRestartInformation() override { return json{}; }
    bool LoadRestartInformation() override { return true; }
    StringList MethodName() const override { return {"TrajectoryAnalysis"}; }
    void ReadControlFile() override {}
    void LoadControlJson() override {}

private:
    /*! \brief Initialize trajectory analysis */
    bool Initialise() override;

    /*! \brief Analyze single timestep */
    void analyzeTimestep(const Molecule& mol, int timestep);

    /*! \brief Calculate time-series statistics */
    void calculateStatistics();

    /*! \brief Perform convergence analysis */
    void analyzeConvergence();

    /*! \brief Calculate autocorrelation functions */
    void calculateAutocorrelations();

    /*! \brief Calculate moving averages */
    void calculateMovingAverages();

    /*! \brief Analyze fluctuations and stability */
    void analyzeFluctuations();

    /*! \brief Export time series data */
    void exportTimeSeries();

    /*! \brief Output results in requested format */
    void outputResults();

    /*! \brief Output results to file */
    void outputToFile(const std::string& filename);

    /*! \brief Calculate basic molecular properties for timestep */
    json calculateBasicProperties(const Molecule& mol, int timestep);

    /*! \brief Calculate geometric properties for timestep */
    json calculateGeometricProperties(const Molecule& mol, int timestep);

    /*! \brief Calculate CG/polymer properties for timestep */
    json calculateCGProperties(const Molecule& mol, int timestep);

    /*! \brief Detect if molecule has CG or polymer structure */
    bool isCGorPolymer(const Molecule& mol);

    // Configuration
    std::string m_filename;
    json m_config;
    json m_results;

    // Analysis parameters
    std::string m_properties;
    std::string m_output_format;
    std::string m_output_file;
    int m_stride;
    int m_start_frame;
    int m_end_frame;
    int m_moving_average_window;
    bool m_correlation_analysis;
    bool m_fluctuation_analysis;
    bool m_convergence_analysis;
    bool m_export_timeseries;
    bool m_center_of_mass;
    bool m_gyration_radius;
    bool m_end_to_end_distance;
    bool m_recenter_structures; // Claude Generated 2025: Center structures at origin before analysis
    int m_verbosity; // Verbosity level (0-3)

    // Claude Generated 2025: PBC logging state
    bool m_pbc_used = false; // Track if PBC notification has been shown

    // Time series data storage
    std::vector<double> m_time_points;
    std::vector<double> m_gyration_radius_series;
    std::vector<double> m_gyration_radius_mass_weighted_series;
    std::vector<double> m_end_to_end_series;
    std::vector<double> m_rout_series;
    std::vector<std::vector<double>> m_center_of_mass_series; // [timestep][x,y,z]
    std::vector<double> m_total_energy_series;
    std::vector<int> m_fragment_count_series;
    std::vector<double> m_molecular_mass_series;

    // Statistics
    struct TimeSeriesStats {
        double mean = 0.0;
        double variance = 0.0;
        double std_dev = 0.0;
        double min_value = 0.0;
        double max_value = 0.0;
        double range = 0.0;
        std::vector<double> moving_average;
        std::vector<double> autocorrelation;
        double equilibration_time = 0.0;
        bool is_converged = false;
    };

    TimeSeriesStats m_gyration_stats;
    TimeSeriesStats m_end_to_end_stats;
    TimeSeriesStats m_rout_stats;
    TimeSeriesStats m_com_displacement_stats;
    TimeSeriesStats m_energy_stats;

    // Analysis helper functions
    TimeSeriesStats calculateTimeSeriesStats(const std::vector<double>& series, const std::string& name);
    std::vector<double> calculateMovingAverage(const std::vector<double>& series, int window);
    std::vector<double> calculateAutocorrelation(const std::vector<double>& series, int max_lag = 100);
    double calculateEquilibrationTime(const std::vector<double>& series);
    bool assessConvergence(const std::vector<double>& series, double threshold = 0.05);
    double calculateCOMDisplacement(int timestep);
};