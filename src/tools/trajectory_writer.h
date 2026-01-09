/* trajectory_writer.h - Claude Generated 2025
 *
 * Unified trajectory output system for Curcuma
 * Consolidates fragmented output code from analysis.cpp, trajectoryanalysis.cpp, rmsdtraj.cpp
 */

#ifndef TRAJECTORY_WRITER_H
#define TRAJECTORY_WRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include "json.hpp"

// Forward declaration
class TrajectoryStatistics;

using json = nlohmann::json;

class TrajectoryWriter {
public:
    enum class Format {
        HumanTable,
        CSV,
        JSON,
        DAT,
        VTF
    };

    enum class StatType {
        None,
        Cumulative,
        Moving,
        Both
    };

    explicit TrajectoryWriter(const json& config = {});
    ~TrajectoryWriter() = default;

    // Core writing methods
    void writeHumanTable(std::ostream& out, const json& trajectory_data) const;
    void writeCSV(std::ostream& out, const json& trajectory_data) const;
    void writeJSON(std::ostream& out, const json& trajectory_data) const;
    void writeDAT(std::ostream& out, const json& data, const std::string& type) const;
    void writeVTF(std::ostream& out, const json& trajectory_data) const;

    // Convenience methods
    void writeToFile(const std::string& filename, Format format, const json& data) const;
    void writeMultipleFiles(const std::string& base_filename, const json& data) const;

    // Scattering per-frame file output helpers - Claude Generated 2026
    bool writeScatteringPerFrameFiles(const json& timestep_data,
                                     const std::string& output_directory,
                                     const std::string& file_prefix) const;

    // Cross-frame scattering statistics aggregation - Claude Generated 2026
    bool writeScatteringStatistics(const json& trajectory_data,
                                   const std::string& output_directory,
                                   const std::string& file_prefix) const;

    // Configuration methods
    void setDefaultFormat(Format format) { m_default_format = format; }
    void setPrecision(int decimal_places);
    void setColumnWidths(const std::vector<int>& widths) { m_column_widths = widths; }
    void setStatisticsConfig(bool enable_cumulative, bool enable_moving);

    // Static helper methods
    static std::string formatHeader(const std::string& metric, StatType stat_type);
    static std::vector<std::string> parseMetricsList(const std::string& metrics_str);
    static std::string bracketNotation(const std::string& metric, const std::string& stat_suffix);

    // Claude Generated (2025): JSON schema conversion helper for geometry commands
    static json convertToTrajectoryJSON(
        const std::vector<double>& values,
        const std::string& metric_name,
        const TrajectoryStatistics& stats,
        const json& additional_data = {});

    // Claude Generated (2025): JSON schema converter for geometry commands
    static json createTrajectoryJSON(
        const std::vector<double>& values,
        const std::string& metric_name,
        const std::string& unit,
        const TrajectoryStatistics& stats,
        const std::vector<json>& atom_info = {},
        const std::vector<int>& frame_numbers = {});

private:
    // Helper struct for per-q-point statistics - Claude Generated 2026
    struct QPointStatistics {
        std::vector<double> P_q_values;  // Store all frames for median
        std::vector<double> S_q_values;  // Store all frames for median
        double P_q_sum = 0.0;            // For mean calculation
        double S_q_sum = 0.0;
        double P_q_M2 = 0.0;             // For std dev (Welford)
        double S_q_M2 = 0.0;
        int count = 0;                   // Number of frames
    };

    json m_config;
    Format m_default_format;
    int m_precision;
    bool m_enable_cumulative;
    bool m_enable_moving;
    std::vector<int> m_column_widths;

    // Helper methods
    void writeHumanHeader(std::ostream& out, const std::vector<std::string>& enabled_metrics) const;
    void writeCSVHeader(std::ostream& out, const std::vector<std::string>& enabled_metrics) const;
    void writeHumanDataLine(std::ostream& out, const json& timestep, size_t structure_num,
                           const std::vector<std::string>& enabled_metrics) const;
    void writeCSVDataLine(std::ostream& out, const json& timestep, size_t structure_num,
                         const std::vector<std::string>& enabled_metrics) const;
    void writeDATStatistics(std::ostream& out, const json& data) const;

    // Format helpers
    std::string formatValue(double value, int width = 9) const;
    std::string formatValueCSV(double value) const;
    std::string getMetricDisplayName(const std::string& metric) const;
    std::vector<int> getMetricColumnCount(const std::vector<std::string>& metrics) const;

    // Safe JSON handling helpers
    double safeGetDouble(const json& j, double default_value = 0.0) const;
    double safeGetNestedDouble(const json& j, const std::string& key, double default_value = 0.0) const;
    bool isValidNumber(const json& j) const;
    bool isValidNestedNumber(const json& j, const std::string& key) const;

    // Scattering-specific helper functions
    double safeGetScatteringValue(const json& scattering, const std::string& key, double default_value = 0.0) const;
    std::vector<double> extractArrayValues(const json& array, const std::vector<int>& indices) const;
    std::vector<double> extractSampledArrayValues(const json& array, int sampling_rate) const;
    std::string formatArrayForHuman(const std::vector<double>& values, int width = 9, int precision = 6) const;
    std::string formatArrayForCSV(const std::vector<double>& values, int precision = 6) const;

    // Scattering per-frame file output helpers - Claude Generated 2026
    std::vector<int> mapQValuesToIndices(const std::vector<double>& user_q_values,
                                        const json& q_values_array) const;
    std::vector<std::pair<double, double>> extractScatteringDataAtIndices(
                                        const json& scattering_data,
                                        const std::vector<int>& indices,
                                        const std::string& data_type) const;
    std::vector<double> parseQValueString(const std::string& q_value_str) const;

    // Array statistics helpers
    double calculateMean(const json& array) const;
    double calculateMin(const json& array) const;
    double calculateMax(const json& array) const;

    // Statistics helpers
    bool metricAvailable(const json& timestep, const std::string& metric) const;
    json getMetricValue(const json& timestep, const std::string& metric) const;
    json getStatisticsValue(const json& timestep, const std::string& metric, const std::string& stat_type) const;
};

// Progress tracking utility
class ProgressTracker {
public:
    explicit ProgressTracker(int report_interval = 5);
    ~ProgressTracker() = default;

    void setReportInterval(int percent_interval) { m_report_interval = percent_interval; }
    void setMessage(const std::string& message) { m_message = message; }
    void setAutoTiming(bool enable) { m_auto_timing = enable; }
    void setCallback(std::function<void(int, int, const std::string&)> callback) { m_callback = callback; }

    void reportProgress(int current, int total, const std::string& custom_message = "");
    void reset();

    // Timing information
    double getElapsedTime() const;
    double getEstimatedTimeRemaining(int current, int total) const;
    std::string formatTime(double seconds) const;

private:
    int m_report_interval;
    std::string m_message;
    bool m_auto_timing;
    std::function<void(int, int, const std::string&)> m_callback;
    std::vector<bool> m_progress_reported;
    std::chrono::steady_clock::time_point m_start_time;
    std::chrono::steady_clock::time_point m_last_report_time;

    void initializeProgressVector(int total);
    bool shouldReport(int current_percent);
    void internalReport(int current, int total);
};

#endif // TRAJECTORY_WRITER_H