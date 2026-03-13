/* trajectory_writer.cpp - Claude Generated 2025
 *
 * Unified trajectory output system implementation for Curcuma
 * Consolidates output code from analysis.cpp, trajectoryanalysis.cpp, rmsdtraj.cpp
 */

#include "trajectory_writer.h"
#include "../capabilities/trajectory_statistics.h"
#include <iomanip>
#include <algorithm>
#include "../core/curcuma_logger.h"

using namespace nlohmann;

// ---------- TrajectoryWriter Implementation ----------

TrajectoryWriter::TrajectoryWriter(const json& config)
    : m_config(config),
      m_default_format(Format::HumanTable),
      m_precision(3),
      m_enable_cumulative(false),
      m_enable_moving(false),
      m_column_widths({5, 9, 11, 9, 10})  // Default: Structure, Data, Cum-mean, Cum-std, Moving
{
    if (config.contains("default_format")) {
        std::string format_str = config["default_format"];
        if (format_str == "HumanTable") m_default_format = Format::HumanTable;
        else if (format_str == "CSV") m_default_format = Format::CSV;
        else if (format_str == "JSON") m_default_format = Format::JSON;
        else if (format_str == "DAT") m_default_format = Format::DAT;
        else if (format_str == "VTF") m_default_format = Format::VTF;
    }

    if (config.contains("precision")) {
        m_precision = config["precision"];
    }

    if (config.contains("enable_cumulative")) {
        m_enable_cumulative = config["enable_cumulative"];
    }

    if (config.contains("enable_moving")) {
        m_enable_moving = config["enable_moving"];
    }
}

void TrajectoryWriter::setPrecision(int decimal_places)
{
    m_precision = decimal_places;
}

void TrajectoryWriter::setStatisticsConfig(bool enable_cumulative, bool enable_moving)
{
    m_enable_cumulative = enable_cumulative;
    m_enable_moving = enable_moving;
}

void TrajectoryWriter::writeHumanTable(std::ostream& out, const json& trajectory_data) const
{
    // Only process trajectories (multi-structure files)
    if (!trajectory_data.contains("timesteps") || trajectory_data["timesteps"].size() <= 1) {
        return;
    }

    const auto& timesteps = trajectory_data["timesteps"];
    if (timesteps.empty()) {
        return;
    }

    // Get configuration from trajectory_data or fallback
    bool has_stats = trajectory_data.contains("statistics_config");
    std::vector<std::string> enabled_metrics;
    bool enable_cumulative = has_stats ? trajectory_data["statistics_config"]["enable_cumulative"].get<bool>() : m_enable_cumulative;
    bool enable_moving = has_stats ? trajectory_data["statistics_config"]["enable_moving"].get<bool>() : m_enable_moving;

    if (has_stats) {
        enabled_metrics = trajectory_data["statistics_config"]["metrics"].get<std::vector<std::string>>();
    } else {
        // Fallback to default metrics
        enabled_metrics = {"gyration", "rout", "end2end"};
    }

    // Write trajectory header
    const auto& first = timesteps[0];
    if (first.contains("basic")) {
        const auto& basic = first["basic"];
        out << "# Trajectory: " << trajectory_data["filename"].get<std::string>()
            << " (Mass: " << std::fixed << std::setprecision(1)
            << basic["mass"].get<double>() << " amu, Atoms: "
            << basic["atom_count"].get<int>() << ")" << std::endl;
    }

    // Write column header
    writeHumanHeader(out, enabled_metrics);

    // Write data lines
    for (size_t i = 0; i < timesteps.size(); ++i) {
        writeHumanDataLine(out, timesteps[i], i + 1, enabled_metrics);
    }
}

void TrajectoryWriter::writeCSV(std::ostream& out, const json& trajectory_data) const
{
    if (!trajectory_data.contains("timesteps") || trajectory_data["timesteps"].empty()) {
        CurcumaLogger::warn("CSV output only available for trajectory data");
        return;
    }

    const auto& timesteps = trajectory_data["timesteps"];

    // Get configuration
    bool has_stats = trajectory_data.contains("statistics_config");
    std::vector<std::string> enabled_metrics;
    bool enable_cumulative = has_stats ? trajectory_data["statistics_config"]["enable_cumulative"].get<bool>() : m_enable_cumulative;
    bool enable_moving = has_stats ? trajectory_data["statistics_config"]["enable_moving"].get<bool>() : m_enable_moving;

    if (has_stats) {
        enabled_metrics = trajectory_data["statistics_config"]["metrics"].get<std::vector<std::string>>();
    } else {
        enabled_metrics = {"gyration", "rout", "end2end"};
    }

    // Write CSV header
    writeCSVHeader(out, enabled_metrics);

    // Write data lines
    for (size_t i = 0; i < timesteps.size(); ++i) {
        writeCSVDataLine(out, timesteps[i], i + 1, enabled_metrics);
    }
}

void TrajectoryWriter::writeJSON(std::ostream& out, const json& trajectory_data) const
{
    out << std::setw(2) << trajectory_data << std::endl;
}

void TrajectoryWriter::writeDAT(std::ostream& out, const json& data, const std::string& type) const
{
    if (type == "rmsd" || type == "summary") {
        writeDATStatistics(out, data);
    } else {
        // Generic DAT format for raw data
        if (data.contains("timesteps")) {
            for (const auto& timestep : data["timesteps"]) {
                if (timestep.contains("rout")) {
                    out << std::fixed << std::setprecision(6) << timestep["rout"].get<double>() << std::endl;
                } else if (timestep.contains("rmsd")) {
                    out << std::fixed << std::setprecision(6) << timestep["rmsd"].get<double>() << std::endl;
                }
            }
        }
    }
}

// Claude Generated 2026: Statistics summary output for consistent Human-readable format
void TrajectoryWriter::writeStatisticsSummary(std::ostream& out, const json& summary_data) const
{
    // Header/title
    if (summary_data.contains("title")) {
        out << "=== " << summary_data["title"].get<std::string>() << " ===" << std::endl;
        out << std::endl;
    }

    // Metadata section
    if (summary_data.contains("metadata")) {
        for (auto& [key, value] : summary_data["metadata"].items()) {
            out << std::setw(20) << std::left << key << ": ";
            if (value.is_string()) {
                out << value.get<std::string>();
            } else if (value.is_number_integer()) {
                out << value.get<int>();
            } else if (value.is_number_float()) {
                out << std::fixed << std::setprecision(m_precision) << value.get<double>();
            } else {
                out << value.dump();
            }
            out << std::endl;
        }
        out << std::endl;
    }

    // Sections with statistics
    if (summary_data.contains("sections")) {
        for (const auto& section : summary_data["sections"]) {
            // Check if section is enabled (default true)
            if (!section.value("enabled", true)) continue;

            std::string name = section.value("name", "Unknown");
            std::string unit = section.value("unit", "");

            out << name << " Statistics:" << std::endl;

            if (section.contains("statistics")) {
                const auto& stats = section["statistics"];

                // Mean ± StdDev
                if (stats.contains("mean")) {
                    out << "  Mean";
                    if (stats.contains("std_dev")) {
                        out << " ± StdDev:    ";
                        out << std::fixed << std::setprecision(m_precision)
                            << stats["mean"].get<double>() << " ± "
                            << stats["std_dev"].get<double>();
                    } else {
                        out << ":            ";
                        out << std::fixed << std::setprecision(m_precision)
                            << stats["mean"].get<double>();
                    }
                    if (!unit.empty()) out << " " << unit;
                    out << std::endl;
                }

                // Range (Min - Max)
                if (stats.contains("min") && stats.contains("max")) {
                    out << "  Range:           "
                        << std::fixed << std::setprecision(m_precision)
                        << stats["min"].get<double>() << " - "
                        << stats["max"].get<double>();
                    if (!unit.empty()) out << " " << unit;
                    out << std::endl;
                }

                // Median
                if (stats.contains("median")) {
                    out << "  Median:          "
                        << std::fixed << std::setprecision(m_precision)
                        << stats["median"].get<double>();
                    if (!unit.empty()) out << " " << unit;
                    out << std::endl;
                }

                // Count
                if (stats.contains("count")) {
                    out << "  Count:           " << stats["count"].get<int>() << std::endl;
                }
            }

            out << std::endl;
        }
    }
}

void TrajectoryWriter::writeVTF(std::ostream& out, const json& trajectory_data) const
{
    // VTF format for VMD - simplified implementation
    out << "# VTF trajectory file generated by Curcuma" << std::endl;

    if (trajectory_data.contains("timesteps")) {
        const auto& timesteps = trajectory_data["timesteps"];
        for (size_t i = 0; i < timesteps.size(); ++i) {
            const auto& ts = timesteps[i];
            if (ts.contains("coordinates")) {
                const auto& coords = ts["coordinates"];
                out << "timestep" << std::endl;
                out << "atom" << std::endl;
                out << coords.size() << std::endl;
                for (size_t j = 0; j < coords.size(); ++j) {
                    const auto& coord = coords[j];
                    out << "C ";
                    for (int k = 0; k < 3; ++k) {
                        out << std::fixed << std::setprecision(6) << coord[k].get<double>() << " ";
                    }
                    out << std::endl;
                }
            }
        }
    }
}

void TrajectoryWriter::writeToFile(const std::string& filename, Format format, const json& data) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        CurcumaLogger::error_fmt("Could not open output file: {}", filename);
        return;
    }

    switch (format) {
        case Format::HumanTable:
            writeHumanTable(file, data);
            break;
        case Format::CSV:
            writeCSV(file, data);
            break;
        case Format::JSON:
            writeJSON(file, data);
            break;
        case Format::DAT:
            writeDAT(file, data, "default");
            break;
        case Format::VTF:
            writeVTF(file, data);
            break;
    }

    file.close();
}

void TrajectoryWriter::writeMultipleFiles(const std::string& base_filename, const json& data) const
{
    // Write human-readable table
    writeToFile(base_filename + ".txt", Format::HumanTable, data);

    // Write CSV for data analysis
    writeToFile(base_filename + ".csv", Format::CSV, data);

    // Write JSON for structured access
    writeToFile(base_filename + ".json", Format::JSON, data);
}

// ---------- Static Helper Methods ----------

std::string TrajectoryWriter::formatHeader(const std::string& metric, StatType stat_type)
{
    if (metric == "gyration") {
        if (stat_type == StatType::Cumulative) return "<Gyr_u>";
        if (stat_type == StatType::Moving) return "<Gyr_u>_w";
    } else if (metric == "rout") {
        if (stat_type == StatType::Cumulative) return "<Rout>";
        if (stat_type == StatType::Moving) return "<Rout>_w";
    } else if (metric == "end2end") {
        if (stat_type == StatType::Cumulative) return "<End2End>";
        if (stat_type == StatType::Moving) return "<End2End>_w";
    }
    return metric;
}

std::vector<std::string> TrajectoryWriter::parseMetricsList(const std::string& metrics_str)
{
    std::vector<std::string> metrics;
    std::string current;
    for (char c : metrics_str) {
        if (c == ',' || c == ' ') {
            if (!current.empty()) {
                metrics.push_back(current);
                current.clear();
            }
        } else {
            current += c;
        }
    }
    if (!current.empty()) {
        metrics.push_back(current);
    }
    return metrics;
}

std::string TrajectoryWriter::bracketNotation(const std::string& metric, const std::string& stat_suffix)
{
    if (stat_suffix == "mean") return "<" + metric + ">";
    if (stat_suffix == "std") return "σ(" + metric + ")";
    if (stat_suffix == "moving_avg") return "<" + metric + ">_w";
    return metric;
}

// ---------- Private Helper Methods ----------

void TrajectoryWriter::writeHumanHeader(std::ostream& out, const std::vector<std::string>& enabled_metrics) const
{
    out << "# Str";
    for (const auto& metric : enabled_metrics) {
        if (metric == "gyration") {
            out << "  Gyr_u";
            if (m_enable_cumulative) out << "   <Gyr_u> σ(Gyr_u)";
            if (m_enable_moving) out << " <Gyr_u>_w";

            out << "  Gyr_m";
            if (m_enable_cumulative) out << "   <Gyr_m> σ(Gyr_m)";
            if (m_enable_moving) out << " <Gyr_m>_w";
        } else if (metric == "rout") {
            out << "   Rout";
            if (m_enable_cumulative) out << "    <Rout>  σ(Rout)";
            if (m_enable_moving) out << "  <Rout>_w";
        } else if (metric == "end2end") {
            out << "  End2End";
            if (m_enable_cumulative) out << "  <End2End> σ(End2End)";
            if (m_enable_moving) out << " <End2End>_w";
        } else if (metric == "mass") {
            out << "   Mass";
            if (m_enable_cumulative) out << "    <Mass>  σ(Mass)";
            if (m_enable_moving) out << "  <Mass>_w";
        } else if (metric == "com") {
            out << "    COM";
            if (m_enable_cumulative) out << "     <COM>   σ(COM)";
            if (m_enable_moving) out << "   <COM>_w";
        } else if (metric == "scattering") {
            // Enhanced header with additional scattering parameters
            out << " P_q_mean  P_q_min  P_q_max  S_q_mean  S_q_min  S_q_max";
            out << " Rg_Guinier    q_min    q_max";
            // Add selected array values (example: first and last few values)
            out << " P_q[0]   P_q[5]  P_q[10]  S_q[0]   S_q[5]  S_q[10]";
        }
    }
    out << std::endl;
}

void TrajectoryWriter::writeCSVHeader(std::ostream& out, const std::vector<std::string>& enabled_metrics) const
{
    out << "Structure";
    for (const auto& metric : enabled_metrics) {
        if (metric == "gyration") {
            out << ",Gyr_u";
            if (m_enable_cumulative) out << ",<Gyr_u>,σ(Gyr_u)";
            if (m_enable_moving) out << ",<Gyr_u>_w";
            out << ",Gyr_m";
            if (m_enable_cumulative) out << ",<Gyr_m>,σ(Gyr_m)";
            if (m_enable_moving) out << ",<Gyr_m>_w";
        } else if (metric == "rout") {
            out << ",Rout";
            if (m_enable_cumulative) out << ",<Rout>,σ(Rout)";
            if (m_enable_moving) out << ",<Rout>_w";
        } else if (metric == "end2end") {
            out << ",End2End";
            if (m_enable_cumulative) out << ",<End2End>,σ(End2End)";
            if (m_enable_moving) out << ",<End2End>_w";
        } else if (metric == "mass") {
            out << ",Mass";
            if (m_enable_cumulative) out << ",<Mass>,σ(Mass)";
            if (m_enable_moving) out << ",<Mass>_w";
        } else if (metric == "com") {
            out << ",COM";
            if (m_enable_cumulative) out << ",<COM>,σ(COM)";
            if (m_enable_moving) out << ",<COM>_w";
        } else if (metric == "scattering") {
            // Enhanced header with additional scattering parameters
            out << ",P_q_mean,P_q_min,P_q_max,S_q_mean,S_q_min,S_q_max";
            out << ",Rg_Guinier,q_min,q_max";
            // Add selected array values (example: first and last few values)
            out << ",P_q[0],P_q[5],P_q[10],S_q[0],S_q[5],S_q[10]";
        }
    }
    out << std::endl;
}

void TrajectoryWriter::writeHumanDataLine(std::ostream& out, const json& timestep, size_t structure_num,
                                        const std::vector<std::string>& enabled_metrics) const
{
    out << std::setw(5) << structure_num;  // Structure number

    for (const auto& metric : enabled_metrics) {
        if (metric == "gyration") {
            // Unweighted gyration
            if (metricAvailable(timestep, "gyration")) {
                const auto& gyr = getMetricValue(timestep, "gyration");
                if (isValidNestedNumber(gyr, "unweighted")) {
                    out << formatValue(gyr["unweighted"].get<double>());
                } else {
                    out << std::setw(9) << "---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "gyration_unweighted", "mean");
                    auto stat_std = getStatisticsValue(timestep, "gyration_unweighted", "std");
                    if (isValidNumber(stat_mean)) {
                        out << formatValue(stat_mean.get<double>(), 11);
                    } else {
                        out << std::setw(11) << "---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << formatValue(stat_std.get<double>());
                    } else {
                        out << std::setw(9) << "---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "gyration_unweighted", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << formatValue(stat_moving.get<double>(), 10);
                    } else {
                        out << std::setw(10) << "---";
                    }
                }

                // Mass-weighted gyration
                if (isValidNestedNumber(gyr, "mass_weighted")) {
                    out << formatValue(gyr["mass_weighted"].get<double>());
                } else {
                    out << std::setw(9) << "---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "gyration_mass", "mean");
                    auto stat_std = getStatisticsValue(timestep, "gyration_mass", "std");
                    if (isValidNumber(stat_mean)) {
                        out << formatValue(stat_mean.get<double>(), 11);
                    } else {
                        out << std::setw(11) << "---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << formatValue(stat_std.get<double>());
                    } else {
                        out << std::setw(9) << "---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "gyration_mass", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << formatValue(stat_moving.get<double>(), 10);
                    } else {
                        out << std::setw(10) << "---";
                    }
                }
            } else {
                // Skip columns with "---"
                int skip_cols = 2 + (m_enable_cumulative ? 4 : 0) + (m_enable_moving ? 2 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "rout") {
            if (metricAvailable(timestep, "rout")) {
                auto rout_value = getMetricValue(timestep, "rout");
                if (isValidNumber(rout_value)) {
                    out << formatValue(rout_value.get<double>());
                } else {
                    out << std::setw(9) << "---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "rout", "mean");
                    auto stat_std = getStatisticsValue(timestep, "rout", "std");
                    if (isValidNumber(stat_mean)) {
                        out << formatValue(stat_mean.get<double>(), 11);
                    } else {
                        out << std::setw(11) << "---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << formatValue(stat_std.get<double>());
                    } else {
                        out << std::setw(9) << "---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "rout", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << formatValue(stat_moving.get<double>(), 10);
                    } else {
                        out << std::setw(10) << "---";
                    }
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "end2end") {
            if (metricAvailable(timestep, "end2end")) {
                auto end2end_value = getMetricValue(timestep, "end_to_end_distance");
                if (isValidNumber(end2end_value)) {
                    out << formatValue(end2end_value.get<double>());
                } else {
                    out << std::setw(9) << "---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "end_to_end_distance", "mean");
                    auto stat_std = getStatisticsValue(timestep, "end_to_end_distance", "std");
                    if (isValidNumber(stat_mean)) {
                        out << formatValue(stat_mean.get<double>(), 11);
                    } else {
                        out << std::setw(11) << "---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << formatValue(stat_std.get<double>());
                    } else {
                        out << std::setw(9) << "---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "end_to_end_distance", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << formatValue(stat_moving.get<double>(), 10);
                    } else {
                        out << std::setw(10) << "---";
                    }
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "mass") {
            if (metricAvailable(timestep, "mass")) {
                auto mass_value = getMetricValue(timestep, "mass");
                if (isValidNumber(mass_value)) {
                    out << formatValue(mass_value.get<double>());
                } else {
                    out << std::setw(9) << "---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "mass", "mean");
                    auto stat_std = getStatisticsValue(timestep, "mass", "std");
                    if (isValidNumber(stat_mean)) {
                        out << formatValue(stat_mean.get<double>(), 11);
                    } else {
                        out << std::setw(11) << "---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << formatValue(stat_std.get<double>());
                    } else {
                        out << std::setw(9) << "---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "mass", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << formatValue(stat_moving.get<double>(), 10);
                    } else {
                        out << std::setw(10) << "---";
                    }
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "com") {
            if (metricAvailable(timestep, "com")) {
                const auto& com = getMetricValue(timestep, "com");
                // Check if com is an array with at least 3 elements and they're all valid numbers
                if (com.is_array() && com.size() >= 3 &&
                    isValidNumber(com[0]) && isValidNumber(com[1]) && isValidNumber(com[2])) {
                    out << std::setw(7) << std::fixed << std::setprecision(3)
                        << "(" << com[0].get<double>() << ","
                        << com[1].get<double>() << ","
                        << com[2].get<double>() << ")";
                } else {
                    out << std::setw(7) << "---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "com", "mean");
                    auto stat_std = getStatisticsValue(timestep, "com", "std");
                    if (isValidNumber(stat_mean)) {
                        out << formatValue(stat_mean.get<double>(), 11);
                    } else {
                        out << std::setw(11) << "---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << formatValue(stat_std.get<double>());
                    } else {
                        out << std::setw(9) << "---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "com", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << formatValue(stat_moving.get<double>(), 10);
                    } else {
                        out << std::setw(10) << "---";
                    }
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "scattering") {
            if (metricAvailable(timestep, "scattering")) {
                const auto& scattering = getMetricValue(timestep, "scattering");
                if (scattering.contains("P_q") && scattering.contains("S_q") &&
                    scattering["P_q"].is_array() && scattering["S_q"].is_array()) {
                    const auto& P_q = scattering["P_q"];
                    const auto& S_q = scattering["S_q"];

                    // Additional safety checks
                    if (!P_q.empty() && !S_q.empty()) {
                        // Original statistical values
                        out << formatValue(calculateMean(P_q))
                            << formatValue(calculateMin(P_q))
                            << formatValue(calculateMax(P_q))
                            << formatValue(calculateMean(S_q))
                            << formatValue(calculateMin(S_q))
                            << formatValue(calculateMax(S_q));

                        // Additional scattering parameters
                        double guinier_rg = safeGetScatteringValue(scattering, "guinier_rg", 0.0);
                        double q_min = safeGetScatteringValue(scattering, "q_min", 0.0);
                        double q_max = safeGetScatteringValue(scattering, "q_max", 0.0);

                        out << formatValue(guinier_rg)
                            << formatValue(q_min)
                            << formatValue(q_max);

                        // Check if q_values array is available
                        if (scattering.contains("q_values") && scattering["q_values"].is_array()) {
                            // We have the q_values array available for mapping
                        }

                        // Selected array values (first few values as example)
                        std::vector<int> selected_indices = {0, 5, 10};
                        auto p_q_values = extractArrayValues(P_q, selected_indices);
                        auto s_q_values = extractArrayValues(S_q, selected_indices);

                        for (double val : p_q_values) {
                            out << formatValue(val);
                        }
                        for (double val : s_q_values) {
                            out << formatValue(val);
                        }
                    } else {
                        // Skip columns with "---" for all enhanced parameters
                        for (int j = 0; j < 15; ++j)  // 6 statistical + 3 additional + 6 array values
                            out << std::setw(9) << "---";
                    }
                } else {
                    // Skip columns with "---" for all enhanced parameters
                    for (int j = 0; j < 15; ++j)  // 6 statistical + 3 additional + 6 array values
                        out << std::setw(9) << "---";
                }
            } else {
                // Skip columns with "---" for all enhanced parameters
                for (int j = 0; j < 15; ++j)  // 6 statistical + 3 additional + 6 array values
                    out << std::setw(9) << "---";
            }
        }
    }
    out << std::endl;
}

void TrajectoryWriter::writeCSVDataLine(std::ostream& out, const json& timestep, size_t structure_num,
                                       const std::vector<std::string>& enabled_metrics) const
{
    out << std::fixed << std::setprecision(6) << structure_num;

    for (const auto& metric : enabled_metrics) {
        if (metric == "gyration") {
            // Unweighted gyration
            if (metricAvailable(timestep, "gyration")) {
                const auto& gyr = getMetricValue(timestep, "gyration");
                if (isValidNestedNumber(gyr, "unweighted")) {
                    out << "," << formatValueCSV(gyr["unweighted"].get<double>());
                } else {
                    out << ",---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "gyration_unweighted", "mean");
                    auto stat_std = getStatisticsValue(timestep, "gyration_unweighted", "std");
                    if (isValidNumber(stat_mean)) {
                        out << "," << formatValueCSV(stat_mean.get<double>());
                    } else {
                        out << ",---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << "," << formatValueCSV(stat_std.get<double>());
                    } else {
                        out << ",---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "gyration_unweighted", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << "," << formatValueCSV(stat_moving.get<double>());
                    } else {
                        out << ",---";
                    }
                }

                // Mass-weighted gyration
                if (isValidNestedNumber(gyr, "mass_weighted")) {
                    out << "," << formatValueCSV(gyr["mass_weighted"].get<double>());
                } else {
                    out << ",---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "gyration_mass", "mean");
                    auto stat_std = getStatisticsValue(timestep, "gyration_mass", "std");
                    if (isValidNumber(stat_mean)) {
                        out << "," << formatValueCSV(stat_mean.get<double>());
                    } else {
                        out << ",---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << "," << formatValueCSV(stat_std.get<double>());
                    } else {
                        out << ",---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "gyration_mass", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << "," << formatValueCSV(stat_moving.get<double>());
                    } else {
                        out << ",---";
                    }
                }
            } else {
                // Skip columns with empty commas
                int skip_cols = 2 + (m_enable_cumulative ? 4 : 0) + (m_enable_moving ? 2 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "rout") {
            if (metricAvailable(timestep, "rout")) {
                auto rout_value = getMetricValue(timestep, "rout");
                if (isValidNumber(rout_value)) {
                    out << "," << formatValueCSV(rout_value.get<double>());
                } else {
                    out << ",---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "rout", "mean");
                    auto stat_std = getStatisticsValue(timestep, "rout", "std");
                    if (isValidNumber(stat_mean)) {
                        out << "," << formatValueCSV(stat_mean.get<double>());
                    } else {
                        out << ",---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << "," << formatValueCSV(stat_std.get<double>());
                    } else {
                        out << ",---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "rout", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << "," << formatValueCSV(stat_moving.get<double>());
                    } else {
                        out << ",---";
                    }
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "end2end") {
            if (metricAvailable(timestep, "end2end")) {
                auto end2end_value = getMetricValue(timestep, "end_to_end_distance");
                if (isValidNumber(end2end_value)) {
                    out << "," << formatValueCSV(end2end_value.get<double>());
                } else {
                    out << ",---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "end_to_end_distance", "mean");
                    auto stat_std = getStatisticsValue(timestep, "end_to_end_distance", "std");
                    if (isValidNumber(stat_mean)) {
                        out << "," << formatValueCSV(stat_mean.get<double>());
                    } else {
                        out << ",---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << "," << formatValueCSV(stat_std.get<double>());
                    } else {
                        out << ",---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "end_to_end_distance", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << "," << formatValueCSV(stat_moving.get<double>());
                    } else {
                        out << ",---";
                    }
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "mass") {
            if (metricAvailable(timestep, "mass")) {
                auto mass_value = getMetricValue(timestep, "mass");
                if (isValidNumber(mass_value)) {
                    out << "," << formatValueCSV(mass_value.get<double>());
                } else {
                    out << ",---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "mass", "mean");
                    auto stat_std = getStatisticsValue(timestep, "mass", "std");
                    if (isValidNumber(stat_mean)) {
                        out << "," << formatValueCSV(stat_mean.get<double>());
                    } else {
                        out << ",---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << "," << formatValueCSV(stat_std.get<double>());
                    } else {
                        out << ",---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "mass", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << "," << formatValueCSV(stat_moving.get<double>());
                    } else {
                        out << ",---";
                    }
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "com") {
            if (metricAvailable(timestep, "com")) {
                const auto& com = getMetricValue(timestep, "com");
                // Check if com is an array with at least 3 elements and they're all valid numbers
                if (com.is_array() && com.size() >= 3 &&
                    isValidNumber(com[0]) && isValidNumber(com[1]) && isValidNumber(com[2])) {
                    out << std::fixed << std::setprecision(3)
                        << ",(" << com[0].get<double>() << ","
                        << com[1].get<double>() << ","
                        << com[2].get<double>() << ")";
                } else {
                    out << ",---";
                }

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    auto stat_mean = getStatisticsValue(timestep, "com", "mean");
                    auto stat_std = getStatisticsValue(timestep, "com", "std");
                    if (isValidNumber(stat_mean)) {
                        out << "," << formatValueCSV(stat_mean.get<double>());
                    } else {
                        out << ",---";
                    }
                    if (isValidNumber(stat_std)) {
                        out << "," << formatValueCSV(stat_std.get<double>());
                    } else {
                        out << ",---";
                    }
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    auto stat_moving = getStatisticsValue(timestep, "com", "moving_avg");
                    if (isValidNumber(stat_moving)) {
                        out << "," << formatValueCSV(stat_moving.get<double>());
                    } else {
                        out << ",---";
                    }
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "scattering") {
            if (metricAvailable(timestep, "scattering")) {
                const auto& scattering = getMetricValue(timestep, "scattering");
                if (scattering.contains("P_q") && scattering.contains("S_q") &&
                    scattering["P_q"].is_array() && scattering["S_q"].is_array()) {
                    const auto& P_q = scattering["P_q"];
                    const auto& S_q = scattering["S_q"];

                    // Additional safety checks
                    if (!P_q.empty() && !S_q.empty()) {
                        // Original statistical values
                        out << "," << formatValueCSV(calculateMean(P_q))
                            << "," << formatValueCSV(calculateMin(P_q))
                            << "," << formatValueCSV(calculateMax(P_q))
                            << "," << formatValueCSV(calculateMean(S_q))
                            << "," << formatValueCSV(calculateMin(S_q))
                            << "," << formatValueCSV(calculateMax(S_q));

                        // Additional scattering parameters
                        double guinier_rg = safeGetScatteringValue(scattering, "guinier_rg", 0.0);
                        double q_min = safeGetScatteringValue(scattering, "q_min", 0.0);
                        double q_max = safeGetScatteringValue(scattering, "q_max", 0.0);

                        out << "," << formatValueCSV(guinier_rg)
                            << "," << formatValueCSV(q_min)
                            << "," << formatValueCSV(q_max);

                        // Check if q_values array is available
                        if (scattering.contains("q_values") && scattering["q_values"].is_array()) {
                            // We have the q_values array available for mapping
                        }

                        // Selected array values (first few values as example)
                        std::vector<int> selected_indices = {0, 5, 10};
                        auto p_q_values = extractArrayValues(P_q, selected_indices);
                        auto s_q_values = extractArrayValues(S_q, selected_indices);

                        for (double val : p_q_values) {
                            out << "," << formatValueCSV(val);
                        }
                        for (double val : s_q_values) {
                            out << "," << formatValueCSV(val);
                        }
                    } else {
                        // Skip columns with empty commas for all enhanced parameters
                        for (int j = 0; j < 15; ++j)  // 6 statistical + 3 additional + 6 array values
                            out << ",";
                    }
                } else {
                    // Skip columns with empty commas for all enhanced parameters
                    for (int j = 0; j < 15; ++j)  // 6 statistical + 3 additional + 6 array values
                        out << ",";
                }
            } else {
                // Skip columns with empty commas for all enhanced parameters
                for (int j = 0; j < 15; ++j)  // 6 statistical + 3 additional + 6 array values
                    out << ",";
            }
        }
    }
    out << std::endl;
}

void TrajectoryWriter::writeDATStatistics(std::ostream& out, const json& data) const
{
    // Write RMSD-style statistics summary with #Mean, #Median, #StdDev, #Shannon
    if (data.contains("statistics")) {
        const auto& stats = data["statistics"];

        // Try to extract RMSD statistics
        std::string rmsd_key = "rmsd";
        if (!stats.contains(rmsd_key) && stats.contains("values")) {
            rmsd_key = "values";  // Fallback key
        }

        if (stats.contains("rmsd") && stats["rmsd"].contains("summary")) {
            const auto& summary = stats["rmsd"]["summary"];
            out << "#Mean\\t" << summary["mean"] << "\\t" << summary["energy_mean"] << std::endl;
            out << "#Median\\t" << summary["median"] << "\\t" << summary["energy_median"] << std::endl;
            out << "#StdDev\\t" << summary["std"] << "\\t" << summary["energy_std"] << std::endl;
            out << "#Shannon\\t" << summary["shannon"] << "\\t" << summary["energy_shannon"] << std::endl;
        }
    }
}

// ---------- Format Helper Methods ----------

std::string TrajectoryWriter::formatValue(double value, int width) const
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(m_precision) << std::setw(width) << value;
    return oss.str();
}

// Helper function to safely extract double from JSON
double TrajectoryWriter::safeGetDouble(const json& j, double default_value) const
{
    if (j.is_number() && !j.is_null()) {
        return j.get<double>();
    }
    return default_value;
}

// Helper function to safely extract double from nested JSON
double TrajectoryWriter::safeGetNestedDouble(const json& j, const std::string& key, double default_value) const
{
    if (j.contains(key) && j[key].is_number() && !j[key].is_null()) {
        return j[key].get<double>();
    }
    return default_value;
}

// Helper function to check if JSON value is valid for formatting
bool TrajectoryWriter::isValidNumber(const json& j) const
{
    return j.is_number() && !j.is_null();
}

// Helper function to check if nested JSON value is valid for formatting
bool TrajectoryWriter::isValidNestedNumber(const json& j, const std::string& key) const
{
    return j.contains(key) && j[key].is_number() && !j[key].is_null();
}

// Scattering-specific helper functions

// Helper function to safely extract double from scattering JSON
double TrajectoryWriter::safeGetScatteringValue(const json& scattering, const std::string& key, double default_value) const
{
    if (scattering.contains(key) && scattering[key].is_number() && !scattering[key].is_null()) {
        return scattering[key].get<double>();
    }
    return default_value;
}

// Helper function to extract values at specific indices from JSON array
std::vector<double> TrajectoryWriter::extractArrayValues(const json& array, const std::vector<int>& indices) const
{
    std::vector<double> result;
    if (!array.is_array()) {
        return result;
    }

    for (int index : indices) {
        if (index >= 0 && index < static_cast<int>(array.size()) &&
            array[index].is_number() && !array[index].is_null()) {
            result.push_back(array[index].get<double>());
        } else {
            result.push_back(0.0); // Default value for invalid indices
        }
    }
    return result;
}

std::string TrajectoryWriter::formatValueCSV(double value) const
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << value;
    return oss.str();
}

double TrajectoryWriter::calculateMean(const json& array) const
{
    if (!array.is_array() || array.empty()) {
        return 0.0;
    }

    double sum = 0.0;
    int count = 0;
    for (const auto& value : array) {
        if (value.is_number() && !value.is_null()) {
            sum += value.get<double>();
            count++;
        }
    }
    return count > 0 ? sum / count : 0.0;
}

double TrajectoryWriter::calculateMin(const json& array) const
{
    if (!array.is_array() || array.empty()) {
        return 0.0;
    }

    double min_val = std::numeric_limits<double>::max();
    bool found_valid = false;
    for (const auto& value : array) {
        if (value.is_number() && !value.is_null()) {
            min_val = std::min(min_val, value.get<double>());
            found_valid = true;
        }
    }
    return found_valid ? min_val : 0.0;
}

double TrajectoryWriter::calculateMax(const json& array) const
{
    if (!array.is_array() || array.empty()) {
        return 0.0;
    }

    double max_val = std::numeric_limits<double>::lowest();
    bool found_valid = false;
    for (const auto& value : array) {
        if (value.is_number() && !value.is_null()) {
            max_val = std::max(max_val, value.get<double>());
            found_valid = true;
        }
    }
    return found_valid ? max_val : 0.0;
}

bool TrajectoryWriter::metricAvailable(const json& timestep, const std::string& metric) const
{
    if (metric == "gyration") {
        return timestep.contains("geometric") &&
               timestep["geometric"].contains("gyration_radius");
    } else if (metric == "rout") {
        return timestep.contains("polymer") &&
               timestep["polymer"].contains("rout");
    } else if (metric == "end2end") {
        return timestep.contains("polymer") &&
               timestep["polymer"].contains("end_to_end_distance");
    } else if (metric == "mass") {
        return timestep.contains("basic") &&
               timestep["basic"].contains("mass");
    } else if (metric == "com") {
        return timestep.contains("geometric") &&
               timestep["geometric"].contains("center_of_mass");
    } else if (metric == "scattering") {
        return timestep.contains("scattering") &&
               timestep["scattering"].contains("P_q") &&
               timestep["scattering"].contains("S_q");
    }
    return false;
}

json TrajectoryWriter::getMetricValue(const json& timestep, const std::string& metric) const
{
    if (metric == "gyration") {
        return timestep["geometric"]["gyration_radius"];
    } else if (metric == "rout") {
        return timestep["polymer"]["rout"];
    } else if (metric == "end2end") {
        return timestep["polymer"]["end_to_end_distance"];
    } else if (metric == "mass") {
        return timestep["basic"]["mass"];
    } else if (metric == "com") {
        return timestep["geometric"]["center_of_mass"];
    } else if (metric == "scattering") {
        return timestep["scattering"];
    }
    return nullptr;
}

json TrajectoryWriter::getStatisticsValue(const json& timestep, const std::string& metric, const std::string& stat_type) const
{
    if (timestep.contains("statistics") && timestep["statistics"].contains(metric)) {
        const auto& stats = timestep["statistics"][metric];
        if (stat_type == "mean") return stats["mean"];
        if (stat_type == "std") return stats["std"];
        if (stat_type == "moving_avg") return stats["moving_avg"];
        if (stat_type == "min") return stats["min"];
        if (stat_type == "max") return stats["max"];
        if (stat_type == "median") return stats["median"];
    }
    return nullptr;
}

// ---------- ProgressTracker Implementation - Simplified KISS version (Claude Generated 2026) ----------

void ProgressTracker::report(int current, int total, const std::string& label)
{
    int percent = (total > 0) ? (current * 100 / total) : 0;
    if (percent != m_last_percent) {
        fmt::print("\r{}: {}%", label, percent);
        std::fflush(stdout);
        m_last_percent = percent;
    }
}

// ---------- JSON Schema Converter Helper Functions ----------
// Claude Generated 2025 - Helper functions for geometry command refactoring

json TrajectoryWriter::convertToTrajectoryJSON(
    const std::vector<double>& values,
    const std::string& metric_name,
    const TrajectoryStatistics& stats,
    const json& additional_data)
{
    json trajectory_data = json::object();

    // Build timesteps array
    json timesteps = json::array();
    for (size_t i = 0; i < values.size(); ++i) {
        json timestep = json::object();
        timestep["step"] = i + 1;
        timestep[metric_name] = values[i];
        timesteps.push_back(timestep);
    }

    trajectory_data["timesteps"] = timesteps;

    // Add statistics configuration
    json config = json::object();
    config["metric"] = metric_name;
    config["statistics"] = {
        {"mean", stats.getMean(metric_name)},
        {"std_dev", stats.getStdDev(metric_name)},
        {"min_value", stats.getMin(metric_name)},
        {"max_value", stats.getMax(metric_name)},
        {"count", static_cast<int>(values.size())}
    };

    // Merge additional data if provided
    if (!additional_data.is_null() && additional_data.is_object()) {
        trajectory_data.merge_patch(additional_data);
    }

    trajectory_data["statistics_config"] = config;

    return trajectory_data;
}

json TrajectoryWriter::createTrajectoryJSON(
    const std::vector<double>& values,
    const std::string& metric_name,
    const std::string& unit,
    const TrajectoryStatistics& stats,
    const std::vector<json>& atom_info,
    const std::vector<int>& frame_numbers)
{
    json trajectory_data = json::object();

    // Build timesteps array with frame numbers and atom info if available
    json timesteps = json::array();
    for (size_t i = 0; i < values.size(); ++i) {
        json timestep = json::object();
        timestep["step"] = frame_numbers.empty() ? (i + 1) : frame_numbers[i];
        timestep[metric_name] = values[i];

        // Add unit if specified
        if (!unit.empty()) {
            timestep[metric_name + "_unit"] = unit;
        }

        // Add atom info if available
        if (i < atom_info.size() && !atom_info[i].is_null()) {
            timestep.merge_patch(atom_info[i]);
        }

        timesteps.push_back(timestep);
    }

    trajectory_data["timesteps"] = timesteps;

    // Add metadata for TrajectoryWriter compatibility
    trajectory_data["total_timesteps"] = static_cast<int>(values.size());

    // Add statistics configuration with unit
    json config = json::object();
    config["metric"] = metric_name;
    config["unit"] = unit;
    config["statistics"] = {
        {"mean", stats.getMean(metric_name)},
        {"std_dev", stats.getStdDev(metric_name)},
        {"min_value", stats.getMin(metric_name)},
        {"max_value", stats.getMax(metric_name)},
        {"count", static_cast<int>(values.size())}
    };

    // Add TrajectoryWriter expected configuration fields
    config["enable_cumulative"] = true;
    config["enable_moving"] = false;

    trajectory_data["statistics_config"] = config;
    trajectory_data["metrics"] = json::array({metric_name});
    trajectory_data["moving_window"] = 10;

    return trajectory_data;
}
// ---------- Scattering Per-Frame File Output Helpers - Claude Generated 2026 ----------

// Helper function to parse comma-separated q-values from string
std::vector<double> TrajectoryWriter::parseQValueString(const std::string& q_value_str) const
{
    std::vector<double> q_values;
    if (q_value_str.empty()) {
        return q_values;
    }

    std::stringstream ss(q_value_str);
    std::string token;

    while (std::getline(ss, token, ',')) {
        // Trim whitespace
        token.erase(0, token.find_first_not_of(" 	"));
        token.erase(token.find_last_not_of(" 	") + 1);

        try {
            double q_val = std::stod(token);
            q_values.push_back(q_val);
        } catch (const std::exception& e) {
            // Skip invalid values
            continue;
        }
    }

    return q_values;
}

// Helper function to map user q-values to array indices
std::vector<int> TrajectoryWriter::mapQValuesToIndices(const std::vector<double>& user_q_values,
                                                       const json& q_values_array) const
{
    std::vector<int> indices;
    if (!q_values_array.is_array() || q_values_array.empty()) {
        return indices;
    }

    // Convert JSON array to vector for easier access
    std::vector<double> available_q_values;
    for (const auto& q_val : q_values_array) {
        if (q_val.is_number() && !q_val.is_null()) {
            available_q_values.push_back(q_val.get<double>());
        }
    }

    if (available_q_values.empty()) {
        return indices;
    }

    // For each user q-value, find the closest matching index
    for (double user_q : user_q_values) {
        int closest_index = 0;
        double min_diff = std::abs(user_q - available_q_values[0]);

        for (size_t i = 1; i < available_q_values.size(); ++i) {
            double diff = std::abs(user_q - available_q_values[i]);
            if (diff < min_diff) {
                min_diff = diff;
                closest_index = static_cast<int>(i);
            }
        }

        indices.push_back(closest_index);
    }

    return indices;
}

// Helper function to extract scattering data at specific indices
std::vector<std::pair<double, double>> TrajectoryWriter::extractScatteringDataAtIndices(
    const json& scattering_data,
    const std::vector<int>& indices,
    const std::string& data_type) const
{
    std::vector<std::pair<double, double>> result;
    if (!scattering_data.contains("q_values") || !scattering_data.contains(data_type) ||
        !scattering_data["q_values"].is_array() || !scattering_data[data_type].is_array()) {
        return result;
    }

    const auto& q_values = scattering_data["q_values"];
    const auto& data_values = scattering_data[data_type];

    for (int index : indices) {
        if (index >= 0 && index < static_cast<int>(q_values.size()) &&
            index < static_cast<int>(data_values.size()) &&
            q_values[index].is_number() && !q_values[index].is_null() &&
            data_values[index].is_number() && !data_values[index].is_null()) {

            double q_val = q_values[index].get<double>();
            double data_val = data_values[index].get<double>();
            result.emplace_back(q_val, data_val);
        }
    }

    return result;
}

// Helper function to write per-frame scattering files
bool TrajectoryWriter::writeScatteringPerFrameFiles(const json& timestep_data,
                                                   const std::string& output_directory,
                                                   const std::string& file_prefix) const
{
    if (!timestep_data.contains("timesteps") || !timestep_data["timesteps"].is_array()) {
        return false;
    }

    const auto& timesteps = timestep_data["timesteps"];
    bool success = true;

    // Check if we have user-specified q-values
    std::string user_q_values_str = "";
    if (timestep_data.contains("config") && timestep_data["config"].contains("scattering_q_values")) {
        user_q_values_str = timestep_data["config"]["scattering_q_values"].get<std::string>();
    }

    std::vector<double> user_q_values = parseQValueString(user_q_values_str);
    bool use_user_q_values = !user_q_values.empty();

    for (size_t i = 0; i < timesteps.size(); ++i) {
        const auto& timestep = timesteps[i];

        // Check if this timestep has scattering data
        if (!timestep.contains("scattering") || !timestep["scattering"].is_object()) {
            continue;
        }

        const auto& scattering = timestep["scattering"];

        // Check if we have the required data
        if (!scattering.contains("q_values") || !scattering.contains("P_q") || !scattering.contains("S_q") ||
            !scattering["q_values"].is_array() || !scattering["P_q"].is_array() || !scattering["S_q"].is_array()) {
            continue;
        }

        // Generate frame number with leading zeros
        std::stringstream frame_ss;
        frame_ss << std::setfill('0') << std::setw(3) << (i + 1);
        std::string frame_num = frame_ss.str();

        // Get q-values and data arrays
        const auto& q_values_array = scattering["q_values"];
        const auto& pq_values_array = scattering["P_q"];
        const auto& sq_values_array = scattering["S_q"];

        // If user specified q-values, map them to indices
        std::vector<std::pair<double, double>> pq_data, sq_data;
        if (use_user_q_values) {
            // Map user q-values to indices in the array
            std::vector<int> indices = mapQValuesToIndices(user_q_values, q_values_array);

            // Extract P(q) and S(q) data at those indices
            pq_data = extractScatteringDataAtIndices(scattering, indices, "P_q");
            sq_data = extractScatteringDataAtIndices(scattering, indices, "S_q");
        } else {
            // Use all q-values
            for (size_t j = 0; j < q_values_array.size() && j < pq_values_array.size() && j < sq_values_array.size(); ++j) {
                if (q_values_array[j].is_number() && !q_values_array[j].is_null() &&
                    pq_values_array[j].is_number() && !pq_values_array[j].is_null() &&
                    sq_values_array[j].is_number() && !sq_values_array[j].is_null()) {
                    double q_val = q_values_array[j].get<double>();
                    pq_data.emplace_back(q_val, pq_values_array[j].get<double>());
                    sq_data.emplace_back(q_val, sq_values_array[j].get<double>());
                }
            }
        }

        // Create P(q) file
        std::string pq_filename = output_directory + "/" + file_prefix + "_" + frame_num + "_Pq.csv";
        std::ofstream pq_file(pq_filename);
        if (pq_file.is_open()) {
            pq_file << "# q (Å⁻¹),P(q)\n";
            for (const auto& data_point : pq_data) {
                pq_file << std::fixed << std::setprecision(6)
                       << data_point.first << ","
                       << data_point.second << "\n";
            }
            pq_file.close();
        } else {
            success = false;
        }

        // Create S(q) file
        std::string sq_filename = output_directory + "/" + file_prefix + "_" + frame_num + "_Sq.csv";
        std::ofstream sq_file(sq_filename);
        if (sq_file.is_open()) {
            sq_file << "# q (Å⁻¹),S(q)\n";
            for (const auto& data_point : sq_data) {
                sq_file << std::fixed << std::setprecision(6)
                       << data_point.first << ","
                       << data_point.second << "\n";
            }
            sq_file.close();
        } else {
            success = false;
        }
    }

    return success;
}

// Cross-frame scattering statistics aggregation - Claude Generated 2026
bool TrajectoryWriter::writeScatteringStatistics(const json& trajectory_data,
                                                const std::string& output_directory,
                                                const std::string& file_prefix,
                                                bool include_median) const
{
    if (!trajectory_data.contains("timesteps") || !trajectory_data["timesteps"].is_array()) {
        return false;
    }

    const auto& timesteps = trajectory_data["timesteps"];
    if (timesteps.empty()) {
        return false;
    }

    // Step 1: Extract common q-grid from first frame with scattering data
    std::vector<double> common_q_values;
    int n_qpoints = 0;

    for (const auto& timestep : timesteps) {
        if (timestep.contains("scattering") && timestep["scattering"].contains("q_values")) {
            const auto& q_arr = timestep["scattering"]["q_values"];
            if (q_arr.is_array() && !q_arr.empty()) {
                n_qpoints = q_arr.size();
                for (const auto& q : q_arr) {
                    if (q.is_number()) {
                        common_q_values.push_back(q.get<double>());
                    }
                }
                break;  // Use first frame's q-grid
            }
        }
    }

    if (common_q_values.empty()) {
        return false;  // No scattering data found
    }

    // Step 2: Initialize statistics storage for each q-point
    std::vector<QPointStatistics> q_stats(n_qpoints);

    // Step 3: Accumulate data across all frames (Welford's algorithm for variance)
    int frame_count = 0;
    for (const auto& timestep : timesteps) {
        if (!timestep.contains("scattering")) continue;

        const auto& scattering = timestep["scattering"];
        if (!scattering.contains("P_q") || !scattering.contains("S_q")) continue;
        if (!scattering["P_q"].is_array() || !scattering["S_q"].is_array()) continue;

        const auto& P_q = scattering["P_q"];
        const auto& S_q = scattering["S_q"];

        // Ensure array sizes match
        if (P_q.size() != n_qpoints || S_q.size() != n_qpoints) continue;

        frame_count++;

        for (int i = 0; i < n_qpoints; ++i) {
            if (!P_q[i].is_number() || !S_q[i].is_number()) continue;

            double p_val = P_q[i].get<double>();
            double s_val = S_q[i].get<double>();

            auto& stats = q_stats[i];

            // Store for median calculation (only if memory-efficient mode disabled)
            if (include_median) {
                stats.P_q_values.push_back(p_val);
                stats.S_q_values.push_back(s_val);
            }

            // Welford's online algorithm for mean and variance
            stats.count++;

            // P(q) statistics
            double delta_P = p_val - stats.P_q_sum / std::max(1, stats.count - 1);
            stats.P_q_sum += p_val;
            double delta2_P = p_val - stats.P_q_sum / stats.count;
            stats.P_q_M2 += delta_P * delta2_P;

            // S(q) statistics
            double delta_S = s_val - stats.S_q_sum / std::max(1, stats.count - 1);
            stats.S_q_sum += s_val;
            double delta2_S = s_val - stats.S_q_sum / stats.count;
            stats.S_q_M2 += delta_S * delta2_S;
        }
    }

    if (frame_count == 0) {
        return false;  // No valid frames
    }

    // Step 4: Calculate median by sorting (requires full series) - only if enabled
    if (include_median) {
        for (auto& stats : q_stats) {
            if (!stats.P_q_values.empty()) {
                std::sort(stats.P_q_values.begin(), stats.P_q_values.end());
                std::sort(stats.S_q_values.begin(), stats.S_q_values.end());
            }
        }
    }

    // Step 5: Write statistics file
    std::string stats_filename = output_directory + "/" + file_prefix + "_statistics.csv";
    std::ofstream stats_file(stats_filename);

    if (!stats_file.is_open()) {
        return false;
    }

    // Header
    std::string mode_note = include_median ? "" : " (memory-efficient, no median)";
    stats_file << "# Cross-frame scattering statistics (N=" << frame_count << " frames)" << mode_note << "\n";
    if (include_median) {
        stats_file << "# q (Å⁻¹),P_avg,P_std,P_median,S_avg,S_std,S_median\n";
    } else {
        stats_file << "# q (Å⁻¹),P_avg,P_std,S_avg,S_std\n";
    }

    // Data rows
    for (int i = 0; i < n_qpoints; ++i) {
        const auto& stats = q_stats[i];

        if (stats.count == 0) continue;

        // Mean
        double P_mean = stats.P_q_sum / stats.count;
        double S_mean = stats.S_q_sum / stats.count;

        // Standard deviation (sample std: n-1 denominator)
        double P_std = (stats.count > 1) ? std::sqrt(stats.P_q_M2 / (stats.count - 1)) : 0.0;
        double S_std = (stats.count > 1) ? std::sqrt(stats.S_q_M2 / (stats.count - 1)) : 0.0;

        // Write CSV row
        stats_file << std::fixed << std::setprecision(6) << common_q_values[i] << "," << P_mean << "," << P_std;

        if (include_median) {
            // Median
            double P_median = 0.0;
            double S_median = 0.0;
            if (!stats.P_q_values.empty()) {
                int mid = stats.P_q_values.size() / 2;
                if (stats.P_q_values.size() % 2 == 0) {
                    P_median = (stats.P_q_values[mid - 1] + stats.P_q_values[mid]) / 2.0;
                    S_median = (stats.S_q_values[mid - 1] + stats.S_q_values[mid]) / 2.0;
                } else {
                    P_median = stats.P_q_values[mid];
                    S_median = stats.S_q_values[mid];
                }
            }
            stats_file << "," << P_median;
        }

        stats_file << "," << S_mean << "," << S_std;

        if (include_median) {
            double S_median = 0.0;
            if (!stats.S_q_values.empty()) {
                int mid = stats.S_q_values.size() / 2;
                if (stats.S_q_values.size() % 2 == 0) {
                    S_median = (stats.S_q_values[mid - 1] + stats.S_q_values[mid]) / 2.0;
                } else {
                    S_median = stats.S_q_values[mid];
                }
            }
            stats_file << "," << S_median;
        }

        stats_file << "\n";
    }

    stats_file.close();
    return true;
}
