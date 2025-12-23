/* trajectory_writer.cpp - Claude Generated 2025
 *
 * Unified trajectory output system implementation for Curcuma
 * Consolidates output code from analysis.cpp, trajectoryanalysis.cpp, rmsdtraj.cpp
 */

#include "trajectory_writer.h"
#include "../capabilities/trajectory_statistics.h"
#include <iomanip>
#include <algorithm>
#include <chrono>
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
                out << formatValue(gyr["unweighted"]);

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "gyration_unweighted", "mean"), 11)
                        << formatValue(getStatisticsValue(timestep, "gyration_unweighted", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "gyration_unweighted", "moving_avg"), 10);
                }

                // Mass-weighted gyration
                out << formatValue(gyr["mass_weighted"]);

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "gyration_mass", "mean"), 11)
                        << formatValue(getStatisticsValue(timestep, "gyration_mass", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "gyration_mass", "moving_avg"), 10);
                }
            } else {
                // Skip columns with "---"
                int skip_cols = 2 + (m_enable_cumulative ? 4 : 0) + (m_enable_moving ? 2 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "rout") {
            if (metricAvailable(timestep, "rout")) {
                out << formatValue(getMetricValue(timestep, "rout"));

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "rout", "mean"), 11)
                        << formatValue(getStatisticsValue(timestep, "rout", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "rout", "moving_avg"), 10);
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "end2end") {
            if (metricAvailable(timestep, "end2end")) {
                out << formatValue(getMetricValue(timestep, "end_to_end_distance"));

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "end_to_end_distance", "mean"), 11)
                        << formatValue(getStatisticsValue(timestep, "end_to_end_distance", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "end_to_end_distance", "moving_avg"), 10);
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "mass") {
            if (metricAvailable(timestep, "mass")) {
                out << formatValue(getMetricValue(timestep, "mass"));

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "mass", "mean"), 11)
                        << formatValue(getStatisticsValue(timestep, "mass", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "mass", "moving_avg"), 10);
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << std::setw(9) << "---";
            }
        } else if (metric == "com") {
            if (metricAvailable(timestep, "com")) {
                const auto& com = getMetricValue(timestep, "com");
                out << std::setw(7) << std::fixed << std::setprecision(3)
                    << "(" << com[0].get<double>() << ","
                    << com[1].get<double>() << ","
                    << com[2].get<double>() << ")";

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "com", "mean"), 11)
                        << formatValue(getStatisticsValue(timestep, "com", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << formatValue(getStatisticsValue(timestep, "com", "moving_avg"), 10);
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
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
                out << formatValueCSV(gyr["unweighted"]);

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "gyration_unweighted", "mean"))
                        << "," << formatValueCSV(getStatisticsValue(timestep, "gyration_unweighted", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "gyration_unweighted", "moving_avg"));
                }

                // Mass-weighted gyration
                out <<  "," << formatValueCSV(gyr["mass_weighted"]);

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "gyration_mass", "mean"))
                        << "," << formatValueCSV(getStatisticsValue(timestep, "gyration_mass", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "gyration_mass", "moving_avg"));
                }
            } else {
                // Skip columns with empty commas
                int skip_cols = 2 + (m_enable_cumulative ? 4 : 0) + (m_enable_moving ? 2 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "rout") {
            if (metricAvailable(timestep, "rout")) {
                out << "," << formatValueCSV(getMetricValue(timestep, "rout"));

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "rout", "mean"))
                        << "," << formatValueCSV(getStatisticsValue(timestep, "rout", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "rout", "moving_avg"));
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "end2end") {
            if (metricAvailable(timestep, "end2end")) {
                out << "," << formatValueCSV(getMetricValue(timestep, "end_to_end_distance"));

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "end_to_end_distance", "mean"))
                        << "," << formatValueCSV(getStatisticsValue(timestep, "end_to_end_distance", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "end_to_end_distance", "moving_avg"));
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "mass") {
            if (metricAvailable(timestep, "mass")) {
                out << "," << formatValueCSV(getMetricValue(timestep, "mass"));

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "mass", "mean"))
                        << "," << formatValueCSV(getStatisticsValue(timestep, "mass", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "mass", "moving_avg"));
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
                    out << ",";
            }
        } else if (metric == "com") {
            if (metricAvailable(timestep, "com")) {
                const auto& com = getMetricValue(timestep, "com");
                out << std::fixed << std::setprecision(3)
                    << ",(" << com[0].get<double>() << ","
                    << com[1].get<double>() << ","
                    << com[2].get<double>() << ")";

                if (m_enable_cumulative && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "com", "mean"))
                        << "," << formatValueCSV(getStatisticsValue(timestep, "com", "std"));
                }
                if (m_enable_moving && timestep.contains("statistics")) {
                    out << "," << formatValueCSV(getStatisticsValue(timestep, "com", "moving_avg"));
                }
            } else {
                int skip_cols = 1 + (m_enable_cumulative ? 2 : 0) + (m_enable_moving ? 1 : 0);
                for (int j = 0; j < skip_cols; ++j)
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

std::string TrajectoryWriter::formatValueCSV(double value) const
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << value;
    return oss.str();
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

// ---------- ProgressTracker Implementation ----------

ProgressTracker::ProgressTracker(int report_interval)
    : m_report_interval(report_interval), m_auto_timing(false)
{
    m_start_time = std::chrono::steady_clock::now();
}

void ProgressTracker::reportProgress(int current, int total, const std::string& custom_message)
{
    if (current <= 0 || total <= 0) return;

    int current_percent = int((double(current) / total) * 100);

    // Initialize progress tracking vector if needed
    if (m_progress_reported.empty()) {
        initializeProgressVector(total);
    }

    int progress_step = current_percent / m_report_interval;

    if (!m_progress_reported[progress_step] && current_percent >= m_report_interval) {
        internalReport(current, total);
        m_progress_reported[progress_step] = true;
    }
}

void ProgressTracker::initializeProgressVector(int total)
{
    size_t vector_size = (100 / m_report_interval) + 1;
    m_progress_reported.assign(vector_size, false);
    m_start_time = std::chrono::steady_clock::now();
}

bool ProgressTracker::shouldReport(int current_percent)
{
    return current_percent >= m_report_interval && (current_percent % m_report_interval == 0);
}

void ProgressTracker::internalReport(int current, int total)
{
    std::string message = m_message;

    if (m_auto_timing) {
        double elapsed = getElapsedTime();
        double eta = getEstimatedTimeRemaining(current, total);
        message += fmt::format(" (Elapsed: {}, ETA: {})", formatTime(elapsed), formatTime(eta));
    }

    if (m_callback) {
        m_callback(current, total, message);
    } else {
        CurcumaLogger::progress(current, total, message);
    }
}

double ProgressTracker::getElapsedTime() const
{
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - m_start_time);
    return elapsed.count();
}

double ProgressTracker::getEstimatedTimeRemaining(int current, int total) const
{
    int processed = current;
    if (processed == 0) return 0.0;

    double elapsed = getElapsedTime();
    return (elapsed * (total - processed)) / processed;
}

std::string ProgressTracker::formatTime(double seconds) const
{
    if (seconds < 60) {
        return fmt::format("{:.0f}s", seconds);
    } else if (seconds < 3600) {
        int minutes = int(seconds / 60);
        int secs = int(seconds) % 60;
        return fmt::format("{}m{:02d}s", minutes, secs);
    } else {
        int hours = int(seconds / 3600);
        int minutes = int((int(seconds) % 3600) / 60);
        return fmt::format("{}h{:02d}m", hours, minutes);
    }
}

void ProgressTracker::reset()
{
    m_progress_reported.clear();
    m_start_time = std::chrono::steady_clock::now();
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