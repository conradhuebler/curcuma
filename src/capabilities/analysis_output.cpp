/*
 * <Analysis output dispatcher - Eliminates code duplication>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "analysis_output.h"
#include "src/core/curcuma_logger.h"

AnalysisOutputDispatcher::AnalysisOutputDispatcher(
    const UnifiedAnalysis::AnalysisConfig& config,
    const TrajectoryWriter& writer)
    : m_config(config)
    , m_writer(writer)
{
}

void AnalysisOutputDispatcher::dispatch(const json& results, std::ostream& out)
{
    // Main output routing based on format
    if (m_config.output_format == "json") {
        m_writer.writeJSON(out, results);
    } else if (m_config.output_format == "csv") {
        m_writer.writeCSV(out, results);
    } else {
        // Default: human-readable table
        m_writer.writeHumanTable(out, results);
    }

    // Analysis-specific outputs (replaces duplicate code)
    handleScatteringOutput(results);
    handleRDFOutput(results);
    handleTopologyOutput(results);
}

void AnalysisOutputDispatcher::dispatchToFile(const json& results, const std::string& filename)
{
    // Determine format from filename or config
    TrajectoryWriter::Format format = TrajectoryWriter::Format::HumanTable;

    bool has_json_ext = (filename.size() >= 5 && filename.substr(filename.size() - 5) == ".json");
    bool has_csv_ext = (filename.size() >= 4 && filename.substr(filename.size() - 4) == ".csv");

    if (m_config.output_format == "json" || has_json_ext) {
        format = TrajectoryWriter::Format::JSON;
    } else if (m_config.output_format == "csv" || has_csv_ext) {
        format = TrajectoryWriter::Format::CSV;
    }

    // Write main output file
    m_writer.writeToFile(filename, format, results);

    // Analysis-specific outputs (replaces duplicate code)
    handleScatteringOutput(results);
    handleRDFOutput(results);
    handleTopologyOutput(results);
}

void AnalysisOutputDispatcher::handleScatteringOutput(const json& results)
{
    // Early exit if scattering not enabled or no data
    if (!m_config.hasScatteringOutput()) return;
    if (!hasScatteringData(results)) return;

    std::string output_dir = m_config.scattering_output_directory;
    std::string file_prefix = m_config.scattering_file_prefix;

    // Per-frame files (replaces lines 1117-1145 and 1221-1246)
    if (m_config.hasPerFrameFiles()) {
        bool success = m_writer.writeScatteringPerFrameFiles(results, output_dir, file_prefix);
        if (success) {
            CurcumaLogger::success_fmt("Per-frame scattering files generated in: {}", output_dir);
        } else {
            CurcumaLogger::error("Failed to generate per-frame scattering files");
        }
    }

    // Cross-frame statistics (replaces lines 1147-1187 and 1248-1272)
    bool success = m_writer.writeScatteringStatistics(
        results, output_dir, file_prefix,
        m_config.scattering_stats_include_median);

    if (success) {
        CurcumaLogger::success_fmt("Cross-frame scattering statistics saved to: {}/{}_statistics.csv",
                                   output_dir, file_prefix);
    }
}

void AnalysisOutputDispatcher::handleRDFOutput(const json& results)
{
    // Placeholder for RDF output
    // TODO: Implement RDF-specific output logic here
    if (!m_config.hasRDFOutput()) return;

    // Future implementation will go here
}

void AnalysisOutputDispatcher::handleTopologyOutput(const json& results)
{
    // Placeholder for topology output
    // TODO: Implement topology-specific output logic here

    // Future implementation will go here
}

bool AnalysisOutputDispatcher::hasScatteringData(const json& results) const
{
    if (!results.contains("timesteps") || !results["timesteps"].is_array()) {
        return false;
    }

    // Check if at least one timestep has scattering data
    for (const auto& timestep : results["timesteps"]) {
        if (timestep.contains("scattering")) {
            return true;
        }
    }

    return false;
}
