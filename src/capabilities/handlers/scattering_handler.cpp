/*
 * <Scattering analysis output handler>
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

#include "scattering_handler.h"
#include "src/core/curcuma_logger.h"
#include <fstream>
#include <cstdio>  // for snprintf

bool ScatteringOutputHandler::isEnabled(const UnifiedAnalysis::AnalysisConfig& config) const
{
    return config.hasScatteringOutput();
}

bool ScatteringOutputHandler::hasData(const json& results) const
{
    // Check if results contain scattering data in at least one timestep
    if (!results.contains("timesteps") || !results["timesteps"].is_array()) {
        return false;
    }

    for (const auto& timestep : results["timesteps"]) {
        if (timestep.contains("scattering")) {
            return true;
        }
    }

    return false;
}

void ScatteringOutputHandler::handleOutput(const json& results,
                                          const UnifiedAnalysis::AnalysisConfig& config,
                                          const TrajectoryWriter& writer,
                                          const std::string& output_filename)
{
    // Early exit if scattering not enabled or no data
    if (!isEnabled(config)) return;
    if (!hasData(results)) return;

    std::string output_dir = config.scattering_output_directory;
    std::string file_prefix = config.scattering_file_prefix;

    // If output_filename is provided, derive file_prefix from it
    if (!output_filename.empty()) {
        // Extract basename: "/path/to/file.json" -> "file"
        size_t last_slash = output_filename.find_last_of("/\\");
        std::string basename = (last_slash != std::string::npos)
            ? output_filename.substr(last_slash + 1)
            : output_filename;

        // Remove extension: "file.json" -> "file"
        size_t last_dot = basename.find_last_of(".");
        if (last_dot != std::string::npos) {
            file_prefix = basename.substr(0, last_dot);
        } else {
            file_prefix = basename;
        }

        // Extract directory from output_filename if present
        if (last_slash != std::string::npos) {
            output_dir = output_filename.substr(0, last_slash);
        }
    }

    // Per-frame files (replaces lines 1117-1145 and 1221-1246)
    // Phase 4: New naming - basename.001.scattering.csv (combined P(q) and S(q))
    if (config.hasPerFrameFiles()) {
        if (!results.contains("timesteps") || !results["timesteps"].is_array()) {
            CurcumaLogger::error("No timesteps data for per-frame scattering output");
            return;
        }

        const auto& timesteps = results["timesteps"];
        int files_written = 0;

        for (size_t i = 0; i < timesteps.size(); ++i) {
            if (!timesteps[i].contains("scattering")) continue;

            const auto& scattering = timesteps[i]["scattering"];
            // Keys are "P_q" and "S_q" (with underscore)
            if (!scattering.contains("P_q") || !scattering.contains("S_q")) continue;

            // New naming: basename.001.scattering.csv
            char frame_file[256];
            snprintf(frame_file, sizeof(frame_file), "%s.%03zu.scattering.csv",
                    file_prefix.c_str(), i + 1);
            std::string full_path = output_dir + "/" + frame_file;

            // Write combined P(q) and S(q) to single CSV
            std::ofstream out(full_path);
            if (!out.is_open()) {
                CurcumaLogger::error_fmt("Failed to open file: {}", full_path);
                continue;
            }

            // Write header
            out << "q,P_q,S_q\n";

            // Write data (P(q) and S(q) should have same q-values)
            const auto& q_values = scattering["q_values"];
            const auto& pq_values = scattering["P_q"];
            const auto& sq_values = scattering["S_q"];

            for (size_t j = 0; j < q_values.size(); ++j) {
                out << q_values[j].get<double>() << ","
                    << pq_values[j].get<double>() << ","
                    << sq_values[j].get<double>() << "\n";
            }

            out.close();
            files_written++;
        }

        if (files_written > 0) {
            CurcumaLogger::success_fmt("Per-frame scattering files generated in: {} ({} files)",
                                       output_dir, files_written);
        } else {
            CurcumaLogger::error("No scattering data found for per-frame output");
        }
    }

    // Cross-frame statistics (replaces lines 1147-1187 and 1248-1272)
    // Phase 4: New naming - basename.scattering instead of basename
    // writeScatteringStatistics will append ".csv" itself
    std::string stats_prefix = file_prefix + ".scattering";
    bool success = writer.writeScatteringStatistics(
        results, output_dir, stats_prefix,
        config.scattering_stats_include_median);

    if (success) {
        CurcumaLogger::success_fmt("Cross-frame scattering statistics saved to: {}/{}.csv",
                                   output_dir, stats_prefix);
    }
}
