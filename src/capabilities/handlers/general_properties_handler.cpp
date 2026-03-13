/*
 * <General properties analysis output handler>
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

#include "general_properties_handler.h"
#include "src/core/curcuma_logger.h"
#include <fstream>

bool GeneralPropertiesHandler::isEnabled(const UnifiedAnalysis::AnalysisConfig& config) const
{
    // General properties are the baseline analysis - always generate output
    return true;
}

bool GeneralPropertiesHandler::hasData(const json& results) const
{
    // Check if results contain timesteps with geometric/polymer data
    if (!results.contains("timesteps") || !results["timesteps"].is_array()) {
        return false;
    }

    // Check if at least one timestep has geometric or polymer data
    for (const auto& timestep : results["timesteps"]) {
        if (timestep.contains("geometric") || timestep.contains("polymer")) {
            return true;
        }
    }

    return false;
}

void GeneralPropertiesHandler::handleOutput(const json& results,
                                           const UnifiedAnalysis::AnalysisConfig& config,
                                           const TrajectoryWriter& writer,
                                           const std::string& output_filename)
{
    if (output_filename.empty()) return;
    if (!hasData(results)) return;

    // Extract basename and directory from output_filename
    std::string basename;
    std::string output_dir = ".";

    // Extract directory: "/path/to/file" → "/path/to"
    size_t last_slash = output_filename.find_last_of("/\\");
    if (last_slash != std::string::npos) {
        output_dir = output_filename.substr(0, last_slash);
        basename = output_filename.substr(last_slash + 1);
    } else {
        basename = output_filename;
    }

    // Remove extension from basename if present: "file.json" → "file"
    size_t last_dot = basename.find_last_of(".");
    if (last_dot != std::string::npos) {
        basename = basename.substr(0, last_dot);
    }

    // Create output filename: basename.general.csv
    std::string general_file = output_dir + "/" + basename + ".general.csv";

    // Generate CSV content manually (TrajectoryWriter uses different format)
    std::ofstream outfile(general_file);
    if (!outfile.is_open()) {
        CurcumaLogger::error_fmt("Failed to open file for writing: {}", general_file);
        return;
    }

    // Write header
    outfile << "Frame,Gyr_u,Gyr_m,Rout,End2End\n";

    // Write data for each timestep
    const auto& timesteps = results["timesteps"];
    for (size_t i = 0; i < timesteps.size(); ++i) {
        const auto& ts = timesteps[i];

        // Extract data (with defaults if missing)
        double gyr_u = 0.0, gyr_m = 0.0, rout = 0.0, end2end = 0.0;

        if (ts.contains("geometric")) {
            const auto& geom = ts["geometric"];
            if (geom.contains("gyration_radius")) {
                const auto& gyr = geom["gyration_radius"];
                if (gyr.contains("unweighted")) gyr_u = gyr["unweighted"].get<double>();
                if (gyr.contains("mass_weighted")) gyr_m = gyr["mass_weighted"].get<double>();
            }
        }

        if (ts.contains("polymer")) {
            const auto& poly = ts["polymer"];
            if (poly.contains("rout")) rout = poly["rout"].get<double>();
            if (poly.contains("end_to_end_distance")) end2end = poly["end_to_end_distance"].get<double>();
        }

        // Write row: Frame number (1-based), Gyr_u, Gyr_m, Rout, End2End
        outfile << (i + 1) << ","
                << gyr_u << ","
                << gyr_m << ","
                << rout << ","
                << end2end << "\n";
    }

    outfile.close();

    CurcumaLogger::success_fmt("General properties saved to: {}", general_file);
}
