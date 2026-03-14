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

    // Claude Generated 2026: Generate gnuplot visualization script
    generateGnuplotScript(output_dir, file_prefix);
}

// Claude Generated 2026: Generate gnuplot script for scattering visualization
void ScatteringOutputHandler::generateGnuplotScript(
    const std::string& output_dir,
    const std::string& file_prefix) const
{
    std::string gnu_file = output_dir + "/" + file_prefix + ".scattering.gnu";
    std::ofstream out(gnu_file);

    if (!out.is_open()) {
        CurcumaLogger::warn_fmt("Failed to create gnuplot script: {}", gnu_file);
        return;
    }

    // Data files (generated by handleOutput)
    // Note: writeScatteringStatistics generates filename with "_statistics.csv" suffix
    std::string stats_file = file_prefix + ".scattering_statistics.csv";
    std::string output_png = file_prefix + ".scattering_plot.png";

    // Gnuplot template (extracted from gnuplot_workflow.sh lines 133-175)
    out << "#!/usr/bin/gnuplot\n";
    out << "# Claude Generated 2026: Scattering analysis visualization\n";
    out << "# Usage: gnuplot " << file_prefix << ".scattering.gnu\n\n";

    out << "set terminal pngcairo size 1200,800 font \"Arial,12\"\n";
    out << "set output '" << output_png << "'\n";
    out << "set datafile separator ','\n";  // CSV format
    out << "set style data lines\n";
    out << "set grid back lt 1 lc rgb \"#e0e0e0\"\n\n";

    out << "set multiplot layout 2,2\n\n";

    // Plot 1: P(q) Log-Log
    out << "# Plot 1: Form Factor P(q) [Log-Log]\n";
    out << "set logscale xy\n";
    out << "set title \"Form Factor P(q) [Log-Log]\"\n";
    out << "set xlabel \"q (Å⁻¹)\"\n";
    out << "set ylabel \"P(q)\"\n";
    out << "plot '" << stats_file << "' using 1:2 with lines lw 2 lc rgb \"#0072B2\" title \"P(q) (mean)\"\n\n";

    // Plot 2: S(q) with ideal gas reference
    out << "# Plot 2: Structure Factor S(q)\n";
    out << "unset logscale y\n";
    out << "set logscale x\n";
    out << "set title \"Structure Factor S(q)\"\n";
    out << "set xlabel \"q (Å⁻¹)\"\n";
    out << "set ylabel \"S(q)\"\n";
    out << "plot '" << stats_file << "' using 1:5 with lines lw 2 lc rgb \"#E69F00\" title \"S(q) (mean)\", \\\n";
    out << "     '" << stats_file << "' using 1:(1.0) with lines lt 2 lc rgb \"#999999\" title \"Ideal Gas\"\n\n";

    // Plot 3: P(q) + S(q) Combined
    out << "# Plot 3: Form Factor vs Structure Factor\n";
    out << "set logscale y\n";
    out << "set title \"Form Factor vs Structure Factor\"\n";
    out << "set xlabel \"q (Å⁻¹)\"\n";
    out << "set ylabel \"P(q), S(q)\"\n";
    out << "plot '" << stats_file << "' using 1:2 with lines lw 2 lc rgb \"#0072B2\" title \"P(q)\", \\\n";
    out << "     '" << stats_file << "' using 1:5 with lines lw 2 lc rgb \"#E69F00\" title \"S(q)\"\n\n";

    // Plot 4: Guinier Analysis
    out << "# Plot 4: Guinier Region: ln(P(q)) vs q²\n";
    out << "unset logscale xy\n";
    out << "set title \"Guinier Region: ln(P(q)) vs q²\"\n";
    out << "set xlabel \"q² (Å⁻²)\"\n";
    out << "set ylabel \"ln(P(q))\"\n";
    out << "set xrange [0:0.01]\n";
    out << "plot '" << stats_file << "' using ($1*$1):(log($2)) with lines lw 2 lc rgb \"#0072B2\" title \"ln(P(q))\"\n\n";

    out << "unset multiplot\n";
    out << "# End of gnuplot script\n";

    out.close();

    CurcumaLogger::success_fmt("Gnuplot script saved: {}", gnu_file);
    CurcumaLogger::info_fmt("To generate plot: gnuplot {}", gnu_file);
}
