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

#pragma once

#include "../analysis_output_handler.h"

/*! \brief Handler for scattering analysis output - Claude Generated 2026
 *
 * Implements the IAnalysisOutputHandler interface for scattering-specific output:
 * - Per-frame CSV files (if enabled)
 * - Cross-frame statistics (P(q) and S(q) aggregation)
 * - Configuration-based conditional execution
 *
 * Extracted from AnalysisOutputDispatcher::handleScatteringOutput() to support
 * registry-based handler dispatch and independent testing.
 */
class ScatteringOutputHandler : public IAnalysisOutputHandler {
public:
    /*! \brief Check if scattering output is enabled in configuration */
    bool isEnabled(const UnifiedAnalysis::AnalysisConfig& config) const override;

    /*! \brief Check if results contain scattering data */
    bool hasData(const json& results) const override;

    /*! \brief Generate scattering output files and statistics */
    void handleOutput(const json& results,
                     const UnifiedAnalysis::AnalysisConfig& config,
                     const TrajectoryWriter& writer,
                     const std::string& output_filename = "") override;

    /*! \brief Identifier for this handler */
    std::string name() const override { return "Scattering"; }

private:
    /*! \brief Generate gnuplot script for scattering visualization - Claude Generated 2026 */
    void generateGnuplotScript(const std::string& output_dir,
                              const std::string& file_prefix) const;
};
