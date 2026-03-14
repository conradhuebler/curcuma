/*
 * <Topology analysis output handler>
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

/*! \brief Handler for topological analysis output (dMatrix, persistence diagrams) - Claude Generated 2026
 *
 * Implements the IAnalysisOutputHandler interface for topology-specific output.
 * Currently a placeholder for future topological data analysis output.
 */
class TopologyOutputHandler : public IAnalysisOutputHandler {
public:
    /*! \brief Check if topology output is enabled in configuration */
    bool isEnabled(const UnifiedAnalysis::AnalysisConfig& config) const override;

    /*! \brief Check if results contain topology data */
    bool hasData(const json& results) const override;

    /*! \brief Generate topology output (placeholder - future implementation) */
    void handleOutput(const json& results,
                     const UnifiedAnalysis::AnalysisConfig& config,
                     const TrajectoryWriter& writer,
                     const std::string& output_filename = "") override;

    /*! \brief Identifier for this handler */
    std::string name() const override { return "Topology"; }
};
