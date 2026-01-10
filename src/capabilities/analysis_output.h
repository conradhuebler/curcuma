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

#pragma once

#include "json.hpp"
#include "src/tools/trajectory_writer.h"
#include "analysis.h"

using json = nlohmann::json;

/*! \brief Unified output dispatcher for analysis results - Claude Generated 2026
 *
 * Eliminates ~160 lines of duplicate code between outputResults() and outputToFile().
 * Centralizes all analysis-specific output handlers (Scattering, RDF, Topology).
 */
class AnalysisOutputDispatcher {
public:
    /*! \brief Constructor takes configuration and writer
     *  @param config Analysis configuration state
     *  @param writer TrajectoryWriter for main output
     */
    AnalysisOutputDispatcher(const UnifiedAnalysis::AnalysisConfig& config,
                            const TrajectoryWriter& writer);

    /*! \brief Dispatch output to stream (for outputResults)
     *  @param results Analysis results JSON
     *  @param out Output stream (typically std::cout)
     */
    void dispatch(const json& results, std::ostream& out);

    /*! \brief Dispatch output to file (for outputToFile)
     *  @param results Analysis results JSON
     *  @param filename Output file path
     */
    void dispatchToFile(const json& results, const std::string& filename);

private:
    const UnifiedAnalysis::AnalysisConfig& m_config;
    TrajectoryWriter m_writer;

    /*! \brief Handle scattering-specific output (per-frame files + statistics)
     *  Replaces duplicate code from lines 1117-1187 and 1221-1272
     */
    void handleScatteringOutput(const json& results);

    /*! \brief Handle RDF-specific output (placeholder for future)
     */
    void handleRDFOutput(const json& results);

    /*! \brief Handle topology-specific output (placeholder for future)
     */
    void handleTopologyOutput(const json& results);

    /*! \brief Check if results contain scattering data
     */
    bool hasScatteringData(const json& results) const;
};
