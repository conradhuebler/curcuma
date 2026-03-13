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

#include <memory>
#include <vector>
#include "json.hpp"
#include "src/tools/trajectory_writer.h"
#include "analysis.h"
#include "analysis_output_handler.h"

using json = nlohmann::json;

/*! \brief Unified output dispatcher for analysis results - Claude Generated 2026
 *
 * **Phase 2 (DEPRECATED):** Eliminates ~160 lines of duplicate code between outputResults() and outputToFile().
 *
 * **Phase 3 (CURRENT):** Registry-based handler dispatch enables pluggable analysis outputs.
 * - Maintains dispatcher role (main output routing)
 * - Delegates analysis-specific output to registered handlers
 * - Allows new analysis types without modifying dispatcher code
 *
 * **Architecture:**
 * 1. dispatch() routes main output (JSON/CSV/Human)
 * 2. Handler loop iterates registered IAnalysisOutputHandler instances
 * 3. Each handler checks: isEnabled() && hasData() → handleOutput()
 * 4. Adding new analysis: Create handler + register in constructor (zero dispatcher changes)
 */
class AnalysisOutputDispatcher {
public:
    /*! \brief Constructor takes configuration and writer
     *  @param config Analysis configuration state
     *  @param writer TrajectoryWriter for main output
     *
     * During construction, all available handlers are registered
     * (Scattering, RDF, Topology, etc.)
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
    std::vector<std::unique_ptr<IAnalysisOutputHandler>> m_handlers;

    /*! \brief Register a handler for analysis-specific output
     *  @param handler Unique pointer to handler instance
     *
     * Called during construction to register all available handlers.
     * Each handler will be invoked if isEnabled() && hasData() are true.
     */
    void registerHandler(std::unique_ptr<IAnalysisOutputHandler> handler);
};
