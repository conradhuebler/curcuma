/*
 * <Analysis output handler interface - Strategy pattern for pluggable handlers>
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

#include <string>
#include <memory>
#include "json.hpp"
#include "analysis.h"
#include "src/tools/trajectory_writer.h"

using json = nlohmann::json;

/*! \brief Abstract interface for analysis-specific output handlers - Claude Generated 2026
 *
 * Implements Strategy pattern to allow pluggable handlers without modifying dispatcher.
 * Each handler is responsible for:
 * - Determining if it should run (based on configuration)
 * - Checking if it has data to output (in results)
 * - Performing analysis-specific output generation
 *
 * This design enables adding new analysis types without touching AnalysisOutputDispatcher:
 * 1. Create new class inheriting IAnalysisOutputHandler
 * 2. Implement 4 virtual methods
 * 3. Register in AnalysisOutputDispatcher constructor
 * Done - zero changes to existing code!
 */
class IAnalysisOutputHandler {
public:
    virtual ~IAnalysisOutputHandler() = default;

    /*! \brief Check if this handler should run based on configuration
     *  @param config Analysis configuration state
     *  @return true if this handler is enabled and should run
     */
    virtual bool isEnabled(const UnifiedAnalysis::AnalysisConfig& config) const = 0;

    /*! \brief Check if results contain data relevant to this handler
     *  @param results Analysis results JSON from all timesteps
     *  @return true if results have data this handler can process
     */
    virtual bool hasData(const json& results) const = 0;

    /*! \brief Execute handler-specific output generation
     *  @param results Analysis results JSON
     *  @param config Analysis configuration state
     *  @param writer TrajectoryWriter for file output operations
     *  @param output_filename Optional main output filename for deriving auxiliary file names
     */
    virtual void handleOutput(const json& results,
                             const UnifiedAnalysis::AnalysisConfig& config,
                             const TrajectoryWriter& writer,
                             const std::string& output_filename = "") = 0;

    /*! \brief Handler identification for debugging and logging
     *  @return Human-readable handler name (e.g., "Scattering", "RDF", "Topology")
     */
    virtual std::string name() const = 0;
};
