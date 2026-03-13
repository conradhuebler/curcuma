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
#include "handlers/general_properties_handler.h"
#include "handlers/scattering_handler.h"
#include "handlers/rdf_handler.h"
#include "handlers/topology_handler.h"

AnalysisOutputDispatcher::AnalysisOutputDispatcher(
    const UnifiedAnalysis::AnalysisConfig& config,
    const TrajectoryWriter& writer)
    : m_config(config)
    , m_writer(writer)
{
    // Register all available handlers - Claude Generated 2026 (Phase 3+4)
    // IMPORTANT: GeneralPropertiesHandler FIRST - baseline analysis for all runs
    // Adding new analysis types: just add registerHandler() call here
    registerHandler(std::make_unique<GeneralPropertiesHandler>());
    registerHandler(std::make_unique<ScatteringOutputHandler>());
    registerHandler(std::make_unique<RDFOutputHandler>());
    registerHandler(std::make_unique<TopologyOutputHandler>());
}

void AnalysisOutputDispatcher::registerHandler(std::unique_ptr<IAnalysisOutputHandler> handler)
{
    m_handlers.push_back(std::move(handler));
}

void AnalysisOutputDispatcher::dispatch(const json& results, std::ostream& out)
{
    // NOTE: Main output is already written by outputResults() before calling this method
    // This method only handles analysis-specific auxiliary files (scattering, RDF, etc.)
    // DO NOT write main output here to avoid duplication!

    // Registry-based handler dispatch - Claude Generated 2026 (Phase 3)
    // NO HARDCODED CALLS - handlers checked and invoked via registry loop
    for (const auto& handler : m_handlers) {
        if (handler->isEnabled(m_config) && handler->hasData(results)) {
            handler->handleOutput(results, m_config, m_writer);
        }
    }
}

void AnalysisOutputDispatcher::dispatchToFile(const json& results, const std::string& filename)
{
    // NOTE: Main output file is already written by outputToFile() before calling this method
    // This method only handles analysis-specific auxiliary files (scattering, RDF, etc.)
    // DO NOT write main output file here to avoid duplication!

    // Registry-based handler dispatch - Claude Generated 2026 (Phase 3)
    // NO HARDCODED CALLS - handlers checked and invoked via registry loop
    // Pass filename to handlers so they can derive auxiliary file names from it
    for (const auto& handler : m_handlers) {
        if (handler->isEnabled(m_config) && handler->hasData(results)) {
            handler->handleOutput(results, m_config, m_writer, filename);
        }
    }
}
