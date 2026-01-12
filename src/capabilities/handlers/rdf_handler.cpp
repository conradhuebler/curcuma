/*
 * <RDF analysis output handler>
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

#include "rdf_handler.h"

bool RDFOutputHandler::isEnabled(const UnifiedAnalysis::AnalysisConfig& config) const
{
    return config.hasRDFOutput();
}

bool RDFOutputHandler::hasData(const json& results) const
{
    // Placeholder: Check if results contain RDF data
    // Will be implemented when RDF analysis is added to results
    if (!results.contains("timesteps") || !results["timesteps"].is_array()) {
        return false;
    }

    for (const auto& timestep : results["timesteps"]) {
        if (timestep.contains("rdf")) {
            return true;
        }
    }

    return false;
}

void RDFOutputHandler::handleOutput(const json& results,
                                   const UnifiedAnalysis::AnalysisConfig& config,
                                   const TrajectoryWriter& writer,
                                   const std::string& output_filename)
{
    // Placeholder for future RDF output implementation
    // When RDF analysis is added to results, implement output logic here
    if (!isEnabled(config)) return;
    if (!hasData(results)) return;

    // TODO: Implement RDF-specific output generation
    // - Generate g(r) CSV files (use output_filename for naming)
    // - Generate coordination number n(r) if enabled
    // - Aggregate statistics across frames
}
