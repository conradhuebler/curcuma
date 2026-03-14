/*
 * <Parallel GFN-FF Parameter Generation Thread>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated (February 2026): Inter-phase parallelisation for GFN-FF
 * parameter generation. Each thread handles one independent generation phase
 * (angles, torsions, inversions, coulomb, repulsion, dispersion) after bonds
 * have been computed sequentially.
 *
 * Thread safety: All generate*() methods only read shared topology data.
 * Each thread writes to its own local JSON result. No locking required.
 */

#pragma once

#include "external/CxxThreadPool/include/CxxThreadPool.hpp"
#include "json.hpp"
#include <functional>
#include <string>

using json = nlohmann::json;

/**
 * @brief CxxThread subclass for parallel GFN-FF parameter generation phases
 *
 * Wraps a parameter generation function (e.g. generateGFNFFTorsions) as a
 * CxxThread task for execution in CxxThreadPool. Each thread stores its
 * result in a thread-local JSON object.
 *
 * Usage:
 *   auto* t = new ParameterGeneratorThread("torsions", [&]() { return gfnff.generateGFNFFTorsions(); });
 *   pool.addThread(t);
 *   pool.StartAndWait();
 *   json result = t->getResult();
 */
class ParameterGeneratorThread : public CxxThread {
public:
    /**
     * @brief Construct with phase name and generation function
     * @param phase_name Human-readable name for timing output (e.g. "torsions")
     * @param generator Lambda/function that returns JSON result
     */
    ParameterGeneratorThread(const std::string& phase_name, std::function<json()> generator)
        : m_phase_name(phase_name)
        , m_generator(std::move(generator))
    {
        // AutoDelete defaults to true - CxxThreadPool destructor handles cleanup
        // Results are read between StartAndWait() and pool destruction
    }

    int execute() override
    {
        m_result = m_generator();
        return 0;
    }

    /// Get the JSON result after execution completes
    const json& getResult() const { return m_result; }

    /// Get phase name for timing output
    const std::string& getPhaseName() const { return m_phase_name; }

private:
    std::string m_phase_name;
    std::function<json()> m_generator;
    json m_result;
};
