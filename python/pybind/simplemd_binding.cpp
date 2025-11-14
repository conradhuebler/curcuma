/*
 * <Python Bindings for Molecular Dynamics - SIMPLIFIED>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated: Simplified Python bindings - SimpleMD is CLI-focused
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "src/capabilities/simplemd.h"
#include "src/core/molecule.h"

namespace py = pybind11;
using namespace curcuma;

/**
 * @brief Simplified Python bindings for molecular dynamics
 *
 * NOTE: SimpleMD is designed for CLI usage and writes results to files.
 * Full programmatic API not yet available.
 *
 * Claude Generated: Minimal bindings matching actual API
 */
void bind_simplemd(py::module& m) {
    // Convenience function placeholder
    m.def("run_molecular_dynamics",
        [](const Molecule& mol, const std::string& method, double temperature,
           int steps, double timestep, int verbosity) {

            py::print("Warning: run_molecular_dynamics() not yet fully implemented");
            py::print("SimpleMD class is CLI-focused and not suitable for Python bindings yet");

            // Return copy of input molecule for now
            return Molecule(mol);
        },
        py::arg("molecule"),
        py::arg("method") = "uff",
        py::arg("temperature") = 300.0,
        py::arg("steps") = 1000,
        py::arg("timestep") = 0.5,
        py::arg("verbosity") = 1,
        R"pbdoc(
            Run molecular dynamics simulation (PLACEHOLDER - not yet implemented)

            WARNING: This function is not yet fully implemented.
            The SimpleMD class is designed for CLI usage and does not
            provide a programmatic API suitable for Python bindings.

            Args:
                molecule: Initial molecular structure
                method: Calculation method (default: "uff")
                temperature: Temperature in Kelvin (default: 300)
                steps: Number of MD steps (default: 1000)
                timestep: Integration timestep in fs (default: 0.5)
                verbosity: Output level 0-3 (default: 1)

            Returns:
                Final molecular structure (currently returns input unchanged)

            Note:
                This is a placeholder. Use curcuma CLI for MD:
                $ curcuma -md input.xyz -method uff -temp 300 -steps 1000
        )pbdoc");
}
