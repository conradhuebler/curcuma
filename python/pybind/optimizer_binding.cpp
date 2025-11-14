/*
 * <Python Bindings for Geometry Optimization - SIMPLIFIED>
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
 * Claude Generated: Simplified Python bindings - CurcumaOpt is CLI-focused
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "src/capabilities/curcumaopt.h"
#include "src/core/molecule.h"

namespace py = pybind11;
using namespace curcuma;

/**
 * @brief Simplified Python bindings for geometry optimization
 *
 * NOTE: CurcumaOpt is designed for CLI usage and writes results to files.
 * Full programmatic API not yet available.
 * Use the convenience function optimize_geometry() for simple use cases.
 *
 * Claude Generated: Minimal bindings matching actual API
 */
void bind_optimizer(py::module& m) {
    // For now, don't bind CurcumaOpt class directly - it's too CLI-focused
    // Only provide the convenience function

    // Convenience function for simple optimization
    m.def("optimize_geometry",
        [](const Molecule& mol, const std::string& method, int max_iter, int verbosity) {
            // CurcumaOpt doesn't have a simple programmatic API
            // This is a workaround using EnergyCalculator + simple gradient descent
            // TODO: Implement proper optimization wrapper

            py::print("Warning: optimize_geometry() not yet fully implemented");
            py::print("CurcumaOpt class is CLI-focused and not suitable for Python bindings yet");

            // Return copy of input molecule for now
            return Molecule(mol);
        },
        py::arg("molecule"),
        py::arg("method") = "uff",
        py::arg("max_iter") = 100,
        py::arg("verbosity") = 1,
        R"pbdoc(
            Optimize molecular geometry (PLACEHOLDER - not yet implemented)

            WARNING: This function is not yet fully implemented.
            The CurcumaOpt class is designed for CLI usage and does not
            provide a programmatic API suitable for Python bindings.

            Args:
                molecule: Molecule to optimize
                method: Calculation method (default: "uff")
                max_iter: Maximum iterations (default: 100)
                verbosity: Output level 0-3 (default: 1)

            Returns:
                Optimized molecule (currently returns input unchanged)

            Note:
                This is a placeholder. Use curcuma CLI for optimization:
                $ curcuma -opt input.xyz -method gfn2
        )pbdoc");
}
