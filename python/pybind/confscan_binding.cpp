/*
 * <Python Bindings for Conformational Scanning - SIMPLIFIED>
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
 * Claude Generated: Simplified Python bindings - ConfScan is CLI-focused
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "src/capabilities/confscan.h"
#include "src/core/molecule.h"

namespace py = pybind11;
using namespace curcuma;

/**
 * @brief Simplified Python bindings for conformational analysis
 *
 * NOTE: ConfScan is designed for CLI usage and writes results to files.
 * Full programmatic API not yet available.
 *
 * Claude Generated: Minimal bindings matching actual API
 */
void bind_confscan(py::module& m) {
    // Convenience function placeholder
    m.def("scan_dihedral",
        [](const Molecule& mol, const std::vector<int>& dihedral_atoms,
           const std::string& method, int steps, int verbosity) {

            py::print("Warning: scan_dihedral() not yet fully implemented");
            py::print("ConfScan class is CLI-focused and not suitable for Python bindings yet");

            // Return empty list for now
            return std::vector<Molecule>();
        },
        py::arg("molecule"),
        py::arg("dihedral_atoms"),
        py::arg("method") = "uff",
        py::arg("steps") = 36,
        py::arg("verbosity") = 1,
        R"pbdoc(
            Scan dihedral angle and generate conformers (PLACEHOLDER - not yet implemented)

            WARNING: This function is not yet fully implemented.
            The ConfScan class is designed for CLI usage and does not
            provide a programmatic API suitable for Python bindings.

            Args:
                molecule: Molecule to scan
                dihedral_atoms: List of 4 atom indices defining dihedral (0-indexed)
                method: Calculation method (default: "uff")
                steps: Number of scan steps (default: 36 = 10° increments)
                verbosity: Output level 0-3 (default: 1)

            Returns:
                List of conformer molecules (currently returns empty list)

            Note:
                This is a placeholder. Use curcuma CLI for conformational scanning:
                $ curcuma -confscan input.xyz -method gfn2
        )pbdoc");
}
