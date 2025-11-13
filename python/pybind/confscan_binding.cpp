/*
 * <Python Bindings for Conformational Scanning>
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
 * Claude Generated: Python bindings for conformational analysis
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "src/capabilities/confscan.h"
#include "src/core/molecule.h"

namespace py = pybind11;
using namespace curcuma;

/**
 * @brief Python bindings for ConfScan (conformational scanning)
 *
 * Exposes conformational analysis capabilities:
 * - Systematic dihedral angle scanning
 * - Energy-based conformer ranking
 * - RMSD-based structure filtering
 * - Trajectory output
 * - Support for all QM/MM methods
 *
 * Educational design: Simple interface for exploring potential energy surfaces
 * and conformational preferences of flexible molecules.
 *
 * Claude Generated: Complete conformational analysis Python interface
 */
void bind_confscan(py::module& m) {
    py::class_<ConfScan>(m, "ConformationalScan", R"pbdoc(
        Conformational scanning and analysis

        Systematically explores conformational space by scanning dihedral angles
        and optimizing structures at each step.

        Example:
        --------
        >>> mol = Molecule("flexible_molecule.xyz")
        >>> scan = ConformationalScan(
        ...     method="gfn2",
        ...     dihedral=[1, 2, 3, 4],  # Atoms defining dihedral
        ...     steps=36  # 10-degree increments
        ... )
        >>> scan.set_molecule(mol)
        >>> scan.run()
        >>> conformers = scan.get_conformers()
    )pbdoc")
        // Constructor
        .def(py::init([](const py::dict& config_dict) {
            json controller = json::object();

            // Convert Python dict to JSON
            for (auto item : config_dict) {
                std::string key = py::str(item.first);
                py::object value = py::reinterpret_borrow<py::object>(item.second);

                if (py::isinstance<py::int_>(value)) {
                    controller[key] = value.cast<int>();
                } else if (py::isinstance<py::float_>(value)) {
                    controller[key] = value.cast<double>();
                } else if (py::isinstance<py::bool_>(value)) {
                    controller[key] = value.cast<bool>();
                } else if (py::isinstance<py::str>(value)) {
                    controller[key] = value.cast<std::string>();
                } else if (py::isinstance<py::list>(value)) {
                    // Handle lists (for dihedral atom indices)
                    controller[key] = json::array();
                    for (auto item_in_list : value) {
                        if (py::isinstance<py::int_>(item_in_list)) {
                            controller[key].push_back(item_in_list.cast<int>());
                        }
                    }
                }
            }

            return new ConfScan(controller);
        }),
             py::arg("config") = py::dict(),
             R"pbdoc(
                Create conformational scanner

                Args:
                    config: Configuration dictionary with options:
                        - method: Calculation method (e.g., "gfn2", "uff")
                        - dihedral: List of 4 atom indices defining dihedral angle
                        - steps: Number of scan steps (default: 36)
                        - optimize: Optimize at each scan point (default: True)
                        - rmsd_threshold: RMSD cutoff for duplicate filtering (default: 0.5 Å)
                        - energy_window: Energy window for conformer selection in kJ/mol (default: 50)
                        - verbosity: Output level 0-3 (default: 1)
             )pbdoc")

        // Main methods
        .def("set_molecule", &ConfScan::setMolecule,
             py::arg("molecule"),
             "Set molecule for conformational scanning")

        .def("run", &ConfScan::start,
             "Run conformational scan")

        .def("get_conformers", &ConfScan::getStructures,
             "Get list of generated conformers")

        .def("write_ensemble", &ConfScan::writeEnsemble,
             py::arg("filename"),
             "Write conformer ensemble to XYZ file")

        // Python special methods
        .def("__repr__", [](const ConfScan&) {
            return "<Conformational Scanner>";
        })
        .def("__str__", [](const ConfScan&) {
            return "Conformational Scanning Engine";
        });

    // Convenience function for simple conformational scans
    m.def("scan_dihedral",
        [](Molecule& mol, const std::vector<int>& dihedral_atoms,
           const std::string& method, int steps, int verbosity) {
            json config;
            config["method"] = method;
            config["steps"] = steps;
            config["verbosity"] = verbosity;
            config["optimize"] = true;

            // Set dihedral atom indices
            config["dihedral"] = json::array();
            for (int atom : dihedral_atoms) {
                config["dihedral"].push_back(atom);
            }

            ConfScan scan(config);
            scan.setMolecule(mol);
            scan.start();

            return scan.getStructures();
        },
        py::arg("molecule"),
        py::arg("dihedral_atoms"),
        py::arg("method") = "uff",
        py::arg("steps") = 36,
        py::arg("verbosity") = 1,
        R"pbdoc(
            Scan dihedral angle and generate conformers

            Args:
                molecule: Molecule to scan
                dihedral_atoms: List of 4 atom indices defining dihedral (0-indexed)
                method: Calculation method (default: "uff")
                steps: Number of scan steps (default: 36 = 10° increments)
                verbosity: Output level 0-3 (default: 1)

            Returns:
                List of conformer molecules

            Example:
                >>> mol = Molecule("butane.xyz")
                >>> conformers = scan_dihedral(
                ...     mol, dihedral_atoms=[0, 1, 2, 3], method="gfn2", steps=72
                ... )
                >>> print(f"Generated {len(conformers)} conformers")
        )pbdoc");
}
