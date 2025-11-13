/*
 * <Python Bindings for Molecular Dynamics>
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
 * Claude Generated: Python bindings for molecular dynamics simulations
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "src/capabilities/simplemd.h"
#include "src/core/molecule.h"

namespace py = pybind11;
using namespace curcuma;

/**
 * @brief Python bindings for SimpleMD (molecular dynamics)
 *
 * Exposes molecular dynamics simulation capabilities:
 * - NVE, NVT, NPT ensembles
 * - Velocity Verlet integration
 * - Thermostat and barostat control
 * - Trajectory output (XYZ, VTF formats)
 * - Energy and temperature monitoring
 *
 * Educational design: Simplified MD interface for learning molecular dynamics
 * fundamentals without complex setup procedures.
 *
 * Claude Generated: Complete molecular dynamics Python interface
 */
void bind_simplemd(py::module& m) {
    py::class_<SimpleMD>(m, "MolecularDynamics", R"pbdoc(
        Molecular dynamics simulation engine

        Runs molecular dynamics simulations using velocity Verlet integration
        with support for different ensembles and thermostats.

        Example:
        --------
        >>> mol = Molecule("molecule.xyz")
        >>> md = MolecularDynamics(
        ...     method="gfn2",
        ...     timestep=0.5,  # fs
        ...     temperature=300,  # K
        ...     steps=1000
        ... )
        >>> md.set_molecule(mol)
        >>> md.run()
        >>> final_structure = md.get_molecule()
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
                }
            }

            return new SimpleMD(controller);
        }),
             py::arg("config") = py::dict(),
             R"pbdoc(
                Create molecular dynamics simulation

                Args:
                    config: Configuration dictionary with options:
                        - method: Calculation method (e.g., "gfn2", "uff")
                        - timestep: Integration timestep in fs (default: 0.5)
                        - temperature: Target temperature in K (default: 300)
                        - steps: Number of MD steps (default: 1000)
                        - ensemble: Ensemble type "nvt" or "nve" (default: "nvt")
                        - dump: Trajectory dump frequency (default: 10)
                        - thermostat: Thermostat type (default: "berendsen")
                        - verbosity: Output level 0-3 (default: 1)
             )pbdoc")

        // Main MD methods
        .def("set_molecule", &SimpleMD::setMolecule,
             py::arg("molecule"),
             "Set initial molecular structure")

        .def("run", &SimpleMD::start,
             "Run molecular dynamics simulation")

        .def("get_molecule", &SimpleMD::getMolecule,
             "Get final molecular structure")

        // Trajectory access
        .def("write_trajectory", &SimpleMD::writeTrajectory,
             py::arg("filename"),
             "Write MD trajectory to file")

        // Python special methods
        .def("__repr__", [](const SimpleMD&) {
            return "<Molecular Dynamics Simulator>";
        })
        .def("__str__", [](const SimpleMD&) {
            return "SimpleMD - Molecular Dynamics Engine";
        });

    // Convenience function for simple MD runs
    m.def("run_molecular_dynamics",
        [](Molecule& mol, const std::string& method, double temperature,
           int steps, double timestep, int verbosity) {
            json config;
            config["method"] = method;
            config["temperature"] = temperature;
            config["steps"] = steps;
            config["timestep"] = timestep;
            config["verbosity"] = verbosity;
            config["ensemble"] = "nvt";
            config["dump"] = 10;

            SimpleMD md(config);
            md.setMolecule(mol);
            md.start();

            return md.getMolecule();
        },
        py::arg("molecule"),
        py::arg("method") = "uff",
        py::arg("temperature") = 300.0,
        py::arg("steps") = 1000,
        py::arg("timestep") = 0.5,
        py::arg("verbosity") = 1,
        R"pbdoc(
            Run molecular dynamics simulation (convenience function)

            Args:
                molecule: Initial molecular structure
                method: Calculation method (default: "uff")
                temperature: Temperature in Kelvin (default: 300)
                steps: Number of MD steps (default: 1000)
                timestep: Integration timestep in fs (default: 0.5)
                verbosity: Output level 0-3 (default: 1)

            Returns:
                Final molecular structure

            Example:
                >>> mol = Molecule("water.xyz")
                >>> final = run_molecular_dynamics(
                ...     mol, method="gfn2", temperature=350, steps=5000
                ... )
                >>> final.write_xyz("md_final.xyz")
        )pbdoc");
}
