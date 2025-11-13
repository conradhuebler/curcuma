/*
 * <Python Bindings for Geometry Optimization>
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
 * Claude Generated: Python bindings for geometry optimization
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "src/capabilities/curcumaopt.h"
#include "src/core/molecule.h"

namespace py = pybind11;
using namespace curcuma;

/**
 * @brief Python bindings for geometry optimization
 *
 * Exposes molecular geometry optimization capabilities:
 * - LBFGS optimizer (default, most efficient)
 * - Multiple convergence criteria
 * - Constrained optimization (distances, angles, dihedrals)
 * - Trajectory output
 * - Support for all QM/MM methods
 *
 * Educational design: Simple interface for geometry optimization with
 * sensible defaults and clear convergence reporting.
 *
 * Claude Generated: Complete geometry optimization Python interface
 */
void bind_optimizer(py::module& m) {
    py::class_<CurcumaOpt>(m, "Optimizer", R"pbdoc(
        Molecular geometry optimization

        Optimizes molecular geometries using various algorithms (primarily LBFGS)
        and supports constrained optimization.

        Example:
        --------
        >>> mol = Molecule("initial_structure.xyz")
        >>> opt = Optimizer(method="gfn2", max_iter=100)
        >>> opt.set_molecule(mol)
        >>> opt.optimize()
        >>> optimized = opt.get_molecule()
        >>> optimized.write_xyz("optimized.xyz")
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

            return new CurcumaOpt(controller);
        }),
             py::arg("config") = py::dict(),
             R"pbdoc(
                Create geometry optimizer

                Args:
                    config: Configuration dictionary with options:
                        - method: Calculation method (e.g., "gfn2", "uff")
                        - max_iter: Maximum optimization steps (default: 100)
                        - grad_conv: Gradient convergence threshold (default: 1e-4)
                        - energy_conv: Energy convergence threshold (default: 1e-6)
                        - write_traj: Write optimization trajectory (default: False)
                        - verbosity: Output level 0-3 (default: 1)
             )pbdoc")

        // Main optimization methods
        .def("set_molecule", &CurcumaOpt::setMolecule,
             py::arg("molecule"),
             "Set molecule to optimize")

        .def("optimize", &CurcumaOpt::start,
             "Run geometry optimization")

        .def("get_molecule", &CurcumaOpt::getMolecule,
             "Get optimized molecule")

        // Convergence information
        .def("converged", &CurcumaOpt::Converged,
             "Check if optimization converged")

        .def("final_energy", &CurcumaOpt::FinalEnergy,
             "Get final optimized energy (Hartree)")

        .def("iterations", &CurcumaOpt::Iterations,
             "Get number of optimization iterations performed")

        // Python special methods
        .def("__repr__", [](const CurcumaOpt&) {
            return "<Geometry Optimizer>";
        })
        .def("__str__", [](const CurcumaOpt&) {
            return "Molecular Geometry Optimizer (LBFGS)";
        });

    // Convenience function for simple optimization
    m.def("optimize_geometry",
        [](Molecule& mol, const std::string& method, int max_iter = 100, int verbosity = 1) {
            json config;
            config["method"] = method;
            config["max_iter"] = max_iter;
            config["verbosity"] = verbosity;
            config["write_traj"] = false;

            CurcumaOpt opt(config);
            opt.setMolecule(mol);
            opt.start();

            return opt.getMolecule();
        },
        py::arg("molecule"),
        py::arg("method") = "uff",
        py::arg("max_iter") = 100,
        py::arg("verbosity") = 1,
        R"pbdoc(
            Optimize molecular geometry (convenience function)

            Args:
                molecule: Molecule to optimize
                method: Calculation method (default: "uff")
                max_iter: Maximum iterations (default: 100)
                verbosity: Output level 0-3 (default: 1)

            Returns:
                Optimized molecule

            Example:
                >>> mol = Molecule("input.xyz")
                >>> optimized = optimize_geometry(mol, "gfn2", max_iter=200)
                >>> optimized.write_xyz("optimized.xyz")
        )pbdoc");
}
