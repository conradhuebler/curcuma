/*
 * <Python Bindings for RMSD Analysis>
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
 * Claude Generated: Python bindings for RMSD calculation and alignment
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "src/capabilities/rmsd.h"
#include "src/core/molecule.h"

namespace py = pybind11;
using namespace curcuma;

/**
 * @brief Python bindings for RMSD (Root Mean Square Deviation) analysis
 *
 * Exposes RMSD calculation and structural alignment capabilities:
 * - Calculate RMSD between two structures
 * - Align molecules to minimize RMSD
 * - Heavy atom RMSD (ignore hydrogens)
 * - Fragment-based RMSD
 * - Reorder atoms to match reference structure
 *
 * Educational design: Simple interface for structural comparison and alignment,
 * essential for conformational analysis and molecular dynamics.
 *
 * Claude Generated: Complete RMSD analysis Python interface
 */
void bind_rmsd(py::module& m) {
    py::class_<RMSDDriver>(m, "RMSD", R"pbdoc(
        RMSD calculation and structural alignment

        Provides tools for calculating root-mean-square deviation between
        molecular structures and performing optimal structural alignment.

        Example:
        --------
        >>> mol1 = Molecule("structure1.xyz")
        >>> mol2 = Molecule("structure2.xyz")
        >>> rmsd = RMSD(mol1, mol2)
        >>> rmsd_value = rmsd.calculate()
        >>> print(f"RMSD: {rmsd_value:.3f} Angstrom")
        >>> aligned = rmsd.get_aligned_molecule()
    )pbdoc")
        // Constructors
        .def(py::init([](const py::dict& config_dict) {
            json controller = json::object();
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

            return new RMSDDriver(controller);
        }),
             py::arg("config") = py::dict(),
             "Create RMSD calculator with configuration")

        .def(py::init([]() {
            json controller = json::object();
            controller["heavy"] = false;
            controller["reorder"] = false;
            controller["check"] = -1;
            controller["silent"] = false;
            return new RMSDDriver(controller);
        }),
             "Create RMSD calculator with default configuration")

        // Main methods
        .def("set_reference", &RMSDDriver::setReference,
             py::arg("molecule"),
             "Set reference structure for RMSD calculation")

        .def("set_target", &RMSDDriver::setTarget,
             py::arg("molecule"),
             "Set target structure for RMSD calculation")

        .def("calculate", &RMSDDriver::CalculateRMSD,
             "Calculate RMSD between reference and target structures (Ångström)")

        .def("rmsd", &RMSDDriver::RMSD,
             "Get calculated RMSD value (Ångström)")

        .def("get_aligned_target", &RMSDDriver::TargetAligned,
             "Get aligned target molecule")

        .def("get_aligned_reference", &RMSDDriver::ReferenceAligned,
             "Get aligned reference molecule")

        // Python-friendly helper methods
        .def("__call__", &RMSDDriver::CalculateRMSD,
             "Calculate RMSD (callable interface)")

        .def("__repr__", [](const RMSDDriver&) {
            return "<RMSD Calculator>";
        })
        .def("__str__", [](const RMSDDriver&) {
            return "RMSD Calculator for structural alignment";
        });

    // Convenience functions for simple RMSD calculations
    m.def("calculate_rmsd",
        [](const Molecule& ref, const Molecule& target, bool heavy_only = false, bool reorder = false) {
            json config;
            config["heavy"] = heavy_only;
            config["reorder"] = reorder;
            config["silent"] = true;

            RMSDDriver rmsd(config);
            rmsd.setReference(ref);
            rmsd.setTarget(target);
            rmsd.CalculateRMSD();
            return rmsd.RMSD();
        },
        py::arg("reference"),
        py::arg("target"),
        py::arg("heavy_only") = false,
        py::arg("reorder") = false,
        R"pbdoc(
            Calculate RMSD between two molecules

            Args:
                reference: Reference molecule
                target: Target molecule to compare
                heavy_only: Only consider heavy atoms (default: False)
                reorder: Allow atom reordering to minimize RMSD (default: False)

            Returns:
                RMSD value in Ångström

            Example:
                >>> mol1 = Molecule("conf1.xyz")
                >>> mol2 = Molecule("conf2.xyz")
                >>> rmsd = calculate_rmsd(mol1, mol2, heavy_only=True)
                >>> print(f"Heavy atom RMSD: {rmsd:.3f} Å")
        )pbdoc");

    m.def("align_molecules",
        [](const Molecule& ref, const Molecule& target) {
            json config;
            config["silent"] = true;
            config["reorder"] = false;

            RMSDDriver rmsd(config);
            rmsd.setReference(ref);
            rmsd.setTarget(target);
            rmsd.CalculateRMSD();  // This performs alignment
            return rmsd.TargetAligned();
        },
        py::arg("reference"),
        py::arg("target"),
        R"pbdoc(
            Align target molecule to reference structure

            Args:
                reference: Reference molecule
                target: Target molecule to align

            Returns:
                Aligned copy of target molecule

            Example:
                >>> ref = Molecule("reference.xyz")
                >>> mol = Molecule("structure.xyz")
                >>> aligned = align_molecules(ref, mol)
                >>> aligned.write_xyz("aligned.xyz")
        )pbdoc");
}
