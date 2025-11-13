/*
 * <Python Bindings for Molecule Class>
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
 * Claude Generated: Python bindings for core Molecule class
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>

#include "src/core/molecule.h"

namespace py = pybind11;
using namespace curcuma;

/**
 * @brief Python bindings for the Molecule class
 *
 * Exposes the core Molecule data structure to Python, allowing:
 * - File I/O (XYZ, MOL2, SDF formats)
 * - Geometry manipulation (coordinates, atoms, bonds)
 * - Property access (charge, mass, energy, fragments)
 * - Structure analysis (distance matrices, connectivity)
 * - JSON import/export
 *
 * Educational design: Methods are exposed with Python-friendly names
 * and intuitive behavior (e.g., __len__ for atom count, __str__ for printing)
 *
 * Claude Generated: Complete Molecule class Python interface
 */
void bind_molecule(py::module& m) {
    py::class_<Molecule>(m, "Molecule", R"pbdoc(
        Molecular structure representation

        The Molecule class is the core data structure for storing and manipulating
        molecular geometries, properties, and metadata. It supports both quantum
        mechanical and force field calculations.

        Constructors:
            Molecule(): Empty molecule for programmatic construction
            Molecule(n, q=0): Pre-allocate n atoms with zero coordinates (for file loading)
            Molecule(filename): Load molecule from file (XYZ, MOL2, SDF)
            Molecule(other): Copy constructor

        Example:
            >>> mol = Molecule("water.xyz")
            >>> print(mol.atom_count())
            3
            >>> print(mol.charge())
            0
            >>> mol.print_geometry()
    )pbdoc")
        // Constructors
        .def(py::init<>(), "Create an empty molecule")
        .def(py::init<int, int>(),
             py::arg("natoms"),
             py::arg("charge") = 0,
             "Pre-allocate molecule with n atoms (zero coordinates)")
        .def(py::init<const std::string&>(),
             py::arg("filename"),
             "Load molecule from file (XYZ, MOL2, SDF)")
        .def(py::init<const Molecule&>(),
             py::arg("other"),
             "Copy constructor")

        // File I/O
        .def("load_geometry", &Molecule::LoadGeometry,
             py::arg("filename"),
             "Load molecular geometry from file")
        .def("append_xyz", &Molecule::appendXYZFile,
             py::arg("filename"),
             "Append XYZ geometry to file")
        .def("write_xyz", &Molecule::writeXYZFile,
             py::arg("filename"),
             "Write molecular geometry to XYZ file")

        // JSON support
        .def("to_json", &Molecule::ExportJson,
             "Export molecule to JSON format")
        .def("from_json", py::overload_cast<const json&>(&Molecule::ImportJson),
             py::arg("json_data"),
             "Import molecule from JSON data")
        .def("write_json", &Molecule::WriteJsonFile,
             py::arg("filename"),
             "Write molecule to JSON file")

        // Basic properties
        .def("atom_count", &Molecule::AtomCount,
             "Get number of atoms")
        .def("charge", &Molecule::Charge,
             "Get molecular charge")
        .def("set_charge", &Molecule::setCharge,
             py::arg("charge"),
             "Set molecular charge")
        .def("mass", &Molecule::Mass,
             "Get molecular mass in amu")
        .def("calculate_mass", &Molecule::CalculateMass,
             "Calculate molecular mass from atomic masses")

        // Energy
        .def("energy", &Molecule::Energy,
             "Get stored energy (Hartree)")
        .def("set_energy", &Molecule::setEnergy,
             py::arg("energy"),
             "Set energy value (Hartree)")

        // Geometry access and manipulation
        .def("geometry", &Molecule::getGeometry,
             "Get geometry matrix (3×N matrix of coordinates in Ångström)")
        .def("set_geometry", py::overload_cast<const Geometry&>(&Molecule::setGeometry),
             py::arg("geometry"),
             "Set geometry from 3×N matrix (Ångström)")
        .def("get_atom", &Molecule::Atom,
             py::arg("index"),
             "Get atom position vector (Ångström)")
        .def("set_atom", &Molecule::setAtom,
             py::arg("position"), py::arg("index"),
             "Set atom position (Ångström)")
        .def("get_element", &Molecule::Atom,
             py::arg("index"),
             "Get element number (atomic number)")

        // Distance calculations
        .def("calculate_distance", &Molecule::CalculateDistance,
             py::arg("i"), py::arg("j"),
             "Calculate distance between atoms i and j (Ångström)")
        .def("distance_matrix", &Molecule::DistanceMatrix,
             "Get full distance matrix (Ångström)")

        // Connectivity and topology
        .def("initialise_connectivity", &Molecule::InitialiseConnectedMass,
             py::arg("scaling") = 1.3, py::arg("protons") = true,
             "Initialize bond connectivity detection")
        .def("print_connectivity", &Molecule::printConnections,
             "Print molecular connectivity information")

        // Fragment analysis
        .def("fragment_count", &Molecule::FragmentCount,
             "Get number of molecular fragments")
        .def("fragment_masses", &Molecule::FragmentMass,
             "Get mass of each fragment")
        .def("print_fragments", &Molecule::printFragmente,
             "Print fragment information")

        // Utility methods
        .def("print_geometry", &Molecule::print_geom,
             py::arg("detailed") = true,
             "Print molecular geometry to console")
        .def("print_atom", &Molecule::printAtom,
             py::arg("index"),
             "Print information about specific atom")

        // Center of mass
        .def("center_of_mass", &Molecule::getCentreOfMass,
             "Get center of mass coordinates (Ångström)")
        .def("center_of_geometry", &Molecule::getCentreOfGeometry,
             "Get geometric center coordinates (Ångström)")

        // Python special methods
        .def("__len__", &Molecule::AtomCount,
             "Number of atoms")
        .def("__repr__", [](const Molecule& mol) {
            return fmt::format("<Molecule: {} atoms, charge={}, mass={:.2f} amu>",
                             mol.AtomCount(), mol.Charge(), mol.Mass());
        })
        .def("__str__", [](const Molecule& mol) {
            return fmt::format("Molecule with {} atoms (charge={}, mass={:.2f} amu)",
                             mol.AtomCount(), mol.Charge(), mol.Mass());
        });
}
