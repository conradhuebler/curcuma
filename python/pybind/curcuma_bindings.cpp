/*
 * <Python Bindings for Curcuma - Main Module>
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
 * Claude Generated: Python interface for Curcuma computational chemistry toolkit
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

// Forward declarations for submodule binding functions
void bind_molecule(py::module& m);
void bind_energy_calculator(py::module& m);
void bind_rmsd(py::module& m);
void bind_optimizer(py::module& m);
void bind_simplemd(py::module& m);
void bind_confscan(py::module& m);

/**
 * @brief Main pybind11 module definition for Curcuma
 *
 * This module exposes the core computational chemistry capabilities of Curcuma
 * to Python, including:
 * - Molecular structure manipulation (Molecule class)
 * - Energy calculations (QM/MM methods via EnergyCalculator)
 * - RMSD analysis and structural alignment
 * - Geometry optimization (multiple algorithms)
 * - Molecular dynamics simulations (SimpleMD)
 * - Conformational analysis (ConfScan)
 *
 * Educational-first design: All bindings prioritize clarity and ease of use
 * over performance optimizations, making them ideal for learning and teaching
 * computational chemistry.
 *
 * Claude Generated: Complete Python interface implementation
 */
PYBIND11_MODULE(curcuma, m) {
    m.doc() = R"pbdoc(
        Curcuma - Molecular Modelling and Simulation Toolkit
        ====================================================

        A Python interface to the Curcuma computational chemistry toolkit,
        providing access to:

        - Quantum mechanical methods (Extended Hückel, GFN1/GFN2, etc.)
        - Force field methods (UFF, GFN-FF, QMDFF)
        - Geometry optimization (LBFGS and other algorithms)
        - Molecular dynamics simulations
        - Conformational analysis and searching
        - RMSD calculations and structural alignment
        - Dispersion corrections (D3, D4)

        Example Usage
        -------------
        >>> import curcuma
        >>> mol = curcuma.Molecule("structure.xyz")
        >>> calc = curcuma.EnergyCalculator("gfn2")
        >>> calc.set_molecule(mol)
        >>> energy = calc.calculate_energy()
        >>> print(f"Energy: {energy} Hartree")

        For more information, see the documentation at:
        https://github.com/conradhuebler/curcuma
    )pbdoc";

    // Bind core molecular data structures
    bind_molecule(m);

    // Bind energy calculation interface
    bind_energy_calculator(m);

    // Bind analysis capabilities
    bind_rmsd(m);

    // Bind optimization capabilities
    bind_optimizer(m);

    // Bind molecular dynamics
    bind_simplemd(m);

    // Bind conformational analysis
    bind_confscan(m);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
