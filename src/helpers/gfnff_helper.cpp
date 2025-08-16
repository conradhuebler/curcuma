/*
 * < GFN-FF External Library Helper >
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 * Claude Generated: External GFN-FF library standalone test program
 */

#ifdef USE_GFNFF
#include "external/gfnff/include/gfnff_interface_c.h"
#endif

#include <iostream>
#include <math.h>
#include <stdio.h>

int main(int argc, char** argv)
{
#ifdef USE_GFNFF
    std::cout << "=== External GFN-FF Library Standalone Test ===" << std::endl;

    // Test molecule: H2 (simple case)
    const int natoms = 2;
    int attyp[2] = { 1, 1 }; // H, H
    const double charge = 0.0;
    const int printlevel = 3; // Maximum verbosity for debugging
    const char* solvent = "none";

    // H2 coordinates in Bohr (0.74 Angstrom = 1.398 Bohr)
    double coord[3 * 2] = {
        0.0, 0.0, 0.0, // H1
        0.0, 0.0, 1.398397 // H2 at 0.74 Angstrom
    };

    double (*coord_ptr)[3] = reinterpret_cast<double (*)[3]>(coord);

    std::cout << "Test setup:" << std::endl;
    std::cout << "  Atoms: " << natoms << " (H2 molecule)" << std::endl;
    std::cout << "  Charge: " << charge << std::endl;
    std::cout << "  Print level: " << printlevel << std::endl;
    std::cout << "  H-H distance: 0.74 Angstrom" << std::endl;
    std::cout << std::endl;

    // Initialize calculator
    std::cout << "1. Initializing GFN-FF calculator..." << std::endl;
    c_gfnff_calculator calculator = c_gfnff_calculator_init(
        natoms,
        attyp,
        coord_ptr,
        charge,
        printlevel,
        solvent);

    if (calculator.ptr == nullptr) {
        std::cerr << "ERROR: Failed to initialize GFN-FF calculator" << std::endl;
        return 1;
    }

    std::cout << "   ✓ GFN-FF calculator initialized successfully" << std::endl;
    std::cout << "   Calculator pointer: " << calculator.ptr << std::endl;
    std::cout << std::endl;

    // Perform single-point calculation
    std::cout << "2. Running single-point calculation..." << std::endl;
    double energy = 0.0;
    double gradient[3 * 2];
    double (*grad_ptr)[3] = reinterpret_cast<double (*)[3]>(gradient);
    int iostat = 0;

    c_gfnff_calculator_singlepoint(&calculator, natoms, attyp, coord_ptr, &energy, grad_ptr, &iostat);

    std::cout << "3. Results:" << std::endl;
    std::cout << "   I/O Status: " << iostat << std::endl;

    if (iostat == 0) {
        std::cout << "   ✓ Calculation completed successfully" << std::endl;
        std::cout << "   Energy: " << energy << " Hartree" << std::endl;
        std::cout << "   Gradient:" << std::endl;
        for (int i = 0; i < natoms; ++i) {
            std::cout << "     Atom " << (i + 1) << ": ["
                      << gradient[3 * i + 0] << ", "
                      << gradient[3 * i + 1] << ", "
                      << gradient[3 * i + 2] << "]" << std::endl;
        }
    } else {
        std::cerr << "   ERROR: Calculation failed with iostat = " << iostat << std::endl;
    }

    // Cleanup
    std::cout << std::endl;
    std::cout << "4. Cleanup..." << std::endl;
    c_gfnff_calculator_deallocate(&calculator);
    std::cout << "   ✓ GFN-FF calculator deallocated" << std::endl;

    std::cout << std::endl;
    std::cout << "=== External GFN-FF Test Completed ===" << std::endl;

    return (iostat == 0) ? 0 : 1;

#else
    std::cerr << "ERROR: External GFN-FF not available - USE_GFNFF not compiled" << std::endl;
    return 1;
#endif
}