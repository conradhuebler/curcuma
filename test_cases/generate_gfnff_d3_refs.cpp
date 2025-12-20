/**
 * GFN-FF D3 Reference Energy Generator
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Purpose: Generates reference D3 dispersion energies for GFN-FF validation tests.
 *
 * Method:
 * 1. Load test molecules
 * 2. Calculate D3 dispersion with GFN-FF parameters (s6=1.0, s8=2.85, a1=0.80, a2=4.60)
 * 3. Output energies as C++ code for copy-paste into test_gfnff_d3.cpp
 *
 * Claude Generated - December 19, 2025
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

// Core
#include "src/core/molecule.h"
#include "src/core/config_manager.h"
#include "src/tools/formats.h"

// D3
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"

#include "json.hpp"
using json = nlohmann::json;

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct MoleculeReference {
    std::string file_path;
    std::string description;
};

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    std::cout << "\n";
    std::cout << "GFN-FF D3 Reference Energy Generator\n";
    std::cout << "================================================================================\n";
    std::cout << "Generating reference D3 energies using GFN-FF parameters:\n";
    std::cout << "  s6  = 1.0\n";
    std::cout << "  s8  = 2.85\n";
    std::cout << "  a1  = 0.80\n";
    std::cout << "  a2  = 4.60 Bohr\n";
    std::cout << "  alp = 14.0\n";
    std::cout << "\n";

    // Test molecules (same as in test_gfnff_d3.cpp)
    std::vector<MoleculeReference> molecules = {
        {"../test_cases/molecules/dimers/HH.xyz", "H2 dimer"},
        {"../test_cases/molecules/dimers/HCl.xyz", "HCl dimer"},
        {"../test_cases/molecules/dimers/OH.xyz", "OH radical"},
        {"../test_cases/molecules/larger/CH4.xyz", "Methane"},
        {"../test_cases/molecules/larger/CH3OH.xyz", "Methanol"},
        {"../test_cases/molecules/larger/CH3OCH3.xyz", "Dimethyl ether"}
    };

    // GFN-FF D3 parameters (Spicher/Grimme, J. Chem. Theory Comput. 2020)
    json d3_config_json = {
        {"d3_s6", 1.0},
        {"d3_s8", 2.85},
        {"d3_a1", 0.80},
        {"d3_a2", 4.60},
        {"d3_alp", 14.0}  // Standard D3 alpha parameter
    };

    std::cout << "Processing " << molecules.size() << " molecules...\n\n";

    int success_count = 0;

    for (const auto& mol_ref : molecules) {
        std::cout << "Processing: " << mol_ref.description << " (" << mol_ref.file_path << ")\n";

        try {
            // Load molecule
            auto mol = Files::LoadFile(mol_ref.file_path);
            if (mol.AtomCount() == 0) {
                std::cerr << "  ERROR: Failed to load molecule\n\n";
                continue;
            }

            std::vector<int> atoms = mol.Atoms();
            Matrix geometry = mol.getGeometry();

            std::cout << "  Atoms: " << atoms.size() << "\n";

            // Calculate D3 dispersion with GFN-FF parameters
            ConfigManager d3_config("d3param", d3_config_json);
            D3ParameterGenerator d3_gen(d3_config);

            d3_gen.GenerateParameters(atoms, geometry);
            double d3_energy = d3_gen.getTotalEnergy();

            std::cout << "  D3 Energy: " << std::scientific << std::setprecision(12)
                      << d3_energy << " Eh\n";
            std::cout << "  ✓ Success\n\n";

            success_count++;

        } catch (const std::exception& e) {
            std::cerr << "  ERROR: " << e.what() << "\n\n";
        }
    }

    std::cout << "================================================================================\n";
    std::cout << "Summary: " << success_count << "/" << molecules.size() << " molecules processed\n";
    std::cout << "\n";

    // Now output C++ code for copy-paste
    if (success_count == molecules.size()) {
        std::cout << "================================================================================\n";
        std::cout << "C++ CODE FOR test_gfnff_d3.cpp (copy-paste into createReferenceData())\n";
        std::cout << "================================================================================\n\n";

        for (const auto& mol_ref : molecules) {
            try {
                auto mol = Files::LoadFile(mol_ref.file_path);
                std::vector<int> atoms = mol.Atoms();
                Matrix geometry = mol.getGeometry();

                ConfigManager d3_config("d3param", d3_config_json);
                D3ParameterGenerator d3_gen(d3_config);
                d3_gen.GenerateParameters(atoms, geometry);
                double d3_energy = d3_gen.getTotalEnergy();

                // Output as C++ struct initializer
                std::cout << "    refs.push_back({\n";
                std::cout << "        \"" << mol_ref.file_path << "\",\n";
                std::cout << "        \"" << mol_ref.description << "\",\n";
                std::cout << "        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,\n";
                std::cout << "        " << std::scientific << std::setprecision(12)
                          << d3_energy << ",  // expected_d3_energy\n";
                std::cout << "        REFERENCE_TOLERANCE, CONSISTENCY_TOLERANCE,\n";
                std::cout << "        -999.0  // xtb_gfnff_dispersion - not available\n";
                std::cout << "    });\n\n";

            } catch (const std::exception& e) {
                std::cerr << "ERROR: " << e.what() << "\n";
            }
        }

        std::cout << "\n";
    }

    return (success_count == molecules.size()) ? 0 : 1;
}
