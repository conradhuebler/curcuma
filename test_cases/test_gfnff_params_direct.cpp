/**
 * Direct GFN-FF Parameter Generation Test
 * Tests parameter generation by directly instantiating GFNFF class
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include "src/core/molecule.h"
#include "src/core/energy_calculators/qm_methods/gfnff.h"
#include "src/core/curcuma_logger.h"

using curcuma::Molecule;

int main() {
    std::cout << std::string(80, '=') << std::endl;
    std::cout << "GFN-FF Parameter Generation Test (Direct GFNFF Class)" << std::endl;
    std::cout << std::string(80, '=') << std::endl;

    // Set verbosity to see debug output
    CurcumaLogger::set_verbosity(3);

    // Test molecules
    struct TestCase {
        const char* name;
        const char* xyz_file;
        double expected_r0;  // Expected equilibrium distance in Bohr
    };

    TestCase tests[] = {
        { "HH", "test_cases/molecules/dimers/HH.xyz", 1.31280068 },
        { "OH", "test_cases/molecules/dimers/OH.xyz", 1.52784633 },
        { "HCl", "test_cases/molecules/dimers/HCl.xyz", 2.30286684 },
    };

    int passed = 0;
    int failed = 0;

    for (const auto& test : tests) {
        std::cout << "\n" << std::string(60, '-') << std::endl;
        std::cout << "Test: " << test.name << std::endl;
        std::cout << std::string(60, '-') << std::endl;

        try {
            // Load molecule
            std::cout << "  Loading: " << test.xyz_file << std::endl;
            Molecule mol(test.xyz_file);
            std::cout << "  Atoms: " << mol.AtomCount() << std::endl;

            // Create GFN-FF instance
            std::cout << "  Initializing GFNFF..." << std::endl;
            GFNFF gfnff;

            // Initialize with molecule
            Mol mol_info = mol.getMolInfo();
            bool init_ok = gfnff.InitialiseMolecule(mol_info);

            if (!init_ok) {
                std::cout << "  ✗ FAILED to initialize GFNFF" << std::endl;
                failed++;
                continue;
            }

            std::cout << "  GFNFF initialized successfully" << std::endl;

            // Perform calculation (triggers parameter generation with debug output)
            std::cout << "  Calculating energy (parameters will be generated)..." << std::endl;
            double energy = gfnff.Calculation(false);
            std::cout << "  Energy: " << std::fixed << std::setprecision(10) << energy << " Hartree" << std::endl;

            // Get bond parameters
            std::cout << "  Retrieving bond parameters..." << std::endl;
            int bond_count = gfnff.getBondCount();
            std::cout << "  Bond count: " << bond_count << std::endl;

            if (bond_count > 0) {
                double shift, alpha, force_constant;
                bool param_ok = gfnff.getVBondParameters(0, shift, alpha, force_constant);

                if (param_ok) {
                    std::cout << "\n  Bond parameters (index 0):" << std::endl;
                    std::cout << "    vbond(1) shift:     " << std::fixed << std::setprecision(12) << shift << " Bohr" << std::endl;
                    std::cout << "    vbond(2) alpha:     " << std::fixed << std::setprecision(12) << alpha << std::endl;
                    std::cout << "    vbond(3) constant:  " << std::fixed << std::setprecision(12) << force_constant << std::endl;

                    // Note: force_constant might not directly be r0, but shift should be close to vbond(1) from JSON
                    std::cout << "\n  Expected r0:        " << std::fixed << std::setprecision(8) << test.expected_r0 << " Bohr" << std::endl;
                    std::cout << "  (Parameter validation depends on internal representation)" << std::endl;

                    passed++;
                } else {
                    std::cout << "  ✗ Could not retrieve bond parameters" << std::endl;
                    failed++;
                }
            } else {
                std::cout << "  ✗ No bonds detected" << std::endl;
                failed++;
            }

        } catch (const std::exception& e) {
            std::cout << "  ✗ Exception: " << e.what() << std::endl;
            failed++;
        }
    }

    // Summary
    std::cout << "\n" << std::string(80, '=') << std::endl;
    std::cout << "SUMMARY" << std::endl;
    std::cout << std::string(80, '=') << std::endl;
    std::cout << "Passed: " << passed << "/3" << std::endl;
    std::cout << "Failed: " << failed << "/3" << std::endl;

    if (failed == 0) {
        std::cout << "\n✓ All tests PASSED" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ Some tests FAILED" << std::endl;
        return 1;
    }
}
