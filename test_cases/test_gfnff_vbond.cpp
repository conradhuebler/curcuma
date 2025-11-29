/**
 * GFN-FF vbond Parameter Verification Test
 * Tests that Curcuma's native cgfnff implementation generates correct vbond parameters
 * that match the reference Fortran GFN-FF implementation.
 *
 * Reference: XTB 6.6.1 (8d0f1dd) GFN-FF topology parameters
 * Purpose: Verify exact parameter generation matches Fortran reference
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "json.hpp"
#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/qm_methods/gfnff.h"

// Using declarations
using curcuma::Molecule;

/**
 * Reference vbond parameters from XTB 6.6.1 GFN-FF
 * Extracted from test_cases/molecules/dimers/HH.json and larger molecules
 *
 * vbond format: [shift, alpha, force_constant]
 */
struct VbondReference {
    std::string molecule_file;     // Path to .xyz file
    std::string xtb_reference;     // Source .json file (for documentation)

    // vbond parameters for a single bond
    double shift;                  // vbond(1,i) - equilibrium distance shift
    double alpha;                  // vbond(2,i) - exponential decay parameter
    double force_constant;         // vbond(3,i) - force constant parameter

    // Tolerance
    double tolerance;

    // For molecules with multiple bonds, we can have multiple references
    int bond_index;                // Index of bond in molecule (0-based)
};

// For molecules with multiple bonds, we can have multiple references per molecule
struct MoleculeVbondReferences {
    std::string molecule_file;     // Path to .xyz file
    std::string xtb_reference;     // Source .json file (for documentation)
    std::vector<VbondReference> bonds;  // All bond references for this molecule
};

class GFNFFVbondTester {
private:
    std::vector<MoleculeVbondReferences> m_references;
    int m_tests_passed = 0;
    int m_tests_failed = 0;

    void setupReferenceValues() {
        // Reference data from XTB 6.6.1 GFN-FF topology parameters
        // File: test_cases/molecules/dimers/HH.json
        MoleculeVbondReferences hh_ref;
        hh_ref.molecule_file = "dimers/HH.xyz";
        hh_ref.xtb_reference = "dimers/HH.json";

        VbondReference hh_bond;
        hh_bond.molecule_file = "dimers/HH.xyz";
        hh_bond.xtb_reference = "dimers/HH.json";
        hh_bond.shift = -0.160000000000000;        // vbond(1,1) shift
        hh_bond.alpha = 0.467792780000000;         // vbond(2,1) alpha
        hh_bond.force_constant = -0.178827447071212; // vbond(3,1) force_constant
        hh_bond.tolerance = 1e-8;                  // very strict tolerance for parameter comparison
        hh_bond.bond_index = 0;                    // first (and only) bond

        hh_ref.bonds.push_back(hh_bond);
        m_references.push_back(hh_ref);

        // Reference data for CH4 from XTB 6.6.1 GFN-FF topology parameters
        // File: test_cases/molecules/larger/CH4.json
        MoleculeVbondReferences ch4_ref;
        ch4_ref.molecule_file = "larger/CH4.xyz";
        ch4_ref.xtb_reference = "larger/CH4.json";

        // CH4 has 4 C-H bonds with identical parameters
        // From CH4.json: vbond = [[-0.182000000000000, 0.482285756225000, -0.167145660211808], ...]
        for (int i = 0; i < 4; i++) {
            VbondReference ch4_bond;
            ch4_bond.molecule_file = "larger/CH4.xyz";
            ch4_bond.xtb_reference = "larger/CH4.json";
            ch4_bond.shift = -0.182000000000000;        // vbond(1,i) shift
            ch4_bond.alpha = 0.482285756225000;         // vbond(2,i) alpha
            ch4_bond.force_constant = -0.167145660211808; // vbond(3,i) force_constant
            ch4_bond.tolerance = 1e-8;                  // very strict tolerance for parameter comparison
            ch4_bond.bond_index = i;                    // bond index

            ch4_ref.bonds.push_back(ch4_bond);
        }

        m_references.push_back(ch4_ref);
    }

    bool fileExists(const std::string& path) {
        std::ifstream f(path);
        return f.good();
    }

    void assertClose(double expected, double actual, double tolerance,
                     const std::string& name, const std::string& component) {
        double error = std::abs(actual - expected);
        double rel_error = (expected != 0.0) ? std::abs((actual - expected) / expected) * 100.0 : 0.0;

        if (error <= tolerance) {
            m_tests_passed++;
            std::cout << "  ✓ " << component << " (" << name << "): "
                      << std::fixed << std::setprecision(12) << actual
                      << " (error: " << std::setprecision(2) << rel_error << "%)"
                      << std::endl;
        } else {
            m_tests_failed++;
            std::cout << "  ✗ " << component << " (" << name << "): "
                      << std::fixed << std::setprecision(12) << actual
                      << " (expected: " << expected << ", "
                      << "error: " << error << " = " << rel_error << "%)"
                      << std::endl;
        }
    }

public:
    GFNFFVbondTester() {
        setupReferenceValues();
    }

    void runTests(const std::string& molecule_base_path) {
        std::cout << "\n====================================================\n"
                  << "   GFN-FF vbond Parameter Verification Test\n"
                  << "====================================================\n";

        int ref_index = 0;
        for (const auto& mol_ref : m_references) {
            ref_index++;
            std::string mol_path = molecule_base_path + "/" + mol_ref.molecule_file;

            std::cout << "\n[Test " << ref_index << "] " << mol_ref.molecule_file
                      << " (XTB: " << mol_ref.xtb_reference << ")\n";
            std::cout << "  ---\n";

            // Check if molecule file exists
            if (!fileExists(mol_path)) {
                std::cout << "  ✗ SKIPPED: Molecule file not found: " << mol_path
                          << std::endl;
                m_tests_failed += mol_ref.bonds.size();
                continue;
            }

            try {
                // Load molecule
                Molecule molecule(mol_path);
                std::cout << "  Atoms: " << molecule.AtomCount() << "\n";

                // Create GFNFF calculator directly to access parameters
                json config = json::object();
                config["verbosity"] = 0;  // Silent mode

                GFNFF gfnff(config);

                // Convert curcuma::Molecule to Mol for GFNFF
                Mol mol = molecule.getMolInfo();
                bool init_success = gfnff.InitialiseMolecule(mol);

                if (!init_success) {
                    std::cout << "  ✗ ERROR: Failed to initialize GFNFF for "
                              << mol_ref.molecule_file << std::endl;
                    m_tests_failed += mol_ref.bonds.size();
                    continue;
                }

                // Calculate energy to ensure parameters are fully generated
                double energy = gfnff.Calculation(false);
                std::cout << "  Energy: " << energy << " Eh\n";

                // Access the vbond parameters for each bond
                int bond_count = gfnff.getBondCount();
                std::cout << "  Bond count: " << bond_count << "\n";

                if (bond_count == 0 && mol_ref.bonds.size() > 0) {
                    std::cout << "  ✗ ERROR: No bonds detected but expected "
                              << mol_ref.bonds.size() << " bonds\n";
                    m_tests_failed += mol_ref.bonds.size();
                    continue;
                }

                // Test each reference bond
                for (const auto& ref : mol_ref.bonds) {
                    if (ref.bond_index >= bond_count) {
                        std::cout << "  ✗ ERROR: Bond index " << ref.bond_index
                                  << " out of range (0-" << (bond_count-1) << ")\n";
                        m_tests_failed++;
                        continue;
                    }

                    double shift, alpha, force_constant;
                    bool param_success = gfnff.getVBondParameters(ref.bond_index, shift, alpha, force_constant);

                    if (param_success) {
                        std::cout << "  Bond " << ref.bond_index << " calculated vbond parameters:\n";
                        std::cout << "    Shift (vbond[1]): " << std::fixed << std::setprecision(12) << shift << "\n";
                        std::cout << "    Alpha (vbond[2]): " << std::fixed << std::setprecision(12) << alpha << "\n";
                        std::cout << "    Force Constant (vbond[3]): " << std::fixed << std::setprecision(12) << force_constant << "\n";

                        // Calculate the ratio between calculated and reference values
                        double ratio = force_constant / ref.force_constant;
                        std::cout << "    Ratio (calculated/reference): " << std::fixed << std::setprecision(12) << ratio << "\n";

                        std::cout << "  Bond " << ref.bond_index << " reference vbond parameters:\n";
                        std::cout << "    Shift (vbond[1]): " << std::fixed << std::setprecision(12) << ref.shift << "\n";
                        std::cout << "    Alpha (vbond[2]): " << std::fixed << std::setprecision(12) << ref.alpha << "\n";
                        std::cout << "    Force Constant (vbond[3]): " << std::fixed << std::setprecision(12) << ref.force_constant << "\n";

                        // Compare parameters
                        assertClose(ref.shift, shift, ref.tolerance,
                                   ref.molecule_file + " bond " + std::to_string(ref.bond_index), "Shift (vbond[1])");
                        assertClose(ref.alpha, alpha, ref.tolerance,
                                   ref.molecule_file + " bond " + std::to_string(ref.bond_index), "Alpha (vbond[2])");
                        assertClose(ref.force_constant, force_constant, ref.tolerance,
                                   ref.molecule_file + " bond " + std::to_string(ref.bond_index), "Force Constant (vbond[3])");
                    } else {
                        std::cout << "  ✗ ERROR: Failed to retrieve vbond parameters for bond "
                                  << ref.bond_index << "\n";
                        m_tests_failed++;
                    }
                }

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                m_tests_failed += mol_ref.bonds.size();
            }
        }

        // Print summary
        std::cout << "\n====================================================\n"
                  << "   Test Summary\n"
                  << "====================================================\n"
                  << "Passed: " << m_tests_passed << "\n"
                  << "Failed: " << m_tests_failed << "\n"
                  << "Total:  " << (m_tests_passed + m_tests_failed) << "\n";

        if (m_tests_failed == 0) {
            std::cout << "\n✓ All tests passed!\n";
        } else {
            std::cout << "\n✗ Some tests failed!\n";
        }
    }

    int getFailureCount() const {
        return m_tests_failed;
    }
};

int main(int argc, char** argv) {
    // Set up logging
    CurcumaLogger::set_verbosity(0);  // Silent for tests

    // Determine molecule base path
    std::string molecule_base = "test_cases/molecules";
    if (argc > 1) {
        molecule_base = argv[1];
    }

    // Run tests
    GFNFFVbondTester tester;
    tester.runTests(molecule_base);

    // Return exit code
    return (tester.getFailureCount() > 0) ? 1 : 0;
}