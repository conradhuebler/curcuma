/**
 * GFN-FF Regression Test Suite
 * Tests native cgfnff implementation against XTB 6.6.1 reference data
 *
 * Reference: XTB 6.6.1 (8d0f1dd) GFN-FF calculations
 * Purpose: Prevent parameter regressions during cgfnff development
 *
 * Test molecules (with XTB reference .out files):
 * - HH.xyz (H2 unoptimized, r=0.47 Å)
 * - HH.opt.xyz (H2 optimized, r=0.78 Å)
 * - OH.xyz (OH radical)
 * - water.xyz (H2O, from trimers/)
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "src/core/energycalculator.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"

// Using declarations
using curcuma::Molecule;

/**
 * Reference energy values from XTB 6.6.1 GFN-FF
 * Extracted from test_cases/molecules/dimers/*.out
 *
 * Each test validates ALL energy components individually:
 * - bond_energy
 * - angle_energy
 * - torsion_energy
 * - repulsion_energy
 * - electrostat_energy
 * - dispersion_energy
 * - total_energy (sum of all components)
 */
struct EnergyReference {
    std::string molecule_file;  // Path to .xyz file
    std::string xtb_reference;  // Source .out file (for documentation)

    // Total energy (sum of all components)
    double total_energy;

    // Individual energy components (Eh)
    double bond_energy;
    double angle_energy;
    double torsion_energy;
    double repulsion_energy;
    double electrostat_energy;  // Coulomb/electrostatic
    double dispersion_energy;
    // Note: HB energy, XB energy, bonded ATM energy, external energy are often zero

    // Tolerance (Eh)
    double tolerance;
};

class GFNFFRegressionTester {
private:
    std::vector<EnergyReference> m_references;
    int m_tests_passed = 0;
    int m_tests_failed = 0;

    void setupReferenceValues() {
        // Reference data from XTB 6.6.1 GFN-FF calculations
        // File: test_cases/molecules/dimers/HH.out
        m_references.push_back({
            "dimers/HH.xyz",
            "dimers/HH.out",
            0.050469788228,        // total energy (Eh)
            -0.164952024621,       // bond energy
            0.0,                   // angle energy
            0.0,                   // torsion energy
            0.215468972820,        // repulsion energy
            0.0,                   // electrostat energy
            -0.000047159971,       // dispersion energy
            1e-6                   // tolerance (very strict for reference)
        });

        // File: test_cases/molecules/dimers/HH_opt.out
        // Note: This is H2 at optimized geometry (r=0.78 Å)
        m_references.push_back({
            "dimers/HH.opt.xyz",
            "dimers/HH_opt.out",
            -0.163854384167,       // total energy (Eh) - NEGATIVE at opt geometry!
            -0.174573883763,       // bond energy
            0.0,                   // angle energy
            0.0,                   // torsion energy
            0.010769311348,        // repulsion energy
            0.0,                   // electrostat energy
            -0.000049811752,       // dispersion energy
            1e-6                   // tolerance
        });

        // File: test_cases/molecules/dimers/OH.out
        m_references.push_back({
            "dimers/OH.xyz",
            "dimers/OH.out",
            -0.233907913429,       // total energy (Eh)
            -0.170910705515,       // bond energy
            0.0,                   // angle energy
            0.0,                   // torsion energy
            0.013572983003,        // repulsion energy
            -0.076511600143,       // electrostat energy (charged radical)
            -0.000058590773,       // dispersion energy
            1e-6                   // tolerance
        });

        // File: test_cases/molecules/trimers/water.out
        m_references.push_back({
            "trimers/water.xyz",
            "trimers/water.out",
            -0.327655852825,       // total energy (Eh)
            -0.268290774250,       // bond energy (2 O-H bonds)
            0.000832786791,        // angle energy (H-O-H angle)
            0.0,                   // torsion energy
            0.027517026985,        // repulsion energy
            -0.087575536936,       // electrostat energy
            -0.000139355416,       // dispersion energy
            1e-6                   // tolerance
        });
    }

    bool fileExists(const std::string& path) {
        std::ifstream f(path);
        return f.good();
    }

    void assertClose(double expected, double actual, double tolerance,
                     const std::string& name, const std::string& component) {
        double error = std::abs(actual - expected);
        double rel_error = std::abs((actual - expected) / expected) * 100.0;

        if (error <= tolerance) {
            m_tests_passed++;
            std::cout << "  ✓ " << component << " (" << name << "): "
                      << actual << " Eh (error: " << rel_error << "%)"
                      << std::endl;
        } else {
            m_tests_failed++;
            std::cout << "  ✗ " << component << " (" << name << "): "
                      << actual << " Eh (expected: " << expected << " Eh, "
                      << "error: " << error << " Eh = " << rel_error << "%)"
                      << std::endl;
        }
    }

public:
    GFNFFRegressionTester() {
        setupReferenceValues();
    }

    void runTests(const std::string& molecule_base_path) {
        std::cout << "\n====================================================\n"
                  << "   GFN-FF Regression Test Suite (XTB Reference)\n"
                  << "====================================================\n";

        int ref_index = 0;
        for (const auto& ref : m_references) {
            ref_index++;
            std::string mol_path = molecule_base_path + "/" + ref.molecule_file;

            std::cout << "\n[Test " << ref_index << "] " << ref.molecule_file
                      << " (XTB: " << ref.xtb_reference << ")\n";
            std::cout << "  ---\n";

            // Check if molecule file exists
            if (!fileExists(mol_path)) {
                std::cout << "  ✗ SKIPPED: Molecule file not found: " << mol_path
                          << std::endl;
                m_tests_failed++;
                continue;
            }

            try {
                // Load molecule
                Molecule molecule(mol_path);
                std::cout << "  Atoms: " << molecule.AtomCount() << "\n";

                // Create EnergyCalculator
                json config = json::object();
                config["method"] = "cgfnff";
                config["verbosity"] = 0;  // Silent mode

                EnergyCalculator ec("cgfnff", config);

                // Convert curcuma::Molecule to Mol for EnergyCalculator
                Mol mol = molecule.getMolInfo();
                ec.setMolecule(mol);

                // Calculate energy (no gradient)
                double energy = ec.CalculateEnergy(false);

                // Extract energy components
                double bond_energy = ec.getBondEnergy();
                double angle_energy = ec.getAngleEnergy();
                double dihedral_energy = ec.getDihedralEnergy();
                double repulsion_energy = ec.getRepulsionEnergy();
                double coulomb_energy = ec.getCoulombEnergy();
                double dispersion_energy = ec.getDispersionEnergy();

                std::cout << "  Energy Component Comparison:\n"
                          << "    XTB Reference              Curcuma              Component\n"
                          << "    " << std::string(60, '-') << "\n";

                // Print reference values and calculated values
                std::cout << "    "
                         << std::setw(18) << ref.bond_energy
                         << std::setw(18) << bond_energy
                         << "  Bond Energy\n"
                         << "    "
                         << std::setw(18) << ref.angle_energy
                         << std::setw(18) << angle_energy
                         << "  Angle Energy\n"
                         << "    "
                         << std::setw(18) << ref.torsion_energy
                         << std::setw(18) << dihedral_energy
                         << "  Torsion Energy\n"
                         << "    "
                         << std::setw(18) << ref.repulsion_energy
                         << std::setw(18) << repulsion_energy
                         << "  Repulsion Energy\n"
                         << "    "
                         << std::setw(18) << ref.electrostat_energy
                         << std::setw(18) << coulomb_energy
                         << "  Electrostat Energy\n"
                         << "    "
                         << std::setw(18) << ref.dispersion_energy
                         << std::setw(18) << dispersion_energy
                         << "  Dispersion Energy\n"
                         << "    "
                         << std::string(60, '-') << "\n"
                         << "    "
                         << std::setw(18) << ref.total_energy
                         << std::setw(18) << energy
                         << "  Total Energy\n";

                // Compare all energy components
                assertClose(ref.bond_energy, bond_energy, ref.tolerance,
                           ref.molecule_file, "Bond Energy");
                assertClose(ref.angle_energy, angle_energy, ref.tolerance,
                           ref.molecule_file, "Angle Energy");
                assertClose(ref.torsion_energy, dihedral_energy, ref.tolerance,
                           ref.molecule_file, "Torsion Energy");
                assertClose(ref.repulsion_energy, repulsion_energy, ref.tolerance,
                           ref.molecule_file, "Repulsion Energy");
                assertClose(ref.electrostat_energy, coulomb_energy, ref.tolerance,
                           ref.molecule_file, "Electrostat Energy");
                assertClose(ref.dispersion_energy, dispersion_energy, ref.tolerance,
                           ref.molecule_file, "Dispersion Energy");
                assertClose(ref.total_energy, energy, ref.tolerance,
                           ref.molecule_file, "Total Energy");

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                m_tests_failed++;
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
    GFNFFRegressionTester tester;
    tester.runTests(molecule_base);

    // Return exit code
    return (tester.getFailureCount() > 0) ? 1 : 0;
}
