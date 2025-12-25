/**
 * Unified GFN-FF Test Suite
 * Combines regression, parameter, and comprehensive validation tests
 *
 * Reference: XTB 6.6.1 (8d0f1dd) GFN-FF calculations
 * Purpose: Complete validation of Curcuma's native cgfnff implementation
 *
 * Test Categories:
 * 1. Energy Component Regression - Total energy breakdown validation
 * 2. Parameter Verification - vbond, vangl, vtor parameter accuracy
 * 3. Comprehensive Validation - CN, charges, all energy terms (CH3OH)
 *
 * Molecules Tested:
 * - HH (H2 dimer) - Basic bond energy
 * - OH (radical) - Electrostatic with charged species
 * - H2O (trimer) - Angle energy terms
 * - CH3OH (methanol) - All energy terms including torsion
 * - CH4 (methane) - Multiple equivalent bonds
 * - CH3OCH3 (dimethyl ether) - Complex torsion terms
 * - C6H6 (benzene) - Aromatic π-system, sp² carbons (D4 dispersion test)
 * - C4H10 (butane) - Flexible alkane chain (D4 dispersion test)
 *
 * Claude Generated - December 25, 2025
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/tools/formats.h"
#include "json.hpp"

using json = nlohmann::json;
using curcuma::Molecule;

// Test configuration
struct TestConfig {
    double energy_tolerance = 1e-6;     // Total energy tolerance (Hartree)
    double component_tolerance = 1e-6; // Individual component tolerance
    double param_tolerance = 1e-8;      // Parameter tolerance (vbond, etc.)
    double cn_tolerance = 0.05;         // Coordination number tolerance
    double charge_tolerance = 0.005;    // EEQ charge tolerance
    int verbosity = 1;                  // Output level (0=silent, 1=results, 2=details)
};

// Energy component reference data
struct EnergyReference {
    std::string molecule_file;
    std::string description;

    // Energy components from XTB 6.6.1 reference
    double total_energy;
    double bond_energy;
    double angle_energy;
    double torsion_energy;
    double repulsion_energy;
    double electrostat_energy;
    double dispersion_energy;
    double batm_energy;
};

// Parameter reference data (vbond format: [shift, alpha, force_constant])
struct VbondReference {
    std::string molecule_file;
    int bond_index;
    std::string bond_description;  // e.g. "H-H bond"
    double shift;
    double alpha;
    double force_constant;
};

// Comprehensive validation data (CH3OH)
struct ComprehensiveReference {
    std::string molecule_file;
    std::string description;

    // Coordination numbers
    std::vector<std::pair<std::string, double>> coordination_numbers;

    // EEQ charges
    std::vector<std::pair<std::string, double>> charges;

    // All energy components
    std::map<std::string, double> energies;
};

class UnifiedGFNFFTest {
private:
    TestConfig config;
    std::vector<EnergyReference> energy_refs;
    std::vector<VbondReference> vbond_refs;
    std::vector<ComprehensiveReference> comprehensive_refs;

    int total_tests = 0;
    int passed_tests = 0;

public:
    UnifiedGFNFFTest(const TestConfig& cfg = TestConfig()) : config(cfg) {
        setup_reference_data();
    }

    void setup_reference_data() {
        // Energy regression references
        energy_refs = {
            {"test_cases/molecules/dimers/HH.xyz", "H2 dimer (unoptimized, r=0.47 Å)",
             0.050469853027, -0.164952020372, 0.0, 0.0, 0.215469033370, 0.0, -0.000047159971, 0.0},
            {"test_cases/molecules/dimers/HH.opt.xyz", "H2 optimized (r=0.78 Å)",
             -0.0192416941, -0.0237366693, 0.0, 0.0, 0.00294998396, -0.000197294282, 0.00195844759, -0.000015962895},
            {"test_cases/molecules/dimers/OH.xyz", "OH radical",
             -0.233907925071, -0.170910719346, 0.0, 0.0, 0.013572992439, -0.076511607391, -0.000058590773, 0.0},
            {"test_cases/molecules/dimers/HCl.xyz", "H-Cl dimer (halogen chemistry)",
             -0.0120384482, -0.0843104989, 0.0, 0.0, 0.080505700573, -0.008093183022, -0.000140466881, 0.0},
            {"test_cases/molecules/trimers/water.xyz", "H2O trimer",
             -0.2130633276, -0.2159236329, 0.00123926595, 0.0, 0.00702065869, -0.00139899241, 0.00480579676, -0.000806423657},
            {"test_cases/molecules/trimers/O3.xyz", "O3 ozone (bent triatomic, π-system)",
             -0.330939597635, -0.354938759395, 0.004615989890, 0.0, 0.023390924422, -0.003608214874, -0.000399537678, 0.0},
            {"test_cases/molecules/larger/CH3OCH3.xyz", "Dimethyl ether (complex torsions)",
             -1.2092092216, -1.2164439418, 0.001779533537, 0.000023390598, 0.053864662977, -0.047825361074, 0.000041946447, -0.000142398327},
            {"test_cases/molecules/larger/CH4.xyz", "Methane - sp3 tetrahedral validation",
             -0.630814967693, -0.656386349843, 0.000068985137, 0.0, 0.027729233574, -0.001576233642, -0.000650602918, 0.0},
            {"test_cases/molecules/larger/C6H6.xyz", "Benzene - aromatic π-system, sp² carbons",
             -2.364775496871, -2.508166283882, 0.000000000758, 0.000000729656, 0.158023153244, -0.004339043253, -0.006294467601, -0.003999585793},
            {"test_cases/validation/butane.xyz", "Butane - flexible alkane chain",
             -1.950514159158, -2.051987951695, 0.000390623328, 0.004144855655, 0.107248927516, -0.004233455377, -0.005848570553, -0.000228588030}
        };

        // vbond parameter references
        vbond_refs = {
            {"test_cases/molecules/dimers/HH.xyz", 0, "H-H bond", 0.0, 1.889726124565, 0.020104827311},
            {"test_cases/molecules/dimers/HH.opt.xyz", 0, "H-H bond (optimized)", 0.0, 1.889726124565, 0.020104827311},
            {"test_cases/molecules/dimers/OH.xyz", 0, "O-H bond", 0.0, 2.101777932193, 0.037073626719},
            {"test_cases/molecules/larger/CH4.xyz", 0, "C-H bond", 0.0, 1.879923356919, 0.009696223926},
            {"test_cases/molecules/larger/CH4.xyz", 1, "C-H bond (equivalent)", 0.0, 1.879923356919, 0.009696223926}
        };

        // Comprehensive validation references
        comprehensive_refs = {
            {"test_cases/molecules/larger/CH3OH.xyz", "Methanol - comprehensive validation",
             {
                 {"C", 3.48}, {"O", 1.91},
                 {"H_methyl", 0.050}, {"H_hydroxyl", 0.050}
             },
             {
                 {"C", 0.048}, {"O", -0.432}, {"H_avg", 0.050}
             },
             {
                 {"E_total", -0.760203302581}, {"E_bond", -0.742357724674},
                 {"E_angle", 0.000931012526}, {"E_torsion", 0.000000954542},
                 {"E_repulsion", 0.043528741758}, {"E_coulomb", -0.061320827255},
                 {"E_dispersion", -0.000981258152}, {"E_batm", -0.000004201327}
             }
            }
        };
    }

    bool test_energy_regression() {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "GFN-FF Energy Component Regression Tests" << std::endl;
        std::cout << std::string(80, '=') << std::endl;

        CurcumaLogger::set_verbosity(config.verbosity);

        for (const auto& ref : energy_refs) {
            std::cout << "\nTesting: " << ref.description << std::endl;
            std::cout << std::string(50, '-') << std::endl;

            try {
                // Load molecule
                Molecule mol(ref.molecule_file);

                // CRITICAL FIX: Set charge and spin like working test_energy_methods.cpp!
                mol.setCharge(0);   // Neutral molecule for GFN-FF
                mol.setSpin(0);     // Doublet = 0 (neutral, closed shell)

                // Debug output to verify molecule state
                if (config.verbosity >= 2) {
                    std::cout << "  DEBUG: Molecule loaded successfully!" << std::endl;
                    std::cout << "  DEBUG: Atoms = " << mol.AtomCount() << std::endl;
                    std::cout << "  DEBUG: Charge = " << mol.Charge() << std::endl;
                    std::cout << "  DEBUG: Spin = " << mol.getMolInfo().m_spin << std::endl;
                }

                Mol mol_info = mol.getMolInfo();

                // PHASE 1 FIX: Direct GFNFF instantiation for component access
                // (Claude Generated 2025-12-13)
                GFNFF gfnff;
                bool init_success = gfnff.InitialiseMolecule(mol_info);
                if (!init_success) {
                    std::cout << "  ✗ ERROR: GFNFF initialization failed" << std::endl;
                    total_tests++;
                    continue;
                }

                double energy = gfnff.Calculation(false);  // false = no gradient

                // NOW we can access all energy components!
                double E_bond = gfnff.BondEnergy();
                double E_angle = gfnff.AngleEnergy();
                double E_torsion = gfnff.DihedralEnergy();
                double E_repulsion = gfnff.RepulsionEnergy();
                double E_coulomb = gfnff.CoulombEnergy();
                double E_dispersion = gfnff.DispersionEnergy();

                // Store in map for later validation (Phase 2)
                std::map<std::string, double> components = {
                    {"bond", E_bond},
                    {"angle", E_angle},
                    {"torsion", E_torsion},
                    {"repulsion", E_repulsion},
                    {"coulomb", E_coulomb},
                    {"dispersion", E_dispersion}
                };

                if (config.verbosity >= 2) {
                    std::cout << "  Total energy: " << std::fixed << std::setprecision(8) << energy << " Hartree" << std::endl;
                    std::cout << "  Energy components:" << std::endl;
                    std::cout << "    Bond:       " << std::setw(12) << E_bond << " Eh" << std::endl;
                    std::cout << "    Angle:      " << std::setw(12) << E_angle << " Eh" << std::endl;
                    std::cout << "    Torsion:    " << std::setw(12) << E_torsion << " Eh" << std::endl;
                    std::cout << "    Repulsion:  " << std::setw(12) << E_repulsion << " Eh" << std::endl;
                    std::cout << "    Coulomb:    " << std::setw(12) << E_coulomb << " Eh" << std::endl;
                    std::cout << "    Dispersion: " << std::setw(12) << E_dispersion << " Eh" << std::endl;
                }

                // Validate total energy
                double energy_error = std::abs(energy - ref.total_energy);
                if (energy_error < config.energy_tolerance) {
                    std::cout << "  ✓ Total energy: error = " << std::scientific << energy_error << " Hartree" << std::endl;
                    passed_tests++;
                } else {
                    std::cout << "  ✗ Total energy: error = " << std::scientific << energy_error << " Hartree (expected " << std::fixed << std::setprecision(8) << ref.total_energy << ")" << std::endl;
                }
                total_tests++;

                // Validate individual energy components
                std::map<std::string, double> ref_components = {
                    {"bond", ref.bond_energy},
                    {"angle", ref.angle_energy},
                    {"torsion", ref.torsion_energy},
                    {"repulsion", ref.repulsion_energy},
                    {"coulomb", ref.electrostat_energy},
                    {"dispersion", ref.dispersion_energy}
                };

                for (const auto& [name, value] : components) {
                    double reference = ref_components[name];
                    double error = std::abs(value - reference);

                    if (error < config.component_tolerance) {
                        std::cout << "  ✓ " << name << " energy: error = "
                                  << std::scientific << error << " Eh" << std::endl;
                        passed_tests++;
                    } else {
                        std::cout << "  ✗ " << name << " energy: error = "
                                  << std::scientific << error
                                  << " Eh (got " << std::fixed << std::setprecision(8) << value
                                  << ", expected " << reference << ")" << std::endl;
                    }
                    total_tests++;
                }

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                total_tests++;
            }
        }

        return (passed_tests == total_tests);
    }

    bool test_parameter_verification() {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "GFN-FF Parameter Verification Tests" << std::endl;
        std::cout << std::string(80, '=') << std::endl;

        CurcumaLogger::set_verbosity(0);  // Silent for cleaner output

        for (const auto& ref : vbond_refs) {
            std::cout << "\nTesting: " << ref.bond_description << std::endl;
            std::cout << std::string(50, '-') << std::endl;

            try {
                // Load molecule and extract GFNFF parameters directly
                Molecule mol(ref.molecule_file);

                // CRITICAL FIX: Set charge and spin like working test_energy_methods.cpp!
                mol.setCharge(0);   // Neutral molecule for GFN-FF
                mol.setSpin(0);     // Doublet = 0 (neutral, closed shell)

                Mol mol_info = mol.getMolInfo();

                // Direct parameter generation via GFNFF class
                GFNFF gfnff;
                gfnff.InitialiseMolecule(mol_info);

                // Get vbond parameters via the correct method signature
                double shift, alpha, force_constant;
                bool success = gfnff.getVBondParameters(ref.bond_index, shift, alpha, force_constant);

                if (success) {
                    double shift_error = std::abs(shift - ref.shift);
                    double alpha_error = std::abs(alpha - ref.alpha);
                    double fc_error = std::abs(force_constant - ref.force_constant);

                    bool shift_ok = shift_error < config.param_tolerance;
                    bool alpha_ok = alpha_error < config.param_tolerance;
                    bool fc_ok = fc_error < config.param_tolerance;

                    if (shift_ok && alpha_ok && fc_ok) {
                        std::cout << "  ✓ vbond parameters: shift=" << std::fixed << std::setprecision(8) << shift
                                  << " alpha=" << alpha << " fc=" << force_constant << std::endl;
                        passed_tests++;
                    } else {
                        std::cout << "  ✗ vbond parameters:" << std::endl;
                        std::cout << "    Expected: shift=" << std::fixed << std::setprecision(8) << ref.shift
                                  << " alpha=" << ref.alpha << " fc=" << ref.force_constant << std::endl;
                        std::cout << "    Actual:   shift=" << shift << " alpha=" << alpha
                                  << " fc=" << force_constant << std::endl;
                    }

                    if (config.verbosity >= 2) {
                        std::cout << "    Errors: shift=" << std::scientific << shift_error
                                  << " alpha=" << alpha_error << " fc=" << fc_error << std::endl;
                    }
                } else {
                    std::cout << "  ✗ ERROR: Could not extract vbond parameters" << std::endl;
                }

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
            }

            total_tests++;
        }

        return true;
    }

    bool test_comprehensive_validation() {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "GFN-FF Comprehensive Validation Tests" << std::endl;
        std::cout << std::string(80, '=') << std::endl;

        for (const auto& ref : comprehensive_refs) {
            std::cout << "\nTesting: " << ref.description << std::endl;
            std::cout << std::string(50, '-') << std::endl;

            try {
                Molecule mol(ref.molecule_file);

                // CRITICAL FIX: Set charge and spin like working test_energy_methods.cpp!
                mol.setCharge(0);   // Neutral molecule for GFN-FF
                mol.setSpin(0);     // Doublet = 0 (neutral, closed shell)

                Mol mol_info = mol.getMolInfo();

                if (config.verbosity >= 2) {
                    std::cout << "  Molecule: " << mol.AtomCount() << " atoms" << std::endl;
                }

                nlohmann::json controller;
                controller["method"] = "cgfnff";
                controller["cgfnff"]["verbosity"] = 0;

                EnergyCalculator calculator("cgfnff", controller);
                calculator.setMolecule(mol_info);
                double total_energy = calculator.CalculateEnergy();

                std::cout << "\n  Total Energy Validation:" << std::endl;
                double energy_error = std::abs(total_energy - ref.energies.at("E_total"));
                if (energy_error < config.energy_tolerance) {
                    std::cout << "  ✓ Total energy: " << std::fixed << std::setprecision(8) << total_energy
                              << " (error=" << std::scientific << energy_error << ")" << std::endl;
                    passed_tests++;
                } else {
                    std::cout << "  ✗ Total energy: " << std::fixed << std::setprecision(8) << total_energy
                              << " (expected=" << ref.energies.at("E_total")
                              << ", error=" << std::scientific << energy_error << ")" << std::endl;
                }
                total_tests++;

                // TODO: Add CN, charges, and individual energy component validation
                // This would require extending GFNFF interface to expose these values

                std::cout << "\n  Note: Detailed CN, charge, and component validation requires GFNFF interface extension" << std::endl;

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                total_tests++;
            }
        }

        return true;
    }

    bool test_parameter_flag_combinations() {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "GFN-FF Parameter Flag Combination Tests" << std::endl;
        std::cout << std::string(80, '=') << std::endl;

        // Test 1: Dispersion disabled
        {
            std::cout << "\nTest 1: Dispersion disabled" << std::endl;
            std::cout << std::string(40, '-') << std::endl;

            try {
                Molecule mol("test_cases/molecules/larger/CH3OH.xyz");
                json params = {
                    {"method", "cgfnff"},
                    {"gradient", 0},
                    {"dispersion", false},
                    {"hbond", true},
                    {"repulsion", true},
                    {"coulomb", true}
                };

                EnergyCalculator calculator("cgfnff", params);
                calculator.setMolecule(mol.getMolInfo());
                double energy = calculator.CalculateEnergy();

                std::cout << "  Energy (no dispersion): " << std::scientific << energy << " Eh" << std::endl;

                // Compare with full calculation - energy should be less negative without dispersion
                bool test_passed = (energy > -1.0); // Basic sanity check
                std::cout << "  " << (test_passed ? "✓ PASSED" : "✗ FAILED") << std::endl;
                total_tests++;
                if (test_passed) passed_tests++;

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                total_tests++;
            }
        }

        // Test 2: Hydrogen bond disabled
        {
            std::cout << "\nTest 2: Hydrogen bond disabled" << std::endl;
            std::cout << std::string(40, '-') << std::endl;

            try {
                Molecule mol("test_cases/molecules/trimers/water.xyz");
                json params = {
                    {"method", "cgfnff"},
                    {"gradient", 0},
                    {"dispersion", true},
                    {"hbond", false},
                    {"repulsion", true},
                    {"coulomb", true}
                };

                EnergyCalculator calculator("cgfnff", params);
                calculator.setMolecule(mol.getMolInfo());
                double energy = calculator.CalculateEnergy();

                std::cout << "  Energy (no hbond): " << std::scientific << energy << " Eh" << std::endl;

                bool test_passed = (energy < -0.1); // Basic sanity check
                std::cout << "  " << (test_passed ? "✓ PASSED" : "✗ FAILED") << std::endl;
                total_tests++;
                if (test_passed) passed_tests++;

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                total_tests++;
            }
        }

        // Test 3: All non-bonded terms disabled (should only have bonded energy)
        {
            std::cout << "\nTest 3: All non-bonded terms disabled" << std::endl;
            std::cout << std::string(40, '-') << std::endl;

            try {
                Molecule mol("test_cases/molecules/larger/CH4.xyz");
                json params = {
                    {"method", "cgfnff"},
                    {"gradient", 0},
                    {"dispersion", false},
                    {"hbond", false},
                    {"repulsion", false},
                    {"coulomb", false}
                };

                EnergyCalculator calculator("cgfnff", params);
                calculator.setMolecule(mol.getMolInfo());
                double energy = calculator.CalculateEnergy();

                std::cout << "  Energy (bonded only): " << std::scientific << energy << " Eh" << std::endl;

                // With only bonded terms, energy should be much less negative
                bool test_passed = (energy > -2.0 && energy < 0.0);
                std::cout << "  " << (test_passed ? "✓ PASSED" : "✗ FAILED") << std::endl;
                total_tests++;
                if (test_passed) passed_tests++;

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                total_tests++;
            }
        }

        // Test 4: Edge case - atoms at cutoff distance
        {
            std::cout << "\nTest 4: Edge case - atoms at cutoff distance" << std::endl;
            std::cout << std::string(40, '-') << std::endl;

            // Create a simple H2 system with atoms exactly at cutoff distance
            try {
                Molecule mol;
                mol.setCharge(0);
                mol.setSpin(1);
                mol.addAtom({1, {0.0, 0.0, 0.0}});  // H at origin
                mol.addAtom({1, {15.0, 0.0, 0.0}});  // H at cutoff distance (15 Bohr)

                json params = {
                    {"method", "cgfnff"},
                    {"gradient", 0},
                    {"dispersion", true},
                    {"hbond", false},
                    {"repulsion", true},
                    {"coulomb", true}
                };

                EnergyCalculator calculator("cgfnff", params);
                calculator.setMolecule(mol.getMolInfo());
                double energy = calculator.CalculateEnergy();

                std::cout << "  Energy (at cutoff): " << std::scientific << energy << " Eh" << std::endl;

                // Should be very close to zero interaction energy at cutoff
                bool test_passed = (std::abs(energy) < 1e-6);
                std::cout << "  " << (test_passed ? "✓ PASSED" : "✗ FAILED") << std::endl;
                total_tests++;
                if (test_passed) passed_tests++;

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                total_tests++;
            }
        }

        // Test 5: Metal-specific correction handling (Fe atom test)
        {
            std::cout << "\nTest 5: Metal-specific correction handling" << std::endl;
            std::cout << std::string(40, '-') << std::endl;

            // Create a single Fe atom
            try {
                Molecule mol;
                mol.setCharge(0);
                mol.setSpin(3);  // Fe has 4 unpaired electrons in ground state
                mol.addAtom({26, {0.0, 0.0, 0.0}});  // Fe (Z=26)

                json params = {
                    {"method", "cgfnff"},
                    {"gradient", 0},
                    {"dispersion", false},
                    {"hbond", false},
                    {"repulsion", false},
                    {"coulomb", false}
                };

                EnergyCalculator calculator("cgfnff", params);
                calculator.setMolecule(mol.getMolInfo());
                double energy = calculator.CalculateEnergy();

                std::cout << "  Energy (single Fe atom): " << std::scientific << energy << " Eh" << std::endl;

                // Single atom should have minimal energy
                bool test_passed = ( std::abs(energy) < 1e-6);
                std::cout << "  " << (test_passed ? "✓ PASSED" : "✗ FAILED") << std::endl;
                total_tests++;
                if (test_passed) passed_tests++;

            } catch (const std::exception& e) {
                std::cout << "  ✗ ERROR: " << e.what() << std::endl;
                total_tests++;
            }
        }

        return true;
    }

    bool run_all_tests() {
        std::cout << "Unified GFN-FF Test Suite" << std::endl;
        std::cout << "========================" << std::endl;
        std::cout << "Configuration:" << std::endl;
        std::cout << "  Energy tolerance: " << std::scientific << config.energy_tolerance << " Hartree" << std::endl;
        std::cout << "  Parameter tolerance: " << std::scientific << config.param_tolerance << std::endl;
        std::cout << "  Verbosity level: " << config.verbosity << std::endl;

        // Run all test categories
        test_energy_regression();
        test_parameter_verification();
        test_comprehensive_validation();
        test_parameter_flag_combinations();

        // Summary
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << "UNIFIED GFN-FF TEST SUITE SUMMARY" << std::endl;
        std::cout << std::string(80, '=') << std::endl;
        std::cout << "Total tests: " << total_tests << std::endl;
        std::cout << "Passed: " << passed_tests << std::endl;
        std::cout << "Failed: " << (total_tests - passed_tests) << std::endl;
        std::cout << "Success rate: " << std::fixed << std::setprecision(1)
                  << (100.0 * passed_tests / total_tests) << "%" << std::endl;

        return (passed_tests == total_tests);
    }

    int get_total_tests() const { return total_tests; }
    int get_passed_tests() const { return passed_tests; }
};

int main(int argc, char** argv) {
    TestConfig config;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--verbose" || arg == "-v") {
            config.verbosity = 2;
        } else if (arg == "--quiet" || arg == "-q") {
            config.verbosity = 0;
        } else if (arg == "--tolerance" && i + 1 < argc) {
            config.energy_tolerance = std::stod(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Unified GFN-FF Test Suite\n";
            std::cout << "Usage: " << argv[0] << " [options]\n";
            std::cout << "Options:\n";
            std::cout << "  -v, --verbose     Enable verbose output\n";
            std::cout << "  -q, --quiet       Silent mode (minimal output)\n";
            std::cout << "  --tolerance VAL   Set energy tolerance (default: 1e-6)\n";
            std::cout << "  -h, --help        Show this help\n";
            return 0;
        }
    }

    UnifiedGFNFFTest test_suite(config);
    bool success = test_suite.run_all_tests();

    return success ? 0 : 1;
}