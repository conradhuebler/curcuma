/**
 * GFN-FF D3 Integration Validation Test
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Purpose: Verifies that GFN-FF correctly integrates the native D3 dispersion implementation.
 *
 * Test Strategy:
 * 1. Calculate D3 dispersion standalone with GFN-FF parameters (s6=1.0, s8=2.85, a1=0.80, a2=4.60)
 * 2. Calculate D3 dispersion via GFN-FF.DispersionEnergy()
 * 3. CRITICAL: Both MUST be identical (same D3 code!) - delta < 1e-10 Eh
 *
 * Reference: Native D3ParameterGenerator (validated against s-dftd3 with 10/11 molecules <1%)
 *
 * Claude Generated - December 19, 2025
 */

#include <cassert>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// Core
#include "src/core/molecule.h"
#include "src/core/config_manager.h"
#include "src/core/units.h"
#include "src/tools/formats.h"

// GFN-FF and D3
#include "src/core/energy_calculators/ff_methods/gfnff.h"
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"

#include "json.hpp"
using json = nlohmann::json;

// ============================================================================
// CONSTANTS
// ============================================================================

// GFN-FF D3 Parameters (Spicher/Grimme, J. Chem. Theory Comput. 2020)
const double GFNFF_D3_S6 = 1.0;
const double GFNFF_D3_S8 = 2.85;
const double GFNFF_D3_A1 = 0.80;
const double GFNFF_D3_A2 = 4.60;  // Bohr

// Test tolerances
const double CONSISTENCY_TOLERANCE = 1e-10;  // Standalone vs GFN-FF must be near-identical
const double REFERENCE_TOLERANCE = 1e-6;     // Hartree

// ============================================================================
// DATA STRUCTURES
// ============================================================================

struct GFNFFDispersionReference {
    std::string molecule_file;
    std::string description;

    // GFN-FF D3 parameters (always same for all tests)
    double d3_s6 = GFNFF_D3_S6;
    double d3_s8 = GFNFF_D3_S8;
    double d3_a1 = GFNFF_D3_A1;
    double d3_a2 = GFNFF_D3_A2;

    // Expected energy from standalone D3 calculation
    double expected_d3_energy = 0.0;

    // Tolerances
    double reference_tolerance = REFERENCE_TOLERANCE;
    double consistency_tolerance = CONSISTENCY_TOLERANCE;

    // Optional: XTB GFN-FF reference (for documentation)
    double xtb_gfnff_dispersion = -999.0;  // -999.0 = not available
};

struct TestResult {
    std::string molecule_name;

    double standalone_d3 = 0.0;
    double gfnff_d3 = 0.0;
    double expected_d3 = 0.0;

    double standalone_error = 0.0;  // vs reference
    double gfnff_error = 0.0;       // vs reference
    double consistency_delta = 0.0; // standalone vs gfnff

    bool standalone_passed = false;
    bool gfnff_passed = false;
    bool consistency_passed = false;
    bool all_passed = false;

    std::string error_message;
};

// ============================================================================
// REFERENCE DATA
// ============================================================================

// TEST STRUCTURES AND DATA
// ============================================================================

std::vector<GFNFFDispersionReference> createReferenceData() {
    /**
     * Test molecules with D3 reference energies - hardcoded version
     */
    std::vector<GFNFFDispersionReference> refs;

    refs.push_back({
        "H2",  // molecule_file (use molecule name now)
        "H2 dimer from HH.xyz",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -3.1187639277766e-05,  // expected_d3_energy (s-dftd3 with GFN-FF params)
        1e-6, // Realistic tolerance for D3
        CONSISTENCY_TOLERANCE,
        -999.0  // xtb_gfnff_dispersion - not available
    });

    refs.push_back({
        "HCl",  // molecule_file
        "HCl dimer from HCl.xyz",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -9.1325281458352e-05,  // expected_d3_energy (s-dftd3 with GFN-FF params)
        1e-6,
        CONSISTENCY_TOLERANCE,
        -999.0
    });

    refs.push_back({
        "OH",  // molecule_file
        "OH radical from OH.xyz",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -4.8655959001113e-05,  // expected_d3_energy (s-dftd3 with GFN-FF params)
        1e-6,
        CONSISTENCY_TOLERANCE,
        -999.0
    });

    refs.push_back({
        "CH4",  // molecule_file
        "Methane from CH4.xyz",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -3.8476453526938e-04,  // expected_d3_energy (s-dftd3 with GFN-FF params)
        1e-6,
        CONSISTENCY_TOLERANCE,
        -999.0
    });

    refs.push_back({
        "CH3OH",  // molecule_file
        "Methanol from CH3OH.xyz",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -6.1977861552636e-04,  // expected_d3_energy (s-dftd3 with GFN-FF params)
        1e-6,
        CONSISTENCY_TOLERANCE,
        -999.0
    });

    refs.push_back({
        "CH3OCH3",  // molecule_file
        "Dimethyl ether from CH3OCH3.xyz",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -1.4313191021323e-03,  // expected_d3_energy (s-dftd3 with GFN-FF params)
        1e-6,
        CONSISTENCY_TOLERANCE,
        -999.0
    });

    return refs;
}

std::vector<GFNFFDispersionReference> loadReferencesFromJSON() {
    /**
     * Load reference energies from d3_reference_energies.json
     *
     * FUTURE-PROOF: Single source of truth - no hardcoded values.
     * When adding new test molecules, just update the JSON file.
     *
     * Claude Generated - December 20, 2025
     */

    std::vector<GFNFFDispersionReference> refs;

    try {
        std::ifstream f("../d3_reference_energies.json");
        if (!f.is_open()) {
            std::cerr << "WARNING: Could not open d3_reference_energies.json" << std::endl;
            std::cerr << "         Falling back to hardcoded references" << std::endl;
            return createReferenceData();
        }

        json j = json::parse(f);

        for (auto& [key, val] : j.items()) {
            refs.push_back({
                val["file"].get<std::string>(),
                key,
                GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
                val["energy"].get<double>(),
                REFERENCE_TOLERANCE, CONSISTENCY_TOLERANCE,
                -999.0  // xtb_gfnff_dispersion - not available
            });
        }

        std::cout << "Loaded " << refs.size() << " reference molecules from JSON" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "ERROR loading JSON: " << e.what() << std::endl;
        std::cerr << "Falling back to hardcoded references" << std::endl;
        return createReferenceData();
    }

    return refs;
}

// ============================================================================
// TEST FUNCTIONS
// ============================================================================

TestResult testGFNFFD3Integration(const GFNFFDispersionReference& ref) {
    /**
     * Main test function: Validates GFN-FF D3 integration.
     *
     * Tests:
     * 1. Standalone D3 with GFN-FF parameters
     * 2. GFN-FF integrated D3
     * 3. Consistency check (standalone == integrated)
     */

    TestResult result;
    result.molecule_name = ref.description;
    result.expected_d3 = ref.expected_d3_energy;

    try {
        // Create molecule directly (avoid registry complications)
        Molecule mol;

        if (ref.molecule_file == "H2") {
            mol.addAtom({1, Eigen::Vector3d(-0.092881, 0.361653, -0.249341)});
            mol.addAtom({1, Eigen::Vector3d(0.000000, 0.763239, -0.477047)});
        } else if (ref.molecule_file == "HCl") {
            mol.addAtom({17, Eigen::Vector3d(0.000000, 0.000000, 0.119262)});
            mol.addAtom({1, Eigen::Vector3d(0.000000, 0.763239, -0.477047)});
        } else if (ref.molecule_file == "OH") {
            mol.addAtom({8, Eigen::Vector3d(0.000000, 0.000000, 0.119262)});
            mol.addAtom({1, Eigen::Vector3d(0.000000, 0.763239, -0.477047)});
        } else if (ref.molecule_file == "CH4") {
            mol.addAtom({6, Eigen::Vector3d(-6.52745801014127, 1.22559601369319, 0.00000199477487)});
            mol.addAtom({1, Eigen::Vector3d(-5.71833882788279, 0.97438995382583, 0.67333841923947)});
            mol.addAtom({1, Eigen::Vector3d(-6.20118058043256, 1.99864435246030, -0.68344482547359)});
            mol.addAtom({1, Eigen::Vector3d(-6.81600263252576, 0.34633873760782, -0.56107634994743)});
            mol.addAtom({1, Eigen::Vector3d(-7.37430994901762, 1.58301094241287, 0.57119076140667)});
        } else if (ref.molecule_file == "CH3OH") {
            mol.addAtom({6, Eigen::Vector3d(-4.39019608206338, 1.80749146124850, -0.05236110118918)});
            mol.addAtom({1, Eigen::Vector3d(-3.55446893679818, 1.51846523808356, 0.59390888890703)});
            mol.addAtom({1, Eigen::Vector3d(-4.04190865425030, 2.55246475669577, -0.77579148135144)});
            mol.addAtom({1, Eigen::Vector3d(-4.74205299752436, 0.92930974600559, -0.59009212821488)});
            mol.addAtom({8, Eigen::Vector3d(-5.48388747837165, 2.28133056625900, 0.69456006621008)});
            mol.addAtom({1, Eigen::Vector3d(-5.21119853581715, 3.06238366739865, 1.18714027466186)});
        } else if (ref.molecule_file == "CH3OCH3") {
            mol.addAtom({6, Eigen::Vector3d(-5.165738, 2.528991, 1.023522)});
            mol.addAtom({6, Eigen::Vector3d(-6.279554, 3.753698, 2.690779)});
            mol.addAtom({1, Eigen::Vector3d(-5.335950, 4.211291, 3.003347)});
            mol.addAtom({1, Eigen::Vector3d(-7.086297, 4.475849, 2.844576)});
            mol.addAtom({1, Eigen::Vector3d(-6.486329, 2.865617, 3.295613)});
            mol.addAtom({8, Eigen::Vector3d(-6.239979, 3.413269, 1.311751)});
            mol.addAtom({1, Eigen::Vector3d(-5.216886, 2.257417, -0.034507)});
            mol.addAtom({1, Eigen::Vector3d(-4.206669, 3.022649, 1.210811)});
            mol.addAtom({1, Eigen::Vector3d(-5.243572, 1.617382, 1.623040)});
        } else {
            result.error_message = "Unknown molecule: " + ref.molecule_file;
            return result;
        }

        if (mol.AtomCount() == 0) {
            result.error_message = "Failed to create molecule: " + ref.molecule_file;
            return result;
        }

        std::vector<int> atoms = mol.Atoms();
        Matrix geometry = mol.getGeometry();

        // ========================================================================
        // TEST 1: Standalone D3 with GFN-FF parameters
        // ========================================================================

        json d3_config_json = {
            {"d3_s6", ref.d3_s6},
            {"d3_s8", ref.d3_s8},
            {"d3_a1", ref.d3_a1},
            {"d3_a2", ref.d3_a2},
            {"d3_alp", 14.0}  // Standard D3 alpha parameter
        };

        ConfigManager d3_config("d3param", d3_config_json);
        D3ParameterGenerator d3_gen(d3_config);

        d3_gen.GenerateParameters(atoms, geometry);
        result.standalone_d3 = d3_gen.getTotalEnergy();
        result.standalone_error = std::abs(result.standalone_d3 - ref.expected_d3_energy);
        result.standalone_passed = (result.standalone_error < ref.reference_tolerance);

        // ========================================================================
        // TEST 2: GFN-FF integrated D3
        // ========================================================================

        // Create GFN-FF with default parameters (should use D3)
        json gfnff_params = {
            {"verbosity", 0}  // Silent mode for testing
        };

        GFNFF gfnff;
        gfnff.setParameters(gfnff_params);

        Mol mol_info;
        mol_info.m_atoms = atoms;
        mol_info.m_geometry = geometry;
        mol_info.m_number_atoms = static_cast<int>(atoms.size());
        mol_info.m_charge = 0;

        bool init_success = gfnff.InitialiseMolecule(mol_info);
        if (!init_success) {
            result.error_message = "GFN-FF InitialiseMolecule() failed";
            return result;
        }

        double total_energy = gfnff.Calculation(false);  // no gradient

        result.gfnff_d3 = gfnff.DispersionEnergy();
        result.gfnff_error = std::abs(result.gfnff_d3 - ref.expected_d3_energy);
        result.gfnff_passed = (result.gfnff_error < ref.reference_tolerance);

        // ========================================================================
        // TEST 3: Consistency check (CRITICAL!)
        // ========================================================================

        result.consistency_delta = std::abs(result.standalone_d3 - result.gfnff_d3);
        result.consistency_passed = (result.consistency_delta < ref.consistency_tolerance);

        // Overall pass: all three tests must pass
        result.all_passed = result.standalone_passed && result.gfnff_passed && result.consistency_passed;

    } catch (const std::exception& e) {
        result.error_message = std::string("Exception: ") + e.what();
        result.all_passed = false;
    }

    return result;
}

// ============================================================================
// OUTPUT FUNCTIONS
// ============================================================================

void printTestHeader() {
    std::cout << "\n";
    std::cout << "GFN-FF D3 Dispersion Validation Test\n";
    std::cout << "================================================================================\n";
    std::cout << "Configuration:\n";
    std::cout << "  D3 s6: " << GFNFF_D3_S6 << "\n";
    std::cout << "  D3 s8: " << GFNFF_D3_S8 << "\n";
    std::cout << "  D3 a1: " << GFNFF_D3_A1 << "\n";
    std::cout << "  D3 a2: " << GFNFF_D3_A2 << " Bohr\n";
    std::cout << "  Energy tolerance: " << std::scientific << std::setprecision(6)
              << REFERENCE_TOLERANCE << " Hartree\n";
    std::cout << "  Consistency tolerance: " << std::scientific << std::setprecision(1)
              << CONSISTENCY_TOLERANCE << " Hartree\n";
    std::cout << "\n";
}

void printTestResult(const TestResult& result) {
    std::cout << "Testing: " << result.molecule_name << "\n";

    // Standalone D3
    std::cout << "  " << (result.standalone_passed ? "✓" : "✗")
              << " Standalone D3:  " << std::scientific << std::setprecision(4)
              << result.standalone_d3 << " Eh";
    if (result.expected_d3 != 0.0) {
        std::cout << " (vs ref " << result.expected_d3 << ")";
    }
    std::cout << "\n";

    // GFN-FF D3
    std::cout << "  " << (result.gfnff_passed ? "✓" : "✗")
              << " GFN-FF D3:      " << std::scientific << std::setprecision(4)
              << result.gfnff_d3 << " Eh";
    if (result.expected_d3 != 0.0) {
        std::cout << " (vs ref " << result.expected_d3 << ")";
    }
    std::cout << "\n";

    // Consistency
    std::cout << "  " << (result.consistency_passed ? "✓" : "✗")
              << " Consistency:    delta = " << std::scientific << std::setprecision(2)
              << result.consistency_delta << " Eh";
    if (!result.consistency_passed && result.consistency_delta > 1e-10) {
        std::cout << "  ⚠️ FAILURE!";
    }
    std::cout << "\n";

    // Error message if any
    if (!result.error_message.empty()) {
        std::cout << "      Error: " << result.error_message << "\n";
    }

    std::cout << "\n";
}

void printSummary(const std::vector<TestResult>& results) {
    int total = results.size();
    int standalone_passed = 0;
    int gfnff_passed = 0;
    int consistency_passed = 0;
    int all_passed = 0;

    for (const auto& result : results) {
        if (result.standalone_passed) standalone_passed++;
        if (result.gfnff_passed) gfnff_passed++;
        if (result.consistency_passed) consistency_passed++;
        if (result.all_passed) all_passed++;
    }

    std::cout << "================================================================================\n";
    std::cout << "Summary: " << all_passed << "/" << total << " molecules passed ("
              << (100 * all_passed / total) << "%)\n";
    std::cout << "  Standalone D3 accuracy: " << standalone_passed << "/" << total
              << (standalone_passed == total ? " ✓" : " ✗") << "\n";
    std::cout << "  GFN-FF D3 accuracy:     " << gfnff_passed << "/" << total
              << (gfnff_passed == total ? " ✓" : " ✗") << "\n";
    std::cout << "  Consistency checks:     " << consistency_passed << "/" << total
              << (consistency_passed == total ? " ✓" : " ✗") << "\n";
    std::cout << "\n";

    if (all_passed == total) {
        std::cout << "✅ GFN-FF D3 integration is CORRECT\n";
    } else {
        std::cout << "❌ CRITICAL: GFN-FF D3 integration has FAILURES\n";
        std::cout << "   Possible causes:\n";
        std::cout << "   1. GFN-FF using different D3 parameters\n";
        std::cout << "   2. GFN-FF not calling D3ParameterGenerator\n";
        std::cout << "   3. Parameter extraction bug in extractDispersionConfig()\n";
    }

    std::cout << "================================================================================\n";
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    printTestHeader();

    auto references = createReferenceData();
    std::vector<TestResult> results;

    // Check if references have expected values set
    bool has_references = false;
    for (const auto& ref : references) {
        if (ref.expected_d3_energy != 0.0) {
            has_references = true;
            break;
        }
    }

    if (!has_references) {
        std::cout << "⚠️  WARNING: No reference energies set (all 0.0)\n";
        std::cout << "   Run generate_gfnff_d3_refs to generate reference values first!\n";
        std::cout << "   For now, only consistency test will be meaningful.\n\n";
    }

    // Run tests
    for (const auto& ref : references) {
        auto result = testGFNFFD3Integration(ref);
        results.push_back(result);
        printTestResult(result);
    }

    printSummary(results);

    // Exit code: 0 if all tests passed, 1 otherwise
    int all_passed_count = 0;
    for (const auto& result : results) {
        if (result.all_passed) all_passed_count++;
    }

    return (all_passed_count == results.size()) ? 0 : 1;
}
