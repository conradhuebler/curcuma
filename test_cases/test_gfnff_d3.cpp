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

std::vector<GFNFFDispersionReference> createReferenceData() {
    /**
     * Reference energies generated with D3ParameterGenerator using GFN-FF parameters.
     *
     * These values serve as the "ground truth" for what GFN-FF's D3 energy should be.
     * The D3ParameterGenerator is already validated against s-dftd3 (10/11 molecules <1% error).
     */

    std::vector<GFNFFDispersionReference> refs;

    // Reference energies generated with D3ParameterGenerator using GFN-FF parameters
    // Generated: December 19, 2025

    refs.push_back({
        "../test_cases/molecules/dimers/HH.xyz",
        "H2 dimer",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -3.117228579267e-05,  // expected_d3_energy
        REFERENCE_TOLERANCE, CONSISTENCY_TOLERANCE,
        -999.0  // xtb_gfnff_dispersion - not available
    });

    refs.push_back({
        "../test_cases/molecules/dimers/HCl.xyz",
        "HCl dimer",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -9.126452453995e-05,  // expected_d3_energy
        REFERENCE_TOLERANCE, CONSISTENCY_TOLERANCE,
        -999.0  // xtb_gfnff_dispersion - not available
    });

    refs.push_back({
        "../test_cases/molecules/dimers/OH.xyz",
        "OH radical",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -4.859150805298e-05,  // expected_d3_energy
        REFERENCE_TOLERANCE, CONSISTENCY_TOLERANCE,
        -999.0  // xtb_gfnff_dispersion - not available
    });

    refs.push_back({
        "../test_cases/molecules/larger/CH4.xyz",
        "Methane",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -3.837717923381e-04,  // expected_d3_energy
        REFERENCE_TOLERANCE, CONSISTENCY_TOLERANCE,
        -999.0  // xtb_gfnff_dispersion - not available
    });

    refs.push_back({
        "../test_cases/molecules/larger/CH3OH.xyz",
        "Methanol",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -6.142434770501e-04,  // expected_d3_energy
        REFERENCE_TOLERANCE, CONSISTENCY_TOLERANCE,
        -999.0  // xtb_gfnff_dispersion - not available
    });

    refs.push_back({
        "../test_cases/molecules/larger/CH3OCH3.xyz",
        "Dimethyl ether",
        GFNFF_D3_S6, GFNFF_D3_S8, GFNFF_D3_A1, GFNFF_D3_A2,
        -1.421175523031e-03,  // expected_d3_energy
        REFERENCE_TOLERANCE, CONSISTENCY_TOLERANCE,
        -999.0  // xtb_gfnff_dispersion - not available
    });

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
        // Load molecule
        auto mol = Files::LoadFile(ref.molecule_file);
        if (mol.AtomCount() == 0) {
            result.error_message = "Failed to load molecule: " + ref.molecule_file;
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
