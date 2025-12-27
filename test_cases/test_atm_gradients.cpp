/**
 * ATM (Axilrod-Teller-Muto) 3-Body Dispersion Gradient Test
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Validates analytical ATM gradients against numerical finite differences
 *
 * Claude Generated - December 2025
 */

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "src/core/molecule.h"
#include "src/core/config_manager.h"
#include "src/core/global.h"
#include "src/core/curcuma_logger.h"
#include "src/tools/formats.h"
#include "core/test_molecule_registry.h"

#include "src/core/energy_calculators/ff_methods/forcefield.h"
#include "src/core/energy_calculators/ff_methods/forcefieldgenerator.h"
#include "src/core/energy_calculators/ff_methods/d3param_generator.h"

using namespace curcuma;
using namespace TestMolecules;

// ============================================================================
// Test Configuration
// ============================================================================

const double GRADIENT_TOLERANCE = 1e-6;  // Hartree/Bohr - analytical vs numerical
const double FINITE_DIFF_STEP = 1e-6;   // Bohr - finite difference step size (reduced from 1e-5 for better accuracy)

// ============================================================================
// Numerical Gradient Calculation
// ============================================================================

Matrix calculateNumericalGradient(ForceField& ff, const Matrix& geometry) {
    const int natoms = geometry.rows();
    Matrix numerical_gradient = Matrix::Zero(natoms, 3);

    Matrix geom_plus = geometry;
    Matrix geom_minus = geometry;

    // Calculate numerical gradient via central finite differences
    for (int i = 0; i < natoms; ++i) {
        for (int coord = 0; coord < 3; ++coord) {
            // E(x + dx)
            geom_plus(i, coord) = geometry(i, coord) + FINITE_DIFF_STEP;
            ff.UpdateGeometry(geom_plus);
            double E_plus = ff.Calculate(false);  // Energy only

            // E(x - dx)
            geom_minus(i, coord) = geometry(i, coord) - FINITE_DIFF_STEP;
            ff.UpdateGeometry(geom_minus);
            double E_minus = ff.Calculate(false);

            // Central difference: dE/dx = (E+ - E-) / (2*dx)
            numerical_gradient(i, coord) = (E_plus - E_minus) / (2.0 * FINITE_DIFF_STEP);

            // Reset geometry
            geom_plus(i, coord) = geometry(i, coord);
            geom_minus(i, coord) = geometry(i, coord);
        }
    }

    return numerical_gradient;
}

// ============================================================================
// Test Functions
// ============================================================================

bool testATMGradients_HCl() {
    std::cout << "========================================" << std::endl;
    std::cout << "Test: HCl ATM Gradient Validation" << std::endl;
    std::cout << "========================================" << std::endl;

    // Get molecule from registry
    Molecule mol = TestMoleculeRegistry::createMolecule("HCl", false);  // Angstrom
    std::vector<int> atoms = mol.Atoms();
    Matrix geometry_bohr = mol.getGeometry() * CurcumaUnit::Length::ANGSTROM_TO_BOHR;

    std::cout << "  Atoms: " << atoms.size() << std::endl;

    // Generate D3-ATM parameters with s9 enabled
    std::cout << "Generating D3-ATM parameters..." << std::endl;

    json d3_params_json = {
        {"d3_s9", 1.0},      // Enable ATM 3-body correction
        {"d3_a1", 0.4145},   // PBE0/BJ damping parameters
        {"d3_a2", 4.8593},
        {"d3_s8", 1.2177},
        {"d3_s6", 1.0}
    };
    ConfigManager d3_config("d3param", d3_params_json);
    D3ParameterGenerator d3_gen(d3_config);
    d3_gen.GenerateParameters(atoms, geometry_bohr);

    json d3_params = d3_gen.getParameters();

    // Check ATM triples were generated
    int natoms = atoms.size();
    if (!d3_params.contains("atm_triples") || d3_params["atm_triples"].size() == 0) {
        std::cerr << "WARNING: No ATM triples generated (need at least 3 atoms for ATM)" << std::endl;
        std::cerr << "         HCl has only " << natoms << " atoms - skipping gradient test" << std::endl;
        return true;  // Not a failure, just not applicable
    }

    int num_triples = d3_params["atm_triples"].size();
    std::cout << "  ATM triples generated: " << num_triples << std::endl;

    // Create ForceField with ATM parameters only (no other terms to isolate ATM)
    json ff_params = json{};
    ff_params["natoms"] = natoms;
    ff_params["method"] = "cgfnff";  // GFN-FF method (string, not int!)
    ff_params["atm_triples"] = d3_params["atm_triples"];

    ForceField ff(ff_params);
    ff.setParameterCaching(false);  // Disable caching for isolated ATM testing
    ff.setAtomTypes(atoms);
    ff.UpdateGeometry(geometry_bohr);
    ff.setParameter(ff_params);

    std::cout << "\nCalculating analytical ATM gradient..." << std::endl;
    double E_analytical = ff.Calculate(true);  // With gradient
    Matrix analytical_gradient = ff.Gradient();

    std::cout << "  ATM Energy: " << std::scientific << std::setprecision(10)
              << E_analytical << " Eh" << std::endl;

    std::cout << "\nCalculating numerical gradient (finite differences)..." << std::endl;
    Matrix numerical_gradient = calculateNumericalGradient(ff, geometry_bohr);

    // Compare gradients
    std::cout << "\nGradient Comparison:" << std::endl;
    std::cout << "  Atom | Coord | Analytical   | Numerical    | Difference   | Status" << std::endl;
    std::cout << "  -----|-------|--------------|--------------|--------------|--------" << std::endl;

    bool all_passed = true;
    double max_error = 0.0;

    for (int i = 0; i < natoms; ++i) {
        for (int coord = 0; coord < 3; ++coord) {
            double analytical = analytical_gradient(i, coord);
            double numerical = numerical_gradient(i, coord);
            double diff = std::abs(analytical - numerical);

            max_error = std::max(max_error, diff);

            bool passed = diff < GRADIENT_TOLERANCE;
            all_passed = all_passed && passed;

            const char* coord_name[] = {"x", "y", "z"};
            std::cout << "  " << std::setw(4) << i
                      << " | " << coord_name[coord]
                      << "     | " << std::scientific << std::setprecision(6) << analytical
                      << " | " << numerical
                      << " | " << diff
                      << " | " << (passed ? "✓ PASS" : "✗ FAIL") << std::endl;
        }
    }

    std::cout << "\n  Maximum error: " << std::scientific << max_error << " Eh/Bohr" << std::endl;
    std::cout << "  Tolerance: " << GRADIENT_TOLERANCE << " Eh/Bohr" << std::endl;

    // Translational invariance check (sum of gradients should be zero)
    Vector gradient_sum = analytical_gradient.colwise().sum();
    double translation_error = gradient_sum.norm();

    std::cout << "\nTranslational Invariance:" << std::endl;
    std::cout << "  Sum of gradients: (" << std::scientific << std::setprecision(6)
              << gradient_sum(0) << ", " << gradient_sum(1) << ", " << gradient_sum(2) << ")" << std::endl;
    std::cout << "  Norm: " << translation_error << " Eh/Bohr" << std::endl;
    std::cout << "  " << (translation_error < 1e-10 ? "✓ PASS" : "✗ FAIL") << std::endl;

    std::cout << "\n========================================" << std::endl;
    if (all_passed && translation_error < 1e-10) {
        std::cout << "✓ HCl ATM Gradient Test PASSED" << std::endl;
    } else {
        std::cout << "✗ HCl ATM Gradient Test FAILED" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    return all_passed && (translation_error < 1e-10);
}

bool testATMGradients_Methane() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test: CH4 ATM Gradient Validation" << std::endl;
    std::cout << "========================================" << std::endl;

    // Get molecule from registry
    Molecule mol = TestMoleculeRegistry::createMolecule("CH4", false);  // Angstrom
    std::vector<int> atoms = mol.Atoms();
    Matrix geometry_bohr = mol.getGeometry() * CurcumaUnit::Length::ANGSTROM_TO_BOHR;

    int natoms = atoms.size();
    std::cout << "  Atoms: " << natoms << std::endl;

    // Generate D3-ATM parameters with s9 enabled
    std::cout << "Generating D3-ATM parameters..." << std::endl;

    json d3_params_json = {
        {"d3_s9", 1.0},      // Enable ATM 3-body correction
        {"d3_a1", 0.4145},   // PBE0/BJ damping parameters
        {"d3_a2", 4.8593},
        {"d3_s8", 1.2177},
        {"d3_s6", 1.0}
    };
    ConfigManager d3_config("d3param", d3_params_json);
    D3ParameterGenerator d3_gen(d3_config);
    d3_gen.GenerateParameters(atoms, geometry_bohr);

    json d3_params = d3_gen.getParameters();

    int num_triples = d3_params["atm_triples"].size();
    int expected_triples = (natoms * (natoms - 1) * (natoms - 2)) / 6;  // C(n,3)
    std::cout << "  ATM triples generated: " << num_triples << std::endl;
    std::cout << "  Expected: " << expected_triples << " (for " << natoms << " atoms: C(" << natoms << ",3))" << std::endl;

    // Create ForceField with ATM parameters
    json ff_params = json{};
    ff_params["natoms"] = natoms;
    ff_params["method"] = "gfnff";  // GFN-FF method (string, not int!)
    ff_params["atm_triples"] = d3_params["atm_triples"];

    ForceField ff(ff_params);
    ff.setParameterCaching(false);  // Disable caching for isolated ATM testing
    ff.setAtomTypes(atoms);
    ff.UpdateGeometry(geometry_bohr);
    ff.setParameter(ff_params);

    std::cout << "\nCalculating analytical ATM gradient..." << std::endl;
    double E_analytical = ff.Calculate(true);
    Matrix analytical_gradient = ff.Gradient();

    std::cout << "  ATM Energy: " << std::scientific << std::setprecision(10)
              << E_analytical << " Eh" << std::endl;

    std::cout << "\nCalculating numerical gradient..." << std::endl;
    Matrix numerical_gradient = calculateNumericalGradient(ff, geometry_bohr);

    // Calculate maximum error
    double max_error = 0.0;
    for (int i = 0; i < natoms; ++i) {
        for (int coord = 0; coord < 3; ++coord) {
            double diff = std::abs(analytical_gradient(i, coord) - numerical_gradient(i, coord));
            max_error = std::max(max_error, diff);
        }
    }

    std::cout << "\nGradient Accuracy:" << std::endl;
    std::cout << "  Maximum error: " << std::scientific << max_error << " Eh/Bohr" << std::endl;
    std::cout << "  Tolerance: " << GRADIENT_TOLERANCE << " Eh/Bohr" << std::endl;

    bool passed = max_error < GRADIENT_TOLERANCE;

    // Translational invariance
    Vector gradient_sum = analytical_gradient.colwise().sum();
    double translation_error = gradient_sum.norm();

    std::cout << "\nTranslational Invariance:" << std::endl;
    std::cout << "  Norm of gradient sum: " << translation_error << " Eh/Bohr" << std::endl;
    std::cout << "  " << (translation_error < 1e-10 ? "✓ PASS" : "✗ FAIL") << std::endl;

    std::cout << "\n========================================" << std::endl;
    if (passed && translation_error < 1e-10) {
        std::cout << "✓ CH4 ATM Gradient Test PASSED" << std::endl;
    } else {
        std::cout << "✗ CH4 ATM Gradient Test FAILED" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    return passed && (translation_error < 1e-10);
}

bool testATMGradients_Monosaccharide() {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test: Monosaccharide ATM Gradient Validation" << std::endl;
    std::cout << "========================================" << std::endl;

    // Get molecule from registry
    Molecule mol = TestMoleculeRegistry::createMolecule("monosaccharide", false);  // Angstrom
    std::vector<int> atoms = mol.Atoms();
    Matrix geometry_bohr = mol.getGeometry() * CurcumaUnit::Length::ANGSTROM_TO_BOHR;

    int natoms = atoms.size();
    std::cout << "  Atoms: " << natoms << std::endl;

    // Generate D3-ATM parameters with s9 enabled
    std::cout << "Generating D3-ATM parameters..." << std::endl;

    json d3_params_json = {
        {"d3_s9", 1.0},      // Enable ATM 3-body correction
        {"d3_a1", 0.4145},   // PBE0/BJ damping parameters
        {"d3_a2", 4.8593},
        {"d3_s8", 1.2177},
        {"d3_s6", 1.0}
    };
    ConfigManager d3_config("d3param", d3_params_json);
    D3ParameterGenerator d3_gen(d3_config);
    d3_gen.GenerateParameters(atoms, geometry_bohr);

    json d3_params = d3_gen.getParameters();

    int num_triples = d3_params["atm_triples"].size();
    int expected_triples = (natoms * (natoms - 1) * (natoms - 2)) / 6;  // C(n,3)
    std::cout << "  ATM triples generated: " << num_triples << std::endl;
    std::cout << "  Expected: " << expected_triples << " (for " << natoms << " atoms: C(" << natoms << ",3))" << std::endl;

    // Verify triple count
    if (num_triples != expected_triples) {
        std::cerr << "ERROR: Triple count mismatch!" << std::endl;
        return false;
    }

    // Create ForceField with ATM parameters
    json ff_params = json{};
    ff_params["natoms"] = natoms;
    ff_params["method"] = "gfnff";  // GFN-FF method (string, not int!)
    ff_params["atm_triples"] = d3_params["atm_triples"];

    ForceField ff(ff_params);
    ff.setParameterCaching(false);  // Disable caching for isolated ATM testing
    ff.setAtomTypes(atoms);
    ff.UpdateGeometry(geometry_bohr);
    ff.setParameter(ff_params);

    std::cout << "\nCalculating analytical ATM gradient..." << std::endl;
    double E_analytical = ff.Calculate(true);
    Matrix analytical_gradient = ff.Gradient();

    std::cout << "  ATM Energy: " << std::scientific << std::setprecision(10)
              << E_analytical << " Eh" << std::endl;

    std::cout << "\nCalculating numerical gradient (finite differences)..." << std::endl;
    std::cout << "  This will take ~2-5 minutes (163 Calculate() calls, 2925 triples each)..." << std::endl;
    Matrix numerical_gradient = calculateNumericalGradient(ff, geometry_bohr);

    // Compare gradients
    std::cout << "\nGradient Comparison:" << std::endl;
    double max_error = 0.0;
    bool all_passed = true;

    for (int i = 0; i < natoms; ++i) {
        for (int coord = 0; coord < 3; ++coord) {
            double analytical = analytical_gradient(i, coord);
            double numerical = numerical_gradient(i, coord);
            double diff = std::abs(analytical - numerical);

            max_error = std::max(max_error, diff);

            if (diff >= GRADIENT_TOLERANCE) {
                all_passed = false;
            }
        }
    }

    std::cout << "  Maximum error: " << std::scientific << max_error << " Eh/Bohr" << std::endl;
    std::cout << "  Tolerance: " << GRADIENT_TOLERANCE << " Eh/Bohr" << std::endl;
    std::cout << "  Gradient accuracy: " << (all_passed ? "✓ PASS" : "✗ FAIL") << std::endl;

    // Translational invariance check (sum of gradients should be zero)
    Vector gradient_sum = analytical_gradient.colwise().sum();
    double translation_error = gradient_sum.norm();

    std::cout << "\nTranslational Invariance:" << std::endl;
    std::cout << "  Sum of gradients: (" << std::scientific << std::setprecision(6)
              << gradient_sum(0) << ", " << gradient_sum(1) << ", " << gradient_sum(2) << ")" << std::endl;
    std::cout << "  Norm: " << translation_error << " Eh/Bohr" << std::endl;
    std::cout << "  " << (translation_error < 1e-10 ? "✓ PASS" : "✗ FAIL") << std::endl;

    std::cout << "\n========================================" << std::endl;
    if (all_passed && translation_error < 1e-10) {
        std::cout << "✓ Monosaccharide ATM Gradient Test PASSED" << std::endl;
    } else {
        std::cout << "✗ Monosaccharide ATM Gradient Test FAILED" << std::endl;
        if (!all_passed) {
            std::cout << "  Reason: Gradient accuracy below tolerance" << std::endl;
        }
        if (translation_error >= 1e-10) {
            std::cout << "  Reason: Translational invariance violated" << std::endl;
        }
    }
    std::cout << "========================================" << std::endl;

    return all_passed && (translation_error < 1e-10);
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    // Initialize CurcumaLogger with environment variable or default to 0
    const char* env_verbosity = std::getenv("CUR_VERBOSITY");
    int verbosity = env_verbosity ? std::atoi(env_verbosity) : 0;
    CurcumaLogger::set_verbosity(verbosity);

    std::cout << "========================================" << std::endl;
    std::cout << "ATM 3-Body Dispersion Gradient Tests" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Validates analytical gradients against numerical finite differences" << std::endl;
    std::cout << "Reference: external/cpp-d4/src/damping/atm.cpp" << std::endl;
    std::cout << std::endl;

    bool all_passed = true;

    // Test 1: HCl (2 atoms - not applicable for ATM, needs 3+)
    all_passed = testATMGradients_HCl() && all_passed;

    // Test 2: Methane (5 atoms, 10 triples)
    all_passed = testATMGradients_Methane() && all_passed;

    // Test 3: Monosaccharide (27 atoms, 2925 triples - large molecule test)
    all_passed = testATMGradients_Monosaccharide() && all_passed;

    std::cout << "\n========================================" << std::endl;
    std::cout << "FINAL RESULT" << std::endl;
    std::cout << "========================================" << std::endl;
    if (all_passed) {
        std::cout << "✓ All ATM gradient tests PASSED" << std::endl;
        std::cout << "========================================" << std::endl;
        return 0;
    } else {
        std::cout << "✗ Some ATM gradient tests FAILED" << std::endl;
        std::cout << "========================================" << std::endl;
        return 1;
    }
}
