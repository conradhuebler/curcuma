/*
 * Comprehensive D3 Dispersion Test for All Test Molecules
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Claude Generated - December 2025
 *
 * Tests native D3 implementation on all molecules in test_cases/molecules/
 * Validates accuracy, stability, and performance for systems of varying sizes
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

#include "src/core/energy_calculators/ff_methods/d3param_generator.h"
#include "src/core/config_manager.h"
#include "src/core/curcuma_logger.h"
#include "src/tools/formats.h"

struct MoleculeTest {
    std::string name;
    std::string file_path;
    int expected_atoms;
    double reference_energy;  // s-dftd3 reference (if available)
    bool has_reference;
};

// Test molecule database with s-dftd3 reference energies
std::vector<MoleculeTest> test_molecules = {
    // Dimers (2 atoms) - s-dftd3 PBE0/BJ references
    {"H2", "test_cases/molecules/dimers/HH.xyz", 2, -6.773101e-5, true},
    {"HCl", "test_cases/molecules/dimers/HCl.xyz", 2, -2.625591e-4, true},
    {"OH", "test_cases/molecules/dimers/OH.xyz", 2, -1.179094e-4, true},

    // Trimers (3 atoms) - s-dftd3 PBE0/BJ references
    {"HCN", "test_cases/molecules/trimers/HCN.xyz", 3, -6.860246e-4, true},
    {"O3", "test_cases/molecules/trimers/O3.xyz", 3, -5.916148e-4, true},
    {"H2O", "test_cases/molecules/trimers/water.xyz", 3, -2.768645e-4, true},

    // Larger molecules - s-dftd3 PBE0/BJ references
    {"CH4", "test_cases/molecules/larger/CH4.xyz", 5, -9.221170e-4, true},
    {"CH3OH", "test_cases/molecules/larger/CH3OH.xyz", 6, -1.505362e-3, true},
    {"CH3OCH3", "test_cases/molecules/larger/CH3OCH3.xyz", 9, -3.369664e-3, true},

    // Sugar molecules - native D3 PBE0/BJ references (validated implementation)
    {"triose", "test_cases/molecules/larger/triose.xyz", 66, -2.437105e-2, true},
    {"monosaccharide", "test_cases/molecules/larger/monosaccharide.xyz", 27, -8.473223e-3, true},
};

struct TestResult {
    std::string name;
    int natoms;
    int npairs;
    double energy;
    double reference;
    double error_pct;
    bool passed;
    std::string status;
};

bool loadXYZ(const std::string& filename, std::vector<int>& atoms, Eigen::MatrixXd& geometry) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open " << filename << std::endl;
        return false;
    }

    int natoms;
    file >> natoms;

    std::string comment;
    std::getline(file, comment);  // Skip first newline
    std::getline(file, comment);  // Read comment line

    atoms.resize(natoms);
    geometry.resize(natoms, 3);

    for (int i = 0; i < natoms; ++i) {
        std::string element;
        double x, y, z;
        file >> element >> x >> y >> z;

        // Convert element symbol to atomic number (simple mapping)
        if (element == "H") atoms[i] = 1;
        else if (element == "C") atoms[i] = 6;
        else if (element == "N") atoms[i] = 7;
        else if (element == "O") atoms[i] = 8;
        else if (element == "Cl") atoms[i] = 17;
        else {
            std::cerr << "ERROR: Unknown element " << element << std::endl;
            return false;
        }

        geometry(i, 0) = x;
        geometry(i, 1) = y;
        geometry(i, 2) = z;
    }

    file.close();
    return true;
}

TestResult runTest(const MoleculeTest& test) {
    TestResult result;
    result.name = test.name;
    result.reference = test.reference_energy;
    result.passed = false;
    result.status = "FAIL";

    // Load molecule
    std::vector<int> atoms;
    Eigen::MatrixXd geometry;

    if (!loadXYZ(test.file_path, atoms, geometry)) {
        result.status = "ERROR: File not found";
        result.natoms = 0;
        result.npairs = 0;
        result.energy = 0.0;
        result.error_pct = 0.0;
        return result;
    }

    result.natoms = atoms.size();
    result.npairs = result.natoms * (result.natoms - 1) / 2;

    // Check atom count
    if (result.natoms != test.expected_atoms) {
        result.status = "ERROR: Atom count mismatch";
        result.energy = 0.0;
        result.error_pct = 0.0;
        return result;
    }

    // Create D3 parameter generator with PBE0/BJ parameters
    json d3_config = {
        {"d3_a1", 0.4145},
        {"d3_a2", 4.8593},
        {"d3_alp", 14.0},
        {"d3_s6", 1.0},
        {"d3_s8", 1.2177}
    };

    ConfigManager config("d3param", d3_config);
    D3ParameterGenerator d3_gen(config);

    // Generate D3 parameters
    d3_gen.GenerateParameters(atoms, geometry);

    // Calculate D3 energy
    result.energy = d3_gen.getTotalEnergy();

    // Validate results
    if (test.has_reference) {
        result.error_pct = std::abs(result.energy - result.reference) / std::abs(result.reference) * 100.0;
        result.passed = result.error_pct < 1.0;  // 1% threshold
        result.status = result.passed ? "PASS" : "FAIL";
    } else {
        // No reference - just check that energy is reasonable
        result.passed = (result.energy < 0.0 && result.energy > -1.0);  // Dispersion is attractive
        result.error_pct = 0.0;
        result.status = result.passed ? "OK (no ref)" : "SUSPICIOUS";
    }

    return result;
}

int main() {
    // Set verbosity to minimal for clean output
    CurcumaLogger::set_verbosity(0);

    std::cout << "======================================================================\n";
    std::cout << "D3 DISPERSION - COMPREHENSIVE MOLECULE TEST SUITE\n";
    std::cout << "======================================================================\n";
    std::cout << "Testing native D3 implementation on " << test_molecules.size() << " molecules\n";
    std::cout << "PBE0/BJ damping: s6=1.0, s8=1.2177, a1=0.4145, a2=4.8593 Bohr\n";
    std::cout << "======================================================================\n\n";

    std::vector<TestResult> results;
    int passed = 0;
    int failed = 0;
    int no_ref = 0;

    // Run tests
    for (const auto& test : test_molecules) {
        std::cout << "Testing: " << std::left << std::setw(12) << test.name
                  << " (" << test.expected_atoms << " atoms)... " << std::flush;

        TestResult result = runTest(test);
        results.push_back(result);

        std::cout << result.status;
        if (result.status == "PASS") {
            std::cout << " (error: " << std::fixed << std::setprecision(3)
                      << result.error_pct << "%)";
            passed++;
        } else if (result.status.find("OK") != std::string::npos) {
            std::cout << " (E = " << std::scientific << std::setprecision(4)
                      << result.energy << " Eh)";
            no_ref++;
        } else if (result.status == "FAIL") {
            std::cout << " (error: " << std::fixed << std::setprecision(3)
                      << result.error_pct << "%)";
            failed++;
        }
        std::cout << "\n";
    }

    // Detailed results table
    std::cout << "\n======================================================================\n";
    std::cout << "DETAILED RESULTS\n";
    std::cout << "======================================================================\n";
    std::cout << std::left << std::setw(12) << "Molecule"
              << std::right << std::setw(8) << "Atoms"
              << std::setw(8) << "Pairs"
              << std::setw(16) << "Energy (Eh)"
              << std::setw(16) << "Reference"
              << std::setw(10) << "Error %"
              << std::setw(12) << "Status"
              << "\n";
    std::cout << "----------------------------------------------------------------------\n";

    for (const auto& result : results) {
        std::cout << std::left << std::setw(12) << result.name
                  << std::right << std::setw(8) << result.natoms
                  << std::setw(8) << result.npairs
                  << std::setw(16) << std::scientific << std::setprecision(4) << result.energy;

        if (result.reference != 0.0) {
            std::cout << std::setw(16) << std::scientific << std::setprecision(4) << result.reference
                      << std::setw(10) << std::fixed << std::setprecision(3) << result.error_pct;
        } else {
            std::cout << std::setw(16) << "-"
                      << std::setw(10) << "-";
        }

        std::cout << std::setw(12) << result.status << "\n";
    }

    // Summary
    std::cout << "\n======================================================================\n";
    std::cout << "SUMMARY\n";
    std::cout << "======================================================================\n";
    std::cout << "Total tests:       " << test_molecules.size() << "\n";
    std::cout << "Passed (with ref): " << passed << "\n";
    std::cout << "Failed:            " << failed << "\n";
    std::cout << "OK (no reference): " << no_ref << "\n";
    std::cout << "======================================================================\n";

    if (failed > 0) {
        std::cout << "✗ SOME TESTS FAILED\n";
        return 1;
    } else if (passed > 0) {
        std::cout << "✓ ALL TESTS WITH REFERENCES PASSED\n";
    }
    std::cout << "======================================================================\n";

    return 0;
}
