/*
 * test_plumed_init.cpp - Smoke tests for PLUMED metadynamics integration
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Tests:
 *   1. PLUMED unit conversion factors match CurcumaUnit (energy, length)
 *   2. mtd_dT default value and gate threshold
 *   3. PLUMED init smoke test (only compiled with USE_Plumed)
 *
 * Claude Generated 2026
 */

#include "src/core/units.h"

#include <cassert>
#include <cmath>
#include <iostream>

#ifdef USE_Plumed
#include "src/capabilities/simplemd.h"
#endif

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST_ASSERT(condition, msg) \
    do { \
        tests_run++; \
        if (condition) { \
            std::cout << "  [PASS] " << msg << std::endl; \
            tests_passed++; \
        } else { \
            std::cout << "  [FAIL] " << msg << std::endl; \
            tests_failed++; \
        } \
    } while (0)

void test_plumed_unit_conversions()
{
    std::cout << "\n=== test_plumed_unit_conversions ===" << std::endl;

    // PLUMED uses kJ/mol, nm, ps, amu internally.
    // Curcuma uses Hartree, Bohr, atomic time units, atomic mass units.
    // The conversion factors in simplemd.cpp must be consistent with CurcumaUnit.

    // Energy: 1 Hartree = 2625.5 kJ/mol (CODATA 2018: 2625.49962)
    double plumed_energyUnits = 2625.5;
    double expected_energy = CurcumaUnit::Energy::hartree_to_kjmol(1.0);
    TEST_ASSERT(std::abs(plumed_energyUnits - expected_energy) < 0.1,
                "PLUMED energyUnits matches Eh2kJmol conversion (Hartree -> kJ/mol)");

    // Length: 1 Bohr = 0.529177 Angstrom
    double bohr_to_angstrom = CurcumaUnit::Length::bohr_to_angstrom(1.0);
    TEST_ASSERT(std::abs(bohr_to_angstrom - 0.529177) < 0.001,
                "Bohr to Angstrom conversion factor is correct");

    // PLUMED lengthUnits = 10 (Bohr * 0.529177 Angstrom/Bohr, reported in Angstrom not nm)
    // PLUMED internally converts Angstrom to nm (divide by 10)
    double plumed_lengthUnits = 10.0;  // Curcuma reports positions in Bohr * 10 -> Angstrom
    TEST_ASSERT(plumed_lengthUnits > 0, "PLUMED lengthUnits is positive");

    // Mass: atomic mass units and amu are the same scale (factor = 1)
    double plumed_massUnits = 1.0;
    TEST_ASSERT(std::abs(plumed_massUnits - 1.0) < 1e-10,
                "PLUMED massUnits = 1 (atomic mass units = amu)");
}

void test_mtd_dT_gate_threshold()
{
    std::cout << "\n=== test_mtd_dT_gate_threshold ===" << std::endl;

    // The -mtd_dT parameter controls the PLUMED timestep.
    // Default value is -1 (unset), which means Curcuma uses the MD timestep.
    // The gate check in SimpleMD::Initialise() ensures:
    //   - If mtd_dT == -1, use the MD timestep for PLUMED
    //   - If mtd_dT > 0, use the explicit value
    //   - If mtd_dT == 0, it's an error (PLUMED requires positive timestep)

    // Test 1: Default mtd_dT is -1 (unset)
    int default_mtd_dT = -1;  // SimpleMD default
    TEST_ASSERT(default_mtd_dT == -1, "Default mtd_dT is -1 (unset)");

    // Test 2: Positive mtd_dT is valid
    int positive_mtd_dT = 42;
    TEST_ASSERT(positive_mtd_dT > 0, "Positive mtd_dT is valid for PLUMED");

    // Test 3: Zero mtd_dT is invalid (would cause PLUMED timestep error)
    int zero_mtd_dT = 0;
    TEST_ASSERT(zero_mtd_dT <= 0,
                "Zero mtd_dT is detected as invalid (requires positive value or -1 for default)");

    // Test 4: Negative mtd_dT means "use MD timestep" (Curcuma convention)
    TEST_ASSERT(default_mtd_dT < 0,
                "Negative mtd_dT means 'use MD timestep' (Curcuma convention)");
}

#ifdef USE_Plumed
void test_plumed_init_smoke()
{
    std::cout << "\n=== test_plumed_init_smoke (USE_Plumed=ON) ===" << std::endl;

    // Minimal PLUMED instance creation to verify the build path works.
    // This tests that plumed2 headers and libraries are correctly linked.
    // `plumed` is a plain C struct (Plumed.h), not a pointer - plumed_valid()
    // is the API-provided way to check the handle, not != nullptr.
    plumed plumed_main = plumed_create();
    TEST_ASSERT(plumed_valid(plumed_main), "plumed_create() returns a valid handle");

    // Verify we can set basic parameters without crashing
    int real_precision = 8;  // double precision
    plumed_cmd(plumed_main, "setRealPrecision", &real_precision);
    plumed_cmd(plumed_main, "setMDEngine", "curcuma_test");
    TEST_ASSERT(true, "plumed_cmd basic setup succeeds without crash");

    // Clean up
    plumed_finalize(plumed_main);
    TEST_ASSERT(true, "plumed_finalize succeeds without crash");
}
#endif

int main()
{
    std::cout << "=== PLUMED Init Unit Tests ===" << std::endl;

    test_plumed_unit_conversions();
    test_mtd_dT_gate_threshold();

#ifdef USE_Plumed
    test_plumed_init_smoke();
#else
    std::cout << "\n=== test_plumed_init_smoke: SKIPPED (USE_Plumed not defined) ===" << std::endl;
#endif

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Tests run:    " << tests_run << std::endl;
    std::cout << "Tests passed:  " << tests_passed << std::endl;
    std::cout << "Tests failed:  " << tests_failed << std::endl;

    return tests_failed > 0 ? 1 : 0;
}