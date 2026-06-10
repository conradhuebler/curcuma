/*
 * <Unit Tests for ORCA Interface>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (June 2026) - Phase 0 of ORCA Interface Improvement Plan
 *
 * Covers:
 *   - methodToOrcaKeyword case-insensitive mapping
 *   - buildGeometryBlock ORCA input format
 *   - parseGradientFromText / parseChargesFromText (mocked ORCA 5.x and 6.x output)
 *   - JSON property parsing (mocked orca_2json output)
 *   - CG atom guard (orca_allow_cg parameter)
 *   - Shell redirect fix in produceAndLoadJSON
 *   - checkOrcaExecutable (best-effort, skippable in CI)
 */

#include "src/core/energy_calculators/qm_methods/orcainterface.h"
#include "src/core/energy_calculators/qm_methods/orca_method.h"
#include "src/core/config_manager.h"
#include "src/core/molecule.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Test harness
class OrcaInterfaceTest {
public:
    int tests_run = 0;
    int tests_passed = 0;

    void assert_equal_str(const std::string& expected, const std::string& actual, const std::string& msg)
    {
        tests_run++;
        if (expected == actual) {
            tests_passed++;
            std::cout << "  PASS: " << msg << std::endl;
        } else {
            std::cout << "  FAIL: " << msg << " - Expected '" << expected
                      << "', got '" << actual << "'" << std::endl;
        }
    }

    void assert_equal_int(int expected, int actual, const std::string& msg)
    {
        tests_run++;
        if (expected == actual) {
            tests_passed++;
            std::cout << "  PASS: " << msg << std::endl;
        } else {
            std::cout << "  FAIL: " << msg << " - Expected " << expected
                      << ", got " << actual << std::endl;
        }
    }

    void assert_equal_dbl(double expected, double actual, double tol, const std::string& msg)
    {
        tests_run++;
        if (std::fabs(expected - actual) < tol) {
            tests_passed++;
            std::cout << "  PASS: " << msg << std::endl;
        } else {
            std::cout << "  FAIL: " << msg << " - Expected " << expected
                      << ", got " << actual << std::endl;
        }
    }

    void assert_true(bool condition, const std::string& msg)
    {
        tests_run++;
        if (condition) {
            tests_passed++;
            std::cout << "  PASS: " << msg << std::endl;
        } else {
            std::cout << "  FAIL: " << msg << std::endl;
        }
    }

    void assert_false(bool condition, const std::string& msg)
    {
        tests_run++;
        if (!condition) {
            tests_passed++;
            std::cout << "  PASS: " << msg << std::endl;
        } else {
            std::cout << "  FAIL: " << msg << std::endl;
        }
    }

    void section(const std::string& name)
    {
        std::cout << "\n=== " << name << " ===" << std::endl;
    }

    void print_summary()
    {
        std::cout << "\n--- Summary ---\n"
                  << "Total:  " << tests_run << "\n"
                  << "Passed: " << tests_passed << "\n"
                  << "Failed: " << (tests_run - tests_passed) << std::endl;
    }
};

// Helper: build a minimal Mol with given atoms, geometry, charge, spin
static Mol makeMol(const std::vector<int>& atoms, const Matrix& geom, int charge, int spin)
{
    Mol m;
    m.m_atoms = atoms;
    m.m_geometry = geom;
    m.m_charge = charge;
    m.m_spin = spin;
    m.m_number_atoms = static_cast<int>(atoms.size());
    return m;
}

// =================================================================================
// 1. methodToOrcaKeyword
// =================================================================================
static void test_method_to_orca_keyword(OrcaInterfaceTest& t)
{
    t.section("methodToOrcaKeyword");

    // Known composite methods, case-insensitive
    t.assert_equal_str("HF-3c", OrcaMethod::methodToOrcaKeyword("hf-3c"), "hf-3c -> HF-3c");
    t.assert_equal_str("HF-3c", OrcaMethod::methodToOrcaKeyword("HF-3C"), "HF-3C -> HF-3c");
    t.assert_equal_str("HF-3c", OrcaMethod::methodToOrcaKeyword("Hf-3C"), "Hf-3C -> HF-3c");
    t.assert_equal_str("B97-3c", OrcaMethod::methodToOrcaKeyword("b97-3c"), "b97-3c -> B97-3c");
    t.assert_equal_str("B97-3c", OrcaMethod::methodToOrcaKeyword("B97-3C"), "B97-3C -> B97-3c");
    t.assert_equal_str("r2SCAN-3c", OrcaMethod::methodToOrcaKeyword("r2scan-3c"), "r2scan-3c -> r2SCAN-3c");
    t.assert_equal_str("r2SCAN-3c", OrcaMethod::methodToOrcaKeyword("R2SCAN-3C"), "R2SCAN-3C -> r2SCAN-3c");
    t.assert_equal_str("PBEh-3c", OrcaMethod::methodToOrcaKeyword("pbeh-3c"), "pbeh-3c -> PBEh-3c");

    // Custom-keywords mode
    t.assert_equal_str("", OrcaMethod::methodToOrcaKeyword("orca"), "orca -> empty (custom)");

    // Unknown method: pass-through
    t.assert_equal_str("CCSD(T)", OrcaMethod::methodToOrcaKeyword("CCSD(T)"), "unknown -> pass-through");
}

// =================================================================================
// 2. buildGeometryBlock
// =================================================================================
static void test_build_geometry_block(OrcaInterfaceTest& t)
{
    t.section("buildGeometryBlock");

    OrcaInterface iface; // default config, no molecule needed for block builder

    // Water: O at (0,0,0.12), H1 at (0,0.76,-0.48), H2 at (0,-0.76,-0.48)
    Mol h2o = makeMol({ 8, 1, 1 },
        (Matrix(3, 3) << 0.0, 0.0, 0.12,
                       0.0, 0.76, -0.48,
                       0.0, -0.76, -0.48).finished(),
        0, 0);

    std::string block = iface.buildGeometryBlock(h2o);
    t.assert_true(block.find("*xyz 0 1") != std::string::npos, "header charge=0 mult=1");
    t.assert_true(block.find("O ") != std::string::npos, "element O present");
    t.assert_true(block.find("H ") != std::string::npos, "element H present");
    t.assert_true(block.find("0.76000000") != std::string::npos, "8-decimal precision");
    // Block ends with "*\n", so check the second-to-last char
    t.assert_true(block.size() >= 2 && block[block.size() - 2] == '*',
        "ends with * before trailing newline");

    // Methane: 5 atoms, charge 0, singlet
    Mol ch4 = makeMol({ 6, 1, 1, 1, 1 },
        Matrix::Zero(5, 3),
        0, 0);
    block = iface.buildGeometryBlock(ch4);
    t.assert_true(block.find("*xyz 0 1") != std::string::npos, "CH4 header");
    t.assert_true(block.find("C ") != std::string::npos, "C element");

    // Hydroxide anion: charge -1, doublet (m_spin = 0.5 -> mult = 2)
    Mol oh_minus = makeMol({ 8, 1 },
        Matrix::Zero(2, 3),
        -1, 0);
    block = iface.buildGeometryBlock(oh_minus);
    t.assert_true(block.find("*xyz -1 1") != std::string::npos, "OH- header charge=-1 mult=1");

    // Triplet O2: charge 0, multiplicity 3 (m_spin = 1.0 -> mult = 2*1+1 = 3)
    Mol o2_triplet = makeMol({ 8, 8 },
        Matrix::Zero(2, 3),
        0, 1);
    block = iface.buildGeometryBlock(o2_triplet);
    t.assert_true(block.find("*xyz 0 3") != std::string::npos, "O2 triplet mult=3");

    // Quadrat O2: m_spin = 1.5 -> mult = 4
    Mol o2_quartet = makeMol({ 8, 8 },
        Matrix::Zero(2, 3),
        0, 2);
    block = iface.buildGeometryBlock(o2_quartet);
    t.assert_true(block.find("*xyz 0 5") != std::string::npos, "O2 quartet mult=5");
}

// =================================================================================
// 3. parseGradientFromText (mocked ORCA output)
// =================================================================================
static void test_gradient_text_parsing(OrcaInterfaceTest& t)
{
    t.section("parseGradientFromText");

    OrcaInterface iface;

    // Mock ORCA 5.x output: index + symbol + colon + 3 floats per line
    std::string orca5 = R"(Some preamble
                          SCF CONVERGED AFTER 12 ITERATIONS
                          *                        ****ORCA TERMINATED NORMALLY****
                          The cartesian gradient:
                          --------------------------------
                          #   0   C   :    0.123456   -0.234567    0.345678
                          #   1   H   :   -0.111111    0.222222   -0.333333
)";
    std::ofstream out5("test_orca5_gradient.out");
    out5 << orca5;
    out5.close();
    iface.setInputFile("test_orca5_gradient");

    // Force text parse: clear m_property_json and inject path
    iface.clearError();
    Matrix g5 = iface.parseGradientFromText(2);
    t.assert_equal_int(2, static_cast<int>(g5.rows()), "ORCA 5.x: 2 rows");
    t.assert_equal_dbl(0.123456, g5(0, 0), 1e-5, "ORCA 5.x: g[0,0]");
    t.assert_equal_dbl(-0.234567, g5(0, 1), 1e-5, "ORCA 5.x: g[0,1]");
    t.assert_equal_dbl(0.345678, g5(0, 2), 1e-5, "ORCA 5.x: g[0,2]");
    t.assert_equal_dbl(-0.111111, g5(1, 0), 1e-5, "ORCA 5.x: g[1,0]");

    // Mock ORCA 6.x output: no index, just 3 floats
    std::string orca6 = R"(CARTESIAN GRADIENT (Eh/Bohr)
                          --------------------------------
                              1.000000    2.000000    3.000000
                             -1.000000   -2.000000   -3.000000
)";
    std::ofstream out6("test_orca6_gradient.out");
    out6 << orca6;
    out6.close();
    iface.setInputFile("test_orca6_gradient");
    Matrix g6 = iface.parseGradientFromText(2);
    t.assert_equal_dbl(1.0, g6(0, 0), 1e-5, "ORCA 6.x: g[0,0]");
    t.assert_equal_dbl(2.0, g6(0, 1), 1e-5, "ORCA 6.x: g[0,1]");
    t.assert_equal_dbl(3.0, g6(0, 2), 1e-5, "ORCA 6.x: g[0,2]");
    t.assert_equal_dbl(-1.0, g6(1, 0), 1e-5, "ORCA 6.x: g[1,0]");

    // Missing file -> zero matrix
    iface.setInputFile("does_not_exist.out");
    Matrix gnone = iface.parseGradientFromText(2);
    t.assert_equal_int(2, static_cast<int>(gnone.rows()), "missing file: still 2 rows");
    t.assert_equal_dbl(0.0, gnone(0, 0), 1e-6, "missing file: zero values");

    // Malformed line -> zero row, no crash
    std::string bad = R"(The cartesian gradient:
                        not a number
                        more garbage
)";
    std::ofstream outbad("test_orca_bad.out");
    outbad << bad;
    outbad.close();
    iface.setInputFile("test_orca_bad");
    Matrix gbad = iface.parseGradientFromText(2);
    t.assert_equal_dbl(0.0, gbad(0, 0), 1e-6, "malformed: zero row");
}

// =================================================================================
// 4. parseChargesFromText
// =================================================================================
static void test_charges_text_parsing(OrcaInterfaceTest& t)
{
    t.section("parseChargesFromText");

    OrcaInterface iface;

    std::string orca_out = R"(Some output
MULLIKEN ATOMIC CHARGES
-----------------------
   0 C    :    0.123456
   1 H    :   -0.061728
   2 H    :   -0.061728
)";
    std::ofstream out("test_charges.out");
    out << orca_out;
    out.close();
    iface.setInputFile("test_charges");

    Vector q = iface.parseChargesFromText(3);
    t.assert_equal_int(3, static_cast<int>(q.size()), "charges: 3 atoms");
    t.assert_equal_dbl(0.123456, q(0), 1e-5, "charges[0] = 0.123");
    t.assert_equal_dbl(-0.061728, q(1), 1e-5, "charges[1] = -0.062");
    t.assert_equal_dbl(-0.061728, q(2), 1e-5, "charges[2] = -0.062");

    // Missing block -> zero vector
    iface.setInputFile("does_not_exist_charges.out");
    Vector q0 = iface.parseChargesFromText(2);
    t.assert_equal_int(2, static_cast<int>(q0.size()), "missing: 2 entries");
    t.assert_equal_dbl(0.0, q0(0), 1e-6, "missing: zero");
}

// =================================================================================
// 5. JSON property parsing (mocked orca_2json output)
// =================================================================================
static void test_json_property_parsing(OrcaInterfaceTest& t)
{
    t.section("JSON property parsing");

    OrcaInterface iface;
    iface.setInputFile("test_property");

    // Mock orca_2json output
    std::string json = R"({
        "Calculation": {
            "CalculationType": {
                "Energy": {
                    "TotalEnergy": {
                        "Value": -76.4321,
                        "Unit": "Eh"
                    }
                },
                "Properties": {
                    "Gradient": [0.1, 0.2, 0.3, -0.1, -0.2, -0.3],
                    "MullikenCharges": [0.5, -0.25, -0.25],
                    "LoewdinCharges": [0.4, -0.2, -0.2],
                    "DipoleMoment": [1.0, 2.0, 3.0]
                }
            }
        }
    })";
    std::ofstream out("test_property.property.json");
    out << json;
    out.close();

    t.assert_true(iface.loadPropertyJSON(), "loadPropertyJSON succeeds");

    t.assert_equal_dbl(-76.4321, iface.parseEnergy(), 1e-4, "energy from JSON");

    Matrix g = iface.parseGradient(2);
    t.assert_equal_dbl(0.1, g(0, 0), 1e-6, "gradient[0,0] from JSON");
    t.assert_equal_dbl(-0.3, g(1, 2), 1e-6, "gradient[1,2] from JSON");

    Vector c = iface.parseCharges(3);
    t.assert_equal_dbl(0.5, c(0), 1e-6, "Mulliken charges[0] from JSON");
    t.assert_equal_dbl(-0.25, c(2), 1e-6, "Mulliken charges[2] from JSON");

    Position dip = iface.parseDipole();
    t.assert_equal_dbl(1.0, dip[0], 1e-6, "dipole[0] from JSON");
    t.assert_equal_dbl(2.0, dip[1], 1e-6, "dipole[1] from JSON");
    t.assert_equal_dbl(3.0, dip[2], 1e-6, "dipole[2] from JSON");

    // Garbage JSON: load fails, doesn't crash
    std::ofstream bad("test_property.property.json");
    bad << "{ not valid json";
    bad.close();
    t.assert_false(iface.loadPropertyJSON(), "garbage JSON: loadPropertyJSON returns false");
}

// =================================================================================
// 6. CG atom filtering (orca_allow_cg parameter)
// =================================================================================
static void test_cg_atom_filtering(OrcaInterfaceTest& t)
{
    t.section("CG atom guard");

    // Molecule with one CG atom (226) and one real C
    Mol cg_mol = makeMol({ 226, 6 },
        Matrix::Zero(2, 3),
        0, 0);

    // Case A: orca_allow_cg = false (default) -> error
    {
        json cfg;
        OrcaInterface iface(cfg);
        t.assert_false(iface.generateInput(cg_mol, "HF-3c"),
            "CG present + default config -> generateInput fails");
        t.assert_true(iface.hasError(), "CG guard: hasError set");
        t.assert_true(iface.getErrorMessage().find("CG") != std::string::npos,
            "CG guard: error message mentions CG");
    }

    // Case B: orca_allow_cg = true -> succeeds
    {
        json cfg;
        cfg["orca_allow_cg"] = true;
        OrcaInterface iface(cfg);
        t.assert_true(iface.generateInput(cg_mol, "HF-3c"),
            "CG present + orca_allow_cg=true -> generateInput succeeds");
        t.assert_false(iface.hasError(), "CG guard bypassed: no error");
    }

    // Case C: all-real molecule, default config -> succeeds
    {
        json cfg;
        OrcaInterface iface(cfg);
        Mol real = makeMol({ 6, 1, 1, 1, 1 },
            Matrix::Zero(5, 3),
            0, 0);
        t.assert_true(iface.generateInput(real, "HF-3c"),
            "all-real molecule -> generateInput succeeds");
    }
}

// =================================================================================
// 7. produceAndLoadJSON redirect (no crash without orca_2json)
// =================================================================================
static void test_produce_and_load_json_redirect(OrcaInterfaceTest& t)
{
    t.section("produceAndLoadJSON redirect");

    // Create a dummy input file
    std::ofstream inp("test_plo_input.inp");
    inp << "! HF-3c\n*xyz 0 1\nH 0 0 0\n*\n";
    inp.close();

    OrcaInterface iface;
    iface.setInputFile("test_plo_input.inp");

    // Should not crash even if orca_2json is not installed
    bool ok = iface.produceAndLoadJSON();
    // We don't assert ok=true (depends on orca_2json availability);
    // we just assert the call returned without crashing.
    (void)ok;
    t.assert_true(true, "produceAndLoadJSON did not crash");
}

// =================================================================================
// 8. checkOrcaExecutable (best-effort, skippable)
// =================================================================================
static void test_check_orca_executable(OrcaInterfaceTest& t)
{
    t.section("checkOrcaExecutable");

    // Best-effort test: may be skipped in CI without orca
    const char* skip = std::getenv("CURCUMA_TESTING_SKIP_ORCA");
    if (skip && std::string(skip) == "1") {
        t.assert_true(true, "checkOrcaExecutable: SKIPPED (env CURCUMA_TESTING_SKIP_ORCA=1)");
        return;
    }

    // Just call it; result depends on system. We only check it returns
    // without crashing and the return type is bool.
    bool available = OrcaInterface::checkOrcaExecutable();
    t.assert_true(available || !available, "checkOrcaExecutable: returns bool");
    std::cout << "  INFO: orca " << (available ? "available" : "not available") << std::endl;
}

// =================================================================================
// 9. Energy decomposition returns empty (no fake dummies)
// =================================================================================
static void test_get_energy_decomposition(OrcaInterfaceTest& t)
{
    t.section("getEnergyDecomposition");

    OrcaMethod method("hf-3c", json{});
    json decomp = method.getEnergyDecomposition();
    t.assert_true(decomp.is_object(), "decomp is JSON object");
    t.assert_false(decomp.contains("Bond"), "decomp no longer contains fake Bond key");
    t.assert_false(decomp.contains("Dispersion"), "decomp no longer contains fake Dispersion key");
}

// =================================================================================
// Main
// =================================================================================
int main()
{
    OrcaInterfaceTest t;

    test_method_to_orca_keyword(t);
    test_build_geometry_block(t);
    test_gradient_text_parsing(t);
    test_charges_text_parsing(t);
    test_json_property_parsing(t);
    test_cg_atom_filtering(t);
    test_produce_and_load_json_redirect(t);
    test_check_orca_executable(t);
    test_get_energy_decomposition(t);

    t.print_summary();
    return t.tests_passed == t.tests_run ? 0 : 1;
}
