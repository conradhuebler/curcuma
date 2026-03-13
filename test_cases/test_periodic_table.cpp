/*
 * <Unit Tests for Periodic Table Utility>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated (Session 5, December 2025)
 */

#include "src/core/periodic_table.h"

#include <iostream>
#include <string>

// Test utilities
class PeriodicTableTest {
public:
    int tests_run = 0;
    int tests_passed = 0;

    void assert_equal(int expected, int actual, const std::string& msg = "")
    {
        tests_run++;
        if (expected == actual) {
            tests_passed++;
            std::cout << "✓ " << msg << std::endl;
        } else {
            std::cout << "✗ FAIL: " << msg << " - Expected: " << expected
                      << ", Got: " << actual << std::endl;
        }
    }

    void assert_true(bool condition, const std::string& msg = "")
    {
        tests_run++;
        if (condition) {
            tests_passed++;
            std::cout << "✓ " << msg << std::endl;
        } else {
            std::cout << "✗ FAIL: " << msg << std::endl;
        }
    }

    void assert_false(bool condition, const std::string& msg = "")
    {
        tests_run++;
        if (!condition) {
            tests_passed++;
            std::cout << "✓ " << msg << std::endl;
        } else {
            std::cout << "✗ FAIL: " << msg << std::endl;
        }
    }

    void print_summary()
    {
        std::cout << "\n============================================" << std::endl;
        std::cout << "Tests: " << tests_passed << "/" << tests_run << " passed" << std::endl;
        std::cout << "============================================" << std::endl;
    }
};

int main()
{
    PeriodicTableTest test;

    // ========== Period 1 ==========
    std::cout << "\n=== Period 1 Elements ===" << std::endl;
    test.assert_equal(1, PeriodicTable::getGroup(1), "H: Group 1");
    test.assert_equal(8, PeriodicTable::getGroup(2), "He: Group 8");
    test.assert_equal(0, PeriodicTable::getMetalType(1), "H: Non-metal");
    test.assert_equal(0, PeriodicTable::getMetalType(2), "He: Non-metal");

    // ========== Period 2 ==========
    std::cout << "\n=== Period 2 Elements (Li-Ne) ===" << std::endl;
    test.assert_equal(1, PeriodicTable::getGroup(3), "Li: Group 1");
    test.assert_equal(2, PeriodicTable::getGroup(4), "Be: Group 2");
    test.assert_equal(3, PeriodicTable::getGroup(5), "B: Group 3");
    test.assert_equal(4, PeriodicTable::getGroup(6), "C: Group 4");
    test.assert_equal(5, PeriodicTable::getGroup(7), "N: Group 5");
    test.assert_equal(6, PeriodicTable::getGroup(8), "O: Group 6");
    test.assert_equal(7, PeriodicTable::getGroup(9), "F: Group 7");
    test.assert_equal(8, PeriodicTable::getGroup(10), "Ne: Group 8");

    test.assert_equal(1, PeriodicTable::getMetalType(3), "Li: Main group metal");
    test.assert_equal(1, PeriodicTable::getMetalType(4), "Be: Main group metal");
    test.assert_equal(0, PeriodicTable::getMetalType(6), "C: Non-metal");
    test.assert_equal(0, PeriodicTable::getMetalType(7), "N: Non-metal");
    test.assert_equal(0, PeriodicTable::getMetalType(8), "O: Non-metal");
    test.assert_equal(0, PeriodicTable::getMetalType(9), "F: Non-metal");

    // ========== Period 3 ==========
    std::cout << "\n=== Period 3 Elements (Na-Ar) ===" << std::endl;
    test.assert_equal(1, PeriodicTable::getGroup(11), "Na: Group 1");
    test.assert_equal(2, PeriodicTable::getGroup(12), "Mg: Group 2");
    test.assert_equal(3, PeriodicTable::getGroup(13), "Al: Group 3");
    test.assert_equal(4, PeriodicTable::getGroup(14), "Si: Group 4");
    test.assert_equal(5, PeriodicTable::getGroup(15), "P: Group 5");
    test.assert_equal(6, PeriodicTable::getGroup(16), "S: Group 6");
    test.assert_equal(7, PeriodicTable::getGroup(17), "Cl: Group 7");
    test.assert_equal(8, PeriodicTable::getGroup(18), "Ar: Group 8");

    test.assert_equal(1, PeriodicTable::getMetalType(11), "Na: Main group metal");
    test.assert_equal(1, PeriodicTable::getMetalType(12), "Mg: Main group metal");
    test.assert_equal(1, PeriodicTable::getMetalType(13), "Al: Main group metal");

    // ========== First Row Transition Metals ==========
    std::cout << "\n=== First Row Transition Metals (K-Kr) ===" << std::endl;
    test.assert_equal(1, PeriodicTable::getGroup(19), "K: Group 1");
    test.assert_equal(2, PeriodicTable::getGroup(20), "Ca: Group 2");
    test.assert_equal(-3, PeriodicTable::getGroup(21), "Sc: Group -3");
    test.assert_equal(-6, PeriodicTable::getGroup(24), "Cr: Group -6");
    test.assert_equal(-8, PeriodicTable::getGroup(26), "Fe: Group -8");
    test.assert_equal(-11, PeriodicTable::getGroup(29), "Cu: Group -11");
    test.assert_equal(-12, PeriodicTable::getGroup(30), "Zn: Group -12");

    test.assert_equal(1, PeriodicTable::getMetalType(19), "K: Main group metal");
    test.assert_equal(1, PeriodicTable::getMetalType(20), "Ca: Main group metal");
    test.assert_equal(2, PeriodicTable::getMetalType(21), "Sc: Transition metal");
    test.assert_equal(2, PeriodicTable::getMetalType(24), "Cr: Transition metal");
    test.assert_equal(2, PeriodicTable::getMetalType(26), "Fe: Transition metal");
    test.assert_equal(2, PeriodicTable::getMetalType(29), "Cu: Transition metal");
    test.assert_equal(2, PeriodicTable::getMetalType(30), "Zn: Transition metal");
    test.assert_equal(1, PeriodicTable::getMetalType(31), "Ga: Main group metal");

    // ========== Lanthanides ==========
    std::cout << "\n=== Lanthanide Series (La-Lu) ===" << std::endl;
    test.assert_equal(-3, PeriodicTable::getGroup(57), "La: Group -3");
    test.assert_equal(-3, PeriodicTable::getGroup(58), "Ce: Group -3");
    test.assert_equal(-3, PeriodicTable::getGroup(71), "Lu: Group -3");

    test.assert_equal(2, PeriodicTable::getMetalType(57), "La: Transition metal (Fortran classification)");
    test.assert_equal(2, PeriodicTable::getMetalType(58), "Ce: Transition metal");
    test.assert_equal(2, PeriodicTable::getMetalType(71), "Lu: Transition metal");

    // ========== Invalid Elements ==========
    std::cout << "\n=== Invalid/Out-of-range Elements ===" << std::endl;
    test.assert_equal(0, PeriodicTable::getGroup(0), "Z=0: Invalid");
    test.assert_equal(0, PeriodicTable::getGroup(87), "Z=87: Above range");
    test.assert_equal(0, PeriodicTable::getGroup(-1), "Z=-1: Negative");
    test.assert_equal(0, PeriodicTable::getMetalType(0), "Z=0: Invalid");
    test.assert_equal(0, PeriodicTable::getMetalType(87), "Z=87: Above range");

    // ========== Halogen Identification ==========
    std::cout << "\n=== Halogen Identification (Group 7) ===" << std::endl;
    test.assert_true(PeriodicTable::isHalogen(9), "F: Halogen");
    test.assert_true(PeriodicTable::isHalogen(17), "Cl: Halogen");
    test.assert_true(PeriodicTable::isHalogen(35), "Br: Halogen");
    test.assert_true(PeriodicTable::isHalogen(53), "I: Halogen");
    test.assert_true(PeriodicTable::isHalogen(85), "At: Halogen");

    test.assert_false(PeriodicTable::isHalogen(1), "H: Not halogen");
    test.assert_false(PeriodicTable::isHalogen(8), "O: Not halogen");
    test.assert_false(PeriodicTable::isHalogen(6), "C: Not halogen");
    test.assert_false(PeriodicTable::isHalogen(26), "Fe: Not halogen");

    // ========== Chalcogen Identification ==========
    std::cout << "\n=== Chalcogen Identification (Group 6) ===" << std::endl;
    test.assert_true(PeriodicTable::isChalcogen(8), "O: Chalcogen");
    test.assert_true(PeriodicTable::isChalcogen(16), "S: Chalcogen");
    test.assert_true(PeriodicTable::isChalcogen(34), "Se: Chalcogen");
    test.assert_true(PeriodicTable::isChalcogen(52), "Te: Chalcogen");
    test.assert_true(PeriodicTable::isChalcogen(84), "Po: Chalcogen");

    test.assert_false(PeriodicTable::isChalcogen(1), "H: Not chalcogen");
    test.assert_false(PeriodicTable::isChalcogen(9), "F: Not chalcogen");
    test.assert_false(PeriodicTable::isChalcogen(6), "C: Not chalcogen");
    test.assert_false(PeriodicTable::isChalcogen(24), "Cr: Not chalcogen");

    // ========== Transition Metal Identification ==========
    std::cout << "\n=== Transition Metal Identification ===" << std::endl;
    test.assert_true(PeriodicTable::isTransitionMetal(21), "Sc: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(24), "Cr: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(26), "Fe: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(29), "Cu: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(30), "Zn: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(39), "Y: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(47), "Ag: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(58), "Ce: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(71), "Lu: TM");
    test.assert_true(PeriodicTable::isTransitionMetal(57), "La: TM (Fortran classification)");

    test.assert_false(PeriodicTable::isTransitionMetal(3), "Li: Not TM");
    test.assert_false(PeriodicTable::isTransitionMetal(6), "C: Not TM");
    test.assert_false(PeriodicTable::isTransitionMetal(8), "O: Not TM");

    // ========== Main Group Metal Identification ==========
    std::cout << "\n=== Main Group Metal Identification ===" << std::endl;
    test.assert_true(PeriodicTable::isMainGroupMetal(3), "Li: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(11), "Na: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(19), "K: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(4), "Be: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(12), "Mg: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(20), "Ca: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(13), "Al: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(31), "Ga: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(50), "Sn: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(82), "Pb: Main group metal");
    test.assert_true(PeriodicTable::isMainGroupMetal(83), "Bi: Main group metal");

    test.assert_false(PeriodicTable::isMainGroupMetal(1), "H: Not main group metal");
    test.assert_false(PeriodicTable::isMainGroupMetal(6), "C: Not main group metal");
    test.assert_false(PeriodicTable::isMainGroupMetal(26), "Fe: Not main group metal");
    test.assert_false(PeriodicTable::isMainGroupMetal(58), "Ce: Not main group metal");
    test.assert_false(PeriodicTable::isMainGroupMetal(57), "La: Not main group metal (classified as TM in Fortran)");

    // Print summary
    test.print_summary();

    // Return non-zero exit code if any tests failed
    return (test.tests_run == test.tests_passed) ? 0 : 1;
}
