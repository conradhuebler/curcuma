/*
 * <Structural Tests for PolymerBuild — Prime Notation Xx Selection>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated: Tests that verify actual polymer structures produced
 * by buildSubchain with different Xx selection patterns. Compares atom
 * counts, element composition, and geometric equivalence.
 */

#include "src/capabilities/polymerbuild.h"
#include "src/core/molecule.h"
#include "src/core/curcuma_logger.h"
#include "src/core/elements.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fmt/format.h>

/// Claude Generated: Friend class to access private methods for structural testing
class TestablePolymerBuild {
public:
    /// Construct with fragment map and disabled optimization for deterministic results
    TestablePolymerBuild(const std::string& fragment_path)
        : m_pb(make_controller(fragment_path), true)
    {
    }

    std::vector<SequenceEntry> parseSequence(const std::string& s) { return m_pb.parseSequence(s); }

    SubchainResult buildSubchain(const std::vector<SequenceEntry>& entries, int offset, int idx)
    {
        return m_pb.buildSubchain(entries, offset, idx);
    }

private:
    static nlohmann::json make_controller(const std::string& fragment_path)
    {
        nlohmann::json ctrl;
        ctrl["fragments"] = {{"A", fragment_path}};
        ctrl["optimize"] = false;      // Disable FF optimization for deterministic geometry
        ctrl["dynamics"] = false;      // No MD
        ctrl["verbosity"] = 0;         // Silent
        ctrl["write_intermediates"] = false;
        return ctrl;
    }

    PolymerBuild m_pb;
};

class StructureTest {
public:
    int tests_run = 0;
    int tests_passed = 0;

    void assert_true(bool condition, const std::string& msg)
    {
        tests_run++;
        if (condition) {
            tests_passed++;
            std::cout << "  ✓ " << msg << std::endl;
        } else {
            std::cout << "  ✗ FAIL: " << msg << std::endl;
        }
    }

    void assert_equal(int expected, int actual, const std::string& msg)
    {
        tests_run++;
        if (expected == actual) {
            tests_passed++;
            std::cout << "  ✓ " << msg << std::endl;
        } else {
            std::cout << "  ✗ FAIL: " << msg
                      << " — expected " << expected << ", got " << actual << std::endl;
        }
    }

    /// Compare two geometries element-by-element with tolerance
    bool geometries_equal(const Molecule& a, const Molecule& b, double tol = 1e-4)
    {
        if (a.AtomCount() != b.AtomCount()) return false;
        for (int i = 0; i < a.AtomCount(); ++i) {
            if (a.Atom(i).first != b.Atom(i).first) return false;
            double dist = (a.Atom(i).second - b.Atom(i).second).norm();
            if (dist > tol) return false;
        }
        return true;
    }

    /// Count atoms of a given element
    int count_element(const Molecule& mol, int elem)
    {
        int count = 0;
        for (int i = 0; i < mol.AtomCount(); ++i)
            if (mol.Atom(i).first == elem) ++count;
        return count;
    }

    /// Count remaining Xx atoms (element 0)
    int count_xx(const Molecule& mol) { return count_element(mol, 0); }

    /// Print element composition
    void print_composition(const Molecule& mol, const std::string& label)
    {
        std::map<int, int> composition;
        for (int i = 0; i < mol.AtomCount(); ++i)
            composition[mol.Atom(i).first]++;
        std::cout << "    " << label << " (" << mol.AtomCount() << " atoms): ";
        for (const auto& [elem, count] : composition)
            std::cout << Elements::ElementAbbr[elem] << count << " ";
        std::cout << std::endl;
    }

    void print_summary()
    {
        std::cout << "\n============================================" << std::endl;
        std::cout << "Tests: " << tests_passed << "/" << tests_run << " passed" << std::endl;
        std::cout << "============================================" << std::endl;
    }
};

int main(int argc, char** argv)
{
    // Determine fragment file path
    // The test is run from build/test_cases/ and the fragment is copied there by CMake
    std::string fragment_path = "polymerbuild/fragment_2xx.xyz";

    // Verify fragment exists
    {
        Molecule test_mol;
        test_mol.LoadMolecule(fragment_path);
        if (test_mol.AtomCount() == 0) {
            std::cerr << "ERROR: Cannot load fragment from " << fragment_path << std::endl;
            return 1;
        }
    }

    CurcumaLogger::set_verbosity(0);

    StructureTest test;

    // ================================================================
    // 1. Single fragment: verify fragment loading
    // ================================================================
    std::cout << "\n=== Single Fragment Loading ===" << std::endl;
    {
        TestablePolymerBuild pb(fragment_path);
        auto seq = pb.parseSequence("A");
        SubchainResult result = pb.buildSubchain(seq, 0, 0);
        test.assert_equal(7, result.molecule.AtomCount(), "Single fragment A: 7 atoms (5 real + 2 Xx)");
        test.assert_equal(2, test.count_xx(result.molecule), "Single fragment A: 2 Xx atoms");
        test.print_composition(result.molecule, "A");
    }

    // ================================================================
    // 2. Two fragments with default Xx (xx=0): (A)2
    // ================================================================
    std::cout << "\n=== Two Fragments: (A)2 ===" << std::endl;
    SubchainResult result_A2;
    {
        TestablePolymerBuild pb(fragment_path);
        auto seq = pb.parseSequence("(A)2");
        result_A2 = pb.buildSubchain(seq, 0, 0);

        // 2 fragments × 7 atoms - 2 removed Xx (1 from each side of interface) = 12 atoms
        // Each connection removes 2 Xx atoms (one from polymer, one from fragment)
        test.assert_true(result_A2.molecule.AtomCount() > 7,
            fmt::format("(A)2: more than 7 atoms (got {})", result_A2.molecule.AtomCount()));
        test.print_composition(result_A2.molecule, "(A)2");

        // The remaining Xx count should be 2 (chain ends)
        int xx_remaining = test.count_xx(result_A2.molecule);
        test.assert_equal(2, xx_remaining,
            fmt::format("(A)2: 2 remaining Xx (chain ends), got {}", xx_remaining));
    }

    // ================================================================
    // 3. Equivalence: (A)2 == AA (identical structures)
    // ================================================================
    std::cout << "\n=== Equivalence: (A)2 == AA ===" << std::endl;
    {
        TestablePolymerBuild pb(fragment_path);
        auto seq = pb.parseSequence("AA");
        SubchainResult result_AA = pb.buildSubchain(seq, 0, 0);

        test.assert_equal(result_A2.molecule.AtomCount(), result_AA.molecule.AtomCount(),
            "(A)2 and AA: same atom count");
        test.assert_true(test.geometries_equal(result_A2.molecule, result_AA.molecule),
            "(A)2 and AA: identical geometry");
    }

    // ================================================================
    // 4. Three fragments: (A)3 — verify linear growth
    // ================================================================
    std::cout << "\n=== Three Fragments: (A)3 ===" << std::endl;
    SubchainResult result_A3;
    {
        TestablePolymerBuild pb(fragment_path);
        auto seq = pb.parseSequence("(A)3");
        result_A3 = pb.buildSubchain(seq, 0, 0);

        test.assert_true(result_A3.molecule.AtomCount() > result_A2.molecule.AtomCount(),
            fmt::format("(A)3 ({} atoms) > (A)2 ({} atoms)",
                result_A3.molecule.AtomCount(), result_A2.molecule.AtomCount()));
        test.assert_equal(2, test.count_xx(result_A3.molecule), "(A)3: 2 Xx at chain ends");
        test.print_composition(result_A3.molecule, "(A)3");

        // Verify linear atom count growth: each added fragment adds (7 - 2) = 5 atoms
        // (fragment atoms minus the 2 removed Xx at interface)
        int atoms_per_addition = result_A3.molecule.AtomCount() - result_A2.molecule.AtomCount();
        int atoms_first_addition = result_A2.molecule.AtomCount() - 7;  // 7 = single fragment
        test.assert_equal(atoms_per_addition, atoms_first_addition,
            fmt::format("Linear growth: each addition adds {} atoms", atoms_per_addition));
    }

    // ================================================================
    // 5. Equivalence: (A)3 == AAA == (AAA)1
    // ================================================================
    std::cout << "\n=== Equivalence: (A)3 == AAA == (AAA)1 ===" << std::endl;
    {
        TestablePolymerBuild pb1(fragment_path);
        auto seq1 = pb1.parseSequence("AAA");
        SubchainResult r1 = pb1.buildSubchain(seq1, 0, 0);

        TestablePolymerBuild pb2(fragment_path);
        auto seq2 = pb2.parseSequence("(AAA)1");
        SubchainResult r2 = pb2.buildSubchain(seq2, 0, 0);

        test.assert_equal(result_A3.molecule.AtomCount(), r1.molecule.AtomCount(),
            "(A)3 and AAA: same atom count");
        test.assert_true(test.geometries_equal(result_A3.molecule, r1.molecule),
            "(A)3 and AAA: identical geometry");
        test.assert_true(test.geometries_equal(result_A3.molecule, r2.molecule),
            "(A)3 and (AAA)1: identical geometry");
    }

    // ================================================================
    // 6. Prime notation: (A')2 — always uses Xx#1
    // ================================================================
    std::cout << "\n=== Prime Notation: (A')2 ===" << std::endl;
    SubchainResult result_Ap2;
    {
        TestablePolymerBuild pb(fragment_path);
        auto seq = pb.parseSequence("(A')2");
        result_Ap2 = pb.buildSubchain(seq, 0, 0);

        test.assert_equal(result_A2.molecule.AtomCount(), result_Ap2.molecule.AtomCount(),
            "(A')2 and (A)2: same atom count (same fragment, different Xx)");
        test.assert_equal(2, test.count_xx(result_Ap2.molecule), "(A')2: 2 Xx at chain ends");
        test.print_composition(result_Ap2.molecule, "(A')2");
    }

    // ================================================================
    // 7. Non-equivalence: (A)2 != (A')2 (different Xx selected → different geometry)
    // ================================================================
    std::cout << "\n=== Non-Equivalence: (A)2 vs (A')2 ===" << std::endl;
    {
        // Same atom count but different geometry because different Xx is consumed
        test.assert_equal(result_A2.molecule.AtomCount(), result_Ap2.molecule.AtomCount(),
            "(A)2 and (A')2: same atom count");
        test.assert_true(!test.geometries_equal(result_A2.molecule, result_Ap2.molecule),
            "(A)2 and (A')2: different geometry (different Xx used)");
    }

    // ================================================================
    // 8. Alternating: (AA')2 vs (A)4 — different structures
    // ================================================================
    std::cout << "\n=== Alternating: (AA')2 vs (A)4 ===" << std::endl;
    {
        TestablePolymerBuild pb1(fragment_path);
        auto seq1 = pb1.parseSequence("(AA')2");
        SubchainResult r_alt = pb1.buildSubchain(seq1, 0, 0);

        TestablePolymerBuild pb2(fragment_path);
        auto seq2 = pb2.parseSequence("(A)4");
        SubchainResult r_homo = pb2.buildSubchain(seq2, 0, 0);

        test.assert_equal(r_alt.molecule.AtomCount(), r_homo.molecule.AtomCount(),
            "(AA')2 and (A)4: same atom count (same fragment, 4 monomers each)");
        test.assert_true(!test.geometries_equal(r_alt.molecule, r_homo.molecule),
            "(AA')2 and (A)4: different geometry (alternating vs uniform Xx)");
        test.print_composition(r_alt.molecule, "(AA')2");
        test.print_composition(r_homo.molecule, "(A)4");
    }

    // ================================================================
    // 9. Alternating equivalence: (AA')2 == AA'AA'
    // ================================================================
    std::cout << "\n=== Alternating Equivalence: (AA')2 == AA'AA' ===" << std::endl;
    {
        TestablePolymerBuild pb1(fragment_path);
        auto seq1 = pb1.parseSequence("(AA')2");
        SubchainResult r1 = pb1.buildSubchain(seq1, 0, 0);

        TestablePolymerBuild pb2(fragment_path);
        auto seq2 = pb2.parseSequence("AA'AA'");
        SubchainResult r2 = pb2.buildSubchain(seq2, 0, 0);

        test.assert_equal(r1.molecule.AtomCount(), r2.molecule.AtomCount(),
            "(AA')2 and AA'AA': same atom count");
        test.assert_true(test.geometries_equal(r1.molecule, r2.molecule),
            "(AA')2 and AA'AA': identical geometry");
    }

    // ================================================================
    // 10. Monomer ID tracking
    // ================================================================
    std::cout << "\n=== Monomer ID Tracking ===" << std::endl;
    {
        TestablePolymerBuild pb(fragment_path);
        auto seq = pb.parseSequence("(A)3");
        SubchainResult result = pb.buildSubchain(seq, 0, 0);

        test.assert_equal(result.molecule.AtomCount(), (int)result.atom_monomer_id.size(),
            "atom_monomer_id size matches atom count");

        // Check that there are exactly 3 distinct monomer IDs (0, 1, 2)
        std::set<int> unique_ids(result.atom_monomer_id.begin(), result.atom_monomer_id.end());
        test.assert_equal(3, (int)unique_ids.size(),
            fmt::format("(A)3: 3 distinct monomer IDs, got {}", unique_ids.size()));
    }

    // ================================================================
    // 11. Interface bonds created
    // ================================================================
    std::cout << "\n=== Interface Bonds ===" << std::endl;
    {
        TestablePolymerBuild pb(fragment_path);
        auto seq = pb.parseSequence("(A)3");
        SubchainResult result = pb.buildSubchain(seq, 0, 0);

        // 3 fragments → 2 interfaces → 2 interface bonds
        test.assert_equal(2, (int)result.interface_bonds.size(),
            "(A)3: 2 interface bonds (connecting 3 fragments)");

        // Interface bond indices must be valid
        for (const auto& [a, b] : result.interface_bonds) {
            test.assert_true(a >= 0 && a < result.molecule.AtomCount()
                          && b >= 0 && b < result.molecule.AtomCount(),
                fmt::format("Interface bond ({},{}) within valid range [0,{})",
                    a, b, result.molecule.AtomCount()));
        }
    }

    // ================================================================
    // 12. No NaN in geometry
    // ================================================================
    std::cout << "\n=== Geometry Validity ===" << std::endl;
    {
        TestablePolymerBuild pb(fragment_path);
        auto seq = pb.parseSequence("(AA')3");
        SubchainResult result = pb.buildSubchain(seq, 0, 0);

        bool has_nan = false;
        for (int i = 0; i < result.molecule.AtomCount(); ++i) {
            Position pos = result.molecule.Atom(i).second;
            if (std::isnan(pos(0)) || std::isnan(pos(1)) || std::isnan(pos(2))) {
                has_nan = true;
                break;
            }
        }
        test.assert_true(!has_nan,
            fmt::format("(AA')3: no NaN in geometry ({} atoms)", result.molecule.AtomCount()));
    }

    // ================================================================
    // Summary
    // ================================================================
    test.print_summary();
    return (test.tests_run == test.tests_passed) ? 0 : 1;
}
