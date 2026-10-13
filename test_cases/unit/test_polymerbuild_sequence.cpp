/*
 * <Unit Tests for PolymerBuild Sequence Parser — Prime Notation>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Claude Generated: Tests for parseSequence/parseTokens with prime notation
 * for Xx connection point selection in polymer assembly.
 */

#include "src/capabilities/polymerbuild.h"

#include <iostream>
#include <string>
#include <vector>

/// Claude Generated: Friend class to access private parsing methods for testing
class TestablePolymerBuild {
public:
    TestablePolymerBuild()
        : m_pb(nlohmann::json{}, true)
    {
    }

    std::vector<SequenceEntry> parseSequence(const std::string& s) { return m_pb.parseSequence(s); }
    std::vector<SequenceEntry> parseTokens(const std::string& s) { return m_pb.parseTokens(s); }

private:
    PolymerBuild m_pb;
};

class SequenceParserTest {
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

    void assert_str(const std::string& expected, const std::string& actual, const std::string& msg)
    {
        tests_run++;
        if (expected == actual) {
            tests_passed++;
            std::cout << "  ✓ " << msg << std::endl;
        } else {
            std::cout << "  ✗ FAIL: " << msg
                      << " — expected \"" << expected << "\", got \"" << actual << "\"" << std::endl;
        }
    }

    /// Helper: format a SequenceEntry vector as human-readable string
    std::string format_seq(const std::vector<SequenceEntry>& seq)
    {
        std::string s;
        for (size_t i = 0; i < seq.size(); ++i) {
            if (i > 0) s += ", ";
            s += seq[i].fragment_name + std::string(seq[i].xx_selection, '\'');
        }
        return "[" + s + "]";
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
    SequenceParserTest test;
    TestablePolymerBuild pb;

    // ================================================================
    // 1. parseTokens — basic token parsing
    // ================================================================
    std::cout << "\n=== parseTokens: Basic Tokens ===" << std::endl;
    {
        auto tokens = pb.parseTokens("A");
        test.assert_equal(1, (int)tokens.size(), "Single token: size=1");
        if (!tokens.empty()) {
            test.assert_str("A", tokens[0].fragment_name, "Single token: name=A");
            test.assert_equal(0, tokens[0].xx_selection, "Single token: xx_selection=0");
        }
    }
    {
        auto tokens = pb.parseTokens("A'");
        test.assert_equal(1, (int)tokens.size(), "Single prime: size=1");
        if (!tokens.empty()) {
            test.assert_str("A", tokens[0].fragment_name, "Single prime: name=A");
            test.assert_equal(1, tokens[0].xx_selection, "Single prime: xx_selection=1");
        }
    }
    {
        auto tokens = pb.parseTokens("B''");
        test.assert_equal(1, (int)tokens.size(), "Double prime: size=1");
        if (!tokens.empty()) {
            test.assert_str("B", tokens[0].fragment_name, "Double prime: name=B");
            test.assert_equal(2, tokens[0].xx_selection, "Double prime: xx_selection=2");
        }
    }

    // ================================================================
    // 2. parseTokens — multi-token strings
    // ================================================================
    std::cout << "\n=== parseTokens: Multi-Token Strings ===" << std::endl;
    {
        auto tokens = pb.parseTokens("AA'");
        test.assert_equal(2, (int)tokens.size(), "AA': size=2");
        if (tokens.size() == 2) {
            test.assert_str("A", tokens[0].fragment_name, "AA'[0]: name=A");
            test.assert_equal(0, tokens[0].xx_selection, "AA'[0]: xx_selection=0");
            test.assert_str("A", tokens[1].fragment_name, "AA'[1]: name=A");
            test.assert_equal(1, tokens[1].xx_selection, "AA'[1]: xx_selection=1");
        }
    }
    {
        auto tokens = pb.parseTokens("A-B");
        test.assert_equal(2, (int)tokens.size(), "A-B with separator: size=2");
        if (tokens.size() == 2) {
            test.assert_str("A", tokens[0].fragment_name, "A-B[0]: name=A");
            test.assert_str("B", tokens[1].fragment_name, "A-B[1]: name=B");
        }
    }
    {
        auto tokens = pb.parseTokens("A'B''C");
        test.assert_equal(3, (int)tokens.size(), "A'B''C: size=3");
        if (tokens.size() == 3) {
            test.assert_equal(1, tokens[0].xx_selection, "A'B''C[0]: xx=1");
            test.assert_equal(2, tokens[1].xx_selection, "A'B''C[1]: xx=2");
            test.assert_equal(0, tokens[2].xx_selection, "A'B''C[2]: xx=0");
        }
    }

    // ================================================================
    // 3. parseSequence — basic sequences (backward compatibility)
    // ================================================================
    std::cout << "\n=== parseSequence: Basic (No Primes) ===" << std::endl;
    {
        auto seq = pb.parseSequence("A");
        test.assert_equal(1, (int)seq.size(), "Single A: size=1");
        if (!seq.empty()) {
            test.assert_str("A", seq[0].fragment_name, "Single A: name=A");
            test.assert_equal(0, seq[0].xx_selection, "Single A: xx=0");
        }
    }
    {
        auto seq = pb.parseSequence("(A)5");
        test.assert_equal(5, (int)seq.size(), "(A)5: size=5");
        for (size_t i = 0; i < seq.size(); ++i) {
            test.assert_str("A", seq[i].fragment_name, "(A)5[" + std::to_string(i) + "]: name=A");
            test.assert_equal(0, seq[i].xx_selection, "(A)5[" + std::to_string(i) + "]: xx=0");
        }
    }
    {
        auto seq = pb.parseSequence("A-B-C");
        test.assert_equal(3, (int)seq.size(), "A-B-C: size=3");
        if (seq.size() == 3) {
            test.assert_str("A", seq[0].fragment_name, "A-B-C[0]: A");
            test.assert_str("B", seq[1].fragment_name, "A-B-C[1]: B");
            test.assert_str("C", seq[2].fragment_name, "A-B-C[2]: C");
        }
    }
    {
        auto seq = pb.parseSequence("(A)3-B-(C)2");
        test.assert_equal(6, (int)seq.size(), "(A)3-B-(C)2: size=6");
        std::cout << "    Parsed: " << test.format_seq(seq) << std::endl;
        if (seq.size() == 6) {
            test.assert_str("A", seq[0].fragment_name, "[0]=A");
            test.assert_str("A", seq[1].fragment_name, "[1]=A");
            test.assert_str("A", seq[2].fragment_name, "[2]=A");
            test.assert_str("B", seq[3].fragment_name, "[3]=B");
            test.assert_str("C", seq[4].fragment_name, "[4]=C");
            test.assert_str("C", seq[5].fragment_name, "[5]=C");
        }
    }

    // ================================================================
    // 4. parseSequence — prime notation (new feature)
    // ================================================================
    std::cout << "\n=== parseSequence: Prime Notation ===" << std::endl;
    {
        // KEY TEST: (AA')2 must expand to [A, A', A, A'] (alternating)
        auto seq = pb.parseSequence("(AA')2");
        test.assert_equal(4, (int)seq.size(), "(AA')2: size=4");
        std::cout << "    Parsed: " << test.format_seq(seq) << std::endl;
        if (seq.size() == 4) {
            test.assert_str("A", seq[0].fragment_name, "(AA')2[0]: name=A");
            test.assert_equal(0, seq[0].xx_selection, "(AA')2[0]: xx=0");
            test.assert_str("A", seq[1].fragment_name, "(AA')2[1]: name=A");
            test.assert_equal(1, seq[1].xx_selection, "(AA')2[1]: xx=1");
            test.assert_str("A", seq[2].fragment_name, "(AA')2[2]: name=A");
            test.assert_equal(0, seq[2].xx_selection, "(AA')2[2]: xx=0");
            test.assert_str("A", seq[3].fragment_name, "(AA')2[3]: name=A");
            test.assert_equal(1, seq[3].xx_selection, "(AA')2[3]: xx=1");
        }
    }
    {
        // (AA')3-B'' → [A, A', A, A', A, A', B'']
        auto seq = pb.parseSequence("(AA')3-B''");
        test.assert_equal(7, (int)seq.size(), "(AA')3-B'': size=7");
        std::cout << "    Parsed: " << test.format_seq(seq) << std::endl;
        if (seq.size() == 7) {
            test.assert_equal(0, seq[0].xx_selection, "[0]: xx=0");
            test.assert_equal(1, seq[1].xx_selection, "[1]: xx=1");
            test.assert_equal(0, seq[2].xx_selection, "[2]: xx=0");
            test.assert_equal(1, seq[3].xx_selection, "[3]: xx=1");
            test.assert_equal(0, seq[4].xx_selection, "[4]: xx=0");
            test.assert_equal(1, seq[5].xx_selection, "[5]: xx=1");
            test.assert_str("B", seq[6].fragment_name, "[6]: name=B");
            test.assert_equal(2, seq[6].xx_selection, "[6]: xx=2");
        }
    }

    // ================================================================
    // 5. Equivalence: (A)4 == AAAA == (AA)2 (all xx_selection=0)
    // ================================================================
    std::cout << "\n=== parseSequence: Equivalence Tests ===" << std::endl;
    {
        auto s1 = pb.parseSequence("(A)4");
        auto s2 = pb.parseSequence("AAAA");
        auto s3 = pb.parseSequence("(AA)2");
        test.assert_equal((int)s1.size(), (int)s2.size(), "(A)4 and AAAA: same size");
        test.assert_equal((int)s1.size(), (int)s3.size(), "(A)4 and (AA)2: same size");

        bool all_match = true;
        for (size_t i = 0; i < s1.size() && i < s2.size(); ++i) {
            if (s1[i].fragment_name != s2[i].fragment_name || s1[i].xx_selection != s2[i].xx_selection)
                all_match = false;
        }
        test.assert_true(all_match, "(A)4 and AAAA: identical entries");

        all_match = true;
        for (size_t i = 0; i < s1.size() && i < s3.size(); ++i) {
            if (s1[i].fragment_name != s3[i].fragment_name || s1[i].xx_selection != s3[i].xx_selection)
                all_match = false;
        }
        test.assert_true(all_match, "(A)4 and (AA)2: identical entries");
    }

    // ================================================================
    // 6. Non-equivalence: (AA')2 != (A)4
    // ================================================================
    std::cout << "\n=== parseSequence: Non-Equivalence Tests ===" << std::endl;
    {
        auto s1 = pb.parseSequence("(AA')2");  // [A, A', A, A']
        auto s2 = pb.parseSequence("(A)4");    // [A, A, A, A]
        test.assert_equal((int)s1.size(), (int)s2.size(), "(AA')2 and (A)4: same size");

        bool has_difference = false;
        for (size_t i = 0; i < s1.size() && i < s2.size(); ++i) {
            if (s1[i].xx_selection != s2[i].xx_selection)
                has_difference = true;
        }
        test.assert_true(has_difference, "(AA')2 and (A)4: different xx_selection values");
    }

    // ================================================================
    // 7. Fragment names with underscores and digits
    // ================================================================
    std::cout << "\n=== parseSequence: Complex Fragment Names ===" << std::endl;
    {
        auto seq = pb.parseSequence("(pdmaema)3");
        test.assert_equal(3, (int)seq.size(), "(pdmaema)3: size=3");
        if (!seq.empty())
            test.assert_str("pdmaema", seq[0].fragment_name, "(pdmaema)3[0]: name=pdmaema");
    }
    {
        auto seq = pb.parseSequence("(frag_1'frag_2)2");
        test.assert_equal(4, (int)seq.size(), "(frag_1'frag_2)2: size=4");
        std::cout << "    Parsed: " << test.format_seq(seq) << std::endl;
        if (seq.size() == 4) {
            test.assert_str("frag_1", seq[0].fragment_name, "[0]: name=frag_1");
            test.assert_equal(1, seq[0].xx_selection, "[0]: xx=1");
            test.assert_str("frag_2", seq[1].fragment_name, "[1]: name=frag_2");
            test.assert_equal(0, seq[1].xx_selection, "[1]: xx=0");
        }
    }

    // ================================================================
    // 8. Edge cases
    // ================================================================
    std::cout << "\n=== parseSequence: Edge Cases ===" << std::endl;
    {
        auto seq = pb.parseSequence("");
        test.assert_equal(0, (int)seq.size(), "Empty string: size=0");
    }
    {
        auto seq = pb.parseSequence("(A)1");
        test.assert_equal(1, (int)seq.size(), "(A)1: size=1");
    }
    {
        auto seq = pb.parseSequence("A'''");
        test.assert_equal(1, (int)seq.size(), "A''': size=1");
        if (!seq.empty()) {
            test.assert_str("A", seq[0].fragment_name, "A''': name=A");
            test.assert_equal(3, seq[0].xx_selection, "A''': xx_selection=3");
        }
    }

    // ================================================================
    // 9. D&C cache key differentiation
    // ================================================================
    std::cout << "\n=== Cache Key Differentiation ===" << std::endl;
    {
        // Simulate cache key generation as done in assemblePolymer
        auto s1 = pb.parseSequence("(A)4");
        auto s2 = pb.parseSequence("(AA')2");

        auto make_key = [](const std::vector<SequenceEntry>& seq) {
            std::string key;
            for (const auto& entry : seq)
                key += entry.fragment_name + std::string(entry.xx_selection, '\'') + ",";
            return key;
        };

        std::string k1 = make_key(s1);
        std::string k2 = make_key(s2);
        std::cout << "    Key for (A)4:   " << k1 << std::endl;
        std::cout << "    Key for (AA')2: " << k2 << std::endl;
        test.assert_true(k1 != k2, "Cache keys differ for (A)4 vs (AA')2");
    }
    {
        // Same pattern must produce same cache key
        auto s1 = pb.parseSequence("(A)4");
        auto s2 = pb.parseSequence("AAAA");

        auto make_key = [](const std::vector<SequenceEntry>& seq) {
            std::string key;
            for (const auto& entry : seq)
                key += entry.fragment_name + std::string(entry.xx_selection, '\'') + ",";
            return key;
        };

        std::string k1 = make_key(s1);
        std::string k2 = make_key(s2);
        test.assert_true(k1 == k2, "Cache keys equal for (A)4 and AAAA");
    }

    // ================================================================
    // Summary
    // ================================================================
    test.print_summary();

    return (test.tests_run == test.tests_passed) ? 0 : 1;
}
