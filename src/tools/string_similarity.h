/*
 * <String similarity utilities for error correction>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 * Claude Generated: Levenshtein distance for method name suggestions
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace StringUtils {

/**
 * @brief Educational implementation of Levenshtein distance algorithm
 *
 * Calculates the minimum number of single-character edits (insertions,
 * deletions, or substitutions) required to change one string into another.
 *
 * Uses dynamic programming with O(m*n) time and space complexity.
 *
 * @param s1 First string
 * @param s2 Second string
 * @return Minimum edit distance between the two strings
 *
 * Example: levenshtein_distance("gfn3", "gfn2") = 1 (one substitution)
 */
inline int levenshtein_distance(const std::string& s1, const std::string& s2)
{
    const size_t m = s1.size();
    const size_t n = s2.size();

    // Dynamic programming table: dp[i][j] = edit distance between s1[0..i-1] and s2[0..j-1]
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1));

    // Initialize base cases: deletion/insertion costs from empty string
    for (size_t i = 0; i <= m; i++)
        dp[i][0] = i; // Delete all characters from s1
    for (size_t j = 0; j <= n; j++)
        dp[0][j] = j; // Insert all characters into empty string

    // Fill table with minimum edit operations
    for (size_t i = 1; i <= m; i++) {
        for (size_t j = 1; j <= n; j++) {
            if (s1[i - 1] == s2[j - 1]) {
                // Characters match - no operation needed
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                // Take minimum of three operations + 1
                dp[i][j] = 1 + std::min({
                    dp[i - 1][j],     // Deletion from s1
                    dp[i][j - 1],     // Insertion into s1
                    dp[i - 1][j - 1]  // Substitution
                });
            }
        }
    }

    return dp[m][n];
}

/**
 * @brief Find closest matching strings based on Levenshtein distance
 *
 * Educational helper for suggesting alternatives when user makes typos.
 * Returns up to max_suggestions candidates with edit distance <= max_distance.
 *
 * @param input User's input string (potentially misspelled)
 * @param candidates List of valid options to compare against
 * @param max_distance Maximum edit distance to consider (default: 3)
 * @param max_suggestions Maximum number of suggestions to return (default: 3)
 * @return Vector of closest matches, sorted by edit distance
 *
 * Example: find_closest_matches("gfn3", {"gfn2", "gfn1", "uff"}, 3, 3)
 *          returns ["gfn2", "gfn1"] (distances 1 and 1)
 */
inline std::vector<std::string> find_closest_matches(
    const std::string& input,
    const std::vector<std::string>& candidates,
    int max_distance = 3,
    int max_suggestions = 3)
{
    std::vector<std::pair<int, std::string>> distances;

    // Calculate edit distance to each candidate
    for (const auto& candidate : candidates) {
        int dist = levenshtein_distance(input, candidate);
        if (dist <= max_distance) {
            distances.push_back({ dist, candidate });
        }
    }

    // Sort by distance (closest first)
    std::sort(distances.begin(), distances.end());

    // Return top N suggestions
    std::vector<std::string> matches;
    for (size_t i = 0; i < std::min(distances.size(), size_t(max_suggestions)); i++) {
        matches.push_back(distances[i].second);
    }

    return matches;
}

} // namespace StringUtils
