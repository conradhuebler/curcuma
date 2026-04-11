/*
 * <Citation Registry — runtime tracking of used computational methods>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Tracks which citations were used during a run, logs each on first use,
 * deduplicates, supports sub-references (e.g. D4 inside GFN-FF),
 * prints a summary at end of run, and writes a BibTeX file.
 */

#pragma once

#include <set>
#include <string>
#include <vector>
#include <utility>

class CitationRegistry {
public:
    /// Register a citation by key (looked up in Citations::database()).
    /// Logs the citation on first use (always visible, verbosity-independent).
    /// parent: non-empty if this is a sub-reference of another method
    ///         (e.g. "d4" is a sub-ref of "gfnff").
    static void cite(const std::string& key, const std::string& parent = "");

    /// Print a human-readable summary of all cited methods.
    static void printSummary();

    /// Write BibTeX file with all cited methods that have bibtex entries.
    /// filename: full path. If empty, derives from basename in CWD.
    static void writeBibTeX(const std::string& output_dir = "",
                             const std::string& basename = "");

    /// Clear registry (for testing)
    static void clear();

    /// Check if a key has been cited
    static bool hasKey(const std::string& key);

private:
    static std::vector<std::string> m_cited_keys;                        // ordered, deduped
    static std::set<std::string> m_seen;                                  // fast lookup
    static std::vector<std::pair<std::string, std::string>> m_subrefs;   // (child, parent)
};