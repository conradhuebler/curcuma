/*
 * <Citation Database — compiled-in reference data for all computational methods>
 * Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * One source of truth for all citation data. Adding a new method citation
 * means adding one entry to the database() map — no other code changes needed.
 */

#pragma once

#include <string>
#include <unordered_map>

namespace Citations {

struct CitationData {
    std::string description;   // Human-readable method name, e.g. "GFN-FF Force Field"
    std::string reference;      // Human-readable reference string
    std::string bibtex_key;    // BibTeX cite key, e.g. "spicher2020gfnff"
    std::string bibtex;        // Full BibTeX entry (can be empty initially)
};

/// Returns the full citation database. Called once, cached internally.
const std::unordered_map<std::string, CitationData>& database();

/// Look up a citation by key. Returns nullptr if not found.
const CitationData* lookup(const std::string& key);

} // namespace Citations