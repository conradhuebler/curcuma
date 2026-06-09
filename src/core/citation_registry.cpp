/*
 * <Citation Registry Implementation>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 */

#include "citation_registry.h"
#include "citation_database.h"
#include "curcuma_logger.h"

#include <fmt/color.h>
#include <fmt/core.h>
#include <fstream>
#include <set>

// Static member definitions
std::vector<std::string> CitationRegistry::m_cited_keys;
std::set<std::string> CitationRegistry::m_seen;
std::vector<std::pair<std::string, std::string>> CitationRegistry::m_subrefs;
std::mutex CitationRegistry::m_mutex;

void CitationRegistry::cite(const std::string& key, const std::string& parent)
{
    // Lock-free fast path (Claude Generated, Jun 2026): cite() is called from every
    // ComputationalMethod::calculateEnergy, i.e. on the hot path — every MD/opt step times every
    // parallel worker (e.g. GFN-FF fires 9 cites per energy eval). The registry below is already
    // mutex-safe, but taking the global lock that often needlessly serialises the workers. A
    // thread_local cache of keys this thread has already registered lets each worker skip the
    // global mutex after the first sighting, so the lock is taken at most once per (thread, key).
    // NOTE: CitationRegistry::clear() (testing only; currently unused) does NOT reset this cache —
    // if clear() is ever used for re-citation, add a generation counter to invalidate the cache.
    thread_local std::set<std::string> tls_seen;
    if (tls_seen.count(key)) return;

    {
        std::lock_guard<std::mutex> lock(m_mutex);

        // Already cited by some thread — record locally so this thread skips the lock next time.
        if (m_seen.count(key)) {
            tls_seen.insert(key);
            return;
        }

        // Look up citation data
        const Citations::CitationData* data = Citations::lookup(key);
        if (!data) {
            CurcumaLogger::warn("CitationRegistry: unknown key '" + key + "'");
            return; // not cached: an unknown key is a bug worth re-warning if it recurs
        }

        // Register
        m_cited_keys.push_back(key);
        m_seen.insert(key);

        // Track sub-reference relationship
        if (!parent.empty()) {
            m_subrefs.push_back({key, parent});
        }
    }
    tls_seen.insert(key);

    // Registration is silent; citations are printed once in the end-of-run summary.
    // (Previously logged immediately, but that cluttered the output during calculations.)
}

void CitationRegistry::printSummary()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_cited_keys.empty()) return;

    fmt::print("\n");
    fmt::print(fg(fmt::color::green), "===============================================================\n");
    fmt::print(fg(fmt::color::green), "  Please cite the following references for methods used in this run:\n");
    fmt::print(fg(fmt::color::green), "===============================================================\n");

    // Collect sub-ref children per parent
    std::set<std::string> children;
    for (const auto& [child, parent] : m_subrefs) {
        children.insert(child);
    }

    // Helper: print one entry (description + reference, two lines)
    auto printEntry = [](const std::string& prefix, const Citations::CitationData* data) {
        fmt::print(fg(fmt::color::green), "  [CITE]  {}\n", prefix + data->description);
        fmt::print(fg(fmt::color::green), "          {}\n", data->reference);
    };

    // Print top-level entries (those that are not children of another)
    for (const auto& key : m_cited_keys) {
        if (children.count(key)) continue;

        const Citations::CitationData* data = Citations::lookup(key);
        if (!data) continue;

        printEntry("", data);

        // Print sub-references of this parent (indented, without raw key)
        for (const auto& [child, parent] : m_subrefs) {
            if (parent != key) continue;
            const Citations::CitationData* child_data = Citations::lookup(child);
            if (!child_data) continue;
            printEntry("  -> ", child_data);
        }
    }

    fmt::print(fg(fmt::color::green), "===============================================================\n");
}

void CitationRegistry::writeBibTeX(const std::string& output_dir, const std::string& basename)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_cited_keys.empty()) return;

    // Derive filename
    std::string filename;
    if (!basename.empty()) {
        std::string stem = basename;
        auto dot = stem.rfind('.');
        if (dot != std::string::npos) stem = stem.substr(0, dot);
        filename = stem + "_citations.bib";
    } else {
        filename = "curcuma_citations.bib";
    }

    if (!output_dir.empty()) {
        filename = output_dir + "/" + filename;
    }

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        CurcumaLogger::warn("CitationRegistry: cannot write BibTeX to '" + filename + "'");
        return;
    }

    ofs << "% Citations generated by Curcuma\n\n";

    for (const auto& key : m_cited_keys) {
        const Citations::CitationData* data = Citations::lookup(key);
        if (!data || data->bibtex.empty()) continue;
        ofs << data->bibtex << "\n\n";
    }

    ofs.close();
    CurcumaLogger::info("BibTeX written to: " + filename);
}

void CitationRegistry::clear()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_cited_keys.clear();
    m_seen.clear();
    m_subrefs.clear();
}

bool CitationRegistry::hasKey(const std::string& key)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_seen.count(key) > 0;
}