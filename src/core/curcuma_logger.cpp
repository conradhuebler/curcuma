/*
 * <Curcuma Logging System Implementation>
 * Copyright (C) 2025 Claude AI - Generated Code
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

#include "curcuma_logger.h"
#include "citation_registry.h"
#include <algorithm>
#include <cstdlib> // for getenv
#include <mutex> // Claude Generated (Jun 2026): thread-safe citation registry (parallel ConfSearch)
#include <set>
#include <unistd.h> // for isatty

// Static member definitions - Claude Generated
int CurcumaLogger::m_verbosity = 1;
bool CurcumaLogger::m_use_colors = true;
bool CurcumaLogger::m_progress_enabled = true; // Claude Generated: global progress-bar switch
CurcumaLogger::OutputFormat CurcumaLogger::m_format = CurcumaLogger::OutputFormat::TERMINAL;
std::chrono::high_resolution_clock::time_point CurcumaLogger::m_start_time;
std::unordered_map<std::string, std::string> CurcumaLogger::m_citation_registry;
std::set<std::string> CurcumaLogger::m_used_citations;

// Claude Generated (Jun 2026): the citation statics above are process-global and are touched
// concurrently when ConfSearch runs MD/opt on a CxxThreadPool (every CurcumaMethod ctor inits
// the registry, every method addCitation()s into m_used_citations). Without serialisation the
// concurrent reassign/rehash corrupts the containers -> bad_alloc / SIGSEGV. The registry text
// is constant, so build it exactly once; guard the used-set inserts/reads with a mutex.
namespace {
std::once_flag g_citation_registry_once;
std::mutex g_citation_mutex;
}

// Unit conversion functions (from global_config.h constants)
inline double hartree_to_kjmol(double eh) { return eh * 2625.4996394798; }
inline double hartree_to_ev(double eh) { return eh * 27.211386245988; }
inline double hartree_to_wavenumber(double eh) { return eh * 219474.6313632; }
inline double bohr_to_angstrom(double bohr) { return bohr * 0.529177210903; }

// ============================================================================
// Public Interface Implementation
// ============================================================================

void CurcumaLogger::initialize(int verbosity, bool auto_detect_colors)
{
    m_verbosity = verbosity;

    bool use_colors = auto_detect_colors;

    // Check environment variables for color override
    if (std::getenv("CURCUMA_NO_COLOR")) {
        use_colors = false;
    } else if (std::getenv("FORCE_COLOR")) {
        use_colors = true;
    }

    // Check if we're outputting to a terminal for color auto-detection
    if (use_colors && !isatty(fileno(stdout))) {
        use_colors = false; // Disable colors when piping to file
    }

    m_use_colors = use_colors;
    m_format = OutputFormat::TERMINAL;

    // Check for plain mode (no colors, no prefixes - like ORCA/Gaussian)
    // Claude Generated
    if (std::getenv("CURCUMA_PLAIN")) {
        m_format = OutputFormat::PLAIN;
        m_use_colors = false;
    }

    m_start_time = std::chrono::high_resolution_clock::now();
}

void CurcumaLogger::error(const std::string& msg)
{
    // Errors always visible
    log_colored(fmt::color::red, "[ERROR] ", msg, true);
}

void CurcumaLogger::warn(const std::string& msg)
{
    if (m_verbosity >= 1) {
        log_colored(fmt::color::orange, "[WARN]  ", msg);
    }
}

void CurcumaLogger::success(const std::string& msg)
{
    if (m_verbosity >= 1) {
        log_colored(fmt::color::lime_green, "[OK]    ", msg);
    }
}

// Claude Generated: Neutral result reporting for non-success/failure outcomes
void CurcumaLogger::result(const std::string& msg)
{
    if (m_verbosity >= 1) {
        log_colored(fmt::color::white, "[RESULT]", msg);
    }
}

void CurcumaLogger::info(const std::string& msg)
{
    if (m_verbosity >= 2) {
        log_plain("        " + msg);
    }
}

void CurcumaLogger::citation(const std::string& key)
{
    CitationRegistry::cite(key);
}

// Claude Generated (April 2026): Citation key registry
void CurcumaLogger::initCitationRegistry()
{
    // Populate the constant registry exactly once for the whole process. call_once blocks any
    // concurrent CurcumaMethod ctor until the first finishes, so workers never read a half-built
    // (rehashing) map. m_used_citations is cleared here too -> it accumulates across the run and
    // every method cited in any phase is printed at the end.
    std::call_once(g_citation_registry_once, []() {
    m_citation_registry = {
        {"gfnff",   "GFN-FF: Spicher & Grimme, Angew. Chem. Int. Ed. 59, 15665 (2020) — doi:10.1002/anie.202004239"},
        {"d4",      "D4: Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019) — doi:10.1063/1.5090222"},
        {"gfn2xtb", "GFN2-xTB: Bannwarth et al., Angew. Chem. Int. Ed. 58, 3628 (2019) — doi:10.1002/anie.201901016"},
        {"tblite",  "TBLite: Caldeweyher et al., J. Chem. Theory Comput. 19, 4466 (2023) — doi:10.1021/acs.jctc.3c00556"},
        {"uff",     "UFF: Rappe et al., J. Am. Chem. Soc. 114, 10024 (1992) — doi:10.1021/ja00079a032"},
        {"eht",     "EHT: Hoffmann, J. Chem. Phys. 39, 1397 (1963) — doi:10.1063/1.1734456"},
        {"lbfgs",   "L-BFGS: Nocedal & Wright, Numerical Optimization (2006), Chapter 7"},
        {"diis",    "DIIS: Pulay, Chem. Phys. Lett. 73, 393 (1980)"},
        {"rfo",     "RFO: Banerjee et al., J. Phys. Chem. 89, 52 (1985)"},
        {"molalign","Molalign: J. Chem. Inf. Model. 2023, 63, 1157 — doi:10.1021/acs.jcim.2c01187"},
        {"d3",      "D3: Grimme et al., J. Chem. Phys. 132, 154104 (2010) — doi:10.1063/1.3382344"},
        {"curcuma", "Curcuma: doi:10.5281/zenodo.4302722"},
    };
    m_used_citations.clear();
    });
}

void CurcumaLogger::addCitation(const std::string& key, const std::string& full_text)
{
    std::string display;
    bool emit = false;
    {
        // Serialise all access to the shared citation containers (called from MD/opt workers).
        std::lock_guard<std::mutex> lock(g_citation_mutex);

        // Deduplicate: only register each key once
        if (m_used_citations.count(key)) return;
        m_used_citations.insert(key);

        // If full_text is provided and key is not in registry, store it directly
        if (!full_text.empty() && m_citation_registry.find(key) == m_citation_registry.end()) {
            m_citation_registry[key] = full_text;
        }

        if (m_verbosity >= 2) {
            auto it = m_citation_registry.find(key);
            display = (it != m_citation_registry.end()) ? it->second : full_text;
            if (display.empty()) display = key;
            emit = true;
        }
    }
    // Immediate feedback at verbosity >= 2 (outside the lock: I/O only, no shared state).
    if (emit)
        log_colored(fmt::color::green, "[CITE]  ", display);
}

void CurcumaLogger::printCitations()
{
    std::lock_guard<std::mutex> lock(g_citation_mutex);
    if (m_used_citations.empty()) return;

    fmt::print("\n");
    log_colored(fmt::color::green, "", "═══════════════════════════════════════════════════════════════");
    log_colored(fmt::color::green, "", "  Please cite the following references for methods used in this run:");
    log_colored(fmt::color::green, "", "═══════════════════════════════════════════════════════════════");
    for (const auto& key : m_used_citations) {
        auto it = m_citation_registry.find(key);
        std::string display = (it != m_citation_registry.end()) ? it->second : key;
        log_colored(fmt::color::green, "  [CITE]  ", display);
    }
    log_colored(fmt::color::green, "", "═══════════════════════════════════════════════════════════════\n");
}

void CurcumaLogger::param(const std::string& key, const std::string& value)
{
    if (m_verbosity >= 2) {
        log_colored(fmt::color::cornflower_blue, "[PARAM] ", key + ": " + value);
    }
}

void CurcumaLogger::param(const std::string& key, int value)
{
    if (m_verbosity >= 2) {
        log_colored(fmt::color::cornflower_blue, "[PARAM] ", key + ": " + std::to_string(value));
    }
}

void CurcumaLogger::param(const std::string& key, double value)
{
    if (m_verbosity >= 2) {
        log_colored(fmt::color::cornflower_blue, "[PARAM] ", fmt::format("{}: {:.6g}", key, value));
    }
}

void CurcumaLogger::param(const std::string& key, bool value)
{
    if (m_verbosity >= 2) {
        log_colored(fmt::color::cornflower_blue, "[PARAM] ", key + ": " + (value ? "true" : "false"));
    }
}

void CurcumaLogger::param_table(const json& parameters, const std::string& title)
{
    if (m_verbosity >= 2) {
        if (!parameters.empty()) {
            log_colored(fmt::color::cyan, "[TABLE] ", title);
            std::string separator(title.length() + 8, '-');
            log_plain("        " + separator);

            // Find the maximum key length for alignment
            size_t max_key_length = 0;
            for (const auto& item : parameters.items()) {
                max_key_length = std::max(max_key_length, item.key().length());
            }

            // Print each parameter
            for (const auto& item : parameters.items()) {
                std::string value_str;
                if (item.value().is_string()) {
                    value_str = item.value().get<std::string>();
                } else if (item.value().is_number_integer()) {
                    value_str = std::to_string(item.value().get<int>());
                } else if (item.value().is_number_float()) {
                    value_str = fmt::format("{:.6g}", item.value().get<double>());
                } else if (item.value().is_boolean()) {
                    value_str = item.value().get<bool>() ? "true" : "false";
                } else {
                    value_str = item.value().dump();
                }

                std::string padded_key = item.key();
                padded_key.resize(max_key_length, ' ');
                log_colored(fmt::color::cornflower_blue, "        ", padded_key + " : " + value_str);
            }
            log_plain("");
        }
    }
}

void CurcumaLogger::param_comparison_table(const json& defaults, const json& controller, const std::string& title)
{
    if (m_verbosity < 1)
        return;

    log_colored(fmt::color::cyan, "[TABLE] ", title);
    std::string separator(title.length() + 8, '-');
    log_plain("        " + separator);

    // Claude Generated (June 2026): Build the rows up front so the key column can be aligned.
    // Iterate only the keys of `defaults` (the module's registered parameters) and read the
    // effective value from `controller`. This keeps the table restricted to registered params -
    // controller-only entries (the input filename, global flags, nested module/"global"/"rmsd"
    // sub-config objects that rendered as ugly "{N keys}") are no longer shown. Object/array
    // values are skipped defensively.
    struct Row {
        std::string key;
        std::string value;
        bool changed;
    };
    std::set<std::string> all_keys;
    for (const auto& item : defaults.items())
        all_keys.insert(item.key());

    std::vector<Row> rows;
    size_t key_width = 0;
    bool any_changed = false;
    for (const auto& key : all_keys) { // std::set iterates in sorted order
        const json* eff = nullptr;
        if (controller.contains(key))
            eff = &controller.at(key);
        else if (defaults.contains(key))
            eff = &defaults.at(key);
        if (eff == nullptr || eff->is_object() || eff->is_array())
            continue; // skip nested sub-configs / non-scalar meta entries

        bool changed = controller.contains(key) && defaults.contains(key)
            && defaults.at(key) != controller.at(key);
        any_changed = any_changed || changed;
        rows.push_back({ key, format_json_value(*eff), changed });
        key_width = std::max(key_width, key.size());
    }

    if (rows.empty()) {
        log_plain("        No parameters found");
        log_plain("");
        return;
    }

    for (const auto& r : rows) {
        std::string padded_key = r.key;
        padded_key.resize(key_width, ' '); // align the ' : ' separators
        std::string line = "        " + padded_key + " : " + r.value;
        if (r.changed) {
            line += "  *";
            log_colored(fmt::color::yellow, "", line);
        } else {
            log_colored(fmt::color::cornflower_blue, "", line);
        }
    }

    if (any_changed)
        log_colored(fmt::color::white, "        ", "Legend: * = modified from default");
    log_plain("");
}

void CurcumaLogger::result_raw(const std::string& data)
{
    // Raw results for scripting - no colors, no prefix
    fmt::print("{}\n", data);
}

void CurcumaLogger::header(const std::string& title)
{
    if (m_verbosity >= 2) {
        std::string separator(title.length() + 4, '=');
        log_colored(fmt::color::cyan, "", separator);
        log_colored(fmt::color::cyan, "", "  " + title);
        log_colored(fmt::color::cyan, "", separator);
    }
}

void CurcumaLogger::progress(int current, int total, const std::string& msg)
{
    if (m_verbosity >= 2) {
        double percent = (100.0 * current) / total;
        log_colored(fmt::color::yellow, "[PROG]  ",
            fmt::format("{} [{:3.0f}%] {}/{}", msg, percent, current, total));
    }
}

// Claude Generated: live in-place ASCII progress bar. Updated via carriage-return on stdout.
// No-op when the global switch is off. ASCII only (terminal-safe). NOTE: gated only on the
// global switch, not on logger verbosity (which sub-objects lower mid-task) - the caller is
// responsible for any verbosity gating.
void CurcumaLogger::progress_bar(int current, int total, const std::string& label)
{
    if (!m_progress_enabled || total <= 0)
        return;
    if (current > total)
        current = total;
    if (current < 0)
        current = 0;
    const int width = 24;
    double frac = static_cast<double>(current) / total;
    int filled = static_cast<int>(frac * width + 0.5);
    if (filled > width)
        filled = width;
    std::string bar(filled, '#');
    bar.append(width - filled, '-');
    // Trailing spaces overwrite any leftover characters from a previous, longer line.
    fmt::print("\r  {} [{}] {:3.0f}% {}/{}    ", label, bar, frac * 100.0, current, total);
    std::fflush(stdout);
}

// Claude Generated: close the progress-bar line with a newline so following output is clean.
// Gated only on the global switch (the caller decides when a bar was actually shown).
void CurcumaLogger::progress_done()
{
    if (!m_progress_enabled)
        return;
    fmt::print("\n");
    std::fflush(stdout);
}

void CurcumaLogger::energy_rel(double value_eh, const std::string& label)
{
    if (m_verbosity >= 1) {
        std::string formatted = format_energy_relative(value_eh);
        log_colored(fmt::color::cornflower_blue, "[ENERGY]", label + ": " + formatted);
    }
}

void CurcumaLogger::energy_abs(double value_eh, const std::string& label)
{
    if (m_verbosity >= 1) {
        log_colored(fmt::color::cornflower_blue, "[ENERGY]",
            fmt::format("{}: {:.8f} Eh", label, value_eh));
    }
}

void CurcumaLogger::length(double value_bohr, const std::string& label)
{
    if (m_verbosity >= 2) {
        log_colored(fmt::color::cornflower_blue, "[LENGTH]",
            fmt::format("{}: {:.4f} Å", label, bohr_to_angstrom(value_bohr)));
    }
}

void CurcumaLogger::time(double value_aut, const std::string& label)
{
    if (m_verbosity >= 2) {
        std::string formatted = format_time(value_aut * 24.188843265857); // AUT_TO_FS
        log_colored(fmt::color::cornflower_blue, "[TIME]  ", label + ": " + formatted);
    }
}

#ifdef CURCUMA_DEBUG
void CurcumaLogger::debug(int level, const std::string& msg)
{
    if (level <= 2) { // CURCUMA_DEBUG_LEVEL from global_config.h
        log_colored(fmt::color::magenta, fmt::format("[DEBUG{}] ", level), msg);
    }
}

void CurcumaLogger::debug_var(const std::string& name, const std::string& value)
{
    if (1 <= 2) { // CURCUMA_DEBUG_LEVEL
        log_colored(fmt::color::magenta, "[DEBUG] ", fmt::format("{} = {}", name, value));
    }
}

void CurcumaLogger::debug_timing(const std::string& label)
{
    if (2 <= 2) { // CURCUMA_DEBUG_LEVEL
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - m_start_time);
        log_colored(fmt::color::magenta, "[TIMING]",
            fmt::format("{}: {:.3f} ms", label, duration.count() / 1000.0));
    }
}
#endif

// ============================================================================
// Private Helper Implementation
// ============================================================================

void CurcumaLogger::log_colored(fmt::color color, const std::string& prefix, const std::string& msg, bool force)
{
    if (m_format == OutputFormat::PLAIN) {
        // Plain mode: content only, no prefixes, no colors (like ORCA/Gaussian)
        fmt::print("{}\n", msg);
    } else if (m_use_colors || force) {
        fmt::print(fmt::fg(color), "{}{}\n", prefix, msg);
    } else {
        fmt::print("{}{}\n", prefix, msg);
    }
}

void CurcumaLogger::log_plain(const std::string& msg)
{
    fmt::print("{}\n", msg);
}

std::string CurcumaLogger::format_json_value(const json& value)
{
    if (value.is_string()) {
        return value.get<std::string>();
    } else if (value.is_number_integer()) {
        return std::to_string(value.get<int>());
    } else if (value.is_number_float()) {
        return fmt::format("{:.6g}", value.get<double>());
    } else if (value.is_boolean()) {
        return value.get<bool>() ? "true" : "false";
    } else if (value.is_array()) {
        return fmt::format("[{} items]", value.size());
    } else if (value.is_object()) {
        return fmt::format("{{{}keys}}", value.size());
    } else {
        return value.dump();
    }
}

std::string CurcumaLogger::format_energy_relative(double eh_diff)
{
    double abs_diff = std::abs(eh_diff);
    if (abs_diff > 0.01) {
        return fmt::format("{:.2f} kJ/mol", hartree_to_kjmol(eh_diff));
    } else if (abs_diff > 1e-6) {
        return fmt::format("{:.1f} meV", hartree_to_ev(eh_diff) * 1000.0);
    } else {
        return fmt::format("{:.0f} cm⁻¹", hartree_to_wavenumber(eh_diff));
    }
}

std::string CurcumaLogger::format_time(double fs)
{
    if (fs < 1000.0) {
        return fmt::format("{:.1f} fs", fs);
    } else if (fs < 1e6) {
        return fmt::format("{:.1f} ps", fs / 1000.0);
    } else {
        return fmt::format("{:.1f} ns", fs / 1e6);
    }
}

std::string CurcumaLogger::get_elapsed_time()
{
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_start_time);
    return fmt::format("{:.3f} s", duration.count() / 1000.0);
}