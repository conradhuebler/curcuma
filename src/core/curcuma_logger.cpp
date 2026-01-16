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
#include <algorithm>
#include <cstdlib> // for getenv
#include <set>
#include <unistd.h> // for isatty

// Static member definitions - Claude Generated
int CurcumaLogger::m_verbosity = 1;
bool CurcumaLogger::m_use_colors = true;
CurcumaLogger::OutputFormat CurcumaLogger::m_format = CurcumaLogger::OutputFormat::TERMINAL;
std::chrono::high_resolution_clock::time_point CurcumaLogger::m_start_time;

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

void CurcumaLogger::citation(const std::string& ref)
{
    if (m_verbosity >= 2) {
        log_colored(fmt::color::green, "[CITE]  ", ref);
    }
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
    if (m_verbosity >= 1) {
        log_colored(fmt::color::cyan, "[TABLE] ", title);
        std::string separator(title.length() + 8, '-');
        log_plain("        " + separator);

        // Collect all unique keys from both JSONs
        std::set<std::string> all_keys;
        for (const auto& item : defaults.items()) {
            all_keys.insert(item.key());
        }
        for (const auto& item : controller.items()) {
            all_keys.insert(item.key());
        }

        if (all_keys.empty()) {
            log_plain("        No parameters found");
            log_plain("");
            return;
        }

        // Sort keys for consistent output
        std::vector<std::string> sorted_keys(all_keys.begin(), all_keys.end());
        std::sort(sorted_keys.begin(), sorted_keys.end());

        // Simple list format for cleaner implementation
        for (const auto& key : sorted_keys) {
            std::string display_value;
            bool is_changed = false;

            // Determine current value and check if changed
            if (controller.contains(key)) {
                const auto& ctrl_val = controller[key];
                display_value = format_json_value(ctrl_val);

                // Check if different from default
                if (defaults.contains(key) && defaults[key] != ctrl_val) {
                    is_changed = true;
                }
            } else if (defaults.contains(key)) {
                display_value = format_json_value(defaults[key]);
            } else {
                display_value = "N/A";
            }

            std::string param_line = "        " + key + " : " + display_value;
            if (is_changed) {
                param_line += "*";
                log_colored(fmt::color::yellow, "", param_line);
            } else {
                log_colored(fmt::color::cornflower_blue, "", param_line);
            }
        }

        if (std::any_of(sorted_keys.begin(), sorted_keys.end(), [&](const auto& key) {
                return controller.contains(key) && defaults.contains(key) && defaults[key] != controller[key];
            })) {
            log_colored(fmt::color::white, "        ", "Legend: * = modified from default");
        }
        log_plain("");
    }
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
    if (m_use_colors || force) {
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