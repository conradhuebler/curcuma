/*
 * <Curcuma Logging System>
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

#pragma once

#include "json.hpp"
#include <chrono>
#include <fmt/color.h>
#include <fmt/core.h>
#include <string>
#include <type_traits>

using json = nlohmann::json;

class CurcumaLogger {
public:
    enum class OutputFormat {
        TERMINAL,
        RAW,
        MARKDOWN
    };

private:
    static int m_verbosity;
    static bool m_use_colors;
    static OutputFormat m_format;
    static std::chrono::high_resolution_clock::time_point m_start_time;

public:
    // Configuration methods
    static void set_verbosity(int level) { m_verbosity = level; }
    static void set_colors(bool enable) { m_use_colors = enable; }
    static void set_format(OutputFormat fmt) { m_format = fmt; }
    static int get_verbosity() { return m_verbosity; }
    static bool colors_enabled() { return m_use_colors; }

    // Initialize logger with environment detection
    static void initialize(int verbosity = 2, bool auto_detect_colors = true);

    // Core logging functions with verbosity control
    static void error(const std::string& msg);
    static void warn(const std::string& msg);
    static void success(const std::string& msg);
    static void info(const std::string& msg);
    static void citation(const std::string& ref);

    // Parameter logging with different types
    static void param(const std::string& key, const std::string& value);
    static void param(const std::string& key, int value);
    static void param(const std::string& key, double value);
    static void param(const std::string& key, bool value);

    // JSON parameter table printing
    static void param_table(const json& parameters, const std::string& title = "Parameters");
    static void param_comparison_table(const json& defaults, const json& controller, const std::string& title = "Parameter Settings");

    // Special formatting functions
    static void result_raw(const std::string& data);
    static void progress(int current, int total, const std::string& msg);
    static void header(const std::string& title);

    // Unit-aware output functions
    static void energy_rel(double value_eh, const std::string& label);
    static void energy_abs(double value_eh, const std::string& label);
    static void length(double value_bohr, const std::string& label);
    static void time(double value_aut, const std::string& label);

// Debug functions (only available in debug builds)
#ifdef CURCUMA_DEBUG
    static void debug(int level, const std::string& msg);
    static void debug_var(const std::string& name, const std::string& value);
    static void debug_timing(const std::string& label);
#endif

    // Modern template-based logging functions - Claude Generated
    template <typename... Args>
    static void error_fmt(fmt::format_string<Args...> format_str, Args&&... args)
    {
        std::string msg = fmt::format(format_str, std::forward<Args>(args)...);
        error(msg);
    }

    template <typename... Args>
    static void warn_fmt(fmt::format_string<Args...> format_str, Args&&... args)
    {
        std::string msg = fmt::format(format_str, std::forward<Args>(args)...);
        warn(msg);
    }

    template <typename... Args>
    static void success_fmt(fmt::format_string<Args...> format_str, Args&&... args)
    {
        std::string msg = fmt::format(format_str, std::forward<Args>(args)...);
        success(msg);
    }

    template <typename... Args>
    static void info_fmt(fmt::format_string<Args...> format_str, Args&&... args)
    {
        std::string msg = fmt::format(format_str, std::forward<Args>(args)...);
        info(msg);
    }

    template <typename... Args>
    static void citation_fmt(fmt::format_string<Args...> format_str, Args&&... args)
    {
        std::string msg = fmt::format(format_str, std::forward<Args>(args)...);
        citation(msg);
    }

    // Template-based param function for any type
    template <typename T>
    static void param_value(const std::string& key, const T& value)
    {
        if constexpr (std::is_same_v<T, std::string>) {
            param(key, value);
        } else if constexpr (std::is_integral_v<T>) {
            param(key, static_cast<int>(value));
        } else if constexpr (std::is_floating_point_v<T>) {
            param(key, static_cast<double>(value));
        } else if constexpr (std::is_same_v<T, bool>) {
            param(key, value);
        } else {
            // Fallback for other types using fmt formatting
            std::string value_str = fmt::format("{}", value);
            param(key, value_str);
        }
    }

private:
    // Helper functions for internal formatting
    static void log_colored(fmt::color color, const std::string& prefix, const std::string& msg, bool force_visible = false);
    static void log_plain(const std::string& msg);
    static std::string format_json_value(const json& value);
    static std::string format_energy_relative(double eh_diff);
    static std::string format_time(double fs);
    static std::string get_elapsed_time();
};