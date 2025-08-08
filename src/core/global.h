/*
 * <Some globale definition for chemical structures.>
 * Copyright (C) 2019 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <iostream>
#include <map>
#include <set>

#include <Eigen/Dense>

#include "src/global_config.h"
#include "src/version.h"

#include "json.hpp"

// for convenience
using json = nlohmann::json;

const double pi = 3.14159265359;
const double au = 0.52917721092; // Angstrom
const double amu2au = 1822.8884850;
const double kb_Eh = 3.166811e-6; // Hartree
const double kb_SI = 1.380649e-23; // SI
const double kb_eV = 8.617333262e-5; // eV
const double fs2amu = 41.34137314;
const double R = 8.31446261815324;
const double atomic_mass = 1.66053906660e-27;
const double T_Eh = 3.1577464e5;
const double eV2Eh = 27.211386245988; // Hartree

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Geometry;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Vector3d Position;
typedef Eigen::Vector4d Vector4d;

typedef Eigen::VectorXd Vector;
typedef std::pair<int, int> IntPair;
typedef std::vector<std::string> StringList;

struct Mol {
    double m_energy;
    double m_spin;

    int m_number_atoms;
    int m_charge;

    std::string m_commentline;
    std::string m_formula; // Claude Generated: Molecular formula

    Geometry m_geometry;
    Vector m_partial_charges;
    
    std::vector<std::pair<int, int>> m_bonds;
    std::vector<int> m_atoms;
};

inline Vector PositionPair2Vector(const std::pair<Position, Position>& pair)
{
    Vector vector = Vector::Zero(6);
    vector(0) = pair.first(0);
    vector(1) = pair.first(1);
    vector(2) = pair.first(2);
    vector(3) = pair.second(0);
    vector(4) = pair.second(1);
    vector(5) = pair.second(2);
    return vector;
}

inline Vector PositionPair2Vector(const Position& first, const Position& second)
{
    Vector vector = Vector::Zero(6);
    vector(0) = first(0);
    vector(1) = first(1);
    vector(2) = first(2);
    vector(3) = second(0);
    vector(4) = second(1);
    vector(5) = second(2);
    return vector;
}

class LimitedStorage {
public:
    inline LimitedStorage(unsigned int size)
        : m_size(size)
    {
    }
    inline ~LimitedStorage()
    {
        m_shelf.clear();
    }
    inline void addItem(double rmsd, const std::vector<int>& vector)
    {
        m_shelf.insert(std::pair<double, std::vector<int>>(rmsd, vector));
        if (m_shelf.size() >= m_size)
            m_shelf.erase(--m_shelf.end());
    }

    inline const std::map<double, std::vector<int>>* data() const { return &m_shelf; }
    inline int size() const { return data()->size(); }

private:
    unsigned int m_size;
    std::map<double, std::vector<int>> m_shelf;
};

inline std::pair<Position, Position> Vector2PositionPair(const Vector& vector)
{
    return std::pair<Position, Position>(Position{ vector(0), vector(1), vector(2) }, Position{ vector(3), vector(4), vector(5) });
}

inline int CompareTopoMatrix(const Matrix& m1, const Matrix& m2)
{
    if (m1.rows() != m2.rows() || m1.cols() != m2.cols() || m1.cols() != m1.rows())
        return -1;
    int result = 0;
    for (int i = 0; i < m1.rows(); ++i)
        for (int j = i + 1; j < m1.cols(); ++j)
            result += (m1(i, j) != m2(i, j)) + (m1(j, i) != m2(j, i));

    return result;
}

inline void CompactTopo(const Matrix& m1)
{
    for (int i = 0; i < m1.rows(); ++i)
        for (int j = i + 1; j < m1.cols(); ++j) {
            if (m1(i, j) == 1)
                std::cout << "  " << i << "   ... " << j << std::endl;
        }
}
/*
inline json CLI2Json(int argc, char** argv)
{
    json controller;
    json key;
    if (argc < 2)
        return controller;
    std::string keyword = argv[1];
    keyword.erase(0, 1);
    for (int i = 2; i < argc; ++i) {
        std::string current = argv[i];
        std::string sub = current.substr(0, 1);
        if ((i + 1) >= argc) {
            current.erase(0, 1);
            if (sub.compare("-") == 0)
                key[current] = true;
            else
                key[current] = false;
        } else {
            if (sub.compare("-") == 0 && ((i + 1) < argc)) {
                std::string next = argv[i + 1];
                std::string next_sub = next.substr(0, 1);
                bool isNumber = true;
                bool isVector = false;
                double number = 0.0;
                // std::size_t found = next.find("|");
                if (next.find("|") != std::string::npos || next.find(",") != std::string::npos || next.find(":") != std::string::npos) {
                    isNumber = false;
                    isVector = true;
                } else {
                    try {
                        number = std::stod(next);
                    } catch (const std::invalid_argument& error) {
                        isNumber = false;
                    }
                }
                if (isNumber) {
                    current.erase(0, 1);
                    key[current] = number;
                } else {
                    if (next_sub.compare("-") == 0 || next.compare("false") == 0) {
                        current.erase(0, 1);
                        key[current] = false;
                        continue;
                    } else if (next_sub.compare("+") == 0 || next.compare("true") == 0) {
                        current.erase(0, 1);
                        key[current] = true;
                        continue;
                    } else {
                        current.erase(0, 1);
                        if (isVector) {
                            key[current] = argv[i + 1];
                        } else {
                            //try {
                            //    key[current] = std::stoi(argv[i + 1]);
                            //} catch (const std::invalid_argument& error) {
                            try {
                                key[current] = std::stod(argv[i + 1]);
                            } catch (const std::invalid_argument& error) {
                                key[current] = argv[i + 1];
                            }
                        }
                        //}

                        ++i;
                    }
                }
            }
        }
    }
    controller[keyword] = key;
    return controller;
}
*/
/* this is the 2nd github copilot version */
/*
inline json CLI2Json(int argc, char** argv)
{
    json controller;
    json key;
    if (argc < 2)
        return controller;

    std::string keyword = argv[1];
    keyword.erase(0, 1);

    for (int i = 2; i < argc; ++i) {
        std::string current = argv[i];
        std::string sub = current.substr(0, 1);

        if (sub == "-") {
            current.erase(0, 1);
            if ((i + 1) >= argc || argv[i + 1][0] == '-' || argv[i + 1] == std::string("true") || argv[i + 1] == std::string("+")) {
                key[current] = true;
            } else if (argv[i + 1] == std::string("false")) {
                key[current] = false;
                ++i;
            } else {
                std::string next = argv[i + 1];
                bool isNumber = true;
                bool isVector = next.find("|") != std::string::npos || next.find(",") != std::string::npos || next.find(":") != std::string::npos;
                bool isRange = next.find("-") != std::string::npos;

                if (!isVector && !isRange) {
                    try {
                        std::stod(next);
                    } catch (const std::invalid_argument&) {
                        isNumber = false;
                    }
                }

                if (isNumber) {
                    key[current] = std::stod(next);
                } else if (isVector || isRange) {
                    key[current] = next;
                } else {
                    key[current] = next;
                }
                ++i;
            }
        }
    }

    controller[keyword] = key;
    return controller;
}*/
/* this is the github copilot version */

template <class T>
inline T Json2KeyWord(const json& controller, std::string name)
{
    T temp;
    bool found = false;
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    for (auto& el : controller.items()) {
        std::string key = el.key();
        transform(key.begin(), key.end(), key.begin(), ::tolower);
        if (key.compare(name) == 0) {
            temp = el.value();
            found = true;
        }
    }
    if (found)
        return temp;
    else
        throw -1;
}

inline json MergeJson(const json& reference, const json& patch)
{
    json result = reference;
    for (const auto& object : patch.items()) {
        bool found = false;
        std::string outer = object.key();
        transform(outer.begin(), outer.end(), outer.begin(), ::tolower);
        for (const auto& local : reference.items()) {
            std::string inner = local.key();
            transform(inner.begin(), inner.end(), inner.begin(), ::tolower);
            if (outer.compare(inner) == 0) {
                result[local.key()] = object.value();
                found = true;
            }
        }
        if (!found) {
            result[outer] = object.value();
        }
    }
    return result;
}

inline json IncludeJson(const json& reference, const json& patch)
{
    json result = reference;
    for (const auto& object : patch.items()) {
        result[object.key()] = object.value();
    }
    return result;
}

inline void PrintController(const json& controller)
{
    for (const auto& entry : controller.items())
        std::cout << "**| " << entry.key() << ": " << entry.value() << " |**      ";
    std::cout << std::endl
              << std::endl;
}

inline int MaxThreads()
{
    int threads = 1;
    const char* val = std::getenv("CurcumaThreads");
    if (val == nullptr) { // invalid to assign nullptr to std::string
    } else {
        try {
            threads = atoi(std::getenv("CurcumaThreads"));
        } catch (const std::invalid_argument& error) {
        }
    }
    return threads;
}

// Unit Conversion Functions - Claude Generated
inline double hartree_to_kjmol(double eh) { return eh * CURCUMA_EH_TO_KJMOL; }
inline double hartree_to_kcalmol(double eh) { return eh * CURCUMA_EH_TO_KCALMOL; }
inline double hartree_to_ev(double eh) { return eh * CURCUMA_EH_TO_EV; }
inline double hartree_to_wavenumber(double eh) { return eh * CURCUMA_EH_TO_WAVENUMBER; }
inline double bohr_to_angstrom(double bohr) { return bohr * CURCUMA_BOHR_TO_ANGSTROM; }
inline double angstrom_to_bohr(double ang) { return ang / CURCUMA_BOHR_TO_ANGSTROM; }

#include <chrono>
#include <fmt/color.h>
#include <fmt/format.h>
#include <unistd.h>

/**
 * @brief Curcuma Logging System - Claude Generated
 *
 * Provides unified logging with verbosity levels, colors, and unit-aware formatting
 * according to curcuma development standards.
 */
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

    // Initialize logger with environment detection
    static void initialize()
    {
        m_verbosity = 1; // Default: Normal Print

        // Auto-detect color support
        m_use_colors = isatty(STDOUT_FILENO) && (std::getenv("CURCUMA_NO_COLOR") == nullptr);

#ifdef CURCUMA_COLOR_DEFAULT
        // Build-time color preference
        if (std::getenv("CURCUMA_NO_COLOR") == nullptr) {
            m_use_colors = true;
        }
#endif

        m_format = OutputFormat::TERMINAL;
        m_start_time = std::chrono::high_resolution_clock::now();
    }

    // Core logging functions with verbosity control
    static void error(const std::string& msg)
    {
        // Errors always visible
        log_colored(fmt::color::red, "[ERROR] ", msg, true);
    }

    static void warn(const std::string& msg)
    {
        if (m_verbosity >= 1) {
            log_colored(fmt::color::orange, "[WARN]  ", msg);
        }
    }

    static void success(const std::string& msg)
    {
        if (m_verbosity >= 1) {
            log_colored(fmt::color::lime_green, "[OK]    ", msg);
        }
    }

    static void info(const std::string& msg)
    {
        if (m_verbosity >= 2) {
            log_plain("        " + msg);
        }
    }

    static void citation(const std::string& ref)
    {
        if (m_verbosity >= 2) {
            log_colored(fmt::color::green, "[CITE]  ", ref);
        }
    }

    static void param(const std::string& key, const std::string& value)
    {
        if (m_verbosity >= 2) {
            log_colored(fmt::color::cornflower_blue, "[PARAM] ", key + ": " + value);
        }
    }

    // Overloaded version for integers - Claude Generated
    static void param(const std::string& key, int value)
    {
        if (m_verbosity >= 2) {
            log_colored(fmt::color::cornflower_blue, "[PARAM] ", key + ": " + std::to_string(value));
        }
    }

    // Overloaded version for doubles - Claude Generated
    static void param(const std::string& key, double value)
    {
        if (m_verbosity >= 2) {
            log_colored(fmt::color::cornflower_blue, "[PARAM] ", fmt::format("{}: {:.6g}", key, value));
        }
    }

    // Overloaded version for booleans - Claude Generated
    static void param(const std::string& key, bool value)
    {
        if (m_verbosity >= 2) {
            log_colored(fmt::color::cornflower_blue, "[PARAM] ", key + ": " + (value ? "true" : "false"));
        }
    }

    // JSON parameter table printing - Claude Generated
    static void param_table(const json& parameters, const std::string& title = "Parameters")
    {
        if (m_verbosity >= 1) {
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

    // Compact 3-column parameter table showing parameter sets side by side - Claude Generated
    static void param_comparison_table(const json& defaults, const json& controller, const std::string& title = "Parameter Settings")
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

            // Separate parameters into short and long value lists
            std::vector<std::string> short_params, long_params;
            std::vector<bool> param_changed;

            for (const auto& key : sorted_keys) {
                std::string display_value;
                bool is_changed = false;

                // Determine current value and check if changed
                if (controller.contains(key)) {
                    const auto& ctrl_val = controller[key];
                    display_value = format_json_value(ctrl_val);

                    // Check if different from default - better comparison logic
                    if (defaults.contains(key)) {
                        try {
                            // Compare JSON values directly for more accurate detection
                            if (defaults[key] != ctrl_val) {
                                is_changed = true;
                            }
                        } catch (...) {
                            // Fallback to string comparison if JSON comparison fails
                            std::string default_str = format_json_value(defaults[key]);
                            if (default_str != display_value) {
                                is_changed = true;
                            }
                        }
                    } else {
                        // New parameter not in defaults
                        is_changed = true;
                    }
                } else if (defaults.contains(key)) {
                    display_value = format_json_value(defaults[key]);
                } else {
                    display_value = "N/A";
                }

                // Separate based on value length
                if (display_value.length() > 25) {
                    long_params.push_back(key);
                    param_changed.push_back(is_changed);
                } else {
                    short_params.push_back(key);
                    param_changed.push_back(is_changed);
                }
            }

            // Display short parameters in 3-column format
            if (!short_params.empty()) {
                size_t total_short = short_params.size();
                size_t rows_needed = (total_short + 2) / 3; // Ceiling division

                // Find maximum parameter name length for alignment
                size_t max_param_length = 0;
                for (const auto& key : short_params) {
                    max_param_length = std::max(max_param_length, key.length());
                }
                max_param_length = std::min(max_param_length, static_cast<size_t>(15)); // Limit width

                // Column width = parameter_name + " : " + value + potential "*"
                size_t col_width = max_param_length + 15; // Space for " : value*"

                // Print header for 3 columns
                std::string header_format = fmt::format("        {{:<{}}} | {{:<{}}} | {{:<{}}}", col_width, col_width, col_width);
                log_colored(fmt::color::white, "", fmt::format(header_format, "Parameter : Value", "Parameter : Value", "Parameter : Value"));

                std::string header_separator = fmt::format("        {0:-<{1}}-+-{0:-<{1}}-+-{0:-<{1}}", "", col_width);
                log_plain(header_separator);

                // Process parameters in groups of 3 (side by side)
                for (size_t row = 0; row < rows_needed; ++row) {
                    std::string col1_content, col2_content, col3_content;
                    bool row_has_changes = false;

                    // Process each column for this row
                    for (int col = 0; col < 3; ++col) {
                        size_t param_index = row + col * rows_needed;

                        if (param_index < total_short) {
                            const std::string& key = short_params[param_index];
                            std::string display_value;
                            bool is_changed = false;

                            // Get the current value and change status
                            size_t original_index = std::find(sorted_keys.begin(), sorted_keys.end(), key) - sorted_keys.begin();
                            is_changed = param_changed[original_index];

                            if (controller.contains(key)) {
                                display_value = format_json_value(controller[key]);
                            } else if (defaults.contains(key)) {
                                display_value = format_json_value(defaults[key]);
                            } else {
                                display_value = "N/A";
                            }

                            // Truncate long values for display
                            if (display_value.length() > 12) {
                                display_value = display_value.substr(0, 9) + "...";
                            }

                            // Truncate long parameter names
                            std::string display_key = key;
                            if (display_key.length() > max_param_length) {
                                display_key = display_key.substr(0, max_param_length - 2) + "..";
                            }

                            // Format the parameter entry with star if changed
                            std::string param_entry = display_key + " : " + display_value;
                            if (is_changed) {
                                param_entry += "*";
                                row_has_changes = true;
                            }

                            // Assign to appropriate column
                            if (col == 0) {
                                col1_content = param_entry;
                            } else if (col == 1) {
                                col2_content = param_entry;
                            } else {
                                col3_content = param_entry;
                            }
                        }
                    }

                    // Pad columns to consistent width
                    if (col1_content.length() > col_width)
                        col1_content = col1_content.substr(0, col_width);
                    if (col2_content.length() > col_width)
                        col2_content = col2_content.substr(0, col_width);
                    if (col3_content.length() > col_width)
                        col3_content = col3_content.substr(0, col_width);

                    // Format and print the row
                    std::string row_format = fmt::format("        {{:<{}}} | {{:<{}}} | {{:<{}}}", col_width, col_width, col_width);
                    std::string formatted_row = fmt::format(row_format, col1_content, col2_content, col3_content);

                    // Color code based on changes in this row
                    if (row_has_changes) {
                        log_colored(fmt::color::yellow, "", formatted_row);
                    } else {
                        log_colored(fmt::color::cornflower_blue, "", formatted_row);
                    }
                }
                log_plain("");
            }

            // Display long parameters separately with full values
            if (!long_params.empty()) {
                log_colored(fmt::color::white, "        ", "Parameters with long values:");
                std::string long_separator(40, '-');
                log_plain("        " + long_separator);

                for (const auto& key : long_params) {
                    std::string full_value;
                    bool is_changed = false;

                    // Get change status
                    size_t original_index = std::find(sorted_keys.begin(), sorted_keys.end(), key) - sorted_keys.begin();
                    is_changed = param_changed[original_index];

                    if (controller.contains(key)) {
                        full_value = format_json_value(controller[key]);
                    } else if (defaults.contains(key)) {
                        full_value = format_json_value(defaults[key]);
                    } else {
                        full_value = "N/A";
                    }

                    std::string param_line = "        " + key + " : " + full_value;
                    if (is_changed) {
                        param_line += "*";
                        log_colored(fmt::color::yellow, "", param_line);
                    } else {
                        log_colored(fmt::color::cornflower_blue, "", param_line);
                    }
                }
                log_plain("");
            }

            // Add legend
            log_colored(fmt::color::white, "        ", "Legend: * = modified from default");
            log_plain("");
        }
    }

    static void result_raw(const std::string& data)
    {
        // Raw results for scripting - no colors, no prefix
        fmt::print("{}\n", data);
    }

    static void header(const std::string& title)
    {
        if (m_verbosity >= 2) {
            std::string separator(title.length() + 4, '=');
            log_colored(fmt::color::cyan, "", separator);
            log_colored(fmt::color::cyan, "", "  " + title);
            log_colored(fmt::color::cyan, "", separator);
        }
    }

    static void progress(int current, int total, const std::string& msg)
    {
        if (m_verbosity >= 2) {
            double percent = (100.0 * current) / total;
            std::string bar = create_progress_bar(percent, 30);
            log_colored(fmt::color::yellow, "[PROG]  ",
                fmt::format("{} [{:3.0f}%] {}/{}", msg, percent, current, total));
        }
    }

    // Unit-aware formatting functions
    static void energy_rel(double value_eh, const std::string& label)
    {
        if (m_verbosity >= 1) {
            std::string formatted = format_energy_relative(value_eh);
            log_colored(fmt::color::cornflower_blue, "[ENERGY]", label + ": " + formatted);
        }
    }

    static void energy_abs(double value_eh, const std::string& label)
    {
        if (m_verbosity >= 1) {
            log_colored(fmt::color::cornflower_blue, "[ENERGY]",
                fmt::format("{}: {:.8f} Eh", label, value_eh));
        }
    }

    static void length(double value_bohr, const std::string& label)
    {
        if (m_verbosity >= 2) {
            log_colored(fmt::color::cornflower_blue, "[LENGTH]",
                fmt::format("{}: {:.4f} Å", label, bohr_to_angstrom(value_bohr)));
        }
    }

    static void time(double value_aut, const std::string& label)
    {
        if (m_verbosity >= 2) {
            std::string formatted = format_time(value_aut * CURCUMA_AUT_TO_FS);
            log_colored(fmt::color::cornflower_blue, "[TIME]  ", label + ": " + formatted);
        }
    }

#ifdef CURCUMA_DEBUG
    // Debug functions (compile-time conditional)
    static void debug(int level, const std::string& msg)
    {
        if (level <= CURCUMA_DEBUG_LEVEL) {
            log_colored(fmt::color::magenta, fmt::format("[DEBUG{}] ", level), msg);
        }
    }

    template <typename T>
    static void debug_var(const std::string& name, const T& value)
    {
        if (1 <= CURCUMA_DEBUG_LEVEL) {
            log_colored(fmt::color::magenta, "[DEBUG] ", fmt::format("{} = {}", name, value));
        }
    }

    static void debug_timing(const std::string& label)
    {
        if (2 <= CURCUMA_DEBUG_LEVEL) {
            auto now = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - m_start_time);
            log_colored(fmt::color::magenta, "[TIMING]",
                fmt::format("{}: {:.3f} ms", label, duration.count() / 1000.0));
        }
    }
#endif

private:
    // Helper functions for parameter table - Claude Generated
    static std::string format_json_value(const json& value)
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

    static std::string get_json_type_name(const json& value)
    {
        if (value.is_string())
            return "string";
        if (value.is_number_integer())
            return "integer";
        if (value.is_number_float())
            return "float";
        if (value.is_boolean())
            return "boolean";
        if (value.is_array())
            return "array";
        if (value.is_object())
            return "object";
        return "unknown";
    }

private:
    static void log_colored(fmt::color color, const std::string& prefix, const std::string& msg, bool force = false)
    {
        if (m_use_colors && !force) {
            fmt::print(fmt::fg(color), "{}{}\n", prefix, msg);
        } else {
            fmt::print("{}{}\n", prefix, msg);
        }
    }

    static void log_plain(const std::string& msg)
    {
        fmt::print("{}\n", msg);
    }

    static std::string format_energy_relative(double eh_diff)
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

    static std::string format_time(double fs)
    {
        if (fs < 1000.0) {
            return fmt::format("{:.1f} fs", fs);
        } else if (fs < 1e6) {
            return fmt::format("{:.1f} ps", fs / 1000.0);
        } else {
            return fmt::format("{:.1f} ns", fs / 1e6);
        }
    }

    static std::string create_progress_bar(double percent, int width)
    {
        int filled = static_cast<int>(percent * width / 100.0);
        std::string bar = "[";
        for (int i = 0; i < width; ++i) {
            bar += (i < filled) ? "█" : "░";
        }
        bar += "]";
        return bar;
    }
};

// Static member definitions moved to curcuma_logger.cpp - Claude Generated
