/*
 * <Configuration Manager for Curcuma Capabilities>
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include <string>
#include <algorithm>
#include "json.hpp"

using json = nlohmann::json;

/*! \brief Modern configuration manager for Curcuma capabilities - Claude Generated 2025
 *
 * Provides type-safe, hierarchical parameter access with automatic merging of
 * defaults from ParameterRegistry and user-provided configuration.
 *
 * **Key Features:**
 * - Automatic default loading from ParameterRegistry
 * - Type-safe parameter access via get<T>()
 * - Hierarchical dot notation: "topological.save_image"
 * - Fallback to flat notation: "topological_save_image"
 * - Case-insensitive parameter names
 * - Default value support
 *
 * **Example Usage:**
 * ```cpp
 * ConfigManager config("analysis", user_json);
 * std::string format = config.get<std::string>("output_format");
 * bool save = config.get<bool>("topological.save_image");  // or "topological_save_image"
 * int window = config.get<int>("window", 10);  // with default
 * ```
 *
 * **Migration from Json2KeyWord:**
 * ```cpp
 * // Old:
 * std::string format = Json2KeyWord<std::string>(m_config, "output_format");
 *
 * // New:
 * std::string format = config.get<std::string>("output_format");
 * ```
 */
class ConfigManager
{
public:
    /*! \brief Constructor - loads defaults and merges with user input
     *
     * @param module Module name (e.g., "analysis", "simplemd")
     * @param user_input User-provided configuration (typically controller["module"])
     */
    ConfigManager(const std::string& module, const json& user_input);

    /*! \brief Type-safe parameter access with exception on missing key
     *
     * @tparam T Parameter type (int, double, std::string, bool)
     * @param key Parameter name (supports dot notation: "topological.save_image")
     * @return Parameter value
     * @throws std::runtime_error if parameter not found
     *
     * Example:
     * ```cpp
     * int max_iter = config.get<int>("max_iterations");
     * bool flag = config.get<bool>("topological.exclude_hydrogen");
     * ```
     */
    template <typename T>
    T get(const std::string& key) const;

    /*! \brief Type-safe parameter access with default value
     *
     * @tparam T Parameter type
     * @param key Parameter name
     * @param default_value Value to return if key not found
     * @return Parameter value or default
     *
     * Example:
     * ```cpp
     * int window = config.get<int>("window", 10);
     * ```
     */
    template <typename T>
    T get(const std::string& key, T default_value) const;

    /*! \brief Check if parameter exists
     *
     * @param key Parameter name
     * @return true if parameter exists, false otherwise
     */
    bool has(const std::string& key) const;

    /*! \brief Export merged configuration as JSON
     *
     * Useful for debugging or passing to external APIs that expect JSON.
     *
     * @return Complete merged configuration
     */
    json exportConfig() const { return m_config; }

    /*! \brief Get module name
     *
     * @return Module name (e.g., "analysis")
     */
    std::string getModule() const { return m_module; }

private:
    std::string m_module;  //!< Module name
    json m_config;         //!< Merged configuration (defaults + user_input)

    /*! \brief Resolve key with dot notation to flat key
     *
     * Supports both notations:
     * - "topological.save_image" → tries both nested and "topological_save_image"
     * - "topological_save_image" → direct flat access
     *
     * @param key Parameter name
     * @return Resolved flat key
     */
    std::string resolveFlatKey(const std::string& key) const;

    /*! \brief Case-insensitive key lookup
     *
     * Maintains compatibility with Json2KeyWord behavior.
     *
     * @param key Parameter name
     * @return Parameter value as JSON
     * @throws std::runtime_error if not found
     */
    json findKey(const std::string& key) const;
};

// Template implementations must be in header for most compilers

template <typename T>
T ConfigManager::get(const std::string& key) const
{
    try {
        json value = findKey(key);
        return value.get<T>();
    } catch (...) {
        throw std::runtime_error("ConfigManager: Parameter '" + key +
                                "' not found in module '" + m_module + "'");
    }
}

template <typename T>
T ConfigManager::get(const std::string& key, T default_value) const
{
    try {
        return get<T>(key);
    } catch (...) {
        return default_value;
    }
}
