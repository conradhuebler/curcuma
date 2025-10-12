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
#include <vector>
#include <map>
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
    /*! \brief Single-Module Constructor - loads defaults and merges with user input
     *
     * @param module Module name (e.g., "analysis", "simplemd")
     * @param user_input User-provided configuration (typically controller["module"])
     */
    ConfigManager(const std::string& module, const json& user_input);

    /*! \brief Multi-Module Constructor - loads defaults from multiple modules - Claude Generated 2025
     *
     * Enables Module Composition pattern: capabilities can import parameters from submodules
     * without parameter duplication. Example: ConfScan imports RMSD parameters.
     *
     * @param modules List of modules to load (first = primary module)
     * @param user_input User-provided configuration with nested module configs
     *
     * Example:
     * ```cpp
     * ConfigManager config({"confscan", "rmsd"}, controller);
     * int confscan_threads = config.get<int>("threads");        // from confscan
     * int rmsd_threads = config.get<int>("rmsd.threads");       // from rmsd module
     * json rmsd_cfg = config.exportModule("rmsd");             // export for RMSDDriver
     * ```
     */
    ConfigManager(const std::vector<std::string>& modules, const json& user_input);

    /*! \brief Type-safe parameter access with exception on missing key
     *
     * @tparam T Parameter type (int, double, std::string, bool)
     * @param key Parameter name (supports dot notation: "topological.save_image" or "rmsd.method")
     * @return Parameter value
     * @throws std::runtime_error if parameter not found
     *
     * Example:
     * ```cpp
     * int max_iter = config.get<int>("max_iterations");           // primary module
     * bool flag = config.get<bool>("topological.exclude_hydrogen");
     * std::string method = config.get<std::string>("rmsd.method");  // from rmsd module
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

    /*! \brief Get primary module name
     *
     * @return Primary module name (e.g., "analysis" or first in multi-module list)
     */
    std::string getModule() const { return m_module; }

    /*! \brief Export configuration for a specific submodule - Claude Generated 2025
     *
     * Returns complete merged configuration (defaults + user overrides) for a submodule,
     * ready to be passed to submodule constructors (e.g., RMSDDriver, PersistentDiagram).
     *
     * @param module Module name to export
     * @return Complete JSON configuration for that module
     *
     * Example:
     * ```cpp
     * json rmsd_config = config.exportModule("rmsd");
     * RMSDDriver driver(rmsd_config);  // Pass complete RMSD config
     * ```
     */
    json exportModule(const std::string& module) const;

private:
    std::string m_module;  //!< Primary module name
    json m_config;         //!< Merged configuration (defaults + user_input) for single-module mode

    // Multi-Module support - Claude Generated 2025
    std::vector<std::string> m_modules;            //!< List of loaded modules (empty in single-module mode)
    std::map<std::string, json> m_module_configs;  //!< Per-module configurations (multi-module mode)
    std::map<std::string, json> m_cross_module_params;  //!< Cross-module alias forwarding (e.g., rmsdmethod → rmsd.method)

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

    /*! \brief Module-specific parameter access - Claude Generated 2025 (PRIVATE to avoid overload ambiguity)
     *
     * Explicitly retrieves parameter from a specific module (for Multi-Module ConfigManager).
     * PRIVATE to prevent ambiguity with get<T>(key, default_value).
     *
     * @tparam T Parameter type
     * @param module Module name (must be in loaded modules list)
     * @param key Parameter name within that module
     * @return Parameter value
     * @throws std::runtime_error if module not loaded or parameter not found
     */
    template <typename T>
    T getFromModule(const std::string& module, const std::string& key) const;
};

// Template implementations must be in header for most compilers

template <typename T>
T ConfigManager::get(const std::string& key) const
{
    // Multi-Module mode: Check for dot notation to identify module
    if (!m_modules.empty()) {
        size_t dot_pos = key.find('.');
        if (dot_pos != std::string::npos) {
            // Dot notation: "rmsd.method" → module="rmsd", param="method"
            std::string module = key.substr(0, dot_pos);
            std::string param = key.substr(dot_pos + 1);

            // Check if this is a loaded module
            if (std::find(m_modules.begin(), m_modules.end(), module) != m_modules.end()) {
                return getFromModule<T>(module, param);  // Use private module-specific getter
            }
        }

        // No dot or unknown module → search in primary module (first in list)
        if (!m_modules.empty()) {
            return getFromModule<T>(m_modules[0], key);
        }
    }

    // Single-Module mode (backward compatible)
    try {
        json value = findKey(key);
        return value.get<T>();
    } catch (...) {
        throw std::runtime_error("ConfigManager: Parameter '" + key +
                                "' not found in module '" + m_module + "'");
    }
}

template <typename T>
T ConfigManager::getFromModule(const std::string& module, const std::string& key) const
{
    // Module-specific getter for Multi-Module mode (PRIVATE to avoid overload ambiguity)
    auto it = m_module_configs.find(module);
    if (it == m_module_configs.end()) {
        throw std::runtime_error("ConfigManager: Module '" + module + "' not loaded");
    }

    const json& module_config = it->second;

    // Case-insensitive search within module config
    std::string key_lower = key;
    std::transform(key_lower.begin(), key_lower.end(), key_lower.begin(), ::tolower);

    for (const auto& item : module_config.items()) {
        std::string config_key = item.key();
        std::string config_key_lower = config_key;
        std::transform(config_key_lower.begin(), config_key_lower.end(),
                      config_key_lower.begin(), ::tolower);

        if (key_lower == config_key_lower) {
            return item.value().get<T>();
        }
    }

    throw std::runtime_error("ConfigManager: Parameter '" + key +
                            "' not found in module '" + module + "'");
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
