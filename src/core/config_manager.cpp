/*
 * <Configuration Manager Implementation>
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

#include "config_manager.h"
#include "parameter_registry.h"
#include <algorithm>
#include <cctype>
#include <stdexcept>

// Claude Generated 2025: Modern configuration manager for Curcuma capabilities

ConfigManager::ConfigManager(const std::string& module, const json& user_input)
    : m_module(module)
{
    // Load defaults from ParameterRegistry
    m_config = ParameterRegistry::getInstance().getDefaultJson(module);

    // DEBUG: Show what we're merging - Claude Generated
    #ifdef DEBUG_CONFIG_MANAGER
    std::cerr << "[ConfigManager] Module: " << module << std::endl;
    std::cerr << "[ConfigManager] Defaults loaded: " << m_config.size() << " parameters" << std::endl;
    std::cerr << "[ConfigManager] User input: " << user_input.dump(2) << std::endl;
    #endif

    // Merge user input (case-insensitive + alias resolution, like MergeJson)
    auto& registry = ParameterRegistry::getInstance();

    for (const auto& item : user_input.items()) {
        std::string user_key = item.key();
        std::string user_key_lower = user_key;
        std::transform(user_key_lower.begin(), user_key_lower.end(),
                      user_key_lower.begin(), ::tolower);

        // First, try alias resolution (case-insensitive)
        std::string resolved_key = registry.resolveAlias(module, user_key);

        // If alias resolved to a parameter name, use that
        if (!resolved_key.empty() && m_config.contains(resolved_key)) {
            m_config[resolved_key] = item.value();
            #ifdef DEBUG_CONFIG_MANAGER
            std::cerr << "[ConfigManager] Alias resolved: " << user_key << " -> " << resolved_key << " = " << item.value() << std::endl;
            #endif
            continue;
        }

        // Try to find matching key in defaults (case-insensitive)
        bool found = false;
        for (const auto& def_item : m_config.items()) {
            std::string def_key = def_item.key();
            std::string def_key_lower = def_key;
            std::transform(def_key_lower.begin(), def_key_lower.end(),
                          def_key_lower.begin(), ::tolower);

            if (user_key_lower == def_key_lower) {
                // Override default with user value
                m_config[def_key] = item.value();
                found = true;
                #ifdef DEBUG_CONFIG_MANAGER
                std::cerr << "[ConfigManager] Merged: " << user_key << " -> " << def_key << " = " << item.value() << std::endl;
                #endif
                break;
            }
        }

        // If not found in defaults, add as new key
        if (!found) {
            m_config[user_key] = item.value();
            #ifdef DEBUG_CONFIG_MANAGER
            std::cerr << "[ConfigManager] Added new: " << user_key << " = " << item.value() << std::endl;
            #endif
        }
    }

    #ifdef DEBUG_CONFIG_MANAGER
    std::cerr << "[ConfigManager] Final config: " << m_config.size() << " parameters" << std::endl;
    #endif
}

bool ConfigManager::has(const std::string& key) const
{
    try {
        findKey(key);
        return true;
    } catch (...) {
        return false;
    }
}

json ConfigManager::findKey(const std::string& key) const
{
    // First try: Exact case-insensitive match
    std::string key_lower = key;
    std::transform(key_lower.begin(), key_lower.end(), key_lower.begin(), ::tolower);

    for (const auto& item : m_config.items()) {
        std::string config_key = item.key();
        std::string config_key_lower = config_key;
        std::transform(config_key_lower.begin(), config_key_lower.end(),
                      config_key_lower.begin(), ::tolower);

        if (key_lower == config_key_lower) {
            return item.value();
        }
    }

    // Second try: Check if key uses dot notation (e.g., "topological.save_image")
    // Try to resolve to flat notation (e.g., "topological_save_image")
    size_t dot_pos = key.find('.');
    if (dot_pos != std::string::npos) {
        std::string flat_key = key;
        std::replace(flat_key.begin(), flat_key.end(), '.', '_');

        // Try flat notation with case-insensitive search
        std::string flat_key_lower = flat_key;
        std::transform(flat_key_lower.begin(), flat_key_lower.end(),
                      flat_key_lower.begin(), ::tolower);

        for (const auto& item : m_config.items()) {
            std::string config_key = item.key();
            std::string config_key_lower = config_key;
            std::transform(config_key_lower.begin(), config_key_lower.end(),
                          config_key_lower.begin(), ::tolower);

            if (flat_key_lower == config_key_lower) {
                return item.value();
            }
        }

        // Third try: Check nested JSON structure
        // e.g., "topological.save_image" → m_config["topological"]["save_image"]
        std::string prefix = key.substr(0, dot_pos);
        std::string suffix = key.substr(dot_pos + 1);

        std::string prefix_lower = prefix;
        std::transform(prefix_lower.begin(), prefix_lower.end(),
                      prefix_lower.begin(), ::tolower);

        for (const auto& item : m_config.items()) {
            std::string config_key = item.key();
            std::string config_key_lower = config_key;
            std::transform(config_key_lower.begin(), config_key_lower.end(),
                          config_key_lower.begin(), ::tolower);

            if (prefix_lower == config_key_lower && item.value().is_object()) {
                // Found nested object, try to find suffix
                const json& nested = item.value();
                std::string suffix_lower = suffix;
                std::transform(suffix_lower.begin(), suffix_lower.end(),
                              suffix_lower.begin(), ::tolower);

                for (const auto& nested_item : nested.items()) {
                    std::string nested_key = nested_item.key();
                    std::string nested_key_lower = nested_key;
                    std::transform(nested_key_lower.begin(), nested_key_lower.end(),
                                  nested_key_lower.begin(), ::tolower);

                    if (suffix_lower == nested_key_lower) {
                        return nested_item.value();
                    }
                }
            }
        }
    }

    // Not found
    throw std::runtime_error("Key not found: " + key);
}

std::string ConfigManager::resolveFlatKey(const std::string& key) const
{
    // Convert dot notation to underscore notation
    std::string flat_key = key;
    std::replace(flat_key.begin(), flat_key.end(), '.', '_');
    return flat_key;
}
