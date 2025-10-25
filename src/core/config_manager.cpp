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
#include "global.h"  // For MergeJson
#include <algorithm>
#include <cctype>
#include <stdexcept>

// Claude Generated 2025: Modern configuration manager for Curcuma capabilities

// Single-Module Constructor
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

// Multi-Module Constructor - Claude Generated 2025
ConfigManager::ConfigManager(const std::vector<std::string>& modules, const json& user_input)
    : m_module(modules.empty() ? "" : modules[0]), m_modules(modules)
{
    if (modules.empty()) {
        throw std::runtime_error("ConfigManager: Multi-Module constructor requires at least one module");
    }

    auto& registry = ParameterRegistry::getInstance();

    #ifdef DEBUG_CONFIG_MANAGER
    std::cerr << "[ConfigManager] Multi-Module mode: " << modules.size() << " modules" << std::endl;
    for (const auto& mod : modules) {
        std::cerr << "  - " << mod << std::endl;
    }
    #endif

    // Load and merge configuration for each module
    for (const auto& module : modules) {
        // Get defaults from registry
        json module_defaults = registry.getDefaultJson(module);

        #ifdef DEBUG_CONFIG_MANAGER
        std::cerr << "[ConfigManager] Module '" << module << "': " << module_defaults.size() << " default params" << std::endl;
        #endif

        // Extract user input for this module - Claude Generated 2025: Fixed nested/flat separation
        json module_user_input;

        // 1. Check for nested module config: controller["rmsd"] = {...}
        if (user_input.contains(module) && user_input[module].is_object()) {
            module_user_input = user_input[module];
            #ifdef DEBUG_CONFIG_MANAGER
            std::cerr << "[ConfigManager] Found nested config for '" << module << "'" << std::endl;
            #endif

            // CRITICAL FIX: If this is the PRIMARY module, filter out nested submodule configs
            // Example: controller["confscan"]["rmsd"] = {...} should be skipped for confscan
            if (module == modules[0]) {
                json filtered = json::object();
                for (const auto& item : module_user_input.items()) {
                    // Skip nested objects that are submodule configs
                    if (item.value().is_object()) {
                        bool is_submodule = std::find(m_modules.begin(), m_modules.end(), item.key()) != m_modules.end();
                        if (is_submodule) {
                            #ifdef DEBUG_CONFIG_MANAGER
                            std::cerr << "[ConfigManager] Filtering nested submodule: " << item.key() << std::endl;
                            #endif
                            continue;  // Skip rmsd: {...} for primary module
                        }
                    }
                    filtered[item.key()] = item.value();
                }
                module_user_input = filtered;
            }

        } else if (module == modules[0]) {
            // Primary module: Extract ONLY flat params, SKIP nested submodule configs
            // This prevents conflict like: defaults has rmsd=0.9, CLI has rmsd={method:"subspace"}
            module_user_input = json::object();

            for (const auto& item : user_input.items()) {
                // Skip nested objects that are submodule configs (e.g., controller["rmsd"] = {...})
                if (item.value().is_object()) {
                    // Check if this key is a known submodule
                    bool is_submodule = std::find(m_modules.begin(), m_modules.end(), item.key()) != m_modules.end();
                    if (is_submodule) {
                        #ifdef DEBUG_CONFIG_MANAGER
                        std::cerr << "[ConfigManager] Skipping nested submodule config: " << item.key() << std::endl;
                        #endif
                        continue;  // Skip controller["rmsd"] = {...} for primary module
                    }
                }

                // Include all non-submodule params for primary module
                module_user_input[item.key()] = item.value();
            }

            #ifdef DEBUG_CONFIG_MANAGER
            std::cerr << "[ConfigManager] Extracted " << module_user_input.size()
                      << " flat params for primary module '" << module << "'" << std::endl;
            #endif
        } else {
            // CRITICAL FIX: Secondary module - check if nested under primary module
            // Example: controller["confscan"]["rmsd"] = {...} for rmsd module
            std::string primary = modules[0];
            if (user_input.contains(primary) && user_input[primary].is_object()) {
                const json& primary_config = user_input[primary];
                if (primary_config.contains(module) && primary_config[module].is_object()) {
                    module_user_input = primary_config[module];
                    #ifdef DEBUG_CONFIG_MANAGER
                    std::cerr << "[ConfigManager] Found nested config for '" << module
                              << "' under primary '" << primary << "'" << std::endl;
                    #endif
                }
            }
        }

        // Merge defaults + user input with alias resolution
        json merged_config = module_defaults;

        for (const auto& item : module_user_input.items()) {
            std::string user_key = item.key();
            std::string user_key_lower = user_key;
            std::transform(user_key_lower.begin(), user_key_lower.end(),
                          user_key_lower.begin(), ::tolower);

            // Try alias resolution
            std::string resolved_key = registry.resolveAlias(module, user_key);

            if (!resolved_key.empty() && merged_config.contains(resolved_key)) {
                merged_config[resolved_key] = item.value();
                #ifdef DEBUG_CONFIG_MANAGER
                std::cerr << "[ConfigManager] Module '" << module << "': Alias " << user_key << " -> " << resolved_key << std::endl;
                #endif
                continue;
            }

            // Case-insensitive match
            bool found = false;
            for (const auto& def_item : merged_config.items()) {
                std::string def_key = def_item.key();
                std::string def_key_lower = def_key;
                std::transform(def_key_lower.begin(), def_key_lower.end(),
                              def_key_lower.begin(), ::tolower);

                if (user_key_lower == def_key_lower) {
                    merged_config[def_key] = item.value();
                    found = true;
                    #ifdef DEBUG_CONFIG_MANAGER
                    std::cerr << "[ConfigManager] Module '" << module << "': Merged " << user_key << " -> " << def_key << std::endl;
                    #endif
                    break;
                }
            }

            // Add as new if not found
            if (!found) {
                merged_config[user_key] = item.value();
                #ifdef DEBUG_CONFIG_MANAGER
                std::cerr << "[ConfigManager] Module '" << module << "': Added new key " << user_key << std::endl;
                #endif
            }
        }

        // Cross-module alias resolution - Claude Generated 2025
        // Check if any params in merged_config are actually aliases for OTHER modules
        std::vector<std::string> keys_to_remove;

        for (const auto& item : merged_config.items()) {
            std::string key = item.key();

            // CRITICAL: Check if this key is a canonical parameter in CURRENT module
            // Example: "cycles" is a parameter in confscan, but also an alias in rmsd
            // We must NOT treat it as a cross-module alias if it's OUR parameter!
            const ParameterDefinition* def = registry.findDefinition(module, key);
            if (def != nullptr && def->name == key) {
                // This is a canonical parameter in current module - keep it!
                #ifdef DEBUG_CONFIG_MANAGER
                std::cerr << "[ConfigManager] Keeping canonical parameter: " << key << " in " << module << std::endl;
                #endif
                continue;
            }

            // Check all OTHER modules to see if this key is an alias there
            for (const auto& other_module : m_modules) {
                if (other_module == module) continue;  // Skip same module

                std::string resolved = registry.resolveAlias(other_module, key);
                // Only treat as cross-module alias if resolved name DIFFERS from input
                // (canonical names map to themselves, so we need to filter those out)
                if (!resolved.empty() && resolved != key) {
                    // This key is an alias for a parameter in a DIFFERENT module!
                    // E.g., "rmsdmethod" in confscan resolves to "method" in rmsd

                    #ifdef DEBUG_CONFIG_MANAGER
                    std::cerr << "[ConfigManager] Cross-module alias: " << key
                              << " (" << module << ") -> " << resolved
                              << " (" << other_module << ")" << std::endl;
                    #endif

                    // Store in cross-module params for later application
                    if (m_cross_module_params.find(other_module) == m_cross_module_params.end()) {
                        m_cross_module_params[other_module] = json::object();
                    }
                    m_cross_module_params[other_module][resolved] = item.value();

                    // Mark for removal from current module
                    keys_to_remove.push_back(key);
                    break;  // Found the target module, stop searching
                }
            }
        }

        // Remove cross-module params from current module config
        for (const auto& key : keys_to_remove) {
            merged_config.erase(key);
        }

        // Apply cross-module params if this module is a target
        if (m_cross_module_params.find(module) != m_cross_module_params.end()) {
            for (const auto& item : m_cross_module_params[module].items()) {
                merged_config[item.key()] = item.value();
                #ifdef DEBUG_CONFIG_MANAGER
                std::cerr << "[ConfigManager] Applied cross-module param to '" << module
                          << "': " << item.key() << std::endl;
                #endif
            }
        }

        // Store merged config for this module
        m_module_configs[module] = merged_config;

        #ifdef DEBUG_CONFIG_MANAGER
        std::cerr << "[ConfigManager] Module '" << module << "' final: " << merged_config.size() << " params" << std::endl;
        #endif
    }

    // For backward compatibility: primary module config also in m_config
    if (!m_modules.empty() && m_module_configs.find(m_modules[0]) != m_module_configs.end()) {
        m_config = m_module_configs[m_modules[0]];
    }
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

    // Second try: Check if key uses dot notation (e.g., "topological.save_image" or "rmsd.method")
    size_t dot_pos = key.find('.');
    if (dot_pos != std::string::npos) {
        std::string prefix = key.substr(0, dot_pos);
        std::string suffix = key.substr(dot_pos + 1);

        // Claude Generated 2025: Multi-Module dot notation support
        // If we're in multi-module mode and prefix is a known submodule, look there
        if (!m_modules.empty()) {
            std::string prefix_lower = prefix;
            std::transform(prefix_lower.begin(), prefix_lower.end(),
                          prefix_lower.begin(), ::tolower);

            // Check if prefix matches any of the loaded modules (case-insensitive)
            for (const auto& module : m_modules) {
                std::string module_lower = module;
                std::transform(module_lower.begin(), module_lower.end(),
                              module_lower.begin(), ::tolower);

                if (prefix_lower == module_lower) {
                    // Found matching module! Look in its config
                    auto it = m_module_configs.find(module);
                    if (it != m_module_configs.end()) {
                        const json& module_config = it->second;

                        // Try to find suffix in this module's config
                        std::string suffix_lower = suffix;
                        std::transform(suffix_lower.begin(), suffix_lower.end(),
                                      suffix_lower.begin(), ::tolower);

                        for (const auto& config_item : module_config.items()) {
                            std::string config_key = config_item.key();
                            std::string config_key_lower = config_key;
                            std::transform(config_key_lower.begin(), config_key_lower.end(),
                                          config_key_lower.begin(), ::tolower);

                            if (suffix_lower == config_key_lower) {
                                #ifdef DEBUG_CONFIG_MANAGER
                                std::cerr << "[ConfigManager] Found dot notation key '" << key
                                          << "' in module '" << module << "'" << std::endl;
                                #endif
                                return config_item.value();
                            }
                        }
                    }
                }
            }
        }

        // Fallback: Try to resolve dot notation to flat notation (e.g., "topological.save_image" → "topological_save_image")
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

        // Third try: Check nested JSON structure in primary module
        // e.g., "topological.save_image" → m_config["topological"]["save_image"]
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

// Export configuration for a specific module - Claude Generated 2025
json ConfigManager::exportModule(const std::string& module) const
{
    // Multi-Module mode
    if (!m_modules.empty()) {
        auto it = m_module_configs.find(module);
        if (it != m_module_configs.end()) {
            #ifdef DEBUG_CONFIG_MANAGER
            std::cerr << "[ConfigManager] Exporting module '" << module << "': " << it->second.size() << " params" << std::endl;
            #endif
            return it->second;
        }

        // Module not found
        #ifdef DEBUG_CONFIG_MANAGER
        std::cerr << "[ConfigManager] Module '" << module << "' not found in loaded modules" << std::endl;
        #endif
        return json{};
    }

    // Single-Module mode
    if (module == m_module) {
        return m_config;
    }

    return json{};
}
