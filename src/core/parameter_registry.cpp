#include "parameter_registry.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include "json.hpp"

using json = nlohmann::json;

ParameterRegistry& ParameterRegistry::getInstance()
{
    static ParameterRegistry instance;
    return instance;
}

void ParameterRegistry::addDefinition(const std::string& module, ParameterDefinition&& def)
{
    std::string canonical_name = def.name;
    registry[module].push_back(std::move(def));

    // Map the canonical name to itself
    alias_to_name_map[module][canonical_name] = canonical_name;
    // Map all aliases to the canonical name
    const auto& added_def = registry[module].back();
    for (const auto& alias : added_def.aliases) {
        alias_to_name_map[module][alias] = canonical_name;
    }

    // Claude Generated: Populate inverse lookup (lowercased name/alias -> owning modules).
    // Used by findOwnerModules to auto-route flat CLI flags.
    auto add_owner = [&](const std::string& key) {
        std::string lower = key;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        // Claude Generated (Jun 2026): canonicalize hyphens to underscores
        lower.erase(std::remove(lower.begin(), lower.end(), '-'), lower.end());
        auto& mods = name_to_modules_map[lower];
        if (std::find(mods.begin(), mods.end(), module) == mods.end()) {
            mods.push_back(module);
        }
    };
    add_owner(canonical_name);
    for (const auto& alias : added_def.aliases) {
        add_owner(alias);
    }
}

const ParameterDefinition* ParameterRegistry::findDefinition(const std::string& module, const std::string& alias) const
{
    auto module_it = alias_to_name_map.find(module);
    if (module_it == alias_to_name_map.end()) {
        return nullptr;
    }

    auto alias_it = module_it->second.find(alias);
    if (alias_it == module_it->second.end()) {
        return nullptr;
    }

    const std::string& canonical_name = alias_it->second;

    auto registry_it = registry.find(module);
    if (registry_it != registry.end()) {
        for (const auto& def : registry_it->second) {
            if (def.name == canonical_name) {
                return &def;
            }
        }
    }

    return nullptr;
}

std::vector<ParameterDefinition> ParameterRegistry::getForModule(const std::string& module) const
{
    auto it = registry.find(module);
    if (it != registry.end()) {
        return it->second;
    }
    return {};
}

// Claude Generated: Print help for a specific module
void ParameterRegistry::printHelp(const std::string& module) const
{
    auto it = registry.find(module);
    if (it == registry.end()) {
        std::cout << "No parameters registered for module: " << module << std::endl;
        return;
    }

    std::cout << "Parameters for module: " << module << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    // Group parameters by category
    std::map<std::string, std::vector<const ParameterDefinition*>> by_category;
    for (const auto& param : it->second) {
        by_category[param.category].push_back(&param);
    }

    // Print each category
    for (const auto& [category, params] : by_category) {
        std::cout << "\n[" << category << "]" << std::endl;

        for (const auto* param : params) {
            std::cout << "  -" << param->name;

            // Show type
            std::cout << " <";
            switch (param->type) {
            case ParamType::String:
                std::cout << "string";
                break;
            case ParamType::Int:
                std::cout << "int";
                break;
            case ParamType::Double:
                std::cout << "double";
                break;
            case ParamType::Bool:
                std::cout << "bool";
                break;
            case ParamType::StringList:
                std::cout << "list";
                break;
            case ParamType::Json:
                std::cout << "json";
                break;
            case ParamType::Selection:
                std::cout << "selection";
                break;
            case ParamType::Path:
                std::cout << "path";
                break;
            }
            std::cout << ">";

            // Show default value
            std::cout << " (default: ";
            try {
                switch (param->type) {
                case ParamType::String:
                    std::cout << std::any_cast<std::string>(param->defaultValue);
                    break;
                case ParamType::Int:
                    std::cout << std::any_cast<int>(param->defaultValue);
                    break;
                case ParamType::Double:
                    std::cout << std::any_cast<double>(param->defaultValue);
                    break;
                case ParamType::Bool:
                    std::cout << (std::any_cast<bool>(param->defaultValue) ? "true" : "false");
                    break;
                case ParamType::StringList:
                case ParamType::Json:
                case ParamType::Selection:
                case ParamType::Path:
                    std::cout << std::any_cast<std::string>(param->defaultValue);
                    break;
                }
            } catch (...) {
                std::cout << "?";
            }
            std::cout << ")";
            // Claude Generated 2026 - What used to live in prose, if it was written
            // down at all: the unit, the permitted values, and when the parameter
            // matters in the first place.
            if (!param->unit.empty())
                std::cout << " [" << param->unit << "]";
            if (!param->allowed.empty()) {
                std::cout << " {";
                for (size_t i = 0; i < param->allowed.size(); ++i)
                    std::cout << (i ? "|" : "") << param->allowed[i];
                std::cout << "}";
            }
            if (param->hasMinimum || param->hasMaximum) {
                std::cout << " range: ";
                if (param->hasMinimum)
                    std::cout << param->minimum;
                std::cout << "..";
                if (param->hasMaximum)
                    std::cout << param->maximum;
            }
            if (!param->relevantWhen.empty())
                std::cout << " (only when " << param->relevantWhen << ")";
            if (param->deprecated) {
                std::cout << " DEPRECATED";
                if (!param->replacedBy.empty())
                    std::cout << ", use " << param->replacedBy;
            }
            std::cout << std::endl;

            // Show help text
            std::cout << "      " << param->helpText << std::endl;

            // Show aliases if any
            if (!param->aliases.empty()) {
                std::cout << "      Aliases: ";
                for (size_t i = 0; i < param->aliases.size(); ++i) {
                    if (i > 0)
                        std::cout << ", ";
                    std::cout << param->aliases[i];
                }
                std::cout << std::endl;
            }
        }
    }
    std::cout << std::endl;
}

// Claude Generated: Print all available modules
void ParameterRegistry::printAllModules() const
{
    std::cout << "Available modules:" << std::endl;
    for (const auto& [module, params] : registry) {
        std::cout << "  " << module << " (" << params.size() << " parameters)" << std::endl;
    }
}

// Claude Generated: Generate default JSON for a module
json ParameterRegistry::getDefaultJson(const std::string& module) const
{
    json result = json::object();  // Initialize as empty object, not null — null causes type_error.306 when .value() is called on it

    auto it = registry.find(module);
    if (it == registry.end()) {
        return result;
    }

    for (const auto& param : it->second) {
        try {
            switch (param.type) {
            case ParamType::String:
                result[param.name] = std::any_cast<std::string>(param.defaultValue);
                break;
            case ParamType::Int:
                result[param.name] = std::any_cast<int>(param.defaultValue);
                break;
            case ParamType::Double:
                result[param.name] = std::any_cast<double>(param.defaultValue);
                break;
            case ParamType::Bool:
                result[param.name] = std::any_cast<bool>(param.defaultValue);
                break;
            // Claude Generated 2026 - Selection and Path are strings that carry a
            // meaning; StringList and Json hold their default as JSON text so a
            // structured parameter can have one at all. temp_regions is the case
            // that made this necessary: an array of objects, read from the
            // controller but never registered, so it had no default and never
            // appeared in -export_config.
            case ParamType::Selection:
            case ParamType::Path:
                result[param.name] = std::any_cast<std::string>(param.defaultValue);
                break;
            case ParamType::StringList:
            case ParamType::Json: {
                const std::string text = std::any_cast<std::string>(param.defaultValue);
                if (text.empty()) {
                    result[param.name] = (param.type == ParamType::StringList)
                        ? nlohmann::json::array()
                        : nlohmann::json::object();
                    break;
                }
                nlohmann::json parsed = nlohmann::json::parse(text, nullptr, false);
                if (parsed.is_discarded()) {
                    std::cerr << "Warning: default of parameter " << param.name
                              << " in module " << module << " is not valid JSON: "
                              << text << std::endl;
                    parsed = nlohmann::json::object();
                }
                result[param.name] = parsed;
                break;
            }
            }
        } catch (const std::bad_any_cast& e) {
            std::cerr << "Warning: Failed to cast default value for parameter "
                      << param.name << " in module " << module << std::endl;
        }
    }

    return result;
}

// Claude Generated: Validate registry for duplicates and type consistency
// Claude Generated 2026 - Modules as concepts, not just as buckets of parameters.
void ParameterRegistry::addModule(ModuleDefinition&& definition)
{
    if (definition.name.empty())
        return;
    module_registry[definition.name] = std::move(definition);
}

const ModuleDefinition* ParameterRegistry::findModule(const std::string& name) const
{
    const auto it = module_registry.find(name);
    return it == module_registry.end() ? nullptr : &it->second;
}

std::vector<ModuleDefinition> ParameterRegistry::modules() const
{
    std::vector<ModuleDefinition> result;
    result.reserve(module_registry.size());
    for (const auto& entry : module_registry)
        result.push_back(entry.second);
    return result;
}

std::vector<std::string> ParameterRegistry::modulesForCommand(const std::string& command) const
{
    std::vector<std::string> result;
    for (const auto& entry : module_registry) {
        const auto& commands = entry.second.commands;
        if (std::find(commands.begin(), commands.end(), command) != commands.end())
            result.push_back(entry.first);
    }
    return result;
}

bool ParameterRegistry::validateRegistry() const
{
    bool valid = true;

    // Check for duplicate parameter names within each module
    for (const auto& [module, params] : registry) {
        std::map<std::string, int> name_counts;

        for (const auto& param : params) {
            name_counts[param.name]++;
            if (name_counts[param.name] > 1) {
                std::cerr << "Error: Duplicate parameter '" << param.name
                          << "' in module '" << module << "'" << std::endl;
                valid = false;
            }

            // Check for alias conflicts
            for (const auto& alias : param.aliases) {
                name_counts[alias]++;
                if (name_counts[alias] > 1) {
                    std::cerr << "Error: Alias '" << alias
                              << "' conflicts with another name/alias in module '"
                              << module << "'" << std::endl;
                    valid = false;
                }
            }
        }
    }

    // Check for shared parameters across modules (type consistency). Claude Generated (June 2026):
    // this only prints developer diagnostics (it never affects `valid`) and previously spammed a
    // block of "Parameter 'X' has different types ..." warnings to stderr on every run. Gate it
    // behind a debug macro (same pattern as DEBUG_CONFIG_MANAGER) so normal output stays clean;
    // these same-name/different-type collisions are handled by the flat-flag routing at runtime.
#ifdef DEBUG_PARAMETER_REGISTRY
    std::map<std::string, std::pair<std::string, ParamType>> first_occurrence;

    for (const auto& [module, params] : registry) {
        for (const auto& param : params) {
            auto it = first_occurrence.find(param.name);
            if (it == first_occurrence.end()) {
                first_occurrence[param.name] = { module, param.type };
            } else {
                // Same parameter name in different modules - check type consistency
                if (it->second.second != param.type) {
                    std::cerr << "Warning: Parameter '" << param.name
                              << "' has different types in modules '"
                              << it->second.first << "' and '" << module << "'" << std::endl;
                }
            }
        }
    }
#endif

    return valid;
}

// Claude Generated: Reverse-lookup a parameter name to its owning module(s).
// Case-insensitive. Returns empty vector for unknown names. Multiple entries indicate
// the name is shared across modules (CLI2Json treats this as ambiguous).
std::vector<std::string> ParameterRegistry::findOwnerModules(const std::string& param_name) const
{
    std::string lower = param_name;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    // Claude Generated (Jun 2026): canonicalize hyphens to underscores for matching
    lower.erase(std::remove(lower.begin(), lower.end(), '-'), lower.end());
    auto it = name_to_modules_map.find(lower);
    if (it == name_to_modules_map.end()) {
        return {};
    }
    return it->second;
}

// Claude Generated: Resolve alias to canonical name (case-insensitive)
std::string ParameterRegistry::resolveAlias(const std::string& module, const std::string& alias) const
{
    auto module_it = alias_to_name_map.find(module);
    if (module_it == alias_to_name_map.end()) {
        return ""; // Module not found
    }

    // Try exact match first (fast path)
    auto alias_it = module_it->second.find(alias);
    if (alias_it != module_it->second.end()) {
        return alias_it->second;
    }

    // Try case-insensitive match (slower path for backward compatibility)
    std::string alias_lower = alias;
    std::transform(alias_lower.begin(), alias_lower.end(), alias_lower.begin(), ::tolower);
    // Claude Generated (Jun 2026): canonicalize hyphens to underscores for matching
    alias_lower.erase(std::remove(alias_lower.begin(), alias_lower.end(), '-'), alias_lower.end());

    for (const auto& entry : module_it->second) {
        std::string entry_key_lower = entry.first;
        std::transform(entry_key_lower.begin(), entry_key_lower.end(), entry_key_lower.begin(), ::tolower);
        entry_key_lower.erase(std::remove(entry_key_lower.begin(), entry_key_lower.end(), '-'), entry_key_lower.end());

        if (alias_lower == entry_key_lower) {
            return entry.second; // Return canonical name
        }
    }

    return ""; // Alias not found
}
