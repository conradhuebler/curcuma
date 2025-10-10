#include "parameter_registry.h"
#include "json.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>

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
                }
            } catch (...) {
                std::cout << "?";
            }
            std::cout << ")" << std::endl;

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
    json result;

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
            }
        } catch (const std::bad_any_cast& e) {
            std::cerr << "Warning: Failed to cast default value for parameter "
                      << param.name << " in module " << module << std::endl;
        }
    }

    return result;
}

// Claude Generated: Validate registry for duplicates and type consistency
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

    // Check for shared parameters across modules (type consistency)
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

    return valid;
}

// Claude Generated: Resolve alias to canonical name
std::string ParameterRegistry::resolveAlias(const std::string& module, const std::string& alias) const
{
    auto module_it = alias_to_name_map.find(module);
    if (module_it == alias_to_name_map.end()) {
        return ""; // Module not found
    }

    auto alias_it = module_it->second.find(alias);
    if (alias_it == module_it->second.end()) {
        return ""; // Alias not found
    }

    return alias_it->second;
}
