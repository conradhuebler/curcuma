#pragma once

#include <any>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

enum class ParamType { String,
    Int,
    Double,
    Bool };

struct ParameterDefinition {
    std::string name; // Kanonischer Name (z.B. "max_iterations")
    std::string module; // Zugehöriges Modul (z.B. "opt", "casino")
    ParamType type; // Datentyp
    std::any defaultValue; // Standardwert
    std::string helpText; // Beschreibung für die Hilfe-Ausgabe
    std::string category = "General"; // Gruppierung für die Hilfe
    std::vector<std::string> aliases; // Alternative Namen/Kurzformen
};

class ParameterRegistry {
public:
    static ParameterRegistry& getInstance();

    void addDefinition(const std::string& module, ParameterDefinition&& def);
    const ParameterDefinition* findDefinition(const std::string& module, const std::string& alias) const;
    std::vector<ParameterDefinition> getForModule(const std::string& module) const;

    // Claude Generated: Auto-generated help system
    void printHelp(const std::string& module) const;
    void printAllModules() const;

    // Claude Generated: Default JSON generation from registry
    nlohmann::json getDefaultJson(const std::string& module) const;

    // Claude Generated: Registry validation
    bool validateRegistry() const;

    // Claude Generated: Alias resolution
    std::string resolveAlias(const std::string& module, const std::string& alias) const;

private:
    ParameterRegistry() = default;
    std::map<std::string, std::vector<ParameterDefinition>> registry;
    std::map<std::string, std::map<std::string, std::string>> alias_to_name_map;
};

// Deklaration der Initialisierungsfunktion, die vom generierten Code kommt
void initialize_generated_registry();
