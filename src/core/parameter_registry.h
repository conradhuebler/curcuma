#pragma once

#include <any>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

// Claude Generated - missing json include (needed when building as submodule without PCH)
#include "external/json.hpp"
using json = nlohmann::json;

enum class ParamType { String,
    Int,
    Double,
    Bool,
    // Claude Generated 2026 - structured and semantic kinds. Their default value is
    // written as JSON text in the macro and parsed by getDefaultJson(); Selection
    // and Path are strings that carry their meaning, so a schema generator can say
    // "atom selection" instead of "string".
    StringList,   ///< default written as a JSON array, e.g. "[]"
    Json,         ///< default written as JSON, e.g. "[]" or "{}" -- see temp_regions
    Selection,    ///< FragString grammar: "1:10,15", "F1", "-1"
    Path };       ///< a file or directory path

/// How prominently a parameter should be offered.
///
/// Claude Generated 2026 - The help category is a *grouping* ("Walls", "SCF"), not
/// a level of exposure, and the two were sharing one field. Anything generating a
/// compact interface -- a tool schema, a dialog's first page -- needs to know which
/// handful of parameters actually matter. Advanced is the default so nothing
/// becomes prominent by accident; being primary has to be stated.
enum class ParamTier { Primary,
    Advanced,
    Expert };

struct ParameterDefinition {
    std::string name; // Kanonischer Name (z.B. "max_iterations")
    std::string module; // Zugehöriges Modul (z.B. "opt", "casino")
    ParamType type; // Datentyp
    std::any defaultValue; // Standardwert
    std::string helpText; // Beschreibung für die Hilfe-Ausgabe
    std::string category = "General"; // Gruppierung für die Hilfe
    std::vector<std::string> aliases; // Alternative Namen/Kurzformen

    // Claude Generated 2026 - Everything below is optional and comes from the
    // PARAM macro's trailing annotation string. The fields sit AFTER aliases and
    // all carry defaults, so the seven-value aggregate initialisation the
    // generator has always emitted keeps working unchanged.
    ParamTier tier = ParamTier::Advanced;
    std::vector<std::string> allowed;   ///< the permitted values, when there is a fixed set
    std::string unit;                   ///< "K", "fs", "A", "Eh", ... -- was prose before
    bool hasMinimum = false;
    double minimum = 0.0;
    bool hasMaximum = false;
    double maximum = 0.0;
    /// When this parameter matters at all, e.g. "wall_type!=none". qurcuma
    /// reimplements exactly this by hand today, writing whole blocks only when the
    /// feature is on so curcuma's own defaults survive otherwise.
    std::string relevantWhen;
    bool deprecated = false;
    std::string replacedBy;
};

/// What a module IS, as opposed to what parameters it has.
///
/// Claude Generated 2026 - `-list_modules` printed "simplemd (78 parameters)" and
/// nothing else, and the mapping from CLI command to module was recorded nowhere:
/// only ten of the thirty-six command names match a module, `md` means `simplemd`,
/// `dock` means `docking`, and a schema generator had no way to find that out.
struct ModuleDefinition {
    std::string name;
    std::string description;
    std::string category;
    std::vector<std::string> commands;  ///< CLI verbs this module configures
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

    // --- Modules ---------------------------------------------------------
    // Claude Generated 2026
    void addModule(ModuleDefinition&& definition);
    const ModuleDefinition* findModule(const std::string& name) const;
    std::vector<ModuleDefinition> modules() const;
    /// Which modules configure @p command ("md" -> {"simplemd"}). Empty when the
    /// command is not claimed by any module -- which is true for most of them.
    std::vector<std::string> modulesForCommand(const std::string& command) const;

    // Claude Generated: Inverse lookup — which modules own a given parameter name or alias?
    // Returns all module names that register the given name (or any alias resolving to it).
    // Case-insensitive. Empty vector means the name is unknown to the registry.
    // Drives flat-CLI auto-routing in main.cpp::CLI2Json.
    std::vector<std::string> findOwnerModules(const std::string& param_name) const;

private:
    ParameterRegistry() = default;
    std::map<std::string, std::vector<ParameterDefinition>> registry;
    std::map<std::string, ModuleDefinition> module_registry;  // Claude Generated 2026
    std::map<std::string, std::map<std::string, std::string>> alias_to_name_map;
    // Lowercased name/alias -> deduped list of owning modules. Built incrementally in addDefinition.
    std::map<std::string, std::vector<std::string>> name_to_modules_map;
};

// Deklaration der Initialisierungsfunktion, die vom generierten Code kommt
void initialize_generated_registry();
