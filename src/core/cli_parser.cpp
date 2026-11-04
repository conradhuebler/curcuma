/*
 * <Command-line parsing: argv -> controller JSON>
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
 */

// Claude Generated (Nov 2026): moved unchanged from main.cpp so the routing can be tested without
// the whole program (docs/debt_2026-11/PARAMS.md, migration step 0).

#include "cli_parser.h"

#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/method_factory.h"
#include "src/core/parameter_registry.h"

#include <cctype>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Claude Generated 2026: Recursive deep-merge — dst wins where it already has a value (CLI > file).
// Unlike nlohmann::merge_patch this never deletes on null and recurses through every nested object,
// fixing the 1-level-only behavior of the old -import_config merge.
void mergeJsonRecursive(json& dst, const json& src)
{
    if (!src.is_object()) {
        return;
    }
    for (auto it = src.begin(); it != src.end(); ++it) {
        if (!dst.contains(it.key())) {
            dst[it.key()] = it.value();
        } else if (dst[it.key()].is_object() && it.value().is_object()) {
            mergeJsonRecursive(dst[it.key()], it.value());
        }
        // else: dst already has a scalar/array at this key — CLI wins, keep dst.
    }
}

// Helper function to set nested JSON values from dot-notation keys - Claude Generated
void setNestedJsonValue(json& target, const std::string& dotKey, const json& value) {
    if (dotKey.find('.') == std::string::npos) {
        // No dots - use as flat key
        target[dotKey] = value;
        return;
    }

    // Split on dots and create nested structure
    std::vector<std::string> keys;
    std::stringstream ss(dotKey);
    std::string key;

    while (std::getline(ss, key, '.')) {
        keys.push_back(key);
    }

    // Navigate/create nested structure
    json* current = &target;
    for (size_t i = 0; i < keys.size() - 1; ++i) {
        if (!current->contains(keys[i]) || !(*current)[keys[i]].is_object()) {
            (*current)[keys[i]] = json::object();
        }
        current = &(*current)[keys[i]];
    }

    // Set the final value
    (*current)[keys.back()] = value;
}

json CLI2Json(int argc, char** argv)
{
    json controller;
    json key = json::object();  // Initialize as empty object (not null) so controller["opt"] is always usable
    if (argc < 2)
        return controller;

    std::string keyword = argv[1];
    keyword.erase(0, 1);

    // Claude Generated (October 2025): Global parameters that should be accessible
    // both at top level (controller[param]) and module level (controller[module][param])
    // ENHANCED: Added "method" to support global energy method specification
    // Claude Generated (March 2026): Added "gpu" for GPU acceleration control.
    //   -gpu <backend>: none (default, CPU) | cuda | rocm | vulkan. ("auto" was removed Oct 26, 2026: warning, CPU.)
    //   Applies to gfn1/gfn2 (USE_*_XTB builds) and gfnff (USE_* builds); an explicit
    //   backend not compiled in warns and falls back to CPU.
    // Claude Generated (Jul 2026): Added "charge","spin" — molecular charge/spin must reach
    // top-level controller["charge"] so -sp/-opt can apply them to the molecule
    // (setCharge reads controller["charge"], not the command-module namespace).
    std::set<std::string> global_params = {
        "verbosity", "threads", "method", "gpu",  // energy_method and gpu apply to all capabilities
        "gpu_device",   // Claude Generated (Sep 2026, multi-GPU): device index for -gpu (see docs/MULTI_GPU.md)
        "gpu_devices", "gpu_workers_per_device",  // batch workers spread over these devices
        "gpu_strict",   // Claude Generated (Sep 2026): any GPU fallback ends the run (exit 3), see src/core/gpu_fallback.h
        "gpu_memory_check",  // native gfn1/gfn2 GPU: refuse a basis that does not fit (default true)
        "gpu_sparse_integrals",  // native gfn1/gfn2 GPU: screened S/H0/multipole storage auto|on|off
        // Claude Generated (Sep 2026, multi-GPU step 3): multi-GPU eigensolve (docs/GPU_TUNING.md)
        "gpu_eigensolver_devices", "gpu_eigensolver_backend", "gpu_eigensolver_block",
        "gpu_eigensolver_min_nao", "gpu_eigensolver_fp32", "gpu_eigensolver_verify",
        "gpu_density_devices", "gpu_density_min_nao",  // split the pattern density (docs/GPU_TUNING.md)
        "charge", "spin",  // molecular charge/spin (top-level, not module-scoped)
        "export_run", // Export current run configuration
        "import_config", // Import custom configuration
        "bak",   // Files to copy back from BMT output directory to CWD
        "no_bmt", // Disable BMT output directory (legacy: write to CWD)
        "noprogress" // Disable live progress bars globally (e.g. when redirecting to a file)
    };

    // Claude Generated (Jul 2026): modules that are sub-scopes of a calculation rather than
    // commands in their own right. Auto-routing a flat flag into one of these is the intended,
    // everyday behaviour ("-sp mol.xyz -cn_cutoff_bohr 5.5" -> controller["gfnff"]), so it is
    // only reported at info level. Moving a flag into another CAPABILITY module means the user
    // almost certainly aimed it at the active command and it silently went elsewhere -> warn.
    // This is the safety net for the class of bug where ConfSearch (unregistered) lost
    // -opt_method to polymerbuild and -thermostat to simplemd without any diagnostic.
    std::set<std::string> scope_modules(MethodFactory::methodParameterScopes().begin(),
                                        MethodFactory::methodParameterScopes().end());
    scope_modules.insert({"ripser", "rmsd", "ancopt", "native_lbfgs"});

    // Claude Generated (October 2025): CLI keyword to module name mapping
    // Maps command-line keywords (e.g., -md) to actual module names (e.g., simplemd)
    std::map<std::string, std::string> keyword_to_module = {
        {"md", "simplemd"},
        {"opt", "opt"},
        {"sp", "opt"},  // single point also uses opt module
        {"interaction", "interaction"},
        {"confscan", "confscan"},
        {"rmsd", "rmsd"},
        {"analysis", "analysis"},
        {"hessian", "hessian"},
        {"casino", "casino"},
        // Claude Generated (Nov 2026): the classes read modules "docking"/"trajectoryanalysis",
        // the CLI keywords are -dock/-traj; without the mapping their flags were never routed.
        {"dock", "docking"},
        {"traj", "trajectoryanalysis"}
    };

    // Get actual module name (for ConfigManager)
    std::string module_name = keyword;
    if (keyword_to_module.count(keyword) > 0) {
        module_name = keyword_to_module[keyword];
    }

    auto isStrictNumber = [](const std::string& s) -> bool {
        try {
            size_t pos = 0;
            std::stod(s, &pos);
            return pos == s.size();
        } catch (const std::exception&) {
            return false;
        }
    };

    // Claude Generated (Jul 2026, F2): a token starting with '-' followed by a digit or '.'
    // is a negative NUMBER value (e.g. "-charge -1", "-threshold -1e-6"), not a flag. The
    // value-consume logic below treats any '-'-prefixed next arg as a flag, which silently
    // drops negative numeric values (charge became boolean true -> get<int>()=1, breaking
    // all negatively-charged species). isFlag() returns false for negative numbers so they
    // flow into the stod() value branch instead.
    auto isFlag = [](const std::string& s) -> bool {
        if (s.empty() || s[0] != '-') return false;
        if (s.size() < 2) return true;             // lone "-"
        char c = s[1];
        // "-3", "-1.5", "-.5", "-1e-6" are numbers, not flags
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') return false;
        return true;
    };

    for (int i = 2; i < argc; ++i) {
        std::string current = argv[i];
        std::string sub = current.substr(0, 1);

        if (sub == "-") {
            // Claude Generated (Oct 2026): a flag token that contains whitespace ("-spin 1" as ONE
            // argv entry, e.g. from a shell variable that was not word-split) used to become the
            // boolean key "spin 1" and was silently ignored - the run then reported a closed-shell
            // energy for a radical. No parameter name contains whitespace, so refuse it.
            if (current.find_first_of(" \t\r\n") != std::string::npos) {
                CurcumaLogger::error("Command-line argument '" + current
                    + "' contains whitespace: flag and value must be separate arguments "
                      "(check shell quoting / variable expansion).");
                exit(1);
            }
            current.erase(0, 1);
            // Claude Generated (Jun 2026): canonicalize hyphens to underscores in CLI
            // parameter names so "-opt-method" matches the JSON key "opt_method".
            std::replace(current.begin(), current.end(), '-', '_');

            // Handle special verbosity shortcuts - Claude Generated
            if (current == "silent" || current == "quiet") {
                key["verbosity"] = 0;
                continue;
            } else if (current == "verbose") {
                key["verbosity"] = 3;
                continue;
            } else if (current == "plain") {
                // Claude Generated: Plain mode - no colors, no prefixes (like ORCA/Gaussian)
                CurcumaLogger::set_plain_mode(true);
                continue;
            } else if (current == "v") {
                // Handle -v N syntax for verbosity level
                if ((i + 1) < argc) {
                    try {
                        int verbosity_level = std::stoi(argv[i + 1]);
                        if (verbosity_level >= 0 && verbosity_level <= 4) {
                            key["verbosity"] = verbosity_level;
                            ++i; // Skip next argument
                            continue;
                        }
                    } catch (const std::exception&) {
                        // Fall through to normal processing
                    }
                }
            }

            // Claude Generated (October 2025 - CRITICAL FIX): Strip redundant keyword prefix from dotted parameters
            // Makes "-md.max_time 10" and "-max_time 10" synonymous within "-md" command context
            // Fixes SimpleMD double-nesting bug: controller["simplemd"]["md"]["max_time"] → controller["simplemd"]["max_time"]
            if (current.find('.') != std::string::npos) {
                size_t dot_pos = current.find('.');
                std::string param_prefix = current.substr(0, dot_pos);
                std::string param_key = current.substr(dot_pos + 1);

                // Resolve prefix to module name using keyword_to_module map (e.g., "md" → "simplemd")
                std::string prefix_module = param_prefix;
                if (keyword_to_module.count(param_prefix) > 0) {
                    prefix_module = keyword_to_module[param_prefix];
                }

                // If prefix matches current command's keyword or module name, strip it
                // Example: "-md input.xyz -md.max_time 10" → keyword="md", module_name="simplemd"
                //          param_prefix="md" matches keyword → strip → current="max_time"
                // Preserves cross-module routing: "-md -rmsd.method subspace" keeps "rmsd.method"
                if (param_prefix == keyword || prefix_module == module_name) {
                    current = param_key;  // Strip prefix: "md.max_time" → "max_time"
                }
            }

            if ((i + 1) >= argc || isFlag(argv[i + 1]) || argv[i + 1] == std::string("true") || argv[i + 1] == std::string("+")) {
                // Claude Generated (May 2026): Multi-value accumulation for -bak flag
                if (current == "bak" && key.contains("bak")) {
                    if (key["bak"].is_array()) {
                        key["bak"].push_back(true);
                    } else {
                        key["bak"] = json::array({key["bak"], true});
                    }
                } else {
                    setNestedJsonValue(key, current, true);
                }
            } else if (argv[i + 1] == std::string("false")) {
                setNestedJsonValue(key, current, false);
                ++i;
            } else {
                std::string next = argv[i + 1];
                //       std::cout << "next: " << next << std::endl;

                bool isNumber = true;
                bool isVector = next.find("|") != std::string::npos || next.find(",") != std::string::npos || next.find(":") != std::string::npos;
                //      std::cout << "isNumber: " << isNumber << std::endl
                //                  << "isVector: " << isVector << std::endl;
                if (isVector) {
                    isNumber = false;
                }
                if (!isVector) {
                    // Claude Generated (Nov 2026): strict number check. std::stod accepts numeric
                    // prefixes ("2026run.json" -> 2026, "1abc.log" -> 1) and throws out_of_range
                    // for "1e999"; only a token consumed completely by stod is a number.
                    if (!isStrictNumber(next)) {
                        isNumber = false;
                        isVector = true;
                    }
                }
                // std::cout << "isNumber: " << isNumber << std::endl
                //             << "isVector: " << isVector << std::endl;

                // Claude Generated (Jul 2026): a numeric-looking token (e.g. "1.0", "2.0")
                // must still be stored as a JSON STRING when the target parameter is
                // declared String (e.g. confscan.slx: "1.0" or "1.0,2.0" per confscan.h).
                // Otherwise setNestedJsonValue() below stores a JSON number, and the later
                // ConfigManager::get<std::string>() -> nlohmann from_json<string> merge
                // throws json::type_error.302 ("type must be string, but is number"),
                // which is uncaught in ConfigManager::getFromModule() and aborts (signal 6).
                // Comma/pipe/colon-separated values ("1.0,2.0") already take the string path
                // via isVector above; this closes the gap for single bare numbers.
                if (isNumber) {
                    std::string mod = module_name, param = current;
                    size_t dot = current.find('.');
                    if (dot != std::string::npos) {
                        std::string pfx = current.substr(0, dot);
                        mod = keyword_to_module.count(pfx) ? keyword_to_module[pfx] : pfx;
                        param = current.substr(dot + 1);
                    }
                    const ParameterDefinition* def = ParameterRegistry::getInstance().findDefinition(mod, param);
                    if (!def) {
                        for (const auto& owner : ParameterRegistry::getInstance().findOwnerModules(param)) {
                            if ((def = ParameterRegistry::getInstance().findDefinition(owner, param)))
                                break;
                        }
                    }
                    if (def && def->type == ParamType::String)
                        isNumber = false;
                }
                if (isNumber) {
                    // Claude Generated (May 2026): Multi-value accumulation for -bak
                    if (current == "bak" && key.contains("bak")) {
                        if (key["bak"].is_array()) {
                            key["bak"].push_back(std::stod(next));
                        } else {
                            key["bak"] = json::array({key["bak"], std::stod(next)});
                        }
                    } else {
                        setNestedJsonValue(key, current, std::stod(next));
                    }
                } else if (isVector) {
                    //        std::cout << next << std::endl;
                    if (current == "bak" && key.contains("bak")) {
                        if (key["bak"].is_array()) {
                            key["bak"].push_back(next);
                        } else {
                            key["bak"] = json::array({key["bak"], next});
                        }
                    } else {
                        setNestedJsonValue(key, current, next);
                    }
                } else {
                    if (current == "bak" && key.contains("bak")) {
                        if (key["bak"].is_array()) {
                            key["bak"].push_back(next);
                        } else {
                            key["bak"] = json::array({key["bak"], next});
                        }
                    } else {
                        setNestedJsonValue(key, current, next);
                    }
                }
                ++i;
            }
        }
    }

    // Claude Generated 2026: Auto-route flat CLI flags to their owning module via ParameterRegistry.
    // Lets "-static_charges true" land in controller["gfnff"]["static_charges"] without requiring
    // dotted form "-gfnff.static_charges true". Same-name-in-command-module wins; truly ambiguous
    // names (multiple owners, none matching the current command) emit a warning and stay put.
    {
        std::vector<std::string> reroute_remove;
        json reroute_targets; // module -> {param_name: value}
        for (auto it = key.begin(); it != key.end(); ++it) {
            const std::string& pname = it.key();
            if (global_params.count(pname) > 0)
                continue;
            if (it.value().is_object())
                continue;
            if (pname.find('.') != std::string::npos)
                continue;

            auto owners = ParameterRegistry::getInstance().findOwnerModules(pname);
            if (owners.empty())
                continue; // legacy / unregistered flag — preserve current behavior

            bool matches_current = false;
            for (const auto& owner : owners) {
                if (owner == module_name) {
                    matches_current = true;
                    break;
                }
            }
            if (matches_current)
                continue;

            // Claude Generated (Oct 16, 2026): a name several modules register (-electronic_temperature:
            // tblite, ulysses, xtb) is routed to the module of the selected -method when that module is
            // one of the owners. Before, it stayed in the command module with a warning and no method
            // ever read it (-electronic_temperature 40000 did not change a native gfn2 energy).
            if (owners.size() > 1 && key.contains("method") && key["method"].is_string()) {
                std::string method_scope = MethodFactory::parameterScopeForMethod(key["method"].get<std::string>());
                // Oct 2026: a shared SCF knob (-max_iterations, -diis_start: module "scf" via aliases, and others)
                // goes to "scf" when the method's own scope does not register it but the method reads "scf".
                if (std::find(owners.begin(), owners.end(), method_scope) == owners.end()
                    && std::find(owners.begin(), owners.end(), std::string("scf")) != owners.end()
                    && MethodFactory::readsSharedScfScope(key["method"].get<std::string>()))
                    method_scope = "scf";
                if (!method_scope.empty()
                    && std::find(owners.begin(), owners.end(), method_scope) != owners.end()) {
                    CurcumaLogger::info("CLI flag -" + pname + " is registered in several modules; routed to -"
                        + method_scope + "." + pname + " for the selected method \"" + key["method"].get<std::string>() + "\".");
                    reroute_targets[method_scope][pname] = it.value();
                    reroute_remove.push_back(pname);
                    continue;
                }
            }

            if (owners.size() > 1) {
                std::string candidates;
                for (size_t i = 0; i < owners.size(); ++i) {
                    if (i)
                        candidates += ", ";
                    candidates += owners[i];
                }
                CurcumaLogger::warn("CLI flag -" + pname + " is registered in multiple modules ("
                    + candidates + ") and none match the active command's module \"" + module_name
                    + "\". Use -<module>." + pname + " to disambiguate. Keeping flag in command module.");
                continue;
            }

            // Claude Generated (Jul 2026): announce the move. Silent rerouting hid the fact that
            // -opt_method / -thermostat / -restart never reached ConfSearch for months.
            {
                const std::string msg = "CLI flag -" + pname + " is not a parameter of the active"
                    " command's module \"" + module_name + "\"; it is owned by \"" + owners[0]
                    + "\" and was routed to -" + owners[0] + "." + pname + ".";
                if (scope_modules.count(owners[0]) > 0)
                    CurcumaLogger::info(msg);
                else
                    CurcumaLogger::warn(msg + " Use the dotted form if that is not what you meant.");
            }
            reroute_targets[owners[0]][pname] = it.value();
            reroute_remove.push_back(pname);
        }
        for (const auto& k : reroute_remove)
            key.erase(k);
        for (auto& [mod, params] : reroute_targets.items()) {
            if (!key.contains(mod) || !key[mod].is_object())
                key[mod] = json::object();
            for (auto& [pname, pval] : params.items())
                key[mod][pname] = pval;
        }
    }

    // Claude Generated 2025: Extract module-specific parameters to top level
    // Parameters like "rmsd.method" should be at controller["rmsd"]["method"], NOT controller["confscan"]["rmsd"]["method"]
    json module_params;
    std::vector<std::string> keys_to_remove;

    for (auto& [param_name, param_value] : key.items()) {
        // Claude Generated 2025: Check if this is a module-specific parameter
        // Could be either:
        // 1. Flat key with dots: "rmsd.method" (stored as param_name contains dot)
        // 2. Nested structure: setNestedJsonValue already created key["rmsd"]["method"]
        //    so param_name = "rmsd" and param_value = {"method": "subspace"}

        bool is_flat_dotted = param_name.find('.') != std::string::npos;
        bool is_nested_object = param_value.is_object();

        if (is_flat_dotted) {
            // Handle flat dot notation: "rmsd.method" or "-md.max_time"
            size_t dot_pos = param_name.find('.');
            std::string module_name = param_name.substr(0, dot_pos);
            std::string param_key = param_name.substr(dot_pos + 1);

            if (module_name != keyword) {
                // Claude Generated (October 2025 - CRITICAL FIX): Map keywords to actual module names
                // e.g., "-md.max_time" has module_name="md" from Punkt-notation,
                // but must be routed to module_params["simplemd"] (the actual module)
                std::string target_module = module_name;
                if (keyword_to_module.count(module_name) > 0) {
                    target_module = keyword_to_module[module_name];
                }
                if (!module_params.contains(target_module)) {
                    module_params[target_module] = json::object();
                }
                setNestedJsonValue(module_params[target_module], param_key, param_value);
                keys_to_remove.push_back(param_name);
            } else {
                // Claude Generated (October 2025 - CRITICAL FIX): Handle "-keyword.param" (e.g., "-md.max_time")
                // These should be FLAT in controller["simplemd"], not nested as key["md"]["max_time"]
                // Extract directly to top level to avoid double-nesting
                setNestedJsonValue(key, param_key, param_value);
                keys_to_remove.push_back(param_name);
            }
        } else if (is_nested_object && param_name != keyword) {
            // Handle nested structure for OTHER modules (not this command's keyword)
            // param_name = "rmsd", param_value = {"method": "subspace"}
            //
            // Claude Generated (Jul 2026): map CLI keywords to module names here too. Because
            // setNestedJsonValue() already split "md.rmsd_mtd" into key["md"]["rmsd_mtd"], this
            // branch (not the flat-dotted one above) is what handles "-md.x" issued from another
            // command -- and it used to store it verbatim under controller["md"], which no module
            // reads. So "curcuma -confsearch mol.xyz -md.rmsd_mtd false" was silently ignored.
            std::string target_module = param_name;
            if (keyword_to_module.count(param_name) > 0)
                target_module = keyword_to_module[param_name];
            if (module_params.contains(target_module) && module_params[target_module].is_object()) {
                for (auto& [k, v] : param_value.items())
                    module_params[target_module][k] = v; // merge, do not clobber a sibling block
            } else {
                module_params[target_module] = param_value;
            }
            keys_to_remove.push_back(param_name);
        }
        // NOTE: We KEEP nested structures where param_name == keyword!
        // E.g., key["md"] = {"max_time": 10} stays in key for backward compat
        // It will be stored in both controller["simplemd"] and controller["md"]
    }

    // Remove extracted module parameters from main command params
    for (const auto& remove_key : keys_to_remove) {
        key.erase(remove_key);
    }

    // Claude Generated (October 2025): Extract global parameters for explicit access
    json global_values;
    for (const auto& param : global_params) {
        if (key.count(param) > 0) {
            global_values[param] = key[param];
        }
    }

    // Build controller with proper structure using actual module name
    // This enables ConfigManager to find parameters under the correct module name
    controller[module_name] = key;  // Changed from keyword to module_name

    // Claude Generated (October 2025): Backward compatibility - keep old keyword-based access
    // Some tests use "-md.max_time" which creates controller["md"],
    // but we now route to controller["simplemd"] via module_name
    // Store under both for compatibility
    if (module_name != keyword) {
        controller[keyword] = key;  // Also keep the old keyword for backward compat
    }

    // Merge module-specific parameters at top level
    for (auto& [mod_name, module_config] : module_params.items()) {
        controller[mod_name] = module_config;  // Direct assignment (may contain nested structure from setNestedJsonValue)
    }

    // Ensure global parameters are always set at top level AND in global section
    for (const auto& [param, value] : global_values.items()) {
        controller[param] = value;  // Top-level access (backward compat)
        controller["global"][param] = value;  // Explicit global namespace
    }
    return controller;
}
