/*
 * <Curcuma main file.>
 * Copyright (C) 2019 - 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *               2024 Gerd Gehrisch <gg27fyla@student.freiberg.de>
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
#include "src/core/energy_calculators/qm_methods/eht.h"
#include "src/core/energy_calculators/qm_methods/orcainterface.h"
#include "src/core/fileiterator.h"
#include "src/core/imagewriter.hpp"
#include "src/core/molecule.h"

#include "src/capabilities/analysenciplot.h"
#include "src/capabilities/analysis.h"
#include "src/capabilities/confscan.h"
#include "src/capabilities/confsearch.h"
#include "src/capabilities/confstat.h"
#include "src/capabilities/curcumaopt.h"
// Modern optimizer system - Claude Generated (simplified)
#include "src/capabilities/docking.h"
#include "src/capabilities/hessian.h"
#include "src/capabilities/casino.h"
#include "src/capabilities/optimisation/modern_optimizer_simple.h"
#include "src/capabilities/nebdocking.h"
#include "src/capabilities/pairmapper.h"
#include "src/capabilities/persistentdiagram.h"
#include "src/capabilities/qmdfffit.h"
#include "src/capabilities/rmsd.h"
#include "src/capabilities/rmsdtraj.h"
#include "src/capabilities/simplemd.h"
#include "src/capabilities/trajectoryanalysis.h"

#include "src/tools/general.h"
#include "src/tools/info.h"

// Claude Generated: Parameter registry system
#include "generated/parameter_registry.h"
#include "src/core/parameter_registry.h"

#include "src/capabilities/optimiser/OptimiseDipoleScaling.h"
#include "src/capabilities/optimisation/modern_optimizer_simple.h"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include <string>
#include <vector>

#ifdef C17
#ifndef _WIN32
#include <filesystem>
namespace fs = std::filesystem;
#endif
#endif
// #include "ripser.h"

#ifndef _WIN32
#if __GNUC__
// Thanks to
// https://stackoverflow.com/questions/77005/how-to-automatically-generate-a-stacktrace-when-my-program-crashes
#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

void bt_handler(int sig)
{
    void* array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Curcuma crashed. Although this is probably unintended, it happened anyway.\n Some kind of backtrace will be printed out!\n\n");
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    fprintf(stderr, "Good-By\n");
#ifdef C17
#ifndef _WIN32
    //  std::filesystem::remove("stop");
#endif
#else
    remove("stop");
#endif
    exit(1);
}

void ctrl_c_handler(int s)
{
    if (std::filesystem::exists("stop")) {
#ifdef C17
#ifndef _WIN32
        std::filesystem::remove("stop");
#endif
#else
        remove("stop");
#endif
        printf("Caught stop signal a second time.\nWill exit now!\n\n");
        exit(1);
    } else {
        printf("Caught stop signal\nWill try to stop current stuff!\n");
        std::ofstream output("stop");
    }
}

#endif
#endif

#include "json.hpp"
using json = nlohmann::json;

void Distance(const Molecule &mol, char **argv)
{
    int donor = stoi(std::string(argv[3]));
    int proton = stoi(std::string(argv[4]));
    int acceptor = stoi(std::string(argv[5]));
    std::cout << "Using atoms " << donor << " " << proton << " " << acceptor << std::endl;
    std::cout << "Donor ";
    mol.printAtom(donor);
    std::cout << std::endl
              << "Proton: ";
    mol.printAtom(proton);
    std::cout << std::endl
              << "Acceptor: ";
    mol.printAtom(acceptor);
    std::cout << std::endl
              << "Hydrogen Bond Angle: " << mol.CalculateAngle(donor - 1, proton - 1, acceptor - 1) << std::endl;
    std::cout << "Hydrogen bond length " << mol.CalculateDistance(proton - 1, acceptor - 1) << std::endl;
}

double DotProduct(const Eigen::Vector3d& pos1, const Eigen::Vector3d& pos2)
{
    return pos1.dot(pos2);
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
    json key;
    if (argc < 2)
        return controller;

    std::string keyword = argv[1];
    keyword.erase(0, 1);

    // Claude Generated (October 2025): Global parameters that should be accessible
    // both at top level (controller[param]) and module level (controller[module][param])
    // ENHANCED: Added "method" to support global energy method specification
    std::set<std::string> global_params = {
        "verbosity", "threads", "method",  // energy_method applies to all capabilities
        "export_run", "export-run", // Export current run configuration
        "import_config", "import-config" // Import custom configuration
    };

    // Claude Generated (October 2025): CLI keyword to module name mapping
    // Maps command-line keywords (e.g., -md) to actual module names (e.g., simplemd)
    std::map<std::string, std::string> keyword_to_module = {
        {"md", "simplemd"},
        {"opt", "opt"},
        {"sp", "opt"},  // single point also uses opt module
        {"confscan", "confscan"},
        {"rmsd", "rmsd"},
        {"analysis", "analysis"},
        {"hessian", "hessian"},
        {"casino", "casino"}
    };

    // Get actual module name (for ConfigManager)
    std::string module_name = keyword;
    if (keyword_to_module.count(keyword) > 0) {
        module_name = keyword_to_module[keyword];
    }

    for (int i = 2; i < argc; ++i) {
        std::string current = argv[i];
        std::string sub = current.substr(0, 1);

        if (sub == "-") {
            current.erase(0, 1);

            // Handle special verbosity shortcuts - Claude Generated
            if (current == "silent" || current == "quiet") {
                key["verbosity"] = 0;
                continue;
            } else if (current == "verbose") {
                key["verbosity"] = 3;
                continue;
            } else if (current == "v") {
                // Handle -v N syntax for verbosity level
                if ((i + 1) < argc) {
                    try {
                        int verbosity_level = std::stoi(argv[i + 1]);
                        if (verbosity_level >= 0 && verbosity_level <= 3) {
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

            if ((i + 1) >= argc || argv[i + 1][0] == '-' || argv[i + 1] == std::string("true") || argv[i + 1] == std::string("+")) {
                setNestedJsonValue(key, current, true);
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
                    try {
                        // std::cout << "stod: " << std::stod(next) << std::endl;
                        std::stod(next);

                    } catch (const std::invalid_argument&) {
                        isNumber = false;
                        isVector = true;
                    }
                }
                // std::cout << "isNumber: " << isNumber << std::endl
                //             << "isVector: " << isVector << std::endl;
                if (isNumber) {
                    setNestedJsonValue(key, current, std::stod(next));
                } else if (isVector) {
                    //        std::cout << next << std::endl;
                    setNestedJsonValue(key, current, next);
                } else {
                    setNestedJsonValue(key, current, next);
                }
                ++i;
            }
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
            module_params[param_name] = param_value;
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

    // Claude Generated: DEBUG - Show what's in key before storing to controller
    std::cerr << "[CLI2Json DEBUG] keyword=" << keyword << ", module_name=" << module_name << std::endl;
    std::cerr << "[CLI2Json DEBUG] key object content:" << std::endl;
    std::cerr << key.dump(2) << std::endl;

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

// Structured Command Dispatch System - Claude Generated
struct CapabilityInfo {
    std::string description;
    std::string category;
    std::vector<std::string> supported_formats;
    std::function<int(const json&, int, char**)> handler;
};

// Capability handler functions - Claude Generated
int executeAnalysis(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        UnifiedAnalysis dummy(json{}, true);
        dummy.printHelp();
        return 0;
    }

    // Extract analysis-specific configuration - Claude Generated
    json analysis_config = controller.contains("analysis") ? controller["analysis"] : controller;

    auto* analysis = new UnifiedAnalysis(analysis_config, false);
    analysis->setFileName(argv[2]);
    analysis->start();
    delete analysis;
    return 0;
}

int executeRMSD(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for rmsd calculation as follows\ncurcuma -rmsd A.xyz B.xyz" << std::endl;
        std::cerr << "Please use curcuma for rmsd calculation as follows\ncurcuma -rmsd AB.xyz" << std::endl;
        return 1;
    }

    Molecule molecule1, molecule2;
    FileIterator file1, file2;
    std::string reffile, tarfile;

    if (std::filesystem::exists(std::string(argv[2]))) {
        file1.setFile(argv[2]);
        molecule1 = file1.Next();
        reffile = file1.Basename();
    }
    if (argc > 3 && std::filesystem::exists(argv[3])) {
        file2.setFile(argv[3]);
        tarfile = file2.Basename();
        molecule2 = file2.Next();
    } else {
        tarfile = file1.Basename();
        molecule2 = file1.Next();
    }

    if (molecule1.AtomCount() == 0 || molecule2.AtomCount() == 0) {
        std::cout << "At least one structure is empty:\n";
        std::cout << argv[2] << " " << molecule1.AtomCount() << " atoms" << std::endl;
        if (argc > 3) std::cout << argv[3] << " " << molecule2.AtomCount() << " atoms" << std::endl;
        return 1;
    }

    json rmsd_config = controller.value("rmsd", json::object());
    auto* driver = new RMSDDriver(rmsd_config, false);
    driver->setReference(molecule1);
    driver->setTarget(molecule2);
    driver->start();
    std::cout << "RMSD for two molecules " << driver->RMSD() << std::endl;

    driver->ReferenceAligned().writeXYZFile(reffile + ".centered.xyz");
    driver->TargetAligned().writeXYZFile(tarfile + ".centered.xyz");
    driver->TargetReorderd().writeXYZFile(tarfile + ".reordered.xyz");

    delete driver;
    return 0;
}

int executeSinglePoint(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for energy calculation as follows:\ncurcuma -sp input.xyz" << std::endl;
        return 1;
    }

    json sp_controller = controller;
    json opt_params = sp_controller.contains("opt") ? sp_controller["opt"] : json{};
    opt_params["single_point"] = true;
    sp_controller["opt"] = opt_params;

    CurcumaOpt opt(sp_controller, false);
    opt.setFileName(argv[2]);
    opt.start();
    return 0;
}

// Additional capability handlers - Claude Generated
int executeOptimization(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        json safe_opt_config = controller.contains("opt") ? controller["opt"] : json{};
        ModernOptimization::ModernOptimizerDispatcher helper(safe_opt_config, false);
        helper.printHelp();
        return 0;
    }

    // Parse optional -optimizer parameter
    std::string optimizer_method = "auto";
    if (argc >= 5 && strcmp(argv[3], "-optimizer") == 0) {
        optimizer_method = argv[4];
        std::cout << "🧪 Using native Curcuma optimizer: " << optimizer_method << std::endl;
    }

    // Check if we should use modern native optimizers
    bool use_modern = (optimizer_method == "native_lbfgs" || optimizer_method == "lbfgs" ||
        optimizer_method == "diis" || optimizer_method == "rfo" || optimizer_method == "auto");

    if (use_modern) {

        try {
            auto molecule = std::make_unique<Molecule>(argv[2]);
            std::string method = controller.value("method", "uff");
            EnergyCalculator energy_calc(method, controller);
            energy_calc.setMolecule(molecule->getMolInfo());

            // Claude Generated (October 2025): Merge with default opt parameters from ParameterRegistry
            json opt_defaults = ParameterRegistry::getInstance().getDefaultJson("opt");
            json opt_config = MergeJson(opt_defaults, controller.contains("opt") ? controller["opt"] : json{});

            auto result = ModernOptimization::ModernOptimizerDispatcher::optimizeStructure(
                molecule.get(), optimizer_method, &energy_calc, opt_config);

            if (result.success) {
                // Claude Generated: Derive output filename from input basename like Legacy CurcumaOpt
                // Extract basename from argv[2] (e.g., "input.xyz" → "input")
                std::string filename(argv[2]);
                std::string basename = filename.size() >= 4 ?
                    filename.substr(0, filename.size() - 4) : filename;  // Remove last 4 chars (.xyz)
                std::string output_file = opt_config.value("output", basename + ".opt.xyz");

                molecule->writeXYZFile(output_file);
                CurcumaLogger::success_fmt("Optimized structure written to: {}", output_file);
                return 0;
            } else {
                // Modern optimizer failed, fall through to legacy
                CurcumaLogger::warn(fmt::format("Modern optimization failed: {}, using legacy optimizer", result.error_message));
            }
        } catch (const std::exception& e) {
            // Fall through to legacy code below
            CurcumaLogger::warn(fmt::format("Modern optimization failed: {}, using legacy optimizer", e.what()));
        }
    }

    // Legacy optimization (fallback for exceptions or non-modern methods)
    // Claude Generated 2025: Apply same parameter merging as modern optimizer to fix JSON null bug
    json legacy_controller = controller;
    json opt_defaults = ParameterRegistry::getInstance().getDefaultJson("opt");
    legacy_controller["opt"] = MergeJson(opt_defaults, controller.contains("opt") ? controller["opt"] : json{});

    CurcumaOpt opt(legacy_controller, false);
    opt.setFileName(argv[2]);
    opt.start();
    return 0;
}

int executeConfScan(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for conformation scan as follows:\ncurcuma -confscan conffile.xyz" << std::endl;
        std::cerr << "Additional arguments are:" << std::endl;
        std::cerr << "-writeXYZ  **** Write results to xyz files!" << std::endl;
        std::cerr << "-rank n    **** Write only the first n results!" << std::endl;
        std::cerr << "-reorder   **** Force reordering of structure!" << std::endl;
        std::cerr << "-heavy     **** Use only heavy atoms for rmsd calculation." << std::endl;
        std::cerr << "-maxenergy **** Maximal energy difference [kJ/mol]." << std::endl;
        return 1;
    }

    // Claude Generated 2025: Simplified - CLI2Json routes parameters, ConfigManager handles defaults
    // ConfScan's ConfigManager (Multi-Module: "confscan", "rmsd") will:
    // 1. Load defaults for both modules from ParameterRegistry
    // 2. Merge user parameters from controller["confscan"] and controller["rmsd"]
    // 3. Handle global parameters (verbosity, threads)
    // No manual merging needed here!

    auto* scan = new ConfScan(controller, false);  // Claude Generated: Explicit false for default verbosity level 1
    scan->setFileName(argv[2]);
    scan->start();
    int accepted = scan->AcceptedCount();
    int reorder_success = scan->ReorderSuccessfull();
    int reuse_count = scan->ReuseCount();
    int skipped_count = scan->ReorderSkippedCount();
    std::cout << accepted << " " << reorder_success << " " << reuse_count << " " << skipped_count << std::endl;
    delete scan;
    return 0;
}

int executeConfStat(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for conformation statistics as follows:\ncurcuma -confstat conffile.xyz" << std::endl;
        return 1;
    }

    auto* stat = new ConfStat(controller, false);  // Claude Generated: Explicit false for default verbosity level 1
    stat->setFileName(argv[2]);
    stat->start();
    delete stat;
    return 0;
}

int executeDocking(const json& controller, int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Please use curcuma for docking as follows:\ncurcuma -dock -host A.xyz -guest B.xyz -Step_x 10 -Step_y 10 -Step_z 10" << std::endl;
        return 1;
    }

    auto* docking = new Docking(controller, false);
    if (!docking->Initialise()) {
        docking->printError();
        delete docking;
        return 1;
    }
    docking->start();
    delete docking;
    return 0;
}

int executeSimpleMD(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Please use curcuma for molecular dynamics as follows:\ncurcuma -md input.xyz" << std::endl;
        SimpleMD help(json::object(), true);
        help.printHelp();
        return 0;
    }

    Molecule mol = Files::LoadFile(argv[2]);
    if (mol.AtomCount() == 0) {
        CurcumaLogger::error_fmt("Could not load molecule from file: {}", argv[2]);
        return 1;
    }

    // Claude Generated 2025: Simplified - SimpleMD's ConfigManager handles defaults
    // SimpleMD constructor uses ConfigManager("simplemd", controller) which:
    // 1. Loads defaults from ParameterRegistry
    // 2. Merges user parameters from controller["md"]
    // 3. Handles alias resolution automatically
    // No manual merging needed here!

    auto* md = new SimpleMD(controller, false);
    md->setFile(argv[2]);  // Claude Generated (October 2025): Set basename for trajectory file naming
    md->setMolecule(mol);
    md->Initialise();
    md->start();
    delete md;
    return 0;
}

int executeCasino(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        Casino dummy(json{}, true);
        dummy.printHelp();
        return 0;
    }
    // Claude Generated 2025: Apply parameter merging to fix JSON null issues
    json casino_defaults = ParameterRegistry::getInstance().getDefaultJson("casino");
    json casino_config = MergeJson(casino_defaults, controller.contains("casino") ? controller["casino"] : json{});

    auto* casino = new Casino(casino_config, false);
    casino->setFileName(argv[2]);
    casino->start();
    delete casino;
    return 0;
}

int executeTrajectoryAnalysis(const json& controller, int argc, char** argv) {
    if (argc < 3) {
        TrajectoryAnalysis dummy(json{}, true);
        dummy.printHelp();
        return 0;
    }

    auto* traj = new TrajectoryAnalysis(controller, false);
    traj->setFileName(argv[2]);
    traj->start();
    delete traj;
    return 0;
}

// Capability registry - Claude Generated
const std::map<std::string, CapabilityInfo> CAPABILITY_REGISTRY = {
    {"analysis", {"Unified molecular analysis (all formats, all properties)", "analysis",
                  {"XYZ", "VTF", "MOL2", "SDF", "PDB"}, executeAnalysis}},
    {"rmsd", {"RMSD calculation between structures", "analysis",
              {"XYZ", "VTF", "MOL2", "SDF"}, executeRMSD}},
    {"sp", {"Single point energy calculation", "calculation",
            {"XYZ", "VTF", "MOL2", "SDF"}, executeSinglePoint}},
    {"opt", {"Geometry optimization with various algorithms", "optimization",
             {"XYZ", "VTF", "MOL2", "SDF"}, executeOptimization}},
    {"confscan", {"Conformational scanning along reaction coordinates", "conformational",
                  {"XYZ", "MOL2", "SDF"}, executeConfScan}},
    {"confstat", {"Conformational statistics and analysis", "conformational",
                  {"XYZ", "MOL2", "SDF"}, executeConfStat}},
    {"dock", {"Molecular docking calculations", "docking",
              {"XYZ", "MOL2", "SDF"}, executeDocking}},
    {"md", {"Molecular dynamics simulation", "dynamics",
            {"XYZ", "VTF", "MOL2", "SDF"}, executeSimpleMD}},
    {"casino", {"Casino Monte Carlo simulation with enhanced sampling", "dynamics",
                {"XYZ", "VTF", "MOL2", "SDF", "PDB"}, executeCasino}},
    {"traj", {"Time-series analysis for molecular trajectories", "analysis",
              {"XYZ.trj", "VTF", "SDF"}, executeTrajectoryAnalysis}},
    // TODO: Add more capabilities here as they are converted
};

void showStructuredHelp(const std::string& category = "") {
    std::cout << "Curcuma - Computational Chemistry Toolkit" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;

    if (category.empty()) {
        // Show overview with enhanced format information
        std::cout << "Modern Molecular Modeling with Universal Format Support" << std::endl;
        std::cout << std::endl;

        // Auto-discover and group capabilities by category
        std::map<std::string, std::vector<std::string>> categories;
        std::set<std::string> all_formats;

        for (const auto& [name, info] : CAPABILITY_REGISTRY) {
            categories[info.category].push_back(name);
            for (const auto& format : info.supported_formats) {
                all_formats.insert(format);
            }
        }

        // Display capabilities by category with enhanced formatting
        for (const auto& [cat, capabilities] : categories) {
            std::string category_name = cat;
            category_name[0] = std::toupper(category_name[0]);

            std::cout << "📊 " << category_name << " Capabilities:" << std::endl;

            for (const auto& cap : capabilities) {
                const auto& info = CAPABILITY_REGISTRY.at(cap);

                // Enhanced formatting with icons
                std::string icon = "🔬";
                if (cat == "dynamics") icon = "⚡";
                else if (cat == "optimization") icon = "🎯";
                else if (cat == "conformational") icon = "🔄";
                else if (cat == "docking") icon = "🔗";
                else if (cat == "calculation") icon = "⚙️";

                std::cout << "  " << icon << " -" << cap
                         << std::string(std::max(1, 18 - static_cast<int>(cap.length())), ' ')
                         << info.description << std::endl;

                // Show supported formats for each capability
                std::cout << "      Formats: ";
                for (size_t i = 0; i < info.supported_formats.size(); ++i) {
                    std::cout << info.supported_formats[i];
                    if (i < info.supported_formats.size() - 1) std::cout << ", ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        // Enhanced format support section with auto-discovery
        std::cout << "🗂️ Supported File Formats:" << std::endl;
        std::cout << "  ┌─────────────────────────────────────────────────────────┐" << std::endl;
        std::cout << "  │ Atomistic:     XYZ, MOL2, SDF, PDB                     │" << std::endl;
        std::cout << "  │ Coarse-Grained: VTF (VMD Trajectory Format)            │" << std::endl;
        std::cout << "  │ Trajectories:   XYZ.trj, VTF multi-timestep           │" << std::endl;
        std::cout << "  │ All formats work transparently with all capabilities   │" << std::endl;
        std::cout << "  └─────────────────────────────────────────────────────────┘" << std::endl;
        std::cout << std::endl;

        // Usage examples
        std::cout << "💡 Quick Start Examples:" << std::endl;
        std::cout << "  curcuma -analysis molecule.xyz             # Complete molecular analysis" << std::endl;
        std::cout << "  curcuma -opt protein.pdb -method gfn2      # Quantum optimization" << std::endl;
        std::cout << "  curcuma -casino polymer.vtf -steps 10000   # Casino Monte Carlo simulation" << std::endl;
        std::cout << "  curcuma -traj md_output.xyz -stride 10     # Trajectory analysis" << std::endl;
        std::cout << std::endl;

        std::cout << "📚 Get detailed help: curcuma -help [category]" << std::endl;
        std::cout << "   Available categories: ";
        bool first = true;
        for (const auto& [cat, capabilities] : categories) {
            if (!first) std::cout << ", ";
            std::cout << cat;
            first = false;
        }
        std::cout << std::endl;

    } else {
        // Enhanced category-specific help
        bool found = false;
        std::vector<std::string> found_capabilities;

        // Collect capabilities in this category
        for (const auto& [name, info] : CAPABILITY_REGISTRY) {
            if (info.category == category) {
                found_capabilities.push_back(name);
                found = true;
            }
        }

        if (found) {
            std::string category_name = category;
            category_name[0] = std::toupper(category_name[0]);

            std::cout << "📋 " << category_name << " Capabilities - Detailed Help" << std::endl;
            std::cout << "══════════════════════════════════════════════════════" << std::endl;
            std::cout << std::endl;

            for (const auto& cap : found_capabilities) {
                const auto& info = CAPABILITY_REGISTRY.at(cap);

                std::cout << "🔹 -" << cap << std::endl;
                std::cout << "   Description: " << info.description << std::endl;
                std::cout << "   Usage: curcuma -" << cap << " input_file [options]" << std::endl;
                std::cout << "   Supported formats: ";

                for (size_t i = 0; i < info.supported_formats.size(); ++i) {
                    std::cout << info.supported_formats[i];
                    if (i < info.supported_formats.size() - 1) std::cout << ", ";
                }
                std::cout << std::endl;

                // Add specific help references
                std::cout << "   Get detailed help: curcuma -" << cap << " (without arguments)" << std::endl;
                std::cout << std::endl;
            }

            // Category-specific examples
            if (category == "analysis") {
                std::cout << "💡 Analysis Examples:" << std::endl;
                std::cout << "  curcuma -analysis molecule.xyz -properties all" << std::endl;
                std::cout << "  curcuma -traj trajectory.vtf -export_timeseries true" << std::endl;
                std::cout << "  curcuma -rmsd ref.xyz target.xyz" << std::endl;
            } else if (category == "dynamics") {
                std::cout << "💡 Dynamics Examples:" << std::endl;
                std::cout << "  curcuma -md system.xyz -temperature 300 -steps 10000" << std::endl;
                std::cout << "  curcuma -casino polymer.vtf -move_type mixed -adaptive_step true" << std::endl;
            } else if (category == "optimization") {
                std::cout << "💡 Optimization Examples:" << std::endl;
                std::cout << "  curcuma -opt molecule.xyz -method uff" << std::endl;
                std::cout << "  curcuma -sp system.pdb -method gfn2" << std::endl;
            }
            std::cout << std::endl;

        } else {
            std::cout << "❌ Unknown category: " << category << std::endl;
            std::cout << std::endl;
            std::cout << "Available categories: ";

            std::set<std::string> available_categories;
            for (const auto& [name, info] : CAPABILITY_REGISTRY) {
                available_categories.insert(info.category);
            }

            bool first = true;
            for (const auto& cat : available_categories) {
                if (!first) std::cout << ", ";
                std::cout << cat;
                first = false;
            }
            std::cout << std::endl;
        }
    }
}

int main(int argc, char **argv) {
#ifndef _WIN32
#if __GNUC__
    signal(SIGINT, ctrl_c_handler);
    signal(SIGSEGV, bt_handler);
    signal(SIGABRT, bt_handler);
#endif
#endif

    General::StartUp(argc, argv);

    // Claude Generated: Initialize parameter registry early
    initialize_generated_registry();
    if (!ParameterRegistry::getInstance().validateRegistry()) {
        std::cerr << "Warning: Parameter registry validation failed" << std::endl;
    }

#ifdef C17
#ifndef _WIN32
    std::filesystem::remove("stop");
#endif
#else
    remove("stop");
#endif
    RunTimer timer(true);

    // Handle no arguments or help requests - Claude Generated
    if(argc < 2) {
        showStructuredHelp();
        exit(1);
    }

    std::string command = argv[1] + 1; // Remove '-' prefix

    // Handle help requests - Claude Generated
    if (command == "help" || command == "h") {
        if (argc >= 3) {
            showStructuredHelp(argv[2]); // Category-specific help
        } else {
            showStructuredHelp();
        }
        return 0;
    }

    // Handle parameter registry commands - Claude Generated (October 2025)
    if (command == "list-modules") {
        ParameterRegistry::getInstance().printAllModules();
        return 0;
    }

    if (command == "export-config") {
        if (argc < 3) {
            std::cerr << "Error: -export-config requires module name" << std::endl;
            std::cerr << "Usage: curcuma -export-config <module>" << std::endl;
            std::cerr << "Example: curcuma -export-config analysis" << std::endl;
            return 1;
        }
        std::string module = argv[2];
        json config = ParameterRegistry::getInstance().getDefaultJson(module);
        if (config.empty()) {
            std::cerr << "Error: Unknown module '" << module << "'" << std::endl;
            std::cerr << "Use: curcuma -list-modules to see available modules" << std::endl;
            return 1;
        }
        std::cout << config.dump(2) << std::endl;
        return 0;
    }

    if (command == "help-module") {
        if (argc < 3) {
            std::cerr << "Error: -help-module requires module name" << std::endl;
            std::cerr << "Usage: curcuma -help-module <module>" << std::endl;
            std::cerr << "Example: curcuma -help-module analysis" << std::endl;
            return 1;
        }
        std::string module = argv[2];
        ParameterRegistry::getInstance().printHelp(module);
        return 0;
    }

    json controller = CLI2Json(argc, argv);

    // Handle config import - Claude Generated (October 2025)
    // Now global parameter, always in controller["import_config"] if present
    if (controller.contains("import_config")) {
        std::string config_file = controller["import_config"];

        std::ifstream file(config_file);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open config file: " << config_file << std::endl;
            return 1;
        }

        json imported_config;
        try {
            file >> imported_config;
        } catch (const std::exception& e) {
            std::cerr << "Error: Failed to parse JSON from " << config_file << std::endl;
            std::cerr << "Details: " << e.what() << std::endl;
            return 1;
        }

        // Merge imported config into controller
        // CLI arguments have higher priority (don't overwrite existing keys)
        for (auto it = imported_config.begin(); it != imported_config.end(); ++it) {
            if (!controller.contains(it.key())) {
                controller[it.key()] = it.value();
            } else if (controller[it.key()].is_object() && it.value().is_object()) {
                // Merge nested objects (e.g., controller["analysis"] with imported["analysis"])
                for (auto nested_it = it.value().begin(); nested_it != it.value().end(); ++nested_it) {
                    if (!controller[it.key()].contains(nested_it.key())) {
                        controller[it.key()][nested_it.key()] = nested_it.value();
                    }
                }
            }
        }

        std::cout << "Loaded configuration from: " << config_file << std::endl;
    }

    // Handle run export - Claude Generated (October 2025)
    // Export current configuration AFTER all merging/importing
    // Now global parameter, always in controller["export_run"] if present
    if (controller.contains("export_run")) {
        std::string export_file = controller["export_run"];

        // Create export config without meta-parameters
        json export_config = controller;

        // Remove all meta-parameters (those used for I/O operations, not actual computation)
        // Both at top level and from module sub-dictionaries
        std::vector<std::string> meta_params = { "export_run", "export-run", "import_config", "import-config" };

        for (const auto& meta : meta_params) {
            export_config.erase(meta);

            // Also remove from all module sub-dictionaries
            for (auto it = export_config.begin(); it != export_config.end(); ++it) {
                if (it.value().is_object()) {
                    it.value().erase(meta);
                }
            }
        }

        std::ofstream outfile(export_file);
        if (!outfile.is_open()) {
            std::cerr << "Error: Could not write to " << export_file << std::endl;
            return 1;
        }

        outfile << export_config.dump(2) << std::endl;
        outfile.close();

        std::cout << "Exported current run configuration to: " << export_file << std::endl;
    }

    // Try structured dispatch first - Claude Generated
    auto it = CAPABILITY_REGISTRY.find(command);
    if (it != CAPABILITY_REGISTRY.end()) {
        return it->second.handler(controller, argc, argv);
    }

    // Fall back to legacy command handling for compatibility - Claude Generated
    if (false) { // Placeholder - will be removed as capabilities are migrated

        // if(strcmp(argv[1], "-rmsd") == 0)
        // {
        //     // MOVED TO STRUCTURED DISPATCH SYSTEM - Claude Generated
        //     [RMSD code moved to executeRMSD function]
        // }

        if (strcmp(argv[1], "-compare") == 0) {
            if (argc < 4) {
                std::cerr << "Please use curcuma to compare structures as follows\ncurcuma -compare A.xyz B.xyz -metric [rmsd, inertia, ripser]" << std::endl;
                exit(1);
            }

            // Parse file names
            std::string file_a = argv[2];
            std::string file_b = argv[3];

            // Parse metric option
            std::string metric = "rmsd"; // Default metric
            for (int i = 4; i < argc; ++i) {
                std::string arg = argv[i];
                if (arg == "-metric" && i + 1 < argc) {
                    metric = argv[i + 1];
                    break;
                }
            }

            // Load molecules
            Molecule molecule1, molecule2;
            FileIterator file1, file2;
            std::string reffile;
            std::string tarfile;

            if (std::filesystem::exists(std::string(argv[2]))) {
                file1.setFile(argv[2]);
                molecule1 = file1.Next();
                reffile = file1.Basename();
            }
            if (std::filesystem::exists(argv[3])) {
                file2.setFile(argv[3]);
                tarfile = file2.Basename();
                molecule2 = file2.Next();
            } else {
                tarfile = file1.Basename();
                molecule2 = file1.Next();
            }

            // Compare structures based on selected metric
            if (metric == "rmsd") {
                json config;
                config["rmsd"] = controller["compare"];
                RMSDDriver driver(config, true);
                driver.setReference(molecule1);
                driver.setTarget(molecule2);
                driver.start();
                std::cout << "RMSD: " << driver.RMSD() << " Å" << std::endl;
            } else if (metric == "inertia") {

                molecule1.CalculateRotationalConstants();
                molecule2.CalculateRotationalConstants();

                std::cout << "Inertia constants for " << reffile << ": " << std::endl;
                std::cout << "Ia: " << molecule1.Ia() << " MHz, Ib: " << molecule1.Ib() << " MHz, Ic: " << molecule1.Ic() << " MHz" << std::endl;
                std::cout << "Inertia constants for " << tarfile << ": " << std::endl;
                std::cout << "Ia: " << molecule2.Ia() << " MHz, Ib: " << molecule2.Ib() << " MHz, Ic: " << molecule2.Ic() << " MHz" << std::endl;
                double Ia_diff = std::abs(molecule1.Ia() - molecule2.Ia());
                double Ib_diff = std::abs(molecule1.Ib() - molecule2.Ib());
                double Ic_diff = std::abs(molecule1.Ic() - molecule2.Ic());
                std::cout << "Inertia differences: " << Ia_diff << ", " << Ib_diff << ", " << Ic_diff << " MHz" << std::endl;
            } else if (metric == "ripser") {
                PersistentDiagram pd(controller["compare"]);
                pd.setDistanceMatrix(molecule1.LowerDistanceVector());
                Eigen::MatrixXd image_a = pd.generateImage(pd.generatePairs());
                pd.setDistanceMatrix(molecule2.LowerDistanceVector());
                Eigen::MatrixXd image_b = pd.generateImage(pd.generatePairs());
                std::cout << "Persistence diagram for " << reffile << ": " << std::endl;
                std::cout << image_a << std::endl;
                std::cout << "Persistence diagram for " << tarfile << ": " << std::endl;
                std::cout << image_b << std::endl;
                double diff = (image_a - image_b).cwiseAbs().sum();
                std::cout << "Persistence diagram difference: " << diff << std::endl;
            } else {
                std::cerr << "Unknown metric: " << metric << std::endl;
                std::cerr << "Supported metrics: rmsd, inertia, ripser" << std::endl;
                exit(1);
            }
        // } else if (strcmp(argv[1], "-analysis") == 0) {
        //     // MOVED TO STRUCTURED DISPATCH SYSTEM - Claude Generated
        //     if (argc < 3) {
        //         UnifiedAnalysis::printHelp();
        //         return 0;
        //     }
        //     auto* analysis = new UnifiedAnalysis(controller, false);
        //     analysis->setFileName(argv[2]);
        //     analysis->start();
        //     delete analysis;
        //     return 0;

        // } else if (strcmp(argv[1], "-dock") == 0) {
        //     // MOVED TO STRUCTURED DISPATCH SYSTEM - Claude Generated
        //     [Docking code moved to executeDocking function]

        } else if (strcmp(argv[1], "-hbonds") == 0) {
            if (argc < 6) {
                std::cerr << "Please use curcuma for hydrogen bond analysis as follows\ncurcuma -hbonds A.xyz index_donor index_proton index_acceptor" << std::endl;
                return -1;
            }
            FileIterator file(argv[2]);
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                if (argc == 6) {
                    if (std::string(argv[1]).find("-hbonds") != std::string::npos) {
                        Distance(mol, argv);
                    }
                } else {
                    mol.print_geom();
                    std::cout << std::endl
                              << std::endl;
                    std::cout << mol.getGeometry() << std::endl;
                }
            }
        // } else if (strcmp(argv[1], "-confscan") == 0) {
        //     // MOVED TO STRUCTURED DISPATCH SYSTEM - Claude Generated
        //     [ConfScan code moved to executeConfScan function]

        // } else if (strcmp(argv[1], "-confstat") == 0) {
        //     // MOVED TO STRUCTURED DISPATCH SYSTEM - Claude Generated
        //     [ConfStat code moved to executeConfStat function]

        } else if (strcmp(argv[1], "-led") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for fragment assignment as follows:\ncurcuma -led input.xyz" << std::endl;
                return 0;
            }

            Molecule mol1 = Files::LoadFile(argv[2]);
            if (!mol1.Atoms().empty())
                mol1.printFragmente();

        } else if (strcmp(argv[1], "-hmap") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for hydrogen bond mapping as follows:\ncurcuma -hmap trajectory.xyz" << std::endl;
                return 0;
            }

            std::vector<std::pair<int, int>> pairs, elements;

            if (argc > 3) {
                for (std::size_t i = 3; i < argc; ++i) {
                    if (strcmp(argv[i], "-pair") == 0) {
                        if (i + 2 < argc) {
                            if (Tools::isInt(argv[i + 1]) && Tools::isInt(argv[i + 2])) {
                                int first = std::stoi(argv[i + 1]) - 1;
                                int second = std::stoi(argv[i + 2]) - 1;
                                ++i;
                                pairs.emplace_back(first, second);
                            } else {
                                int first = Elements::String2Element(argv[i + 1]);
                                int second = Elements::String2Element(argv[i + 2]);
                                ++i;
                                elements.emplace_back(first, second);
                            }
                        }
                    }
                    if (strcmp(argv[i], "-pairfile") == 0) {
                        if (i + 1 < argc) {
                            std::ifstream input(argv[i + 1]);
                            for (std::string line; getline(input, line);) {
                                std::vector<std::string> numbers = Tools::SplitString(line);
                                if (numbers.size() == 2) {
                                    if (Tools::isInt(numbers[0]) && Tools::isInt(numbers[1])) {
                                        int first = std::stoi(numbers[0]) - 1;
                                        int second = std::stoi(numbers[1]) - 1;
                                        pairs.emplace_back(first, second);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            PairMapper mapper;
            mapper.setFile(argv[2]);
            for (const std::pair<int, int>& pair : pairs)
                mapper.addPair(pair);
            for (const std::pair<int, int>& pair : elements)
                mapper.addElementPair(pair);

            mapper.FindPairs();

        } else if (strcmp(argv[1], "-nci") == 0) {
            if (argc < 4) {
                std::cerr << "Please use curcuma to post-process two RDG vs rho plots from NCIPLOT as follows:\ncurcuma -nci file1.dat file2.dat" << std::endl;
                std::cerr << "Additonal arguments are:" << std::endl;
                std::cerr << "-bins             **** Number of bins during indexing the file!" << std::endl;
                std::cerr << "-scale_d1         **** Scale minimal distance for file1.dat!" << std::endl;
                std::cerr << "-scale_d2         **** Scale minimal distance for file2.dat!" << std::endl;
                std::cerr << "-local_distance   **** Recalculate distance for every bin (false = default)" << std::endl;
                return 0;
            }

            AnalyseNCIPlot analyse(controller);
            analyse.setFiles(argv[2], argv[3]);
            analyse.start();

        // } else if (strcmp(argv[1], "-opt") == 0) {
        //     // MOVED TO STRUCTURED DISPATCH SYSTEM - Claude Generated
        //     [Optimization code moved to executeOptimization function]

        } else if (strcmp(argv[1], "-modern-opt") == 0) {
            // Claude Generated - Modern optimizer system showcase
            if (argc < 3) {
                std::cerr << "Please use curcuma for modern optimization as follows:\ncurcuma -modern-opt input.xyz [method]" << std::endl;
                std::cerr << "Available methods: lbfgspp, lbfgs, diis, rfo, auto" << std::endl;
                return 0;
            }

            try {
                // Load molecule
                Molecule molecule = Files::LoadFile(argv[2]);
                CurcumaLogger::header("Modern Optimization System - Claude Generated");
                CurcumaLogger::info_fmt("Input file: {}", argv[2]);
                CurcumaLogger::info_fmt("Atoms: {}", molecule.AtomCount());

                // Determine optimization method
                std::string method = "lbfgspp"; // Default
                if (argc >= 4) {
                    method = argv[3];
                }
                CurcumaLogger::info_fmt("Selected method: {}", method);

                // Create energy calculator (using controller settings)
                EnergyCalculator energy_calc("uff", controller); // Use UFF as default force field
                energy_calc.setMolecule(molecule.getMolInfo());

                // Show available optimizers
                auto available = ModernOptimization::ModernOptimizerDispatcher::getAvailableOptimizers();
                CurcumaLogger::info("Available optimizers:");
                for (const auto& pair : available) {
                    CurcumaLogger::info_fmt("  {} - {}", pair.first, pair.second);
                }

                // Demonstrate modern features
                ModernOptimization::ModernOptimizerDispatcher::demonstrateModernFeatures(molecule);

                /*
                 * ARCHITECTURAL DECISION RECORD: Modern Optimization Dispatcher
                 *
                 * CONTEXT: Multiple optimization algorithms with different capabilities
                 * - LBFGSPP: External LBFGSpp library (working implementation)
                 * - NATIVE_LBFGS: Native L-BFGS two-loop recursion (Claude Generated)
                 * - NATIVE_DIIS: Native DIIS acceleration method (Claude Generated)
                 * - NATIVE_RFO: Native RFO eigenvector following (Claude Generated)
                 * - INTERNAL: Legacy internal LBFGS (placeholder)
                 *
                 * DECISION: Strategy pattern with OptimizerType enum selection
                 * - Type-safe method selection: no magic strings, clear enum values
                 * - Auto-selection logic: >200 atoms → LBFGSpp, otherwise native methods
                 * - Unified SimpleOptimizationResult interface for all algorithms
                 *
                 * IMPLEMENTATION CHAIN:
                 * 1. main.cpp:692 → ModernOptimizerDispatcher::optimizeStructure()
                 * 2. modern_optimizer_simple.cpp:parseOptimizerType() → enum conversion
                 * 3. optimizeWithLBFGSpp() / optimizeWithNativeLBFGS() method selection
                 * 4. Individual optimizer execution with EnergyCalculator integration
                 *
                 * RUNTIME BEHAVIOR:
                 * - method="lbfgspp" → External LBFGSpp library wrapper
                 * - method="native_lbfgs" → Curcuma's native L-BFGS implementation
                 * - method="diis" → Native DIIS acceleration (Pulay method)
                 * - method="rfo" → Native RFO eigenvector following (Banerjee method)
                 * - method="auto" → Automatic selection based on system size
                 *
                 * DEBUGGING ENTRY POINTS:
                 * - Set verbosity ≥ 2 to see algorithm selection and mathematical progress
                 * - Each optimizer logs literature citations and scientific parameters
                 * - Energy calculator integration with gradient validation and step monitoring
                 */
                // Perform optimization using modern system with documented architecture above
                auto result = ModernOptimization::ModernOptimizerDispatcher::optimizeStructure(
                    &molecule, method, &energy_calc, controller["opt"]);

                if (result.success) {
                    CurcumaLogger::success("Optimization completed successfully!");
                    CurcumaLogger::energy_abs(result.final_energy, "Final energy");
                    CurcumaLogger::info_fmt("Iterations: {}", result.iterations_performed);
                    CurcumaLogger::info_fmt("Time: {:.3f} seconds", result.optimization_time_seconds);

                    // Write optimized structure
                    std::string outfile = argv[2];
                    size_t dot_pos = outfile.find_last_of('.');
                    if (dot_pos != std::string::npos) {
                        outfile = outfile.substr(0, dot_pos);
                    }
                    outfile += ".opt.xyz";

                    std::ofstream output(outfile);
                    output << result.final_molecule.XYZString();
                    CurcumaLogger::success_fmt("Optimized structure written to: {}", outfile);

                } else {
                    CurcumaLogger::error_fmt("Optimization failed: {}", result.error_message);
                    return 1;
                }

            } catch (const std::exception& e) {
                CurcumaLogger::error_fmt("Modern optimization error: {}", e.what());
                return 1;
            }
            return 0;

        // } else if (strcmp(argv[1], "-sp") == 0) {
        //     // MOVED TO STRUCTURED DISPATCH SYSTEM - Claude Generated
        //     if (argc < 3) {
        //         std::cerr << "Please use curcuma for energy calculation as follows:\ncurcuma -sp input.xyz" << std::endl;
        //         return 0;
        //     }
        //     json sp = controller["sp"];
        //     sp["SinglePoint"] = true;
        //     controller["sp"] = sp;
        //     CurcumaOpt opt(controller, false);
        //     opt.setFileName(argv[2]);
        //     opt.start();
        //     return 0;

        } else if (strcmp(argv[1], "-block") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma to split a file with many structures (trajectories) into several smaller:\ncurcuma block input.xyz X" << std::endl;
                std::cerr << "With X the number of files to produce!" << std::endl;

                return 0;
            }
            int blocks = std::stoi(argv[3]);
            std::string outfile = std::string(argv[2]);
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();
            FileIterator file(argv[2]);
            int mols = file.MaxMolecules();
            std::multimap<double, Molecule> results;
            int block = mols / blocks;
            int index = 0;
            int i = 0;
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                mol.appendXYZFile(outfile + "_" + std::to_string(index + 1) + ".xyz");
                i++;
                if (i > block) {
                    index++;
                    i = 0;
                }
            }

            return 0;

        // } else if (strcmp(argv[1], "-md") == 0) {
        //     // MOVED TO STRUCTURED DISPATCH SYSTEM - Claude Generated
        //     [SimpleMD code moved to executeSimpleMD function]

        } else if (strcmp(argv[1], "-confsearch") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for conformational search as follows:\ncurcuma -confsearch input.xyz" << std::endl;
                return 0;
            }

            ConfSearch confsearch(controller, false);
            confsearch.setFile(argv[2]);
            confsearch.start();

        } else if (strcmp(argv[1], "-rmsdtraj") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for rmsd analysis of trajectories as follows:\ncurcuma -rmsdtraj input.xyz" << std::endl;
                std::cerr << "Additonal arguments are:" << std::endl;
                std::cerr << "-write        **** Write unique conformers!" << std::endl;
                std::cerr << "-rmsd d       **** Set rmsd threshold to d ( default = 1.0)!" << std::endl;
                std::cerr << "-fragment n   **** Set fragment to n." << std::endl;
                std::cerr << "-reference    **** Add different xyz structure as reference." << std::endl;
                std::cerr << "-second       **** Add second trajectory." << std::endl;
                std::cerr << "-heavy        **** Check only heavy atoms. Do not use with -write." << std::endl;
                RMSDTraj traj(controller, false);
                traj.printHelp();
                return 0;
            }

            RMSDTraj traj(controller, false);
            traj.setFile(argv[2]);
            traj.Initialise();
            traj.start();

        } else if (strcmp(argv[1], "-nebprep") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for geometry preparation for nudge-elastic-band calculation follows:\ncurcuma -nebprep first.xyz second.xyz" << std::endl;
                return 0;
            }
            int pt = 0;
            for (std::size_t i = 3; i < argc; ++i) {
                if (strcmp(argv[i], "-pt") == 0) {
                    if (i + 1 < argc) {
                        pt = std::stoi(argv[i + 1]);
                        ++i;
                    }
                }
            }

            Molecule mol1 = Files::LoadFile(argv[2]);
            Molecule mol2 = Files::LoadFile(argv[3]);
            auto* nebdock = new NEBDocking;
            nebdock->setStructures(mol1, mol2);
            nebdock->setProtonTransfer(pt);
            nebdock->Prepare();
            delete nebdock;

        } else if (strcmp(argv[1], "-centroid") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma for centroid calculation of user definable fragments:\ncurcuma -centroid first.xyz" << std::endl;
                return 0;
            }

            std::cout << controller << std::endl;
            std::vector<int> frag_list, atom_list;
            if (controller["centroid"].contains("addfragment"))
                frag_list = Tools::CreateList(controller["centroid"]["addfragment"].get<std::string>());

            if (controller["centroid"].contains("addatoms")) {
                if (controller["centroid"]["addatoms"].is_number())
                    atom_list.push_back(controller["centroid"]["addatoms"]);
                else
                    atom_list = Tools::CreateList(controller["centroid"]["addatoms"].get<std::string>());
                for (int i : atom_list)
                    std::cout << i << " ";
                std::cout << std::endl;
            }
            if (!frag_list.empty() && !atom_list.empty()) {
                std::cout << "Having both, fragments and atoms, added is for now mutale exclusive. Might be changed someday ...";
                exit(1);
            }
            int pt = 0, fragment = 0;
            /*
            std::vector<int> frag;
            for (std::size_t i = 3; i < argc; ++i) {
                if (strcmp(argv[i], "-fragment") == 0) {
                    if (i + 1 < argc) {
                        fragment = std::stoi(argv[i + 1]);
                        ++i;
                    }
                    // continue;
                }
                if (strcmp(argv[i], "-addfragment") == 0) {
                    bool loop = true;
                    while (i + 1 < argc && loop) {
                        StringList list = Tools::SplitString(argv[i + 1]);
                        if (list.size() > 1) {
                            for (auto element : list)
                                frag.push_back(std::stoi(element) - 1);
                        }
                        if (Tools::isInt(argv[i + 1]))
                            frag.push_back(std::stoi(argv[i + 1]) - 1);
                        else
                            loop = false;
                        ++i;
                    }
                    // continue;
                }
            }
            */
            if (!frag_list.empty()) {
                std::cout << "Using fragment of atoms :";
                for (int atom : frag_list)
                    std::cout << atom + 1 << " ";
                std::cout << std::endl;
                std::cout << "to calculate centroid!" << std::endl;

                std::ofstream result_file;
                result_file.open("centroids.dat");
                FileIterator file(argv[2]);
                while (!file.AtEnd()) {
                    Molecule mol = file.Next();
                    if (!frag_list.empty()) {
                        result_file << GeometryTools::Centroid(mol.getGeometry(frag_list)).transpose() << std::endl;
                        std::cout << mol.getGeometry(frag_list) << std::endl;
                    } else {
                        mol.GetFragments(1.2);
                        result_file << GeometryTools::Centroid(mol.getGeometryByFragment(fragment)).transpose() << std::endl;
                    }
                }
            }
            if (!atom_list.empty()) {
                std::ofstream result_file;
                result_file.open("centroids.dat");
                FileIterator file(argv[2]);
                while (!file.AtEnd()) {
                    Molecule mol = file.Next();
                    Molecule tmp;
                    // for (int atom : atom_list)
                    //     tmp.addPair(mol.Atom(atom - 1));
                    result_file << GeometryTools::Centroid(mol.getGeometry(atom_list)).transpose() << std::endl;
                }
            }

        } else if (strcmp(argv[1], "-split") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma to split supramolecular structures as follows:\ncurcuma -split molecule.xyz" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);
            int index = 1;
            std::string outfile = argv[2];
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                mol.setScaling(1.2);
                std::cout << file.MaxMolecules() << std::endl;
                if (file.MaxMolecules() <= 1)
                    mol.writeXYZFragments(outfile);
                else
                    mol.writeXYZFragments(outfile + "_M" + std::to_string(index));
                index++;
            }

        } else if (strcmp(argv[1], "-distance") == 0) {
            if (argc < 4) {
                std::cerr << "Please use curcuma to calculate distances as follows:\ncurcuma -distance molecule.xyz indexA indexB" << std::endl;
                return 0;
            }
            std::string atomsA = std::string(argv[3]);
            std::string atomsB = std::string(argv[4]);

            /*
            int indexA = 0, indexB = 0;
            try {
                indexA = std::stoi(argv[3]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate distances as follows:\ncurcuma -distance molecule.xyz indexA indexB" << std::endl;
                return 0;
            }
            try {
                indexB = std::stoi(argv[4]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate distances as follows:\ncurcuma -distance molecule.xyz indexA indexB" << std::endl;
                return 0;
            }*/
            std::vector<int> A = Tools::CreateList(atomsA);
            std::vector<int> B = Tools::CreateList(atomsB);
            /*    if(A.size() == 1 && B.size() == 1)
                {
                    int indexA = A[0], indexB = B[0];
                    FileIterator file(argv[2]);
                    outfile = argv[2];

                    while (!file.AtEnd()) {
                        Molecule mol = file.Next();
                        std::cout << ":: " << mol.CalculateDistance(indexA - 1, indexB - 1) << "::" << std::endl;
                    }
                }else
                */
            {
                for (int i : A)
                    std::cout << i << " ";
                std::cout << std::endl;

                for (int i : B)
                    std::cout << i << " ";
                std::cout << std::endl;

                FileIterator file(argv[2]);
                std::string outfile = argv[2];

                while (!file.AtEnd()) {
                    Molecule mol = file.Next();
                    Molecule tmpA, tmpB;
                    /* for (int atom : A)
                         tmpA.addPair(mol.Atom(atom - 1));
                     for (int atom : B)
                         tmpB.addPair(mol.Atom(atom - 1));*/
                    // std::cout << mol.getGeometry(A) << std::endl << mol.getGeometry(B) << std::endl;
                    auto cA = GeometryTools::Centroid(mol.getGeometry(A));
                    auto cB = GeometryTools::Centroid(mol.getGeometry(B));
                    std::cout << ":: " << sqrt((((cA[0] - cB[0]) * (cA[0] - cB[0])) + ((cA[1] - cB[1]) * (cA[1] - cB[1])) + ((cA[2] - cB[2]) * (cA[2] - cB[2])))) << "::" << std::endl;
                }
            }
        // TODO: This command is part of the legacy system and should be migrated to a modern capability handler.
        // During migration, replace Tools::CreateList with Tools::ParseStringToVector and adjust the input string format (',' to ';', ':' to '-').
        } else if (strcmp(argv[1], "-angle") == 0) {
            if (argc < 6) {
                std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -angle molecule.xyz indexA indexB indexC" << std::endl;
                return 0;
            }

            /*
            int indexA = 0, indexB = 0, indexC = 0;
            try {
                indexA = std::stoi(argv[3]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -angle molecule.xyz indexA indexB indexC" << std::endl;
                return 0;
            }
            try {
                indexB = std::stoi(argv[4]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -angle molecule.xyz indexA indexB indexC" << std::endl;
                return 0;
            }
            try {
                indexC = std::stoi(argv[5]);
            } catch (const std::invalid_argument& arg) {
                std::cerr << "Please use curcuma to calculate angles as follows:\ncurcuma -angle molecule.xyz indexA indexB indexC" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);

            printf("\n  Angle\t\tr(%u,%u)\tr(%u,%u)\tr(%u,%u)\n", indexA - 1, indexB - 1, indexA - 1, indexC - 1, indexC - 1, indexB - 1);

            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                printf(":: %8.4f\t%8.4f\t%8.4f\t%8.4f ::\n", mol.CalculateAngle(indexA - 1, indexB - 1, indexC - 1), mol.CalculateDistance(indexA - 1, indexB - 1), mol.CalculateDistance(indexA - 1, indexC - 1), mol.CalculateDistance(indexC - 1, indexB - 1));
            }
            printf("\n\n");
            */
            std::string atomsA = std::string(argv[3]);
            std::string atomsB = std::string(argv[4]);
            std::string atomsC = std::string(argv[5]);

            std::vector<int> A = Tools::CreateList(atomsA);
            std::vector<int> B = Tools::CreateList(atomsB);
            std::vector<int> C = Tools::CreateList(atomsC);

            /* if(A.size() == 1 && B.size() == 1 && C.size() == 1)
             {
                 int indexA = A[0], indexB = B[0], indexC = C[0];

                 FileIterator file(argv[2]);

                 printf("\n  Angle\t\tr(%u,%u)\tr(%u,%u)\tr(%u,%u)\n", indexA - 1, indexB - 1, indexA - 1, indexC - 1, indexC - 1, indexB - 1);

                 while (!file.AtEnd()) {
                     Molecule mol = file.Next();
                     printf(":: %8.4f\t%8.4f\t%8.4f\t%8.4f ::\n", mol.CalculateAngle(indexA - 1, indexB - 1, indexC - 1), mol.CalculateDistance(indexA - 1, indexB - 1), mol.CalculateDistance(indexA - 1, indexC - 1), mol.CalculateDistance(indexC - 1, indexB - 1));
                 }
                 printf("\n\n");
             }else
            */
            {
                for (int i : A)
                    std::cout << i << " ";
                std::cout << std::endl;

                for (int i : B)
                    std::cout << i << " ";
                std::cout << std::endl;

                for (int i : C)
                    std::cout << i << " ";
                std::cout << std::endl;
                FileIterator file(argv[2]);

                // printf("\n  Angle\t\tr(%u,%u)\tr(%u,%u)\tr(%u,%u)\n", indexA - 1, indexB - 1, indexA - 1, indexC - 1, indexC - 1, indexB - 1);

                while (!file.AtEnd()) {
                    Molecule mol = file.Next();
                    /* Molecule tmpA, tmpB, tmpC;
                     for (int atom : A)
                         tmpA.addPair(mol.Atom(atom));
                     for (int atom : B)
                         tmpB.addPair(mol.Atom(atom));
                     for (int atom : C)
                         tmpC.addPair(mol.Atom(atom));*/
                    auto cA = GeometryTools::Centroid(mol.getGeometry(A));
                    auto cB = GeometryTools::Centroid(mol.getGeometry(B));
                    auto cC = GeometryTools::Centroid(mol.getGeometry(C));
                    auto cAB = cA - cB;
                    auto cCB = cC - cB;
                    double angle = acos(DotProduct(cAB, cCB) / (sqrt(DotProduct(cAB, cAB) * DotProduct(cCB, cCB)))) * 360 / 2.0 / pi;
                    printf(":: %8.4f\t%8.4f\t%8.4f\t%8.4f ::\n",
                        angle,
                        sqrt((((cA[0] - cB[0]) * (cA[0] - cB[0])) + ((cA[1] - cB[1]) * (cA[1] - cB[1])) + ((cA[2] - cB[2]) * (cA[2] - cB[2])))),
                        sqrt((((cA[0] - cC[0]) * (cA[0] - cC[0])) + ((cA[1] - cC[1]) * (cA[1] - cC[1])) + ((cA[2] - cC[2]) * (cA[2] - cC[2])))),
                        sqrt((((cB[0] - cC[0]) * (cB[0] - cC[0])) + ((cB[1] - cC[1]) * (cB[1] - cC[1])) + ((cB[2] - cC[2]) * (cB[2] - cC[2])))));
                }
            }
        // TODO: This command is part of the legacy system and should be migrated to a modern capability handler.
        } else if (strcmp(argv[1], "-dMatrix") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma to calculate a distance matrix for a molecule as follows:\ncurcuma -dMatrix molecule.xyz [options]" << std::endl;
                std::cerr << "\nFile output options:" << std::endl;
                std::cerr << "  -exclude_bonds  Exclude bonds from distance matrix." << std::endl;
                std::cerr << "  -print_elements  Print elements in distance matrix." << std::endl;
                std::cerr << "  -print_energy    Print energy in distance matrix." << std::endl;
                std::cerr << "  -stride          Process every N-th image (default: 1)." << std::endl;
                std::cerr << "  -save_dmat       Save the distance matrix file (.dMat)." << std::endl;
                std::cerr << "  -save_pairs      Save the persistence pairs file (.pairs)." << std::endl;
                std::cerr << "  -save_pd_text    Save the persistence diagram as text (.PD)." << std::endl;
                std::cerr << "  -save_pd_image   Save the persistence diagram as an image (default: true)." << std::endl;
                std::cerr << "  -stride          Process every N-th image (default: 1)." << std::endl;
                std::cerr << "  -save_pi_text    Save the persistence image as text (.PI)." << std::endl;
                std::cerr << "  -save_pi_image   Save the persistence image as an image." << std::endl;
                std::cerr << "\nImage options:" << std::endl;
                std::cerr << "  -format <format> Image format (png, jpg, bmp, tga, default: png)." << std::endl;
                std::cerr << "  -colormap <map>  Colormap (grayscale, jet, hot, viridis, coolwarm, default: hot)." << std::endl;
                std::cerr << "  -resolution <w>x<h> Image resolution (default: 800x800)." << std::endl;
                std::cerr << "  Post-processing options:" << std::endl;
                std::cerr << "  -post_processing <type> Post-processing type (none, adaptive, ring_focused, default: none)." << std::endl;
                std::cerr << "  -temperature <f>      Temperature for adaptive scaling (default: 2.0)." << std::endl;
                std::cerr << "  -damping <f>          Damping for adaptive scaling (default: 1.5)." << std::endl;
                std::cerr << "  -preserve_structure <true|false> Preserve structure in adaptive scaling (default: true)." << std::endl;
                std::cerr << "\nPersistence diagram options:" << std::endl;
                std::cerr << "  -ripser_xmax <f>  Max x-value for persistence diagram (default: 4.0)." << std::endl;
                std::cerr << "  -ripser_xmin <f>  Min x-value for persistence diagram (default: 0.0)." << std::endl;
                std::cerr << "  -ripser_ymax <f>  Max y-value for persistence diagram (default: 4.0)." << std::endl;
                std::cerr << "  -ripser_ymin <f>  Min y-value for persistence diagram (default: 0.0)." << std::endl;
                std::cerr << "  -ripser_bins <n>  Number of bins for persistence image (default: 10)." << std::endl;
                std::cerr << "  -ripser_scaling <f> Scaling factor for persistence image (default: 0.1)." << std::endl;
                std::cerr << "  -ripser_stdx <f>  Standard deviation for x-axis in persistence image (default: 10.0)." << std::endl;
                std::cerr << "  -ripser_stdy <f>  Standard deviation for y-axis in persistence image (default: 10.0)." << std::endl;
                std::cerr << "  -ripser_ratio <f> Ratio for ripser calculation (default: 1.0)." << std::endl;
                std::cerr << "  -ripser_dimension <n> Dimension for ripser calculation (default: 2)." << std::endl;
                std::cerr << "  -ripser_epsilon <f> Epsilon for ripser calculation (default: 0.4)." << std::endl;
                std::cerr << "  -ripser_max <f>  Max value for ripser calculation (default: 0.0)." << std::endl;
                return 0;
            }
            bool exclude_bonds = false;
            bool print_elements = false;
            bool print_energy = false;
            bool exclude_hydrogen = true;
            std::string format = "png";
            EigenImageWriter::ColorMap colormap = EigenImageWriter::HOT;
            int width = 800, height = 800;
            int stride = 1;
            double min = 0;
            bool save_dmat = false, save_pairs = false, save_pd_text = false, save_pd_image = true, save_pi_text = false, save_pi_image = false, save_pd_average_image = false, save_pd_stddev_image = false, save_pi_average_image = false, save_pi_stddev_image = false;
            EigenImageWriter::PostProcessing post_processing = EigenImageWriter::NONE;
            double temperature = 2.0, damping = 1.5;
            bool preserve_structure = true;

            if (controller.contains("dMatrix")) {
                json dMatrix_opts = controller["dMatrix"];
                if (dMatrix_opts.contains("stride"))
                    stride = dMatrix_opts["stride"];
                if (dMatrix_opts.contains("exclude_bonds"))
                    exclude_bonds = dMatrix_opts["exclude_bonds"];
                if (dMatrix_opts.contains("print_elements"))
                    print_elements = dMatrix_opts["print_elements"];
                if (dMatrix_opts.contains("print_energy"))
                    print_energy = dMatrix_opts["print_energy"];
                if (dMatrix_opts.contains("save_dmat"))
                    save_dmat = dMatrix_opts["save_dmat"];
                if (dMatrix_opts.contains("save_pairs"))
                    save_pairs = dMatrix_opts["save_pairs"];
                if (dMatrix_opts.contains("save_pd_text"))
                    save_pd_text = dMatrix_opts["save_pd_text"];
                if (dMatrix_opts.contains("save_pd_image"))
                    save_pd_image = dMatrix_opts["save_pd_image"];
                if (dMatrix_opts.contains("save_pi_text"))
                    save_pi_text = dMatrix_opts["save_pi_text"];
                if (dMatrix_opts.contains("save_pi_image"))
                    save_pi_image = dMatrix_opts["save_pi_image"];
                if (dMatrix_opts.contains("save_pd_average_image"))
                    save_pd_average_image = dMatrix_opts["save_pd_average_image"];
                if (dMatrix_opts.contains("save_pd_stddev_image"))
                    save_pd_stddev_image = dMatrix_opts["save_pd_stddev_image"];
                if (dMatrix_opts.contains("save_pi_average_image"))
                    save_pi_average_image = dMatrix_opts["save_pi_average_image"];
                if (dMatrix_opts.contains("save_pi_stddev_image"))
                    save_pi_stddev_image = dMatrix_opts["save_pi_stddev_image"];
                if (dMatrix_opts.contains("format"))
                    format = dMatrix_opts["format"].get<std::string>();
                if (dMatrix_opts.contains("colormap")) {
                    std::string cm = dMatrix_opts["colormap"].get<std::string>();
                    if (cm == "jet") colormap = EigenImageWriter::JET;
                    else if (cm == "hot") colormap = EigenImageWriter::HOT;
                    else if (cm == "viridis") colormap = EigenImageWriter::VIRIDIS;
                    else if (cm == "coolwarm") colormap = EigenImageWriter::COOLWARM;
                    else if (cm == "grayscale")
                        colormap = EigenImageWriter::GRAYSCALE;
                }
                if (dMatrix_opts.contains("resolution")) {
                    std::string res = dMatrix_opts["resolution"].get<std::string>();
                    size_t pos = res.find('x');
                    if (pos != std::string::npos) {
                        width = std::stoi(res.substr(0, pos));
                        height = std::stoi(res.substr(pos + 1));
                    } else {
                        width = height = std::stoi(res);
                    }
                }
                if (dMatrix_opts.contains("post_processing")) {
                    std::string pp = dMatrix_opts["post_processing"].get<std::string>();
                    if (pp == "adaptive")
                        post_processing = EigenImageWriter::ADAPTIVE;
                    else if (pp == "ring_focused")
                        post_processing = EigenImageWriter::RING_FOCUSED;
                }
                if (dMatrix_opts.contains("temperature"))
                    temperature = dMatrix_opts["temperature"];
                if (dMatrix_opts.contains("damping"))
                    damping = dMatrix_opts["damping"];
                if (dMatrix_opts.contains("preserve_structure"))
                    preserve_structure = dMatrix_opts["preserve_structure"];
            }
            std::cout << "Excluding bonds: " << exclude_bonds << std::endl;
            std::cout << "Printing elements: " << print_elements << std::endl;
            std::cout << "Printing energy: " << print_energy << std::endl;
            FileIterator file(argv[2]);
            json dMatrix = controller["dMatrix"];
            fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "\nPlease cite the follow research report!\nTownsend, J., Micucci, C.P., Hymel, J.H. et al. Representation of molecular structures with persistent homology for machine learning applications in chemistry. Nat Commun 11, 3230 (2020). https://doi.org/10.1038/s41467-020-17035-5\n\n");

            std::string outfile = argv[2];
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();

            int index = 0;
            Eigen::MatrixXd total_pd_image_sum;
            Eigen::MatrixXd total_pi_image_sum;
            Eigen::MatrixXd total_pd_image_sq_sum;
            Eigen::MatrixXd total_pi_image_sq_sum;
            int num_images = 0;
            while (!file.AtEnd()) {
                if (stride > 1 && index % stride != 0) {
                    file.Next();
                    index++;
                    continue;
                }
                std::cout << "Processing image " << index + 1 << " of " << file.MaxMolecules() << std::endl;
                Molecule mol = file.Next();

                std::ofstream input;

                if (save_dmat) {
                    input.open(outfile + "_" + std::to_string(index) + ".dMat", std::ios::out);
                    if (print_energy)
                        input << std::setprecision(10) << mol.Energy() << std::endl;
                    input << mol.DistanceMatrixString(exclude_bonds, print_elements);
                    input.close();
                }

                auto vector = mol.LowerDistanceVector(exclude_hydrogen);

                PersistentDiagram diagram(controller["dMatrix"]);
                diagram.setDistanceMatrix(vector);
                diagram.setENScaling(mol.DeltaEN());
                {
                    auto l_pd = diagram.generatePairs();
                    if (save_pairs) {
                        input.open(outfile + "_" + std::to_string(index) + ".pairs", std::ios::out);
                        for (const auto& r : l_pd) {
                            input << r.first << " " << r.second << std::endl;
                        }
                        input.close();
                    }
                    if (save_pd_text) {
                        std::cout << "Writing Persistence diagram as " + outfile + "_" + std::to_string(index) + ".PD" << std::endl;
                        input.open(outfile + "_" + std::to_string(index) + ".PD", std::ios::out);
                        input << diagram.generateImage(l_pd);
                        input.close();
                    }
                    if (save_pd_image) {
                        Eigen::MatrixXd pd_image_matrix = diagram.generateImage(l_pd);
                        EigenImageWriter::saveMatrix(pd_image_matrix, outfile + "_" + std::to_string(index) + ".PD." + format, colormap, 90, false, width, height, post_processing, temperature, damping, preserve_structure);
                    }
                    Eigen::MatrixXd pd_image_matrix_for_acc = diagram.generateImage(l_pd);
                    if (num_images == 0) {
                        total_pd_image_sum = pd_image_matrix_for_acc;
                        total_pd_image_sq_sum = pd_image_matrix_for_acc.array().square().matrix();
                    } else {
                        total_pd_image_sum += pd_image_matrix_for_acc;
                        total_pd_image_sq_sum += pd_image_matrix_for_acc.array().square().matrix();
                    }
                }
                diagram.setDistanceMatrix(vector);
                {
                    std::cout << "Writing Persistence Image (EN scaled bond topology) as " + outfile + "_" + std::to_string(index) + ".PI" << std::endl;
                    auto l_pi = diagram.generateTriples();
                    if (save_pi_text) {
                        input.open(outfile + "_" + std::to_string(index) + ".PI", std::ios::out);
                        input << diagram.generateImage(l_pi);
                        input.close();
                    }
                    if (save_pi_image) {
                        Eigen::MatrixXd pi_image_matrix = diagram.generateImage(l_pi);
                        EigenImageWriter::saveMatrix(pi_image_matrix, outfile + "_" + std::to_string(index) + ".PI." + format, colormap, 90, false, width, height, post_processing, temperature, damping, preserve_structure);
                    }
                    Eigen::MatrixXd pi_image_matrix_for_acc = diagram.generateImage(l_pi);
                    if (num_images == 0) {
                        total_pi_image_sum = pi_image_matrix_for_acc;
                        total_pi_image_sq_sum = pi_image_matrix_for_acc.array().square().matrix();
                    } else {
                        total_pi_image_sum += pi_image_matrix_for_acc;
                        total_pi_image_sq_sum += pi_image_matrix_for_acc.array().square().matrix();
                    }
                }
                num_images++;
                index++;
            }
            if (num_images > 0) {
                if (save_pd_average_image) {
                    Eigen::MatrixXd average_pd_image = total_pd_image_sum / num_images;
                    EigenImageWriter::saveMatrix(average_pd_image, outfile + "_average.PD." + format, colormap, 90, false, width, height, post_processing, temperature, damping, preserve_structure);
                }
                if (save_pd_stddev_image) {
                    Eigen::MatrixXd average_pd_image = total_pd_image_sum / num_images;
                    Eigen::MatrixXd stddev_pd_image = ((total_pd_image_sq_sum / num_images) - average_pd_image.array().square().matrix()).cwiseSqrt();
                    EigenImageWriter::saveMatrix(stddev_pd_image, outfile + "_stddev.PD." + format, colormap, 90, false, width, height, post_processing, temperature, damping, preserve_structure);
                    Eigen::MatrixXd mean2std_matrix = average_pd_image.cwiseProduct(stddev_pd_image);
                    EigenImageWriter::saveMatrix(mean2std_matrix, outfile + "_s2n.PI." + format, colormap, 90, false, width, height, post_processing, temperature, damping, preserve_structure);
                }
                if (save_pi_average_image) {
                    Eigen::MatrixXd average_pi_image = total_pi_image_sum / num_images;
                    EigenImageWriter::saveMatrix(average_pi_image, outfile + "_average.PI." + format, colormap, 90, false, width, height, post_processing, temperature, damping, preserve_structure);
                }
                if (save_pi_stddev_image) {
                    Eigen::MatrixXd average_pi_image = total_pi_image_sum / num_images;
                    Eigen::MatrixXd stddev_pi_image = ((total_pi_image_sq_sum / num_images) - average_pi_image.array().square().matrix()).cwiseSqrt();
                    EigenImageWriter::saveMatrix(stddev_pi_image, outfile + "_stddev.PI." + format, colormap, 90, false, width, height, post_processing, temperature, damping, preserve_structure);
                    Eigen::MatrixXd mean2std_matrix = average_pi_image.cwiseProduct(stddev_pi_image);
                    EigenImageWriter::saveMatrix(mean2std_matrix, outfile + "_s2n.PI." + format, colormap, 90, false, width, height, post_processing, temperature, damping, preserve_structure);
                }
            }
        } else if (strcmp(argv[1], "-center") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma to center a structure as follows:\ncurcuma -center molecule.xyz" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);
            int index = 1;
            std::string outfile = argv[2];
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                std::cout << mol.Centroid() << std::endl;
                // mol.setGeometry(GeometryTools::TranslateGeometry(mol.getGeometry(), GeometryTools::Centroid(mol.getGeometry()), Position{ 0, 0, 0 }));
                mol.Center();
                std::cout << mol.Centroid() << std::endl;

                if (file.MaxMolecules() <= 1)
                    mol.writeXYZFragments(outfile);
                else
                    mol.writeXYZFragments(outfile + "_M" + std::to_string(index));
                index++;
            }
        } else if (strcmp(argv[1], "-reorder") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma to center a structure as follows:\ncurcuma -center molecule.xyz" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);
            int index = 1;
            std::string outfile = argv[2];
            for (int i = 0; i < 4; ++i)
                outfile.pop_back();
            outfile += ".random.xyz";
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                mol.writeXYZFile(outfile, Tools::RandomVector(0, mol.AtomCount()));
            }
        } else if (strcmp(argv[1], "-hessian") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma to analyse hessians of molecules as follow:\ncurcuma -hessian -hess_read_xyz molecule.xyz -hess_read_file hessian" << std::endl;
                std::cerr << "or" << std::endl;
                std::cerr << "Please use curcuma to analyse hessians of molecules as follow:\ncurcuma -hessian -hess_read_file hessian.json" << std::endl;
                return 0;
            }
            Molecule mol1 = Files::LoadFile(argv[2]);

            Hessian hessian(controller["hessian"]);
            hessian.setMolecule(mol1);

            hessian.start();
        } else if (strcmp(argv[1], "-qmdfffit") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma to analyse hessians of molecules as follow:\ncurcuma -hessian -hess_read_xyz molecule.xyz -hess_read_file hessian" << std::endl;
                std::cerr << "or" << std::endl;
                std::cerr << "Please use curcuma to analyse hessians of molecules as follow:\ncurcuma -hessian -hess_read_file hessian.json" << std::endl;
                return 0;
            }
            Molecule mol1 = Files::LoadFile(argv[2]);

            QMDFFFit qmdfffit(controller["qmdfffit"]);
            qmdfffit.setMolecule(mol1);
            qmdfffit.start();
        } else if (strcmp(argv[1], "-eht") == 0) {
            if (argc < 3) {
                return 0;
            }
            FileIterator file(argv[2]);
            while (!file.AtEnd()) {
                /*
                Molecule mol1 = file.Next();

                EHT eht;
                eht.InitialiseMolecule(mol1.getMolInfo());
                eht.Calculation();
                */
            }

        } else if (strcmp(argv[1], "-gyration") == 0) {
            FileIterator file(argv[2]);
            int count = 1;
            double sum = 0, sum_mass = 0, sqrt_sum = 0, sqrt_sum_mass = 0, hmass = 1;
            for (std::size_t i = 2; i < argc; ++i) {

                if (strcmp(argv[i], "-hmass") == 0) {
                    if (i + 1 < argc)
                        hmass = std::stoi(argv[i + 1]);
                }
            }
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                std::pair<double, double> gyr = mol.GyrationRadius(hmass);
                if (std::isnan(gyr.first) || std::isnan(gyr.second))
                    continue;
                sum += gyr.first;
                sum_mass += gyr.second;
                sqrt_sum += sqrt(gyr.first);
                sqrt_sum_mass += sqrt(gyr.second);
                std::cout << ":: " << gyr.first << " " << sum / static_cast<double>(count) << " " << gyr.second << " " << sum_mass / static_cast<double>(count) << " " << sqrt(gyr.first) << " " << sqrt_sum / static_cast<double>(count) << " " << sqrt(gyr.second) << " " << sqrt_sum_mass / static_cast<double>(count) << std::endl;
                count++;
            }
        } else if (strcmp(argv[1], "-dipole") == 0) {
            if (argc < 3) {
                std::cerr << "Please use curcuma to optimise the dipole of molecules as follow:\ncurcuma -dipole molecule.xyz" << std::endl;
                return 0;
            }

            FileIterator file(argv[2]);
            auto lm_basename = file.Basename();
            // TODO: Add additional arguments...
            // TODO: if -md then do first molecular dynamic to generate some conformers, then fit
            // TODO: if -scale <array of int> then calc dipole with scalingfactor and xtb2
            // TODO: if -methode <String> then change methode

            const json blob = controller["dipole"]; // declare blob as json, const why not used for now


            std::vector<Molecule> conformers;
            while (!file.AtEnd()) { // calculation and output dipole moment
                Molecule mol = file.Next(); // load Molecule
                mol.Center(false); //sets the Centroid to the origin
                EnergyCalculator interface("gfn2", blob); // set method to gfn2-xtb
                interface.setMolecule(mol.getMolInfo()); // set molecule
                interface.CalculateEnergy(false); // calc energy and Wave function
                mol.setPartialCharges(interface.Charges()); // calc partial Charges and set it to mol
                mol.setDipole(interface.Dipole() * au); //calc dipole moments and set it to mol in eA
                conformers.push_back(mol);
            }
            Molecule mol = conformers.at(0); // maybe molecule in ground state
            //Calculation of the scaling vector linear and nonlinear
            const auto linear_vector = DipoleScalingCalculation(conformers); //linear
            const auto nonlinear_vector = OptimiseDipoleScaling(conformers, linear_vector); //nonlinear
            // output scaling vector as JSON
            std::vector<double> vec_linear_scaling(linear_vector.data(), linear_vector.data() + linear_vector.rows() * linear_vector.cols());
            std::vector<double> vec_nonlinear_scaling(nonlinear_vector.data(), nonlinear_vector.data() + nonlinear_vector.rows() * nonlinear_vector.cols());
            json scaling_vector;
            scaling_vector["scaling_vector_linear"] = Tools::DoubleVector2String(vec_linear_scaling);
            scaling_vector["scaling_vector_nonlinear"] = Tools::DoubleVector2String(vec_nonlinear_scaling);
            std::ofstream out(lm_basename + "_scaling_vector.json");
            out << scaling_vector << std::endl;

            double mean_dipole_gfn2 = 0;
            double mean_dipole_nonlinear = 0;
            double mean_dipole_linear = 0;
            double r2_lin = 0;
            double r2_nlin = 0;
            double r2_lin_diffofnorm = 0;
            double r2_nlin_diffofnorm = 0;
            //output Dipole moments + Calculation of Mean Dipole
            std::ofstream file_dipole;
            file_dipole.open(lm_basename + "_dipole.out", std::ios_base::app);
            file_dipole << "linear Dipole (x y z magn.); nonlinear Dipole (x y z magn.); gfn2 Dipoles (x y z magn.)" << std::endl;
            for (const auto& conf : conformers){
                const auto dipole_lin = conf.CalculateDipoleMoment(vec_linear_scaling);
                const auto dipole_nlin = conf.CalculateDipoleMoment(vec_nonlinear_scaling);
                const auto dipole_gfn2 = conf.getDipole();
                mean_dipole_linear += dipole_lin.norm()/ conformers.size();
                mean_dipole_nonlinear += dipole_nlin.norm()/ conformers.size();
                mean_dipole_gfn2 += dipole_gfn2.norm()/ conformers.size();
                file_dipole << dipole_lin[0] << " " << dipole_lin[1] << " " << dipole_lin[2] << " " << dipole_lin.norm() << "; ";
                file_dipole << dipole_nlin[0] << " " << dipole_nlin[1] << " " << dipole_nlin[2] << " " << dipole_nlin.norm() << "; ";
                file_dipole << dipole_gfn2[0] << " " << dipole_gfn2[1] << " " << dipole_gfn2[2] << " " << dipole_gfn2.norm() << std::endl;
                const double residual = (conf.CalculateDipoleMoment(linear_vector) - conf.getDipole()).norm();
                r2_lin += residual * residual;
                const double residual_1 = (conf.CalculateDipoleMoment(nonlinear_vector) - conf.getDipole()).norm();
                r2_nlin += residual_1 * residual_1;
                const double residual_2 = conf.CalculateDipoleMoment(linear_vector).norm() - conf.getDipole().norm();
                r2_lin_diffofnorm += residual_2 * residual_2;
                const double residual_3 = conf.CalculateDipoleMoment(nonlinear_vector).norm() - conf.getDipole().norm();
                r2_nlin_diffofnorm += residual_3 * residual_3;
            };
            file_dipole.close();

            std::cout << "\nMean xtb2-Dipole: " << mean_dipole_gfn2 << " [eA], " << mean_dipole_gfn2/0.2082 << " [D]" << std::endl;
            std::cout << "Mean linear Dipole: " << mean_dipole_linear << " [eA], " << mean_dipole_linear/0.2082 << " [D]" << std::endl
                      << "Mean nonlinear Dipole: " << mean_dipole_nonlinear << " [eA], " << mean_dipole_nonlinear/0.2082 << " [D]" << std::endl
                      << std::endl;

            std::cout << "linear Scaling vector:\n"
                      << linear_vector << "\n"
                      << "nonlinear Scaling vector:\n"
                      << nonlinear_vector << "\n"
                      << std::endl;

            std::cout << "Square Sum of Residuals of Components:" << std::endl
            << "linear: " << r2_lin << std::endl
            << "nonlinear " << r2_nlin << std::endl;
            std::cout << "Square Sum of Residuals of Magnitudes" << std::endl
            << "linear: " << r2_lin_diffofnorm << std::endl
            << "nonlinear: " << r2_nlin_diffofnorm << std::endl;

        } else if (strcmp(argv[1], "-dipole_calc") == 0) {
            if (argc < 5) {
                std::cerr << "Please use curcuma to optimise the dipole of molecules as follow:\ncurcuma -dipole molecule.xyz -scaling_json scaling_vector.json" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);
            auto lm_basename = file.Basename();
            const json blob = controller["dipole_calc"]; // declare blob as json, const why not used for now
            int m_natoms;
            Molecule mol;
            while (!file.AtEnd()) { // calculation and output dipole moment
                mol = file.Next(); // load Molecule
                mol.Center(false); // sets the Centroid to the origin
                EnergyCalculator interface("gfn2", blob); // set method to gfn2-xtb
                interface.setMolecule(mol.getMolInfo()); // set molecule
                interface.CalculateEnergy(false); // calc energy and Wave function
                mol.setPartialCharges(interface.Charges()); // calc partial Charges and set it to mol
                mol.setDipole(interface.Dipole() * au); // calc dipole moments and set it to mol in eA
                m_natoms = mol.AtomCount();
            }

            std::vector<double> scaling_vector_linear = std::vector<double>(m_natoms, 1);
            std::vector<double> scaling_vector_nonlinear = std::vector<double>(m_natoms, 1);
            if (strcmp(argv[4], "none") != 0) {
                json scaling;
                std::ifstream file1(argv[4]);
                try {
                    file1 >> scaling;
                } catch ([[maybe_unused]] nlohmann::json::type_error& e) {
                    throw 404;
                } catch ([[maybe_unused]] nlohmann::json::parse_error& e) {
                    throw 404;
                }
                std::string str1, str2;
                try {
                    str1 = scaling["scaling_vector_linear"];
                    str2 = scaling["scaling_vector_nonlinear"];
                } catch ([[maybe_unused]] json::type_error& e) {
                }
                if (!str1.empty()) {
                    scaling_vector_linear = Tools::String2DoubleVec(str1, "|");
                }
                if (!str2.empty()) {
                    scaling_vector_nonlinear = Tools::String2DoubleVec(str2, "|");
                }
            }
            std::cout << "scaling_vector_linear:\n"
                      << scaling_vector_linear[0] << std::endl;
            std::cout << "scaling_vector_nonlinear:\n"
                      << scaling_vector_nonlinear[0] << std::endl;

            auto dipole_lin = mol.CalculateDipoleMoment(scaling_vector_linear);
            auto dipole_nlin = mol.CalculateDipoleMoment(scaling_vector_nonlinear);

            std::cout << "Dipole form xtb2: "
                      << mol.getDipole().norm() << " [eA] " << mol.getDipole().norm() * 4.803 << " [D] " << std::endl;
            std::cout << "Dipole form partial Charges and lin. Scaling: "
                      << dipole_lin.norm() << " [eA] " << dipole_lin.norm() * 4.803 << " [D] " << std::endl;
            std::cout << "Dipole form partial Charges and nonlin. Scaling: "
                      << dipole_nlin.norm() << " [eA] " << dipole_nlin.norm() * 4.803 << " [D] " << std::endl;

        } else if (strcmp(argv[1], "-dipole_calc") == 0) {
            if (argc < 5) {
                std::cerr << "Please use curcuma to optimise the dipole of molecules as follow:\ncurcuma -dipole molecule.xyz -scaling_json scaling_vector.json" << std::endl;
                return 0;
            }
            FileIterator file(argv[2]);
            auto lm_basename = file.Basename();
            const json blob = controller["dipole_calc"]; // declare blob as json, const why not used for now
            int m_natoms;
            Molecule mol;
            while (!file.AtEnd()) { // calculation and output dipole moment
                mol = file.Next(); // load Molecule
                mol.Center(false); //sets the Centroid to the origin
                EnergyCalculator interface("gfn2", blob); // set method to gfn2-xtb
                interface.setMolecule(mol.getMolInfo()); // set molecule
                interface.CalculateEnergy(false); // calc energy and Wave function
                mol.setPartialCharges(interface.Charges()); // calc partial Charges and set it to mol
                mol.setDipole(interface.Dipole() * au); //calc dipole moments and set it to mol in eA
                m_natoms = mol.AtomCount();
            }

            std::vector<double> scaling_vector_linear = std::vector<double>(m_natoms, 1);
            std::vector<double> scaling_vector_nonlinear = std::vector<double>(m_natoms, 1);
            if (strcmp(argv[4], "none") != 0) {
                json scaling;
                std::ifstream file1(argv[4]);
                try {
                    file1 >> scaling;
                } catch ([[maybe_unused]] nlohmann::json::type_error& e) {
                    throw 404;
                } catch ([[maybe_unused]] nlohmann::json::parse_error& e) {
                    throw 404;
                }
                std::string str1, str2;
                try {
                    str1 = scaling["scaling_vector_linear"];
                    str2 = scaling["scaling_vector_nonlinear"];
                } catch ([[maybe_unused]] json::type_error& e) {
                }
                if (!str1.empty()) {
                    scaling_vector_linear = Tools::String2DoubleVec(str1, "|");
                }
                if (!str2.empty()) {
                    scaling_vector_nonlinear = Tools::String2DoubleVec(str2, "|");
                }
            }
            std::cout << "scaling_vector_linear:\n" << scaling_vector_linear[0] << std::endl;
            std::cout << "scaling_vector_nonlinear:\n" << scaling_vector_nonlinear[0] << std::endl;

            auto dipole_lin = mol.CalculateDipoleMoment(scaling_vector_linear);
            auto dipole_nlin = mol.CalculateDipoleMoment(scaling_vector_nonlinear);

            std::cout << "Dipole form xtb2: "
                      << mol.getDipole().norm() << " [eA] " << mol.getDipole().norm()*4.803 << " [D] " << std::endl;
            std::cout << "Dipole form partial Charges and lin. Scaling: "
                      << dipole_lin.norm() << " [eA] " << dipole_lin.norm()*4.803 << " [D] " << std::endl;
            std::cout << "Dipole form partial Charges and nonlin. Scaling: "
                      << dipole_nlin.norm() << " [eA] " << dipole_nlin.norm()*4.803 << " [D] " << std::endl;

        } else if (strcmp(argv[1], "-orca") == 0) {

            if (argc < 3) {
                std::cerr << "Please use curcuma as follows:\ncurcuma -orca input" << std::endl;
                return -1;
            }

            OrcaInterface orca;
            // Eingabedatei zuweisen
            orca.setInputFile(argv[2]);

            // ORCA ausführen
            if (!orca.runOrca()) {
                return -1;
            }
            // ORCA-Ausgabe lesen
            orca.getOrcaJSON();
            orca.readOrcaJSON();

        } else if (strcmp(argv[1], "-stride") == 0) {

            if (argc < 4) {
                std::cerr << "Please use curcuma to keep only every nth structure as follows:\ncurcuma -stride trjectory.xyz 100" << std::endl;
                return 0;
            }
            int stride = std::stoi(argv[3]);

            FileIterator file(argv[2]);
            int index = 1;
            while (!file.AtEnd()) {
                Molecule mol = file.Next();
                if (index % stride == 0)
                    mol.appendXYZFile(std::string("blob.xyz"));
                index++;
            }

        // } else if (argc >= 2) {
        //     // LEGACY CATCH-ALL FILE PROCESSING - COMMENTED OUT FOR STRUCTURED DISPATCH - Claude Generated
        //     [Legacy file processing code that was interfering with structured dispatch]
        } else {
            // Claude Generated - Show help for unknown options
            std::cout << "Curcuma - Computational Chemistry Toolkit" << std::endl;
            std::cout << "==========================================" << std::endl;
            std::cout << std::endl;
            std::cout << "Usage: curcuma [option] file.xyz [parameters]" << std::endl;
            std::cout << std::endl;
            std::cout << "Optimization Options:" << std::endl;
            std::cout << "  -opt file.xyz           Traditional optimization (legacy system)" << std::endl;
            std::cout << "  -modern-opt file.xyz    Modern optimization with Strategy Pattern" << std::endl;
            std::cout << "                          Available methods: lbfgspp, lbfgs, diis, rfo, auto" << std::endl;
            std::cout << "                          Example: curcuma -modern-opt benzene.xyz lbfgspp" << std::endl;
            std::cout << "  -sp file.xyz            Single point energy calculation" << std::endl;
            std::cout << std::endl;
            std::cout << "Analysis Options:" << std::endl;
            std::cout << "  -analysis file.xyz      Unified molecular analysis (all formats, all properties)" << std::endl;
            std::cout << "  -confscan file.xyz      Conformational scanning" << std::endl;
            std::cout << "  -confsearch file.xyz    Conformational searching" << std::endl;
            std::cout << "  -confstat file.xyz      Conformational statistics" << std::endl;
            std::cout << "  -rmsd file1.xyz file2.xyz  RMSD calculation" << std::endl;
            std::cout << "  -rmsdtraj file.trj.xyz  Trajectory RMSD analysis" << std::endl;
            std::cout << "  -md file.xyz            Molecular dynamics simulation" << std::endl;
            std::cout << std::endl;
            std::cout << "Other Options:" << std::endl;
            std::cout << "  -fragment file.xyz      Fragment molecule" << std::endl;
            std::cout << "  -block file.xyz N       Split trajectory into N files" << std::endl;
            std::cout << "  -distance file.xyz i j  Calculate distance between atoms i and j" << std::endl;
            std::cout << std::endl;
            std::cout << "Modern Optimizer Features:" << std::endl;
            std::cout << "  • Strategy Pattern architecture for extensible optimization methods" << std::endl;
            std::cout << "  • Integrated CurcumaLogger with verbosity control and color output" << std::endl;
            std::cout << "  • Unit-aware output (energies in kJ/mol, distances in Å, etc.)" << std::endl;
            std::cout << "  • Comprehensive parameter validation and bounds checking" << std::endl;
            std::cout << "  • Real-time progress reporting with timing analysis" << std::endl;
            std::cout << "  • Type-safe method selection (no more magic numbers)" << std::endl;
            std::cout << std::endl;
        }
    }
#ifdef C17
#ifndef _WIN32
    std::filesystem::remove("stop");
#endif
#else
    remove("stop");
#endif

    return 0;
}
