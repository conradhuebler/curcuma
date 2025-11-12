/*
 * < Method Factory Implementation >
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
#include "src/tools/general.h"

#include "method_factory.h"
#include "src/core/curcuma_logger.h"

// Method implementations
#include "ff_methods/forcefield_method.h"
#include "qm_methods/dispersion_method.h"
#include "qm_methods/eht_method.h"
#include "qm_methods/external_gfnff_method.h"
#include "qm_methods/gfn1_method.h"
#include "qm_methods/gfn2_method.h"
#include "qm_methods/gfnff_method.h"
#include "qm_methods/tblite_method.h"
#include "qm_methods/ulysses_method.h"
#include "qm_methods/xtb_method.h"

#include <iostream>
#include <algorithm>
#include <fmt/format.h>

using namespace std;

// =================================================================================
// Method Lists (from original EnergyCalculator)
// =================================================================================

const std::vector<std::string> MethodFactory::m_ff_methods = { 
    "uff", "uff-d3", "qmdff", "cgfnff" 
};

const std::vector<std::string> MethodFactory::m_tblite_methods = { 
    "ipea1", "gfn1", "gfn2" 
};

const std::vector<std::string> MethodFactory::m_xtb_methods = { 
    "gfnff", "xtb-gfn1", "xtb-gfn2" 
};

const std::vector<std::string> MethodFactory::m_ulysses_methods = {
    // Claude Generated: Complete Ulysses methods with correction mode support
    "ugfn2", "pm6", "am1", "pm3", "mndo", "mndod", "rm1", "pm3pddg", "mndopddg", "pm3bp",
    // D3H4X correction modes
    "pm6-d3h4x", "am1-d3h4x", "pm3-d3h4x", "mndo-d3h4x", "mndod-d3h4x",
    "rm1-d3h4x", "pm3pddg-d3h4x", "mndopddg-d3h4x", "pm3bp-d3h4x",
    // D3H+ correction modes
    "pm6-d3h+", "am1-d3h+", "pm3-d3h+", "mndo-d3h+", "mndod-d3h+",
    "rm1-d3h+", "pm3pddg-d3h+", "mndopddg-d3h+", "pm3bp-d3h+"
};

const std::vector<std::string> MethodFactory::m_d3_methods = { "d3" };
const std::vector<std::string> MethodFactory::m_d4_methods = { "d4" };

// =================================================================================
// Compilation Flag Checks
// =================================================================================

bool MethodFactory::hasTBLite() {
#ifdef USE_TBLITE
    return true;
#else
    return false;
#endif
}

bool MethodFactory::hasXTB() {
#ifdef USE_XTB
    return true;
#else
    return false;
#endif
}

bool MethodFactory::hasUlysses() {
#ifdef USE_ULYSSES
    return true;
#else
    return false;
#endif
}

bool MethodFactory::hasD3() {
#ifdef USE_D3
    return true;
#else
    return false;
#endif
}

bool MethodFactory::hasGFNFF()
{
#ifdef USE_GFNFF
    return true;
#else
    return false;
#endif
}

bool MethodFactory::hasD4() {
#ifdef USE_D4
    return true;
#else
    return false;
#endif
}

// =================================================================================
// Method Registry System
// =================================================================================

const std::vector<MethodFactory::MethodPriority>& MethodFactory::getPriorityMethods() {
    static const std::vector<MethodPriority> priority_methods = {
        // GFN2: Priority TBLite > Ulysses > XTB > Native
        {
            "gfn2",
            {
                { "TBLite", [](const json& config) { return hasTBLite() ? std::make_unique<TBLiteMethod>("gfn2", config) : nullptr; } },
                { "Ulysses", [](const json& config) { return hasUlysses() ? std::make_unique<UlyssesMethod>("ugfn2", config) : nullptr; } },
                { "XTB", [](const json& config) { return hasXTB() ? std::make_unique<XTBMethod>("gfn2", config) : nullptr; } },
                { "Native", [](const json& config) { return std::make_unique<GFN2Method>(config); } }
            } },

        // GFN1: Priority TBLite > XTB > Native
        {
            "gfn1",
            {
                { "TBLite", [](const json& config) { return hasTBLite() ? std::make_unique<TBLiteMethod>("gfn1", config) : nullptr; } },
                { "XTB", [](const json& config) { return hasXTB() ? std::make_unique<XTBMethod>("xtb-gfn1", config) : nullptr; } },
                { "Native", [](const json& config) { return std::make_unique<GFN1Method>(config); } }
            } },

        // IPEA1: TBLite only
        {
            "ipea1",
            { { "TBLite", [](const json& config) { return hasTBLite() ? std::make_unique<TBLiteMethod>("ipea1", config) : nullptr; } } } },

        // GFN-FF: Priority External GFN-FF > XTB > Native cgfnff
        {
            "gfnff",
            { { "External GFN-FF", [](const json& config) { return hasGFNFF() ? std::make_unique<ExternalGFNFFMethod>(config) : nullptr; } },
                { "XTB", [](const json& config) { return hasXTB() ? std::make_unique<XTBMethod>("gfnff", config) : nullptr; } },
                { "Native", [](const json& config) { return std::make_unique<GFNFFMethod>(config); } } } }
    };

    return priority_methods;
}

const std::vector<MethodFactory::ExplicitMethod>& MethodFactory::getExplicitMethods() {
    static const std::vector<ExplicitMethod> explicit_methods = {
        // Native methods (always available)
        { "eht", [](const json& config) { return std::make_unique<EHTMethod>(config); }, "Native", false, "" },
        { "cgfnff", [](const json& config) { return std::make_unique<GFNFFMethod>(config); }, "Native", false, "" },

        // Force field methods (always available)
        { "uff", [](const json& config) { return std::make_unique<ForceFieldMethod>("uff", config); }, "ForceField", false, "" },
        { "uff-d3", [](const json& config) { return std::make_unique<ForceFieldMethod>("uff-d3", config); }, "ForceField", false, "" },
        { "qmdff", [](const json& config) { return std::make_unique<ForceFieldMethod>("qmdff", config); }, "ForceField", false, "" },

        // XTB-specific methods
        { "gfnff", [](const json& config) { return std::make_unique<XTBMethod>("gfnff", config); }, "XTB", true, "USE_XTB" },
        { "xtb-gfn1", [](const json& config) { return std::make_unique<XTBMethod>("xtb-gfn1", config); }, "XTB", true, "USE_XTB" },
        { "xtb-gfn2", [](const json& config) { return std::make_unique<XTBMethod>("xtb-gfn2", config); }, "XTB", true, "USE_XTB" },

        // Ulysses-specific methods (27 total: 9 base methods × 3 correction modes)
        // Base methods without corrections
        { "ugfn2", [](const json& config) { return std::make_unique<UlyssesMethod>("ugfn2", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm6", [](const json& config) { return std::make_unique<UlyssesMethod>("pm6", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "am1", [](const json& config) { return std::make_unique<UlyssesMethod>("am1", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndo", [](const json& config) { return std::make_unique<UlyssesMethod>("mndo", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndod", [](const json& config) { return std::make_unique<UlyssesMethod>("mndod", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "rm1", [](const json& config) { return std::make_unique<UlyssesMethod>("rm1", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3pddg", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3pddg", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndopddg", [](const json& config) { return std::make_unique<UlyssesMethod>("mndopddg", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3bp", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3bp", config); }, "Ulysses", true, "USE_ULYSSES" },

        // Methods with D3H4X corrections
        { "pm6-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("pm6-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "am1-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("am1-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndo-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("mndo-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndod-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("mndod-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "rm1-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("rm1-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3pddg-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3pddg-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndopddg-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("mndopddg-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3bp-d3h4x", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3bp-d3h4x", config); }, "Ulysses", true, "USE_ULYSSES" },

        // Methods with D3H+ corrections
        { "pm6-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("pm6-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "am1-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("am1-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndo-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("mndo-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndod-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("mndod-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "rm1-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("rm1-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3pddg-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3pddg-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "mndopddg-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("mndopddg-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },
        { "pm3bp-d3h+", [](const json& config) { return std::make_unique<UlyssesMethod>("pm3bp-d3h+", config); }, "Ulysses", true, "USE_ULYSSES" },

        // Dispersion corrections
        { "d3", [](const json& config) { return std::make_unique<DispersionMethod>("d3", config); }, "DFT-D3", true, "USE_D3" },
        { "d4", [](const json& config) { return std::make_unique<DispersionMethod>("d4", config); }, "DFT-D4", true, "USE_D4" }
    };

    return explicit_methods;
}

// =================================================================================
// Main Factory Method
// =================================================================================

/*
 * ARCHITECTURAL DECISION RECORD: Computational Method Factory
 *
 * CONTEXT: Multiple computational method providers with overlapping capabilities
 * - TBLite: Modern, fastest GFN methods (gfn1, gfn2, ipea1)
 * - XTB: Established, stable implementation (gfn1, gfn2, gfnff)
 * - Ulysses: Legacy semi-empirical methods (PM3, AM1, ugfn2)
 * - Native: Educational implementations (EHT, cgfnff)
 * - Force Fields: UFF, QMDFF with performance optimizations
 *
 * DECISION: Priority-based factory with hierarchical fallbacks
 * - Method resolution: explicit names > hierarchical priorities > error
 * - Educational focus: clear method resolution visible in debug output
 * - API preservation: maintains EnergyCalculator compatibility
 *
 * IMPLEMENTATION CHAIN:
 * 1. src/core/energycalculator.cpp:109 → MethodFactory::create()
 * 2. method_factory.cpp:198 → method name normalization and resolution
 * 3. method_factory.cpp:250+ → priority-based creation (createGFN2, createGFN1)
 * 4. qm_methods/*_method.cpp → specific method wrapper construction
 * 5. interface/*.cpp → actual computational library initialization
 *
 * RUNTIME BEHAVIOR:
 * - "gfn2" → createGFN2() tries TBLite > Ulysses > XTB fallback chain
 * - "eht" → createEHT() → direct EHTMethod wrapper creation
 * - "uff" → createForceField() → ForceFieldMethod with threading support
 * - "cgfnff" → createCGFNFF() → native GFN-FF implementation
 *
 * DEBUGGING ENTRY POINTS:
 * - Set verbosity ≥ 2 to see method resolution process and fallback attempts
 * - Failed method creation logs provider availability and compilation flags
 * - Each method wrapper reports initialization status via CurcumaLogger
 * - Check method_factory.cpp hasTBLite()/hasXTB() etc. for compilation flag status
 */
std::unique_ptr<ComputationalMethod> MethodFactory::create(const std::string& method_name, const json& config) {
    CurcumaLogger::info("MethodFactory::create called");
    CurcumaLogger::param("requested_method", method_name);
    
    // Convert to lowercase for consistent matching
    std::string method_lower = method_name;
    std::transform(method_lower.begin(), method_lower.end(), method_lower.begin(), ::tolower);
    CurcumaLogger::param("normalized_method", method_lower);
    
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::info("=== MethodFactory Resolution Process ===");
        CurcumaLogger::param("input_method", method_name);
        CurcumaLogger::param("normalized_method", method_lower);
    }
    /*std::cout << hasTBLite() << " "
              << hasXTB() << " "
              << hasUlysses() << " "
              << hasGFNFF() << " "
              << hasD3() << " "
              << hasD4() << std::endl;*/
    // 1. Check priority methods first (shared methods with fallback hierarchy)
    CurcumaLogger::info("Step 1: Checking priority methods (shared methods with fallback hierarchy)");
    for (const auto& priority_method : getPriorityMethods()) {
        if (method_lower == priority_method.method_name) {
            CurcumaLogger::success("✓ Found priority method: " + priority_method.method_name);
            CurcumaLogger::info("Priority resolution order: " + std::to_string(priority_method.priorities.size()) + " providers");

            for (size_t i = 0; i < priority_method.priorities.size(); ++i) {
                const auto& priority = priority_method.priorities[i];
                CurcumaLogger::info("→ Trying provider " + std::to_string(i + 1) + "/" + std::to_string(priority_method.priorities.size()) + ": " + priority.first);

                try {
                    auto method = priority.second(config);  // Try to create method
                    if (method != nullptr) {
                        CurcumaLogger::success("✓ SUCCESS: Method '" + method_name + "' resolved to: " + priority.first + " (priority " + std::to_string(i + 1) + ")");
                        if (CurcumaLogger::get_verbosity() >= 2 && i > 0) {
                            CurcumaLogger::warn("Note: Higher priority providers failed, using fallback #" + std::to_string(i + 1));
                        }
                        return method;
                    } else {
                        CurcumaLogger::warn("✗ Provider " + priority.first + " returned nullptr (not available/failed initialization)");
                    }
                } catch (const std::exception& e) {
                    CurcumaLogger::error("✗ Provider " + priority.first + " threw exception: " + std::string(e.what()));
                } catch (...) {
                    CurcumaLogger::error("✗ Provider " + priority.first + " threw unknown exception");
                }
            }

            CurcumaLogger::error("✗ ALL PROVIDERS FAILED for method: " + priority_method.method_name);
            // If we get here, no provider was available
            throw MethodCreationException(fmt::format(
                "Method '{}' not available - no compiled providers (need TBLite, XTB, or Ulysses)", 
                method_name));
        }
    }
    
    // 2. Check explicit methods (single provider)
    CurcumaLogger::info("Step 2: Checking explicit methods (single provider)");
    for (const auto& explicit_method : getExplicitMethods()) {
        if (method_lower == explicit_method.method_name) {
            CurcumaLogger::info("Found explicit method: " + explicit_method.method_name);
            CurcumaLogger::param("provider", explicit_method.provider_name);
            CurcumaLogger::param("requires_compilation", explicit_method.requires_compilation_flag);

            // Check compilation requirement
            if (explicit_method.requires_compilation_flag) {
                bool available = false;
                if (explicit_method.compilation_flag == "USE_TBLITE") available = hasTBLite();
                else if (explicit_method.compilation_flag == "USE_XTB") available = hasXTB();
                else if (explicit_method.compilation_flag == "USE_ULYSSES") available = hasUlysses();
                else if (explicit_method.compilation_flag == "USE_GFNFF")
                    available = hasGFNFF();
                else if (explicit_method.compilation_flag == "USE_D3") available = hasD3();
                else if (explicit_method.compilation_flag == "USE_D4") available = hasD4();
                
                CurcumaLogger::param("compilation_flag", explicit_method.compilation_flag);
                CurcumaLogger::param("flag_available", available);
                
                if (!available) {
                    CurcumaLogger::error("Method not available - compilation flag not set");
                    throw MethodCreationException(fmt::format(
                        "Method '{}' not available - {} not compiled", 
                        method_name, explicit_method.provider_name));
                }
            }
            
            CurcumaLogger::info("Creating method using explicit mapping...");
            auto method = explicit_method.creator(config);
            if (method != nullptr) {
                CurcumaLogger::success("Method '" + method_name + "' resolved to: " + explicit_method.provider_name + " (explicit mapping)");
                return method;
            } else {
                CurcumaLogger::error("Failed to create method using explicit mapping");
                throw MethodCreationException(fmt::format(
                    "Failed to create method '{}' using {}", 
                    method_name, explicit_method.provider_name));
            }
        }
    }
    
    // 3. Fallback: Check if method matches any pattern in method lists (library-specific)
    // This handles cases where method names partially match known patterns
    
    // Check force field methods
    if (matchesMethodList(method_lower, m_ff_methods)) {
        CurcumaLogger::success("Method '" + method_name + "' resolved to: ForceField (pattern match)");
        return std::make_unique<ForceFieldMethod>(method_lower, config);
    }
    
    // Check TBLite methods
    if (matchesMethodList(method_lower, m_tblite_methods) && hasTBLite()) {
        CurcumaLogger::success("Method '" + method_name + "' resolved to: TBLite (pattern match)");
        return std::make_unique<TBLiteMethod>(method_lower, config);
    }
    
    // Check XTB methods  
    if (matchesMethodList(method_lower, m_xtb_methods) && hasXTB()) {
        CurcumaLogger::success("Method '" + method_name + "' resolved to: XTB (pattern match)");
        return std::make_unique<XTBMethod>(method_lower, config);
    }
    
    // Check Ulysses methods
    if (matchesMethodList(method_lower, m_ulysses_methods) && hasUlysses()) {
        CurcumaLogger::success("Method '" + method_name + "' resolved to: Ulysses (pattern match)");
        return std::make_unique<UlyssesMethod>(method_lower, config);
    }
    
    // 4. Default fallback to ForceField (as in original EnergyCalculator)
    CurcumaLogger::warn("Unknown method: '" + method_name + "', using default ForceField");
    return std::make_unique<ForceFieldMethod>("uff", config);
}

// =================================================================================
// Utility Methods
// =================================================================================

std::vector<std::string> MethodFactory::getAvailableMethods() {
    std::vector<std::string> available;
    
    // Add priority methods that have at least one available provider
    for (const auto& priority_method : getPriorityMethods()) {
        for (const auto& priority : priority_method.priorities) {
            // Test if any provider can create this method (without actually creating it)
            try {
                auto test_method = priority.second(json{});
                if (test_method != nullptr) {
                    available.push_back(priority_method.method_name);
                    break; // Found one provider, move to next method
                }
            } catch (...) {
                // Provider not available, try next
            }
        }
    }
    
    // Add explicit methods that are available
    for (const auto& explicit_method : getExplicitMethods()) {
        bool method_available = true;
        if (explicit_method.requires_compilation_flag) {
            if (explicit_method.compilation_flag == "USE_TBLITE") method_available = hasTBLite();
            else if (explicit_method.compilation_flag == "USE_XTB") method_available = hasXTB();
            else if (explicit_method.compilation_flag == "USE_ULYSSES") method_available = hasUlysses();
            else if (explicit_method.compilation_flag == "USE_GFNFF")
                method_available = hasGFNFF();
            else if (explicit_method.compilation_flag == "USE_D3") method_available = hasD3();
            else if (explicit_method.compilation_flag == "USE_D4") method_available = hasD4();
        }
        
        if (method_available) {
            available.push_back(explicit_method.method_name);
        }
    }
    
    return available;
}

bool MethodFactory::isMethodAvailable(const std::string& method_name) {
    auto available = getAvailableMethods();
    return std::find(available.begin(), available.end(), method_name) != available.end();
}

json MethodFactory::getMethodInfo(const std::string& method_name) {
    json info;
    info["method"] = method_name;
    info["available"] = isMethodAvailable(method_name);
    info["type"] = "unknown";
    info["providers"] = json::array();
    
    // Check priority methods
    for (const auto& priority_method : getPriorityMethods()) {
        if (method_name == priority_method.method_name) {
            info["type"] = "priority_based";
            for (const auto& priority : priority_method.priorities) {
                json provider_info;
                provider_info["name"] = priority.first;
                provider_info["available"] = (priority.second(json{}) != nullptr);
                info["providers"].push_back(provider_info);
            }
            return info;
        }
    }
    
    // Check explicit methods
    for (const auto& explicit_method : getExplicitMethods()) {
        if (method_name == explicit_method.method_name) {
            info["type"] = "explicit";
            json provider_info;
            provider_info["name"] = explicit_method.provider_name;
            provider_info["available"] = !explicit_method.requires_compilation_flag || (explicit_method.compilation_flag == "USE_TBLITE" && hasTBLite()) || (explicit_method.compilation_flag == "USE_XTB" && hasXTB()) || (explicit_method.compilation_flag == "USE_ULYSSES" && hasUlysses()) || (explicit_method.compilation_flag == "USE_GFNFF" && hasGFNFF()) || (explicit_method.compilation_flag == "USE_D3" && hasD3()) || (explicit_method.compilation_flag == "USE_D4" && hasD4());
            info["providers"].push_back(provider_info);
            return info;
        }
    }
    
    return info;
}

void MethodFactory::printAvailableMethods() {
    fmt::print("=== Available Curcuma Methods ===\n");
    
    fmt::print("Core Methods (always available):\n");
    fmt::print("  - EHT: Extended Hückel Theory\n");
    fmt::print("  - cgfnff: Native GFN-FF implementation\n");
    fmt::print("  - ForceField: UFF, QMDFF methods\n");
    
    fmt::print("\nOptional QM Methods:\n");
    fmt::print("  - TBLite: {}\n", hasTBLite() ? "YES" : "NO");
    fmt::print("  - XTB: {}\n", hasXTB() ? "YES" : "NO");
    fmt::print("  - Ulysses: {}\n", hasUlysses() ? "YES" : "NO");
    fmt::print("  - External GFN-FF: {}\n", hasGFNFF() ? "YES" : "NO");
    fmt::print("  - DFT-D3: {}\n", hasD3() ? "YES" : "NO");
    fmt::print("  - DFT-D4: {}\n", hasD4() ? "YES" : "NO");
    
    fmt::print("\nShared Methods (priority-based resolution):\n");
    for (const auto& priority_method : getPriorityMethods()) {
        fmt::print("  - {}: ", priority_method.method_name);
        std::vector<std::string> available_providers;
        for (const auto& priority : priority_method.priorities) {
            if (priority.second(json{}) != nullptr) {
                available_providers.push_back(priority.first);
            }
        }
        if (!available_providers.empty()) {
            fmt::print("{} ", fmt::join(available_providers, " > "));
        } else {
            fmt::print("UNAVAILABLE");
        }
        fmt::print("\n");
    }
    
    fmt::print("===================================\n");
}

bool MethodFactory::matchesMethodList(const std::string& method, const std::vector<std::string>& method_list) {
    for (const auto& list_method : method_list) {
        if (method.find(list_method) != std::string::npos) {
            return true;
        }
    }
    return false;
}

