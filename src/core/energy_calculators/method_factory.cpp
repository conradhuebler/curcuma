/*
 * < Method Factory Implementation >
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "src/tools/string_similarity.h"

#include "method_factory.h"
#include "src/core/curcuma_logger.h"

// Method implementations
#include "ff_methods/forcefield_method.h"
#include "ff_methods/gfnff_method.h"
#include "qm_methods/dispersion_method.h"
#include "qm_methods/eht_method.h"
#include "qm_methods/external_gfnff_method.h"
#include "qm_methods/gfnff_method.h"
#ifdef USE_TBLITE
#include "qm_methods/tblite_method.h"
#endif
#include "qm_methods/ulysses_method.h"
#include "qm_methods/xtb_method.h"

#include <iostream>
#include <algorithm>
#include <fmt/format.h>

using namespace std;

// =================================================================================
// Ulysses Method List (9 base methods x 3 correction modes = 27)
// =================================================================================

const std::vector<std::string> MethodFactory::m_ff_methods = {
    "uff", "uff-d3", "d3", "qmdff", "gfnff", "gfnff-d3"
};

const std::vector<std::string> MethodFactory::m_tblite_methods = {
    "ipea1", "gfn1", "gfn2"
};

const std::vector<std::string> MethodFactory::m_xtb_methods = {
    "xtb-gfnff", "xtb-gfn1", "xtb-gfn2"
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

bool MethodFactory::hasGFNFF() {
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

// Claude Generated: Deduplicated compilation flag check
bool MethodFactory::checkCompilationFlag(const std::string& flag) {
    if (flag == "USE_TBLITE") return hasTBLite();
    if (flag == "USE_XTB") return hasXTB();
    if (flag == "USE_ULYSSES") return hasUlysses();
    if (flag == "USE_GFNFF") return hasGFNFF();
    if (flag == "USE_D3") return hasD3();
    if (flag == "USE_D4") return hasD4();
    return false;
}

// =================================================================================
// Method Classification Helpers
// =================================================================================

bool MethodFactory::isUlyssesMethod(const std::string& method) {
    return std::find(m_ulysses_methods.begin(), m_ulysses_methods.end(), method)
        != m_ulysses_methods.end();
}

// =================================================================================
// Priority-Based Method Creation (multiple providers with fallback chains)
// =================================================================================

// Claude Generated: Direct implementation replacing registry lambdas
std::unique_ptr<ComputationalMethod> MethodFactory::createGFN2(const json& config) {
    // Priority: TBLite > Ulysses > XTB
#ifdef USE_TBLITE
    CurcumaLogger::info("GFN2: trying TBLite (priority 1)");
    try {
        auto method = std::make_unique<TBLiteMethod>("gfn2", config);
        CurcumaLogger::success("GFN2 resolved to TBLite");
        return method;
    } catch (const std::exception& e) {
        CurcumaLogger::warn("TBLite failed: " + std::string(e.what()));
    }
#endif
    if (hasUlysses()) {
        CurcumaLogger::info("GFN2: trying Ulysses (fallback)");
        try {
            auto method = std::make_unique<UlyssesMethod>("ugfn2", config);
            if (method) {
                CurcumaLogger::success("GFN2 resolved to Ulysses");
                return method;
            }
        } catch (const std::exception& e) {
            CurcumaLogger::warn("Ulysses failed: " + std::string(e.what()));
        }
    }
    if (hasXTB()) {
        CurcumaLogger::info("GFN2: trying XTB (fallback)");
        try {
            auto method = std::make_unique<XTBMethod>("gfn2", config);
            if (method) {
                CurcumaLogger::success("GFN2 resolved to XTB");
                return method;
            }
        } catch (const std::exception& e) {
            CurcumaLogger::warn("XTB failed: " + std::string(e.what()));
        }
    }
    throw MethodCreationException(
        "Method 'gfn2' not available - no compiled providers (need TBLite, XTB, or Ulysses)");
}

std::unique_ptr<ComputationalMethod> MethodFactory::createGFN1(const json& config) {
    // Priority: TBLite > XTB
#ifdef USE_TBLITE
    CurcumaLogger::info("GFN1: trying TBLite (priority 1)");
    try {
        auto method = std::make_unique<TBLiteMethod>("gfn1", config);
        CurcumaLogger::success("GFN1 resolved to TBLite");
        return method;
    } catch (const std::exception& e) {
        CurcumaLogger::warn("TBLite failed: " + std::string(e.what()));
    }
#endif
    if (hasXTB()) {
        CurcumaLogger::info("GFN1: trying XTB (fallback)");
        try {
            auto method = std::make_unique<XTBMethod>("xtb-gfn1", config);
            if (method) {
                CurcumaLogger::success("GFN1 resolved to XTB");
                return method;
            }
        } catch (const std::exception& e) {
            CurcumaLogger::warn("XTB failed: " + std::string(e.what()));
        }
    }
    throw MethodCreationException(
        "Method 'gfn1' not available - no compiled providers (need TBLite or XTB)");
}

std::unique_ptr<ComputationalMethod> MethodFactory::createIPEA1(const json& config) {
    // TBLite only
#ifdef USE_TBLITE
    CurcumaLogger::info("IPEA1: using TBLite");
    return std::make_unique<TBLiteMethod>("ipea1", config);
#else
    throw MethodCreationException(
        "Method 'ipea1' requires TBLite which was not compiled (cmake .. -DUSE_TBLITE=ON)");
#endif
}

std::unique_ptr<ComputationalMethod> MethodFactory::createGFNFF(const json& config) {
    // Priority: External GFN-FF > XTB > Native cgfnff
#ifdef USE_GFNFF
    CurcumaLogger::info("GFN-FF: trying External GFN-FF (priority 1)");
    try {
        auto method = std::make_unique<ExternalGFNFFMethod>(config);
        if (method) {
            CurcumaLogger::success("GFN-FF resolved to External GFN-FF");
            return method;
        }
    } catch (const std::exception& e) {
        CurcumaLogger::warn("External GFN-FF failed: " + std::string(e.what()));
    }
#endif
    if (hasXTB()) {
        CurcumaLogger::info("GFN-FF: trying XTB (fallback)");
        try {
            auto method = std::make_unique<XTBMethod>("gfnff", config);
            if (method) {
                CurcumaLogger::success("GFN-FF resolved to XTB");
                return method;
            }
        } catch (const std::exception& e) {
            CurcumaLogger::warn("XTB failed: " + std::string(e.what()));
        }
    }
    CurcumaLogger::info("GFN-FF: using native implementation");
    return std::make_unique<GFNFFComputationalMethod>("gfnff", config);
}

// =================================================================================
// Explicit Method Creation (single provider)
// =================================================================================

std::unique_ptr<ComputationalMethod> MethodFactory::createXTBExplicit(const std::string& method, const json& config) {
    if (!hasXTB()) {
        CurcumaLogger::error_fmt("Method '{}' requires XTB which was not compiled in this build", method);
        CurcumaLogger::error("To enable, reconfigure CMake: cmake .. -DUSE_XTB=ON");
        throw MethodCreationException(fmt::format(
            "Method '{}' unavailable - XTB not compiled (flag: USE_XTB)", method));
    }
    return std::make_unique<XTBMethod>(method, config);
}

std::unique_ptr<ComputationalMethod> MethodFactory::createUlyssesExplicit(const std::string& method, const json& config) {
    if (!hasUlysses()) {
        CurcumaLogger::error_fmt("Method '{}' requires Ulysses which was not compiled in this build", method);
        CurcumaLogger::error("To enable, reconfigure CMake: cmake .. -DUSE_ULYSSES=ON");
        throw MethodCreationException(fmt::format(
            "Method '{}' unavailable - Ulysses not compiled (flag: USE_ULYSSES)", method));
    }
    return std::make_unique<UlyssesMethod>(method, config);
}

std::unique_ptr<ComputationalMethod> MethodFactory::createDFTD3(const json& config) {
    if (!hasD3()) {
        CurcumaLogger::error("Method 'd3' requires DFT-D3 which was not compiled in this build");
        CurcumaLogger::error("To enable, reconfigure CMake: cmake .. -DUSE_D3=ON");
        throw MethodCreationException("Method 'd3' unavailable - DFT-D3 not compiled (flag: USE_D3)");
    }
    return std::make_unique<DispersionMethod>("d3", config);
}

std::unique_ptr<ComputationalMethod> MethodFactory::createDFTD4(const json& config) {
    if (!hasD4()) {
        CurcumaLogger::error("Method 'd4' requires DFT-D4 which was not compiled in this build");
        CurcumaLogger::error("To enable, reconfigure CMake: cmake .. -DUSE_D4=ON");
        throw MethodCreationException("Method 'd4' unavailable - DFT-D4 not compiled (flag: USE_D4)");
    }
    return std::make_unique<DispersionMethod>("d4", config);
}

// =================================================================================
// Main Factory Method - Direct if/else dispatch
// =================================================================================

/*
 * ARCHITECTURAL DECISION RECORD: Computational Method Factory
 *
 * CONTEXT: Multiple computational method providers with overlapping capabilities
 * - TBLite: Modern, fastest GFN methods (gfn1, gfn2, ipea1)
 * - XTB: Established, stable implementation (gfn1, gfn2, gfnff)
 * - Ulysses: Legacy semi-empirical methods (PM3, AM1, ugfn2)
 * - Native: Educational implementations (EHT, gfnff)
 * - Force Fields: UFF, QMDFF with performance optimizations
 *
 * DECISION: Priority-based factory with hierarchical fallbacks
 * - Method resolution: explicit names > hierarchical priorities > error
 * - Educational focus: clear method resolution visible in debug output
 * - API preservation: maintains EnergyCalculator compatibility
 *
 * RUNTIME BEHAVIOR:
 * - "gfn2" → createGFN2() tries TBLite > Ulysses > XTB fallback chain
 * - "eht" → createEHT() → direct EHTMethod wrapper creation
 * - "uff" → createForceField() → ForceFieldMethod with threading support
 * - "gfnff" → GFNFFComputationalMethod → native GFN-FF implementation (always available)
 * - "xtb-gfnff" → ExternalGFNFF or XTBMethod("gfnff") → Fortran/XTB GFN-FF
 */
std::unique_ptr<ComputationalMethod> MethodFactory::create(const std::string& method_name, const json& config) {
    // Convert to lowercase for consistent matching
    std::string method = method_name;
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);

    CurcumaLogger::info("MethodFactory::create called");
    CurcumaLogger::param("requested_method", method_name);
    CurcumaLogger::param("normalized_method", method);

    // Priority methods (multiple providers, fallback chains)
    if (method == "gfn2") return createGFN2(config);
    if (method == "gfn1") return createGFN1(config);
    if (method == "ipea1") return createIPEA1(config);

    // Native methods (always available)
    if (method == "eht") {
        CurcumaLogger::success("Method 'eht' resolved to native EHT");
        return std::make_unique<EHTMethod>(config);
    }

    // Native GFN-FF (always available, Curcuma's own implementation)
    if (method == "gfnff") {
        CurcumaLogger::success("Method 'gfnff' resolved to native GFN-FF");
        return std::make_unique<GFNFFComputationalMethod>("gfnff", config);
    }

    // Force field methods (always available)
    if (method == "uff" || method == "uff-d3" || method == "qmdff") {
        CurcumaLogger::success("Method '" + method_name + "' resolved to ForceField");
        return std::make_unique<ForceFieldMethod>(method, config);
    }

    // XTB-specific methods
    if (method == "xtb-gfn1" || method == "xtb-gfn2") {
        return createXTBExplicit(method, config);
    }

    // External GFN-FF (External GFNFF > XTB, no native fallback)
    if (method == "xtb-gfnff") return createGFNFF(config);

    // Ulysses methods (all 27 with one check)
    if (isUlyssesMethod(method)) {
        return createUlyssesExplicit(method, config);
    }

    // Dispersion corrections
    if (method == "d3") return createDFTD3(config);
    if (method == "d4") return createDFTD4(config);

    // Unknown method - suggest alternatives
    auto available = getAvailableMethods();
    auto suggestions = StringUtils::find_closest_matches(method_name, available, 3, 3);

    if (!suggestions.empty()) {
        CurcumaLogger::error_fmt("Unknown computational method: '{}'", method_name);
        CurcumaLogger::error("Did you mean one of these?");
        for (const auto& suggestion : suggestions) {
            CurcumaLogger::error_fmt("  - {}", suggestion);
        }
    } else {
        CurcumaLogger::error_fmt("Unknown computational method: '{}'", method_name);
        CurcumaLogger::error("Run 'curcuma --methods' to see available methods");
    }

    throw MethodCreationException(fmt::format(
        "Unknown computational method: '{}'", method_name));
}

// =================================================================================
// Utility Methods
// =================================================================================

// Claude Generated: No-instantiation available methods check
std::vector<std::string> MethodFactory::getAvailableMethods() {
    std::vector<std::string> available;

    // Always available: native methods and force fields
    available.insert(available.end(), {"eht", "gfnff", "uff", "uff-d3", "qmdff"});

    // Priority methods: add if ANY provider is available
    if (hasTBLite() || hasUlysses() || hasXTB())
        available.push_back("gfn2");
    if (hasTBLite() || hasXTB())
        available.push_back("gfn1");
    if (hasTBLite())
        available.push_back("ipea1");

    // XTB-specific
    if (hasXTB()) {
        available.push_back("xtb-gfn1");
        available.push_back("xtb-gfn2");
    }

    // External GFN-FF (available if external GFNFF or XTB present)
    if (hasGFNFF() || hasXTB())
        available.push_back("xtb-gfnff");

    // Ulysses methods
    if (hasUlysses()) {
        available.insert(available.end(), m_ulysses_methods.begin(), m_ulysses_methods.end());
    }

    // Dispersion corrections
    if (hasD3()) available.push_back("d3");
    if (hasD4()) available.push_back("d4");

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

    // Priority methods
    if (method_name == "gfn2") {
        info["type"] = "priority_based";
        if (hasTBLite()) info["providers"].push_back({{"name", "TBLite"}, {"available", true}});
        if (hasUlysses()) info["providers"].push_back({{"name", "Ulysses"}, {"available", true}});
        if (hasXTB()) info["providers"].push_back({{"name", "XTB"}, {"available", true}});
        return info;
    }
    if (method_name == "gfn1") {
        info["type"] = "priority_based";
        if (hasTBLite()) info["providers"].push_back({{"name", "TBLite"}, {"available", true}});
        if (hasXTB()) info["providers"].push_back({{"name", "XTB"}, {"available", true}});
        return info;
    }
    if (method_name == "ipea1") {
        info["type"] = "priority_based";
        info["providers"].push_back({{"name", "TBLite"}, {"available", hasTBLite()}});
        return info;
    }
    if (method_name == "gfnff") {
        info["type"] = "priority_based";
        if (hasGFNFF()) info["providers"].push_back({{"name", "External GFN-FF"}, {"available", true}});
        if (hasXTB()) info["providers"].push_back({{"name", "XTB"}, {"available", true}});
        info["providers"].push_back({{"name", "Native"}, {"available", true}});
        return info;
    }

    // Explicit methods
    if (method_name == "eht" || method_name == "cgfnff") {
        info["type"] = "explicit";
        info["providers"].push_back({{"name", "Native"}, {"available", true}});
        return info;
    }
    if (method_name == "uff" || method_name == "uff-d3" || method_name == "qmdff") {
        info["type"] = "explicit";
        info["providers"].push_back({{"name", "ForceField"}, {"available", true}});
        return info;
    }
    if (method_name == "xtb-gfn1" || method_name == "xtb-gfn2") {
        info["type"] = "explicit";
        info["providers"].push_back({{"name", "XTB"}, {"available", hasXTB()}});
        return info;
    }
    if (isUlyssesMethod(method_name)) {
        info["type"] = "explicit";
        info["providers"].push_back({{"name", "Ulysses"}, {"available", hasUlysses()}});
        return info;
    }
    if (method_name == "d3") {
        info["type"] = "explicit";
        info["providers"].push_back({{"name", "DFT-D3"}, {"available", hasD3()}});
        return info;
    }
    if (method_name == "d4") {
        info["type"] = "explicit";
        info["providers"].push_back({{"name", "DFT-D4"}, {"available", hasD4()}});
        return info;
    }

    return info;
}

void MethodFactory::printAvailableMethods() {
    fmt::print("=== Available Curcuma Methods ===\n");

    fmt::print("Core Methods (always available):\n");
    fmt::print("  - EHT: Extended Hückel Theory\n");
    fmt::print("  - gfnff: Native C++ GFN-FF implementation\n");
    fmt::print("  - ForceField: UFF, QMDFF methods\n");

    fmt::print("\nOptional QM Methods:\n");
    fmt::print("  - TBLite: {}\n", hasTBLite() ? "YES" : "NO");
    fmt::print("  - XTB: {}\n", hasXTB() ? "YES" : "NO");
    fmt::print("  - Ulysses: {}\n", hasUlysses() ? "YES" : "NO");
    fmt::print("  - External GFN-FF: {}\n", hasGFNFF() ? "YES" : "NO");
    fmt::print("  - DFT-D3: {}\n", hasD3() ? "YES" : "NO");
    fmt::print("  - DFT-D4: {}\n", hasD4() ? "YES" : "NO");

    fmt::print("\nShared Methods (priority-based resolution):\n");

    // GFN2
    fmt::print("  - gfn2: ");
    std::vector<std::string> gfn2_providers;
    if (hasTBLite()) gfn2_providers.push_back("TBLite");
    if (hasUlysses()) gfn2_providers.push_back("Ulysses");
    if (hasXTB()) gfn2_providers.push_back("XTB");
    fmt::print("{}\n", gfn2_providers.empty() ? "UNAVAILABLE" : fmt::format("{}", fmt::join(gfn2_providers, " > ")));

    // GFN1
    fmt::print("  - gfn1: ");
    std::vector<std::string> gfn1_providers;
    if (hasTBLite()) gfn1_providers.push_back("TBLite");
    if (hasXTB()) gfn1_providers.push_back("XTB");
    fmt::print("{}\n", gfn1_providers.empty() ? "UNAVAILABLE" : fmt::format("{}", fmt::join(gfn1_providers, " > ")));

    // IPEA1
    fmt::print("  - ipea1: {}\n", hasTBLite() ? "TBLite" : "UNAVAILABLE");

    // GFN-FF
    fmt::print("  - gfnff: ");
    std::vector<std::string> gfnff_providers;
    if (hasGFNFF()) gfnff_providers.push_back("External GFN-FF");
    if (hasXTB()) gfnff_providers.push_back("XTB");
    gfnff_providers.push_back("Native");
    fmt::print("{}\n", fmt::format("{}", fmt::join(gfnff_providers, " > ")));

    fmt::print("===================================\n");
}
