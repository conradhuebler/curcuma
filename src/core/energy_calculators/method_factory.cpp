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
#include "qm_methods/gfn1_method.h"
#include "qm_methods/gfn2_method.h"
#include "qm_methods/gfnff_method.h"
#include "qm_methods/nddo_method.h"
#ifdef USE_CUDA
#include "qm_methods/gfnff_gpu_method.h"
#endif
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
// Method Lists
// =================================================================================

const std::vector<std::string> MethodFactory::m_ff_methods = {
    "uff", "uff-d3", "d3", "qmdff", "gfnff", "gfnff-d3"
};

const std::vector<std::string> MethodFactory::m_tblite_methods = {
    "ipea1"
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

const std::vector<std::string> MethodFactory::m_d3_methods = {
    "d3"
};

const std::vector<std::string> MethodFactory::m_d4_methods = {
    "d4"
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

// AP3 (2026-04-25): Native xTB is now the canonical gfn2 provider.
// For other providers use explicit names: "ipea1" (TBLite), "ugfn2" (Ulysses), "xtb-gfn2" (XTB).
std::unique_ptr<ComputationalMethod> MethodFactory::createGFN2(const json& config) {
    CurcumaLogger::info("GFN2: using native xTB implementation");
    return std::make_unique<GFN2Method>(config);
}

// AP3 (2026-04-25): Native xTB is now the canonical gfn1 provider.
// For other providers use explicit names: "xtb-gfn1" (XTB), "ipea1" (TBLite).
std::unique_ptr<ComputationalMethod> MethodFactory::createGFN1(const json& config) {
    CurcumaLogger::info("GFN1: using native xTB implementation");
    return std::make_unique<GFN1Method>(config);
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
    // Priority for "xtb-gfnff" method: External GFN-FF > XTB
    // Note: This is for the "xtb-gfnff" method name only.
    // The "gfnff" method is handled in create() with GPU dispatch support.
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
// Main Factory Method - Direct if/else dispatch (GCC 15 compatible, no lambdas)
// =================================================================================

/*
 * ARCHITECTURAL DECISION RECORD: Computational Method Factory
 *
 * CONTEXT: Multiple computational method providers with overlapping capabilities
 * - TBLite: Modern, fastest GFN methods (gfn1, gfn2, ipea1)
 * - XTB: Established, stable implementation (gfn1, gfn2, gfnff)
 * - Ulysses: Legacy semi-empirical methods (PM3, AM1, ugfn2)
 * - Native: Educational implementations (EHT, GFN1, GFN2, PM3, MNDO, AM1, PM6, gfnff)
 * - Force Fields: UFF, QMDFF with performance optimizations
 *
 * DECISION: Direct if/else dispatch with priority-based fallback functions
 * - No registry lambdas (avoids GCC 15 brace-init with std::function issue)
 * - Method resolution: explicit names > hierarchical priorities > error
 * - Educational focus: clear method resolution visible in debug output
 * - API preservation: maintains EnergyCalculator compatibility
 *
 * RUNTIME BEHAVIOR:
 * - "gfn2" → createGFN2() tries TBLite > Ulysses > XTB > Native fallback chain
 * - "gfn1" → createGFN1() tries TBLite > XTB > Native fallback chain
 * - "eht"  → direct EHTMethod creation (always available)
 * - "pm3"  → native PM3Method (always available, no external deps)
 * - "mndo" → native MNDOMethod (always available, no external deps)
 * - "am1"  → native AM1Method (always available, no external deps)
 * - "pm6"  → native PM6Method (always available, no external deps)
 * - "ngfn2"→ native GFN2Method directly (bypass priority chain)
 * - "uff"  → ForceFieldMethod with threading support
 * - "gfnff"→ createGFNFF() tries External > XTB > Native chain
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

    // Native methods (always available, no external dependencies)
    if (method == "eht") {
        CurcumaLogger::success("Method 'eht' resolved to native EHT");
        return std::make_unique<EHTMethod>(config);
    }

    if (method == "pm3") {
        CurcumaLogger::success("Method 'pm3' resolved to native PM3");
        return std::make_unique<NDDOMethod>(NDDOMethodType::PM3, config);
    }

    if (method == "mndo") {
        CurcumaLogger::success("Method 'mndo' resolved to native MNDO");
        return std::make_unique<NDDOMethod>(NDDOMethodType::MNDO, config);
    }

    if (method == "am1") {
        CurcumaLogger::success("Method 'am1' resolved to native AM1");
        return std::make_unique<NDDOMethod>(NDDOMethodType::AM1, config);
    }

    if (method == "pm6") {
        CurcumaLogger::success("Method 'pm6' resolved to native PM6");
        return std::make_unique<NDDOMethod>(NDDOMethodType::PM6, config);
    }

    // Direct native GFN2 access (bypasses priority chain, always uses native implementation)
    if (method == "ngfn2") {
        CurcumaLogger::success("Method 'ngfn2' resolved to native GFN2");
        return std::make_unique<GFN2Method>(config);
    }

    // Direct native GFN1 access (bypasses priority chain, always uses native implementation)
    if (method == "ngfn1") {
        CurcumaLogger::success("Method 'ngfn1' resolved to native GFN1");
        return std::make_unique<GFN1Method>(config);
    }

    // Native GFN-FF (always available, Curcuma's own implementation)
    // GPU acceleration via -gpu cuda flag
    if (method == "gfnff") {
        std::string gpu_mode = config.value("gpu", "none");
        std::transform(gpu_mode.begin(), gpu_mode.end(), gpu_mode.begin(), ::tolower);

        // Auto-detect: use GPU if compiled with CUDA
        if (gpu_mode == "auto") {
#ifdef USE_CUDA
            gpu_mode = "cuda";
#else
            gpu_mode = "none";
#endif
        }

        if (gpu_mode == "cuda") {
#ifdef USE_CUDA
            CurcumaLogger::info("GFN-FF: using GPU acceleration (CUDA)");
            return std::make_unique<GFNFFGPUComputationalMethod>("gfnff", config);
#else
            throw MethodCreationException(
                "GPU acceleration requested (--gpu cuda) but Curcuma was compiled "
                "without CUDA support. Recompile with: cmake -DUSE_CUDA=ON");
#endif
        }

        CurcumaLogger::info("GFN-FF: using CPU implementation");
        return std::make_unique<GFNFFComputationalMethod>("gfnff", config);
    }

    // Force field methods (always available)
    if (method == "uff" || method == "uff-d3" || method == "qmdff") {
        CurcumaLogger::success("Method '" + method_name + "' resolved to ForceField");
        return std::make_unique<ForceFieldMethod>(method, config);
    }

    // XTB-specific methods (explicit XTB selection)
    if (method == "xtb-gfn1" || method == "xtb-gfn2") {
        return createXTBExplicit(method, config);
    }

    // External GFN-FF (External GFNFF > XTB, no native fallback)
    if (method == "xtb-gfnff") return createGFNFF(config);

    // Ulysses methods (all 27 variants with one check)
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

    // Always available: native methods, force fields, and native xTB (gfn1/gfn2/ngfn1/ngfn2)
    available.insert(available.end(), {"eht", "pm3", "mndo", "am1", "pm6",
                                       "gfn1", "gfn2", "ngfn1", "ngfn2",
                                       "gfnff", "uff", "uff-d3", "qmdff"});

    if (hasTBLite())
        available.push_back("ipea1");

    // XTB-specific explicit methods
    if (hasXTB()) {
        available.push_back("xtb-gfn1");
        available.push_back("xtb-gfn2");
    }

    // External GFN-FF (available if external GFNFF or XTB present)
    if (hasGFNFF() || hasXTB())
        available.push_back("xtb-gfnff");

    // Ulysses methods (27 variants)
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

    // Native xTB methods (AP3: canonical providers since 2026-04-25)
    if (method_name == "gfn2" || method_name == "ngfn2") {
        info["type"] = "explicit";
        info["providers"].push_back({{"name", "Native xTB"}, {"available", true}});
        return info;
    }
    if (method_name == "gfn1" || method_name == "ngfn1") {
        info["type"] = "explicit";
        info["providers"].push_back({{"name", "Native xTB"}, {"available", true}});
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
#ifdef USE_CUDA
        info["providers"].push_back({{"name", "Native+GPU"}, {"available", true}});
        info["gpu_support"] = true;
#else
        info["gpu_support"] = false;
#endif
        return info;
    }

    // Native methods (always available)
    if (method_name == "eht" || method_name == "cgfnff" ||
        method_name == "pm3" || method_name == "mndo" ||
        method_name == "am1" || method_name == "pm6") {
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

    fmt::print("Core Methods (always available, no external dependencies):\n");
    fmt::print("  - eht:   Extended Hückel Theory\n");
    fmt::print("  - pm3:   Native PM3 semi-empirical\n");
    fmt::print("  - mndo:  Native MNDO semi-empirical\n");
    fmt::print("  - am1:   Native AM1 semi-empirical\n");
    fmt::print("  - pm6:   Native PM6 semi-empirical\n");
    fmt::print("  - gfn2:  Native GFN2-xTB (canonical, via curcuma::xtb::XTB) — alias ngfn2\n");
    fmt::print("  - gfn1:  Native GFN1-xTB (canonical, via curcuma::xtb::XTB) — alias ngfn1\n");
    fmt::print("  - ngfn2: Alias for gfn2 (backward compatibility)\n");
    fmt::print("  - ngfn1: Alias for gfn1 (backward compatibility)\n");
    fmt::print("  - gfnff: Native C++ GFN-FF implementation\n");
#ifdef USE_CUDA
    fmt::print("    (GPU acceleration: use '-gpu cuda' or '-gpu auto')\n");
#endif
    fmt::print("  - uff, uff-d3, qmdff: Force field methods\n");

    fmt::print("\nOptional External Libraries:\n");
    fmt::print("  - TBLite: {}\n", hasTBLite() ? "YES" : "NO");
    fmt::print("  - XTB: {}\n", hasXTB() ? "YES" : "NO");
    fmt::print("  - Ulysses: {}\n", hasUlysses() ? "YES" : "NO");
    fmt::print("  - External GFN-FF: {}\n", hasGFNFF() ? "YES" : "NO");
    fmt::print("  - DFT-D3: {}\n", hasD3() ? "YES" : "NO");
    fmt::print("  - DFT-D4: {}\n", hasD4() ? "YES" : "NO");
#ifdef USE_CUDA
    fmt::print("  - CUDA GPU: YES (gfnff -gpu cuda)\n");
#else
    fmt::print("  - CUDA GPU: NO\n");
#endif

    fmt::print("\nExternal Provider Methods:\n");

    // IPEA1
    fmt::print("  - ipea1: {}\n", hasTBLite() ? "TBLite" : "UNAVAILABLE");

    // GFN-FF
    fmt::print("  - gfnff: ");
    std::vector<std::string> gfnff_providers;
    if (hasGFNFF()) gfnff_providers.push_back("External GFN-FF");
    if (hasXTB()) gfnff_providers.push_back("XTB");
    gfnff_providers.push_back("Native");
#ifdef USE_CUDA
    gfnff_providers.push_back("Native+GPU");
#endif
    fmt::print("{}\n", fmt::format("{}", fmt::join(gfnff_providers, " > ")));

    fmt::print("===================================\n");
}
