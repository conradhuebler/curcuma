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
#include "qm_methods/native_xtb_method.h"
#include "qm_methods/gfnff_method.h"
#include "qm_methods/nddo_method.h"
// CUDA GPU methods are loaded at runtime from libcurcuma_cuda.so via gpu_plugin (so the
// cuBLAS/cuSOLVER runtime never touches the CPU startup path). No CUDA headers here.
#include "gpu_plugin.h"
#ifdef USE_ROCM
#include "qm_methods/gfnff_hip_method.h"
#endif
#if defined(USE_ROCM)
#include "qm_methods/xtb_hip_method.h"
#endif
#if defined(USE_VULKAN)
#include "qm_methods/xtb_vulkan_method.h"
#endif
#ifdef USE_TBLITE
#include "qm_methods/tblite_method.h"
#endif
#include "qm_methods/ulysses_method.h"
#include "qm_methods/xtb_method.h"
#include "qm_methods/orca_method.h"

#include <iostream>
#include <algorithm>
#include <fmt/format.h>

using namespace std;

// =================================================================================
// Method Lists
// =================================================================================

const std::vector<std::string> MethodFactory::m_ff_methods = {
    "uff", "uff-d3", "d3", "qmdff", "gfnff", "gfnff-fast"
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
// Claude Generated (2026-06): resolve the -gpu mode for the native xTB path across all
// compiled GPU backends. Returns the backend name ("cuda"/"rocm"/"vulkan") only when it
// was requested AND this build supports it; "auto" picks the first compiled backend
// (priority cuda > rocm > vulkan); otherwise warns and returns "none" (CPU). The gfnff
// dispatch in create() follows the same scheme.
static std::string resolveNativeXtbGpuMode(const json& config, const char* label) {
    std::string gpu_mode = config.value("gpu", std::string("none"));
    std::transform(gpu_mode.begin(), gpu_mode.end(), gpu_mode.begin(), ::tolower);
    if (gpu_mode.empty() || gpu_mode == "none" || gpu_mode == "cpu")
        return "none";

    // Compile-time availability of each native-xTB GPU backend.
#if defined(USE_CUDA)
    constexpr bool has_cuda = true;
#else
    constexpr bool has_cuda = false;
#endif
#if defined(USE_ROCM)
    constexpr bool has_rocm = true;
#else
    constexpr bool has_rocm = false;
#endif
#if defined(USE_VULKAN)
    constexpr bool has_vulkan = true;
#else
    constexpr bool has_vulkan = false;
#endif

    if (gpu_mode == "auto") {
        if (has_cuda)   return "cuda";
        if (has_rocm)   return "rocm";
        if (has_vulkan) return "vulkan";
        return "none";  // no GPU backend compiled: silent CPU (auto = best-effort)
    }
    if (gpu_mode == "cuda") {
        if (has_cuda) return "cuda";
    } else if (gpu_mode == "rocm") {
        if (has_rocm) return "rocm";
    } else if (gpu_mode == "vulkan") {
        if (has_vulkan) return "vulkan";
    } else {
        CurcumaLogger::warn(std::string(label) + ": unknown -gpu value '" + gpu_mode
            + "' (use cuda|rocm|vulkan|auto|none). Falling back to CPU.");
        return "none";
    }
    // A specific backend was requested that this build does not provide.
    std::string upper = gpu_mode;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
    CurcumaLogger::warn(std::string(label) + ": GPU acceleration requested (-gpu " + gpu_mode
        + ") but native xTB " + gpu_mode + " support was not compiled. Falling back to CPU.");
    CurcumaLogger::warn("To enable it, recompile with: cmake -DUSE_" + upper
        + "=ON -DUSE_" + upper + "_XTB=ON");
    return "none";
}

std::unique_ptr<ComputationalMethod> MethodFactory::createGFN2(const json& config) {
    const std::string gpu = resolveNativeXtbGpuMode(config, "GFN2");
    (void)gpu;  // only read inside the GPU-backend #ifdefs below
#if defined(USE_CUDA)
    if (gpu == "cuda") {
        CurcumaLogger::info("GFN2: using native xTB on GPU (CUDA)");
        if (auto m = gpu_plugin::createNativeXtb("cuda",
                static_cast<int>(curcuma::xtb::MethodType::GFN2), config))
            return m;   // else the plugin was unavailable -> fall through to CPU below
    }
#endif
#if defined(USE_ROCM)
    if (gpu == "rocm") {
        CurcumaLogger::info("GFN2: using native xTB on GPU (ROCm/HIP)");
        return std::make_unique<XtbHipComputationalMethod>(curcuma::xtb::MethodType::GFN2, config);
    }
#endif
#if defined(USE_VULKAN)
    if (gpu == "vulkan") {
        CurcumaLogger::info("GFN2: using native xTB on GPU (Vulkan)");
        return std::make_unique<XtbVulkanComputationalMethod>(curcuma::xtb::MethodType::GFN2, config);
    }
#endif
    CurcumaLogger::info("GFN2: using native xTB implementation");
    return std::make_unique<NativeXtbMethod>(curcuma::xtb::MethodType::GFN2, config);
}

// AP3 (2026-04-25): Native xTB is now the canonical gfn1 provider.
// For other providers use explicit names: "xtb-gfn1" (XTB), "ipea1" (TBLite).
std::unique_ptr<ComputationalMethod> MethodFactory::createGFN1(const json& config) {
    const std::string gpu = resolveNativeXtbGpuMode(config, "GFN1");
    (void)gpu;  // only read inside the GPU-backend #ifdefs below
#if defined(USE_CUDA)
    if (gpu == "cuda") {
        CurcumaLogger::info("GFN1: using native xTB on GPU (CUDA)");
        if (auto m = gpu_plugin::createNativeXtb("cuda",
                static_cast<int>(curcuma::xtb::MethodType::GFN1), config))
            return m;   // else the plugin was unavailable -> fall through to CPU below
    }
#endif
#if defined(USE_ROCM)
    if (gpu == "rocm") {
        CurcumaLogger::info("GFN1: using native xTB on GPU (ROCm/HIP)");
        return std::make_unique<XtbHipComputationalMethod>(curcuma::xtb::MethodType::GFN1, config);
    }
#endif
#if defined(USE_VULKAN)
    if (gpu == "vulkan") {
        CurcumaLogger::info("GFN1: using native xTB on GPU (Vulkan)");
        return std::make_unique<XtbVulkanComputationalMethod>(curcuma::xtb::MethodType::GFN1, config);
    }
#endif
    CurcumaLogger::info("GFN1: using native xTB implementation");
    return std::make_unique<NativeXtbMethod>(curcuma::xtb::MethodType::GFN1, config);
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
// ORCA Method Creation (runtime availability check)
// =================================================================================

std::unique_ptr<ComputationalMethod> MethodFactory::createOrca(const std::string& method, const json& config) {
    if (!hasOrca()) {
        CurcumaLogger::error("Method '" + method + "' requires ORCA which was not found in PATH");
        CurcumaLogger::error("Please install ORCA and ensure the 'orca' executable is accessible.");
        throw MethodCreationException(
            "Method '" + method + "' unavailable - ORCA executable not found in PATH");
    }
    CurcumaLogger::success("Method '" + method + "' resolved to ORCA");
    return std::make_unique<OrcaMethod>(method, config);
}

bool MethodFactory::hasOrca() {
    return OrcaMethod::isAvailable();
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
 * - "gfn1"/"gfn2" → native curcuma xTB (NativeXtbMethod), the canonical provider
 * - "xtb-gfn1"/"xtb-gfn2" → external GFN: TBLite > external XTB binary (like xtb-gfnff)
 * - "tblite-gfn1"/"tblite-gfn2" → TBLite explicitly; "ipea1" → TBLite iPEA1
 * - "eht"  → direct EHTMethod creation (always available)
 * - "pm3"  → native PM3Method (always available, no external deps)
 * - "mndo" → native MNDOMethod (always available, no external deps)
 * - "am1"  → native AM1Method (always available, no external deps)
 * - "pm6"  → native PM6Method (always available, no external deps)
 * - "uff"  → ForceFieldMethod with threading support
 * - "gfnff"→ createGFNFF() tries External > XTB > Native chain
 */
// =================================================================================
// Method table (Claude Generated, Sep 2026)
//
// One row per method family. create(), getAvailableMethods(), getMethodInfo() and
// printAvailableMethods() are all driven by this table, so adding a method means
// adding ONE row here (plus the ComputationalMethod subclass and its PARAM block).
// Providers are listed in priority order; a method is "available" when at least one
// provider is compiled in / reachable.
// =================================================================================

namespace {

// Native GFN-FF with the optional GPU back-ends (-gpu cuda|rocm|vulkan|auto).
std::unique_ptr<ComputationalMethod> createNativeGfnff(const std::string& method, const json& config)
{
    json gfnff_config = config;
    if (method == "gfnff-fast") {
        gfnff_config["static_charges"] = true;
        gfnff_config["static_cn"] = true;
        CurcumaLogger::warn("Method 'gfnff-fast': NON-POLARIZING fast GFN-FF — EEQ charges and "
            "CN/D4 are frozen after the first geometry for speed. Valid for equilibrium dynamics "
            "/ relaxation of pre-equilibrated systems only; NOT for charge transfer, ionic "
            "dynamics, large conformational change, or reactions. See docs/GFNFF_FAST_WP.md.");
    }
    std::string gpu_mode = gfnff_config.value("gpu", "none");
    std::transform(gpu_mode.begin(), gpu_mode.end(), gpu_mode.begin(), ::tolower);
    if (gpu_mode == "auto") {
#if defined(USE_CUDA)
        gpu_mode = "cuda";
#elif defined(USE_ROCM)
        gpu_mode = "rocm";
#elif defined(USE_VULKAN)
        gpu_mode = "vulkan";
#else
        gpu_mode = "none";
#endif
    }
    if (gpu_mode == "cuda") {
#ifdef USE_CUDA
        CurcumaLogger::info("GFN-FF: using GPU acceleration (CUDA)");
        if (auto m = gpu_plugin::createGfnff("cuda", gfnff_config))
            return m;   // else the plugin was unavailable -> fall through to CPU below
#else
        CurcumaLogger::warn("GPU acceleration requested (-gpu cuda) but Curcuma was compiled "
            "without CUDA support. Falling back to CPU.");
        CurcumaLogger::warn("To enable CUDA, recompile with: cmake -DUSE_CUDA=ON");
#endif
    } else if (gpu_mode == "rocm") {
#ifdef USE_ROCM
        CurcumaLogger::info("GFN-FF: using GPU acceleration (ROCm/HIP)");
        return std::make_unique<GFNFFHipComputationalMethod>("gfnff", gfnff_config);
#else
        CurcumaLogger::warn("GPU acceleration requested (-gpu rocm) but Curcuma was compiled "
            "without ROCm GFN-FF support. Falling back to CPU.");
        CurcumaLogger::warn("To enable ROCm GFN-FF, recompile with: cmake -DUSE_ROCM=ON");
#endif
    } else if (gpu_mode == "vulkan") {
        CurcumaLogger::warn("GFN-FF -gpu vulkan: compute shaders not yet ported; using CPU.");
    } else if (gpu_mode != "none" && gpu_mode != "cpu" && !gpu_mode.empty()) {
        CurcumaLogger::warn("GFN-FF: unknown -gpu value '" + gpu_mode
            + "' (use cuda|rocm|vulkan|auto|none). Using CPU.");
    }
    CurcumaLogger::info("GFN-FF: using CPU implementation");
    return std::make_unique<GFNFFComputationalMethod>("gfnff", gfnff_config);
}

// External GFN1/GFN2 through TBLite (preferred) or the xtb binary.
std::unique_ptr<ComputationalMethod> createExternalGfn(const std::string& method, const json& config)
{
    const std::string gfn = (method == "xtb-gfn2") ? "gfn2" : "gfn1";
    if (MethodFactory::hasTBLite()) {
        CurcumaLogger::success("Method '" + method + "' resolved to TBLite " + gfn);
#ifdef USE_TBLITE
        return std::make_unique<TBLiteMethod>(gfn, config);
#endif
    }
    if (MethodFactory::hasXTB()) {
        CurcumaLogger::success("Method '" + method + "' resolved to external XTB " + gfn);
        return MethodFactory::createXTBExplicit(method, config);
    }
    throw MethodCreationException("Method '" + method + "' requires TBLite or XTB "
        "(cmake .. -DUSE_TBLITE=ON or -DUSE_XTB=ON)");
}

std::unique_ptr<ComputationalMethod> createTbliteExplicit(const std::string& method, const json& config)
{
#ifdef USE_TBLITE
    const std::string gfn = (method == "tblite-gfn2") ? "gfn2" : "gfn1";
    CurcumaLogger::info("TBLite " + gfn + ": using TBLite explicitly");
    return std::make_unique<TBLiteMethod>(gfn, config);
#else
    (void)config;
    throw MethodCreationException("Method '" + method
        + "' requires TBLite which was not compiled (cmake .. -DUSE_TBLITE=ON)");
#endif
}

auto always = []() { return true; };

} // namespace

const std::vector<MethodDescriptor>& MethodFactory::methodTable()
{
    static const std::vector<MethodDescriptor> table = {
        // ---- native QM (no external dependencies) ----
        { {"gfn2"}, "Quantum Methods (native)", "GFN2-xTB, native SCF (canonical; -gpu cuda|rocm|vulkan)",
          always, {{"Native xTB", always}},
          [](const std::string&, const json& c) { return createGFN2(c); } },
        { {"gfn1"}, "Quantum Methods (native)", "GFN1-xTB, native SCF (canonical; -gpu cuda|rocm|vulkan)",
          always, {{"Native xTB", always}},
          [](const std::string&, const json& c) { return createGFN1(c); } },
        { {"eht"}, "Quantum Methods (native)", "Extended Hueckel theory",
          always, {{"Native", always}},
          [](const std::string&, const json& c) -> std::unique_ptr<ComputationalMethod> { return std::make_unique<EHTMethod>(c); } },
        { {"pm3"}, "Quantum Methods (native)", "PM3 (NDDO)", always, {{"Native", always}},
          [](const std::string&, const json& c) -> std::unique_ptr<ComputationalMethod> { return std::make_unique<NDDOMethod>(NDDOMethodType::PM3, c); } },
        { {"mndo"}, "Quantum Methods (native)", "MNDO (NDDO)", always, {{"Native", always}},
          [](const std::string&, const json& c) -> std::unique_ptr<ComputationalMethod> { return std::make_unique<NDDOMethod>(NDDOMethodType::MNDO, c); } },
        { {"am1"}, "Quantum Methods (native)", "AM1 (NDDO)", always, {{"Native", always}},
          [](const std::string&, const json& c) -> std::unique_ptr<ComputationalMethod> { return std::make_unique<NDDOMethod>(NDDOMethodType::AM1, c); } },
        { {"pm6"}, "Quantum Methods (native)", "PM6 (NDDO)", always, {{"Native", always}},
          [](const std::string&, const json& c) -> std::unique_ptr<ComputationalMethod> { return std::make_unique<NDDOMethod>(NDDOMethodType::PM6, c); } },
        // ---- native force fields ----
        { {"gfnff", "gfnff-fast"}, "Force Fields (native)", "GFN-FF, native (gfnff-fast: frozen charges/CN; -gpu cuda|rocm)",
          always, {{"Native", always},
#ifdef USE_CUDA
                   {"Native+GPU", always},
#endif
          },
          [](const std::string& m, const json& c) { return createNativeGfnff(m, c); } },
        { {"uff", "uff-d3", "qmdff"}, "Force Fields (native)", "UFF / UFF-D3 / QMDFF (ForceField engine)",
          always, {{"ForceField", always}},
          [](const std::string& m, const json& c) -> std::unique_ptr<ComputationalMethod> {
              CurcumaLogger::success("Method '" + m + "' resolved to ForceField");
              return std::make_unique<ForceFieldMethod>(m, c); } },
        // ---- external QM providers ----
        { {"ipea1"}, "Quantum Methods (external providers)", "iPEA1-xTB via TBLite",
          hasTBLite, {{"TBLite", hasTBLite}},
          [](const std::string&, const json& c) { return createIPEA1(c); } },
        { {"tblite-gfn1", "tblite-gfn2"}, "Quantum Methods (external providers)", "GFN1/GFN2-xTB via TBLite explicitly",
          hasTBLite, {{"TBLite", hasTBLite}},
          [](const std::string& m, const json& c) { return createTbliteExplicit(m, c); } },
        { {"xtb-gfn1", "xtb-gfn2"}, "Quantum Methods (external providers)", "GFN1/GFN2-xTB: TBLite > xtb binary",
          []() { return hasTBLite() || hasXTB(); }, {{"TBLite", hasTBLite}, {"XTB", hasXTB}},
          [](const std::string& m, const json& c) { return createExternalGfn(m, c); } },
        { {"xtb-gfnff"}, "Quantum Methods (external providers)", "GFN-FF: external Fortran GFN-FF > xtb binary > native",
          always, {{"External GFN-FF", hasGFNFF}, {"XTB", hasXTB}, {"Native", always}},
          [](const std::string&, const json& c) { return createGFNFF(c); } },
        { m_ulysses_methods, "Quantum Methods (external providers)", "Ulysses semi-empirical methods (PM6, AM1, ..., ugfn2; -d3h4x / -d3h+ variants)",
          hasUlysses, {{"Ulysses", hasUlysses}},
          [](const std::string& m, const json& c) { return createUlyssesExplicit(m, c); } },
        // ---- dispersion-only ----
        { {"d3"}, "Dispersion Corrections", "DFT-D3 (external s-dftd3)", hasD3, {{"DFT-D3", hasD3}},
          [](const std::string&, const json& c) { return createDFTD3(c); } },
        { {"d4"}, "Dispersion Corrections", "DFT-D4 (external cpp-d4)", hasD4, {{"DFT-D4", hasD4}},
          [](const std::string&, const json& c) { return createDFTD4(c); } },
        // ---- ORCA (external process) ----
        { {"hf-3c", "b97-3c", "r2scan-3c", "pbeh-3c", "orca"}, "ORCA (external process)", "ORCA composite methods / custom input (-orca_input)",
          hasOrca, {{"ORCA", hasOrca}},
          [](const std::string& m, const json& c) { return createOrca(m, c); } },
    };
    return table;
}

const MethodDescriptor* MethodFactory::findMethod(const std::string& lower_name)
{
    for (const auto& d : methodTable())
        for (const auto& n : d.names)
            if (n == lower_name) return &d;
    return nullptr;
}

const std::vector<std::string>& MethodFactory::methodParameterScopes()
{
    // JSON sub-scopes that carry method-specific parameters (controller["gfnff"], ...).
    // Union of the registry modules of the energy methods plus the legacy dotted-flag
    // scopes; consumed by EnergyCalculator, the opt/sp driver, SimpleMD and ConfSearch.
    static const std::vector<std::string> scopes = {
        "gfnff", "eeq_solver", "gfnff_external", "forcefield", "uff", "qmdff",
        "xtb", "tblite", "ulysses", "eht", "orca",
        "d3", "d4", "dftd3", "dftd4", "d3param", "d4param"
    };
    return scopes;
}

std::unique_ptr<ComputationalMethod> MethodFactory::create(const std::string& method_name, const json& config) {
    std::string method = method_name;
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);

    CurcumaLogger::info("MethodFactory::create called");
    CurcumaLogger::param("requested_method", method_name);
    CurcumaLogger::param("normalized_method", method);

    if (const MethodDescriptor* d = findMethod(method)) {
        if (!d->available()) {
            std::string providers;
            for (const auto& p : d->providers) providers += (providers.empty() ? "" : ", ") + p.first;
            throw MethodCreationException(fmt::format(
                "Method '{}' is not available in this build (needs one of: {}). "
                "Run 'curcuma --methods' for the compiled-in providers.", method, providers));
        }
        return d->create(method, config);
    }

    auto available = getAvailableMethods();
    auto suggestions = StringUtils::find_closest_matches(method_name, available, 3, 3);
    CurcumaLogger::error_fmt("Unknown computational method: '{}'", method_name);
    if (!suggestions.empty()) {
        CurcumaLogger::error("Did you mean one of these?");
        for (const auto& suggestion : suggestions)
            CurcumaLogger::error_fmt("  - {}", suggestion);
    } else {
        CurcumaLogger::error("Run 'curcuma --methods' to see available methods");
    }
    throw MethodCreationException(fmt::format("Unknown computational method: '{}'", method_name));
}

std::vector<std::string> MethodFactory::getAvailableMethods() {
    std::vector<std::string> available;
    for (const auto& d : methodTable())
        if (d.available())
            available.insert(available.end(), d.names.begin(), d.names.end());
    return available;
}

bool MethodFactory::isMethodAvailable(const std::string& method_name) {
    std::string m = method_name;
    std::transform(m.begin(), m.end(), m.begin(), ::tolower);
    const MethodDescriptor* d = findMethod(m);
    return d && d->available();
}

json MethodFactory::getMethodInfo(const std::string& method_name) {
    json info;
    info["method"] = method_name;
    info["available"] = isMethodAvailable(method_name);
    info["type"] = "unknown";
    info["providers"] = json::array();
    std::string m = method_name;
    std::transform(m.begin(), m.end(), m.begin(), ::tolower);
    if (const MethodDescriptor* d = findMethod(m)) {
        info["type"] = (d->providers.size() > 1) ? "priority_based" : "explicit";
        info["family"] = d->family;
        info["description"] = d->description;
        for (const auto& p : d->providers)
            info["providers"].push_back({{"name", p.first}, {"available", p.second()}});
        if (m == "gfnff" || m == "gfnff-fast") {
#ifdef USE_CUDA
            info["gpu_support"] = true;
#else
            info["gpu_support"] = false;
#endif
        }
    }
    return info;
}

void MethodFactory::printAvailableMethods() {
    fmt::print("=== Available Curcuma Methods ===\n");
    std::string current_family;
    for (const auto& d : methodTable()) {
        if (d.family != current_family) {
            current_family = d.family;
            fmt::print("\n{}:\n", current_family);
        }
        std::string names;
        if (d.names.size() <= 5) {
            for (const auto& n : d.names) names += (names.empty() ? "" : " / ") + n;
        } else {
            names = d.names.front() + " (+" + std::to_string(d.names.size() - 1) + " more)";
        }
        std::string providers;
        for (const auto& p : d.providers)
            if (p.second()) providers += (providers.empty() ? "" : " > ") + p.first;
        fmt::print("  - {:<26} {}  [{}]\n", names + ":", d.description,
                   providers.empty() ? "UNAVAILABLE" : providers);
    }
#ifdef USE_CUDA
    fmt::print("\nCUDA GPU: yes (-gpu cuda for gfn1/gfn2/gfnff)\n");
#endif
#ifdef USE_ROCM
    fmt::print("ROCm GPU: yes (-gpu rocm for gfn1/gfn2/gfnff)\n");
#endif
#ifdef USE_VULKAN
    fmt::print("Vulkan GPU: yes (-gpu vulkan for gfn1/gfn2)\n");
#endif
    fmt::print("===================================\n");
}
