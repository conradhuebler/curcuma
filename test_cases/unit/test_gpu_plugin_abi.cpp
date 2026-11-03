/*
 * Copyright (C) 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 */
/**
 * GPU plugin ABI version check (S0) - Claude Generated (Nov 2026).
 *
 * The real loader (gpu_plugin::available / rejectionReason, the same pluginHandle() the create* calls use) is pointed
 * at stand-in plugins through $CURCUMA_PLUGIN_DIR. Runs without a GPU. Stand-ins (gpu_plugin_abi_stub/stub_plugin.cpp):
 *   libcurcuma_abitest_none.so   no curcuma_abitest_none_abi_version symbol (plugin built before S0) -> refused
 *   libcurcuma_abitest_wrong.so  version CURCUMA_GPU_PLUGIN_ABI_VERSION + 1                         -> refused
 *   libcurcuma_abitest_ok.so     version CURCUMA_GPU_PLUGIN_ABI_VERSION                             -> accepted
 *   libcurcuma_abitest_align.so  right version, Eigen alignment 8 (never a real value)              -> refused
 * Checked: the decision, the reason text (path, found and expected version, "rebuild the GPU plugin from the same
 * source tree"), and checkPluginCompatibility() directly on a dlopen handle.
 * The CLI tests cli_gpu_plugin_abi_* cover the method factory (counted fallback, CPU energy, -gpu_strict exit 3).
 */
#include "src/core/curcuma_logger.h"
#include "src/core/energy_calculators/gpu_plugin.h"
#include "src/core/energy_calculators/gpu_plugin_abi.h"

#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <string>

static int g_failures = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok)
        ++g_failures;
}
static bool has(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::printf("usage: test_gpu_plugin_abi <stub plugin dir>\n");
        return 2;
    }
    const std::string dir = argv[1];
    setenv("CURCUMA_PLUGIN_DIR", dir.c_str(), 1);
    CurcumaLogger::set_verbosity(1);   // the refusal warnings are part of what is shown
    const std::string expected = std::to_string(CURCUMA_GPU_PLUGIN_ABI_VERSION);

    // (1) no version symbol
    check(!gpu_plugin::available("abitest_none"), "(1) plugin without abi_version symbol is refused");
    const std::string r1 = gpu_plugin::rejectionReason("abitest_none");
    check(has(r1, dir + "/libcurcuma_abitest_none.so"), "(1) reason names the plugin path: " + r1);
    check(has(r1, "curcuma_abitest_none_abi_version missing"), "(1) reason names the missing symbol");
    check(has(r1, "expects " + expected), "(1) reason names the expected version");
    check(has(r1, "rebuild the GPU plugin from the same source tree"), "(1) reason says to rebuild");
    check(gpu_plugin::resolveEntryPoint("abitest_none", "curcuma_abitest_none_device_count") == nullptr,
        "(1) no entry point of a refused plugin is reachable");
    check(gpu_plugin::deviceCount("abitest_none") == 0, "(1) refused plugin reports 0 devices");

    // (2) wrong version
    check(!gpu_plugin::available("abitest_wrong"), "(2) plugin with another ABI version is refused");
    const std::string r2 = gpu_plugin::rejectionReason("abitest_wrong");
    check(has(r2, "ABI version " + std::to_string(CURCUMA_GPU_PLUGIN_ABI_VERSION + 1) + ", this curcuma expects " + expected),
        "(2) reason names found and expected version: " + r2);

    // (3) right version
    check(gpu_plugin::available("abitest_ok"), "(3) plugin with the right ABI version is accepted");
    check(gpu_plugin::rejectionReason("abitest_ok").empty(), "(3) no rejection reason");
    check(gpu_plugin::resolveEntryPoint("abitest_ok", "curcuma_abitest_ok_device_count") != nullptr,
        "(3) entry point of the accepted plugin is reachable");

    // (4) right version, wrong Eigen alignment (the older variant check, still in force)
    check(!gpu_plugin::available("abitest_align"), "(4) plugin with another Eigen alignment is refused");
    check(has(gpu_plugin::rejectionReason("abitest_align"), "Eigen alignment 8"), "(4) reason names the alignment");

    // (5) absent plugin: not a rejection
    check(!gpu_plugin::available("abitest_absent") && gpu_plugin::rejectionReason("abitest_absent").empty(),
        "(5) absent plugin: unavailable, no rejection reason");

    // (6) the check function on its own
    for (const char* b : { "none", "wrong", "ok" }) {
        const std::string backend = std::string("abitest_") + b;
        const std::string path = dir + "/libcurcuma_" + backend + ".so";
        void* h = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
        check(h != nullptr, std::string("(6) dlopen ") + path);
        if (!h) continue;
        const bool accepted = gpu_plugin::checkPluginCompatibility(h, backend, path).empty();
        check(accepted == (std::string(b) == "ok"), std::string("(6) checkPluginCompatibility on ") + b);
        dlclose(h);
    }

    std::printf(g_failures == 0 ? "gpu_plugin_abi: all checks passed\n" : "gpu_plugin_abi: %d failure(s)\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
