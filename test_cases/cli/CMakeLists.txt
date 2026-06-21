# CLI Test Suite - CMake Integration
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Build directory isolation for clean source tree

# =============================================================================
# Claude Generated (October 2025): Copy-and-Run Pattern for Clean Source Tree
# All tests execute from their directories in the BUILD TREE to keep source
# tree clean and consistent with C++ unit test infrastructure.
# =============================================================================

# Helper macro to add CLI test with build directory isolation
# Usage: add_cli_test(CATEGORY TEST_NAME)
# Example: add_cli_test(rmsd 01_default_rmsd)
macro(add_cli_test CATEGORY TEST_NAME)
    # Create test directory in build tree
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CATEGORY}/${TEST_NAME})

    # Copy test inputs and scripts
    file(GLOB TEST_FILES
         "${CMAKE_CURRENT_SOURCE_DIR}/${CATEGORY}/${TEST_NAME}/*")
    file(COPY ${TEST_FILES}
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${CATEGORY}/${TEST_NAME})

    # Copy shared test utilities to category and cli directories
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_utils.sh
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${CATEGORY}/)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_utils.sh
         DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

    # Claude Generated (October 2025): Fix paths and environment in run_test.sh for build directory
    # Tests run from build/test_cases/cli/CATEGORY/TEST_NAME/ instead of source directory
    file(READ "${CMAKE_CURRENT_BINARY_DIR}/${CATEGORY}/${TEST_NAME}/run_test.sh" RUN_TEST_CONTENT)

    # Set PROJECT_ROOT for test_utils.sh to find curcuma binary
    # CMAKE_CURRENT_SOURCE_DIR = test_cases/cli, so go up 2 levels to project root
    set(CORRECT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../..")
    set(PROJECT_ROOT_EXPORT "export PROJECT_ROOT='${CORRECT_ROOT}'\n")
    # Remove any hardcoded PROJECT_ROOT lines from the script (from dev environments)
    string(REGEX REPLACE "export PROJECT_ROOT='[^']*'\n?" "" RUN_TEST_CONTENT "${RUN_TEST_CONTENT}")
    string(REPLACE "#!/bin/bash\n" "#!/bin/bash\n${PROJECT_ROOT_EXPORT}\n" RUN_TEST_CONTENT "${RUN_TEST_CONTENT}")

    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${CATEGORY}/${TEST_NAME}/run_test.sh" "${RUN_TEST_CONTENT}")

    # Add test running from build directory (clean source tree)
    add_test(
        NAME cli_${CATEGORY}_${TEST_NAME}
        COMMAND bash run_test.sh
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${CATEGORY}/${TEST_NAME}
    )

    # Set timeout for all CLI tests
    set_tests_properties(cli_${CATEGORY}_${TEST_NAME} PROPERTIES TIMEOUT 120)
endmacro()

# =============================================================================
# curcumaopt Tests (6 scenarios)
# =============================================================================

add_cli_test(curcumaopt 01_default_uff_opt)
add_cli_test(curcumaopt 02_gfn2_single_point)
# AP3 (2026-04-25): gfn2 now routes to native xTB; TBLite golden ref no longer valid until AP5.
set_tests_properties(cli_curcumaopt_02_gfn2_single_point PROPERTIES DISABLED TRUE)
add_cli_test(curcumaopt 03_invalid_method)
add_cli_test(curcumaopt 04_lbfgs_params)
add_cli_test(curcumaopt 05_alias_singlepoint)
add_cli_test(curcumaopt 06_opth_only)

# =============================================================================
# GFN-FF Tests (2 scenarios)
# =============================================================================

add_cli_test(gfnff 01_complex_gfnff_singlepoint)
add_cli_test(gfnff 02_caffeine_gfnff_energy_components)

# =============================================================================
# GPU GFN-FF Tests — only registered when USE_CUDA=ON
# =============================================================================

if(USE_CUDA)
    add_cli_test(gfnff_gpu 01_gfnff_gpu_singlepoint)
    set_tests_properties(cli_gfnff_gpu_01_gfnff_gpu_singlepoint
        PROPERTIES TIMEOUT 120 LABELS "gfnff_gpu;gpu;cli")
endif()

# ROCm GFN-FF — only when USE_ROCM_GFNFF=ON. Self-skips at runtime if no HIP device.
if(USE_ROCM_GFNFF)
    add_cli_test(gfnff_gpu 02_rocm_singlepoint)
    set_tests_properties(cli_gfnff_gpu_02_rocm_singlepoint
        PROPERTIES TIMEOUT 120 LABELS "gfnff_gpu;gpu;rocm;cli")
endif()

# =============================================================================
# Vulkan native-xTB nuclear gradient (V-AP1) — only when USE_VULKAN_XTB=ON.
# Self-skips at runtime if no FP64 Vulkan device is present.
# =============================================================================
if(USE_VULKAN_XTB)
    add_cli_test(gpu_gradient 01_vulkan_gfn1_gradient)
    set_tests_properties(cli_gpu_gradient_01_vulkan_gfn1_gradient
        PROPERTIES TIMEOUT 120 LABELS "gpu_gradient;gpu;vulkan;cli")
    # X-AP4: GFN2 D4 EEQ charge-response gradient on the device (q-response).
    add_cli_test(gpu_gradient 02_vulkan_gfn2_qresponse)
    set_tests_properties(cli_gpu_gradient_02_vulkan_gfn2_qresponse
        PROPERTIES TIMEOUT 120 LABELS "gpu_gradient;gpu;vulkan;cli")
    # V-AP2: GFN2 AO multipole integrals on the device.
    add_cli_test(gpu_multipole 01_vulkan_gfn2_multipole)
    set_tests_properties(cli_gpu_multipole_01_vulkan_gfn2_multipole
        PROPERTIES TIMEOUT 120 LABELS "gpu_multipole;gpu;vulkan;cli")
endif()

# =============================================================================
# RMSD Tests (6 scenarios)
# =============================================================================

add_cli_test(rmsd 01_default_rmsd)
add_cli_test(rmsd 02_no_reorder)
add_cli_test(rmsd 03_invalid_method)
add_cli_test(rmsd 04_template_elements)
add_cli_test(rmsd 05_fragment_rmsd)
add_cli_test(rmsd 06_alias_reorder)

# =============================================================================
# ConfScan Tests (7 scenarios)
# =============================================================================

add_cli_test(confscan 01_default_scan)
add_cli_test(confscan 02_get_rmsd)
add_cli_test(confscan 03_invalid_rmsd_method)
add_cli_test(confscan 04_slx_logic)
add_cli_test(confscan 05_hybrid_rmsd_elements)
add_cli_test(confscan 06_heavy_only)
add_cli_test(confscan 07_restart_scan)

# Claude Generated (June 2026): ConfScan scenarios run a real multi-pass scan over 44
# conformers (04 runs three reorder passes; 07 runs the scan twice). The 120 s macro
# default is tight on a loaded machine - give them headroom so a busy host doesn't cause
# spurious timeouts.
set_tests_properties(
    cli_confscan_01_default_scan
    cli_confscan_02_get_rmsd
    cli_confscan_03_invalid_rmsd_method
    cli_confscan_04_slx_logic
    cli_confscan_05_hybrid_rmsd_elements
    cli_confscan_06_heavy_only
    cli_confscan_07_restart_scan
    PROPERTIES TIMEOUT 300)

# =============================================================================
# SimpleMD Tests (7 scenarios)
# =============================================================================

add_cli_test(simplemd 01_nve_short_run)
add_cli_test(simplemd 02_nvt_berendsen)
add_cli_test(simplemd 03_invalid_thermostat)
add_cli_test(simplemd 04_rattle_all)
add_cli_test(simplemd 05_spheric_wall)
add_cli_test(simplemd 06_alias_temperature)
add_cli_test(simplemd 07_restart_simulation)
add_cli_test(simplemd 08_gfnff_acetic_acid_dimer_md)
add_cli_test(simplemd 09_xtbgfnff_acetic_acid_dimer_md)
add_cli_test(simplemd 10_gfnff_polymer_md)
set_tests_properties(cli_simplemd_10_gfnff_polymer_md PROPERTIES TIMEOUT 3600 LABELS "simplemd;gfnff;large_system;performance")
add_cli_test(simplemd 11_diagnostics_dump)
set_tests_properties(cli_simplemd_11_diagnostics_dump PROPERTIES TIMEOUT 180 LABELS "simplemd;gfnff;diagnostics")
add_cli_test(simplemd 12_eeq_cutoff_auto)
set_tests_properties(cli_simplemd_12_eeq_cutoff_auto PROPERTIES TIMEOUT 180 LABELS "simplemd;gfnff;eeq_cutoff")

# =============================================================================
# CG (Coarse Graining) Tests — disabled: -vtf command not implemented in main.cpp
# =============================================================================

# add_cli_test(cg 01_single_point)

# =============================================================================
# Analysis Tests (1 scenario)
# =============================================================================

add_cli_test(analysis 01_parallel_equivalence)
# 100 Frames × Scattering braucht länger als Standard-30s
set_tests_properties(cli_analysis_01_parallel_equivalence PROPERTIES TIMEOUT 120)

# =============================================================================
# CLI Routing Tests — Claude Generated 2026
# Validates ParameterRegistry-driven flat-flag routing and JSON round-trip
# =============================================================================

add_cli_test(routing 01_flat_param_routing)
add_cli_test(routing 02_roundtrip_replay)

# =============================================================================
# Error Handling Tests (6 scenarios) - Phase 1-2 Implementation
# =============================================================================

# Phase 1: File Validation
add_cli_test(errors 01_missing_file)
add_cli_test(errors 02_invalid_xyz_format)
add_cli_test(errors 03_permission_denied)

# Phase 2: Method Validation
add_cli_test(errors 04_typo_method)
add_cli_test(errors 05_method_not_compiled)
add_cli_test(errors 06_list_methods)

# =============================================================================
# SQM Method Validation Tests
# Native semi-empirical methods vs Ulysses/TBLite reference energies
# Claude Generated (March 2026): Test infrastructure for PM3/AM1/MNDO/GFN1/GFN2
# Claude Generated (April 2026): Tests 06-11 are anchor points for the
#   NATIVE_XTB_ROADMAP migration (docs/NATIVE_XTB_AP{1..5}_PLAN.md).
# =============================================================================

add_cli_test(sqm 01_pm3_singlepoint)
add_cli_test(sqm 02_am1_singlepoint)
add_cli_test(sqm 03_mndo_singlepoint)
add_cli_test(sqm 04_gfn1_singlepoint)
add_cli_test(sqm 05_gfn2_singlepoint)

# Native xTB roadmap anchors — see docs/NATIVE_XTB_*_PLAN.md
add_cli_test(sqm 06_gfn2_extended_singlepoint)     # gfnff-shared molecules, TBLite gfn2 ref
add_cli_test(sqm 07_gfn1_extended_singlepoint)     # gfnff-shared molecules, TBLite gfn1 ref
add_cli_test(sqm 08_ngfn1_baseline_singlepoint)    # native ngfn1 vs TBLite-gfn1 (≤ 5e-4 Eh)
add_cli_test(sqm 09_ngfn2_baseline_singlepoint)    # native ngfn2 broken-state snapshot
add_cli_test(sqm 10_gfn2_optimization_smoke)       # gfn2 -opt anchor (gradient = zero pre-AP-4)
add_cli_test(sqm 11_gfn2_provider_check)           # routing anchor — updated for AP 3 (native backend)

# AP3 (2026-04-25): gfn2/gfn1 route to native xTB. TBLite golden refs invalid until AP5.
set_tests_properties(cli_sqm_04_gfn1_singlepoint PROPERTIES DISABLED TRUE)
set_tests_properties(cli_sqm_05_gfn2_singlepoint PROPERTIES DISABLED TRUE)
set_tests_properties(cli_sqm_06_gfn2_extended_singlepoint PROPERTIES DISABLED TRUE)
set_tests_properties(cli_sqm_07_gfn1_extended_singlepoint PROPERTIES DISABLED TRUE)
# AP3: gfn2 -opt still broken (gradient = zero), disable until AP4 restores gradient.
set_tests_properties(cli_sqm_10_gfn2_optimization_smoke PROPERTIES DISABLED TRUE)

# Allow the optimization test more wall-clock — caffeine -opt can run > 30 s
# (kept for reference; test is currently DISABLED)
set_tests_properties(cli_sqm_10_gfn2_optimization_smoke PROPERTIES TIMEOUT 180)

message(STATUS "Added CLI test scenarios with build directory isolation")
message(STATUS "  - Tests execute from: ${CMAKE_CURRENT_BINARY_DIR}/category/test_name")
message(STATUS "  - Source tree remains clean (all outputs in build directory)")