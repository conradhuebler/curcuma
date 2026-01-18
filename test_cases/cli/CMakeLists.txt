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
    set(PROJECT_ROOT_EXPORT "export PROJECT_ROOT='${CMAKE_CURRENT_SOURCE_DIR}/../..'\n")
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

# =============================================================================
# CG (Coarse Graining) Tests (1 scenario)
# =============================================================================

add_cli_test(cg 01_single_point)

message(STATUS "Added 29 CLI test scenarios with build directory isolation")
message(STATUS "  - Tests execute from: ${CMAKE_CURRENT_BINARY_DIR}/category/test_name")
message(STATUS "  - Source tree remains clean (all outputs in build directory)")