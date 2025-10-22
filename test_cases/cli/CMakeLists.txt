# CLI Test Suite - CMake Integration
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Comprehensive CLI argument validation tests

# All tests execute from their own directories to access input files

# =============================================================================
# curcumaopt Tests (6 scenarios)
# =============================================================================

add_test(
    NAME cli_curcumaopt_01_default_uff
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/01_default_uff_opt/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/01_default_uff_opt
)

add_test(
    NAME cli_curcumaopt_02_gfn2_single_point
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/02_gfn2_single_point/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/02_gfn2_single_point
)

add_test(
    NAME cli_curcumaopt_03_invalid_method
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/03_invalid_method/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/03_invalid_method
)

add_test(
    NAME cli_curcumaopt_04_lbfgs_params
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/04_lbfgs_params/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/04_lbfgs_params
)

add_test(
    NAME cli_curcumaopt_05_alias_singlepoint
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/05_alias_singlepoint/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/05_alias_singlepoint
)

add_test(
    NAME cli_curcumaopt_06_opth_only
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/06_opth_only/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/curcumaopt/06_opth_only
)

# =============================================================================
# RMSD Tests (6 scenarios)
# =============================================================================

add_test(
    NAME cli_rmsd_01_default
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/01_default_rmsd/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/01_default_rmsd
)

add_test(
    NAME cli_rmsd_02_no_reorder
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/02_no_reorder/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/02_no_reorder
)

add_test(
    NAME cli_rmsd_03_invalid_method
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/03_invalid_method/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/03_invalid_method
)

add_test(
    NAME cli_rmsd_04_template_elements
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/04_template_elements/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/04_template_elements
)

add_test(
    NAME cli_rmsd_05_fragment
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/05_fragment_rmsd/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/05_fragment_rmsd
)

add_test(
    NAME cli_rmsd_06_alias_reorder
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/06_alias_reorder/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rmsd/06_alias_reorder
)

# =============================================================================
# ConfScan Tests (7 scenarios)
# =============================================================================

add_test(
    NAME cli_confscan_01_default
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/confscan/01_default_scan/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/confscan/01_default_scan
)

add_test(
    NAME cli_confscan_02_get_rmsd
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/confscan/02_get_rmsd/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/confscan/02_get_rmsd
)

add_test(
    NAME cli_confscan_03_invalid_method
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/confscan/03_invalid_rmsd_method/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/confscan/03_invalid_rmsd_method
)

add_test(
    NAME cli_confscan_04_slx_logic
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/confscan/04_slx_logic/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/confscan/04_slx_logic
)

add_test(
    NAME cli_confscan_05_hybrid_elements
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/confscan/05_hybrid_rmsd_elements/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/confscan/05_hybrid_rmsd_elements
)

add_test(
    NAME cli_confscan_06_heavy_only
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/confscan/06_heavy_only/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/confscan/06_heavy_only
)

add_test(
    NAME cli_confscan_07_restart
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/confscan/07_restart_scan/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/confscan/07_restart_scan
)

# =============================================================================
# SimpleMD Tests (7 scenarios)
# =============================================================================

add_test(
    NAME cli_simplemd_01_nve_short
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/01_nve_short_run/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/01_nve_short_run
)

add_test(
    NAME cli_simplemd_02_nvt_berendsen
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/02_nvt_berendsen/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/02_nvt_berendsen
)

add_test(
    NAME cli_simplemd_03_invalid_thermostat
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/03_invalid_thermostat/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/03_invalid_thermostat
)

add_test(
    NAME cli_simplemd_04_rattle
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/04_rattle_all/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/04_rattle_all
)

add_test(
    NAME cli_simplemd_05_wall
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/05_spheric_wall/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/05_spheric_wall
)

add_test(
    NAME cli_simplemd_06_alias_temperature
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/06_alias_temperature/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/06_alias_temperature
)

add_test(
    NAME cli_simplemd_07_restart
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/07_restart_simulation/run_test.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simplemd/07_restart_simulation
)

# =============================================================================
# Test Properties
# =============================================================================

# Set timeout for all CLI tests (30 seconds should be enough for short runs)
set_tests_properties(
    cli_curcumaopt_01_default_uff
    cli_curcumaopt_02_gfn2_single_point
    cli_curcumaopt_03_invalid_method
    cli_curcumaopt_04_lbfgs_params
    cli_curcumaopt_05_alias_singlepoint
    cli_curcumaopt_06_opth_only
    cli_rmsd_01_default
    cli_rmsd_02_no_reorder
    cli_rmsd_03_invalid_method
    cli_rmsd_04_template_elements
    cli_rmsd_05_fragment
    cli_rmsd_06_alias_reorder
    cli_confscan_01_default
    cli_confscan_02_get_rmsd
    cli_confscan_03_invalid_method
    cli_confscan_04_slx_logic
    cli_confscan_05_hybrid_elements
    cli_confscan_06_heavy_only
    cli_confscan_07_restart
    cli_simplemd_01_nve_short
    cli_simplemd_02_nvt_berendsen
    cli_simplemd_03_invalid_thermostat
    cli_simplemd_04_rattle
    cli_simplemd_05_wall
    cli_simplemd_06_alias_temperature
    cli_simplemd_07_restart
    PROPERTIES TIMEOUT 30
)

message(STATUS "Added 26 CLI test scenarios (curcumaopt: 6, rmsd: 6, confscan: 7, simplemd: 7)")
