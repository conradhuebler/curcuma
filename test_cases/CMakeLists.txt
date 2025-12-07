file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input_aa.xyz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input_ab.xyz
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/AAA-bGal/A.xyz
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/AAA-bGal/B.xyz
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(rmsd_test
        rmsd/main.cpp)

target_link_libraries(rmsd_test curcuma_core)

add_executable(reorder_test
        reorder/main.cpp)

    add_executable(AAAbGal
            AAAbGal.cpp)
# Claude Generated: Test executable for Molecule class
add_executable(test_molecule
        test_molecule.cpp)

target_link_libraries(AAAbGal curcuma_core)
target_link_libraries(reorder_test curcuma_core)
target_link_libraries(test_molecule curcuma_core)
if(USE_PCH)
        target_precompile_headers(AAAbGal REUSE_FROM curcuma_core)
        target_precompile_headers(reorder_test REUSE_FROM curcuma_core)
        target_precompile_headers(test_molecule REUSE_FROM curcuma_core)
endif()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/confscan/input.xyz
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(confscan_test
        confscan.cpp)

target_link_libraries(confscan_test curcuma_core)

#target_link_libraries(curcuma_tests curcuma_core)￼

add_executable(confstat_tests
    test_confstat.cpp
)

target_link_libraries(confstat_tests
        curcuma_core
)

# Claude Generated: Energy calculation test suite
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/AAA-bGal/A.xyz
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_energy_methods
        test_energy_methods.cpp)

target_link_libraries(test_energy_methods curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_energy_methods REUSE_FROM curcuma_core)
endif()

# Claude Generated: CG Potentials unit test suite
add_executable(test_cg_potentials
        test_cg_potentials.cpp)

target_link_libraries(test_cg_potentials curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_cg_potentials REUSE_FROM curcuma_core)
endif()

# Note: test_cg_forcefield_integration removed due to complex linking with curcuma_D3
# CG functionality will be tested via CLI integration test instead

# Claude Generated: GFN-FF Regression Test Suite
# Tests native cgfnff against XTB 6.6.1 reference data (Nov 2025)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/HH.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/HH.opt.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/OH.xyz
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/molecules/dimers/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/molecules/trimers/water.xyz
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/molecules/trimers/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/molecules/larger/CH4.xyz
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/)

add_executable(test_gfnff_regression
        test_gfnff_regression.cpp)

target_link_libraries(test_gfnff_regression curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_regression REUSE_FROM curcuma_core)
endif()

# Claude Generated: GFN-FF vbond Parameter Verification Test
# Tests that Curcuma's native cgfnff implementation generates correct vbond parameters (Nov 2025)
add_executable(test_gfnff_vbond
        test_gfnff_vbond.cpp)

target_link_libraries(test_gfnff_vbond curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_vbond REUSE_FROM curcuma_core)
endif()

# Claude Generated: Direct GFN-FF Parameter Generation Test
# Tests GFNFF class parameter generation directly with debug output (Nov 2025)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/HH.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/OH.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/HCl.xyz
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/molecules/dimers/)

add_executable(test_gfnff_params_direct
        test_gfnff_params_direct.cpp)

target_link_libraries(test_gfnff_params_direct curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_params_direct REUSE_FROM curcuma_core)
endif()

# Register with CTest
add_test(NAME test_gfnff_regression
         COMMAND test_gfnff_regression ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_gfnff_vbond
         COMMAND test_gfnff_vbond ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_gfnff_params_direct
         COMMAND test_gfnff_params_direct)

# Claude Generated: GFN-FF CH3OH Comprehensive Validation Test (Dec 2025)
# Validates Curcuma native GFN-FF against XTB 6.6.1 Fortran reference
# Tests: Coordination Numbers, EEQ charges, Total energy
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/molecules/larger/CH3OH.xyz
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/)

add_executable(test_gfnff_ch3oh
        test_gfnff_ch3oh.cpp)

target_link_libraries(test_gfnff_ch3oh curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_ch3oh REUSE_FROM curcuma_core)
endif()

add_test(NAME test_gfnff_ch3oh
         COMMAND test_gfnff_ch3oh)

# Claude Generated (Session 5): Periodic Table utility tests
add_executable(test_periodic_table
        test_periodic_table.cpp)

target_link_libraries(test_periodic_table curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_periodic_table REUSE_FROM curcuma_core)
endif()

add_test(NAME test_periodic_table
         COMMAND test_periodic_table)

# Claude Generated (Session 5): Two-Phase EEQ System Tests
# Tests Phase 1 topology charges and Phase 2 environment corrections
add_executable(test_eeq_two_phase
        test_eeq_two_phase.cpp)

target_link_libraries(test_eeq_two_phase curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_eeq_two_phase REUSE_FROM curcuma_core)
endif()

add_test(NAME test_eeq_two_phase
         COMMAND test_eeq_two_phase)

# todo restore

# Claude Generated: Debug Ulysses methods test
#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/../release/debug_ulysses_test.cpp
#        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#add_executable(debug_ulysses_test
#        ${CMAKE_CURRENT_BINARY_DIR}/debug_ulysses_test.cpp)

#target_link_libraries(debug_ulysses_test curcuma_core)

#if(USE_PCH)
#        target_precompile_headers(debug_ulysses_test REUSE_FROM curcuma_core)
#endif()

# Claude Generated: CLI argument validation test suite (October 2025)
# Include comprehensive Bash-based CLI tests for all capabilities
add_subdirectory(cli)
