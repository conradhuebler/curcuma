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

# Claude Generated: Unified GFN-FF Test Suite
# Combines regression, parameter, and comprehensive validation tests (Dec 2025)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/HH.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/HH.opt.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/OH.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/dimers/HCl.xyz
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/molecules/dimers/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/molecules/trimers/water.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/trimers/O3.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/trimers/HCN.xyz
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/molecules/trimers/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/molecules/larger/CH4.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/larger/CH3OCH3.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/larger/CH3OH.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/molecules/larger/CH3OCH3.log
          DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/)

add_executable(test_gfnff_unified
        test_gfnff_unified.cpp)

target_link_libraries(test_gfnff_unified curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_unified REUSE_FROM curcuma_core)
endif()

# Register with CTest
add_test(NAME test_gfnff_unified
         COMMAND test_gfnff_unified
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_cases)

# Configure test timeout and grouping
set_tests_properties(test_gfnff_unified PROPERTIES TIMEOUT 60 LABELS "gfnff")

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

# Claude Generated (December 2025): DFT-D3/D4 Dispersion Energy Test Suite
# Tests dispersion corrections with various functionals and parameter combinations

# Copy test molecules to build directory
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/validation/HH.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/validation/benzene.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/validation/ethene.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/validation/butane.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/validation/HCl.xyz
          ${CMAKE_CURRENT_SOURCE_DIR}/validation/OH.xyz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_cases/validation/)

# Create test executable
add_executable(test_dispersion
        test_dispersion.cpp)

# Link against core library
target_link_libraries(test_dispersion curcuma_core)

# Precompiled headers (performance optimization)
if(USE_PCH)
    target_precompile_headers(test_dispersion REUSE_FROM curcuma_core)
endif()

# Register with CTest
add_test(NAME test_dispersion
         COMMAND test_dispersion
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Configure test properties
set_tests_properties(test_dispersion PROPERTIES
    TIMEOUT 120          # 2 minutes (36 tests)
    LABELS "dispersion;D3;D4")

# Conditional test execution (skip if neither D3 nor D4 compiled)
if(NOT USE_D3 AND NOT USE_D4)
    set_tests_properties(test_dispersion PROPERTIES DISABLED TRUE)
endif()

# Claude Generated (December 2025): D3 Debug Test
# Simple test for debugging D3 parameter generation
add_executable(test_d3_debug
        test_d3_debug.cpp)

target_link_libraries(test_d3_debug curcuma_core)

if(USE_PCH)
    target_precompile_headers(test_d3_debug REUSE_FROM curcuma_core)
endif()

add_test(NAME test_d3_debug
         COMMAND test_d3_debug
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Claude Generated (December 2025): D3 Diagnostic Test
# Compares Curcuma D3 implementation against validated formula
add_executable(test_d3_diagnostic
        test_d3_diagnostic.cpp)

target_link_libraries(test_d3_diagnostic curcuma_core)

if(USE_PCH)
    target_precompile_headers(test_d3_diagnostic REUSE_FROM curcuma_core)
endif()

add_test(NAME test_d3_diagnostic
         COMMAND test_d3_diagnostic
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Claude Generated (December 2025): GFN-FF D3 Reference Generator
# Temporary utility to generate reference energies for test_gfnff_d3.cpp
add_executable(generate_gfnff_d3_refs
        generate_gfnff_d3_refs.cpp)

target_link_libraries(generate_gfnff_d3_refs curcuma_core)

if(USE_PCH)
    target_precompile_headers(generate_gfnff_d3_refs REUSE_FROM curcuma_core)
endif()

# Claude Generated (December 2025): GFN-FF D3 Integration Validation Test
# Verifies that GFN-FF correctly uses the native D3 implementation
# Tests standalone D3 vs GFN-FF integrated D3 consistency

# Note: Test molecules already copied by test_gfnff_unified (lines 82-95)
# Molecules needed: dimers/HH.xyz, dimers/HCl.xyz, dimers/OH.xyz,
#                   larger/CH4.xyz, larger/CH3OH.xyz, larger/CH3OCH3.xyz

# Create test executable
add_executable(test_gfnff_d3
        test_gfnff_d3.cpp)

# Link against core library
target_link_libraries(test_gfnff_d3 curcuma_core)

# Precompiled headers (performance optimization)
if(USE_PCH)
    target_precompile_headers(test_gfnff_d3 REUSE_FROM curcuma_core)
endif()

# Register with CTest
add_test(NAME test_gfnff_d3
         COMMAND test_gfnff_d3
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Configure test properties
set_tests_properties(test_gfnff_d3 PROPERTIES
    TIMEOUT 60           # 1 minute (6 molecules, 3 tests each)
    LABELS "gfnff;D3;dispersion")

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

# Add debug executable
#add_executable(debug_molecule
#        ../debug_molecule.cpp)

#target_link_libraries(debug_molecule curcuma_core)

#if(USE_PCH)
#        target_precompile_headers(debug_molecule REUSE_FROM curcuma_core)
#endif()
