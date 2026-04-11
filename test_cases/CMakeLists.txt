file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input_aa.xyz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/input_ab.xyz
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/AAA-bGal/A.xyz
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/AAA-bGal/B.xyz
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Claude Generated: Test Molecule Registry - Central hardcoded molecule database
# Provides zero-dependency molecule access for all test cases
add_library(test_molecule_registry core/test_molecule_registry.cpp)

target_link_libraries(test_molecule_registry curcuma_core)

if(USE_PCH)
    target_precompile_headers(test_molecule_registry REUSE_FROM curcuma_core)
endif()

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

# Claude Generated (February 2026): Unified GFN-FF Validation Runner
# Purpose: High-precision validation against standardized analyzer reference data
add_executable(test_gfnff_validation
        test_gfnff_validation.cpp)

target_link_libraries(test_gfnff_validation curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_validation REUSE_FROM curcuma_core)
endif()

# Copy reference data to build tree so tests can find them
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/reference_data
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Register all reference data files as individual CTests
file(GLOB GFNFF_REF_FILES "${CMAKE_CURRENT_SOURCE_DIR}/reference_data/*.ref.json")
foreach(REF_FILE ${GFNFF_REF_FILES})
    get_filename_component(MOL_NAME ${REF_FILE} NAME_WE)
    set(TEST_NAME "gfnff_val_${MOL_NAME}")

    # Path to the reference file in the build directory
    set(BUILD_REF_FILE "${CMAKE_CURRENT_BINARY_DIR}/reference_data/${MOL_NAME}.ref.json")

    add_test(NAME ${TEST_NAME}
             COMMAND test_gfnff_validation ${BUILD_REF_FILE}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 60 LABELS "gfnff;validation")
endforeach()

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

# Claude Generated (March 2026): GPU GFN-FF Validation Test (gfnff GPU vs CPU)
# Only built and registered when USE_CUDA=ON
if(USE_CUDA)
    add_executable(test_gfnff_gpu test_gfnff_gpu.cpp)
    target_link_libraries(test_gfnff_gpu curcuma_core CUDA::cudart)

    if(USE_PCH)
        target_precompile_headers(test_gfnff_gpu REUSE_FROM curcuma_core)
    endif()

    # Built-in molecule test (water, methane, ethanol) — no arguments
    add_test(NAME gfnff_gpu_builtin_molecules
             COMMAND test_gfnff_gpu
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(gfnff_gpu_builtin_molecules
        PROPERTIES TIMEOUT 120 LABELS "gfnff_gpu;gpu;validation")

    # Per-molecule tests using the same reference data as CPU GFN-FF validation
    foreach(REF_FILE ${GFNFF_REF_FILES})
        get_filename_component(MOL_NAME ${REF_FILE} NAME_WE)
        set(BUILD_REF_FILE "${CMAKE_CURRENT_BINARY_DIR}/reference_data/${MOL_NAME}.ref.json")
        set(GPU_TEST_NAME "gfnff_gpu_vs_cpu_${MOL_NAME}")

        add_test(NAME ${GPU_TEST_NAME}
                 COMMAND test_gfnff_gpu ${BUILD_REF_FILE}
                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(${GPU_TEST_NAME}
            PROPERTIES TIMEOUT 120 LABELS "gfnff_gpu;gpu;validation")
    endforeach()

    # Claude Generated (March 2026): GPU analytical vs numerical gradient test
    add_executable(test_gpu_numgrad test_gpu_numgrad.cpp)
    target_link_libraries(test_gpu_numgrad curcuma_core CUDA::cudart)
    if(USE_PCH)
        target_precompile_headers(test_gpu_numgrad REUSE_FROM curcuma_core)
    endif()
endif()


# Claude Generated: CLI argument validation test suite (October 2025)
# Include comprehensive Bash-based CLI tests for all capabilities
add_subdirectory(cli)
