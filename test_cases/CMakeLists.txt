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
