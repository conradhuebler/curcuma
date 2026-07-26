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

# Claude Generated 2026: BMT Utils unit test suite
add_executable(test_bmt_utils
        test_bmt_utils.cpp)

target_link_libraries(test_bmt_utils curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_bmt_utils REUSE_FROM curcuma_core)
endif()

add_test(NAME test_bmt_utils COMMAND test_bmt_utils)
set_tests_properties(test_bmt_utils PROPERTIES TIMEOUT 10 LABELS "bmt;unit")

# Claude Generated 2026: PLUMED init unit test suite (only built with USE_Plumed)
if(USE_Plumed)
    add_executable(test_plumed_init
            test_plumed_init.cpp)
    target_link_libraries(test_plumed_init curcuma_core)
    if(USE_PCH)
        target_precompile_headers(test_plumed_init REUSE_FROM curcuma_core)
    endif()
    add_test(NAME test_plumed_init COMMAND test_plumed_init)
    set_tests_properties(test_plumed_init PROPERTIES TIMEOUT 10 LABELS "plumed;unit")
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

# Claude Generated (June 2026): per-frame CPU-vs-GPU gradient diagnostic on a fixed
# trajectory — the clean force-backend metric (replaces the misleading MD heat-exchange).
# Builds unconditionally; the gpu config falls back to CPU when CUDA is not compiled in.
add_executable(test_gfnff_grad_traj test_gfnff_grad_traj.cpp)
target_link_libraries(test_gfnff_grad_traj curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_gfnff_grad_traj REUSE_FROM curcuma_core)
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

# Claude Generated (June 2026, Phase 0 ORCA improvement plan):
# Unit tests for ORCA interface helpers, parsers, and CG-atom guard.
add_executable(test_orca_interface
        test_orca_interface.cpp)

target_link_libraries(test_orca_interface curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_orca_interface REUSE_FROM curcuma_core)
endif()

add_test(NAME test_orca_interface
         COMMAND test_orca_interface)
set_tests_properties(test_orca_interface PROPERTIES TIMEOUT 30 LABELS "orca;unit")

# Claude Generated (May 2026, WP-D Stage B): fast_exp_neg_sq_block unit test
add_executable(test_fast_exp
        test_fast_exp.cpp)

target_link_libraries(test_fast_exp curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_fast_exp REUSE_FROM curcuma_core)
endif()

# Claude Generated (June 2026): opt-in multi-step EEQ PCG warm-start extrapolation.
# Drives GFN-FF (forced PCG) over a smooth trajectory of a multi-fragment system;
# asserts aspc/gauss energies == none (PCG result-neutral) and that the extrapolation
# path fired. Run via: ctest -L eeq_extrapolation
add_executable(test_eeq_extrapolation
        test_eeq_extrapolation.cpp)
target_link_libraries(test_eeq_extrapolation curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_eeq_extrapolation REUSE_FROM curcuma_core)
endif()
add_test(NAME test_eeq_extrapolation
         COMMAND test_eeq_extrapolation
                 ${CMAKE_CURRENT_SOURCE_DIR}/molecules/larger/complex.xyz 12)
set_tests_properties(test_eeq_extrapolation PROPERTIES TIMEOUT 180 LABELS "eeq_extrapolation")

add_test(NAME test_fast_exp
         COMMAND test_fast_exp)

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

# Claude Generated (Jul 2026): GFN-FF topology re-entrancy regression guard.
# calculateTopologyInfoOnce() must give identical charges when called twice on
# identical input, both plainly and across an EEQ cache invalidation (the MD/opt
# rebuild path). Guards the Phase-1/Phase-2 Cholesky cache collision.
add_executable(test_gfnff_topology_reentrancy
        test_gfnff_topology_reentrancy.cpp)

target_link_libraries(test_gfnff_topology_reentrancy curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_topology_reentrancy REUSE_FROM curcuma_core)
endif()

add_test(NAME test_gfnff_topology_reentrancy
         COMMAND test_gfnff_topology_reentrancy ${CMAKE_SOURCE_DIR}/test_cases)

# Claude Generated: PolymerBuild sequence parser tests (prime notation)
add_executable(test_polymerbuild_sequence
        test_polymerbuild_sequence.cpp)

target_link_libraries(test_polymerbuild_sequence curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_polymerbuild_sequence REUSE_FROM curcuma_core)
endif()

add_test(NAME test_polymerbuild_sequence
         COMMAND test_polymerbuild_sequence)

# Claude Generated: PolymerBuild structural tests (Xx selection, geometry comparison)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/polymerbuild/fragment_2xx.xyz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/polymerbuild)

add_executable(test_polymerbuild_structure
        test_polymerbuild_structure.cpp)

target_link_libraries(test_polymerbuild_structure curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_polymerbuild_structure REUSE_FROM curcuma_core)
endif()

add_test(NAME test_polymerbuild_structure
         COMMAND test_polymerbuild_structure
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(test_polymerbuild_structure PROPERTIES TIMEOUT 120)

# Claude Generated (April 2026): CPU GFN-FF analytical vs numerical gradient test
# Tests gradient correctness for all built-in molecules (H2O, CH4, CH3OH, CH3OCH3).
# gfnff_numgrad_builtin compares the analytical Gradient() against the FULL finite
# difference (NumGrad, EEQ + CN recomputed each step) — the physically meaningful
# total derivative an optimiser/MD consumes.
add_executable(test_gfnff_numgrad test_gfnff_numgrad.cpp)
target_link_libraries(test_gfnff_numgrad curcuma_core test_molecule_registry)

if(USE_PCH)
    target_precompile_headers(test_gfnff_numgrad REUSE_FROM curcuma_core)
endif()

add_test(NAME gfnff_numgrad_builtin
         COMMAND test_gfnff_numgrad
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(gfnff_numgrad_builtin
    PROPERTIES TIMEOUT 120 LABELS "gfnff;gradient;cpu")

# NOTE (Claude Generated, June 2026): the former gfnff_numgrad_fixed_charges ctest was
# removed. It compared the analytical gradient against the FIXED-charge finite difference
# (--fixed-charges, EEQ frozen) on the premise that the analytical gradient was a
# frozen-charge approximation (no dq/dx). That premise is obsolete: since the Feb-2026
# Coulomb Term-1b / dynamic-charge gradient work the analytical Gradient() is the TOTAL
# derivative (it includes the charge response), so it correctly disagrees with the
# fixed-charge partial derivative — for polar molecules the charge-response term is large
# and opposes the partial geometry term (H2O: total -1.2e-2 vs partial +2.9e-2 Eh/Bohr).
# The full-EEQ gfnff_numgrad_builtin validates the gradient; the --fixed-charges flag is
# retained on the binary for manual formula diagnostics only (not a pass/fail gate).

# Claude Generated (WP5, June 2026): GFN-FF self-consistent ALPB solvation gradient.
# Analytical Gradient() vs NumGrad() (full FD incl. solvated EEQ re-solve) in water.
# Same 5e-3 Eh/Bohr full-EEQ tolerance as the gas-phase test (the reference gfnff
# gradient is itself frozen-charge — see docs/SQM_SOLVATION_WP.md WP5/S5).
add_test(NAME gfnff_numgrad_alpb_water
         COMMAND test_gfnff_numgrad --solvent water
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(gfnff_numgrad_alpb_water
    PROPERTIES TIMEOUT 180 LABELS "gfnff;gradient;solvation;cpu")

# Claude Generated (June 2026): native GFN1/GFN2 ALPB-solvation gradient test.
# Analytical XTB::Gradient() vs central finite difference of the total energy,
# GFN1 (CM5) + GFN2 (Mulliken) in water. See docs/SQM_SOLVATION_WP.md (WP1.4).
add_executable(test_xtb_solvation_numgrad test_xtb_solvation_numgrad.cpp)
target_link_libraries(test_xtb_solvation_numgrad curcuma_core test_molecule_registry)
target_include_directories(test_xtb_solvation_numgrad PRIVATE ${CMAKE_SOURCE_DIR})
if(USE_PCH)
    target_precompile_headers(test_xtb_solvation_numgrad REUSE_FROM curcuma_core)
endif()
add_test(NAME xtb_solvation_numgrad
         COMMAND test_xtb_solvation_numgrad
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(xtb_solvation_numgrad
    PROPERTIES TIMEOUT 180 LABELS "sqm;solvation;gradient;cpu")

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


# Claude Generated (AP5, Apr 2026): GFN1/GFN2 analytical vs numerical gradient test
add_executable(test_xtb_gradient test_xtb_gradient.cpp)
target_link_libraries(test_xtb_gradient curcuma_core)
if(USE_PCH)
    target_precompile_headers(test_xtb_gradient REUSE_FROM curcuma_core)
endif()

set(_GRAD_MOL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/sqm_reference/molecules)
foreach(_mol H2O CH4 NH3)
    add_test(NAME xtb_gradient_${_mol}
             COMMAND test_xtb_gradient ${_GRAD_MOL_DIR}/${_mol}.xyz
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(xtb_gradient_${_mol}
        PROPERTIES TIMEOUT 120 LABELS "gradient;native_xtb")
endforeach()

# Claude Generated (WP4, Jun 2026): native self-contained symmetric eigensolver
# (the `-eigensolver native` backend) vs Eigen on random symmetric matrices.
add_executable(test_native_eigensolver test_native_eigensolver.cpp)
target_link_libraries(test_native_eigensolver curcuma_core)
if(USE_PCH)
    target_precompile_headers(test_native_eigensolver REUSE_FROM curcuma_core)
endif()
add_test(NAME native_eigensolver COMMAND test_native_eigensolver
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(native_eigensolver PROPERTIES TIMEOUT 60 LABELS "native_xtb;eigensolver")

# Same binary, the alternative tridiagonalization backends (CURCUMA_EIG_TRED2):
# 'blocked' = our MKL-free BLAS-3 panel reduction with compact-WY block-Q accumulation
# (the GPU-portable kernel), 'scalar' = the unvectorized reference baseline. The default
# run above covers Eigen's vectorized reduction. Claude Generated.
add_test(NAME native_eigensolver_blocked COMMAND test_native_eigensolver
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(native_eigensolver_blocked PROPERTIES TIMEOUT 60
    LABELS "native_xtb;eigensolver" ENVIRONMENT "CURCUMA_EIG_TRED2=blocked")
add_test(NAME native_eigensolver_scalar COMMAND test_native_eigensolver
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(native_eigensolver_scalar PROPERTIES TIMEOUT 60
    LABELS "native_xtb;eigensolver" ENVIRONMENT "CURCUMA_EIG_TRED2=scalar")

# Claude Generated (Jun 2026): 0 K density-matrix purification kernel (the GEMM-only
# `-eigensolver purify` density path) vs the occupied-subspace projector from Eigen.
add_executable(test_purification test_purification.cpp)
target_link_libraries(test_purification curcuma_core)
if(USE_PCH)
    target_precompile_headers(test_purification REUSE_FROM curcuma_core)
endif()
add_test(NAME purification COMMAND test_purification
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(purification PROPERTIES TIMEOUT 60 LABELS "native_xtb;eigensolver;purification")

# Claude Generated (Jun 2026): seeded block LOBPCG kernel (the iterative `-eigensolver lobpcg`
# path) vs Eigen's lowest-k eigenpairs, incl. a warm-start (subspace-recycling) case.
add_executable(test_lobpcg test_lobpcg.cpp)
target_link_libraries(test_lobpcg curcuma_core)
if(USE_PCH)
    target_precompile_headers(test_lobpcg REUSE_FROM curcuma_core)
endif()
add_test(NAME lobpcg COMMAND test_lobpcg
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(lobpcg PROPERTIES TIMEOUT 60 LABELS "native_xtb;eigensolver;lobpcg")

# Claude Generated (2026, AP ∂q/∂x — Phase 1): D4 dE/dq analytical-vs-FD test
add_executable(test_d4_dedq test_d4_dedq.cpp)
target_link_libraries(test_d4_dedq curcuma_core test_molecule_registry)
if(USE_PCH)
    target_precompile_headers(test_d4_dedq REUSE_FROM curcuma_core)
endif()
add_test(NAME d4_dedq COMMAND test_d4_dedq
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(d4_dedq PROPERTIES TIMEOUT 120 LABELS "gradient;dispersion")

# Claude Generated (Phase 3b-2, May 2026): CPSCF orbital-Hessian + Z-vector gate test
add_executable(test_xtb_cpscf test_xtb_cpscf.cpp)
target_link_libraries(test_xtb_cpscf curcuma_core test_molecule_registry)
if(USE_PCH)
    target_precompile_headers(test_xtb_cpscf REUSE_FROM curcuma_core)
endif()
add_test(NAME xtb_cpscf COMMAND test_xtb_cpscf
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(xtb_cpscf PROPERTIES TIMEOUT 120 LABELS "cpscf;native_xtb")

# Claude Generated (large_system_mode, Jun 2026): native-GFN large-system modes
# (-large_system_mode fragments|dc|sparse) fast-proxy regression on caffeine —
# each mode vs its exact/known dense limit + the default path unchanged, plus
# the -eigensolver propagation through fragments/dc. See docs/SQM_LARGE_SYSTEMS.md.
add_test(NAME large_system_modes
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/large_system_modes.sh
                 $<TARGET_FILE:curcuma> ${CMAKE_CURRENT_SOURCE_DIR}/molecules/larger
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(large_system_modes PROPERTIES TIMEOUT 180 LABELS "large_system;native_xtb")

# Claude Generated: CLI argument validation test suite (October 2025)
# Include comprehensive Bash-based CLI tests for all capabilities
add_subdirectory(cli)

# Claude Generated (Phase 0, Apr 2026): reference SCC dumps from tblite
# for native GFN1/GFN2 validation (see scripts/extract_xtb_params.py + xtb_native.*)
add_subdirectory(sqm_reference)
