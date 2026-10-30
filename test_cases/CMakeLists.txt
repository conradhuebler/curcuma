# Structures come from the library test_cases/structures (rules in its README), staged byte for byte into the build tree.
include(${CMAKE_CURRENT_SOURCE_DIR}/structures/structures.cmake)
curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR}
    input_aa=aaa-host.rmsd-a input_ab=aaa-host.rmsd-b
    A=aaa-bglc.conf76 B=aaa-bglc.conf96)

# The former molecules/ directories, filled from the library under their old file names
curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR}/molecules/dimers
    HH=dihydrogen.r0.47 HCl=hydrogen-chloride OH=hydroxyl)
curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR}/molecules/trimers
    HCN=hydrogen-cyanide O3=ozone water=water)
curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR}/molecules/ionic NaCl_cluster=na2cl2-isolated-ions)
curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR}/molecules/larger
    acetic_acid_dimer=acetic-acid-dimer C6H6=benzene C6H5COOH=benzoic-acid caffeine=caffeine
    CH3OCH3=dimethyl-ether CH3OH=methanol CH4=methane complex=macrocycle-bgal H2O=water.alt1
    HCl=hydrogen-chloride.alt1 mixture2=urea400-water1000 polymer=peo201-chain
    polymer_2x=peo201x2-water1500 polymer_2x_gfn2_opt=peo201x2-water1500.gfn2-relaxed
    polymer_2x_gfnff_opt=peo201x2-water1500.gfnff-min triose=trisaccharide-c18h32o16
    water8_cluster=water8-cluster)
curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR}/molecules/staged cl2m_ea25=cl2-anion)
# sqm_reference/molecules, same file names as before
curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR}/sqm_reference/molecules
    acetic_acid_dimer=acetic-acid-dimer C6H6=benzene.ideal-d6h caffeine=caffeine caffeine.opt=caffeine.opt
    CH4=methane.ideal-td complex=macrocycle-bgal FI_pyr=fi-pyridine H2=dihydrogen.r0.741 H2O=water.ideal-c2v
    H2S=hydrogen-sulfide.d-shell-validation HCl=hydrogen-chloride.d-shell-validation HCN=hydrogen-cyanide
    He2=helium-dimer.r3.0 LiH=lithium-hydride.r1.595 NH3=ammonia.ideal-c3v PH3=phosphine.d-shell-validation
    SiH4=silane.d-shell-validation triose=trisaccharide-c18h32o16)

# Claude Generated: Test Molecule Registry - Central hardcoded molecule database
# Provides zero-dependency molecule access for all test cases
add_library(test_molecule_registry unit/test_molecule_registry.cpp)

target_link_libraries(test_molecule_registry curcuma_core)
# the registry reads its geometries from the structure library (test_cases/structures)
target_compile_definitions(test_molecule_registry PRIVATE
    CURCUMA_STRUCTURE_LIBRARY_DIR="${CURCUMA_STRUCTURES_DIR}")

if(USE_PCH)
    target_precompile_headers(test_molecule_registry REUSE_FROM curcuma_core)
endif()

add_executable(rmsd_test
        rmsd/main.cpp)

target_link_libraries(rmsd_test curcuma_core)

add_executable(reorder_test
        reorder/main.cpp)

    add_executable(AAAbGlc
            unit/aaabglc_test.cpp)
# Claude Generated: Test executable for Molecule class
add_executable(test_molecule
        unit/test_molecule.cpp)

target_link_libraries(AAAbGlc curcuma_core)
target_link_libraries(reorder_test curcuma_core)
target_link_libraries(test_molecule curcuma_core)
if(USE_PCH)
        target_precompile_headers(AAAbGlc REUSE_FROM curcuma_core)
        target_precompile_headers(reorder_test REUSE_FROM curcuma_core)
        target_precompile_headers(test_molecule REUSE_FROM curcuma_core)
endif()

curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR} input=c33h45n9-agal-frames44)

add_executable(confscan_test
        unit/confscan_test.cpp)

target_link_libraries(confscan_test curcuma_core)

#target_link_libraries(curcuma_tests curcuma_core)￼

add_executable(confstat_tests
    unit/test_confstat.cpp
)

target_link_libraries(confstat_tests
        curcuma_core
)

# Claude Generated: Energy calculation test suite
add_executable(test_energy_methods
        unit/test_energy_methods.cpp)

target_link_libraries(test_energy_methods curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_energy_methods REUSE_FROM curcuma_core)
endif()

# Claude Generated 2026: BMT Utils unit test suite
add_executable(test_bmt_utils
        unit/test_bmt_utils.cpp)

target_link_libraries(test_bmt_utils curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_bmt_utils REUSE_FROM curcuma_core)
endif()

add_test(NAME test_bmt_utils COMMAND test_bmt_utils)
set_tests_properties(test_bmt_utils PROPERTIES TIMEOUT 10 LABELS "bmt;unit")

# Claude Generated 2026: PLUMED init unit test suite (only built with USE_Plumed)
if(USE_Plumed)
    add_executable(test_plumed_init
            unit/test_plumed_init.cpp)
    target_link_libraries(test_plumed_init curcuma_core)
    if(USE_PCH)
        target_precompile_headers(test_plumed_init REUSE_FROM curcuma_core)
    endif()
    add_test(NAME test_plumed_init COMMAND test_plumed_init)
    set_tests_properties(test_plumed_init PROPERTIES TIMEOUT 10 LABELS "plumed;unit")
endif()

# Claude Generated: CG Potentials unit test suite
add_executable(test_cg_potentials
        unit/test_cg_potentials.cpp)

target_link_libraries(test_cg_potentials curcuma_core)
add_test(NAME cg_potentials COMMAND test_cg_potentials)

if(USE_PCH)
        target_precompile_headers(test_cg_potentials REUSE_FROM curcuma_core)
endif()

# Note: test_cg_forcefield_integration removed due to complex linking with curcuma_D3
# CG functionality will be tested via CLI integration test instead

# Claude Generated (February 2026): Unified GFN-FF Validation Runner
# Purpose: High-precision validation against standardized analyzer reference data
add_executable(test_gfnff_validation
        unit/test_gfnff_validation.cpp)

target_link_libraries(test_gfnff_validation curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_validation REUSE_FROM curcuma_core)
endif()

# Claude Generated (June 2026): per-frame CPU-vs-GPU gradient diagnostic on a fixed
# trajectory — the clean force-backend metric (replaces the misleading MD heat-exchange).
# Builds unconditionally; the gpu config falls back to CPU when CUDA is not compiled in.
add_executable(test_gfnff_grad_traj unit/test_gfnff_grad_traj.cpp)
target_link_libraries(test_gfnff_grad_traj curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_gfnff_grad_traj REUSE_FROM curcuma_core)
endif()

# Claude Generated (Aug 2026): analytic vs finite-difference gradient across a
# react-topology rebuild event (topology_mode=react bond formation).
add_executable(test_gfnff_react_fd unit/test_gfnff_react_fd.cpp)
target_link_libraries(test_gfnff_react_fd curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_gfnff_react_fd REUSE_FROM curcuma_core)
endif()
add_test(NAME gfnff_react_fd_gradient
         COMMAND test_gfnff_react_fd
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Claude Generated (Sep 2026): GeometryTools angle/dihedral/radius of gyration.
# External potentials: parsing and selection resolution. Claude Generated 2026.
add_executable(test_external_potentials test_external_potentials.cpp)
target_link_libraries(test_external_potentials curcuma_cap curcuma_core)
add_test(NAME external_potentials COMMAND test_external_potentials)
set_tests_properties(external_potentials PROPERTIES TIMEOUT 300)

add_executable(test_geometry_tools test_geometry_tools.cpp)
target_link_libraries(test_geometry_tools curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_geometry_tools REUSE_FROM curcuma_core)
endif()
add_test(NAME geometry_tools COMMAND test_geometry_tools)
set_tests_properties(geometry_tools PROPERTIES TIMEOUT 60)

# Claude Generated (Sep 2026): react-scan bookkeeping — refractory count, authoritative
# empty bond set, exchange-resolution sigma counts, event record, bond orders.
add_executable(test_gfnff_react_filters test_gfnff_react_filters.cpp)
target_link_libraries(test_gfnff_react_filters curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_gfnff_react_filters REUSE_FROM curcuma_core)
endif()
add_test(NAME gfnff_react_filters
         COMMAND test_gfnff_react_filters
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(gfnff_react_filters PROPERTIES TIMEOUT 120 LABELS "gfnff;react")
# Claude Generated (Sep 27, 2026): frag_charge_model ensemble must not depend on the geometry
# the topology was built at (Known Issue #31; check 5 of cli_gfnff_05_frag_charge_ensemble,
# which needs a reused calculator instance the -sp CLI cannot express).
add_executable(test_gfnff_frag_charge_history unit/test_gfnff_frag_charge_history.cpp)
target_link_libraries(test_gfnff_frag_charge_history curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_gfnff_frag_charge_history REUSE_FROM curcuma_core)
endif()
add_test(NAME gfnff_frag_charge_history
         COMMAND test_gfnff_frag_charge_history
                 ${CMAKE_CURRENT_BINARY_DIR}/molecules/staged/cl2m_ea25.xyz
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(gfnff_frag_charge_history PROPERTIES TIMEOUT 120 LABELS "gfnff;frag_charge")

# Claude Generated (Oct 2026): the bounded placement selection of frag_charge_model ensemble
# (GFNFF::fragSelectPlacements) returns the same placements in the same order as the former full
# enumeration (20000 random cases), and scales to hundreds of fragment groups (docs/FRAG_CHARGE_MODEL.md).
add_executable(test_frag_charge_selection unit/test_frag_charge_selection.cpp)
target_link_libraries(test_frag_charge_selection curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_frag_charge_selection REUSE_FROM curcuma_core)
endif()
add_test(NAME gfnff_frag_charge_selection COMMAND test_frag_charge_selection 20000)
set_tests_properties(gfnff_frag_charge_selection PROPERTIES TIMEOUT 300 LABELS "gfnff;frag_charge")

# Claude Generated (Oct 2026): -gpu backend vs EEQ solver names (rocm: auto/cholesky only; vulkan: all) - name logic, no device.
add_executable(test_eeq_gpu_supports unit/test_eeq_gpu_supports.cpp)
target_link_libraries(test_eeq_gpu_supports curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_eeq_gpu_supports REUSE_FROM curcuma_core)
endif()
add_test(NAME eeq_gpu_supports COMMAND test_eeq_gpu_supports)
set_tests_properties(eeq_gpu_supports PROPERTIES TIMEOUT 60 LABELS "gfnff;eeq")

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
        unit/test_periodic_table.cpp)

target_link_libraries(test_periodic_table curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_periodic_table REUSE_FROM curcuma_core)
endif()

add_test(NAME test_periodic_table
         COMMAND test_periodic_table)

# Claude Generated (June 2026, Phase 0 ORCA improvement plan):
# Unit tests for ORCA interface helpers, parsers, and CG-atom guard.
add_executable(test_orca_interface
        unit/test_orca_interface.cpp)

target_link_libraries(test_orca_interface curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_orca_interface REUSE_FROM curcuma_core)
endif()

add_test(NAME test_orca_interface
         COMMAND test_orca_interface)
set_tests_properties(test_orca_interface PROPERTIES TIMEOUT 30 LABELS "orca;unit")

# Claude Generated (May 2026, WP-D Stage B): fast_exp_neg_sq_block unit test
add_executable(test_fast_exp
        unit/test_fast_exp.cpp)

target_link_libraries(test_fast_exp curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_fast_exp REUSE_FROM curcuma_core)
endif()

# Claude Generated (June 2026): opt-in multi-step EEQ PCG warm-start extrapolation.
# Drives GFN-FF (forced PCG) over a smooth trajectory of a multi-fragment system;
# asserts aspc/gauss energies == none (PCG result-neutral) and that the extrapolation
# path fired. Run via: ctest -L eeq_extrapolation
add_executable(test_eeq_extrapolation
        unit/test_eeq_extrapolation.cpp)
target_link_libraries(test_eeq_extrapolation curcuma_core)
if(USE_PCH)
        target_precompile_headers(test_eeq_extrapolation REUSE_FROM curcuma_core)
endif()
add_test(NAME test_eeq_extrapolation
         COMMAND test_eeq_extrapolation
                 ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/complex.xyz 12)
set_tests_properties(test_eeq_extrapolation PROPERTIES TIMEOUT 180 LABELS "eeq_extrapolation")

add_test(NAME test_fast_exp
         COMMAND test_fast_exp)

# Claude Generated (Sep 2026): vendored fdlibm erf/acos/exp/log (CURCUMA_PORTABLE_MATH
# layer, docs/PORTABLE_ERF.md) pinned against the host libm + special values. Header-only
# functions, so no library link is needed. Run via: ctest -R test_portable_math
add_executable(test_portable_math unit/test_portable_math.cpp)
target_include_directories(test_portable_math PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME test_portable_math COMMAND test_portable_math)
set_tests_properties(test_portable_math PROPERTIES LABELS "portable_math;unit")

# Claude Generated (Session 5): Two-Phase EEQ System Tests
# Tests Phase 1 topology charges and Phase 2 environment corrections
add_executable(test_eeq_two_phase
        unit/test_eeq_two_phase.cpp)

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
        unit/test_gfnff_topology_reentrancy.cpp)

target_link_libraries(test_gfnff_topology_reentrancy curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_gfnff_topology_reentrancy REUSE_FROM curcuma_core)
endif()

add_test(NAME test_gfnff_topology_reentrancy
         COMMAND test_gfnff_topology_reentrancy ${CMAKE_CURRENT_BINARY_DIR})

# Claude Generated (Sep 2026, Known Issue #32): energy-only calls on a reused GFN-FF
# calculator must use their own CN (Coulomb chi(CN)) and D4 C6, not the last gradient
# call's. Checks FD-vs-analytic via NumGrad (reused energy-only calls), energy-only vs
# gradient call at one geometry, and reused vs fresh Coulomb / total / D4 dispersion
# energies (the last two fail when only the CN fix, not the C6 refresh, is present).
add_executable(test_gfnff_stale_cn unit/test_gfnff_stale_cn.cpp)
target_link_libraries(test_gfnff_stale_cn curcuma_core test_molecule_registry)
if(USE_PCH)
        target_precompile_headers(test_gfnff_stale_cn REUSE_FROM curcuma_core)
endif()
add_test(NAME gfnff_stale_cn_energy_only
         COMMAND test_gfnff_stale_cn
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(gfnff_stale_cn_energy_only
    PROPERTIES TIMEOUT 300 LABELS "gfnff;stale_cn")

# Claude Generated: PolymerBuild sequence parser tests (prime notation)
add_executable(test_polymerbuild_sequence
        unit/test_polymerbuild_sequence.cpp)

target_link_libraries(test_polymerbuild_sequence curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_polymerbuild_sequence REUSE_FROM curcuma_core)
endif()

add_test(NAME test_polymerbuild_sequence
         COMMAND test_polymerbuild_sequence)

# Claude Generated: PolymerBuild structural tests (Xx selection, geometry comparison)
curcuma_stage_structures(${CMAKE_CURRENT_BINARY_DIR}/polymerbuild fragment_2xx=polymerbuild-fragment-2xx)

add_executable(test_polymerbuild_structure
        unit/test_polymerbuild_structure.cpp)

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
add_executable(test_gfnff_numgrad unit/test_gfnff_numgrad.cpp)
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
add_executable(test_xtb_solvation_numgrad unit/test_xtb_solvation_numgrad.cpp)
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
    add_executable(test_gfnff_gpu unit/test_gfnff_gpu.cpp)
    # Sep 2026: USE_CUDA no longer propagates from curcuma_core (the backend is a
    # plugin); the test links the plugin directly for the GFNFFGPUComputationalMethod
    # typeinfo it dynamic_casts to.
    target_link_libraries(test_gfnff_gpu curcuma_core curcuma_cuda CUDA::cudart)
    target_compile_definitions(test_gfnff_gpu PRIVATE USE_CUDA)

    if(USE_PCH)
        target_precompile_headers(test_gfnff_gpu REUSE_FROM curcuma_core)
    endif()

    # Built-in molecule test (water, methane, ethanol) — no arguments
    add_test(NAME gfnff_gpu_builtin_molecules
             COMMAND test_gfnff_gpu
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(gfnff_gpu_builtin_molecules
        PROPERTIES TIMEOUT 120 RESOURCE_LOCK "gpu" LABELS "gfnff_gpu;gpu;validation")

    # Per-molecule tests using the same reference data as CPU GFN-FF validation
    foreach(REF_FILE ${GFNFF_REF_FILES})
        get_filename_component(MOL_NAME ${REF_FILE} NAME_WE)
        set(BUILD_REF_FILE "${CMAKE_CURRENT_BINARY_DIR}/reference_data/${MOL_NAME}.ref.json")
        set(GPU_TEST_NAME "gfnff_gpu_vs_cpu_${MOL_NAME}")

        add_test(NAME ${GPU_TEST_NAME}
                 COMMAND test_gfnff_gpu ${BUILD_REF_FILE}
                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(${GPU_TEST_NAME}
            PROPERTIES TIMEOUT 120 RESOURCE_LOCK "gpu" LABELS "gfnff_gpu;gpu;validation")
    endforeach()

    # Claude Generated (March 2026): GPU analytical vs numerical gradient test
    add_executable(test_gpu_numgrad unit/test_gpu_numgrad.cpp)
    target_link_libraries(test_gpu_numgrad curcuma_core CUDA::cudart)
    target_compile_definitions(test_gpu_numgrad PRIVATE USE_CUDA)
    if(USE_PCH)
        target_precompile_headers(test_gpu_numgrad REUSE_FROM curcuma_core)
    endif()

endif()

# Claude Generated (Sep 2026, MULTI_GPU_GAPS F-11): an energy-only call at a new geometry must
# equal a gradient call there (CPU: stale Coulomb CN; GPU: stale D4 C6).
add_executable(test_gfnff_energy_only_geometry unit/test_gfnff_energy_only_geometry.cpp)
target_link_libraries(test_gfnff_energy_only_geometry curcuma_core)
# The GPU plugin is looked up next to the executable (gpu_plugin.cpp, /proc/self/exe).
set_target_properties(test_gfnff_energy_only_geometry PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
if(USE_PCH)
    target_precompile_headers(test_gfnff_energy_only_geometry REUSE_FROM curcuma_core)
endif()
foreach(EO_MOL triose caffeine)
    add_test(NAME gfnff_energy_only_geometry_${EO_MOL}
             COMMAND test_gfnff_energy_only_geometry
                     ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/${EO_MOL}.xyz none
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(gfnff_energy_only_geometry_${EO_MOL}
        PROPERTIES TIMEOUT 120 LABELS "gfnff")
    if(USE_CUDA)
        add_test(NAME gfnff_gpu_energy_only_geometry_${EO_MOL}
                 COMMAND test_gfnff_energy_only_geometry
                         ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/${EO_MOL}.xyz cuda
                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(gfnff_gpu_energy_only_geometry_${EO_MOL}
            PROPERTIES TIMEOUT 120 RESOURCE_LOCK "gpu" LABELS "gfnff_gpu;gpu")
    endif()
    # Claude Generated (Oct 2026): the same check on ROCm, with the fused dispersion gather (default) and with
    # CURCUMA_ROCM_DISP_FUSED=0 (per-pair C6/dC6dCN kernel). The fused path leaves the per-pair C6 stale in gradient
    # calls and refills it before the next energy-only call; this test is the only one that runs that refill.
    if(USE_ROCM)
        add_test(NAME gfnff_rocm_energy_only_geometry_${EO_MOL}
                 COMMAND test_gfnff_energy_only_geometry
                         ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/${EO_MOL}.xyz rocm
                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(gfnff_rocm_energy_only_geometry_${EO_MOL}
            PROPERTIES TIMEOUT 120 RESOURCE_LOCK "gpu" LABELS "gfnff_gpu;gpu;rocm")
        add_test(NAME gfnff_rocm_energy_only_geometry_unfused_${EO_MOL}
                 COMMAND test_gfnff_energy_only_geometry
                         ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/${EO_MOL}.xyz rocm
                 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(gfnff_rocm_energy_only_geometry_unfused_${EO_MOL}
            PROPERTIES TIMEOUT 120 RESOURCE_LOCK "gpu" LABELS "gfnff_gpu;gpu;rocm"
                       ENVIRONMENT "CURCUMA_ROCM_DISP_FUSED=0")
    endif()
endforeach()


# Claude Generated (AP5, Apr 2026): GFN1/GFN2 analytical vs numerical gradient test
add_executable(test_xtb_gradient unit/test_xtb_gradient.cpp)
target_link_libraries(test_xtb_gradient curcuma_core)
if(USE_PCH)
    target_precompile_headers(test_xtb_gradient REUSE_FROM curcuma_core)
endif()

set(_GRAD_MOL_DIR ${CMAKE_CURRENT_BINARY_DIR}/sqm_reference/molecules)
foreach(_mol H2O CH4 NH3)
    add_test(NAME xtb_gradient_${_mol}
             COMMAND test_xtb_gradient ${_GRAD_MOL_DIR}/${_mol}.xyz
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(xtb_gradient_${_mol}
        PROPERTIES TIMEOUT 120 LABELS "gradient;native_xtb")
endforeach()

# Claude Generated (WP4, Jun 2026): native self-contained symmetric eigensolver
# (the `-eigensolver native` backend) vs Eigen on random symmetric matrices.
add_executable(test_native_eigensolver unit/test_native_eigensolver.cpp)
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
add_executable(test_purification unit/test_purification.cpp)
target_link_libraries(test_purification curcuma_core)
if(USE_PCH)
    target_precompile_headers(test_purification REUSE_FROM curcuma_core)
endif()
add_test(NAME purification COMMAND test_purification
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(purification PROPERTIES TIMEOUT 60 LABELS "native_xtb;eigensolver;purification")

# Claude Generated (Jun 2026): seeded block LOBPCG kernel (the iterative `-eigensolver lobpcg`
# path) vs Eigen's lowest-k eigenpairs, incl. a warm-start (subspace-recycling) case.
# The measurement capability. Claude Generated 2026 - checked against geometry
# worked out by hand, not against a previous run of the same code.
# The CIF reader and supercell builder. Claude Generated 2026 - the face-centring
# case is checked against positions worked out on paper, because a reader that
# ignores symmetry returns a well-formed molecule with the wrong atom count.
add_executable(test_cif test_cif.cpp)
target_link_libraries(test_cif curcuma_core)
add_test(NAME cif COMMAND test_cif)
set_tests_properties(cif PROPERTIES TIMEOUT 60)

add_executable(test_measurement test_measurement.cpp)
target_link_libraries(test_measurement curcuma_cap curcuma_core)
add_test(NAME measurement_capability COMMAND test_measurement)
set_tests_properties(measurement_capability PROPERTIES TIMEOUT 120)

add_executable(test_lobpcg unit/test_lobpcg.cpp)
target_link_libraries(test_lobpcg curcuma_core)
if(USE_PCH)
    target_precompile_headers(test_lobpcg REUSE_FROM curcuma_core)
endif()
add_test(NAME lobpcg COMMAND test_lobpcg
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(lobpcg PROPERTIES TIMEOUT 60 LABELS "native_xtb;eigensolver;lobpcg")

# Claude Generated (2026, AP ∂q/∂x — Phase 1): D4 dE/dq analytical-vs-FD test
add_executable(test_d4_dedq unit/test_d4_dedq.cpp)
target_link_libraries(test_d4_dedq curcuma_core test_molecule_registry)
if(USE_PCH)
    target_precompile_headers(test_d4_dedq REUSE_FROM curcuma_core)
endif()
add_test(NAME d4_dedq COMMAND test_d4_dedq
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(d4_dedq PROPERTIES TIMEOUT 120 LABELS "gradient;dispersion")

# Claude Generated (Phase 3b-2, May 2026): CPSCF orbital-Hessian + Z-vector gate test
add_executable(test_xtb_cpscf unit/test_xtb_cpscf.cpp)
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
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/scripts/large_system_modes.sh
                 $<TARGET_FILE:curcuma> ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(large_system_modes PROPERTIES TIMEOUT 180 LABELS "large_system;native_xtb")

# Claude Generated: CLI argument validation test suite (October 2025)
# Include comprehensive Bash-based CLI tests for all capabilities
add_subdirectory(cli)

# Claude Generated (Phase 0, Apr 2026): reference SCC dumps from tblite
# for native GFN1/GFN2 validation (see scripts/extract_xtb_params.py + xtb_native.*)
add_subdirectory(sqm_reference)

# Claude Generated 2026: Path collective variables (Branduardi s,z) validation
add_executable(test_path_cv
        test_path_cv.cpp)

target_link_libraries(test_path_cv curcuma_core)

if(USE_PCH)
        target_precompile_headers(test_path_cv REUSE_FROM curcuma_core)
endif()

add_test(NAME test_path_cv COMMAND test_path_cv)
set_tests_properties(test_path_cv PROPERTIES TIMEOUT 60 LABELS "pathcv;unit")

# Claude Generated 2026: WHAM solver vs a known double-well profile (header-only, no core link)
add_executable(test_wham test_wham.cpp)
target_include_directories(test_wham PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME test_wham COMMAND test_wham)
set_tests_properties(test_wham PROPERTIES TIMEOUT 120 LABELS "pathcv;unit")

# Claude Generated 2026: SHAKE/RATTLE scalar constraint, NVE energy conservation (distance CV)
add_executable(test_scalar_constraint test_scalar_constraint.cpp)
target_link_libraries(test_scalar_constraint curcuma_cap curcuma_core)
target_compile_definitions(test_scalar_constraint PRIVATE CURCUMA_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
add_test(NAME test_scalar_constraint COMMAND test_scalar_constraint)
set_tests_properties(test_scalar_constraint PROPERTIES TIMEOUT 300 LABELS "pathcv;unit")

# Claude Generated 2026: SimpleMD umbrella_q restraint (gradient vs FD, parser)
add_executable(test_umbrella_q test_umbrella_q.cpp)
target_link_libraries(test_umbrella_q curcuma_cap curcuma_core)
add_test(NAME test_simplemd_umbrella_q COMMAND test_umbrella_q)
set_tests_properties(test_simplemd_umbrella_q PROPERTIES TIMEOUT 60 LABELS "pathcv;unit")

# Claude Generated 2026: pure string/path helpers (Kabsch, reparametrisation, double-well FTS)
add_executable(test_string_path test_string_path.cpp)
target_link_libraries(test_string_path curcuma_core)
add_test(NAME test_string_path COMMAND test_string_path)
set_tests_properties(test_string_path PROPERTIES TIMEOUT 60 LABELS "pathcv;unit")

# Claude Generated 2026: automatic z threshold (population separation) for -nebmd.cv_histogram
add_executable(test_z_filter test_z_filter.cpp)
target_include_directories(test_z_filter PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME test_z_filter COMMAND test_z_filter)
set_tests_properties(test_z_filter PROPERTIES TIMEOUT 60 LABELS "pathcv;unit")

# Claude Generated 2026: grid-free (natural cubic spline) barrier of the mean-force profile
add_executable(test_fe_profile test_fe_profile.cpp)
target_include_directories(test_fe_profile PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME test_fe_profile COMMAND test_fe_profile)
set_tests_properties(test_fe_profile PROPERTIES TIMEOUT 60 LABELS "pathcv;unit")

# Claude Generated 2026: absolute TST rate / Eyring free energy across s* (fe_profile.h tstRate)
add_executable(test_tst_rate test_tst_rate.cpp)
target_include_directories(test_tst_rate PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME test_tst_rate COMMAND test_tst_rate)
set_tests_properties(test_tst_rate PROPERTIES TIMEOUT 60 LABELS "pathcv;unit")

# Claude Generated (Oct 2026): property tests generated from the structure library (test_cases/struct_test_policy.json):
# translation/rotation/permutation invariance, sum of forces and torque, finite-difference gradient, repeat and thread
# independence, CPU against CUDA. The selection is made by scripts/gen_struct_tests.py at configure time.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(CURCUMA_STRUCT_PYTHON "${Python3_EXECUTABLE}")
    get_filename_component(CURCUMA_STRUCT_REPO "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
    execute_process(COMMAND "${Python3_EXECUTABLE}" "${CURCUMA_STRUCT_REPO}/scripts/gen_struct_tests.py"
                            --out "${CMAKE_CURRENT_BINARY_DIR}/generated_struct_tests.cmake" --cuda "${USE_CUDA}" --rocm "${USE_ROCM}"
                    RESULT_VARIABLE _struct_gen_rc ERROR_VARIABLE _struct_gen_err)
    if(_struct_gen_rc EQUAL 0)
        include("${CMAKE_CURRENT_BINARY_DIR}/generated_struct_tests.cmake")
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${CURCUMA_STRUCTURES_DIR}/manifest.json" "${CMAKE_CURRENT_SOURCE_DIR}/struct_test_policy.json"
            "${CURCUMA_STRUCT_REPO}/scripts/gen_struct_tests.py")
    else()
        message(WARNING "structure property tests not generated: ${_struct_gen_err}")
    endif()
else()
    message(STATUS "Python3 not found - structure property tests are not generated")
endif()

# Claude Generated (Sep 2026): interface-level gradient + frequency contract.
# Guards two things nothing else covered:
#   1. EnergyCalculator::Gradient() must be Eh/Angstrom for EVERY method — SimpleMD
#      integrates it against Angstrom coordinates and the FD Hessian differences it
#      against Angstrom displacements. GFN-FF returned Eh/Bohr here until Sep 2026;
#      gfnff_numgrad_builtin could not see it because it compares GFN-FF's INTERNAL
#      gradient against an internal FD, both in Bohr.
#   2. Vibrational frequencies vs committed xtb 6.7.1 values (no xtb needed at run
#      time). Catches the Eh/Ang^2-vs-Eh/Bohr^2 Hessian bug and the empirical
#      "+47.349 cm-1" offset that used to mask it.
# Verified to FAIL on a pre-fix binary (gfnff gradient ratio 0.5287 = au; frequencies
# 92 % / 92 % / 40 % high for gfn1 / gfn2 / gfnff).
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME gradient_unit_contract
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_gradient_units.py
                     $<TARGET_FILE:curcuma>)
    set_tests_properties(gradient_unit_contract
        PROPERTIES TIMEOUT 300 LABELS "gradient;frequencies;native_xtb;gfnff")

    # Claude Generated (Sep 2026, docs/MULTI_GPU_GAPS.md F-3): the GFN-FF multi-GPU split of one
    # molecule (Coulomb tiles + projected-PCG EEQ) must reproduce the single-device run.
    # Skipped (77) with fewer than two CUDA devices.
    if(USE_CUDA)
        add_test(NAME gfnff_gpu_split_equals_single
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_gfnff_gpu_split.py
                         $<TARGET_FILE:curcuma>
                         ${CMAKE_CURRENT_BINARY_DIR}/molecules/larger/polymer_2x.xyz)
        set_tests_properties(gfnff_gpu_split_equals_single
            PROPERTIES TIMEOUT 600 RESOURCE_LOCK "gpu" LABELS "gfnff_gpu;gpu" SKIP_RETURN_CODE 77)
    endif()

    # Claude Generated (Sep 2026): adaptive (step-rejecting) integrator.
    #   1. adaptive_step is OFF by default and an explicit false is identical to it.
    #   2. The failure mode it exists for reproduces on ONE water: 1750 K, dt = 1 fs, NVE
    #      runs away by +8.78 Eh because the GFN-FF O-H bond stiffens under compression.
    #   3. With the feature the same run conserves energy (-0.003 Eh), with no constraint
    #      and no mass modification.
    # See docs/MD_LARGE_SYSTEMS.md.
    add_test(NAME md_adaptive_step
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_adaptive_step.py
                     $<TARGET_FILE:curcuma>)
    set_tests_properties(md_adaptive_step
        PROPERTIES TIMEOUT 300 LABELS "md;integrator;gfnff")

    # Claude Generated (Sep 2026): the MD clock. SimpleMD integrates in Angstrom/amu/Hartree,
    # whose time unit is 1.9516 fs, so a step given in femtoseconds has to be converted before
    # it multiplies a velocity. Until Sep 2026 it was not, and every MD ran 1.9516x too fast
    # while energies, forces, temperature and energy conservation all stayed correct - which is
    # why nothing noticed. This test ties the reported time to a vibrational period taken from
    # the (independently validated) Hessian path. See docs/MD_LARGE_SYSTEMS.md.
    add_test(NAME md_time_axis
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_md_time_axis.py
                     $<TARGET_FILE:curcuma>)
    set_tests_properties(md_time_axis
        PROPERTIES TIMEOUT 300 LABELS "md;integrator;gfnff;units")

    # Claude Generated (Oct 2026): CPU curcuma single points (gfnff, gfn2) against the stored
    # xtb 6.7.1 energies of the structure library. No xtb run; see scripts/check_xtb_library_energies.py.
    add_test(NAME xtb_library_energies
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_xtb_library_energies.py
                     --curcuma $<TARGET_FILE:curcuma>
                     --reference ${CMAKE_CURRENT_SOURCE_DIR}/reference_data/xtb_library_energies.json
                     --structures ${CMAKE_CURRENT_SOURCE_DIR}/structures)
    set_tests_properties(xtb_library_energies
        PROPERTIES TIMEOUT 900 LABELS "gfnff;sqm;xtb_reference")
endif()
