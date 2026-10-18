# Native local PNO-RI-MP2 (WP15, method lmp2). Claude Generated (Oct 2026).
#
# check_lmp2.py runs `curcuma -sp <xyz> -method lmp2 -verbosity 2` (def2-SVP, aux def2-universal-JKFIT).
#   qm_lmp2_exact_*     complete pair spaces, no PNO cut, all pairs and Fock couplings, -qm.lmp2_check true:
#                       the iterated LMP2 must equal curcuma's canonical RI-MP2 (same aux) within 1e-9 Eh
#                       (RI-MP2 itself is checked against PySCF in ctest -L qm_mp2)
#   qm_lmp2_default_*   default thresholds against reference.json, produced by scripts/lmp2_reference.py with
#                       curcuma itself (regression reference, tolerance 1e-7 Eh; not an external validation)
#   qm_lmp2_opt_error   -opt with lmp2 must stop with the "no analytic gradient" error
#   qm_lmp2_mem_error   a too small -qm.lmp2_max_memory_mb must stop with an error naming the size
# Structures from the library (test_cases/structures). Run:  ctest -L qm_lmp2 --output-on-failure

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_S ${CMAKE_SOURCE_DIR}/test_cases/structures)
    set(_EX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_lmp2.py --curcuma $<TARGET_FILE:curcuma> --exact)
    add_test(NAME qm_lmp2_exact_butane_fc COMMAND ${_EX} --fc --xyz ${_S}/medium/butane.anti.xyz)
    add_test(NAME qm_lmp2_exact_butane COMMAND ${_EX} --xyz ${_S}/medium/butane.anti.xyz)
    add_test(NAME qm_lmp2_exact_boys_water_dimer COMMAND ${_EX} --xyz ${_S}/clusters/water-dimer.roo2.8.xyz
        "--extra=-qm.lmp2_localizer boys")
    add_test(NAME qm_lmp2_exact_boys_benzene_fc COMMAND ${_EX} --fc --xyz ${_S}/small/benzene.ideal-d6h.xyz
        "--extra=-qm.lmp2_localizer boys")
    add_test(NAME qm_lmp2_exact_rijk_water_dimer_fc COMMAND ${_EX} --fc --xyz ${_S}/clusters/water-dimer.roo2.8.xyz
        "--extra=-qm.ri jk")
    add_test(NAME qm_lmp2_exact_jfit_water_dimer_fc COMMAND ${_EX} --fc --xyz ${_S}/clusters/water-dimer.roo2.8.xyz
        "--extra=-qm.ri_mp2_aux def2-universal-JFIT")
    add_test(NAME qm_lmp2_exact_threads1_acetic_acid_dimer_fc COMMAND ${_EX} --fc --threads 1
        --xyz ${_S}/clusters/acetic-acid-dimer.cyclic-hb.xyz)

    set(_DF ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_lmp2.py --curcuma $<TARGET_FILE:curcuma>
        --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    foreach(_c butane_fc butane_boys_fc benzene_fc acetic_acid_dimer_fc)
        add_test(NAME qm_lmp2_default_${_c} COMMAND ${_DF} --case ${_c})
    endforeach()
    add_test(NAME qm_lmp2_default_threads1_butane_fc COMMAND ${_DF} --case butane_fc --threads 1)
    set_tests_properties(qm_lmp2_exact_butane_fc qm_lmp2_exact_butane qm_lmp2_exact_boys_water_dimer
        qm_lmp2_exact_boys_benzene_fc qm_lmp2_exact_rijk_water_dimer_fc qm_lmp2_exact_jfit_water_dimer_fc
        qm_lmp2_exact_threads1_acetic_acid_dimer_fc qm_lmp2_default_butane_fc qm_lmp2_default_butane_boys_fc
        qm_lmp2_default_benzene_fc qm_lmp2_default_acetic_acid_dimer_fc qm_lmp2_default_threads1_butane_fc
        PROPERTIES TIMEOUT 900 LABELS "qm_lmp2;qm")

    add_test(NAME qm_lmp2_opt_error
        COMMAND ${CMAKE_COMMAND} -E env sh -c
            "$<TARGET_FILE:curcuma> -opt ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e/H2O.xyz -method lmp2 -no_bmt 2>&1 | grep -q 'lmp2 has no analytic gradient'")
    add_test(NAME qm_lmp2_mem_error
        COMMAND ${CMAKE_COMMAND} -E env sh -c
            "$<TARGET_FILE:curcuma> -sp ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e/H2O.xyz -method lmp2 -qm.lmp2_max_memory_mb 0.001 -no_bmt 2>&1 | grep -q 'LMP2: the factors B'")
    set_tests_properties(qm_lmp2_opt_error qm_lmp2_mem_error PROPERTIES TIMEOUT 300 LABELS "qm_lmp2;qm")
endif()
