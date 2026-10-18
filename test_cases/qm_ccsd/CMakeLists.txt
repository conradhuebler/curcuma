# Native canonical closed-shell CCSD (WP13) vs PySCF. Claude Generated (Oct 2026).
#
# check_ccsd.py runs `curcuma -sp <xyz> -method ccsd -verbosity 2` and compares E(HF), the iteration-0
# (MP2) and the CCSD correlation energy, the total and the T1 diagnostic with reference.json, produced by
# scripts/ccsd_reference.py (PySCF 2.14 RHF conv_tol 1e-12 + cc.CCSD conv_tol 1e-12, def2-SVP spherical,
# frozen=None or the chemical frozen core). Tolerance 1e-8 Eh. ORCA 6.1 values (`! RHF CCSD def2-SVP
# ExtremeSCF [No]FrozenCore`, STol 1e-10) are stored in the same file and printed, not tested.
#   qm_ccsd_<mol>        all electrons,     qm_ccsd_<mol>_fc   frozen core
#   qm_ccsd_direct_*     integral-direct SCF: the AO tensor is rebuilt for CCSD
#   qm_ccsd_threads1_*   one thread
#   qm_ccsd_opt_error    -opt with ccsd must stop with the "no analytic gradient" error
#   qm_ccsd_ri_error     -qm.ri jk with ccsd must stop with an error
#   qm_ccsd_mem_error    a too small -qm.cc_max_memory_mb must stop with an error naming the size
#
# Run:  ctest -L qm_ccsd --output-on-failure

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_CC ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_ccsd.py
        --curcuma $<TARGET_FILE:curcuma> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    set(_CC_TESTS)
    foreach(_mol H2O NH3 CH4 formaldehyde ethane water_dimer)
        add_test(NAME qm_ccsd_${_mol} COMMAND ${_CC} --mol ${_mol})
        add_test(NAME qm_ccsd_${_mol}_fc COMMAND ${_CC} --mol ${_mol} --fc)
        list(APPEND _CC_TESTS qm_ccsd_${_mol} qm_ccsd_${_mol}_fc)
    endforeach()
    add_test(NAME qm_ccsd_benzene_fc COMMAND ${_CC} --mol benzene --fc --threads 8)
    add_test(NAME qm_ccsd_direct_formaldehyde COMMAND ${_CC} --mol formaldehyde "--extra=-qm.scf_direct on")
    add_test(NAME qm_ccsd_threads1_H2O_fc COMMAND ${_CC} --mol H2O --fc --threads 1)
    list(APPEND _CC_TESTS qm_ccsd_benzene_fc qm_ccsd_direct_formaldehyde qm_ccsd_threads1_H2O_fc)
    set_tests_properties(${_CC_TESTS} PROPERTIES TIMEOUT 900 LABELS "qm_ccsd;qm")

    add_test(NAME qm_ccsd_opt_error
        COMMAND ${CMAKE_COMMAND} -E env sh -c
            "$<TARGET_FILE:curcuma> -opt ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e/H2O.xyz -method ccsd -no_bmt 2>&1 | grep -q 'ccsd has no analytic gradient'")
    add_test(NAME qm_ccsd_ri_error
        COMMAND ${CMAKE_COMMAND} -E env sh -c
            "$<TARGET_FILE:curcuma> -sp ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e/H2O.xyz -method ccsd -qm.ri jk -no_bmt 2>&1 | grep -q 'is not supported with ccsd'")
    add_test(NAME qm_ccsd_mem_error
        COMMAND ${CMAKE_COMMAND} -E env sh -c
            "$<TARGET_FILE:curcuma> -sp ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e/H2O.xyz -method ccsd -qm.cc_max_memory_mb 1 -no_bmt 2>&1 | grep -q 'CCSD: needs about'")
    set_tests_properties(qm_ccsd_opt_error qm_ccsd_ri_error qm_ccsd_mem_error PROPERTIES TIMEOUT 300 LABELS "qm_ccsd;qm")
endif()
