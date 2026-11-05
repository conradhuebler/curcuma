# Seminumerical exchange (-qm.exchange_method cosx, COSX / RIJCOSX). Claude Generated (Nov 2026).
#
# check_qm_cosx.py runs curcuma -sp with the COSX option and compares with reference.json (PySCF RHF/RKS with the same
# COSX expression on curcuma's own grid, 1e-8 Eh); the grid error against the exact exchange is a measurement in
# docs/QM_COSX.md, the tests only bound it (--bound, about 10x the measured value).
#   qm_cosx_<case>            energy vs the same-grid PySCF COSX reference; |E_cosx - E_exact| <= bound
#   qm_cosx_fine_h2o          no overlap fit on a fine grid (75x302) approaches the exact exchange
#   qm_cosx_rij_h2o           RIJCOSX: COSX error with RI-J equals that with exact J (1e-6 Eh)
#   qm_cosx_threads_h2o       1 and 4 threads bitwise equal (dump_qm_energy)
#   qm_cosx_off_h2o           -qm.exchange_method exact = option absent, 17 digits incl. gradient (dump_qm_energy)
#   qm_cosx_refuse_*          gradient, open shell, -qm.ri jk, ECP, MP2 must fail loud
#
# Run:  ctest -L qm_cosx --output-on-failure
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_COSX ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_cosx.py
        --curcuma $<TARGET_FILE:curcuma> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    add_test(NAME qm_cosx_h2o_hf_svp COMMAND ${_COSX} --case h2o_hf_svp --bound 3e-5)
    add_test(NAME qm_cosx_nh3_b3lyp_svp COMMAND ${_COSX} --case nh3_b3lyp_svp --bound 3e-5)
    add_test(NAME qm_cosx_ch2o_hf_svp_nofit COMMAND ${_COSX} --case ch2o_hf_svp_nofit --bound 4e-3)
    add_test(NAME qm_cosx_ch4_hf_tzvp COMMAND ${_COSX} --case ch4_hf_tzvp --bound 3e-5)
    add_test(NAME qm_cosx_fine_h2o COMMAND ${_COSX} --case h2o_hf_svp --check fine --bound 1e-6)
    add_test(NAME qm_cosx_rij_h2o COMMAND ${_COSX} --case h2o_hf_svp --check rij --tol 1e-6)
    add_test(NAME qm_cosx_threads_h2o COMMAND ${_COSX} --case h2o_hf_svp --check threads --dump $<TARGET_FILE:dump_qm_energy>)
    add_test(NAME qm_cosx_off_h2o COMMAND ${_COSX} --case h2o_hf_svp --check off --dump $<TARGET_FILE:dump_qm_energy>)
    add_test(NAME qm_cosx_refuse_gradient COMMAND ${_COSX} --case h2o_hf_svp --check refuse "--extra=-gradient"
        "--expect=analytic gradient not implemented for -qm.exchange_method cosx")
    add_test(NAME qm_cosx_refuse_open_shell COMMAND ${_COSX} --case h2o_hf_svp --check refuse "--extra=-charge 1 -spin 1"
        "--expect=closed shell only")
    add_test(NAME qm_cosx_refuse_rijk COMMAND ${_COSX} --case h2o_hf_svp --check refuse "--extra=-qm.ri jk"
        "--expect=-qm.ri jk fits the exchange already")
    add_test(NAME qm_cosx_refuse_ecp COMMAND ${_COSX} --case h2o_hf_svp --check refuse --xyz ../structures/small/hydrogen-iodide.xyz
        "--expect=effective core potentials")
    add_test(NAME qm_cosx_refuse_mp2 COMMAND ${_COSX} --case h2o_hf_svp --check refuse "--extra=-method mp2"
        "--expect=MP2 is not implemented with -qm.exchange_method cosx")
    set_tests_properties(qm_cosx_h2o_hf_svp qm_cosx_nh3_b3lyp_svp qm_cosx_ch2o_hf_svp_nofit qm_cosx_ch4_hf_tzvp
        qm_cosx_fine_h2o qm_cosx_rij_h2o qm_cosx_threads_h2o qm_cosx_off_h2o qm_cosx_refuse_gradient
        qm_cosx_refuse_open_shell qm_cosx_refuse_rijk qm_cosx_refuse_ecp qm_cosx_refuse_mp2
        PROPERTIES TIMEOUT 600 LABELS "qm_cosx")
endif()
