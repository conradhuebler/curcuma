# Native meta-GGA analytic gradient (WP11 stage 2): tpss, scan on the CPU. Claude Generated (Oct 2026).
#
# Reuses the WP8 driver test_qm_dft_grad (target of ../qm_dft_grad, --functional tpss|scan) and its check script
# check_qm_dft_grad.py; reference.json here is written by
#   ~/.venvs/qmref/bin/python scripts/qm_dft_grad_pyscf_reference.py --set mgga   (PySCF 2.14, libxc 7.0.0, scipy BFGS)
#   qm_mgga_grad_pyscf_<case>   without / with grid response vs PySCF grid_response False / True, same grid, tol 1e-9
#   qm_mgga_grad_fd_<case>_<g>  analytic vs 4-point FD of curcuma's energy: moving grid (h 5e-4 Bohr) = with
#                               response, tol 2e-9; frozen grid (h 1e-4 Bohr) = without response, tol 1e-8 (at
#                               h 2e-3 the frozen-grid FD of the heavy atom's on-axis component is off by up to 1.4e-5,
#                               see docs/DFT_ROADMAP/WP11_meta_gga.md)
#   qm_mgga_grad_{translation,rotation,threads,update}_*   as in qm_dft_grad
#   qm_mgga_opt_<mol>_tpss      -opt reaches the PySCF/scipy minimum (sg1, grid response on)
# Geometries: xyz files of qm_1e, qm_hf3c and qm_grad (no new structure files).
# Run:  ctest -L qm_mgga_grad --output-on-failure

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_MG_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../qm_dft_grad/check_qm_dft_grad.py
        --exe $<TARGET_FILE:test_qm_dft_grad> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    set(_MG_PYSCF_CASES)
    foreach(_f tpss scan)
        foreach(_m H2O NH3 CH4 HF benzene)
            list(APPEND _MG_PYSCF_CASES ${_m}_${_f}_sg1)
        endforeach()
        list(APPEND _MG_PYSCF_CASES H2O_${_f}_sg3 H2O_${_f}_rij_sg1 NH3_${_f}_150x974)
    endforeach()
    list(APPEND _MG_PYSCF_CASES benzene_tpss_rijk_sg1)
    foreach(_case ${_MG_PYSCF_CASES})
        add_test(NAME qm_mgga_grad_pyscf_${_case} COMMAND ${_MG_CHECK} --check pyscf --case ${_case})
        set_tests_properties(qm_mgga_grad_pyscf_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_mgga_grad")
    endforeach()
    foreach(_case H2O_tpss H2O_scan NH3_tpss NH3_scan)
        add_test(NAME qm_mgga_grad_fd_${_case}_moving
                 COMMAND ${_MG_CHECK} --check fd --case ${_case} --fd moving --h 5e-4 --tol 2e-9)
        add_test(NAME qm_mgga_grad_fd_${_case}_frozen
                 COMMAND ${_MG_CHECK} --check fd --case ${_case} --fd frozen --h 1e-4 --tol 1e-8)
        set_tests_properties(qm_mgga_grad_fd_${_case}_moving qm_mgga_grad_fd_${_case}_frozen
                             PROPERTIES TIMEOUT 900 LABELS "qm_mgga_grad")
    endforeach()
    foreach(_case NH3_tpss H2O_scan)
        add_test(NAME qm_mgga_grad_translation_${_case} COMMAND ${_MG_CHECK} --check translation --case ${_case})
        set_tests_properties(qm_mgga_grad_translation_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_mgga_grad")
    endforeach()
    foreach(_case NH3_tpss H2O_tpss H2O_scan_sg3)
        add_test(NAME qm_mgga_grad_rotation_${_case} COMMAND ${_MG_CHECK} --check rotation --case ${_case})
        set_tests_properties(qm_mgga_grad_rotation_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_mgga_grad")
    endforeach()
    foreach(_case NH3_tpss NH3_scan benzene_tpss)
        add_test(NAME qm_mgga_grad_threads_${_case} COMMAND ${_MG_CHECK} --check threads --case ${_case})
        set_tests_properties(qm_mgga_grad_threads_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_mgga_grad")
    endforeach()
    foreach(_case H2O_tpss H2O_scan)
        add_test(NAME qm_mgga_grad_update_${_case} COMMAND ${_MG_CHECK} --check update --case ${_case})
        set_tests_properties(qm_mgga_grad_update_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_mgga_grad")
    endforeach()
    foreach(_case H2O_distorted_tpss NH3_tpss formaldehyde_distorted_tpss)
        add_test(NAME qm_mgga_opt_${_case}
            COMMAND ${_MG_CHECK} --check opt --case ${_case} --curcuma $<TARGET_FILE:curcuma>)
        set_tests_properties(qm_mgga_opt_${_case} PROPERTIES TIMEOUT 1800 LABELS "qm_mgga_grad")
    endforeach()
endif()
