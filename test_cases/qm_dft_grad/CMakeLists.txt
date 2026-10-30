# Native Kohn-Sham analytic gradient (WP8 DFT): lda/pbe/b3lyp. Claude Generated (Oct 2026).
#
# test_qm_dft_grad (C++) prints JSON; check_qm_dft_grad.py compares it with reference.json, written by
# scripts/qm_dft_grad_pyscf_reference.py (PySCF 2.14, scipy BFGS; venv ~/.venvs/qmref). No numpy needed here.
#   qm_dft_grad_pyscf_<case>     gradient without / with grid response vs PySCF grid_response False / True, same grid
#   qm_dft_grad_fd_<case>_<grid> analytic vs 4-point FD of curcuma's energy: moving grid = with response,
#                                frozen grid (fixed in space) = without response
#   qm_dft_grad_translation_*, _rotation_*, _threads_*, _update_*   invariances, 1 vs 4 threads bitwise, UpdateMolecule
#   qm_dft_opt_<mol>_<fn>        -opt reaches the PySCF minimum (same grid, grid response on)
# Geometries: xyz files of qm_1e, qm_hf3c and qm_grad (no new structure files, as qm_lda/qm_pbe/qm_b3lyp).
# Run:  ctest -L qm_dft_grad --output-on-failure

add_executable(test_qm_dft_grad test_qm_dft_grad.cpp)
target_link_libraries(test_qm_dft_grad PRIVATE curcuma_core)
target_include_directories(test_qm_dft_grad PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_dft_grad PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_DG_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_dft_grad.py
        --exe $<TARGET_FILE:test_qm_dft_grad> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    set(_DG_PYSCF_CASES)
    foreach(_m H2O NH3 CH4 HF benzene)
        list(APPEND _DG_PYSCF_CASES ${_m}_lda_vwn5_sg1 ${_m}_pbe_vwn5_sg1 ${_m}_b3lyp_vwn5_sg1 ${_m}_b3lyp_vwn_rpa_sg1)
    endforeach()
    list(APPEND _DG_PYSCF_CASES H2O_lda_vwn_rpa_sg1 H2O_b3lyp_vwn5_sg3 H2O_b3lyp_vwn5_rij_sg1 H2O_b3lyp_vwn5_rijk_sg1
         NH3_pbe_vwn5_rijk_sg1 NH3_lda_vwn5_rij_sg1 benzene_b3lyp_vwn_rpa_rijk_sg1)
    foreach(_case ${_DG_PYSCF_CASES})
        add_test(NAME qm_dft_grad_pyscf_${_case} COMMAND ${_DG_CHECK} --check pyscf --case ${_case})
        set_tests_properties(qm_dft_grad_pyscf_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dft_grad")
    endforeach()
    foreach(_case H2O_lda H2O_pbe H2O_b3lyp NH3_b3lyp_vwn_rpa H2O_b3lyp_rij H2O_pbe_rijk)
        foreach(_fd moving frozen)
            add_test(NAME qm_dft_grad_fd_${_case}_${_fd} COMMAND ${_DG_CHECK} --check fd --case ${_case} --fd ${_fd})
            set_tests_properties(qm_dft_grad_fd_${_case}_${_fd} PROPERTIES TIMEOUT 900 LABELS "qm_dft_grad")
        endforeach()
    endforeach()
    foreach(_chk translation rotation)
        foreach(_case NH3_b3lyp H2O_pbe)
            add_test(NAME qm_dft_grad_${_chk}_${_case} COMMAND ${_DG_CHECK} --check ${_chk} --case ${_case})
            set_tests_properties(qm_dft_grad_${_chk}_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dft_grad")
        endforeach()
    endforeach()
    foreach(_case NH3_lda NH3_pbe NH3_b3lyp benzene_b3lyp)
        add_test(NAME qm_dft_grad_threads_${_case} COMMAND ${_DG_CHECK} --check threads --case ${_case})
        set_tests_properties(qm_dft_grad_threads_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dft_grad")
    endforeach()
    foreach(_case H2O_b3lyp H2O_pbe)
        add_test(NAME qm_dft_grad_update_${_case} COMMAND ${_DG_CHECK} --check update --case ${_case})
        set_tests_properties(qm_dft_grad_update_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dft_grad")
    endforeach()
    foreach(_mol H2O_distorted formaldehyde_distorted)
        foreach(_fn lda pbe b3lyp)
            add_test(NAME qm_dft_opt_${_mol}_${_fn}
                COMMAND ${_DG_CHECK} --check opt --case ${_mol}_${_fn} --curcuma $<TARGET_FILE:curcuma>)
            set_tests_properties(qm_dft_opt_${_mol}_${_fn} PROPERTIES TIMEOUT 1800 LABELS "qm_dft_grad")
        endforeach()
    endforeach()
endif()
