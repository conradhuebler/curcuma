# Native QM with f and g shells (Nov 2026). Claude Generated.
#
# check_qm_fg.py runs dump_qm_energy (target of qm_basis_embed, built-in basis sets only) and compares with
# reference.json, written by scripts/qm_fg_reference.py (PySCF 2.14 / libxc 7.0.0 and ORCA 6.1.0, settings there):
#   qm_fg_ref_<case>   energy vs PySCF 1e-8 Eh (and ORCA 1e-8 Eh where stored, RI-MP2 2e-8), gradient vs PySCF 1e-8 Eh/Bohr
#   qm_fg_fd_<case>    analytic gradient vs a 4-point central difference of curcuma's own energies (1e-8 Eh/Bohr)
#   qm_fg_gpu_<case>   -gpu cuda run finishes with the CPU energy; 77 = skip without CUDA
#   qm_gpu_fg_*        stage 6b (label qm_gpu_fg): dump_qm_energy --gpu cuda vs the CPU (energy 1e-9 Eh, gradient 1e-10
#                      Eh/Bohr) and FD gate / PySCF reference on the device path; fail on any counted fallback of the f/g
#                      device kernels (qm_eri, qm_direct, qm_grad, qm_grad1e, qm_mp2_grad, qm_ri_grad); 77 = skip
# Geometries: the xyz files of qm_1e (no new structure files).
# Run:  ctest -L qm_fg --output-on-failure
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_FG_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_fg.py --dump $<TARGET_FILE:dump_qm_energy>)
    foreach(_case hf_H2O_def2-TZVP hf_H2O_def2-TZVPP hf_H2O_def2-TZVPD hf_H2O_def2-QZVP hf_NH3_def2-TZVPP
                  hf_CH4_def2-TZVPPD hf_HF_def2-QZVPPD hfcart_H2O_def2-TZVP hfcart_H2O_def2-QZVP
                  lda_H2O_def2-TZVP_sg1 pbe_H2O_def2-TZVP_sg1 b3lyp_H2O_def2-TZVP_sg1 tpss_H2O_def2-TZVP_sg1
                  b3lyp_H2O_def2-QZVP_sg1 tpss_H2O_def2-QZVP_sg1
                  hfgrad_H2O_def2-TZVP hfgrad_NH3_def2-TZVPD b3lypgrad_H2O_def2-TZVP_sg1 tpssgrad_H2O_def2-TZVP_sg1
                  pbegrad_H2O_def2-TZVP_sg1 rijk_H2O_def2-TZVP pberij_H2O_def2-TZVP_sg1 mp2_H2O_def2-TZVP
                  rimp2_H2O_def2-TZVP rimp2_H2O_def2-QZVPP ccsd_H2O_def2-TZVP b2plyp_H2O_def2-TZVP_sg1)
        add_test(NAME qm_fg_ref_${_case} COMMAND ${_FG_CHECK} --check ref --case ${_case})
        set_tests_properties(qm_fg_ref_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_fg")
    endforeach()
    foreach(_case hfgrad_H2O_def2-TZVP b3lypgrad_H2O_def2-TZVP_sg1 tpssgrad_H2O_def2-TZVP_sg1 rijk_H2O_def2-TZVP
                  mp2_H2O_def2-TZVP b2plyp_H2O_def2-TZVP_sg1)
        add_test(NAME qm_fg_fd_${_case} COMMAND ${_FG_CHECK} --check fd --case ${_case})
        set_tests_properties(qm_fg_fd_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_fg")
    endforeach()
    foreach(_case hf_H2O_def2-TZVP b3lyp_H2O_def2-TZVP_sg1)
        add_test(NAME qm_fg_gpu_${_case} COMMAND ${_FG_CHECK} --check gpu --case ${_case} --curcuma $<TARGET_FILE:curcuma>)
        set_tests_properties(qm_fg_gpu_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_fg" SKIP_RETURN_CODE 77)
    endforeach()
    # stage 6b (Nov 2026): f/g shells in the CUDA ERI, J/K and gradient kernels (no counted fallback allowed)
    set(_FG_GPU ${_FG_CHECK} --check gpudev --gtol 1e-10)
    foreach(_case hfgrad_H2O_def2-TZVP hfgrad_H2O_def2-QZVPP hfgrad_NH3_def2-TZVPD b3lypgrad_H2O_def2-QZVP_sg1
                  pbegrad_H2O_def2-TZVP_sg1 rijk_H2O_def2-TZVP pberij_H2O_def2-TZVP_sg1)
        add_test(NAME qm_gpu_fg_${_case} COMMAND ${_FG_GPU} --case ${_case})
    endforeach()
    add_test(NAME qm_gpu_fg_direct_hfgrad_H2O_def2-QZVPP COMMAND ${_FG_GPU} --case hfgrad_H2O_def2-QZVPP --qm scf_direct=on)
    add_test(NAME qm_gpu_fg_cart_H2O_def2-QZVP COMMAND ${_FG_GPU} --case hfcart_H2O_def2-QZVP --gradient)
    add_test(NAME qm_gpu_fg_mp2grad_H2O_def2-TZVP COMMAND ${_FG_GPU} --case mp2_H2O_def2-TZVP --gradient)
    add_test(NAME qm_gpu_fg_b2plypgrad_H2O_def2-TZVP COMMAND ${_FG_GPU} --case b2plyp_H2O_def2-TZVP_sg1 --gradient
             --allow qm_dh_grad_pt2)
    add_test(NAME qm_gpu_fg_ecp_hfgrad_hit_def2-TZVP COMMAND ${_FG_GPU} --case hfgrad_hit_def2-TZVP
             --reference ${CMAKE_CURRENT_SOURCE_DIR}/../qm_ecp/reference.json)
    foreach(_case hfgrad_H2O_def2-QZVPP b3lypgrad_H2O_def2-QZVP_sg1)
        add_test(NAME qm_gpu_fg_ref_${_case} COMMAND ${_FG_CHECK} --check ref --gpu --case ${_case})
    endforeach()
    add_test(NAME qm_gpu_fg_fd_hfgrad_H2O_def2-TZVP COMMAND ${_FG_CHECK} --check fd --gpu --case hfgrad_H2O_def2-TZVP)
    add_test(NAME qm_gpu_fg_fd_b3lypgrad_H2O_def2-QZVP_sg1 COMMAND ${_FG_CHECK} --check fd --gpu --case b3lypgrad_H2O_def2-QZVP_sg1)
    get_property(_fg_tests DIRECTORY PROPERTY TESTS)
    foreach(_t ${_fg_tests})
        if(_t MATCHES "^qm_gpu_fg_")
            set_tests_properties(${_t} PROPERTIES TIMEOUT 900 LABELS "qm_gpu_fg" SKIP_RETURN_CODE 77)
        endif()
    endforeach()
endif()
