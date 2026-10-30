# Native B3LYP (WP7): B88/LYP kernels vs libxc, V_xc vs FD, hybrid SCF energies vs PySCF RKS on the identical grid.
# Claude Generated (Oct 2026).
#
# test_qm_b3lyp (C++) prints JSON; check_qm_b3lyp.py compares it with reference.json, written by
# scripts/qm_b3lyp_pyscf_reference.py (PySCF 2.14 / libxc 7.0.0, venv ~/.venvs/qmref; settings in the script).
#   qm_b3lyp_kernel            B88, LYP eps, v_rho, v_sigma vs libxc on 472 (rho, sigma) pairs (rho 1e-12..1e2,
#                              sigma 0..8e12); v_rho, v_sigma vs finite differences; no NaN/inf down to rho 1e-20
#   qm_b3lyp_vxcfd_<mol>       V_xc vs finite difference of E_xc in the density matrix; 1 vs 4 threads bitwise
#   qm_b3lyp_energy_<case>     E, E_xc, -(a_x/4) Tr(PK) vs PySCF RKS (same grid) to 1e-8 Eh, both VWN variants
#                              (vwn5 = libxc B3LYP5, vwn_rpa = libxc B3LYP); *_rij: J fitted + exact K (PySCF
#                              only_dfj), *_rijk: J and K fitted, same JFIT/JKFIT files
#   qm_b3lyp_gradient          a gradient request for b3lyp succeeds (WP8 DFT, values in qm_dft_grad)
#   qm_b3lyp_invariance_<mol>  atom order and translation 1e-9 Eh; rotation (grid fixed in space) 1e-4 Eh at sg1
#   qm_b3lyp_update_H2O, qm_b3lyp_threads_H2O, qm_b3lyp_direct_H2O, qm_b3lyp_rij_H2O
# Geometries: the xyz files of qm_1e and qm_hf3c (no new structure files, as qm_lda/qm_pbe).
# Run:  ctest -L qm_b3lyp --output-on-failure

add_executable(test_qm_b3lyp test_qm_b3lyp.cpp)
target_link_libraries(test_qm_b3lyp PRIVATE curcuma_core)
target_include_directories(test_qm_b3lyp PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_b3lyp PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_B3_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_b3lyp.py
        --exe $<TARGET_FILE:test_qm_b3lyp> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    add_test(NAME qm_b3lyp_kernel COMMAND ${_B3_CHECK} --check kernel)
    set_tests_properties(qm_b3lyp_kernel PROPERTIES TIMEOUT 60 LABELS "qm_b3lyp")
    foreach(_mol H2O benzene)
        add_test(NAME qm_b3lyp_vxcfd_${_mol} COMMAND ${_B3_CHECK} --check vxcfd --case ${_mol})
        set_tests_properties(qm_b3lyp_vxcfd_${_mol} PROPERTIES TIMEOUT 600 LABELS "qm_b3lyp")
    endforeach()
    foreach(_case H2O_vwn5_sg1 CH4_vwn5_sg1 NH3_vwn5_sg1 HF_vwn5_sg1 LiH_vwn5_sg1 Ne_vwn5_sg1 benzene_vwn5_sg1
                  H2O_vwn_rpa_sg1 CH4_vwn_rpa_sg1 NH3_vwn_rpa_sg1 HF_vwn_rpa_sg1 LiH_vwn_rpa_sg1 Ne_vwn_rpa_sg1
                  benzene_vwn_rpa_sg1 H2O_vwn5_sg3 NH3_vwn_rpa_sg2 H2O_vwn5_rij_sg1 benzene_vwn5_rijk_sg1
                  H2O_vwn_rpa_rijk_sg1)
        add_test(NAME qm_b3lyp_energy_${_case} COMMAND ${_B3_CHECK} --check energy --case ${_case})
        set_tests_properties(qm_b3lyp_energy_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_b3lyp")
    endforeach()
    add_test(NAME qm_b3lyp_gradient COMMAND ${_B3_CHECK} --check gradient)
    set_tests_properties(qm_b3lyp_gradient PROPERTIES TIMEOUT 300 LABELS "qm_b3lyp")
    foreach(_mol H2O NH3)
        add_test(NAME qm_b3lyp_invariance_${_mol} COMMAND ${_B3_CHECK} --check invariance --case ${_mol})
        set_tests_properties(qm_b3lyp_invariance_${_mol} PROPERTIES TIMEOUT 600 LABELS "qm_b3lyp")
    endforeach()
    foreach(_chk update threads direct rij)
        add_test(NAME qm_b3lyp_${_chk}_H2O COMMAND ${_B3_CHECK} --check ${_chk} --case H2O)
        set_tests_properties(qm_b3lyp_${_chk}_H2O PROPERTIES TIMEOUT 600 LABELS "qm_b3lyp")
    endforeach()
endif()
