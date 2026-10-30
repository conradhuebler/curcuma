# Native PBE (WP6): XC kernels vs libxc, V_xc vs finite differences, SCF energies vs PySCF RKS on the identical grid.
# Claude Generated (Oct 2026).
#
# test_qm_pbe (C++) prints JSON; check_qm_pbe.py compares it with reference.json, written by
# scripts/qm_pbe_pyscf_reference.py (PySCF 2.14 / libxc 7.0.0, venv ~/.venvs/qmref; settings in the script).
#   qm_pbe_kernel             PBE-X, PBE-C, PW92 eps, v_rho, v_sigma vs libxc on 464 (rho, sigma) pairs
#                             (rho 1e-12..1e2, sigma 0..1e10 and s 0.01..1e3); v_rho, v_sigma vs finite differences
#   qm_pbe_vxcfd_<mol>        V_xc vs finite difference of E_xc in the density matrix; 1 vs 4 threads bitwise
#   qm_pbe_energy_<case>      E and E_xc vs PySCF RKS (same grid: Treutler x Lebedev, Becke, no pruning) to 1e-8 Eh;
#                             *_rij/_rijk: -qm.ri j|jk vs PySCF density_fit with the same JFIT/JKFIT file (J only)
#   qm_pbe_gradient           a gradient request for pbe succeeds (WP8 DFT, values in qm_dft_grad)
#   qm_pbe_invariance_<mol>   atom order and translation 1e-9 Eh; rotation (grid fixed in space) 1e-4 Eh at sg1
#   qm_pbe_update_H2O, qm_pbe_threads_H2O, qm_pbe_direct_H2O, qm_pbe_rij_H2O
# Geometries: the xyz files of qm_1e and qm_hf3c (no new structure files, as qm_lda).
# Run:  ctest -L qm_pbe --output-on-failure

add_executable(test_qm_pbe test_qm_pbe.cpp)
target_link_libraries(test_qm_pbe PRIVATE curcuma_core)
target_include_directories(test_qm_pbe PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_pbe PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_PBE_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_pbe.py
        --exe $<TARGET_FILE:test_qm_pbe> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    add_test(NAME qm_pbe_kernel COMMAND ${_PBE_CHECK} --check kernel)
    set_tests_properties(qm_pbe_kernel PROPERTIES TIMEOUT 60 LABELS "qm_pbe")
    foreach(_mol H2O benzene)
        add_test(NAME qm_pbe_vxcfd_${_mol} COMMAND ${_PBE_CHECK} --check vxcfd --case ${_mol})
        set_tests_properties(qm_pbe_vxcfd_${_mol} PROPERTIES TIMEOUT 600 LABELS "qm_pbe")
    endforeach()
    foreach(_case H2O_sg1 CH4_sg1 NH3_sg1 HF_sg1 LiH_sg1 Ne_sg1 benzene_sg1 H2O_sg3 NH3_sg2
                  H2O_rij_sg1 benzene_rijk_sg1)
        add_test(NAME qm_pbe_energy_${_case} COMMAND ${_PBE_CHECK} --check energy --case ${_case})
        set_tests_properties(qm_pbe_energy_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_pbe")
    endforeach()
    add_test(NAME qm_pbe_gradient COMMAND ${_PBE_CHECK} --check gradient)
    set_tests_properties(qm_pbe_gradient PROPERTIES TIMEOUT 300 LABELS "qm_pbe")
    foreach(_mol H2O NH3)
        add_test(NAME qm_pbe_invariance_${_mol} COMMAND ${_PBE_CHECK} --check invariance --case ${_mol})
        set_tests_properties(qm_pbe_invariance_${_mol} PROPERTIES TIMEOUT 600 LABELS "qm_pbe")
    endforeach()
    foreach(_chk update threads direct rij)
        add_test(NAME qm_pbe_${_chk}_H2O COMMAND ${_PBE_CHECK} --check ${_chk} --case H2O)
        set_tests_properties(qm_pbe_${_chk}_H2O PROPERTIES TIMEOUT 600 LABELS "qm_pbe")
    endforeach()
endif()
