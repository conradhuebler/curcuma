# Native LDA (WP5): XC kernels vs libxc, SCF energies vs PySCF RKS on the identical grid, invariances, caches.
# Claude Generated (Oct 2026).
#
# test_qm_lda (C++) prints JSON; check_qm_lda.py compares it with reference.json, written by
# scripts/qm_lda_pyscf_reference.py (PySCF 2.14 / libxc 7.0.0, venv ~/.venvs/qmref; settings in the script).
#   qm_lda_kernel             Slater, VWN5, VWN-RPA eps and v vs libxc (rho 1e-12..1e2) to 1e-12 relative, v vs FD
#   qm_lda_energy_<case>      E and E_xc vs PySCF RKS (same grid: Treutler x Lebedev, Becke, no pruning) to 1e-8 Eh;
#                             *_rij/_rijk: -qm.ri j|jk vs PySCF density_fit with the same JFIT/JKFIT file (J only)
#   qm_lda_gradient           a gradient request for lda succeeds (WP8 DFT, values in qm_dft_grad)
#   qm_lda_invariance_<mol>   atom order and translation 1e-9 Eh; rotation (grid fixed in space) 1e-4 Eh at sg1
#   qm_lda_update_H2O         UpdateMolecule invalidates the grid and the V_xc cache
#   qm_lda_threads_H2O, qm_lda_direct_H2O, qm_lda_rij_H2O
# Geometries: the xyz files of qm_1e and qm_hf3c (no new structure files, as qm_grid).
# Run:  ctest -L qm_lda --output-on-failure

add_executable(test_qm_lda test_qm_lda.cpp)
target_link_libraries(test_qm_lda PRIVATE curcuma_core)
target_include_directories(test_qm_lda PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_lda PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_LDA_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_lda.py
        --exe $<TARGET_FILE:test_qm_lda> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    add_test(NAME qm_lda_kernel COMMAND ${_LDA_CHECK} --check kernel)
    set_tests_properties(qm_lda_kernel PROPERTIES TIMEOUT 60 LABELS "qm_lda")
    foreach(_case H2O_sg1 CH4_sg1 NH3_sg1 HF_sg1 LiH_sg1 Ne_sg1 benzene_sg1 H2O_sg3 H2O_rpa_sg1 NH3_rpa_sg2
                  H2O_rij_sg1 benzene_rijk_sg1)
        add_test(NAME qm_lda_energy_${_case} COMMAND ${_LDA_CHECK} --check energy --case ${_case})
        set_tests_properties(qm_lda_energy_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_lda")
    endforeach()
    add_test(NAME qm_lda_gradient COMMAND ${_LDA_CHECK} --check gradient)
    set_tests_properties(qm_lda_gradient PROPERTIES TIMEOUT 300 LABELS "qm_lda")
    foreach(_mol H2O NH3)
        add_test(NAME qm_lda_invariance_${_mol} COMMAND ${_LDA_CHECK} --check invariance --case ${_mol})
        set_tests_properties(qm_lda_invariance_${_mol} PROPERTIES TIMEOUT 600 LABELS "qm_lda")
    endforeach()
    foreach(_chk update threads direct rij)
        add_test(NAME qm_lda_${_chk}_H2O COMMAND ${_LDA_CHECK} --check ${_chk} --case H2O)
        set_tests_properties(qm_lda_${_chk}_H2O PROPERTIES TIMEOUT 600 LABELS "qm_lda")
    endforeach()
endif()
