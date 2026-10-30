# Native meta-GGAs TPSS and SCAN (WP11 stage 1): kernels vs libxc, V_xc vs FD, energies vs PySCF on the same grid.
# Claude Generated (Oct 2026).
#
# test_qm_mgga (C++) prints JSON; check_qm_mgga.py compares it with reference.json, written by
# scripts/qm_mgga_pyscf_reference.py (PySCF 2.14 / libxc 7.0.0, venv ~/.venvs/qmref; settings in the script).
#   qm_mgga_kernel               TPSS-X/C, SCAN-X/C eps, v_rho, v_sigma, v_tau vs libxc on 1760 (rho, sigma, tau)
#                                points (rho 1e-12..1e2, s 0..30, alpha 0..1e3, SCAN cut-offs, tau < tau_W clamp);
#                                v_* vs finite differences
#   qm_mgga_vxcfd_<mol>_<fn>     V_xc vs finite difference of E_xc in P; 1 vs 4 threads bitwise (unscreened, screened)
#   qm_mgga_energy_<case>        E and E_xc vs PySCF RKS (same grid: Treutler x Lebedev, Becke, no pruning) to 1e-8 Eh;
#                                *_rij/_rijk: -qm.ri j|jk vs PySCF density_fit with the same JFIT/JKFIT file (J only)
#   qm_mgga_gradient_<fn>        a gradient request succeeds (stage 2; values: test_cases/qm_mgga_grad)
#   qm_mgga_invariance_<case>    atom order and translation 1e-9 Eh; rotation (grid fixed in space) at sg1 and sg3
#   qm_mgga_{update,threads,screen}_H2O_<fn>
#   qm_mgga_gpu_H2O_<fn>         -gpu cuda: counted CPU fallback of the meta-GGA XC part, same energy; -gpu_strict exits
# Geometries: the xyz files of qm_1e and qm_hf3c (no new structure files, as qm_pbe).
# Run:  ctest -L qm_mgga --output-on-failure

add_executable(test_qm_mgga test_qm_mgga.cpp)
target_link_libraries(test_qm_mgga PRIVATE curcuma_core)
target_include_directories(test_qm_mgga PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_mgga PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_MGGA_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_mgga.py
        --exe $<TARGET_FILE:test_qm_mgga> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    add_test(NAME qm_mgga_kernel COMMAND ${_MGGA_CHECK} --check kernel)
    set_tests_properties(qm_mgga_kernel PROPERTIES TIMEOUT 120 LABELS "qm_mgga")
    foreach(_case H2O_tpss H2O_scan H2O_tpss_screen H2O_scan_screen benzene_scan_screen)
        add_test(NAME qm_mgga_vxcfd_${_case} COMMAND ${_MGGA_CHECK} --check vxcfd --case ${_case})
        set_tests_properties(qm_mgga_vxcfd_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_mgga")
    endforeach()
    foreach(_fn tpss scan)
        foreach(_case H2O_${_fn}_sg1 CH4_${_fn}_sg1 NH3_${_fn}_sg1 HF_${_fn}_sg1 LiH_${_fn}_sg1 Ne_${_fn}_sg1
                      benzene_${_fn}_sg1 H2O_${_fn}_sg2 H2O_${_fn}_sg3 NH3_${_fn}_sg2 NH3_${_fn}_sg3
                      benzene_${_fn}_sg3 H2O_${_fn}_rij_sg1 benzene_${_fn}_rijk_sg1)
            add_test(NAME qm_mgga_energy_${_case} COMMAND ${_MGGA_CHECK} --check energy --case ${_case})
            set_tests_properties(qm_mgga_energy_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_mgga")
        endforeach()
        add_test(NAME qm_mgga_gradient_${_fn} COMMAND ${_MGGA_CHECK} --check gradient --case H2O_${_fn})
        set_tests_properties(qm_mgga_gradient_${_fn} PROPERTIES TIMEOUT 300 LABELS "qm_mgga")
        foreach(_case H2O_${_fn}_sg1 NH3_${_fn}_sg1 H2O_${_fn}_sg3)
            add_test(NAME qm_mgga_invariance_${_case} COMMAND ${_MGGA_CHECK} --check invariance --case ${_case})
            set_tests_properties(qm_mgga_invariance_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_mgga")
        endforeach()
        foreach(_chk update threads screen)
            add_test(NAME qm_mgga_${_chk}_H2O_${_fn} COMMAND ${_MGGA_CHECK} --check ${_chk} --case H2O_${_fn})
            set_tests_properties(qm_mgga_${_chk}_H2O_${_fn} PROPERTIES TIMEOUT 600 LABELS "qm_mgga")
        endforeach()
        add_test(NAME qm_mgga_gpu_H2O_${_fn} COMMAND ${_MGGA_CHECK} --check gpu --case H2O_${_fn}
                 --curcuma $<TARGET_FILE:curcuma>)
        set_tests_properties(qm_mgga_gpu_H2O_${_fn} PROPERTIES TIMEOUT 600 LABELS "qm_mgga")
    endforeach()
endif()
