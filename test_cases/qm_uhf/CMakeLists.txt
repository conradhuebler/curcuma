# Open-shell native QM (stage 5b, Nov 2026): UHF and UKS. Claude Generated.
#
# check_qm_uhf.py runs dump_qm_energy (target of qm_basis_embed, built-in basis sets only, --charge/--spin) and
# compares with reference.json, written by scripts/qm_uhf_reference.py (PySCF 2.14 / libxc 7.0.0 UHF/UKS, ORCA 6.1.0 UHF
# VeryTightSCF; settings in the script). Structures: test_cases/structures (atoms/, small/).
#   qm_uhf_ref_<case>     energy vs PySCF 1e-8 Eh and ORCA 1e-8 Eh, <S^2> vs PySCF 1e-6 (ORCA 1e-5), gradient 1e-8 Eh/Bohr
#   qm_uhf_fd_<case>      analytic UHF gradient vs a 4-point central difference of curcuma's own energies (1e-8 Eh/Bohr)
#   qm_uhf_rhf_hf         -qm.reference uhf on a singlet = RHF (energy 1e-9 Eh, gradient 1e-8 Eh/Bohr, <S^2> = 0)
#   qm_uhf_direct_<case>  integral-direct UHF = stored tensor (1e-9 Eh, gradient 1e-8 Eh/Bohr)
#   qm_uhf_reject_<what>  closed-shell-only methods (mp2, ri-mp2, ccsd, lmp2, b2plyp), meta-GGAs (tpss, scan) and a
#                         contradicting spin state (odd electrons with -spin 0, -qm.reference rhf with -spin 1) are errors
#   qm_gpu_uhf_*          stage 6a, label qm_gpu_uhf: -gpu cuda (dump_qm_energy --gpu) vs the CPU (energy 1e-9 Eh, <S^2>
#                         1e-8, gradient 1e-8 Eh/Bohr), vs PySCF, and the FD gate of the device path; only the counted
#                         qm_uhf_scf fallback (host SCF loop) is allowed; 77 = skip without CUDA plugin/device
#   qm_uhf_opt_hf         curcuma -opt of a distorted CH3 radical
#   qm_uks_kernel         spin-polarised kernels vs libxc (spin = 1), own FD, closed-shell limit
#   qm_uks_ref_<case>     UKS energy / <S^2> / gradient vs PySCF on the same grid (sg1)
#   qm_uks_fd_<case>      UKS gradient (grid response) vs FD; qm_uks_rhf_<f>: -qm.reference uhf on a singlet = RKS
# Run:  ctest -L "qm_uhf|qm_uks" --output-on-failure;  GPU: ctest -L qm_gpu_uhf
add_executable(test_qm_uks_kernel test_qm_uks_kernel.cpp)
target_link_libraries(test_qm_uks_kernel PRIVATE curcuma_core)
target_include_directories(test_qm_uks_kernel PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external
                           )
target_compile_features(test_qm_uks_kernel PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_UHF_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_uhf.py --dump $<TARGET_FILE:dump_qm_energy>)
    foreach(_case uhf_oh_svp uhf_oh_tzvp uhf_nh2_svp uhf_ch3_svp uhf_o2_svp uhf_o2_tzvp uhf_no_svp uhf_cn_svp uhf_ch2_svp
                  uhf_h2op_svp uhf_li_svp uhf_na_svp uhf_k_svp uhf_cl_tzvp uhf_i_svp uhf_cu_svp uhf_oh_svp_rij
                  uhf_oh_svp_rijk uhfgrad_ch3p_svp uhfgrad_ch2d_svp uhfgrad_ch2d_svp_rijk uhfgrad_ch2d_svp_rij
                  uhfgrad_h2op_svp)
        add_test(NAME qm_uhf_ref_${_case} COMMAND ${_UHF_CHECK} --check ref --case ${_case})
        set_tests_properties(qm_uhf_ref_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_uhf")
    endforeach()
    # CH2 is planar (z = 0): coordinates in the plane; CH3 (distorted) has no symmetry
    foreach(_case uhfgrad_ch2d_svp uhfgrad_ch2d_svp_rijk)
        add_test(NAME qm_uhf_fd_${_case} COMMAND ${_UHF_CHECK} --check fd --case ${_case} --coords 0:x 1:y 2:x)
        set_tests_properties(qm_uhf_fd_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_uhf")
    endforeach()
    add_test(NAME qm_uhf_fd_hf3cgrad_ch3p COMMAND ${_UHF_CHECK} --check fd --case hf3cgrad_ch3p)
    set_tests_properties(qm_uhf_fd_hf3cgrad_ch3p PROPERTIES TIMEOUT 900 LABELS "qm_uhf")
    add_test(NAME qm_uhf_rhf_hf COMMAND ${_UHF_CHECK} --check rhf --case hf)
    set_tests_properties(qm_uhf_rhf_hf PROPERTIES TIMEOUT 600 LABELS "qm_uhf")
    add_test(NAME qm_uhf_direct_uhfgrad_ch2d_svp COMMAND ${_UHF_CHECK} --check direct --case uhfgrad_ch2d_svp)
    set_tests_properties(qm_uhf_direct_uhfgrad_ch2d_svp PROPERTIES TIMEOUT 600 LABELS "qm_uhf")
    foreach(_what mp2 ri-mp2 ccsd lmp2 b2plyp tpss scan parity rohf)
        add_test(NAME qm_uhf_reject_${_what} COMMAND ${_UHF_CHECK} --check reject --case ${_what})
        set_tests_properties(qm_uhf_reject_${_what} PROPERTIES TIMEOUT 300 LABELS "qm_uhf")
    endforeach()
    # stage 6a (Nov 2026): open shell on the CUDA backend (device J/K, V_xc, gradient parts; host SCF loop)
    set(_UHF_GPU ${_UHF_CHECK} --plugin-dir $<TARGET_FILE_DIR:curcuma>)
    foreach(_case uhfgrad_ch2d_svp uhfgrad_ch3p_svp uhfgrad_ch2d_svp_rijk uhf_o2_svp uksgrad_pbe_ch2d_svp
                  uksgrad_b3lyp_ch3p_svp uksgrad_lda_h2op_svp uks_b3lyp_vwn5_nh2_svp_rijk hf3cgrad_ch3p)
        add_test(NAME qm_gpu_uhf_${_case} COMMAND ${_UHF_GPU} --check gpu --case ${_case})
        set_tests_properties(qm_gpu_uhf_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_gpu_uhf" SKIP_RETURN_CODE 77)
    endforeach()
    add_test(NAME qm_gpu_uhf_direct_uhfgrad_ch2d_svp COMMAND ${_UHF_GPU} --check gpu --case uhfgrad_ch2d_svp
             --extra-qm scf_direct=on)
    add_test(NAME qm_gpu_uhf_direct_uksgrad_b3lyp_ch3p_svp COMMAND ${_UHF_GPU} --check gpu --case uksgrad_b3lyp_ch3p_svp
             --extra-qm scf_direct=on)
    set_tests_properties(qm_gpu_uhf_direct_uhfgrad_ch2d_svp qm_gpu_uhf_direct_uksgrad_b3lyp_ch3p_svp PROPERTIES
                         TIMEOUT 600 LABELS "qm_gpu_uhf" SKIP_RETURN_CODE 77)
    foreach(_case uhfgrad_ch2d_svp uhf_k_svp uksgrad_b3lyp_ch3p_svp uksgrad_pbe_ch2d_svp)
        add_test(NAME qm_gpu_uhf_ref_${_case} COMMAND ${_UHF_GPU} --check ref --gpu --case ${_case})
        set_tests_properties(qm_gpu_uhf_ref_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_gpu_uhf" SKIP_RETURN_CODE 77)
    endforeach()
    add_test(NAME qm_gpu_uhf_fd_uhfgrad_ch2d_svp COMMAND ${_UHF_GPU} --check fd --gpu --case uhfgrad_ch2d_svp --coords 0:x 1:y 2:x)
    add_test(NAME qm_gpu_uhf_fd_uksgrad_pbe_ch2d_svp COMMAND ${_UHF_GPU} --check fd --gpu --case uksgrad_pbe_ch2d_svp
             --coords 0:x 1:y 2:x)
    add_test(NAME qm_gpu_uhf_fd_uksgrad_b3lyp_ch3p_svp COMMAND ${_UHF_GPU} --check fd --gpu --case uksgrad_b3lyp_ch3p_svp)
    set_tests_properties(qm_gpu_uhf_fd_uhfgrad_ch2d_svp qm_gpu_uhf_fd_uksgrad_pbe_ch2d_svp qm_gpu_uhf_fd_uksgrad_b3lyp_ch3p_svp
                         PROPERTIES TIMEOUT 900 LABELS "qm_gpu_uhf" SKIP_RETURN_CODE 77)
    add_test(NAME qm_uhf_opt_hf COMMAND ${_UHF_CHECK} --check opt --case hf --curcuma $<TARGET_FILE:curcuma>)
    set_tests_properties(qm_uhf_opt_hf PROPERTIES TIMEOUT 900 LABELS "qm_uhf")

    add_test(NAME qm_uks_kernel COMMAND ${_UHF_CHECK} --check kernel --exe $<TARGET_FILE:test_qm_uks_kernel>)
    set_tests_properties(qm_uks_kernel PROPERTIES TIMEOUT 60 LABELS "qm_uks")
    foreach(_case uks_lda_vwn5_nh2_svp uks_pbe_nh2_svp uks_b3lyp_vwn5_nh2_svp uks_lda_vwn5_o2_svp uks_pbe_o2_svp
                  uks_b3lyp_vwn5_o2_svp uks_lda_vwn5_ch3_svp uks_pbe_ch3_svp uks_b3lyp_vwn5_ch3_svp uks_lda_vwnrpa_nh2_svp
                  uks_b3lyp_vwnrpa_o2_svp uks_b3lyp_vwn5_na_tzvp uks_pbe_li_svp uks_b3lyp_vwn5_nh2_svp_rijk
                  uksgrad_pbe_ch2d_svp uksgrad_b3lyp_ch3p_svp uksgrad_lda_h2op_svp)
        add_test(NAME qm_uks_ref_${_case} COMMAND ${_UHF_CHECK} --check ref --case ${_case})
        set_tests_properties(qm_uks_ref_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_uks")
    endforeach()
    add_test(NAME qm_uks_fd_uksgrad_pbe_ch2d_svp COMMAND ${_UHF_CHECK} --check fd --case uksgrad_pbe_ch2d_svp --coords 0:x 1:y 2:x)
    set_tests_properties(qm_uks_fd_uksgrad_pbe_ch2d_svp PROPERTIES TIMEOUT 900 LABELS "qm_uks")
    add_test(NAME qm_uks_fd_uksgrad_b3lyp_ch3p_svp COMMAND ${_UHF_CHECK} --check fd --case uksgrad_b3lyp_ch3p_svp)
    set_tests_properties(qm_uks_fd_uksgrad_b3lyp_ch3p_svp PROPERTIES TIMEOUT 900 LABELS "qm_uks")
    foreach(_f pbe b3lyp)
        add_test(NAME qm_uks_rhf_${_f} COMMAND ${_UHF_CHECK} --check rhf --case ${_f})
        set_tests_properties(qm_uks_rhf_${_f} PROPERTIES TIMEOUT 600 LABELS "qm_uks")
    endforeach()
endif()
