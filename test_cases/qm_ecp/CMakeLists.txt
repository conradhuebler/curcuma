# Native QM beyond Ne (Nov 2026). Claude Generated.
#
# Same checker as qm_fg (../qm_fg/check_qm_fg.py with --reference reference.json); reference.json is written by
# scripts/qm_ecp_reference.py (PySCF 2.14 / libxc 7.0.0 and ORCA 6.1.0, settings there). Structures: test_cases/structures.
#   qm_z_ref_<case>     all-electron Na-Kr (def2-SVP/TZVP): energy vs PySCF 1e-8 Eh (and ORCA 1e-8 Eh), gradient vs PySCF 1e-8 Eh/Bohr
#   qm_ecp_ref_<case>   def2-ECP elements (I, Cs, Xe, Hg, Pb, Ce): same tolerances, ORCA per case "orca_tol" (Xe 5e-8, Hg 1e-7)
#   qm_*_fd_<case>      analytic gradient vs a 4-point central difference of curcuma's own energies (1e-8 Eh/Bohr)
#   qm_ecp_gpu_<case>   -gpu cuda run finishes with the CPU energy; 77 = skip without CUDA
# Run:  ctest -L "qm_z|qm_ecp" --output-on-failure
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_ECP_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../qm_fg/check_qm_fg.py
        --reference ${CMAKE_CURRENT_SOURCE_DIR}/reference.json --dump $<TARGET_FILE:dump_qm_energy>)
    foreach(_case hf_nacl_def2-SVP hf_zncl2_def2-SVP hf_hbr_def2-TZVP b3lyp_so2_def2-SVP_sg1 pbe_kcl_def2-TZVP_sg1
                  hfgrad_h2s_def2-SVP b3lypgrad_so2_def2-SVP_sg1 mp2fc_nacl_def2-SVP)
        add_test(NAME qm_z_ref_${_case} COMMAND ${_ECP_CHECK} --check ref --case ${_case})
        set_tests_properties(qm_z_ref_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_z")
    endforeach()
    add_test(NAME qm_z_fd_hfgrad_h2s_def2-SVP COMMAND ${_ECP_CHECK} --check fd --case hfgrad_h2s_def2-SVP --coords 0:z 1:x 1:z)
    set_tests_properties(qm_z_fd_hfgrad_h2s_def2-SVP PROPERTIES TIMEOUT 600 LABELS "qm_z")
    foreach(_case hf_hi_def2-SVP hf_hi_def2-TZVP hf_ch3i_def2-SVP hf_csi_def2-SVP hf_xef2_def2-SVP hf_hgcl2_def2-SVP
                  hf_pbh4_def2-SVP hf_cef4_def2-SVP b3lyp_hi_def2-SVP_sg1 hfgrad_hit_def2-SVP hfgrad_hit_def2-TZVP
                  b3lypgrad_hit_def2-SVP_sg1 rijkgrad_hit_def2-SVP mp2fc_hi_def2-SVP rimp2fc_hi_def2-SVP)
        add_test(NAME qm_ecp_ref_${_case} COMMAND ${_ECP_CHECK} --check ref --case ${_case})
        set_tests_properties(qm_ecp_ref_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_ecp")
    endforeach()
    foreach(_case hfgrad_hit_def2-SVP b3lypgrad_hit_def2-SVP_sg1 mp2grad_hit_def2-SVP b2plypgrad_hit_def2-SVP_sg1
                  rijkgrad_hit_def2-SVP)
        add_test(NAME qm_ecp_fd_${_case} COMMAND ${_ECP_CHECK} --check fd --case ${_case} --coords 0:x 0:z 1:z)
        set_tests_properties(qm_ecp_fd_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_ecp")
    endforeach()
    foreach(_case hf_ch3i_def2-SVP b3lyp_hi_def2-SVP_sg1)
        add_test(NAME qm_ecp_gpu_${_case} COMMAND ${_ECP_CHECK} --check gpu --case ${_case} --curcuma $<TARGET_FILE:curcuma>)
        set_tests_properties(qm_ecp_gpu_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_ecp" SKIP_RETURN_CODE 77)
    endforeach()
endif()
