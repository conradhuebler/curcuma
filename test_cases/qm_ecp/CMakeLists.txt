# Native QM beyond Ne (Nov 2026). Claude Generated.
#
# Same checker as qm_fg (../qm_fg/check_qm_fg.py with --reference reference.json); reference.json is written by
# scripts/qm_ecp_reference.py (PySCF 2.14 / libxc 7.0.0 and ORCA 6.1.0, settings there). Structures: test_cases/structures.
#   qm_z_ref_<case>     all-electron Na-Kr (def2-SVP/TZVP): energy vs PySCF 1e-8 Eh (and ORCA 1e-8 Eh), gradient vs PySCF 1e-8 Eh/Bohr
#   qm_ecp_ref_<case>   def2-ECP elements (I, Cs, Xe, Hg, Pb, Ce): same tolerances, ORCA per case "orca_tol" (Xe 5e-8, Hg 1e-7)
#   qm_*_fd_<case>      analytic gradient vs a 4-point central difference of curcuma's own energies (1e-8 Eh/Bohr)
#   qm_ecp_gpu_<case>   -gpu cuda run finishes with the CPU energy; 77 = skip without CUDA
#   qm_gpu_ecp_*        ECP matrix and gradient on the device (label qm_gpu_ecp, block at the end)
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
    # ECP matrix and gradient on the device (Nov 2026, label qm_gpu_ecp; 77 = skip without CUDA plugin/device):
    #   qm_gpu_ecp_<case>     check_qm_gpu_ecp.py: device vs CPU of the same binary (E 1e-10 Eh, gradient 1e-10 Eh/Bohr),
    #                         ECP matrix element-wise 1e-13 Eh, equal radial orders, ECP gradient 1e-11 Eh/Bohr; fails on
    #                         a counted qm_ecp fallback
    #   qm_gpu_ecp_fd_<case>  FD gate of the assembled device gradient (check_qm_fg.py --check fd --gpu, 1e-8 Eh/Bohr)
    set(_GPU_ECP ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_ecp.py --dump $<TARGET_FILE:dump_qm_energy>)
    foreach(_spec "hit_tzvp_hf:hfgrad_hit_def2-TZVP"
                  "hgcl2_tzvp_hf:hf_hgcl2_def2-SVP:basis=def2-TZVP"
                  "hit_svp_b3lyp:b3lypgrad_hit_def2-SVP_sg1"
                  "ch3i_svp_hf:hf_ch3i_def2-SVP" "csi_svp_hf:hf_csi_def2-SVP" "xef2_svp_hf:hf_xef2_def2-SVP"
                  "pbh4_svp_hf:hf_pbh4_def2-SVP" "cef4_svp_hf:hf_cef4_def2-SVP")
        string(REPLACE ":" ";" _p "${_spec}")
        list(GET _p 0 _name)
        list(GET _p 1 _case)
        set(_extra "")
        list(LENGTH _p _np)
        if(_np GREATER 2)
            list(GET _p 2 _q)
            set(_extra --qm ${_q})
        endif()
        add_test(NAME qm_gpu_ecp_${_name} COMMAND ${_GPU_ECP} --case ${_case} ${_extra})
        set_tests_properties(qm_gpu_ecp_${_name} PROPERTIES TIMEOUT 600 LABELS "qm_gpu_ecp" SKIP_RETURN_CODE 77)
    endforeach()
    # open shell: HI+ doublet (UHF, the SCF loop stays on the host: qm_uhf_scf is counted by design)
    add_test(NAME qm_gpu_ecp_hit_cation_uhf COMMAND ${_GPU_ECP} --case hfgrad_hit_def2-SVP --charge 1 --spin 1
             --allow qm_uhf_scf)
    set_tests_properties(qm_gpu_ecp_hit_cation_uhf PROPERTIES TIMEOUT 600 LABELS "qm_gpu_ecp" SKIP_RETURN_CODE 77)
    set(_GPU_FD ${_ECP_CHECK} --check fd --gpu --coords 0:x 1:z)
    add_test(NAME qm_gpu_ecp_fd_hit_tzvp_hf COMMAND ${_GPU_FD} --case hfgrad_hit_def2-TZVP)
    add_test(NAME qm_gpu_ecp_fd_hit_tzvp_b3lyp COMMAND ${_GPU_FD} --case b3lypgrad_hit_def2-SVP_sg1 --qm basis=def2-TZVP)
    add_test(NAME qm_gpu_ecp_fd_hgcl2_tzvp_hf COMMAND ${_GPU_FD} --case hf_hgcl2_def2-SVP --qm basis=def2-TZVP --coords 1:z 2:z)
    add_test(NAME qm_gpu_ecp_fd_hgcl2_tzvp_b3lyp COMMAND ${_GPU_FD} --case hf_hgcl2_def2-SVP --method b3lyp
             --qm basis=def2-TZVP grid=sg1 --coords 1:z 2:z)
    add_test(NAME qm_gpu_ecp_fd_hit_cation_uhf COMMAND ${_GPU_FD} --case hfgrad_hit_def2-SVP --charge 1 --spin 1
             --allow qm_uhf_scf)
    set_tests_properties(qm_gpu_ecp_fd_hit_tzvp_hf qm_gpu_ecp_fd_hit_tzvp_b3lyp qm_gpu_ecp_fd_hgcl2_tzvp_hf
                         qm_gpu_ecp_fd_hgcl2_tzvp_b3lyp qm_gpu_ecp_fd_hit_cation_uhf
                         PROPERTIES TIMEOUT 900 LABELS "qm_gpu_ecp" SKIP_RETURN_CODE 77)
endif()
