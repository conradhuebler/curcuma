# Native double-hybrid analytic gradient (WP12 stage 2, CPU). Claude Generated (Oct 2026).
#
# test_qm_dh_grad (C++) prints JSON; check_qm_dh_grad.py applies the gates. reference.json:
# scripts/qm_dh_grad_pyscf_reference.py (PySCF 2.14, libxc 7.0.0); reference_orca.json: scripts/qm_dh_grad_orca_compare.py
# (ORCA 6.1.0).
#   qm_dh_grad_kernel              jet derivatives (qm_xc_fxc_kernels.h) vs production v, libxc second derivatives, FD
#   qm_dh_grad_fxc_<case>          F_xc[D] and the directional XC gradient vs FD of V_xc / of the XC gradient
#   qm_dh_grad_fd_<case>_<grid>    analytic vs 4-point FD (h 1e-3 Bohr) of curcuma's energy: moving grid = with grid
#                                  response, frozen grid = without; tol 1e-8 Eh/Bohr
#   qm_dh_grad_pyscf_<case>        vs FD of the PySCF double-hybrid energy (same grid, moving)
#   qm_dh_grad_orca_<case>         vs ORCA 6.1 EnGrad (matched settings, 250x1454, frozen core), stored gradient
#   qm_dh_grad_{translation,rotation,threads,update}_<case>, qm_dh_opt_<case> (-opt reaches a stationary point)
#   qm_dh_grad_direct_<case>       forced integral-direct J/K (eri_screening 1e-10) vs stored tensor: Z-vector converges
# Geometries: xyz files of qm_1e, qm_hf3c and qm_grad (no new structure files).
# Run:  ctest -L qm_dh_grad --output-on-failure
add_executable(test_qm_dh_grad test_qm_dh_grad.cpp)
target_link_libraries(test_qm_dh_grad PRIVATE curcuma_core)
target_include_directories(test_qm_dh_grad PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_dh_grad PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_DG ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_dh_grad.py --exe $<TARGET_FILE:test_qm_dh_grad>
        --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    add_test(NAME qm_dh_grad_kernel COMMAND ${_DG} --check kernel)
    set_tests_properties(qm_dh_grad_kernel PROPERTIES TIMEOUT 120 LABELS "qm_dh_grad")
    foreach(_case H2O_b2plyp H2O_dsd-pbep86 H2O_dsd-pbeb95-d3bj)
        add_test(NAME qm_dh_grad_fxc_${_case} COMMAND ${_DG} --check fxc --case ${_case})
        set_tests_properties(qm_dh_grad_fxc_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_dh_grad")
    endforeach()
    # every table row on H2O (RI, moving and frozen grid), plus canonical / frozen-core / other molecules. The P86 rows
    # run with -qm.p86_lda vwn5 (_vwn5, the curcuma default since Oct 31 2026; the driver default is pz81): the PZ81 local part (libxc LDA_C_PZ) is discontinuous at r_s = 1, so the
    # FD of a PZ81-P86 energy picks up jumps of ~1e-9 Eh per grid point crossing the seam (3e-5 Eh/Bohr at h 1e-3,
    # WP12 doc stage 2); the PZ81 second derivatives are checked by qm_dh_grad_kernel and _fxc_H2O_dsd-pbep86.
    set(_FD_CASES)
    foreach(_fn b2plyp b2plyp-d3bj b2gp-plyp b2gp-plyp-d3bj dsd-blyp dsd-blyp-d3bj dsd-pbep86_vwn5 dsd-pbep86-d3bj_vwn5
                dsd-pbeb95-d3bj revdsd-pbep86_vwn5)
        list(APPEND _FD_CASES H2O_${_fn})
    endforeach()
    list(APPEND _FD_CASES H2O_b2plyp_canonical H2O_b2plyp_fc H2O_dsd-pbeb95-d3bj_canonical_fc NH3_dsd-pbep86_fc_vwn5
         HF_revdsd-pbep86_canonical_fc_vwn5 CH4_dsd-blyp_fc NH3_dsd-pbeb95-d3bj_fc CH4_dsd-pbep86-d3bj_canonical_vwn5)
    foreach(_case ${_FD_CASES})
        foreach(_g moving frozen)
            add_test(NAME qm_dh_grad_fd_${_case}_${_g} COMMAND ${_DG} --check fd --case ${_case} --fd ${_g} --tol 1e-8)
            set_tests_properties(qm_dh_grad_fd_${_case}_${_g} PROPERTIES TIMEOUT 900 LABELS "qm_dh_grad")
        endforeach()
    endforeach()
    # (the PZ81-P86 cases H2O_dsd-pbep86_ri, HF_revdsd-pbep86_ri_fc of reference.json are not gated: PySCF's energy has
    # the same r_s = 1 jumps, measured 3.3e-5 / 2.1e-6 Eh/Bohr, equal to curcuma's own FD)
    foreach(_case H2O_b2plyp_ri H2O_b2plyp_ri_fc H2O_b2plyp_canonical NH3_dsd-blyp_ri H2O_dsd-pbep86_vwn5_ri
                  HF_revdsd-pbep86_vwn5_ri_fc H2O_dsd-pbeb95-d3bj_ri CH4_b2gp-plyp_canonical_fc)
        add_test(NAME qm_dh_grad_pyscf_${_case} COMMAND ${_DG} --check pyscf --case ${_case} --tol 1e-8)
        set_tests_properties(qm_dh_grad_pyscf_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_dh_grad")
    endforeach()
    foreach(_case H2O_b2plyp H2O_dsd-pbep86)
        add_test(NAME qm_dh_grad_orca_${_case} COMMAND ${_DG} --check orca --case ${_case} --tol 3e-6
                 --orca ${CMAKE_CURRENT_SOURCE_DIR}/reference_orca.json)
        set_tests_properties(qm_dh_grad_orca_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh_grad")
    endforeach()
    foreach(_case NH3_b2plyp H2O_dsd-pbeb95-d3bj)
        add_test(NAME qm_dh_grad_translation_${_case} COMMAND ${_DG} --check translation --case ${_case})
        set_tests_properties(qm_dh_grad_translation_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh_grad")
    endforeach()
    foreach(_case NH3_dsd-pbep86 H2O_b2plyp_canonical)
        add_test(NAME qm_dh_grad_rotation_${_case} COMMAND ${_DG} --check rotation --case ${_case})
        set_tests_properties(qm_dh_grad_rotation_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh_grad")
    endforeach()
    foreach(_case NH3_b2plyp NH3_dsd-pbeb95-d3bj_canonical_fc benzene_dsd-pbep86-d3bj)
        add_test(NAME qm_dh_grad_threads_${_case} COMMAND ${_DG} --check threads --case ${_case})
        set_tests_properties(qm_dh_grad_threads_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh_grad")
    endforeach()
    # integral-direct J/K in the Z-vector (power-of-two scaling in dhG4, Oct 31 2026): eri_screening 1e-10 makes the
    # absolute direct screen bite on H2O/NH3 (before the fix: 100 iterations, residual 4e-9 / 1e-8)
    foreach(_case H2O_dsd-pbep86-d3bj_vwn5 NH3_b2plyp)
        add_test(NAME qm_dh_grad_direct_${_case} COMMAND ${_DG} --check direct --case ${_case} --screening 1e-10 --tol 1e-10)
        set_tests_properties(qm_dh_grad_direct_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_dh_grad")
    endforeach()
    foreach(_case H2O_b2plyp H2O_dsd-pbeb95-d3bj_fc)
        add_test(NAME qm_dh_grad_update_${_case} COMMAND ${_DG} --check update --case ${_case})
        set_tests_properties(qm_dh_grad_update_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh_grad")
    endforeach()
    # formaldehyde with dsd-pbep86-d3bj on PZ81-P86 stalls at |g| ~ 3e-5 Eh/Bohr (r_s = 1 jumps of PZ81, the line search
    # fails, exit 1; WP12 doc stage 2): run on VWN5-P86. The opt cases pass -qm.p86_lda / -qm.frozen_core explicitly
    # (driver convention: pz81, no fc); _default runs curcuma's defaults (since Oct 31 2026 vwn5 + frozen core) and
    # must converge with exit 0.
    foreach(_case H2O_distorted_b2plyp NH3_b2plyp formaldehyde_distorted_b2plyp H2O_distorted_dsd-pbep86-d3bj
                  NH3_dsd-pbep86-d3bj formaldehyde_distorted_dsd-pbep86-d3bj_vwn5 H2O_distorted_dsd-pbeb95-d3bj
                  formaldehyde_distorted_dsd-pbep86-d3bj_default)
        add_test(NAME qm_dh_opt_${_case} COMMAND ${_DG} --check opt --case ${_case} --curcuma $<TARGET_FILE:curcuma>)
        set_tests_properties(qm_dh_opt_${_case} PROPERTIES TIMEOUT 1800 LABELS "qm_dh_grad")
    endforeach()
endif()
