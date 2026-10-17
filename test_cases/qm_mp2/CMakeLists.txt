# Native canonical closed-shell MP2 (WP11) vs PySCF. Claude Generated (Oct 2026).
#
# check_mp2.py runs `curcuma -sp <xyz> -method mp2|scs-mp2|sos-mp2 -verbosity 2` and compares
# E(HF), E_OS, E_SS, E_corr, the SCS/SOS values and the total with reference.json, produced by
# scripts/mp2_reference.py (PySCF 2.14 RHF conv_tol 1e-12 + mp.MP2, def2-SVP spherical, frozen=None
# or the chemical frozen core). Tolerance 1e-8 Eh. ORCA 6.1 (`! RHF MP2 def2-SVP ExtremeSCF
# NoFrozenCore|FrozenCore`) was compared by hand, see docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md.
#   qm_mp2_<mol>         full MP2,             qm_mp2_<mol>_fc       frozen core
#   qm_mp2_scs/sos_*     the scaled variants as method names
#   qm_mp2_direct_*      integral-direct SCF: the tensor is rebuilt for MP2
#   qm_mp2_batch1_*      -qm.mp2_batch_mb tiny: one occupied orbital per batch
#   qm_mp2_threads1_*    one thread
#   qm_mp2_opt_<mol>     -opt -method mp2 reaches the PySCF MP2 minimum (replaces qm_mp2_opt_error, which checked
#                        the "no analytic gradient" error before the MP2 gradient existed)
#   qm_mp2_grad_*        analytic MP2 gradient (check_mp2_gradient.py, dump_qm_gradient) against reference_grad.json
#                        from scripts/mp2_gradient_reference.py (PySCF 2.14, ~/.venvs/qmref: RHF conv_tol 1e-13,
#                        conv_tol_grad 1e-11; four-point central difference of the PySCF MP2 energy, h = 1e-3 Bohr,
#                        tol 1e-8 Eh/Bohr; PySCF's analytic mp.MP2(mf).nuc_grad_method() with CPHF tol 1e-13 within
#                        2.5e-7, because it differs from its own FD by up to 1.1e-7); *_fd also against the FD of
#                        curcuma's own energy (1e-7); scs/sos only against the own FD (PySCF has no SCS gradient)
#   qm_gpu_mp2_grad_*    the same with -gpu cuda (2e derivative part on the device; label qm_gpu, exit 77 = skip)
#   qm_gpu_mp2_*         -gpu cuda (label qm_gpu, exit 77 = skip without plugin/device)
# RI-MP2 (methods ri-mp2, scs-ri-mp2, sos-ri-mp2; Claude Generated, Oct 2026) against reference_ri.json, produced by
# scripts/rimp2_reference.py (PySCF 2.14 with the auxiliary basis read from curcuma's own
# def2-universal-{JKFIT,JFIT}.turbomole files; manual Cholesky evaluation, cross-checked with mp.dfmp2.DFMP2 to
# <=3e-15). ORCA 6.1 (`! RHF RI-MP2 def2-SVP ExtremeSCF [Frozen|NoFrozen]Core` + `%basis AuxC "def2/JK" end`) was
# compared by hand, see the WP10 doc. Tolerance 1e-8 Eh.
#   qm_rimp2_<mol>[_fc]  aux def2-universal-JKFIT (default), exact SCF
#   qm_rimp2_jfit_*      -qm.ri_mp2_aux def2-universal-JFIT
#   qm_rimp2_rijk_*      -qm.ri jk: RI-JK SCF, its tensor reused for the correlation part
#   qm_rimp2_eigen_*     -qm.ri_metric eigen;  qm_rimp2_threads1_*  one thread
#   qm_rimp2_opt_error   -opt with ri-mp2 must stop with the "no analytic gradient" error
#   qm_gpu_rimp2_*       -gpu cuda, device pair energies (label qm_gpu, exit 77 = skip)
#
# Run:  ctest -L qm_mp2 --output-on-failure

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_MP2 ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mp2.py
        --curcuma $<TARGET_FILE:curcuma> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    set(_MP2_TESTS)
    foreach(_mol H2O NH3 CH4 HF H2O_distorted formaldehyde_distorted formaldehyde benzene)
        add_test(NAME qm_mp2_${_mol} COMMAND ${_MP2} --mol ${_mol})
        add_test(NAME qm_mp2_${_mol}_fc COMMAND ${_MP2} --mol ${_mol} --fc)
        list(APPEND _MP2_TESTS qm_mp2_${_mol} qm_mp2_${_mol}_fc)
    endforeach()
    add_test(NAME qm_mp2_scs_formaldehyde COMMAND ${_MP2} --mol formaldehyde --method scs-mp2)
    add_test(NAME qm_mp2_sos_NH3_fc COMMAND ${_MP2} --mol NH3 --method sos-mp2 --fc)
    add_test(NAME qm_mp2_direct_formaldehyde_distorted COMMAND ${_MP2} --mol formaldehyde_distorted
        "--extra=-qm.scf_direct on")
    add_test(NAME qm_mp2_batch1_formaldehyde_fc COMMAND ${_MP2} --mol formaldehyde --fc
        "--extra=-qm.mp2_batch_mb 0.001")
    add_test(NAME qm_mp2_threads1_CH4 COMMAND ${_MP2} --mol CH4 "--extra=-qm.threads 1")
    list(APPEND _MP2_TESTS qm_mp2_scs_formaldehyde qm_mp2_sos_NH3_fc qm_mp2_direct_formaldehyde_distorted
        qm_mp2_batch1_formaldehyde_fc qm_mp2_threads1_CH4)
    set_tests_properties(${_MP2_TESTS} PROPERTIES TIMEOUT 900 LABELS "qm_mp2")

    # MP2 gradient (Claude Generated, Oct 2026)
    set(_MP2G ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mp2_gradient.py
        --dump $<TARGET_FILE:dump_qm_gradient> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference_grad.json)
    set(_MP2G_TESTS)
    foreach(_mol H2O NH3 CH4 HF H2O_distorted formaldehyde_distorted formaldehyde benzene)
        add_test(NAME qm_mp2_grad_${_mol} COMMAND ${_MP2G} --mol ${_mol})
        add_test(NAME qm_mp2_grad_${_mol}_fc COMMAND ${_MP2G} --mol ${_mol} --fc)
        list(APPEND _MP2G_TESTS qm_mp2_grad_${_mol} qm_mp2_grad_${_mol}_fc)
    endforeach()
    add_test(NAME qm_mp2_grad_fd_H2O_distorted COMMAND ${_MP2G} --mol H2O_distorted --fd)
    add_test(NAME qm_mp2_grad_fd_NH3_fc COMMAND ${_MP2G} --mol NH3 --fc --fd)
    add_test(NAME qm_mp2_grad_scs_formaldehyde_distorted COMMAND ${_MP2G} --mol formaldehyde_distorted --method scs-mp2 --fd)
    add_test(NAME qm_mp2_grad_sos_NH3_fc COMMAND ${_MP2G} --mol NH3 --fc --method sos-mp2 --fd)
    add_test(NAME qm_mp2_grad_direct_formaldehyde_distorted COMMAND ${_MP2G} --mol formaldehyde_distorted
        "--extra=--qm scf_direct on")
    add_test(NAME qm_mp2_grad_batch1_formaldehyde_fc COMMAND ${_MP2G} --mol formaldehyde --fc
        "--extra=--qm mp2_batch_mb 0.001")
    add_test(NAME qm_mp2_grad_threads1_H2O_distorted COMMAND ${_MP2G} --mol H2O_distorted --threads 1)
    list(APPEND _MP2G_TESTS qm_mp2_grad_fd_H2O_distorted qm_mp2_grad_fd_NH3_fc qm_mp2_grad_scs_formaldehyde_distorted
        qm_mp2_grad_sos_NH3_fc qm_mp2_grad_direct_formaldehyde_distorted qm_mp2_grad_batch1_formaldehyde_fc
        qm_mp2_grad_threads1_H2O_distorted)
    set_tests_properties(${_MP2G_TESTS} PROPERTIES TIMEOUT 1800 LABELS "qm_mp2")

    # -opt -method mp2 must reach the PySCF minimum (until Oct 2026 this was qm_mp2_opt_error: the
    # "no analytic gradient" error, now obsolete because the gradient exists)
    foreach(_mol H2O formaldehyde)
        add_test(NAME qm_mp2_opt_${_mol}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mp2_opt.py --curcuma $<TARGET_FILE:curcuma>
                --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference_grad.json --mol ${_mol})
        set_tests_properties(qm_mp2_opt_${_mol} PROPERTIES TIMEOUT 2500 LABELS "qm_mp2")
    endforeach()

    foreach(_mol formaldehyde_distorted benzene)
        add_test(NAME qm_gpu_mp2_grad_${_mol} COMMAND ${_MP2G} --mol ${_mol} --gpu --curcuma $<TARGET_FILE:curcuma>)
        add_test(NAME qm_gpu_mp2_grad_${_mol}_fc COMMAND ${_MP2G} --mol ${_mol} --fc --gpu --curcuma $<TARGET_FILE:curcuma>)
        set_tests_properties(qm_gpu_mp2_grad_${_mol} qm_gpu_mp2_grad_${_mol}_fc
            PROPERTIES TIMEOUT 1800 LABELS "qm_gpu;gpu;qm_mp2" SKIP_RETURN_CODE 77)
    endforeach()

    # Device path: stored device tensor of the SCF (default), integral-direct SCF on the device
    # (tensor rebuilt on the host and uploaded), with the CPU result as second reference.
    foreach(_mol H2O formaldehyde benzene)
        add_test(NAME qm_gpu_mp2_${_mol} COMMAND ${_MP2} --mol ${_mol} --gpu --cpu-check)
        add_test(NAME qm_gpu_mp2_${_mol}_fc COMMAND ${_MP2} --mol ${_mol} --fc --gpu)
        set_tests_properties(qm_gpu_mp2_${_mol} qm_gpu_mp2_${_mol}_fc
            PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu;qm_mp2" SKIP_RETURN_CODE 77)
    endforeach()
    add_test(NAME qm_gpu_mp2_formaldehyde_direct COMMAND ${_MP2} --mol formaldehyde_distorted --gpu --cpu-check
        "--extra=-qm.gpu_direct on")
    add_test(NAME qm_gpu_mp2_benzene_batch1 COMMAND ${_MP2} --mol benzene --fc --gpu "--extra=-qm.mp2_batch_mb 30")
    set_tests_properties(qm_gpu_mp2_formaldehyde_direct qm_gpu_mp2_benzene_batch1
        PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu;qm_mp2" SKIP_RETURN_CODE 77)

    # ---- RI-MP2 ----
    set(_RIMP2 ${_MP2} --ri-ref ${CMAKE_CURRENT_SOURCE_DIR}/reference_ri.json --method ri-mp2)
    set(_RIMP2_TESTS)
    foreach(_mol H2O NH3 CH4 HF H2O_distorted formaldehyde_distorted formaldehyde benzene)
        add_test(NAME qm_rimp2_${_mol} COMMAND ${_RIMP2} --mol ${_mol})
        add_test(NAME qm_rimp2_${_mol}_fc COMMAND ${_RIMP2} --mol ${_mol} --fc)
        list(APPEND _RIMP2_TESTS qm_rimp2_${_mol} qm_rimp2_${_mol}_fc)
    endforeach()
    add_test(NAME qm_rimp2_jfit_H2O COMMAND ${_RIMP2} --mol H2O --aux def2-universal-JFIT)
    add_test(NAME qm_rimp2_jfit_benzene_fc COMMAND ${_RIMP2} --mol benzene --fc --aux def2-universal-JFIT)
    add_test(NAME qm_rimp2_rijk_formaldehyde COMMAND ${_RIMP2} --mol formaldehyde --scf rijk)
    add_test(NAME qm_rimp2_rijk_benzene_fc COMMAND ${_RIMP2} --mol benzene --fc --scf rijk)
    add_test(NAME qm_rimp2_scs_formaldehyde COMMAND ${_MP2} --ri-ref ${CMAKE_CURRENT_SOURCE_DIR}/reference_ri.json
        --method scs-ri-mp2 --mol formaldehyde)
    add_test(NAME qm_rimp2_sos_NH3_fc COMMAND ${_MP2} --ri-ref ${CMAKE_CURRENT_SOURCE_DIR}/reference_ri.json
        --method sos-ri-mp2 --mol NH3 --fc)
    add_test(NAME qm_rimp2_eigen_formaldehyde_distorted COMMAND ${_RIMP2} --mol formaldehyde_distorted
        "--extra=-qm.ri_metric eigen")
    add_test(NAME qm_rimp2_threads1_CH4 COMMAND ${_RIMP2} --mol CH4 "--extra=-qm.threads 1")
    list(APPEND _RIMP2_TESTS qm_rimp2_jfit_H2O qm_rimp2_jfit_benzene_fc qm_rimp2_rijk_formaldehyde
        qm_rimp2_rijk_benzene_fc qm_rimp2_scs_formaldehyde qm_rimp2_sos_NH3_fc qm_rimp2_eigen_formaldehyde_distorted
        qm_rimp2_threads1_CH4)
    set_tests_properties(${_RIMP2_TESTS} PROPERTIES TIMEOUT 900 LABELS "qm_mp2")
    add_test(NAME qm_rimp2_opt_error
        COMMAND ${CMAKE_COMMAND} -E env sh -c
            "$<TARGET_FILE:curcuma> -opt ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e/H2O.xyz -method ri-mp2 -no_bmt 2>&1 | grep -q 'ri-mp2 has no analytic gradient'")
    set_tests_properties(qm_rimp2_opt_error PROPERTIES TIMEOUT 300 LABELS "qm_mp2")
    add_test(NAME qm_gpu_rimp2_benzene COMMAND ${_RIMP2} --mol benzene --gpu --cpu-check)
    add_test(NAME qm_gpu_rimp2_formaldehyde_fc COMMAND ${_RIMP2} --mol formaldehyde --fc --gpu --cpu-check)
    add_test(NAME qm_gpu_rimp2_benzene_batch1 COMMAND ${_RIMP2} --mol benzene --fc --gpu --cpu-check
        "--extra=-qm.mp2_batch_mb 0.001")
    set_tests_properties(qm_gpu_rimp2_benzene qm_gpu_rimp2_formaldehyde_fc qm_gpu_rimp2_benzene_batch1
        PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu;qm_mp2" SKIP_RETURN_CODE 77)
endif()
