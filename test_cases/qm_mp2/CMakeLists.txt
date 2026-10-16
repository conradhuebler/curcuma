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
#   qm_mp2_opt_error     -opt with mp2 must stop with the "no analytic gradient" error
#   qm_gpu_mp2_*         -gpu cuda (label qm_gpu, exit 77 = skip without plugin/device)
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

    # -opt needs a gradient: must fail with the clear error, not optimise on a wrong gradient
    add_test(NAME qm_mp2_opt_error
        COMMAND ${CMAKE_COMMAND} -E env sh -c
            "$<TARGET_FILE:curcuma> -opt ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e/H2O.xyz -method mp2 -no_bmt 2>&1 | grep -q 'mp2 has no analytic gradient'")
    set_tests_properties(qm_mp2_opt_error PROPERTIES TIMEOUT 300 LABELS "qm_mp2")

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
endif()
