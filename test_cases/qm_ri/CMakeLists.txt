# Native HF with density fitting (WP10: RI-J, RI-JK) vs PySCF and ORCA 6.1. Claude Generated (Oct 2026).
#
# check_qm_ri.py runs `curcuma -sp <xyz> -method hf -qm.ri j|jk` and compares the total energy with
# reference.json (PySCF with the same auxiliary basis files, 1e-8 Eh; ORCA RIJONX/RIJK, 1e-7 Eh; the
# generating commands are in the file). The RI error against the exact path is a measurement in
# docs/DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md, not a test criterion.
#   qm_ri_<mol>_<j|jk>      spherical def2-SVP, Cholesky metric
#   qm_ri_H2O_cartesian_*   -qm.cartesian_d true (orbital and auxiliary functions cartesian)
#   qm_ri_*_eigen           -qm.ri_metric eigen (J^-1/2, nothing dropped at the default cutoff)
#   qm_ri_grad_<mol>_<j|jk> analytic RI gradient (-sp -gradient) vs PySCF density_fit().nuc_grad_method() with the
#                           same aux files (reference_grad.json, scripts/ri_gradient_reference.py), 5e-10 Eh/Bohr
#   qm_ri_grad_fd_*         analytic RI gradient vs the finite difference of curcuma's own RI energy (dump_qm_gradient)
#   qm_ri_opt_*_jk          -opt -qm.ri jk reaches the PySCF RI-JK minimum (scripts/ri_opt_reference.py)
#   qm_ri_gradient_error    a gradient that cannot be formed (memory limit) must fail loud (until Oct 2026: every
#                           RI gradient request had to fail; the analytic RI gradient replaced that)
#   qm_ri_hf3c_ignored      hf-3c ignores -qm.ri with a warning
#   qm_gpu_ri_benzene_jk    device contractions (cuBLAS) vs the CPU, exit 77 without device
#   qm_gpu_ri_grad_benzene_* device 3-index RI gradient (cuda/qm_gpu_ri_grad.cu) vs the CPU (1e-10), exit 77 without device
#
# Run:  ctest -L qm_ri --output-on-failure
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_RI_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_ri.py
        --curcuma $<TARGET_FILE:curcuma> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    foreach(_mol H2O NH3 formaldehyde benzene)
        foreach(_ri j jk)
            add_test(NAME qm_ri_${_mol}_${_ri} COMMAND ${_RI_CHECK} --mol ${_mol} --ri ${_ri} "--extra=-qm.threads 4")
            set_tests_properties(qm_ri_${_mol}_${_ri} PROPERTIES TIMEOUT 600 LABELS "qm_ri")
        endforeach()
    endforeach()
    foreach(_ri j jk)
        add_test(NAME qm_ri_H2O_cartesian_${_ri}
            COMMAND ${_RI_CHECK} --mol H2O --ri ${_ri} --set pyscf_cartesian "--extra=-qm.cartesian_d true")
        set_tests_properties(qm_ri_H2O_cartesian_${_ri} PROPERTIES TIMEOUT 600 LABELS "qm_ri")
    endforeach()
    add_test(NAME qm_ri_formaldehyde_jk_eigen
        COMMAND ${_RI_CHECK} --mol formaldehyde --ri jk "--extra=-qm.ri_metric eigen")
    add_test(NAME qm_ri_gradient_error COMMAND ${_RI_CHECK} --mol H2O --ri jk --check gradient_error)
    add_test(NAME qm_ri_hf3c_ignored COMMAND ${_RI_CHECK} --mol H2O --ri jk --check hf3c)
    set_tests_properties(qm_ri_formaldehyde_jk_eigen qm_ri_gradient_error qm_ri_hf3c_ignored
        PROPERTIES TIMEOUT 600 LABELS "qm_ri")
    # RI gradient (Claude Generated, Oct 2026)
    foreach(_mol H2O NH3 formaldehyde benzene)
        foreach(_ri j jk)
            add_test(NAME qm_ri_grad_${_mol}_${_ri} COMMAND ${_RI_CHECK} --mol ${_mol} --ri ${_ri} --check gradient "--extra=-qm.threads 4")
            set_tests_properties(qm_ri_grad_${_mol}_${_ri} PROPERTIES TIMEOUT 900 LABELS "qm_ri")
        endforeach()
    endforeach()
    foreach(_ri j jk)
        add_test(NAME qm_ri_grad_H2O_cartesian_${_ri} COMMAND ${_RI_CHECK} --mol H2O --ri ${_ri} --check gradient --set pyscf_cartesian)
        set_tests_properties(qm_ri_grad_H2O_cartesian_${_ri} PROPERTIES TIMEOUT 600 LABELS "qm_ri")
    endforeach()
    add_test(NAME qm_ri_grad_formaldehyde_jk_eigen
        COMMAND ${_RI_CHECK} --mol formaldehyde --ri jk --check gradient "--extra=-qm.ri_metric eigen")
    add_test(NAME qm_ri_grad_fd_H2O_jk COMMAND ${_RI_CHECK} --mol H2O --ri jk --check fd --dump $<TARGET_FILE:dump_qm_gradient>)
    add_test(NAME qm_ri_grad_fd_NH3_j COMMAND ${_RI_CHECK} --mol NH3 --ri j --check fd --dump $<TARGET_FILE:dump_qm_gradient>)
    foreach(_mol H2O formaldehyde)
        add_test(NAME qm_ri_opt_${_mol}_jk COMMAND ${_RI_CHECK} --mol ${_mol} --ri jk --check opt "--extra=-qm.threads 4")
    endforeach()
    set_tests_properties(qm_ri_grad_formaldehyde_jk_eigen qm_ri_grad_fd_H2O_jk qm_ri_grad_fd_NH3_j qm_ri_opt_H2O_jk
        qm_ri_opt_formaldehyde_jk PROPERTIES TIMEOUT 1500 LABELS "qm_ri")
    add_test(NAME qm_gpu_ri_benzene_jk COMMAND ${_RI_CHECK} --mol benzene --ri jk --check gpu)
    set_tests_properties(qm_gpu_ri_benzene_jk PROPERTIES TIMEOUT 600 LABELS "qm_ri;qm_gpu;gpu" SKIP_RETURN_CODE 77)
    foreach(_ri j jk)
        add_test(NAME qm_gpu_ri_grad_benzene_${_ri} COMMAND ${_RI_CHECK} --mol benzene --ri ${_ri} --check gpu_grad)
        set_tests_properties(qm_gpu_ri_grad_benzene_${_ri} PROPERTIES TIMEOUT 900 LABELS "qm_ri;qm_gpu;gpu" SKIP_RETURN_CODE 77)
    endforeach()
endif()
