# Native QM analytic RHF Hessian (Nov 2026). Claude Generated.
#
# dump_qm_hessian prints the analytic Hessian of `hf` (QMEngine::computeHessian) and, with --fd, the central
# finite difference of curcuma's own analytic gradient. check_qm_hessian.py gates
#   (a) analytic vs FD of the fully assembled gradient (Eh/Bohr^2),
#   (b) symmetry, translational and rotational invariance (H (L r) = L g at every geometry),
#   (c) analytic vs PySCF 2.14 pyscf.hessian.rhf at the same geometry (reference.json, scripts/hf_hessian_reference.py),
#   (d) with --pieces: the derivative integrals against finite differences at fixed density.
# Run:  ctest -L qm_hess --output-on-failure

add_executable(dump_qm_hessian dump_qm_hessian.cpp)
target_link_libraries(dump_qm_hessian PRIVATE curcuma_core)
target_include_directories(dump_qm_hessian PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(dump_qm_hessian PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_QMHESS_CASES H2O_svp NH3_svp CH4_svp HF_svp H2O_distorted_svp formaldehyde_distorted_svp)
    foreach(_case ${_QMHESS_CASES})
        add_test(NAME qm_hess_${_case}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_hessian.py
                --dump $<TARGET_FILE:dump_qm_hessian>
                --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json
                --case ${_case})
        set_tests_properties(qm_hess_${_case} PROPERTIES TIMEOUT 1800 LABELS "qm_hess")
    endforeach()
    # piecewise integral checks (fixed density) on the distorted water
    add_test(NAME qm_hess_pieces_H2O_distorted
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_hessian.py
            --dump $<TARGET_FILE:dump_qm_hessian>
            --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json
            --case H2O_distorted_svp --pieces --no-fd)
    set_tests_properties(qm_hess_pieces_H2O_distorted PROPERTIES TIMEOUT 1800 LABELS "qm_hess")
    # integral-direct SCF: the CPHF J/K go through the engine's buildJK per density instead of the batched GEMMs
    add_test(NAME qm_hess_direct_formaldehyde_distorted
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_hessian.py
            --dump $<TARGET_FILE:dump_qm_hessian>
            --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json
            --case formaldehyde_distorted_svp --qm scf_direct on)
    set_tests_properties(qm_hess_direct_formaldehyde_distorted PROPERTIES TIMEOUT 1800 LABELS "qm_hess")
    # -hessian capability: analytic frequencies equal the numerical (seminumerical) ones
    add_test(NAME qm_hess_cli_frequencies_H2O
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_hessian_cli.py
            --curcuma $<TARGET_FILE:curcuma>
            --xyz ${CMAKE_SOURCE_DIR}/test_cases/qm_1e/H2O.xyz)
    set_tests_properties(qm_hess_cli_frequencies_H2O PROPERTIES TIMEOUT 1800 LABELS "qm_hess")
    # refused case (density fitting): warning and the numerical Hessian, no wrong analytic result
    add_test(NAME qm_hess_cli_fallback_ri_H2O
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_hessian_cli.py
            --curcuma $<TARGET_FILE:curcuma>
            --xyz ${CMAKE_SOURCE_DIR}/test_cases/qm_1e/H2O.xyz "--fallback=-qm.ri jk")
    set_tests_properties(qm_hess_cli_fallback_ri_H2O PROPERTIES TIMEOUT 1800 LABELS "qm_hess")
endif()
