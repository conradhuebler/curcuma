# Native QM analytic gradient (WP8, Sep 2026). Claude Generated.
#
# dump_qm_gradient prints the analytic gradient of `hf` / `hf-3c` and the central
# finite difference of the same energy. check_qm_gradient.py gates
#   (a) analytic vs FD of curcuma's own energy (self-consistency of the derivative),
#   (b) analytic vs PySCF's analytic RHF gradient (+ simple-dftd3 D3/gCP gradients
#       for hf-3c), stored in reference.json by scripts/hf_gradient_reference.py.
# Run:  ctest -L qm_grad --output-on-failure

add_executable(dump_qm_gradient dump_qm_gradient.cpp)
target_link_libraries(dump_qm_gradient PRIVATE curcuma_core)
target_include_directories(dump_qm_gradient PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(dump_qm_gradient PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_QMGRAD_CASES H2O_svp NH3_svp CH4_svp HF_svp LiH_svp H2O_distorted_svp
        formaldehyde_distorted_svp HCN_svp H2O_hf3c NH3_hf3c H2O_distorted_hf3c
        formaldehyde_distorted_hf3c water_dimer_hf3c BF3_hf3c LiF_hf3c He_CH4_hf3c benzene_hf3c)
    foreach(_case ${_QMGRAD_CASES})
        add_test(NAME qm_grad_${_case}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gradient.py
                --dump $<TARGET_FILE:dump_qm_gradient>
                --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json
                --case ${_case})
        set_tests_properties(qm_grad_${_case} PROPERTIES TIMEOUT 3000 LABELS "qm_grad")
    endforeach()

    # -opt with hf-3c must reach the independent (PySCF + simple-dftd3) minimum.
    foreach(_mol H2O formaldehyde)
        add_test(NAME qm_opt_hf3c_${_mol}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_opt.py
                --curcuma $<TARGET_FILE:curcuma>
                --start ${CMAKE_CURRENT_SOURCE_DIR}/${_mol}_distorted.xyz
                --ref ${CMAKE_CURRENT_SOURCE_DIR}/${_mol}_hf3c_min_reference.xyz)
        set_tests_properties(qm_opt_hf3c_${_mol} PROPERTIES TIMEOUT 1500 LABELS "qm_grad")
    endforeach()
endif()
