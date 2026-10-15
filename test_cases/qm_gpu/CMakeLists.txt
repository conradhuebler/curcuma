# Native QM on the GPU (Q-G1, Oct 2026): CPU-vs-GPU energy. Claude Generated.
# Skipped (exit 77) when libcurcuma_cuda.so or a device is missing.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    foreach(_case H2O:../qm_grad/H2O_distorted.xyz formaldehyde:../qm_grad/formaldehyde_distorted.xyz benzene:../qm_hf3c/benzene.xyz)
        string(REPLACE ":" ";" _p ${_case})
        list(GET _p 0 _name)
        list(GET _p 1 _xyz)
        add_test(NAME qm_gpu_jk_${_name}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu.py
                --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/${_xyz})
        set_tests_properties(qm_gpu_jk_${_name} PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
    endforeach()
endif()
