# Native QM on the GPU, CPU-vs-GPU. Claude Generated (Oct 2026).
# Skipped (exit 77) when libcurcuma_cuda.so or a device is missing.
#   qm_gpu_jk_*     energy with -gpu cuda (Q-G2 device ERI is the default since Oct 2026)
#   qm_gpu_jk_upload_H2O  the Q-G1 path (CPU tensor uploaded, -qm.gpu_eri false)
#   qm_gpu_eri_*    device ERI tensor element-wise vs the CPU tensor (<= 1e-12) + energy (<= 1e-9)
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
        add_test(NAME qm_gpu_eri_${_name}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_eri.py
                --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/${_xyz})
        set_tests_properties(qm_gpu_eri_${_name} PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
    endforeach()
    add_test(NAME qm_gpu_jk_upload_H2O
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz
            "--extra=-qm.gpu_eri false")
    set_tests_properties(qm_gpu_jk_upload_H2O PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
    # cartesian 6d (no transform) and unscreened (every quartet) variants
    add_test(NAME qm_gpu_eri_H2O_cartesian
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_eri.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz
            "--extra=-qm.cartesian_d true")
    add_test(NAME qm_gpu_eri_formaldehyde_noscreen
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_eri.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/formaldehyde_distorted.xyz
            "--extra=-qm.eri_screening 0")
    # HF-3c (MINIX basis) through the same device path
    add_test(NAME qm_gpu_eri_H2O_hf3c
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_eri.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz
            --method hf-3c --basis MINIX)
    set_tests_properties(qm_gpu_eri_H2O_cartesian qm_gpu_eri_formaldehyde_noscreen qm_gpu_eri_H2O_hf3c
        PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
endif()
