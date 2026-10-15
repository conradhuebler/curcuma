# Native QM on the GPU, CPU-vs-GPU. Claude Generated (Oct 2026).
# Skipped (exit 77) when libcurcuma_cuda.so or a device is missing.
#   qm_gpu_jk_*     energy with -gpu cuda (Q-G2 device ERI is the default since Oct 2026)
#   qm_gpu_jk_upload_H2O  the Q-G1 path (CPU tensor uploaded, -qm.gpu_eri false)
#   qm_gpu_eri_*    device ERI tensor element-wise vs the CPU tensor (<= 1e-12) + energy (<= 1e-9)
#   qm_gpu_direct_* integral-direct J/K on the device (Q-G3): J/K of one build vs CPU DirectJK,
#                   stored device and stored CPU tensor (<= 1e-10), quartet counts, E_elec (<= 1e-9)
#   qm_gpu_scf_*    device-resident SCF loop (Q-G4) vs the host loop and the CPU: E_elec (<= 1e-9),
#                   same iteration count as the host loop
#   qm_gpu_grad_*   two-electron gradient on the device (Q-G7) vs the CPU 2e gradient (<= 1e-9 Eh/Bohr),
#                   its translational invariance, total gradient vs the CPU run, FD on H2O
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
    # Q-G3 (Claude Generated, Oct 2026): forced device-direct (-qm.gpu_direct on)
    foreach(_case H2O:../qm_grad/H2O_distorted.xyz formaldehyde:../qm_grad/formaldehyde_distorted.xyz benzene:../qm_hf3c/benzene.xyz)
        string(REPLACE ":" ";" _p ${_case})
        list(GET _p 0 _name)
        list(GET _p 1 _xyz)
        add_test(NAME qm_gpu_direct_${_name}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_direct.py
                --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/${_xyz}
                "--extra=-qm.gpu_direct on")
        set_tests_properties(qm_gpu_direct_${_name} PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
    endforeach()
    # stored tensor "does not fit" (device limited to 100 MB): scf_direct auto must go direct on the
    # device, without a counted fallback; CPU reference integral-direct as well
    add_test(NAME qm_gpu_direct_benzene_nofit
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_direct.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_hf3c/benzene.xyz
            "--extra=-qm.gpu_max_memory_mb 100 -qm.eri_max_memory_mb 100" "--cpu-extra=-qm.scf_direct on")
    # cartesian 6d, unscreened, HF-3c/MINIX through the device-direct path
    add_test(NAME qm_gpu_direct_H2O_cartesian
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_direct.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz
            "--extra=-qm.gpu_direct on -qm.cartesian_d true" "--cpu-extra=-qm.cartesian_d true")
    add_test(NAME qm_gpu_direct_formaldehyde_noscreen
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_direct.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/formaldehyde_distorted.xyz
            "--extra=-qm.gpu_direct on -qm.eri_screening 0" "--cpu-extra=-qm.eri_screening 0")
    add_test(NAME qm_gpu_direct_H2O_hf3c
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_direct.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz
            --method hf-3c --basis MINIX "--extra=-qm.gpu_direct on")
    set_tests_properties(qm_gpu_direct_benzene_nofit qm_gpu_direct_H2O_cartesian qm_gpu_direct_formaldehyde_noscreen
        qm_gpu_direct_H2O_hf3c PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
    # Q-G4 (Claude Generated, Oct 2026): device-resident SCF loop, stored device tensor (default path)
    foreach(_case H2O:../qm_grad/H2O_distorted.xyz formaldehyde:../qm_grad/formaldehyde_distorted.xyz benzene:../qm_hf3c/benzene.xyz)
        string(REPLACE ":" ";" _p ${_case})
        list(GET _p 0 _name)
        list(GET _p 1 _xyz)
        add_test(NAME qm_gpu_scf_${_name}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_scf.py
                --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/${_xyz})
        set_tests_properties(qm_gpu_scf_${_name} PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
    endforeach()
    # device integral-direct J/K (incremental Fock builds), cartesian 6d, HF-3c/MINIX, ADIIS, plain damping
    add_test(NAME qm_gpu_scf_benzene_direct
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_scf.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_hf3c/benzene.xyz
            "--extra=-qm.gpu_direct on" "--cpu-extra=-qm.scf_direct on")
    add_test(NAME qm_gpu_scf_H2O_cartesian
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_scf.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz
            "--extra=-qm.cartesian_d true" "--cpu-extra=-qm.cartesian_d true")
    add_test(NAME qm_gpu_scf_H2O_hf3c
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_scf.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz
            --method hf-3c --basis MINIX)
    add_test(NAME qm_gpu_scf_formaldehyde_adiis
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_scf.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/formaldehyde_distorted.xyz
            "--extra=-qm.scf_mode adiis" "--cpu-extra=-qm.scf_mode adiis")
    add_test(NAME qm_gpu_scf_H2O_plain
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_scf.py
            --curcuma $<TARGET_FILE:curcuma> --xyz ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz
            "--extra=-qm.scf_mode plain" "--cpu-extra=-qm.scf_mode plain")
    set_tests_properties(qm_gpu_scf_benzene_direct qm_gpu_scf_H2O_cartesian qm_gpu_scf_H2O_hf3c
        qm_gpu_scf_formaldehyde_adiis qm_gpu_scf_H2O_plain
        PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
    # Q-G7 (Claude Generated, Oct 2026): two-electron gradient on the device vs the CPU 2e gradient from
    # the same density (<= 1e-9 Eh/Bohr), translational invariance of the device 2e part, total
    # gradient of the device run vs the CPU run; SCF paths: stored device tensor (default), device-
    # direct J/K, device-resident SCF loop; cartesian 6d, unscreened, HF-3c/MINIX; FD on H2O
    set(_QMGG ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_grad.py)
    set(_QMGG_H2O ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz)
    set(_QMGG_FORM ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/formaldehyde_distorted.xyz)
    set(_QMGG_BENZ ${CMAKE_CURRENT_SOURCE_DIR}/../qm_hf3c/benzene.xyz)
    add_test(NAME qm_gpu_grad_H2O COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma> --xyz ${_QMGG_H2O} --fd)
    add_test(NAME qm_gpu_grad_formaldehyde COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma> --xyz ${_QMGG_FORM})
    add_test(NAME qm_gpu_grad_benzene COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma> --xyz ${_QMGG_BENZ})
    add_test(NAME qm_gpu_grad_benzene_direct COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QMGG_BENZ} "--dev-extra=-qm.gpu_direct on -qm.scf_threshold 1e-8" --tol-total 1e-7)
    add_test(NAME qm_gpu_grad_benzene_scf COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QMGG_BENZ} "--dev-extra=-qm.gpu_scf on")
    add_test(NAME qm_gpu_grad_H2O_upload COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QMGG_H2O} "--dev-extra=-qm.gpu_eri false")
    add_test(NAME qm_gpu_grad_H2O_cartesian COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QMGG_H2O} "--extra=-qm.cartesian_d true")
    add_test(NAME qm_gpu_grad_formaldehyde_noscreen COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QMGG_FORM} "--extra=-qm.eri_screening 0" --tol-2e 1e-12)
    add_test(NAME qm_gpu_grad_H2O_hf3c COMMAND ${Python3_EXECUTABLE} ${_QMGG} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QMGG_H2O} --method hf-3c --basis MINIX)
    set_tests_properties(qm_gpu_grad_H2O qm_gpu_grad_formaldehyde qm_gpu_grad_benzene qm_gpu_grad_benzene_direct
        qm_gpu_grad_benzene_scf qm_gpu_grad_H2O_upload qm_gpu_grad_H2O_cartesian qm_gpu_grad_formaldehyde_noscreen
        qm_gpu_grad_H2O_hf3c PROPERTIES TIMEOUT 1500 LABELS "qm_gpu;gpu" SKIP_RETURN_CODE 77)
endif()
