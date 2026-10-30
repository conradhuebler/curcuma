# Kohn-Sham DFT on the GPU vs the CPU of the same binary. Claude Generated (Oct 2026, GPU KS-DFT G1-G4).
# Skipped (exit 77) when libcurcuma_cuda.so or a device is missing. Label qm_gpu_dft (and gpu).
#   qm_gpu_dft_xc_*       G1: V_xc / E_xc on the device vs the CPU from the same density (rel <= 1e-11 per build)
#   qm_gpu_dft_scf_*      G2: SCF energy vs the CPU (<= 1e-9 Eh), iterations +-1; device-resident loop, device-direct
#                         J/K, RI-J / RI-JK (host loop with device V_xc), streamed small blocks
#   qm_gpu_dft_grad_*     G3: XC gradient parts and the 2e part vs the CPU from the same density (<= 1e-9 Eh/Bohr),
#                         total gradient vs the CPU run (<= 1e-8 Eh/Bohr)
#   qm_gpu_dft_opt_*      G3: -opt with -gpu cuda reaches the minimum of the CPU -opt (<= 1e-8 Eh, max |dr_ij| <= 1e-4 A)
#   qm_gpu_dft_fallback_* G4: too little device memory for one block -> counted qm_xc fallback, CPU energy;
#                         -gpu_strict true -> exit code 3
#   qm_gpu_dft_*screen*   AO screening (-qm.xc_ao_threshold, default on): unscreened device path and small groups
#   qm_gpu_dft_wgrad_*    Becke weight derivative device vs CPU, incl. 1410 atoms (test_qm_gpu_wgrad.cpp)
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_QD ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_gpu_dft.py)
    set(_QD_H2O ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/H2O_distorted.xyz)
    set(_QD_NH3 ${CMAKE_CURRENT_SOURCE_DIR}/../structures/small/ammonia.ideal-c3v.xyz)
    set(_QD_BENZ ${CMAKE_CURRENT_SOURCE_DIR}/../qm_hf3c/benzene.xyz)
    set(_QD_FORM ${CMAKE_CURRENT_SOURCE_DIR}/../qm_grad/formaldehyde_distorted.xyz)
    set(_qd_tests "")
    # G1: every functional and VWN variant on H2O, NH3, benzene
    foreach(_mol H2O NH3 BENZ)
        foreach(_f lda:vwn5 lda:vwn_rpa pbe:vwn5 b3lyp:vwn5 b3lyp:vwn_rpa)
            string(REPLACE ":" ";" _p ${_f})
            list(GET _p 0 _m)
            list(GET _p 1 _v)
            string(TOLOWER ${_mol} _ml)
            set(_name qm_gpu_dft_xc_${_ml}_${_m}_${_v})
            add_test(NAME ${_name} COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma> --xyz ${_QD_${_mol}}
                --mode xc --method ${_m} "--extra=-qm.vwn_variant ${_v}")
            list(APPEND _qd_tests ${_name})
        endforeach()
    endforeach()
    # G2: SCF energies
    add_test(NAME qm_gpu_dft_scf_benzene_pbe_devloop COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_BENZ} --mode scf --method pbe --loop device "--dev-extra=-qm.gpu_scf on")
    add_test(NAME qm_gpu_dft_scf_benzene_b3lyp_devloop COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_BENZ} --mode scf --method b3lyp --loop device "--dev-extra=-qm.gpu_scf on")
    add_test(NAME qm_gpu_dft_scf_h2o_lda_devloop COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_H2O} --mode scf --method lda --loop device "--dev-extra=-qm.gpu_scf on")
    add_test(NAME qm_gpu_dft_scf_benzene_b3lyp_direct COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_BENZ} --mode scf --method b3lyp --loop device "--dev-extra=-qm.gpu_scf on -qm.gpu_direct on"
        "--cpu-extra=-qm.scf_direct on"
        # +-2 iterations: device direct J/K uses atomics, so the count varies from run to run at the same energy
        --tol-it 2)
    add_test(NAME qm_gpu_dft_scf_benzene_pbe_rij COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_BENZ} --mode scf --method pbe --loop host "--extra=-qm.ri j"
        # +-2 iterations: the CPU RI-J SCF of benzene ends at dP ~ 1e-10, on the convergence threshold, so the screened
        # V_xc (change ~1e-13) moves it by two iterations (docs/QM_GPU_ROADMAP.md section 4.1); energies agree.
        --tol-it 2)
    add_test(NAME qm_gpu_dft_scf_benzene_b3lyp_rijk COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_BENZ} --mode scf --method b3lyp --loop host "--extra=-qm.ri jk")
    add_test(NAME qm_gpu_dft_scf_formaldehyde_pbe_blocks COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_FORM} --mode scf --method pbe --loop device "--dev-extra=-qm.gpu_scf on -qm.gpu_xc_block_points 96")
    list(APPEND _qd_tests qm_gpu_dft_scf_benzene_pbe_devloop qm_gpu_dft_scf_benzene_b3lyp_devloop qm_gpu_dft_scf_h2o_lda_devloop
        qm_gpu_dft_scf_benzene_b3lyp_direct qm_gpu_dft_scf_benzene_pbe_rij qm_gpu_dft_scf_benzene_b3lyp_rijk
        qm_gpu_dft_scf_formaldehyde_pbe_blocks)
    # G3: gradients (LDA/PBE: kx = 0, B3LYP: kx = 0.2), without grid response, small blocks
    foreach(_case h2o_lda:H2O:lda h2o_pbe:H2O:pbe h2o_b3lyp:H2O:b3lyp formaldehyde_pbe:FORM:pbe benzene_b3lyp:BENZ:b3lyp)
        string(REPLACE ":" ";" _p ${_case})
        list(GET _p 0 _n)
        list(GET _p 1 _mol)
        list(GET _p 2 _m)
        add_test(NAME qm_gpu_dft_grad_${_n} COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_${_mol}} --mode grad --method ${_m})
        list(APPEND _qd_tests qm_gpu_dft_grad_${_n})
    endforeach()
    add_test(NAME qm_gpu_dft_grad_formaldehyde_b3lyp_noresponse COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_FORM} --mode grad --method b3lyp "--extra=-qm.grid_response false")
    add_test(NAME qm_gpu_dft_grad_h2o_pbe_blocks COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_H2O} --mode grad --method pbe "--dev-extra=-qm.gpu_xc_block_points 64")
    list(APPEND _qd_tests qm_gpu_dft_grad_formaldehyde_b3lyp_noresponse qm_gpu_dft_grad_h2o_pbe_blocks)
    # G3: -opt with the device gradient reaches the CPU minimum
    foreach(_case h2o_pbe:H2O:pbe h2o_b3lyp:H2O:b3lyp formaldehyde_b3lyp:FORM:b3lyp)
        string(REPLACE ":" ";" _p ${_case})
        list(GET _p 0 _n)
        list(GET _p 1 _mol)
        list(GET _p 2 _m)
        add_test(NAME qm_gpu_dft_opt_${_n} COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_${_mol}} --mode opt --method ${_m})
        list(APPEND _qd_tests qm_gpu_dft_opt_${_n})
    endforeach()
    # G4: not enough device memory for one block of points -> counted CPU fallback (and -gpu_strict)
    add_test(NAME qm_gpu_dft_fallback_h2o_memory COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_H2O} --mode fallback --method pbe "--dev-extra=-qm.gpu_xc_max_memory_mb 0.01")
    add_test(NAME qm_gpu_dft_fallback_h2o_strict COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_H2O} --mode fallback --method pbe --strict "--dev-extra=-qm.gpu_xc_max_memory_mb 0.01")
    list(APPEND _qd_tests qm_gpu_dft_fallback_h2o_memory qm_gpu_dft_fallback_h2o_strict)
    # AO screening (Oct 2026, -qm.xc_ao_threshold): unscreened device path (threshold 0) and screened groups of
    # 1024 points vs the CPU of the same density / the CPU run
    add_test(NAME qm_gpu_dft_xc_benzene_pbe_noscreen COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_BENZ} --mode xc --method pbe "--extra=-qm.xc_ao_threshold 0")
    add_test(NAME qm_gpu_dft_grad_benzene_pbe_screen_groups COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_BENZ} --mode grad --method pbe "--dev-extra=-qm.gpu_xc_block_points 1024")
    add_test(NAME qm_gpu_dft_grad_benzene_b3lyp_noscreen COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
        --xyz ${_QD_BENZ} --mode grad --method b3lyp "--extra=-qm.xc_ao_threshold 0")
    list(APPEND _qd_tests qm_gpu_dft_xc_benzene_pbe_noscreen qm_gpu_dft_grad_benzene_pbe_screen_groups
        qm_gpu_dft_grad_benzene_b3lyp_noscreen)
    set_tests_properties(${_qd_tests} PROPERTIES TIMEOUT 1500 LABELS "qm_gpu_dft;gpu" SKIP_RETURN_CODE 77)
endif()

# Becke weight derivative on the device vs qmgrid::gridWeightGradient, also above the former 700-atom limit
# (Claude Generated, Oct 2026). Structures from the library; small grids keep the O(npts natoms^2) CPU reference
# short (peo201-chain 1410 atoms, 4x6 grid: 33840 points, f nonzero on every 16th = 2115 points). rel = max|g_dev - g_cpu| / max|g_cpu| <= 1e-12.
add_executable(test_qm_gpu_wgrad test_qm_gpu_wgrad.cpp)
# whole archive: the CUDA plugin binds its undefined core symbols to the executable at dlopen time (as for curcuma,
# built -rdynamic); the archive has header-defined non-inline functions in two objects (rmsd/munkres.h), which a
# whole-archive link sees twice, hence --allow-multiple-definition for this test executable only
target_link_libraries(test_qm_gpu_wgrad PRIVATE "$<LINK_LIBRARY:WHOLE_ARCHIVE,curcuma_core>")
target_link_options(test_qm_gpu_wgrad PRIVATE "LINKER:--allow-multiple-definition")
target_include_directories(test_qm_gpu_wgrad PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_gpu_wgrad PRIVATE cxx_std_17)
set(_WG_PEO ${CMAKE_CURRENT_SOURCE_DIR}/../structures/bulk/peo201-chain.xyz)
set(_WG_CAF ${CMAKE_CURRENT_SOURCE_DIR}/../structures/medium/caffeine.xyz)
add_test(NAME qm_gpu_dft_wgrad_caffeine_sg0 COMMAND test_qm_gpu_wgrad --xyz ${_WG_CAF} --grid sg0 --threads 4)
add_test(NAME qm_gpu_dft_wgrad_peo201_becke COMMAND test_qm_gpu_wgrad --xyz ${_WG_PEO} --grid 4x6 --stride 16 --threads 8)
add_test(NAME qm_gpu_dft_wgrad_peo201_stratmann COMMAND test_qm_gpu_wgrad --xyz ${_WG_PEO} --grid 4x6 --partition stratmann
    --stride 16 --threads 8)
# the plugin is found next to the executable or on the loader path: point the loader at the curcuma binary's directory
set_tests_properties(qm_gpu_dft_wgrad_caffeine_sg0 qm_gpu_dft_wgrad_peo201_becke qm_gpu_dft_wgrad_peo201_stratmann
    PROPERTIES TIMEOUT 1500 LABELS "qm_gpu_dft;gpu" SKIP_RETURN_CODE 77
    ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:curcuma>:$ENV{LD_LIBRARY_PATH}")
