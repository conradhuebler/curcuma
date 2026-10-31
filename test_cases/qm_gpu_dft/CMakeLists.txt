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
#   *_tpss*, *_scan*      meta-GGA on the device (WP11 stage 3), only registered with USE_CUDA
#   *_dh_*, DH names      double hybrids on the device (WP12 stage 3, USE_CUDA only): P86/B95 V_xc, DH energies with
#                         -gpu_strict true, gradient parts incl. F_xc / XC response / PT2 2e (qm_dh_grad_pt2 allowed)
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
    # Meta-GGA TPSS/SCAN on the device (WP11 stage 3, Claude Generated, Oct 2026); registered only with a CUDA build.
    # Same gates as LDA/GGA/B3LYP; SCAN is grid-sensitive, so every comparison is CPU vs GPU on the SAME grid.
    if(USE_CUDA)
        set(_qm_tests "")
        foreach(_mol H2O NH3 BENZ)
            foreach(_m tpss scan)
                string(TOLOWER ${_mol} _ml)
                add_test(NAME qm_gpu_dft_xc_${_ml}_${_m} COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
                    --xyz ${_QD_${_mol}} --mode xc --method ${_m})
                list(APPEND _qm_tests qm_gpu_dft_xc_${_ml}_${_m})
            endforeach()
        endforeach()
        add_test(NAME qm_gpu_dft_xc_h2o_scan_fine COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_H2O} --mode xc --method scan "--extra=-qm.grid 150x974")
        add_test(NAME qm_gpu_dft_xc_benzene_tpss_noscreen COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_BENZ} --mode xc --method tpss "--extra=-qm.xc_ao_threshold 0")
        list(APPEND _qm_tests qm_gpu_dft_xc_h2o_scan_fine qm_gpu_dft_xc_benzene_tpss_noscreen)
        # SCF energies, device-resident loop; -gpu_strict true must pass (no counted fallback left for tpss/scan)
        foreach(_case benzene_tpss:BENZ:tpss benzene_scan:BENZ:scan nh3_scan_fine:NH3:scan)
            string(REPLACE ":" ";" _p ${_case})
            list(GET _p 0 _n)
            list(GET _p 1 _mol)
            list(GET _p 2 _m)
            set(_grid "")
            if(_n MATCHES "_fine")
                set(_grid "-qm.grid 150x974")
            endif()
            add_test(NAME qm_gpu_dft_scf_${_n}_devloop COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
                --xyz ${_QD_${_mol}} --mode scf --method ${_m} --loop device "--extra=${_grid}"
                "--dev-extra=-qm.gpu_scf on -gpu_strict true")
            list(APPEND _qm_tests qm_gpu_dft_scf_${_n}_devloop)
        endforeach()
        add_test(NAME qm_gpu_dft_scf_benzene_tpss_rijk COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_BENZ} --mode scf --method tpss --loop host "--extra=-qm.ri jk" "--dev-extra=-gpu_strict true")
        list(APPEND _qm_tests qm_gpu_dft_scf_benzene_tpss_rijk)
        # Gradients (tau term, grid response, small blocks, screened groups, unscreened path)
        foreach(_case h2o_tpss:H2O:tpss h2o_scan:H2O:scan formaldehyde_tpss:FORM:tpss benzene_scan:BENZ:scan)
            string(REPLACE ":" ";" _p ${_case})
            list(GET _p 0 _n)
            list(GET _p 1 _mol)
            list(GET _p 2 _m)
            add_test(NAME qm_gpu_dft_grad_${_n} COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
                --xyz ${_QD_${_mol}} --mode grad --method ${_m})
            list(APPEND _qm_tests qm_gpu_dft_grad_${_n})
        endforeach()
        add_test(NAME qm_gpu_dft_grad_h2o_scan_noresponse COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_H2O} --mode grad --method scan "--extra=-qm.grid_response false")
        add_test(NAME qm_gpu_dft_grad_h2o_tpss_blocks COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_H2O} --mode grad --method tpss "--dev-extra=-qm.gpu_xc_block_points 64")
        add_test(NAME qm_gpu_dft_grad_benzene_tpss_screen_groups COMMAND ${Python3_EXECUTABLE} ${_QD}
            --curcuma $<TARGET_FILE:curcuma> --xyz ${_QD_BENZ} --mode grad --method tpss "--dev-extra=-qm.gpu_xc_block_points 1024")
        add_test(NAME qm_gpu_dft_grad_benzene_tpss_noscreen COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_BENZ} --mode grad --method tpss "--extra=-qm.xc_ao_threshold 0")
        add_test(NAME qm_gpu_dft_opt_h2o_tpss COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_H2O} --mode opt --method tpss)
        add_test(NAME qm_gpu_dft_fallback_h2o_scan_memory COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_H2O} --mode fallback --method scan "--dev-extra=-qm.gpu_xc_max_memory_mb 0.01")
        list(APPEND _qm_tests qm_gpu_dft_grad_h2o_scan_noresponse qm_gpu_dft_grad_h2o_tpss_blocks
            qm_gpu_dft_grad_benzene_tpss_screen_groups qm_gpu_dft_grad_benzene_tpss_noscreen qm_gpu_dft_opt_h2o_tpss
            qm_gpu_dft_fallback_h2o_scan_memory)
        list(APPEND _qd_tests ${_qm_tests})
        # Double hybrids on the device (WP12 stage 3, Claude Generated, Oct 31 2026): P86 (VWN5 default, PZ81) through the
        # GGA path, B95 through the meta-GGA path, energies with -gpu_strict true; gradient parts F_xc (Z-vector), XC
        # response gradient, PT2 1e/2e/RI-3-index on the device. The PT2 intermediates stay on the host as the counted
        # qm_dh_grad_pt2 fallback, the only one the gradient cases allow; -gpu_strict true with a gradient exits 3 on it.
        set(_dh_tests "")
        foreach(_case h2o_dsd-pbep86:H2O:dsd-pbep86: h2o_dsd-pbep86_pz81:H2O:dsd-pbep86:-qm.p86_lda=pz81
                      h2o_dsd-pbeb95-d3bj:H2O:dsd-pbeb95-d3bj: benzene_dsd-pbep86-d3bj:BENZ:dsd-pbep86-d3bj:
                      benzene_dsd-pbeb95-d3bj:BENZ:dsd-pbeb95-d3bj:)
            string(REPLACE ":" ";" _p "${_case}")
            list(GET _p 0 _n)
            list(GET _p 1 _mol)
            list(GET _p 2 _m)
            list(LENGTH _p _len)
            set(_x "")
            if(_len GREATER 3)
                list(GET _p 3 _x)
                string(REPLACE "=" " " _x "${_x}")
            endif()
            add_test(NAME qm_gpu_dft_xc_${_n} COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
                --xyz ${_QD_${_mol}} --mode xc --method ${_m} "--extra=${_x}")
            list(APPEND _dh_tests qm_gpu_dft_xc_${_n})
        endforeach()
        foreach(_case benzene_dsd-pbep86-d3bj:BENZ:dsd-pbep86-d3bj benzene_dsd-pbeb95-d3bj:BENZ:dsd-pbeb95-d3bj
                      nh3_b2plyp:NH3:b2plyp)
            string(REPLACE ":" ";" _p ${_case})
            list(GET _p 0 _n)
            list(GET _p 1 _mol)
            list(GET _p 2 _m)
            add_test(NAME qm_gpu_dft_scf_${_n}_devloop COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
                --xyz ${_QD_${_mol}} --mode scf --method ${_m} --loop device --dh "--dev-extra=-qm.gpu_scf on -gpu_strict true")
            list(APPEND _dh_tests qm_gpu_dft_scf_${_n}_devloop)
        endforeach()
        foreach(_case h2o_dsd-pbep86-d3bj:H2O:dsd-pbep86-d3bj:
                      h2o_dsd-pbeb95-d3bj_canonical:H2O:dsd-pbeb95-d3bj:-qm.dh_pt2=canonical
                      formaldehyde_b2plyp_noresponse:FORM:b2plyp:-qm.grid_response=false
                      nh3_revdsd-pbep86_nofc:NH3:revdsd-pbep86:-qm.frozen_core=false
                      benzene_dsd-pbeb95-d3bj:BENZ:dsd-pbeb95-d3bj:)
            string(REPLACE ":" ";" _p "${_case}")
            list(GET _p 0 _n)
            list(GET _p 1 _mol)
            list(GET _p 2 _m)
            list(LENGTH _p _len)
            set(_x "")
            if(_len GREATER 3)
                list(GET _p 3 _x)
                string(REPLACE "=" " " _x "${_x}")
            endif()
            add_test(NAME qm_gpu_dft_grad_dh_${_n} COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
                --xyz ${_QD_${_mol}} --mode grad --method ${_m} --dh --allow-fallback qm_dh_grad_pt2 "--extra=${_x}")
            list(APPEND _dh_tests qm_gpu_dft_grad_dh_${_n})
        endforeach()
        add_test(NAME qm_gpu_dft_opt_dh_h2o_dsd-pbep86-d3bj COMMAND ${Python3_EXECUTABLE} ${_QD} --curcuma $<TARGET_FILE:curcuma>
            --xyz ${_QD_H2O} --mode opt --method dsd-pbep86-d3bj --allow-fallback qm_dh_grad_pt2)
        add_test(NAME qm_gpu_dft_fallback_dh_h2o_b2plyp_grad_strict COMMAND ${Python3_EXECUTABLE} ${_QD}
            --curcuma $<TARGET_FILE:curcuma> --xyz ${_QD_H2O} --mode fallback --method b2plyp --strict --dh
            --fallback-name qm_dh_grad_pt2)
        list(APPEND _dh_tests qm_gpu_dft_opt_dh_h2o_dsd-pbep86-d3bj qm_gpu_dft_fallback_dh_h2o_b2plyp_grad_strict)
        list(APPEND _qd_tests ${_dh_tests})
    endif()
    set_tests_properties(${_qd_tests} PROPERTIES TIMEOUT 1500 LABELS "qm_gpu_dft;gpu" SKIP_RETURN_CODE 77)
    if(USE_CUDA)
        set_tests_properties(${_dh_tests} PROPERTIES LABELS "qm_gpu_dft;qm_dh;gpu")
    endif()
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
