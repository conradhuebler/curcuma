# Native double hybrids (WP12 stage 1): P86/B95 kernels vs libxc, KS part and PT2 components vs PySCF on the same
# grid, D3(BJ) vs ORCA 6.1, CLI behaviour. Claude Generated (Oct 2026).
#
# test_qm_dh (C++) prints JSON; check_qm_dh.py compares it with reference.json (scripts/qm_dh_pyscf_reference.py,
# PySCF 2.14 / libxc 7.0.0) and reference_orca.json (ORCA 6.1 D3(BJ) energies, scripts/qm_dh_orca_compare.py).
#   qm_dh_kernel                 P86 (PZ81 and VWN5 local), B95: eps, v_rho, v_sigma (v_tau) vs libxc; v vs FD
#   qm_dh_energy_<case>          E_KS, E_xc, E_os, E_ss (RI, def2-universal-JKFIT) vs PySCF on the same grid
#   qm_dh_fc_<case>              the same with -qm.frozen_core true
#   qm_dh_canonical_<case>       -qm.dh_pt2 canonical vs pyscf.mp.MP2 on the KS orbitals
#   qm_dh_d3_<case>              D3(BJ) term of the -d3bj names vs ORCA 6.1 (two-body)
#   qm_dh_{gradient,update,threads}_<fn>   gradient request is an error; UpdateMolecule; 1 vs 4 threads bitwise
#   qm_dh_cli_*                  curcuma -sp prints the energy, -opt stops with an error, D4 request is an error
#   qm_dh_gpu_H2O                -gpu cuda: counted qm_xc fallback, same energy, -gpu_strict exits; skip (77) without
#                                a CUDA plugin
# Run:  ctest -L qm_dh --output-on-failure

add_executable(test_qm_dh test_qm_dh.cpp)
target_link_libraries(test_qm_dh PRIVATE curcuma_core)
target_include_directories(test_qm_dh PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_dh PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_DH_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_dh.py
        --exe $<TARGET_FILE:test_qm_dh> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json
        --orca ${CMAKE_CURRENT_SOURCE_DIR}/reference_orca.json --curcuma $<TARGET_FILE:curcuma>)
    add_test(NAME qm_dh_kernel COMMAND ${_DH_CHECK} --check kernel)
    set_tests_properties(qm_dh_kernel PROPERTIES TIMEOUT 120 LABELS "qm_dh")
    foreach(_fn b2plyp b2gp-plyp dsd-blyp dsd-pbep86 dsd-pbeb95-d3bj revdsd-pbep86)
        foreach(_mol H2O CH4 NH3 HF LiH Ne benzene)
            if(_fn STREQUAL "dsd-pbeb95-d3bj" AND _mol STREQUAL "LiH")
                continue()  # no PySCF reference: its SCF stalls at |g| ~ 6e-5 with MGGA_C_BC95 (WP12 doc)
            endif()
            add_test(NAME qm_dh_energy_${_mol}_${_fn}_sg1 COMMAND ${_DH_CHECK} --check energy --case ${_mol}_${_fn}_sg1)
            set_tests_properties(qm_dh_energy_${_mol}_${_fn}_sg1 PROPERTIES TIMEOUT 900 LABELS "qm_dh")
        endforeach()
    endforeach()
    foreach(_fn dsd-blyp-d3bj dsd-pbep86-d3bj dsd-pbep86@vwn5)
        foreach(_mol H2O NH3 benzene)
            add_test(NAME qm_dh_energy_${_mol}_${_fn}_sg1 COMMAND ${_DH_CHECK} --check energy --case ${_mol}_${_fn}_sg1)
            set_tests_properties(qm_dh_energy_${_mol}_${_fn}_sg1 PROPERTIES TIMEOUT 900 LABELS "qm_dh")
        endforeach()
    endforeach()
    foreach(_case H2O_b2plyp_sg3 H2O_dsd-pbeb95-d3bj_sg3)
        add_test(NAME qm_dh_energy_${_case} COMMAND ${_DH_CHECK} --check energy --case ${_case})
        set_tests_properties(qm_dh_energy_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh")
    endforeach()
    foreach(_case H2O_b2plyp_sg1 NH3_dsd-pbep86_sg1 benzene_dsd-blyp_sg1 benzene_dsd-pbeb95-d3bj_sg1 LiH_b2gp-plyp_sg1)
        add_test(NAME qm_dh_fc_${_case} COMMAND ${_DH_CHECK} --check fc --case ${_case})
        set_tests_properties(qm_dh_fc_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh")
    endforeach()
    foreach(_case H2O_b2plyp_sg1 benzene_dsd-pbep86_sg1 H2O_dsd-pbeb95-d3bj_sg1)
        add_test(NAME qm_dh_canonical_${_case} COMMAND ${_DH_CHECK} --check canonical --case ${_case})
        set_tests_properties(qm_dh_canonical_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh")
    endforeach()
    foreach(_case H2O_b2plyp-d3bj H2O_b2gp-plyp-d3bj H2O_dsd-blyp-d3bj benzene_dsd-pbep86-d3bj
                  benzene_dsd-pbeb95-d3bj NH3_b2plyp-d3bj)
        add_test(NAME qm_dh_d3_${_case} COMMAND ${_DH_CHECK} --check d3 --case ${_case})
        set_tests_properties(qm_dh_d3_${_case} PROPERTIES TIMEOUT 900 LABELS "qm_dh")
    endforeach()
    foreach(_fn b2plyp dsd-pbeb95-d3bj)
        foreach(_chk gradient update threads)
            add_test(NAME qm_dh_${_chk}_${_fn} COMMAND ${_DH_CHECK} --check ${_chk} --case H2O_${_fn})
            set_tests_properties(qm_dh_${_chk}_${_fn} PROPERTIES TIMEOUT 600 LABELS "qm_dh")
        endforeach()
    endforeach()
    foreach(_chk cli_sp cli_opt cli_d4 cli_unknown_option)
        add_test(NAME qm_dh_${_chk} COMMAND ${_DH_CHECK} --check ${_chk} --case H2O_b2plyp)
        set_tests_properties(qm_dh_${_chk} PROPERTIES TIMEOUT 600 LABELS "qm_dh")
    endforeach()
    add_test(NAME qm_dh_gpu_H2O COMMAND ${_DH_CHECK} --check gpu --case H2O_dsd-pbeb95-d3bj)
    set_tests_properties(qm_dh_gpu_H2O PROPERTIES TIMEOUT 600 LABELS "qm_dh" SKIP_RETURN_CODE 77)
endif()
