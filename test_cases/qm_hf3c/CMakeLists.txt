# Native HF-3c validation (Sep 2026). Claude Generated.
#
# check_hf3c.py runs `curcuma -sp <xyz> -method hf-3c` and compares every term
# (HF/MINIX, D3(BJ), gCP, SRB, total) with reference.json, which holds values
# from programs that share no code with curcuma: PySCF for HF/MINIX and the
# simple-dftd3 Python API for D3(BJ) and gCP/SRB (scripts/hf3c_reference.py).
# H2O is additionally checked against ORCA 6.1's `! HF-3c` total.
#
# Run:  ctest -L qm_hf3c --output-on-failure

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_HF3C_MOLS H2 He LiH BeH2 BH CH4 NH3 H2O HF Ne
                   benzene water_dimer formaldehyde HCN BF3 LiF He_CH4)
    foreach(_mol ${_HF3C_MOLS})
        add_test(NAME qm_hf3c_${_mol}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_hf3c.py
                --curcuma $<TARGET_FILE:curcuma>
                --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json
                --mol ${_mol})
        set_tests_properties(qm_hf3c_${_mol} PROPERTIES TIMEOUT 600 LABELS "qm_hf3c")
    endforeach()
endif()

# updateGeometry() must invalidate the integral/SCF caches (opt/scan/MD path).
add_executable(test_qm_update_geometry test_qm_update_geometry.cpp)
target_link_libraries(test_qm_update_geometry PRIVATE curcuma_core)
target_include_directories(test_qm_update_geometry PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_update_geometry PRIVATE cxx_std_17)
add_test(NAME qm_update_geometry COMMAND test_qm_update_geometry)
set_tests_properties(qm_update_geometry PROPERTIES TIMEOUT 120 LABELS "qm_hf3c")

# Basis-set citations of the native QM (Oct 2026): the literature of the basis sets a run loads is registered with
# CitationRegistry and printed at the end (qm_methods/qm_citations.h, entries in src/core/citation_database.cpp). H2O keeps it fast.
set(_QM_CITE_XYZ ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e/H2O.xyz)
add_test(NAME qm_cite_def2 COMMAND curcuma -sp ${_QM_CITE_XYZ} -method hf -verbosity 1)
set_tests_properties(qm_cite_def2 PROPERTIES TIMEOUT 120 LABELS "qm_hf3c;qm_cite"
    PASS_REGULAR_EXPRESSION "Weigend & Ahlrichs 2005.*Basis Set Exchange|Basis Set Exchange.*Weigend & Ahlrichs 2005")
add_test(NAME qm_cite_jkfit COMMAND curcuma -sp ${_QM_CITE_XYZ} -method hf -qm.ri jk -verbosity 1)
set_tests_properties(qm_cite_jkfit PROPERTIES TIMEOUT 120 LABELS "qm_hf3c;qm_cite"
    PASS_REGULAR_EXPRESSION "def2-universal-JKFIT; Weigend 2008")
add_test(NAME qm_cite_hf3c COMMAND curcuma -sp ${_QM_CITE_XYZ} -method hf-3c -verbosity 1)
set_tests_properties(qm_cite_hf3c PROPERTIES TIMEOUT 120 LABELS "qm_hf3c;qm_cite"
    PASS_REGULAR_EXPRESSION "Sure & Grimme 2013")
# Method literature (Oct 2026): MP2/RI-MP2 variants, CCSD, local MP2, integral-direct SCF.
foreach(_c "mp2;-method mp2;Moeller-Plesset perturbation theory \\(1934\\)"
           "scsrimp2;-method scs-ri-mp2;Spin-component-scaled MP2 \\(Grimme 2003\\)"
           "ccsd;-method ccsd;Purvis & Bartlett 1982"
           "lmp2;-method lmp2;Pseudonatural orbitals \\(Edmiston & Krauss 1966\\)"
           "direct;-method hf -qm.scf_direct on;Integral-direct SCF \\(Almloef")
    list(GET _c 0 _n)
    list(GET _c 1 _a)
    list(GET _c 2 _re)
    separate_arguments(_args UNIX_COMMAND "${_a}")
    add_test(NAME qm_cite_${_n} COMMAND curcuma -sp ${_QM_CITE_XYZ} ${_args} -verbosity 1)
    set_tests_properties(qm_cite_${_n} PROPERTIES TIMEOUT 300 LABELS "qm_hf3c;qm_cite" PASS_REGULAR_EXPRESSION "${_re}")
endforeach()
