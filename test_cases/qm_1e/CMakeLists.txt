# Native QM 1-electron integral validation suite (WP1).
# Claude Generated (WP1, June 2026).
#
# dump_qm_1e builds the curcuma native-QM basis + 1e integrals (S, T, V,
# Hc=T+V) and emits them as JSON. diff_qm_1e.py orchestrates the comparison:
#   (a) element-wise S/T/V vs the independent Python OS witness (cartesian,
#       same AO order) -- the kernel gate, no external deps;
#   (b) internal consistency: symmetry, H=T+V, S_ii=1, Tr(P*S)=1, nbf;
#   (c) optional generalized-eigenvalue MO energies vs an ORCA reference JSON
#       (<mol>.orca_ref.json next to the xyz, from scripts/qm_1e_reference.py)
#       -- auto-discovered by the script, so a checkout without ORCA still
#       runs (a)+(b).
#
# Run:  ctest -L qm_1e --output-on-failure

add_executable(dump_qm_1e dump_qm_1e.cpp)
target_link_libraries(dump_qm_1e PRIVATE curcuma_core)
target_include_directories(dump_qm_1e PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(dump_qm_1e PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_QM1E_PY    ${CMAKE_CURRENT_SOURCE_DIR}/diff_qm_1e.py)
    set(_QM1E_WIT   ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/qm_1e_python_ints.py)
    set(_QM1E_MOLS  H2 He LiH BeH2 BH CH4 NH3 H2O HF Ne)
    set(_QM1E_TOL   1e-10)

    foreach(_mol ${_QM1E_MOLS})
        add_test(NAME qm_1e_${_mol}
            COMMAND ${Python3_EXECUTABLE} ${_QM1E_PY}
                --dump $<TARGET_FILE:dump_qm_1e>
                --python ${_QM1E_WIT}
                --xyz ${CMAKE_CURRENT_SOURCE_DIR}/${_mol}.xyz
                --basis def2-SVP
                --tol ${_QM1E_TOL}
                --quiet)
        set_tests_properties(qm_1e_${_mol}
            PROPERTIES TIMEOUT 120 LABELS "qm_1e")
    endforeach()
endif()