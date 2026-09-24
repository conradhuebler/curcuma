# Native QM 2-electron integral validation suite (WP2).
# Claude Generated (WP2, June 2026).
#
# dump_qm_2e builds the curcuma native-QM basis + 4-centre ERI (cartesian 6d,
# chemists' (mu nu | lam sig), McMurchie-Davidson) via QMEngine::cartesianERI(), and
# the Coulomb J / exchange K matrices from a dummy closed-shell density P. It
# emits them as JSON. diff_qm_2e.py orchestrates the comparison:
#   (a) element-wise ERI vs the independent Python MD witness (cartesian, same
#       AO order) -- the kernel gate, no external deps;
#   (b) internal consistency: 8-fold ERI symmetry, J/K symmetric, Tr(P*S)=2,
#       Tr(P*J)==Tr(P*K) single-orbital identity, J/K element-wise vs witness;
#   (c) optional xcDFT He-VDZ ERI cross-check (gated by --xcDFT-ref, auto-skipped
#       if the reference file is absent, so a checkout without it still runs
#       (a)+(b)).
#
# The 10 molecules and .xyz geometries are shared with qm_1e (no new inputs).
# Run:  ctest -L qm_2e --output-on-failure

add_executable(dump_qm_2e dump_qm_2e.cpp)
target_link_libraries(dump_qm_2e PRIVATE curcuma_core)
target_include_directories(dump_qm_2e PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(dump_qm_2e PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_QM2E_PY    ${CMAKE_CURRENT_SOURCE_DIR}/diff_qm_2e.py)
    set(_QM2E_WIT   ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/qm_2e_python_ints.py)
    set(_QM2E_XYZ   ${CMAKE_CURRENT_SOURCE_DIR}/../qm_1e)
    set(_QM2E_MOLS  H2 He LiH BeH2 BH CH4 NH3 H2O HF Ne)
    set(_QM2E_TOL   1e-10)
    # Optional xcDFT He-VDZ ERI reference (auto-skipped if absent).
    set(_QM2E_XCDFT ${CMAKE_CURRENT_SOURCE_DIR}/ERI.He.VDZ.dat)

    foreach(_mol ${_QM2E_MOLS})
        add_test(NAME qm_2e_${_mol}
            COMMAND ${Python3_EXECUTABLE} ${_QM2E_PY}
                --dump $<TARGET_FILE:dump_qm_2e>
                --python ${_QM2E_WIT}
                --xyz ${_QM2E_XYZ}/${_mol}.xyz
                --basis def2-SVP
                --tol ${_QM2E_TOL}
                --xcDFT-ref ${_QM2E_XCDFT}
                --quiet)
        set_tests_properties(qm_2e_${_mol}
            PROPERTIES TIMEOUT 600 LABELS "qm_2e")
    endforeach()
endif()
# Timing only (not a ctest): old vs shell-blocked ERI, J/K and spherical transform.
add_executable(bench_qm_integrals bench_qm_integrals.cpp)
target_link_libraries(bench_qm_integrals PRIVATE curcuma_core)
target_include_directories(bench_qm_integrals PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(bench_qm_integrals PRIVATE cxx_std_17)
