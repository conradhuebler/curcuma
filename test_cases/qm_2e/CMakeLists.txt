# Native DFT 2-electron integral validation suite (WP2).
# Claude Generated (WP2, June 2026).
#
# dump_dft_2e builds the curcuma native-DFT basis + 4-centre ERI (cartesian 6d,
# chemists' (mu nu | lam sig), McMurchie-Davidson) via DFT::cartesianERI(), and
# the Coulomb J / exchange K matrices from a dummy closed-shell density P. It
# emits them as JSON. diff_dft_2e.py orchestrates the comparison:
#   (a) element-wise ERI vs the independent Python MD witness (cartesian, same
#       AO order) -- the kernel gate, no external deps;
#   (b) internal consistency: 8-fold ERI symmetry, J/K symmetric, Tr(P*S)=2,
#       Tr(P*J)==Tr(P*K) single-orbital identity, J/K element-wise vs witness;
#   (c) optional xcDFT He-VDZ ERI cross-check (gated by --xcDFT-ref, auto-skipped
#       if the reference file is absent, so a checkout without it still runs
#       (a)+(b)).
#
# The 10 molecules and .xyz geometries are shared with dft_1e (no new inputs).
# Run:  ctest -L dft_2e --output-on-failure

add_executable(dump_dft_2e dump_dft_2e.cpp)
target_link_libraries(dump_dft_2e PRIVATE curcuma_core)
target_include_directories(dump_dft_2e PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(dump_dft_2e PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_DFT2E_PY    ${CMAKE_CURRENT_SOURCE_DIR}/diff_dft_2e.py)
    set(_DFT2E_WIT   ${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/dft_2e_python_ints.py)
    set(_DFT2E_XYZ   ${CMAKE_CURRENT_SOURCE_DIR}/../dft_1e)
    set(_DFT2E_MOLS  H2 He LiH BeH2 BH CH4 NH3 H2O HF Ne)
    set(_DFT2E_TOL   1e-10)
    # Optional xcDFT He-VDZ ERI reference (auto-skipped if absent).
    set(_DFT2E_XCDFT ${CMAKE_CURRENT_SOURCE_DIR}/ERI.He.VDZ.dat)

    foreach(_mol ${_DFT2E_MOLS})
        add_test(NAME dft_2e_${_mol}
            COMMAND ${Python3_EXECUTABLE} ${_DFT2E_PY}
                --dump $<TARGET_FILE:dump_dft_2e>
                --python ${_DFT2E_WIT}
                --xyz ${_DFT2E_XYZ}/${_mol}.xyz
                --basis def2-SVP
                --tol ${_DFT2E_TOL}
                --xcDFT-ref ${_DFT2E_XCDFT}
                --quiet)
        set_tests_properties(dft_2e_${_mol}
            PROPERTIES TIMEOUT 600 LABELS "dft_2e")
    endforeach()
endif()