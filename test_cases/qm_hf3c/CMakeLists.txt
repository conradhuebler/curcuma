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
