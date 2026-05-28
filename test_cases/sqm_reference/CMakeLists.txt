# Phase 0 reference dumper for native xTB validation.
# Claude Generated (Phase 0, Apr 2026).
#
# The four test_xtb_* executables build unconditionally (header-only).
# The two dump executables require tblite.

# === Phase 3.1 overlap diagnostic ===
add_executable(test_xtb_overlap test_xtb_overlap.cpp)
target_include_directories(test_xtb_overlap PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(test_xtb_overlap PRIVATE cxx_std_17)

# === Phase 3.2 bare-Hamiltonian diagnostic ===
add_executable(test_xtb_h0 test_xtb_h0.cpp)
target_include_directories(test_xtb_h0 PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(test_xtb_h0 PRIVATE cxx_std_17)

# === Phase 3.3 isotropic Coulomb diagnostic ===
add_executable(test_xtb_coulomb test_xtb_coulomb.cpp)
target_include_directories(test_xtb_coulomb PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(test_xtb_coulomb PRIVATE cxx_std_17)

# === Phase 3.6 SCF-from-snapshot driver ===
add_executable(test_xtb_scf_snapshot test_xtb_scf_snapshot.cpp)
target_include_directories(test_xtb_scf_snapshot PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(test_xtb_scf_snapshot PRIVATE cxx_std_17)

# === Step 3a diagnostic: curcuma D4 SCF potential vs (vat_tblite − gf2.vat_extra) ===
add_executable(diag_curcuma_d4_potential diag_curcuma_d4_potential.cpp)
target_link_libraries(diag_curcuma_d4_potential PRIVATE curcuma_core)
target_include_directories(diag_curcuma_d4_potential PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(diag_curcuma_d4_potential PRIVATE cxx_std_17)

# === tblite-dependent executables (dump + multipole) ===
if(USE_TBLITE)
    # Phase 0 reference dumper
    add_executable(dump_tblite dump_tblite.cpp)
    target_include_directories(dump_tblite PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_SOURCE_DIR}/external/tblite/include
    )
    target_link_libraries(dump_tblite PRIVATE tblite gfortran)

    # Phase 3.5 multipole dumper
    add_executable(dump_tblite_multipole dump_tblite_multipole.cpp)
    target_include_directories(dump_tblite_multipole PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_SOURCE_DIR}/external/tblite/include
    )
    target_link_libraries(dump_tblite_multipole PRIVATE tblite gfortran)
    target_compile_features(dump_tblite_multipole PRIVATE cxx_std_17)

    # Molecule set for Phase 0
    set(_SQM_MOLS
        H2
        He2
        LiH
        H2O
        CH4
        NH3
        C6H6
    )

    set(_SQM_MOL_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/molecules)
    set(_SQM_DUMP_DIR  ${CMAKE_CURRENT_BINARY_DIR}/dumps)
    file(MAKE_DIRECTORY ${_SQM_DUMP_DIR})

    foreach(_mol ${_SQM_MOLS})
        foreach(_meth gfn1 gfn2)
            add_test(
                NAME sqm_reference_dump_${_mol}_${_meth}
                COMMAND dump_tblite ${_meth}
                    ${_SQM_MOL_DIR}/${_mol}.xyz
                    ${_SQM_DUMP_DIR}/${_mol}_${_meth}.json
            )
            set_tests_properties(sqm_reference_dump_${_mol}_${_meth}
                PROPERTIES TIMEOUT 60)
            add_test(
                NAME sqm_overlap_${_mol}_${_meth}
                COMMAND test_xtb_overlap ${_SQM_DUMP_DIR}/${_mol}_${_meth}.json
            )
            set_tests_properties(sqm_overlap_${_mol}_${_meth}
                PROPERTIES TIMEOUT 30
                           DEPENDS sqm_reference_dump_${_mol}_${_meth})
            add_test(
                NAME sqm_h0_${_mol}_${_meth}
                COMMAND test_xtb_h0 ${_SQM_DUMP_DIR}/${_mol}_${_meth}.json
            )
            set_tests_properties(sqm_h0_${_mol}_${_meth}
                PROPERTIES TIMEOUT 30
                           DEPENDS sqm_reference_dump_${_mol}_${_meth})
            add_test(
                NAME sqm_coulomb_${_mol}_${_meth}
                COMMAND test_xtb_coulomb ${_SQM_DUMP_DIR}/${_mol}_${_meth}.json
            )
            set_tests_properties(sqm_coulomb_${_mol}_${_meth}
                PROPERTIES TIMEOUT 30
                           DEPENDS sqm_reference_dump_${_mol}_${_meth})
            add_test(
                NAME sqm_scf_${_mol}_${_meth}
                COMMAND test_xtb_scf_snapshot ${_SQM_DUMP_DIR}/${_mol}_${_meth}.json
            )
            set_tests_properties(sqm_scf_${_mol}_${_meth}
                PROPERTIES TIMEOUT 60
                           DEPENDS sqm_reference_dump_${_mol}_${_meth})
        endforeach()
    endforeach()
endif()
