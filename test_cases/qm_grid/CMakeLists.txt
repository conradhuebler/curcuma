# Native DFT quadrature grid (WP4): Lebedev order, grid vs PySCF, electron count. Claude Generated (Oct 2026).
#
# test_qm_grid (C++) prints JSON; check_qm_grid.py compares it with reference.json, written by
# scripts/qm_grid_pyscf_reference.py (PySCF 2.14, venv ~/.venvs/qmref; the generating settings are in the file).
#   qm_grid_lebedev           every Lebedev grid integrates all real Y_lm with l <= degree to 1e-13 and is NOT exact
#                             at l = degree + 1; grids <= 590 points: orthonormality of Y_lm, l <= degree/2, to 1e-13
#   qm_grid_pts_<case>        points and weights vs PySCF gen_grid (same radial scheme, counts, partition, size
#                             adjustment, no pruning): every 211th point: coordinates and partition factor w/w_quad to 1e-12, totals and test integrals 1e-12 relative
#   qm_grid_rho_<case>        HF/def2-SVP density: integral of rho on the grid vs PySCF (1e-8) and vs N_el
#                             (reported), rho and grad rho at sampled points vs PySCF (1e-7 of the local scale)
#   qm_grid_aofd_<mol>        AO gradients vs 5-point finite differences of the AO values (1e-6 relative)
# Geometries: the xyz files of qm_1e and qm_hf3c (no new structure files).
# Run:  ctest -L qm_grid --output-on-failure

add_executable(test_qm_grid test_qm_grid.cpp)
target_link_libraries(test_qm_grid PRIVATE curcuma_core)
target_include_directories(test_qm_grid PRIVATE ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/external)
target_compile_features(test_qm_grid PRIVATE cxx_std_17)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_GRID_CHECK ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_grid.py
        --exe $<TARGET_FILE:test_qm_grid> --ref ${CMAKE_CURRENT_SOURCE_DIR}/reference.json)
    add_test(NAME qm_grid_lebedev COMMAND ${_GRID_CHECK} --check lebedev)
    set_tests_properties(qm_grid_lebedev PROPERTIES TIMEOUT 300 LABELS "qm_grid")
    foreach(_case H2O_sg1 H2O_sg3 H2O_stratmann H2O_adjust_becke H2O_adjust_none H2O_euler CH4_sg1 benzene_sg1
                  LiH_sg2 Ne_sg0)
        add_test(NAME qm_grid_pts_${_case} COMMAND ${_GRID_CHECK} --check points --case ${_case})
        set_tests_properties(qm_grid_pts_${_case} PROPERTIES TIMEOUT 300 LABELS "qm_grid")
    endforeach()
    foreach(_case H2O_sg1 CH4_sg1 benzene_sg1 LiH_sg2 Ne_sg0)
        add_test(NAME qm_grid_rho_${_case} COMMAND ${_GRID_CHECK} --check density --case ${_case})
        set_tests_properties(qm_grid_rho_${_case} PROPERTIES TIMEOUT 600 LABELS "qm_grid")
    endforeach()
    foreach(_mol H2O CH4)
        add_test(NAME qm_grid_aofd_${_mol} COMMAND ${_GRID_CHECK} --check aofd --case ${_mol}_sg1)
        set_tests_properties(qm_grid_aofd_${_mol} PROPERTIES TIMEOUT 300 LABELS "qm_grid")
    endforeach()
    add_test(NAME qm_grid_aofd_H2O_cartesian COMMAND ${_GRID_CHECK} --check aofd --case H2O_sg1 --cartesian)
    set_tests_properties(qm_grid_aofd_H2O_cartesian PROPERTIES TIMEOUT 300 LABELS "qm_grid")
endif()
