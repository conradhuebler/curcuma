# Native QM SCF robustness (Nov 2026): level shift (-qm.scf_level_shift) and wavefunction stability analysis
# (-qm.scf_stability). Claude Generated.
#
# check_qm_stability.py runs dump_qm_energy (target of qm_basis_embed) and compares with reference.json, written by
# scripts/qm_stability_reference.py (PySCF 2.14: energies, dense (A + B) lowest eigenvalues, followed solutions).
#   qm_levelshift_<case>          energy (and analytic gradient) with -qm.scf_level_shift = without (1e-9 Eh, 1e-8 Eh/Bohr)
#   qm_levelshift_rescue_<case>   undamped Roothaan (scf_mode plain, scf_damping 1) fails without and converges with the
#                                 shift to the PySCF energy (1e-8 Eh)
# Run:  ctest -L qm_stability --output-on-failure
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_STAB ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_qm_stability.py --dump $<TARGET_FILE:dump_qm_energy>)
    # level shift: converged result independent of the shift (RHF, RKS, UHF; gradients = tight SCF)
    add_test(NAME qm_levelshift_h2o_hf COMMAND ${_STAB} --check shift --case h2o_hf --shift 0.5 --gradient)
    add_test(NAME qm_levelshift_h2o_hf_plain COMMAND ${_STAB} --check shift --case h2o_hf --shift 0.3
             --extra scf_mode=plain)
    add_test(NAME qm_levelshift_h2_b3lyp COMMAND ${_STAB} --check shift --case h2_r2.0_b3lyp --shift 0.5 --gradient)
    add_test(NAME qm_levelshift_ch3_uhf COMMAND ${_STAB} --check shift --case ch3_uhf --shift 0.5 --gradient)
    add_test(NAME qm_levelshift_o3_hf COMMAND ${_STAB} --check shift --case o3_hf --shift 0.3
             --extra scf_level_shift_off=1e-3)
    # level shift rescues undamped Roothaan iterations that oscillate without it
    add_test(NAME qm_levelshift_rescue_o3_hf COMMAND ${_STAB} --check shift --case o3_hf --shift 0.3
             --expect-unconverged-without --extra scf_mode=plain scf_damping=1.0 scf_max_iterations=150)
    add_test(NAME qm_levelshift_rescue_n2_r2.0_hf COMMAND ${_STAB} --check shift --case n2_r2.0_hf --shift 1.0
             --expect-unconverged-without --extra scf_mode=plain scf_damping=1.0 scf_max_iterations=150)
    set_tests_properties(qm_levelshift_h2o_hf qm_levelshift_h2o_hf_plain qm_levelshift_h2_b3lyp qm_levelshift_ch3_uhf
                         qm_levelshift_o3_hf qm_levelshift_rescue_o3_hf qm_levelshift_rescue_n2_r2.0_hf
                         PROPERTIES TIMEOUT 600 LABELS "qm_stability")
endif()
