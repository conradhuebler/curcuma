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
# Reads tblite dumps from ${D4_DIAG_DUMP_DIR} (default: <source>/release_tblite/dumps).
# Dumps must be generated separately via dump_tblite_multipole (USE_TBLITE build).
add_executable(diag_curcuma_d4_potential diag_curcuma_d4_potential.cpp)
target_link_libraries(diag_curcuma_d4_potential PRIVATE curcuma_core)
target_include_directories(diag_curcuma_d4_potential PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(diag_curcuma_d4_potential PRIVATE cxx_std_17)

# === GFN2 component audit (Claude Generated): per-container energy diff vs tblite ===
# Reads tblite dumps containing the e_components_tblite block (added with the
# AP-component-audit tblite patch). Injects tblite's converged density into
# curcuma's XTB and prints per-container |E_curcuma - E_tblite|.
add_executable(diag_curcuma_energy_components diag_curcuma_energy_components.cpp)
target_link_libraries(diag_curcuma_energy_components PRIVATE curcuma_core)
target_include_directories(diag_curcuma_energy_components PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(diag_curcuma_energy_components PRIVATE cxx_std_17)

# === GFN1 D3 C6 probe (Claude Generated, WP2): curcuma interpolated D3 C6 ===
# Dumps curcuma's per-pair interpolated GFN1 D3 C6 for diffing against s-dftd3
# (dump_dftd3_atomic_c6). Localises the GFN1 dispersion residual to the C6 table.
add_executable(diag_d3_c6 diag_d3_c6.cpp)
target_link_libraries(diag_d3_c6 PRIVATE curcuma_core)
target_include_directories(diag_d3_c6 PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(diag_d3_c6 PRIVATE cxx_std_17)

# === GFN2-D4 reference-C6 audit (Claude Generated): curcuma cache vs dftd4 model%c6 ===
# Reads a c6 dump from dump_dftd4_c6 (Fortran, dftd4 new_d4_model ref=gfn2) and
# diffs against curcuma's D4ParameterGenerator::getC6FlatCache() per
# (Z_i, Z_j, iref, jref). Localises whether the C-path residual lives in the
# constructed reference C6 itself. C++-only (reads a JSON file); the Fortran
# dumper is built in the USE_TBLITE block below.
add_executable(diag_curcuma_d4_c6 diag_curcuma_d4_c6.cpp)
target_link_libraries(diag_curcuma_d4_c6 PRIVATE curcuma_core)
target_include_directories(diag_curcuma_d4_c6 PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(diag_curcuma_d4_c6 PRIVATE cxx_std_17)

# === GFN2-D4 WEIGHTED per-atom C6 audit (Claude Generated) ===
# Diffs curcuma's weightedC6Gfn2 (the per-atom C6 entering the energy) against
# tblite's get_atomic_c6 (from dump_dftd4_atomic_c6) at tblite's converged
# charges, plus a D4-CN comparison. Localises the residual to the CN/charge
# weighting (gwvec) vs ruling D4 out entirely. C++-only; Fortran dumper below.
add_executable(diag_curcuma_atomic_c6 diag_curcuma_atomic_c6.cpp)
target_link_libraries(diag_curcuma_atomic_c6 PRIVATE curcuma_core)
target_include_directories(diag_curcuma_atomic_c6 PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(diag_curcuma_atomic_c6 PRIVATE cxx_std_17)

# === Multi-step SCC extrapolation (Claude Generated, June 2026) ===========
# Self-contained (reads committed xyz, no tblite dumps): drives the native
# GFN1/GFN2 SCF along a smooth interpolated trajectory and checks that the
# opt-in aspc/gauss SCC extrapolation (a) lands on the same fixpoint as the
# 1-step warm-start and (b) does not cost more SCF iterations. Run via:
#   ctest -L scf_extrapolation
add_executable(test_xtb_scf_extrapolation test_xtb_scf_extrapolation.cpp)
target_link_libraries(test_xtb_scf_extrapolation PRIVATE curcuma_core)
target_include_directories(test_xtb_scf_extrapolation PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/external
)
target_compile_features(test_xtb_scf_extrapolation PRIVATE cxx_std_17)
add_test(
    NAME scf_extrapolation_caffeine
    COMMAND test_xtb_scf_extrapolation
        ${CMAKE_CURRENT_SOURCE_DIR}/molecules/caffeine.xyz
        ${CMAKE_CURRENT_SOURCE_DIR}/molecules/caffeine.opt.xyz
)
set_tests_properties(scf_extrapolation_caffeine PROPERTIES TIMEOUT 120 LABELS "scf_extrapolation")

# Per-molecule ctest entries — track the GFN2 D4 SCF potential vs tblite. Each
# tolerance is set just above the current state so the suite is green now, and
# tightens as the D4 model gets fixed. Run via:  ctest -R "d4_diag"
set(D4_DIAG_DUMP_DIR "${CMAKE_SOURCE_DIR}/release_tblite/dumps"
    CACHE PATH "Directory containing tblite GFN2 dumps (*_gfn2.json) for diag_d4 ctest")
set(_D4_DIAG_MOLS_TOL
    "H2:1.5e-6"      # symmetric, refh=0 → exact-ish
    "H2O:1.0e-6"     # sub-µ after CN+α-zeta
    "NH3:1.5e-6"     # sub-µ after CN+α-zeta
    "CH4:1.5e-5"     # C-path overcorrection — to be fixed
    "triose:1.5e-4"  # 66-atom sum of C-path residual — to be fixed
)
foreach(_mt ${_D4_DIAG_MOLS_TOL})
    string(REPLACE ":" ";" _mt_list ${_mt})
    list(GET _mt_list 0 _mol)
    list(GET _mt_list 1 _tol)
    add_test(
        NAME d4_diag_${_mol}
        COMMAND diag_curcuma_d4_potential --quiet --tol ${_tol}
            ${D4_DIAG_DUMP_DIR}/${_mol}_gfn2.json
            ${CMAKE_CURRENT_SOURCE_DIR}/molecules/${_mol}.xyz
    )
    set_tests_properties(d4_diag_${_mol} PROPERTIES TIMEOUT 60 LABELS "d4_diag")
endforeach()

# Per-molecule ctest entries for the GFN2 component audit. Reads the same
# dumps as d4_diag (must contain the new e_components_tblite block from the
# refreshed dump_tblite_multipole). Tolerances are set just above the current
# combined |disp+elec| residual so the suite is green now and tightens as
# the D4-C-path / triose residual gets fixed. Per-container max|diff| is NOT
# the gating number — see the tool's header for why (tblite splits D4 across
# dispersion/eelec; curcuma puts it all in m_E_dispersion). Run via:
# ctest -L gfn2_align
set(_ECOMP_MOLS_TOL
    "H2:5.0e-10"     # bit-identical-ish
    "H2O:1.0e-8"     # 2.7e-10
    "NH3:1.0e-8"     # 1.4e-10
    "CH4:1.0e-8"     # 2.4e-11 (C-path resolved by refcovcn 2-body fix)
    "triose:1.0e-6"  # 5.9e-9 (3.2 mEh -> 5.9 nEh after refcovcn 2-body + ATM 3-body)
)
foreach(_mt ${_ECOMP_MOLS_TOL})
    string(REPLACE ":" ";" _mt_list ${_mt})
    list(GET _mt_list 0 _mol)
    list(GET _mt_list 1 _tol)
    add_test(
        NAME ecomp_${_mol}
        COMMAND diag_curcuma_energy_components --quiet --tol ${_tol}
            ${D4_DIAG_DUMP_DIR}/${_mol}_gfn2.json
            ${CMAKE_CURRENT_SOURCE_DIR}/molecules/${_mol}.xyz
    )
    set_tests_properties(ecomp_${_mol} PROPERTIES TIMEOUT 60 LABELS "gfn2_align")
endforeach()

# Per-molecule ctest entries for the GFN1 component audit (label gfn1_align).
# Same injected-density per-container diff as ecomp_* above, but reads the
# *_gfn1.json forensic dumps (dump_tblite_multipole gfn1 ...). Repulsion is
# bit-clean. Three WP2 fixes brought the injected-density combined |disp+elec|
# from ~1e-2 (triose) down to <1e-8 on every molecule:
#   2026-05-30  D3 reference-CN table regenerated (was scrambled; C6 wrong).
#   2026-05-31  D3 C8/C6 exact s-dftd3 form (was empirical 10.72*r4r2, +0.06%/pair).
#   2026-05-31  D3 CN real-space cutoff 25 Bohr to match tblite (tblite/disp/d3.f90;
#               curcuma had summed all pairs -> CN too high on large molecules).
# The worst container is now "electronic" (GFN1 electronic floor), not dispersion.
# Tolerances are set just above the current values. Like the GFN2 align/d4_diag
# tests these read generated (not committed) dumps from ${D4_DIAG_DUMP_DIR}; a
# clean checkout must run dump_tblite_multipole first (incl. the 231-atom complex,
# ~33 MB; its diag runs in <1 s). Run: ctest -L gfn1_align
set(_ECOMP_MOLS_TOL_GFN1
    "H2:1.0e-9"      # 7.8e-11 (electronic)
    "H2O:1.0e-9"     # 2.5e-10
    "NH3:1.0e-9"     # 2.1e-10
    "CH4:1.0e-9"     # 2.2e-10
    "triose:1.0e-8"  # 5.3e-9
    "complex:1.0e-8" # 7.4e-9 (231 atoms; worst term electronic, dispersion now exact)
)
foreach(_mt ${_ECOMP_MOLS_TOL_GFN1})
    string(REPLACE ":" ";" _mt_list ${_mt})
    list(GET _mt_list 0 _mol)
    list(GET _mt_list 1 _tol)
    add_test(
        NAME ecomp_${_mol}_gfn1
        COMMAND diag_curcuma_energy_components --quiet --tol ${_tol}
            ${D4_DIAG_DUMP_DIR}/${_mol}_gfn1.json
            ${CMAKE_CURRENT_SOURCE_DIR}/molecules/${_mol}.xyz
    )
    set_tests_properties(ecomp_${_mol}_gfn1 PROPERTIES TIMEOUT 60 LABELS "gfn1_align")
endforeach()

# === Native GFN1/GFN2 vs committed tblite references (Claude Generated, 2026-05) ===
# End-to-end validation: runs the real `curcuma` binary (native ngfn1/ngfn2) and
# compares the single-point energy against the committed tblite reference JSONs in
# reference_data/ (produced by dump_tblite_reference in a USE_TBLITE build). Needs
# NO tblite build -- the references are committed (like the GFN-FF suite).
#
# HONEST GATE: every test gates total energy at 1e-8 Eh (the target agreement with
# tblite). Molecules that do NOT yet reach 1e-8 are marked WILL_FAIL (ctest
# expected-failure) with their MEASURED residual in the list below -- the gate is
# never loosened to mask a gap. When a method reaches 1e-8 for an xfail molecule,
# ctest reports it as a FAILURE ("passed but expected to fail"), forcing removal of
# the xfail: progress is self-documenting and a closed gap cannot silently re-open.
# Baseline measured 2026-05-30 (release/, USE_TBLITE=OFF). See docs/SQM_VALIDATION.md
# and docs/SQM_VALIDATION_ROADMAP.md. Run:  ctest -L gfn2_validation / gfn1_validation
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    set(_SQM_VAL_PY   ${CMAKE_CURRENT_SOURCE_DIR}/validate_sqm.py)
    set(_SQM_VAL_REF  ${CMAKE_CURRENT_SOURCE_DIR}/reference_data)
    set(_SQM_VAL_MOL  ${CMAKE_CURRENT_SOURCE_DIR}/molecules)
    set(_SQM_VAL_TOL  1e-8)
    set(_SQM_MOLS H2 He2 LiH H2O CH4 NH3 C6H6 HCN acetic_acid_dimer caffeine triose complex)

    # Molecules NOT yet at 1e-8 -> tracked expected-failures (WILL_FAIL). Measured
    # |dE| = |E_native - E_tblite| (Eh), 2026-05-31 (after the GFN1 third-order
    # double-count fix, the D3 C8/C6 exact-form fix, and the D3 CN 25-Bohr cutoff
    # fix that matches tblite):
    #   GFN2: only complex (6.95e-5). All others <= 1e-8 (bit-identical).
    #   GFN1: 10/12 now <= 1e-8. Two remain at the same ~1.5e-8 numerical floor:
    #         He2 (dispersion-only) and complex (231 atoms). complex was 2.6e-7
    #         until the D3 CN cutoff fix removed the dispersion gap (now exact,
    #         <1e-9 at fixed density); its residual is now the GFN1 electronic
    #         floor (~7e-9 at fixed density) plus SCF density relaxation. See
    #         docs/SQM_WP2_gfn1_accuracy.md.
    set(_GFN2_XFAIL complex)
    set(_GFN1_XFAIL He2 complex)

    foreach(_meth gfn1 gfn2)
        if(_meth STREQUAL gfn1)
            set(_xfail ${_GFN1_XFAIL})
        else()
            set(_xfail ${_GFN2_XFAIL})
        endif()
        foreach(_mol ${_SQM_MOLS})
            add_test(NAME sqm_val_${_mol}_${_meth}
                COMMAND ${Python3_EXECUTABLE} ${_SQM_VAL_PY}
                    --curcuma $<TARGET_FILE:curcuma> --quiet --tol-energy ${_SQM_VAL_TOL}
                    --ref ${_SQM_VAL_REF}/${_mol}_${_meth}.ref.json
                    --xyz ${_SQM_VAL_MOL}/${_mol}.xyz)
            set_tests_properties(sqm_val_${_mol}_${_meth}
                PROPERTIES TIMEOUT 120 LABELS "${_meth}_validation")
            if(${_mol} IN_LIST _xfail)
                set_tests_properties(sqm_val_${_mol}_${_meth} PROPERTIES WILL_FAIL TRUE)
            endif()
        endforeach()
    endforeach()

    # ── Native ALPB solvation validation (Claude Generated, June 2026) ──
    # Compares native GFN1/GFN2 + ALPB total energy vs the committed tblite ALPB
    # references (reference_data/<mol>_<meth>_<solvent>.ref.json, generated by
    # dump_tblite_reference --solvent). validate_sqm.py reads the "solvent" field
    # from the ref and passes -xtb.solvent/-xtb.solvent_model 3. Same 1e-8 gate;
    # all currently pass (worst ~6e-9 Eh). Labels: gfn{1,2}_solvation.
    # See docs/SQM_SOLVATION_WP.md (WP1.4).
    set(_SQM_SOLV_MOLS H2O NH3 CH4 HCN C6H6 acetic_acid_dimer caffeine)
    set(_SQM_SOLV_SOLVENTS water dmso acetone chcl3)
    foreach(_meth gfn1 gfn2)
        foreach(_mol ${_SQM_SOLV_MOLS})
            foreach(_solv ${_SQM_SOLV_SOLVENTS})
                set(_ref ${_SQM_VAL_REF}/${_mol}_${_meth}_${_solv}.ref.json)
                if(EXISTS ${_ref})
                    add_test(NAME sqm_solv_${_mol}_${_meth}_${_solv}
                        COMMAND ${Python3_EXECUTABLE} ${_SQM_VAL_PY}
                            --curcuma $<TARGET_FILE:curcuma> --quiet --tol-energy ${_SQM_VAL_TOL}
                            --ref ${_ref}
                            --xyz ${_SQM_VAL_MOL}/${_mol}.xyz)
                    set_tests_properties(sqm_solv_${_mol}_${_meth}_${_solv}
                        PROPERTIES TIMEOUT 120 LABELS "${_meth}_solvation")
                endif()
            endforeach()
        endforeach()
    endforeach()

    # GPU native xTB validation (Claude Generated, GPU port): identical 1e-8 gate
    # vs the committed tblite references, but with -gpu cuda so the eigensolve (and,
    # in later stages, the whole SCF/gradient) runs on the device. The GPU energies
    # match the CPU path to <1e-8, so the SAME xfail set applies. Without a CUDA
    # device the wrapper falls back to CPU, so these never false-fail off-GPU.
    # Labels: gpu_gfn1_validation / gpu_gfn2_validation. Only built with USE_CUDA_XTB.
    if(USE_CUDA_XTB)
        foreach(_meth gfn1 gfn2)
            if(_meth STREQUAL gfn1)
                set(_xfail ${_GFN1_XFAIL})
            else()
                set(_xfail ${_GFN2_XFAIL})
            endif()
            foreach(_mol ${_SQM_MOLS})
                add_test(NAME sqm_val_gpu_${_mol}_${_meth}
                    COMMAND ${Python3_EXECUTABLE} ${_SQM_VAL_PY}
                        --curcuma $<TARGET_FILE:curcuma> --quiet --tol-energy ${_SQM_VAL_TOL}
                        --gpu cuda
                        --ref ${_SQM_VAL_REF}/${_mol}_${_meth}.ref.json
                        --xyz ${_SQM_VAL_MOL}/${_mol}.xyz)
                set_tests_properties(sqm_val_gpu_${_mol}_${_meth}
                    PROPERTIES TIMEOUT 120 LABELS "gpu_${_meth}_validation")
                if(${_mol} IN_LIST _xfail)
                    set_tests_properties(sqm_val_gpu_${_mol}_${_meth} PROPERTIES WILL_FAIL TRUE)
                endif()
            endforeach()
        endforeach()
    endif()
endif()

# === Stage 3 GPU integral component tests (Claude Generated, GPU port) ======
# Each compares a device integral kernel to its CPU reference elementwise @1e-9.
# Built only with USE_CUDA_XTB; the exes runtime-skip (exit 0) without a device,
# so they never false-fail off-GPU. Linked against curcuma_core (CPU reference
# + the compiled device kernels).
if(USE_CUDA_XTB)
    set(_SQM_CUDA_MOL ${CMAKE_CURRENT_SOURCE_DIR}/molecules)
    set(_SQM_CUDA_COMP_MOLS H2O CH4 NH3 C6H6 HCN triose)

    # 3a: coordination numbers + CN-shifted self-energies.
    add_executable(test_xtb_cuda_cn test_xtb_cuda_cn.cpp)
    target_link_libraries(test_xtb_cuda_cn PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_cn PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_cn PRIVATE cxx_std_17)
    foreach(_meth gfn1 gfn2)
        foreach(_mol ${_SQM_CUDA_COMP_MOLS})
            add_test(NAME sqm_cuda_cn_${_mol}_${_meth}
                COMMAND test_xtb_cuda_cn --quiet --tol 1e-9 ${_meth}
                        ${_SQM_CUDA_MOL}/${_mol}.xyz)
            set_tests_properties(sqm_cuda_cn_${_mol}_${_meth}
                PROPERTIES TIMEOUT 120 LABELS "gpu_integrals")
        endforeach()
    endforeach()

    # 3b: overlap S + bare Hamiltonian H0 + device Cholesky L.
    add_executable(test_xtb_cuda_h0 test_xtb_cuda_h0.cpp)
    target_link_libraries(test_xtb_cuda_h0 PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_h0 PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_h0 PRIVATE cxx_std_17)
    foreach(_meth gfn1 gfn2)
        foreach(_mol ${_SQM_CUDA_COMP_MOLS})
            add_test(NAME sqm_cuda_h0_${_mol}_${_meth}
                COMMAND test_xtb_cuda_h0 --quiet --tol 1e-9 ${_meth}
                        ${_SQM_CUDA_MOL}/${_mol}.xyz)
            set_tests_properties(sqm_cuda_h0_${_mol}_${_meth}
                PROPERTIES TIMEOUT 120 LABELS "gpu_integrals")
        endforeach()
    endforeach()

    # 3c: Coulomb gamma matrix.
    add_executable(test_xtb_cuda_gamma test_xtb_cuda_gamma.cpp)
    target_link_libraries(test_xtb_cuda_gamma PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_gamma PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_gamma PRIVATE cxx_std_17)
    foreach(_meth gfn1 gfn2)
        foreach(_mol ${_SQM_CUDA_COMP_MOLS})
            add_test(NAME sqm_cuda_gamma_${_mol}_${_meth}
                COMMAND test_xtb_cuda_gamma --quiet --tol 1e-9 ${_meth}
                        ${_SQM_CUDA_MOL}/${_mol}.xyz)
            set_tests_properties(sqm_cuda_gamma_${_mol}_${_meth}
                PROPERTIES TIMEOUT 120 LABELS "gpu_integrals")
        endforeach()
    endforeach()

    # 3d: GFN2 multipole integrals (dipole + traceless quadrupole). GFN2 only.
    add_executable(test_xtb_cuda_multipole test_xtb_cuda_multipole.cpp)
    target_link_libraries(test_xtb_cuda_multipole PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_multipole PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_multipole PRIVATE cxx_std_17)
    foreach(_mol ${_SQM_CUDA_COMP_MOLS})
        add_test(NAME sqm_cuda_multipole_${_mol}
            COMMAND test_xtb_cuda_multipole --quiet --tol 1e-9
                    ${_SQM_CUDA_MOL}/${_mol}.xyz)
        set_tests_properties(sqm_cuda_multipole_${_mol}
            PROPERTIES TIMEOUT 120 LABELS "gpu_integrals")
    endforeach()

    # 4 (primitive): overlap derivative dS/dR — crux of the nuclear gradient.
    add_executable(test_xtb_cuda_overlap_grad test_xtb_cuda_overlap_grad.cpp)
    target_link_libraries(test_xtb_cuda_overlap_grad PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_overlap_grad PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_overlap_grad PRIVATE cxx_std_17)
    foreach(_meth gfn1 gfn2)
        foreach(_mol ${_SQM_CUDA_COMP_MOLS})
            add_test(NAME sqm_cuda_overlapgrad_${_mol}_${_meth}
                COMMAND test_xtb_cuda_overlap_grad --quiet --tol 1e-9 ${_meth}
                        ${_SQM_CUDA_MOL}/${_mol}.xyz)
            set_tests_properties(sqm_cuda_overlapgrad_${_mol}_${_meth}
                PROPERTIES TIMEOUT 120 LABELS "gpu_gradient")
        endforeach()
    endforeach()

    # 4a: full GFN1 nuclear gradient (device sections 1/2/3 + host 3b/4) vs CPU.
    add_executable(test_xtb_cuda_gradient test_xtb_cuda_gradient.cpp)
    target_link_libraries(test_xtb_cuda_gradient PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_gradient PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_gradient PRIVATE cxx_std_17)
    foreach(_meth gfn1 gfn2)
        foreach(_mol ${_SQM_CUDA_COMP_MOLS})
            add_test(NAME sqm_cuda_gradient_${_mol}_${_meth}
                COMMAND test_xtb_cuda_gradient --quiet --tol 1e-7 ${_meth}
                        ${_SQM_CUDA_MOL}/${_mol}.xyz)
            set_tests_properties(sqm_cuda_gradient_${_mol}_${_meth}
                PROPERTIES TIMEOUT 180 LABELS "gpu_gradient")
        endforeach()
    endforeach()

    # 5 (Part A): single-shot D4 EEQ charges + dq/dx charge-response on the device
    # vs D4ChargeModel CPU. Charges FP64 LU (~1e-10); response has reductions → 1e-9.
    add_executable(test_xtb_cuda_eeq test_xtb_cuda_eeq.cpp)
    target_link_libraries(test_xtb_cuda_eeq PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_eeq PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_eeq PRIVATE cxx_std_17)
    foreach(_meth gfn1 gfn2)
        foreach(_mol ${_SQM_CUDA_COMP_MOLS})
            add_test(NAME sqm_cuda_eeq_${_mol}_${_meth}
                COMMAND test_xtb_cuda_eeq --quiet --tol 1e-9 ${_meth}
                        ${_SQM_CUDA_MOL}/${_mol}.xyz)
            set_tests_properties(sqm_cuda_eeq_${_mol}_${_meth}
                PROPERTIES TIMEOUT 120 LABELS "gpu_integrals")
        endforeach()
    endforeach()

    # 5 (Part B1): device atomic Mulliken charges from the resident density vs the
    # host updatePopulations (GFN2 only — the resident ao2at map + in-SCF D4 q_at).
    add_executable(test_xtb_cuda_qat test_xtb_cuda_qat.cpp)
    target_link_libraries(test_xtb_cuda_qat PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_qat PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_qat PRIVATE cxx_std_17)
    foreach(_mol ${_SQM_CUDA_COMP_MOLS})
        add_test(NAME sqm_cuda_qat_${_mol}
            COMMAND test_xtb_cuda_qat --quiet --tol 1e-9 gfn2
                    ${_SQM_CUDA_MOL}/${_mol}.xyz)
        set_tests_properties(sqm_cuda_qat_${_mol}
            PROPERTIES TIMEOUT 120 LABELS "gpu_integrals")
    endforeach()

    # 5 (Part B2): in-SCF GFN2 D4 atom-potential dE_D4/dq on the device vs the host
    # D4Evaluator at a frozen charge set (GFN2 only).
    add_executable(test_xtb_cuda_d4pot test_xtb_cuda_d4pot.cpp)
    target_link_libraries(test_xtb_cuda_d4pot PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_d4pot PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_d4pot PRIVATE cxx_std_17)
    foreach(_mol ${_SQM_CUDA_COMP_MOLS})
        add_test(NAME sqm_cuda_d4pot_${_mol}
            COMMAND test_xtb_cuda_d4pot --quiet --tol 1e-9 gfn2
                    ${_SQM_CUDA_MOL}/${_mol}.xyz)
        set_tests_properties(sqm_cuda_d4pot_${_mol}
            PROPERTIES TIMEOUT 120 LABELS "gpu_integrals")
    endforeach()

    # === Stage 6 device-resident SCF component tests (Claude Generated) ======
    # The fused per-iteration loop body broken into device-vs-host kernel checks.
    # All runtime-skip (exit 0) without a CUDA device. Label: gpu_scf.

    # 6 (S6.1): device occupations (Fermi / integer fill) vs host occupationsFromEps
    # at a converged eps spectrum. Method-independent kernel — gfn1 + gfn2.
    add_executable(test_xtb_cuda_occ test_xtb_cuda_occ.cpp)
    target_link_libraries(test_xtb_cuda_occ PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_occ PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_occ PRIVATE cxx_std_17)
    foreach(_meth gfn1 gfn2)
        foreach(_mol ${_SQM_CUDA_COMP_MOLS})
            add_test(NAME sqm_cuda_occ_${_mol}_${_meth}
                COMMAND test_xtb_cuda_occ --quiet --tol 1e-9 ${_meth}
                        ${_SQM_CUDA_MOL}/${_mol}.xyz)
            set_tests_properties(sqm_cuda_occ_${_mol}_${_meth}
                PROPERTIES TIMEOUT 120 LABELS "gpu_scf")
        endforeach()
    endforeach()

    # 6 (S6.2): device resident shell Mulliken charges vs host getShellCharges (gfn2).
    add_executable(test_xtb_cuda_qsh test_xtb_cuda_qsh.cpp)
    target_link_libraries(test_xtb_cuda_qsh PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_qsh PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_qsh PRIVATE cxx_std_17)
    foreach(_mol ${_SQM_CUDA_COMP_MOLS})
        add_test(NAME sqm_cuda_qsh_${_mol}
            COMMAND test_xtb_cuda_qsh --quiet --tol 1e-9 gfn2
                    ${_SQM_CUDA_MOL}/${_mol}.xyz)
        set_tests_properties(sqm_cuda_qsh_${_mol}
            PROPERTIES TIMEOUT 120 LABELS "gpu_scf")
    endforeach()

    # 6 (S6.2b): device D4 reference weights W/dWq rebuilt from q vs host buildRefWFlat
    # at a frozen charge set (gfn2).
    add_executable(test_xtb_cuda_d4refw test_xtb_cuda_d4refw.cpp)
    target_link_libraries(test_xtb_cuda_d4refw PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_d4refw PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_d4refw PRIVATE cxx_std_17)
    foreach(_mol ${_SQM_CUDA_COMP_MOLS})
        add_test(NAME sqm_cuda_d4refw_${_mol}
            COMMAND test_xtb_cuda_d4refw --quiet --tol 1e-9 gfn2
                    ${_SQM_CUDA_MOL}/${_mol}.xyz)
        set_tests_properties(sqm_cuda_d4refw_${_mol}
            PROPERTIES TIMEOUT 120 LABELS "gpu_scf")
    endforeach()

    # 6 (S6.3): device SCC energy (coulomb / third-order / multipole) vs the host
    # evaluateComponentsAtFixedDensity at the identical resident state (gfn2).
    add_executable(test_xtb_cuda_scc_energy test_xtb_cuda_scc_energy.cpp)
    target_link_libraries(test_xtb_cuda_scc_energy PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_scc_energy PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_scc_energy PRIVATE cxx_std_17)
    foreach(_mol ${_SQM_CUDA_COMP_MOLS})
        add_test(NAME sqm_cuda_scc_energy_${_mol}
            COMMAND test_xtb_cuda_scc_energy --quiet --tol 1e-9 gfn2
                    ${_SQM_CUDA_MOL}/${_mol}.xyz)
        set_tests_properties(sqm_cuda_scc_energy_${_mol}
            PROPERTIES TIMEOUT 180 LABELS "gpu_scf")
    endforeach()

    # 6 (S6.4): device modified-Broyden mixer vs host BroydenMixer over a synthetic
    # contractive sequence (standalone numeric, no molecule).
    add_executable(test_xtb_cuda_broyden test_xtb_cuda_broyden.cpp)
    target_link_libraries(test_xtb_cuda_broyden PRIVATE curcuma_core)
    target_include_directories(test_xtb_cuda_broyden PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(test_xtb_cuda_broyden PRIVATE cxx_std_17)
    add_test(NAME sqm_cuda_broyden
        COMMAND test_xtb_cuda_broyden --quiet --tol 1e-10)
    set_tests_properties(sqm_cuda_broyden
        PROPERTIES TIMEOUT 120 LABELS "gpu_scf")
endif()

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

    # Phase 3.5 multipole dumper. Needs a patched tblite exposing the extended
    # result C-API (atomic multipoles, multipole integrals, per-container energy
    # components). Default OFF so `make all` builds against the bundled upstream
    # tblite (a stub main); set -DTBLITE_EXT_RESULT_API=ON with a patched tblite
    # to build the real dumper. See dump_tblite_multipole.cpp header.
    option(TBLITE_EXT_RESULT_API
        "Build dump_tblite_multipole against a patched tblite with the extended result C-API" OFF)
    add_executable(dump_tblite_multipole dump_tblite_multipole.cpp)
    target_include_directories(dump_tblite_multipole PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_SOURCE_DIR}/external/tblite/include
    )
    target_link_libraries(dump_tblite_multipole PRIVATE tblite gfortran)
    target_compile_features(dump_tblite_multipole PRIVATE cxx_std_17)
    if(TBLITE_EXT_RESULT_API)
        target_compile_definitions(dump_tblite_multipole PRIVATE TBLITE_EXT_RESULT_API)
    endif()

    # Compact reference dumper for the committed SQM validation suite.
    # Emits only validation-relevant fields (small, git-friendly *.ref.json).
    add_executable(dump_tblite_reference dump_tblite_reference.cpp)
    target_include_directories(dump_tblite_reference PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/external
        ${CMAKE_SOURCE_DIR}/external/tblite/include
    )
    target_link_libraries(dump_tblite_reference PRIVATE tblite gfortran)
    target_compile_features(dump_tblite_reference PRIVATE cxx_std_17)

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

    # === GFN2-D4 reference-C6 dumper (Claude Generated) ===
    # Standalone Fortran exe: dftd4 new_d4_model(ref=gfn2) -> model%c6 JSON.
    # Links dftd4 + mctc-lib directly. tblite's FetchContent compiled Fortran
    # but only enabled the language in its own subtree; enable it here so this
    # target gets Fortran compile/link rules. Reuses the already-detected
    # CMAKE_Fortran_COMPILER. Module (.mod) dirs come from the build trees.
    enable_language(Fortran)
    add_executable(dump_dftd4_c6 dump_dftd4_c6.f90)
    target_include_directories(dump_dftd4_c6 PRIVATE
        ${CMAKE_BINARY_DIR}/_deps/dftd4-build/include
        ${CMAKE_BINARY_DIR}/_deps/mctc-lib-build/include
    )
    target_link_libraries(dump_dftd4_c6 PRIVATE dftd4 mctc-lib multicharge)

    # Weighted per-atom C6 dumper (tblite get_atomic_c6 at given charges).
    add_executable(dump_dftd4_atomic_c6 dump_dftd4_atomic_c6.f90)
    target_include_directories(dump_dftd4_atomic_c6 PRIVATE
        ${CMAKE_BINARY_DIR}/_deps/dftd4-build/include
        ${CMAKE_BINARY_DIR}/_deps/mctc-lib-build/include
    )
    target_link_libraries(dump_dftd4_atomic_c6 PRIVATE dftd4 mctc-lib multicharge)

    # ATM three-body energy probe (tblite get_dispersion3 at q=0 C6).
    add_executable(dump_dftd4_atm dump_dftd4_atm.f90)
    target_include_directories(dump_dftd4_atm PRIVATE
        ${CMAKE_BINARY_DIR}/_deps/dftd4-build/include
        ${CMAKE_BINARY_DIR}/_deps/mctc-lib-build/include
    )
    target_link_libraries(dump_dftd4_atm PRIVATE dftd4 mctc-lib multicharge)

    # === simple-dftd3 atomic-C6 dumper (Claude Generated, WP2) ===
    # Dumps s-dftd3's weighted per-atom C6 (the exact term GFN1 D3 uses, since
    # tblite GFN1 dispersion IS s-dftd3 with s9=0). Diffed against curcuma's
    # interpolated D3 C6 to localise the GFN1 dispersion residual to the C6 table.
    add_executable(dump_dftd3_atomic_c6 dump_dftd3_atomic_c6.f90)
    target_include_directories(dump_dftd3_atomic_c6 PRIVATE
        ${CMAKE_BINARY_DIR}/_deps/s-dftd3-build/include
        ${CMAKE_BINARY_DIR}/_deps/mctc-lib-build/include
    )
    target_link_libraries(dump_dftd3_atomic_c6 PRIVATE s-dftd3 mctc-lib)

    # c6 audit: dump dftd4 model%c6, then diff against curcuma's cache.
    # Element coverage: H2O (H,O), CH4 (C,H), NH3 (N,H), C6H6 (C,H),
    # triose (C,H,O). Tolerance is on max relative |diff| per (Zi,Zj,iref,jref).
    # After the GenerateParameters cache-rebuild fix (2026-05-28) the reference
    # C6 is bit-identical to dftd4 (~1e-10). NOTE: this fix is energy-neutral for
    # H/C/N/O at typical geometries (the corrected references carry ~0 CN/charge
    # weight) — it removes a latent reference-C6 bug, but the CH4/triose energy
    # residual is NOT here. See docs/GFN2_D4_STATUS.md.
    set(_C6_MOLS_TOL
        "H2O:1.0e-9"
        "NH3:1.0e-9"
        "CH4:1.0e-9"
        "C6H6:1.0e-9"
        "triose:1.0e-9"
    )
    foreach(_mt ${_C6_MOLS_TOL})
        string(REPLACE ":" ";" _mt_list ${_mt})
        list(GET _mt_list 0 _mol)
        list(GET _mt_list 1 _tol)
        add_test(
            NAME d4c6_dump_${_mol}
            COMMAND dump_dftd4_c6
                ${_SQM_MOL_DIR}/${_mol}.xyz
                ${_SQM_DUMP_DIR}/${_mol}_c6.json
        )
        set_tests_properties(d4c6_dump_${_mol} PROPERTIES TIMEOUT 60 LABELS "d4_c6")
        add_test(
            NAME d4c6_diff_${_mol}
            COMMAND diag_curcuma_d4_c6 --quiet --tol ${_tol}
                ${_SQM_DUMP_DIR}/${_mol}_c6.json
                ${_SQM_MOL_DIR}/${_mol}.xyz
        )
        set_tests_properties(d4c6_diff_${_mol}
            PROPERTIES TIMEOUT 60 LABELS "d4_c6"
                       DEPENDS d4c6_dump_${_mol})
    endforeach()
endif()
