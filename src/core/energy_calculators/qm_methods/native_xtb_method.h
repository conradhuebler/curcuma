/*
 * <Native xTB Method Wrapper — unified GFN1 / GFN2>
 * Copyright (C) 2025 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software under GPL-3.0.
 *
 * Claude Generated: a single ComputationalMethod adapter for the native
 * curcuma::xtb::XTB solver, parametrized by MethodType. It replaces the former
 * GFN1Method / GFN2Method wrappers, which were near-identical delegation to the
 * same XTB instance — GFN2 additionally carried try/catch error handling,
 * element validation, and orbital accessors. Those are kept here for both
 * methods (all backed by XTB, which supports GFN1 and GFN2 alike). The only
 * per-method differences are the default config (D3 vs D4) and the method name.
 */

#pragma once

#include "../computational_method.h"
#include "xtb_native.h"
#include "xtb_fragment_scf.h"
#include "src/core/molecule.h"
#include "src/core/parameter_macros.h"

#include <memory>
#include <string>

#include "src/core/solvation/cpcmx.h"

// Parameters of the NATIVE GFN1/GFN2 solver (Claude Generated, moved here Sep 2026
// from xtbinterface.h, which wraps the EXTERNAL xtb/tblite binaries and has nothing
// to do with this code path). The registry module stays "xtb" (-xtb.<name>); the
// definitions live next to the solver they configure.
BEGIN_PARAMETER_DEFINITION(xtb)
    PARAM(d4_charge_source, String, "mulliken", "GFN2 D4 charge response in the gradient: 'mulliken' (exact, variational, as tblite), 'eeq' (approximate, single-shot EEQ) or 'cpscf' (Z-vector solve, ~100x slower). Only 'mulliken' matches finite differences to ~1e-8 Eh/A; 'eeq' and 'cpscf' are about 1e-4 Eh/A off on small molecules (test xtb_cpscf, Nov 2026). The energy is the same for all three.", "Dispersion", {})
    // Implicit solvation (Claude Generated, June 2026): self-consistent ALPB coupled
    // into the native GFN1/GFN2 SCF. Use the dotted form (-xtb.solvent water) because
    // the flat -solvent flag is ambiguous (also registered by tblite/ulysses/gfnff_external).
    PARAM(solvent, String, "none", "Implicit solvent for native GFN1/GFN2 (e.g. 'water', 'dmso', 'chloroform'); 'none' runs in the gas phase. Names are case-insensitive and accept the xtb aliases (dcm, dichloromethane, diethylether, thf, h2o, ...); an unknown solvent is an error. With solvent_model cosmo or cpcm a plain number is the dielectric constant (as xtb --cosmo 78.4, tblite --cpcm 78.4); cpcm takes epsilon of a name from tblite's solvent table (water 78.36). Use the dotted form -xtb.solvent, the flat -solvent is shared with other providers.", "Solvation", {})
    PARAM(solvent_model, String, "none", "Native solvation model: 'alpb', 'gbsa', 'cosmo' (ddCOSMO as xtb --cosmo), 'cpcm' (as tblite --cpcm: same ddCOSMO solver, keps / (1 + 0.571412), 74 points, eta 2, no surface term) or 'cpcmx' (CPCM-X as xtb --cpcmx: conductor COSMO SCF plus the post-SCF COSMO-RS/SMD step; single points only, no gradient; needs the CPCM-X solvent database, see cpcmx_db); numeric codes 0-3 accepted for none/cpcm/gbsa/alpb; 'none' with a solvent set means ALPB.", "Solvation", {})
    // CPCM-X (Claude Generated, Nov 2026): port of CPCM-X 1.1.0 as called by xtb 6.7.1 --cpcmx (src/core/solvation/cpcmx.h).
    PARAM(cpcmx_db, String, "", "Directory with the CPCM-X solvent database (the files DB/xtb/<solvent>.cosmo of github.com/grimme-lab/CPCM-X, LGPL-3.0; not shipped with curcuma). Empty: $CURCUMA_CPCMX_DB, then <executable dir>/cpcmx_db, then <executable dir>/../share/curcuma/cpcmx_db.", "Solvation", {})
    PARAM(cpcmx_smd_oxygen_radius, String, "atom8", "CPCM-X SMD term for solvents with Abraham acidity alpha < 0.43 (not water): which radius gets the SMD oxygen enlargement 1.8 (0.43 - alpha), added in Bohr units as in the Fortran sources. 'atom8' (default) reproduces xtb 6.7.1, whose CPCM-X enlarges the radius of the 8th atom whatever its element (nothing below 8 atoms, so the result depends on the atom order); 'oxygen' every oxygen atom, as CPCM-X 1.1.0; 'none' no enlargement.", "Solvation", {})
    PARAM(cpcmx_reference_rounding, Bool, true, "CPCM-X: round the solute segment data to the decimals of xtb's xtb.cosmo file (positions, charges, areas, potentials f14.9, coordinates f19.14, energy f21.10) before the post-SCF step, as xtb --cpcmx does by passing the data through that file. false uses full precision (measured effect <= 2e-10 Eh on four molecules, Nov 2026).", "Solvation", {})
    PARAM(solvent_state, String, "gsolv", "Reference state of the solvation free energy, as xtb --alpb water [state]: 'gsolv' (default), 'bar1M' (1 bar gas to 1 mol/L solution, +1.9 kcal/mol) or 'reference' (1 bar gas to the pure solvent). Case-insensitive; anything else is an error. COSMO and CPCM allow gsolv only (xtb adds no state shift to COSMO, tblite rejects one for CPCM).", "Solvation", {})
    PARAM(solvent_epsilon, Double, -1.0, "Dielectric constant for solvent_model cosmo or cpcm; -1 takes it from the solvent name. With no solvent name it alone switches COSMO/CPCM on (COSMO: surface tensions of xtb's default COSMO set). No effect on ALPB and GBSA, which take it from the solvent parameters.", "Solvation", {})
    // Native ddCOSMO discretisation (Claude Generated, Nov 2026); defaults are those of xtb 6.7.1 (cosmo.f90).
    PARAM(cosmo_grid, Int, 230, "COSMO: Lebedev points per atom of the surface-area grid; the cavity grid is the Lebedev grid nearest to half of it (xtb: 230 and 110). Larger grids are smoother and slower. The cosmo_* settings do not apply to solvent_model cpcm, which keeps tblite's fixed defaults (74 points, lmax 6, conv 1e-8, eta 2).", "Solvation", {})
    PARAM(cosmo_lmax, Int, 6, "COSMO: highest angular momentum of the spherical-harmonics expansion of the surface charge on each sphere (xtb: 6).", "Solvation", {})
    PARAM(cosmo_conv, Double, 1.0e-8, "COSMO: convergence threshold of the iterative ddCOSMO solver (norm of the increment; xtb: 1e-8). The gradient re-solves the adjoint problem to 1e-3 of it.", "Solvation", {})
    PARAM(cosmo_eta, Double, 0.2, "COSMO: width of the smooth switching region between buried and exposed surface points, relative to the sphere radius (xtb: 0.2).", "Solvation", {})
    // SCF knobs shared with other methods (scf_threshold, scf_max_iterations, scf_damping,
    // scf_allow_unconverged, scf_mode, scf_diis_start, scf_diis_subspace) are in module "scf"
    // (scf_parameters.h); -xtb.<name> still overrides them for gfn1/gfn2. Claude Generated (Oct 2026).
    PARAM(scf_guess, String, "eeq", "Initial charges: 'eeq' (EEQ charges, fewer iterations on large or polar systems), 'h0' (bare Hamiltonian) or 'fragments' (converge each molecule alone first; for large neutral clusters of small molecules).", "SCF", {})
    PARAM(sto6g_legacy_4sp, Bool, false, "Use xtb's older STO-6G tables for the 4s/4p shells instead of tblite's. Affects K, Ca and Ge-Kr only; set true to reproduce the xtb binary on these elements.", "Reference", {})
    PARAM(d4_pair_cutoff, Double, 0.0, "Real-space cutoff (Bohr) of the GFN2 D4 two-body term (energy, gradient and SCF potential alike; CPU and CUDA, a unit cell always uses 50). 0 = automatic: 50 (tblite), or 72 with scf_criterion xtb. 72 reproduces the energy of the xtb 6.7.1 binary (its dispersion matrix, build_wdispmat; xtb takes its gradient at 60). A value > 0 always wins. Matters only for systems extending past 50 Bohr.", "Reference", {})
    PARAM(d4_atm_cutoff, Double, 0.0, "Real-space cutoff (Bohr) of the GFN2 D4 three-body (ATM) term. 0 = automatic: 25 (tblite), or 40 with scf_criterion xtb (not with a unit cell). 40 reproduces the xtb binary on molecules larger than ~25 Bohr, at higher cost. A value > 0 always wins.", "Reference", {})
    PARAM(pbc_integral_cutoff_bohr, Double, 0.0, "Periodic cell only: radius (Bohr) of the lattice-image sum of the overlap and H0 integrals. 0 = tblite rule min(sqrt(50/min_alpha), 40). No effect without a unit cell.", "Periodic", {})
    PARAM(pbc_ewald_tolerance, Double, 1.0e-10, "Periodic cell only: tolerance of the Ewald sum of the isotropic Coulomb (and GFN2 multipole) matrices; sets alpha and the k-space cutoff from the real-space cutoff. No effect without a unit cell.", "Periodic", {})
    PARAM(pbc_ewald_r_cut_bohr, Double, 0.0, "Periodic cell only: real-space cutoff (Bohr) of the Ewald sum, minimum image, must stay below L_min/2; 0 = L_min/2. No effect without a unit cell.", "Periodic", {})
    PARAM(pbc_coulomb_switch_bohr, Double, 10.0, "Periodic cell only: distance (Bohr) at which the Klopman-Ohno kernel is switched to 1/r (quintic switch over the last Bohr, tblite rcut = 10). The lattice sum of the Klopman-Ohno remainder diverges without this switch; another value gives no tblite agreement. No effect without a unit cell.", "Periodic", {})
    PARAM(pbc_multipole_tblite_factors, Bool, false, "Periodic cell, GFN2 only: use tblite 0.4.0's reciprocal-space factors of the multipole Ewald (charge-dipole 2, charge-quadrupole +1) instead of the Ewald-consistent 1 and -1/3. true reproduces tblite's periodic GFN2 energy (alpha dependent, up to 1e-3 Eh at 8 A); false (default) is alpha independent. No effect without a unit cell. Diagnostic; see docs/WP-PERIODIC-XTB.md P2b.", "Periodic", {})
    PARAM(spin_polarization, Bool, false, "Spin polarisation of native GFN1/GFN2 (tblite spin-polarization): separate alpha and beta Fock matrices with the shell spin constants W; acts only with -spin N > 0 (closed shell is unchanged). Broyden SCF on the CPU and on -gpu cuda (host-driven device loop, FP64 or mixed precision); -gpu rocm|vulkan, -scf_gpu_partial_diag, -large_system_mode, -scf_mode other than broyden, -scf_precond and non-default -eigensolver give an error. AI-generated, machine-tested only; see SPIN_POLARISATION.md.", "SCF", {})
    PARAM(spin_polarization_scale, Double, 1.0, "Factor on all spin constants W of spin_polarization (tblite wscale); 1.0 is the tblite model, other values are for calibration only.", "SCF", {})
    PARAM(level_shift, Double, 0.2, "Virtual-orbital shift (Eh) for scf_mode level-shift. It is faded out near convergence, so the converged result is unshifted.", "SCF", {})
    PARAM(warm_start, Bool, true, "Start each SCF from the converged charges of the previous geometry (optimisation, MD). false always starts from scf_guess.", "SCF", {})
    PARAM(keep_diis, Bool, false, "Keep the DIIS/Broyden history across geometry steps instead of resetting it. Experimental: may help on smooth MD trajectories and hurt on large steps.", "SCF", {})
    PARAM(eigensolver, String, "mkl", "Eigensolver: 'mkl' (LAPACK dsyevd, fastest on CPU), 'native' (built-in divide-and-conquer), 'purify' (density purification, needs -electronic_temperature 0 and a gap) or 'lobpcg' (iterative, experimental).", "Eigensolver", {})
    PARAM(scf_mixed_precision, Bool, true, "Solve the eigenproblem in FP32 until max|dq| < scf_fp32_threshold, then in FP64; the converged result is FP64. CUDA GPUs with full-rate FP64 and Vulkan default to false; ROCm always uses it.", "Eigensolver", {"mixed_precision"})
    PARAM(scf_fp32_threshold, Double, 1.0e-3, "max|dq| below which the eigensolve switches from FP32 to FP64; keep it above scf_threshold. Smaller stays longer in FP32 (faster on FP64-weak GPUs, less safe); CUDA uses 1e-5 unless set.", "Eigensolver", {"fp32_threshold"})
    PARAM(scf_gpu_gfn1_resident, Bool, true, "GFN1 on the GPU: run the whole SCF iteration on the device, as for GFN2. false keeps the host-driven GFN1 loop with one device eigensolve per iteration.", "GPU", {})
    PARAM(gpu_eigen_strategy, String, "auto", "CUDA, periodic cell only: how the full FP64 eigensolves of the SCF run. auto = from the device attributes (FP64-weak device with cuBLAS FP64 emulation, compute capability >= 12: fp32_refine; otherwise fp64_direct), fp64_direct = Cholesky reduction + dsyevd in FP64, fp32_refine = FP32 solve + FP64 iterative refinement (Ogita-Aishima) with emulated FP64 GEMMs, falls back to fp64_direct when it misses the FP64 acceptance, emulated = fp64_direct with the cuBLAS/cuSOLVER FP64 emulation. Measured on an RTX 5080 only. No effect without a unit cell or on ROCm/Vulkan.", "GPU", {})
    PARAM(scf_gpu_partial_diag, Bool, false, "CUDA SCF: solve only the occupied eigenpairs plus a small buffer instead of the full spectrum; results are unchanged. Opt-in: the dense tridiagonalisation dominates, so the gain is small.", "GPU", {})
    PARAM(scf_pseudo_diag, Bool, true, "GPU SCF: in FP32 steps, rotate the previous orbitals once (first-order Jacobi) instead of a full eigensolve, within gap and step-size gates. On by default for CUDA; ROCm only when given.", "GPU", {})
    PARAM(scf_pseudo_diag_cpu, Bool, false, "CPU SCF: replace the eigensolve of middle iterations by one first-order orbital rotation (closed shell, same gates as scf_pseudo_diag); the final steps are full solves. Opt-in, experimental.", "Eigensolver", {})
    PARAM(scf_pseudo_diag_cpu_min_nao, Int, 1500, "Smallest number of basis functions for which scf_pseudo_diag_cpu is used; below it the extra full solve costs more than the rotation steps save.", "Eigensolver", {})
    PARAM(scf_pseudo_diag_max_dq, Double, 0.3, "Largest max|dq| of the previous step at which a pseudo-diagonalisation step may replace the eigensolve (GPU and CPU). Smaller values start the rotation steps later.", "Eigensolver", {})
    PARAM(scf_pseudo_diag_max_steps, Int, 12, "Most consecutive pseudo-diagonalisation steps before a full eigensolve refreshes the orbitals (GPU and CPU).", "Eigensolver", {})
    PARAM(scf_pseudo_diag_fp64, Bool, false, "With scf_pseudo_diag: use the rotation in FP64 SCF steps as well (same gates). Experimental: the gradient assumes exact eigenvectors. See GFN2_GPU_COST_PLAN.md.", "GPU", {})
    // Performance knobs that used to be environment variables or hard-coded heuristics
    // (Claude Generated, Sep 2026). Every default below is the measured one, so a plain run is
    // unchanged; they exist so a machine can be tuned without a rebuild and so
    // scripts/tuning_sweep.py can scan them. See docs/GPU_TUNING.md and docs/SQM_PERFORMANCE.md.
    PARAM(eigensolver_max_threads, Int, 0, "Thread cap for the eigensolve alone; 0 uses -threads, CURCUMA_EIG_MAX_THREADS overrides. The eigensolve is memory-bound, so fewer threads can be faster on machines with few memory channels.", "Eigensolver", {})
    PARAM(scf_reduce, String, "auto", "Reduction of the Fock matrix with the Cholesky factor of S: 'sygst' (fewer flops), 'trsm' (scales better with threads) or 'auto' (trsm from scf_reduce_threads up). Results agree to rounding.", "Eigensolver", {})
    PARAM(scf_reduce_threads, Int, 8, "Thread count from which scf_reduce auto uses trsm instead of sygst.", "Eigensolver", {})
    PARAM(scf_fp32_stall_patience, Int, 3, "Mixed precision: after this many FP32 iterations without 30 % progress in max|dq|, the rest of the SCF runs in FP64. 0 disables the guard.", "Eigensolver", {})
    PARAM(scf_fp32_false_fixpoint_factor, Double, 10.0, "Mixed precision: if an FP64 step shows a residual this many times larger than the last FP32 one, FP32 found a false fixed point and FP64 takes over. 0 disables the check.", "Eigensolver", {})
    PARAM(gpu_multipole_otf, String, "auto", "GPU GFN2: 'auto' stores the multipole interaction matrices up to ~1 GB (~2700 atoms; ROCm ~5770) and rebuilds them on the fly above; 'on' always rebuilds, 'off' always stores. Energies are identical.", "GPU", {})
    // ROCm resident-SCF solver switches (Claude Generated, Oct 2026), read by xtb_hip_method.cpp.
    PARAM(rocm_cached_cholesky, Bool, false, "ROCm SCF: factor the overlap once per geometry and solve each step by sygst/syevd/trsm instead of sygvd. Each factor is checked; on failure the step uses sygvd.", "GPU", {})
    PARAM(rocm_density_syrk, Bool, false, "ROCm SCF: build the density with a symmetric rank-k update (dsyrk) instead of a GEMM.", "GPU", {})
    PARAM(rocm_explicit_linv, Bool, false, "ROCm SCF with rocm_cached_cholesky: form L^-1 once per geometry and transform each step with trmm. Falls back to sygst/trsm if this fails.", "GPU", {})
    // Claude Generated (Oct 2026)
    PARAM(rocm_device_multipole, Bool, false, "ROCm GFN2, host-driven loop: compute the multipole potential and energy on the device, up to ~5770 atoms. No effect while rocm_device_potential is on.", "GPU", {})
    PARAM(rocm_device_potential, Bool, true, "ROCm GFN2, host-driven loop: build the whole SCF potential and the Fock matrix on the device, up to ~5770 atoms. false builds the potential on the host.", "GPU", {})
    // Claude Generated (Oct 2026), ROCm WP2
    PARAM(rocm_resident_loop, Bool, true, "ROCm GFN1/GFN2: run the whole SCF iteration, including the Broyden mixer, on the device. false keeps the host-driven loop; after a failure the SCF reruns there.", "GPU", {})
    PARAM(rocm_resident_open_shell, Bool, true, "ROCm with rocm_resident_loop: keep open-shell systems (n_alpha != n_beta) on the device loop. false runs them on the host-driven loop.", "GPU", {})
    PARAM(rocm_resident_precond, Bool, true, "ROCm with rocm_resident_loop: run -scf_precond frag on the device loop. false moves it to the host-driven loop, where the preconditioner is ignored with a warning.", "GPU", {})
    PARAM(rocm_device_solvation, Bool, true, "ROCm with native ALPB/GBSA: add the solvation reaction field on the device, so solvated runs keep the device potential and loop. Solvation energy and gradient stay on the host.", "GPU", {})
    PARAM(rocm_d4_postscf_device, Bool, true, "ROCm GFN2: compute the post-SCF D4 dispersion (energy, gradient, three-body term) on the device. Needs the resident loop; false uses the host code.", "GPU", {})
    PARAM(rocm_d3_device, Bool, true, "ROCm GFN1: compute the D3(BJ) dispersion energy and gradient on the device. false uses the host code.", "GPU", {})
    PARAM(rocm_cn_device, Bool, true, "ROCm: take the coordination numbers from the device (Hamiltonian CN; GFN2 also the D4 CN) instead of recomputing them on the host.", "GPU", {})
    PARAM(rocm_grad_resident, Bool, true, "ROCm after the resident loop (gas phase; GFN2 with d4_charge_source mulliken): build the gradient potentials on the device. false rebuilds them on the host.", "GPU", {})
    PARAM(rocm_grad_tail_device, Bool, true, "ROCm: compute the O(nat^2) part of the gradient (GFN2 multipole term, CN chain rule) and the repulsion energy on the device. false uses the host loops.", "GPU", {})
    PARAM(rocm_fixed_point, Bool, true, "ROCm: sum the device gradient and host-loop moments in 64-bit fixed point, so the result does not depend on thread order. false uses FP64 atomics.", "GPU", {})
    PARAM(rocm_otf_multipole, Bool, true, "ROCm GFN2 with rocm_device_potential: allow the on-the-fly multipole interaction (see gpu_multipole_otf). false builds the potential on the host above ~5770 atoms.", "GPU", {})
    PARAM(rocm_defer_host_integrals, Bool, true, "ROCm resident loop: download host copies of the device integrals only when a host step needs them. Not with d4_charge_source cpscf; false always builds them on the host.", "GPU", {})
    // Broyden seed preconditioner (Claude Generated, Oct 2026): opt-in, docs/SCF_PRECOND.md,
    // xtb_precond.h. Default 'none' leaves the SCF path unchanged.
    PARAM(scf_precond, String, "none", "Broyden preconditioner: 'none' (scalar step) or 'frag' (inverse SCC Jacobian of each molecule). Same result, fewer iterations on molecular clusters; neutral closed shell only. See SCF_PRECOND.md.", "Preconditioner", {})
    PARAM(scf_precond_guard, String, "g2:4:best", "Rule that drops the scf_precond matrix when it slows the SCF, e.g. 'g2:4:best'; 'off' disables it. Syntax in SCF_PRECOND.md; CURCUMA_SCF_PRECOND_GUARD overrides.", "Preconditioner", {})
    PARAM(scf_precond_rebuild_rms, Double, 0.0, "Rebuild the scf_precond frag matrix once the atoms moved more than this RMS (Angstrom) from where it was built. 0 never rebuilds for the same molecule.", "Preconditioner", {})
    PARAM(scf_precond_max_fragment_atoms, Int, 100, "Largest molecule (atoms) that gets its own scf_precond frag block; the build cost grows roughly as N^3.4. 0 = no limit; larger molecules use the scalar step or segments.", "Preconditioner", {})
    PARAM(scf_precond_segment_atoms, Int, 0, "With scf_precond frag: cut molecules above scf_precond_max_fragment_atoms into segments of at least this many heavy atoms (rings, C=O stay whole). 0 = off; about 3-12 is useful.", "Preconditioner", {})
    PARAM(scf_precond_segment_radius, Double, 2.0, "Radius (Angstrom) of the local cluster around a segment that its Jacobian block is computed on. Larger is more accurate and costlier; used when scf_precond_segment_atoms > 0.", "Preconditioner", {})
    PARAM(scf_precond_jacobian, String, "fd", "How scf_precond frag computes each Jacobian block: 'fd' (central differences of the SCC map) or 'response' (linear response, one eigensolve).", "Preconditioner", {})
    // Multi-step SCC extrapolation (Claude Generated, June 2026): opt-in generalisation
    // of the 1-step warm-start. Predicts the new-geometry SCC vector ([q_sh] for GFN1,
    // [q_sh; dp_at; qp_at] for GFN2) from the history of previous converged steps so the
    // SCF starts closer to the new fixpoint (fewer iterations, esp. in MD). Default
    // 'none' = the existing 1-step warm-start (unchanged). See docs/SQM_SCF_EXTRAPOLATION.md.
    PARAM(scf_criterion, String, "default", "SCF convergence rule of native gfn1/gfn2 (CPU loop). 'default': max|dq| < scf_threshold and |dE| < 100*scf_threshold. 'xtb': the rule of xtb 6.7.1, RMS of the full SCC residual (shell charges, dipoles, quadrupoles) divided by sqrt(number of atoms) < 1e-4*scf_accuracy (GFN2) or 2e-5*scf_accuracy (GFN1) and |dE| < 1e-6*scf_accuracy, followed by one more iteration that must satisfy it again; scf_threshold is ignored. 'xtb' also sets the GFN2 D4 cutoffs to xtb's 40 (ATM) and 72 (pair) Bohr unless d4_atm_cutoff / d4_pair_cutoff are given. Not implemented in the device-resident GPU loop (a warning is printed, the default rule runs).", "SCF", {},
        "enum=default|xtb")
    PARAM(scf_accuracy, Double, 1.0, "Factor on the xtb convergence limits for scf_criterion xtb (xtb 'acc'; xtb MD uses sccacc 2.0, i.e. -scf_accuracy 2). Larger is looser.", "SCF", {},
        "min=0")
    PARAM(scf_extrapolation, String, "none", "Predict the starting charges of a new geometry from previous steps: 'none' (last step only), 'aspc' (fixed-step MD) or 'gauss' (polynomial fit, irregular steps). The converged result is unchanged.", "Extrapolation", {})
    PARAM(scf_extrapolation_order, Int, 3, "Order of scf_extrapolation: for 'aspc' the predictor order k (uses k+2 steps), for 'gauss' the polynomial degree. Higher follows smooth MD better but can overshoot in optimisations.", "Extrapolation", {})
    PARAM(scf_extrapolation_apply, String, "guess", "'guess' seeds the SCF with the prediction; 'xlbomd' propagates the charges as extended-Lagrangian MD variables (MD only, with aspc; experimental). The SCF converges in both.", "Extrapolation", {})
    PARAM(scf_xlbomd_correctors, Int, 1, "XL-BOMD: minimum number of SCF cycles per MD step before convergence is accepted. 1 = no effect; experimental.", "Extrapolation", {})
    // Native GFN large-system fragmentation (Claude Generated, June 2026): opt-in
    // approximate scaling to >1000 atoms by exploiting locality. Default 'none' = the
    // exact dense path (unchanged). See docs/SQM_LARGE_SYSTEMS.md.
    //
    // Combination with -eigensolver: 'none' uses the chosen -eigensolver as usual.
    // 'fragments' PROPAGATEs -eigensolver to every sub-system (each fragment is a
    // small dense SCF; -eigensolver=purify forces T=0 per fragment).
    // 'dc' propagates -eigensolver per sub-block, BUT purify and lobpcg fall back
    // to the dense GES because DC needs the full spectrum for Fermi occupation /
    // chemical-potential bisection (purify gives projector eigenvalues 0/1, lobpcg
    // gives only a partial spectrum). Only 'native'/'mkl' apply directly.
    // 'sparse' IS the 0-K density-purification path — -eigensolver is ignored.
    // For 'fragments' or 'dc' with -eigensolver=purify, -electronic_temperature
    // MUST be 0 (hard error otherwise: purification is 0-K only by construction).
    PARAM(large_system_mode, String, "none", "Approximate large-system path: 'none' (exact dense), 'fragments' (one SCF per disconnected molecule), 'dc' (divide-and-conquer) or 'sparse' (O(N) purification, -electronic_temperature 0). See SQM_LARGE_SYSTEMS.md.", "LargeSystem", {})
    PARAM(large_system_buffer_bohr, Double, 10.0, "dc mode: buffer radius (Bohr) around the atoms of each cell. Larger converges to the dense energy at higher cost.", "LargeSystem", {})
    PARAM(large_system_cell_bohr, Double, 12.0, "dc mode: edge length (Bohr) of the spatial cells; smaller cells give more, smaller fragments.", "LargeSystem", {})
    PARAM(large_system_sparse_threshold, Double, 1.0e-6, "sparse mode: matrix elements below this magnitude are dropped. Smaller converges to the dense result at higher cost.", "LargeSystem", {})
END_PARAMETER_DEFINITION

/**
 * @brief Unified ComputationalMethod wrapper for native GFN1-xTB / GFN2-xTB.
 *
 * Construct with the MethodType; the factory passes GFN1 or GFN2. All
 * ComputationalMethod calls delegate to the wrapped curcuma::xtb::XTB.
 *
 * Claude Generated.
 */
class NativeXtbMethod : public ComputationalMethod {
public:
    explicit NativeXtbMethod(curcuma::xtb::MethodType method, const json& config = json{});
    ~NativeXtbMethod() override = default;

    // ---- ComputationalMethod interface ------------------------------------
    bool setMolecule(const Mol& mol) override;
    bool updateGeometry(const Matrix& geometry) override;
    double calculateEnergy(bool gradient = false) override;

    Matrix getGradient() const override;
    Vector getCharges() const override;
    Vector getBondOrders() const override { return Vector{}; }     // not computed
    Position getDipole() const override { return Position{0.0, 0.0, 0.0}; }  // TODO
    /// CPCM-X (solvent_model cpcmx) has no analytic gradient, as in xtb; every other setting has one.
    bool hasGradient() const override { return !m_cpcmx; }
    /// WP-PERIODIC-XTB (Nov 2026): true once the engine holds a cell (every term of the energy is then an image
    /// sum or an Ewald sum); SimpleMD then wraps fragments and removes only the translation (A6/A6b rules).
    bool hasPeriodicEnergy() const override { return m_xtb && m_xtb->hasPeriodicCell(); }

    std::string getMethodName() const override;
    bool isThreadSafe() const override { return true; }
    // Forward the global -threads budget to the solver as its intra-molecule
    // thread count. The solver auto-gates it (serial when run under molecule-level
    // parallelism or for small systems), so this is safe to set unconditionally.
    void setThreadCount(int threads) override
    {
        m_thread_count = std::max(1, threads);
        if (m_xtb) m_xtb->setIntraThreads(m_thread_count);
        if (m_xtb_conductor) m_xtb_conductor->setIntraThreads(m_thread_count);
        if (m_c1_driver) m_c1_driver->setIntraThreads(m_thread_count);
    }

    void setParameters(const json& params) override;
    json getParameters() const override { return m_parameters; }

    bool hasError() const override { return m_has_error; }
    void clearError() override { m_has_error = false; m_error_message.clear(); }
    std::string getErrorMessage() const override { return m_error_message; }

    // Optional advanced (QM) features — delegate to the native solver.
    Vector getOrbitalEnergies() const override;
    int getNumElectrons() const override;
    json getEnergyDecomposition() const override;
    bool saveToFile(const std::string& filename) const override;

    // ---- SCF warm-start / iterative-mode (Claude Generated) ---------------
    void setWarmStart(bool on) override;
    void resetWarmStart() override;
    void setIterativeMode(bool on) override;

    // Access the underlying dense native solver (null when a large_system_mode
    // driver is active). Used by the GPU wrapper to install an external
    // eigensolver hook. Claude Generated (GPU port).
    curcuma::xtb::XTB* solver() { return m_xtb.get(); }
    /// WP-PERIODIC-XTB P4/G1: EnergyCalculator reaches the periodic virial and the cell update (XTB::setWantVirial/
    /// virial/virialValid/updateCell) through this virtual, for the CPU wrapper and the GPU adapters alike.
    curcuma::xtb::XTB* xtbInstance() const override { return m_xtb.get(); }

    // ---- Native xTB extras (not in the base interface) --------------------
    Matrix getMolecularOrbitals() const;
    double getHOMOLUMOGap() const;
    double getHOMOEnergy() const;
    double getLUMOEnergy() const;
    Vector getCoordinationNumbers() const;

private:
    curcuma::xtb::MethodType m_method;
    std::unique_ptr<curcuma::xtb::XTB> m_xtb;
    // large_system_mode driver (opt-in, mode != none). When active, all
    // calculate/gradient/charge calls delegate here instead of m_xtb. Built in
    // setMolecule once the mode is known. Claude Generated.
    std::unique_ptr<curcuma::xtb::FragmentScfDriver> m_c1_driver;
    // CPCM-X (solvent_model cpcmx, Claude Generated Nov 2026): m_xtb runs the gas-phase SCF (E_gas, and all reported
    // properties, as xtb prints the gas-phase results last), m_xtb_conductor the SCF in xtb's conductor COSMO set whose
    // segments feed the post-SCF step. The energy is E_gas + dG_solv(CPCM-X).
    bool m_cpcmx = false;
    std::unique_ptr<curcuma::xtb::XTB> m_xtb_conductor;
    std::string m_cpcmx_solvent;                     ///< normalised CPCM-X solvent name
    std::unique_ptr<curcuma::solvation::cpcmx::Molecule> m_cpcmx_solvent_data;
    curcuma::solvation::cpcmx::Result m_cpcmx_result;
    double m_cpcmx_conductor_energy = 0.0, m_cpcmx_gas_energy = 0.0;
    double calculateCpcmx();
    std::string scopedString(const char* key, const std::string& fallback) const;
    bool scopedBool(const char* key, bool fallback) const;

    Mol m_molecule;
    bool m_calculation_done = false;
    double m_last_energy = 0.0;

    // Read large_system_mode / eigensolver / electronic_temperature from the
    // controller ("xtb" scope, top-level fallback). Used by the T=0 hard-error
    // on -eigensolver=purify when large_system_mode=fragments|dc.
    std::string largeSystemModeString() const;
    std::string eigensolverString() const;
    double      electronicTemperature() const;
    bool        spinPolarizationRequested() const;   // -xtb.spin_polarization (Claude Generated, Nov 2026)
    std::string scfPrecondString() const;   // -scf_precond, for the large_system_mode warning

    // Push the controller settings into the native solver: the D4 charge-response
    // source (harmless for GFN1, which uses D3) and the SCF-convergence settings
    // (mode/guess/damping/DIIS/level-shift) via curcuma::xtb::applyXtbScfConfig.
    void applyConfig();
    void handleError(const std::string& operation);

    static json getDefaultConfig(curcuma::xtb::MethodType method);
    static bool isMoleculeSupported(const Mol& mol);  // native xTB params cover Z = 1..86
};
