# Native RHF: SCF convergence test and iteration count (Oct 25, 2026)

> 🤖 AI-generated, ⚙️ machine-tested. One machine (Core Ultra 9 285K, RTX 5090), n = 1 per number, not reviewed by
> the operator. Human production testing is pending.

**Default change** (`-qm.*`, native `hf`, `hf-3c`, `mp2`, `ri-mp2`, `lmp2`, `ccsd`): the SCF no longer stops on
`||P_new - P||_F < 1e-6` alone. New default `-qm.scf_convergence auto`: all of

| condition | parameter | default |
|---|---|---|
| energy change between two iterations, \|E_n - E_(n-1)\| | `-qm.scf_energy_threshold` | 1e-8 Eh |
| orbital gradient \|\|X^T (F P S - S P F) X\|\|_F (the DIIS error, Lowdin basis) | `-qm.scf_gradient_threshold` | 1e-5 |
| the same, when an analytic gradient or a correlation method follows | `-qm.scf_gradient_threshold_tight` | 3e-7 |
| density change \|\|P_new - P\|\|_F | `-qm.scf_threshold` | 1e-4 (was 1e-6, the only test) |

E_n = 1/2 Tr(P_n (H + F(P_n))) is the energy of the density whose Fock matrix was just built (variational);
the converged SCF still finishes with one more diagonalisation step and a full J/K build of that density, as before.
The old rule: `-qm.scf_convergence density -qm.scf_threshold 1e-6`. All four conditions are combined with AND, so a
set `-qm.scf_threshold` (all ctest checks use 1e-8..1e-10) still forces at least that tight an SCF.

Implementation: `ScfConvergence` in `qm_scf.cpp` (host loop and device-resident loop), device scalars from
`QMGpuSCF::scfMeasures` (`cuda/qm_gpu_scf.cu`). An energy-only SCF that is followed by a gradient request on the same
geometry is continued from its density with the tight threshold (`QMEngine::Calculation`).

## Why the old test cost iterations

The HF energy is stationary in the orbitals: its error is second order in the orbital gradient g. A density
criterion of 1e-7 (timing study) or 1e-6 (default) on the Frobenius norm keeps iterating long after the energy is
exact to 1e-10 Eh. Anthracene, 24 threads, per iteration of one tight run (`-qm.scf_convergence density
-qm.scf_threshold 1e-10`): E_n reaches its converged value to 1e-10 Eh at iteration 11 (g = 1.6e-5, dP = 2.3e-5),
the 1e-7 criterion stops at 19.

Final energies when the run is stopped early (the final energy always comes from a full build one step later):

| system | nbf | reference E_elec (Eh), tight (dP 1e-10, old binary for benzene/anthracene) | new default | g < 3e-5 | g < 1e-4 | old 1e-7 |
|---|---|---|---|---|---|---|
| benzene (stored tensor) | 114 | -433.6725435404 | 0 | - | - | 0 |
| anthracene | 246 | -1307.2299507207 | 0 | +1e-10 | +1e-10 | 0 |
| C20H42 | 490 | -2163.0355541443 | 0 (1T run: +1e-10) | 0 | 0 | 0 |
| C40H82 | 970 | -5029.1769178829 | 0 | 0 | - | 0 |
| C20H42, device loop | 490 | -2163.0355541443 | 0 | - | - | 0 |
| C80H162, device loop | 1930 | -11508.2026479164 (old 1e-7) | 0 | - | - | 0 |

Deviations in Eh at the 10 printed decimals (resolution 1e-10). The task's acceptance limit was 1e-9 Eh.

## Quantities that are first order in the orbital error

Analytic gradients and correlation energies (MP2, RI-MP2, LMP2, CCSD) are not stationary in the HF orbitals. Anthracene,
deviation from the tight reference:

| orbital gradient threshold | RI-MP2 E_corr (Eh) | HF gradient max dev. (Eh/Å) | iterations |
|---|---|---|---|
| 1e-6 | 1.9e-9 | 1.8e-8 | 14 |
| 3e-7 (`scf_gradient_threshold_tight` default) | 1.7e-9 | 6.1e-9 | 15 |
| 1e-7 | 5.2e-10 | 2.7e-9 | 16 |
| old default (dP 1e-6) | 1.7e-9 | 6.1e-9 | 15 |
| old 1e-7 (timing study) | 7e-11 | - | 19 |

The tight default reproduces the accuracy of the old default for these quantities (same iteration count); it is
not tighter. Below ~5e-10 Eh in E_corr the incremental integral-direct Fock builds set a floor (their screening
noise in the orbitals) unless the run reaches the full-build switch (`DirectFockCadence`).

## Energy threshold and incremental-build noise

E_n from an incremental (integral-direct) Fock build carries the screening noise of those builds: on C20H42 and
C40H82 it drifts by up to 3e-9 and 9e-9 Eh between iterations once g < 1e-6. With
`scf_energy_threshold 1e-9` C40H82 needed 19 instead of 11 iterations (it waited for the full-build switch);
1e-8 (ORCA TightSCF TolE) avoids that. The accuracy guarantee comes from the orbital gradient; the energy
condition is a safeguard.

## Convergence acceleration and guess (checked, unchanged)

At the new criterion (24 threads, iterations):

| variant | anthracene | C20H42 |
|---|---|---|
| default (atom guess, DIIS 8, start 1, early screening) | 12 | 11 |
| `-qm.scf_diis_subspace` 6 / 12 / 16 | 12 / 12 / 12 | 11 / 11 / 11 |
| `-qm.scf_mode adiis` | 12 | 11 |
| `-qm.scf_diis_start 2` | 12 | 11 |
| `-qm.scf_guess sad` | 16 | 15 |
| `-qm.scf_guess h0` | 17 | 30 |
| `-qm.scf_early_screening 0` | 12 (2.3 s vs 2.0 s) | 11 (9.1 s vs 7.0 s) |

The per-iteration contraction of curcuma's DIIS matches PySCF's (anthracene: dP 1.1e-4, 2.3e-5, 7e-6 at iterations
10-12 vs PySCF |ddm| 1.1e-4, 2.4e-5, 3.4e-6 at cycles 10-12, both with an atomic guess), so the remaining
difference to PySCF/ORCA iteration counts is the stopping rule, not the extrapolation. The atomic-RHF guess (default
since Oct 19) is the best of the three guesses; no SAP guess was implemented. Early screening and incremental builds
cost no iteration at the new criterion: the runs converge before the incremental-build noise floor, so the stall
switch does not fire. The energy condition is reset when the early screening or the FP32 phase ends (one more
iteration at most; in all runs above the gradient condition was still unmet at that point).

## Iterations and wall time, old binary vs new default

Old = build of `c60957c8` (before this change) with `-qm.scf_threshold 1e-7` (setting of the timing study) and with
its default 1e-6; new = this change, no SCF flags. `-sp -method hf -verbosity 2 -no_bmt`, def2-SVP, 1 thread pinned
with `taskset -c 2`, 24 threads unpinned; device loop: `-gpu cuda -qm.ri off -qm.gpu_direct on -qm.gpu_scf on`
(FP32 early phase on). Iterations = SCF loop iterations (each run adds one final full J/K build). n = 1.

| system | nbf | threads | old 1e-7: it / wall (s) | old 1e-6: it / wall (s) | new: it / wall (s) | speed-up vs 1e-7 | ORCA nori (s, study) | new / ORCA |
|---|---|---|---|---|---|---|---|---|
| anthracene | 246 | 1 | 19 / 48.2 | - | 12 / 33.2 | 1.45 | 42.7 | 0.78 |
| C20H42 | 490 | 1 | 16 / 129.3 | - | 11 / 107.3 | 1.21 | 113 | 0.95 |
| C40H82 | 970 | 1 | 16 / 548.2 | - | 11 / 472.3 | 1.16 | 488 | 0.97 |
| anthracene | 246 | 24 | 19 / 2.9 | 15 / 2.4 | 12 / 2.1 | 1.38 | 12.7 | 0.17 |
| C20H42 | 490 | 24 | 16 / 8.2 | 13 / 7.7 | 11 / 7.0 | 1.17 | 24.8 | 0.28 |
| C40H82 | 970 | 24 | 16 / 38.4 | 13 / 39.6 | 11 / 36.2 | 1.06 | 101 | 0.36 |
| C20H42 GPU | 490 | 24 | 16 / 2.4 | 14 / 2.2 | 12 / 2.0 | 1.20 | - | - |
| C80H162 GPU | 1930 | 24 | 17 / 30.0 | 14 / 27.3 | 12 / 25.4 | 1.18 | 471 | 0.05 |

The iterations drop by 30-37 %, the wall time by 6-31 %: the iterations that are no longer run are the cheap ones
(incremental builds on a small dP, where the density-weighted screening skips most quartets); the first eight or
nine iterations carry most of the J/K cost. At 1 thread curcuma is now below the ORCA 6.1 wall time of the study
(0.78-0.97x; ORCA not rerun), at 24 threads 0.17-0.36x. 24-thread runs vary by 5-10 % between repeats.

A looser early screening (`-qm.scf_early_screening 1e-9`: C20H42 6.9 -> 6.5 s, anthracene 2.3 -> 2.1 s, same
iterations; 1e-8: anthracene 12 -> 20 iterations) or a later switch (`-qm.scf_early_switch 1e-4`: +0 to +2
iterations) was measured once and not adopted.

## What was tested / not tested

- Tested: benzene, anthracene, C20H42, C40H82 (CPU 1 and 24 threads), C20H42 and C80H162 (device loop, RTX 5090),
  RI-MP2 and the HF gradient on anthracene, `-opt` hf and hf-3c on distorted H2O/formaldehyde at the defaults
  (same minima as before), `ctest -L qm` (274/276, the two `qm_opt_hf3c_*` need numpy and pass with the qmref venv).
- Not tested: small-gap or near-degenerate systems (the energy error is g^2 over the orbital Hessian; a small gap
  makes it larger), charged molecules, other bases than def2-SVP/MINIX, RI-JK mixed precision with the new rule
  (stage-2 exact builds are triggered by `scf_threshold`, now 1e-4), ROCm/Vulkan (no device loop there),
  MD energy conservation with the tight threshold.
- Not implemented: skipping the final full build when the last loop build was full (would save one build; done for energy-only SCFs on Oct 27, 2026, `-qm.scf_verified_stop`, [QM_SHARED_SCF_INFRA_2026-10-27.md](QM_SHARED_SCF_INFRA_2026-10-27.md)), SAP guess,
  second-order SCF.
