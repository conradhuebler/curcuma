# Spin polarisation of native GFN1/GFN2 (`-xtb.spin_polarization`)

Status: 🤖 AI-implemented, ⚙️ machine-tested (Nov 21, 2026). Not ✅ TESTED: human production testing is pending until the operator removes this line.

## What it is

`-spin N` alone (see [SCF_MODES.md](SCF_MODES.md), "Open shell") fills one orbital set with two Fermi fillings (n_alpha, n_beta). That is the
tblite model *without* spin polarisation: the triplet O2 lies 1.65 kcal/mol above the singlet. With `-xtb.spin_polarization true` (flat
`-spin_polarization true` is routed to the module `xtb`, as is `-xtb.spin_polarization_scale`) the SCF runs the tblite `spin-polarization`
container with two Fock matrices and the spin constants W (tblite `src/tblite/spin.f90`, nspin = 2):

- Shell magnetisation `m_l = sum_{mu in l, nu} (Pa - Pb)_{mu nu} S_{nu mu}` (alpha minus beta Mulliken population of shell l).
- `E_spin = 1/2 sum_A sum_{l,l' on A} W^A_{ll'} m_l m_l'` (all shell pairs of one atom; W < 0), added to the electronic energy.
- `u_l = sum_{l'} W_{ll'} m_{l'}`; `F^sigma = F + sigma * 1/2 S_{mu nu} (u_mu + u_nu)`, sigma = +1 (alpha), -1 (beta); `F` is the charge Fock matrix.
- Both spins have their own orbitals and Fermi filling (`fermiFillingChannel`); `P = Pa + Pb` feeds every other term unchanged. `m_sh` is
  appended to the Broyden vector (`[q_sh | dp | qp | m_sh]`) and to the warm-start state.
- Gradient: W has no geometry dependence. Two changes only (`xtb_gradient.cpp`): the energy-weighted density sums both orbital sets, and the
  overlap-gradient coefficient gets `+ Pm_{mu nu} (u_a + u_b)` with `Pm = Pa - Pb` (opposite in sign to the charge term `- P (v_a + v_b)`).
- Acts only with `-spin N > 0`; with `-spin 0` nothing changes. Without the option all results are bit-identical to before.
- `spin_polarization_scale` multiplies all W (tblite `wscale`); calibration only.

Code: `XTB::solveEigenUhf`, `spinPotential`, `energySpin`, `buildSpinFocks` (`xtb_scf.cpp`); constants `parameters/spin_constants.hpp`.

## Constants: origin and open points

86 elements x 6 values (ss sp pp sd pd dd), copied value for value from tblite **v0.8.0** (`src/tblite/data/spin.f90`), indexed by angular momentum, not
by shell position (transition metals are ordered d,s,p). 85 of the 86 rows differ from the tblite commit pinned in `CMakeLists.txt` (9ca469a, Dec 2024).
- **TODO operator: licence.** tblite is LGPL-3.0-or-later; whether the table may be copied into this GPL-3.0 tree is not decided.
- **TODO operator: literature reference.** The tblite source gives none (only the PR title "Updated spGFN2 spin constants (#155)"); not looked up, not guessed.
- The tblite pin in `CMakeLists.txt` still gives the old constants for `USE_TBLITE` runs with spin polarisation; raise it before comparing with `tblite-gfn2`.

## Not supported (hard error, no silent fallback)

`-gpu vulkan`, `-scf_gpu_partial_diag` with `-gpu cuda|rocm`, `-large_system_mode fragments|dc|sparse`, `-scf_mode` other than `broyden`, `-scf_precond frag`
(CPU and GPU), `-eigensolver` other than `mkl`, `-scf_extrapolation_apply xlbomd`. The analytic response (`xtb_response.cpp`) already rejects open shells. `-spin` stays global (no spin per fragment).
Not implemented: UHF without the spin constants, spin contamination analysis, orbital output for the beta set (`C`/`eps` hold the alpha set).

## What was tested (Nov 21, 2026, release build, CPU, `-threads 1`)

Reference: tblite 0.8.0 Python wheel, `add("spin-polarization", 1.0)`, same electronic temperature.

- **Energies** against tblite, GFN2 and GFN1, 300 K and 2000 K: H, O, OH (uhf 1, 3), H2 (0, 2), O2 (0, 2, 4), Pt7 (0, 2, 4, 6, 8): all differences <= 2.8e-8 Eh
  (most < 1e-8; the 2000 K Pt7 values 1e-8 to 2.8e-8 at the default `-scf_threshold`). O2 triplet GFN2: -7.931969836872 vs -7.931969836906 Eh.
- **Gradients** against tblite (Eh/Bohr): O2, O2 tilted, OH, H2, Pt7 (uhf 4, 6, GFN1 and GFN2), Pt7 + 2 H2 + O2 (uhf 2, 4, 2000 K): max difference 4e-8, mostly < 1e-8.
- **FD gate** (`-dump_gradient` against central differences of `-sp`, Angstrom steps, `-scf_threshold 1e-10`): O2 (uhf 2, 4), O2 GFN1, OH, H2 (uhf 2), Pt7 uhf 4, 6, GFN1 uhf 4: error scales
  with h^2 (O2: 2.1e-6 at h = 1e-3, 5.3e-7 at 5e-4); smallest 5e-8 Eh/A (Pt7, h = 5e-4). Pt7 uhf 2 at 300 K is bistable: displaced geometries land on another SCF solution (3.5 mEh apart); only
  the 12 of 21 components that stay on one branch are comparable (6e-8 Eh/A); at 2000 K all 21 are (2.4e-7).
- **Regression**: without the option 11 cases (water, O2 spin 0 and 2, OH, H2, Pt7, Pt7 + 2 H2 + O2; GFN1 and GFN2) give bit-identical energies and `-dump_gradient` files against the binary of `267dc325`;
  `ctest -R "cli_sqm_|xtb_gradient|gradient_unit_contract|sqm_"`: 343 of 343; `ctest -L structinv`: 168 of 172, the 4 failures (`reactions_doebereiner_control_gfn2`) fail identically with the old binary (SCF not converged).
- **ctest** `cli_sqm_16_spin_polarization`: triplet O2 energy and gradient (GFN1, GFN2) against tblite, option off, `-spin 0`, three excluded modes.
- **Pt7 + 2 H2 + O2** (structure `reactions-doebereiner-pt7`, GFN2, 2000 K): see the spin ladder below.

## Spin ladder of the Doebereiner Pt7 + 2 H2 + O2 start geometry (GFN2, `-sp`, spin_polarization true)

| uhf | native 2000 K | tblite 2000 K | native 300 K | tblite 300 K |
|---|---|---|---|---|
| 0 | -41.94941114 | -41.94941116 | not converged (500 it) | not converged |
| 2 | -41.95057830 | -41.95057830 | -41.88990742 | not converged |
| 4 | -41.95253772 | -41.95253773 | -41.90091402 | not converged |
| 6 | -41.96038415 | -41.96038416 | -41.90283926 | -41.89334208 (other solution) |
| 8 | -41.96143396 | -41.96143398 | -41.90507124 | -41.90507123 |
| 10 | -41.95432148 | -41.95432150 | -41.90141416 | -41.90141416 |

The minimum is uhf 8 at both temperatures (not 4 or 6), 0.66 kcal/mol below uhf 6 at 2000 K. Short MD (cold Pt7 input, spin 8, 1000 steps, 1 thread, machine under foreign load 25): 38.3 s with the option
against 30.1 s without (26 against 33 steps/s, not half), 20.4 against 20.0 SCF iterations per step (200 steps), no SCF error.

## Known limits

- The SCF has several solutions (different magnetisation patterns). The result depends on guess and damping: GFN1 full system uhf 4, 2000 K: `-scf_guess eeq` (default) -42.50317, `-scf_guess h0` -42.51162 (tblite -42.51162).
  Multi-solution cases were not resolved; compare energies only where both programs converge to the same solution. At 300 K the 13-atom system often needs more than 500 Broyden cycles (also in tblite).
- Not tested: MD beyond 1000 steps, `-opt` of open-shell metal clusters (only a smoke run on O2), solvation, `-hessian`, periodic or charged systems, elements other than H, O, Pt, GFN1 d-shell systems other than Pt.
- GPU: CUDA and ROCm (sections below); Vulkan is a hard error.
- Solvation (checked Nov 22, 2026 after the merge with the native ddCOSMO/CPCM/CPCM-X work, n = 1, OH radical, GFN2, `-spin 1`, water): ALPB, cosmo and cpcm change the energy with spin polarisation as without it
  (cpcm -5.1e-4 Eh both, cosmo +1.1e-3 against +1.2e-3, alpb -9.9e-3 against -9.8e-3); no external reference exists for the combination (tblite 0.4.0 has no spin polarisation) and no test covers it.

## GPU (CUDA, `-gpu cuda`)

Status: 🤖 AI-implemented, ⚙️ machine-tested on one RTX 5080 (CUDA 13, Nov 22, 2026). Human production testing pending.

- Runs on the host-driven device loop: the host builds the charge potential v and the spin potential u (as on the CPU), the device solves
  F^alpha = F(v - u) and F^beta = F(v + u) (the spin shift +- 1/2 S (u_mu + u_nu) has the form of the isotropic potential term, so the existing
  Fock builders are reused), keeps both orbital sets (`dCa`, `dC`), builds Pa, Pb, P = Pa + Pb, Pm = Pa - Pb and the Mulliken populations of P
  and Pm; the host does the two Fermi fillings, m_sh, E_spin and the Broyden step on `[q_sh | dp | qp | m_sh]`.
- Gradient on the device: W' = Wa + Wb - 1/2 Pm o (u_mu + u_nu) (`XtbGpuContext::spinWeightedDensity`), read by the unchanged gradient kernels
  in place of W, which reproduces the CPU term + Pm (u_a + u_b). If the device gradient fails, the host gradient runs with W of both spins.
- Mixed precision (`-scf_mixed_precision`, default on for FP64-weak GPUs) is supported: both solves of a step use the same precision and no FP32
  step is accepted as converged (unchanged rule).
- Not used with spin polarisation (automatically, GPU still used): the GFN2 device potential build, the fully device-resident loops (GFN1 and
  GFN2), pseudo-diagonalisation (resident loop only). Screened integral storage works (`-gpu_sparse_integrals on`, checked on the 13-atom system
  only, where the pattern holds 100 % of the pairs). Multi-GPU eigensolver/density split: not run (one GPU).
- Code: spin branch of the SCF loop in `XTB::Calculation` (`xtb_native.cpp`), `GpuScfBackend::solveSpin/densitySpin/finalizeSpin/spinWeightedDensity`
  (`xtb_native.h`, CUDA in `xtb_gpu_method.cpp`), `XtbGpuContext::residentSpinKeepAlpha/residentDensitySpin/residentFinalizeSpin/spinWeightedDensity`
  (`cuda/xtb_gpu_context.cu`). ctest `cli_sqm_22_spin_polarization_cuda`.

Measured (Nov 22, 2026, GPU against the CPU of the same build, `-scf_threshold 1e-8`, `-threads 1` on the CPU):
- O2 (uhf 2, 4; GFN2, GFN1 uhf 2), OH (uhf 1, GFN2 and GFN1), H2 (uhf 2), Pt7 (`pt-cluster-7`, uhf 4 and 6, 2000 K, GFN2 and GFN1), Pt7 + 2 H2 + O2
  (`reactions-doebereiner-pt7`, 2000 K, GFN2 uhf 4 and 8, GFN1 uhf 4): FP64 energy <= 1.7e-14 Eh, gradient <= 5.7e-12 Eh/A; mixed precision energy
  <= 3.6e-14 Eh, gradient <= 4.2e-9 Eh/A. The Doebereiner system at uhf 2 does not converge in 500 cycles on the CPU nor in FP64 on the GPU at
  `-scf_threshold 1e-8` (the mixed-precision run converged to -41.95068297, another solution), and uhf 6 converges to -41.95410536 on both (the
  ladder below lists -41.96038415): multi-solution case, not compared.
- FD gate (GPU `-dump_gradient` against central differences of GPU energies, FP64, `-scf_threshold 1e-10`): O2 uhf 2 GFN2 2.1e-6 (h = 1e-3) and
  5.3e-7 Eh/A (h = 5e-4); Pt7 uhf 4 2000 K GFN2 1.6e-7 / 3.3e-8, GFN1 6.3e-8 (h = 5e-4); Pt7 + 2 H2 + O2 uhf 8 GFN2 3.1e-7 Eh/A (h = 5e-4, 39 components).
- Without the option: 10 cases (8 GPU, 2 CPU; caffeine GFN1/GFN2, O2 spin 2, OH spin 1, Pt7 spin 8, Pt7 + 2 H2 + O2 spin 4, trisaccharide GFN2 in mixed
  precision and FP64, GFN1 CPU) give the same energies (8 decimals) and bit-identical `-dump_gradient` files as the binary of f8b0fa3f.

Timing (Nov 22, 2026, GFN2, RTX 5080 + 32 CPU threads, n = 2 each, the two runs within 2 % except CPU 16 threads on the small systems; foreign load
average 8 to 13 from another session's curcuma runs during the measurement, so CPU numbers may be pessimistic). GPU = `-gpu cuda -threads 16`, default
mixed precision; FP64 = `-scf_mixed_precision false`. SCF = the "SCF converged" time of `-sp`; MD = wall time of 200 steps (`-dt 0.5 -max_time 100`).

| System | atoms, uhf | SCF CPU 1 thr | SCF CPU 16 thr | SCF GPU | SCF GPU FP64 | MD CPU 1 thr | MD CPU 16 thr | MD GPU | MD GPU FP64 |
|---|---|---|---|---|---|---|---|---|---|
| Pt7 + 2 H2 + O2, 2000 K | 13, 8 | 31 ms | 21-28 ms | 213 ms | 191 ms | 5.3 s | 5.0-6.1 s | 15.2 s | 20.7 s |
| trisaccharide minus one H | 65, 1 | 100 ms | 64-96 ms | 224 ms | 254 ms | - | - | - | - |
| pep20 peptide minus one H | 308, 1 | 4.65 s | 1.92 s | 1.07 s | 1.41 s | 421 s | 128 s | 58 s | 102 s |

The GPU is slower than the CPU up to at least 65 atoms (per-iteration launch and transfer cost of the host-driven loop) and 1.8x (SCF) / 2.2x (MD)
faster than 16 CPU threads at 308 atoms. The same 308-atom radical without spin polarisation (fully resident loop): SCF 665 ms on the GPU against
1.08 s on 16 CPU threads, i.e. the spin branch costs 1.6x the closed-loop GPU time per SCF (two eigensolves, host potential).

## GPU (ROCm, `-gpu rocm`)

Status: 🤖 AI-implemented, ⚙️ machine-tested on one RX 9060 XT (gfx1200, Nov 24, 2026). Human production testing pending.

- Port of the CUDA path above, same split between host and device: `HipScfBackend::supportsSpinPolarization/solveSpin/densitySpin/finalizeSpin/spinWeightedDensity`
  (`xtb_hip_method.cpp`), `XtbHipContext::residentSpinKeepAlpha/residentDensitySpin/residentFinalizeSpin/spinWeightedDensity` with the kernels
  `k_spin_combine`, `k_spin_w_shift` (`rocm/xtb_hip_context.hip`); the device gradient reads W' once (`w_spin_override`), as on CUDA.
- Not used with spin polarisation (automatically, GPU still used): the device potential build (`-xtb.rocm_device_potential`), the device multipole
  potential, both fully device-resident loops, the device gradient potentials (`-xtb.rocm_grad_resident`), pseudo-diagonalisation; post-SCF D4 runs on
  the host (no device D4 reference tables on the host-driven loop). `-xtb.rocm_density_syrk` is not applied to the two spin densities (GEMM).
- Screened integral storage (`-gpu_sparse_integrals on|auto`) with spin polarisation on ROCm (since Nov 25, 2026, SPE; 🤖 ⚙️): P = Pa + Pb,
  Pm = Pa - Pb and W' live on the screened pattern (`rocm/xtb_hip_sparse_impl.hiph`, `*Sp` twins of the four spin functions, kernel
  `k_spin_w_combine_sp`); CUDA keeps P, Pm and W dense next to the screened S/H0. Screened vs dense: <= 1.4e-12 Eh, <= 2.0e-13 Eh/A (O2, pep20 and
  peo201-chain triplets, n = 1 each); geometry optimisation, MD, solvation and `auto` with spin untested ([ROCM_PARITY_SPE.md](ROCM_PARITY_SPE.md)).
- Measured (Nov 24, 2026, GPU against the CPU of the same build, `-scf_threshold 1e-10`, `-threads 1`): O2 (uhf 2, 4; GFN2, GFN1 uhf 2), OH (uhf 1, GFN2
  and GFN1), H2 (uhf 2), Pt7 (`pt-cluster-7`, uhf 4 and 6, 2000 K, GFN2 and GFN1), 10 cases: FP64 energy <= 2.5e-14 Eh, gradient <= 1.4e-12 Eh/A;
  mixed precision energy <= 4.7e-14 Eh, gradient <= 7.6e-12 Eh/A. O2 uhf 2 with ALPB water (GFN2, GFN1): energy <= 3.6e-15 Eh, gradient <= 2.9e-15 Eh/A.
- FD gate (GPU `-dump_gradient` against central differences of GPU energies, FP64, `-scf_threshold 1e-10`): O2 uhf 2 GFN2 5.3e-7 (h = 5e-4) and 1.3e-7 Eh/A
  (h = 2.5e-4), GFN1 6.1e-7 / 1.5e-7 (h^2 scaling, truncation error); Pt7 uhf 4 2000 K GFN2 3.3e-8, GFN1 6.3e-8 Eh/A (h = 5e-4, 21 components each).
- Two runs md5-identical (`-dump_gradient`, default mixed precision; O2 GFN2, Pt7 GFN2 and GFN1). Without the option: 10 cases (8 GPU, 2 CPU; caffeine
  GFN1/GFN2 in mixed precision and FP64, O2 spin 2, OH spin 1, Pt7 spin 8 and 4, option explicitly off) give bit-identical `-dump_gradient` files as the
  binary without the port.
- Not tested: systems above 13 atoms, MD, `-opt`, d-shell main-group elements, multi-GPU, gfx other than gfx1200, timing. ctest
  `cli_sqm_23_spin_polarization_rocm` (also checks that `-scf_gpu_partial_diag`, `-scf_precond frag` and `-large_system_mode fragments` stop with an error).

