# GPU run-to-run reproducibility (CUDA GFN1/GFN2 and GFN-FF)

Status Oct 17, 2026, AI-generated, machine-tested. Known Issue #36 is the short form.

## Question and method

Two identical GPU runs should give the same result. Measured on PR18 (MOR41, GFN2, `-threads 1`, one device),
`CURCUMA_PSEUDO_TRACE=1` prints the SCF decision quantities with 17 digits (max|dq| of the previous step, HOMO/LUMO
gap, rotation angle of every pseudo-diagonalisation step), which exposes differences far below the 10 decimals of the
SCF table. Load was varied with 8 concurrent processes on 4 GPUs (two per device); the sequential case is the control.

## Findings (before the change)

| Setting | Runs | Result |
|---|---:|---|
| sequential, one process at a time | 40 | all 25 iterations identical to 10 decimals; gradient vs run 1: median 4e-15, max 4e-12 Eh/A |
| 8 concurrent processes | 32 | 9-10 runs (28-31 %) end in a different SCF state: the energy of iteration 22 differs by 9e-9 Eh, the gradient by 1.5e-7 (4.2e-7 in another sweep) Eh/A; 22 of 32 runs already differ in digit 16 of `last_dq` from iteration 5 on |
| same, host code pieces against host code pieces | 6 sweeps | the jump (4.2e-7) appears in 1 sweep: not caused by the device-side host-work changes of Oct 16-17 |

**Mechanism.** The traces of the diverging runs have the same structure (same pseudo steps taken, same gates), but the
FP32 quantities of one late pseudo step differ at 1e-9 relative (HOMO -0.35866244084708188 vs -0.35866244076341036).
The cause is FP64 `atomicAdd` in the loop, whose order depends on the load of the device: the shell and atom charges
(`k_qsh_scatter`, `k_qat_scatter`, `k_qsh_to_qat`), the atomic multipole moments (`k_multipole_moments*`) and the block
sums of the SCC energies (`k_energy_*`). That is 1e-16 noise in the charges. The mixed-precision path rounds the Fock
matrix to float; the noise occasionally moves an element over a float rounding boundary (1 ulp = 6e-8 relative), the
FP32 iterations carry that to ~1e-9 and the SCF ends in a neighbouring state, still inside the SCF threshold (1e-5)
but not the same state. Every outcome is a valid converged SCF; it just is not the same one.

## Change

The charge and moment sums use a fixed order now: CSR offsets of the AO -> shell, AO -> atom and shell -> atom maps
(built in `beginBasis`, used when the maps are sorted), one thread per shell/atom (`k_qsh_gather`, `k_qat_gather`,
`k_qat_from_qsh`), per-AO moment partials plus `k_mom_gather`, and per-block energy partials summed in index order
(`k_sum_part`). `CURCUMA_GPU_ATOMIC_SUMS=1` restores the atomic kernels (energy partials stay).

| Setting after the change | Runs | Result |
|---|---:|---|
| 8 concurrent processes, PR18 | 32 | one distinct 17-digit trace (32/32 identical); all final energies -49.26678083; gradient spread among the runs max 1.1e-16 |
| 4 GPUs, distributed eigensolver (peo201-water500, GFN2) | 5 | one distinct trace; energies identical; gradient spread 3-5e-16 |

Cost: none measurable (peo201-water500 GFN2 TOTAL 10.9 s, GFN1 18.0 s; before 11.3 / 18.4 s). MOR41 against the
all-host run: energies identical, gradient max 2.8e-14 (GFN2).

## What is still not bitwise

- **GFN-FF list builds** still take their slot from an atomic counter (`k_disp_pairs_build`, `k_generate_cn_pairs_write`);
  the order of a list no longer matters because every sum over it is an integer sum.
- **ROCm, Vulkan**: not looked at.
- CPU runs at different thread counts differ in the last bit (Known Issue #36).

## GFN-FF on the GPU: fixed-point accumulation (Oct 17, 2026, AI/machine-tested)

Integer addition is associative, so a sum of 64-bit integers does not depend on the arrival order of the threads.
`fix_add` (`gfnff_kernels.cu`) rounds a contribution to a multiple of 2^-52 (forces, dE/dCN, zz_hb) or 2^-56 (CN) and adds
it as `unsigned long long` into the 8 bytes of the double buffer; `k_fix_to_double` converts in place when the
accumulation is over (`k_fix_copy` for the diagnostic snapshots).

- CN pairs and HB CN (`k_cn_compute_pairs`, `k_hb_cn_per_atom`) -> converted before the log squashing / bond map.
- Fragment sums of the EEQ PCG: block per fragment over a CSR atom list (`k_frag_sum_csr`, built where `d_atom_frag` is
  uploaded); `k_eeq_reduce_sums_det` (one block, fixed tree); `k_eeq_reduce_fragment_sums_csr`.
- Forces: `add_grad` (all 62 sites), dispersion/bond dE/dCN, `zz_hb`, the shared-memory partner gradients of
  `k_coulomb_tiles`; the partial buffers of the helper devices (`k_add_coulomb_parts`, `k_add_disp_parts`) are added as
  integers. `d_dEdcn` is converted once before `k_coulomb_postprocess` / `k_cn_chainrule`, `d_grad` after the CN chain
  rule and before the download.

Measured (RTX A4500, `-method gfnff -gpu cuda`; AI-run, not human-reviewed):

| Check | Result |
|---|---|
| Gradient GPU vs CPU, polymer (1410 atoms) / caffeine (24) | max 1.7e-13 / 1.5e-13 Eh/A, energies equal to 8 decimals |
| MD polymer, 5000 steps, 8 concurrent processes on 4 GPUs vs 1 solo run | trajectory file (1000 frames) md5-identical in 9 of 9; with only the charge path fixed (forces still atomic) the same test differed from the solo run at printed rows 283-322 of 1001 in 8 of 8; the unchanged original binary was not run in this test |
| MD polymer_2x (7320 atoms), 300 steps, `-gpu_split_devices all` (4 GPUs) twice and `none` | trajectory file (6 printed digits) md5-identical in 3 of 3; full-precision final JSON identical for the two split runs, but split vs none differs at 1e-12 relative (Epot, Ekin): the EEQ matvec partition changes the charges in the last bits, so split and unsplit are not bitwise equal to each other |
| `ctest -R "gfnff\|gpu\|cuda\|eeq"` | 278 of 278 passed |
| Cost, polymer MD 5000 steps solo | 20 s before, 20-21 s after (3 runs) |

Limits: one system class per MD test (peptide-like polymer, polymer_2x), one GPU model, HB/XB host stages were not
varied, periodic systems and the EEQ split path with more than one fragment are not covered by the bitwise check. The
range of the accumulator is +-2048 Eh/Bohr per atom component (+-128 for CN-like terms, never reached); a clash beyond
that would wrap silently. ROCm (`rocm/gfnff_kernels_hip.h`) has its own kernels and is unchanged.

## GFN1/GFN2 gradient: same scheme (Oct 17, 2026, AI/machine-tested)

`k_grad_repulsion`, `k_grad_cn_onsite*`, `k_grad_h0_pulay*` and `k_grad_coulomb` add with `gfix_add` (2^-52, 64-bit
integer into the double buffer, `xtb_gpu_context.cu`); `k_grad_fix_to_double` converts `dGrad` and `dEdcn` before the
download.

| Check | Result |
|---|---|
| PR18 gradient (GFN2), 32 runs in 4 rounds of 8 concurrent processes on 4 GPUs + 1 solo run | dump file (energy and gradient, 17 digits) md5-identical in 33 of 33; before the change the gradient varied by 1e-16 to 5e-16 between runs |
| MOR41 (95 structures), GFN2, `-scf_threshold 1e-10`, GPU vs CPU, 1 thread | energies equal to 8 decimals, max gradient difference 1.08e-10 Eh/A (bounded by the SCF threshold; the run before the change was not repeated) |
| `ctest -R "sqm\|gpu\|cuda\|gfn"` | 430 of 431 passed; `cli_sqm_12` fails as before (BLAS, `docs/BROKEN_TESTS_2026-10-14.md`) |

Not covered: the post-SCF tails that already ran fixed-order (D4, D3, multipole), the multi-GPU cuSOLVERMp path, GFN1 on
the check set above (same kernels, not run separately).

## GFN-FF energies: fixed-point too (Oct 17, 2026, AI/machine-tested)

`blockReduceAddEnergy` rounds every thread value to an integer (2^-48 Eh; dispersion 2^-52, one rounded value per pair),
sums integers through warp, block and atomic, `k_coulomb_postprocess` and the helper `k_add_*_parts` add integers; the
14 slots are converted before the download (`k_fix_to_double`). The cause found on the way: the dispersion pair list is
built in atomic-counter order, so the thread-to-pair map changed between runs and the floating-point thread sums with it.

| Check | Result |
|---|---|
| polymer `-sp` gradient dump (energy 18 digits + gradient), 32 runs in 4 rounds of 8 concurrent processes, `-threads 1` | md5-identical 32 of 32; with the gradient already fixed-point, but FP energy blocks: 8 distinct energies (spread 1.8e-13 Eh) |
| Energy GPU vs CPU | polymer 4.6e-12 Eh (2^-44 first try: 1e-10), caffeine 2.5e-14 Eh |
| MD polymer 5000 steps, 8 concurrent + 1 solo | trajectory identical; final JSON (Epot, Etot, Ekin, geometry, velocities, gradient) one distinct value in 9 of 9 |
| `ctest -R "gfnff\|gpu\|cuda\|eeq"` | 278 of 278 |

Cost: energies carry a rounding error of up to 3.6e-15 Eh per thread value (random, grows with the square root of the
number of contributions). Range +-32768 Eh per term, dispersion +-2048 Eh; beyond that the sum wraps silently.

## Options for the rest (not done)

1. ROCm and Vulkan: not looked at.
2. Split vs unsplit GFN-FF (EEQ matvec partition) is not bitwise equal; same-configuration runs are.
