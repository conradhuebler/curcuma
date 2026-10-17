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

- **Gradient (GFN1/GFN2)**: the gradient kernels accumulate with FP64 `atomicAdd` (`k_grad_repulsion`,
  `k_grad_cn_onsite*`, `k_grad_h0_pulay*`, `k_grad_coulomb`): 1e-16 to 5e-16 Eh/A between runs. No decision depends on it.
- **GFN-FF on the GPU**: 38 atomic sites in `gfnff_kernels.cu` (force accumulation of the bonded terms, Coulomb tiles,
  dispersion, CN pairs), the EEQ reductions `k_eeq_reduce_sums` / `k_eeq_reduce_fragment_sums` (they set the
  charges), and pair-list builds whose slot order comes from an atomic counter (`k_disp_pairs_build`,
  `k_generate_cn_pairs_write`). MD on polymer (1410 atoms, 1000 steps): 3 sequential runs identical to the printed
  10 decimals at 300 steps on 1 GPU and on 4 GPUs; with 8 concurrent runs the Epot of the runs first differs from a
  solo run at print index 414 or 611 (of 1001, one print per step in units of 1 fs) and stays within 1e-10 up to step 1000.
- **ROCm, Vulkan**: not looked at.
- CPU runs at different thread counts differ in the last bit (Known Issue #36).

## Options for the rest (not done)

1. GFN-FF EEQ reductions and pair-list builds deterministic (small; they feed the charges and the summation order).
2. GFN-FF and GFN1/GFN2 gradient accumulation: per-term partial buffers with a fixed-order reduction per atom, or
   64-bit fixed-point atomics (integer addition is associative; resolution/range trade-off 2^-56 / +-128 Eh/Bohr).
   The first costs memory traffic, the second changes the numbers at 1e-17 per term.
3. A `-gpu_deterministic` switch that bundles 1 and 2 and disables nothing else.
