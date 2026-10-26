# ROCm vs CPU vs CUDA: evaluation of GFN-FF and GFN2 on gfx1200 (Oct 13, 2026)

Status: AI-measured, machine-tested only. Human review pending. Not a validation of the ROCm path.

## Setup

- Host: Intel Core Ultra 5 225 (10 threads, `-threads 10`), AMD Radeon RX 9060 XT (gfx1200, 16 GB), ROCm in `/opt/rocm`.
- Binary: `release/curcuma`, `USE_ROCM=ON`, `ROCM_GPU_ARCH=gfx1200`, source commit `9d65c6b1`, clean working tree.
- CPU run: `-gpu none`; GPU run: `-gpu rocm`; same structure files from `test_cases/structures/`, fresh directory per run, `-no_bmt`.
- No CUDA hardware on this host: nothing in this document is a measured CUDA number. CUDA statements come from code and docs.
- Scripts and logs: session scratchpad (not kept in the repository). Commands are listed under "Reproduce".

## Result 1: energies equal CPU

`-sp`, energies agree to all 8 printed decimals (Eh):

| Method | Structure | Atoms | CPU = ROCm |
|---|---|---|---|
| gfnff | caffeine | 24 | -4.67273706 |
| gfnff | peo201-chain | 1410 | -203.51759212 |
| gfnff | peo201-water500 | 2910 | -367.39771445 |
| gfnff | peo201x2-water1500.gfnff-min | 7320 | -917.51616997 |
| gfnff | urea400-water1000 | 6200 | -857.07822827 |
| gfn2 | caffeine | 24 | -42.14723025 |
| gfn2 | trisaccharide-c18h32o16 | 66 | -119.90699419 |
| gfn2 | host-guest frame (frame 1 of c33h45n9-agal-frames44) | 114 | -166.61385997 |
| gfn2 | peo201-chain | 1410 | -2088.25340678 |

MD step 0 (`-md`, gfnff): Epot, Ekin and Etot identical on peo201-chain and urea400-water1000.
Gradients and trajectories beyond step 0 were not compared in this evaluation.

## Result 2: speed, ROCm is not faster than the CPU here

Wall time of the whole process (includes plugin load and HIP init). `-sp` single runs unless noted.

| Case | Atoms | CPU | ROCm |
|---|---|---|---|
| gfnff `-sp` peo201-chain | 1410 | 0.12 s | 0.29 s |
| gfnff `-sp` peo201-water500 | 2910 | 0.34 s | 0.69 s |
| gfnff `-sp` peo201x2-water1500 | 7320 | 1.57 s | 2.48 s |
| gfnff `-sp` urea400-water1000 | 6200 | 1.58 s | 3.38 s |
| gfnff MD, peo201x2-water1500 | 7320 | ~150 ms/step | ~208 ms/step |
| gfn2 `-sp` frame 1 | 114 | 0.11 s | 0.68 s |
| gfn2 `-sp` trisaccharide | 66 | 0.05 s | 0.48 s |
| gfn2 `-opt` to convergence, frame 1 (n=3 each) | 114 | 22.0 s | 62.1-63.8 s |
| gfn2 `-opt` to convergence, trisaccharide (n=3) | 66 | 6.3-6.4 s | 20.5-23.1 s |
| gfn2 `-sp` peo201-chain | 1410 | 18 s | 16 s |

Caveats:
- MD ms/step is the difference of a 10 and a 40 `-MaxTime` run divided by the added time, one run each: coarse.
- gfn2 CPU time at 1410 atoms probably uses the large-system mode (not checked), so it is not a like-for-like SCF comparison.
- Small systems are dominated by process start and device init; they say little about the per-step cost in a long run.
- The `-opt` rows are full optimisations to convergence. `-opt.max_iter 5` (and `-MaxIter 5`) has no effect (same wall time with and without, checked afterwards on trisaccharide gfn2: 6 s each); an earlier version of this document wrongly called them 5-iteration runs. The step count was not recorded and may differ between CPU and ROCm, and the energy regex of the script failed, so the energies of these runs were not compared (the `-sp` energies above were).
- The docs report ROCm gfn2 about 6-8x and ROCm gfnff MD about 2.3x faster than CPU (Radeon 890M, 16-thread class CPU, `SQM_ROCM.md`). Not reproduced on gfx1200 with this CPU. No cause established; the gfx1200 FP64 rate and the host CPU differ from the documented machine.

## Result 3: ROCm is not optimised like CUDA (from the code)

Counts by `grep -c __global__` / unique `k_*` names (counted once, 2026-10-13):

| Item | CUDA | ROCm |
|---|---|---|
| GFN-FF unique kernel names | 73 | 54 |
| xTB context file `__global__` count | 70 (`xtb_gpu_context.cu`) | 26 (`xtb_hip_context.hip`) |
| GFN-FF device source | 3 large files (workspace, EEQ solver, kernels) | one hipcc TU `gfnff_rocm.hip` (7238 lines) + kernels/EEQ headers |

Present on CUDA only (kernel names): EEQ iterative solvers (`k_ppcg_*`, `k_eeq_block_jacobi_apply`, `k_eeq_diag_inv`, `k_li_*` local-inverse preconditioner), fragment projection (`k_frag_*`), `k_coulomb_tiles`, `k_coulomb_implicit`.
Present on ROCm only: `k_coulomb_gather`, `k_coulomb_dense`, `k_dispersion_gather` (per-atom gather instead of atomics).

Functional differences (details in `SQM_ROCM.md`, `GPU_TUNING.md`, `MULTI_GPU_GAPS.md`):
- EEQ on ROCm: `-eeq_solver.solve_method` accepts `auto`, `cholesky` and `ppcg` (Oct 2026; the Oct 13 text said only `auto`/`cholesky`, see `EEQSolver::gpuSupports`). The CPU solver runs only for 16 or more fragments with an explicitly set `-gfnff.eeq_rocm_cpu_fragment_threshold` or when the solver does not resolve to ppcg.
- The four opt-in GFN-FF tuning switches of `GPU_TUNING.md` are ignored on ROCm.
- xTB: no distributed eigensolver, so no split of one molecule over several GPUs (G2-1). GFN2 potential, the multipole SD/DD/SQ gradient and CN chain rule stay on the host; the GFN2 D4 pair gradient (`k_d4_grad`, `k_d4_atm`) runs on the device since Oct 2026, the GFN1 D3 gradient stays on the host (list: ROCM_HOST_PARTS_2026-10-17.md).
- GFN2 SCF on ROCm is bound by the FP64 eigensolve; `-scf_mixed_precision false` is ignored (G2-13); warp reductions assume wave32 (F-17).

## Result 4: documented defect F-1 is fixed, docs are stale

`TODO.md`, `MULTI_GPU_GAPS.md` (F-1) and Known Issue #35 say ROCm GFN-FF computes no Coulomb term by default. Commit `743cffba` ("Fix ROCm GFN-FF: implicit Coulomb had no HIP consumer") added the consumer (`k_coulomb_dense`, `gfnff_rocm.hip:6303`). The equal energies on 1410 and 6200 atoms (Result 1) would not hold without the Coulomb term. The three entries were not edited in this change.

## Not tested

- Gradient and trajectory parity beyond MD step 0; long MD stability.
- Any CUDA measurement; any other AMD GPU; multi-GPU.
- Charged or open-shell systems, transition metals, solvation on ROCm.
- Per-phase GPU timing of gfn2 (the cause of the 3-6x slowdown on small systems is not established).

## Reproduce

```
release/curcuma -sp a.xyz -method gfnff|gfn2 -gpu none|rocm -no_bmt -threads 10
release/curcuma -md a.xyz -method gfnff -gpu none|rocm -MaxTime 10|40 -dt 1 -no_bmt -threads 10
release/curcuma -opt a.xyz -method gfn2 -gpu none|rocm -no_bmt -threads 10
```
