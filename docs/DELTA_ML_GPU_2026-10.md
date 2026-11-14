# Delta-ML on the GPU: descriptor, forces and fit on the device (Oct 27, 2026)

Status: AI-generated, machine-tested. Human production testing pending.

Branch `feature/delta-ml`, base commit `483e3804`, one RTX 5080 (sm_120, CUDA 13.3), build
`build_deltaml_cuda` (Release, MKL, `-DUSE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=120`). Models: peo8
(58 atoms, 274 labelled frames, M = 1000) and the 90-atom host M1 (353 conformers, M = 2000), both
fitted with `-activelearn -dataset ... -method gfnff -ref_method gfn2`.

## Profile before the port

Per MD step (`-md`, dt 0.5 fs, mean over 501 steps; env-gated scoped timers, removed afterwards;
load average about 14 on 32 cores from two unrelated MD jobs):

| ms per step | peo8 CPU | peo8 old GPU | M1 CPU | M1 old GPU |
|---|---|---|---|---|
| descriptor alone | 0.58 | 0.57 | 1.17 | 1.11 |
| descriptor + Jacobian (host) | 1.09 | 1.11 | 2.71 | 2.62 |
| host-device transfer | - | 0.06 | - | 0.08 |
| kernel on device | - | 0.09 | - | 0.16 |
| W-Jacobian contraction (host) | 0.46 incl. kernel | 0.12 | 2.14 incl. kernel | 0.65 |

Per fit (M1, 282 training structures, M = 2000, old GPU path, total 13.3 s): host design block
(kernel derivative x Jacobian GEMMs) 8.8 s, host descriptor + Jacobian 0.93 s, H2D 0.28 s, DSYRK
0.68 s, evaluation of the training set (CPU prediction) 1.6 s, descriptor pass 0.30 s, solve 0.22 s.

Decision: the host part was 85-90 % of the old GPU step and two thirds of the fit, so both were ported.

## What runs on the device now (`src/core/ml/cuda/ml_model_cuda.cu`)

- Prediction: coordinates and element ids in; neighbour lists, standardised ACSF G2/G4, kernel against
  the resident centres, adjoint W, contraction into per-centre slots and an ordered gather into the
  gradient; energy, indicator and gradient out. No descriptor matrix or Jacobian crosses the bus.
- Fit: same descriptor and kernel, the compact Jacobian on the device, the gradient rows of the design
  block from it, then cuBLAS DSYRK/DGEMV as before. The solve (LDLT) stays on the CPU. The fitted model
  keeps the device context, so the post-fit `evaluate()` runs there too.
- Definitions as in `descriptor.cpp` (cutoff, C1 cutoff function, element-resolved shells, unknown
  elements contribute nothing). The topology anchor acts only on the base method's bond graph, not on the
  descriptor, so the device path needs no anchor logic.
- No floating-point atomics; every sum in fixed order, the gradient gathered over centres in ascending
  index. Integer zeta by repeated multiplication, per-neighbour cutoff values computed once.
- Refusal instead of divergence: grid beyond 32 radial / 8 zeta / 8 lambda, or more than 512 neighbours
  inside the cutoff -> error message and the CPU path. Unavailable plugin or device -> one warning.
- Removed: `GpuAccumulator`, `GpuInference` and their plugin entry points (replaced, no caller left).

## Acceptance

1. FD gate, fully assembled surface (`-method deltaml:<model> -gpu cuda`, GFN-FF on the GPU as well),
   all 3N coordinates, h = 1e-4 Angstrom: max |g - FD| 3.5e-8 Eh/A (peo8), 1.1e-8 Eh/A (M1). The CPU path
   gives 3.5e-8 / 1.1e-8 on the same structures, so the residual is the base surface, not the correction.
2. GPU vs CPU, 25 distinct frames per molecule, base overridden to UFF (no GPU path; base checked
   bit-identical with and without `-gpu`), so differences are the correction alone: max |dE| 8.5e-14 /
   4.3e-14 Eh (correction incl. element baseline up to 76 / 103 Eh), max |dg| 1.5e-15 / 1.0e-14 Eh/A
   (largest correction gradient component 0.04 / 0.10 Eh/A), relative 5.6e-11 / 1.2e-11 on components
   above 1e-6.
3. Fit: held-out (55 / 71 structures) relative energy MAE, gradient RMSE and Spearman of CPU, old-GPU and
   new-GPU fits identical to every printed digit (peo8 2.3662 kcal/mol, 0.006208 Eh/A, 0.9164; M1 0.8329,
   0.003668, 0.4035). Model to model on 5 frames: new GPU fit vs CPU fit 1.3e-8 Eh / 2.6e-8 Eh/A; the old
   GPU fit shows the same 1.4e-8 / 2.7e-8 (summation order of the normal equations). Energy-only fit
   (9-atom fixture): 1.2e-11 Eh.
4. Timing, alternating, n = 3, load average about 1 on 32 cores, 2000 MD steps, `-threads 4`; correction
   = Delta-ML run minus base-only run of the same mode (median / mean, ms per step):

   | | CPU | old GPU (HEAD) | new GPU |
   |---|---|---|---|
   | peo8 correction | 1.27 / 1.28 | 1.09 / 1.08 | **0.17 / 0.16** |
   | peo8 whole step | 1.44 / 1.44 | 1.76 / 1.76 | 0.84 / 0.84 |
   | M1 correction | 3.13 / 3.15 | 2.15 / 2.15 | **0.29 / 0.28** |
   | M1 whole step | 3.43 / 3.44 | 3.22 / 3.22 | 1.35 / 1.35 |
   | fit peo8 (train_seconds) | 2.5 / 2.5 | 2.1 / 2.1 | 0.8 / 0.8 |
   | fit M1 (train_seconds) | 11.7 / 11.7 | 7.5 / 7.5 | 2.1 / 2.1 |

   The whole GPU step is now dominated by GFN-FF on the GPU (0.67 / 1.06 ms against 0.16 / 0.29 ms on
   the CPU at these sizes); `-gpu` is one global switch, so that placement cannot be split.
5. `ctest -R 'deltaml|cli_simplemd'`: 26/26 in `build_deltaml_cuda`, 26/26 in `build_deltaml` (CPU-only;
   `cli_deltaml_06_gpu_path` skips there with exit 77 after the "unavailable" warning).

`cli_deltaml_06_gpu_path` repeats 1-3 in small (9 atoms, 3 frames, 5 FD coordinates, fit agreement).

## Not tested

ROCm and Vulkan (no entry points: they report unavailable), multi-GPU, charged or open-shell systems,
elements beyond C/H/N/O, M above 2000, more than 90 atoms, more than 512 neighbours (refused, path not
exercised), concurrent use of one model by several threads (calls are serialised by a mutex, not measured).
