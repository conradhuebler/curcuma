# ROCm parity A6: multi-GPU on ROCm (Nov 24, 2026)

Item A6 of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md). 🤖 AI-generated, ⚙️ machine-tested on ONE GPU.
Machine: one RX 9060 XT (gfx1200), ROCm 7.2.4, rocSOLVER 0.7. **No run with more than one physical device was possible;
every real multi-device path below is untested.**

## 1. GFN-FF: projected-PCG EEQ split (ported, opt-in)

What was ported (port of `cuda/eeq_solver_gpu.cu`, F-3 of MULTI_GPU_GAPS.md):

| ROCm (file, function) | CUDA reference |
|---|---|
| `rocm/eeq_solver_hip.hiph`: `EEQSolverHipImpl::SplitPart`, `releaseSplit` | `EEQSolverGPUImpl::SplitPart`, `releaseSplit` |
| `EEQSolverHip::setSplitDevices`, `setSplitMixedPrecision`, `splitStatus` | same members |
| `splitWantedHip`, `splitPrepareHip`, `splitBuildHip`, `splitMatvecHip`, `splitMatvecF32Hip` | `splitWanted`, `splitPrepare`, `splitBuild`, `splitMatvec`, `splitMatvecF32` |
| kernels `k_eeq_build_cols_hip`, `k_split_coldot_f64_hip`, `k_split_coldot_f32_hip`, `k_eeq_diag_inv_hip` | `k_eeq_build_cols`, `k_split_coldot_f64/f32`, `k_eeq_diag_inv` |
| `solveWithDeviceRHSAndGPUProjectedPCG`: split build, diagonal, dense matrix only for the local-inverse rebuild; `eeqMatvecHip`/`eeqMatvecF32Hip` dispatch to the split | `solveWithDeviceRHSAndGPUProjectedPCG`, `eeqMatvec` |
| `rocm/gfnff_rocm.hip`: `FFWorkspaceHip::setCoulombSplitDevices` (no-op), `coulombSplitStatus` | Coulomb-tile and D4-pair split (not ported) |
| `gfnff_gpu_method_impl.h`: split block compiled for `has_multi_gpu_split \|\| has_multi_gpu_split_opt_in`, opt-in gate | unchanged for CUDA (`has_multi_gpu_split_opt_in = false`) |

Only HIP runtime peer copies (`hipMemcpyPeerAsync`, `hipDeviceEnablePeerAccess`); no RCCL (not installed here, not
needed: CUDA uses peer copies too).

Switches:
- `-gfnff.rocm_multi_gpu_split` (PARAM in `gfnff.h`, default false, wrapper fallback false). Without it an explicit
  `-gpu_split_devices <list>|all` naming more than one device is refused: warning, counted GPU fallback
  (`-gpu_strict true` exits 3); `auto` stays silent as before. With it the EEQ split runs and the first calculation
  prints "unverified, never run with more than one device".
- `has_multi_gpu_split` stays **false** on ROCm; the new trait `has_multi_gpu_split_opt_in` is true on ROCm, false on CUDA.
- Test hook `CURCUMA_ROCM_SPLIT_SAME_DEVICE=1`: helper entries may name the calculation's own device ("0,0" = two
  participants on device 0, "0,0,0" = three), device-local copies instead of peer copies.

Measured (RX 9060 XT, `-md -maxtime 4 -dt 1 -thermostat none -seed 7 -md_diagnostics true`, 4 MD records each,
charges, every energy term, gradient norm; all runs inside `flock /tmp/claude-1000/rocm-gpu.lock`):

| structure | comparison | max abs dq (e) | max abs dE term (Eh) | max abs dgrad |
|---|---|---|---|---|
| polymer_2x (7320 atoms, 1502 frag., localinv device build, 4 PPCG solves of 16-18 it.) | base build vs new build, no switch, n = 3 runs | md5 identical (bitwise) | | |
| polymer_2x | split 2 participants (FP32 matvec, default) vs single device (FP32 default) | 1.6e-11 | 8.0e-13 | 9.7e-13 |
| polymer_2x | split 2 vs split 3 participants (FP32), and split 2 run twice | 0 (bitwise) | 0 | 0 |
| polymer_2x | split 2 FP64 (`-gfnff.gpu_eeq_mixed_precision false`) vs single FP64 (`CURCUMA_ROCM_EEQ_FP32=0`) | 8.4e-15 | 3.4e-13 | 1.5e-14 |
| polymer_2x | split 2 vs split 3 participants (FP64) | 0 (bitwise) | 0 | 0 |
| polymer (1410 atoms, 1 fragment, `-gfnff.gpu_split_eeq_min_atoms 0`) | split 2 vs single | 5.4e-12 | 6.0e-13 | 5.2e-13 |
| polymer | split 2 vs split 3 | 0 (bitwise) | 0 | 0 |

For scale: single device FP32 vs FP64 matvec differ by 4.3e-12 e / 1.6e-12 Eh on polymer_2x; the CUDA test limits
(`check_gfnff_gpu_split.py`) are 1e-9 e, 1e-9 Eh, 1e-8.

New ctest `gfnff_rocm_split_one_device` (`test_cases/scripts/check_gfnff_rocm_split.py`, labels `gfnff_gpu;gpu;rocm`):
split(2) vs single within the CUDA limits, split(2) == split(3) bitwise, the status line present, and the refusal without
the opt-in (`-gpu_split_devices 0,1`, exit 3). Passed.

Not tested:
- any run over two or more physical GPUs: peer access, P2P copy bandwidth, cross-device event waits, the stream
  priority on a helper, helper allocation failures, `-gpu_devices` pool interaction, multi-GPU timing;
- devices without peer access (`hipDeviceCanAccessPeer` false; how `hipMemcpyPeerAsync` behaves then was not checked);
- wave64 devices (CDNA/MI): the coldot kernels use a fixed 256-thread block and LDS tree, independent of the wave width
  by construction, but never run there;
- a split that fails during set-up or build (status "unavailable: ...", counted fallback) was not injected.

Not implemented relative to CUDA:
- the Coulomb-tile split and the D4-pair split (`prepareCoulombSplit`, `prepareDispSplit`): both terms stay on the
  calculation device; the status line says so;
- the CUDA mixed-precision scheme (defect-correction rounds with `gpu_eeq_mixed_rho`): the ROCm split keeps the ROCm
  reliable-update FP32 scheme of `runProjectedPCGHip`; `gpu_eeq_mixed_rho` is stored and not used;
- no memory relief on the calculation device: `d_A` is allocated for max_natoms^2 in the `EEQSolverHip` constructor
  and is filled when the local inverse is rebuilt on the device (CUDA does the same rebuild, but allocates on demand).

## 2. GFN1/GFN2: multi-GPU eigensolver and density split (refused)

Assessment (ROCm 7.2.4 on this machine, headers read, no measurement):
- rocSOLVER offers `syev`, `syevd`, `syevdx`, `syevj`, `syevdj`, `syevx`, `sygvd` and their batched / strided-batched
  forms: all single-device. There is no counterpart of cusolverMg or cuSOLVERMp; hipSOLVER wraps the same single-device
  API; rocBLAS has no multi-device (XT) BLAS.
- Distributed eigensolvers with AMD GPU support exist outside ROCm (SLATE, ELPA, MAGMA multi-GPU); none is installed or
  a dependency of curcuma, and any of them would be a new build dependency whose multi-device path cannot be run here.
- Spectrum slicing with `rocsolver_dsyevdx` (each device the full matrix, its own index range) needs no new
  library, but every device repeats the tridiagonal reduction, the dominant O(n^3) part, so it cannot be faster than one
  device; eigenvectors of clustered eigenvalues from different devices are not guaranteed orthogonal.
- An own block-cyclic distributed solver (reduction to tridiagonal form over devices, back-transformation) is a
  separate project and could not be validated with one device.

Decision: refusal with a message, no port. `XtbHipComputationalMethod::refuseMultiGpuSplit` (`xtb_hip_method.cpp`):
`-gpu_eigensolver_devices` or `-gpu_density_devices` naming two or more devices (or `all`/`auto` with more than one
visible device, `XtbHipContext::deviceCount()`) gives a warning "... not available on ROCm ... uses one device" and a
counted GPU fallback (`-gpu_strict true` exits 3). Before, both keys were ignored on ROCm without a message (checked
with the base build: `-gpu_eigensolver_devices 0,1 -gpu_strict true`, exit 0, no message). Without the keys nothing
changes; CUDA's default split above 4000 basis functions has no ROCm counterpart and stays unmentioned.
The pattern-density split (`-gpu_density_devices`) needs only peer copies and a daxpy and could be ported like the EEQ
split; it is refused here as well (not part of A6, would also be untestable beyond one device).

New ctest `xtb_rocm_multi_gpu_refusal` (`test_cases/scripts/check_xtb_rocm_multi_gpu_refusal.py`, caffeine gfn2):
both keys with "0,1" exit 3 with the message, `none`/`none` exits 0. Passed.

## 3. Regression

- Gate "without the new switch bitwise equal to the unchanged build": polymer_2x GFN-FF MD diagnostics md5-identical
  between the base build (worktree HEAD + base.patch) and the new build, 2 + 3 runs, also with the default `auto`.
- `ctest -L rocm -j1` (117 tests, inside the GPU lock): 112 passed, 5 not run (test executables not built); after building
  `test_gfnff_energy_only_geometry` and `test_hbond_energies_gpu` the 5 passed. The full `ctest -j1` was not run (wave C).
- No finite-difference check was run: the split changes only the summation order of the EEQ matvec, and energies and
  gradients equal the single-device run to 1e-12 (table above), whose surface the existing ROCm gates cover.
- No new compiler warnings (warning counts of `eeq_solver_hip.hiph`, `xtb_hip_method.cpp`, `gfnff_gpu_method_impl.h`
  equal in the base and the new build log).

## 4. Proposed lines for the central files

- `scripts/gpu_parity_known.json`, traits: `"has_multi_gpu_split_opt_in": "ROCm: EEQ split compiled in, opt-in -gfnff.rocm_multi_gpu_split, never run with >1 device (ROCM_PARITY_A6.md); CUDA split is released"`; keep `has_multi_gpu_split` with the reason "ROCm split only opt-in, unverified on >1 device; Coulomb/D4 split not ported (ROCM_PARITY_A6.md)".
- `docs/GPU_PARITY_CUDA_ROCM.md` table A, row "Multi-GPU eigensolver (cuSOLVERMg/Mp)": ROCm "refused with a counted fallback (`-gpu_eigensolver_devices`/`-gpu_density_devices` with >1 device), no rocSOLVER multi-device eigensolver (ROCM_PARITY_A6.md)". GFN-FF line 30: "`has_multi_gpu_split` (ROCm false; EEQ split opt-in `-gfnff.rocm_multi_gpu_split`, checked on one GPU only)".
- `TODO.md`: "ROCm GFN-FF EEQ split (`-gfnff.rocm_multi_gpu_split`): run on two real GPUs (split vs single 1e-9, md5 of two runs, timing), then decide the trait; Coulomb/D4 split not ported (ROCM_PARITY_A6.md)."
- `AIChangelog.md`: "2026-11-24 ROCm A6: GFN-FF projected-PCG EEQ split ported to ROCm behind -gfnff.rocm_multi_gpu_split (one-GPU checks only, CURCUMA_ROCM_SPLIT_SAME_DEVICE); GFN1/GFN2 -gpu_eigensolver_devices/-gpu_density_devices refused on ROCm with a counted fallback (docs/ROCM_PARITY_A6.md)."
- `src/core/energy_calculators/ff_methods/CLAUDE.md`, GPU bullet: "ROCm multi-GPU: EEQ split only, opt-in `-gfnff.rocm_multi_gpu_split` (unverified on >1 device, ROCM_PARITY_A6.md)."
- `docs/MULTI_GPU_GAPS.md` G2-1: "ROCm refuses an explicit eigensolver/density split with a counted fallback (Nov 24, 2026)".

## 5. Plugin ABI

No change needed: `ComputationalMethod`, the plugin C entry points and the class layouts shared across the plugin
boundary are unchanged. `EEQSolverHipImpl` (pimpl, plugin-internal) and the new static `XtbHipContext::deviceCount()`
live inside `libcurcuma_rocm.so`; `FFWorkspaceGPU`'s declaration is unchanged (two members of it are now defined on ROCm).
The CUDA trait struct gained a `static constexpr` member (no layout change); CUDA was not compiled here.
