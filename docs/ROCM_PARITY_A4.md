# ROCm parity A4: spin polarisation on `-gpu rocm` (Nov 24, 2026)

Status: 🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT (gfx1200). Human production testing pending.
Item A4 of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md). CUDA reference: host-driven device loop of
[SPIN_POLARISATION.md](SPIN_POLARISATION.md), section "GPU (CUDA)".

## What was ported

| File | Function | CUDA counterpart |
|---|---|---|
| `src/core/energy_calculators/qm_methods/xtb_hip_method.cpp` | `HipScfBackend::supportsSpinPolarization`, `solveSpin`, `densitySpin`, `finalizeSpin`, `spinWeightedDensity` | `XtbGpuScfBackend` in `xtb_gpu_method.cpp` |
| `src/core/energy_calculators/qm_methods/rocm/xtb_hip_context.{h,hip}` | `XtbHipContext::residentSpinKeepAlpha`, `residentDensitySpin`, `residentFinalizeSpin`, `spinWeightedDensity`; `Impl::ensureSpin/freeSpin/resSpinKeepAlpha/spinChannel/resDensitySpin/resFinalizeSpin/spinW`; kernels `k_spin_combine`, `k_spin_w_shift` (verbatim) | `cuda/xtb_gpu_context.cu` functions and kernels of the same names |
| same `.hip`, `Impl::gradient` | one-shot `w_spin_override`: W' from `spinW` replaces the W build (also the WP7 `w_from_goo` branch); flag consumed at function entry | `computeGradient` `w_spin_override` |
| `src/core/energy_calculators/qm_methods/xtb_native.cpp` | error texts only (ROCm no longer refused; Vulkan still is) | - |
| `src/core/energy_calculators/qm_methods/native_xtb_method.h` | PARAM help of `spin_polarization` (cuda\|rocm) | - |
| `test_cases/cli/sqm/23_spin_polarization_rocm/`, `test_cases/cli/CMakeLists.txt` | ctest `cli_sqm_23_spin_polarization_rocm` (USE_ROCM only, skip 77 without device) | `cli_sqm_22_spin_polarization_cuda` |

The core (`xtb_native.cpp`) needed no logic change: it already keeps the device potential build, the device multipole potential, both resident loops,
partial diagonalisation and the device gradient potentials off with spin polarisation, and only `supportsSpinPolarization()` gated ROCm.
`parameters/spin_constants.hpp` is unchanged.

## Measured

Build `release/` (USE_ROCM, MKL, AVX2, gfx1200), every GPU run under `flock /tmp/claude-1000/rocm-gpu.lock`, `-threads 1`, `-scf_threshold 1e-10`,
`-gradient -dump_gradient` (energy 15 decimals, gradient 17 digits). Script: scratchpad `a4/gates.py` (not in the repo); structures from `test_cases/structures`.

- GPU against CPU of the same build, 10 cases (O2 `dioxygen.bent-triplet` uhf 2 and 4 GFN2, uhf 2 GFN1; OH `hydroxyl` uhf 1 GFN2, GFN1; H2 `dihydrogen.r0.741` uhf 2 GFN2;
  Pt7 `pt-cluster-7` 2000 K uhf 4 and 6, GFN2 and GFN1): FP64 (`-scf_mixed_precision false`) |dE| <= 2.5e-14 Eh, max |dg| <= 1.4e-12 Eh/A; mixed precision
  |dE| <= 4.7e-14 Eh, max |dg| <= 7.6e-12 Eh/A. Every GPU run printed the spin-branch log line and ran with `-gpu_strict true`; the device gradient ran (gradient tail log line).
- O2 uhf 2 with `-solvent water` (ALPB), GFN2 and GFN1: |dE| <= 3.6e-15 Eh, max |dg| <= 2.9e-15 Eh/A (n = 2).
- FD gate of the assembled GPU surface (central differences of GPU `-sp` energies, FP64): O2 uhf 2 GFN2 5.3e-7 (h = 5e-4 A), 1.3e-7 Eh/A (h = 2.5e-4 A); GFN1
  6.1e-7 / 1.5e-7 (ratio 4.0 for h/2: truncation error); Pt7 uhf 4 GFN2 3.3e-8, GFN1 6.3e-8 Eh/A (h = 5e-4 A, 21 components each).
- Determinism: two runs md5-identical (`-dump_gradient`, default mixed precision), O2 GFN2, Pt7 GFN2 and GFN1 (n = 3).
- Without spin polarisation, new binary against the binary built from the worktree HEAD + base patch before this change: 10 of 10 `-dump_gradient` files
  bit-identical (caffeine GFN2/GFN1 GPU mixed, caffeine GFN2 GPU FP64, O2 spin 2 GFN2 GPU, OH spin 1 GFN1 GPU, Pt7 spin 8 GFN2 GPU, Pt7 spin 4 GFN1 GPU,
  O2 spin 2 with `-xtb.spin_polarization false` GPU, caffeine GFN2 CPU, O2 spin 2 spin-polarised CPU).
- ctest: `cli_sqm_23_spin_polarization_rocm` 27/27 checks pass; the same script against the base binary fails 18 checks (negative control);
  `cli_sqm_16_spin_polarization` and the 32 `sqm_val_rocm_*` gates pass.

## Not tested

Systems above 13 atoms, MD, `-opt`, `-hessian`, d-shell main-group elements (S, P, Cl) with spin polarisation on ROCm, charged systems, solvation models other than ALPB,
`-gpu_sparse_integrals`, multi-GPU, any gfx other than gfx1200, timing (no measurement made). The full `ctest -j1` was not run (only the targets
`curcuma` and `curcuma_rocm` were built).

## Not implemented relative to CUDA

- Nothing on the CUDA spin path is missing. ROCm-only paths that are not used with spin polarisation (automatic, GPU stays in use): `-xtb.rocm_device_potential`,
  `-xtb.rocm_device_multipole`, `-xtb.rocm_resident_loop` (GFN1/GFN2), `-xtb.rocm_grad_resident`, pseudo-diagonalisation; post-SCF D4 runs on the host
  (no device D4 reference tables on the host-driven loop); `-xtb.rocm_density_syrk` is not applied to the two spin densities (GEMM).
- Still refused with an error (tested in `cli_sqm_23`): `-scf_gpu_partial_diag`, `-scf_precond frag`, `-large_system_mode fragments`; `-gpu vulkan` stays refused.

## Proposed lines (for the coordinator)

- AIChangelog: `Nov 24, 2026: ROCm spin polarisation (-xtb.spin_polarization, -gpu rocm): host-driven device loop ported from CUDA (HipScfBackend, XtbHipContext spin functions); GPU = CPU <= 2.5e-14 Eh / 1.4e-12 Eh/A (FP64, 10 cases), ctest cli_sqm_23 (docs/ROCM_PARITY_A4.md)`
- TODO.md: remove the ROCm part of any "spin polarisation ROCm/Vulkan" item; Vulkan remains open.
- GPU_PARITY_CUDA_ROCM.md: row spin polarisation: CUDA yes, ROCm yes (host-driven loop, Nov 24, 2026), Vulkan no.
- Root CLAUDE.md "Spin polarisation" bullet: "`-gpu cuda` on the host-driven device loop" -> "`-gpu cuda|rocm` on the host-driven device loop (ROCm Nov 24, 2026, one RX 9060 XT)", "ROCm/Vulkan and large-system modes are errors" -> "Vulkan and large-system modes are errors".
- `scripts/gpu_parity_known.json`: entries for `supportsSpinPolarization`, `k_spin_combine`, `k_spin_w_shift`, `residentSpinKeepAlpha`, `residentDensitySpin`, `residentFinalizeSpin`, `spinWeightedDensity` (if listed as CUDA-only) can be removed; not checked whether the script lists them.

## Plugin ABI

No change needed: `GpuScfBackend` already had the five spin virtuals (CUDA), `HipScfBackend` lives inside the plugin, and `XtbHipContext` is compiled
only into the ROCm plugin; no class layout shared between `curcuma` and the plugin changed. The base binary and the new plugin were not mixed in any test.

## Remarks

- The CUDA test `cli_sqm_22` greps `spin_polarization.*<flag>` over the whole log, which the echoed command line already satisfies; `cli_sqm_23` requires the `[ERROR]` line. Same weakness left in `cli_sqm_22` (not changed, out of scope).
- The worktree has an untracked symlink `external -> /home/conrad/src/curcuma/external` (needed for the build, not part of the change).
