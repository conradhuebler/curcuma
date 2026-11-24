# ROCm parity SPL: GFN-FF Coulomb-tile and D4 multi-GPU split (Nov 25, 2026)

Follow-up of A6 ([ROCM_PARITY_A6.md](ROCM_PARITY_A6.md)) in [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md).
🤖 AI-generated, ⚙️ machine-tested on ONE GPU. Machine: one RX 9060 XT (gfx1200), ROCm 7.2.4. **No run with more than one
physical device was possible; every real multi-device path (peer copies between cards, cross-device events, helper
streams on another device) is compiled but untested.**

## 1. What was ported

Same opt-in as the EEQ split of A6: `-gfnff.rocm_multi_gpu_split true` with `-gpu_split_devices <list>` (trait
`has_multi_gpu_split_opt_in`); without it an explicit list stays refused with a counted fallback (`-gpu_strict` exit 3).
Thresholds as CUDA: Coulomb from `gpu_split_min_atoms` atoms (default 1000), D4 in addition from 200000 pairs;
`CURCUMA_GFNFF_DISP_SPLIT=0` switches the D4 split off; both are off with a periodic cell.

| ROCm (`rocm/gfnff_rocm.hip`) | CUDA reference (`cuda/ff_workspace_gpu.cu`, `cuda/gfnff_kernels.cu`) |
|---|---|
| `FFWorkspaceHipImpl::CoulSplitPartHip`, `DispSplitPartHip`, `releaseCoulSplit`, `releaseDispSplit` | `FFWorkspaceGPUImpl::CoulSplitPart`, `DispSplitPart`, release members |
| `FFWorkspaceHip::setCoulombSplitDevices` (was a no-op), `coulombSplitStatus`, `dispSplitStatus`, `prepareCoulombSplit` | same members |
| `launchCoulombSplitHip`, `launchCoulombTilesRowsHip` | `launchCoulombSplit` |
| `prepareDispSplitHip`, `launchDispSplitHelpersHip`, `appendDispSplitPartsHip` | `prepareDispSplit`, `launchDispSplitHelpers`, `appendDispSplitParts` |
| new kernels `k_dispersion_gather_uv_rows`, `k_split_edge_params_hip`, `k_add_split_parts_hip` | `k_dispersion` on a pair range, `k_add_coulomb_parts`, `k_add_disp_parts` |
| call sites: implicit Coulomb in `launchChargeDependentAndFinish`, fused gather in `launchDispersionGatherHip`, energy-only `k_dispersion` in `prepareAndLaunchChargeIndependent` | `k_coulomb_tiles` / dispersion call sites |

How the work is divided (ROCm-specific where the ROCm kernels differ from CUDA):
- Coulomb: the tile ROWS of the 2-D `k_coulomb_tiles_hip` grid are split into ranges of equal pair count. A range is
  launched with every atom pointer shifted by `row0 * TI` and the tile shape of the full N, so each block is the block
  of the unsplit launch with the same pairs in the same order; the former kernel is reused unchanged.
- D4, gradient: ROCm uses the per-atom gather `k_dispersion_gather_uv` (C6 from `d_disp_uv` per edge), not the CUDA
  pair kernel. Helpers take contiguous ATOM ranges of equal gather-CSR edge count, with a row-local CSR and per-edge
  pair parameters (built on the own device, copied once per pair list, `disp_gen`); per step coordinates, `d_gw` and
  `d_disp_uv` are copied. Only the fused path is split; the unfused gather (`CURCUMA_ROCM_DISP_FUSED=0`, tabulated
  dc6dcn, no per-step C6 refresh) stays on the own device for that call and is counted in the status line.
- D4, energy only: contiguous PAIR ranges of `k_dispersion` as on CUDA; the per-pair C6 slice is copied per call.
- Parts come back into buffers on the own device; `k_add_split_parts_hip` adds them (integer sums with
  `-gfnff.rocm_fixed_point`, the default; FP64 sums otherwise). The fixed-point switch is uploaded to every helper device.
- No split while a stream records a graph (`ffStreamCapturingHip`).
- Status: `GFN-FF multi-GPU: Coulomb tiles on devices [...]` and `dispersion on devices [...]` at verbosity 2 (the wrapper
  prints the D4 line for the opt-in backend only, CUDA output unchanged); a single point prints both lines when the
  workspace is destroyed (the wrapper reports after the second calculation only). Set-up failures give
  "unavailable: ..." and a counted fallback as in A6.

Test hook `CURCUMA_ROCM_SPLIT_SAME_DEVICE`: `1` keeps helper entries that name the own device; such a helper has no
stream of its own, its partition runs on the calling stream before the own partition (partitions one after the other).
`streams` gives a same-device helper its own stream and event, the code path of a real helper device.

## 2. Measured (RX 9060 XT, every GPU run inside `flock /tmp/claude-1000/rocm-gpu.lock`)

| check | structure, n | result |
|---|---|---|
| ISA of `gfnff_rocm.hip.o` (`a2/isa.py`), base vs new | all kernels | 142/142 former functions identical, 3 new kernels |
| hipcc warnings in `gfnff_rocm.hip` | base vs new build log | 210 / 210 |
| no switch, base vs new: md5 of `-sp -dump_gradient` | caffeine, triose, complex, polymer; polymer with refused `-gpu_split_devices 0,1` | 5/5 identical |
| no switch, base vs new: md5 of MD diagnostics (4 records) | polymer_2x (7320 atoms) | identical |
| second run of the new build | the 6 runs above | 6/6 md5-identical |
| Coulomb + D4 split only (EEQ split off), fixed point: 2 and 4 participants vs single | polymer_2x MD, 4 records (charges, every energy term, gradient norm) | bitwise identical (2/2) |
| same, 2, 4 and 8 participants, md5 of the diagnostics file | polymer_2x MD | identical to `none` (3/3) |
| same, FP64 atomics (`-gfnff.rocm_fixed_point false`), 2 and 4 participants vs single | polymer_2x MD | max dq 1.3e-12 e, max dE term 3.4e-13 Eh, max dgrad-norm 1.3e-13 |
| helper on its own stream (`streams`) vs sequential, 2 participants | polymer_2x MD | bitwise identical |
| `-sp -dump_gradient`, 3 participants vs single, fixed point | polymer (1410 atoms) | dE 0, max dG 0 (bitwise) |
| energy-only `-sp` (pair-kernel split, status "1 energy calls split"), 3 participants vs single, both accumulations | polymer | Dispersion, Coulomb terms and total equal to the 10 printed decimals |
| EEQ + Coulomb + D4 split (A6 checks, now with all three split) | polymer_2x MD | split(2) vs single 1.6e-11 e / 8.0e-13 Eh / 9.7e-13; split(2) == split(3) bitwise |
| FD of the assembled surface, split active (3 participants), h 2e-4 / 1e-4 A, central | polymer, atoms 0, 100, 700, 1400 (12 components) | 2.0e-8 / 6.6e-9 Eh/A, Richardson 4.9e-9 Eh/A |
| ctest `gfnff_rocm_split_one_device` (extended, see below) | | passed |
| `ctest -j1 -R "gfnff.*rocm\|rocm.*gfnff"` (28 tests: `gfnff_rocm_*`, structinv GFN-FF, `cli_gfnff_gpu_*_rocm`) | final build | 28/28 passed (`gfnff_rocm_charge_groups` needs numpy on `PYTHONPATH`) |

Stall found with the `streams` mode: 4 participants (6 extra streams on the device) made the third MD step take
58.6 s instead of 0.08 s (also with `GPU_MAX_HW_QUEUES=8`: 65.5 s); `timeout` (SIGTERM) did not end the process, it
needed SIGKILL. Not root-caused; with 2 participants and own streams no stall. This is why the default hook mode runs
same-device partitions on the calling stream. Whether helper streams on a real second device can stall the same way
(several helper streams per device: Coulomb, D4, EEQ) is open.

Extended ctest `gfnff_rocm_split_one_device` (`test_cases/scripts/check_gfnff_rocm_split.py`, now also given
`polymer.xyz`): the A6 checks plus Coulomb/D4-only split with 2 and 4 participants bitwise equal to single (fixed
point), within 1e-8 with FP64 atomics, the `streams` helper bitwise equal to the sequential one, and a component-wise
`-sp` gradient comparison with 3 participants. Full script 74 s on this machine.

## 3. Not tested

- Any run with two or more physical GPUs: peer access, `hipMemcpyPeerAsync` between cards, cross-device event waits,
  helper stream priorities, the fixed-point flag upload to another device, allocation failures on a helper, timing.
- A range overflow or non-finite value on a helper device: `d_gfnff_fix_bad` of a helper is not read back, so such a
  step is not turned into NaN as on the own device.
- Devices without peer access, wave64 devices (CDNA/MI), systems with a cell (refused), other structures than the
  polymer pair (polymer, polymer_2x), the unfused gather path with the split configured (only counted, never measured).
- No timing claim: on one device the partitions run one after the other, so no speed-up is possible or measured.

## 4. Not implemented relative to CUDA

- CUDA splits the D4 pairs by a cost scan (`disp_pair_cost_scan`) and computes dC6/dCN on the helper; ROCm splits the
  gather rows by edge count and copies `d_disp_uv` per step. The energy-only pair ranges are equal-sized (no cost scan).
- The unfused gather and the tabulated dc6dcn path are not split (CUDA splits its tabulated path).
- CUDA's "switched off after start" state of the D4 split has no counterpart: ROCm helpers hold no per-step state, so the
  split is decided per call.

## 5. Proposed lines for the central files

- `scripts/gpu_parity_known.json`, trait `has_multi_gpu_split_opt_in`: "ROCm: EEQ, Coulomb-tile and D4 split compiled in, opt-in -gfnff.rocm_multi_gpu_split, never run with >1 device (ROCM_PARITY_A6.md, ROCM_PARITY_SPL.md)"; `has_multi_gpu_split` reason: "ROCm split only opt-in, unverified on >1 device".
- `docs/GPU_PARITY_CUDA_ROCM.md`, GFN-FF line 30: "`has_multi_gpu_split` (ROCm false; EEQ, Coulomb and D4 split opt-in `-gfnff.rocm_multi_gpu_split`, checked on one GPU only)".
- `docs/MULTI_GPU.md`, ROCm block: replace "Not ported: the Coulomb-tile and D4-pair splits" by "Coulomb tile rows and D4 (gather rows, energy-only pair ranges) split behind the same switch, bitwise equal to one device with fixed point on one GPU with up to 8 participants ([ROCM_PARITY_SPL.md](ROCM_PARITY_SPL.md))".
- `TODO.md`: "ROCm GFN-FF multi-GPU split (EEQ, Coulomb, D4; `-gfnff.rocm_multi_gpu_split`): run on two real GPUs (split vs single, md5 of two runs, timing, stall check with several helper streams per device), read back `d_gfnff_fix_bad` of helpers (ROCM_PARITY_SPL.md)."
- `AIChangelog.md`: "2026-11-25 ROCm SPL: GFN-FF Coulomb-tile and D4 multi-GPU split ported behind -gfnff.rocm_multi_gpu_split (one-GPU checks with up to 8 participants, bitwise with fixed point; docs/ROCM_PARITY_SPL.md)."
- `src/core/energy_calculators/ff_methods/CLAUDE.md`, GPU bullet: "ROCm multi-GPU: EEQ, Coulomb and D4 split opt-in `-gfnff.rocm_multi_gpu_split` (never run on two devices, ROCM_PARITY_A6.md, ROCM_PARITY_SPL.md)."

## 6. Plugin ABI

No change needed: `cuda/ff_workspace_gpu.h` (the declaration shared by both plugins) is unchanged; the ROCm build now
defines members it already declared (`prepareCoulombSplit`, `dispSplitStatus`). New state lives in the pimpl
`FFWorkspaceHipImpl` inside `libcurcuma_rocm.so`. `ComputationalMethod`, the plugin C entry points and shared class
layouts are unchanged. The wrapper change in `gfnff_gpu_method_impl.h` is a message and one extra status line under
`if constexpr (!Backend::has_multi_gpu_split)` (not compiled for CUDA; CUDA was not built here).

## 7. Other changed files

`qm_methods/gfnff_gpu_method_impl.h` (messages, D4 status line), `ff_methods/gfnff.h` (help texts of
`gpu_split_devices`, `rocm_multi_gpu_split`), `qm_methods/gfnff_hip_method.h` (comment), `test_cases/CMakeLists.txt`
(second argument of `gfnff_rocm_split_one_device`), `test_cases/scripts/check_gfnff_rocm_split.py`.
