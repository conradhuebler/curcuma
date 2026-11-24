# ROCm parity A1: GFN-FF charge groups and `cn_zero_guard` on the device (Nov 24, 2026)

Status: 🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT (gfx1200, ROCm). Human production testing pending.
Plan: [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md), item A1. CUDA reference: commit `9e8a609c`.

## What was ported

- No kernel or solver code was missing: the ROCm side of `9e8a609c` already carried the `cn_zero_guard` flag of
  `k_coulomb_postprocess` (`rocm/gfnff_rocm.hip`, launched from `FFWorkspaceHip::launchChargeDependentAndFinish` with
  `m_cn_zero_guard`), `GFNFFGpuMethodImpl::initGPUWorkspace` sets it for both backends (`cnZeroGuardEffective()`), and the
  groups reach the ROCm device EEQ through the shared `prepareEEQParametersForGPU` (fraglist -> `uploadFragmentTopology`,
  qfrag -> `rhs_constraints`). The ROCm solvers take the constraint charges from there (`eeqSmallSolveHip`: `qfrag`;
  dense Cholesky + host Schur; projected PCG `frag_delta` from `rhs_constraints`).
- Change: trait `GFNFFRocmBackend::has_charge_groups_cn_guard` false -> true (`qm_methods/gfnff_hip_method.h`), so the
  wrapper no longer requests the CPU method for `-gfnff.charge_groups` / `-gfnff.cn_zero_guard` on `-gpu rocm`.
  Comments updated in `gfnff_gpu_method_impl.h` and `rocm/gfnff_rocm.hip` (comment only, same line count).
- Test: `test_cases/scripts/check_gpu_charge_groups.py` takes an optional backend argument (`cuda` default, `rocm`);
  ctest `gfnff_rocm_charge_groups` registered under `USE_ROCM` (`test_cases/CMakeLists.txt`).
- Docs brought in line: `ff_methods/CLAUDE.md`, PARAM help of `charge_groups` and `cn_zero_guard` (`gfnff.h`),
  `docs/GFNFF_STATUS.md`, `docs/POLYMERBUILD_CHARGE.md`, `docs/KNOWN_ISSUES_ARCHIVE.md` (#43 switch table and text).

## What was measured (gfx1200, all GPU runs `-gpu rocm -gpu_strict true`, GPU against CPU of the same build)

Gate scripts (pure Python, numpy is not installed on this machine): scratchpad `a1_gates.py`, `a1_extra.py`; each GPU
call under `flock /tmp/claude-1000/rocm-gpu.lock`. n = 1 per case.

| Gate | Structure | Result |
|---|---|---|
| Baseline refuses | caffeine, `cn_zero_guard true`, unchanged build | exit 3 (`-gpu_strict`, counted fallback) |
| Guard, gradient | Cl- ... H2O (H-Cl 3.2 A, CN(Cl) = 0, Cl first, `frag_charge_model reference`, q = -1) | GPU vs CPU 6.4e-15 (off), 5.3e-15 (on) Eh/A; switch changes GPU gradient by 3.56e-4 Eh/A; energy identical (-1.30732893 Eh) |
| FD of the assembled surface | same, switch on, 12 components | max \|FD - analytic\| 5.05e-6 Eh/A (h = 2e-3 A); switch off: Cl 3.57e-4 Eh/A (reference defect visible) |
| Guard is a no-op at CN > 0 | Cl- ... H2O 2.2 A, water dimer | GPU gradient bit-identical on/off; GPU vs CPU 4.6e-14 / 3.9e-15 Eh/A |
| 200 MD steps (CSVR, seed 303, dt 0.5 fs) | (AB)2 PDMAEMA chain (108 atoms, charge 2, polymerbuild group file / `auto` / file + guard), NH4+ Cl- `ionpairs` | charges GPU vs CPU first 10 frames <= 8.7e-15 e, last frame <= 1.9e-13 e; mean Epot <= 1.2e-13 Eh; group sums <= 9.4e-15; total charge exact; EEQ line `[GPU ROCm, WP2 dense Cholesky + CPU Schur]`, nfrag 3 / 2 |
| Second run | the four 200-step GPU runs | trajectory, diag, final.json, md.csv identical (md.csv without the wall-clock columns `remaining`, `dt`) |
| FD with groups | chain + group file (+ guard), 12 random components | 2.76e-6 Eh/A; sp GPU vs CPU dE 0, dG 3.9e-14 Eh/A |
| Small vs dense device EEQ | chain + group file, 50 steps, `CURCUMA_ROCM_EEQ_SMALL=0` vs default | default (one-workgroup solve) vs dense 5.7e-15 e; dense vs CPU 7.4e-15 e; group sums 6.2e-15 |
| Projected PCG | chain + group file, `-gfnff.solve_method ppcg`, 100 steps | `WP7-E projected PCG, precond=localinv gpu build=device`; charges vs CPU ppcg 2.3e-11 e (first 10), group sums 4.4e-16 |
| Without groups bitwise | caffeine (24), polystyrene20 (325 atoms): `-sp` gradient dump, 20 MD steps | new build == unchanged build (md5; md.csv without wall-clock columns) |

## Not tested

- The ctests `gfnff_rocm_charge_groups`, `gfnff_gpu_charge_groups`, `gfnff_charge_groups`, `gfnff_cn_zero_guard` were not run
  through ctest: their scripts import numpy, which this machine lacks; the gate scripts above re-implement their checks.
- Systems above 108 atoms with groups, `"rest": "per_fragment"` files with many neutral solvent groups (nfrag >= 16, the
  ROCm many-fragment route and `eeq_rocm_cpu_fragment_threshold`), mixed-precision EEQ, `sytrf`, multi-GPU, gfx1150.
- Timing (no per-step timing was measured).

## Not implemented relative to CUDA

- Nothing found for this item. Unchanged on both backends: a charged multi-fragment system with groups under the default
  `frag_charge_model ensemble` still takes the ensemble CPU fallback.
- `k_subtract_qtmp` (unlaunched in both backends) has no guard flag, as on CUDA.
- The `-gpu_strict` message of the unchanged build named the ensemble reason for a guard-only fallback (caffeine, charge 0);
  not changed here, the path is gone on ROCm.

## Proposed lines for central files (coordinator)

- TODO.md: remove the item "ROCm: Ladungsgruppen und `cn_zero_guard` bleiben auf der CPU".
- AIChangelog.md: `- **GFN-FF charge groups and cn_zero_guard on ROCm without CPU fallback (Nov 24, 2026, 🤖 ⚙️ gfx1200)**: trait has_charge_groups_cn_guard true; GPU vs CPU gradient 6e-15 Eh/A, FD 5e-6 Eh/A, 200 MD steps charges <= 2e-13 e, ctest gfnff_rocm_charge_groups; docs/ROCM_PARITY_A1.md.`
- docs/GPU_PARITY_CUDA_ROCM.md, row "GFN-FF charge groups and `cn_zero_guard` on the device": ROCm column "yes (Nov 24, 2026, gfx1200, ROCM_PARITY_A1.md)".
- docs/GPU_PERIODIC_GAPS.md B7 and `scripts/gpu_parity_known.json`: `has_charge_groups_cn_guard` is no longer a difference.
- Root CLAUDE.md: nothing.

## Plugin ABI

No change needed: only a `static constexpr` trait of the ROCm backend (inside the plugin) and comments changed; no
`ComputationalMethod`, entry point or shared class layout changed (`gfnff.h` changed PARAM help strings only).
