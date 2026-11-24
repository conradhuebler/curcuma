# ROCm parity plan (Nov 24, 2026)

Goal: every GPU path that CUDA has also runs with `-gpu rocm`, or is refused with a clear message and listed in
[GPU_PARITY_CUDA_ROCM.md](GPU_PARITY_CUDA_ROCM.md). CUDA is the specification: the ROCm code mirrors the CUDA code, same
expressions, same defaults, same fallbacks (counted, `-gpu_strict true` stops). Source of the gap list: `scripts/gpu_parity_check.py`,
[GPU_PERIODIC_GAPS.md](GPU_PERIODIC_GAPS.md) part B, tables A and C of the parity document.

State at the start: GFN-FF XB gradient (B1) and CN range (B10) are ported (`feature/rocm-parity`, one RX 9060 XT gfx1200 on the
development machine; no second GPU, no CUDA device).

## Items

| ID | Item | CUDA reference | Wave |
|---|---|---|---|
| A1 | GFN-FF charge groups and `cn_zero_guard` on the device (trait `has_charge_groups_cn_guard` true; TODO.md "ROCm: Ladungsgruppen") | commit `9e8a609c`, `k_coulomb_postprocess`, `prepareEEQParametersForGPU` | 1 |
| A2 | GFN-FF periodic foundation (GPU_PERIODIC_GAPS B2): minimum image in all device kernels, cell upload, `template<bool PBC>` local terms, device D4/repulsion pair builders, device HB detection with a cell, D4 tail, trait `has_periodic_local_terms` true | commit `8ee4f638` (S1), `gfnffUploadCell`, `micPlace`, `k_hb_detect_pbc`, `k_d4_tail_*` | 1 |
| A3 | Periodic native GFN1/GFN2 on ROCm (B4): `setPeriodicDeviceScf` for the HIP backend, device stages (CN, D3/D4 with cell, ATM), virial, barostat | `xtb_gpu_method.cpp`, `k_cn_pbc`, `k_d3_*_pbc`, `k_d4_*_pbc`, `k_atm_pbc_nb`, [WP-PERIODIC-XTB.md](WP-PERIODIC-XTB.md) | 1 |
| A4 | Spin polarisation in the GFN1/GFN2 SCF on ROCm (`supportsSpinPolarization`, today a hard error) | host-driven device loop, [SPIN_POLARISATION.md](SPIN_POLARISATION.md) section GPU | 1 |
| A5 | Remaining GFN1/GFN2 ROCm gaps: `-scf_precond frag` on the ROCm resident loop (today ignored with a warning), sparse integrals / `-large_system_mode` on the device (`k_build_fock_iso`, `k_build_fock_sp`, `k_density_sp_sym`) | [SCF_PRECOND.md](SCF_PRECOND.md), [SQM_LARGE_SYSTEMS.md](SQM_LARGE_SYSTEMS.md), `xtb_gpu_context.cu` | 1 |
| A6 | Multi-GPU on ROCm: `has_multi_gpu_split` (GFN-FF mixed-precision split EEQ), multi-GPU eigensolver feasibility (cuSOLVERMg has no rocSOLVER counterpart: decide port or documented refusal). Only one GPU is available: no multi-device run is possible, state it | [MULTI_GPU.md](MULTI_GPU.md), `eeq_solver_gpu.cu` | 1 |
| B1 | Periodic GFN-FF hybrid on ROCm: device Ewald (S2/S3: `GfnffEwaldGPU` counterpart, pcg operator), host-charge hybrid, `hasPeriodicEnergy()` | `cuda/ewald_gpu.cu`, ctests `gfnff_gpu_pbc_hybrid`, `gfnff_gpu_pbc_ewald` | 2 (after A2) |
| B2 | Periodic GFN-FF pressure on ROCm: device virial (A1-GPU, B3), HB/XB crossing kernels, V-LRC, `updateCell` for the barostat (NPT with `-gpu rocm`), L1 fast list rebuild | `*_vir` kernels, `k_hb_crossing`, `k_xb_crossing`, `cuda/gfnff_virial_slots.h`, ctest `gfnff_gpu_pbc_virial` | 2 (after A2) |
| C | Parity sweep: re-run `scripts/gpu_parity_check.py`, every remaining CUDA-only kernel or trait either ported or recorded in `scripts/gpu_parity_known.json` with a reason; `ctest -L structinv`, full `ctest -j1`, determinism (md5 of two runs), plugin ABI version, update the parity document, GPU_PERIODIC_GAPS, TODO, AIChangelog | whole tree | 3 |

Order: wave 1 in parallel (A2 first in priority, B1 and B2 depend on it), wave 2 after A2 is merged into `feature/rocm-parity`, wave 3 last.

## Rules for every agent

- Work in the own worktree. Build in `<worktree>/release`: `cmake -DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON -DUSE_AVX2=ON -DCMAKE_BUILD_TYPE=Release ..`,
  `make -j3` (10 cores shared by all agents); only the targets needed (`curcuma`, the ROCm plugin).
- One GPU is shared: every run that touches the GPU goes through `flock /tmp/claude-1000/rocm-gpu.lock <command>`; timing and hash
  comparisons only inside such a call.
- Redirect build and run output into files, read with `tail -n`/`grep -c`; never put logs into the context.
- Change only what the item needs. Do not edit `TODO.md`, `AIChangelog.md`, `CLAUDE.md`, `docs/GPU_PARITY_CUDA_ROCM.md`, `docs/GPU_PERIODIC_GAPS.md`,
  `scripts/gpu_parity_known.json` or `CURCUMA_GPU_PLUGIN_ABI_VERSION` (central files, merged by the coordinator); put the proposed lines into the status file.
- No `git commit`, no `git add -A`, no push, no branch switch. Leave the changes uncommitted in the worktree.
- Gates (from GPU_PERIODIC_GAPS "Gates that apply to every item"): finite-difference check of the assembled surface, GPU against CPU
  (energy 1e-8 Eh, gradient 1e-8 Eh/A, virial 1e-7 Eh unless the existing test sets another limit), second run md5-identical, without a
  cell/new switch bitwise equal to the unchanged build (build the base once from the worktree HEAD for the comparison), new device
  path falls back to the old one with a counted fallback, unsupported combinations refuse with a message.
- New code is marked `Claude Generated`, with literature or CUDA reference in a comment. Labels: only 🤖 and ⚙️, never ✅.
- Status file `docs/ROCM_PARITY_<ID>.md` (dated, English): what was ported (file, function), what was measured (structure, n, numbers), what was
  **not** tested, what is **not** implemented relative to CUDA, proposed lines for TODO/AIChangelog/parity document, plugin ABI change needed (yes/no, why).
- Final report to the coordinator: at most 25 lines, numbers not dumps.
