# GPU gaps for periodic simulations with pressure (CUDA here, ROCm elsewhere)

Claude Generated (Nov 22, 2026). 🤖 AI-generated list; every entry names its source: [read] = read in the code by the author of
this list or by a pre-merge audit (`08_premerge` audits A and B of Nov 22, 2026), [reported] = stated in an agent report and not
re-checked, [open] = not examined. Nothing here is human-tested.

## Goal and order of work

**Goal (operator decision, Nov 22, 2026):** GFN-FF and native GFN1/GFN2 run complete periodic simulations **with pressure**
(NPT: virial, barostat) on every supported backend. Order: first close every CUDA gap on the CUDA machine (RTX 5080 development
machine), then the ROCm gaps on a machine with a ROCm GPU. A gap is closed only when its gate below passes.

## Where we stand (Nov 22, 2026)

| Method | CPU | CUDA | ROCm | Vulkan |
|---|---|---|---|---|
| GFN-FF periodic energy, gradient | complete (MIC, Ewald, EEQ) | hybrid S1/S2: local terms on the device, EEQ and Ewald on the host [read] | refused with a clear error, counted CPU fallback [read] | not built |
| GFN-FF virial, pressure, barostat | complete (N1, N2, V-LRC) | complete for the hybrid S2 (A1 closed Nov 23, 2026: device virial + host Ewald, `-compute_pressure`, `-barostat`) [measured, WP-PERIODIC-NONBONDED A1-GPU] | missing | not built |
| GFN1/GFN2 periodic energy, gradient | complete (Gamma point, images, Ewald) | complete, device-resident SCF (closed shell, neutral) [reported] | refused, counted fallback [read] | not built |
| GFN1/GFN2 virial, pressure, barostat | complete (P4) | complete [reported: parity 1e-9 Eh, NPT smoke test] | missing | not built |

Since A1 (Nov 23, 2026) GFN-FF NPT runs on CUDA through the S2 hybrid; what is left for GFN-FF on CUDA is speed (A2: host EEQ
and Ewald) and the approximate pressure terms (A3). GFN2 is complete on CUDA for closed-shell neutral cells.

## Part A: CUDA gaps (close on the CUDA machine)

Priority P0 = needed for the goal, P1 = correctness or honesty, P2 = performance, P3 = later.

| # | P | Gap | Source | Gate |
|---|---|---|---|---|
| A1 | P0 | **closed Nov 23, 2026** (`feature/gfnff-gpu-virial`): device virial of every local term and both chain rules, host Ewald virial, cell update followed by the device; GPU vs CPU <= 5e-11 Eh per component (432/864/1728 H2O, mixture, 5 small cells), FD of the hybrid energy under isotropic and axial cell scaling <= 5.6e-9 Eh, second run bitwise, no-cell CPU and CUDA 38/38 bitwise to the base, CLI NPT 432 H2O 20 steps CPU = GPU in every printed digit; MD step 1.7x faster (432-1728 H2O) | [measured] WP-PERIODIC-NONBONDED A1-GPU, ctest `gfnff_gpu_pbc_virial` | met |
| A2 | P0/P2 | GFN-FF stages S3/S4: device Ewald and device EEQ with a cell (today the host solves EEQ and Ewald: 107 of 123 ms per step at 5748 atoms) | [reported] WP-PERIODIC-NONBONDED S2 | charges and energy vs host <= 1e-9; speed table 3456 H2O |
| A3 | P0 | Pressure terms of GFN-FF that are still approximate or open: HB/XB and repulsion pairs crossing the cutoff (open in `TODO.md`), HB search radii not converged in the pressure (x0.8: +87 bar, x1.2: -29 bar) | [read] `TODO.md`, Versuchsplan P1 | pressure vs a larger-cutoff reference; document the residual |
| A4 | P1 | Silent host fallback of the periodic xTB device stages: not counted, visible only at verbosity 3, `-gpu_strict` does not trigger (e.g. `xtb_periodic.cpp:163-174`) | [read] audit B | `reportGpuFallback` for every stage; test with fault injection |
| A5 | P1 | `-large_system_mode` with a cell returns a vacuum energy (`native_xtb_method.cpp:184-223` returns before the cell is set); docs contradict each other (`WP-PERIODIC-XTB.md` lines 45, 167, 1376) | [read] | hard error with a clear message; test |
| A6 | P1 | Untested combinations on the device: open shell and charged cells (allowed on the host-driven loop), spin polarisation together with the periodic device stages (merged Nov 22, 2026; excluded by `!m_pbc && !spin_pol` and a hard error with a cell, never run together), multi-GPU with a cell (never ran; FP32 refinement only for one device) | [read] audit B, merge | one test per combination or an explicit refusal |
| A7 | P1 | Fallback paths of the new device stages were never triggered (failed device potential, failed resident step, failed device matrices) | [reported] | fault-injection test |
| A8 | P1 | Plugin ABI 25: every older CUDA plugin archive is incompatible (also the ROCm one); release notes must say so | [read] `gpu_plugin_abi.h` | note in `docs/RELEASE_BINARIES.md` |
| A9 | P1 | Tests: `pbc_xtb_gpu_parity` has no `RESOURCE_LOCK "gpu"` and reports SKIP as passed without a GPU (`test_cases/CMakeLists.txt` ~807); no ctest for the ROCm/Vulkan refusal with a cell (runs without a GPU) | [read] audit B | lock added; SKIP shown as skipped |
| A10 | P1 | `-xtb.gpu_eigen_strategy auto` was measured on the RTX 5080 only; datacenter GPUs (H200, A100) and older cards (A4500, GTX 1660) use rules derived from device attributes, unmeasured; cuBLAS FP64 emulation guarded by `CUBLAS_VER_MAJOR >= 13` only | [reported] | run on one datacenter GPU and one older card, record in the table of `WP-PERIODIC-XTB.md` |
| A11 | P1 | GFN-FF diagnostic switches (`hb_scale`, `*_scale`) fall back to the CPU under `-gpu`; the fallback was never run (`TODO.md:70`) | [read] | one test |
| A12 | P2 | GFN2 host remainder: ATM (~186 ms at 432 H2O), FP64 acceptance solve, first SCF iteration, small O(N) items; GFN1 S/H0 integrals (208 ms per call) and H0/Pulay kernel (202 ms) | [reported] | timing table |
| A13 | P2 | Periodic D3 pair term scales badly on the CPU (276 -> 179 ms for 1 -> 8 threads); Ewald point-charge gradient and KO gradient stay serial (~100 ms at 432 H2O) | [reported] | timing table |
| A14 | P3 | Device-resident MD (integrator, thermostat, barostat on the device) | [read] vault note "GPU-residente MD" | after S3/S4 |
| A15 | P3 | Halogen bond with a cell in xTB (GFN1): Br and I have d shells, which a cell refuses; no XB parity cell | [reported] | needs periodic d-shell support or a documented limit |

## Part B: ROCm gaps (close on a machine with a ROCm GPU)

ROCm was only ever run on a Radeon 890M (gfx1150). Nothing below can be checked on the CUDA machine.

| # | Gap | Source |
|---|---|---|
| B1 | GFN-FF XB gradient fix (`fix/gfnff-xb-gradient`): CPU and CUDA fixed, ROCm NOT (`gfnff_rocm.hip` ~2498-2505 and ~2630-2637 unchanged) | [read] audit A |
| B2 | Periodic GFN-FF local terms (S1/S2 equivalent: minimum image in all device kernels, `cellFromMolecule`, FMA-safe `template<bool PBC>`); today ROCm refuses a cell | [read] |
| B3 | GFN-FF device virial (A1 equivalent) so that NPT runs with `-gpu rocm` | derived from A1 |
| B4 | Periodic xTB on ROCm (`xtb_hip_method.cpp` never calls `setPeriodicDeviceScf`): device stages, virial, barostat; the refusal path has no test | [read] audit B |
| B5 | Plugin ABI 25: the ROCm plugin must be rebuilt and loaded once; the new core symbol `EEQSolver::hostSytrsUpper` is resolved at load time | [read] merge |
| B6 | The HB-CN-after-refresh change sits in the shared template (`gfnff_gpu_method_impl.h`) and changes GPU MD trajectories at the refresh steps; compiled for ROCm, never run | [read] audit A |
| B7 | Parity script: two differences without a reason (`has_charge_groups_cn_guard` from `9e8a609c`, `has_periodic_local_terms` from `8ee4f638`); add them to `docs/GPU_PARITY_CUDA_ROCM.md` and `scripts/gpu_parity_known.json` with the reason or close them | [read] `scripts/gpu_parity_check.py` |
| B8 | Bitwise reproducibility of the periodic paths on ROCm (`-gfnff.rocm_fixed_point`, `-xtb.rocm_fixed_point`) | [open] |
| B9 | Vulkan: code and `USE_VULKAN` exist but nothing is built or shipped; either drop it or give it the same refusal tests | [read] profile |

## Gates that apply to every item

1. FD gate of the assembled surface (energy under displacement, and virial under cell scaling at fixed fractional coordinates).
2. Parity against the CPU path on the cells of `pbc_xtb_gpu_parity` and the GFN-FF periodic test set (energy 1e-8 Eh, gradient 1e-8 Eh/A, virial 1e-7 Eh).
3. A second run md5-identical (deterministic reductions, no atomic adds in fixed-order sums; `docs/GPU_REPRODUCIBILITY.md`).
4. Without a cell bitwise equal to the base build (build the base with the same compiler and flags; gcc and CUDA minor versions change the last bits).
5. A new device path falls back to the previous path on a failure, and the fallback is counted (`-gpu_strict` stops at the first).
6. Full `ctest -j1` in the CPU and the CUDA build; unsupported combinations refuse with a clear message.

## Related documents

[WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md) (GFN-FF periodic, S0-S2), [WP-PERIODIC-XTB.md](WP-PERIODIC-XTB.md) (xTB periodic, CUDA stages), [NPT.md](NPT.md), [GPU_PARITY_CUDA_ROCM.md](GPU_PARITY_CUDA_ROCM.md), [GPU_REPRODUCIBILITY.md](GPU_REPRODUCIBILITY.md).
