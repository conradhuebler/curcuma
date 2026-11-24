# ROCm parity B1: periodic GFN-FF hybrid on ROCm (Nov 24, 2026)

Claude Generated. 🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT (gfx1200); human production testing pending.
Item B1 of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md). CUDA reference: S2 (host-charge hybrid) and S3
(device Ewald, `cuda/ewald_gpu.cu`) of [WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md). Base: `feature/rocm-parity`
`d4dbdfc5` plus the coordinator patch `base_a2.patch` (build fix, B1/B10 XB and CN, plan, A2). Nothing committed.

## What was ported

S2, host-charge hybrid (`-gpu rocm -gfnff.gpu_periodic true`, `coulomb_ewald` true = default):

- `rocm/gfnff_rocm.hip`, `FFWorkspaceHip::launchChargeDependentAndFinish()`: the two pair Coulomb launches
  (`launchCoulombTilesHip`, `k_coulomb_gather`) run only with `m_coulomb_pair_on_device` (CUDA S2 text); the self
  term and TERM 1b (`k_coulomb_postprocess`) stay. Before B1 the flag had no reader in the HIP unit.
- `qm_methods/gfnff_hip_method.h`: `has_periodic_hybrid = true`. The hybrid paths of the shared wrapper
  (`gfnff_gpu_method_impl.h`: host Ewald set-up check, `setCoulombPairOnDevice(false)` in `setMolecule()` and
  `initGPUWorkspace()`, charges through `eeqCpuFragmentRoute()` -> `calculateFinalChargesPeriodic`,
  `FFWorkspace::evaluateEwaldCoulomb()` added to energy and gradient, `hasPeriodicEnergy()` true) are the CUDA ones,
  unchanged.
- New trait `has_periodic_virial` (CUDA true, ROCm false), split from `has_periodic_hybrid`: it gates
  `setWantVirial`/`setHBCrossing` of the device workspace and `assembleHybridVirial()` (B2 territory, not defined in the
  HIP unit). On ROCm the wrapper calls `GFNFF::setBackendVirialUnavailable()` (new inline setter in `gfnff.h`, stored
  in `m_parameters["_backend_virial_unavailable"]`, no new member), so `GFNFF::virialUnavailableReason()` returns the
  reason and `setWantVirial(true)` (`-compute_pressure`, barostat) refuses at set-up; a second guard in
  `calculateEnergy()` prints an error once if a virial is requested anyway.
- Messages: wrapper refusal texts and `MethodFactory::periodicCellNeedsCpu()` say "CUDA and ROCm".

S3, device Ewald:

- `rocm/ewald_hip.h`, `rocm/ewald_hip.hip` (own hipcc translation unit, CMake target `curcuma_hip_ewald`, linked into
  `libcurcuma_rocm.so`): `GfnffEwaldHip` implements `curcuma::ewald::EwaldDeviceOps`, generated statement for
  statement from `cuda/ewald_gpu.cu` (`k_ew_wrapf`, `k_ew_count`, `k_ew_fill`, `k_ew_kernel`, `k_ew_phase`, `k_ew_sf`,
  `k_ew_kpass<GRAD>`, `k_ew_real_pot`, `k_ew_real_eg<GRAD, VIR>`, `k_ew_sum_items<NV>`, `k_ew_kitems`; pair list, kernel,
  phase table, pcg operator, energy/gradient/virial at fixed charges). Differences: 64-bit shuffle masks with width 32,
  the ballot of `k_ew_fill` shifted to the 32-lane group (wave64-safe), `hipcub::DeviceScan`, `hipHostMalloc`.
  Fixed-order sums, no atomics. The pair virial (`k_ew_real_eg<true,true>`, `k_ew_kitems`) is ported with it (same
  kernels as CUDA) and tested on the object (below), not through the hybrid.
- `GFNFFRocmBackend::createEwaldDevice()` (`gfnff_hip_method.cpp`) returns it (nullptr without a HIP stream: counted
  fallback `gfnff_periodic_ewald_device`, host Ewald). `-gfnff.gpu_ewald_device false` = S2 in the same binary; fault
  injection `CURCUMA_GFNFF_EWALD_DEVICE_FAIL=prepare|potentials|evaluate`, profile `CURCUMA_GFNFF_EWALD_PROFILE=1`.

Tests (same sources as CUDA, ROCm variants by a define, 77 without a HIP device):
`test_gfnff_rocm_pbc_hybrid` (`GFNFF_PBC_HYBRID_ROCM`; env `GPU_HYBRID_EWALD_DEVICE=0` runs it on the S2 host Ewald),
`test_gfnff_rocm_pbc_ewald` (`GFNFF_PBC_EWALD_ROCM`; checks (a)-(d) without the hybrid virial, plus (e) the device
Ewald object against the CPU `EwaldSum` range functions incl. the pair virial; (e) is compiled for ROCm only),
`test_gfnff_rocm_pbc_refuse` (`GFNFF_PBC_REFUSE_ROCM`). ctests `gfnff_rocm_pbc_hybrid`, `gfnff_rocm_pbc_ewald`,
`gfnff_rocm_pbc_refuse` (labels `gfnff_gpu;gpu;pbc;rocm`, lock `gpu`).

## Measurements

One RX 9060 XT (gfx1200), ROCm 7.2.4, build `-DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON -DUSE_AVX2=ON`, every GPU run under
`flock /tmp/claude-1000/rocm-gpu.lock`. Base = this worktree with the base patch, built before the B1 edits (binary,
plugin and `gfnff_rocm.hip.o` saved). Thresholds are the ones of the CUDA tests.

| gate | result |
|---|---|
| S2: `gfnff_rocm_pbc_hybrid` with `GPU_HYBRID_EWALD_DEVICE=0` (host Ewald): 50 H2O L 44 (dense), H3O+(H2O)2 + 18 H2O Q = +1, 970 H2O L 58 (pcg, 4 threads), input and moved geometry, default tol and `eeq_ewald_tol` 1e-14 | all checks passed; GPU vs CPU charges <= 1.0e-15 e, energy <= 3.2e-12 Eh (limit 1e-9; tight 1e-11), gradient <= 4.0e-14 Eh/A (limit 1e-12) |
| S2: FD gate (same test), 6 atoms nearest a face, h 2e-4 / 1e-4 A | 20 H2O 8.3e-8 -> 2.1e-8, Richardson 6.3e-10 Eh/A; Q = +1 7.9e-8 -> 2.0e-8, Richardson 8.2e-10 (limit 1e-8) |
| S2: invariances (translation, torn, initialised torn, 2x1x1 supercell) and two instances bit for bit (1 and 4 threads) | passed (charges <= 1.7e-15 e, energy <= 4.4e-15 Eh, gradient <= 3.1e-15 Eh/A in the printed torn lines) |
| S3: `gfnff_rocm_pbc_hybrid` default (device Ewald) | all checks passed; charges <= 4.3e-15 e, energy <= 3.1e-12 Eh, gradient <= 4.0e-14 Eh/A; FD Richardson 6.3e-10 / 8.0e-10 Eh/A |
| S3: `gfnff_rocm_pbc_ewald` (a): 20 H2O (FD settings), 60 H2O, Q = +1, 970 H2O (pcg), Na2Cl2 + 380 H2O L 48 (pcg), input and moved geometry | 137 checks pass, 0 failures; device vs CPU charges <= 3.8e-15 e, energy <= 3.0e-12 Eh, gradient <= 7.4e-14 Eh/A; device vs host-Ewald GPU energy <= 6.0e-13 Eh, gradient <= 1.8e-15 Eh/A; device used (7 operator calls per solve at 2910 and 1144 atoms), no fallback counted |
| S3 (b) FD gate of the device-Ewald hybrid, dense and forced pcg (558 device operator calls) | gradient <= 9.6e-10 Eh/A (limit 1e-8); virial FD not run (no device virial, B2) |
| S3 (c) second instance, 970 H2O, two steps | energy, gradient, charges bit for bit (virial not compared: B2) |
| S3 (d) fault injection prepare / potentials / evaluate | energy = host-Ewald hybrid within 0 / 2.8e-13 / 5.7e-14 Eh, fallback counted 2 / 1 / 1 |
| S3 (e) `GfnffEwaldHip` vs CPU `EwaldSum` (stored-kernel cache, range functions), random neutral charges; 20 H2O, Q = +1 box, 970 H2O | real 5.6e-14 Eh, reciprocal 1.8e-14 Eh, gradient 9.6e-16 Eh/Bohr, potentials 9.9e-15 Eh/e, pair virial 4.0e-14 Eh (limits 1e-10 / 1e-9); repeated evaluate bit for bit |
| `gfnff_rocm_pbc_refuse` | all checks passed ((d) hybrid accepted, `hasPeriodicEnergy()` true, host Ewald part non-zero) |
| CLI `-sp -dump_gradient`, 4 threads: caffeine in a 25 A cube; 1728 H2O (2x2x2 replica of `water216-box.sc-s1`, L 37.25 A, 5184 atoms) | device Ewald vs CPU: caffeine dE 2.8e-14 Eh, max dG 1.6e-13 Eh/A; 1728 H2O dE 1.2e-11 Eh, dG 2.9e-14 Eh/A; host Ewald vs CPU 8.6e-12 Eh / 2.9e-14 Eh/A; two device runs md5-identical; `-gpu_strict true`, no fallback; `-gpu rocm` without `gpu_periodic`: counted CPU fallback, md5 = CPU |
| CLI MD 1728 H2O, 100 steps (dt 0.5 fs, seed 42, CSVR, 4 threads): device Ewald twice, host Ewald, CPU | status columns (Epot, Ekin, Etot, T at 6 decimals, 22 rows) identical in all four; trajectory md5 identical in all four; final state md5 identical between the two device runs; no fallback |
| CLI MD with `-compute_pressure true` on the ROCm hybrid (1728 H2O) | refused at set-up: "the atomic virial is not available: the GPU periodic hybrid on ROCm has no device virial ...", "MD setup failed", exit 1 (an intermediate build without `setBackendVirialUnavailable` ended in `std::terminate` from `SimpleMD::readVirial`, exit 134); CPU control (caffeine cell) runs with the pressure columns |
| no cell, 6 structures (water, caffeine, acetic acid dimer, trisaccharide, macrocycle-bgal, peo201-chain; `-sp -dump_gradient -gpu rocm -gpu_strict true -threads 1`), base vs B1 | 6/6 md5-identical; second B1 run 6/6 identical; base twice 6/6 identical (run on the S3 build and again on the final build) |
| device code of `gfnff_rocm.hip.o`, base vs B1 (`.hip_fatbin` section extracted with objcopy) | byte-identical (1140480 bytes): the S2 edit is host code only |
| final build, `ctest -R "gfnff_rocm_pbc|cli_gfnff_gpu_0.*_rocm"` (refuse, terms, hybrid, ewald, three CLI) | 7/7 pass |
| `scripts/gpu_parity_check.py` | `has_periodic_hybrid` equal (true/true); one unexplained difference: the new `has_periodic_virial` (CUDA true, ROCm false); the `k_ew_*` kernels are no longer CUDA-only |

Timing (one MD run each, n = 1, machine shared with other agents, load average about 17): 100 MD steps of 1728 H2O took
9.6 s wall with the device Ewald, 21.3 s with the host Ewald, 26.7 s on the CPU (4 threads), start-up included. Not a
performance measurement; no profile was taken.

Not run: full `ctest -j1`, `ctest -L structinv` (only the targets the item needs were built), `-threads 1` vs 4 bitwise
comparison of the device-Ewald hybrid beyond the test's two-instance checks.

## Not tested

- CUDA build of the changed shared files (no CUDA toolkit here): the CUDA trait `has_periodic_virial` is true, so every
  CUDA `if constexpr` resolves as before; `test_gfnff_gpu_pbc_ewald.cpp`/`_hybrid.cpp`/`_refuse.cpp` changed only through
  aliases and `kHybridVirial` (true on CUDA: same checks), not compiled with USE_CUDA here.
- wave64 hardware (the HIP Ewald kernels are written for 32-lane groups on either wavefront size; only gfx1200, wave32,
  was run), gfx1150, ROCm `rocm_fixed_point false`, mixed precision, multi-stream mode with a cell above 500 atoms in
  the hybrid ctests (the CLI 5184-atom box ran it), charged boxes under the ensemble charge model (CPU fallback as on
  CUDA), react/`topology_mode auto` with a cell, NVE energy conservation of the ROCm hybrid, more than 5184 atoms,
  Windows. Human production testing pending.

## Not implemented relative to CUDA

- Device virial of the hybrid and the HB/XB crossing kernels (B2): with a cell on ROCm `-compute_pressure` and the
  barostat are refused (`has_periodic_virial` false).
- Multi-GPU with a cell: off on both backends.

## Interfaces for B2 (virial, pressure, NPT)

- Trait `GFNFFRocmBackend::has_periodic_virial` (`gfnff_hip_method.h`): set it to true when `FFWorkspaceHip::setWantVirial`,
  `downloadVirial` (declared in `cuda/ff_workspace_gpu.h`, slots `cuda/gfnff_virial_slots.h`, `GFNFF_VIR_BUFFER`,
  `GFNFF_VIR_TAIL_E8`) and the crossing passes behind `setHBCrossing` exist in `rocm/gfnff_rocm.hip`. The trait gates, in
  `gfnff_gpu_method_impl.h`: the `setWantVirial`/`setHBCrossing` block in `calculateEnergy()`, `assembleHybridVirial()`,
  the `setBackendVirialUnavailable()` call in `setMolecule()` and the once-per-instance error in `calculateEnergy()`
  (both disappear with the trait true).
- The Ewald part of the virial is ready: `GfnffEwaldHip::evaluate(q, terms, grad, &virial)` returns the real + reciprocal
  pair virial (test (e), 4.0e-14 Eh); `FFWorkspace::evaluateEwaldCoulomb` asks for it when the virial is wanted, and
  `assembleHybridVirial()` reads it from `ws->virialComponents()` (coulomb, coulomb_bg), as on CUDA.
- Tests to switch on for B2: `kHybridVirial` in `test_gfnff_gpu_pbc_ewald.cpp` (ROCm branch false today: model virial in
  (a), virial FD in (b), virial in (c)); a ROCm variant of `test_gfnff_gpu_pbc_virial.cpp` in the same pattern.
- Barostat: `followHostCellChanges()` (shared) already runs on ROCm; with a cell change the device Ewald re-keys on the
  cell (`prepare()` compares L), so no extra hook is needed for the Ewald object.

## Proposed central lines (coordinator)

- TODO.md: replace "ROCm periodic hybrid (B1) and virial (B2)" by "ROCm periodic virial and NPT (B2): docs/ROCM_PARITY_B1.md
  section 'Interfaces for B2'".
- AIChangelog.md: `2026-11-24 ROCm B1: periodic GFN-FF hybrid on ROCm (host periodic EEQ + Ewald, device pair Coulomb off) and
  device Ewald GfnffEwaldHip (HIP port of GfnffEwaldGPU); virial refused until B2 (docs/ROCM_PARITY_B1.md)`.
- GPU_PERIODIC_GAPS.md: B1/B3 hybrid and device Ewald on ROCm done Nov 24, 2026 (gfx1200); virial/NPT = B2.
- GPU_PARITY_CUDA_ROCM.md / gpu_parity_known.json: `has_periodic_hybrid` now equal; new trait `has_periodic_virial` differs
  (CUDA true, ROCm false; reason "B2 of ROCM_PARITY_PLAN_2026-11").
- ff_methods/CLAUDE.md: in "Periodic cell, CPU" replace "ROCm refuses" by "ROCm: hybrid and device Ewald since Nov 24, 2026 (B1),
  no virial"; in "Device Ewald of the CUDA hybrid" add "ROCm: `GfnffEwaldHip` (`rocm/ewald_hip.hip`)".
- `SimpleMD::readVirial` message names only CUDA for the GPU hybrid; with the set-up refusal it is not reached on ROCm.

## Plugin ABI

No change needed: no `ComputationalMethod` member, no plugin entry point, no shared class layout changed. `GFNFF` gained
an inline member function that writes a JSON key (`setBackendVirialUnavailable`), no data member; `GfnffEwaldHip` and
the wrapper template are plugin-private; the HIP workspace layout is unchanged (the S2 flag existed).
