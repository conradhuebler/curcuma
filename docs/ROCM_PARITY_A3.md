# ROCm parity A3: periodic native GFN1/GFN2 on ROCm (Nov 24, 2026)

Item A3 of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md), gap B4 of [GPU_PERIODIC_GAPS.md](GPU_PERIODIC_GAPS.md).
🤖 AI-generated, ⚙️ machine-tested on one AMD RX 9060 XT (gfx1200, ROCm 7.2) only; human production testing pending.

## State before

- `EnergyCalculator::setMolecule` swapped `-gpu rocm` + cell to the CPU method (counted fallback "native xTB: periodic
  cell set, rocm path has no lattice images"); `XTB::Calculation` refused a cell with a backend that did not set
  `setPeriodicDeviceScf` (hard error). `xtb_hip_method.cpp` never set it.
- The refusal had a test: `pbc_xtb_cell_refusals` checked `-gpu rocm|vulkan` with a cell for "error or counted CPU
  fallback with the CPU energy" (it now accepts the ROCm device energy within 1e-8 Eh, see below).

## What was ported (file, function)

- `src/core/energycalculator.cpp` `setMolecule`: the P5 swap no longer applies to `-gpu rocm` (Vulkan keeps it).
- `xtb_hip_method.cpp`: `XtbHipComputationalMethod` sets `XTB::setPeriodicDeviceScf(true)`; `HipScfBackend`
  ("BACKEND-SPECIFIC 12") forwards the CUDA periodic overrides verbatim: `cnPeriodic`, `d4CovalentCnPeriodic`,
  `d4CnFoldPeriodic`, `cnChainGradientPeriodic`, `repulsionPeriodic`, `d3DispersionPeriodic`, `d4PairsPeriodic`,
  `dispersionATMPeriodic`, `beginDispersionImages`/`dispersionDedqPeriodic`, `ewaldGammaPeriodic`,
  `multipoleMatricesPeriodic`, `multipoleGradientPeriodic`, `ewaldGradientPeriodic`, `koGradientPeriodic`,
  `h0PulayGradientPeriodic`, `energyWeightedDensity`; new `beginMultipole` (upload of the host-built GFN2 multipole
  integrals; without it GFN2 + cell ran the CPU SCF, because the shared template returned false).
- `rocm/xtb_hip_periodic_kernels.hiph` (new): the CUDA kernels of these stages, extracted from `cuda/xtb_gpu_context.cu`
  with mechanical changes only (`gfix_add` -> `pfix_add` with the same 2^52 fixed point and chunking, wave32
  `__shfl_down`/`__ballot` instead of the `_sync` forms, `c_pauling`/`c_atomic_rad` constant tables -> kernel
  arguments `dpauling`/`darad`), namespace `pbc_dev`; `PbcBuf`/`PbcState` buffers (member `Impl::pbc`).
- `rocm/xtb_hip_periodic_host.hiph` (new): the CUDA host entry points with cuda -> hip, cuBLAS -> rocBLAS
  (`rocblas_dsyrk`, `rocblas_dgemm`), the beginDispersion data mapped to the HIP `Impl` members; plus
  `uploadMultipoleInts` and `residentEnergyWeightedDensity`.
- `xtb_native.h`: new virtual `GpuScfBackend::periodicStageOnDevice(token)` (default true). HIP returns false for
  `ints` and `devpot`; `xtb_native.cpp` (ints, devpot) and `xtb_periodic.cpp` (h0grad) then take the host path without
  a fallback count, as every host-only part of the periodic path. A device call that fails is still counted.
- Tests: `test_pbc_xtb_gpu_parity` and `test_pbc_xtb_gpu_fallback` take the backend as argument, `test_pbc_xtb_combinations`
  gained the part `rocm`; ctests `pbc_xtb_gpu_parity_rocm`, `pbc_xtb_gpu_fallback_rocm`, `pbc_xtb_gpu_combinations_rocm`
  (only with `USE_ROCM`). Docs: WP-PERIODIC-XTB status lines, `qm_methods/CLAUDE.md`.

## What was measured

All runs `release/` of this worktree, `scf_threshold` 1e-9 unless noted, GPU runs under the shared lock.

- `pbc_xtb_gpu_parity_rocm` (10 cases: water dimer 8 and 12 A, water trimer 6.5 A, 216 H2O 18.6 A with 4 threads,
  caffeine 14 A; gfn1 and gfn2): ROCm vs CPU worst |dE| 4.3e-12 Eh, max |dG| 2.1e-10 Eh/A, max |dXi| 4.6e-10 Eh
  (limits 1e-8 / 1e-8 / 1e-7); no stage fallback. With the final code the test ran 4 times: passed once; the three
  failures were only the run-to-run check (216 H2O box in all three, gfn1 and/or gfn2, dE 5e-13 to 5e-12 Eh between two
  calculators in one process; caffeine gfn2 once, 7e-15 Eh), see "Open".
- `pbc_xtb_gpu_fallback_rocm`: 24 stage tokens (gfn2 14, gfn1 10) injected: each counted, energy/gradient/virial equal to
  the host variant (worst 5e-16 Eh/A); `-gpu_strict` exits 3; a normal run counts 0 stage fallbacks. Passed.
- `pbc_xtb_gpu_combinations_rocm`: open shell (OH, spin 1, 8 A) and charged cell (hydronium_water2, +1, 10 A), gfn1/gfn2
  on the host-driven loop: worst |dE| 1.6e-14 Eh, |dG| 1.2e-10 Eh/A, |dXi| 1.2e-10 Eh, second run bitwise, 0 fallbacks.
- FD gate of the assembled surface on ROCm (scratch program, central differences of `-gpu rocm` energies, h 1e-4 A;
  strain h 1e-5, Xi_aa = -dE/d eps_aa and the isotropic strain): water dimer 8 A gfn1 / gfn2 8.1e-9 / 8.9e-9 Eh/A,
  virial 5.8e-10 / 9.1e-10 Eh; caffeine 14 A gfn1 / gfn2 2.9e-8 / 2.9e-8 Eh/A, virial 8.6e-9 / 9.2e-9 Eh. The same
  program on the CPU gives the identical caffeine numbers (2.90e-8 / 9.23e-9), so that part is FD truncation.
- Run-to-run in separate processes: 216 H2O gfn1 (CPU 4 threads, ROCm 4 and 1 threads) 3/3 md5-identical each;
  caffeine 14 A gfn2 8/8 identical with and without `-scf_mixed_precision`; in one process without a preceding CPU
  calculation 3/3 to 6/6 identical (216 H2O, caffeine).
- Without a cell, CPU: 12/12 md5-identical to the base build (water dimer, caffeine, acetic acid dimer x gfn1/gfn2 x
  1/4 threads). ROCm: 7/12 identical; the base build itself is not run-to-run reproducible there (caffeine gfn1: 3
  distinct hashes in 3 base runs), and with FP64 only (`-scf_mixed_precision false`, caffeine gfn1, 8 runs each) the
  most frequent hash is the same in both builds (base 7/8, new 6/8; the outliers differ per run). The vacuum ROCm path
  is changed only by the new `beginMultipole`, which is reached only when the device integral build fails.
- NPT smoke: 216 H2O gfn1, `-gpu rocm -barostat berendsen -compute_pressure true`, 20 steps of 0.5 fs, seeded: runs,
  no fallback, cell 18.6258 -> 18.6226 A (n = 1, no physics claim).

## Not tested

- Any GPU other than one gfx1200; wave64 hardware (refused at configure time); multi-GPU with a cell.
- Larger cells than 216 H2O, other elements than H/C/N/O, MD energy conservation with a cell, NPT beyond the 20-step
  smoke, `-opt` with a cell on ROCm, thread counts other than 1 and 4.

## Not implemented relative to CUDA

- G1b device periodic integrals (`beginComputedPeriodic`, `k_overlap_h0_pbc`, `k_multipole_ints_pbc`): the host builds
  S/H0/L and the multipole integrals and uploads them (`ints` host by design, not counted).
- R2-R5 resident periodic potential and loop (`beginPotentialPeriodic`, device copies of the periodic gamma and multipole
  matrices, `periodicGammaResident`, `downloadConvergedPotentials`): ROCm keeps the host-driven device loop for every
  periodic SCF (`devpot` host by design; `resloop`/`potdl`/`mpdl` unreachable).
- A1r FP64 refinement of the FP32 solve, G1c host-loop pseudo-diagonalisation, E1 density complement: the HIP backend
  keeps the default no-ops of these switches (plain ROCm solves).
- Superseded CUDA comparison kernels (`triview`, `mppair`, `dedq1`) were ported along (same INACTIVE status).

## Open

- In-process nondeterminism of the ROCm SCF: after a CPU calculation in the same process, repeated `-gpu rocm`
  calculations of the 216 H2O box differ at 1e-12 Eh (3 distinct of 3; also with every periodic stage on the host
  (`CURCUMA_PBC_GPU_HOST=all`), with FP64 only, and without a cell). Not caused by the periodic stages; cause not found.
  It makes `pbc_xtb_gpu_parity_rocm` fail intermittently on its run-to-run check.

## Proposed lines (central files, for the coordinator)

- TODO.md: "ROCm: SCF nicht reproduzierbar im selben Prozess nach einer CPU-Rechnung (216 H2O, auch ohne Zelle, auch
  FP64), Ursache offen; `pbc_xtb_gpu_parity_rocm` scheitert daran zeitweise (docs/ROCM_PARITY_A3.md)." and "ROCm periodisch:
  G1b (Geraete-Integrale) und R2-R5 (residente periodische Schleife) fehlen (docs/ROCM_PARITY_A3.md)."
- AIChangelog.md: "2026-11-24 ROCm parity A3: periodic native gfn1/gfn2 with `-gpu rocm` (host-driven device loop, CUDA
  stages G2-G5/H1-H5/B1/C1/D1 ported; vs CPU <= 4.3e-12 Eh, 2.1e-10 Eh/A, 4.6e-10 Eh, 10 cells) -> docs/ROCM_PARITY_A3.md"
- GPU_PERIODIC_GAPS.md B4: "partly closed Nov 24, 2026: device stages, virial, barostat (smoke); open: G1b, R2-R5,
  in-process reproducibility (docs/ROCM_PARITY_A3.md)".
- GPU_PARITY_CUDA_ROCM.md / `gpu_parity_known.json`: CUDA-only kernels `k_overlap_h0_pbc`, `k_multipole_ints_pbc` (G1b not
  ported on ROCm); seam `periodicStageOnDevice` new (HIP false for ints, devpot).
- Root CLAUDE.md "GPU" bullet: "periodic gfn1/gfn2 on ROCm (host-driven device loop, no G1b/R2-R5): docs/ROCM_PARITY_A3.md".

## Plugin ABI

Change needed: yes. `GpuScfBackend` (xtb_native.h, shared by core and plugins) gained the virtual
`periodicStageOnDevice`, so a plugin built before this change has a different vtable; raise
`CURCUMA_GPU_PLUGIN_ABI_VERSION` once (not edited here, central file). The CUDA plugin needs no code change.
