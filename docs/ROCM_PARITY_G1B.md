# ROCm parity G1B: device periodic integrals (Nov 25, 2026)

🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT (gfx1200). Human production testing pending.
Worktree start `dd4812c4`; changes uncommitted.

## What was ported

CUDA reference: `cuda/xtb_gpu_context.cu` (`k_overlap_h0_pbc`, `k_multipole_ints_pbc`,
`XtbGpuContext::computeIntegralsPeriodic`) and `xtb_gpu_method.cpp` (`beginComputedPeriodic`), WP-PERIODIC-XTB G1b.

- `rocm/xtb_hip_periodic_kernels.hiph`: `pbc_dev::k_overlap_h0_pbc`, `pbc_dev::k_multipole_ints_pbc`. Bodies verbatim;
  the CUDA `__constant__` tables `c_pauling`/`c_atomic_rad` become kernel arguments (the HIP device tables `dpauling`/`darad`,
  as in the A3 port of `k_grad_h0_pulay_pbc`). New kernels only; no existing kernel changed. `PbcState` gained
  `dPbcIntTrans` and `ints_error`.
- `rocm/xtb_hip_periodic_host.hiph`: `XtbHipContext::computeIntegralsPeriodic` (+ `periodicIntegralError`). Writes the
  resident dense `dS0`/`dH0`/`dL` (and `dDpInt`/`dQpInt` for GFN2), uploads xyz and the host self-energies into
  `dxyz`/`dse`, L = chol(S) by `rocsolver_dpotrf` (lower, on a copy, info checked). Differences to CUDA: the resident
  buffers exist from `begBasis` (no allocation, no memory estimate here); a screened storage is turned dense
  (`spEnsureDense`, CUDA sets `I.sparse = false`); the Cholesky probe of CUDA (`probeCholeskyFactor`) is the HIP
  `factorCheck` at the first use of `dL` (`dL_ver = s_version`), as for the vacuum `computeIntegrals`; `cn_ok` stays false
  (the periodic CN is the host one).
- `xtb_hip_method.cpp`: `HipScfBackend::beginComputedPeriodic` override; `periodicStageOnDevice` now false only for
  `devpot`. Failure reason via `lastError()`; `xtb_native.cpp` then runs the host build and counts the fallback (token
  `ints`, `-gpu_strict` exits 3, injection `CURCUMA_PBC_GPU_FAIL=ints`, host switch `CURCUMA_PBC_GPU_HOST=ints`).
- `xtb_native.cpp` (shared host code): `m_pbc_mp_device_only` (R4) additionally requires
  `periodicStageOnDevice("devpot")`. Without it the first ROCm G1b run counted an `mpmat` fallback (R4 asked for
  device-only periodic multipole matrices, which ROCm does not keep). CUDA returns true there: no change on CUDA.
- Tests: `test_pbc_xtb_gpu_fallback` has `ints` in the ROCm token lists (gfn1 and gfn2); comments of it and of
  `test_pbc_xtb_gpu_parity` updated.

## What was measured

Build `release/` of this worktree (`-DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON -DUSE_AVX2=ON`); base = the same worktree at
`dd4812c4` built before any edit (binary + plugin copied to the scratchpad). Every GPU run under
`flock /tmp/claude-1000/rocm-gpu.lock`, one job at a time. Cells: water dimer (`water-dimer.roo2.8`) 8 A, caffeine 14 A,
216 H2O (`water216-box.sc-s1`) 18.6258 A (4 threads); `-scf_threshold 1e-9`, `-dump_gradient`.

| Gate | Result |
|---|---|
| ISA of `xtb_hip_context.hip.o` (scratch `a2/isa.py`), base vs new | 145/145 former functions identical, 2 new kernels |
| Device vs host S/H0 (`CURCUMA_PBC_INTS_CHECK`, water dimer 8 A gfn2) | max dS 1.7e-18, max dH0 1.3e-18 |
| ROCm vs CPU, CLI, 3 cells x gfn1/gfn2 | worst dE 3.4e-12 Eh (216 H2O gfn2), max dG 6.2e-11 Eh/A |
| `pbc_xtb_gpu_parity_rocm` (10 cases) | passed; worst dE 3.4e-12 Eh, dG 2.1e-10 Eh/A, dXi 3.3e-10 Eh; 10/10 second run bitwise |
| Second run, separate process (CLI, 6 cases) | 6/6 dump md5-identical |
| `CURCUMA_PBC_GPU_HOST=ints` (new) vs base, 6 cases | 6/6 dump md5-identical |
| No cell, new vs base (caffeine, water dimer x gfn1/gfn2), `-gpu rocm` and CPU | 4/4 and 4/4 md5-identical |
| New vs base with a cell (device vs host integrals) | dE <= 1.1e-12 Eh, dG <= 7.1e-15 Eh/A |
| FD of the assembled `-gpu rocm` surface (h 1e-4 A, central) | water dimer 8 A: all 18 components, gfn1 8.1e-9, gfn2 8.9e-9 Eh/A; caffeine 14 A: atoms 0, 5, 12 (9 components), 1.7e-8 Eh/A both methods; the same caffeine FD on the CPU gives 1.7e-8 (FD truncation) |
| `pbc_xtb_gpu_fallback_rocm` | passed; 26 tokens counted (24 + `ints` for gfn1 and gfn2), `ints` fallback = host variant to 0.0; `-gpu_strict` exit 3 |
| CLI `CURCUMA_PBC_GPU_FAIL=ints -gpu_strict true` | exit 3, message names stage `ints` |
| `pbc_xtb_gpu_combinations_rocm` | passed (open shell, charged cell; 0 stage fallbacks, second run bitwise) |

Timing, `-sp` with gradient, `CURCUMA_GPU_PROFILE=1`, n = 3 per variant, base and new alternating (medians):

| Cell | Method | nao / integral images | setup + pre-SCF base -> new | SCF per iteration | wall base -> new |
|---|---|---|---|---|---|
| 216 H2O, 648 atoms | gfn1 | 1728 / 125 | 484 -> 687 ms | 268 -> 268 ms (14 it) | 4.74 -> 4.93 s |
| 216 H2O, 648 atoms | gfn2 | 1296 / 27 | 419 -> 332 ms (post-SCF E 197 -> 170 ms) | 170 -> 170 ms (17 it) | 3.71 -> 3.61 s |
| caffeine 14 A, 24 atoms | gfn1 | 76 | 99 -> 100 ms | unchanged | 0.386 -> 0.385 s |
| caffeine 14 A, 24 atoms | gfn2 | 66 | 95 -> 100 ms | unchanged | 0.391 -> 0.398 s |

Per SCF iteration nothing changes (the integrals are built once per geometry). GFN1 on the 216 H2O box is slower: a
one-off trace (n = 1; temporary stream-synchronised timers in `computeIntegralsPeriodic`, removed again) gave `k_overlap_h0_pbc` 497 ms for 1728^2 AO pairs x 125 images and the device Cholesky 97 ms,
against 272 ms for the host build on 4 threads. GFN2 (27 images) gains from the skipped host multipole integral build
(kernel 53 ms, multipole 23 ms, Cholesky 85 ms). The kernels are one thread per AO pair with the full image loop, as on
CUDA; on this GPU the FP64 rate makes that slower than the host for many images.

## Not tested

- Any other GPU than one gfx1200; cells larger than 216 H2O; non-cubic orthorhombic cells; elements other than H/C/N/O.
- `-gpu_sparse_integrals auto|on` with a cell (code turns the storage dense; not run), spin polarisation with a cell on
  ROCm, MD/opt/NPT with a cell (geometry updates through the new path), thread counts other than 1 and 4.
- FD of the virial on ROCm (only GPU vs CPU virial in the parity test); the device Cholesky probe on an ill-conditioned S.

## Not implemented relative to CUDA

- No device memory estimate before the build (CUDA `estimateStorageBytes`); the buffers are those of `begBasis`.
- No probe right after the device Cholesky (CUDA `probeCholeskyFactor`); the HIP `factorCheck` runs at the first solve.
- No `profMark` buckets inside the build (HIP has no profiler marks there).
- R2-R5 (resident periodic loop), A1r, G1c, E1 remain unported (unchanged from A3).

## Proposed lines (central files, for the coordinator)

- TODO.md: "ROCm periodisch: G1b (Geräte-Integrale) portiert (docs/ROCM_PARITY_G1B.md); GFN1 mit vielen Bildern (216 H2O,
  125 Bilder) auf gfx1200 langsamer als der Host-Aufbau (687 statt 484 ms Setup), Kernel pro AO-Paar; offen: Default
  für ROCm GFN1 oder Schalenpaar-Kernel." Replace the A3 line "G1b ... fehlen" by "R2-R5 (residente periodische Schleife)
  fehlen".
- AIChangelog.md: "2026-11-25 ROCm parity G1B: periodic S/H0/L + GFN2 multipole integrals on the ROCm device
  (`k_overlap_h0_pbc`, `k_multipole_ints_pbc`, `beginComputedPeriodic`); vs CPU <= 3.4e-12 Eh, 2.1e-10 Eh/A, 3.3e-10 Eh
  (10 cells), ISA of the old kernels identical -> docs/ROCM_PARITY_G1B.md"
- GPU_PARITY_CUDA_ROCM.md / `gpu_parity_known.json`: drop `k_overlap_h0_pbc`, `k_multipole_ints_pbc` from the CUDA-only
  list; seam `periodicStageOnDevice` HIP false only for `devpot`.
- GPU_PERIODIC_GAPS.md B4: "G1b closed on ROCm Nov 25, 2026 (docs/ROCM_PARITY_G1B.md); open: R2-R5".
- `qm_methods/CLAUDE.md`, ROCm periodic bullet: "`periodicStageOnDevice` false for `devpot` only; G1b device integrals
  ported (docs/ROCM_PARITY_G1B.md); not ported: R2-R5, A1r, G1c, E1".

## Plugin ABI change needed

No. `GpuScfBackend::beginComputedPeriodic` exists since ABI 14; the new `XtbHipContext` members are non-virtual functions
behind the pimpl (class layout unchanged), `PbcState` lives inside the plugin. An old plugin with the new core keeps the
host integrals (`ints` false there); a new plugin with an old core runs G1b but counts an `mpmat` fallback (correct
numbers, host matrices), because the old core lacks the `devpot` condition of R4.
