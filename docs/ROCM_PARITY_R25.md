# ROCm parity R25: resident periodic potential and SCF loop on ROCm (Nov 24-25, 2026)

Port of the CUDA WP-PERIODIC-XTB items R2-R5 ([WP-PERIODIC-XTB.md](WP-PERIODIC-XTB.md)) to ROCm, closing the
"R2-R5" line of [ROCM_PARITY_A3.md](ROCM_PARITY_A3.md) "Not implemented relative to CUDA".
🤖 AI-generated, ⚙️ machine-tested on one AMD RX 9060 XT (gfx1200, ROCm 7.2) only; human production testing pending.

## State before

ROCm ran every periodic SCF on the host-driven device loop: `HipScfBackend::periodicStageOnDevice("devpot")` returned
false (host potential, not counted), `ewaldGammaPeriodic` / `multipoleMatricesPeriodic` kept no device copy, so
`periodicGammaResident()` was false and the tokens `resloop`/`potdl`/`mpdl` were unreachable.

## What was ported (file, function)

- `rocm/xtb_hip_periodic_resident.hiph` (new, included at the end of `xtb_hip_context.hip`): `k_mp_pbc_self_r25` (CUDA
  `k_mp_pbc_self`), `XtbHipContext::pbcKeepGammaResident`, `pbcKeepMultipoleResident` (private), `periodicGammaResident`,
  `beginPotentialPeriodic`, `residentLoopDims`, `downloadConvergedPotentials`.
- `rocm/xtb_hip_periodic_host.hiph`: `ewaldGammaPeriodic` keeps the periodic gamma in the resident `dgamma` of the basis
  (CUDA: `dGamma`); `multipoleMatricesPeriodic` adds the self terms on the device and keeps the 15 blocks of the
  `k_mp_potential_stored` layout in `dMpAmat` while the 18 blocks stay below 2 GB (CUDA limit), and accepts
  `out18 == nullptr` (R4, no host copy).
- `rocm/xtb_hip_context.hip` (`Impl`): flag `pot_pbc` (cleared by `potBeginCommon`, set by `beginPotentialPeriodic`);
  `potBuildCore` takes the periodic D4 gather `pbc_dev::k_d4_dedq_pbc_b` with `pot_pbc`; every cell-free writer of
  `dMpAmat` (`potBegin`, `potBeginStored`, `potBeginOtf`, `mpPotUpload`, `freeMpPot`) clears `pbc.mp_pbc_resident`,
  `beginBasis`/`beginComputed` clear `pbc.gamma_pbc_resident`.
- `xtb_hip_method.cpp`: "BACKEND-SPECIFIC 13" overrides `beginPotentialPeriodic`, `periodicGammaResident`,
  `downloadConvergedPotentials` (CUDA wording; not routed through `notePotentialBegin`, so a periodic failure does not
  switch off the cell-free device potential of later calculations); `periodicStageOnDevice("devpot")` true.
- No host change: `xtb_native.cpp`/`xtb_periodic.cpp` already drive R2-R5 through the `GpuScfBackend` virtuals (ABI 22/23).
- `rocm/xtb_hip_periodic_kernels.hiph`: comments of `PbcState::gamma_pbc_resident`/`mp_pbc_resident`.
- Tests: `test_pbc_xtb_gpu_fallback` ROCm GFN2 token list gained `devpot`, `resloop`, `resstep`, `potdl`; header
  comments of that test and of `test_pbc_xtb_gpu_parity`.

Form differences to CUDA, same expressions: multipole potential and energy by the existing block-per-atom kernels
`k_mp_potential_stored`/`k_mp_energy_stored` on the 15-block layout (dd packed, the periodic builder is exactly
symmetric in a, b, as the CUDA `_cw` kernels assume) instead of new warp-per-atom kernels.

## What was measured

All runs: `release/` of this worktree (`USE_ROCM`, MKL, AVX2), one RX 9060 XT, every GPU job alone under
`/tmp/claude-1000/rocm-gpu.lock`, `-scf_threshold 1e-9`, 4 threads, `-gpu_strict true`. Base = the same worktree at
`dd4812c4` built before the change (main binary md5-identical to the new one; only the plugin differs). Cells:
water dimer 8 A, caffeine 14 A, 216 H2O 18.626 A (`water216-box.sc-s1`), 432 H2O 37.252 x 18.626 x 18.626 A (the
216 box doubled along x). Energies and gradients from `-sp -dump_gradient` (17 digits).

| Check | Result |
|---|---|
| Path (verbosity 2) with a cell, GFN2 | "periodic GFN2 potential ... on the device" and "fully device-resident loop" in 4/4 cells; base: host-driven loop |
| ROCm vs CPU, GFN2, 4 cells | worst \|dE\| 1.1e-11 Eh (432 H2O), max \|dG\| 2.9e-11 Eh/A |
| `pbc_xtb_gpu_parity_rocm` (10 cases, gfn1 + gfn2, incl. virial) | passed; worst \|dE\| 3.6e-12 Eh, \|dG\| 2.1e-10 Eh/A, \|dXi\| 3.3e-10 Eh; second run bitwise 10/10 |
| `pbc_xtb_gpu_fallback_rocm` | passed; 28 injected stage tokens counted (gfn2 19 incl. the new devpot, resloop, resstep, potdl; gfn1 9), each equal to its host token to 0.0 (E, G, Xi); normal runs 0 fallbacks; `-gpu_strict` exit 3 |
| `pbc_xtb_gpu_combinations_rocm`, `pbc_xtb_cell_refusals` | passed (open shell / charged cells stay on the host-driven loop) |
| Second run, new build, 8 cases (4 cells x gfn1/gfn2) + 2 later repeats | md5-identical 10/10 |
| `CURCUMA_PBC_GPU_HOST=devpot` (new) vs base, 8 cases | md5-identical 8/8 |
| GFN1 with a cell, new vs base | md5-identical 4/4 (GFN1 keeps the host-driven loop, see "Not implemented") |
| No cell, ROCm, new vs base (water dimer, caffeine, acetic-acid dimer x gfn1/gfn2 x 1/4 threads) | md5-identical 12/12 |
| ISA of `xtb_hip_context.hip.o` (gfx1200), base vs new (`a2/isa.py`) | 145/145 former functions identical, 1 new kernel (`k_mp_pbc_self_r25`) |
| FD of the assembled ROCm surface, GFN2, h 1e-4 A, 6 components each, `-scf_threshold 1e-11` | water dimer 8 A 2.3e-9 Eh/A, caffeine 14 A 1.7e-8 Eh/A (A3 measured 2.9e-8 for caffeine on the CPU as well: FD truncation) |
| NVE MD, 216 H2O gfn2, 100 steps of 0.5 fs, seed 7, T0 300 K, new vs base | all 101 status lines and the last trajectory frame identical as printed; max \|Etot - Etot0\| 5.53e-3 Eh in both (hot unrelaxed box, T 300 -> 642 K) |
| NPT MD, same plus csvr and Berendsen (`-pcoupl_interval 10 -tau_p 100`) | identical as printed in both builds (cell 18.626 -> 19.189 A), 0 fallbacks |
| GFN1 resident loop with a cell (R5) through a temporary host gate, not kept (see below) | 4/4 cells on the device loop, vs CPU worst \|dE\| 4.9e-11 Eh, \|dG\| 6.2e-11 Eh/A, second run md5-identical 4/4; SCF time not lower than the host-driven loop (216 H2O 3829 vs 3753 ms, 432 H2O 18000 vs 17886 ms, n = 1) |

Time, GFN2 single point, SCF wall time printed by the SCF (median of n = 4 runs each, range in brackets):

| Cell | iterations | base (host-driven loop) | new (resident loop) | per iteration base / new |
|---|---|---|---|---|
| caffeine 14 A (24 atoms) | 20 | 202 ms [190-213] | 230 ms [212-238] | 4.3 / 4.5 ms |
| 216 H2O (648 atoms) | 17 | 2859 ms [2830-2879] | 2820 ms [2773-2985] | 196.5 / 192.6 ms |
| 432 H2O (1296 atoms) | 17 | 11120 ms [10927-11142] | 10688 ms [10659-10700] | 836 / 815 ms |

MD, 216 H2O gfn2, 100 steps incl. gradient (n = 1 each): NVE 270.8 s base, 263.7 s new (-2.6 %); NPT 278.3 s base,
271.6 s new (-2.4 %). The rocSOLVER generalized eigensolve dominates each iteration, so moving the potential,
mixing and energies onto the device saves 2-3 % at these sizes. The first iteration of a small cell is 30-40 ms
slower (water dimer 139 vs 98 ms, n = 1; caffeine total +28 ms median); cause not analysed (candidates: first launch
of the resident-loop kernels in the process, the D4 weight-table upload).

## Not tested

- Any GPU other than one gfx1200; cells larger than 432 H2O; the 2 GB keep limit of the 18 multipole blocks (about
  3700 atoms) and the counted fallback above it; elements other than H/C/N/O; triclinic cells (refused anyway).
- `out18 == nullptr` (R4 device-only multipole matrices, token `mpdl`): reachable only with device-built periodic
  integrals (`m_pbc_mp_device_only` needs `pbc_dev_ints`, G1b), so not run.
- FD of the virial (only GPU vs CPU of the virial in `pbc_xtb_gpu_parity_rocm`).
- MD energy conservation of a relaxed system; long MD; NPT with C-rescale; `-opt` with a cell.
- The new code with `-gpu_sparse_integrals auto|on`, `-scf_precond frag`, pseudo-diagonalisation, `-gpu_multipole_otf on`
  together with a cell.
- Several calculators in one process alternating periodic and cell-free geometries beyond what the parity and
  fallback tests do (the resident flags are cleared by every cell-free writer of `dgamma`/`dMpAmat`, checked by
  reading the code only).

## Not implemented relative to CUDA

- R5 GFN1 resident loop with a cell is ported on the device side (`periodicGammaResident` true), but the host gate in
  `xtb_native.cpp` also requires the device-built integrals (`gpu_computed`, G1b), which ROCm does not have yet; GFN1
  with a cell therefore stays on the host-driven loop until G1b lands. The temporary host gate used for the
  measurement above (`gpu_computed || (m_pbc && env)`) was reverted; no host change is part of this item.
- R3 complement density (E1, `setHostLoopDensityComplement`): no HIP port; the periodic resident loop keeps the GEMM
  density (same numbers, no speed-up from E1).
- CUDA's warp-per-atom multipole kernels (`k_multipole_potential_cw`, `k_energy_multipole_cw`): ROCm uses its existing
  block-per-atom kernels on the 15-block layout (same expressions).
- Open shell and charged cells keep the host-driven loop, as on CUDA (host gate `pbc_neutral_closed`).

## Merge notes (G1B works in the same files)

- Comments that still say R2-R5 is not ported and need the merge to reconcile: `rocm/xtb_hip_context.h` block comment
  above `beginDispersionImages` ("Not ported: ... beginPotentialPeriodic, R2-R5; periodicGammaResident() is therefore
  always false"), `xtb_hip_method.cpp` BACKEND-SPECIFIC 12 header ("Not ported ... beginPotentialPeriodic and the
  resident periodic loop (R2-R5)") and the comment above `periodicStageOnDevice` ("devpot ... host by design"). Left
  unchanged here because G1B rewrites the same sentences for `ints`.
- `periodicStageOnDevice`: R25 inserted one line before the existing `return t != "ints" && t != "devpot";`; when G1B
  removes `ints` there, the whole function can become `return true` (or be deleted).
- `rocm/xtb_hip_context.hip`: the include of `xtb_hip_periodic_resident.hiph` sits right after the include of
  `xtb_hip_periodic_host.hiph`.
- If G1B's `beginComputedPeriodic` port writes the resident `dgamma` (CUDA's does not), it must run before
  `ewaldGammaPeriodic` of the same geometry or clear `pbc.gamma_pbc_resident`; `beginComputed` clears it already.

## Proposed lines (central files, for the coordinator)

- TODO.md: replace "ROCm periodisch: G1b ... und R2-R5 (residente periodische Schleife) fehlen" by "ROCm periodisch:
  G1b (Geraete-Integrale) fehlt; damit auch die residente GFN1-Schleife mit Zelle (R5 geraeteseitig portiert,
  docs/ROCM_PARITY_R25.md)."
- AIChangelog.md: "2026-11-25 ROCm parity R25: periodic GFN2 SCF device-resident on ROCm (CUDA R2-R4: periodic device
  potential, resident loop, converged potentials; R5 device side) -> docs/ROCM_PARITY_R25.md"
- GPU_PARITY_CUDA_ROCM.md / GPU_PERIODIC_GAPS.md B4: "R2-R4 ported to ROCm Nov 25, 2026 (GFN2, closed shell, neutral;
  vs CPU <= 1.1e-11 Eh, 2.9e-11 Eh/A, 4 cells; 2-3 % per SCF iteration at 648-1296 atoms); R5 waits for G1b; E1 and
  the warp-per-atom kernels not ported."
- `qm_methods/CLAUDE.md` ROCm periodic line: "`periodicStageOnDevice` false for `ints` only; R2-R4 resident GFN2 loop
  with a cell (`rocm/xtb_hip_periodic_resident.hiph`); not ported: G1b, R5 host gate (needs G1b), A1r, G1c, E1".
- ROCM_PARITY_A3.md "Not implemented relative to CUDA": strike the R2-R5 line, point to this file.

## Plugin ABI

No change needed: no new `GpuScfBackend` virtual and no shared class layout change (`XtbHipContext` is plugin-internal);
the virtuals used exist since ABI 22/23.
