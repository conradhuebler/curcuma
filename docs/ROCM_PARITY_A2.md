# ROCm parity A2: periodic GFN-FF local terms on ROCm (Nov 24, 2026)

Claude Generated. 🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT (gfx1200); human production testing pending.
Item A2 of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md) (GPU_PERIODIC_GAPS B2). CUDA reference: S1
(`8ee4f638`, WP-PERIODIC-NONBONDED S1). Base: `feature/rocm-parity` `d4dbdfc5` plus the coordinator patch (build fix,
B1, B10, plan). Nothing committed.

## What was ported

`src/core/energy_calculators/ff_methods/rocm/gfnff_rocm.hip` (single translation unit):

- `micPlace()` / `micChainOfFour()`: the CPU `placeMinimumImage()` arithmetic, copied from `cuda/gfnff_kernels.cu`.
- `k_angles`, `k_dihedrals`, `k_inversions`, `k_storsions`, `k_xbonds`, `k_hbonds` are `template <bool PBC>`;
  `<true>` builds the image frame (HB frame on H, neighbours of B and the acceptor parent C next to B, D next to C;
  XB frame on X), `<false>` keeps the former statements (`if constexpr` at each divergence point, the CUDA S1 text).
  Launch through `GFNFF_HIP_PBC_LAUNCH(kernel, impl.pbc_active, ...)`; the declarations in `rocm/gfnff_kernels_hip.h`
  are templates as well (they carry the `__launch_bounds__`).
- `k_xbonds_mixed` (FP32, no frame) runs only without a cell; the XB launch follows `hbond` as on CUDA/CPU.
- Minimum image in `k_check_displacement`, `k_check_displacement_pair`, the host test `hostDisplacementExceeded()`
  (single-stream mode; host copy of the cell in `FFWorkspaceHipImpl::h_lat/h_lat_inv`), `k_disp_pairs_count/_build`,
  `dispPairInList()` (`k_disp_rows_count/_build`), `repPairInListHip()` (`k_rep_rows_count/_build`), `d_hb_dist2()`,
  `d_hb_case1_strength()` (frame on H). The FP32 pre-filters (`d_hb_strength_decision_f`, the float queue filter of
  `k_hb_detect_rows`) are skipped with a cell; the double tests decide.
- D4 tail: `k_d4_tail_W`, `k_d4_tail_atoms` (CUDA text; dE/dCN through `acc_add`, i.e. fixed point or FP64 as the
  rest of the HIP step), energy slot `E_DISP_TAIL`, `FFWorkspaceHip::setDispersionTail()`, launched in
  `launchChargeDependentAndFinish()` before the Coulomb post-processing and the CN chain rule; `disp_tail` in
  `energyComponents()`.
- `FFWorkspaceHip::setUnitCell()`: per-device cache of the uploaded state (as CUDA), status checked through
  `gfnffUploadCellHip()`; `pbc_active` flag. The CUDA trap (copy through an `extern __constant__` of another unit)
  cannot occur in the single HIP unit; the status is checked anyway. `FFWorkspaceHip::hbondCount()` added.
- Already periodic before A2 (applyMIC): CN (`k_cn_compute`, `_pairs`, chain rule), CN pair list build, dispersion
  and repulsion energy kernels, HB-CN, Coulomb pair kernels (minimum-image pair sum, S1 semantics).

Wrapper (`qm_methods/`):

- `gfnff_hip_method.h`: `has_periodic_local_terms = true`, new trait `has_periodic_hybrid = false`.
- `gfnff_gpu_method.h` (CUDA): `has_periodic_hybrid = true` (new trait, CUDA code paths unchanged).
- `gfnff_gpu_method_impl.h`: the hybrid-only `if constexpr` blocks (host Ewald set-up, `setCoulombPairOnDevice`,
  `setWantVirial`/`setHBCrossing`, `assembleHybridVirial`) are gated by `has_periodic_hybrid`; the local-term blocks
  (D4 tail upload, pair build with the evaluation radius) stay on `has_periodic_local_terms`. `setMolecule()` refuses a
  cell with `coulomb_ewald true` on a backend without the hybrid (message names the backend and `-gfnff.coulomb_ewald
  false`). ROCm with a cell therefore runs only with `-gfnff.gpu_periodic true -gfnff.coulomb_ewald false` (S1 mode:
  minimum-image Coulomb, non-periodic EEQ, `hasPeriodicEnergy()` false).

Test: `test_cases/unit/test_gfnff_gpu_pbc_terms.cpp` builds a second target `test_gfnff_rocm_pbc_terms`
(`GFNFF_PBC_TERMS_ROCM`, HIP method class), ctest `gfnff_rocm_pbc_terms` (labels `gfnff_gpu;gpu;pbc;rocm`, lock
`gpu`, 77 when no HIP device answers).

## Measurements

All on one RX 9060 XT (gfx1200), ROCm build `-DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON -DUSE_AVX2=ON`, every GPU run
under `flock /tmp/claude-1000/rocm-gpu.lock`. Base = the same worktree and options before the A2 edits (binary and
plugin saved, `gfnff_rocm.hip.o` saved for the ISA comparison).

| gate | result |
|---|---|
| ctest `gfnff_rocm_pbc_terms` (mirror of `gfnff_gpu_pbc_terms`: 8 systems x 5 states + torn-initialised, dispersion/repulsion off, supercell, FD, wrap, bitwise; thresholds of the CUDA test, 1e-11 Eh per term, 1e-12 Eh/A) | all checks passed; GPU vs CPU worst term 3.2e-13 Eh (nonbonded_rep; disp_tail 2.2e-13, bond 7.1e-14, HB 3.1e-14, XB 2.6e-15, sTors 0: no structure has sTors), worst gradient 1.5e-13 Eh/A; isolated dispersion/repulsion gradients <= 9.4e-16 Eh/A; HB triple counts GPU = CPU (caffeine 358, composite 550, 60 waters 1468, the CUDA numbers) |
| same test: supercell 2x1x1 (7 neutral systems) | every term = 2 x cell within 5.7e-14 Eh, per-atom gradient 3.3e-14 Eh/A, HB counts doubled |
| same test: FD of the GPU energy on torn geometries (h 2e-4 / 1e-4 A) | caffeine 8.5e-8 -> 2.1e-8, Richardson 3.2e-10; composite 3.1e-8 -> 7.9e-9, 3.7e-10; 40 waters 9.8e-8 -> 2.4e-8, 1.7e-9 Eh/A (limits 5e-8 / 2e-9) |
| same test: wrap of one molecule by +L | 0 D4 and 0 repulsion rebuilds, dE 3.6e-15 Eh, dG 2.5e-15 Eh/A; a 1.5 A move rebuilds (1) |
| same test: two device instances | energy and gradient bit for bit (60 waters, composite) |
| no cell, 6 structures (water, caffeine, acetic acid dimer, trisaccharide, macrocycle-bgal, peo201-chain; `-sp -dump_gradient -gpu rocm -gpu_strict true -threads 1`), base vs A2 | 6/6 md5-identical; second A2 run 6/6 md5-identical; base run twice 6/6 identical |
| no cell, MD 200 fs (caffeine, trisaccharide; seed 42, `-gpu rocm -gpu_strict true`), base vs A2 vs A2 again | trajectory file and final state md5-identical (3/3 per structure) |
| ISA (`llvm-objdump` of the gfx1200 code object, pc-relative literals and trailing padding normalised), base vs A2 | 100 of 110 base kernels instruction-identical, among them `k_angles/k_dihedrals/k_inversions/k_storsions/k_xbonds/k_hbonds<false>` against the former kernels; the 10 that differ are the decision kernels with an added minimum-image branch (`k_check_displacement(_pair)`, `k_disp_pairs_count/_build`, `k_disp_rows_count/_build`, `k_rep_rows_count/_build`, `k_hb_ab_rows`, `k_hb_detect_rows`); 8 new kernels (six `<true>`, two tail) |
| CLI, caffeine in a 25 A cube | `-gpu rocm` without opt-in: counted CPU fallback, gradient md5 = CPU run; `-gfnff.gpu_periodic true` (coulomb_ewald true): error message naming ROCm and `coulomb_ewald false`, exit 1 (no CPU fallback); `gpu_periodic true coulomb_ewald false`: GPU vs CPU dE 1.5e-14 Eh, max dG 1.5e-13 Eh/A, two GPU runs md5-identical, no fallback (`-gpu_strict true`) |
| ctest `cli_gfnff_gpu_0*_rocm` (3) + `gfnff_rocm_pbc_terms` | 4/4 pass |
| `scripts/gpu_parity_check.py` | `has_periodic_local_terms` equal; one unexplained difference: the new `has_periodic_hybrid` (CUDA true, ROCm false) |

Not run: full `ctest -j1` (only `curcuma`, the plugin and the new test were built, per the plan's build rule),
`ctest -L structinv`, `-threads 4` probe.

## Not tested

- CUDA build of the wrapper change (no CUDA device or toolkit here): the CUDA trait `has_periodic_hybrid` is `true`, so
  every CUDA `if constexpr` resolves as before; compiled only in the ROCm build.
- ROCm with `rocm_fixed_point false`, mixed precision, multi-stream mode above 500 atoms with a cell (the test boxes
  have at most 180 atoms), sTors and case-4 H-bonds under a cell (no structure), `topology_mode auto`/react with a
  cell, MD with a cell on ROCm, gfx1150, Windows.

## Not implemented relative to CUDA (needed by wave 2)

B1 (periodic hybrid, S2/S3):
- `has_periodic_hybrid` true for ROCm and the wrapper paths it unlocks: `FFWorkspaceHip::setCoulombPairOnDevice`
  is the header inline (exists), but the HIP unit never reads `m_coulomb_pair_on_device` (0 uses; CUDA 3):
  `launchChargeDependentAndFinish()` must skip the pair Coulomb kernels (`launchCoulombTilesHip`, `k_coulomb_gather`)
  when it is false and keep the self term and TERM 1b (`k_coulomb_postprocess`), or the hybrid counts Coulomb twice.
- Device Ewald: a HIP `curcuma::ewald::EwaldDeviceOps` (`prepare`, `potentials`, `evaluate`, `lastError`,
  `lastPrepareMs`; ewald.h) mirroring `GfnffEwaldGPU` (`cuda/ewald_gpu.cu`: `k_ew_count`, `k_ew_fill`, `k_ew_kernel`,
  `k_ew_phase`, `k_ew_kitems`, `k_ew_kpass`, `k_ew_sf`, `k_ew_real_pot`, `k_ew_real_eg`, `k_ew_sum_items`,
  `k_ew_wrapf`) and `GFNFFRocmBackend::createEwaldDevice()` returning it (today nullptr: the wrapper then counts a
  fallback and keeps the host Ewald, which is enough for a first hybrid).
- ctests to mirror: `gfnff_gpu_pbc_hybrid`, `gfnff_gpu_pbc_ewald`, `gfnff_gpu_pbc_refuse` (ROCm variant).

B2 (pressure, NPT):
- Virial: `FFWorkspaceHip::setWantVirial`, `downloadVirial` (declared in the shared header, not defined in the HIP unit:
  a ROCm build that calls them fails at dlopen), the `*_body<PBC, VIR>` split of the six local-term kernels with
  `k_angles_vir`, `k_dihedrals_vir`, `k_inversions_vir`, `k_storsions_vir`, `k_xbonds_vir`, `k_hbonds_vir`, `k_batm_vir`,
  and `k_bonds_virial`, `k_dispersion_virial`, `k_repulsion_virial`, `k_cn_chainrule_virial`, `k_hb_alpha_virial`,
  slots `cuda/gfnff_virial_slots.h`, the tail C8 energy `GFNFF_VIR_TAIL_E8`. Trap A1-GPU: new kernel parameters go into
  new kernels around a shared body, never into the `<false>` kernels (they are ISA-checked against this state).
- Crossing: `k_hb_crossing`, `k_xb_crossing`, `FFWorkspaceHip` side of `setHBCrossing` (inline in the header, needs the
  device passes).
- Cell update: the wrapper's `followHostCellChanges()` is shared and calls `uploadDispersionTailToDevice()` (works on
  ROCm since A2) and `invalidateGraph()`; the HIP CN-derivative pair list rebuild after a cell change is untested.
- L1 fast rebuilds: `k_hb_detect_pbc`, `k_hb_hnb` (neighbour-limited HB detection with float minimum-image filters),
  `list_fast_rebuild` row-block builds; on ROCm the D4/repulsion builds are row builds already (`k_disp_rows_*`,
  `k_rep_rows_*`), the HB detection with a cell queues every hydrogen (no float filter): correct, slower.

## Proposed central lines (coordinator)

- TODO.md: remove "ROCm: periodic GFN-FF local terms (B2)"; add "ROCm periodic hybrid (B1) and virial (B2): see
  docs/ROCM_PARITY_A2.md section 'Not implemented'".
- AIChangelog.md: `2026-11-24 ROCm A2: periodic GFN-FF local terms on ROCm (template<bool PBC> kernels, MIC builders,
  device HB detection, D4 tail), -gfnff.gpu_periodic true with coulomb_ewald false; hybrid refused (docs/ROCM_PARITY_A2.md)`.
- GPU_PERIODIC_GAPS.md B2: "local terms done Nov 24, 2026 (S1 mode, gfx1200); hybrid = B1 of the parity plan".
- GPU_PARITY_CUDA_ROCM.md / gpu_parity_known.json: `has_periodic_local_terms` now equal; new trait
  `has_periodic_hybrid` differs (CUDA true, ROCm false: reason "B1/B2 of ROCM_PARITY_PLAN_2026-11").
- ff_methods/CLAUDE.md "Periodic cell, CUDA (S1)": add "ROCm since Nov 24, 2026 (A2), S1 mode only".
- `MethodFactory::periodicCellNeedsCpu` warning says "CUDA only"; with A2 the opt-in also reaches ROCm (S1 mode).
  Decision for the coordinator: ROCm + `gpu_periodic true` + `coulomb_ewald true` now ends the run with the wrapper's
  error (exit 1, as on CUDA before S2); a counted CPU fallback instead would need a backend check in
  `periodicCellNeedsCpu()` (core file, not changed here). B1 removes the case.

## Plugin ABI

No change needed: no shared class layout changed (`FFWorkspaceGPUImpl`/`FFWorkspaceHipImpl` are plugin-private Pimpls;
the wrapper template and its traits are compiled into each plugin), no `ComputationalMethod` or entry point change.
