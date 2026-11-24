# ROCm parity L1: fast list rebuilds with a cell on ROCm (Nov 24, 2026)

Claude Generated. 🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT (gfx1200); human production testing pending.
Item L1 (open after [ROCM_PARITY_B2.md](ROCM_PARITY_B2.md) "Not implemented") of
[ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md). CUDA reference: commit `ae154fe3`
([WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md) section L1). Base: `feature/rocm-parity` `ebd1b357`. Nothing committed.

## What was ported

`src/core/energy_calculators/ff_methods/rocm/gfnff_rocm.hip`, all behind `-gfnff.list_fast_rebuild` (default true; the
wrapper already called `setListFastRebuild` for the HIP backend, `has_cn_cutoff` is true there) and a cell
(`impl.pbc_active`); `false` or no cell runs the former kernels.

- Float helpers (CUDA text): `hbWrapF`, `RowCellF`/`rowCellF()` (orthorhombic, L <= 400 Bohr, off-diagonal < 1e-9 Bohr,
  else off), `rowBandDecide` (r^2 outside cutoff^2 (1 +- 1e-4) decided in float, inside or NaN by the double test).
- D4 pair list: `k_disp_rows_count_pbc`, `k_disp_rows_build_pbc` (predicate `dispPairInListF`). ROCm already had the
  ordered row builds (`k_disp_rows_*`, Oct 2026) with the gather CSR on the device, so only the float band filter was
  missing; `k_disp_pairs_static`, `k_disp_csr_lower` unchanged.
- Repulsion list: `k_rep_rows_count_pbc`, `k_rep_rows_build_pbc` (predicate `repPairInListHipF`), same scheme.
- HB detection: `k_hb_hnb` (count and write pass, one thread per AB atom; stores hyd_list positions, which the ROCm
  queue uses, where CUDA stores atom indices) and `k_hb_detect_rows_pbc` (hydrogens near the donor only, float
  minimum-image pre-filter of the two distance tests, `d_hb_strength_decision_fv` on minimum-image displacements);
  host side in `FFWorkspaceHip::detectHBondsOnDevice()` (buffers `hnb_cnt/off/idx`, device scan).
- The cell-list D4 CN of the generator (host code, `CNCalculator::calculateGFNFFCNPeriodicLocal`) is shared and was
  already active with `-gpu rocm` before this item.
- Diagnostics (no numerical effect): `CURCUMA_GFNFF_LISTHASH=1` prints `d4pairs` (the CUDA columns), `d4csr` (ROCm
  gather CSR, extra), `reppairs` (as a sorted set) and `hbtriples`; `CURCUMA_GFNFF_STEPPHASE=1` prints `HBDEV`,
  `DISPBUILD`, `REPBUILD` timers (they synchronise the stream).
- Deviation from CUDA (measured reason): CUDA builds `k_hb_detect` and `k_hb_detect_pbc` from one template body. On
  gfx1200 the template bodies changed the ISA of all 5 former kernels (e.g. `k_hb_detect_rows` 3579 -> 2175
  instructions). So the former kernels keep their source verbatim and the `_pbc` kernels are separate copies. A
  second caller of `d_hb_case1_strength` alone changed its inlining into `k_hb_detect_rows` (one FMA operand order,
  `__forceinline__` did not restore it), so `k_hb_detect_rows_pbc` calls a copy `d_hb_case1_strength_pbc`. Any later
  change to these five kernels or to `d_hb_case1_strength` has to be made in both places.
- `cuda/ff_workspace_gpu.h`: comment of `setListFastRebuild` only (no layout change).

## Measured

Structure: the 1728-H2O box of B2 (5184 atoms, L = 37.25 Angstrom), NPT CSVR 300 K, Berendsen every 2 steps, 4 host
threads, `-gpu rocm -gfnff.gpu_periodic true -gpu_strict true -gfnff.frag_charge_model reference`, RX 9060 XT. Base =
HEAD `ebd1b357` built from this worktree; "old" = the new binary with `-gfnff.list_fast_rebuild false`.

| check (thresholds fixed beforehand: equality) | result |
|---|---|
| ISA of the plugin object `gfnff_rocm.hip.o` (scratch tool `a2/isa.py`), base vs new | 136/136 former functions identical, 6 new kernels |
| list hashes `CURCUMA_GFNFF_LISTHASH=1`, new vs old, 100 steps (51 D4 CN, 51 D4 lists, 51 D4 CSR, 60 HB lists, 50 repulsion lists) | 263/263 identical |
| trajectory md5 and md.csv (14 numeric columns) of those runs, new = old = base | identical (3 binaries/paths) |
| no cell: 6 x `-sp -dump_gradient` (water, caffeine, acetic-acid dimer, trisaccharide, macrocycle-bgal, peo201-chain) + 100-step NVT MD of macrocycle-bgal, base vs new, second new run | 7/7 identical, second run 7/7 identical |
| fallbacks in the cell runs | 0 (`-gpu_strict true`) |

Time per list build (`CURCUMA_GFNFF_STEPPHASE=1`, ms, median and range over the builds of a 50-step run without the
set-up build; n = 29 HB detections, 25 D4 builds, 24 repulsion builds; old = former kernels = base):

| part | old | new |
|---|---|---|
| HB detection kernel (`k_hb_detect_rows` / `k_hb_hnb` + `k_hb_detect_rows_pbc`) | 101.9 (100.2-103.2) | 2.4 + 9.9 (9.7-10.1) |
| D4 device list (count, build, static data, CSR) | 25.3 (25.0-28.7) | 17.9 (15.6-18.5) |
| repulsion device list | 12.1 (10.5-13.4) | 4.6 (4.5-5.9) |
| D4 generator CN on the host (shared code, already fast in the base) | 101.4 (88.6-158.8) with `false` | 5.6 (2.2-6.3) |

Time per MD step (md.csv `dt`, 1 ms resolution, steps 10-60 of a 60-step run, median per run, then median and range
over n = 3 runs; inside one GPU lock, base and new runs interleaved; trajectories md5-identical):

| binary | step without barostat | step with barostat (list rebuild) | mean over all steps |
|---|---|---|---|
| base `ebd1b357` | 77 (77-77) | 225 (224-226) | 162 (161-163) |
| new | 77 (77-77) | 129 (128-130) | 104 (104-105) |

An earlier session (other agents' runs on the machine, plain steps 84-87 ms) gave new 147 (146-147) against
`list_fast_rebuild false` 373 (369-383) ms per barostat step; that `false` path also switches off the host cell-list
D4 CN, so it is slower than the base.

ctests (`ctest -j1` under the GPU lock): `gfnff_rocm_pbc_refuse`, `_terms`, `_hybrid`, `_virial`, `_ewald`, the 14
`structinv_rocm_*_gfnff`, `cli_gfnff_gpu_02/04/05_rocm`: 22/22 passed. `gfnff_rocm_charge_groups` not run (its script
needs numpy, absent on the machine; no cell involved).

## Not tested

- Non-orthorhombic cells (the float filters switch off; GFN-FF refuses them anyway), boxes above L = 400 Bohr, systems
  other than water (the 5184-atom box of B2), long runs, wave64 hardware, gfx1150, Windows, multi-GPU (off with a cell).
- `-gfnff.gpu_hb_detect_on_device false` and `CURCUMA_ROCM_HB_GATHER=0` / `CURCUMA_GPU_DISP_SORT=0` with a cell
  (the fast HB path is in the device detection only; the D4 fast path only in the ordered build).

## Not implemented relative to CUDA

- Nothing of L1. CUDA-only parts of L1 that ROCm did not need: the row-block replacement of the atomic/key-sort D4 and
  repulsion builds (ROCm had ordered row builds before) and `k_disp_atom_zeta` (ROCm computes zeta in
  `k_disp_pairs_static`).

## Proposed central lines (coordinator)

- TODO.md: remove "ROCm L1 fast list rebuilds with a cell (k_hb_detect_pbc, list_fast_rebuild)".
- AIChangelog.md: `2026-11-24 ROCm L1: list_fast_rebuild on ROCm (float minimum-image band filters in the D4/repulsion
  row builds, k_hb_hnb + k_hb_detect_rows_pbc), same lists (LISTHASH 263/263), barostat step 5184 atoms 225 -> 129 ms
  on gfx1200 (docs/ROCM_PARITY_L1.md)`.
- ff_methods/CLAUDE.md, L1 bullet: "CUDA row-block D4/repulsion builds and `k_hb_detect_pbc`" -> add "; ROCm since
  Nov 24, 2026: float band filters in `k_disp_rows_*_pbc`/`k_rep_rows_*_pbc`, `k_hb_hnb` + `k_hb_detect_rows_pbc`
  (separate copies, former kernels ISA-identical)".
- GPU_PARITY_CUDA_ROCM.md / GPU_PERIODIC_GAPS.md: L1 fast list rebuilds on ROCm done Nov 24, 2026 (gfx1200).
- gpu_parity_known.json: the CUDA kernel names `k_hb_detect_pbc` and `k_disp_atom_zeta` have no ROCm namesake
  (ROCm: `k_hb_detect_rows_pbc`; zeta inside `k_disp_pairs_static`); `k_disp_rows_*`, `k_rep_rows_*`, `k_hb_hnb` exist.

## Plugin ABI

No change needed: no shared class layout or entry point changed (`m_list_fast_rebuild` and `setListFastRebuild` were in
the shared header since `ae154fe3`; new buffers live in the HIP Pimpl `FFWorkspaceHipImpl`).
