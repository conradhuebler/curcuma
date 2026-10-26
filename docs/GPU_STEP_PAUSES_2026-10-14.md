# CUDA GFN-FF MD: periodic step pauses (Oct 14, 2026)

Status: AI-generated changes, machine-measured on one GPU. Not human-tested. Branch `feature/gpu-step-phase-timing`.
Machine: RTX 5090 (32 GB), 24 threads, run with `-threads 12`. Origin: a GUI client saw single MD steps of several hundred
ms at a regular rhythm on a CUDA server (6200 atoms); this note reproduces and reduces them without the client.

## Where the long steps came from

`curcuma -md urea400-water1000.xyz -method gfnff -gpu cuda -dt 1 -md_diagnostics true -md_diagnostics_timing true`, 200 to 250 steps,
single runs (n = 1) unless a count is given. The median step is 15.4 ms; the mean was 36.3 ms with defaults, 22 steps longer than 100 ms.

| Pause | Interval | Cost before | Phase field |
|---|---|---|---|
| HB/XB list re-detection (`hb_update_force_every` 10) | every 10 steps | 99 to 165 ms | `gpu_hbxb_update` |
| D4 pair list rebuild (displacement > half skin) | about every 80 to 100 steps | 169 to 217 ms | `gpu_disp_update` |
| Non-bonded repulsion list rebuild (skin 2 Bohr) | every 6 to 8 steps | 19 to 25 ms | `gpu_rep_update` |

`hb_update_force_every 0` removes the 10-step rhythm (0 of 200 steps with `gpu_hbxb_update` > 50 ms), the other two stay.
New fields in the `timing_ms` block of the diagnostics file, all host wall clock: `gpu_pre`, `gpu_rep_update`, `gpu_rep_upload`,
`gpu_rep_pairs` (count), `gpu_disp_update`, `gpu_coul_update`, `gpu_hb_detect`, `gpu_hb_rebuild`, `gpu_hb_upload`. With them
the named phases add up to `total` for every step above twice the median.

## Changes

- **HB/XB refresh** (no change to the detection result): the GPU wrapper no longer copies the lists into the CPU workspace
  (`setSkipHostHBWorkspaceUpdate`, every backend now); the device upload of the H-bond list uses the staged pinned upload
  (`stagedHBondUpload`, `-gfnff.gpu_staged_list_upload`); the offsets of the packed neighbour arrays and the bond-HB cross
  reference (`rebuildBondHBData`) are computed in parallel chunks joined in list order. 98.8 -> 35.4 ms per re-detection
  (median of 11): detection 50 -> 22, `rebuildBondHBData` 9.4 -> 2.5, upload 37.9 -> 10.7.
- **Repulsion list**: the storage of the previous list and the per-chunk buffers are reused, the chunks are joined in parallel,
  and the CPU workspace copy is skipped (`setSkipHostRepulsionWorkspaceUpdate`). The list is the serial one element for element.
  `gpu_rep_update` 21.9 -> 12.5 ms (median of 33), CPU runs bit-identical (60 MD steps, all energy terms).
- **`-gfnff.gpu_disp_pairs_on_device`** (opt-in, default still false), two defects fixed: the setup build ran before the first
  `setGeometry()`, so the 60 Bohr cutoff kept every pair (19.2 M instead of 13.7 M, +1.7 ms per step); the pair order depended on
  the arrival at the atomic counter, now sorted by `i*N+j` (host order; switch `CURCUMA_GPU_DISP_SORT=0` keeps the arrival
  order). Rebuild 169 to 217 -> 23 ms, median step 15.3 -> 15.4 ms (default off: 15.4).

Result of the host-side changes alone (200 steps): median 15.4 -> 15.4 ms, mean 36.3 -> 23.2 ms, steps > 100 ms 22 -> 2 (the two D4 rebuilds); with the two device switches (now the CUDA default) see below.
What remains per step (mean, estimated from the intervals above): HB/XB about 3.5 ms, repulsion 1.8 ms, D4 1.9 ms (0.2 with the switch).

## Device-side H-bond detection (`-gfnff.gpu_hb_detect_on_device`, CUDA, opt-in)

The (A, H, B) triples are enumerated on the device instead of by the host detection plus upload: `k_hb_ab_pairs` compacts the
AB pairs, `k_hb_detect` (one warp per AB pair, lanes over the HB hydrogens, float pre-filter, three-valued float strength
decision with the double formula for the 1 % band) emits the triples, they are sorted by (A, H, B) and `k_hb_fill_columns` /
`k_hb_fill_neighbors` write the H-bond SoA. Only (A, H, B) go back to the host for `rebuildBondHBDataFromTriples()`. Halogen
bonds stay on the host. Host inputs: `GFNFF::prepareHBDetectInputs()` (a copy of the O(N) part of the host detection, kept in
step by the verification mode), trigger `GFNFF::hbxbRefreshDue()`. Case 1 triples have A < B by atom index (the host orientation
follows its cell order, the energy is symmetric); below `nb_cell_list_min_atoms` the pruning of weak case-1 triples follows the
host rule (only its parallel legacy branch prunes).

- **Same list as the host**: `-gfnff.gpu_hb_detect_verify true` runs the host detection at every refresh and prints the set
  difference (`HBVERIFY` lines). 21 structures from 3 to 7320 atoms, 4 refreshes each, and 31 refreshes along 300 MD steps of
  urea400-water1000, peo201x2-water1500, acetic-acid dimer and pep20: 0 differences, 0 duplicates (582k to 605k triples on urea).
- **Time per re-detection** on urea400-water1000 (median of 19, RTX 5090, 12 threads): `gpu_hbxb_update` 36.2 -> 8.1 ms
  (device enumeration 2.3 ms, sort 0.3, SoA 1, download 1.7, host inputs 0.7, `rebuildBondHBData` 2.0); mean MD step 23.2 -> 19.9 ms,
  median 15.4 -> 14.9 ms. peo201x2-water1500: 6.8 -> 3.1 ms. Step 100 energy equal to 1e-10 Eh.
- **Single points** (16 structures, 12 to 7320 atoms): device list against host list ΔE <= 4.6e-13 Eh, Δgradient <= 5.6e-16 Eh/A.
- First version of the kernel took 78 ms (one thread per pair, 0.4 % of the 19 M threads active, each looping over 3600
  hydrogens, and double-precision `pow`/`exp`/division in the strength test on a GPU with 1/64 double rate).

Not checked: gradients over long runs, `-opt` beyond caffeine (38 steps, trajectory equal) and a 399-step unconverged cluster
optimisation (scatter between identical runs larger than the on/off difference), charged systems with several fragments (CPU
fallback), memory of the device buffers (about 20 B per candidate entry plus the sort), MOR41/GMTKN55. ROCm: own section below.

## `gpu_disp_pairs_on_device` checked further (same day)

27 structures from 2 to 7320 atoms (`-sp`, energy and gradient, off / on / CPU): on against off ΔE <= 3e-13 Eh, Δgradient <=
5e-16 Eh/A. MD (2 alternating runs each): median step equal for caffeine (0.35 ms), trisaccharide (0.71), urea cluster (1.86),
peo201-chain (2.9), peo201-water500 (5.5), peo201x2-water1500 (16.0 ms); mean and longest step lower with the switch on the three
bulk systems (peo201x2-water1500: mean 19.9 -> 18.4 ms, longest step 174 -> 36 ms). Energies in step 50 equal on 9 digits.
Both `gpu_disp_pairs_on_device` and `gpu_hb_detect_on_device` default to **true on CUDA since Oct 14, 2026** (operator decision after
these checks). Both switches are mirrored on ROCm (next two sections). With the defaults
(urea400-water1000, 200 steps, n = 1): mean step 18.1 ms (36.3 ms at the start of the session), longest step 48 ms (234 ms), median
14.9 ms; each switch set to `false` restores its former cost (H-bond refresh 8.0 -> 35.9 ms, D4 rebuild 23 -> 167-218 ms), E(step 100) equal in
all four combinations.

## ROCm mirror of the device D4 build (same day, RX 9060 XT, `-threads 10`)

Before, `gfnff_rocm.hip` built the list with an atomic slot (arrival order) and then rebuilt the per-atom gather CSR on the
host. Phase split of one rebuild on urea400-water1000 (13.7 M pairs, 3 rebuilds, ms): count 3.9, build 21-35 (incl. nine
`hipFree`/`hipMalloc`), index download 13-26, host CSR loops 123-137, CSR upload 10; host CN and Gaussian weights 16.
Now (`FFWorkspaceHip::generateDispersionPairListOnGPU`): one workgroup per atom counts its partners below and above its index
(`k_disp_rows_count`), the host scans the two N-vectors, `k_disp_rows_build` writes the pairs in ascending partner order
(host order, key i*N+j, no sort needed) and the CSR rows, `k_disp_pairs_static` contracts C6 per pair, `k_disp_csr_lower` finds
the pair index of the lower half by binary search. Buffers grow only (`ensure`), no index download.
`CURCUMA_GPU_DISP_SORT=0` keeps the former path. Default on for ROCm as well; `-gfnff.gpu_disp_pairs_on_device false` is the way off.

- Same list: a temporary comparison against the former kernel, sorted on the host, and against the host CSR loop gave 0
  differences in pairs, order, the five per-pair values (bitwise), CSR and offsets for 16 single points (2 to 7320 atoms) and
  9 MD rebuilds (urea400-water1000, peo201x2-water1500).
- `-sp` on 16 structures from 2 to 7320 atoms (on / off / `-gpu none`): on against off ΔE <= 1.1e-13 Eh, Δgradient <=
  1.7e-16 Eh/A; on against CPU ΔE <= 1.8e-11 Eh, Δgradient <= 1.6e-13 Eh/A (off against CPU the same).
- 100 MD steps on / off: E(last step) differs by 3.4e-13 Eh on urea400-water1000 and on peo201x2-water1500; last frame equal at
  the written precision (1e-6 A).
- 200 MD steps, 3 alternating runs per binary, default before (host build) / after (device build): D4 rebuild median 427 -> 47 ms
  (max 464 -> 58, n = 6) on urea400-water1000 and 306 -> 42 ms (max 340 -> 42, n = 9) on peo201x2-water1500; median step
  unchanged (36.2 and 36.6 ms), mean step 52.8 -> 49.0 and 47.0 -> 42.6 ms. Of the 47 ms, 16 are the host CN and Gaussian
  weights, about 13 the C6 contraction.

## ROCm mirror of the device H-bond detection (same day, RX 9060 XT gfx1200, `-threads 10`)

`FFWorkspaceHip::detectHBondsOnDevice()` (`gfnff_rocm.hip`, trait `has_device_hb_detect` true, default on,
`-gfnff.gpu_hb_detect_on_device false` is the way off). Same tests as the CUDA kernels; the build links no rocPRIM/hipCUB, so
the order is made without a key sort: `k_hb_detect_rows` counts every triple into its row (case, A) while emitting it, a device
scan gives the row offsets, the triples are scattered into their rows and `k_hb_row_sort` (one workgroup per row, bitonic
in shared memory sized by the longest row) orders each row by (H, B). The result is the list sorted by (A, H, B) and grouped
by case type, which is what `fillHBondSoA()` with `hbCasePermutation()` uploads; the gather-slot tables (`buildHBSlotTables()`)
are built on the device the same way (rows = atoms, key = slot index). AB pairs are enumerated per row (`k_hb_ab_rows`), the
case-1 test runs one workgroup per AB pair with a shared-memory queue of the float pre-filter survivors. No wavefront
intrinsic is used (only `__syncthreads` and atomics), so wave32 and wave64 need no separate path. `CURCUMA_ROCM_HB_GATHER=0`
gives rows by A only and no slot tables, as in `updateHBonds()`.

- **Same list as the host** (`-gfnff.gpu_hb_detect_verify true`): 28 structures from 3 to 7320 atoms (halogens: fi-pyridine,
  na2cl2; charged: caffeine +1, hydronium-water2 +1 with `-gfnff.frag_charge_model reference`), 5 refreshes each, and 31
  refreshes along 300 MD steps of acetic-acid dimer, pep20, urea400-water1000 (582k to 605k triples) and peo201x2-water1500:
  0 differences, 0 duplicates in 264 refreshes.
- **Same device arrays**: a temporary comparison of every uploaded array (SoA columns, neighbour lists, slot offsets, per-atom
  slot lists) against `fillHBondSoA()` of the same triples sorted on the host: 0 differences in 239 lists, also with the
  shared-memory cap of the row sort forced to 16 entries (global-memory path, 35 lists) and with `CURCUMA_ROCM_HB_GATHER=0`
  (30 lists).
- **Time per refresh**, 3 alternating runs of 200 MD steps per binary (n = 60 refreshes each), median before -> after:
  urea400-water1000 `gpu_hbxb_update` 48.1 -> 13.0 ms (detection 28.1 -> 9.5, upload 17.4 -> 1.0, `rebuildBondHBData`
  2.4 -> 2.2); peo201x2-water1500 8.9 -> 3.9 ms. MD step (n = 597): median 43.16 -> 43.10 and 44.85 -> 44.81 ms, mean
  55.86 -> 51.99 and 50.90 -> 50.11 ms, longest 129.8 -> 99.9 and 86.7 -> 87.4 ms. Energy of step 100: before and after differ
  by at most 1.1e-12 Eh (three runs each).
- `rebuildBondHBDataCore()` now uses 16 chunks per thread instead of 4: in the case-grouped order all hits (cases 2-4) are in
  the last quarter of the list and were scheduled last (3.1-4.0 ms instead of 2.0-2.5 ms at 590k triples); the joined result
  does not depend on the chunk count.
- **Single points** (16 structures, 2 to 7320 atoms, on / off / `-gpu none`): on against off dE <= 3.4e-13 Eh, dgradient <=
  2.2e-16 Eh/A; on against CPU dE <= 1.8e-11 Eh, dgradient <= 1.6e-13 Eh/A (off against CPU the same).
- **Device memory** on urea400-water1000: candidate buffers 22.5 MB (five ints per entry, capacity 2^20, plus the order),
  slot scratch 12.9 MB (2.4 M slots), rows 0.3 MB; pinned host stage of the triple download about 8 MB.
- `ctest -R "gpu|rocm|gfnff" -j1`: 110 tests, 109 passed, `cli_simplemd_14_gfnff_react_h_recombination` failed (fails before
  this change as well).

Not checked: wave64 hardware and gfx1150 (only gfx1200 was run), trajectories beyond 300 steps, `-opt`, MOR41/GMTKN55,
charged systems with several fragments under the ensemble charge model (CPU fallback), the CUDA build of the shared changes
(PARAM text, `rebuildBondHBDataCore()` chunk count; not compiled here).

## ROCm repulsion list upload (same day, RX 9060 XT gfx1200, `-threads 10`)

Phase split of the old upload (`gpu_rep_upload`, 25.4 ms median on urea400-water1000, 1.47 M non-bonded pairs; temporary
timers, n = 28 refreshes, medians): SoA conversion 9.1 ms (fresh host vectors 4.4, serial pack 3.3, pageable copy 1.4, device
reallocation 0 in the median), gather CSR 13.3 ms (degree count 4.6, fresh vectors 1.8, serial scatter 6.1, copy 0.8), bonded
list 0.4 ms. The copies were not the cost; the serial host work into freshly allocated vectors was.
`FFWorkspaceHip::updateRepulsion()` now calls `stagedRepulsionUploadHip()` when `-gfnff.gpu_staged_list_upload` is true (the
default; trait `has_staged_list_upload` is now true on ROCm, where the switch selects this path only): both SoA lists and the
gather CSR are written by the `-threads` host threads into one pinned block kept across refreshes (1/8 headroom, about 79 MB at
1.47 M pairs by the layout), each thread range counts its own atom degrees and places its edges after those of the lower ranges
(row order as in the serial build), copies are `hipMemcpyAsync` on the main stream followed by one synchronisation, device
buffers grow with 1/8 headroom. `false`, or no pinned memory, keeps the serial path.

- **Device arrays** (temporary byte comparison against the serial build, removed afterwards): 0 differences in all five SoA
  columns of both lists, `n`, and the CSR (`off`, `partner`, `pairidx`, `N`) over 15 + 15 refreshes on urea400-water1000
  (`-threads` 10 and 3), 13 + 13 on peo201x2-water1500, 9 with `CURCUMA_ROCM_REP_GATHER=0`, 23 on peo201-chain.
- **Time**, 3 alternating runs of 200 MD steps per binary against a copy of the previous build, median (max) before -> after:
  urea400-water1000 (n = 81 refreshes) `gpu_rep_update` 35.7 (54.5) -> 15.1 (28.1) ms, `gpu_rep_upload` 25.4 (28.6) -> 4.6
  (13.0) ms; peo201x2-water1500 (n = 75, 0.92 M pairs) 19.3 (27.6) -> 9.0 (18.5) ms and 12.9 (17.5) -> 3.0 (9.3) ms. The maximum
  is the first refresh of each run (pinned and device allocation); later ones stay at or below 5.3 ms. MD step (n = 600): median
  43.21 -> 43.26 and 45.08 -> 45.12 ms (the refresh steps are not the median), mean 52.23 -> 49.56 and 50.43 -> 49.21 ms.
- **Energies**: Epot of step 100 over the 6 runs per system within 8e-13 Eh (urea) and 1e-13 Eh (peo); `-sp` on 12 structures
  from 2 to 7320 atoms against `-gpu none`: dE <= 1.8e-11 Eh, dgradient <= 1.5e-13 Eh/A (a single point does not refresh the
  list, so this checks only that nothing else moved).
- `ctest -R "gpu|rocm|gfnff" -j1`: 110 tests, 109 passed, `cli_simplemd_14_gfnff_react_h_recombination` failed (fails before as
  well, see above).

Not checked: the CUDA build of the shared changes (trait removal `lean_hb_refresh`, comments; not compiled here), wave64 and
gfx1150, runs longer than 200 steps, `-threads` other than 3 and 10 for the time, pinned memory exhaustion (the fallback path
was not forced), memory use beyond the layout estimate. The remaining 10.5 ms of the urea400-water1000 refresh are the host list
generation in `gfnff_method.cpp`, unchanged here.

## Checked

- Energies in steps 0, 10 and 100 equal to the previous binary within 1e-12 Eh for every change above (GPU runs are not bitwise
  reproducible, run-to-run scatter about 1e-13 Eh).
- `ctest -R "gpu|gfnff|simplemd"`: 164 tests, 160 passed, 1 skipped, 3 failed (`cli_gfnff_06_eeq_solver_control`,
  `cli_gfnff_07_eeq_localinv`, `cli_simplemd_14_gfnff_react_h_recombination`); the same three fail on the unchanged master.
- No ctest exercises `gpu_disp_pairs_on_device`.

## Not checked

- Gradients and long trajectories (energy only, 250 steps); other systems than urea400-water1000 (polymer 7320 atoms, small,
  charged ones); several runs for the scatter of the medians.
- The device HB list was not compared byte for byte (equal energies only); memory use (pinned stage 2 x 32 MB, repulsion chunk
  cache about 70 MB at 6200 atoms, the D4 sort needs about 275 MB temporarily).
- ROCm: the HB/XB and repulsion changes reach it through the shared wrapper only.
- MOR41/GMTKN55 regression (`scripts/refset_regression.py`): detection code is unchanged, not run.

## Not implemented

- D4 and repulsion rebuild stay on the step that triggers them (no overlap with device work).

Human production testing is pending.
