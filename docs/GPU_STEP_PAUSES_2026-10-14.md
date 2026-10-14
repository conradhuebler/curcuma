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

Result with default parameters (200 steps): median 15.4 -> 15.4 ms, mean 36.3 -> 23.2 ms, steps > 100 ms 22 -> 2 (the two D4 rebuilds).
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
fallback), memory of the device buffers (about 20 B per candidate entry plus the sort), ROCm (not compilable here; trait
`has_device_hb_detect` is false there), MOR41/GMTKN55.

## `gpu_disp_pairs_on_device` checked further (same day)

27 structures from 2 to 7320 atoms (`-sp`, energy and gradient, off / on / CPU): on against off ΔE <= 3e-13 Eh, Δgradient <=
5e-16 Eh/A. MD (2 alternating runs each): median step equal for caffeine (0.35 ms), trisaccharide (0.71), urea cluster (1.86),
peo201-chain (2.9), peo201-water500 (5.5), peo201x2-water1500 (16.0 ms); mean and longest step lower with the switch on the three
bulk systems (peo201x2-water1500: mean 19.9 -> 18.4 ms, longest step 174 -> 36 ms). Energies in step 50 equal on 9 digits.
The default stays `false` until the ROCm side has the same setup-order and sort fixes.

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
- ROCm (not compilable here): `gfnff_rocm.hip` has the same unsorted device D4 build and the same setup-before-`setGeometry`
  order; the HB/XB and repulsion changes reach it through the shared wrapper only.
- MOR41/GMTKN55 regression (`scripts/refset_regression.py`): detection code is unchanged, not run.

## Not implemented

- ROCm mirror of the device D4 fixes and of the device H-bond detection.
- `gpu_disp_pairs_on_device` as default; D4 and repulsion rebuild stay on the step that triggers them (no overlap with device work).

Human production testing is pending.
