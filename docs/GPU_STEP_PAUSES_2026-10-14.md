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

- HB/XB detection kernel on the device (the 35 ms per re-detection would shrink to the upload; large, must return the host list).
- `gpu_disp_pairs_on_device` as default; D4 and repulsion rebuild stay on the step that triggers them (no overlap with device work).

Human production testing is pending.
