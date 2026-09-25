# Second merge of `origin/feature/multi-gpu` into the reactff branch (Sep 25, 2026)

Claude Generated. 🤖 AI-merged, machine-tested (CPU only). Not pushed. Nothing on `reactff2-llm`
itself or in the main working tree was touched.

## Where it is

| | |
|---|---|
| worktree | `/home/conrad/src/curcuma_branches/curcuma/.claude/worktrees/agent-a339a20dbf71e1a21` |
| branch | `reactff2-llm-merge-multigpu` |
| merge commit | `4c39e4a16b8efab4476deeea9c697ff86d05cae2` (parents `9c69d40e` = reactff2-llm as instructed, `7d2cceb7` = origin/feature/multi-gpu after `git fetch`) |
| HEAD | this status file on top of the merge commit |

**Important — `reactff2-llm` moved while this ran.** At the start its tip was `9c69d40e`
(as the task said); by the time of the merge it is `63a3e4de`, five commits further
(`18aa7503`, `ad7f754c`, `837266d8`, `d4482d2f`, `63a3e4de` — the packages 14-32 stage-2 work).
This merge is based on `9c69d40e` and does **not** contain them. Reconciling means merging
`reactff2-llm` (63a3e4de) into this branch, or redoing the merge on the new tip. Expect a
**semantic** overlap there that git will not flag: `837266d8` adds a "Stale-CN fix B" per-step
D4 C6 refresh (`refreshC6WeightsForCN` in `prepareCNAndEEQ`), and the remote's `f51f5200` adds
its own per-step C6 refresh (`GFNFF::refreshDispersionC6`, PARAM `dispersion_c6_update`). Both
fix the same frozen-C6 defect; only one should survive.

## What was incoming

`merge-base(reactff2-llm@9c69d40e, origin/feature/multi-gpu) = 43e77c92` (the point of the
Sep 18 package-2 merge). **44 new commits**, 71 files, +28.7k/-0.5k lines (≈9k of it three
polymer_2x xyz files). By area:

- **GPU / multi-GPU (CUDA, ROCm)**: verified distributed eigensolve, optional multi-GPU library
  finder, FP32 stall guard + false-fixpoint detection, mixed precision off on full-rate-FP64
  GPUs, GPU gradient 38 -> 7 s at 7320 atoms (W build distributed), Coulomb-gradient kernel,
  GFN-FF EEQ projected-PCG on GPU (WP7-E), every tuning knob a CLI flag + `scripts/tuning_sweep.py`.
- **MD**: its own MD-clock fix (`ef462fcf`, same bug and same constant as this branch's
  `1b7e5ff0`), step-rejecting integrator `-adaptive_step` (+ local hottest-atom criterion,
  default OFF), cell-list crash guard for diverging MD.
- **GFN-FF physics / defaults (these change numbers)**:
  - `ba0319dd`/`f51f5200`: non-bonded repulsion, D4, explicit-Coulomb pair lists refreshed
    during MD; per-step D4 C6 refresh `dispersion_c6_update` **default true**.
  - `7d2cceb7`: `hh_repulsion_bpair` **default true**, `dispersion_atm` (bonded-triple ATM)
    **default false**, opt-in `amideh_acidity_order_bug`.
  - `f4828b06`: `coulomb_r_cut` PARAM exposed (default 100 Bohr unchanged).
  - `2dd151d0`/`526fcb8a`/`1902a19c`: sparse `bpair`/`topo_distances` (`SparseTopoTable`),
    cell-list `nb_hc`/`nb_nometal`, HB/XB energy-magnitude pruning.
- **Tests/tools**: `cli_curcumaopt_07` repaired (runtime reference), `md_adaptive_step` ctest,
  `scripts/refset_regression.py`, `scripts/scan_convergence.py`.

## Conflicts and resolutions

| file | resolution | confidence |
|---|---|---|
| `src/core/units.h` | both added `MD_TIME_UNIT_FS`/`FS_TO_MD_TIME` with the identical value; kept ours, appended the remote's O-H-period verification sentence to the comment | high |
| `src/capabilities/simplemd.h` | kept our `wrapIntoContainer()` AND the remote's `IntegratorStep()/adaptiveStepTolerance()/hottestAtomRatio()`; `m_dt2` removal: ours (with comment); `hydrogen_mass` kept ours (`min=1`) + the remote's 9 `adaptive_step*` PARAMs | high |
| `src/capabilities/simplemd.cpp` | `Verlet()`/`Rattle()`/`NoseHover()`: ours (identical code; taking both would redeclare `dt`). Exactly three `FS_TO_MD_TIME` conversions remain (checked). Final report: our `CurcumaLogger::raw` line + the remote's adaptive-step summary | high |
| `test_cases/check_md_time_axis.py` (add/add) | ours (covers gfn2 + gfnff; it is the one the branch's test history refers to) | high |
| `test_cases/CMakeLists.txt` (auto-merged, but **broken**) | both sides registered `md_time_axis` -> duplicate `add_test` NAME; removed the remote's registration, left a comment | high |
| `src/core/energy_calculators/ff_methods/gfnff.h` | kept both PARAM blocks; dropped the remote's second `coulomb_implicit` (ours already default true, with the Sep 18 note) | high |
| `src/core/energy_calculators/ff_methods/gfnff_method.cpp` | (1) our react-jump diagnostics + the remote's `updateNonbondedRepulsionIfNeeded/updateDispersionPairsIfNeeded/refreshDispersionC6/updateCoulombPairsIfNeeded`: both. (2) `Calculation()`: our reuse-topology check first, then the remote's D4 skin check (topology before pair lists). (3) `generateRepulsionPairsNative()`: the remote's `hh_table.get(i,j)`; **plus two non-conflict compile fixes** in our rev-gfnff code, which still indexed `topo_distances[i][j]` (now a `SparseTopoTable`): -> `.get(i,j)`; and the `hh_table` declaration hoisted above our `nonbonded_set` lambda so the rev blend partner set uses the same H...H classification as the plain list | medium - see "Unsure" |
| `src/core/energy_calculators/qm_methods/xtb_scf.cpp` | the remote's `-scf_reduce` knob, but inside our `CURCUMA_XTB_HAVE_LAPACK_SYEVD` guard: the remote called `dsygst_` unconditionally again, the same BLAS-less build break the Sep 18 merge fixed | high |
| `ff_workspace.h` (auto-merged, but **broken**) | identical `bonds()` accessor added on both sides -> overload error; kept one | high |
| `AIChangelog.md`, `CLAUDE.md`, `docs/REV_GFNFF_TODO.md` | both sides kept. Numbering collisions fixed: remote REV_GFNFF_TODO #11/#12/#13 -> **#16/#17/#18**, remote CLAUDE.md Known Issues #31/#32 -> **#34/#35**; cross-references updated in AIChangelog, CLAUDE.md, TODO.md, docs/GFNFF_STATUS.md, ff_methods/CLAUDE.md | high |

## Build and test (CPU only, `release/`, same flags as the main tree's release: BLAS, AVX2/512, march=native, no TBLITE/XTB/ULYSSES/CUDA/ROCm/Vulkan)

External deps were copied from the main tree's `external/` (read-only copy, no network fetch).
Baseline = `9c69d40e` built in the same worktree before the merge (binary md5 `f45065c9…`),
merged binary md5 `25a953b4…`. Both builds exit 0. `make GenerateParams` reports no warnings.
`ctest -j16 --timeout 1500`:

| | total | failed |
|---|---:|---:|
| baseline `9c69d40e` | 304 | **12** |
| merged | 305 (+`md_adaptive_step`) | **13** |

Failing in both (11): `confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`,
`cli_confscan_01..07` (7), `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`.

- **Fixed by the merge (1)**: `cli_curcumaopt_07_opt_multixyz` (the remote's `4e40fbb2`).
- **New failures (2)**, both traced to the remote's intentional default changes, not to the
  conflict resolution:
  - `cli_gfnff_04_rev_well_form`: plain gfnff and explicit gauss energies are 1.2195e-8 Eh off
    their pinned values. With `-gfnff.dispersion_atm true` the merged binary reproduces both
    pins **exactly** (-4.672737068614 / -4.673521653477). Cause: ATM term now off by default
    (`7d2cceb7`). Fix = re-pin the two values, or decide whether rev-gfnff should keep ATM.
  - `cli_simplemd_20_gfnff_rev_h_budget`: the negative control (budget fix off) no longer
    violates the bounds (29.5 kJ/mol instead of 1234.3 on baseline). With
    `-gfnff.dispersion_atm true -gfnff.hh_repulsion_bpair false -gfnff.dispersion_c6_update false`
    the merged binary reproduces the baseline control **bit-for-bit (1234.29 kJ/mol)**; ATM
    alone (26.2) or C6-refresh alone (27.6) is not enough, hh_bpair alone does nothing here.
    This is the chaotic single-trajectory fragility the test header already documents
    ("a rare event of ONE chaotic trajectory"); the test is an operator decision, not re-pinned.
  Neither test was modified.

## Unsure / not verified

- **rev-gfnff 2^k corner blending vs the per-step pair-list refresh.** The remote's refresh
  functions rebuild the repulsion / D4 / Coulomb lists and the D4 C6 of the workspace's *slot*
  corner only. During an active transition the other corners (`FFWorkspace::TopologyState`)
  keep the lists and C6 they were built with at transition start. Below 800 atoms (every rev
  test system) the repulsion rebuild is skipped, but the C6 refresh is not, so during a blend
  the corners' dispersion may be evaluated with C6 from different CNs. The react/rev ctests
  that pass (13-17, 19, 21+, `gfnff_react_*`) do not show a problem, but none targets this.
  Worth an explicit check (or refreshing C6 per corner) before relying on blended energies.
- `getCachedBondList()`/`getCachedTopology()` are reused by the per-step repulsion rebuild;
  in plain (non-react) auto mode, a geometry-tracker re-perception could make that list
  differ from the bonded terms. Pre-existing in the remote, not merge-specific.
- `hh_repulsion_bpair` now also governs the H...H factor of the rev-gfnff blend partner set
  (my hoist). Consistent with the plain list, but a behaviour change for rev mode.
- **GPU code: not compiled, not run.** No CUDA/ROCm/Vulkan SDK is configured in this build and
  no GPU was used. Everything in `cuda/`, `rocm/`, the GPU plugins, multi-GPU eigensolve,
  `eeq_solver_gpu` (WP7-E), `bench_syevd_mg.cpp` is merged as-is from the remote (no
  conflicts there) and untested here.
- Threads > 1 and the `-adaptive_step` integrator on rev systems were not exercised beyond
  ctest.
