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

# Reconciliation onto `reactff2-llm` (Sep 25, 2026, second agent)

Claude Generated. CPU only, not pushed. Done in the main working tree on `reactff2-llm`, whose tip was
`556490cd` at the start. Method: `git merge --no-ff reactff2-llm-merge-multigpu`, which reuses the
first agent's conflict resolutions from `4c39e4a1`. A fresh merge of `origin/feature/multi-gpu`
(`7d2cceb7`, the same commit) would have needed all of them redone.

## Conflicts on top of the first merge

- `CLAUDE.md` / `AIChangelog.md`: both sides kept. **The numbers collided again.** Today's commits
  already use Known Issues #34 (frag_charge) and #35 (stale CN), so the first agent's #34/#35 for
  the multi-gpu entries are now **#36** (pair-list refresh) and **#37** (pprcht deviations). Their
  cross-references in `CLAUDE.md`, `AIChangelog.md` and `TODO.md` were updated to match.
- `test_cases/cli/curcumaopt/07_opt_multixyz/golden_energies.txt`: our side modified the file, the
  remote side deleted it. Deleted, because multi-gpu's `run_test.sh` computes the reference at run
  time and no longer reads it. The test passes.

## Duplicate C6-refresh fix (operator decision 1): equivalent, so multi-gpu's was removed

- **Code comparison.** Both set `d.C6 = getChargeWeightedC6(...)` from Gaussian weights computed at
  the current CN, on the CPU only, and both are skipped when `gpu_only`/`reuse_cn` is set. The
  differences:
  - (a) Ours also refreshes every stored rev-gfnff corner D4 list (`forEachD4PairList`); multi-gpu's
    refreshes only the slot list. So ours covers more.
  - (b) Ours reads the weights through `refreshC6WeightsForCN`, which has no CN-change cache.
    Multi-gpu's goes through `updateCNValuesForGradient` and therefore also forces
    `d4_cn_cache_threshold = 0`.
  - (c) Multi-gpu's skips the setup geometry.
  - (d) Multi-gpu's has an on/off PARAM, `dispersion_c6_update`, which also switches the GPU
    device-side refresh.
- **Numerical check.**
  - The merged binary with only our fix gives trajectories and optimised geometries **byte-identical**
    to the merged binary with both fixes: 200 fs gfnff MD at 800 K, and gfnff `-opt`, each on triose
    and caffeine (`in.trj.xyz` / `in.opt.xyz` compared with `cmp`).
  - At the setup geometry, recomputing the C6 changes nothing: the dispersion term is identical to
    10 digits with `dispersion_c6_update` true and false (triose -0.0181251182, caffeine
    -0.0562085554).
- **What was kept.**
  - The PARAM `dispersion_c6_update` now gates our block on the CPU; the GPU refresh stays as it was.
  - The forced `d4_cn_cache_threshold = 0` is kept (with a new comment): it keeps the gradient's
    dC6/dCN at the same CN as the energy's C6.
- **What was removed:** `refreshDispersionC6`, `dispersionC6Stale`, `m_disp_c6_geometry` and
  `FFWorkspace::d4DispersionsForC6Refresh`.

## Tests

- **`cli_gfnff_04` (decision 2):** re-pinned for ATM-off: gauss -4.673521653477 -> -4.673521641282,
  gfnff -4.672737068614 -> -4.672737056419. With `-gfnff.dispersion_atm true` the old pins come back
  exactly. The comment asking for a re-evaluation at stage 3 (WP5) is next to `PARAM(dispersion_atm`
  in `gfnff.h`.
- **`cli_gfnff_05` (not in the brief, same cause):** this test did not exist on the first merge's base.
  Its H3OpH2O2 identity value moved 0.466545744259 -> 0.466545743025 (-1.2e-9 Eh), and
  `-dispersion_atm true` alone restores it, so I re-pinned it under decision 2.
- **`cli_simplemd_20` (decision 3), NOT green.**
  - The control arm is now informational and no longer gates the test. It no longer explodes; it
    already did not explode on the pre-merge binary (22.79 kJ/mol, Fix B), and on the merged binary
    it gives 25.03.
  - **But now the SHIPPED-default arm violates**: 224.44 kJ/mol against the 150 bound. Total energy
    is not conserved at one step, t = 1.3797 ps (+165 kJ/mol), right after an H-H topology event.
  - Restoring three multi-gpu defaults brings back the pre-merge trajectory bit-for-bit
    (43.35 kJ/mol): `dispersion_atm true`, `hh_repulsion_bpair false` and
    `d4_cn_cache_threshold 0.01`.
  - Replicates, 12 runs each (T = 1975..2030 K in 5 K steps; the seed does not change the start
    velocities). Violations: pre-merge **1/12**, merged **2/12**, merged with cache 0.01 **3/12**.
    That is the known rare event of one chaotic trajectory, and n=12 cannot tell these apart. The
    test's T = 2000 K simply draws a violation now.
  - The bounds were left unchanged and the test is left failing on purpose. **Operator question.**

## Verification (CPU, same configuration as `build_rev`, which is the package-32 baseline, md5 `999b8f90`)

- **Build:** `make` exit 0; merged binary md5 `3c08e62b`.
- **`ctest`** (with `CURCUMA=<merged binary>`): **295/307 passed, 12 failed** vs baseline
  **294/306, 12 failed**. The set of failing tests is identical:
  - `confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`
  - `cli_confscan_01..07`
  - `cli_simplemd_18`, `cli_simplemd_20`

  The one extra test is the new passing `md_adaptive_step`.
  - Pitfall: without `CURCUMA`, the CLI tests pick up `<repo>/release/curcuma` (an old Sep 18 binary)
    and `cli_gfnff_03/04/05/06` fail spuriously.
  - Pitfall: the test scripts are copied at configure time, so rerun `cmake .` after editing a
    `run_test.sh`.
- **Old binary vs new binary, gfnff** (`refset_regression.py`, and a two-binary script for
  S30L-CI). With the multi-gpu defaults:

  | set | structures | changed | MAD (kcal/mol) | max (kcal/mol) |
  |---|---:|---:|---:|---:|
  | GMTKN55 | 2462 | 1782 | 1.2e-4 | 0.066 (`AL2X6/al2me5`) |
  | MOR41 | 285 | 267 | 5.2e-4 | 2.7e-3 |
  | S30L-CI | 90 | 89 | 2.0e-4 | 1.3e-3 |

  Adding `-gfnff.dispersion_atm true -gfnff.hh_repulsion_bpair false -gfnff.hb_min_pair_energy_eh 0
  -gfnff.xb_min_pair_energy_eh 0` makes **all three sets bit-identical** to the old binary. So every
  energy change comes from those four multi-gpu defaults: ATM off, bpair H...H repulsion, and the
  HB/XB energy pruning. The first status report listed all of these as number-moving.
- **GMTKN55 gfnff vs xtb** (`gmtkn55_compare.py`, copied to scratch, curcuma keys recomputed, xtb
  cache reused): MAD 0.860 -> 0.859, max 131.608 unchanged, n=2460.
- **gfn1/gfn2, old binary vs new binary:** MOR41 gfn1 285/285 and gfn2 285/285 bit-identical,
  GMTKN55 gfn2 2462/2462 bit-identical (MAD 0.000).
- **Not verified:** any GPU code (CUDA/ROCm/Vulkan, multi-GPU eigensolve, WP7-E EEQ); no SDK or
  GPU here. That includes the GPU C6 refresh, which `dispersion_c6_update` still switches.

## Commits (not pushed)

1. `df839027`: merge commit, conflict resolution and Known Issue renumbering.
2. `228f5d55`: removes the duplicate C6 refresh; the source code change only.
3. The commit that adds this section: the test re-pins (04, 05), the re-scoped control arm of
   test 20, the ATM re-evaluation comment, and the doc notes in `CLAUDE.md` #35/#36 and
   `AIChangelog.md`.

Left uncommitted on purpose: the orchestrator's concurrent package-33 notes in `AGENT_STATE.md`,
`WORK_STATUS.md` and `docs/REV_GFNFF_STAGE2.md`. They are not part of this reconciliation.
