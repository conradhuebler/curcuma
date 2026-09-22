# rev-gfnff: agent state and handover index

**Read this file first after a context clear or in a fresh session.** It is the one place where
delegated work, git state and pending decisions are recorded so nothing is lost. In-process
subagents do NOT survive a clear; their output files do. Rule: every delegated agent gets a fixed
status-file path at launch, writes progress there (numbers, not dumps), and is listed below.
**This file holds ORIENTATION, not archive** — one line per finished package, pointing at
`WORK_STATUS.md` (the full chronological record, package by package) or a package's own
`*_STATUS.md` file for every number. Keep this file under ~250 lines; if it grows past that,
compress old package rows further before adding new ones, following the pattern below (this file
itself grew to 713 lines by 2026-09-22 before a full compression — the pre-package-1 history that
used to fill 300+ lines here now lives only in `docs/REV_GFNFF_ROADMAP.md` WP0-WP5 and git history).

**Agent tiers (operator rule, 2026-09-12).** Lower tier launched without asking (Sonnet = routine:
campaigns, measurements, tests, scripts to spec; Opus = bounded coding). Higher tier (**Fable**) is
an operator decision: propose the task, then ask before launching. rev-gfnff needs Fable for
method/chemistry judgement calls. See memory `agent-division-of-labour`.

## Current state (2026-09-22, no agent running)

Branch `reactff2-llm`, **106 commits ahead of `origin/reactff2-llm`, nothing pushed**, working
tree clean. rev-gfnff defaults, all operator-decided and orchestrator-verified independently of
the delegated agent that built them (see the package index below for each):

| setting | value | since |
|---|---|---|
| `-gfnff.rev_share_form` | `conserving` (+ `rev_share_donor_rule true`) | package 6 |
| `-gfnff.rev_budget_fix_h` | `true` (a no-op under `conserving`, H's cap is 0 by element there) | package 6 |
| `-gfnff.rev_well_form` | `mg3` (bond-order-resolved, continuous key, no discrete switch) | package 12 |
| `-gfnff.coulomb_implicit` | `true` (from the `feature/multi-gpu` merge, package 2) | package 2 |
| `-gfnff.rev_dt_cap` | `0.25` fs — **outside the measured-safe 0.0625 fs band, operator decision open** | pre-existing |
| plain `-method gfnff` | bit-identical throughout every package above | — |

`ctest -R "gfnff|sqm_val|react|cli_simplemd_|cli_gfnff_"`: **111/113**. The two failures are
diagnosed, neither is a code defect, both await an operator call on the TEST rather than the code:

- **`cli_simplemd_18`**: its `Etot(t)` OLS-slope statistic is not well-posed on this 12-H2/8000K
  NVE bath — the trajectory is a saturating step (one H2 dissociates in ~1 ps, then flat for 9 ps),
  so the slope measures a one-time plateau height, not a rate; `mg3`'s deeper H-H well makes that
  plateau bigger, which reads as "more dissipative" but the actual injection rate is LOWER for
  `mg3` than `mg` (package 13, `MG3_DISSIPATION_STATUS.md`). Needs a better statistic or a
  differently-chosen bath, not a threshold tweak — the same 2.5e-3 floor also fails plain `gauss`
  at small dt and is non-monotone in dt for every arm.
- **`cli_simplemd_20`**: negative control no longer violates since the MD clock fix (package 10) —
  awaits recalibration to the true-fs regime.

**Also open, not urgent:** `CurcumaUnit::Constants::ATOMIC_TIME_TO_FS` / `Time::ATOMIC_TIME_TO_FS`
in `units.h` are really aut->**atto**seconds (off by 1000x), zero use sites, deliberately not
touched (package 10). Stage 3b is only partly exploited (`mg3` splits by bond order, not yet by
element environment; `mg2`'s free-curvature fit is not combined with the bond-order table). Stage
2 (split-charge model) has not been started. Further `feature/multi-gpu` material (adaptive
step-rejecting integrator, large-system MD work, GPU eigensolver stack) deliberately not imported
— it was reviewed (see package 10's background) but never requested beyond the MD clock fix.
Merge to `master` vs. keep developing on `reactff2-llm`: not yet decided.

## Package index (full detail in `WORK_STATUS.md`, one section per package)

| # | what it did | headline result | status file |
|---:|---|---|---|
| 1 | `rev_budget_fix_h` default on + new 20-cell baseline | H keeps one valence; removed the original runaway | `WORK_STATUS.md` §1 |
| 2 | merged `origin/feature/multi-gpu` | `coulomb_implicit` adopted; identity proven at `-threads 1` | `WORK_STATUS.md` §2 |
| 3 | built the `conserving` valence share | adducts −87..−107 → +2..+32 kcal/mol; dative/ylide neutrals regressed 73-110 kcal (donor rule needed) | `WORK_STATUS.md` §3 |
| 4 | built `mg`/`erfmorse` wells, curvature-pinned | class-A rms 24.5 → 19.5/19.2; `mg` cheaper, chosen | `WORK_STATUS.md` §4 |
| 5 | calibration, 3 new falsifier tests, docs | `cli_simplemd_16/18/19/20`, `cli_gfnff_03/04` | `WORK_STATUS.md` §5 |
| 6 | donor rule + BOTH flips shipped | `conserving`+donor and `mg` become default; dative/ylide regression → 0.00 | `WORK_STATUS.md` §6 |
| 7 | corrected package 6's "mg×conserving interacts" | 130-cell sweep: it's `conserving` alone, on an untested cell; `mg` only shifts frequency; resolution-failure mechanism (fixed-width transition window under-resolved at high T) | `WORK_STATUS.md` §7 |
| 8 | shipped a runtime warning (no default change) | react-mode MD + `conserving` above 0.125 fs (later corrected to 0.0625, package 10) warns once at start | `WORK_STATUS.md` §8 |
| 9a/9b | built `mg2` (free curvature+r0) and `mg3` (+ bond-order), both opt-in | class-A rms → 15.8/13.2; guard cost 1.044→1.054/1.055; bond-order key is continuous, verified no hidden switch | `WORK_STATUS.md` §9 |
| 10 | fixed a GENERAL curcuma bug: MD clock ran 1.9516144204x too fast | every pre-Sep-22 "fs"/"ps" number describes the wrong timescale (relative comparisons unaffected); re-derived the react warning/`rev_dt_cap` band in true fs | `WORK_STATUS.md` §10 |
| 11 | clean paired-replicate re-measurement of the mg/mg2/mg3 tail | 11 700 trajectories: NOT distinguishable at any dt — both the "mg3 worse" and "mg3 is best" prior claims were single-trajectory artefacts | `WORK_STATUS.md` §11, `tail_remeasure.csv` |
| 12 | made `mg3` the `rev_well_form` default | class-A rms 22.2→13.2, bond length error 0.025→0.004 Å, guard +0.011 kcal/mol; `gfnff` untouched | `WORK_STATUS.md` §12 |
| 13 | diagnosed `cli_simplemd_18`'s "mg3 more dissipative" reading | a saturating-trajectory statistic artefact, not a defect (see Current state above) | `WORK_STATUS.md` §13, `MG3_DISSIPATION_STATUS.md` |

**Standalone diagnostic agents from the 2026-09-14/15 runaway investigation** (before "package"
numbering started; each is its own `*_STATUS.md` under this directory, referenced from
`WORK_STATUS.md`'s early sections and `docs/REV_GFNFF_ROADMAP.md`): `r0fix`/`wall-etot`/
`nve-test-2`/`poly-jump`/`orca-ref`/`outliers`/`eeq-drift`/`refpaths` (data-basis and stage-1/3a(i)
work), `valfix`/`proxy`/`fable-bondstate`/`qp-bondstate` (the exhausted-weight-space investigation
into a smooth valence share, closed — see [[Valenzanteil im reaktiven GFN-FF]] in the vault),
`h-budget`/`verbosity-traj`/`runaway`/`scan-cadence`/`break-tail` (the hydrogen-budget runaway root
cause), `harness-kept`/`baseline-head`/`guards`/`classa-harness`/`pairs-blend`/`uks-inspect`/`cij`
(the class-A/guard measurement apparatus this whole project stands on).

## Recurring traps, worth reading before writing a new harness

- **A `cmd.txt` is not provenance.** Fingerprint an MD arm by its trajectory (rebuild count, max
  jump), not by what its launch command claims. Memory `revgfnff-run-provenance-fingerprint`.
- **`curcuma -version` embeds `git describe`, so raw binary md5 is not a same-source identity test
  across commits** (a doc-only commit changes it). Use a behavioural fingerprint; record md5 only
  to say which build produced a number.
- **A react-MD tail/rate number from ONE trajectory is one sample of a chaotic or saturating
  process.** Use paired perturbation replicates (`scripts/revgfnff_tail_sweep.py`) and check the
  trajectory SHAPE (early- vs late-window slope) before trusting a single derived statistic — the
  orchestrator itself was caught by this twice (packages 7 and 13). Memory
  `revgfnff-tail-needs-replicates`.
- **Shell quoting**: never bundle a flag and its value in one shell variable (zsh does not
  word-split); write separate literal argv tokens, one variable per token.
- **`*.topo.json` caches silently reuse a stale topology across geometries with the same element
  list — delete it before every re-measurement, fresh directory per structure.**
- **Count rebuilds as `max(REACT rebuild #N)`, not `grep -c`** — from calculator verbosity 2 on,
  the line is printed twice (GFN-FF + SimpleMD flush).
- **`SimpleMD` runs the calculator one verbosity level below the run** (Known Issue #31) —
  `CURCUMA_BLENDDUMP`/`CURCUMA_SHAREDUMP` etc. need `-verbosity` one higher than you'd expect.
- **`test_cases/cli/test_utils.sh` picks `release/curcuma` before `build_rev`** — export
  `CURCUMA=$PWD/build_rev/curcuma` before any `ctest -R cli_...`; and `ctest` runs the COPIES of
  scripts in the build tree, so run `cmake .` there after editing a `run_test.sh`.
- **A fresh single point cannot reproduce a mid-transition corner-blend state** (`w_a` latches to
  the coordinate at transition start) — an FD check "at a mid-transition frame" via a fresh SP is
  not well-posed; use dt-scaling instead to bound a chain-rule term.

## Still-open items not superseded by any package above

- **H-H repulsion blend**: `rev_bo5_center` 1.3 → 2.5 measured to fix the residual +16.28 → +2.81
  kcal/mol (D_e unchanged) but never built in as a default — small, self-contained, still pending.
- **WP0 leftovers**: the default hysteresis 1.6/2.6 fires from ordinary vibrations at 2500 K;
  react-mode exchange resolutions are folded into "broken" in the event record, losing
  comparability with the old 115/248 count.
- **Hypothesis, untested**: the class-H UKS "state instability" in the r2SCAN-3c reference data may
  be ORCA's `$new_job` chain artefact (rkt02's chain was 11.4 kcal/mol off), not SCF
  nondeterminism.
- **Methodological note from `nve-test-2` (pre-dates the 12-H2 bath)**: 2 H2 in NVE produces zero
  react events at any wall radius/temperature tried — a colder, non-reactive system is needed for
  a sensitive NVE-drift discriminator; and an active confinement wall removes energy through the
  rotation projection unless `-md.rm_COM 0 -md.rmrottrans 0` is set explicitly.
- **Stage 2 plan** (next major stage, not started): localised integer `q0` per fragment instead of
  the uniform rule, `kappa_Z` fitted on closed-shell charged NCI (AHB21/CHB6/IL16, n=43) + anionic
  SN2, guarded by S66 and the 285-reaction conformer set; Cl2-/F2- stay report-only (open-shell).
  Design doc: `docs/REV_GFNFF_STAGE2.md`.
