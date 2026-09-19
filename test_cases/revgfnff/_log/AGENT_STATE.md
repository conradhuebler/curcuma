# rev-gfnff: agent state and handover index

**Read this file first after a context clear or in a fresh session.** It is the one place where
delegated work, git state and pending decisions are recorded so nothing is lost. In-process
subagents do NOT survive a clear; their output files do. Rule: every delegated agent gets a fixed
status-file path at launch, writes progress there (numbers, not dumps), and is listed below. Keep
this file under ~160 lines; move finished detail into docs/.

**Agent tiers (operator rule, 2026-09-12).** Lower tier launched without asking (Sonnet = routine:
campaigns, measurements, tests, scripts to spec; Opus = bounded coding). Higher tier (**Fable**) is
an operator decision: propose the task, then ask before launching. rev-gfnff needs Fable for
method/chemistry judgement calls. See memory `agent-division-of-labour`.

## Git state (branch `reactff2-llm`, verified 2026-09-14 night)

Commits since 2026-09-13 (newest last):

| commit | summary |
|---|---|
| `e36d9925` | **stage 3a (i)**: drop the stretched pair's own `cn_ij` from its own r0 |
| `aef7aaf0` | Restate the equilibrium guard as a relative criterion (roadmap decision 8) |
| `c66ffa13` | Attribute the post-r0-fix bond outliers to their owning stage (WP5) |
| `36dfe72f` | Restate two acceptance criteria that measured the deliberate stage-1 deviations |
| `ceb2d160` | **stage 3a (ii)**: valence-share factor `c_ij` (rkt06 rms 15.79 -> 2.71) |
| `2ea9f8cf` | Four missing class-A pairs + the class-S topology-frame fix |
| `bb6079b8` | classa harness `--mode kept` opt-out (the topo-reuse fix changes its meaning) |
| `283ed23d` | Track `test_cases/revgfnff/_log` (was untracked, 24 files) |
| `11b1baea` | Smooth 1,3 proxy, opt-in default off (`-gfnff.rev_share_onethree`) |
| `0bd3c80a` | **Cherry-picked** the batch topology-reuse fix (`76e7f83a` from `fix/topo-cache-reuse`) |
| `95d879a2` | Fable consultation recorded (its .md files were only added in `a474be89`) + full-ctest verification |
| `a474be89` | **QP gate not cleared** + `rev_bo13_ordinary_join` (default off) + `shareD` dump; **3331 kJ event re-attributed to the proxy**, FABLE_BOND_STATE/PROXY_STATUS corrected |
| `1b7e5ff0` | **MD time-step unit fix (general curcuma)**: `-md.time_step 1.0` integrated 1.9516 fs; new `MD_TIME_UNIT_FS`, ctest `md_time_axis` |
| `cd8c64d9` | **rev-gfnff re-derived in true fs**: warning threshold 0.125 -> 0.0625, `rev_dt_cap` help, mg/mg2/mg3 ordering falsified, `cli_simplemd_20` flagged |

**FULL `ctest` at HEAD (2026-09-14 night, release/ rebuilt after the cherry-pick): 14 failures out
of the whole suite, none caused by the cherry-pick — orchestrator-verified each one.**

- **4 pre-existing, unrelated**: `confscan_dtemplate` (flaky), `test_orca_interface`, `xtb_cpscf`,
  `cli_curcumaopt_07_opt_multixyz`.
- **7 newly-surfaced but pre-existing**: `cli_confscan_01/02/03/04/05/06/07` all fail on a stale
  assertion string (`"All 44 input structures were read"`) that has **never existed** in
  `confscan.cpp`'s history — the actual output is `"Processed: N / N"`. A test-string bug, not a
  code regression; nobody had run the full confscan suite in a while.
- **3 newly-surfaced, real but pre-existing**: `cli_simplemd_16/18/19` (all `gfnff react rev`)
  fail on their MD-behaviour thresholds (test 16: T_max 7946 K vs a 6400 K cap; test 18: rebuild
  counts 8-13x the calibrated 46-52 and slopes ~2-3 orders above the 1.6e-3 Eh/ps floor; test 19:
  the join mechanism itself still works exactly as designed — 0 formations under `order`, the
  `weight` negative control still fires — but its energy-agreement tolerance is now exceeded by
  -0.3995 vs 0.01 kcal/mol). **Orchestrator-verified all three reproduce identically on
  `build_rev/curcuma`, which does NOT contain the cherry-pick** — so the cause is the accumulated
  stage 3a (i) r0-fix + stage 3a (ii) valence-share physics change, and these three tests were
  simply never re-calibrated against it. This is the first full-suite run since both changes
  landed together in one binary. **Needs a recalibration pass** (like the original `nve-test-2`
  task) before these three can be trusted again; not attempted here.

**Not yet decided**: merge to `master`, or keep developing on `reactff2-llm`.

## Live agents

**2026-09-20: `package-10` DONE — a general curcuma bug fixed (the MD clock), and every rev-gfnff
"fs" re-derived.** No agent running. Two local commits, nothing pushed. Detail `WORK_STATUS.md`
package 10 (`Packages done: 10/10`).

- **The core fix (general curcuma, not rev-gfnff)**: `SimpleMD` integrates in Angstrom/amu/Hartree,
  whose implied time unit is `sqrt(amu*A^2/Eh) = 1.9516144204 fs` (derived three ways from
  curcuma's own CODATA-2018 constants; new `CurcumaUnit::Constants::MD_TIME_UNIT_FS`). The user's
  step was handed over unconverted, so **every curcuma MD ran 1.9516144x faster than requested**
  and `-MaxTime` was stretched by the same factor — every method, every caller (ConfSearch's
  exploration MD and polymerbuild included). Re-implemented here rather than cherry-picked, because
  `ef462fcf`'s diff also touches an adaptive integrator this branch does not have. Converted only
  in `Verlet()` / `Rattle()` / `NoseHover()`; everything that counts or schedules time stays in fs.
- **Evidence**: water local-mode period against the Hessian (a path validated to 0.13 % vs xtb
  6.7.1) — gfnff / gfn2 / gfn1 ratios **1.9536 / 1.9520 / 1.9511** before, **-0.09 / -0.01 /
  +0.02 %** after. The scheme is untouched, proven exactly: pre-fix `-dt 0.25` and post-fix
  `-dt 0.4879036051` give **bit-identical** trajectories (202 frames, max |dx| = 0.000e+00).
  Static paths byte-identical (9/9 `-dump_gradient` files, 12-digit energies, `diff -r` clean).
  NVE still scales dt^2 (CH4 rms ratios 4.01 gfnff / 4.17 gfn2). New ctest **`md_time_axis`**,
  verified to fail on the pre-fix binary.
- **`ctest` blast radius: exactly ONE new failure.** Full suite 12/304 both before and after —
  before = 11 pre-existing + `md_time_axis` (failing by design), after = the same 11 +
  **`cli_simplemd_20_gfnff_rev_h_budget`**. `cli_simplemd_16/18/19` pass unchanged (their slope /
  mean / count criteria turned out to be time-scale robust).
- **Test 20 is FLAGGED, not recalibrated — operator decision.** Its negative control
  (`-gfnff.rev_budget_fix_h false`) went 2593.60 kJ / 0.559 a0 -> 32.68 / 2.120 and no longer
  violates. Not the fix's fault: restoring the OLD physical regime on the FIXED binary (every
  dt-derived setting x 1.9516144204) reproduces all three arms **exactly** (59.26 / 2593.60 /
  70.77 kJ, 1.770 / 0.559 / 1.536 a0). And there is no dt to recalibrate to — with the discrete
  dynamics held fixed, the violation appears only at 0.4879 and 0.60 fs while 0.45, 0.50 and 0.55
  stay inside, so it is one chaotic trajectory, not a threshold in dt (`-md.seed` is inert on this
  path, 8 seeds bit-identical). The falsifier for `rev_budget_fix_h` now exists only ABOVE the
  shipped 0.25 fs cap; it needs a new cell/temperature or a re-scoped assertion. A dated header
  block in `run_test.sh` records all of this; no threshold or flag was touched.
- **True-fs re-characterization** (130 cells, 9.758 ps each = package 7's actual physical exposure;
  the protocol was anchored first by reproducing package 7's shipped-default row exactly — 9708
  rebuilds, median 49.4, p90 214.2, max 800.1, 27 cells >100 kJ, T_max 28392): per-step |dEpot|
  median / cells above 100 kJ/mol = 97.8/64 at true 1.0 fs, 51.8/25 at 0.5, **27.5/11 at 0.25**,
  13.1/10 at 0.125, **6.3/0 at 0.0625**. Warning threshold `rev_react_dt_advice` **0.125 ->
  0.0625** — the first true step at which NO cell exceeds 100 kJ/mol, which also happens to equal
  the old value rescaled (0.125/1.9516 = 0.0640). Warning text rewritten with the table and an
  explicit "these are real femtoseconds" note; gate re-verified (fires at 0.25 and 0.125, silent at
  0.0625 and 0.05, silent for `delivered` and for plain `gfnff`).
- **`rev_dt_cap` 0.25 now means a genuine 0.25 fs** (it used to integrate 0.488). That alone halved
  every robust tail statistic for free — median 51.8 -> 27.5, cells >100 kJ 25 -> 11, T_max 18928
  -> 11948 K — but 0.25 fs is still **outside** the band where the whole sample is bounded, so the
  warning still fires for a plain react run. **Default NOT changed; operator decides.**
- **mg/mg2/mg3 tail claim: DISPUTED — orchestrator spot-check contradicts the "mg3 is vindicated"
  punchline on the single most-scrutinised cell.** The package's own 130-cell AGGREGATE numbers
  (rebuild `dE_jump` events >= 50 kJ per 1000 rebuilds: old clock mg 1.65 < mg2 2.87 < mg3 5.78;
  true 0.25 fs mg 0.62 < mg3 1.36 < mg2 1.78) are not disputed as raw numbers, but the summary drawn
  from them ("mg3 is the BEST of the three... should not weigh against mg3 in the adoption
  decision") does NOT hold on direct re-check. **Orchestrator, true dt = 0.25 fs on the current
  `build_rev` binary (no conversion needed, the fix makes `-md.time_step 0.25` genuinely 0.25 fs),
  on `c2h6/T2000_f0`** (the single cell every prior package used as the worst-case reference,
  package 7's original finding): **mg 28.3 kJ / 0 events, mg2 32.2 kJ / 0 events, mg3 256.8 kJ / 12
  events.** `mg3` is 8-9x worse than `mg`/`mg2` on exactly the cell that matters most, not
  comparable or better. **Do not adopt `mg3` on the strength of package 10's aggregate summary
  alone — the aggregate and the single-cell picture disagree, and a proper full re-sweep (not just
  one summary statistic) is needed before this decision is safe to make.** The conformer/S66 guard
  (1.0543/1.0547) and the 0.0212 A equilibrium shift are unaffected by any of this (they are
  time-step independent, static-energy quantities) and stand as previously measured.
- **Not done**: no default value changed; nothing else imported from `origin/feature/multi-gpu`
  (no adaptive integrator, no GPU eigensolver, no `MD_LARGE_SYSTEMS.md`); `test_cases/revgfnff/_log/*.md`
  other than `WORK_STATUS`/`AGENT_STATE` still quote the OLD time scale.

**2026-09-20: `package-9` DONE — 3a(iii) step 2 and 3b BOTH delivered, both OPT-IN, no default
flipped.** No agent running. Four local commits, nothing pushed. Binary `89a587bb` — note that
`build_rev/curcuma` embeds `git describe`, so the md5 changes on every commit with no source
change (`89a587bb` -> `7f2e2504`, verified identical: 5/5 well forms to 12 digits, 9/9 react-MD
cells). Use the trajectory fingerprint, not the md5, to compare builds across commits here.

- `-gfnff.rev_well_form mg2` (free curvature + r0 re-solve) and `mg3` (bond-order-resolved),
  table `src/core/energy_calculators/ff_methods/rev_well_table_v2.h`.
- **Guardrail verified at the level that can see it**: `gauss`/`mg`/`erfmorse` are bit-identical
  to the pre-session binary `4c323b80` — 12 digits on the four reference single points AND
  8 of 8 react-MD cells reproduce rebuild count / per-step max / n>=50 / T_max exactly.
- Class-A harness rms 24.68 (gauss) / 22.15 (mg) -> **15.83 (mg2) / 13.22 (mg3)**; median dev D_e
  -25.79 / -14.32 -> -16.17 / **-6.86**; median dev k -169 / -158 -> **-44** / -82; median
  |b_model - b_r2SCAN-3c| over 32 class-A bonds 0.0293 / 0.0246 -> **0.0060 / 0.0044** A.
  Class D 4.980 / 5.152 -> **4.590 / 4.555** (grad_RMS 14.6 / 14.4 -> **11.4 / 11.3**).
  rkt06 2.384 / 2.379 -> **2.267**. Adducts under `conserving` unchanged. FD gradient 1.14e-07.
  `ctest` 113/113.
- **Two costs, flagged**: the conformer/S66 guard opens 1.0341 -> **1.0543 / 1.0547**, i.e. **2x**
  the package-4 precedent (1.0341 -> 1.0439); and `mg3`'s rebuild `dE_jump` tail over 130 cells is
  **456.3 kJ / 89 events >= 50** against `mg`'s 272.3 / 14. The tail was checked against the
  design risk and is NOT the order dimension (0 of the 10 largest jumps over 259 rebuilds carries
  any order change; the forced-order response is linear to 1 part in 8000).
- Measured and REJECTED alternatives, both recorded: a fitted curvature (runs to a flat quartic
  bottom on 8 of 32 bonds) and an r0 solved on the rigid scan (optimised water O-H 0.9179 A
  against a reference 0.9644). A `dr0 = 0` variant gives a gentler guard (1.0508/1.0490) and a
  worse equilibrium (0.020 A from the reference) — also recorded.
- **Bug found and fixed in `scripts/revgfnff_wellfit.py`**: its scan did not force
  `-gfnff.rev_well_form gauss`, so after the Sep 19 default flip the Gaussian recovery of
  (r0, alpha, k_b) was being run against an MG curve. With the fix the script regenerates the
  COMMITTED `mg` table bit-for-bit.
- **Operator decision open**: whether to flip `rev_well_form` to `mg2` or `mg3`. Detail:
  `WORK_STATUS.md` packages 9a/9b, `docs/REV_GFNFF_STAGE3A.md` section 2.3.
- **Correction (orchestrator, own rebuild md5 `7f2e2504`)**: the agent's hand-back message (not
  written into `WORK_STATUS.md`, so no file needed fixing) claimed "package 7's 130-cell tail does
  not reproduce on the current HEAD binary" (`gauss+conserving` on `c2h6/T2000_f7` giving 96.6 kJ /
  4983 K instead of the recorded 1071.3 / 26265). **That claim is false** — re-run independently on
  the current binary, same protocol: **1071.3 kJ / 172 rebuilds / 309 events / T_max 26265**,
  matching package 7 exactly. Likely a probe slip on the agent's side (wrong flags/binary), not a
  regression — the actual guardrail claim (shipped default + `gfnff`/`revgfnff` single points
  untouched) WAS independently reconfirmed separately: `c2h6/T2000_f0` shipped default still gives
  800.1 kJ / 76 rebuilds / 62 events / T_max 28392, caffeine energies unchanged.

**2026-09-20: `package-8` DONE — package 7's option (ii) shipped as a WARNING, no default changed.**
No agent running. One local commit `84c540d5`, nothing pushed. Detail `WORK_STATUS.md` 8.1-8.5
(`Packages done: 8/8`), `docs/REV_GFNFF_STAGE3A.md` new section 2.2 "Recommended MD settings".
**Orchestrator-independently verified** (own rebuild, md5 `4c323b80`): the `rev_dt_cap` PARAM
(default 0.25, `simplemd.h:779`) pre-dates this work and does clamp `-method revgfnff`'s time
step, confirming the "effective default is 0.25 fs" claim; a plain no-flag react run prints the
warning exactly once with the reported text; `ctest -R "gfnff|sqm_val|react|cli_simplemd_|
cli_gfnff_"` reconfirmed **113/113**.

- **Operator decision (2026-09-20)**: adopt `-md.time_step <= 0.125 fs` for react-mode MD with the
  conserving share as an **operating recommendation** — warn at run start and document it — and do
  NOT change the default time step. The structural fix (transition window in distance instead of
  bond order, package 7 option 1) stays **deferred**; the window code was not touched.
- **The effective default time step is 0.25 fs, not 1.0**: `-md.time_step` defaults to 1.0 but
  `-method revgfnff` clamps it to `-md.rev_dt_cap` (default 0.25, stage 1). So the warning fires
  for a plain no-flag react run — 0.25 fs is exactly where package 7's tail lives.
- **Warning gate** (`SimpleMD::Initialise`, `CurcumaLogger::warn`, once per run, verbosity >= 1):
  `revgfnff|gfnff-rev` AND `topology_mode == react` AND `rev_share_form == conserving` AND
  `dt > 0.125`. The **well form is deliberately not in the gate** (package 7.5: `mg` changes the
  frequency, not the mechanism). Verified on 10 cases: fires for conserving at dt 0.25 (explicit,
  clamped-default, `gauss` well, flat `-rev_share_form conserving`), silent for dt 0.1, for
  `delivered` with either well form, for `topology_mode auto` and for plain `gfnff`.
- **Inert**: `gfnff` caffeine -4.6727370686 unchanged, `ctest -R "gfnff|sqm_val|react|cli_simplemd_|cli_gfnff_"`
  **113/113** (tests 16/18/20 run the warning's own trigger combination and pass unchanged).
  Binary `4c323b80`.

**Earlier, 2026-09-19 night: `package-7` DONE — CORRECTS package 6's framing.** No agent running. One
local commit `7d028e8a` (diagnostic only, no default changed), parent `df928122`. Full detail
`WORK_STATUS.md` 7.1-7.10, `docs/REV_GFNFF_STAGE3A.md` 2.1 rewritten.

**Package 6's "mg+conserving interact" was a 20-cell sampling artefact — the tail belongs to
`conserving` alone, `mg` only changes how often it is visited.** Over 130 cells (same 3 systems,
every frame, both T) `gauss+conserving` alone reaches **1071.3 kJ/mol, T_max 46601 K** on
`c2h6/T2000_f7` — a cell outside the 20-cell grid, worse than the shipped default's 800.1/28392.
**Orchestrator-independently reproduced** (own rebuild, md5 `fa09bcac`): `c2h6/T2000_f7` under
`gauss+conserving` gives exactly **1071.3 kJ / 309 events / T_max 26265** (own run); the shipped
arm on the SAME cell gives only 299.6/132/9571, i.e. `mg` does not create this cell's extreme,
`conserving` alone already has it. Mechanism: a transient geminal H2 (both H still on one carbon)
gets `c = f_H^2 = 0.2531` of its well under `conserving` where the delivered left-over rule gives
exactly 0 — a forbidden configuration becomes cheap, visited 1.5-2.8x more often. **It is a
resolution failure, not a discontinuity**: the corner weight `s` moves 0.000 -> 0.516 in ONE 0.25
fs step (the transition window is ~0.175 a0 wide in distance, ~2 MD steps at 2000 K); halving `dt`
collapses it. **Orchestrator-reproduced on `c2h6/T2000_f0`**: dt 0.25 -> 0.125 fs gives
**800.1 -> 38.7 kJ**, T_max **28392 -> 5579** — matches the agent's report to the digit. Both
share-mechanism hypotheses from the package-7 brief were tested and falsified (the min-width
smoothing parameter has no effect; `gauss+conserving`, which keeps `w` on its well, is the WORST
arm, so a missing `w` on the new well forms is not the cause). No minimal defect found — **this is
a design tension between the conserving share and the fixed-width bo3 transition window**, not a
bug. Options (none built): (i) widen/redefine the transition window in distance not bond order —
fixes it, re-calibrates ctests 13-20; (ii) **`-md.time_step <= 0.125 fs` with `conserving`** —
measured, zero falsifier cost, 2x wall time, cheapest honest mitigation; (iii) multiply `c` by the
settled weight `sig_p` — likely breaks rkt06 (a half-formed TS bond also has sig=0 by
construction); (iv) fix the H-H perception itself (admitted at r/rcov 1.6) — largest
re-validation; (v) do nothing (context: `rev_valence_share false` is 15x worse, median 744 kJ).
Also measured and rejected: `rev_share_onethree true` is catastrophic (median 973, max 22642).
**Operator decision pending**: adopt (ii) as a documented/default react-mode MD time-step
recommendation with `conserving`, or leave as characterized and move on.

**Earlier, superseded framing (2026-09-19 evening): `package-6` DONE — the donor rule and BOTH default flips are in.** No agent is
running in this tree. Final default state of rev-gfnff:
`rev_budget_fix_h true` (now a NO-OP under conserving), **`rev_share_form conserving`**,
**`rev_share_donor_rule true`** (new), **`rev_well_form mg`**. Five local commits
`118b289c` / `a4712e83` / `1327aaee` / `33199345` / `e688254b`; nothing pushed. Full detail:
`WORK_STATUS.md` package 6 (`Packages done: 6/6`). **Orchestrator-verified independently** (own
rebuild, md5 `8f03e0c5`; `ctest` 113/113 reconfirmed): caffeine revgfnff **-4.54694390** / gfnff
-4.67273707 Eh (the absolute-energy-shift claim); the 20-cell grid aggregate reproduces **to the
digit** — 1169 rebuilds, step max **800.12** kJ/mol, 271 events >= 50 kJ, T_max **28391.5** K,
worst cell `c2h6/T2000_f0`.

Headlines: (1) the **donor rule** (a group-13 partner, or a partner with fewer partners than its
own nominal valence, grants X_i >= 1) takes the dative/ylide regression H3N-BH3 / H3N-O / H3N-CH2
from +94.4 / +73.4 / +109.5 kcal/mol to **bit-identical with the share-off arm**, and moves NO
other falsifier; sulfoxides/phosphine oxides never needed it. (2) `conserving` is the default:
class-C adducts -87..-107 -> -1.5..+0.0 kcal/mol. (3) `mg` is the default: class-A median rms
24.49 -> 19.50, dev r90 -0.330 -> -0.058, guard 1.0341 -> 1.0439. (4) `ctest` 113/113.

**2026-09-19: `package-7` DONE — the "grid interaction" was a 20-cell sampling artefact, and the
real tail is root-caused.** One local commit (diagnostic only), nothing pushed. `WORK_STATUS.md`
package 7 (`Packages done: 7/7`), `docs/REV_GFNFF_STAGE3A.md` 2.1 rewritten.

- **The tail belongs to `conserving`, not to the combination.** Over **130 cells** (the same three
  systems, every frame, both T, one frozen binary `24a57b1c`) `gauss + conserving` reaches
  **1071.3 kJ/mol** per step and **T_max 46 601 K** — worse in the maximum than the shipped
  default's 800.1 / 28 392. Its worst cell `c2h6/T2000_f7` is simply not in the 20-cell grid. MG
  adds FREQUENCY (27 vs 16 cells above 100 kJ; on c2h6 at 2000 K 20/25 vs 9/25), not mechanism
  (0.6 kJ/mol on the corner gap; the H-H well is within 2-3.5 % of the Gaussian at every r).
- **Mechanism**: the conserving share keeps a transient geminal H2 (both H still on the same C) at
  `c = f_H^2 = 0.2531` of a 0.204 Eh well = **-136 kJ/mol**, where the delivered left-over rule
  gives that pair exactly 0. Same geometry, bond term of the bridged corner minus the unbridged
  one: **+124.4 / +128.6 (delivered) vs -10.6 / -11.2 (conserving) kJ/mol** — free instead of
  forbidden, so the molecule visits it 1.5-2.8x as often.
- **The jump is a blend-window resolution failure, not a step in the potential**: the new
  `CURCUMA_BLENDDUMP=1` shows the break transition moving `s` 0.000000 -> 0.515777 in ONE 0.25 fs
  step (window ~0.175 a0 wide in distance for H-H, ~2 steps at 2000 K; corner gap median 153, max
  390 kJ/mol). **dt 0.25 -> 0.125 -> 0.0625 fs gives 800.1 -> 38.7 -> 13.5 kJ and T_max 28 392 ->
  5 579 -> 5 441 K** on the worst cell; a discontinuity would survive a smaller step.
- **Verdict: design tension, no minimal defect, nothing changed.** Five costed options in
  `WORK_STATUS.md` 7.8; the cheapest honest mitigation is `dt <= 0.125 fs` with `conserving`.
  Both brief-suggested hypotheses were falsified by measurement (share `min` width has no effect
  over 130 cells; `gauss + conserving` keeps `w` on the well and is still worst).
- The only source change is an env-gated `blendD` log line in `FFWorkspace::updateTransitions()`.
  Verified inert: `gfnff` caffeine -4.6727370686 / benzene -2.3627255262 and `revgfnff` caffeine
  -4.5469438980 unchanged, the 20-cell grid x 4 arms identical in **80/80** cells, `ctest` 113/113.
- **Correction to package 6.4**: the `step max` of the two DELIVERED arms (391.44 / 396.88) came
  from an older binary; on the current one the same cells give **222.53 / 190.23** with the event
  counts matching to 2. The two conserving rows reproduce exactly.

Two smaller findings from package 6 still open: `ncl3_N-Cl` is the single class-A bond type the
conserving share makes worse (rms 20.10 -> 23.57, 30 of 32 bit-identical), and an MG well moves
`revgfnff`'s ABSOLUTE energy away from `gfnff` by construction (-133.6 kcal/mol on the acetic-acid
dimer, -0.56 before) while relative energies do not move.

**2026-09-18 evening: `work-packages` DONE — 5/5 packages, 12 local commits `378a13cc..64f0109f`,
`WORK_STATUS.md` (601 lines). Orchestrator-verified on the final binary (md5 bd82dff3): c2h6 cell
default 59.3 kJ/step (136 rebuilds), `rev_budget_fix_h false` 2593.6 (90), `rev_share_form
conserving` 54.4, `rev_well_form mg` 51.2; caffeine SPs unchanged.** Delivered: (1) H fix default
ON, 20-cell baseline, adduct falsifier; (2) multi-gpu merged, identity at T=1 proven,
`coulomb_implicit true` adopted (new dump_params md5s 4013d6fc / 6c3a87c8), a remote build break
fixed (`dsygst_` outside its BLAS guard), GPU code merged but not compiled here; (3) `conserving`
share built, default off — adducts −87..−107 → +2..+32, ch4_H f10 388 → 34 kJ, but dative/ylide
neutrals +73..+110 kcal/mol (H3N-BH3, amine oxide, N-ylide) → recommendation: adopt only with a
donor rule; (4) `rev_well_form gauss|mg|erfmorse`, default gauss bit-identical — class-A median
rms 24.5 → 19.5/19.2, r90 dev −0.33 → −0.06, guard 1.034 → 1.044 (holds), r_eq shifts up to 0.0067
A, MG and erf-Morse indistinguishable, MG cheaper (erf-Morse bisection cost 1.42x before caching)
→ recommendation MG if any, but the pair table (no C-C/C=C/C#C split) is why 19.5 not 2.1; (5)
tests 16/18/19 recalibrated, three new two-armed tests (20, gfnff_03, gfnff_04),
`docs/REV_GFNFF_STAGE3A.md`. **Harness finding**: `test_utils.sh` prefers `release/curcuma`; run
`CURCUMA=build_rev/curcuma ctest ...` — then 111/113. **The two failures `cli_simplemd_08/09`
(plain gfnff, acetic-acid dimer, 0.4465 Eh NVE drift) are a BLAS-less-build defect, not a branch
regression**: orchestrator rebuilt `release/` (USE_BLAS ON) from HEAD → 3/3 pass (md5 831bb4fb);
`build_rev` is USE_BLAS OFF. Method note from the agent: an algebraically equal re-association in
the bond gradient cost 1 ulp and 77 grid rebuilds at bit-identical single points — the trajectory
fingerprint catches what the energy identity cannot.

**CORRECTION (same evening): the "BLAS-less-build defect" above was the orchestrator's wrong
inference** — `release` and `build_rev` differ in FOUR options, not one, and a `build_blas/`
(= build_rev + USE_BLAS ON) fails too. `blasless-drift` (Opus, `BLASLESS_DRIFT_STATUS.md`) showed:
gradients of the two builds agree to 1e-15 Eh/A and both match FD; the NVE drift shrinks as dt^2
(0.67 / 8.7e-4 / 1.7e-4 / 6.4e-5 Eh at dt 1.0 / 0.5 / 0.25 / 0.125 fs, identical in both builds);
at dt = 1 fs the free O-H stretch (`rattle_12 false`) sits at the Verlet stability edge since the
gradient-unit fix (Known Issue #28) and the 10 ps run diverges for 4 of 30 (build_rev) resp. 6 of
30 (release) perturbed starts — a coin flip decided by 1 ulp, not a code defect. Test 09 runs test
08's trajectory bit-for-bit (`xtb-gfnff` falls back to native gfnff here), so it was n = 1.
Orchestrator confirmed release == build_rev on the trajectory to the printed digit, then set
`cli_simplemd_08/09` to `-md.time_step 0.5` (0 of 30 fail there, 17x margin): both pass against
both builds. Note: `ctest` runs the COPIES of the scripts in the build tree — re-run `cmake .`
after editing a `run_test.sh`, or the old copy is tested.

~~**Operator decisions pending**: `rev_share_form conserving`, `rev_well_form mg`, whether to
commission the donor rule.~~ **All three DECIDED and DELIVERED 2026-09-19 — see the block at the
top of this section.** The donor rule was built, both defaults were flipped. Nothing pushed.

**Earlier the same day: `work-packages` (Opus), the ONLY agent in the main tree, builds in
`build_rev/`, measures with frozen copies, commits locally per logical change, never pushes.**
Operator decisions 2026-09-18: H fix default ON; build the charge-based budget rule; build BOTH
well forms (MG and erf-Morse) and compare in real use; merges in the reviewed order; the agent
works through without hand-backs. `origin/reactff2-llm` already merged by the orchestrator
(`6530631f`; caffeine SPs and the c2h6 trajectory bit-identical). Packages, status in
`WORK_STATUS.md` (`Packages done: n/5`): (1) `rev_budget_fix_h` default true + new baseline on the
20 live cells with the three smoothness columns, all falsifiers, and the four class-C radical
adduct scans as a new falsifier; (2) merge `origin/feature/multi-gpu` with `coulomb_implicit`
pinned false, prove identity, flip it as a separate commit; (3) share mode `conserving`
(`f_i = min(1, Val_i/S_i)`, `c = f_i f_j`, excess budget by charge), default unchanged; (4)
`rev_well_form gauss|mg|erfmorse`, curvature-pinned two-parameter step, `w` off the well, default
`gauss`; (5) recalibrate `cli_simplemd_16/18/19`, new ctests, docs. It stops early only on a
failed acceptance. **A fresh session resumes it from `WORK_STATUS.md`.**

**Before that (2026-09-17).** fable-review-2 is DONE — `FABLE_REVIEW_2.md`, 360 lines, 9/9 sections;
headline results (orchestrator spot-checked the two corrections to the briefing, both hold):
(1) H-budget mechanism and "no falsifier moves" reproduced to the digit -> make `rev_budget_fix_h`
the default. (2) **The residual after the H fix is the same defect on CARBON**, not thermal: ch4_H
f10 starts as CH4 with an H in a face, Val(C) 4.53 -> 4.99, all five C-H shares 0.77 -> 1.00 in the
first step (orchestrator: c 0.773 -> 0.945 -> 0.999, Val(C) 4.98); 438 of 487 grid events >= 50 kJ
sit in the two f10 cells, the other 49 are thermal. (3) **The carbon budget is NOT the class-B
root** (corrects the 2026-09-15 briefing): on `ref/P/rkt03` Val(C) = 4.00 throughout (orchestrator:
4.000-4.012 on 4 frames), TS error -13.5 kcal/mol, capping C changes 0.00; the "78 kcal/mol" was a
pre-share estimate. The budget's real damage is on class-C approach scans: radical adducts 54-107
kcal/mol too deep (CH4+H, NH3+H, H2O+H, N2H4+H, n = 4). (4) **Proposed rule**: valence-conserving
share `f_i = min(1, Val_i/S_i)`, `c = f_i f_j`, plus excess budget by CHARGE not element (H, F
never; B/Al +1; period >= 3 groups 15-17 as today; C/N/O by positive topological group charge) —
offline: adducts +2..+32 of the reference, rkt06 2.71 -> 2.8, NH4+/H3O+ +0.17/+0.05 (0 if built
from w), ClO4-/BF4- 0.00; not covered: neutral dative N, H5O2+, N2H7+; MD effect inferred, not
run. (5) **erf wells (operator's question)**: only erf-Morse `E = -D(2y - y^2)`, `y =
erfc((r-R)/sigma)/erfc((r0-R)/sigma)` satisfies E'(r0) = 0; on all 32 class-A curves median rms
delivered 19.2 -> MG 1.32 / erf-Morse 1.36 (never > 0.12 apart), curvature-pinned 2-parameter
variant 2.1 for both; no unification with the share (fitted erfc centre lies inside r0); fit depth
on charge-frozen rests. Verdict: build 3a(iii) as the curvature-pinned step first; MG or erf-Morse
is the operator's taste. (6) multi-gpu: no OpenMP reduction in `ff_methods`, thread-count
independent; `coulomb_implicit` serves the rev corners but breaks the `dump_params` md5 yardstick;
merge order `origin/reactff2-llm` first, then decide the H fix and re-baseline, then multi-gpu with
`coulomb_implicit` pinned false. Also: the grid has 20 live cells, not 22.

Launch record:

| task | tier | what it does | status file | tree |
|---|---|---|---|---|
| fable-review-2 (DONE) | **Fable** (operator-requested 2026-09-17) | Independent review seeded ONLY with the 2026-09-15 briefing, free to read everything, read-only on the repo. (A) check the diagnosis chain link by link and recommend the rule for which elements may grow a share budget (H fix, carbon 4.95 in CH4 + H / class B, residual 222 kJ per step). (B) **the operator's question: can the bond wells be realised through erf functions** — candidates, equilibrium constraints, parameters, relation to the erf bond order and the share, offline test against class-A curves; verdict for stage 3a(iii). (C) what `origin/feature/multi-gpu` (16 commits, 3 in the GFN-FF core: on-the-fly CPU Coulomb pairs, parallel topology loops, setup speed-up; dry-run merge conflicts in `gfnff_method.cpp`, `gfnff_gpu_method_impl.h`, `main.cpp`, `AIChangelog.md`, `.gitignore`) means for rev-gfnff and what to re-measure after a merge | `FABLE_REVIEW_2.md` | frozen binary `scratchpad/fable2/curcuma_frozen` (md5 f7a37866 = HEAD 4d5287bc) |

Also noted 2026-09-17: `origin/reactff2-llm` has 3 commits this checkout lacks (`ca71d8c0` addPair,
`3f60ea31` CIF, `13eb736d` ANCOpt hold); local is 37 ahead and unpushed. No merge done yet.

**Before that (2026-09-15 evening).** Five agents closed (qp-bondstate, break-tail,
scan-cadence, verbosity-traj, h-budget) — Finished table below. Worktree `curcuma-head/` is on
branch `fix/h-valence-budget` with a clean tree (its patch is in the main tree as `35d7bb52`);
its `build/curcuma` (md5 a4b9de6e) is the h-budget agent's binary and now differs from the
worktree source, rebuild before reuse.

**Operator decisions pending**: (1) make `rev_budget_fix_h` the default; (2) which elements may
grow a share budget at all (the carbon of CH4 + H reaches 4.95 — likely the class-B root);
(3) adopt max |dEpot| per MD step outside rebuilds + hard-swap count as the smoothness
falsifier; (4) recalibrate `cli_simplemd_16/18/19` after (1).

qp-bondstate is DONE (gate not cleared; see the Finished table); its tree is committed as
`a474be89` together with the attribution correction below.

**Attribution correction (2026-09-15)**: the "3331.6 kJ event occurs bit-identically with the proxy
OFF" statement that this file, FABLE_BOND_STATE 2.3 and the vault note carried since 2026-09-14
night is **false** — it was read from run directories whose `cmd.txt` says proxy-off but whose
trajectory is the proxy-ON one (86-rebuild fingerprint). Corrected in FABLE_BOND_STATE (header),
PROXY_STATUS section 5, the fable-bondstate row below, and the vault. Rule: a `cmd.txt` is not
provenance; record the binary md5 per run and fingerprint the arm before quoting its number.

**PENDING**: `cli_simplemd_16/18/19` need a recalibration pass (their MD-behaviour thresholds
predate the stage 3a (i)+(ii) changes; see the Git-state section above for the numbers).

All six earlier agents closed on 2026-09-14; their results are in the Finished table below.

**PENDING, and deliberately deferred — the main-tree cherry-pick of the topology fix.** The
cherry-pick was attempted in the main tree and **refused** because the `proxy` agent has
`ff_methods/gfnff.h` and `gfnff_method.cpp` dirty while it works. Stashing another agent's live edits
would be harmful, so the combined state was built in an **isolated worktree** instead
(`curcuma-head/`, detached at HEAD `ceb2d160` + `76e7f83a` -> commit **`002ad100`**) and the baseline
re-measurement runs against that. **The main-tree cherry-pick still has to happen** once `proxy`
finishes: `git cherry-pick 76e7f83a`, then a full ctest (the worktree counts, 60/65 and 17/22, are
environment-limited — the worktree lacks gitignored data and the same 5+5 failures reproduce on the
pre-change binary there), and the main tree last measured 65/65 + 22/22.

**Binary discipline, and a CORRECTION to my own instruction**: the four preparation agents were told
to measure on `release/curcuma`, which I labelled "the pre-`c_ij` committed build". **That label was
wrong in a material way** — the binary is also **pre-r0-fix** (stage 3a (i), `e36d9925`).
Behaviourally proven, twice by the agents and once by me: it returns caffeine `revgfnff`
**-4.67273522** Eh, the pre-fix value, against -4.67352165 post-fix (equilibrium offset +0.0012 vs
-0.4923 kcal/mol against `gfnff`), and both agents found it reproduces the recorded *pre*-fix
class-A excesses (+15.8/+19.1 at 1.4/1.6 r_eq for ch4 C-H, where post-fix is +2.6/+5.9). Four
`src/` commits sit between its build and HEAD: `b4a3e0a4`, `4601be27`, `4c56d42a`, `e36d9925`.

**Consequences**: the `gfnff` columns are valid (gfnff is bit-identical throughout). The
**`revgfnff` columns and every class-A react-mode column are pre-r0-fix** and must be re-measured on
a HEAD build before they serve as a yardstick — and the relative-energy metrics (guard MADs, D_e)
are affected too, because the r0 fix changes the curve *shape*, which is what those metrics measure.

**Second trap found**: the embedded version string is baked at **cmake configure** time, so a `make`
rebuild after new commits keeps reporting the old `GIT_COMMIT_HASH` — `release/curcuma` self-reports
`0.0.263-66-g4ef8d30c` while containing later code. Do not trust the version string for provenance;
mtime plus a known energy value is the reliable pair. The `pairs-blend` agent handled this correctly
on its own, using a frozen post-3a(i) copy as its primary model column.


## Finished 2026-09-13/14 (committed where a hash is given)

| task | result | commit |
|---|---|---|
| r0fix | **stage 3a (i)**: `CN_i' = CN_i - cn_ij + 1` in the rev r0. 0 new parameters, the gradient chain-rule term added and FD-validated. Feedback at 1.4 r_eq over 28 class-A types **median 11.72 -> 1.24**, max 54.06 -> 36.15; X-H at 1.4 r_eq now +2.6 / -1.4 / +4.1 (were +15..+16). `gfnff` bit-identical, ctest 65/65. **Consequence: an absolute revgfnff-vs-gfnff offset of ~-0.02 kcal/mol per bond at equilibrium** — relative energies preserved to 0.001 kcal/mol, hence the guard restatement in `aef7aaf0` | `e36d9925` |
| wall-etot | `Etot = Epot + Ekin + Wall` at all five sites, `Epot` stays pure. The wall work was **19.4 % of the conserved total** in the measured row; drift over 10 ps 37.7 mEh -> 2.4e-5 Eh (orchestrator reproduced both). Wall-free runs bit-identical in every column | `4c56d42a` |
| nve-test-2 | test 18 now a 12 H2 bath, wall-free, 11500 K (46-52 rebuilds per dt, floor 1.6e-3 Eh/ps, new >=20-rebuild assertion); new `cli_simplemd_19_gfnff_rev_form_refuses_hbond` (verified to FAIL under `rev_form_switch weight`). ctest `cli_simplemd_` 22/22, `gfnff` 65/65 | `74a8a98b` |
| poly-jump | polyatomic jump baseline, 11527 rebuilds; the tail is on the BREAK side (see corrections) | `48aea7c4` |
| orca-ref | **three reference campaigns, ~4 h ORCA / 1.3 h elapsed, 16 cores.** (1) **Hirshfeld charges at all class-A points** (new class H, 63 series, 1086/1260 points, verified bit-identical coordinates/energies to class A): the reference's bonded-atom charge change r_eq -> 3.5 r_eq is **<= 0.078 e for every one of the eight drifting bonds** (C=C 0.016, O-O 0.012, N-O 0.044, C=N 0.078, N-N 0.012, N=N 0.015, O-Cl 0.040, C-O 0.058; larger 0.13-0.21 only on the halides and CO) — so the reference charges barely move while our Coulomb drifts 11-35 kcal/mol. (2) contact scans (class S, 4 systems x 20 pts: water 3.00 A/-3.84, HF 2.90/-2.74, NH3...H2O 3.25/-1.92, CH4...H2O 3.80/-0.51 kcal/mol). (3) the four missing pairs NF3/NCl3/OF2/ClF (RKS 20/20 each, r_eq and D_e in `ORCA_REF_STATUS.md`). **10 of 32 class-H UKS series failed SCF** (all 31 RKS complete); `--slowconv --uks-inside-out` recovered of2/clf 1/20 -> 16/20 | uncommitted (`_log/HIRSHFELD_CHARGES.md`, `ORCA_REF_STATUS.md`) |
| outliers | attribution of the remaining X-H residuals — section below | `c66ffa13` |
| eeq-drift | **the chi(CN) falsifier FAILED for all eight drifting bonds** — |dq_EEQ - dq_ref| = C=C 0.267, N-O 0.388, C=N 0.262, N-N 0.199, N=N 0.118, O-O 0.120, O-Cl 0.114, C-O 0.152 e (threshold 0.1; C-H passes at 0.066). EEQ moves **1.8-17.7x** the physical charge and for C=C/N-O/N-N/O-O the **wrong way**. `gfnff-fast` gives <= 5.6 kcal/mol of geometry-only Coulomb change, so **92-146 % of the drift is charge-driven**: **stage 2 owns it**, and fitting the bond D to it would fit a charge-model error. Also fixed `scripts/revgfnff_fit.py`'s `_topology_index`/`_order_points` for class S (verified: all four `ref/S/*` now select the separated frame; a 230-system comparison shows exactly 8 differences, all class S) | `EEQ_DRIFT_STATUS.md` |
| refpaths | **the acceptance gap CLOSED: 18/18 BH76 RKT hydrogen-transfer reactions now have a relaxed r2SCAN-3c NEB path** (was ~6; RKT22 excluded, it is a C5H8 isomerisation). 16 new + rkt06/rkt14 redone, one NEB attempt each, all converged, ~35 min elapsed. Layout `ref/P/<rkt>/` with per-point energies, `xi`, `s2`, band + isolated-rerun branches. **Orchestrator-verified on rkt06**: reactant -1.6666, image_5 -1.6626 -> barrier **2.51** kcal/mol at image 5 (agent: 2.52), benchmark TS 0.0001 Eh above image 5, `<S^2>` rising 0.750 -> 0.762 at the TS. **10 reactions have an interior barrier** (rkt02 1.47, rkt03 9.22, rkt04 0.84, rkt06 2.52, rkt11 4.09, rkt14 3.35, rkt18 5.98, rkt19 6.99, rkt20 6.81, rkt21 9.69); **8 are barrierless at r2SCAN-3c** (rkt01/07/08/09/10/12/16/17, maximum at the reactant endpoint) — verified as the reference surface, not the driver: OptTS at the RKT01 benchmark geometry converges in place as a genuine saddle (one imaginary mode, -1016 cm-1) while RKT10's is not even a saddle. **Those 8 support only the rms half of the acceptance criterion.** Also: the ORCA `$new_job` chain is unreliable on these radicals — every geometry is now its own job, and rkt02's chain was **11.4 kcal/mol** high on two images (its first-reported 12.6 barrier is 1.47). **Hypothesis for future work**: the class-H UKS "state instability" recorded below may be the same chain artefact, not SCF nondeterminism | `REF_PATHS_STATUS.md` |
| valfix (iter 2) | **the conflict is structural, and it is now measured.** Per-pair: the F...F contact in BF4- (1.867 A) has tight bond order **0.4739** (r/R2 1.0077, r/r0 1.246, w 0.9991) and rkt06's migrating pair **0.4985** (r/R2 1.0005, r/r0 1.260, w 0.9963) — **5 % apart in every continuous quantity the model has**, yet rkt06 needs the partner to consume a valence and BF4- needs it not to. Only the topology separates them, and the falsifiers forbid a discrete topological test. Variants built and measured: delivered (wide weight) BF4- +569.7 / rkt06 2.71 / smoothness 4 events, max 479.4; **settled weight (my literal suggestion) is NOT a fix** — the contacts stop consuming valence but then have none themselves and collapse to c ~ 0.0006, still **+473** kcal/mol; **topological 1,3 mask fixes BF4- exactly (dE = 0)** but is the discrete switch the smoothness falsifier forbids — 50 events >= 50 kJ, max **61140** kJ/mol, one cell to 8306 K — correctly **not delivered**. Remaining candidate: a **smooth** 1,3 proxy (bond-order leak onto the atom's other neighbours) needing a three-body chain rule, not yet built. **Corrections to my own claims**: BF4- at the **experimental** 1.394 A is fixed exactly (dE = 0) — my 1.143 A counterexample was an unphysical geometry; and my "settled weight" proposal was wrong. Kept: hypervalency fixed for NH4+/H3O+/ClO4-/CH5+, rkt06 2.71, bit-identity 20/20, ctest 65/65 + 22/22, and the `dcdw` gradient fix (FD 1.44e-4 -> 1.19e-8) | `VALFIX_STATUS.md` |
| valfix (iter 1) | **hypervalent valence fixed for NH4+/H3O+/ClO4- but NOT BF4-** (orchestrator-verified: +0.0006 / +0.0005 / +0.0000 vs **+569.7** kcal/mol). Form `Val_i = Val_Z + softplus_50(N_i − Val_Z)`, `N_i` = sum over **settled** bonds via `shareClip(2b−1)` — smooth, C1, 1 new global constant, 0 element-wise. rkt06 rms unchanged at **2.71**. Bit-identity 20/20 toggle combinations, `gfnff` identical. **Real bug fixed on the way**: `calcBonds`'s `dcdw` was the "sums free" form while its comment said "fixed sums", so the share's gradient was 0.4 % short — FD of the total energy at rkt06 pt10 **1.44e−4 → 1.19e−8**; worth backporting. Also found `rev_h_not_sp` registered but read nowhere (now read; `CURCUMA_REVDUMP=1` added to expose never-read PARAMs). ctest 65/65 + 22/22. **Smoothness REGRESSED: max \|dE_jump\| 10.7 → 479.4 kJ/mol, 0 → 4 events ≥ 50 kJ** (isolated to the valence rule; the `dcdw` fix alone keeps 10.7/0) — reported honestly, not hidden, and iteration 2 must either fix it or leave it standing | `VALFIX_STATUS.md` |
| topo-reuse | **the batch topology-reuse defect FIXED** (branch `fix/topo-cache-reuse`, `76e7f83a`). **Cause was not the cache**: `-batch_reuse_topology true` reuses one `GFNFF` object, but only `TopologyInfo` was ever invalidated — the force-field interaction lists (bonds/angles/torsions, bonded-vs-nonbonded partition, EEQ fragments) are built once in `initializeForceField()` and were **never rebuilt**, so every frame ran on frame 0's bond graph. Fix: `Calculation()` compares the graph the lists were built from against `perceiveGeometricBonds()` and calls `rebuildForceFieldForCurrentGeometry()` on mismatch — per-frame, not per-0.5 Bohr (H2 gains and loses a bond between two ticks of that trigger). `react`/`constant` excluded. New PARAM `-gfnff.reuse_topology_check` (default false = old), enabled automatically by `-batch_reuse_topology true`; `false` is the trust-frame-0 opt-out. **Falsifier verified by the orchestrator on the agent's binary**: the 1.4 r_eq frame is **-0.44046797 Eh in both seeding orders** (was -0.57481746 vs -0.35147217 with the opt-out, a 140 kcal/mol order dependence and 84.31 kcal/mol from the correct value). Homogeneous caffeine batch (240 frames) bit-identical, wall 0.0802 vs 0.0809 s; `gfnff` bit-identical. **Scoping was mandatory**: run globally the check also moved plain-gfnff MD. **Part B**: the react drop follows `rev_bo3_center` exactly (1.6->2.0 moves 1.80->2.25 r_eq) and `rev_tr_prebreak`; `rev_bo_break` and `react_bond_break_factor` change nothing. The drop itself is **energy-neutral**; the well is truncated by the bond-weight decay *before* it -> **stage 3a (iii)'s join radius, not extra hysteresis** | `TOPO_REUSE_STATUS.md` |
| h-budget | **DONE** (2026-09-15) — **`-gfnff.rev_budget_fix_h` (hydrogen keeps Val = 1 in the 3a(ii) share budget, derivative channel zero; DEFAULT OFF, operator decides) removes both runaways and moves no falsifier.** Per-step max abs dEpot outside rebuilds, off -> on: c2h6/T2000_f16 **2593.6 -> 59.3 kJ/mol**, T_max 62 098 -> 5 086 K, min r(H-H) 0.559 -> 1.770 a0; ch4_H/T2000_f10 **1089.0 -> 222.5** (the 391.4 at step 1 is in both arms), T_max 1.6e8 -> 8 306 K, min r(H-H) 0.124 -> 1.594 a0 — orchestrator-reproduced on the worktree binary (2593.6/22/62098 vs 59.3/17/5086) and again on the main-tree build. 22-cell grid on-arm: rebuild dE_jump max **471.1 -> 48.6**, events >= 50 kJ **3 -> 0**, hard swaps **5/478 -> 0/591**, T_max 1.6e8 -> 8 306 K; per-step max over the grid 2593.6 -> 222.5, count >= 50 kJ 143 -> 88 outside the self-destroying off-arm cell (322 -> 485 including it, because the off arm stops doing chemistry after it blows up); 12 of 22 cells identical between arms. Falsifiers: NH4+/H3O+/CH5+/ClO4- bit-identical, BF4- 1.143/1.394 A bit-identical, rkt06 rms 2.714 both with max dE 0, equilibrium 20/20 toggle 0.000000000, gfnff identical; FD on: rkt06 pt10 1.19e-8, c2h6 runaway frame 1.21e-7 Eh/A (off-arm 7.0e-6 there is FD truncation across the beta-50 softplus, dx^2 scaling shown). ctest gfnff 60/65, simplemd 17/22, same 5 in both arms (16/18/19 known; 08/09 run plain gfnff, bit-identical, pre-existing in that build). Two harness defects found and fixed by the agent (verbosity-3 filter dropped status rows with an ANSI prefix; a parallel relaunch overwrote a run dir — caught by rebuild-count fingerprint). Patch applied to the main tree | `HBUDGET_STATUS.md` |
| verbosity-traj | **DONE, FIXED** (2026-09-15) — **a numerical path was keyed on the print level.** `D4ParameterGenerator::getChargeWeightedC6` (`dispersion/d4param_generator.cpp:1414`) took its Lever-3 half-contraction fast path only at `get_verbosity() < 3`, falling through to the flat double loop at higher verbosity so a `C6_DEBUG` print could fire; the two summation orders differ by 1 ulp in C6, which react-mode MD amplifies (bit-identical for 4776 steps, 1 ulp at step 4777, topology sequence diverges at ~10495) until the cell's outcome differs (90 rebuilds / +471.1 vs 136 / +0.6). Proven by exclusion: disabling every `>= 3` gate in 12 FF files did not remove it; this inverted `< 3` gate is the only numerical one in `src/`. Fix: fast path unconditional, debug fall-through on `CURCUMA_C6DEBUG`. Default verbosity already used the fast path, so **no default result changes**; verbosity 1-4 now give the identical trajectory. Binary 79bb76bb -> e7b6c62a; ctest gfnff/sqm_val/react 95/98 with the same 3 pre-existing failures byte-identical on the unmodified source. Also from the report: count rebuilds as `max(#N)`, not `grep -c` — from calculator level 2 on the line is printed twice (GFN-FF + SimpleMD flush), which is why grep gives 135 for this cell | `VERBOSITY_TRAJ_STATUS.md` |
| runaway (orchestrator) | **DONE** (2026-09-15) — **the explosions are the hydrogen valence budget, not the scan.** Per-step status rows show Epot swinging 0.35-1.0 Eh per step and Ekin reaching 1.77 Eh (62 000 K) 1.5 fs BEFORE the +471.1 hard swap, with no REACT event in those steps; the pair moves 1.22 a0 in one 0.25 fs step. `CURCUMA_SHAREDUMP` per step: Val(H4) 1.058 -> 1.928 between two adjacent steps as the transient H-H tight order crosses the settled window, so c(C-H) 0.505 -> 0.992 and c(H-H) 0.010 -> 0.985 — all three wells full, -0.47 Eh in one step, r(H-H) collapses to 0.56 a0 (brep +1.08 Eh). ch4_H identical (Val 1.15 -> 2.00, r -> 0.78 a0). 2/2 events. A bridging H must share ONE valence; the old fixed `Val = revValence(Z)` did, hence its 10.7/0 (mechanism, not exposure — corrects BREAK_TAIL's reading). New metric needed: max abs dEpot per step outside rebuilds | `RUNAWAY_STATUS.md` |
| scan-cadence | **DONE** (2026-09-15) — **cadence is not the lever**: `revgfnff` already forces `react_check_every 1` (method_factory.cpp:387), so arms default / every-1 / every-1+disp-0.05 are bit-identical (960 / 471.1 / 3 / 5 of 478); every-2 (coarser) 1012 / 392.5 / 3 / 4 of 506 with a new +392.5 hard swap. Frozen binary (worktree build, md5 6febc06f) passed provenance (caffeine SP + arm A bit-exact vs QP_STATUS). First attempt was invalid: it shared `build_rev/` with the verbosity-traj agent, which rebuilt it mid-run — caught by the required md5; rule recorded in memory | `SCAN_CADENCE_STATUS.md` |
| break-tail | **DONE** (2026-09-15) — **the default's three large jumps are hard swaps paying the neighbours' static force constants; share and budget contribute exactly 0.** Per bond, at the step before/after each `begin_break` (all `CURCUMA_SHAREDUMP` share columns bit-identical across the swap; the swapped pair's own well already 0): c2h6/T2000_f16 bond +451.0 = C1-H4 +240.4, C1-H5 +188.8, C1-H3 +20.9 (`fc` 0.2726 -> 0.1667, ratio 1.6352); ch4_H +276.4 likewise. **Named cause**: while the transient H-H bond is in the list both hydrogens are 2-coordinate -> hybridised sp -> their C-H bonds take `bsmat[sp][sp3] = 1.3234`, C-H-H is perceived as a 3-ring -> `ringf 1.18`, and the 3-ring `fxh 1.05` fires on every C-H of that carbon; product x fqq = 1.6353, measured 1.635163 / 1.636648 (5 digits, both molecules). So the corner difference is ~64 % of each neighbouring C-H well, and it is paid at once because the pair crossed its whole break window between two scans (1.88 -> 3.00 a0 in 2 fs; the pair had been re-formed by a revert 2.0 fs earlier, 39 fs after its original formation). ch4_H's pair was formed 4.7 fs before its +310.4 break; the +51.6 is a bookkeeping break at r = 24.9 a0 in the already-exploded system. **`-gfnff.rev_valence_share false` does NOT remove the class** (+383.3 / +453.6 on the two cells, n = 2) — so the 10.7/0 of the old valence arm was exposure, not mechanism, and "0 events over 22 cells" is a weak statistic for a rare event; the mechanism-based metric is the hard-swap count. No src change, binary 79bb76bb. **Side finding, open**: `-verbosity 4` changes the c2h6 cell's trajectory (136 rebuilds / 0.6 kJ vs 90 / 471.1 at 1-3) | `BREAK_TAIL_STATUS.md` |
| qp-bondstate | **DONE** (2026-09-15) — **gate NOT cleared, Stage 2 not started; and it corrected Fable's attribution.** (A) The 3331.6 kJ event exists only with the proxy ON (`ch3nh2/T2000_f8`: off 222 rebuilds / max 2.1 kJ / 0 of 111 `begin_*` at s >= 0.99; on 86 / 3331.6 / 1 of 43 — the log Fable read was a proxy-ON run mislabelled by its `cmd.txt`; orchestrator-verified: every ad-hoc "off" dir of 14:37-14:39 carries the 86-rebuild fingerprint, the 5x5 `rep_off_*` replicates are sound, and the operator reproduced 2.1/222 vs 3331.6/86 independently). The new switch `-gfnff.rev_bo13_ordinary_join` (default off) is an **exact no-op on the default** over 22 cells (960 reb / 471.1 / 3 events / 5 of 478 hard swaps, every column identical) and in the proxy arm moves 19 -> 17 hard swaps, max 3331.6 -> 3146.0. **The delivered default's whole >= 50 kJ tail is three `begin_break` well dumps** (c2h6/T2000_f16 +471.1, ch4_H/T2000_f10 +310.4 and +51.6 kJ/mol) — the VALFIX section 6 class, not the E_over-paying ring closure Fable's 3.6 targets. (B) Offline QP over the 11 rkt06 points (constraints to 2.8e-14): path rms 2.71 (delivered) -> 3.14 at beta 0.01, 3.38 at 0.05, 3.88 at 0.1, 14.1 at 1.0 — **worse at every beta**, both TS points move further below the reference, and the QP cannot act at points 3-4/6-7 at all (one live bond, budget slack, x = 1 by the box bound) although those are the delivered path's actual worst points (+2.55 / +1.62). BF4- compressed: per-pair criterion met for beta <= 0.2 (x_BF = 1, x_FF = 0.0046), total +470.8 vs the 10-bond share-off and **+788.9 vs the pinned 4-bond evaluation** — Fable's predictions (~473 / ~790) confirmed to < 3 kcal, i.e. the total is not fixable by any share (the excess is the mis-topology's angle/Coulomb/repulsion, not the wells). One `CURCUMA_SHAREDUMP`-gated `shareD` line (per-pair well depth) added to `calcBonds`. Default path bit-identical (22-cell grid column for column; caffeine/benzene 12 digits both methods); `ctest -R gfnff` 62/65 and simplemd 19/22 = the known 16/18/19 set | `QP_STATUS.md` |
| proxy | **DONE** (2026-09-14) — **the smooth 1,3 proxy does not resolve it, and that is now a measured result.** Form: sigma = shareClip(2b-1) from the tight bond order, leak t_p = sum over the shared neighbours of sigma_ik·sigma_jk, g_p = shareClip(1-t_p), c_p = 1 - g_p(1 - (f_i+f_j)/2); continuous, C1, no threshold or count, with a three-body chain rule. **It meets every functional requirement except smoothness**: BF4- compressed **+569.7 -> +0.0000** kcal/mol (all 10 pairs read c = 1.0000), rkt06 rms **2.71** with all 11 points identical, equilibrium bit-identical 20/20 toggle combinations, `gfnff` bit-identical, ctest 65/65 + 22/22. **But smoothness is much WORSE: max |dE_jump| 3331.6 kJ/mol with 24 events >= 50 kJ** (T_max 14255 K), against 471.1 / 3 for the delivered variant and 10.7 / 0 for the build without the valence rule. Reason: the BF4- requirement forces c = 1 on EVERY pair including the six F...F 1,3 contacts, so each contact gets the full well of its pair where the plain share suppresses it; in hot react MD those wells appear and vanish. **The two falsifiers are in direct contradiction** -- BF4- demands what smoothness forbids. Delivered as **opt-in, default off** (`-gfnff.rev_share_onethree`). Four variants have now been built, each satisfying one side (wide weight / settled weight / discrete mask / this proxy), so the weight-space is exhausted and **a bond-existence state variable is the next step, not a fifth weight**. Also from this pass: the FD check found **two pre-existing chain-rule defects** (the sum channel was not scaled by g; `dcdw` was missing g) -- inert at g = 1, so the default path is unchanged -- and `reduce()` never aggregated `dEdshare` across threads (now 1.4e-17 Eh/Å agreement). The H-H handover and the contact scans are unchanged on/off (the H-H case is a 2-atom system with no 1,3 pair; at d = 2.30 Å g = 1.0000 for every pair) | `PROXY_STATUS.md` |
| fable-bondstate | **DONE** (2026-09-14), 825 lines, 6/6 sections, `FABLE_BOND_STATE.md`. **Corrects the BF4- falsifier itself**: the "share-off" reference every one of the four variants above was measured against (10-bond fresh perception at the compressed geometry, -0.78960808 Eh) is **+318 kcal/mol above the same force field's own 4-bond evaluation** pinned from the realistic 1.394 A geometry (revgfnff -1.29638269 Eh, gfnff -1.37154674 Eh) — **orchestrator-verified independently, all four numbers reproduce to the last printed digit**, including the 318 kcal/mol gap. The excess is not in the wells (the six extra F...F wells LOWER the energy) but in the mis-topology's Angle/Coulomb/lost-repulsion re-parametrisation. So variants 3/4 "pass" by reproducing a wrong reference; v2 and the new QP proposal (below) land ~+473 above it for the structurally right reason (a 1,3 contact between saturated atoms must pay two valence budgets) and that is closer to correct, not a regression. **Prior-art verdict**: no reactive FF has solved this; Tersoff/Brenner/REBO pass by tabulated coordination splines (not a real fix), MS-EVB/SCC-DFTB have the missing ingredient (a conserved quantity apportioned by a global solve, not a pairwise function of geometry). **Proposal**: replace `c_ij` with a per-corner convex QP — `x_p* = argmin sum_p D_p[-x_p - (beta/2)x_p(1-x_p)]` s.t. each atom's `sum x_p <= B_i` (one global `beta`, no graph input, no per-element data) — exact at equilibrium (proven, not measured), gradient via the envelope theorem (**no three-body chain rule at all**, `dE/dr` only needs `dD_p/dr` and `-lambda_i dB_i/dr`), 1,3 exemption falls out of budget competition. **Stress test**: rkt06 predicted <=10 (direction: TS 2.7-3.4 kcal lower at beta=0.1); BF4- per-pair correct (B-F full share, F...F exactly zero for beta<=0.2) but total still +473 above the (now-corrected) reference for the structural reason above; equilibria exact by construction. **Root cause of variants 3/4's blowup, closed form (new)**: every 394 logged events >=50 kJ in this project are hard `REACT rebuild`s, never smooth transitions; a graph-dependent exemption has amplitude ~1 full well and fan-out ~degree^2 per edge event, and smoothing only divides the force by a window-width ratio — it cannot remove the fan-out. ~~The 3331.6 kJ v4 headline event is a hard swap that occurs bit-identically with the proxy OFF~~ — **RETRACTED 2026-09-15**: read from mislabelled proxy-ON run dirs; with the proxy honoured off that cell has max 2.1 kJ and no hard swap (QP_STATUS A.1, operator-reproduced). The event is the proxy's; only the mechanism description (1,3-closure window overrun paid by E_over, `over +2976 kJ`) stands, in the proxy-ON arm. **Recommendation, decisive**: (1) operator must first decide whether to keep the literal BF4- criterion (unmeetable by anything smooth) or restate it per-pair + against the pinned 4-bond total (which the QP already satisfies structurally); (2) two cheap diagnostics before any code — replay the 3331 kJ event with only the 1,3-closure scan criterion changed (no src edit) and an offline Python QP over the existing `CURCUMA_SHAREDUMP` tables for the 11 rkt06 points at beta in {0.05,0.1,0.2}; (3) only then implement as a switchable `rev_valence_share qp` mode; (4) do NOT build a fifth clip/weight variant and do NOT turn v4 on — the reweight-existing-quantities space is provably exhausted (2.4's fan-out argument). Also flags two structural debts for later: `E_over`'s argument should read an apportioned order, not a raw topology-gated sum (it paid the +2976 kJ in one step); the settled-count budget `sig=clip(2b-1)` is itself a second geometric switch and already decides FHF- on its own. Orchestrator spot-checked one side-note (a claimed CLI flag-order sensitivity) and could **not** reproduce it — both orders give the same, correct value | `FABLE_BOND_STATE.md` |
| harness-kept | **DONE** (2026-09-14), committed `bb6079b8`. `--mode kept` now passes `-gfnff.reuse_topology_check false` (unless `--extra` names the flag, so a caller can override); `fresh` untouched; the reason is in the module docstring. **Verified**: the `gfnff` class-A column differed on **32/32 rows before** the repair (max per-cell 2770 %, median row-max 447 %, e.g. c2h6_C-C dev@1.4 47.42 against the recorded 9.12) and on **0/32 after** (max 0.0 %), byte-identical to an explicit-opt-out run; `revgfnff` + react is **byte-identical before and after** (worst delta 0.0000 over 32 bonds x 11 cells), so the primary yardstick is unaffected; `fresh` byte-identical; the 20 QUALITY exclusions unchanged. **Side finding**: `BASELINE_HEAD.md` §1c (plain `revgfnff`, no react) was recorded while `kept` was silently fresh — median dev D_e -22.254, rms 32.26 — and measures **+9.664 / 22.86** under the repaired kept protocol; corrected in place. **§1 and §2 (react) remain the valid yardstick** | `CLASSA_HARNESS.md` (appended) |
| baseline-head | **DONE** (2026-09-14) — `BASELINE_HEAD.md` (137 lines), binary = commit **`002ad100`** (HEAD `ceb2d160` incl. stage 3a (ii) + cherry-picked topo fix), md5 `ad306214…`, provenance **PASS** (caffeine revgfnff -4.673521653477, gfnff -4.672737068614). **Guard**: revgfnff moves ≤2.5 % per set vs the stale baseline; pooled conformers+S66 1.0345 → **1.0341**. **class-A** (revgfnff+react, 32 bonds): median dev D_e -25.27 (was -25.80), median dev r90 -0.318 (was -0.402), median rms 24.68 (was 26.1); worst co_CTO -138.8, o2_ODO +73.2, hcn_CTN +68.6. **Movers are entirely the r0 fix** — 26/26 comparable bonds match `R0_FIX_STATUS.md` §4's post-r0-fix values to ≤0.09 kcal, so **the 3a (ii) share is inert on this yardstick** (clean separation: the yardstick measures the bond form, not the share). ch4_C-H dev@1.4 **+15.82 → +2.59**, h2o_O-H +16.24 → -1.41, h2_H-H +13.95 → -5.93. **Class D with rev: dE_MAD 4.775 → 4.980 (+4.3 %), grad_RMS 16.305 → 16.599 (+1.8 %)** — the r0 fix buys the stretch-curve correction at ~0.2 kcal/mol cost on the class-D relative energies; the trade is clearly favourable but it is a real cost and is recorded as such. **PROTOCOL FINDING (the most consequential result here)**: the cherry-picked topo fix **changes what "topology kept across a batch" means** — a batch now re-perceives the graph per frame, i.e. `kept` silently becomes `fresh`. Orchestrator-verified on the 1.4 r_eq frame of the hcn scan: default gives **-0.44046797 in both seeding orders** (re-perception, order-independent), `-gfnff.reuse_topology_check false` gives -0.57481746 vs -0.35147217 (the old kept-frame-0 behaviour). That is why the **gfnff control failed 31/32 class-A rows** (median per-cell 49 %, max 2572 %, e.g. c2h6_C-C dev@1.4 +9.12 → +47.42): the protocol changed under it, not the force field. With the opt-out it is 0/32. **The class-A react column is unaffected** (react excludes the check), so the table in `BASELINE_HEAD.md` is valid; only `gfnff` batch measurements in `kept` mode need the opt-out. Controls that passed: 316/316 guard energies bit-identical, class-D 20/20 rows to <5e-4, class-A rev+react 31/32 (the one is a rounding artefact) | `BASELINE_HEAD.md` |
| guards | **DONE — the yardstick now exists** (2026-09-14). `GUARD_BASELINE.md`, 632 single points + 500 class-D frames, 0 errors. Per set (MAD vs published reference, kcal/mol, n): ACONF 0.155/15, ICONF 3.310/17, MCONF 0.589/51, PCONF21 1.648/18, S66 0.825/66, **all 8 conformer sets 1.492/285**. The recorded values were **cited and independently reproduced**, not re-derived (1.49 -> 1.4924; 0.83 -> 0.8252). **class D measured for the first time ever** (500 frames, 20 systems): gfnff dE_MAD **4.775**, dE_RMS 7.391, **grad_RMS 16.305 kcal/mol/A**, max frame dE 39.40; per system 0.819 (h2o_1000K) to 14.206 (h2co_2000K), grad_RMS 5.11 to 32.31; also for `-batch_reuse_topology true` (4.807 / 16.466). **revgfnff == gfnff to <=2e-4 kcal/mol on every set** — the near-constant r0 offset cancels in a relative energy, and that near-zero gap is the property future changes must not open. Caveat: the `revgfnff` column is pre-r0-fix | `GUARD_BASELINE.md` |
| classa-harness | **DONE** (2026-09-14) — `scripts/revgfnff_classa.py` (468 lines), `--binary/--method/--systems/--mode/--json/--extra`. **Verification: the metric definitions were recovered from the recorded numbers rather than chosen.** ch4_C-H and h2o_O-H **agreed exactly**; over the whole 28-bond table D_e 28/28, k 28/28, r50 27/28, r90 25/28 (28/28 on all four with the QUALITY exclusions off, so the three misses are exclusion bracket shifts). Reference side of `OUTLIER_STATUS.md` 21/21, model `fast` 7/7, `fresh` 7/7 — **model `react` 0/7**, cause = the provenance error above; it independently reproduced `R0_FIX_STATUS.md`'s *old* column 16/16. **Fresh baseline** (`revgfnff -gfnff.topology_mode react`, 32 types, own-min convention): signed D_e deviation median **-26.1**, worst **co_CTO -138.8**, then o2_ODO +73.2, hcn_CTN +68.6, h2o_O-H -63.4, ch3oh_C-O -63.1; r90 deviation median -0.402, worst -0.950 (h2_H-H); curve rms median 26.1, worst 62.4 (o2_ODO). QUALITY exclusions removed **20 points** over 32 types; no excluded radius is the largest grid point, so D_e and k are unchanged everywhere. **The react column must be re-run on a HEAD build** | `CLASSA_HARNESS.md` |
| pairs-blend | **DONE** (2026-09-14). **A.** the pair table is **32 pairs, not 21** (28 existing + NF3/NCl3/OF2/ClF + `o2_ODO`); coverage 31/32 RKS at 20/20, 26/32 UKS; **12/32 rows flagged**, every flag transcribed from `QUALITY.md`; worst signed ΔD_e: **C#O -138.75**, **O=O +73.21**, **C#N +68.60**. Only unusable reference: `of2_O-F`/`clf_F-Cl` D_e (pointwise min biased +53.1/+51.5). **B. the H-H blend case is fully diagnosed, with a lever.** The pair sits on the **non-bonded** list at 1.6 r_eq, so **`rev_bo5_center` owns it** (bo5 1.0->2.5 moves RepulsionNonbonded +13.53 -> +0.06; bo4 moves only RepulsionBonded); **neither blend switch times the handover** — that is `rev_bo3_center`, whose onset moves 1:1 with its value while bo5's centre is already saturated. Of the +13.53, **+11.40 is pure relabelling** out of RepulsionBonded (Bond bit-identical) and **+4.23 is real**; it is a **smooth ramp**, no jump, showing as a residual only because the non-bonded set is 4.2 kcal harder. **The handover cost is exactly 0.000000 for 30 of 32 bond types** at their own 1.6 r_eq — non-zero only for H-H (+13.53) and H2O2 O-O (+2.69, new). **Sweep: `rev_bo5_center` 1.3 -> 2.5 takes the H-H residual +16.28 -> +2.81** (half at 1.9), D_e unchanged at 104.88, so the effect lives only in the 1.4-1.8 r_eq band; **bo4 is not a lever** | `PAIR_TABLE.md` |
| uks-inspect | dossier for the operator's ruling. **All ten series agree with class A within 0.052 mEh and dS² <= 7.2e-5 — no state disagreement; the `--slowconv` retries are the same state.** Usable as-is: `ch3cl`, `hocl`, `ncl3`, `o2`. Not usable: `co` (no UKS state exists). Only the 16 retry points: `of2`, `clf`. For D_e only: `h2co`, `h2o2`. Far region only: `cl2`. **`of2`/`clf` pre-retry class-A trees were overwritten and are unrecoverable**, and their retry changed two things at once. `D_e` bias is **upward** (missing far UKS -> higher RKS substituted) but only for `of2` **+53.1** and `clf` **+51.5** and `co`; the four systems I flagged are unbiased. **New: five series OUTSIDE the list are state-unstable** — `c2h2_CTC`, `f2_F-F`, `hcn_CTN`, `n2_NTN`, `n2h2_NDN`, 1-4 points each where two runs of the same input land on different BS solutions (up to 55 kcal/mol, dS² up to 1.09); `ncl3` proves the mechanism (A and H `job.inp` byte-identical except the Hirshfeld line, A 20/20, H 2/20). **`revgfnff_hirshfeld.py:pick()` accepts non-converged points**, so the published UKS near-r_eq column for `h2o2`/`hocl`/`clf`/`ch3cl` is the same point as "far" (change exactly 0.000); RKS rows and the drift conclusion are unaffected | `UKS_INSPECTION.md` |
| cij | **stage 3a (ii) implemented** (5 src files, uncommitted). `c_ij = 1/2(shareClip(g_i)+shareClip(g_j))`, `g_i = (Val_i - sum_i + w)/w`, `shareClip` = smoothstep on [0,1], S = the existing `revWeight`; plus `rev_h_not_sp` (the bridging-H rule no longer fires for H in rev mode). **Falsifier met**: rkt06 path rms **15.79 -> 2.71** (target <= 10), barrier +39.80 -> +3.41 (ref +2.57). **Equilibrium and `gfnff` bit-identical** (c == 1 everywhere; caffeine -4.673521653477 Eh). Smoothness: max |dE_jump| 641.3 -> **10.70**, 47 -> **0** events >= 50 kJ. ctest gfnff 65/65. **Blocker (verified by the orchestrator): hypervalent equilibria collapse** — NH4+ **+217.6** kcal/mol and H3O+ **+130.0** kcal/mol with the share on, because nominal `Val_N` = 3 < the 4 bonds. Under repair by `valfix`.
**(WITHDRAWN 2026-09-14) a second "blocker" I reported — that rev PARAMs read into `m_rev_settings` are inert in `-sp` — was MY OWN MEASUREMENT ARTEFACT, not a defect.** The Bash tool's shell is zsh and does not word-split unquoted variables: my probe passed `-gfnff.rev_valence_share false` bundled in one variable, so flag and value travelled as a **single argv element** and the parser discarded it silently; my comparison run had them as separate literals. Verified on one binary and one path: with two argv elements `-gfnff.rev_over_preset fit2026-09-12` gives 0.87153521 Eh, bundled into one it gives the default 0.82889307. `setupRevSettings()` is reached through `loadParameterOverrides()` (`gfnff_method.cpp:3318`, called at 565 and 3278), i.e. from the constructor **and** after the CLI merge — no call-order defect exists. Independently confirmed by `sp-audit`, which classified every 2026-09-12/13 verification: **no committed conclusion changes.** What IS real and smaller: `rev_h_not_sp` was registered but read nowhere (genuinely dead), now being fixed by `valfix`; `rev_valence_share` is read at line 12062. **Lesson, and it is in my own memory file already (`curcuma-shell-and-bench-gotchas`, "zsh no word-split"): never bundle a flag and its value in a shell variable — write both as separate literal tokens.** Agent also reports: across 11 other class-B systems neutral-to-worse (hfhts 11.93 -> 22.51), and the rkt06 barrier sits at pt4 against the reference's pt10 | `CIJ_STATUS.md` |

Two findings from `nve-test-2` worth keeping: **2 H2 in NVE produces zero react events at every wall
radius 1.6-8.0 A and every temperature 1000-8000 K** (the old test system is unreachable, not merely
cold), and **the spherical wall dominated the NVE slope by ~100x** before the wall-into-Etot fix.
**Still unverified:** the react *formation* path for a genuinely fresh non-bonded pair at
non-dissociative temperature — every event observed at every setting was a break or a re-formation.

## Outlier attribution (measured 2026-09-13, `OUTLIER_STATUS.md`, recorded in roadmap WP5)

Remaining X-H residuals at 1.6 r_eq, owner assigned from the **per-term decomposition** (not the
well-shape table):

| bond | 1.6 residual | owner | mechanism -> stage |
|---|---:|---|---|
| HC-H (HCN) | +19.3 | Bond | D_e +8 % too deep **and** r90 -19 % too narrow -> **well form (iii)** |
| H-H | +16.3 | Bond **+ RepulsionNonbonded +13.5** (zero in every other mode) | the stage-1 repulsion blend hands the pair to the non-bonded branch at 1.6 r_eq -> **stage-1 switch, NOT the well** |
| H-F | -10.5 | Bond (+ Coulomb +10.5) | D_e -28 %, k -20 % -> **depth/curvature (iii)** |
| H-Cl | -24.5 | Bond | D_e -43 %, k -49 %; identical in react/rtopo/fast, i.e. the static fc -> **(iii)** |
| C-H | +5.9 | Bond | D_e -12 % yet the bond term is +9.2 over the reference -> **(iii)** |
| O-H | -2.2 | Bond/Coulomb | fine at 1.6; the react bond drop at ~1.9 r_eq truncates the well (D_e 58.2 vs 121.5) -> **tail + join radius**, part of the (iii) package |
| N-H | +7.9 | Bond | curvature/width plus the same drop at ~2.0 r_eq -> **(iii)** |

**(iii) owns five of the seven outright** and contributes to the other two. Two riders for (iii):
fix the **depth** as well as the width, and judge acceptance on the **break** side.

**Protocol caveat, quantified**: the react/fast columns depend on which frame seeded the bond graph
— 18-117 kcal/mol — because a scan whose frame 0 is already stretched does not perceive the pair as
a bond at all. The r_eq-seeded (kept-topology) number is the one a stretch scan means to measure.
Same hazard as the stale `*.topo.json` of Known Issue #11.

**One inconsistency to ignore**: the agent's own summary sentence reads "stage 3a(iii) does not own
hcn/ch4/nh3/hcl", which contradicts its measurement table D (Bond-dominated with width/depth
signatures). Table D is the measured content and agrees with the ceiling analysis; read the sentence
as a typo for "(i) does not own them".

## Results already on disk (survive everything)

- `_log/FABLE_ROADMAP_REVIEW.md` — the 2026-09-12 design review (276 lines) that set the staging
  order. Read before touching stage 2 or 3.
- `_log/CLASSA_FROZENCN.md` — frozen-CN decomposition, all class-A types (the pre-(i) baseline).
- `_log/R0_FIX_STATUS.md` — stage 3a (i) full before/after tables.
- `_log/POLY_JUMP_BASELINE.md` — polyatomic jump baseline (11527 rebuilds).
- `_log/PRESET_STATUS.md`, `REACT_JOIN_STATUS.md`, `NVE_TEST_STATUS.md`, `WALL_ETOT_STATUS.md`.
- `ref/_log/WP2_STATUS.md` — ORCA campaign (classes A 55 / C 11 / D 20 complete, E partial, B 6/15).
- `fit_work/barriers/reactions_{gfnff,revgfnff_default,revgfnff_fitted}.csv` — stage-1 barrier
  acceptance: NOT met and unreachable by stage 1 alone (the reason the staging changed).
- `fit_work/wp3_fit1/`, `wp3_fit2/` — E_over fits; `wp3_fit2` is the opt-in preset (`b4a3e0a4`).
- `fit_work/wellshape/*.md` — per-bond well shape; superseded as diagnosis by `CLASSA_FROZENCN.md`
  and the review (the "too narrow" reading only holds beyond ~1.8 r_eq).
- `jump_stats/{eeq,sqe}_t{0.95,1.0,1.05}/summary.md` — diatomic MD jumps (the NaN is explained:
  rebuild #1 has no previous-step state; reporting fixed in `4601be27`).
- `params/` — the stage-1a and 2026-09-12-fit override JSON files.

## Corrections to numbers that were circulating

- **"Topology-rebuild jumps are negligible" holds for diatomics only, and the tail is on the BREAK
  side.** Full polyatomic baseline: pooled median **0.00 kJ/mol**, 96.6 % below 1 kJ, max **423.5**
  kJ/mol against 0.7-45.3 diatomic. 58 of the 61 events >= 50 kJ are bond **BREAKS**; the three worst
  are H-H breaks in ethane (3000 K +423.5 and +394.4, 2000 K +310.3).
- **Correction to my own earlier number**: the "median 44.5 / max 69.2 kJ/mol at 1000 K" entry came
  from **n = 3 events**; over 84 events at the same temperature the median is 0.00 and the max 0.00
  (re-verified twice). It had propagated into `docs/REV_GFNFF_STAGE1.md` and is corrected there.
- **The r0(CN) feedback is larger on heavy-heavy bonds than on X-H** (pre-fix: heavy-heavy median
  12.13 / max 35.21 at 1.4 r_eq against X-H median 17.12 / max 21.57; the three largest were O=O,
  N-O, C-O). The polar heavy-heavy bonds are also where the EEQ drifts -20..-35 kcal/mol on
  dissociation, so the r0 fix alone does not straighten them — that is the charge model's.

## Open, characterised, not assigned

- **NVE drift is not measurable on 2 H2 / 3000 K in 10 ps**: all four fitted slopes lie between
  -1.7e-4 and +1.2e-4 Eh/ps with standard errors of 5-8e-5. `cli_simplemd_18` is a sanity bound, not
  a sensitive discriminator. A sensitive test needs a colder, non-reactive system.
- **The rotation projection is not conserved under confinement**: with `rm_COM`/`rmrottrans` at
  their defaults an active wall still removes ~2e-2 Eh/ps (projection off: -2.5e-7 +- 3.4e-7).
  Pre-existing and unrelated to the Etot accounting — the trajectories are bit-identical before and
  after `4c56d42a` — but NVE-with-wall measurements must set `-md.rm_COM 0 -md.rmrottrans 0`.
- **The H-H stage-1 repulsion blend** (see the outlier table): worth +13.5 kcal/mol at 1.6 r_eq and
  zero in every other bond studied, so a switch issue rather than a well issue. Not root-caused.

## Next steps, in order

1. **stage 3a (ii)** — running (`cij`). Then refit E_over (`wp3_fit2` config, ~18 s): if `p_O` stays
   ~1 Eh or `valence_O` stays ~2.5, (ii) missed something, since those fitted values are compensation.
2. **stage 3a (iii)** — the well form (MG: `phi = a x + beta x^2`, `D = s|k_b|`, `a = 1.26 sqrt(alpha)`),
   aimed per the attribution table: fix depth **and** width, move the join radius out (~2.6x) with the
   real tail, judge on the break side. Then the element-table refit (3b).
3. **stage 2** — localised integer q0 per fragment instead of the uniform rule, kappa_Z fitted on
   closed-shell charged NCI (AHB21/CHB6/IL16, n=43) + anionic SN2, S66 and the 285 conformer
   reactions as guards. Cl2-/F2- stay report-only (open-shell). Owns the Coulomb residuals.
4. **The H-H repulsion blend** — small, self-contained, the only outlier (iii) does not own.
