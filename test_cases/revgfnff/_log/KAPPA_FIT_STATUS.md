# Stage-2 SQE kappa_Z fit — run status (2026-09-22)

Binary: `build_rev/curcuma`, `-version` -> `ci-feature-multi-gpu-145-g25d89fcd` (commit `25d89fcd`,
matches current HEAD). Config: `test_cases/revgfnff/fit_work/stage2_kappa_config.json` (unmodified).
Method: `lm`, `--max-iter 40`, `--jobs 16`.

**LM iterations actually run: 1 of 40 requested.** n_evaluations=18, wall=11.6 s. The optimizer
took one damped step (lambda climbed from 1e-2 to 3.33e7 before a step was accepted at all), then
its own `relative loss change < 1e-4` check fired immediately afterwards and it stopped.

## A script bug was hit and worked around (not fixed)

The literal command in the brief (`--workdir test_cases/revgfnff/fit_work/stage2_kappa_run`,
relative) made every one of the 8 systems fail (8/8 failed, 0 evals moved anything, loss printed
as `0`). Cause: `revgfnff_fit.py` never resolves `--workdir`/`--out` to absolute paths (`type=Path`
only), but builds per-file arguments as `workdir / name` and then runs
`subprocess.run(cmd, cwd=str(workdir))`. With a relative workdir, that file argument gets the
workdir prefix applied *twice* once the child's cwd is already `workdir` — curcuma aborts with
`File not found` (confirmed by direct reproduction, SIGABRT). stdout/stderr are `DEVNULL`'d in the
script so this was invisible except as "8/8 failed, loss=0". Not touched per instructions; instead
re-ran with **absolute** `--workdir`/`--out` (same target directories), which is a correct,
already-supported CLI usage and not a script/config edit. Worth fixing in the script for next time
(resolve `args.workdir`/`args.out` early).

## Parameters (name: p0 -> final)

| element (Z) | p0 | final |
|---|---|---|
| H (1) | 0.0 | 0.0 |
| C (6) | 0.0 | 0.0 |
| N (7) | 0.0 | 6.05e-07 |
| O (8) | 0.0 | 0.0 |
| F (9) | 0.0 | 0.0 |
| Cl (17) | 0.85 | 0.849999 |

All six parameters are, for practical purposes, **unchanged from p0**.

## Dataset stats (before -> after; kcal/mol unless noted)

| dataset | n | before | after |
|---|---:|---|---|
| class E, rms_dE | 117 | 533.9568 | 533.9568 |
| class E, rms_grad (kcal/mol/A) | 117 | 1326.173 | 1326.173 |
| AHB21 MAD / RMS | 21 | 13.994 / 25.508 | 13.994 / 25.508 |
| CHB6 MAD / RMS | 6 | 1550.630 / 2161.314 | 1550.630 / 2161.314 (bit-identical) |
| IL16 MAD / RMS | 16 | 75.914 / 82.370 | 75.914 / 82.370 |
| PX13 MAD / RMS | 13 | 383.720 / 487.825 | 383.720 / 487.825 |
| BH76 MAD / RMS (not fitted, separate check) | 76 | 48.07 / 69.33 (plain EEQ) | 39.38 / 56.85 (SQE @ fitted kappa) |

Loss: 843030.1889 -> 843030.1887 (Delta = 2.4e-4, relative 2.8e-10).

Class-D guard (dE-RMS, kcal/mol): baseline 11.9002, before 11.8990, after 11.8990, limit 13.0903,
**ok=True both before and after.**

## Read

The fit did not move: loss is unchanged to 10 significant figures and every parameter sits at
(or within FD noise of) p0. CHB6's p0 MAD (1550.6 kcal/mol) matches the orchestrator's own prior
estimate (~1550), so the p0 evaluation itself is sane. The BH76 "improvement" (48.07 -> 39.38 MAD)
is **not** attributable to this campaign's fit — the "after" kappa vector is statistically p0, so
that number reflects switching the charge model from plain EEQ to SQE-with-the-pre-existing-guess
(kappa_Cl=0.85, everything else ~0), not any calibration done here.

**This is flagged as suspicious, not glossed over.** I checked whether kappa_H/kappa_N are simply
unwired for non-Cl elements: they are not — a direct single-point sweep on `nh4_nh3_pt` (class E,
NH4+/NH3 proton transfer) shows kappa_N=0 -> Single Point Energy 0.25181761 Eh vs kappa_N=2.0 ->
0.32227762 Eh, a ~44 kcal/mol swing, and kappa_H shows the identical swing (both parameters enter
the same N-H hardness sum). So there is real, sizeable sensitivity in the data the LM never
exploited. By contrast, `cl2m_Cl-Cl-` and `f2m_F-F-` (2 of the 8 class-E systems, ~34% of its 117
points) are exactly symmetric homonuclear diatomic anions: kappa_F and kappa_Cl sweeps (0.0 to 2.0)
on their equilibrium frame return **bit-identical** energies at every value tested — expected by
symmetry (no chi imbalance to drive charge transfer) but worth knowing when reading "class E"
numbers, since a third of that class is structurally insensitive to any kappa.
The likely proximate cause of the stall is that CHB6's own residual (RMS 2161 kcal/mol, ~59% of
the total loss) and class E's (34%) dominate the sum of squares and the LM's forward-difference
Jacobian on this loss surface produced a direction so poorly scaled that lambda saturated to 3.33e7
(near-pure gradient descent at a tiny step) on iteration 1, took an essentially-zero step, and then
tripped its own convergence check. This looks like an LM-implementation/tuning issue (initial
damping, step scaling relative to the CHB6-dominated loss magnitude), not evidence that kappa
calibration is infeasible — but per instructions I did not touch the optimizer, config weights, or
bounds. `override_fitted.json` was written and is kept, but it is numerically p0 and should not be
treated as a calibrated result.

## Nelder-Mead re-run (2026-09-22, same config, `--method nm --max-iter 300 --jobs 16`)

Same binary/commit/config as the LM run above (unmodified). Absolute `--workdir`/`--out`
(`stage2_kappa_run_nm` / `stage2_kappa_out_nm`) used from the start, per the known relative-path bug.

**Headline finding: NM also failed to move meaningfully off p0 — and this run gives a cleaner,
stronger diagnostic than the LM stall did.** `loss_history` in `fit_result.json` is
`[827323.5526538984, 827323.5526538984]` — **bit-identical to the last printed digit** between the
best of the initial simplex probe and the best after one full NM iteration (reflection + failed
contraction + full 6-point shrink, 8 extra evaluations). The outer convergence check
(`relative change < 1e-4`) then correctly detected zero progress and stopped after **1 of 300**
requested iterations, `n_evaluations=15`, wall 4.2 s.

## Parameters (name: p0 -> final)

| element (Z) | p0 | final |
|---|---|---|
| H (1) | 0.0 | **0.2** |
| C (6) | 0.0 | 0.0 |
| N (7) | 0.0 | 0.0 |
| O (8) | 0.0 | 0.0 |
| F (9) | 0.0 | 0.0 |
| Cl (17) | 0.85 | 0.85 |

5 of 6 parameters are bit-identical to p0. The 6th, kappa_H = 0.2, is not evidence of a directed
step: 0.2 is exactly the `scale` used to build the initial simplex's H-axis probe vertex
(`x0 + scale_H`), and the loss at that vertex is exactly the same 827323.55... that the whole
of "iteration 1" (reflection/contraction/shrink) could not improve on. So the returned optimum is
literally one corner of the initial exploratory simplex, not a point found by search.

## Dataset stats (before -> after)

| dataset | n | before | after |
|---|---:|---|---|
| class E, rms_dE | 117 | 533.9568 | 521.9934 |
| class E, rms_grad (kcal/mol/A) | 117 | 1326.173 | 1497.636 |
| AHB21 MAD / RMS | 21 | 13.994 / 25.508 | 15.965 / 26.352 |
| CHB6 MAD / RMS | 6 | 1550.630 / 2161.314 | 1549.915 / 2161.440 |
| IL16 MAD / RMS | 16 | 75.914 / 82.370 | 75.546 / 82.272 |
| PX13 MAD / RMS | 13 | 383.720 / 487.825 | 371.964 / 473.740 |

Loss: 843030.19 -> 827323.55 (the whole drop comes from the single kappa_H=0.2 probe vertex, not
from any subsequent search step). Class-D guard ok=True before and after (11.899 / 11.887, limit
13.090). Note AHB21 got slightly *worse* (13.99 -> 15.96 MAD) at that vertex — a further sign this
is an arbitrary simplex corner, not a calibrated optimum.

## BH76 before/after, redone with the NM result (same method as the LM report)

Re-verified all three points fresh (new `--evaluate-only` workdirs, no shared `.topo.json` cache):

| config | MAD (kcal/mol) | RMS (kcal/mol) |
|---|---:|---:|
| before (plain EEQ) | 48.07 | 69.33 |
| after, LM kappa (~p0) | 39.38 | 56.85 |
| after, NM kappa (~p0) | **38.87** | **56.14** |

NM's number is marginally better than LM's, but both "after" kappa vectors are statistically p0
(only kappa_Cl=0.85, the pre-set initial guess, actually differs from zero in either). The BH76
improvement (48 -> ~39 MAD) is attributable to switching `charge_model` from plain EEQ to SQE
*at the untouched initial guess*, not to anything either optimizer calibrated.

## Read: LM and NM agree

Two independently-implemented optimizers (Jacobian-based damped least squares vs. derivative-free
simplex) both returned essentially p0 on the same 6-parameter kappa_Z problem, via two different
stopping mechanisms (LM: damping saturation on an ill-conditioned trial step; NM: zero-improvement
outer check after iteration 1). Per the task's own criterion this is the **stronger signal of a
genuinely flat-for-these-optimizers loss landscape**, not two isolated optimizer bugs. Plausible
contributor (not disentangled here): CHB6's RMS ~2161 kcal/mol dominates the sum-of-squares loss by
~2 orders of magnitude over the kappa-sensitive channels (~44 kcal/mol swing on `nh4_nh3_pt`,
measured in the LM section), so a 6-D step that doesn't also fix CHB6 looks negligible to both
optimizers' convergence checks. Both `override_fitted.json` files are numerically p0-equivalent and
neither should be treated as calibrated. No optimizer code, config, or test file was modified; no
commit was made.

## Third attempt (after the alkali/alkaline-earth OverCoord fix, 2026-09-22)

Same config (`stage2_kappa_config.json`, unmodified), same method (`lm --max-iter 40 --jobs 16`),
absolute `--workdir`/`--out` (`stage2_kappa_run2` / `stage2_kappa_out2`). Binary: `build_rev/curcuma`
rebuilt from an uncommitted working-tree fix to `revValence()` (alkali/alkaline-earth over-coordination);
`-version` still reports `ci-feature-multi-gpu-145-g25d89fcd` because there is no new commit (`git
describe` is commit-derived, not content-derived) — this is expected and was not re-verified here per
the orchestrator's explicit instruction; it was already confirmed behaviorally (CHB6 sweep MAD
1550.6 -> 47.6 kcal/mol pre-fit).

**LM iterations run: 2 of 40 requested.** `n_evaluations=23`, wall time 13.9 s. Stopped on the
optimizer's own `relative loss change < 1e-4` check after iteration 2 (lambda climbed 3.33e4 ->
1.11e5, both steps accepted).

### Parameters (name: p0 -> final)

| element (Z) | p0 | final |
|---|---|---|
| H (1) | 0.0 | 0.0 |
| C (6) | 0.0 | 0.0420467 |
| N (7) | 0.0 | 0.00279194 |
| O (8) | 0.0 | 0.0 |
| F (9) | 0.0 | 0.00057574 |
| Cl (17) | 0.85 | 0.848997 |

### Dataset stats (before -> after; kcal/mol unless noted)

| dataset | n | before | after |
|---|---:|---|---|
| class E, rms_dE | 117 | 533.957 | 533.901 |
| class E, rms_grad (kcal/mol/A) | 117 | 1326.173 | 1329.920 |
| AHB21 MAD / RMS | 21 | 13.994 / 25.508 | 14.017 / 25.518 |
| CHB6 MAD / RMS | 6 | 47.634 / 66.101 | 47.558 / 66.017 |
| IL16 MAD / RMS | 16 | 75.914 / 82.370 | 75.811 / 82.316 |
| PX13 MAD / RMS | 13 | 383.720 / 487.825 | 383.688 / 487.780 |
| BH76 MAD / RMS (not fitted, separate `--evaluate-only` check) | 76 | 48.07 / 69.33 (plain EEQ) | 39.43 / 56.77 (SQE @ final kappa) |

Loss: 343004.318 -> 342931.316 (Delta = 73.0, relative 2.1e-4 — this is the criterion that
triggered the stop, since it fell just under the next iteration's threshold).

Class-D guard (dE-RMS, kcal/mol): baseline 11.9002 (same reference as the first two attempts,
i.e. this guard system is not sensitive to the OverCoord fix at these geometries), before 11.8990,
after 11.9097, limit 13.0903, **ok=True both before and after.**

### Read

**CHB6 is fixed and the fit correctly finds nothing further to do there.** Its p0 MAD is 47.63
kcal/mol, matching the orchestrator's independent post-fix confirmation (47.6) to within rounding —
this is a second, independent corroboration of the same number via the full 6-reaction sweep the
fit script runs internally, not a re-verification I performed on purpose. CHB6 barely moves under
the fit (47.63 -> 47.56), consistent with it no longer being the dominant residual.

**Did it converge to something meaningfully different from p0 this time?** Marginally, but not in a
way that matters. Unlike the first two attempts (LM: parameters bit-identical to p0 to 10 significant
figures; NM: 5/6 bit-identical, the 6th an artifact of the initial simplex probe), this run's kappa_C
moved to 0.0420467 and kappa_N/kappa_F to ~0.001-0.003 — real, LM-computed, non-trivial steps rather
than floating-point noise. But relative to the parameter bounds (0-2, scale 0.2) these are still
tiny (kappa_C at ~21% of one scale unit), the loss dropped by only 0.02%, and **every barrier/class
dataset changed by less than 0.3% MAD** (CHB6 -0.16%, AHB21 +0.16%, IL16 -0.14%, PX13 -0.01%). So
this is a small, real, but chemically negligible step, not a calibration.

**Why it stopped after 2 iterations rather than progressing:** the loss is now dominated by PX13
(MAD 384 / RMS 488 kcal/mol on n=13) and IL16 (MAD 76 / RMS 82 on n=16), not by CHB6 anymore — CHB6
went from ~59% of the total loss (before the OverCoord fix) to a minor contributor. PX13 and IL16
show no meaningful response to any of the 6 kappa_Z parameters in this run (sub-0.1% MAD change),
so the LM's `relative loss change < 1e-4` outer check fires almost immediately once the
kappa-sensitive channels (class E, the ~44 kcal/mol N/H swing documented in the first attempt) are
exhausted. This is the same structural pattern as before — a small number of huge, kappa-insensitive
residuals swamping the sum-of-squares — just with PX13/IL16 now playing CHB6's former role.

**Guard**: ok=True before and after, unchanged from the first two attempts (11.90 vs limit 13.09) —
not affected by either the OverCoord fix or this fit.

**BH76 cross-check**: MAD 48.07 -> 39.43 kcal/mol (RMS 69.33 -> 56.77), a genuine ~18% improvement
that is *larger* than any change seen in the fit's own datasets. This closely reproduces the first
two attempts' BH76 numbers (39.38 with the earlier, ~p0 LM kappa; 38.87 with NM) even though this
run's kappa vector is no longer strictly p0 (kappa_C = 0.042 now, vs ~0 before) — reinforcing the
earlier reading that the BH76 gain comes from switching `charge_model` from plain EEQ to SQE at
essentially the initial guess, not from anything calibrated in any of the three fit attempts.

**Still suspicious, flagged plainly**: PX13 (MAD 384 kcal/mol) and IL16 (MAD 76 kcal/mol) are large
enough that they — not CHB6 anymore — are now the dominant unexplained residuals in this fit's own
loss, and none of the three attempts (LM x2, NM) has moved either one by more than a fraction of a
percent. Per the task's instructions this was not investigated further (no reweighting, no bounds
change, no new optimizer) — it is reported here as the next open question for whoever picks up
stage-2 kappa calibration: PX13/IL16 look like they need either a structural fix (analogous to the
OverCoord fix that resolved CHB6) or a large weight increase before kappa_Z has any leverage over
them, not more optimizer iterations on the current loss balance.

## Fourth attempt (corrected dataset: BH76_anionic replaces PX13, 2026-09-22)

Config: `stage2_kappa_config.json` as already updated (package 16) — fitted on
`{"classes": ["E"]}` + `{"barriers": ["AHB21","CHB6","IL16","BH76_anionic"], "weight_R": 1.0}`;
`{"barriers": ["PX13","BH76"], "weight_R": 0.0}` kept as report-only. `--jobs 16 --max-iter 40
--method lm`, binary `build_rev/curcuma` (uncommitted packages 15+16 rebuilt in place). Run
directory `stage2_kappa_run4`, output `stage2_kappa_out4`, full log `KAPPA_FIT_RUN4.log`.

**3 LM iterations of 40 allowed, n_evaluations=32, wall time 20.5 s.** Per-iteration relative loss
change shrank fast — iter1 4.9e-4, iter2 1.8e-4, iter3 7.3e-6 — and the run stopped on the `relative
loss change < 1e-4` criterion at iter 3, not on max_iter or a lambda cap.

| kappa_Z (element) | p0 | final |
|---|---:|---:|
| rev.sqe_kappa.1 (H) | 0 | 0 |
| rev.sqe_kappa.6 (C) | 0 | 0.408695 |
| rev.sqe_kappa.7 (N) | 0 | 0.011652 |
| rev.sqe_kappa.8 (O) | 0 | 0 |
| rev.sqe_kappa.9 (F) | 0 | 0.000167 |
| rev.sqe_kappa.17 (Cl) | 0.85 | 0.837064 |

| dataset | n | MAD/RMS before | MAD/RMS after | fitted? |
|---|---:|---|---|---|
| class E, rms_dE (kcal/mol) | 117 | 533.957 | 533.850 | yes (weight_E=1.0) |
| class E, rms_grad (kcal/mol/A) | 117 | 1326.173 | 1373.706 | no (weight_G=0.0, reported only) |
| AHB21 | 21 | 13.994 / 25.508 | 14.143 / 25.581 | yes |
| CHB6 | 6 | 47.634 / 66.101 | 46.939 / 65.319 | yes |
| IL16 | 16 | 75.914 / 82.370 | 75.043 / 81.897 | yes |
| BH76_anionic | 16 | 69.187 / 83.667 | 70.098 / 82.519 | yes |
| PX13 (report-only) | 13 | 383.720 / 487.825 | 383.688 / 487.804 | no (weight_R=0) |
| BH76, full set (report-only) | 76 | 39.376 / 56.846 | 39.745 / 56.746 | no (weight_R=0) |

Class-D guard (dE-RMS, kcal/mol): baseline 11.9002, before 11.8990, after 11.9690, limit 13.0903,
**ok=True both before and after.**

Loss: 289524.112 -> 289328.158 (Delta = 195.95, relative 6.8e-4 cumulative over all 3 iterations).

### Read

**Did kappa_Z move meaningfully this time?** Yes, more than in any of the first three attempts —
kappa_C went from 0 to 0.409 (roughly 2x its own scale unit of 0.2, and ~20% of the 0-2 bound), a
real LM step rather than the bit-identical-to-p0 result of attempts 1/2 or the tiny 0.042 of attempt
3. kappa_N and kappa_F moved to ~0.01 and ~0.0002 respectively (negligible); kappa_H, kappa_O stayed
pinned at 0; kappa_Cl barely moved off its 0.85 starting guess (0.837).

**Does that movement look chemically sane on the two targets the design doc actually cares about?**
Not cleanly. **Class E (r2SCAN-3c Cl2-/F2- dissociation curves), the intended headline target for
kappa_Z, is essentially unchanged in the fitted quantity** (rms_dE 533.957 -> 533.850, a 0.02%
change — noise-level) **and got 3.6% worse in the unfitted gradient metric** (rms_grad 1326 -> 1374
kcal/mol/A, weight_G=0 so nothing penalized this). **BH76_anionic, the other genuinely
charge-sensitive target, got slightly WORSE by MAD** (69.19 -> 70.10, +1.3%) though marginally better
by RMS (83.67 -> 82.52, -1.4%) — a mixed, not a clean, signal. Only CHB6 (-1.4% MAD, n=6) and IL16
(-1.1% MAD, n=16) improved, both small enough in magnitude and sample size that this could plausibly
be noise rather than a calibrated effect; AHB21 and the report-only full BH76 got slightly worse.
PX13 stayed flat as expected — confirms the scoping fix worked (it no longer responds to kappa_Z,
consistent with it being charge-neutral and correctly demoted to report-only).

**Plain read: this is a fourth stall, with a different character than the first three.** The
optimizer did leave p0 substantially this time (kappa_C in particular), but that movement does not
correspond to a coherent improvement on the design doc's two intended targets — if anything it
trades a flat class-E energy for a worse class-E gradient and a slightly worse BH76_anionic MAD, in
exchange for a ~1% gain on CHB6/IL16 that may just be noise. The loss landscape near p0 along the
kappa_Z directions looks very shallow for the currently weighted objective (dominated by class E's
~534 kcal/mol rms_dE, which itself barely responds to kappa_Z here) — LM's relative-loss-change
criterion is triggering because the loss genuinely stopped moving, not because it found a good
minimum.

**Guard**: ok=True before and after (11.90-11.97 vs limit 13.09), same as every prior attempt — not
a concern in this run.

**Not acted on, per the task's instructions** (no reweighting, no bounds change, no new optimizer):
the recurring pattern across all four attempts is that class E's rms_dE — the term meant to carry
the fit — has not moved by more than 0.1% in any attempt to date, regardless of which barrier
subsets were in or out of scope. Before a fifth attempt, whoever owns this next should check whether
class E's dissociation-curve residual actually responds to kappa_Z at all in the 0-2 bound, or
whether the sensitivity requires a different parameterization/starting point — reweighting or
swapping subsets again is unlikely to help if the headline term itself is flat.

## Fifth attempt (Layer-A scoring fix + guards wired, 2026-09-22)

Implements `FABLE_REVIEW_3.md`'s Q1.4 Layer A (items 1-5, class-E scoring) and Q4 (guards),
entirely inside `scripts/revgfnff_fit.py` + `stage2_kappa_config.json` + a one-field addition to
`ref/E/cl2m_Cl-Cl-/energies.json` and `ref/E/f2m_F-F-/energies.json`. No change to
`gfnff_method.cpp` or `test_gfnff_sqe.cpp`. No git commit made.

### Item 1 -- reference-energy cap (VERIFIED against the raw data before writing any code)

Computed directly from `ref/E/*/energies.json` (script in the implementation, numbers below from
a one-off check): with `max_ref_dE_kcal=100`, **exactly 7 points drop, all in
`ahb21_21_stretch`** (r_FH = 1.8069 .. 3.5134 A, ref dE-from-that-system's-own-minimum = 145.65 ..
5129.22 kcal/mol), **and nothing else in class E** -- every other system's own maximum is under
62 kcal/mol (`f2m_F-F-` 61.09, `nh3_2_transit` 57.86, `cl2m_Cl-Cl-` 39.77, `h2o2_transit` 43.92,
`hf2_transit` 36.87, `nh4_nh3_pt` 2.81, `fch3f_umbrella` 14.49). Matches the spec's own numbers
exactly. The running fit's own printed report confirms the same count live: `class E cap:
max_ref_dE_kcal=100.0` with `ahb21_21_stretch: dropped_cap=7` and `dropped_cap=0` on all seven
other systems (see the p0 evaluation block quoted under "Acceptance checks" below).

### Item 2 -- diatomic range cap + monotonicity guard

Implemented as specified: `max_scan={"cl2m_Cl-Cl-":3.3,"f2m_F-F-":2.5}` (config, verbatim from the
task), points beyond it excluded from the fitted anchor and instead feed
`sqrt(1e2) * sum_k max(0, E_k - E_{k+1})` over the excluded tail in increasing-r order (own
`_process_class_e` in `revgfnff_fit.py`). At p0: cl2m 9 included / 11 excluded, monotonicity
guard excess exactly 0.0 kcal/mol (the excluded tail IS monotonic there); f2m 10 included / 10
excluded, guard excess 0.79 kcal/mol (there is one small non-monotonic dip in the tail -- a real,
if small, finding about the model's own long-range behaviour, not a reference comparison).

### Item 3 -- fragment anchor (diatomics) + TS residual (the rest)

**Fragment energies.** Model side: verified empirically, not assumed, that an isolated atom's
energy under `-method revgfnff -gfnff.rev_charge_model sqe` is exactly kappa-independent (Cl
atom/anion and F atom/anion bit-identical at kappa in {0, 1.5} and under plain `eeq`) -- an
isolated atom has zero bond pairs, so the SQE hardness term structurally cannot act on it. Cached
once at `FitContext` construction (`_compute_fragment_energies_model`), not per evaluation.
Values: E(Cl)=0.0, E(Cl-)=-0.96970056 Eh (matches FABLE_REVIEW_3's number exactly);
E(F)=0.0, E(F-)=-0.88314547 Eh.

**Reference side.** Cl/Cl- r2SCAN-3c single points already existed on disk
(`test_cases/revgfnff/ref/L/cl_radical`, `cl_minus`, from the TODO#3 work) -- read directly, not
recomputed. F/F- did **not** exist anywhere in the repo (grepped for the values and for
`fragment_energies`/`F-`/`f_minus` -- nothing). `orca` IS on this machine's PATH
(`/opt/orca/orca`, version 6.0.0), so per the task's own stated fallback order ("if you can find
them... if not, and orca is on PATH, compute them") **two new r2SCAN-3c single points were run**
(F atom, UKS doublet, charge 0; F-, closed-shell singlet, charge -1; same `! r2SCAN-3c TightSCF
EnGrad` keyword line as the existing Cl job.inp files, ~10 s wall each). Results: E(F)=
-99.726374651072 Eh (Hirshfeld spin 0.999999, charge 0.000007 -- clean), E(F-)=-99.829167088777
Eh (Hirshfeld charge -0.999982 -- clean). Both energies and the full provenance note are written
into `ref/E/f2m_F-F-/energies.json`'s new `fragment_energies_eh` field (and the pre-existing Cl
numbers into `ref/E/cl2m_Cl-Cl-/energies.json`'s, for symmetry/self-containedness, sourced from
`ref/L` rather than recomputed). **Sanity check against the design's own target**: E(Cl2-,
r=2.7282, this system's own reference minimum) - E(Cl) - E(Cl-) = -0.06612466 Eh = **-41.49
kcal/mol**, matching `docs/REV_GFNFF_TODO.md`'s **-41.5** to the first decimal -- this is an
independent confirmation that the fragment-anchor mechanism computes the same physical quantity
the design doc's acceptance criterion 3 specifies, not just plausible-looking code. F2- has no
prior published target in this repo; its own value at r=2.0160 (this system's reference minimum)
comes out **-49.52 kcal/mol** (recorded here as the new reference point for that side).

**TS residual, judgment call.** "the frame with the largest reference dE" is read relative to
that system's own reference MINIMUM (same anchor as the item-1 cap), not relative to frame 0 --
frame 0 is the shortest-scan-value endpoint (`_topology_index`'s rule for class A/E/L), not the
system's actual minimum, so reading it relative to frame 0 would make the residual trivially zero
whenever frame 0 is already, by construction, the frame of maximum dE relative to itself (true of
every frame relative to itself). Model and reference are both re-anchored at the SAME frame (where
the REFERENCE achieves its minimum) so this is one genuine barrier-height comparison, not a
mixed-anchor quantity. Full reasoning in `_process_class_e`'s docstring in the code. At p0 this
already does real work: `ahb21_21_stretch`'s TS residual is model=16.78 vs ref=42.20 (-25.42
kcal/mol, i.e. the model is 25 kcal too shallow at its own highest physical barrier point,
`r=0.7529`) -- exactly the kind of signal class E was supposed to carry and never did in attempts
1-4.

### Item 4 -- per-system normalisation, judgment call

`n_systems` (the denominator) counts only systems that end up CONTRIBUTING a non-zero residual
(`weight_E * system_weight_E > 0` and data available) -- **5** here (ahb21_21_stretch, cl2m,
f2m, fch3f_umbrella, nh4_nh3_pt), not all 8 loaded. Reasoning: counting the 3 report-only transits
in the denominator would shrink the 5 fitted systems' combined weight to 5/8 of `weight_E` instead
of the full `weight_E` split equally among them, contradicting the stated intent ("each system
counts equally"). Printed every evaluation as `n_systems_contributing_to_fit=5`.

### Item 5 -- per-system weight override

`system_weight_E={"h2o2_transit":0.0,"nh3_2_transit":0.0,"hf2_transit":0.0}`. Confirmed in the p0
report: those three print `weight_E_eff=0` and carry NO `TS@...` line (the append is gated on
`if w and n_sys`), while their `curve_rms_kcal` is still computed and printed (report-only, as
specified) -- 27.88 / 288.36 / 47.83 kcal/mol respectively at p0.

### Item 6 -- guards (S66 / conformers / charged_nci), report-only

`charged_nci` (AHB21+CHB6+IL16) is scored for free from the SAME per-eval fitted-barrier run --
zero extra structures, since those subsets are already loaded as a fitted dataset. `s66` and
`conformers` (8 subsets, 524 structures total, 351 reactions) are a SEPARATE reaction/group set
(`extra_guard_*`), deliberately kept OUT of the per-evaluate hot path.

**Cadence, judgment call (documented in code as well).** The class-D guard's own existing cadence
is "every `evaluate()` call, including every Jacobian FD probe" -- a pre-existing cost, left
untouched. `FABLE_REVIEW_3` Q4 explicitly asks for something cheaper for these new guards ("only
at accepted LM/NM steps and in `--evaluate-only`"), and the task text allows this as long as it
does not restructure the whole `evaluate()` call graph. Implemented via one added parameter
(`evaluate(x, want_guards=False)`) plus one extra call site inside the LM accept branch and one at
p0/final in `main()` -- NM's per-iteration simplex operations do **not** get guards (only its
final returned point does, via the p0/final call in `main()`), which is a real, stated scope
limitation of this implementation, not an oversight: the actual run uses `--method lm`, where every
accepted step DOES get guards (see the LM log lines below). A guard here can never feed back into
`ev.loss` retroactively (base residuals/loss for a given x are already finalised and cached by the
time `want_guards=True` can be requested) -- arming a real `factor`-driven penalty would need
guards computed INSIDE `evaluate()` before the loss is summed, which is the "restructure the whole
call graph" the task scoped out. All three guards in this run use `factor: null` (report-only) so
this limitation has no effect on the requested behaviour, but it is a real ceiling on what this
implementation can do if a future change wants to arm a penalty.

**Collision bug caught and fixed while implementing.** `BarrierGroup.name` was
`f"barrier_c{charge}_m{mult}"` with no other input; both `self.barrier_groups` (AHB21/CHB6/IL16/
BH76_anionic, mostly neutral singlet buckets) and the new guard-only groups (S66 + 8 conformer
subsets, ALL neutral singlet) would collide on the exact same `barrier_c0_m1.xyz`/`.jsonl`
filenames in the same workdir, run by two different `ThreadPoolExecutor` waves -- a real, silent
data-corruption risk (last-write-wins on the xyz, and two subprocess calls racing on the same
output path), not a hypothetical. Fixed by adding a `prefix` field to `BarrierGroup` (default `""`,
kept backward compatible); guard-only groups get `prefix="guard_"`.

### Acceptance checks (`--evaluate-only`, fresh absolute workdir, `--jobs 16`)

Ran clean, no crash (`test_cases/revgfnff/_log/KAPPA_FIT_CHECK5.log`):

```
loaded 8 reference systems (117 points)
loaded 148 barrier reactions (5 charge/spin batches, 234 structures)
loaded 351 guard-only reactions (2 charge/spin batches, 524 structures)
loaded 20 class-D guard systems (500 points, static topology)

p0 evaluation: loss=11586.4  failed systems=0/8  failed frames=0
    class E: n_E=  89 rms_dE= 120.3419 kcal/mol   n_G= 117 rms_grad=1326.1727 kcal/mol/A
    class E cap: max_ref_dE_kcal=100.0  n_systems_contributing_to_fit=5
      ahb21_21_stretch: dropped_cap=7  weight_E_eff=1  n_kept=13  curve_rms=32.578 kcal/mol  TS@r=0.7529: model=16.78 ref=42.20 residual=-25.42 kcal/mol
      cl2m_Cl-Cl-: dropped_cap=0  weight_E_eff=1  n_included=9/n_excluded=11 (max_scan=3.3)  anchored_rms=100.293 kcal/mol  mono_guard_excess=0.0000 kcal/mol
      f2m_F-F-: dropped_cap=0  weight_E_eff=1  n_included=10/n_excluded=10 (max_scan=2.5)  anchored_rms=153.155 kcal/mol  mono_guard_excess=0.7862 kcal/mol
      fch3f_umbrella: dropped_cap=0  weight_E_eff=1  n_kept=11  curve_rms=7.596 kcal/mol  TS@FCH=75: model=4.14 ref=14.49 residual=-10.35 kcal/mol
      h2o2_transit: dropped_cap=0  weight_E_eff=0  n_kept=11  curve_rms=27.879 kcal/mol
      hf2_transit: dropped_cap=0  weight_E_eff=0  n_kept=11  curve_rms=288.355 kcal/mol
      nh3_2_transit: dropped_cap=0  weight_E_eff=0  n_kept=11  curve_rms=47.831 kcal/mol
      nh4_nh3_pt: dropped_cap=0  weight_E_eff=1  n_kept=13  curve_rms=8.934 kcal/mol  TS@zbridge=-0.354: model=-6.36 ref=2.81 residual=-9.18 kcal/mol
    barrier AHB21     : n=  21 MAD=   13.99 kcal/mol   RMS=   25.51 kcal/mol
    barrier BH76      : n=  76 MAD=   39.38 kcal/mol   RMS=   56.85 kcal/mol
    barrier BH76_anionic: n=  16 MAD=   69.19 kcal/mol   RMS=   83.67 kcal/mol
    barrier CHB6      : n=   6 MAD=   47.63 kcal/mol   RMS=   66.10 kcal/mol
    barrier IL16      : n=  16 MAD=   75.91 kcal/mol   RMS=   82.37 kcal/mol
    barrier PX13      : n=  13 MAD=  383.72 kcal/mol   RMS=  487.83 kcal/mol
    guard (class D dE-RMS, penalized): current 11.8990 kcal/mol, baseline 11.9002, limit 13.0903, ok=True
    guard (charged_nci, reported only): n=  43 current MAD=41.7283 kcal/mol, baseline=41.0867, limit=38.9000, ok=False
    guard (conformers, reported only): n= 285 current MAD=1.5052 kcal/mol, baseline=1.5052, limit=1.6000, ok=True
    guard (s66, reported only): n=  66 current MAD=0.8238 kcal/mol, baseline=0.8238, limit=8.0000, ok=True

evals=1  wall=0.7s
```

Guard MADs are plausible: `conformers` baseline 1.505 matches the roadmap's `gfnff` 1.49 almost
exactly. `s66` baseline 0.824 looks low against the task's quoted "~7.8" -- that "~7.8" is the
roadmap's FULL `nci` category (9 subsets: ADIM6/CARBHB12/HAL59/HEAVY28/PNICO23/RG18/S22/S66/
WATER27), not S66 alone; our guard is deliberately `"barriers": ["S66"]` only, per
`FABLE_REVIEW_3`'s own sketch, and S66 (dispersion-bound small dimers) is easier than the NCI
average (halogen bonds, heavy hydrides). Not a wiring bug -- checked by re-reading exactly what
subset list the sketch specifies for this guard. `charged_nci` baseline 41.09 vs the quoted ~39 is
a small, plausible gap (revgfnff-at-kappa=0 vs plain `gfnff` are close but not identical codepaths).
`h2o2_transit`/`nh3_2_transit`/`hf2_transit` print `weight_E_eff=0` and contribute no `TS@...`
line -- confirmed not entering the loss (their curve_rms is report-only).

### A finding from item 2+3 together, verified directly against the two diatomics' own energies (IMPORTANT, read before trusting any kappa_Cl/kappa_F number below)

Before running the real fit, kappa leverage was checked directly (bypassing the fit script's
aggregate RMS, comparing curcuma's own per-frame energies at kappa=0 vs kappa=1.5 for every point
of both diatomics, react topology, exactly as `run_system()` issues it):

| system | frame | diff (kappa 1.5 - kappa 0), kcal/mol | included by max_scan? |
|---|---|---:|---|
| cl2m_Cl-Cl- | r=3.5466 | 93.15 | **NO** (excluded, > 3.3) |
| cl2m_Cl-Cl- | r=3.8194 | 111.03 | **NO** (excluded, > 3.3) |
| cl2m_Cl-Cl- | every other frame (18/20) | 0.0 (bit-identical) | -- |
| f2m_F-F- | r=2.4960 | 208.97 | **YES** (included, <= 2.5) |
| f2m_F-F- | every other frame (19/20) | 0.0 (bit-identical) | -- |

This exactly reproduces `FABLE_REVIEW_3` Q1.2's own probe table ("react, short->long: 18/20 zero
leverage, kappa acts only at r=3.55, 3.82"). The consequence, not previously spelled out this
concretely: **the `max_scan` cutoffs given in the task spec (3.3 / 2.5, chosen in Q1.4 item 2 for
a different reason -- cutting the DFT delocalisation-error plateau, not with kappa-leverage in
mind) happen to EXCLUDE both of cl2m's kappa-sensitive frames while, by a hair, INCLUDING f2m's
one kappa-sensitive frame (r=2.496 <= 2.5)**. Confirmed live in the fit script itself: sweeping
`rev.sqe_kappa.17` (Cl) from 0.0 to 1.5 leaves `cl2m_Cl-Cl-`'s `anchored_rms` bit-identical at
100.293 kcal/mol in both `--evaluate-only` runs; sweeping `rev.sqe_kappa.9` (F) the same way moves
`f2m_F-F-`'s `anchored_rms` from 153.155 to 145.245 kcal/mol. **This is not a bug in this
implementation and not a contradiction of the task -- `FABLE_REVIEW_3` explicitly predicts exactly
this as Layer A's known limitation** ("Layer A makes the loss honest but cannot give kappa_Cl/
kappa_F a lever on the bonded region of Cl2-/F2- ... a design change [Layer B], out of scope for a
fit campaign"). It is recorded here, verified rather than assumed, because it directly predicts
what the actual fit below can and cannot do: **kappa_Cl's only source of gradient in this whole
campaign is the barrier subsets (same as every prior attempt); kappa_F additionally gets one real
frame of class-E signal that no prior attempt had.** Per CLAUDE.md's rule on the anchoring trap,
this is flagged before, not after, reading the fit's outcome.

### Fit run -- completed

`--jobs 16 --max-iter 40 --method lm`, same six `rev.sqe_kappa.{1,6,7,8,9,17}` parameters and
`fixed_override: {"rev": {"charge_model": "sqe"}}` as attempts 1-4, updated config above.
Absolute `--workdir test_cases/revgfnff/fit_work/stage2_kappa_run5` / `--out
.../stage2_kappa_out5`. Full log: `test_cases/revgfnff/_log/KAPPA_FIT_RUN5.log` (quoted in full
above under "Acceptance checks" for p0; the fit itself below). A 2-iteration smoke test run
earlier (`--max-iter 2`, scratch dir, not kept) previewed the same direction of movement seen
below and is not reported separately here.

**8 of 40 requested LM iterations ran, all 8 accepted, n_evaluations=67, wall time 20.7 s.**
Stopped on the optimizer's own `relative loss change < 1e-4` criterion, same stopping mechanism
as every prior attempt -- but this time after real, monotone progress across 8 iterations, not
after 1-3.

**Loss: 11586.4 -> 10385.2 (Delta = 1201.2, relative 10.4%).** For comparison, the largest
relative loss change in any of the first four attempts was attempt 3's 2.1e-4 (0.02%) and attempt
4's 6.8e-4 (0.07%) -- this run's loss moved **~150-500x further**, and it did so over 8 real
iterations rather than converging trivially at iteration 1-3.

#### Parameters (name: p0 -> final)

| element (Z) | p0 | final |
|---|---:|---:|
| H (1) | 0.0 | 0.112713 |
| C (6) | 0.0 | 0.226373 |
| N (7) | 0.0 | 0.0 |
| O (8) | 0.0 | 0.113586 |
| F (9) | 0.0 | 0.00726141 |
| Cl (17) | 0.85 | 0.00586669 |

Four of six parameters moved by a real, non-trivial amount for the first time in this whole
campaign (H, C, O, and -- see the dedicated read below -- Cl); N stayed pinned at 0 (as in every
prior attempt); F moved but only slightly.

#### Dataset stats (before -> after; kcal/mol unless noted)

| dataset | n | before | after | fitted? |
|---|---:|---|---|---|
| class E, rms_dE | 89 (was 117; -7 cap, -21 diatomic range-excluded) | 120.342 | 117.505 | yes (per-system, see below) |
| class E, rms_grad (kcal/mol/A) | 117 | 1326.173 | **12005.179** | no (weight_G=0, reported only) |
| `ahb21_21_stretch` curve_rms | 13 (post-cap) | 32.578 | 12.733 | yes |
| `ahb21_21_stretch` TS residual (model/ref/residual) | -- | 16.78/42.20/**-25.42** | 44.54/42.20/**+2.34** | yes |
| `cl2m_Cl-Cl-` anchored_rms | 9 included | 100.293 | **100.293 (bit-identical)** | yes |
| `f2m_F-F-` anchored_rms | 10 included | 153.155 | 144.450 | yes |
| `fch3f_umbrella` curve_rms / TS residual | 11 / -- | 7.596 / -10.35 | 7.653 / -10.46 | yes |
| `nh4_nh3_pt` curve_rms / TS residual | 13 / -- | 8.934 / -9.18 | 11.418 / **-2.34** | yes |
| `h2o2_transit` / `nh3_2_transit` / `hf2_transit` curve_rms | 11 each | 27.9 / 47.8 / 288.4 | 20.7 / 45.6 / 285.7 | **no** (weight_E_eff=0, report-only) |
| AHB21 MAD / RMS | 21 | 13.99 / 25.51 | 15.99 / 27.19 | yes (worse) |
| CHB6 MAD / RMS | 6 | 47.63 / 66.10 | 46.88 / 65.59 | yes (~flat) |
| IL16 MAD / RMS | 16 | 75.91 / 82.37 | 67.14 / 77.39 | yes (better, -11.6%) |
| BH76_anionic MAD / RMS | 16 | 69.19 / 83.67 | 66.02 / 79.62 | yes (better, -4.6%) |
| PX13 MAD / RMS (report-only) | 13 | 383.72 / 487.83 | 373.89 / 477.37 | no (weight_R=0) |
| BH76, full set (report-only) | 76 | 39.38 / 56.85 | 38.33 / 55.25 | no (weight_R=0) |

Class-D guard (dE-RMS, kcal/mol): baseline 11.9002, before 11.8990, after **11.9357**, limit
13.0903, **ok=True both before and after.**

#### Guards (report-only; `factor: null`, never penalized)

| guard | n | baseline (kappa=0) | before (p0) | after (final) | limit | ok before / after |
|---|---:|---:|---:|---:|---:|---|
| charged_nci (AHB21+CHB6+IL16) | 43 | 41.087 | 41.728 | 39.338 | 38.900 | False / False (closer, still over) |
| conformers (8 subsets) | 285 | 1.505 | 1.505 | **1.640** | 1.600 | **True / False** |
| s66 | 66 | 0.824 | 0.824 | 1.333 | 8.000 | True / True |

**The conformers guard flips from `ok=True` at p0 to `ok=False` at the fitted vector** -- exactly
the neutral-chemistry side effect Q1.3/Q4 predicted and the class-D guard cannot see (class-D's
own dE-RMS, 11.90 -> 11.94, stayed comfortably under its 13.09 limit throughout, confirming again
that it is the wrong instrument for this). This is the guards' first real catch, not a null
result -- worth treating as a genuine finding, not noise: 1.640 vs the roadmap's 1.6 threshold is
a small absolute miss (2.5%), but it is a miss in a system this campaign was never testing before.

### Did kappa_Cl / kappa_F get a real gradient this time? (the task's explicit question)

**kappa_F: yes, a real but modest movement (0 -> 0.0073), and it visibly moves the ONE class-E
signal it has** -- `f2m_F-F-`'s anchored_rms tracks it (153.155 -> 144.450), consistent with the
single-frame leverage (r=2.496, 209 kcal/mol swing at kappa=1.5) measured directly before the fit
(see the finding above). This is new: none of attempts 1-4 could show this, because the old
frame-0-relative scoring gave the diatomics no meaningful metric at all.

**kappa_Cl: it moved MUCH more than in any prior attempt (0.85 -> 0.0059, a ~99.3% reduction,
vs. attempt 4's 0.85 -> 0.837, a 1.5% reduction) -- but NOT through any class-E gradient.**
Verified directly, not inferred: `cl2m_Cl-Cl-`'s anchored_rms is **bit-identical** (100.293...)
before and after the fit, to the same precision as the pre-fit kappa-sweep probe that found zero
leverage on 18/20 of its frames. So kappa_Cl's entire movement in this run comes from the fitted
barrier subsets (AHB21/CHB6/IL16/BH76_anionic -- the anionic SN2 / charged-NCI reactions), which
apparently prefer LESS Cl charge-hardness than the arbitrary 0.85 starting guess, not more. This
is a real, LM-computed signal, and a legitimate result of a barrier-driven fit -- but it says
**nothing about whether kappa_Cl=0.0059 is compatible with the Cl2- dissociation-energy target**
docs/REV_GFNFF_STAGE2.md's acceptance criterion 3 actually asks for, because that target has no
representation anywhere in this loss (Layer A structurally cannot give it one, as predicted).
Reading "kappa_Cl -> 0" as "chlorine doesn't need charge-hardness" would be the anchoring mistake
CLAUDE.md warns about -- the correct reading is "the only lever pulling on kappa_Cl in this
config is the SN2 barriers, and 0.0059 is where THAT lever alone comes to rest."

**Also worth noting**: attempts 1-4's headline positive result ("BH76 MAD 48 -> ~39 from
switching charge model to SQE at kappa_Cl~0.85") reproduces here too (BH76 39.38 -> 38.33) --
but now at kappa_Cl approx 0, not approx 0.85, with H/C/O doing more of the work instead. That
number was never actually attributable to kappa_Cl specifically; this run is further evidence of
that (already suspected in attempt 3/4's notes, now reinforced from a different starting point).

### Does the diatomic anchor score move sensibly with kappa? (the task's explicit question)

Yes for f2m (kappa_F increases -> anchored_rms decreases, i.e. moves toward the reference, in
both the isolated probe and inside the real fit). No detectable movement for cl2m at any kappa_Cl
value tried (0, 0.85, 1.5, and whatever the LM's forward-difference Jacobian probed during the
real run) -- confirmed structurally (zero leverage on 18/20 of its frames, and the 2 leverage
frames fall just outside `max_scan=3.3`), not just an absence of evidence.

### My own plain read: is this fifth attempt more trustworthy than the first four?

**Yes, but only for specific, verifiable reasons -- and it does not reach the design's stated
goal.** What is genuinely better, all independently checked above: (1) the loss actually moved by
an order of magnitude more than any prior attempt, over 8 real iterations, not 1-3; (2) the
TS-residual mechanism (item 3) produced real, interpretable, large corrections at p0
(`ahb21_21_stretch` -25.4 -> +2.3 kcal/mol, `nh4_nh3_pt` -9.2 -> -2.3) -- class E is finally
carrying chemistry-shaped signal instead of noise from fusion frames; (3) the guards caught a
real regression (conformers crossing its threshold) that none of the first four attempts' own
instrumentation (class-D guard, barrier MADs) could see, which is exactly what they were built
for and is itself evidence the new instrumentation works, independent of whether the fit's
outcome is good.

**What is NOT fixed, and was not expected to be** (per FABLE_REVIEW_3's own Layer A / Layer B
split, confirmed rather than assumed here): kappa_Cl still has no connection to the Cl2-
dissociation-energy target that motivated it; its fitted value is an artifact of which OTHER
terms happened to be in the loss, not a calibration against acceptance criterion 3. Anyone reading
`rev.sqe_kappa.17 = 0.0059` as "the calibrated chlorine hardness" would be repeating the exact
anchoring mistake this file's own governing CLAUDE.md warns about. G2 of Q3's promotion gate
("Cl2-/F2- curve SHAPE... rms <= 5 kcal/mol") is **not met and cannot be assessed from this run
at all** -- the anchored_rms numbers above (100.3 / 144.5 kcal/mol) are the honest, now-visible
size of that gap, not evidence toward or against a Layer-B fix.

**A new caveat this run surfaced that none of the first four flagged**: class E's rms_grad
(unfitted, weight_G=0, report-only) went from 1326 to **12005** kcal/mol/A -- a 9x degradation in
force quality alongside the energy improvement. Attempt 4 saw a much smaller version of the same
pattern (1326 -> 1374, +3.6%) and it was not called out prominently there; here the effect is an
order of magnitude larger and is flagged explicitly: **this fitted vector should not be used for
any gradient-consuming application (MD, optimisation) without separately checking forces** --
nothing in this task's scope (energies-only fit) constrains them.

**Net assessment**: this attempt is a real methodological improvement (the scoring is now honest,
in FABLE_REVIEW_3's own words) and it is the first attempt where the fit visibly DID something
rather than stalling -- but its output is not a calibrated kappa_Z vector and should not be
promoted. This *reinforces* Q3's existing recommendation (keep `rev_charge_model=eeq` default,
SQE opt-in) with a concrete new data point (the conformers-guard regression) rather than
overturning it. The open item is unchanged from FABLE_REVIEW_3's own conclusion: kappa_Cl/kappa_F
need Layer B (a q0-placement model change, `EEQSolver`/`GFNFF` code, an operator decision) before
a fit campaign can do anything more for them -- no amount of further reweighting inside
`revgfnff_fit.py` will manufacture a Cl2--bonded-region gradient that structurally does not exist
in the current model.

`override_fitted.json`/`fit_result.json` are written to `stage2_kappa_out5/` and kept, but per the
above should be read as "what an energies-only fit does under Layer-A scoring", not as a
production-ready parameter set.
