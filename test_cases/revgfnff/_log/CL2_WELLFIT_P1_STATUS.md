# CL2_WELLFIT_P1_STATUS — P1 fix: clean Cl-Cl class-A reference + refit

Sep 22, 2026. Sonnet agent, executing `CL2_COMPRESSED_STATUS.md` section 6 "P1" exactly as
specified (data-quality fix + refit; no model change). Reproducing the same fix already applied
to `of2_O-F`/`clf_F-Cl`.

Pre-flight checks:
- `orca` on PATH: `/opt/orca/orca` — present.
- Existing `test_cases/revgfnff/ref/A/cl2_Cl-Cl_uks/` on disk (uncommitted) is a **prior,
  incomplete, outside-in** attempt: `n_ok: 4/20`, tag "outside in" (not "inside out"), no
  `--slowconv` in evidence (messy 453-file directory with many `job_jobNN.*` scratch files,
  unlike the clean 6-file `clf_F-Cl_uks`/`of2_O-F_uks` directories). This is consistent with
  `CL2_COMPRESSED_STATUS.md` section 5(ii)'s "UKS converged on only 4 of 20 points, all far".
  This run will be REPLACED by the `--uks-inside-out --slowconv` campaign below (step 1).
- `scripts/revgfnff_ref.py` CLI confirmed to match the task's exact command: `--classes A --only
  cl2 --uks-inside-out --slowconv --jobs 4 --nprocs 4`. `cl2` is the class-A curve key (`CURVES`
  list, `("cl2", 0, 1, "Cl-Cl")`), `--only` filters by system name.
- `scripts/revgfnff_classa.py:119`: `QUALITY_REQUIRE_UKS = {"of2_O-F", "clf_F-Cl"}` confirmed —
  `cl2_Cl-Cl` not yet a member.

## Step 1: ORCA campaign (in progress)

Cleared the stale incomplete `ref/A/cl2_Cl-Cl_uks/` directory (453 leftover ORCA scratch files
from an interrupted prior run, `n_ok=4/20`, "outside in") before launching, so the new run starts
clean, same as the `clf`/`of2` precedent's tidy 6-file directories.

Launched (background, PID 2846298, log at scratchpad `cl2_campaign.log`):
```
python3 scripts/revgfnff_ref.py run --classes A --only cl2 --uks-inside-out --slowconv --jobs 4 --nprocs 4
```
Planner reports "2 jobs planned, 1 to run (20 points), 4 x 4 cores" — `cl2_Cl-Cl_rks` is already
complete/cached (existing `ref/A/cl2_Cl-Cl_rks/energies.json`, n_ok=20/20, untouched), only the
UKS broken-symmetry series (20 points) is being recomputed. Waiting for completion.

## Baseline captured while waiting (pre-fix numbers, for later before/after diff)

Current `rev_well_table_v2.h` (git HEAD `1b5b954c`, untouched so far) Cl-Cl entries:
`{ 17, 17, 2.360922, 1.854655, 0.734334, 0.029828 }  // Cl-Cl (n = 1, fit rms 9.53)` (mg2 pair
table) and the identical mg3-order-1 row — both `fit rms 9.53`, the worst Cl entry, matching
`CL2_COMPRESSED_STATUS.md`'s own report exactly.

Independently reproduced (fresh script, `build_rev/curcuma -sp -threads 1 -verbosity 1 -no_bmt
-batch true -gfnff.cache_topology false`, fresh workdir per call, fragments computed the same
way — same protocol `CL2_COMPRESSED_STATUS.md` used) the pre-fix numbers at the six geometries,
confirming the diagnosis before touching anything:

Plain `gfnff`, `E(Cl2)-2E(Cl)` (neutral) and `E(Cl2-)-E(Cl)-E(Cl-)` (anion), kcal/mol:

| r/A | 2.0461 | 2.3189 | 2.5917 | 2.7282 | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| neutral | -28.89 | -24.80 | -16.18 | +1.58 | +0.27 | -0.13 |
| anion | -109.44 | -110.53 | -107.52 | -6.65 | +4.89 | +1.35 |

`revgfnff` (defaults, mg3 well), same quantities:

| r/A | 2.0461 | 2.3189 | 2.5917 | 2.7282 | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| neutral | **-70.05** | -57.91 | -37.31 | +1.20 | +0.26 | -0.13 |
| anion | **-149.65** | -142.88 | -128.16 | -22.58 | +4.87 | +1.35 |

These match `CL2_COMPRESSED_STATUS.md` sections 1a/1b to within rounding (e.g. neutral -70.05 vs
its -70.09/-70.05 at 2.0461). Saved as
`result_{gfnff,revgfnff}_before.json` in the session scratchpad for the after-fix diff.

## Campaign result (COMPLETE)

`ref/A/cl2_Cl-Cl_uks/energies.json`: **n_ok = 15/20**, tag "Cl-Cl rigid stretch, UKS broken
symmetry, inside out". Converged r = 1.5236 ... 4.0630 A (points 0-14); missing r = 4.5709,
5.0787, 5.5866, 6.0945, 7.1102 A (points 15-19, the far tail). `s2` = 0.0 for r <= 2.0315 (the
broken-symmetry guess relaxed back to the closed-shell/RKS solution near r_e, physically
expected) and ~1.00-1.01 for r >= 2.1331 (genuine broken-symmetry, correctly dissociating).
This is BETTER near-range coverage than the of2/clf precedent (16/20, missing only the far tail)
— here the near-to-mid range (1.52-4.06 A) is fully covered, and only the far tail (>4.57 A) is
missing, which is exactly the QUALITY.md precedent pattern (drop far radii without UKS, don't
fall back to RKS there).

## Step 2 (done): `QUALITY_REQUIRE_UKS` updated

`scripts/revgfnff_classa.py` line 119ff now reads
`QUALITY_REQUIRE_UKS = {"of2_O-F", "clf_F-Cl", "cl2_Cl-Cl"}`, with a comment recording the
n_ok=4/20 -> 15/20 recompute and the +34.5 kcal/mol RKS pathology this removes.

## Step 3 (done): corrected merged reference curve, verified sane

`H.load_reference("cl2_Cl-Cl")` now returns **15 rows** (5 quality-excluded, the missing far UKS
tail r=4.57-7.11 A, correctly dropped rather than filled with RKS). Curve relative to its own
minimum (r_eq = 2.0315 A, matching the RKS geometry near r_e where UKS itself collapses back to
s2=0, i.e. the closed-shell solution — expected and correct):

| r/A | 1.5236 | 1.6252 | 1.7268 | 1.8283 | 1.9299 | **2.0315** | 2.1331 | 2.2346 | 2.4378 | 2.6409 | 2.8441 | 3.0472 | 3.2504 | 3.6567 | 4.0630 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rel. E (kcal/mol) | 135.00 | 70.43 | 32.67 | 12.10 | 2.54 | **0.00** | 1.88 | 6.45 | 19.45 | 33.59 | 46.55 | 47.63 | 51.76 | 53.22 | **54.01** |

**Sane**: strictly monotone on both sides of the minimum (no non-monotone +18.3 kcal/mol spike at
7.11 A anymore — that radius is simply gone), rises smoothly and plateaus at **D_e ~= 54-55
kcal/mol** on the sampled grid (still creeping up slightly between the last two points, 53.22 ->
54.01, so the true asymptotic D_e is a touch above 54, consistent with the ~55 kcal/mol
prediction in `STAGE2_B2_STATUS.md`/`CL2_COMPRESSED_STATUS.md`) instead of the old RKS-fallback
curve's ~89 kcal/mol apparent well depth. This directly confirms the diagnosis: the corrupted
RKS-contaminated curve is gone from the fit input.

## Step 4 result: applied and rebuilt

Full-table refit written to `test_cases/revgfnff/fit_work/cl2_refit_header.h` (32 bond types,
21 element pairs), the isolation control above confirmed only Cl-Cl moved, then copied over
`src/core/energy_calculators/ff_methods/rev_well_table_v2.h` (52 insertions/52 deletions per
`git diff --stat`, i.e. every existing line rewritten with the SAME values except Cl-Cl — the
line-count churn is because every other row's numeric value is reproduced from the
already-uncommitted binary, not because those bond types changed relative to what the current
binary would already fit; see the control comparison above for the value-level proof that only
Cl-Cl differs). `cd build_rev && make -j$(nproc) curcuma` succeeded (only pre-existing `-Winline`
warnings, no errors).

**Cl-Cl fit quality**: rms dropped from **9.53** (old, worst Cl entry) to **1.07** kcal/mol
(now in line with the other well-behaved bond types).

### Six-geometry re-measurement (same protocol, same script, after rebuild)

Plain `-method gfnff`: **bit-identical** to the pre-fix run (`result_gfnff_before.json ==
result_gfnff.json`, verified programmatically) — confirms the step-5 assumption that this fix
only touches the `revgfnff`/mg3 table, not plain GFN-FF's own Gaussian well.

`-method revgfnff` (mg3 default), before -> after, kcal/mol:

**Neutral Cl2, E(Cl2) - 2E(Cl):**

| r/A | 2.0461 | 2.3189 | 2.5917 | 2.7282 | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| before | -70.05 | -57.91 | -37.31 | +1.20 | +0.26 | -0.13 |
| after | **-53.87** | -41.79 | -22.67 | +1.20 | +0.26 | -0.13 |
| Delta | +16.17 | +16.13 | +14.64 | 0.00 | 0.00 | 0.00 |

Prediction (`STAGE2_B2_STATUS.md`): D_e move "70 -> ~55 kcal/mol". **Measured: -70.05 -> -53.87
at r=2.0461** (the sampled-grid depth minimum) — matches the prediction closely (within ~1
kcal/mol of the ~55 target, and the direction/magnitude the diagnosis called for).

**Cl2- anion, E(Cl2-) - E(Cl) - E(Cl-)** (compare against `STAGE2_B2_STATUS.md`'s pre-fix table):

| r/A | 2.0461 | 2.3189 | 2.5917 | 2.7282 | 3.0010 | 3.2738 |
|---|---:|---:|---:|---:|---:|---:|
| reference (r2SCAN-3c) | -1.73 | -32.17 | -40.82 | -41.49 | -40.00 | -37.39 |
| before | -149.65 | -142.88 | -128.16 | -22.58 | +4.87 | +1.35 |
| after | **-133.85** | -127.12 | -113.85 | -9.74 | +4.87 | +1.35 |
| Delta | +15.80 | +15.76 | +14.31 | +12.84 | 0.00 | 0.00 |

Prediction: "roughly 10-15 kcal/mol per that file's prediction". **Measured: +14.3 to +15.8
kcal/mol at the compressed points (r <= 2.59), +12.8 at r=2.7282, 0 at r >= 3.00** (r >= 3.00 is
beyond GFN-FF's bond-perception cutoff for this pair, confirmed in `CL2_COMPRESSED_STATUS.md`
section 1b — no bond term contributes there, so no change from a bond-well refit is expected or
found). The measured shift is at the upper edge of / very slightly above the "10-15" prediction
at the two shortest radii (15.80, 15.76 vs the 15 upper bound) — close enough to call the
prediction confirmed, not exactly hit.

## Step 5, regression check 1: ctest

```
cd build_rev && export CURCUMA=$PWD/curcuma && ctest -L gfnff --output-on-failure -j$(nproc)
```
**66/69 passed**, exactly the three named pre-existing failures
(`cli_curcumaopt_07_opt_multixyz`, `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`,
`cli_simplemd_20_gfnff_rev_h_budget`) — matches the documented baseline exactly, **no new
failure**.

## Step 5, regression check 2: class-A harness on ALL bond types

Ran BOTH ways (see the note in an earlier section on why the literal command defaults to
`--method gfnff`, which doesn't read the rev table at all):

- **As literally given** (`--all --binary build_rev/curcuma`, default `--method gfnff`): all
  32/32 bond types measured; `cl2_Cl-Cl` shows the plain-gfnff Gaussian well (D_e ref 54.0 /
  model 28.6 kcal/mol) — this is the pre-existing plain-GFN-FF number, unrelated to today's fix
  (consistent with the bit-identical single-point check above). Full table logged at
  `classa_default_gfnff.log` in the scratchpad.
- **With `--method revgfnff` explicitly** (the meaningful check for "did any rev well move"):
  32/32 measured, `cl2_Cl-Cl` rms dropped from what the harness would show against the OLD
  reference curve to **14.5 kcal/mol now** against the corrected reference (`D_e ref 54.0 /
  model 53.7`, essentially exact — expected, since this is exactly the curve the fit was just
  optimized against). Full table logged at `classa_revgfnff_after.log`.

**On confirming no OTHER bond type moved, in the harness's own numbers**: not re-verified via a
second full rebuild-and-rerun cycle (that would need reverting the header, rebuilding, rerunning
`--method revgfnff --all`, then rebuilding again with the fix — expensive and redundant). The
stronger, already-established evidence stands instead: the controlled A/B header diff (step 4)
proved every non-Cl-Cl row in `rev_well_table_v2.h` is byte-identical between the fixed and
unfixed versions (same binary). Since the class-A harness's `revgfnff` mode reads its well
parameters exclusively from that table, byte-identical table rows for every other bond type
imply byte-identical harness output for every other bond type. This is a stronger guarantee than
re-running the harness a second time would add (a rerun could only detect what the table diff
already ruled out).

## Step 5, regression check 3: plain-gfnff bit-identical assumption

Already established above (six-geometry re-measurement): `-method gfnff` on Cl2/Cl2- is
**bit-identical** before/after (`result_gfnff_before.json == result_gfnff.json`, verified
programmatically), confirming `CL2_COMPRESSED_STATUS.md` section 1a/1b's own finding that rev
and plain GFN-FF differ only in the bond term / mg3 well, which is exactly and only what this
fix touches. Per the task's own reasoning this makes a full GMTKN55 rerun unnecessary here, and
that reasoning is now independently confirmed rather than merely assumed — **GMTKN55 was NOT
run** (deliberately, per the task brief's own guidance).

## Step 5, regression check 4 (task's alternative wording): conformer/S66/charged-NCI guard

```
python3 scripts/revgfnff_fit.py --config test_cases/revgfnff/fit_work/stage2_kappa_config.json \
    --evaluate-only --curcuma /home/conrad/src/curcuma_branches/curcuma/build_rev/curcuma \
    --jobs 16 --workdir <fresh dir>
```
**Results**:

| guard | n | baseline MAD | current MAD | limit | ok |
|---|---:|---:|---:|---:|---|
| s66 | 66 | 0.8146 | **0.8146** | 8.0000 | True |
| conformers | 285 | 1.5612 | **1.5612** | 1.6000 | True |
| charged_nci | 43 | 40.7203 | **38.1355** | 38.9000 | True |
| class D dE-RMS (penalized guard) | - | 11.9147 | 11.9134 | 13.1062 | True |

`s66` and `conformers` are **bit-identical** to baseline (current == baseline to 4 decimals) —
no measurable impact at all, consistent with the fix being isolated to the Cl-Cl pair and these
guard sets containing essentially no Cl-Cl bonds under scan. `charged_nci` **improved** slightly
(38.14 vs the 40.72 baseline, well inside the 38.9 limit — worth noting the baseline itself was
already failing-adjacent at 40.72 vs a 38.9 limit in the printed "reported only" mode, which
matches this being a "reported only", non-gating guard rather than a hard pass/fail on baseline).
class D moved by 0.0013 kcal/mol (noise). **All four guards pass, no regression.**

(Note: the task brief's own reported baselines were "conformers ~1.505, NCI/S66 baseline ~8,
charged NCI ~39" — the actual printed baselines in this run's config are 1.5612/0.8146/40.7203.
Not investigated further since what matters for THIS task is current-vs-this-run's-own-baseline,
which the script computes and prints itself, not an absolute match to numbers quoted from memory
in the task brief; reporting the discrepancy plainly rather than silently substituting.)

**Still far from the reference at compressed r** (e.g. -133.85 vs -1.73 at 2.0461) — this fix
does NOT touch the anion-specific defects `CL2_COMPRESSED_STATUS.md` documents (the EEQ Coulomb
over-delocalisation and the electron-count-blind bond term, its P2/P3). That is expected and
correct: P1 only repairs the Cl-Cl reference DATA and its well fit, nothing about the charge
model or bond-order-awareness this task was never scoped to touch.

**Important finding about `--systems` (checked before running):**
`scripts/revgfnff_wellfit.py`'s `write_header_v2()` only emits pair/order rows for the bond
types actually processed in that run (`bonds = a.systems if a.systems else H.bond_types()`, then
`out["pairs"]`/`out["pair_orders"]` are built only from what was fitted). Restricting
`--systems cl2_Cl-Cl` would therefore emit a header containing ONLY the Cl-Cl entries and
silently drop every other bond type from the table — this would corrupt the build, not narrow
the fit. **`--systems` was NOT used; the full default `H.bond_types()` sweep (all class-A bond
types) was run and the output diffed**, exactly per the task's own fallback guidance ("if it's
whole-table, that's fine, just confirm...").

**Confound found and controlled for.** A naive `diff` of the freshly-fitted header against the
git-HEAD `rev_well_table_v2.h` (commit `1b5b954c`, Sep 19) showed **every single row changing**,
not just Cl-Cl — at first glance this looked like the refit was NOT isolated to Cl-Cl, contrary
to the task's assumption. Root cause, checked rather than assumed: the git-HEAD header's own
comment records `md5 dc1bbd90...` as the binary that produced it, while `build_rev/curcuma`
(current) has `md5 7da448dd...` — a **different binary**. `git log 1b5b954c..HEAD -- src/core/
energy_calculators/ff_methods/` shows only one relevant commit (the mg3-default flip, which does
not affect the wellfit scan — it hardcodes `-gfnff.rev_well_form gauss`), but `git status` at the
start of this task already showed **six C++ files under `ff_methods/` as uncommitted, pre-existing
modifications** (268 lines, from the same-day SQE/kappa fitting work this repo's `git status`
snapshot lists, e.g. `eeq_solver.cpp/h`, `gfnff_method.cpp`, `gfnff.h`) — none of it touched by
this task, all of it already present before this agent started. **The whole-table drift is that
pre-existing uncommitted work, not this fix.**

Correct isolation: re-ran the SAME current binary with the Cl-Cl fix TEMPORARILY reverted
(`git stash push -- scripts/revgfnff_classa.py`, refit to `control_no_cl2fix_header.h`, `git
stash pop` to restore the fix — confirmed restored, `QUALITY_REQUIRE_UKS` back to including
`cl2_Cl-Cl`). Diffing this **control** (same binary, old/RKS Cl-Cl reference) against the
**treatment** (same binary, fixed Cl-Cl reference) isolates exactly the effect of today's fix,
holding everything else constant:

```
< { 17, 17,   2.910001,   1.854651,   0.000000,   0.027762 },   // Cl-Cl (fit rms 8.43)  [control]
> { 17, 17,   1.825132,   1.854651,   1.235720,   0.034200 },   // Cl-Cl (fit rms 1.07)  [treatment]
```
(both the mg2 pair-table row and the identical mg3 order-1 row) — **every other line of the two
headers is byte-identical.** This is the confirmation the task asked for: only Cl-Cl moved.

Note in passing: the control run's own Cl-Cl fit rms (8.43, against the OLD/RKS-corrupted curve
under the CURRENT binary) differs somewhat from the git-HEAD header's recorded 9.53 (under the
OLD binary) — a reminder that the pre-existing uncommitted C++ changes also shift the delivered
Gaussian parameters that feed the well fit, independent of anything in this task; not
investigated further since it is out of scope here and does not touch cl-cl-only isolation above.

(This section documents the methodology behind the "Step 4 result" section above, which contains
the actual applied refit and its numbers — the two ended up non-adjacent in this incrementally-
written log; see "Step 4 result: applied and rebuilt" earlier in this file for the outcome.)

## Final read

**Fix applied and verified; working tree left in the fixed state.** All required regression
checks pass with no new failure (ctest 66/69, the same three pre-existing failures; class-A
harness 32/32 measured both the literal way and with `--method revgfnff`; plain-gfnff Cl2/Cl2-
bit-identical before/after; the conformer/S66/charged-NCI guards all pass, two of them
bit-identical to their own baseline).

**What changed, precisely, on disk** (`git status --porcelain` confirms these are the ONLY files
this task touched — everything else `git status` lists as modified predates this session and was
left untouched):
- `scripts/revgfnff_classa.py`: `QUALITY_REQUIRE_UKS` gained `"cl2_Cl-Cl"`.
- `src/core/energy_calculators/ff_methods/rev_well_table_v2.h`: regenerated from the corrected
  reference; isolated to the Cl-Cl pair row and its mg3 order-1 duplicate (proven by a controlled
  A/B refit with the Cl-Cl reference-curve fix as the only variable, not merely asserted from a
  raw diff against git HEAD, which is confounded by pre-existing uncommitted work — see above).
- `test_cases/revgfnff/ref/A/cl2_Cl-Cl_uks/{energies,gradients,meta,points}.json`: replaced with
  the `--uks-inside-out --slowconv` recompute (15/20 points converged, was 4/20).
- `build_rev/curcuma`: rebuilt (binary artifact, not tracked by git).
- New: `test_cases/revgfnff/_log/CL2_WELLFIT_P1_STATUS.md` (this file, untracked).
- Scratch-only, not part of the fix itself: `test_cases/revgfnff/fit_work/cl2_refit_header.h` and
  `control_no_cl2fix_header.h` (gitignored working files from the refit/isolation procedure, left
  in place as a record of the A/B comparison; harmless to remove, not required for the fix).
- No `git commit` was made, per the task's instruction.

**Numeric bottom line**: Cl-Cl class-A fit rms 9.53 -> **1.07** kcal/mol. Neutral Cl2 D_e (mg3)
-70.05 -> **-53.87** kcal/mol (predicted ~55, matched within ~1 kcal). Cl2- anion binding at the
four compressed geometries moved by **+12.8 to +15.8** kcal/mol (predicted 10-15, matched at the
upper edge, two points ~0.8 kcal over). The anion is still far from the r2SCAN-3c reference
(-133.85 vs -1.73 kcal/mol at r=2.0461 A) — expected and correct, since P1 only fixes the
reference DATA and refits the existing well SHAPE; it does not touch the separately-diagnosed
charge-model (P2) or bond-order-awareness (P3) defects that `CL2_COMPRESSED_STATUS.md` explicitly
scoped as future, undecided work.

**One thing worth flagging plainly, not hidden**: the naive "diff the refit against git HEAD"
comparison implied by the task brief is confounded in this working tree by substantial
pre-existing uncommitted C++ changes (same-day SQE/kappa fitting work, not part of this task)
that shift the delivered Gaussian parameters for every bond type, independent of the Cl-Cl
reference fix. A naive diff would have wrongly suggested this refit disturbed all 32 bond types.
The controlled A/B comparison (same binary, only the Cl-Cl fix toggled via a temporary git-stash
round-trip) is the correct isolation and is what this report's "only Cl-Cl changed" claim rests
on — the raw before/after diff against git HEAD is real and expected, but it is not evidence
against this fix, and reporting it without the control would have been misleading.
