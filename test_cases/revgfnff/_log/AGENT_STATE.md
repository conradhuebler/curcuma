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

## Current state (2026-09-24, plain-GFN-FF fragment-charge default just flipped)

Stage 3a/3b (bond term) is complete and shipped (see the default table below); work has moved to
**stage 2, the split-charge charge model** (`docs/REV_GFNFF_STAGE2.md`). Package 14 fixed the
uniform-q0 limitation the design doc had already diagnosed but not applied, and wired
`scripts/revgfnff_fit.py` to fit `rev_sqe_kappa_{H,C,N,O,F,Cl}` (a `fixed_override` config key for
the non-fitted `rev.charge_model: "sqe"` setting was the only gap — class E and the GMTKN55
`barriers` mechanism for AHB21/CHB6/IL16/PX13 were already generically wired). Both verified
end-to-end at p0 (mechanics only, not yet fitted) — see `WORK_STATUS.md` §14.

**Five kappa_Z fit attempts, a Fable design review, and an implemented model change ("B2") all
ran this session — no calibrated kappa_Z exists, but the open question has moved from stage 2 to
stage 1/3a.** Attempts 1-2 stalled on a real stage-1 bug (package 15: bare alkali/alkaline-earth
cations mis-read as grossly over-coordinated in cation-pi contacts — fixed, CHB6 MAD
1550.6 -> 47.6 kcal/mol). Attempt 3 stalled on a wrongly-scoped dataset (package 16: PX13 is
neutral proton-transfer, not anionic SN2 — split out `BH76_anionic`, the real 16-reaction target).
A Fable review (`FABLE_REVIEW_3.md`) found class-E scoring itself broken two ways (94% of the
loss was 4 unphysical fused-atom frames; separately Cl2-/F2- had near-zero kappa leverage under
`uniform` q0). A Sonnet agent's Layer-A fix + guards gave attempt 5 real movement (loss -10.4%,
package 17) but kappa_Cl still couldn't reach its target. **Operator chose B2** (package 18): an
Opus agent added a second q0 rule (`rev_sqe_q0_rule=mu`, new default — places charge by chemical
potential, gives kappa a real lever where `uniform` gave none) and a second `kappa(b)` form
(`rev_sqe_kappa_form=vanishing`, opt-in — the only form that protects delocalised systems like
carboxylates, at the cost of giving back most of the Cl2- gain). Result: the point-target IS hit
(-41.49 vs -41.5) but shown to be largely coincidental; the curve-shape target is proven
**structurally unreachable by any kappa_Z** — the remaining ~104 of ~148 kcal/mol error at
compressed Cl2- cannot be in the charge model (the hardness term can cancel at most ~44 kcal/mol).
**This redirects the open question to a stage-1/3a term decomposition of Cl2-, not a sixth fit.**
One real defect found, not yet fixed: `mu`'s hard argmin produces a force cusp (2.4e-3 Eh/A,
formate probe) — must be fixed before any kappa>0 MD. `docs/REV_GFNFF_STAGE2.md` fully reconciled
with all of this. SQE stays opt-in (Fable's Q3), five-point gate for reconsidering unchanged.

**Packages 19-22 (2026-09-22/23, compressed — fully superseded by 23-27 below): diagnosing the
Cl2-/F2- problem before building the fix.** §19 found 38 of a 106 kcal/mol compressed-geometry
error IS still Coulomb (kappa can't reach it) and the rest an electron-count-blind bond term;
separately, the r2SCAN-3c reference every kappa number this session chased has a large
self-interaction-error tail, and the neutral Cl-Cl well-fit data was RKS-contaminated — 3
proposals (P1 data-fix/P2 charge self-consistency/P3 bond-term redesign), none shipped yet. §20
executed P1 (clean UKS Cl-Cl refit; neutral Cl2 D_e -70.1->-53.9 kcal/mol). §21 ran a 46-job
DLPNO-CCSD(T) campaign that confirmed and quantified the SIE: true D_e(Cl2-)=-28.4, D_e(F2-)=
-26.8 kcal/mol, both 32-46% shallower than the r2SCAN-3c target — this became every later
package's reference. §22 re-measured P1+B2 against it (full-grid rms 87.5/60-65, the baseline
P2+P3 had to beat). Full detail: `CL2_COMPRESSED_STATUS.md`, `CL2_WELLFIT_P1_STATUS.md`,
`CL2F2_CCSDT_STATUS.md`, `WORK_STATUS.md` §§19-22.

**Package 23 (2026-09-23): P2+P3 implemented — the double-counting is resolved for the tested
case.** Opt-in flags `-gfnff.rev_sqe_phase1`/`-gfnff.rev_excess_electron` (both off by default,
bit-identical then) move Cl2-/F2- against the DLPNO-CCSD(T) target from rms 87.6/64.9 to
**11.5/11.7 kcal/mol full-grid** (2.0/1.0 on bonded points), at kappa_Z=0 — no fit needed.
Design: the 2c-3e resonance now lives ONLY in the bond well (reasoned, not assumed — EEQ's
delocalisation has the wrong element AND r-trend); cost: X2- charges become broken-symmetry
(-1,0). Zero regression at kappa=0 (1379 frames + 32 bond types checked). **Real open issues**:
react mode breaks at kappa_Cl=0 (-100 kcal/mol collapse past r=3.5 A, needs kappa_Cl>0 or a
proper fix); F2-'s double-counting is only *partly* resolved (well absorbs an unrelated genuine
Coulomb term too); **a NEW, unrelated bug found**: energy-only calls use a stale cached CN in the
Coulomb self-energy — the long-standing "1.77e-2 Eh/A Cl2- gradient residual" turns out to be
mostly this bug (true value 2.41e-5), not an inherent 2c-3e difficulty; not fixed, own follow-up
needed. One outdated ctest assertion (targeting the now-refuted -41.49 r2SCAN target) retired;
ctest back to 66/69. Docs folded into `docs/REV_GFNFF_STAGE2.md` by the orchestrator. Full
detail: `P2P3_STATUS.md`, `WORK_STATUS.md` §23.

**Package 24 (2026-09-23/24): the operator flagged package 23's design as dangerous — analysed,
alternatives tested, partly fixed.** Putting the whole 2c-3e resonance in the well and none in
Coulomb (`kappa_x=100`) forces an isolated X2- to broken-symmetry charges (-1,0). Confirmed real,
two separate failure modes: (1) **react-mode collapse** of -107/-91 kcal/mol (Cl2- break/form)
and -203/-169 (F2-), independent of `kappa_x` value and scan step — a bookkeeping-consistency
bug, not a `kappa_x`-scale one; **fixed** via `-gfnff.rev_excess_react_consistent` (default TRUE
within P3), react rms 28->1.8 / 48->2.7 kcal/mol, static results bit-identical. (2)
**Broken-symmetry charges mislead a real neighbour**: a water probe near one end of X2- vs. the
other differs 8.9-14.8 kcal/mol vs. DLPNO-CCSD(T); the shipped design's mean label-gap (11.84) is
actually WORSE than plain `sqe` (3.05) or `eeq` (6.15) baselines — **NOT fixed**. Two
alternatives tested and rejected: moderate `kappa_x` (no sweet spot, react mode still worse), a
symmetric `frac` correction (fails — only removes half the delocalisation P2 already fixed
elsewhere, and amplifies real asymmetry by `1/(1-c)`). **Recommendation, endorsed**: keep
`kappa_x=100` with the react repair, keep P3 opt-in, do not use P3 for any X2- expected to
interact with its environment — the real fix (a non-self-consistent Harris-like correction) is
sketched, not built. `ctest -L gfnff` 66/69 unchanged. Full detail:
`P2P3_ALTERNATIVES_STATUS.md`, `WORK_STATUS.md` §24.

**Package 25 (2026-09-24): the section-6 "real fix" built as `-gfnff.rev_excess_mode harris`.**
Charges relax freely (bit-identical to `flat` at `kappa_x=0`); a new additive, non-self-
consistent term `x*g(r)` (own header `rev_harris_table.h`, `rev_well_table_v2.h` untouched)
corrects the energy without touching charges. **Fixes** the danger wherever P3 previously forced
charges on an otherwise-free system (compressed side, below the pass-1 fragment split): water-
probe label gap 8.9/14.3 -> **0.00** kcal/mol there; mean gap over 4 probe geometries 11.84 ->
5.06; react-mode collapse gone (needed 2 new repairs), matches repaired `flat` to 0.1 kcal/mol;
static full-grid rms basically unchanged (11.58/11.75 vs 11.52/11.65). **Does NOT fix, and
introduces a new risk, both stated plainly**: (1) the label gap AT THE REFERENCE MINIMA (the
pass-1 fragment split, Cl 2.64/F 1.92 A) remains 3.7/16.5 kcal/mol — for F2- as bad as `flat`'s
worst case, root cause is the discrete Phase-1 qa placement there, unreachable by any (r,x)
function; (2) a NEW topology-history dependence (up-vs-down scan differs up to 21.7 kcal/mol for
F2-, was 0.4 under `flat` — `flat`'s charge-forcing had been MASKING this pre-existing property of
free charges). Falsifiers all pass (1379/1379 frames bit-identical elsewhere, gradient exact,
ctest 66/69). **The implementing agent recommends switching P3's default mode to `harris` — NOT
adopted by the orchestrator without operator sign-off**, since this is exactly a static-danger-
vs-MD-history-danger trade-off the project reserves for an explicit decision, and the
improvement is near-zero at F2-'s actual equilibrium geometry. `harris` sits in the tree, opt-in,
next to `flat`/`frac`; nothing switches by default. Full detail: `P2P3_HARRIS_STATUS.md`,
`WORK_STATUS.md` §25.

**Package 26 (2026-09-24): the GENERAL fix, in PLAIN GFN-FF — `-gfnff.frag_charge_model
ensemble`.** Operator explicitly chose the broad (not rev-gfnff-scoped) fix after correcting the
orchestrator's own imprecise framing (the old "whole charge on fragment 0" rule is only right for
heterolytic-type separation). Opt-in, default `reference` = today's rule, proven bit-identical on
2462 GMTKN55 + 185 MOR41/S30L-CI structures. Two parts: chemistry-based carrier selection (fixes a
REAL index bug — WATER27 reaction MAD 58.6->21.4 alone) + a continuous window across the
topology-perception threshold (closes BOTH package-25 residuals: label gap ->0.00 everywhere,
up-vs-down history ->=<0.02 kcal/mol). **Honest cost**: the window only helps energetically when
paired with `harris` (P3's one-fragment-side fix) — s_max=1.2 + harris gives **8.50/8.00** kcal/mol
Cl2-/F2- rms (best result yet); plain GFN-FF ALONE gets WORSE with the window on (its own
one-fragment X2- is independently 80-220 kcal/mol too deep, a pre-existing defect this doesn't
fix). GMTKN55 WTMAD-2 94.09->91.50 but SIE4x4/BH76RC get worse (a separate pre-existing
ion-energetics inconsistency the old buggy rule was accidentally masking). MD through the window
needs dt<=0.05 fs. Falsifiers: full `ctest` (not just gfnff-labelled) same 13 pre-existing
failures +1 new test; gradient exact, 2 adversarial corruptions both caught. **Not adopted as any
default by the orchestrator** — genuinely a multi-part decision (adopt carrier fix? adopt window,
and at what s_max, for plain GFN-FF vs. rev-gfnff+harris only?) reserved for the operator. Full
detail: `FRAG_CHARGE_STATUS.md`, `WORK_STATUS.md` §26.

**Package 27 (2026-09-24): the fragment-charge default flip shipped and verified.** Operator
decision (AskUserQuestion): adopt chemistry-aware carrier selection as the new plain-GFN-FF
default, continuous window OFF by default (reserved as a documented recommendation for
`harris`+`s_max=1.2`, not a default). `frag_charge_model` reference->**ensemble**,
`frag_charge_s_max` 1.1->**1.0** are now the shipped defaults. **A real deployment bug found and
fixed along the way**: editing the PARAM macro defaults in `gfnff.h` alone had ZERO runtime
effect — no CLI path merges full registry defaults into the controller, so the matching
hardcoded fallbacks in `gfnff_method.cpp` had to be updated too (see Recurring traps below).
Verified to match package 26's own table exactly: GMTKN55 26/2462 structures move, WTMAD-2
92.595 (target 92.60); MOR41+S30L-CI (185, all neutral) bit-identical; full `ctest` (312 tests)
same 13 pre-existing failures, no new ones; `cli_gfnff_05` fixed to test explicitly rather than
relying on the (now-changed) default, still 9/9. Binary md5 `d9908823`. **Orchestrator's docs
folding (GFNFF_STATUS.md, REV_GFNFF_STAGE2.md, top-level CLAUDE.md) still pending.** Full detail:
`FRAG_CHARGE_STATUS.md` §17, `WORK_STATUS.md` §27.

**Package 28 (2026-09-24): n=2 scope assessment (Br2-/I2-/O2-/SN2-TS), no new compute spent.**
Mechanism generalizes structurally for free (fires correctly for Br2-/I2-/ClF-/BrCl-/HO-OH-, no
hardcoded element names). Br2- confirmed in the exact pre-fix broken state Cl2-/F2- were in
(well -105.6 vs true ~-25..-28 kcal/mol, same double-counting); needs TWO new well rows (order-1
AND half-order) not one. Native GFN2 confirmed unusable as a cheap sanity reference (same
long-range failure now shown for Cl2- too, not just F2-). **O2-/S2- are an architectural dead
end**, not a data gap — their extra electron sits in a pi* orbital the conserving-share
perception cannot see; needs a code change. A latent SN2-TS foot-gun noted (harmless today, would
bite if a C-X half-order row is ever added). **A NEW danger found, affecting Cl2- too**: with the
CURRENTLY RECOMMENDED settings, charges are asymmetric again just past the bond cutoff — a
distance window package 26 never tested — misleading a nearby water by up to 7.8 kcal/mol
(Cl2-)/6.4 (Br2-). NOT fixed. Proposed a ~63-job, <1h Br2- DLPNO-CCSD(T) campaign; **explicitly
did not launch it, asks the operator first**. Full detail: `X2_SCOPE_STATUS.md`, `WORK_STATUS.md`
§28.

**Package 29 (2026-09-24): stale-CN bug fixed (plain GFN-FF), a second related bug found.** Fix
A (applied): energy-only calls on a REUSED calculator used a stale cached CN in the Coulomb
chi(CN) term (gradient calls were always correct) — reproduces the on-record falsifier exactly
(Cl2- 1.77e-2 -> 2.41e-5 Eh/A); zero regression on GMTKN55(2462)/MOR41(95)/S30L-CI(90) single
points (bug invisible to per-structure benchmarks, only bites reused instances); real wins found
(native LBFGS optimizer: caffeine 5000-iter non-convergence -> 44 steps; H2O Hessian freq error
81->0.3 cm-1). Fix B (found while verifying A, NOT applied, kept as a reviewable patch): D4
pairwise C6 also never refreshes — drops remaining FD residuals to ~1e-11 but costs 2 more
`cli_curcumaopt_07` golden-value failures + one shifted optimization minimum (UPU23/2h, n=1) +
~20-40% energy-only runtime on large systems — a real trade-off, awaiting operator review.
**Important operational finding**: 3 Opus agents ran in parallel sharing ONE working tree/build
dir (an orchestrator coordination gap, no worktree isolation) — this agent caught a concrete
consequence (a concurrent agent's built binary silently diverged from the reverted source tree)
and flagged it directly to the other agents. Also: `/tmp` filled up completely mid-session;
anything measured in that window by any agent needs re-confirming. Full detail:
`STALE_CN_STATUS.md`, `WORK_STATUS.md` §29.

**Package 31 (2026-09-24): Br2- campaign done (n=2 -> n=3 proven), and a real fidelity-invariant
break found in the "best result yet" config.** 23 DLPNO-CCSD(T) + 40 r2SCAN-3c jobs, Br2- D_e=
28.45 kcal/mol @2.78 A; fitted rows give full-grid rms 88.6->9.55, matching Cl2-/F2-. Investigating
package 28's past-cutoff finding found the real cause: the continuous window's merged corner has
NO charge path between its two atoms past the bond cutoff, silently defaulting to the lower-
indexed atom — **breaking the kappa=0-equals-EEQ invariant by up to 108 kcal/mol** (this session's
most-repeated fidelity check). Fixed for 12/15 cases by new opt-in `rev_sqe_virtual_pairs`
(label gap ->0.00 everywhere it applies); 3 anionic-SN2-TS cases have a separate, pre-existing
leak, explicitly NOT fixed (would touch the shipped Cl2-/F2- fits). **Cost: package 26's reported
8.50/8.00 kcal/mol was ~1.5 kcal/mol better than it should've been** (riding on the bug); honest
number ~9.7/9.7 — corrected in `docs/REV_GFNFF_STAGE2.md`. Orchestrator independently confirmed
via a from-scratch clean rebuild (md5 `aa7d9cde`, matches the agent's own) that the source tree is
now consistent; full `ctest` **293/306**, exactly the 13 documented pre-existing failures (the
14th, `gfnff_sqe`'s stale `B2/3d` threshold from package 30's mu-fix, fixed directly by the
orchestrator — converted to informational, same treatment as package 23's precedent). All three
Sep-24 follow-up tasks (mu-cusp, stale-CN, n=2-scope) are now complete. Full detail:
`X2_SCOPE_STATUS.md` §8-18, `MU_CUSP_STATUS.md`, `STALE_CN_STATUS.md`, `WORK_STATUS.md` §§29-31.

**A concurrency lesson, resolved**: 3 Opus agents shared one working tree/build dir this session
with no isolation, causing real interference (a corrupted `build_rev`, a transiently-missing
external dependency file, agents attributing the same new test failure to each other). All
resolved by direct investigation + a clean rebuild; going forward use `isolation: "worktree"` for
concurrent source-editing agents (see Recurring traps below).

**Package 32 (2026-09-25): stale-CN Fix B applied (operator: "ja, übernehmen").** Clean rebuild
(md5 `999b8f90`); `cli_curcumaopt_07`'s 17 golden values regenerated from the corrected binary —
uncovered 2 of them had encoded a REAL ~100 kcal/mol optimizer misconvergence, not noise; now
20/20 PASS. **Full `ctest`: 294/306, 12 failures (was 13)** — this is now the current baseline
failure count for the whole suite, superseding "13" in earlier package notes above. Zero
regression (GMTKN55/MOR41/S30L-CI bit-identical, checked via a with/without-Fix-B binary diff).
Aside found, not caused by this fix: this WIP branch's vs-xtb MAD is higher than CLAUDE.md's
documented baseline (0.86/16.0 vs ~0.26/~0.00 GMTKN55/MOR41) — present with or without Fix B,
flagged in Known Issue #35, not root-caused. Docs folded: CLAUDE.md Known Issue #35,
`AIChangelog.md`, `docs/GFNFF_STATUS.md`; Known Issue #34's stale "8.50/8.00" text also corrected
in place. `rev_sqe_virtual_pairs` was already documented as the package-31 recommendation, no
further action needed for that half of the operator's decision. Full detail:
`STALE_CN_STATUS.md` §11, `WORK_STATUS.md` §32.

Full detail: `docs/REV_GFNFF_STAGE1.md` (package 15), `FABLE_REVIEW_3.md` (Q1-Q4),
`STAGE2_B2_STATUS.md` (package 18), `P2P3_STATUS.md` (package 23),
`P2P3_ALTERNATIVES_STATUS.md` (package 24), `P2P3_HARRIS_STATUS.md` (package 25),
`FRAG_CHARGE_STATUS.md` (packages 26-27), `WORK_STATUS.md` §§15-27, `KAPPA_FIT_STATUS.md` (all
five fit attempts), `docs/REV_GFNFF_STAGE2.md` (packages 24-27 not yet folded in),
`docs/GFNFF_STATUS.md` (package 26/27 belongs here too, plain-GFN-FF scope — not yet written).

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
| plain `-method gfnff` | bit-identical through package 25; **package 27 changed its default** (`frag_charge_model=ensemble`, `frag_charge_s_max=1.0` — carrier-selection fix only, window off) | package 27 |
| `-gfnff.frag_charge_model`/`_s_max` | `ensemble`/`1.0` (plain default); recommend `ensemble`+`1.2` alongside `harris` for rev-gfnff X2- (not itself a default) | package 27 |

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
| 14 | stage-2 resumed: q0 affine-shift fix + `revgfnff_fit.py` SQE wiring | q0 fix verified (fidelity/static-gradient bit-identical; a small pre-existing react-corner gradient gap newly exposed, documented, non-blocking); fit tooling validated end-to-end at p0, kappa_Z fit itself not yet run | `WORK_STATUS.md` §14, `docs/REV_GFNFF_STAGE2.md` |
| 15 | kappa_Z fit stalled (LM+NM) on CHB6; root-caused + fixed a stage-1 bug | bare alkali/alkaline-earth cations (cation-pi, e.g. Li+-benzene) mis-read as grossly over-coordinated (+4.97 Eh spurious); fixed, CHB6 MAD 1550.6->47.6 kcal/mol, no ctest regression | `WORK_STATUS.md` §15, `docs/REV_GFNFF_STAGE1.md` |
| 16 | 3rd fit attempt stalled on PX13; found it's neutral, not anionic SN2 | added a `BH76_anionic` virtual subset (the real 16-reaction anionic-SN2 target); PX13/full BH76 kept report-only | `WORK_STATUS.md` §16 |
| 17 | Fable review (Q1-Q4) + Layer-A scoring fix + 4th/5th fit attempts | class-E metric was broken 2 ways (94% loss = 4 fusion frames; Cl2-/F2- had near-zero kappa leverage); Layer-A fix gives real movement (loss -10.4%) but kappa_Cl still can't reach its target (range cap excludes its one react-mode leverage frame) — open question is now B1 vs B2 (model-side), not another fit | `WORK_STATUS.md` §17, `FABLE_REVIEW_3.md`, `KAPPA_FIT_STATUS.md` |
| 18 | B2 implemented (q0 by chemical potential + a 2nd kappa(b) form) | real model improvement, verified by 5 falsifiers; Cl2- curve target proven structurally unreachable by ANY kappa_Z (~104 of ~148 kcal/mol error is not in the charge model) — redirects the open question to stage 1/3a; one new gradient-cusp defect found, not fixed | `WORK_STATUS.md` §18, `STAGE2_B2_STATUS.md`, `docs/REV_GFNFF_STAGE2.md` (reconciled) |
| 19 | compressed Cl2-/F2- decomposed term-by-term; corrects §18 | 38 of 106 kcal/mol IS still Coulomb (qa frozen at Phase-1, kappa can't reach); bond term is electron-count-blind (wrong r0); r2SCAN-3c reference itself has a large SIE tail — the day's -41.5 calibration target may be inflated; mg3's Cl-Cl well-fit data is separately RKS-contaminated; 3 proposals, none shipped | `WORK_STATUS.md` §19, `CL2_COMPRESSED_STATUS.md` |
| 20 | P1 executed: clean Cl-Cl UKS data + refit, applied | Cl-Cl fit rms 9.53->1.07; neutral Cl2 D_e -70.1->-53.9 kcal/mol; anion moved +13..+16 (still far from reference, as expected); all regressions pass; P2+P3 next, to be planned as one decision | `WORK_STATUS.md` §20, `CL2_WELLFIT_P1_STATUS.md` |
| 21 | DLPNO-CCSD(T) reference for Cl2-/F2-: the r2SCAN-3c target itself was wrong | true D_e(Cl2-) -28.4 (was targeting -41.5, 32% too deep), D_e(F2-) -26.8 (was -49.5, 46% too deep); long-range SIE confirmed 2 orders of magnitude; GFN2 found worse than r2SCAN-3c for F2- at long range (new); STAGE2.md acceptance-3 corrected | `WORK_STATUS.md` §21, `CL2F2_CCSDT_STATUS.md` |
| 22 | quick re-eval of P1+B2 against the corrected DLPNO target | on the original 6-pt grid, actually BETTER than the old (wrong) target suggested (Cl2- rms 55 vs 63, F2- rms 14); on the full grid (incl. never-before-tested deep compression) markedly worse (Cl2- 87.5, F2- 60-65) — the honest baseline P2/P3 had to beat | `WORK_STATUS.md` §22 |
| 23 | P2+P3 implemented — double-counting resolved for the tested case | Cl2-/F2- full-grid rms 87.6/64.9 -> 11.5/11.7 kcal/mol at kappa=0, no fit; resonance moved entirely into the bond well; zero regression; react-mode breakage + a newly found stale-CN bug (unrelated) both flagged, not fixed; one outdated ctest assertion retired, ctest back to 66/69 | `WORK_STATUS.md` §23, `P2P3_STATUS.md` |
| 24 | operator-flagged danger in P2+P3's design analysed, alternatives tested | react collapse (-107..-203 kcal/mol) fixed via a corner-consistency repair (default on); broken-symmetry charges mislead a 3rd-body probe by 9-15 kcal/mol, NOT fixed (2 alternatives tried, both rejected); P3 stays opt-in, recommended off for any interacting X2- | `WORK_STATUS.md` §24, `P2P3_ALTERNATIVES_STATUS.md` |
| 25 | built package 24's section-6 sketch: `rev_excess_mode harris` (free charges + non-self-consistent x*g(r) correction) | fixes the label-gap danger on the compressed side (8.9/14.3 -> 0.00 kcal/mol); does NOT fix it at the reference minima (3.7/16.5 remains); introduces a NEW topology-history dependence (up to 21.7 kcal/mol) `flat` didn't have; agent recommends switching the default, orchestrator did NOT adopt it pending operator decision | `WORK_STATUS.md` §25, `P2P3_HARRIS_STATUS.md` |
| 26 | GFN-FF-wide fix: chemistry-aware + continuous fragment-charge placement (`frag_charge_model ensemble`) | closes BOTH package-25 residuals everywhere (label gap ->0.00, history ->=<0.02); fixes a real index bug (WATER27 MAD 58.6->21.4); costs: plain-GFN-FF-alone gets WORSE in the window (pre-existing over-binding not fixed here), SIE4x4/BH76RC worsen (unmasked pre-existing defect), MD needs dt<=0.05 fs; harris+s_max1.2 gives the best Cl2-/F2- result yet (8.50/8.00); opt-in, no default changed | `WORK_STATUS.md` §26, `FRAG_CHARGE_STATUS.md` |
| 27 | operator decision shipped: carrier-selection fix -> plain-GFN-FF default, window stays opt-in | `frag_charge_model`/`_s_max` defaults ensemble/1.0; found+fixed a PARAM-macro-default-has-no-runtime-effect deployment bug along the way; verified exactly against package 26's table (WTMAD-2 92.595, 26/2462 moved, MOR41/S30L-CI bit-identical, full ctest clean) | `WORK_STATUS.md` §27, `FRAG_CHARGE_STATUS.md` §17 |

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
- **Concurrent agents editing the SAME source tree/build dir can silently diverge from each
  other** — a rebuild by agent B can change the binary an already-running agent A is measuring
  against, without A knowing (found package 29: agent A's own built binary md5 no longer matched
  what the shared source tree contained after B reverted a change). Launch concurrent agents that
  edit shared C++ source with `isolation: "worktree"`, not a bare shared tree; if that wasn't done,
  message the other agents directly and have each do a final clean rebuild + re-verify right
  before finalizing, not trust an earlier build.
- **Editing a PARAM macro's default value in `gfnff.h` alone does NOT change runtime behaviour.**
  `GFNFF::GFNFF(const json&)` reads each setting via `m_parameters.value(key, HARDCODED_FALLBACK)`
  in `gfnff_method.cpp`, and no normal CLI path merges the full ParameterRegistry defaults into
  the controller (only `-export_run`'s dump does). Any future default change needs the matching
  hardcoded fallback updated too, or the change silently does nothing (package 27 caught this
  only because it verified against known numbers before concluding success).

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
