# FABLE_REVIEW_3 — stage-2 (SQE kappa_Z) method review (2026-09-22)
Sections written: 4/4 (complete)

AI-generated review, read-only on the repository. Own measurements: `build_rev/curcuma`
(md5 `a7c8c85165f0912c6bc7c160da036a97`, the binary the four fit attempts used), `-threads 1`,
fresh directory per run, no `*.topo.json` reuse. Harness and raw tables in
`scratchpad/kprobe/` (`probe.json`, `assoc/`, `static/`). Labels: **MEASURED** (own run),
**INFERRED**, **PROPOSED**.

---

## Summary (for the operator, read this first)

1. **Q1 — Yes, the class-E metric measures the wrong thing, but not for the reason the
   orchestrator guessed.** 94.3 % of class E's sum of squares is four frames of ONE system
   (`ahb21_21_stretch`, r_FH >= 2.0 A) in which the bridging hydrogen sits 0.21–0.49 A from the
   oxygen — nuclear-repulsion geometries that r2SCAN-3c itself puts 540–5100 kcal/mol above the
   minimum. The fit's loss was 98 % "how well does a force field reproduce two nuclei fused
   together". That is what attempt 4's kappa_C = 0.41 was earned on (it shaved 17 kcal off a
   9605 kcal residual). Separately, and worse: **the Cl2-/F2- channel — kappa_Cl/kappa_F's
   only design target — has exactly zero kappa leverage on 18 of 20 frames in every mode the fit
   can run** (react either direction, static), because the q0 rule spreads the integer charge
   uniformly over a perceived fragment and a symmetric uniform q0 admits no charge flow for
   kappa to penalise. The design doc's Cl2- table was a single-frame regime (static perception
   with nfrag = 2 at exactly 2.73 A) that the fit never enters. Under the current q0 rule the
   model's Cl2- well is **-150 kcal/mol at 2.05 A at every kappa** against the reference's +40
   relative to r_eq; acceptance criterion 3 ("-41.5 at r_eq") is satisfiable only as a
   one-point coincidence on a curve with the wrong shape. Concrete replacement scoring in Q1.
2. **Q2 — The 2.36e-4 Eh/A react-corner gap is not demonstrably an SQE defect and can be
   deferred.** Read line by line, the SQE gradient is variationally complete (frozen q0, one
   bond-order switch for solve and kernel, envelope on p, blend term carries `sqe_hardness`);
   the "Term 1b assumes fragment-Lagrange stationarity" hypothesis does not hold up (Term 1b is
   dE/dCN at fixed q and is the same expression under either stationarity condition). A CLI
   reproduction of case 2c gives a gap of 4.63e-5 on the incoming H that is **bit-for-bit
   identical for `eeq` and for kappa = 0, 0.25, 0.5, 1.0, 2.0** — a pre-existing react-corner
   residual, kappa-independent to < 1e-7 — while the entire SQE-dependent force on that atom
   is only ~2e-4 Eh/A. The test's 2.36e-4 could not be reproduced and most plausibly depends
   on the order in which its FD sweep visits the atoms (30 react scans happen before it
   reaches the incoming H). Three concrete checks for an Opus agent and the explicit
   "stop deferring when" criteria are in Q2.
3. **Q3 — No. SQE stays opt-in; do not promote any current kappa_Z.** No calibrated
   kappa_Z exists (every fitted value was earned on the fusion frames), the two headline
   targets are unmet, kappa > 0 moves *neutral* proton-transfer barriers by tens of kcal/mol
   at kappa ~ 1 and the class-D guard cannot see that, and the `mg3` precedent (measured
   class-A rms 22 -> 13, 11 700 paired trajectories with no tail cost) has no counterpart here.
   The one real gain, BH76 MAD 48 -> 39 at the untouched p0, comes from kappa_Cl = 0.85 acting
   through the one-frame perception regime of Q1.2 and is a reason to keep the option, not to
   default it. A five-point gate for flipping the default is in Q3.
4. **Q4 — Real gap, wire it now as report-only, arm it as a penalty when Q3's gate is being
   assessed.** S66 is already on disk inside the GMTKN55 checkout and loadable through the
   existing `barriers` mechanism with zero new code (`load_reactions(["S66"])`); the conformer
   set is `CHEM_CLASS["conformers"]` in `gmtkn55_reactions.py`; the roadmap's validation
   matrix already fixes the absolute thresholds (conformers <= 1.6, NCI <= 8, charged NCI
   <= 38.9). A `guards` config entry reusing the class-D `excess -> sqrt(1e3)*excess` pattern,
   evaluated once per accepted LM step rather than in every Jacobian probe, is a ~40-line
   change; sketch in Q4. It matters precisely because Q1.3 showed the neutral-molecule side
   effect the class-D guard is blind to.

**What the operator should decide next, in order:** (1) approve Layer A of Q1.4 (data
hygiene, no model change) and a fifth fit on the cleaned set — cheap, and it will move
kappa_H/O/C/F for the first time on real chemistry; (2) decide between B1 (fit-side, thin
signal via the perception accident) and B2 (q0 localised by mu + steeper kappa(b), a design
change) for the Cl2-/F2- channel — this is the actual stage-2 design question and it is not
a fit-campaign decision; (3) leave SQE opt-in until the Q3 gate is met; (4) let a Sonnet
agent wire the S66/conformer guard report-only in the same pass as Layer A.

---

## Q1. Is the class-E loss measuring the wrong thing? — YES, on two independent counts (MEASURED)

**Finding, stated plainly.** The hypothesis "the aggregate RMS is dominated by a flat tail kappa
cannot influence" is *refuted in its mechanism and confirmed in its conclusion*. The tail of
the diatomics is indeed kappa-invariant, but it is only ~3 % of the loss. What actually
dominates is a data defect (Q1.1), and what actually blocks kappa_Cl/kappa_F is a q0-rule
property of the model, not the sampling (Q1.2). Attempt 4's "kappa moved, targets did not" is
fully explained by the two together (Q1.3). The recommendation is in Q1.4.

### Q1.1 The loss is 94 % one system's four unphysical frames (MEASURED)

Per-point residuals recomputed from the attempt-4 run directory
(`fit_work/stage2_kappa_run4/*.jsonl` = the last evaluation, i.e. the fitted kappa) against
`ref/E/*/energies.json`, using the fit's own definition (dE relative to frame 0 after the react
ordering shortest -> longest; `revgfnff_fit.py:573-585`). This reproduces the reported
`rms_dE = 533.85` to the digit, so the decomposition is the fit's own number, not a proxy:

| system | n | sum of squares | share of class E | rms of the system |
|---|---:|---:|---:|---:|
| `ahb21_21_stretch` | 20 | 3.14e7 | **94.3 %** | 1253.6 |
| `hf2_transit` | 11 | 9.15e5 | 2.7 % | 288.4 |
| `f2m_F-F-` | 20 | 5.77e5 | 1.7 % | 169.8 |
| `cl2m_Cl-Cl-` | 20 | 3.86e5 | 1.2 % | 139.0 |
| `nh3_2_transit` | 11 | 2.49e4 | 0.1 % | 47.6 |
| `h2o2_transit` / `nh4_nh3_pt` / `fch3f_umbrella` | 11/13/11 | < 1e4 each | ~0 | 27.9 / 9.0 / 7.7 |

Inside `ahb21_21_stretch` the residual is 1–30 kcal/mol for r_FH <= 1.6 A and then 41 / 356 /
**2660 / 4520 / 1928** / 234 / 135 for r_FH = 1.81 … 3.51 A. The scan holds O…F fixed at
2.460 A and moves the H along the F->O line, so r_FH >= 2.0 A puts the H *at or through* the
oxygen: min O–H = 0.488, 0.280, **0.211**, 0.369, 0.595 A for r_FH = 2.01 … 3.01 A
(`points.xyz`, own geometry read-back). r2SCAN-3c itself puts those frames 541 / 2775 / 5087 /
1569 / 456 kcal/mol above the minimum. These are not proton-transfer geometries; they are
nuclear-fusion geometries produced by a scan that ran ~1 A past the point where the coordinate
stops meaning anything. No force-field parameter can or should fit them, and in a
sum-of-squares they outweigh the 13 physical frames of the same system by 400:1.

Since the class-E chunk enters the loss as `sum (m-r)^2 / n_E` with `weight_E = 1` and the
barrier chunks as `sum (m-r)^2 / n_R`, class E = 533.85^2 = 285 000 of the 289 328 total loss:
**98.5 % of the objective was these four frames.**

### Q1.2 Kappa has zero leverage on the diatomic anions in every mode the fit can run (MEASURED)

Direct per-point sensitivity, `-method revgfnff`, react mode exactly as `run_system()` issues
it, five charge settings: `eeq`, and `sqe` with all six kappa_Z = 0 / the attempt-4 vector /
1.0 / 2.0 (the bound). "Leverage" = max − min over the four sqe settings, kcal/mol, per point:

| system | points with leverage = 0.0 | where kappa acts | leverage there |
|---|---:|---|---:|
| `cl2m_Cl-Cl-` (react, short->long, as fitted) | **18 / 20** | r = 3.55, 3.82 A only (the blend window of the breaking transition) | 93, 111 |
| `f2m_F-F-` (same) | **19 / 20** | r = 2.50 A only | 209 |
| `cl2m` in the **association** direction (long->short) | **20 / 20** | nowhere | 0 |
| `f2m` in the association direction | **20 / 20** | nowhere | 0 |
| `cl2m` **static**, one perception per frame | 7 / 8 tested | r = 2.728 A only (perceived nfrag = 2, pair still present) | 103 (k0 -125.8 -> k2 -40.9) |

The reason is the q0 rule, not the sampling. `docs/REV_GFNFF_STAGE2.md` ("Reference charges
q0"): the integer fragment charge is *spread uniformly over the fragment's atoms*. For a
symmetric homonuclear anion perceived as one fragment that is q0 = (−0.5, −0.5), which is
already the EEQ minimum; the SQE increments p are then zero for every kappa, and the hardness
term `1/2 kappa0/b p^2` is identically zero. kappa can only penalise charge that *flows away
from q0*, and there is nothing to flow. This holds:

- in react short->long (frame 0 at 2.05 A is one fragment; the frozen q0 is uniform; the
  rounding rule at the breaking transition produces (−1, 0) only inside the blend window, after
  which the pair is gone);
- in react long->short (frame 0 is (−1, 0), but the formation transition's completion rebuilds
  the slot and re-derives q0 from the merged perception — the probe shows E(2.05 A) = −156.0
  kcal/mol relative to separately computed Cl + Cl- for eeq, k0, k0.85, k1 and k2 alike);
- in static mode for every frame the perception calls one fragment (r <= 2.59 A: −150 / −136
  / −128 at every kappa) and for every frame with no pair (r >= 2.86 A: +7.9 / +4.9 / +1.4 at
  every kappa).

The design doc's five-row Cl2- table (`STAGE2.md:161-167`) therefore contains exactly ONE
kappa-sensitive row, the r = 2.73 A one, and it is sensitive only because plain perception
carries nfrag = 2 from the q-loop's first pass there (Known Issue #17) while the rev pair
(b = 0.6) still exists. That is a one-frame regime. "kappa_Cl ~ 0.85 hits −41.5" was a
calibration on one point of a curve whose minimum is −150 at 2.05 A regardless of kappa
(reference: +39.8 above r_eq there). **kappa_Cl was never actually tested by any of the four
fits**: it entered at 0.85, had no gradient from class E, and drifted to 0.837 on the
Cl-containing barrier reactions.

Note also what the reference curves are past ~1.2 r_eq: r2SCAN-3c UKS Cl2- at 9.55 A is only
3.7 kcal/mol above the 2.73 A minimum, with Hirshfeld charges (−0.4997, −0.4997); F2- at
6.72 A *is* the curve's minimum (0.0), 0.8 below the bonded 2.0–2.1 A region. That is the
DFT delocalisation error (symmetric dissociation to two half-charged atoms), not chemistry.
The −41.5 the design targets was obtained, correctly, with *separately computed* Cl and Cl-
(`docs/REV_GFNFF_TODO.md:58`), not from the curve tail. So the 11+10 tail points of the two
diatomics must never be fitted against — they would teach the model that Cl2- does not
dissociate.

### Q1.3 Why attempt 4 moved kappa_C to 0.41 without touching either target (INFERRED from Q1.1/Q1.2)

- The LM's forward-difference Jacobian was dominated by the four fusion frames: at r_FH = 2.51 A
  the model residual is 4520 kcal/mol and kappa moves it (k0 9621.9, kfit 9604.8, k1 9247.6 —
  leverage 374). Any kappa direction that lowers those by 0.2 % is worth more to the loss than
  fixing every chemistry target in the set. `ahb21_21_stretch` contains H, C, O, F; kappa_C
  is the one of those with no competing barrier signal, so it is the one that moved.
- Class E's headline systems (`cl2m`, `f2m`) contributed no gradient at all (Q1.2).
- `BH76_anionic` and `CHB6` moved by ±1 kcal/mol because they sit at 1.5 % of the loss.
- The report-only PX13 "staying flat" was read as confirming that neutral systems are
  kappa-insensitive. **That reading is wrong.** The neutral transits in class E have large
  leverage — `h2o2_transit` up to 79 kcal/mol at the TS frame (eeq +76.9, k1 +16.2, k2 −2.4;
  ref +43.9), `nh3_2_transit` up to 28, `hf2_transit` up to 38 — through kappa_H/O/N/F, which
  the fit never stepped. Package 16's premise ("no net charge anywhere, so a charge-hardness
  parameter cannot move it") does not hold for SQE: the hardness damps *every* bond-charge
  increment, in neutral molecules too, and at kappa_Z ~ 1 it changes neutral proton-transfer
  barriers by tens of kcal/mol. Demoting PX13 to report-only is still defensible (its residual
  is a stage-1/3 bond-term problem — `hf2_transit`'s TS frame is +877 kcal/mol at eeq and +840
  at k1, i.e. the bulk of it is not charge), but the stated reason is false and must not be
  reused as a rule.

### Q1.4 Replacement scoring — concrete, in two layers (PROPOSED)

**Layer A — data hygiene, config/script-level, no model change, a Sonnet agent can do it in one
pass.** All of this goes into `scripts/revgfnff_fit.py`'s class-E handling and
`stage2_kappa_config.json`; none of it touches `gfnff_method.cpp`.

1. *Reference-energy cap per system.* New config key `"max_ref_dE_kcal": 100.0` (dataset-level,
   applies to `classes`): drop every point whose reference energy lies more than that above the
   system's reference minimum, before `class_pts` is built. Removes exactly the 7 frames of
   `ahb21_21_stretch` with r_FH >= 1.81 A (145.6 … 5129 above the minimum) and nothing else in
   class E (the largest surviving value is `nh3_2_transit`'s 57.9). Rationale: 100 kcal/mol is
   about one bond energy; a force field is not asked to reproduce the wall beyond that.
   Belt-and-braces alternative that also fixes the data on disk: regenerate the AHB21/21 scan
   with r_FH in [0.75, 1.70] A (`scripts/revgfnff_ref.py`, ~11 ORCA points, ~2 min) — but the cap
   should exist anyway, so the next scan that overshoots is caught automatically.
2. *Range cap for the two diatomic anions.* Per-system `"max_scan": {"cl2m_Cl-Cl-": 3.3,
   "f2m_F-F-": 2.5}` (i.e. r <= ~1.2 r_eq,ref): beyond it the r2SCAN-3c curve is the
   delocalisation-error plateau (Q1.2, last paragraph). Keep the excluded tail as a
   *reference-free shape guard*: the model's own E(r) must be monotonically non-decreasing for
   r > r_eq,ref (that is acceptance 3's "stays monotonic beyond it", and it needs no reference).
   Implement as one extra residual per system, `sqrt(1e2) * sum_k max(0, E_k - E_{k+1})` over
   the excluded frames, same pattern as the class-D guard's `sqrt(1e3) * excess`.
3. *Anchor.* For the diatomics, replace "dE relative to frame 0" by "E(r) − [E(X) + E(X-)]" on
   both sides — the model's fragment energies from two one-atom single points (already trivially
   available: Cl + Cl- = −0.96970056 Eh, F + F- = −0.88314547 Eh with this binary), the
   reference's from the separately computed r2SCAN-3c atoms that produced the −41.5 (add them to
   `ref/E/<system>/energies.json` as `"fragment_energies_eh"`, or to a new `ref/L` entry; the
   loader gains one optional field). This is the quantity acceptance 3 is written in, and it
   removes the compressed frame 0 as the zero. For every other class-E system keep frame 0 but
   also add an explicit *TS/midpoint residual*: `(m − r)` at the frame of maximum reference dE
   (the barrier), weighted `sqrt(n_sys)` so that one point counts like the whole curve. That is
   the point-wise target the design doc asks for, without giving up the curve.
4. *Per-system normalisation.* `sqrt(w_E / (n_systems * n_sys))` instead of `sqrt(w_E / n_E)`
   (`revgfnff_fit.py:622-623`), so a 20-frame scan does not outweigh an 11-frame one by 1.8x.
5. *Weight rebalancing after 1–4.* With the fusion frames gone class E's rms drops from 534 to
   roughly the diatomic-tail level (~100 before item 2, ~30 after); the barrier chunks (MAD 14–75)
   then carry comparable weight without any explicit reweighting. Do NOT raise `weight_E` to
   compensate — the whole point is that class E stops dominating.
6. *Promote the neutral transits to fitted, with a cap, once `hf2_transit`'s +840 kcal/mol TS
   frame is understood* (it is not charge: eeq and k2 differ by 37 there). Until then leave
   `h2o2_transit`/`nh3_2_transit`/`hf2_transit` at `weight_E = 0` (report-only), exactly as PX13,
   and for the same *correct* reason: their residual is dominated by a non-charge term and would
   make kappa_H/F absorb it. `nh4_nh3_pt` stays fitted (leverage 23–52, ref flat within 3, k1
   already overshoots to −29 — it will pin kappa_H/kappa_N at a small value, likely 0.1–0.2).

**Layer B — the model, an operator decision, NOT a fit-config change.** Layer A makes the
loss honest but cannot give kappa_Cl/kappa_F a lever on the bonded region of Cl2-/F2- (Q1.2:
zero leverage at r <= 2.59 A in every mode). Two ways out, in increasing invasiveness:

- *B1, fit-side only:* score the diatomics in **static** mode with one perception per frame
  (`-batch_reuse_topology false` for those two systems), which recovers the design doc's regime
  on the frames where plain perception gives nfrag = 2 while the rev pair still exists — with this
  binary that is the r = 2.73 A frame alone for Cl2-; for F2- it would have to be measured
  (frames 2.0–2.5 A). One or two frames per system is a thin but honest signal for the
  −41.5 / −(F2- analogue) targets, and the r < r_eq frames become a *guard* ("must not get
  deeper than the reference by more than X") rather than a fit target. Cheap; recommended as
  the immediate step. Its weakness is that the signal exists only because of a perception
  accident (Known Issue #17), which a later perception fix could remove.
- *B2, model-side:* change `revSqeQ0Fragments` (the static initialisation rule) so that a
  charged fragment's integer charge is placed on the fragment atom with the lowest EEQ chemical
  potential mu (`EEQSolver::calculateChemicalPotential`, already there for the corner rule)
  instead of uniformly. kappa0 -> 0 still reproduces EEQ exactly (the design's own argument:
  placement inside a connected fragment is immaterial at kappa = 0), so the fidelity acceptance
  is untouched, but for kappa > 0 the whole bonded region becomes kappa-governed: the delocalised
  (−0.5, −0.5) state costs `1/2 kappa0/b p^2` with p = 0.5, i.e. `kappa0/8` per pair at b ~ 1,
  and the Cl2- well would be `-150 + 627.5 * kappa_Cl/8` kcal/mol at 2.05 A — kappa_Cl ~ 1.3
  would land the compressed well near the reference. The price: every charged fragment's
  intramolecular delocalisation (carboxylate resonance, the two O of nitro anions) is also
  damped by the same kappa, since kappa(b) = kappa0/b barely distinguishes a b = 0.97 bond from
  a b = 0.6 one (factor 1.6). If B2 is chosen, the b-dependence of kappa must be made steeper at
  the same time (e.g. `kappa0 * (1 - b)/b` or `kappa0 / b^n`, n >= 3), otherwise AHB21/CHB6/IL16
  will pay for Cl2-. That is a design change with its own falsifiers and is out of scope for a
  fit campaign; it is the actual open design question of stage 2.

**What would refute Q1.4's premise** (answering the CLAUDE.md gegenfrage): if after Layer A
the class-E rms still does not respond to kappa in the 0–2 bound. The per-point table in
`scratchpad/kprobe/probe.json` says it will — `ahb21_21_stretch`'s 13 physical frames have
leverage 11–116 kcal/mol with the reference strictly inside the swept range at every frame
(e.g. r = 1.20 A: ref −37.5, k0 +1.1, kfit −16.7, k1 −95.2), and `nh4_nh3_pt` 23–52 — so a fit
on the cleaned set will move kappa_H/O/C/F for certain. Whether it moves them to *chemically
right* values is what Q3 is about.

Sanity note on my own numbers: the react-mode react-direction dependence is real hysteresis
(E(2.05 A) relative to fragments is −150 in dissociation order, −156 in association order,
−149.7 static); it is a stage-1b property (frozen topology corners), not a measurement artefact,
and it means class-E dE values depend on the scan direction at the 5 kcal/mol level. Any
scoring scheme must fix the direction per system and say so in the config.

---

## Q2. The react-corner SQE gradient gap (~2.36e-4 Eh/A) — defer; no SQE-specific defect is demonstrated

**Finding, stated plainly.** By reading, the SQE gradient is complete; by measurement, the
residual at the case-2c geometry is kappa-independent and equal to the plain `eeq` react-corner
residual, and the SQE terms are too small on that atom to hide a 2.4e-4 error. The number in the
test is real but most likely a property of the react scan's state during the FD sweep, not of
the SQE formula. It does not block kappa > 0 MD or optimisation at the accuracy the model has
anyway; the criteria for revisiting are at the end of this section.

### Q2.1 What the code actually does (INFERRED from reading; file:line for the follow-up agent)

The energy is `E(p; r) = E_EEQ(q0 + B p; r) + 1/2 sum kappa0_ij / b_ij(r) p_ij^2`, minimised
over p without constraint. Its exact derivative is the envelope expression
`dE/dr = dE_EEQ/dr|_q + sum_ij (-1/2 p^2 kappa0/b^2) db/dr`, and every piece of it is present:

| piece | where | status |
|---|---|---|
| q0 frozen at corner creation only | `captureCornerEEQ(corner_generation=true)` at `gfnff_method.cpp:2649` (transition begin) and `:12973` (corner generation); revert keeps the old q0 explicitly (`:13093-13096`); the slot reuses `m_rev_corner_eeq.back().q0` while a transition is in flight (`revSlotCorner`, `:12848`) | correct — q0 is a constant of the FD, so no dE/dq0 term is owed |
| one bond-order switch for solve and kernel | solve: `RevGFNFF::bondOrder(r, rv.R2(i,j), rv.bo2_width)` at the CURRENT geometry (`:12878`); kernel: `revOrder()` = the same call (`ff_workspace.h:799-801`, used at `ff_workspace_gfnff.cpp:2901`) | consistent — p is optimal for the b the kernel differentiates |
| E_EEQ at fixed q: Coulomb kernel (Term 1) + chi(CN) chain (Term 1b) | `calcCoulomb` uses `m_eeq_charges` = the corner's SQE q (`ff_workspace_gfnff.cpp:1244/1290`); `postProcess` builds `qtmp = q cnf/(2 sqrt cn)` from the same vector (`ff_workspace.cpp:577-590`) | correct for SQE as written |
| hardness kernel `-1/2 p^2 kappa0/b^2 db/dr` | `calcSqeHardness` (`ff_workspace_gfnff.cpp:2884-2916`) | correct; clamped pairs contribute no force, and the solver already forced p = 0 there |
| blend: `sum_c w_c grad E_c + sum_t (sum_c +-w E_c) ds/dw dw/dr` | `ff_workspace.cpp:849-873`; `E_c` = `calculateSingle()` total, which includes `sqe_hardness` (`ff_workspace.h:103`, `addScaledComponents:800`) | complete |

On the hypothesis recorded in `test_gfnff_sqe.cpp:226-232` and `STAGE2.md:199-202` ("Term 1b
may assume the standard fragment-Lagrange-multiplier EEQ stationarity structure"): that is not
how the term works. Term 1b is `dE/dCN_i|_q = -q_i dchi_i/dCN_i`, an explicit derivative at
fixed charges; it is the same expression whichever stationarity condition produced q. The
stationarity structure only matters for the *implicit* term `sum_i dE/dq_i dq_i/dr`, which the
gradient omits — and it is allowed to omit it in BOTH models: for EEQ because `dE/dq_i = lambda_f`
and the fragment sums are constant, for SQE because `dE/dp = 0` and q0 is constant. So the
kappa = 0 coincidence the note cites is not evidence for that mechanism.

### Q2.2 The residual at the case-2c geometry is kappa-independent (MEASURED)

Case 2c reproduced through the batch CLI with the test's settings (`react_check_every 1`,
`react_valence_cap false`, `react_refractory_scans 0`, `react_exchange_scans 0`, walk
2.2 -> 1.8 -> 1.5 -> 1.3 -> 1.25 A, then +-1e-5 A displacements as further frames, all in one
persistent calculator exactly like the fit's `run_system()`; `scratchpad/q2/`):

| charge model | gap on H5 x / y / z (Eh/A) | gap on C x | gap on H1 x | FD of Coulomb on H5 x | FD of SqeHardness on H5 x |
|---|---:|---:|---:|---:|---:|
| `eeq` | +4.63e-5 (all three) | -4.64e-5 | -3.5e-6 | +1.2e-4 | — |
| `sqe`, kappa 0 | +4.63e-5 | -4.64e-5 | -3.5e-6 | +1.2e-4 | 0 |
| `sqe`, kappa 0.25 / 0.5 / 1.0 / 2.0 | +4.63e-5 (each) | -4.64e-5 | -3.5e-6 | +1.4e-4 / 1.6e-4 / 1.9e-4 / 2.1e-4 | +3e-5 / 3e-5 / 3e-5 / 2e-5 |

Three things follow. (a) The gap is the same to three digits with the charge model off, and
does not change by more than 1e-7 across the whole kappa range — it is the pre-existing
react-corner residual of stage 1b (the same class as `test_gfnff_rev_fd`'s 1.0e-4 on "forming
H2 1.7->1.05->0.95"), not something the hardness term or the SQE charges add. (b) The *entire*
force the SQE-dependent terms exert on the incoming H is ~2e-4 Eh/A (Coulomb + SqeHardness FD);
a defect confined to those terms could not produce a kappa-independent 2.4e-4 gap. (c) I could
not reproduce 2.36e-4. The one material difference to the test is the FD *order*: `fdResidual`
sweeps atoms 0..5 x/y/z in sequence, so by the time it reaches the incoming H (atom 5, last) the
calculator has done 30 more react scans with `check_every 1` and `refractory 0`; my run displaced
H5 first. That the residual was "h-independent from 1e-3 to 1e-7" is consistent with a
state-dependent offset in the react bookkeeping just as much as with a missing analytic term.
INFERRED, not verified — it is the first thing the follow-up agent should test (below).

### Q2.3 Severity, for MD and optimisation (INFERRED from recorded numbers)

- Relative size: 2.36e-4 on |g| = 0.31 Eh/A is 0.08 %. The model's own recorded react-mode
  residuals at comparable geometries: 1.1e-4 (CH4+H 1.2 A), 4.6e-4 (1.6 A), 3.3e-3 (N2H2
  stretched) — all in `test_gfnff_rev_fd` — and up to 1e-1 Eh/A for the H-H dynamic-r0 CN chain
  (`ff_methods/CLAUDE.md`, "Recorded pre-existing gradient residual"). The SQE number is at or
  below the smallest of these.
- MD consequence: a 2.4e-4 Eh/A non-conservative force on one atom inside a blend window is
  ~0.13 kJ/mol per Angstrom travelled; the stage-1 tail already tolerates 5–15 kJ/mol well-join
  events (`STAGE2.md:236-238`). It cannot be seen against that.
- Optimisation: a uniform 0.08 % gradient error does not move a stationary point measurably and
  only changes when a gradient-norm criterion fires.

### Q2.4 Recommendation: defer, with three cheap checks (PROPOSED, for an Opus agent, ~1 h)

1. In `test_gfnff_sqe.cpp` case 2c, run `fdResidual` with the atom order permuted (incoming H
   first, then reversed) and print the residual per ordering. If it changes with the ordering,
   record the number as react-scan state and close the "SQE gradient gap" as such. This is a
   test-only change.
2. Add an envelope audit behind an env gate (pattern: `CURCUMA_MP_GRAD_AUDIT` in gfn2):
   `CURCUMA_SQE_GRAD_AUDIT=1` makes the corner-prepare callback (`installCornerPrepare`,
   `gfnff_method.cpp:12987`) reuse the p vector of the previous call instead of re-solving.
   FD with frozen p must then agree with the analytic gradient to ~1e-8 (only explicit r
   terms remain); if it does and the re-solved FD does not, the defect is in the solve's
   variationality (convergence of the 1e-12-ridge + 2-step refinement path, `eeq_solver.cpp`
   `calculateSplitCharges`) and the fix is there; if both agree, there is no SQE gap at all.
3. Put the kappa scan of Q2.2 into the test (kappa 0 vs 0.5 vs 2.0 on the same walk) and assert
   the gap is kappa-independent to 1e-6. That is the regression test that would catch a real
   SQE gradient defect in the future, which the current absolute tolerance cannot.

**When this stops being deferrable**: (a) if check 2 shows a kappa-dependent gap; or (b) before
stage 2's react-MD acceptance (item 4) is re-run with kappa > 0 — measure the NVE drift slope
over >= 10 ps (Known Issue #32 protocol) for the `sqe` arm against the `eeq` arm on the same
seeds; a slope difference is the practical criterion that a gradient inconsistency matters.
Until one of those fires, kappa > 0 MD and `-opt` are no worse than what stage 1b already
ships.

---

## Q3. Should any kappa_Z become the default now? — No. SQE stays opt-in.

**Recommendation.** Keep `rev_charge_model = eeq` as the default. Do not promote the attempt-4
vector, the p0 vector (kappa_Cl = 0.85, rest 0), or anything a fit on the current class-E
scoring produces. Document `-gfnff.rev_charge_model sqe -gfnff.rev_sqe_kappa_Cl 0.85` as the
*experimental* anionic-SN2 setting with its one measured effect (BH76 MAD 48.1 -> 39.4 at p0,
`KAPPA_FIT_STATUS.md`) and its known limits (Q1.2).

**Reasoning, in order of weight.**

1. *There is no calibrated kappa_Z.* Attempt 4's kappa_C = 0.41 was earned on four
   nuclear-fusion frames (Q1.1/Q1.3); kappa_Cl never received a gradient from its design target
   (Q1.2); kappa_H/O stayed at 0 because nothing in the loss pulled them. Promoting that vector
   would ship a parameter whose value has no chemical provenance.
2. *The design's own acceptance is unmet on both headline targets.* Acceptance 3 (Cl2- curve)
   holds at one frame and is wrong in shape everywhere else (Q1.2: −150 kcal/mol at 2.05 A at
   every kappa). BH76_anionic got slightly worse by MAD in attempt 4 (69.2 -> 70.1). Charged
   NCI (AHB21/CHB6/IL16, the roadmap's decision-#6 headline) moved by −1 to +1 %.
3. *kappa > 0 is a global change to the electrostatics, not a charged-species patch.* At
   kappa_Z ~ 1 the neutral `h2o2_transit` TS moves by 60 kcal/mol, `nh3_2_transit` by 24
   (Q1.3). GFN-FF's chi, J, alpha, dxi and every charge-dependent bonded factor (fqq) were fitted
   with EEQ charges; SQE at fitted kappa changes every molecule's charges. That is a
   re-parametrisation and needs the full validation matrix (`REV_GFNFF_ROADMAP.md`,
   "Validation matrix"), not one guard.
4. *The class-D guard is blind to exactly this.* It is rms(dE) over 1000/2000 K snapshots of
   ten small neutral molecules with a 1.10 factor on 11.9 kcal/mol. It stayed `ok=True` at every
   kappa tried, including the attempt-4 vector, while the neutral proton-transfer TS frames were
   moving by tens of kcal/mol under kappa. "Guard passed" therefore carries no information about
   the property that matters for a default.
5. *The `mg3` precedent does not transfer.* `mg3` was promoted on a measured, monotone gain on
   its own falsifier (class-A rms 22.15 -> 13.22, median D_e error −14 -> −7 kcal/mol) and on
   11 700 paired trajectories showing no tail cost (`ROADMAP.md:228-245`). SQE has no equivalent
   measurement — its only positive number is a p0 side effect on one subset.
6. *The one real gain is worth keeping as an option.* BH76 MAD 48 -> 39 at p0 is genuine and
   reproducible (three independent re-measurements in `KAPPA_FIT_STATUS.md`), and it is exactly
   what an opt-in preset is for while the model question of Q1.4-B is open.

**Gate for flipping the default** (all five, measured, numbers recorded in the status file):

| # | criterion | number |
|---|---|---|
| G1 | kappa_Z from a fit on the cleaned class-E set (Q1.4 Layer A) is reproducible from two starting vectors | each kappa_Z within 0.1 Eh between the two runs |
| G2 | Cl2- and F2- curve SHAPE, not one point: E(r) − [E(X)+E(X-)] vs the r2SCAN-3c fragment-anchored values for r <= 1.2 r_eq | rms <= 5 kcal/mol, minimum within 0.1 A of the reference |
| G3 | BH76_anionic improves and charged NCI does not pay for it | BH76_anionic MAD >= 20 % below `eeq`; AHB21/CHB6/IL16 MAD each not worse than `eeq` by more than 1 kcal/mol |
| G4 | the roadmap's own matrix on neutral chemistry | GMTKN55 conformers <= 1.6, NCI <= 8 kcal/mol MAD (`gfnff` 1.49 / 7.76) |
| G5 | react-MD acceptance 4 with the fitted kappa | jump statistics inside the stage-1 numbers, n >= 3 per cell, and NVE slope (Q2.4) not worse than `eeq` |

G2 is the one that cannot be met without deciding Q1.4-B; that is why the model decision comes
before any further fit campaign.

---

## Q4. Missing guards (S66, conformer set) — real gap; wire it now as report-only, arm it with Q3

**Finding.** The omission is real and it is not academic: Q1.3 showed the one side effect of
kappa that a charged-species campaign will not see — neutral molecules' charges and barriers
move — and the class-D guard cannot detect it (Q3, point 4). S66 (neutral NCI) and the
conformer subsets are the two cheapest sets that would. But a *penalty* against a metric the
fit does not yet trust would only fight noise on the current loss (Q1); so: wire the numbers
into every evaluation's report now, arm the penalty when the Q3 gate is being assessed.

**Where the data already are (checked).**
- S66 is a GMTKN55 subset and is on disk: `test_cases/GMTKN55-testset/S66/{01, 01A, 01B, …}`
  (66 dimers, 198 structures). It is not in `fetch_testset.py`'s registry as its own entry and
  does not need to be — the `gmtkn55` entry brought it. `gmtkn55_reactions.py` lists it under
  `CHEM_CLASS["nci"]` (`:52`), and `load_reactions(["S66"])` scores each dimer as a
  reaction `AB − A − B`, which is what the existing `barriers` mechanism in `revgfnff_fit.py`
  already consumes ("scores ANY named subset as a general reaction/stoichiometry residual",
  WORK_STATUS 14.2). Zero new loader code.
- The conformer set is `CHEM_CLASS["conformers"]` = ACONF, Amino20x4, BUT14DIOL, ICONF, MCONF,
  PCONF21, SCONF, UPU23 (`gmtkn55_reactions.py:50`); the roadmap's "285-reaction" figure is
  that class (305 reactions if all eight are taken; the exact count depends on which of the
  eight the WP0 baseline used — take the eight and record n). Same loader.
- The thresholds already exist and are absolute, not relative: the validation matrix in
  `REV_GFNFF_ROADMAP.md` fixes `gfnff` at conformers 1.49 / NCI 7.76 / charged NCI 38.9 and
  the `revgfnff` targets at <= 1.6 / <= 8 / <= 38.9 kcal/mol.

**Plug-in sketch (PROPOSED, ~40 lines in `revgfnff_fit.py`, Sonnet-executable).**

```json
"guards": [
  {"name": "s66",        "barriers": ["S66"],                                        "max_mad": 8.0,  "factor": 1.10},
  {"name": "conformers", "barriers": ["ACONF","Amino20x4","BUT14DIOL","ICONF","MCONF","PCONF21","SCONF","UPU23"],
                                                                                    "max_mad": 1.6,  "factor": 1.10},
  {"name": "charged_nci","barriers": ["AHB21","CHB6","IL16"],                        "max_mad": 38.9, "factor": 1.10}
],
"guard_every_eval": false
```

1. `load_barrier_data()` gets the union of fitted and guard subsets, so the structures join the
   existing per-(charge, spin) batch groups — one more batch per group, no new run path.
2. In `Evaluation.evaluate()`, after `barrier_subset_stats` (`:608`): for each guard, MAD over
   its subsets; `limit = min(factor * baseline, max_mad)` with `baseline` = the MAD at p0 under
   `fixed_override` (i.e. `sqe` with kappa = 0 == `eeq`, computed once like
   `_compute_d_guard_baseline`, `:525`); `excess = max(0, MAD − limit)`; append
   `sqrt(1e3) * excess` to `residual_chunks` exactly as the class-D guard does (`:635-644`), and
   report `{baseline, current, limit, ok}` per guard next to `guard_d`.
3. `guard_every_eval: false` evaluates the guards only at accepted LM/NM steps and in
   `--evaluate-only`, not in the 6 Jacobian probes per iteration — the conformer set is ~600
   structures of up to ~70 atoms (order 10 s per evaluation at 16 jobs, against ~0.6 s for the
   current loss), so this keeps the iteration cost within 2x. A guard is a gate, not a
   direction; it does not need a derivative.
4. Report-only first: ship with `"factor": null` meaning "print, do not penalise" (the code
   path already exists for the old class-A gradient guard, `:646-653`). Switch the three to
   penalties in the same commit that starts the Q3 gate assessment.

**Why the relative factor alone is not enough.** With `gfnff` at 1.49 / 7.76, a 1.10 factor
would tolerate 1.64 / 8.5 — looser than the matrix's 1.6 / 8. Taking `min(factor*baseline,
max_mad)` makes the guard exactly the matrix, and the relative part only matters if the
baseline is already below the absolute target by more than 10 %.

**Order of work.** Layer A of Q1.4 and this guard wiring are the same kind of change in the same
file and should go to one Sonnet agent in one pass, with the guard in report-only mode; the
fifth fit then runs with the guards visible. That way the first honest kappa fit already
carries the numbers Q3's gate needs, instead of discovering a conformer regression afterwards.

---

## Measurement provenance

- Binary: `build_rev/curcuma`, md5 `a7c8c85165f0912c6bc7c160da036a97` (= the attempt-4
  binary; uncommitted packages 15/16 in the working tree, HEAD `9c69d40e`).
- Q1.1: `scratchpad/kprobe` decomposition of `fit_work/stage2_kappa_run4/*.jsonl` vs
  `ref/E/*/energies.json`; reproduces `rms_dE = 533.8496`.
- Q1.2: `scratchpad/kprobe/probe.json` (react, fit order; 8 systems x 5 charge settings),
  `scratchpad/kprobe/assoc/` (react, reversed frames, 2 systems x 5), `scratchpad/kprobe/static/`
  (8 single frames x 4). Fragment anchors: Cl + Cl- = −0.96970056 Eh, F + F- = −0.88314547 Eh
  (`revgfnff`, one-atom single points; kappa-independent by construction).
- Q2.2: `scratchpad/q2/walk.xyz` + six run directories; h = 1e-5 A, central differences.
- Test binaries as built: `test_gfnff_sqe` PASS (2c residual 2.357e-4, tolerance 3e-4),
  `test_gfnff_rev_fd` PASS (numbers quoted in Q2.3).
- No file under `src/`, `scripts/` or `test_cases/` was modified; nothing committed.
