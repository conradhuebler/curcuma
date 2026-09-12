# rev-gfnff stage 1: continuous bond order, blended repulsion, over-coordination, multi-transition blending

AI-generated (Sep 2026), machine-tested only. Method name `revgfnff` (alias `gfnff-rev`);
`gfnff` is untouched (golden `gfnff_val_*` ctests bit-identical). CPU only. Roadmap:
`docs/REV_GFNFF_ROADMAP.md` WP3; data: `docs/REV_GFNFF_DATA_BASIS.md`.

## Stage 1a: switches and terms (`rev_bond_order.h`, kernels in `ff_workspace_gfnff.cpp`)

All switches have the erf form `b(r) = 1/2 (1 + erf(k (r - R)/R))`, `R = f (rcov_i + rcov_j) fat_i fat_j`
(the react scan's threshold expression). Five of them, with different roles:

| switch | f | k | used for | why |
|---|---|---|---|---|
| term weight `w` | 2.0 | -7.5 | bond well and angle/torsion/inversion damping x the SHIFTED weight `(w - w_join)/(1 - w_join)` (w_join = `rev_bo_form` 0.05); react scan join (w > 0.05) and well removal (w < 0.02) on the raw w | wide, so the Gaussian well keeps its tail; with the shift a pair joins and leaves the lists at exactly zero weight (the raw join cost 1-5 kJ/mol and was the last visible jump) |
| bond order `b2` | 1.4 | -8 | over-coordination sums (topological 1,3 pairs excluded); formation threshold of a 1,3 pair (ring closure, `rev_bo13_form` 0.1) | 1 at the bond, ~0.03 at a heavy-atom 1,3 distance; a fit may soften it |
| repulsion blend, bonded list `b4` | 1.7 | -12 | `b4 E_bonded + (1-b4) E_nonbonded` for pairs in the bond list | 1 - 1e-9 at 1.1x and 0.9993 at 1.38x (the turning point of a hot X-H bond): at 1.5x/-10 the 4.5 % of the huge non-bonded repulsion there blew up 0.5 fs MD; 0.5 at 1.7x where a bond really breaks |
| repulsion blend, non-bonded list `b5` | 1.3 | -12 | the same blend for pairs NOT in the bond list (topological 1,3 excluded) | ~0 at 1,4 and H-bond distances, so nothing leaks at equilibrium; a pair that approaches joins the bond list at 2.3x, where both switches read ~0, so the list move is continuous |
| transition coordinate `c` | 1.6 | -8 | stage 1b: the blend coordinate s of a topology change (below) | plateau (> 0.9) out to 1.4x, so a hot X-H bond does not chatter; 0.02 at 1.9x |

Over-coordination: `E_over,i = p_Z softplus_k(sum_j b2_ij - Val_Z - 0.5)^2`, sigma partners only
(k = 10, p = 0.3 Eh, per-element overrides via the `rev` section of `-gfnff.param_file`).

**Why the switches were softened (Sep 12, 2026).** The first version used `b2` at 1.3x/-16 for
the repulsion blend and E_over. Every `revgfnff` MD run then exploded within 1 ps (T = 1e9 K),
also with a static topology, while `gfnff` on the same input was fine and dt = 0.25 fs cured it.
That is an integrator instability, not a bookkeeping error: a hot H2 at 3000 K reaches 1.32-1.38x
its covalent sum, i.e. the centre of that switch, so at every outer turning point 60 % of the
non-bonded repulsion and the E_over slope switched on with curvatures of 10-20 Eh/Bohr^2 (the H2
bond itself: 0.37). Turning single terms off made it worse (they partly cancel). Repulsion blend
and E_over now run on 1.4x/-8 (curvature 4x the H2 bond instead of 20x). That is still marginal
at 0.5 fs (T spikes to 30 kK on 2 H2 / N2+3H2 even with a static topology), so **`revgfnff` runs
at dt = 0.25 fs**; there all three demo systems stay at their thermostat temperature.

**X-H stability at 0.5 fs (Sep 12, 2026, second pass).** With the switches above a static
`revgfnff` topology is stable at 0.5 fs (2 H2 3.8 kK, N2 + 3 H2 3.7 kK, was 1e9). Reactive runs
at 0.5 fs still blow up occasionally (1 of 3 runs on 2 H2 and on the 4 H square): a
re-parametrisation of several hundred kJ/mol crossed within one step catapults a hydrogen,
which then crosses other windows in one step. Tried and measured: a softer transition switch
(k = -5) is LESS stable (more overlapping transitions); a softer E_over (softplus k = 5) leaks
0.3 kcal/mol at equilibrium and makes the 4 H square explode; forcing a two-coordinate H to
hyb = 0 (no sp-sp well doubling) makes the jumps larger because the angle term at that H then
takes a tetrahedral theta0 - the consistent treatment is stage 3. What did help: every
transition window starts at the coordinate the scan actually saw (s = 0 by construction even
when c moved 0.15 in one step) and s is a smoothstep in the window (C1 at both ends). So
`revgfnff` MD is capped at **dt = 0.25 fs** (`-md.rev_dt_cap`, plain-printed warning + clamp, 0
disables). Final numbers (3 x 3 runs): at 0.25 fs 2 H2 78-97 % of the jumps below 1 kJ/mol and
100 % below 5 (max 1-3), N2 + 3 H2 96-99 % below 1, 99-100 % below 5 (max 4, one run 43); at 0.5 fs
with the cap off 2 H2 and N2 + 3 H2 no longer explode (T_max 9.7-12.5 / 3.7-6.3 kK) but only
80-93 % of the jumps stay below 5 kJ/mol, and the 4 H square still blows up in 1 of 3 runs.

## Stage 1b: multi-transition blending over topology corners

The discontinuity of the reactive mode is not the well of the pair that changes (that is ~5
kJ/mol at w = 0.05) but the **re-selected parameters of the neighbours** (an H with two partners
is perceived sp and its H-H bond takes `bsmat[1][1] = 1.98` instead of 1.00; a nitrogen going sp
-> sp2 swaps the N-N bond strength). Stage 1b makes every topology change a continuous blend:

- A **transition** t is a pair (i, j) forming or breaking, with a coordinate `s_t in [0, 1]`
  (smoothstep) over a window of the transition switch: formation from the coordinate seen at
  detection (at least `rev_tr_begin` 0.02) to `rev_tr_end` 0.8, break from the coordinate seen
  (at most `rev_tr_prebreak` 0.5) to 0.02; a formed bond "sticks", so it does not chatter; a
  1,3 ring closure runs on `b2: 0.1 -> 0.9`.
- Up to `rev_max_transitions` (4) transitions run at once; the force field is evaluated on every
  **corner** of the topology hypercube (2^k complete parameter sets: bonded and non-bonded lists,
  Coulomb self-energy, e0, and an own EEQ solve per corner per step) and blended,
  `E = sum_b W_b E_b`, `W_b = prod_t (b_t ? s_t : 1 - s_t)`, with the gradient
  `sum_b W_b grad E_b + sum_t (dE/ds_t) (ds_t/dr_t) r_t`. Each corner is an `FFWorkspace::TopologyState`
  swapped into the evaluation slot in O(1) (`swapState`), so no kernel knows about the blend.
- The **well of a transition pair is in every corner** (a forming pair's Bond entry is copied
  into the corners without the bond, a broken bond's entry stays as a *fading well* until
  w < 0.02), so the blend only carries the re-parametrisation. Wells of in-flight formations are
  also appended to corners generated later (`appendFadingWells`).
- Ends: `s = 1` completes (keep the corners with the bit set), `s = 0` plus hysteresis reverts
  (`rev_tr_revert` 0.75 for a break; w < 0.025 for a formation). Both are exactly continuous
  (measured 0.0 kJ/mol). When a fifth event arrives the transition closest to either end is
  snapped there (measured, reported as `promote`/`demote`).
- The 1,3 test of the scan uses only **settled** bonds (transitions below s = 0.5 do not count),
  otherwise a pre-bond at w = 0.05 turns a third atom's ordinary well join into a "ring closure"
  that waits for the tight switch and then joins with the full well.
- The scan never holds a formation back (a delayed well join *is* the hard swap); a demoted
  formation restarts at once (`rev_demote_cooldown` 0).

Everything is switchable: `rev_blend false` gives stage 1a (measured for comparison).

## Measured (dt = 0.25 fs, 2-5 ps, CSVR, spherical wall; `scripts/revgfnff_jump_stats.py`)

Per-event jump `dE_jump` = E(after) - E(before) at the same geometry at every topology event,
split per term (`REACT jump terms`, verbosity 3):

| run | median abs dE_jump [kJ/mol] | share < 1 kJ/mol | share < 5 | max |
|---|---:|---:|---:|---:|
| 2 H2, 3000 K, no blend (stage 1a) | 20.5 | 0.38 | 0.45 | 1657 |
| 2 H2, 3000 K, stage 1b | 0.4 | 0.33-0.62 | 0.90-1.00 | 5-8 |
| 4 H, 2000 K, no blend | 79 | 0.18 | 0.29 | 413 |
| 4 H, 2000 K, stage 1b | 0.1 | 0.40-0.43 | 0.77-0.88 | 93 (start artefact, below) |
| N2 + 3 H2, 3500 K, no blend | 192 | 0.02 | 0.12 | 1020 |
| N2 + 3 H2, 3500 K, stage 1b | 0.1 | 0.89-0.98 | 0.98-1.00 | 4-9 |

Ranges are three runs at 0.95/1.0/1.05 of the temperature (`-md.seed` does not change the
initial velocities, so the temperature is the only cheap lever for statistics). Temperatures
stay at the thermostat value (2 H2 / N2+3H2: T_max 6-9 / 3.8-4.3 kK). Final numbers (all fixes,
3 x 3 runs): N2 + 3 H2 95-96 % of the jumps below 1 kJ/mol, 99-100 % below 5, max 4-7 kJ/mol;
2 H2 87-90 % below 5, max 6 (one run has no event at all in 2 ps once the repulsion leak that
used to dissociate hot H2 is gone). Completions,
reverts and demotes cost exactly 0.0 kJ/mol; the Coulomb jump of a fragment merge (up to -50
kJ/mol at s = 0 before) is gone because the EEQ is solved per corner.

**Remaining tail.** (1) The 4 H square start geometry (1.3 A, all pairs at w = 0.5) begins with
six wells joining at half weight - a property of the test input (seeding the initial topology with
those pairs was tried and rejected: it also seeds their hybridisation). (2) The well join itself at w = 0.05 costs 1-5 kJ/mol
(`rev_bo_form` could go to 0.02 at the price of more candidates); no jump above 40 kJ/mol
remains on 2 H2 and N2 + 3 H2 (3 x 3 runs). (3) Chatter: hot bonds cross `c = 0.5` and come back (100-200 reverts per
3 ps), each cycle two parameter generations and 2^k evaluations; harmless for the energy, a
cost for large systems. (4) The blend force `dE/ds ds/dr` of a re-parametrisation worth
300-600 kJ/mol over a ~0.5 Bohr window is 0.2-0.3 Eh/Bohr; that is what makes dt = 0.25 fs
necessary and what a physically sensible re-parametrisation (an H is never sp; stage 3) will
shrink.

**Stage 2 (split charges) on the same two systems.** `docs/REV_GFNFF_STAGE2.md` repeats this table
for `-gfnff.rev_charge_model sqe` (kappa_H = kappa_N = 0.2 Eh) against an `eeq` baseline of the
same binary and seeds, n = 3 each: both stay inside the numbers above, and the single ~43-45
kJ/mol N2 + 3 H2 outlier that this table records as "one run 43" appears once in three runs for
BOTH charge models.

- **Gradient**: `test_gfnff_rev_fd` (H2 scans, CH4 + H, stretched N2H2, and four react
  sequences walked into their transition windows) agrees with central finite differences to the
  same residual as plain `gfnff`.
- **NVE**: `revgfnff` and `gfnff` drift identically on 2 H2 at 3000 K (2.1e-4 / 8.4e-5 / 5.2e-5
  Eh per ps at dt 0.25 / 0.125 / 0.0625 with no events) - but that drift is LINEAR in dt for
  plain `gfnff` too, unlike the dt^2 scaling Known Issue #28 measured on CH4/H2O/NH3. Not
  root-caused (not a rev issue; the wall is ~0 there). CLI test 18 therefore compares the two.
- **Equilibrium**: caffeine `revgfnff` vs `gfnff`: total **+0.002 kcal/mol** (-4.67273413 vs
  -4.67273707 Eh; over-coordination 2e-6 Eh, bonded repulsion +2e-6). Two leaks were found and
  closed on the way: the torsion weights were taken on the stored quartet's (i,j) and (k,l)
  pairs, which are the 1,3 distances (the outer atoms hang crosswise, i on k and l on j) and
  therefore read 0.5-0.96 at equilibrium, halving the torsion term; and the repulsion blend
  switch (see the table). Cache hazard found doing this: the `.topo.json` cache replayed a
  `gfnff` parameter set for `revgfnff` (and vice versa); the fingerprint now carries the rev
  settings. Also fixed: `-gfnff.rev_*` and `-gfnff.param_file` were ignored in `-sp` runs (the
  gfnff sub-config arrives after construction; `setParameters` now re-reads it).
- Cost: 2^k force evaluations and EEQ solves per step while k transitions are in flight; a
  locality optimisation (only terms touched by a transition differ between corners) is the
  obvious next step for systems beyond the demo size.

## Parameters (all `-gfnff.<name>`, defaults in `gfnff.h`)

`rev_bo_center/width` (2.0/-7.5), `rev_bo2_center/width` (1.4/-6), `rev_bo3_center/width`
(1.6/-8), `rev_bo4_center/width` (1.7/-12), `rev_bo5_center/width` (1.3/-12), `-md.rev_dt_cap` (0.25 fs), `rev_bo_form/break` (0.05/0.02), `rev_tr_begin/end/prebreak/revert`
(0.02/0.8/0.5/0.75), `rev_bo13_form` (0.1), `rev_max_transitions` (4), `rev_demote_cooldown`
(0), `rev_blend` (true), `rev_over_p/k/shift`. Diagnostics: `CURCUMA_REACTSCAN=1` prints every
stretched bond and every candidate pair per scan (call, r, w, c, fading/1,3/in-transition
flags); verbosity 3 prints the per-term jump of every event.
