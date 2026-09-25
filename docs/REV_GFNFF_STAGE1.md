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

**FIXED (Sep 2026, AI/machine-tested): a bare alkali/alkaline-earth cation read as grossly
over-coordinated.** `revValence()` (`gfnff_method.cpp`) gave Li/Na/K (Z=3/11/19) a nominal
valence of 1.0 and Be/Mg/Ca (Z=4/12/20) 2.0 - aimed at simple molecular compounds (LiH, MgO) - and
`revOverP()` falls back to the generic default `p_Z = 0.3 Eh` for any element without its own
override, which these never got (the WP3 fit only ever touched H/C/N/O). Found via GMTKN55 CHB6
(charged H-bond / cation-pi complexes) during the stage-2 kappa_Z fit's dataset check: a bare Li+
sitting non-covalently above a benzene ring accumulates a non-trivial continuous bond order
`b2_ij` to all six ring carbons simultaneously (none of them a real sigma bond), so `sum_j b2_ij`
easily exceeds `Val_Z = 1`, and the term contributed **+4.97 Eh** where the correct total energy
is -1.16 Eh (confirmed against both plain `gfnff`, unaffected, and `xtb --gfnff`, -1.155727 vs
-1.155726 Eh) - the entire ~3000 kcal/mol CHB6 residual on the three cation-pi reactions. Fixed by
letting Z=3/4/11/12/19/20 fall through to the same `default: return 6.0` ("hypervalence-capable
main group and metals: effectively no penalty") every other metal already gets - nothing here was
ever calibrated against alkali/alkaline-earth chemistry either, so nothing was traded away; only
the failure mode is gone. CHB6 (all 6 reactions) MAD 1550.6 -> 47.6 kcal/mol, RMS 2161.3 -> 66.1
(worst single reaction now -145.3, the Li+-benzene case - genuinely hard for any classical force
field, not a bug). `ctest -L gfnff`: same 66/69 as before the fix (the pre-existing
`cli_curcumaopt_07`/`cli_simplemd_18`/`cli_simplemd_20`, unrelated) - plain `gfnff` is untouched
since `revValence()` is only reachable through `-method revgfnff`. See
`test_cases/revgfnff/_log/WORK_STATUS.md` package 15.

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
disables).

> **Time-scale note (Sep 2026).** Every "fs" on this page predates the MD time-step unit fix
> (`CurcumaUnit::Constants::MD_TIME_UNIT_FS`, `AIChangelog`), so it must be multiplied by
> **1.9516144** to read what was actually integrated: the cap's 0.25 was really 0.488 fs and the
> "0.5 fs with the cap off" runs were 0.976 fs. The **PARAM value was not changed**, so the cap
> now clamps to a genuine 0.25 fs — a real tightening by that factor, and every stability number
> below is therefore conservative. Re-measured tail statistics at true fs:
> `docs/REV_GFNFF_STAGE3A.md` section 2.2.

Final numbers (3 x 3 runs): at 0.25 fs 2 H2 78-97 % of the jumps below 1 kJ/mol and
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

## Formation criterion: which switch joins a new bond (`rev_form_switch`, Sep 12, 2026)

Until Sep 12, 2026 a non-bonded pair joined the bond list once the **wide term-weight switch**
(`rev_bo_center` 2.0, width -7.5) rose above `rev_bo_form` = 0.05. That switch crosses 0.05 at
**2.31x the covalent sum**, which is inside hydrogen-bond and van-der-Waals contact range, so the
criterion joined contacts that are not bonds. Measured on the Cs water dimer (O-O 2.910 A, donor
H...O 1.952 A), `revgfnff -gfnff.topology_mode react`, 300 K, 1 ps, dt 0.25 fs, no thermostat:

| switch at | H...O 1.952 A (1.954x) | O...O 2.910 A (2.185x) | value at a bond (0.96x) |
|---|---:|---:|---:|
| wide term weight `rev_bo_*` (2.0/-7.5) | **0.5956** | **0.1631** | 1.000000 |
| narrow bond order `rev_bo2_*` (1.4/-6) | **3.9e-04** | **9.7e-07** | 0.9961 |

The wide switch reads the hydrogen bond as 60 % of a bond: both contacts joined at t = 0 and the
run produced **18 formations / 15 breaks / 33 rebuilds per ps** on a molecule that does not react,
with a mean Epot **4.03 kcal/mol below** the same run in static topology. The narrow switch reads
3.9e-04 there - a factor 260 below even a 0.1 threshold.

**Default since Sep 12, 2026: `rev_form_switch = order`** - a fresh pair joins once its NARROW bond
order exceeds `rev_bo2_form` (0.1, crossing at **1.611x** the covalent sum). Three independent
radii coincide there, which is why 0.1 was chosen rather than a tighter or looser value:

- the non-rev react mode has formed at `react_form_factor` = **1.6x** since it existed (tests 13-15),
- a bond starts **breaking** at transition coordinate 0.5, which is exactly **1.600x**, so formation
  and break are now symmetric about the same radius,
- it is the value `rev_bo13_form` already uses for a 1,3 ring closure on the same switch.

The result is insensitive to the exact number: 0.02 (1.739x) and 0.05 (1.671x) also give zero
formations on the water dimer, as does the default up to 800 K.

Two things follow from joining at a radius where the well is no longer ~0, and both are part of the
mode rather than of the threshold:

1. **The transition runs on the narrow switch too** (`tr.tight`, as a 1,3 closure always has),
   window `rev_bo2_form` -> 0.9, i.e. 1.611x -> 1.189x. On the wider bo3 coordinate the join would
   already sit at c = 0.47 against a window end of 0.8 - 0.06x the covalent sum for H-H, and the FD
   gradient of an active transition degraded to 6.9e-04 Eh/A there (2.1e-05 with the narrow window).
2. **The forming pair's own well is not copied into the old corners** (`RevTransition::well_blend`):
   at 1.611x the wide term weight is already 0.982, so copying it would switch the full well on at
   the event - which is exactly the non-rev hard swap (measured |dE_jump| up to 0.158 Eh on the H4
   test). With the well left to the blend the event is energy-neutral by construction (s = 0).
   A **fading well** is exempt: it re-forms on the bo3 coordinate and is already in every corner,
   so it keeps the old treatment (dropping it on a revert cost +190 kJ/mol until that was fixed).

`-gfnff.rev_form_switch weight` restores the previous behaviour bit-for-bit (verified: water dimer
18/15/33 and H4 94/94/240 with identical jump statistics).

Measured effect, 4 free H atoms, 6000 K, wall 2.5 A, dt 0.25 fs, 5 ps, `revgfnff`:

| | formations | rebuilds | median \|dE_jump\| | max \|dE_jump\| |
|---|---:|---:|---:|---:|
| `weight` (previous) | 94 | 240 | 0.00 Eh | 0.0335 Eh |
| `order` (new default) | 17 | 50 | 0.00 Eh | **3e-05 Eh** |

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

**Scope of these numbers: diatomics only, and the tail does not transfer.** Every cell above is
2 H2, 4 H or N2 + 3 H2, where a bond event cannot re-parametrise a molecular fragment. The full
polyatomic baseline is now measured (5 systems spanning C-O / C-N / C-C / N-N / O-O bonds,
1000/2000/3000 K, 3 start frames each, 5 ps, 11527 rebuilds; full tables in
`test_cases/revgfnff/_log/POLY_JUMP_BASELINE.md`):

| | n | median | max | < 1 kJ | < 5 kJ |
|---|---:|---:|---:|---:|---:|
| polyatomic, all (order criterion) | 11527 | 0.00 | **423.5** | 0.966 | 0.986 |
| ... formations only | 5617 | 0.00 | 78.7 | 0.996 | - |
| ... breaks only | 5630 | 0.00 | **423.5** | 0.934 | 0.974 |
| diatomic (cells above, 12 runs, older binary) | - | 0.0-0.4 | 0.7-45.3 | - | - |

So the **bulk matches** (median 0.00 kJ/mol, 96.6 % below 1) and the **tail is ~10x larger**:
423.5 vs 45.3 kJ/mol. The sharp result is that the roughness is almost entirely on the
**break** side - 58 of the 61 events >= 50 kJ/mol are bond breaks, formations are effectively
smooth (3 of 5617 above 50 kJ, 99.6 % below 1), and the three worst events are all H-H breaks in
ethane (3000 K: +423.5 at t = 1849 fs and +394.4 at t = 2011 fs; 2000 K: +310.3), where the bond
term carries +452 kJ/mol and the bonded/non-bonded repulsion partially cancels it. The worst
temperature dependence is equally clear: 1000 K max 2.9, 2000 K max 310.3, 3000 K max 423.5
kJ/mol. **Stage 3 must therefore measure the break side, not the formation side.**

Two earlier numbers in this file were wrong and are superseded: a methanol run reported "median
44.5 / max 69.2 kJ/mol at 1000 K" from **n = 3 events**, which a 5 ps re-measurement over 84
events at the same temperature (median 0.00, max 0.00) shows was unrepresentative - a median from
three events is not a statistic. The class-A bond scans remain valid as *signed* per-term
information (8 polyatomic bonds at -25..-85 kJ/mol per rebuild, the new topology lower, i.e.
spurious heating) but they sample single stretched bonds, not MD event distributions.
`-method revgfnff` without `-gfnff.topology_mode react` produces **zero** rebuilds - the jump
statistics exist only in react mode.

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
(0.02/0.8/0.5/0.75), `rev_bo13_form` (0.1), `rev_form_switch` (order), `rev_bo2_form` (0.1), `rev_max_transitions` (4), `rev_demote_cooldown`
(0), `rev_blend` (true), `rev_over_p/k/shift`. Diagnostics: `CURCUMA_REACTSCAN=1` prints every
stretched bond and every candidate pair per scan (call, r, w, c, fading/1,3/in-transition
flags); verbosity 3 prints the per-term jump of every event.
