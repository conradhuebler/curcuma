# rev-gfnff react join tightened + dE_jump "n/a" (Sep 12, 2026, branch reactff2-llm, build/, uncommitted)

## A. Change
- `-gfnff.rev_form_switch` (String, default **order**; `weight` = previous behaviour, bit-for-bit).
- `-gfnff.rev_bo2_form` (Double, **0.1**): threshold on the NARROW bond order (rev_bo2_*).
- order mode, fresh pair only (fading wells and 1,3 pairs keep their rules): join at `bo2 > rev_bo2_form`;
  transition runs on that same switch (`tr.tight`), window `rev_bo2_form..0.9`; the pair's well is left to
  the blend (`RevTransition::well_blend`; `beginTransition` no longer copies it into the old corners);
  revert at `c < rev_tr_begin`.
- unmeasurable `dE_jump` prints `n/a`, event keeps NaN, statistics skip it (never 0.0).

## B. Switch values  (thr = (rcov_i+rcov_j)*fat_i*fat_j; wide = rev_bo_* 2.0/-7.5, narrow = rev_bo2_* 1.4/-6)
| geometry | r | r/thr | wide | narrow |
|---|---:|---:|---:|---:|
| water dimer H3...O4 (H bond) | 1.952 A | 1.954 | **0.5956** | **3.9e-04** |
| water dimer O1...O4 | 2.910 A | 2.185 | **0.1631** | **9.7e-07** |
| H-H formation radius of the non-rev react mode (test 14) | 1.6x thr | 1.600 | 0.983 | 0.113 |
| O-H / C-C at equilibrium | 0.96 / 1.39 A | 0.96 / 0.93 | 1.000000 | 0.996 / 0.998 |
wide crosses 0.05 at **2.310x** thr (inside H-bond range); narrow crosses 0.1 at **1.611x**, which coincides
with (a) `react_form_factor` 1.6 (non-rev react, tests 13-15), (b) the break radius (coordinate 0.5 == 1.600x),
(c) `rev_bo13_form` 0.1. Insensitive: 0.02 (1.739x) and 0.05 (1.671x) also give 0 water-dimer formations.

## C. Water dimer, revgfnff react, 1 ps, dt 0.25, no thermostat
`-md.seed` does not change the initial velocities here (42/43/44 byte-identical) -> 3 runs at 300/310/320 K.
| T | before f/b/rebuilds | after | <Epot> after | <Epot> gfnff static | diff |
|---|---|---|---:|---:|---:|
| 300 K | 18 / 15 / 33 | **0/0/0** | -0.662305 | -0.662307 | 0.0013 kcal/mol |
| 310 K | 15 / 12 / 27 | **0/0/0** | -0.662163 | -0.662163 | 0.0000 |
| 320 K | 13 / 11 / 24 | **0/0/0** | -0.662059 | -0.662059 | 0.0000 |
before at 300 K: <Epot> -0.668734 Eh = **-4.03 kcal/mol** off static. 0 formations also at 400/600/800 K.
## D. Genuine formation and tests
| system | before | after |
|---|---|---|
| test 14 (`-method gfnff`, legacy criterion, untouched) | 3 formed / 3 broken / 6 rebuilds | **identical**, same jumps |
| H4 6000 K wall 2.5 A, 5 ps, revgfnff | 94/94/240, med 0, max abs 0.0335 Eh | 17/17/50, med 0, **max 3e-05 Eh** |
| H4 same with `rev_form_switch weight` | - | 94/94/240, jump stats **bit-identical to before** |
| test 16 (2 H2 4000 K CSVR) | 66 events, med 0.0, max 0.8 kJ/mol | **identical** (its events are fading-well re-forms) |
| test 17 (N2+3H2 3500 K) | 220 ev, med 0.0, max 4.0, T_max 4602 K | 172 ev, med 0.0, max 7.5, frac<5 0.994, T_max 3986 K, PASS |
| test 18 system (2 H2 3000 K NVE dt 0.25) | 2 formed + 2 broken | **0 events** (see G1) |

`ctest -R gfnff` in build/: **63/64**, only `cli_simplemd_18_gfnff_rev_nve_vs_gfnff` fails (re-specified by
another agent; its failure line is the old endpoint criterion, ratio 101.3). `gfnff_rev_fd` PASSES (it failed
transiently while the transition still ran on the bo3 coordinate: 6.9e-04 -> 2.1e-05 Eh/A).

## E. Guards (before vs after)
- `gfnff` static SP: C6H6 -2.36477552, caffeine -4.67273707, acetic acid dimer -2.47120432 Eh: IDENTICAL;
  full `-verbosity 3` dump (1887 / 2644 lines, per-bond factors to 12 digits, terms to 10 decimals) BYTE-IDENTICAL.
- `revgfnff` static SP: C6H6 -2.36477475, caffeine -4.67273522 Eh: IDENTICAL; dumps (1890 / 2647 lines) byte-identical.

## F. dE_jump n/a
`REACT rebuild #1: 5 bonds, dE_jump = n/a (no previous-step state)` (was `nan Eh (+0.0 kJ/mol)`); the literal
`dE_jump` is kept (test 14 counts rebuild vs dE_jump lines). `ReactEvent::de_jump_eh` stays NaN
(`finishTransition` no longer writes 0.0); `Results()["react"]["dE_jump_kJ"]` holds only measured values plus
`["dE_jump_unavailable"]`; summary prints `energy jumps (32 measured, 1 n/a)`. `revgfnff_jump_stats.py` and
`react_baseline.py` updated; `revgfnff_jump_by_event.py` needs none (its `w_scan` regex still matches the
formation line, which now also prints `o_scan` and `r/rcov`).
**Cause**: the scan runs at the top of `Calculation()`; at an event in the first step the workspace has not yet
been given this step's CN (`setD3CN`) and EEQ charges (`setEEQCharges`) - both are set later in the same call.

## G. Open
1. **test 18's system no longer reacts at 3000 K** (4 events -> 0); its drift ratio rev/gfnff at dt 0.25 is then
   1.06 (-1.289e-04 vs -1.214e-04 Eh/ps, 10 ps linear fit). A hotter setting is needed if it must show events.
2. A **snapped/demoted** transition now costs s*E_well (the well is no longer in the old corners). Not seen in the
   final runs; it gave one +190 kJ/mol event per run while fading wells were not yet exempted. Lever: `rev_max_transitions`.
3. `test_gfnff_rev_fd.cpp`: the two `1.7->1.3->...` rows no longer start a transition in order mode; added
   `1.7->1.05->0.95 A` (s ~ 0.33). Its residual 1.03e-04 Eh/A is EXACTLY plain gfnff's own FD residual at
   r = 0.95 A (0.90 A: 4.024e-04 for both, and it does not shrink at h = 1e-6) - GFN-FF's short-H-H gradient
   residual, not a blend defect.
