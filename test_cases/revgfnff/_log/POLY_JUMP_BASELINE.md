# Polyatomic topology-rebuild jump baseline (rev-gfnff react MD)

AI-generated (measurement job, 2026-09-13). **No src/ change, no build, no git state change.**
One frozen binary copy: `release/curcuma` at HEAD `08373ee7` (branch `reactff2-llm`), md5
`58512a18c523456e8db7ca826ad84897`. Re-run with that same copy; `release/` may be rebuilt.

**Grid per cell:** `-md <frame>.xyz -method revgfnff -gfnff.topology_mode react -temperature T
-maxtime 5000 -md.time_step 0.25 -md.thermostat csvr -md.coupling 10 -md.rattle_12 false
-md.no_restart -md.seed 42 -threads 1 -verbosity 1 -no_bmt` -- 5 ps, no wall, fresh dir per run,
`*.topo.json` deleted first. 1000/2000/3000 K x **3 replicates from 3 START FRAMES** (0/8/16 of
`fit_work/<sys>_1000K.xyz`, 0/5/10 of `n2h4_H.xyz`, 0/60/120 of a self-made 2 ps h2o2 seed traj).
**`-md.seed` is inert** (seeds 42/43/99 byte-identical, n=128; a repeat run is identical too).
Systems: **ch3oh** C-O, **ch3nh2** C-N, **c2h6** C-C, **n2h4_H** (N2H4+H) N-N, **h2o2** O-O --
all 4-8 atom polyatomics with existing multi-frame trajectories, covering the five bond types and
separating a heavy-atom dissociator (c2h6) from hydrogen-exchange-only systems.

## 1. Headline: polyatomic vs diatomic rebuild jumps (kJ/mol)

| set | cells/runs | events | pooled median | pooled max | <1 kJ | <5 kJ | worst |
|---|---|---|---|---|---|---|---|
| polyatomic, all T (order) | 45 cells | 11527 | 0.00 | 423.50 | 0.966 | 0.986 | 423.5 |
| ... breaks only | 45 cells | 5630 | 0.00 | 423.50 | 0.934 | 0.974 | 423.5 |
| ... formations only | 45 cells | 5617 | 0.00 | 78.70 | 0.996 | 0.998 | 78.7 |
| diatomic (2 H2 / N2+3H2, 6 tags) | 12 runs | n/a | 0.0-0.4 (per run) | 0.7-45.3 (per run) | 0.78-1.00 | 0.99-1.00 | 45.3 |

The diatomic block is read from `test_cases/revgfnff/jump_stats/{eeq,sqe}_t{0.95,1.0,1.05}/summary.md`
(12 runs, old binary `build/curcuma`, before the `order` criterion commit) -- not re-run. Those
files report per-run median/max only, so no pooled median exists for them.

**Reading:** the bulk is as smooth as the diatomics (96.6 % of 11527 rebuilds < 1 kJ/mol) but the
TAIL is 10x larger (423.5 vs 45.3 kJ/mol), and **58 of the 61 events >= 50 kJ are BREAKS**; the 3
formation-side ones are all c2h6 3000 K C-H formations at +67..+79 (formations: max 78.7 overall).
The old methanol figure (median 44.5
/ max 69.2, 1000 K) rested on n=3 events; at 1000 K over 629 events it is median 0.00 / max 2.9.

## 2. Per-cell results (order criterion; 5 ps each)

| system | T | frame | rebuilds | form/break | med \|dE\| | max \|dE\| | <1 kJ | <5 kJ | n/a | T_max |
|---|---|---|---|---|---|---|---|---|---|---|
| c2h6 | 1000 | 0 | 0 | 0/0 | - | - | - | - | 0 | 1304 |
| c2h6 | 1000 | 16 | 4 | 2/2 | 0.00 | 0.00 | 1.000 | 1.000 | 0 | 1375 |
| c2h6 | 1000 | 8 | 0 | 0/0 | - | - | - | - | 0 | 1740 |
| c2h6 | 2000 | 0 | 40 | 20/20 | 0.00 | 0.10 | 1.000 | 1.000 | 0 | 3174 |
| c2h6 | 2000 | 16 | 100 | 49/49 | 0.00 | 43.20 | 0.990 | 0.990 | 0 | 3468 |
| c2h6 | 2000 | 8 | 102 | 47/48 | 0.00 | 310.30 | 0.931 | 0.941 | 0 | 2244 |
| c2h6 | 3000 | 0 | 257 | 120/120 | 0.00 | 143.80 | 0.887 | 0.887 | 0 | 4019 |
| c2h6 | 3000 | 16 | 544 | 248/250 | 0.00 | 394.40 | 0.932 | 0.939 | 0 | 3687 |
| c2h6 | 3000 | 8 | 442 | 202/204 | 0.00 | 423.50 | 0.928 | 0.930 | 0 | 4395 |
| ch3nh2 | 1000 | 0 | 12 | 6/6 | 0.00 | 0.70 | 1.000 | 1.000 | 0 | 1839 |
| ch3nh2 | 1000 | 16 | 36 | 18/18 | 0.00 | 0.40 | 1.000 | 1.000 | 0 | 1678 |
| ch3nh2 | 1000 | 8 | 28 | 14/14 | 0.00 | 0.50 | 1.000 | 1.000 | 0 | 1756 |
| ch3nh2 | 2000 | 0 | 220 | 110/110 | 0.00 | 2.80 | 0.991 | 1.000 | 0 | 3431 |
| ch3nh2 | 2000 | 16 | 246 | 123/123 | 0.00 | 1.80 | 0.980 | 1.000 | 0 | 3337 |
| ch3nh2 | 2000 | 8 | 304 | 152/152 | 0.00 | 4.10 | 0.977 | 1.000 | 0 | 2578 |
| ch3nh2 | 3000 | 0 | 412 | 206/206 | 0.00 | 37.60 | 0.964 | 0.995 | 0 | 4648 |
| ch3nh2 | 3000 | 16 | 427 | 211/211 | 0.00 | 31.60 | 0.981 | 0.995 | 0 | 5285 |
| ch3nh2 | 3000 | 8 | 706 | 351/352 | 0.00 | 14.40 | 0.956 | 0.996 | 0 | 3601 |
| ch3oh | 1000 | 0 | 84 | 42/42 | 0.00 | 0.00 | 1.000 | 1.000 | 0 | 1695 |
| ch3oh | 1000 | 16 | 44 | 22/22 | 0.00 | 0.00 | 1.000 | 1.000 | 0 | 1424 |
| ch3oh | 1000 | 8 | 40 | 20/20 | 0.00 | 0.00 | 1.000 | 1.000 | 0 | 1841 |
| ch3oh | 2000 | 0 | 128 | 64/64 | 0.00 | 0.50 | 1.000 | 1.000 | 0 | 3467 |
| ch3oh | 2000 | 16 | 213 | 106/105 | 0.00 | 18.80 | 0.995 | 0.995 | 0 | 3411 |
| ch3oh | 2000 | 8 | 524 | 261/261 | 0.00 | 0.90 | 1.000 | 1.000 | 0 | 4907 |
| ch3oh | 3000 | 0 | 560 | 265/265 | 0.00 | 42.60 | 0.964 | 0.977 | 0 | 3686 |
| ch3oh | 3000 | 16 | 804 | 385/385 | 0.00 | 96.20 | 0.970 | 0.984 | 0 | 4471 |
| ch3oh | 3000 | 8 | 700 | 340/340 | 0.00 | 38.40 | 0.983 | 0.991 | 0 | 5685 |
| h2o2 | 1000 | 0 | 12 | 6/6 | 0.00 | 0.10 | 1.000 | 1.000 | 0 | 2104 |
| h2o2 | 1000 | 120 | 28 | 14/14 | 0.00 | 0.10 | 1.000 | 1.000 | 0 | 1862 |
| h2o2 | 1000 | 60 | 34 | 17/17 | 0.00 | 0.10 | 1.000 | 1.000 | 0 | 1810 |
| h2o2 | 2000 | 0 | 196 | 98/98 | 0.00 | 1.20 | 0.995 | 1.000 | 0 | 6468 |
| h2o2 | 2000 | 120 | 198 | 99/99 | 0.00 | 1.50 | 0.995 | 1.000 | 0 | 2481 |
| h2o2 | 2000 | 60 | 208 | 104/104 | 0.00 | 0.30 | 1.000 | 1.000 | 0 | 4716 |
| h2o2 | 3000 | 0 | 280 | 140/140 | 0.00 | 1.70 | 0.996 | 1.000 | 0 | 6443 |
| h2o2 | 3000 | 120 | 322 | 161/161 | 0.00 | 0.40 | 1.000 | 1.000 | 0 | 6157 |
| h2o2 | 3000 | 60 | 324 | 162/162 | 0.00 | 0.40 | 1.000 | 1.000 | 0 | 4333 |
| n2h4 | 1000 | 0 | 160 | 80/80 | 0.00 | 2.90 | 0.975 | 1.000 | 0 | 1732 |
| n2h4 | 1000 | 10 | 84 | 42/42 | 0.00 | 1.50 | 0.988 | 1.000 | 1 | 1393 |
| n2h4 | 1000 | 5 | 64 | 32/32 | 0.00 | 1.60 | 0.984 | 1.000 | 0 | 1770 |
| n2h4 | 2000 | 0 | 288 | 144/144 | 0.00 | 1.80 | 0.962 | 1.000 | 0 | 3459 |
| n2h4 | 2000 | 10 | 244 | 122/122 | 0.00 | 1.50 | 0.988 | 1.000 | 1 | 2834 |
| n2h4 | 2000 | 5 | 372 | 186/186 | 0.00 | 2.30 | 0.968 | 1.000 | 0 | 3593 |
| n2h4 | 3000 | 0 | 715 | 358/357 | 0.00 | 13.70 | 0.929 | 0.992 | 0 | 5037 |
| n2h4 | 3000 | 10 | 610 | 305/305 | 0.00 | 11.70 | 0.941 | 0.993 | 1 | 4501 |
| n2h4 | 3000 | 5 | 414 | 206/207 | 0.00 | 10.10 | 0.908 | 0.978 | 0 | 3358 |

- cells with 0 rebuilds in 5 ps: c2h6 1000 K frames 0 and 8 (and the c2h6 1000 K weight cell).
- `n/a` jumps (rebuild whose jump cannot be measured yet): 3, all in `n2h4` frame 10 (one per T).
- no NaN and no instability message in any of the 55 runs; the `Remaining -nan` in the first
  status row of every run is a cosmetic, pre-existing print bug (Epot/Ekin/Etot/T columns are clean).

## 3. Break vs formation, and the temperature scaling

| rebuild kind (bond-count delta) | events | median | p99 | max | <1 kJ | <5 kJ | >=50 kJ |
|---|---|---|---|---|---|---|---|
| BREAK (-1) | 5630 | 0.00 | 69.0 | 423.5 | 0.934 | 0.974 | 58 |
| FORM (+1) | 5617 | 0.00 | 0.4 | 78.7 | 0.996 | 0.998 | 3 |
| no count change (0) | 221 | 0.00 | 0.2 | 27.6 | 0.995 | 0.995 | 0 |

| T | break events | median | max | <5 kJ | >=50 kJ |
|---|---|---|---|---|---|
| 1000 | 314 | 0.00 | 2.9 | 1.000 | 0 |
| 2000 | 1684 | 0.00 | 310.3 | 0.996 | 4 |
| 3000 | 3632 | 0.00 | 423.5 | 0.962 | 54 |

Composition of the 11527 rebuilds: 5630 BREAK, 5617 FORM, 221 no-count-change, 16 multi-bond
(|delta| >= 2) and 43 FIRST rebuilds of a run (median 0.00, max 0.70 -- the case the old binary
printed as NaN/`n/a`; only 3 `n/a` remain, all n2h4 frame 10, i.e. `n2h4` rebuild #1 there).

Dominant term at the large breaks (verbosity-3 reruns of the 3 worst cells, `REACT jump terms`):
the jump is the BOND term (+452.2, +446.9, +366.1 kJ/mol) partly cancelled by `brep` (-38.8, -26.6,
-33.2) and `nbrep` (+15.7, +8.6, +14.1); every other term < 5. `begin_form` transitions are exactly
0.0 in every term except `coul` (max 0.10) -- the well-blend join is energy-neutral as designed.

## 4. New criterion (`order`, default) vs previous (`weight`), matched cells

| system | T | criterion | rebuilds | events | median | max | <1 kJ | <5 kJ |
|---|---|---|---|---|---|---|---|---|
| c2h6 | 1000 | order | 0 | 0 | - | - | - | - |
| c2h6 | 1000 | weight | 0 | 0 | - | - | - | - |
| c2h6 | 2000 | order | 40 | 40 | 0.00 | 0.10 | 1.000 | 1.000 |
| c2h6 | 2000 | weight | 94 | 94 | 0.00 | 144.60 | 0.830 | 0.840 |
| ch3nh2 | 1000 | order | 12 | 12 | 0.00 | 0.70 | 1.000 | 1.000 |
| ch3nh2 | 1000 | weight | 4 | 4 | 34.35 | 68.90 | 0.500 | 0.500 |
| ch3nh2 | 2000 | order | 220 | 220 | 0.00 | 2.80 | 0.991 | 1.000 |
| ch3nh2 | 2000 | weight | 66 | 66 | 0.00 | 329.10 | 0.773 | 0.788 |
| ch3oh | 1000 | order | 84 | 84 | 0.00 | 0.00 | 1.000 | 1.000 |
| ch3oh | 1000 | weight | 4 | 4 | 20.10 | 68.70 | 0.500 | 0.500 |
| ch3oh | 2000 | order | 128 | 128 | 0.00 | 0.50 | 1.000 | 1.000 |
| ch3oh | 2000 | weight | 94 | 94 | 0.00 | 67.20 | 0.926 | 0.947 |
| h2o2 | 1000 | order | 12 | 12 | 0.00 | 0.10 | 1.000 | 1.000 |
| h2o2 | 1000 | weight | 2 | 2 | 40.75 | 41.50 | 0.000 | 0.000 |
| h2o2 | 2000 | order | 196 | 196 | 0.00 | 1.20 | 0.995 | 1.000 |
| h2o2 | 2000 | weight | 94 | 94 | 0.00 | 41.60 | 0.904 | 0.926 |
| n2h4 | 1000 | order | 160 | 160 | 0.00 | 2.90 | 0.975 | 1.000 |
| n2h4 | 1000 | weight | 6 | 6 | 55.65 | 56.80 | 0.333 | 0.333 |
| n2h4 | 2000 | order | 288 | 288 | 0.00 | 1.80 | 0.962 | 1.000 |
| n2h4 | 2000 | weight | 120 | 120 | 0.00 | 57.20 | 0.950 | 0.967 |

- weight pooled over all kinds: n=484, median 0.00, max 329.10
- weight FORMS: n=173, <1 kJ 0.763, >=50 kJ 23, max 146.7  |  order FORMS: n=5617, <1 kJ 0.996, >=50 kJ 3, max 78.7
- weight BREAKS: n=164, <1 kJ 0.963, max 278.5  |  order BREAKS: n=5630, <1 kJ 0.934, max 423.5
- the weight cells are 10 runs vs 45 order cells, so pooled maxima are not comparable; only the
  matched rows above are. The weight max (329.1) sits in a matched cell where order gives <= 0.5 kJ.

## 5. Three worst single events (all cells, order criterion)

| # | system | T (K) | frame | t (fs) | rebuild # | event | dE_jump (kJ/mol) | dominant term |
|---|---|---|---|---|---|---|---|---|
| 1 | c2h6 | 3000 | 8 | 1849.2 | 147 | H4-H5 broken | +423.5 | bond +452.2, brep -38.8, nbrep +15.7 |
| 2 | c2h6 | 3000 | 16 | 2010.8 | 213 | H3-H5 broken | +394.4 | bond +446.9, angle -33.8, brep -26.6 |
| 3 | c2h6 | 2000 | 8 | 2792.8 | 78 | H3-H6 broken | +310.3 | bond +366.1, angle -34.6, brep -33.2 |

- verified real, not a reporting artefact: in `c2h6/T3000_f8` the 250-fs status rows bracket the
  event with Epot -1.004495 Eh (1750 fs) -> -0.875435 Eh (2000 fs), a +0.129 Eh rise containing the
  +0.161 Eh jump; the energy really does step.
- events 1 and 3 belong to a flurry of 9-10 form/break events within ~2 fs (rebuilds #143-#152),
  i.e. many pairs in transition at once; the single worst is NOT a lone event.
- 15 of the 61 events >= 50 kJ/mol are H-H contacts, the rest are C-H (~110 kJ/mol cluster) and C-C.

## 6. Method notes for whoever continues this

- Total wall 59.4 s summed over the 55 runs (max single run 5.89 s); grid wall ~25 s at 16-way
  parallel, whole job ~7 min. More replicates are cheap; the 1000 K cells (0-160 rebuilds, spread
  0/0/4 across frames for c2h6) are the thin ones.
- `-method revgfnff` WITHOUT `-gfnff.topology_mode react` gives ZERO rebuilds in 5 ps -- a first
  pilot was silently useless that way. Always assert the rebuild count is non-zero.
- Runs are deterministic and `-md.seed` is inert. Replicate = new start frame. One directory per
  run, `*.topo.json` deleted (GFN-FF caches the topology next to the input).
- Raw data: `/tmp/claude-1000/.../06e80755-.../scratchpad/polyjump/` -- `runs/<system>/<cell>/
  {input.xyz,cmd.txt,run.log,wall.txt}`, `cells.json`, `analyse.py`, `final_tables.py`, `STATUS.md`.

## 7. Still open / next

- The largest jump is a bond REMOVAL (`begin_break`), so stage 3 must measure break-side
  smoothness, not only the formation join.
- c2h6 dominates (9 of the 12 worst events); h2o2 stays <= 1.7 kJ/mol even at 3000 K.
- 3000 K cells overshoot (T_max up to 5685 K) and are past dissociation -- 1000/2000 K is the
  regime in which the blend is actually being tested.
