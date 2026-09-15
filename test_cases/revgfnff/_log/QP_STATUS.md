# rev-gfnff stage 3a(ii): the apportioned bond share `x_p` (QP) — staged diagnostics (2026-09-14)

AI-generated measurement job, machine-tested. Branch `reactff2-llm`, build dir `build_rev/`.
Pre-change binary md5 `f45bf28a376106c381f29b22eb499548` (= PROXY_STATUS's delivered binary);
binary after the Diagnostic-A switch `087c9efae7d8f8a65ca2bf0cc710cb21`, final binary
(+ the Diagnostic-B dump line) `79bb76bb6720278c9da8b35007ed50ee` (both MAKE_EXIT=0, no new
warning in any touched line). Harness copied from `proxy/jumpP` into
`.../ac8185e8-.../scratchpad/qp/jump/` (`run_one.sh`, `seq.sh`, `analyse.py`, `joblist.txt`).

**Task gate: Stage 1 only. Stage 2 (the QP implementation) was NOT started — Diagnostic A failed
its gate (numbers below).**

## Stage 1 / Diagnostic A — replay of the 3331.6 kJ event + the 1,3 ordinary-join switch

### A.0 Harness calibration (this environment reproduces PROXY_STATUS)

22-cell grid, same `joblist.txt` (c2h6/ch3nh2/ch4_H x 1000/2000 K x 3 frames, 5 ps, dt 0.25,
`-threads 1`, seed 42, sequential). 2 of the 22 cells exit rc=1 with 0 rebuilds in EVERY arm
(`ch4_H/T1000_f16`, `ch4_H/T2000_f16`: the ch4_H trajectory has fewer than 17 frames, so
`extract_frame.py` aborts) — a pre-existing harness defect, identical in all arms, 20 live cells.

| arm | binary | rebuilds | median | max | p99 | <1 kJ | >=50 kJ | T_max |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| PROXY_STATUS "delivered (= default now)" | f45bf28a | 958 | 0.00 | 471.10 | 2.10 | 0.981 | 3 | 8306 |
| this session, default arm | f45bf28a | 960 | 0.00 | **471.10** | 1.80 | 0.981 | **3** | 2037 |
| this session, default arm | 087c9efa | 960 | 0.00 | **471.10** | 1.80 | 0.981 | **3** | 2037 |

max, event count and `<1 kJ` reproduce exactly; the rebuild count is within 0.2 %. (T_max differs
because this `analyse.py` reads the filtered log, where no non-status line can be mis-parsed as a
status row; it is self-consistent across all arms below.) **The new binary reproduces the old one
column for column, i.e. the new PARAM's default path is unchanged.**

### A.1 The 3331.6 kJ event: mechanism confirmed, ATTRIBUTION CORRECTED

Cell `ch3nh2 / T2000 / frame 8`, verbosity 3 (`REACT jump terms` is `CurcumaLogger::info`, and
SimpleMD runs the calculator one level below the run — Known Issue #31 — so verbosity 2 prints
nothing). The trajectory is identical at verbosity 1/2/3 (same rebuild count, same max).

| arm | rebuilds | max \|dE_jump\| | >=50 kJ | begin_* | hard swaps (s >= 0.99) |
|---|---:|---:|---:|---:|---:|
| `-gfnff.rev_share_onethree false` (= default) | 222 | **2.1** | 0 | 111 | **0 / 111** |
| default, no flag at all | 222 | 2.1 | 0 | 111 | 0 / 111 |
| `-gfnff.rev_share_onethree true` (proxy on) | 86 | **3331.6** | 2 | 43 | **1 / 43** |

The event's mechanism is confirmed verbatim in the proxy-ON arm:
`jump terms (begin_form): bond +54.0 angle +446.9 tors -20.3 brep +70.3 nbrep -192.7 coul -2.8
over +2976.2 | s 1.00 w_ws 1.0000 c_ws 0.9867 r_ws 1.3010 ij 2 4` — a `begin_form` admitted at
transition coordinate s = 1.00 (window 0.8 -> 0.9, `c_now` already 0.9867), paid by E_over.

**But Fable's attribution does not reproduce.** FABLE_BOND_STATE section 2.3 states the event
"occurs with the identical value 1.268935 Eh in the proxy-OFF and the share-OFF arms". Here the
proxy-OFF arm of that cell has max 2.1 kJ, 0 events >= 50 kJ and **zero** `begin_*` at s >= 0.99
out of 111. Two independent invocations of the proxy-OFF arm (old and new binary) give identical
numbers, so this is not the invocation-environment sensitivity PROXY_STATUS section 5 warns about.
The "1 / 43" baseline count is reproduced exactly — but it is the **proxy-ON** cell.

**Where the mis-attribution came from, proven**: the log Fable read
(`.../06e80755-.../scratchpad/proxy/jumpP/runs_dbg/ch3nh2/T2000_f8/run.log`, whose `cmd.txt`
carries no `rev_share_onethree` flag and was therefore taken to be the default) has its 86
`REACT jump terms` lines **identical line for line** to my `-gfnff.rev_share_onethree true`
replay, and shares not a single line with the default replay (which has 111 such events). That
log is a proxy-ON run mislabelled by its `cmd.txt`. Consequence: the 3331.6 kJ headline event IS
caused by the proxy, and the delivered default has no 1,3-closure hard swap at all in this cell.

### A.2 The companion change: `-gfnff.rev_bo13_ordinary_join` (new, DEFAULT OFF)

Implemented per FABLE_BOND_STATE 3.6, second bullet, as its own switch, independent of anything in
Stage 2. `gfnff_method.cpp` scan: the "shares a SETTLED neighbour" bit is still computed (and still
printed by the scan trace) but no longer selects the formation criterion or the transition window,
so a 1,3 pair joins on the ordinary `rev_bo2_form` criterion with the ordinary window.
Files: `gfnff_method.cpp` (scan `one_three_raw`/`one_three`, `setupRevSettings` read, flags dump),
`gfnff.h` (PARAM + member). 3 files, 1 new PARAM, no other behaviour touched.

**Result — default (proxy-off) arm, 22 cells: the switch is a NO-OP, every column identical.**

| arm (delivered default + …) | rebuilds | median | max | p99 | <1 kJ | >=50 kJ | T_max | hard swaps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 960 | 0.00 | 471.10 | 1.80 | 0.9812 | 3 | 2037 | **5 / 478** |
| `rev_bo13_ordinary_join true` | 960 | 0.00 | 471.10 | 1.80 | 0.9812 | 3 | 2037 | **5 / 478** |

**Result — proxy-ON arm, 22 cells (the arm in which 1,3 formations actually occur):**

| arm (proxy on + …) | rebuilds | median | max | p99 | <1 kJ | >=50 kJ | T_max | hard swaps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| baseline | 414 | 0.00 | 3331.6 | 470.9 | 0.8641 | 24 | 2037 | **19 / 205** |
| `rev_bo13_ordinary_join true` | 434 | 0.00 | **3146.0** | 414.6 | 0.8773 | 23 | 2037 | **17 / 215** |

and on the single debug cell `ch3nh2/T2000_f8` with the proxy on: 86 reb / max 3331.6 / 2 events
/ 1 of 43 hard -> 80 reb / max **401.1** / 3 events / **2 of 40** hard. The +2976 kJ `begin_form`
at s = 1.00 is gone from that cell; the two remaining hard swaps there are `begin_break`.

**What the switch does NOT do**: the delivered default is unaffected (0 % change on 22 cells), and
in the proxy arm the hard-swap population falls only 19/205 -> 17/215 (9.3 % -> 7.9 % of `begin_*`)
while the grid maximum stays in the same class (3331.6 -> 3146.0 kJ) and the >=50 kJ count falls
24 -> 23. Nothing approaches the pre-valence-share level (10.7 kJ, 0 events >= 50).

**Why the default arm cannot move** (source, `gfnff_method.cpp` scan): with `rev_form_switch =
order` the 1,3 branch and the ordinary branch are already numerically identical — `rev_bo13_form`
= `rev_bo2_form` = 0.1 and both windows end at 0.9 (`tr.tight` is true either way). The
classification therefore changes behaviour for exactly ONE case: a **fading** 1,3 pair (a
recently-broken bond whose well is still carried), where `one_three` overrides the fading branch
(`o_ij > 0.1` on the bo2 switch, window 0.1..0.9 with `tight = true`) instead of `c_ij >
rev_tr_begin = 0.02` on the bo3 switch, window 0.02..0.8. That case occurs in the proxy-ON
trajectories and does not occur anywhere in the 22 default cells.

**Hard swaps in the delivered default are breaks, not 1,3 closures.** All 5 of the 478 `begin_*`
events at s >= 0.99 are `begin_break`, and only one of them is energetically significant
(`c2h6/T2000_f16`, `bond +451.0 angle +28.9 brep -12.1 nbrep +3.7 | s 1.00 w_ws 0.0206`) — it IS
the 471.1 kJ grid maximum. The other four (all `ch4_H/T2000_f10`) are +276, +56, -4 and +0 kJ.
The whole >= 50 kJ tail of the default arm is three events: `c2h6/T2000_f16` +471.1,
`ch4_H/T2000_f10` +310.4 and +51.6 kJ/mol — bond wells dumped at a forced break, the class
VALFIX section 6 already identified, not the E_over-paying ring closure that Fable's 3.6 targets.
Neither the scan switch measured here nor the QP of section 3 addresses that class.

### A.3 Diagnostic A verdict

**FAILED.** Required: "a real reduction in hard-swap population/smoothness cost from the
ordinary-join switch". Measured: exactly zero change on the delivered default (every one of the 8
reported columns identical over 22 cells), and 19 -> 17 hard swaps / 24 -> 23 events >= 50 kJ /
max 3331.6 -> 3146.0 kJ in the proxy arm. The switch is delivered anyway, DEFAULT OFF, because it
is a faithful implementation of the proposed companion change and its null result is the finding.

## Stage 1 / Diagnostic B — offline QP over the rkt06_h_h2 path, no src physics change

### B.0 The one src addition this needed, and why

`CURCUMA_SHAREDUMP=1` did **not** print `D_p`. Its table (`prepareValenceShare`) has w, tight b,
sig, g, sum, Val, u, f, c — the well depth is formed later, in `calcBonds`, because the dynamic r0
(and hence the Gaussian) is built there. One env-gated line was added to `calcBonds` under the
SAME `CURCUMA_SHAREDUMP` gate (`ff_workspace_gfnff.cpp`, `s_share_dump`):

    shareD  <idx> <i>-<j> r <r> D <-k_b e^{-a dr^2} w> w <w> c <c> E <energy>

Diagnostic check on every path point: `sum_p E_p` equals the reported `Bond` component to the
dump's 8 decimals. No physics touched; the print is off unless the variable is set.

### B.1 The solver

`.../scratchpad/qp/pathqp.py`, standalone, no external QP dependency. Per point: dual Gauss-Seidel
over the atoms, each atom a bisection on the monotone
`h_i(lambda_i) = sum_{p at i} x_p(lambda_i + lambda_{j_p}) - B_i`, with the closed form (3.5)
`x_p = clip(1/2 + 1/beta - Lambda_p/(beta D_p), 0, 1)`. Max constraint violation over all points
and all beta: **<= 2.8e-14**. The Bond term is rebuilt as `F(x*) = sum_p D_p[-x_p -
(beta/2) x_p(1-x_p)]` and substituted for the delivered Bond term in each point's total, so the
well shape `D_p(r)` is held exactly at what the delivered code computes.

### B.2 Result — 11 points, kcal/mol relative to the reactant (r2SCAN-3c reference)

Harness calibration: the delivered path reproduces **rms = 2.71** (CIJ_STATUS's value) exactly.

| beta | 0 (delivered) | 0.01 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 | 1.0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| path rms | **2.71** | 3.14 | 3.19 | **3.38** | **3.88** | 4.98 | 8.38 | 14.12 |

| pt | ref | delivered | b=0.05 | b=0.1 | b=0.2 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 1 | -0.05 | -0.04 | -0.04 | -0.04 | -0.04 |
| 2 | -0.18 | -0.14 | -0.14 | -0.14 | -0.14 |
| 3 | -0.35 | -0.07 | **-0.07** | **-0.07** | **-0.07** |
| 4 | 0.86 | 3.41 | **3.41** | **3.41** | **3.41** |
| 5 | 2.53 | -3.06 | -5.03 | -6.07 | -8.62 |
| 6 | 0.51 | 2.13 | **2.13** | **2.13** | **2.13** |
| 7 | -0.32 | -0.17 | **-0.17** | **-0.17** | **-0.17** |
| 8 | -0.06 | -0.06 | -0.06 | -0.06 | -0.06 |
| 9 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 10 | 2.57 | -3.79 | -5.15 | -6.50 | -9.21 |

**Points 3-4 and 6-7 are bit-identical to the delivered path at every beta** — and the reason is
structural, not a coincidence of the crossover width: at those geometries the perception carries
exactly ONE bond, so the budget is slack, `x = 1` by the box bound, and the QP's candidate test
(3.4 step 2) never even runs. The QP can only act at points 5 and 10, the two points with two live
pairs. So "no visible misfit at 3-4/6-7" is satisfied trivially and carries no information about
the QP; it also means the QP cannot repair the delivered path's actual worst points (pt 4 +2.55,
pt 6 +1.62 kcal off the reference).

Per-pair shares at the two points that do move:

| pt | r (Bohr) | D_p (Eh) | delivered c | x(0.05) | x(0.1) | x(0.2) |
|---|---|---|---|---|---|---|
| 5 | 1.851 / 1.687 | 0.16928 / 0.17512 | 0.500 / 0.500 | 0.161 / 0.839 | 0.331 / 0.669 | 0.415 / 0.585 |
| 10 | 1.762 / 1.762 | 0.17271 / 0.17272 | 0.500 / 0.500 | 0.500 / 0.500 | 0.500 / 0.500 | 0.500 / 0.500 |

Fable's 4.2 per-point prediction is confirmed quantitatively: at pt 10 the symmetric split gives
`E = -D(1 + beta/4)`, i.e. **-2.71 kcal/mol at beta = 0.1** (measured -3.79 -> -6.50 = -2.71). His
*estimated* rms of 5-6 was too pessimistic (the measured in-H3 depth ratio at pt 5 is 1.0345, not
the 1.075 he read off the lone-H2 wells), so the real numbers are 3.38 / 3.88.

### B.3 Diagnostic B verdict

**Literal threshold met, substance negative.** The stated bar — "rms <= ~4 kcal/mol at at least
one beta with no visible misfit at points 3-4/6-7" — is met at beta = 0.01…0.1 (3.14…3.88), and
points 3-4/6-7 are untouched. But the QP is **worse than the delivered share at every beta tested**
(2.71 -> 3.14 at the beta -> 0 limit, monotonically rising to 14.12 at beta = 1), and both transition
states move further below the reference, i.e. in the wrong direction. There is no beta at which the
QP improves this falsifier; `beta` is a one-sided knob whose best value on rkt06 is the limit in
which the sharing term vanishes.

### B.4 Bonus (offline, same script): the compressed BF4- per-pair verdict

Inputs measured from the dump at B-F = 1.143 A, charge -1 (10-bond corner): `D_BF = 0.07696004`,
`D_FF = 0.12565950` Eh, `B_B = 4.0000`, `B_F = 1.0139` — identical to FABLE 4.1.

| beta | x_BF | x_FF | E_bond | E_tot | vs 10-bond share-off (-0.78960808) | vs pinned 4-bond (-1.29638269) |
|---|---:|---:|---:|---:|---:|---:|
| 0.05 / 0.1 / 0.2 | **1.000000** | **0.004633** | -0.31142 / -0.31151 / -0.31168 | -0.03923 / -0.03932 / -0.03949 | **+470.8** kcal/mol | **+788.9** kcal/mol |
| 0.25 | 0.866 | 0.049 | -0.31264 | -0.04045 | +470.2 | +788.1 |
| 0.30 | 0.770 | 0.081 | -0.31496 | -0.04277 | +468.8 | +786.7 |

The operator's restated per-pair criterion (i) is met exactly for beta <= 0.202, as Fable derived:
B-F at full share, F...F claiming essentially nothing. The residual `x_FF = 0.004633` is not noise
— it is the fluorine's own budget slack `(B_F - 1)/3 = 0.0139/3` distributed over its three F...F
pairs, i.e. the softplus floor of the settled-count valence, and it would vanish with a hard
`B_F = 1`. Criterion (ii): the measured gap is **+470.8 kcal/mol against the 10-bond share-off
reference (Fable predicted ~+473) and +788.9 against the pinned 4-bond evaluation (Fable predicted
~+790)**. Both predictions are confirmed; the deviation from Fable is <3 kcal in each case, i.e.
nowhere near the "large deviation from the prediction" that would have been the red flag.

## Gate

**NOT CLEARED.** Diagnostic A failed outright (no reduction). Diagnostic B met its literal rms bar
but the QP is worse than the delivered share at every beta on the falsifier it was supposed to
defend. Stage 2 (the `-gfnff.rev_valence_share qp` implementation) was therefore **not started**.

## What is left in the working tree (uncommitted)

1. `gfnff.h` / `gfnff_method.cpp`: PARAM `rev_bo13_ordinary_join` (Bool, **default false**) + its
   member, read, flags-dump entry and the scan use; the 1,3 classification is still computed and
   still printed by the scan trace.
2. `ff_workspace_gfnff.cpp`: one `CURCUMA_SHAREDUMP`-gated `shareD` line in `calcBonds` (well depth
   per pair) + the `s_share_dump` static that gates it.

Nothing else. Final binary md5 `79bb76bb6720278c9da8b35007ed50ee`.

### Default-path verification of those two changes

- 22-cell react-MD grid: the new binary reproduces the pre-change binary column for column
  (960 rebuilds, max 471.10, p99 1.80, `<1 kJ` 0.9812, 3 events >= 50, 5/478 hard swaps).
- 12-digit single points (`-batch`, gradient on): revgfnff caffeine **-4.673521653477**, benzene
  **-2.363224128930**; gfnff caffeine **-4.672737068614**, benzene **-2.362725526194** — all four
  equal to the CIJ_STATUS / PROXY_STATUS records.
- `ctest -R gfnff` in `build_rev/`: **62/65, 87.12 s**. `ctest -R "cli_simplemd_"`: **19/22,
  84.63 s**. The same three failures in both: `cli_simplemd_16_gfnff_rev_h2_smooth`,
  `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`, `cli_simplemd_19_gfnff_rev_form_refuses_hbond` — the
  known pre-existing set the task names (MD-behaviour calibration predating the stage 3a changes;
  PROXY_STATUS's 65/65 was measured at commit `11b1baea`, before them). None of the three sets the
  new PARAM, and the default path is bit-identical as shown above. No test file and no
  `test_cases/cli/CMakeLists.txt` touched.

