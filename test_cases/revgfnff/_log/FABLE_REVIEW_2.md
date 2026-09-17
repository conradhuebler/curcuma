# FABLE_REVIEW_2 — independent method review of the rev-gfnff 3a(ii) state (2026-09-17)
Sections written: 9/9

AI-generated review, read-only on the repository. All own measurements: frozen binary
`scratchpad/fable2/curcuma_frozen`, md5 `f7a37866cf88f9054060c86721616688` (= HEAD `4d5287bc`),
`-threads 1`, fresh directory per structure, no `*.topo.json` reuse. Harness and raw tables in
`scratchpad/fable2/` (`run_md.sh`, `wfilt.py`, `perstep.py`, ...). Labels: **MEASURED** (own run),
**INFERRED**, **PROPOSED**.

## A.1 The H-budget mechanism and the per-step numbers — reproduced (MEASURED)

Four MD runs (2 cells x flag off/on), 20 001 status rows each, `CURCUMA_SHAREDUMP=1` piped through a
window filter. Per-step max |dEpot| outside rebuild intervals:

| cell | arm | max / kJ | 2nd | n >= 50 kJ | T_max / K | rebuild-line fingerprint |
|---|---|---:|---:|---:|---:|---:|
| c2h6/T2000_f16 | off | 2593.60 @3.20825 ps | 2088.36 | 22 | 62 098 | 135 |
| c2h6/T2000_f16 | on | 59.26 @4.02625 | 57.91 | 17 | 5 086 | 204 |
| ch4_H/T2000_f10 | off | 1088.96 @1.42425 | 1076.55 | 180 | 1.594e8 | 27 |
| ch4_H/T2000_f10 | on | 391.44 @0.00025, then 222.53 @0.02125 | 192.55 @0.01725 | 398 | 8 306 | 126 |

Every number of HBUDGET_STATUS section 2 reproduces to the last digit on an independent binary
copy and harness. The briefing's links (3) and (4) stand as far as the hydrogen is concerned.

**But the briefing's reading of the residual is wrong, and this is the main finding of part A**
(detail in A.3): the 391 kJ first step and the 222/193 kJ steps of ch4_H/T2000_f10 all sit in the
first 21 fs of the run and are the SAME defect on the carbon, not thermal amplitude.

## A.2 "No falsifier moves with `rev_budget_fix_h true`" — spot-checked, holds (MEASURED)

- rkt06 (`ref/B/rkt06_h_h2`, 11 points, fresh SP each): rms **2.714** off and on, barrier +3.41 @pt4
  (ref +2.57 @pt10), max |E_on - E_off| = **0.000e+00 Eh** over the 11 points.
- NH4+ / H3O+ / CH5+ / ClO4- / BF4- 1.394 A / BF4- compressed: off = on to all 8 printed digits
  (0.82730341 / 1.20721990 / 0.78269212 / 0.03959777 / -1.47018405 / 0.11826903 Eh, n = 6).
- Not in the project's falsifier list but the place an H budget WOULD matter, so I ran them:
  class-E proton transfers `nh4_nh3_pt` (13 pts), `nh3_2_transit`, `h2o2_transit` (11 each):
  off = on at every point; `hf2_transit` identical on 10/11 points (last point 739.1 vs 779.1
  kcal/mol, a geometry that is already +700 off the reference); class-C `h2_Hp` (H3+) -65.58 vs
  -65.56 kcal/mol. **Verdict: make `rev_budget_fix_h` the default.** No measured cost in 9 species
  / 5 paths; it removes both runaways.
- Side finding on those paths (not the share's, stage 2's): `nh4_nh3_pt` is a mirror-symmetric
  scan (ref 0.0 / 0.0 at both ends) and the model gives 0.0 / **+43.8** kcal/mol, with +146.0 at the
  midpoint (ref -1.7) in all three arms incl. share off. Consistent with the "whole charge on
  fragment 0" rule (Known Issue #13) flipping which fragment is the cation. Any proton-transfer
  acceptance has to wait for stage 2's localised fragment charge.

## A.3 The residual after the H fix is the SAME defect on carbon, not thermal amplitude (MEASURED)

The grid has **20 cells, not 22**: `ch4_H.xyz` has 15 frames, so `ch4_H/T1000_f16` and
`T2000_f16` exit with rc = 1 and an empty log — in my runs and in `hb/runs_ON` (`wall.txt`:
`rc=1 wall=0.01`). Every "22-cell" statement since CIJ_STATUS is n = 20.

Per-step events >= 50 kJ/mol over the 20 cells, flag on, with a compact per-step share record
(`cfilt.py`: Val per atom, min c, first-corner bond energy): **487 events** (HBUDGET: 485), of which

| where | events | max / kJ | max Val (any atom) | events with a same-step budget change >= 0.05 |
|---|---:|---:|---:|---:|
| ch4_H/T2000_f10 + ch4_H/T1000_f10 | **438** (398 + 40) | 391.4 / 388.4 | **5.00** (carbon) | **387** |
| the other 18 cells | 49 (c2h6 33, ch3nh2 16, all T = 2000 K) | 71.4 | 4.01 | **0** |

`f10` is the class-C frame d(C-H) = 1.40 A: the sixth atom sits in a face of CH4, carbon has five
partners. Share dump of the first three steps (flag on):

| t / fs | r(C-H6) / a0 | sig(C-H6) | Val(C) | c (all five C-H) | sum of the five wells / Eh |
|---|---:|---:|---:|---:|---:|
| 0.00 | 2.646 | 0.534 | 4.531 | 0.773 | -0.6083 |
| 0.25 | 2.527 | 0.798 | 4.794 | 0.945 | -0.7521 |
| 0.50 | 2.043 | 0.9995 | 4.990 | 1.000 | — |

-0.1438 Eh of the -0.1491 Eh (391.4 kJ) first-step jump is the bond term through c: the carbon's
budget crosses the settled window, ALL FIVE wells go from 0.77 to 1.00 together, H6 is pulled in by
0.6 a0 in two steps and then rattles between 1.58 and 2.6 a0 in a fifth full C-H well; every
re-crossing of the window modulates five wells in phase (the 222.5 and 192.6 kJ steps at 21.25 and
17.25 fs). This is link (3) of the briefing with C in place of H. HBUDGET section 6's "within the
thermal amplitude" is wrong for these two cells. For the other 18 cells it is right: their 49
events have no budget change, Val <= 4.01, and their per-step |dEpot + dEkin| (median 6.8 / 8.8
kJ) sits inside the all-step distribution (p99 20.9 / 19.4) — ordinary Epot<->Ekin exchange of an
8-atom molecule at 2000-5000 K and dt = 0.25 fs. INFERRED: the 50 kJ threshold is at the thermal
ceiling there; scale the threshold with 3N k T or use |dEtot|.

## A.4 Is the carbon budget the class-B root ("H + CH4 TS 78 kcal/mol too low")? No. (MEASURED)

1. The "78" is not a measurement of any delivered model: it is FABLE_ROADMAP_REVIEW lines 153-154,
   a pre-share estimate from two isolated pair curves. An anchor.
2. `ref/P/rkt03` IS H + CH4 -> CH3 + H2 (14 r2SCAN-3c points). Delivered model, relative to
   CH4 + H: TS **-6.0** (ref +7.6), CH3 + H2 -26.8 (ref -1.6), path rms 15.3 (12.5 against the
   file's own zero). Val(C) = **4.00-4.01 at all 14 points** (collinear abstraction: carbon never
   has a fifth partner), so re-evaluating the path offline with the carbon held at 4 gives the
   identical curve to 0.00 kcal. The rkt03 error is the -25 kcal/mol reaction energy, i.e. the
   C-H vs H-H well depths (class A / 3a(iii)), not the budget.
3. Where the carbon budget does decide: the class-C approach scans (rigid X...H, 15 points each).
   kcal/mol relative to d = 3.0 A, single-corner points only (exact offline re-evaluation; the
   offline "delivered" reproduces the binary to <= 2e-5 Eh):

| scan, d / A | ref | delivered | same formula, Val = nominal | proportional form, nominal |
|---|---:|---:|---:|---:|
| CH4 + H 1.0 / 1.2 / 1.3 | +97.9 / +60.4 / +47.3 | **+10.9 / -17.6 / -21.2** | +269.8 / +237.8 / +221.7 | +114.5 / +84.4 / +70.6 |
| NH3 + H 1.0 / 1.2 / 1.3 | +50.5 / +35.4 / +27.6 | **-56.5 / -56.5 / -30.5** | +161.5 / +153.7 / +143.2 | +52.5 / +48.1 / +40.5 |
| H2O + H 1.0 / 1.2 | +42.5 / +29.8 | **-11.7 / -23.9** | +118.2 / +95.6 | +74.9 / +53.7 |
| N2H4 + H 1.0 / 1.2 / 1.3 | +34.2 / +24.2 / +19.6 | **-56.2 / -58.1 / -32.9** | +149.7 / +139.9 / +129.8 | +46.8 / +40.4 / +33.2 |
| NH3 + H+ 1.0 / 1.2 | -146.5 / -135.8 | -87.0 / -87.1 | +131.0 / +123.1 | +22.0 (nominal) / **-86.9 with Val + q** |

   and CH4 + H delivered has +69.1 at 1.40 A next to -21.2 at 1.30 A: a 90 kcal/mol drop over 0.1 A
   = the budget snap of A.3 seen statically. So: **a radical H that hits a saturated C, N or O
   finds an artificial adduct 54-107 kcal/mol below the reference** (n = 4 systems; HF + H stays
   repulsive). That is the real class-B/C consequence of the growing budget, and simply capping
   the element (column 3) over-corrects by +100..+170 kcal/mol. Reason, in the formula:
   `f_i = clip((Val_i - sum_{k != j} w_ik)/w_ij)` is a LEFT-OVER rule and every listed pair has
   w ~ 1, so one partner too many makes every pair of that atom see "nothing left": five
   partners on Val 4 get f_C = 0 on all five, the atom hands out 0 of its 4 valences. The share
   does not conserve valence; the growing budget exists to paper over exactly that.

## A.5 Recommendation — which elements may grow a budget, and by what rule (PROPOSED)

One rule, two parts, both needed.

**(1) Replace the left-over share by a valence-conserving one**: `f_i = min(1, Val_i / S_i)`,
`S_i = sum_k w_ik` (smooth min), `c_ij = f_i f_j`. Then `sum_j f_i w_ij = min(Val_i, S_i)` exactly:
an over-coordinated atom spreads what it has instead of forfeiting it. Same inputs, same chain
rule structure (Lambda pass over S_i), no new parameter. Measured offline (column 5 above): the
four radical adducts turn from -54..-107 to **+2..+32** kcal/mol of the reference (n = 4 systems,
11 single-corner points); rkt06 rms 2.71 -> **2.8** (bridging H: f = 1/2, c = 1/2 x 1, as now);
rkt03 unchanged (0.00); equilibria c = 1 exactly as now (S_i <= Val_i). The product instead of
the mean matters: with the mean the rkt06 TS pair gets c = 0.75 and the path breaks.

**(2) Excess budget is granted by CHARGE, not by element.** What separates NH4+ (4 full bonds)
from NH3 + H (no bond; ref +50) is not the element and not the geometry — the two are identical
in both — it is the electron count, and the only electron-count information a force field has is
the charge. So: `Val_i = Val_Z + min(G(N_i - Val_Z), X_i)` with
- H, F: X = 0 always (H = the delivered fix; F is never hypervalent).
- group 13 (B, Al): X = 1 unconditionally (empty orbital: BF4-, BH4-, AlCl4-, H3N-BH3).
- period >= 3, groups 15-17: X = 6 - Val_Z as the table already does for P and S; extend to
  Cl/Br/I (ClO4-, ClF3). Costs nothing on H + HCl (Cl has one partner along rkt01).
- everything else incl. C, N, O: `X_i = clip(Q_i)`, Q_i = positive topological (phase-1 EEQ)
  charge of atom i plus its settled H partners, per corner — a constant inside a corner, carried
  across swaps by the existing s-blend, no chain rule. Isolated NH4+, H3O+, CH5+: Q = +1 exactly,
  so c = 1 and all three falsifiers stay bit-identical; NH3 + H, CH4 + H, H2O + H, N2H4 + H:
  Q = 0, nominal valence, column 5.

Falsifier consequences, offline from the binary's own share dump (E(share) - E(share off),
kcal/mol; delivered / proposed (1)+(2) / (1) with nominal valences): NH4+ 0.00 / **+0.17** /
+108.3; H3O+ 0.00 / **+0.05** / +86.6; ClO4- 0.00 / **0.00** / +209.0; BF4- 1.394 A 0.00 /
**0.00** / +49.2; BF4- 1.143 A +569.7 / +588.4 / +600.4 (unchanged in kind: a perception
question, FABLE_BOND_STATE's structural argument stands); CH5+ not evaluable offline (two
corners, weights not dumped). The +0.17 / +0.05 are the settled count's sig = 0.9993 < 1 showing
through `Val/S`; build the excess from the same wide w as S (`G(S_i - Val_Z)`) and they are 0
exactly. rkt06 2.71 -> 2.8, rkt03 0.00 change. The grid: the two f10 cells lose the CH5 well
(INFERRED from the static curve, +71 instead of -21 kcal/mol at 1.3 A; to be measured after a
build). Part (1) and part (2) are both needed: (1) alone fails every hypervalent ion (last
column), (2) alone on the left-over formula is column 4 of the A.4 table.
**Not covered and must be measured before adoption**: neutral dative N (H3N-BH3, amine oxides,
ylides) where Q_i < 1 — they keep c = 1 under the delivered rule and would lose up to
(1 - Val/4) of four wells here; H5O2+ / N2H7+ where the charge is split over two groups. If the
operator wants the minimal step first: (1) alone plus X = 0 for H, C, F and the delivered
growth for N, O, B, period >= 3 — it fixes CH4 + H and the f10 cells and leaves NH3 + H / H2O + H
adducts (-57 / -24) as a known electron-count limit.

## A.6 The still-discrete perception rules inside a reactive FF (judgement)

BREAK_TAIL's numbers are not re-measured here. INFERRED from them: the x1.635 neighbour force
constant is the product of three rules that are all keyed on a TRANSIENT bond (2-coordinate H ->
"sp", C-H-H -> 3-ring `ringf` 1.18, `fxh` 1.05). With the H budget fixed the transient H-H that
triggers them no longer collapses (hard swaps 5 -> 0 on 20 cells), so the symptom is gone, the
rules are not. They are not acceptable as they stand: a hybridisation or ring flag that a pair in
flight can set re-parametrises wells it does not belong to, and the s-blend only hides that while
the transition is slow. Minimal treatment, no new physics: in rev mode derive `hyb` and ring
membership from SETTLED bonds only (sigma_p above ~0.5, or pairs not in flight), so that a
transient contact cannot change its neighbours' `bsmat`/`ringf`/`fxh`; `rev_h_not_sp` is already
the first instance of that policy, generalise it rather than adding per-rule exemptions.

## B. Can the bond wells be realised through erf functions? (operator's question)

### B.1 What the family can and cannot do (analysis)

`d/dr erf(sqrt(alpha)(r - r0))` IS the GFN-FF Gaussian well: the delivered well is already the
derivative of an erf switch centred at r0. The stage-1 bond order `b = 1/2 erfc(|k|(r - R)/R)` has a
Gaussian as ITS derivative. Same family, one integration apart. Candidates, against the project's
guard (r_min = r0 from the existing dynamic-r0 pipeline, curvature 2 alpha |k_b|, depth from k_b):

| candidate | E'(r0) = 0 | curvature reproducible | new params / bond type | verdict |
|---|---|---|---|---|
| C1 `-(D/2) erfc((r - r_c)/sigma)` vs the existing repulsion | **never** (E' > 0 everywhere) | only together with a refit repulsion | 3 | r_eq becomes a balance of two terms; every `rabshift`/CN/pi shift of GFN-FF loses its meaning. **Reject** |
| C2 two-erf `-(A/2)[erf(a(x - x1)) - erf(b(x - x2))]` | by solving x2 (needs b > a e^{-a^2 x1^2}) | yes, one more equation | 4 free after the constraint | flat-bottom box when x2 - x1 >> 1/a (k -> 0 at r0); measured no better than the 3-parameter forms (all-point fit, n = 32: median rms 6.2-6.4 vs 5.8-7.0). **Reject** |
| C3 erf envelope on the Gaussian | yes | yes | 2 | is the delivered `G w`: an envelope can only cut, never add the missing tail or depth. **Already there, not a fix** |
| **C4 erf-Morse** `E = -D (2y - y^2)`, `y = erfc((r - R)/sigma) / erfc((r0 - R)/sigma)` | **identically**, for any R, sigma (y = 1) | `E'' = 2 D h(z0)^2 / sigma^2`, h = (2/sqrt(pi)) e^{-z0^2}/erfc(z0), z0 = (r0 - R)/sigma; solve z0 by bisection once per bond type | D, R, sigma (3) or D, sigma with the curvature pinned (2) | **the one viable erf form**, tested below |

C4's properties: depth exactly D at r0; y depends on r - r0 only (R = r0 + u with constant u), so
the CN chain rule is the Gaussian's unchanged (`dE/dr0 = -dE/dr`); gradient
`dE/dr = -2D(1 - y) y'`, `y' = -(2/(sqrt(pi) sigma)) e^{-z^2}/erfc(z0)` — one exp and one erfc,
the vendored portable `curcuma_erf` applies; inner side bounded (y <= 1/erfc-ratio, so
E >= -D(2y_max - y_max^2)): the wall stays the repulsion's job, as now; tail ~ e^{-z^2}/z, faster
than Morse — the direction the roadmap review measured the reference to need ("r90 too far by
+0.1..+0.4 for Morse in 27/28"). Expanding `-ln y = a x + beta x^2 + ...` gives `a = h(z0)/sigma`,
`beta = h(z0)(h(z0) - 2 z0)/(2 sigma^2) > 0` (h(z) > 2z for every z, the Mills-ratio inequality): **C4 is the planned MG form
(`phi = a x + beta x^2`) with beta >= 0 built in** — a re-parametrisation, not a new shape.

### B.2 Offline test against all 32 class-A curves (MEASURED; binary f7a37866, `wellscan.py`/`wellfitB.py`)

Protocol: `revgfnff`, static topology kept from the r_eq frame (no react, so the pair stays in the
list and there is one corner), per-frame term decomposition from the batch JSONL and the pair's own
`shareD` line (D, w, c). `E_rest = E_total - E_pair` frozen; candidate = `E_rest + E_cand(r - r0)`
with r0 the delivered dynamic r0 (recovered from ln(D/w), residual <= 6e-3). Objective and rms on
the BREAK side (r >= r_eq - one grid point), relative to the value at the reference minimum;
metrics by the harness definitions. "q-frozen" = the Coulomb term held at its r_eq value.
NB this protocol's delivered numbers (median dD_e +6.5, dr90 -0.124) are not the react
yardstick's (-25.3 / -0.318; I reproduced ch4_C-H there: rms 5.87, D_e dev -13.32, r90 dev
-0.351): in react mode the bond drop cuts D_e, here it does not.

| n = 32, medians | params | break rms (raw rest) | break rms (q-frozen) | abs dD_e | abs dr90 | n(rms < 3) |
|---|---:|---:|---:|---:|---:|---:|
| delivered Gaussian x w | 0 | **19.24** | — | 6.5 (signed +6.5) | 0.124 (signed -0.124) | — |
| Gaussian, D refit | 1 | 5.71 | 4.83 | 3.2 | 0.145 | 7 |
| MG, free (D, a, beta) | 3 | 1.67 | **1.32** | 1.0 | 0.024 | 25 |
| erf-Morse C4, free (D, R, sigma) | 3 | 1.70 | **1.36** | 1.0 | 0.024 | 25 |
| MG, curvature pinned to 2 alpha k_b (D, beta) | 2 | 2.10 | **2.07** | 0.9 | 0.065 | 19 |
| erf-Morse, curvature pinned (D, sigma) | 2 | 2.17 | **2.13** | 1.0 | 0.088 | 19 |

Per bond (q-frozen; ref D_e / r90, then delivered, then erf-Morse pinned | free, D_e / r90 / rms):

| bond | ref | delivered (rms) | erf-Morse pinned | erf-Morse free | MG free |
|---|---|---|---|---|---|
| ch4 C-H | 115.2 / 2.224 | 103.8 / 1.889 (5.52) | 114.3 / 2.198 / 0.69 | 114.7 / 2.219 / 0.56 | 114.8 / 2.223 / 0.58 |
| h2 H-H | 107.0 / 2.633 | 104.9 / 1.703 (13.04) | 106.9 / 2.225 / 4.33 | 109.7 / 2.567 / 3.11 | 109.7 / 2.570 / 3.11 |
| c2h6 C-C | 108.2 / 1.877 | 90.7 / 1.688 (10.03) | 107.9 / 1.870 / 0.92 | 108.3 / 1.901 / 0.70 | 108.3 / 1.904 / 0.72 |
| ch3oh C-O | 98.0 / 1.723 | 99.9 / 1.636 (6.59) | 98.1 / 1.776 / 2.22 | 98.0 / 1.774 / 2.22 | 98.1 / 1.778 / 2.25 |
| c2h4 C=C | 188.5 / 2.099 | 201.7 / 1.770 (19.91) | 187.8 / 1.988 / 4.65 | 191.2 / 2.132 / 3.01 | 191.2 / 2.133 / 3.01 |
| n2 N#N | 216.4 / 1.787 | 278.1 / 1.909 (39.06) | 211.2 / 1.680 / 4.61 | 214.2 / 1.758 / 1.69 | 214.3 / 1.758 / 1.62 |
| h2o O-H | 121.5 / 2.082 | 105.0 / 1.844 (10.27) | 123.3 / 1.928 / 3.90 | 124.1 / 2.054 / 2.80 | 124.3 / 2.060 / 2.83 |
| hcl H-Cl | 104.4 / 1.878 | 59.4 / 1.694 (30.03) | 104.6 / 1.782 / 2.55 | 105.4 / 1.931 / 1.84 | 105.5 / 1.937 / 1.91 |

- **erf-Morse == MG on every one of the 32 curves**: max |rms difference| 0.12 (q-frozen), 0.24
  (raw); 0 bonds where either wins by > 0.5. The data cannot tell them apart (B.1 says why).
- The curvature-pinned two-parameter variant keeps r_min = r0 AND E''(r0) = the Gaussian's, i.e.
  it changes the equilibrium region only at third order plus a constant depth offset per bond,
  and still takes the break-side rms 19.2 -> 2.1, abs dD_e to 0.9, abs dr90 to 0.07-0.09.
  Freeing the curvature buys 2.1 -> 1.3 and costs the guard (fitted K_pair / K_Gauss: C-H 1.07-1.10,
  C-C 1.16-1.19, O-H 1.83, H-Cl 1.60-1.66, C=C 1.58, N#N 1.96-2.03, H-H 1.97).
- Depth scale s = D/|k_b| (pinned, q-frozen): median 1.02, range 0.46-2.36; X-H (n = 8) median
  1.21, heavy-heavy (n = 24) 0.91. One global s does not exist; this is the 3b table.
- Worst after the fit (pinned): h2o2 O-O 12.8, cl2 9.5, h2co C=O 9.4, c2h2 8.5, hcn C#N 8.1 —
  the QUALITY-flagged / state-unstable references and the multiple bonds; H-H stays at 3.1-4.3
  with a +0.05 A r_min shift and K x 1.97: that is the stage-1 repulsion hand-over (PAIR_TABLE
  B), not the well — do not let the well fit it.
- **Join radius**: |E_pair| at the last grid point (3-3.5 r_eq), pinned forms: median 0.000, max
  0.57-0.69 kcal/mol (CO). The well ends by itself; it does not need to be multiplied by `w`.
  TOPO_REUSE part B's truncation disappears if `w` stays on the angle/torsion terms only and the
  list drop sits where the well's own tail is below tolerance (closed form from R, sigma or a, beta).
- **The Coulomb drift, quantified for this purpose**: 9 of 32 bonds have |dCoulomb(r_eq -> last)|
  >= 10 kcal/mol in this static protocol, with POSITIVE sign (+11 C-O, +16 C=C, +21 H-F, +22 O-H,
  +41 C=O), where EEQ_DRIFT_STATUS measured -11..-36 in the react protocol: the sign depends on
  the protocol, the reference Hirshfeld charges move <= 0.078 e either way. Fitted D raw -> q-frozen:
  C-O 84 -> 94, C=C 188 -> 202, O-H 105 -> 129, H-F 132 -> 152, C=O 227 -> 266, C-F 100 -> 113
  (8-24 %). **Fit s on the q-frozen rest (or after stage 2); a depth fitted on the raw rest absorbs
  a charge-model error whose sign flips with the topology protocol.** Median rms is lower q-frozen
  (1.32 vs 1.67), consistent with the drift being the model's, not the reference's.

### B.3 Does it sit with the erf bond order and the share? No — and that is the useful negative

The hope is that the well's own erf becomes the bond order the share apportions. The fitted
erfc does not sit where a switch sits: z0 = (r0 - R)/sigma has median **+0.38** over 32 bonds
(range -0.73..+5.9), i.e. the erfc's centre R lies INSIDE r0 and `1/2 erfc(z0)` = 0.30 at the
equilibrium distance. The well lives on the erfc's far flank, where erfc is a quasi-exponential
with a Gaussian roll-off; it is not "1 at r_eq, 0 at 2 r_eq". Forcing a switch-like placement
(z0 <= -1.5) gives h(z0) <= 0.06, and the curvature constraint then demands sigma ~ 0.05 A — a
step. So one erf cannot be both the share's bond order and the well; the families coincide, the
variables do not. For the H-transfer partial-well problem a wider, deeper well makes the pair sum
at the TS worse, not better (roadmap decision 5) — the share has to carry that alone, which is one
more reason for A.5 (1).

### B.4 Verdict for stage 3a(iii)

An erf-based well is a sound basis **only as the erf-Morse form, and there it is the MG form in
other coordinates** — identical quality on 32/32 curves, identical parameter count, identical
sourcing (class-A per pair, then the 3b element factorisation). It brings beta >= 0 for free and
reuses `erfc`; MG brings closed forms (k_e = 2 D a^2, x50/x90) where erf-Morse needs a bisection
for z0 at setup. Neither unifies with the share. Recommendation: build 3a(iii) in two steps,
(1) the curvature-pinned two-parameter variant (s, tail) — equilibrium-safe by construction,
rms 19 -> 2, to be confirmed on the guard harness (pooled 1.0341) because the cubic term moves
r_eq at the 0.003-A level (INFERRED, my r_min interpolation is not finer than 0.01 A);
(2) only then free the curvature with the r0 re-solve the roadmap review already specifies.
Take MG unless the operator values "one function family" for teaching — then erf-Morse costs
nothing. Remove `w` from the well in either case; fit depth on q-frozen rests.

## C. `origin/feature/multi-gpu` and `origin/reactff2-llm` (read-only; `git show`, `git merge-tree` to stdout)

**43e77c92 — on-the-fly CPU Coulomb, `coulomb_implicit` DEFAULT TRUE.** `generateGFNFFParameterSet()`
leaves `params.coulombs` empty and sets `params.coulomb_implicit`; `setInteractionLists()` copies
the flag, `partition()` gives every partition an ATOM range balanced by pair count, `calcCoulomb()`
loops i < j over `m_eeq_charges` and `m_coul_alp`. For rev-gfnff (READ, not run):
- The corner machinery is served: every corner is built by `beginTransition()` ->
  `rebuildInteractionLists()` = `setInteractionLists()` + `partition()`, the two functions the
  commit patches, and `TopologyState` already swaps `eeq_charges`, `coul_alp` and `partitions`, so
  per-corner EEQ / per-corner alpeeq / the 2^k blend read the right per-atom data. The jump-term
  measurement still gets its `coul` column (same accumulator). `rebuildReactiveTopology()` gets
  cheaper (no N^2/2 list per corner per rebuild).
- Two behavioural differences to the stored list: a NaN charge now SKIPS the pair (the stored path
  fell back to the static charge), and gamma_ij is recomputed from `m_coul_alp` instead of read —
  same expression, so INFERRED bit-identical at `-threads 1`; for T > 1 the partition boundaries
  move (atom ranges instead of pair ranges), so the reduction order and the last ulp change.
- **It breaks the project's bit-identity yardsticks by construction**: `-gfnff.dump_params` no
  longer contains a Coulomb list, so the recorded md5 (d297bc3b / 77c134bb) cannot survive; the
  12-digit energies probably do. Re-baseline both, or pin `-gfnff.coulomb_implicit false` in the
  guard scripts.

**298aed77 / ab6e3f5e — parallel topology loops, OpenMP budget.** Every new `omp parallel for`
writes row-local data and merges in row order (bond list, b3list, nbondmat levels, BFS rows, erf
matrix rows); there is **no `reduction` clause anywhere under `ff_methods/` on that branch**
(`git grep`), the pre-existing angle block merges under `omp critical` and then sorts by a total
order (j, i, k), and the CN loop accumulates per row in serial j order. So no reordered
floating-point reduction: topology and parameters are thread-count independent by design. What
DOES change: ab6e3f5e opens the OpenMP budget (`ScopedBlasThreads(m_threads)`), so loops that
"silently ran serial even with -threads 36" now really run parallel — at `-threads 1` (every
smoothness measurement) nothing changes, at T > 1 the react rebuild path executes them per
rebuild. The Known-Issue-#33 class (1 ulp amplified by react MD) is therefore not opened by these
two commits; the only ulp-level change is the Coulomb partitioning above, and only for T > 1.
One thing to check by eye at merge: 298aed77 gates an O(N^2) debug scan on `verbosity >= 3` —
harmless (print only), but it is the pattern #33 forbids for anything numerical.

**The three source conflicts.** (1) `gfnff_method.cpp`, 63 lines: this branch replaced the inline
geometric bond loop by `bonds = perceiveGeometricBonds();` (needed by the topo-reuse check
`76e7f83a` and by react), the remote parallelised the inline loop. Resolve: keep the call, move
the row-buffer parallelisation INTO `perceiveGeometricBonds()`; the list order must stay (i, j)
ascending because react and the reuse check compare graphs. (2) `gfnff_gpu_method_impl.h`: both
sides add an accessor (`getGFNFF()/gfnffInstance()` vs `gpuDevice()`), keep both. (3) `main.cpp`
is the semantic one: this branch's explicit `-batch true` JSONL mode (the class-A / guard / path
harnesses all stand on it, incl. `-batch_reuse_topology` = kept topology) against the remote's
AUTOMATIC batch for any multi-frame `-sp` input, spread over `-threads` workers. Keep ours first
(it returns), theirs after it; then a multi-frame `-sp` WITHOUT `-batch` changes meaning (was: one
structure; becomes: a threaded batch with a fresh calculator per worker) — grep the rev scripts
for that usage before merging.

**Order.** (a) `origin/reactff2-llm` first: 3 commits (ANCOpt atom holding, CIF reader,
`Molecule::addPair`), 9 files, none under `ff_methods/`, one trivial `AIChangelog.md` conflict; it
does not touch this work. (b) Decide `rev_budget_fix_h` and re-baseline the 20-cell grid BEFORE
the multi-gpu merge, so a trajectory change can be attributed. (c) Merge multi-gpu with
`coulomb_implicit` initially pinned false in the rev harness, confirm everything identical, then
flip it and measure the flip alone. **Re-measure after (c)**, all at `-threads 1` plus one T = 4
control: caffeine/benzene 12-digit energies for `gfnff` and `revgfnff`; FD gradient at rkt06 pt10
(1.19e-8) — it exercises Coulomb inside a blended corner; the hypervalent six; rkt06 rms 2.714;
the class-A kept protocol on ch4_C-H (rms 5.87, D_e dev -13.32: the `main.cpp` batch path);
`ctest -R "gfnff|react|cli_simplemd_1[3-9]"`; and the rebuild-count fingerprints of two cells
(c2h6/T2000_f16 135|204, ch4_H/T2000_f10 27|126 off|on) as the 1-ulp detector.

## Summary — where this review disagrees with the 2026-09-15 briefing

| briefing statement | verdict | basis |
|---|---|---|
| (3)/(4) H-budget mechanism; fix moves 2594 -> 59 and 1089 -> 223, T_max 62 098 -> 5 086 K | **confirmed to the digit** | own runs, A.1 |
| "Kein Falsifikator bewegt sich" | **confirmed** (rkt06 2.714 / max dE 0; six ions to 8 digits) and extended to 4 proton-transfer paths + H3+ | A.2 |
| residual 222 kJ/step and the 391 kJ first step are thermal / "in both arms" | **wrong**: carbon budget 4.53 -> 4.99 in two steps, five wells 0.77 -> 1.00; 438 of 487 grid events are these two cells | A.3 |
| "22-Zellen-Gitter" | **20 cells**; two `ch4_H` f16 cells never ran (rc = 1) | A.3 |
| carbon budget 4.95 is "vermutlich die Wurzel der Klasse-B-Frage" (H + CH4 TS 78 too low) | **no**: the 78 is a pre-share estimate; on rkt03 Val(C) = 4.00 throughout, TS error -13.5 not -78, capping C changes 0.00. The budget's real damage is the class-C adducts (-54..-107 kcal/mol, n = 4) | A.4 |
| open decision "which elements get a growing budget" | the element list is second order: the left-over share formula forfeits ALL valence at one partner too many; fix the formula (conserving, product) and grant excess by charge | A.5 |
| (1) proxy / QP closed, (2) break-tail = neighbour force constants x1.64 | **not re-measured here**; accepted from QP_STATUS / BREAK_TAIL (5-digit factor agreement), judged in A.6 | — |
| erf wells | possible only as erf-Morse = MG in other coordinates; no unification with the share | B |

Not done / limits of this review: no source build, so A.5 is offline arithmetic on the binary's
own share dump (exact on single-corner points, not evaluable where two corners blend — CH5+, and
the d = 1.4-1.75 A points of the class-C scans); the MD consequence of A.5 on the two f10 cells is
INFERRED; B is a frozen-rest fit on one static protocol; C is read, not merged or run.
