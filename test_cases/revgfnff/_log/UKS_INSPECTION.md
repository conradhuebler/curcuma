# UKS reference inspection -- class A vs class H (2026-09-13)

Test: each class-H point against the class-A point at the **label r** (not the index -- 4 series
differ in point order). Match = `|dE| <= 1e-6 Eh` and `|dS2| <= 1e-4`. Sources:
`ref/A/<sys>/energies.json`, `ref/H/<sys>/energies.json`.

## What is on disk

- **Class A is the reference.** 27 of 32 UKS series are the original 09-11 outside-in `TightSCF`
  runs; four are re-runs: `hcn_CTN_uks` (SlowConv+inside-out, 20/20), `o2_ODO_uks` (inside-out, no
  SlowConv, 20/20), `of2_O-F_uks`/`clf_F-Cl_uks` (SlowConv+inside-out, 16/20). **The pre-retry of2
  (1/20) and clf (0/20) class-A data is gone** -- the retry wrote into `ref/A/<sys>/`; no other copy
  exists (untracked in git, absent from the curcuma-topo worktree).
- **Class H** = 09-13 outside-in `TightSCF` + Hirshfeld, all 32 UKS series, never re-run. `cl2`,
  `co`, `h2co`, `h2o2` class-A trees carry a stale `job.inp` (SlowConv, 09-11 22:17) from a discarded
  recovery experiment; `energies.json` is the original run, not overwritten.
- Class H is **not** byte-identical to A: its input carries A's coordinates **rounded to 5 decimals**
  (max deviation 2.8e-6 A, ~1e-8 Eh). Only 6/32 series are byte-identical -- including the two the
  campaign report verified (`h2_H-H_rks`, `hf_H-F_uks`), which is why that check passed.

## The 10 affected series

| series | class A on disk | A missing at r | H conv. | H==A | H differs | H null | verdict |
|---|---|---|---|---|---|---|---|
| ch3cl_C-Cl | 20/20 complete | -- | 0/20 | -- | -- | 20 | **A usable as-is**; H unusable |
| co_CTO | 0/20 | all | 0/20 | -- | -- | 20 | **not usable** -- no UKS state at all |
| of2_O-F | 16/20 SlowConv retry | 3.52-4.93 | 0/20 | -- | -- | 20 | **A usable, 16 pts, no far point** |
| hocl_O-Cl | 20/20 complete | -- | 1/20 | 1 | 0 | 19 | **A usable as-is**; H unusable |
| h2co_CDO | 1/20 (far pt only) | 0.90-3.61 | 1/20 | 1 | 0 | 19 | **A usable for D_e only** (1 pt) |
| h2o2_O-O | 1/20 (far pt only) | 1.10-4.41 | 1/20 | 1 | 0 | 19 | **A usable for D_e only** (1 pt) |
| clf_F-Cl | 16/20 SlowConv retry | 4.14-5.79 | 1/20 | 0* | 0 | 19 | **A usable, 16 pts, no far point** |
| ncl3_N-Cl | 20/20 complete | -- | 2/20 | 2 | 0 | 18 | **A usable as-is** |
| cl2_Cl-Cl | 4/20 (the 4 far pts) | 1.52-4.57 | 4/20 | 2 | 0 | 16 | **A usable, 4 far pts only** |
| o2_ODO | 20/20 complete | -- | 16/20 | 12 | 0 | 4 | **A usable as-is** (no RKS series) |

\* clf's one H point (r=5.7947, E=-559.8403487, S^2=1.006252) is at an r the retry never reached --
no overlap, not a disagreement.

**No point in the 10 series differs materially.** Every H point with an A counterpart agrees within
0.052 mEh (0.033 kcal/mol) and `dS2 <= 7.2e-5`; the only sub-1e-6 misses are 2 cl2 and 4 o2 points
(max 0.052 / 0.0104 mEh). Where both classes converged, they agree in state.

**Which failures are robust.** `co`, `h2co`, `h2o2`, `cl2` failed identically in both classes, on the
same points, energies matching -> real. `ch3cl`, `hocl`, `ncl3`, `o2` failed in H while the A run of
the same grid converged 20/20; for `ncl3` the two `job.inp` files are **identical except the
`Print[P_Hirshfeld]` line** (byte-identical `points.xyz` too), yet A gave 20/20 and H 20 -> 2. Those
four H failures are run-dependent SCF instability, not a property of the point.

## Two findings outside the 10

1. **Five further UKS series disagree materially between A and H** (both 20/20; nobody had compared
   them): `c2h2_CTC` 1 pt (r=2.401, +26 kcal, dS2 -1.09), `f2_F-F` 3 pts (to +47 kcal, +1.00),
   `hcn_CTN` 1 pt (+23, -0.66), `n2_NTN` 2 pts (to -55, +0.02), `n2h2_NDN` 4 pts (to +41, -0.14).
   Different broken-symmetry solutions -> **class H cannot stand in for class A anywhere**.
2. **The class-H charge consumer picks the wrong geometry where UKS is incomplete.**
   `revgfnff_hirshfeld.py:pick()` requires only `hirshfeld != null`, not `energy_eh != null`, and
   falls back to the nearest available point. In `HIRSHFELD_CHARGES.md`: `h2o2_O-O_uks`,
   `hocl_O-Cl_uks`, `clf_F-Cl_uks`, `ch3cl_C-Cl_uks` report "near-r_eq" and "far" as **the same
   single point** (change exactly 0.000); `h2co_CDO_uks`, `ncl3`, `cl2` take "near-r_eq" from a point
   2.4-3.6 A away; `o2` takes "far" from r=2.72 instead of 4.23. The **RKS** rows (31/31 complete)
   and the `ORCA_REF_STATUS` drift conclusion, which rests on RKS, are unaffected.

## 1. Verdicts

Usable as-is: `ch3cl`, `hocl`, `ncl3`, `o2` (A complete; H adds nothing). Usable for the points that
have A: `cl2` (4 far points), `of2`, `clf` (the 16 recovered). Usable for D_e only, not for a curve:
`h2co`, `h2o2` (one far point each). Not usable: `co` (0 UKS points; the reference is pure RKS).
Not usable at all: the class-H series of `ch3cl`/`of2`/`co` (0 points) and the class-H near-r_eq
column of `h2o2`/`hocl`/`clf`/`ch3cl`/`h2co`/`ncl3`/`cl2`/`o2` (finding 2).

## 2. Does the missing UKS move the reference D_e?

`revgfnff_curves.py` keys `min(RKS,UKS)` on r and fills a missing UKS with the **RKS** value, which
at the far end is 45-60 kcal/mol above the BS state. D_e = `E_ref(r_max) - E_ref(min)`.

| system | D_e as-is | D_e with UKS-missing r dropped | bias | far-point source | RKS-only counterfactual |
|---|---:|---:|---:|---|---:|
| h2co_CDO | **215.1** | 215.1 | **0** | UKS (converged) | 263.7 |
| h2o2_O-O | **48.2** | 48.2 | **0** | UKS (converged) | 108.6 |
| hocl_O-Cl | **50.8** | 50.8 | **0** | UKS (complete) | -- |
| o2_ODO | **142.6** | 142.6 | **0** | UKS only, complete | (no RKS series) |
| cl2_Cl-Cl | 73.0 | 73.0 | 0 | UKS (the 4 far pts) | 92.4 |
| ch3cl / ncl3 | 86.9 / 32.8 | same | 0 | UKS (complete) | -- |
| **of2_O-F** | **92.9** | **39.8** | **+53.1** | **RKS** | -- |
| **clf_F-Cl** | **107.4** | **56.0** | **+51.5** | **RKS** | -- |
| **co_CTO** | **351.2** | **n/a (no UKS)** | **n/a** | **RKS only** | 351.2 |

- The four systems you named are **unbiased**: `h2co`/`h2o2` because their one converged UKS point
  *is* the far end, `hocl`/`o2` because their A series are complete; `cl2` likewise.
- **Direction where the bias exists: upward.** A missing far UKS point makes `E_ref(r_max)` the
  higher RKS value, so the reference well is too **deep**, i.e. the "model's well is too shallow"
  conclusion is **strengthened spuriously**. Size: 53.1 (of2), 51.5 (clf) kcal/mol; `co` has no UKS
  point at all, so its 351.2 kcal/mol is 100 % the ionic-limit RKS curve.
- The missing *inner* points (h2co 0.9-3.6 A, h2o2 1.1-4.4, cl2 1.5-4.6) do not touch D_e but do
  distort the tail: substituting RKS between the RKS/UKS crossing (~2 r_eq) and the far point lifts
  `E_ref(r)` by up to ~50 kcal/mol, so `rms_tail`/`max_tail` for those bonds are measured against a
  reference that is too steep there.

## 3. Is the SlowConv retry the same state?

**Only indirectly answerable for of2/clf, and the evidence says yes.**
- No common r to compare: of2's A tree was overwritten by the retry and H has 0 points; clf's one H
  point sits at a label the retry never reached.
- clf tail test: fitting `E = Einf + C/r^n` to the retry's last 4-8 converged points predicts
  E(5.7947 A) = -559.8385..-559.8442 Eh. The class-H point (plain `TightSCF`, outside-in, **no**
  SlowConv) is **-559.8403487** -- inside that band, 0.03 kcal from the best fit -- while RKS there is
  50.6 kcal higher and its `S^2 = 1.006252` continues the retry's own trend (1.006245 at r=3.73).
  Same state family, not a damped closed-shell artefact.
- Precedent: `hcn_CTN_uks` is the one series with both a SlowConv+inside-out chain (A, 20/20) and a
  plain outside-in chain (H, 20/20). They agree to 1e-6 Eh and 1e-4 in S^2 at 19/20 points; the
  single 22.8 kcal miss at r=3.452 is a different BS solution *in both*, i.e. SCF path sensitivity
  that is equally present without SlowConv.
- Real caveat: the retry changed **two things at once** (`--slowconv` *and* `--uks-inside-out`), and
  the direction change is the likelier cause of recovery -- walking outward hands each hard far point
  a physical BS guess, which is why the 4 remaining failures are the farthest points. The *state*
  question is settled; which keyword did it is not.

## 4. What would be needed, and cost

- **co** (0/20): needs a BS guess for C#O itself. First try the retry that fixed of2/clf
  (`--slowconv --uks-inside-out`, ~4 min at 4 cores); else a fractional-occupation / `GuessMode`
  start from a pre-converged stretched solution. Until then do not use the 351.2 kcal/mol RKS D_e.
- **of2 / clf far points** (4 each): retry with a finer far-end step, or restart from the last
  converged point's `.gbw`. ~5 min per system. Meanwhile quote D_e as `<= 39.8` / `<= 56.0` kcal/mol
  or drop the pair -- do **not** quote 92.9 / 107.4.
- **h2co / h2o2 / cl2 inner points**: only the tail shape is missing; accept the RKS-filled tail with
  a note, or run the inside-out chain (which recovers exactly these). ~5 min each; D_e sound anyway.
- **Robustness, not coverage, is the real gap**: 5 more series have points where two runs of the same
  input land on different states. A reference needs each point pinned to a definite solution (re-run
  from the other solution's `.gbw`, or make `S^2` a checked run output) -- otherwise "converged"
  means "converged to whatever it found". Worth doing before this reference is consumed.
- **Not determined**: why two byte-identical inputs diverge (SCF nondeterminism vs environment) -- no
  re-run was possible in this read-only pass; and the pre-retry of2/clf class-A numbers, unrecoverable.
