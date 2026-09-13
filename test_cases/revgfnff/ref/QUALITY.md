# Reference-data quality record (created 2026-09-13)

What in `ref/` is trustworthy, what is not, and which points a consumer must exclude.
Derived from `_log/UKS_INSPECTION.md` (per-series inspection, 2026-09-13) with the unstable
points re-extracted from the raw `energies.json` by the orchestrator.

**Operator decision (2026-09-13): label the weaknesses, do not recompute.** Nothing here is
silently usable — read this file before fitting anything to class A/H UKS data.

## 1. RKS is complete and unaffected

All **31 class-H RKS series are 20/20** and the **class-A RKS series are complete**. Everything
that rests on RKS — the `HIRSHFELD_CHARGES.md` drift table and the conclusion drawn from it — is
not affected by anything below. `ref/A/*_uks` is the older campaign and is complete except where
noted.

## 2. The ten class-H UKS series that failed or needed a retry (per-series verdict)

| series | class A | class H | verdict |
|---|---|---|---|
| `ch3cl_C-Cl` | 20/20 complete | 0/20 | **A usable as-is**; H unusable |
| `co_CTO` | 0/20 | 0/20 | **no UKS state exists** — the reference is pure RKS |
| `of2_O-F` | 16/20 (SlowConv retry) | 0/20 | **A usable for those 16 points only**; no far point |
| `hocl_O-Cl` | 20/20 complete | 1/20 | **A usable as-is** |
| `h2co_CDO` | 1/20 (that point IS the far end) | 1/20 | **A usable for D_e only** |
| `h2o2_O-O` | 1/20 (far point) | 1/20 | **A usable for D_e only** |
| `clf_F-Cl` | 16/20 (SlowConv retry) | 1/20 | **A usable for those 16 points only**; no far point |
| `ncl3_N-Cl` | 20/20 complete | 2/20 | **A usable as-is** |
| `cl2_Cl-Cl` | 4/20 (the four far points) | 4/20 | **A usable, far region only** |
| `o2_ODO` | 20/20 complete | 16/20 | **A usable as-is** (there is no RKS series) |

Where both classes converged they **agree in state**: every class-H point with a class-A
counterpart matches within **0.052 mEh** and `dS^2 <= 7.2e-5`. The `--slowconv` retries of
`of2`/`clf` are therefore **the same electronic state**, not a different one.

`co`, `h2co`, `h2o2`, `cl2` failed identically in both classes on the same points, energies
matching — those failures are a property of the point. `ch3cl`, `hocl`, `ncl3`, `o2` failed in H
while the class-A run of the same grid converged 20/20, and for `ncl3` the two `job.inp` files are
**byte-identical except the `Print[P_Hirshfeld]` line** — those are run-dependent SCF instability,
not a property of the point.

## 3. EXCLUDE THESE POINTS — class A and class H land on different broken-symmetry solutions

Five UKS series (all 20/20 in both classes, nobody had compared them before) have points where two
runs of the same input give **different BS solutions**. A consumer taking `min(RKS, UKS)` or
reading a class-H point as a stand-in for class A must drop these radii:

| series | label r | dE (class H - A) | dS^2 |
|---|---:|---:|---:|
| `c2h2_CTC_uks` | 2.4008 | +25.569 kcal/mol | -1.0868 |
| `f2_F-F_uks` | 1.5400 | +47.245 | +1.0038 |
| `f2_F-F_uks` | 1.8200 | -15.481 | -0.3781 |
| `f2_F-F_uks` | 1.9600 | +6.869 | +0.2013 |
| `f2_F-F_uks` | 2.8000 | +0.952 | +0.0002 (marginal; drops below the 0.5 kcal threshold if that is raised) |
| `hcn_CTN_uks` | 3.4520 | +22.802 | -0.6571 |
| `n2_NTN_uks` | 3.2822 | -55.165 | +0.0212 |
| `n2_NTN_uks` | 1.6411 | +15.273 | -0.5475 |
| `n2h2_NDN_uks` | 1.7335 | +36.172 | +0.3384 |
| `n2h2_NDN_uks` | 1.8573 | -37.605 | -0.1679 |
| `n2h2_NDN_uks` | 1.9811 | +40.547 | -0.0761 |
| `n2h2_NDN_uks` | 2.4764 | +20.638 | -0.1361 |

**Consequence: class H cannot stand in for class A on these five series** — for these radii they
are not the same calculation. For the other class-H series the two classes agree.

## 4. Two data weaknesses to know about

- **`of2` and `clf`: the class-A trees were overwritten by the retry.** The pre-retry numbers exist
  nowhere else and are not in git, and the retry changed two things at once (`--slowconv` *and*
  `--uks-inside-out`), so the recovery cannot be attributed to either. If an O-F or F-Cl reference
  `D_e` is ever needed, it must be recomputed — the `min(RKS, UKS)` convention would otherwise
  substitute the higher RKS value at the missing far point, making the reference well too deep by
  about **+53.1 (`of2`) / +51.5 (`clf`) kcal/mol** and spuriously strengthening a "the model's well
  is too shallow" conclusion.
- **`revgfnff_hirshfeld.py:pick()` accepts non-converged points** (it tests `hirshfeld != null`,
  not `energy_eh != null`, and falls back to the nearest available point). In `HIRSHFELD_CHARGES.md`
  the UKS rows for `h2o2_O-O`, `hocl_O-Cl`, `clf_F-Cl`, `ch3cl_C-Cl` therefore report "near-r_eq"
  and "far" as **the same single point** (change exactly 0.000), and `h2co_CDO`, `ncl3`, `cl2`
  take "near-r_eq" from a point 2.4-3.6 A away, `o2` takes "far" from r=2.72 instead of 4.23. The
  **RKS rows are correct** (31/31 complete) and the drift conclusion rests on them.

## 5. What to do when consuming this data

1. Use the **RKS** rows for charge and drift questions.
2. For a `min(RKS, UKS)` well-shape number, exclude the radii in section 3, and drop `of2`/`clf`
   far points rather than substituting RKS.
3. Do not treat a class-H point as a class-A stand-in on the five series of section 3.
