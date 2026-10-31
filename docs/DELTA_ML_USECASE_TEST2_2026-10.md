# Delta-ML application test, second round: consistent labels, wide energy range, minima (Oct 2026)

Status: 🤖 AI-generated scripts, ⚙️ machine-tested, two molecules, not human-reviewed. Measurements by Sonnet agents
(`scripts/deltaml_usecase2/`, `scripts/deltaml_usecase3/`; absolute paths of the original working directories are
hard-coded). **This round supersedes the Tests 1 and 2 of [DELTA_ML_USECASE_TEST_2026-10.md](DELTA_ML_USECASE_TEST_2026-10.md).**

## Why a second round: two label defects of `-activelearn`

The first round used datasets written by the internal labelling of `-activelearn`. Two defects were found afterwards:

1. **Gradient labels 1.89 times too large** (`gradient_base`, `gradient_ref`): `LabelThread` multiplied
   `EnergyCalculator::Gradient()` by `ANGSTROM_TO_BOHR`, but the contract is Eh/Angstrom since the merge with master
   (Known Issue #28). Dataset label against `-sp -dump_gradient`: ratio 1.8906 (n=1 structure) before, 1.0000 (max deviation
   2e-16) after the fix in commit bc7742e6. The energy rows and the gradient rows of the fit contradicted each other.
2. **Stored base energy differs from a fresh single point** by -1.21 +- 0.70 kcal/mol (353 conformers, min -3.05),
   -1.13 +- 0.44 (M1 thermal set) and -5.07 +- 1.98 (M2, min -13.6): training and deployment used different base
   energies. Cause under investigation; also, the datasets come out in thread-chunk order, not in input order.

Everything measured with the internal labelling (the joint-model parameter scan and its energy/force trade-off, the first
application test, the held-out numbers in DELTA_ML.md) is affected; the energy-only model is not affected by defect 1.
The round below labels every structure externally: `-sp -method gfnff -gfnff.cache_topology false -dump_gradient` and
`-sp -method gfn2 -dump_gradient`, one fresh directory per structure, dataset JSON built by script, training with
`-activelearn.dataset`. Check in every run: total energy of `-sp -method deltaml:<model>` minus `deltaml_correction`
equals the fresh GFN-FF energy to 5e-9 Eh (all 706 + 802 + 353 + 102 predictions).

Models: A = default (1000 centres), C = energy-only model, S = joint with lambda 1e-3. Model B of the first round was not rerun.

## Part I: the narrow set again (353 conformers of M1, GFN2 span 2.2 kcal/mol; relabelled)

5-fold CV x 2 repeats (n = 10 folds, not independent), fresh GFN-FF as baseline:

| | Spearman | pair accuracy (|dE| > 1) | recall@5 | regret [kcal/mol] | centred MAE [kcal/mol] |
|---|---|---|---|---|---|
| GFN-FF | 0.52 +- 0.09 | 0.86 | 0.30 | 0.63 | 0.62 |
| A | 0.68 +- 0.07 | 0.95 | 0.48 | 0.29 | 0.38 |
| C | 0.81 +- 0.04 | 0.99 | 0.72 | 0.06 | 0.24 |
| S | 0.44 +- 0.06 | 0.81 | 0.32 | 0.60 | 0.52 |

Paired Spearman against GFN-FF: A +0.161 +- 0.113 (9 of 10 folds better), C +0.289 +- 0.067 (10/10), S -0.079 +- 0.088 (2/10);
recall@5 A +0.18 (6/10), C +0.42 (9/10). The first round (inconsistent labels) gave A -0.179 +- 0.136 (2/10) and
C +0.155 +- 0.078 (9/10) in the deployment frame. The range is unchanged, so the failure of the first round came from the
labels, not from the narrow range; which of the two defects (gradient scale or base offset) matters was not separated.
Blocked folds on this set (5 blocks in file order, 5-frame gap, n=5) give the same pattern (A +0.134, C +0.285, S -0.082);
whether the file order is a time series was not checked.

## Wide range: 400 K GFN-FF MD snapshots (401 frames per molecule)

M1 (90 atoms): GFN2 span 79.4 kcal/mol, std 12.5; M2 (trisaccharide, 66 atoms): span 62.2, std 11.2. GFN-FF against GFN2:
Spearman 0.890 / 0.775, RMSE after linear fit 5.49 / 6.95 kcal/mol. Nearest-neighbour RMSD between frames 0.66 / 0.39 A.

| paired against GFN-FF, n = 10 folds, Spearman | A | C | S |
|---|---|---|---|
| M1 | +0.064 +- 0.020 (10/10) | +0.041 +- 0.015 (10/10) | +0.031 +- 0.013 (10/10) |
| M2 | +0.180 +- 0.051 (10/10) | +0.152 +- 0.047 (10/10) | +0.134 +- 0.041 (10/10) |

Pair accuracy (|dE| > 1): M1 0.87 -> 0.93 (A); M2 0.80 -> 0.92. Spearman in the lowest 25 % by GFN2: M1 0.65 -> 0.81 (A),
M2 0.44 -> 0.76. Centred MAE: M1 4.50 -> 2.80, M2 5.57 -> 2.81 kcal/mol. The top end is not reliably better: recall@5
M1 0.56 -> 0.68 (+0.12 +- 0.27, 5/10 folds), regret M1 -1.39 +- 2.54 (5/10).

**Leakage check (Part II):** blocked folds (5 contiguous time blocks of 80 frames, 5 frames dropped at each block edge, n = 5):
A +0.066 +- 0.025 (M1, 5/5; random folds +0.064) and +0.132 +- 0.024 (M2, 5/5; random folds +0.180). Interpolation between
neighbouring frames explains about 0 % (M1) to 27 % (M2) of the gain; revisited basins in other blocks still leak, so this
bounds the leakage from above only for the time-neighbour part. At the top end the blocked result is weaker: regret of A on M1
+1.24 +- 2.40 kcal/mol (0/5 better).

## Pre-optimisation from the MD frames (20 frames per molecule, A and S trained on the other 381)

The optimiser `lbfgs` fails for large start gradients (about > 0.3 Eh/Bohr, which these frames have): zero first step, exit 1, no
error text, also for GFN2 and GFN-FF alone (converged runs M1: GFN2 4/20, GFN-FF 1/20, A 4/20, S 0/20). Results below
use the default optimiser (all 160 runs exit 0 except one stalled run, M2/A frame 164).

| | RMSD to GFN2 optimum [A] | GFN2 energy at the result minus GFN2 optimum [kcal/mol] |
|---|---|---|
| M1 GFN-FF / A / S | 1.20 / 1.36 / 1.27 | 11.20 / 7.55 / 7.23 |
| M2 GFN-FF / A / S | 0.84 / 0.70 / 0.48 | 23.80 / 12.65 / 12.38 |

Paired against GFN-FF, GFN2 energy at the result: M1 A -3.64 +- 3.16 (18 of 20 better), S -3.96 +- 1.52 (20/20);
M2 A -11.24 +- 7.76 (17 of 19), S -11.42 +- 2.09 (20/20). RMSD: M1 not better (A +0.16 +- 0.53, S +0.07 +- 0.31), M2 S -0.36 +- 0.35
(17/20). The model's own energy at the result is below GFN2 for A (M1 -2.4, M2 -10.8 kcal/mol, mild spurious lowering) and above it
for S (+15.6, +7.4). The energies are relative to GFN2's optimum from the same start, not to the global minimum.

## Part III: ranking of conformer minima (the actual use case; one split per molecule)

60 evenly spaced MD frames per molecule optimised with GFN-FF, duplicates removed (RMSD < 0.3 A and |dE_GFN-FF| < 0.2 kcal/mol),
single points of GFN-FF and GFN2 at the minima, models trained on the 341 other frames. M1: 60 minima (GFN2 span 8.2 kcal/mol);
M2: 42 minima (58 of 60 optimisations converged, 16 duplicates; span 15.9).

| Spearman (pair accuracy) | GFN-FF | A | C | S |
|---|---|---|---|---|
| M1, n = 60 | 0.152 (0.52) | 0.531 (0.80) | 0.473 (0.78) | 0.331 (0.71) |
| M2, n = 42 | 0.445 (0.67) | 0.791 (0.85) | 0.690 (0.80) | 0.685 (0.80) |

Paired Spearman difference with 95 % bootstrap interval over the minima: M1 A +0.379 [-0.003, +0.712], C +0.321 [-0.105, +0.698],
S +0.180 [-0.235, +0.569]; M2 A +0.346 [+0.016, +0.730], C +0.245 [-0.142, +0.657], S +0.240 [-0.072, +0.588]. Pair accuracy
is higher for A, C and S on both molecules (the pairs are not independent). The lowest minimum is hardly ever found: top-1 hit
only for A on M1; recall@3 between 0 and 0.67; regret of A 0.00 against 2.85 (M1) and 0.77 against 3.85 kcal/mol (M2).

## Conclusion (two molecules, trajectory-based sets)

With consistent labels the Delta-ML orders conformers better than GFN-FF on every set tried: the narrow M1 set (C +0.289,
A +0.161 Spearman), the thermal sets (A +0.064 / +0.180, kept to 73 % to 100 % with blocked folds) and, as a tendency, the
GFN-FF minima (pair accuracy, Spearman; the intervals mostly include 0). The choice of the lowest minimum is not shown to
improve. Pre-optimisation lowers the GFN2 energy at the end point (S in 20 of 20 cases on both molecules) but not the geometry on M1.
Which model is best depends on the set: C on the narrow set, A on the wide and minima sets; S (lambda 1e-3) is worse than A at
ranking. Cost (first round, not affected by labels): 11 to 25 times GFN-FF per step.

Not tested: other molecules, a DFT reference, charged or open-shell systems, a minima set built by a real conformer search
(ensembles with several independent MD runs and different basins), longer-range generalisation (blocks of whole basins),
model B, the force accuracy and the energy/force trade-off with the corrected gradient labels (the "4x better forces" of the
earlier documents need a new measurement), hyper-parameters after the label fix (all scans were run with the defective labels).
Caveats: folds share frames and seeds; one MD run per molecule; Part III uses one split; no human review.
