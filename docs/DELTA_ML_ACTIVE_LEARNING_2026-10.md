# Delta-ML active learning with prospective validation: measurement (Oct 2026)

Status: 🤖 AI-generated scripts (`scripts/deltaml_usecase9/`), ⚙️ machine-tested, peo8 only (58 atoms, one start structure), not human-reviewed. Measured
by a Sonnet agent with the binary of commit f585879f. Follows [DELTA_ML_TRAINING_DATA_2026-10.md](DELTA_ML_TRAINING_DATA_2026-10.md).

What changed in f585879f: reproducible sampling MD for a given `-activelearn.seed` and thread count (per-instance thermostat engines; replaced at the merge of devel/methods on Oct 31, 2026 by the instance engine `m_rng`, see DELTA_ML.md), 10 validation frames per
cycle scored with the model that drove the cycle (before the refit), and `al_criterion` uncertainty | validation | both with `al_target_energy` 1.0 kcal/mol,
`al_target_hole` 2.0 kcal/mol, `al_patience` 2. Setup of the runs as in the training-data experiment (300 K, 4 MD runs of 10 ps per cycle, 15 candidates, 15 initial, indicator
threshold 0.3, 8 threads, at most 20 cycles); every final model then ran 4 MD seeds of 40 ps (NVT 300 K) against the pooled GFN2 seeds 2 to 4 (noise floor of the torsion
JSD 0.198, standard error of one JSD 0.04 to 0.07).

## Results

| run | cycles | converged | labels | wall [s] | torsion JSD (4 MD seeds) | model minus GFN2 on own frames: mean / deepest [kcal/mol] |
|---|---|---|---|---|---|---|
| N0 defaults (both, 1.0 / 2.0, patience 2), seeds 1 / 2 / 3 | 20 / 20 / 20 | no / no / no | 379 / 376 / 401 | 495 / 496 / 509 | 0.198 / 0.211 / 0.253 (mean 0.221) | -0.7 / -9.1, -0.8 / -11.1, -1.1 / -11.2 |
| N1 relaxed targets (2.0 / 4.0) | identical to N0 (curve md5 equal): never met twice in a row | no | - | - | - | - |
| N2 = N0 with 4000 centres, seeds 1 / 2 | 20 / 20 | no / no | 486 / 500 | 1507 / 1518 | 0.213 / 0.217 | -0.7 / -9.0, -0.7 / -8.5 |
| N3 `uncertainty` only (old stop rule 0.3), seeds 1 / 2 / 3 | 17 / 16 / 14 | yes | 330 / 307 / 285 | 432 / 405 / 356 | 0.212 / 0.221 / 0.272 (mean 0.235) | -0.7 / -9.1, -1.4 / -9.7, -1.4 / -10.5 |
| old active learning (before f585879f), 3 runs | 12 | no | 163 to 180 | about 260 | 0.236 / 0.256 / 0.239 (mean 0.244) | -3.0 / -12.8, -2.5 / -15.1, -1.3 / -11.7 |
| reference: V0 B2 (400 K GFN-FF data), GFN-FF | - | - | 351 / 0 | - | 0.321 / 0.271 | -6.6 / -15.8 / - |

The N3 curves are byte-for-byte prefixes of the N0 curves (md5 checked): the same runs stopped early. Half of the labels of N0 are validation frames (10 per cycle).

- **Reproducibility:** N0 seed 1 repeated with the same seed and 8 threads (run alongside other jobs): learning curve, dataset and model JSON have identical md5 (3 of 3 files). One seed and one thread count.
- **Quality:** N0 torsion JSD 0.221 against the old 0.244 (-0.023, not resolved), against V0 B2 -0.10 (closer to the reference in 11 of 12 seed pairs), against GFN-FF in 9 of 12; one of three N0 runs
  reaches the noise floor (0.198, equal not below). The holes are flatter on average (mean -0.9 against -2.3 kcal/mol), single frames still lie 9 to 11 kcal/mol below GFN2.
  4000 centres do not change the torsion JSD (0.213, 0.217) but halve the energy distance (W1 of E_GFN2 2.9 and 3.4 against 5.9 kcal/mol) at three times the wall time; the uncertainty indicator never reaches 0.3 with them.
- **The new convergence criterion never fired (0 of 8 runs, 20 cycles):** with the default targets the validation was never met (smallest MAE 1.0 to 1.4 against 1.0, smallest hole 2.8 to 3.3 against 2.0); with
  the relaxed targets it was met in 3 to 4 cycles per run, never twice in a row. A plausible cause that was not tested: ten validation frames let the MAE of one and the same model fluctuate
  by 1.4 to 3.2 kcal/mol between cycles (peo8, N0 seed 1: 1.35 to 2.29 in the last ten cycles), the hole is the extremum of ten frames (about -4 to -5.5 kcal/mol at a typical MAE of 1.7) and contains a constant offset
  of the model. The targets 1.0 (MAE) and 2.0 (hole) are not consistent with each other at this scale. The question whether validation-converged runs are better than uncertainty-stopped runs could not be answered.
- **Cost:** N0 uses 2.1 to 2.4 times the labels of the old procedure and about 1.9 times the wall time (the time is the 10 ps sampling MD of each cycle, not the labels). Plain `uncertainty` gives models
  of the same quality with 285 to 330 labels.

## Conclusion

Reproducible sampling and the prospective validation numbers (energy error and hole depth per cycle, from -19 kcal/mol in cycle 1 to about -4 at the plateau) are delivered and useful as diagnostics. The
validation-based stop is not usable with the defaults: it did not trigger once, and the quality of the final models does not distinguish runs that would have stopped by uncertainty from those that did not.
The torsion populations improve over the first model and over GFN-FF but the improvement over the old active learning is not resolved. Not tested: other molecules or start structures, `al_onthefly` with validation, more
than 20 cycles, pooled or smoothed validation (several cycles, centred errors, a quantile instead of the extremum).

## Addendum (27 Oct 2026): pooled, centred validation criterion

Status: 🤖 AI-generated (`src/capabilities/activelearn.{h,cpp}`, scripts `scripts/deltaml_usecase10/`), ⚙️ machine-tested, peo8 only, not human-reviewed.
Setup as above (binary of HEAD 75d03431 plus this change, 8 threads, 3 processes at once).

What changed: the criterion pools the validation frames of the last `al_window` (3) cycles (30 frames) and re-predicts them with the model entering the
cycle; pool frames stay out of training until they leave the window (otherwise 20 of 30 frames would be training data of the model being judged), cycle-0 frames
go straight into training, the criterion is only evaluated on a full pool. Errors are centred (offset reported separately), the hole is the 10 % quantile
(`al_hole_quantile` 0.9) of the centred signed errors. `both` now means validation in `al_patience` consecutive cycles and the indicator in the stopping cycle:
in the tuning runs the indicator was <= 0.3 only in single cycles, never twice in a row, so the old `both` could not fire with any target (0 of 3).
`uncertainty` keeps its learning bit-identical (model and dataset md5 equal to HEAD, dimethyl ether, 6 cycles, one seed; the curve file gains columns 12-15
and a line for the converged cycle).

Targets from tuning runs with seeds 4, 5, 6 (unreachable target, 20 cycles; replaying the rule on the recorded curve is exact because the trajectory does not
depend on the targets until the stop): pooled centred MAE over cycles 14-20 median 2.12, 75th percentile 2.47; hole quantile median -3.20, 25th percentile
-3.57 kcal/mol (n = 21 correlated cycle values, 3 runs). Chosen `al_target_energy` 2.5, `al_target_hole` 3.5: met in 15 of 21 plateau cycles and in 1 of 15
cycles 4-8. Pooled values still fluctuate (one run: MAE 1.15 to 2.45 within cycles 14-20).

| run (new defaults, seeds as N0) | stop | labels | wall [s] | torsion JSD (4 MD seeds) | model minus GFN2 own frames: mean / deepest |
|---|---|---|---|---|---|
| F seed 1 / 2 / 3 | cycle 18 / none in 20 / cycle 19 | 373 / 386 / 362 | 612 / 684 / 654 | 0.223 / 0.174 / 0.266 (mean 0.221) | -0.9 / -7.8, -0.9 / -9.5, -1.6 / -9.9 |
| N0 (fixed 20 cycles, above) | none | 376 to 401 | 495 to 509 | 0.198 / 0.211 / 0.253 (mean 0.221) | |
| N3 (`uncertainty`, above) | 17 / 16 / 14 | 285 to 330 | 356 to 432 | 0.212 / 0.221 / 0.272 (mean 0.235) | |

Pooled curve at the stop (centred MAE / offset / hole quantile, kcal/mol): seed 1 cycle 17-18 2.09/-1.03/-2.90, 1.97/-0.46/-3.23; seed 3 cycle 16-19
2.09/-1.57/-2.56 to 1.69/-0.63/-3.04; seed 2 met the validation from cycle 18 on but the indicator stayed 0.32 to 0.50. Full per-cycle curves:
`scripts/deltaml_usecase10/curves10.py`. The plain `validation` criterion would have stopped at cycles 7 / 19 / 9 on the same runs (replayed).

- The criterion fires in 2 of 3 runs, late (18, 19), and the timing is set by the indicator, not by the validation (validation held from cycle 17 and 16).
- Final quality equals N0 (mean JSD 0.221 both; between-run spread 0.046, SE of one JSD 0.04 to 0.07): no gain and no loss resolved, n = 3.
- Wall time is 20 to 38 % above N0; not root-caused (N0 used the older binary f585879f; the tuning runs took 746 to 751 s with the same setup).
- Not tested: other molecules, start structures, thread counts, `al_onthefly` with validation, `validation` alone as default, whether a model stopped at cycle 7
  to 9 by `validation` is as good, other window lengths and quantiles.
