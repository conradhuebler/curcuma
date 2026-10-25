# Delta-ML active learning with prospective validation: measurement (Oct 2026)

Status: 🤖 AI-generated scripts (`scripts/deltaml_usecase9/`), ⚙️ machine-tested, peo8 only (58 atoms, one start structure), not human-reviewed. Measured
by a Sonnet agent with the binary of commit f585879f. Follows [DELTA_ML_TRAINING_DATA_2026-10.md](DELTA_ML_TRAINING_DATA_2026-10.md).

What changed in f585879f: reproducible sampling MD for a given `-activelearn.seed` and thread count (per-instance thermostat engines), 10 validation frames per
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
