# Delta-ML sampling: training data or model class? (Oct 2026)

Status: 🤖 AI-generated scripts (`scripts/deltaml_usecase8/`), ⚙️ machine-tested, peo8 only (58 atoms, one start structure), not human-reviewed. Measured
by a Sonnet agent with the binary that contains the label repairs and the `-md.seed` fix. Follows [DELTA_ML_SAMPLING_2026-10.md](DELTA_ML_SAMPLING_2026-10.md), where the
MD on the corrected surface reached the GFN2 energy level but not the GFN2 torsion populations (JSD GFN-FF 0.268, A 0.291, B2 0.305, noise 0.198).

Question: do the models miss the torsion populations because they were trained on 400 K GFN-FF frames (not the target distribution), or because of the model class?

Design: NVT 300 K (csvr, 1 fs, 40 ps, MD seeds 1 to 4, first 5 ps dropped), 21 backbone torsions, every row against the same reference, the pooled GFN2 seeds 2 to 4 (the oracle
trains on seed 1). Noise floor (one GFN2 seed against the pool of the other three of seeds 1 to 4): JSD 0.198 +- 0.049 (the stricter definition seed 2 to 4 against the other two: 0.237 +- 0.052).
A = 1000 centres, B2 = 4000 centres; label check of V1 and V2 against `-sp -dump_gradient`: |dE| <= 7e-14 Eh, |dG| <= 1e-14 Eh/A.

| variant | labels | torsion JSD | seeds <= GFN-FF | W1 E_GFN2 [kcal/mol] | model energy minus GFN2 on own MD frames [kcal/mol] |
|---|---|---|---|---|---|
| GFN-FF | 0 | 0.271 +- 0.052 | 2/4 | 17.4 | - |
| V0 A (400 K GFN-FF data) | 351 | 0.301 +- 0.068 | 1/4 | 6.5 | -9 |
| V0 B2 | 351 | 0.321 +- 0.080 | 1/4 | 3.7 | -6.7 |
| V1 A (300 K GFN-FF data) | 351 | 0.299 +- 0.020 | 0/4 | 9.3 | -12.8 |
| V1 B2 | 351 | 0.273 +- 0.035 | 3/4 | 3.8 | -4.8 |
| V2 A (oracle: GFN2-MD frames) | 351 | 0.360 +- 0.072 | 1/4 | 25.1 | -41.5 |
| V2 B2 (oracle) | 351 | 0.269 +- 0.096 | 3/4 | 6.5 | -12.7 |
| V3a (active learning, default stop) | 87 | 0.318 +- 0.119 | 1/4 | 8.3 | -4.3 |
| V3b (active learning, 12 cycles, threshold 0.3) | 172 | 0.236 +- 0.028 | 4/4 | 6.6 | -3.0 |
| V3c (replicate) | 163 | 0.256 +- 0.028 | 3/4 | 5.8 | -2.5 |
| V3d (replicate) | 180 | 0.239 +- 0.034 | 3/4 | 6.7 | -1.3 |

## Findings (four MD seeds per model, one training draw for V0 to V2; the standard error of a JSD is about 0.04 to 0.07, differences below 0.05 are not resolved)

1. **Target-temperature data (V1) does not help:** A -0.002 +- 0.037 against V0 (2/4 seeds closer), B2 -0.048 +- 0.048 (3/4); V1 B2 (0.273) equals GFN-FF and is 1.4 times the noise floor.
2. **The oracle does not reach the noise floor:** V2 A is worse than GFN-FF (0.360, W1 E_GFN2 25 kcal/mol, model 41 kcal/mol below GFN2 on its own frames), V2 B2 is 0.269 (3/4 seeds below
   0.271, spread 0.185 to 0.402). Perfect coverage of the GFN2 frames therefore does not give the GFN2 torsion populations. Held-out error on GFN2 frames is small (relative-energy MAE
   1.6 to 1.9 kcal/mol against 2.82 for GFN-FF), but the MD on the learned surface runs into energy holes outside the training frames (model 13 to 41 kcal/mol below GFN2): the error is in the
   region that the sampler visits, not in the target distribution. At fixed data more centres shrink the holes (B2 -13 against A -41).
3. **Active learning helps, not to the noise floor:** with the default stop (V3a, 87 labels, "converged" at an uncertainty of 0.48) there is no gain (JSD 0.318): the uncertainty measure certifies coverage of the sampled region, not
   the distribution. With threshold 0.3 and 12 cycles (163 to 180 labels, about half of 351, three independent runs) the JSD is 0.236 / 0.256 / 0.239 (mean 0.244), each below GFN-FF and V0 B2 (paired
   against V0 B2 -0.066 to -0.085, closer in 3/4, 3/4, 4/4 seeds; against GFN-FF -0.015 to -0.035), still above the noise floor (0.198, 0.9 standard deviations of the noise). The holes are shallow (-1 to -4 kcal/mol), but none of the three runs converged and
   the uncertainty is not monotone (V3b 0.48 in cycle 10, 0.78 in 11, 1.00 in 12). The 1000-centre class is thus not the limit by itself. W1 of the energy (5.8 to 6.7) is worse than for V0 B2 and V1 B2 (3.7, 3.8).

**Conclusion:** the wrong torsion populations are not explained by the temperature or the origin of the training frames (the oracle fails); they follow the model energy holes in regions the corrected-surface MD itself visits, and data
collected where that MD goes (active learning) reduces them. A safe use of Delta-ML in MD therefore needs training data from the model's own trajectories, and a convergence criterion better than the present uncertainty indicator.

Observations on the tool: `-activelearn.seed 1` does not fix the sampling MD (same parameters give different uncertainty sequences, e.g. cycle 2: 7.19 against 5.55 sigma), the active learning is not reproducible; the "training set" MAE in
the log only covers the training data (0.02 rising to 0.4 kcal/mol with more data), a held-out measure is missing. The threshold used for the judgement (noise + 1 standard deviation = 0.247) and the choice of V3b to V3d (after V3a stopped early)
were made after seeing the data. Not tested: other molecules, longer active learning, B2-sized active-learning models, several start structures, enhanced sampling. The statement of the agent that a GFN2 label costs 5 ms is not
consistent with 30 ms per GFN2 MD step and is not used here.
