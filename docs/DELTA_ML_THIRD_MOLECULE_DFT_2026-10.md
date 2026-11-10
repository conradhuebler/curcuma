# Delta-ML: third molecule and DFT reference (Oct 2026)

Status: 🤖 AI-generated scripts (`scripts/deltaml_usecase5/`), ⚙️ machine-tested, four molecules overall, not human-reviewed.
Measured by a Sonnet agent with the binary that contains the label and GFN-FF repairs. Continues
[DELTA_ML_RETEST_2026-10.md](DELTA_ML_RETEST_2026-10.md) (M1, M2, GFN2 reference).

Label check (5 structures per data set, fresh single point against internal label): GFN-FF and GFN2 labels bit-identical; ORCA labels
differ by at most 4.9e-8 Eh in energy and 3.2e-5 Eh/A in the gradient (fresh against fresh the ORCA result is bit-identical, so the label
depends on the SCF start, probably the `.gbw` of the previous structure in the worker directory; about 0.3 % of the held-out gradient
RMSE, not isolated). `-ref_method r2scan-3c` works directly (EnGrad TightSCF, charge 0, multiplicity 1; `/opt/orca` must be in PATH).
One caffeine single point with gradient takes 15.6 s at 8 threads; labelling 150 caffeine frames 19 min, 80 peo8 frames 18 min.

## Part A: third molecule peo8 (C16H34O8, 58 atoms, flexible ether chain, GFN2 reference)

351 frames of a 400 K MD, 70 held out, split seeds 1 to 5 (model against GFN-FF, 5 of 5 seeds better in every row):

| | model / GFN-FF | factor |
|---|---|---|
| relative energy MAE, default (1000 centres) | 2.47 / 4.62 kcal/mol | 1.9x |
| gradient RMSE, default | 7.04e-3 / 11.76e-3 Eh/A | 1.67x |
| Spearman, default | 0.940 / 0.786 | +0.154 |
| relative energy MAE, 4000 centres | 1.77 / 4.62 | 2.6x |
| gradient RMSE, 4000 centres | 5.34e-3 / 11.76e-3 | 2.20x |
| energy-only model: energy MAE / gradient RMSE | 2.61 / 4.62; 1.03e-2 / 1.18e-2 | 1.8x; 1.14x |

Minima in new basins (3 splits, 71 / 78 / 74 minima, none close to a training minimum), mean over splits:

| | Spearman | pair accuracy | top-1 (of 3) | recall@3 | regret [kcal/mol] | paired dSpearman |
|---|---|---|---|---|---|---|
| GFN-FF | 0.489 | 0.708 | 0 | 0.00 | 2.22 | - |
| A (default) | 0.704 | 0.813 | 1 | 0.22 | 1.77 | +0.215 (3/3) |
| C (energy model) | 0.737 | 0.835 | 1 | 0.33 | 1.28 | +0.248 (3/3) |
| S (lambda 1e-3) | 0.744 | 0.840 | 0 | 0.22 | 2.39 | +0.255 (3/3) |
| B2 (4000 centres) | 0.704 | 0.815 | 2 | 0.33 | 0.55 | +0.215 (3/3) |

The bootstrap interval of dSpearman excludes 0 in 1 of 3 splits for A and in 2 of 3 for C, S and B2; the sign is always positive.
Pre-optimisation (20 frames, one split, default optimiser): GFN2 energy at the end point minus the GFN2 optimum 14.17 (GFN-FF), 8.73 (A),
5.05 (S), 6.34 (B2) kcal/mol, better than GFN-FF in 20 of 20 frames each; the model energy at its own optimum is more than 2 kcal/mol
below GFN2 for A and B2 (19 of 20, minimum -15.5), never for S (but 8.5 too high). RMSD to the GFN2 optimum is not better than GFN-FF
(1.20 against 1.62 / 1.29 / 1.54 A). Same picture as M1 and M2.

## Part B: DFT reference r2SCAN-3c (gradients measured against the DFT gradients)

Caffeine, 150 frames from a 500 K MD (120 training, 30 test), split seeds 1 to 5:

| model | energy MAE model / GFN-FF [kcal/mol] | gradient RMSE model / GFN-FF [Eh/A] | Spearman model / GFN-FF |
|---|---|---|---|
| joint, default | 4.12 / 3.33 (worse, 0/5 better) | 1.16e-2 / 2.31e-2 (1.99x, 5/5) | 0.873 / 0.912 (worse, 0/5) |
| joint, 4000 centres (= all 2880) | 3.72 / 3.33 (1/5) | 0.82e-2 / 2.31e-2 (2.80x, 5/5) | 0.923 / 0.912 (3/5) |
| energy-only model | 2.36 / 3.33 (1.41x, 5/5) | 2.49e-2 / 2.31e-2 (worse) | 0.948 / 0.912 (+0.036, 5/5) |

Weighting the energy more (`sigma_energy` 1e-5) gives energy MAE 3.63 (1/5 better) at gradient factor 1.64x; at 1e-6 energy MAE 2.52 (5/5)
but the gradient factor falls to 1.00: an energy/force trade-off that the GFN2 thermal sets did not show.

With more data (501 caffeine frames from a second start, 400 training, 100 test, folds share frames in an autocorrelated MD): joint default
energy 3.12 / 3.75 (1.20x, 5/5), gradient 1.85x, Spearman 0.935 / 0.866; joint 4000 centres energy 3.05, gradient 2.59x, Spearman 0.950;
energy-only model energy 1.82 (2.06x, 5/5), Spearman 0.960. Independent MD run (one model and one seed per direction, no interval):
trained on 150 and tested on 351 frames the joint energy is 3.75 against 3.61 (no gain), gradient 1.64x, energy-only model 1.70x; trained
on 351 and tested on 150: 3.02 against 3.42 (1.13x), gradient 1.97x, energy-only model 2.21x.

peo8 with r2SCAN-3c (80 frames, 64 training and only 16 test per split, noisy): joint default energy 5.14 / 6.04 (1.18x, 4/5), gradient
9.2e-3 / 16.6e-3 (1.81x, 5/5), Spearman 0.721 / 0.627 (4/5); joint 4000 centres (= all 3712) energy 4.55, gradient 2.44x, Spearman 0.792;
energy-only model energy 4.85 (1.24x), gradient worse than GFN-FF; `sigma_energy` 1e-5 energy 4.15 (1.45x), gradient 1.40x.

## Conclusion

- **Third molecule (GFN2 reference):** the four findings of M1 and M2 hold for peo8: forces 1.67x (2.2x with 4000 centres), Spearman +0.154
  (5/5 seeds), minima of unseen basins better ordered (paired dSpearman +0.215, 3/3 splits), model energy at its own optimum too low for A and B2.
  Not shown for peo8 either: that the lowest conformer is found reliably (top-1 1 to 2 of 3 against 0 of 3), and the geometry after pre-optimisation
  is not closer to the GFN2 optimum.
- **DFT reference:** the forces transfer: 2.0x (default) to 2.8x (4000 centres) better than GFN-FF against the DFT gradient on caffeine and
  1.8x to 2.4x on peo8, in all 20 splits of the four DFT sets. Energy and order transfer with the energy-only model (caffeine 1.4x with 120 frames,
  2.1x with 400; peo8 1.2x) or an energy-weighted joint model; the default joint model is worse than GFN-FF in energy and order with 120 caffeine
  frames and needs about 400 frames to reach 1.2x. No setting gave energy and force better by about 2x at the same time. The earlier figure of
  2.3x (energy) and 2.1x (gradient) for caffeine in DELTA_ML.md (96 frames, measured with defective labels) is not reproduced for the joint model.
- **Not tested:** charged or open-shell systems, other functionals, long MD, independent MD seeds (`-md.seed` is ineffective), more than one MD run per DFT
  set. The DFT sets are small (30 and 16 test frames per split) and the seeds share frames; 4000 centres means all available environments there.
