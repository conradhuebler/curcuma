# Delta-ML re-measured with repaired labels: forces, trade-off, spurious minima, new basins (Oct 2026)

Status: 🤖 AI-generated scripts (`scripts/deltaml_usecase4/`), ⚙️ machine-tested, two molecules, not human-reviewed. Measured by a
Sonnet agent with the binary that contains the label and GFN-FF repairs (commits bc7742e6, 373b1c13, 43786e04). This
replaces the numbers in DELTA_ML.md that were measured with the defective internal labelling.

Checks: internal labelling equals the externally computed fresh labels bit for bit (20 frames of M2: geometry, both
energies, both gradients, max difference 0.0, input order); model energy minus correction equals a fresh GFN-FF single
point to 9.6e-9 Eh at 892 minima. Sets: E1 = 353 conformers of M1 (GFN2 span 2.2 kcal/mol), M1 and M2 = 401 frames of a 400 K
MD each (M1 90 atoms, M2 trisaccharide 66 atoms).

## Part 1: held-out forces and energies (20 % held out, split seeds 1 to 5, joint model, default settings)

| set | rel. energy MAE model / GFN-FF [kcal/mol] | gradient RMSE model / GFN-FF [Eh/A] | factor | Spearman model / GFN-FF |
|---|---|---|---|---|
| E1 | 0.346 / 0.672 | 2.25e-3 / 9.53e-3 | 4.2x | 0.742 / 0.506 |
| M1 | 2.75 / 4.36 | 7.30e-3 / 12.3e-3 | 1.7x | 0.957 / 0.905 |
| M2 | 2.83 / 5.40 | 8.08e-3 / 14.7e-3 | 1.8x | 0.926 / 0.741 |

The model beats GFN-FF in energy, gradient and Spearman in 5 of 5 seeds on every set; fresh seeds 11 to 15 repeat it
(gradient factor 4.24 / 1.68 / 1.80). **The "4x better forces" holds on the narrow set only; on the MD sets it is 1.7 to
1.8x.** The energy-only model has the gradient of the base (factor 1.00 / 1.16 / 1.16): only the joint model improves
forces. Its energies are good (held-out MAE E1 0.234, M1 3.50, M2 3.36).

Parameter scan (one factor at a time, paired to the default, n = 5 seeds):
- **n_centers** is monotone and beyond the seed noise on all sets: 500 / 2000 / 4000 centres give MAE E1 0.488 / 0.287 /
  0.242 (default 0.346), M1 2.99 / 2.56 / 2.44, M2 3.23 / 2.64 / 2.37; fit time 12 s (1000) to 76 s (4000). On fresh seeds
  11 to 15 4000 centres beat the default in 5 of 5 seeds on all sets (gradient factor against GFN-FF 5.4 / 1.85 / 2.2).
- **lambda** 1e-4 and 1e-2 are worse in MAE and gradient on all sets; 1e-8 is within noise.
- **sigma_energy / sigma_force:** only the ratio acts (swapped pairs are bit-identical). **Trade-off:** on E1 it exists
  (energy weighted up: MAE 0.151, gradient 4.4e-3; force weighted up: MAE 1.08, worse than GFN-FF, gradient 1.6e-3); on
  the thermal sets there is none in the scanned range, the default ratio is best for MAE and both directions cost 0.27 to
  0.93 kcal/mol.
- **ACSF parameters:** effects below 0.17 kcal/mol and 3e-4 Eh/A, mostly within noise. kernel_sigma 0.5x and 2x are worse.

## Part 2: pre-optimisation and spurious minima (20 frames per molecule, one split, default optimiser)

GFN2 energy at the end point minus the GFN2 optimum from the same start (GFN-FF: 11.2 / 23.8 kcal/mol):

| model | M1 | M2 | paired to GFN-FF (better / n) | model energy below GFN2 by > 2 kcal/mol |
|---|---|---|---|---|
| A (lambda 1e-6) | 7.55 | 12.7 | 18/20, 17/19 | 7/20, 19/19 |
| S (lambda 1e-3) | 7.23 | 12.4 | 20/20, 20/20 | 0/20, 0/20 (model energy 7 to 16 kcal/mol too high) |
| B2 (4000 centres) | 5.70 | 8.16 | 20/20, 20/20 | 7/20, 20/20 |

All models reach lower GFN2 energies than GFN-FF. The model energy at its own optimum is not reliable for A and B2
(spurious lowering of up to 34 kcal/mol); lambda 1e-3 removes the lowering and biases the energy high instead. The uncertainty
abort (`-deltaml.uncertainty_abort` 0.5 and 0.9) does not remove the spurious minima (the indicator is small at the results,
median 0.08 to 0.12) and aborts some A runs. RMSD to the GFN2 optimum: M1 not better than GFN-FF (7 to 9 of 20), M2 better
(13 to 17 of 19 to 20).

## Part 3: new basins (6 MD runs per molecule from 6 different start frames, 3 splits of 4 training and 2 test runs)

67 to 80 GFN-FF minima per split, none close to a training minimum. Mean over 3 splits (M1 / M2):

| | Spearman | pair accuracy | regret [kcal/mol] | top-1 hits (of 3) |
|---|---|---|---|---|
| GFN-FF | 0.12 / 0.67 | 0.54 / 0.77 | 3.05 / 3.91 | 0 / 0 |
| A | 0.75 / 0.82 | 0.88 / 0.85 | 1.02 / 3.07 | 1 / 0 |
| B2 (4000 centres) | 0.81 / 0.82 | 0.91 / 0.85 | 0.78 / 2.64 | 1 / 1 |
| S | 0.38 / 0.75 | 0.68 / 0.82 | 2.20 / 4.36 | 0 / 0 |

Paired Spearman gain over GFN-FF in 3 of 3 splits for A and B2 on both molecules (M1 +0.63 and +0.69, M2 +0.14 and +0.15); regret
lower in 3 of 3 splits. The bootstrap interval over the minima excludes 0 in 3 of 3 splits on M1 and in 1 of 3 on M2. The lowest
conformer is found rarely (1 of 3 splits against 0 of 3); it is not shown to be found more reliably.

## Defect found on the way

`-md.seed` does not reach the initial velocities (`InitVelocities` uses an unseeded static `default_random_engine`): six runs with
seeds 21 to 26 from one start were bit-identical, and identical to an earlier run with seed 11. Part 3 therefore uses six different
start frames of an earlier MD, so the runs are not seed-independent.

## Conclusion (two molecules)

With repaired labels the joint model improves held-out energies, forces (1.7x to 4.2x, more with 4000 centres) and the order of
conformers over GFN-FF, also on minima from basins it was not trained on. Not shown: that it picks the lowest conformer more
reliably. Pre-optimisation lowers the GFN2 energy at the end point in the large majority of frames, and the model energy at its own
optimum must not be trusted for the default and the 4000-centre model. Cost (earlier measurement): 11 to 25 times GFN-FF per
step. Not tested: other molecules, a DFT reference, charged or open-shell systems, independent MD seeds, long MD with the
corrected labels, the energy/force trade-off at ratios beyond the scanned range. The seeds of Part 1 share frames, so the paired
n = 5 is only a screen; Part 2 uses one split per molecule; Part 3 splits share runs.
