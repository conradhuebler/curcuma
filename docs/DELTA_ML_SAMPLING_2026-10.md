# Delta-ML as a sampler: how close is the MD distribution to GFN2? (Oct 2026)

Status: 🤖 AI-generated scripts (`scripts/deltaml_usecase7/`), ⚙️ machine-tested, peo8 (58 atoms) and a shortened M1 (90 atoms), one start
structure per molecule, not human-reviewed. Measured by a Sonnet agent with the binary that contains the label repairs and the `-md.seed` fix.
Continues [DELTA_ML_PIPELINE_2026-10.md](DELTA_ML_PIPELINE_2026-10.md), where the energy conservation of the corrected surface was shown (n = 1).

Setup: `-md.seed N` (N > 0) seeds the initial velocities now (checked: seed 1 twice gives the same trajectory md5, seed 2 another one). Models A (1000 centres) and
B2 (4000 centres) trained on all 351 frames of the 400 K GFN-FF data set of peo8 (fit 8 s and 54 s; M1: 401 frames, 15 s and 97 s). GFN2 single points of
GFN2-MD frames equal the MD energy (4 frames, 2.5e-4 kcal/mol). The models were trained on 400 K GFN-FF data, the runs below are at 300 K and 500 K (not tested how far this
leaves the training range).

## Part A: peo8, NVT 300 K (csvr, 1 fs, 40 ps, seeds 1 to 4, first 5 ps dropped, 21 backbone torsions O-C-C-O and C-O-C-C)

Every trajectory against the pooled GFN2 reference (4 seeds); the GFN2 row is one GFN2 seed against the other three (noise floor). Mean +- std over 4 seeds. JSD =
Jensen-Shannon divergence of the torsion distributions (30 degree bins, bits, averaged over the torsions); W1 = Wasserstein-1 distance of the 1D distributions.

| | JSD torsions | W1 E_GFN2 [kcal/mol] | W1 Rg [A] | W1 end-to-end [A] | wall s per ps |
|---|---|---|---|---|---|
| GFN2 (noise) | 0.198 | 0.78 | 0.22 | 1.21 | 30.6 |
| GFN-FF | 0.268 | 17.3 | 0.49 | 2.31 | 0.2 |
| Delta-ML A | 0.291 | 6.4 | 0.50 | 1.90 | 2.0 |
| Delta-ML B2 | 0.305 | 3.6 | 0.17 | 1.49 | 3.5 |

Paired by seed, seeds out of 4 in which Delta-ML is closer to GFN2 than GFN-FF: E_GFN2 A 4/4, B2 4/4; Rg A 2/4, B2 4/4; end-to-end A 2/4, B2 3/4 (difference below 1 SE);
torsion JSD A 1/4, B2 2/4. A control with equal halves (2 seeds against the other 2, 6 splits) gives the torsion JSD 0.151 (GFN2), 0.213 (GFN-FF), 0.207 (A), 0.225 (B2).
No bond broke (largest stretch 1.14 to 1.17 times the start).

## Part B: peo8 at 500 K, minima found per wall time

40 ps, 4 seeds for GFN-FF and B2, one 10 ps run for GFN2; every picosecond quenched with GFN-FF, GFN2 single points of the minima (energies relative to the lowest of
all 330 minima):

| | frames | distinct minima | lowest GFN2 energy [kcal/mol] | minima within 2 kcal/mol | sampling wall time per minimum |
|---|---|---|---|---|---|
| GFN2 | 10 | 10 | 3.19 | 0 | 29.9 s |
| GFN-FF | 160 | 158 | 1.55 | 1 | 0.17 s |
| Delta-ML B2 | 160 | 160 | 0.00 | 2 | 3.4 s |

B2 leads in the lowest minimum only because of one seed (B2 minus GFN-FF -0.95 +- 1.25 kcal/mol, better in 2 of 4 seeds). Of the 330 quenched structures only 2 pairs meet the duplicate
criterion (RMSD < 0.3 A and |dE| < 0.2 kcal/mol; coarser criteria change nothing), so "distinct minima" is the number of frames and the minima are not saturated.

## Part C: M1 (90 atoms), NVT 300 K, 10 ps, 2 seeds only

W1 of E_GFN2: GFN-FF 21.0, A 12.0, B2 5.1 against a noise floor of 1.5 kcal/mol (A and B2 closer in 2 of 2 seeds); Rg closer for A and B2 in 2 of 2 seeds; torsions: A closer in
2 of 2 seeds, B2 worse than GFN-FF in 2 of 2 (the noise floor comes from one seed against the other).

## Conclusion

The corrected surface puts the MD at the GFN2 energy level (E_GFN2 of the frames 4/4 seeds closer than GFN-FF; B2 about 4.6 times the GFN2 noise floor) and, for B2, reproduces the
radius of gyration. It does **not** reproduce the conformer populations: torsion populations and the end-to-end distance are not closer to GFN2 than those of GFN-FF (A 1/4 and 2/4,
B2 2/4 and 3/4 of the seeds), and the sampling at 500 K does not find more or lower minima per wall time than GFN-FF (about 20 times the sampling time per minimum, comparable
lowest minima). The claim "Delta-ML MD samples closer to GFN2 than GFN-FF" holds for the energy level only. Cost per ps: B2 11 per cent and A 6.5 per cent of a GFN2 step on peo8
(M1 6 and 4 per cent), measured under load; the labelling of the 351 frames is not included.

Follow-up: [DELTA_ML_TRAINING_DATA_2026-10.md](DELTA_ML_TRAINING_DATA_2026-10.md) (not the training temperature or the origin of the frames; energy holes in regions visited by the model MD, active learning reduces them).

Not tested: longer runs (the end-to-end distance of GFN2 varies by about 0.5 A between seeds, 40 ps may not be equilibrated), several start structures, models trained on
data at the sampling temperature, enhanced sampling (metadynamics), charged and open-shell systems, GPU.
