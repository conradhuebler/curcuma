# src/core/ml - Delta-Machine-Learning Core

Learns the energy difference between a cheap base method and an expensive reference:
`E_ref ~ E_base + Delta(R)`. Intended for conformer search and pre-optimisation near reference quality;
measured on one molecule: forces 4x better on held-out conformers, but ranking and pre-optimisation no better than GFN-FF, cost 11 to 25x GFN-FF (application test, docs/DELTA_ML_USECASE_TEST_2026-10.md). See [docs/DELTA_ML.md](../../../docs/DELTA_ML.md).

## Files

- `descriptor.h/.cpp` - ACSF (Behler G2 radial + G4 angular), element-resolved, with
  analytic position derivatives in compact per-atom Jacobian blocks
- `kernel.h` - Gaussian kernel; distances via GEMM so a GPU port is a BLAS swap
- `ml_dataset.h/.cpp` - `MLStructure` (geometry + base/reference energy and gradient),
  `MLDataset` with JSON persistence and a deterministic train/test split
- `krr_model.h/.cpp` - `DeltaMLModel`: sparse kernel ridge (subset of regressors, GAP-style)

## Key design points

- **Units**: geometry Angstrom, energy Hartree, gradient Eh/Angstrom (curcuma convention)
- **`n_centers` (M) is the accuracy/cost lever** - prediction costs O(M), independent of
  dataset size; whether large M pays off depends on the reference/base cost ratio
- **Energies and gradients are trained jointly** - one reference point yields 1 energy but
  3N gradient components, which is what makes a few hundred points enough
- **Element-blocked kernel** - environments of different elements never mix
- **Element-count baseline** absorbs the constant offset between two energy zeros; it has
  zero gradient and therefore never perturbs forces
- **C1 cutoff function** is mandatory, not cosmetic: a kink in the descriptor becomes a
  force discontinuity and shows up as NVE energy drift
- **ACSF grid is configurable** via `ACSFSettings` (`DeltaMLTrainingOptions::acsf`, CLI
  `-activelearn.acsf_*`); default 6 shells, eta_r 4 since Oct 2026 (was 8 / 8), keep `ACSFSettings` and the PARAMs equal
- **Model files carry their own descriptor block**: an old model reloads with its old grid, unaffected by default changes
- **`role`** in the model file: `joint` (energies + gradients) or `energy` (`use_forces` false); missing = `joint`;
  `DeltaMLMethod` warns once per process when gradients are taken from an `energy` model
- **Graceful extrapolation**: far from training data the Gaussian kernel decays, so both
  Delta and its gradient vanish and the model falls back to the pure base method

## Verification status (AI-generated, machine-tested)

Automated: `ctest -R deltaml` (unit test `test_deltaml` + CLI tests 01-04; 03 is the FD gate of the energy model, 04 checks labels against fresh `-sp`).

- Descriptor Jacobian vs. finite differences: h^2 convergence to the roundoff floor
  (1.5e-5 at h=1e-3 down to 6e-10 at h=1e-6); translation invariance 2e-16, rotation 1e-15
- Model gradient vs. finite differences: 1.2e-11 Eh/Angstrom
- Application test (one molecule, Oct 2026, docs/DELTA_ML_USECASE_TEST_2026-10.md): joint models worse than GFN-FF in ranking
  (paired Spearman -0.146 +- 0.095 and -0.099 +- 0.112, n=10 folds) and pre-optimisation (GFN2 energy at the result +59 and +75 kcal/mol, n=20);
  cost 11x and 25x GFN-FF; the stored base energy differs from fresh GFN-FF by 1.2 +- 0.7 kcal/mol (n=353), so numbers from
  `-activelearn` and from `-sp` are not in the same frame
- Labels (fixed Oct 2026): fresh calculators per structure, anchor bonds as in `DeltaMLMethod`, dataset in input order;
  stored energy/gradient labels equal `-sp -dump_gradient` exactly on 2 x 401 frames, `-threads 1` == `-threads 4` byte for byte.
  Cause of the old offset: a reused GFNFF keeps stale EEQ distances (docs/DELTA_ML.md, "Open: the anchor freezes more than it should").
  Datasets labelled before this fix, and every number measured on them, carry the offset; the seed now gives a different split
- Energy from the gradient path and the plain path agree bit-for-bit
- JSON round-trip reproduces energies and gradients exactly
- Energy-only fit interpolates its training energies to 1.5e-4 kcal/mol, which separates a
  regression defect from a mere energy/gradient weighting trade-off
- Real system (353 conformers, 90 atoms, C/H/N, GFN-FF -> GFN2), held out (71 structures, split seeds):
  old grid (8 shells, eta 8, M 1000, lambda 1e-6): rel. energy MAE 1.154 +- 0.054 vs base 0.894 +- 0.050
  kcal/mol (n=5, paired +0.260 +- 0.059, i.e. worse than the base), gradient RMSE 4.50e-3 vs 18.0e-3 Eh/A
- New grid (6 shells, eta 4): energy paired -0.215 +- 0.100 vs old grid (n=10, 10/10 better), still
  +0.110 +- 0.095 above the base; energy-only fit (`energy` model, lambda 1e-4): 0.623 +- 0.028, base -0.191 +- 0.050 (n=10)
- NOT REPRODUCED: an earlier line here ("relative energy MAE 0.647 -> 0.167, Spearman 0.41 -> 0.89");
  its protocol is unknown and the held-out numbers above do not confirm it
- **Not tested**: other molecule classes, charged or open-shell systems, transition metals,
  transferability across molecules, long-time MD stability, any DFT reference.
  Human production testing pending.

## Open

- Uncertainty is a descriptor-space distance, not an error estimate - adequate for steering
  active learning, not an accuracy guarantee
- `sigma_energy` / `sigma_force` set a genuine trade-off, not a formality: on a small rough
  set the 3N gradient rows outvote the single energy row and energies get traded away
  (measured: sigma_energy 1e-4 -> 3.01 kcal/mol at gradient 0.0185, 1e-6 -> 0.0039 at 0.0351)
- Energy model: smaller relative energy error, rank correlation not better than GFN-FF (Spearman 0.467 vs ~0.50, n=10 seeds)
- Training cost grows as M^2; centres are capped at 12000 with an explicit error
- GPU port of the kernel and descriptor evaluation not started
