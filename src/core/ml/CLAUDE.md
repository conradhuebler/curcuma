# src/core/ml - Delta-Machine-Learning Core

Learns the energy difference between a cheap base method and an expensive reference:
`E_ref ~ E_base + Delta(R)`. Used to run conformer search and pre-optimisation at
GFN-FF cost with reference-level ranking. See [docs/DELTA_ML.md](../../../docs/DELTA_ML.md).

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
- **Graceful extrapolation**: far from training data the Gaussian kernel decays, so both
  Delta and its gradient vanish and the model falls back to the pure base method

## Verification status (AI-generated, machine-tested)

Automated: `ctest -R deltaml` (20 unit checks + 13 CLI assertions).

- Descriptor Jacobian vs. finite differences: h^2 convergence to the roundoff floor
  (1.5e-5 at h=1e-3 down to 6e-10 at h=1e-6); translation invariance 2e-16, rotation 1e-15
- Model gradient vs. finite differences: 1.2e-11 Eh/Angstrom
- Energy from the gradient path and the plain path agree bit-for-bit
- JSON round-trip reproduces energies and gradients exactly
- Energy-only fit interpolates its training energies to 1.5e-4 kcal/mol, which separates a
  regression defect from a mere energy/gradient weighting trade-off
- Real system (353 conformers, 90 atoms, C/H/N, GFN-FF -> GFN2): relative energy MAE
  0.647 -> 0.167 kcal/mol, gradient RMSE 5.04e-3 -> 0.90e-3 Eh/A, Spearman 0.41 -> 0.89
- **Not tested**: other molecule classes, charged or open-shell systems, transition metals,
  transferability across molecules, long-time MD stability, any DFT reference.
  Human production testing pending.

## Open

- Uncertainty is a descriptor-space distance, not an error estimate - adequate for steering
  active learning, not an accuracy guarantee
- `sigma_energy` / `sigma_force` set a genuine trade-off, not a formality: on a small rough
  set the 3N gradient rows outvote the single energy row and energies get traded away
  (measured: sigma_energy 1e-4 -> 3.01 kcal/mol at gradient 0.0185, 1e-6 -> 0.0039 at 0.0351)
- Training cost grows as M^2; centres are capped at 12000 with an explicit error
- GPU port of the kernel and descriptor evaluation not started
