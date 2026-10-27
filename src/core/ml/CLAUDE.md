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
- `ml_gpu.h/.cpp`, `ml_gpu_abi.h`, `cuda/ml_model_cuda.cu` - `GpuModel`: descriptor, kernel, forces and the fit's
  design block on the device (plugin `curcuma_cuda`, entry points `curcuma_cuda_ml_model_*`); [docs/DELTA_ML_GPU_2026-10.md](../../../docs/DELTA_ML_GPU_2026-10.md)

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
- **GPU path keeps the CPU definitions**: a change to `descriptor.cpp` (cutoff, grid, species handling) must be
  repeated in `cuda/ml_model_cuda.cu`; test `cli_deltaml_06_gpu_path` (skips with 77 without a device) checks parity
- **GPU never diverges silently**: unavailable backend -> one warning; device failure or > 512 neighbours -> error, CPU path

## Verification status (AI-generated, machine-tested)

Automated: `ctest -R deltaml` (unit test `test_deltaml` + CLI tests 01-06; 03 is the FD gate of the energy model, 04 checks labels against fresh `-sp`,
05 the `learn` loop: same seed -> byte-identical curve/dataset, validation columns, `al_criterion`; 06 GPU parity, FD gate and fit on the device).

- GPU (Oct 2026, RTX 5080, 58/90 atoms): correction 0.17/0.29 ms per MD step (old GPU 1.09/2.15, CPU 1.27/3.13), fit (90 atoms, M 2000) 2.1 s vs 7.5 s old GPU;
  vs CPU 25 frames each max |dE| 8.5e-14 Eh, |dg| 1.0e-14 Eh/A; FD gate 3.5e-8 Eh/A; not tested: ROCm/Vulkan, multi-GPU, charged systems, M > 2000

- `learn` reproducible since Oct 2026 (per-instance thermostat engine in SimpleMD for seed > 0, sorted trajectory reading): byte-identical at -threads 1/4/8 on 9 atoms and 4/8 on peo8
- `DeltaMLErrors::energy_mean_signed` / `energy_min_signed` (model minus reference): the latter is the "hole" of the `learn` validation (`al_validation_frames`, `al_criterion`, docs/DELTA_ML.md)
- Not tested: validation frames with `al_onthefly`, criterion defaults on molecules other than peo8 and dimethyl ether; at defaults peo8 did not converge in 12 cycles

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
- With `-gpu cuda` GFN-FF also moves to the GPU and is slower there at 58-90 atoms (0.67 vs 0.16 ms); `-gpu` is one global switch
