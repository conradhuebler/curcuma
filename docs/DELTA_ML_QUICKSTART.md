# Delta-ML quickstart (experimental)

Status: 🤖 AI-generated, ⚙️ machine-tested, **experimental**; human production testing pending. Details: [DELTA_ML.md](DELTA_ML.md).

**What it is.** `-activelearn` fits a sparse kernel model (ACSF descriptor) to `Delta = E_ref - E_base` for one molecule, from
energies and gradients of a cheap base (GFN-FF) and a reference (GFN2 or DFT). `-method deltaml:<model.json>` then adds the
correction wherever a method is accepted (`-sp`, `-opt`, `-md`). It is meant for one molecule at fixed topology, never for reactions.

## What the measurements say (GFN2 reference unless noted; four molecules: M1 90 atoms, M2 trisaccharide 66, peo8 58, caffeine 24)

| claim | verdict | evidence |
|---|---|---|
| held-out forces 1.7x-4.2x better than GFN-FF (1.8x-2.8x against r2SCAN-3c) | shown, 5/5 split seeds per set | [RETEST Part 1](DELTA_ML_RETEST_2026-10.md), [DFT Part B](DELTA_ML_THIRD_MOLECULE_DFT_2026-10.md) |
| conformer minima (also from unseen basins) ordered better than by GFN-FF | shown, Spearman gain in 3/3 splits, M1, M2, peo8 | [RETEST Part 3](DELTA_ML_RETEST_2026-10.md), [THIRD_MOLECULE Part A](DELTA_ML_THIRD_MOLECULE_DFT_2026-10.md) |
| NVE energy conservation as good as GFN-FF | shown, n=1 trajectory per molecule (M1, peo8) | [PIPELINE Part 2](DELTA_ML_PIPELINE_2026-10.md) |
| MD frames at the GFN2 energy level | shown, 4/4 MD seeds (peo8) | [SAMPLING Part A](DELTA_ML_SAMPLING_2026-10.md) |
| lowest conformer found more reliably | not shown (top-1 in 1-2 of 3 splits vs 0 of 3) | [RETEST Part 3](DELTA_ML_RETEST_2026-10.md) |
| pre-optimisation before GFN2 beats a GFN-FF pre-optimisation | not shown (GFN-FF as fast or faster on 2 of 3 molecules) | [PIPELINE Part 1](DELTA_ML_PIPELINE_2026-10.md) |
| MD torsion populations closer to GFN2 than GFN-FF | not shown (JSD A 1/4, B2 2/4 seeds) | [SAMPLING](DELTA_ML_SAMPLING_2026-10.md), [TRAINING_DATA](DELTA_ML_TRAINING_DATA_2026-10.md) |
| model energy at its own optimum is reliable | not shown (up to 34 kcal/mol spurious lowering) | [RETEST Part 2](DELTA_ML_RETEST_2026-10.md) |
| reactions, topology changes | not tested (argued to fail) | [DELTA_ML.md, "When this applies"](DELTA_ML.md) |
| charged or open-shell systems | not tested | [RETEST, conclusion](DELTA_ML_RETEST_2026-10.md) |
| molecules other than the four above, transfer between molecules | not tested | [THIRD_MOLECULE](DELTA_ML_THIRD_MOLECULE_DFT_2026-10.md) |
| GPU beyond fit/step timing | not tested | [DELTA_ML.md, "GPU"](DELTA_ML.md) |

## Worked example: n-butane, GFN-FF -> GFN2 (run 27 Oct 2026, n=1)

Binary `build_deltaml/curcuma` (source = HEAD 483e3804; banner prints 75d03431, built from the then uncommitted tree), 32 cores, input
`test_cases/structures/medium/butane.anti.xyz` (CCCC dihedral 116.4 deg, so only an MD start). Training: 30 ps GFN-FF MD at 500 K, a frame per 100 fs
(302 frames, 167 with |CCCC| > 120 deg); fresh test: second MD, seed 2, 10 ps (102 frames, 95 with |CCCC| > 120 deg). GFN2 is cheap at 14 atoms.

```sh
curcuma -md butane.xyz -method gfnff -max_time 30000 -time_step 1 -temperature 500 \
        -dump_frequency 100 -md.seed 1 -no_bmt                                   # 1.4 s
cp butane.snapshots/butane.trj.xyz frames.xyz                                    # (test set: seed 2, 10 ps)
curcuma -activelearn frames.xyz -method gfnff -ref_method gfn2 -threads 8 -no_bmt  # fit, 1.5 s, exit 0
curcuma -activelearn test_frames.xyz -al_mode test -model frames.deltaml.json \
        -method gfnff -ref_method gfn2 -threads 8 -no_bmt                        # 0.2 s, exit 0
curcuma -sp  x.xyz -method deltaml:frames.deltaml.json -no_bmt                     # 0.03-0.05 s
curcuma -opt x.xyz -method deltaml:frames.deltaml.json -threads 1 -no_bmt          # 0.05 s
```

Fit log (held-out, 60 of 302 frames) and `-al_mode test` on the 102 fresh frames (joint model):

```
[RESULT]  held-out set (joint model: energies + gradients, for -opt/-md)  (60 structures)
[RESULT]    relative energy MAE : model    0.3304 kcal/mol   base    1.5313 kcal/mol
[RESULT]    gradient RMSE       : model  0.003306 Eh/A       base  0.009492 Eh/A
--- -al_mode test ---
[RESULT]  test set  (102 structures)
[RESULT]    relative energy MAE : model    0.3982 kcal/mol   base    1.6161 kcal/mol
[RESULT]    gradient RMSE       : model  0.003313 Eh/A       base  0.009295 Eh/A
[RESULT]    rank correlation    : model    0.9950            base    0.9477
```

`-opt` from one anti and one gauche test frame with each method, then a GFN2 `-sp` at every end point (kcal/mol):

| | GFN-FF | Delta-ML joint | Delta-ML energy model | GFN2 |
|---|---|---|---|---|
| fresh test frames: rel. energy MAE / gradient RMSE [Eh/A] | 1.616 / 9.30e-3 | 0.398 / 3.31e-3 | 0.604 / 6.80e-3 | - |
| E(gauche) - E(anti) at the method's own optima | 0.45 | 0.57 | 1.07 (at the GFN2 optima) | 0.60 |
| GFN2 energy at end point above GFN2 optimum, anti / gauche | 0.49 / 0.52 | 0.004 / 0.085 | - | 0 |
| gauche CCCC dihedral [deg] | 65.6 | 61.9 | - | 67.5 |
| model energy minus GFN2 at its own optimum, anti / gauche | - | -0.16 / -0.26 | - | - |

Read plainly: on this molecule the joint model lifts forces 2.8x and lands closer to the GFN2 minimum in energy, but its gauche
dihedral is further from GFN2 than GFN-FF's, and the energy model gives a worse gauche/anti gap than GFN-FF (one structure pair). Test
frames come from the same kind of MD as the training frames (in distribution). One molecule, one seed per MD: n=1.

## Single-molecule active learning (`-al_mode learn`)

```sh
curcuma -activelearn mol.xyz -al_mode learn -method gfnff -ref_method gfn2 -threads 8
```

Defaults in `src/capabilities/activelearn.h`: 5 cycles (`al_cycles`), each 4 MD runs (`al_md_runs`) of 2000 fs (`al_md_time`) at 500 K
(`al_temperature`) on the current corrected surface, every 50th frame a candidate (`al_dump_frequency`), 15 labels per cycle
(`al_candidates`, `al_initial` in cycle 0) plus 10 random validation frames (`al_validation_frames`), on-the-fly off (`al_onthefly`).
Stop rule `al_criterion`:
- `uncertainty`: no sampled frame above the indicator threshold `al_uncertainty` (0.5). The indicator is a distance in descriptor space, not an error.
- `validation`: the validation frames of the last `al_window` (3) cycles, re-predicted by the model entering the cycle and kept out of
  training while pooled, have a centred energy MAE <= `al_target_energy` (2.5 kcal/mol) and a hole (the `1 - al_hole_quantile` = lowest
  10 % of centred signed errors) not below `-al_target_hole` (-3.5 kcal/mol), in `al_patience` (2) consecutive cycles.
- `both` (**default**): `validation` as above and the indicator in the stopping cycle. Changed in 483e3804; before it, `both` needed the
  indicator in `al_patience` consecutive cycles too, which never happened on peo8 (0 of 3 runs).

Measured on peo8 only (n=3 per variant, [ACTIVE_LEARNING](DELTA_ML_ACTIVE_LEARNING_2026-10.md)): the new `both` stopped in 2 of 3 runs at
cycles 18 and 19 with 362 to 386 labels; `uncertainty` (threshold 0.3, not the default 0.5) stopped at 285 to 330 labels with a mean
torsion JSD of 0.235 against 0.221, a difference not resolved: equal quality at lower cost. The default was not changed. Butane
run with all defaults (n=1, 4.6 s): not converged after 5 cycles (validation held in cycles 4-5, indicator 0.84 > 0.5 in cycle 5),
85 labels; on the 102 fresh frames 0.510 kcal/mol / 3.80e-3 Eh/A, i.e. worse than the 302-frame fit above.

## Common failures

- **Dispersion-free reference**: a DFT reference without D3/D4 puts long-range dispersion into `Delta`, which a local descriptor
  cannot learn; results look plausible and are wrong ([DELTA_ML.md](DELTA_ML.md)). r2SCAN-3c carries D4.
- **Topology change** (bond breaking, proton transfer): the base keeps the anchor bond graph, the reference does not; `Delta` becomes steep.
- **Stale binary in the ctests**: `test_cases/cli/test_utils.sh` takes `release/curcuma` first. On 27 Oct 2026 `release/curcuma`
  (14 Oct) did not contain `al_hole_quantile`; run e.g. `CURCUMA=$PWD/build_deltaml/curcuma ctest -R deltaml`.
- **Units**: `getGradient()` is Eh/Angstrom; labels once multiplied it by Angstrom-to-Bohr (1.89x, fixed in bc7742e6). Results made
  with older internal labels are void ([USECASE_TEST2](DELTA_ML_USECASE_TEST2_2026-10.md)). Check new gradients with an FD of `-sp` energies.
- **One model per molecule**: a model applied to another molecule is untested; `-deltaml.base` different from the stored base gives
  silently wrong energies. Check every model on fresh frames with `-al_mode test` before use.
- **Joint vs energy model**: `.deltaml.json` for `-opt`/`-md`; `.deltaml.energy.json` only for `-sp` energies, its forces are no better than the base.

## Which document answers what

| question | document | scripts |
|---|---|---|
| first ranking/pre-opt test (WITHDRAWN in part: Tests 1-2 used defective labels) | [USECASE_TEST](DELTA_ML_USECASE_TEST_2026-10.md) | `scripts/deltaml_usecase/` |
| same with consistent labels, two molecules | [USECASE_TEST2](DELTA_ML_USECASE_TEST2_2026-10.md) | `scripts/deltaml_usecase2/`, `scripts/deltaml_usecase3/` |
| forces, parameters, spurious minima, new basins | [RETEST](DELTA_ML_RETEST_2026-10.md) | `scripts/deltaml_usecase4/` |
| third molecule, DFT reference | [THIRD_MOLECULE_DFT](DELTA_ML_THIRD_MOLECULE_DFT_2026-10.md) | `scripts/deltaml_usecase5/` |
| pre-optimisation before GFN2, NVE stability | [PIPELINE](DELTA_ML_PIPELINE_2026-10.md) | `scripts/deltaml_usecase6/` |
| MD distribution vs GFN2 | [SAMPLING](DELTA_ML_SAMPLING_2026-10.md) | `scripts/deltaml_usecase7/` |
| training data or model class? | [TRAINING_DATA](DELTA_ML_TRAINING_DATA_2026-10.md) | `scripts/deltaml_usecase8/` |
| active learning with validation, stop criteria | [ACTIVE_LEARNING](DELTA_ML_ACTIVE_LEARNING_2026-10.md) | `scripts/deltaml_usecase9/`, `scripts/deltaml_usecase10/` |
| barcode descriptor instead of ACSF | [BARCODE_ABLATION](DELTA_ML_BARCODE_ABLATION_2026-10.md) | `src/helpers/ml_barcode_ablation.cpp` |

Call order of the scripts: [scripts/DELTA_ML_SCRIPTS.md](../scripts/DELTA_ML_SCRIPTS.md). The scripts have hard-coded paths of the original runs.
