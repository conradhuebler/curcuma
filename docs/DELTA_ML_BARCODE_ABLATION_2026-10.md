# Delta-ML: do persistent-homology barcodes add energy information beyond ACSF? (Oct 2026)

Status: 🤖 AI-generated tool, ⚙️ machine-tested. One molecule, energy only. Not human-reviewed.

## Question

Can a local Vietoris-Rips barcode embedding replace or complement the ACSF descriptor of Delta-ML
(`E_GFN2 - E_GFN-FF`)? Barcodes are permutation- and rotation-invariant but only C0 in the coordinates
(birth and death values are individual pair distances; the critical edge changes), so forces are not
usable. The test is therefore energy-only for every variant.

## Method (`src/helpers/ml_barcode_ablation.cpp`, target `ml_barcode_ablation`)

- Data: 353 conformers of one 90-atom C/H/N molecule, dataset JSON written by `-activelearn`
  (`conf.deltaml.dataset.json`); GFN-FF base, GFN2 reference.
- Per atom: patch of all atoms within 5 A; compressed distance matrix; Ripser (H0, H1) through
  `PersistentDiagram::generateBarcodes()`. Features: H0 death histogram (16 bins), H1 persistence image
  (8 x 8, weight = persistence), and a control with the same dimension (80): histogram of all pair distances in
  the patch ("RDF control"). A null variant shuffles the barcode rows over the conformers.
- Model for every variant: dense kernel ridge on structures with the summed element-blocked atomic Gaussian
  kernel (the dense, force-free limit of `DeltaMLModel`); features standardised per element and column;
  hyper-parameters (sigma scale, lambda) selected by inner 4-fold CV inside each training fold
  (5-fold outer, 3 repeats).
- Metric: centred MAE of the out-of-fold prediction error [kcal/mol]; Spearman of `E_base + prediction`
  against `E_ref`; 2000-fold bootstrap over conformers (paired for differences).

## Result (n = 353, one molecule)

| Variant | MAE | difference to ACSF [95 %] | Spearman |
|---|---|---|---|
| no model (E_base only) | 0.821 | - | 0.490 |
| ACSF | 0.615 | - | 0.509 |
| barcode H0+H1 | 0.654 | +0.039 [+0.010, +0.068] | 0.508 |
| pair-distance histogram (control) | 0.659 | +0.045 [+0.017, +0.071] | 0.467 |
| ACSF + barcode (equal weight) | 0.656 | +0.042 [+0.020, +0.064] | 0.506 |
| ACSF + 0.25 barcode | 0.631 | +0.017 [+0.005, +0.028] | 0.501 |
| ACSF + 0.25 histogram | 0.619 | +0.004 [-0.007, +0.015] | 0.499 |
| ACSF + barcode, rows shuffled (null) | 0.707 | +0.092 [+0.058, +0.127] | 0.457 |

Reading: barcodes are not empty (better than the shuffled null) but not better than the pair-distance
histogram, and no admixture improves ACSF. The hypothesis "barcodes add energy information" is refuted for
this design.

## Not tested

Other molecules or transferability, patch radii other than 5 A, H2, element-resolved filtrations, charged
or open-shell systems. Bootstrap intervals treat the conformers as independent; they are correlated, so the
intervals are too narrow. No force test: the barcode features have no analytic gradient.

## Reproduce

```
curcuma -activelearn conf.xyz -method gfnff -ref_method gfn2 -no_bmt     # writes conf.deltaml.dataset.json
ml_barcode_ablation conf.deltaml.dataset.json --threads 24 --json result.json
```

## Side findings

- `external/ripser/ripser.h` has four function-local `static` enumerators bound to the first ripser object;
  concurrent calls crash (SIGSEGV, gdb), serial calls work only by stack-address luck. `PersistentDiagram` now
  serialises the Ripser call with a mutex and a fixed construction address (`runRipser()`), which is safe but
  gives no speed-up. A real fix removes the `static` in the Ripser fork.
- Fable's design review (lead: smooth Hodge-Laplace spectra instead of barcodes if topology is wanted at all)
  was not pursued: it estimated a 15 to 25 % chance of a robust gain, and the measurement above is negative.
