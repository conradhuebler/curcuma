# Delta-ML: cheap method with a learned correction

**Status**: 🤖 AI-generated, ⚙️ machine-tested. Human production testing pending.


Learns the difference between a cheap base method and an expensive reference:

```
E_ref(R) ~ E_base(R) + Delta(R)
```

The gain is statistical. `Var(Delta)` is far smaller than `Var(E_ref)` and `Delta` is a
smoother function of geometry, so a few hundred reference points suffice where direct
learning would need tens of thousands. The intended purpose is conformer work: run conformer
search and pre-optimisation at GFN-FF cost while ranking conformers at the reference level.
**This purpose is not met on the one molecule tested** (ranking and pre-optimisation no better than
GFN-FF, cost 11 to 25x GFN-FF): see [DELTA_ML_USECASE_TEST_2026-10.md](DELTA_ML_USECASE_TEST_2026-10.md).

## When this applies, and when it does not

GFN-FF freezes a topology at setup: bond list, hybridisation, pi system, EEQ charges. Within
one topology `E_GFNFF(R)` is smooth and so is `Delta(R)` -- that is the regime this is built
for. As soon as the topology would change (bond breaking, proton transfer, rearrangement)
GFN-FF does not follow while the reference does, `Delta` becomes steep and effectively
unlearnable. **Conformers at fixed topology: yes. Reactions: no** -- use GFN2 as the base
method there instead.

A model is fitted **per molecule**. The descriptors are element-resolved, so nothing in the
architecture prevents transfer across molecules, but that is not validated.

## The reference must be dispersion-corrected

This is a hard requirement, not a preference. GFN-FF and GFN2 both carry D4 internally; a
DFT reference must have D3/D4 added. Then the dispersion contribution to `Delta` is a
difference of two D4 terms with *identical* C6 coefficients and only different damping. With
BJ damping both go to `-C6/r^6` at long range, so their difference decays as **r^-12** and
sits comfortably inside the descriptor cutoff.

Pair a dispersion-free reference with a D4-corrected base and the entire long-range
dispersion term lands in `Delta`, where a local descriptor cannot represent it. The result
will look plausible and be wrong.

The residual long-range term is electrostatic (EEQ point charges against the true density,
monopole differences decaying as r^-1), not dispersive. Small for a neutral molecule at
fixed topology, but it is the honest remainder.

## Usage

```bash
# 1a. Fit a model against an existing conformer set or trajectory
curcuma -activelearn conf.xyz -method gfnff -ref_method gfn2 -threads 8

# 1b. Or let it generate its own structures (see the caveats below)
curcuma -activelearn mol.xyz -al_mode learn -method gfnff -ref_method gfn2 \
        -al_cycles 6 -al_temperature 500 -threads 8

# 2. Use it anywhere a method is accepted
curcuma -sp        mol.xyz -method deltaml:conf.deltaml.json
curcuma -opt       mol.xyz -method deltaml:conf.deltaml.json
curcuma -confsearch mol.xyz -method deltaml:conf.deltaml.json

# Measure accuracy against evaluation cost before committing to a model size
curcuma -activelearn conf.xyz -center_sweep 200,500,1000,2000,4000

# Check an existing model against fresh reference calculations
curcuma -activelearn other.xyz -al_mode test -model conf.deltaml.json
```

Two preset configurations ship with the source and carry the measured reasoning in their
comments; the command line still overrides anything in them:

```bash
curcuma -activelearn conformers.xyz -import_config presets/deltaml_conformer.json
curcuma -activelearn start.xyz      -import_config presets/deltaml_trajectory.json
```

`-method deltaml -deltaml.model <file>` is equivalent to the `deltaml:<file>` shorthand. The
shorthand exists because several capabilities rebuild their sub-configuration from the
method name alone and would otherwise drop the `deltaml` block on the way down.

## How the model works

**Descriptor** -- ACSF (Behler & Parrinello): radial G2 functions are smeared distance
histograms per neighbour element, angular G4 functions do the same for angles. Invariant
under translation, rotation and permutation by construction. The cutoff function
`fc(R) = 0.5(cos(pi R/Rc) + 1)` vanishes together with its first derivative, which is
required rather than decorative: a kink in the descriptor becomes a force discontinuity and
shows up as energy drift in NVE dynamics.

**Regression** -- sparse kernel ridge (subset of regressors, as in GAP):

```
Delta(A) = sum_Z n_Z c_Z  +  sum_{i in A} sum_{m=1..M} beta_m k(x_i, z_m)
```

The first term is an element-count baseline absorbing the large constant offset between two
methods with different energy zeros; it has zero gradient and never touches forces. The
second is the learned part over `M` representative atomic environments picked by
farthest-point sampling. The kernel is element-blocked -- environments of different elements
never mix.

**M is the accuracy/cost lever.** A prediction costs O(M) regardless of how much training
data exists. Whether a large M pays off depends only on how expensive the reference is
relative to the base method, which is why `-center_sweep` exists rather than a hard-coded
default.

**Energies and gradients are trained jointly.** One reference calculation yields a single
energy but 3N gradient components; that ratio is what makes a few hundred points enough,
and gradients are required for optimisation and dynamics anyway.

**Extrapolation is graceful.** Far outside the training region the Gaussian kernel decays,
so `Delta` and its gradient both go to zero and the model falls back continuously to the bare
base method. A directly-learned force field has no such safety net.

## Active learning (`-al_mode learn`)

Each cycle runs molecular dynamics at `al_temperature` on the currently corrected surface,
picks the least familiar sampled structures, computes references for them, and refits.
Selection is greedy farthest-point in descriptor space among the candidates above the
uncertainty threshold: maximising uncertainty alone piles every pick onto one pathological
region and buys almost no new information for the same reference cost.

**On-the-fly (`-al_onthefly`).** Each sampling run is given the extrapolation threshold as
an abort condition, so it stops at the first structure that leaves the trained region and
that structure -- rather than a fixed-length batch -- becomes the candidate. The trajectory
nominates its own labels, and the model is refitted before sampling continues. Measured on
a 9-atom test: with the threshold at 0.5 a cycle collected 4 frames instead of 44, three of
which had tripped the threshold.

**The known weakness is structural.** Sampling happens on the surface being learned, so
cycle 0 only visits what the base method considers relevant. If the reference surface has a
minimum elsewhere, that region is never sampled and never learned, and the model confirms
itself. Elevated temperature and metadynamics mitigate this; they do not remove it. The
alternative — triggering a reference calculation the moment the uncertainty rises during
the dynamics, rather than in batches afterwards — is more robust and not implemented here.

**Status: the loop runs, its benefit is not demonstrated.** On a 9-atom test system it
completed its cycles, grew the dataset and refitted each time, but with only a handful of
held-out structures the error estimates are noise. The offline path (above) is the one with
measured results.

## Uncertainty

The reported indicator is `max_i (1 - max_m k(x_i, z_m))`: zero when every atomic environment
sits on a kernel centre, one when the structure is entirely outside the covered region.

Mind the scale. With a Gaussian kernel a nearest centre at distance `d` gives
`1 - exp(-d^2/2 sigma^2)`, so the indicator saturates fast: `d = 1.18 sigma` already reads
0.5, and `d = 2.1 sigma` reads 0.9. Thresholds below roughly 0.2 are effectively unreachable.
The log therefore also prints the corresponding distance in units of the kernel width, which
stays informative where the bounded indicator has flattened out.

It measures **distance in descriptor space, not error**. It answers "have I been here
before", never "am I right here". Adequate for steering active learning and for flagging
suspicious predictions; useless as an accuracy guarantee.

## Parameters

Module `activelearn` (fitting) and `deltaml` (evaluation). `curcuma -activelearn` without
arguments prints the summary; the full list comes from the parameter registry.

> The defaults of `lambda`, `kernel_sigma`, `sigma_*`, `energy_model_lambda` and the ACSF grid were
> chosen from scans on the defective labels described under "Held-out accuracy"; they are to be
> re-measured ([DELTA_ML_USECASE_TEST2_2026-10.md](DELTA_ML_USECASE_TEST2_2026-10.md)).

| Parameter | Default | Effect |
|---|---|---|
| `-ref_method` | `gfn2` | the level being learned; must match the base method's dispersion treatment |
| `-n_centers` | 1000 | kernel centres M; the accuracy/cost lever |
| `-kernel_sigma` | 0 (auto) | Gaussian width; strong lever, too narrow cannot interpolate, too wide oversmooths |
| `-activelearn.lambda` | 1e-6 | regulariser relative to the data term; raise to 1e-2 if NVE drifts. The flat `-lambda` does not route -- use the dotted form |
| `-jitter_relative` | 1e-8 | conditioning ridge on the normal equations; see below |
| `-use_forces` | true | train on gradients as well as energies |
| `-test_fraction` | 0.2 | held-out share, split randomly because neighbouring frames are correlated |
| `-center_sweep` | -- | comma-separated M values to fit and compare in one run |
| `-activelearn.dataset` | -- | train on a stored `<basename>.deltaml.dataset.json` instead of labelling again (`train` mode); the input file still names the outputs and supplies the anchor |
| `-activelearn.acsf_n_shell` | 6 (8 before Oct 2026) | radial shells, equidistant from `acsf_r_first` (0.9) to `acsf_r_last` (5.2 A); changes the feature count |
| `-activelearn.acsf_radial_eta` | 4.0 (8.0 before Oct 2026) | radial Gaussian width parameter [1/A^2] |
| `-activelearn.acsf_rcut_radial` / `acsf_rcut_angular` | 6.0 / 4.0 | cutoffs [A] |
| `-activelearn.acsf_angular_eta` | 0.05 | angular radial damping [1/A^2] |
| `-activelearn.acsf_zeta` | `1,4,16` | angular exponents, each with lambda = +1 and -1 |
| `-activelearn.energy_model` | true | also fit an energy-only model (`use_forces` false, same data, split and descriptor) and write `<basename>.deltaml.energy.json` with `"role": "energy"`; `train` mode only (also with `-activelearn.dataset`), not in `learn` mode or the centre sweep |
| `-activelearn.energy_model_n_centers` | 0 | centres of the energy model; 0 = `n_centers` |
| `-activelearn.energy_model_lambda` | 1e-4 | regulariser of the energy model (see "Held-out accuracy") |
| `-activelearn.energy_model_kernel_sigma` | 0 (auto) | kernel width of the energy model; 0 = median heuristic |

The `acsf_*` values apply to every fit of a run (all cycles in `learn` mode, the centre sweep)
and are stored in the model's `descriptor` block, which `-method deltaml:` reads back. A
model written with the pre-October-2026 grid (8 shells, eta 8) therefore keeps it: the old
`conf.deltaml.json` of the 90-atom set gave the same correction (-102.9052269203 Eh) before and
after the default change. Hyperparameter study without new reference
calculations: `curcuma -activelearn conf.xyz -activelearn.dataset conf.deltaml.dataset.json
-activelearn.acsf_n_shell 12 -no_bmt`. Use the dotted form; a flat flag that does not route is
ignored silently.

## Held-out accuracy (October 2026)

> **Measured with defective labels; to be re-measured.** Every number in this section (and the
> hyperparameter scans, "forces 4x better" and the energy/force trade-off elsewhere in this file)
> was fitted on datasets whose gradient labels were 1.89 times too large (fixed in bc7742e6) and
> whose stored base energies were not the GFN-FF energies of the structures (fixed Oct 2026, see
> "Open: the anchor freezes more than it should"). Kept for the record; current measurements:
> [DELTA_ML_USECASE_TEST2_2026-10.md](DELTA_ML_USECASE_TEST2_2026-10.md).

353 conformers of the 90-atom C/H/N molecule below, base GFN-FF, reference GFN2, 282 training /
71 held-out structures (`-test_fraction 0.2`), several split seeds (`-activelearn.seed`), fitted
from the stored dataset. Relative energy MAE in kcal/mol, gradient RMSE in Eh/A, mean +- standard
deviation over the seeds. The seeds share the same 353 conformers, so the intervals are rather
too narrow than too wide.

| fit | rel. energy MAE | paired vs base | gradient RMSE | n |
|---|---|---|---|---|
| base GFN-FF | 0.894 +- 0.050 | -- | 18.0e-3 | 5 |
| joint, old grid (M 1000, lambda 1e-6, 8 shells, eta 8) | 1.154 +- 0.054 | **+0.260 +- 0.059** | 4.50e-3 | 5 |
| joint, 6 shells, eta 4 (new default) | paired -0.215 +- 0.100 vs old grid, 10/10 better | +0.110 +- 0.095 | about as old grid | 10 |
| joint, M 4000, lambda 1e-7, 6 shells, eta 4, rcut 7, zeta 1,2,4,8,16 | 0.833 (old grid 1.139) | +0.019 +- 0.076 | 3.05e-3 (4.41e-3) | 10 |
| energy only, lambda 1e-6 | 0.698 +- 0.094 | -0.196 +- 0.119 | 18e-3 | 5 |
| energy only, lambda 1e-4 (energy-model default) | 0.623 +- 0.028 | -0.191 +- 0.050, 10/10 better | ~ base | 10 |

Read plainly:

- **The joint model is not better than GFN-FF in relative energies on unseen conformers**, with
  the old grid it is worse. Its gain is in the forces: gradient RMSE about 4x lower.
- **Energy and force trade against each other** (`lambda`, `sigma_energy`/`sigma_force`,
  `use_forces`). The energy-only fit has the smallest relative energy error and forces no better
  than the base method.
- The best joint combination measured (row 4) brings the energy to the base level and the
  gradient to 3.05e-3, but the fit takes 95 s instead of 10.7 s, and it is not better than the base
  in energy. Each ACSF parameter alone moved the energy MAE by at most 0.16 kcal/mol.
- **The energy model lowers the relative energy error but does not improve the ranking**:
  Spearman 0.467 against about 0.50 for the base (lambda 1e-4, paired -0.031 +- 0.099, 5 of 10
  seeds better, n=10). For the energy model, the centre count and the ACSF parameters were
  within the seed noise; lambda <= 1e-8 or a kernel width <= 0.75x the median made it clearly
  worse; 3x the median width alone was as good (0.622), combined with lambda 1e-4 and an ACSF
  variant worse (0.746 +- 0.061, n=10).
- One molecule only. The earlier single-split numbers below (M = 2000) and an older note of
  "relative energy MAE 0.647 -> 0.167, Spearman 0.41 -> 0.89" were not reproduced held-out.

## Two models: which one for what

`-activelearn` in `train` mode writes both: `<basename>.deltaml.json` (`"role": "joint"`,
energies + gradients) and `<basename>.deltaml.energy.json` (`"role": "energy"`). The log prints
the held-out block of each, the joint model's last. Model files without `role` read as joint.

- **Optimisation and dynamics** (`-opt`, `-md`): the joint model. Using the energy model there
  works but prints a one-time warning: its forces are no more accurate than the base method's.
- **Relative energies of given structures** (`-sp`, also on a multi-XYZ file, which `-sp`
  evaluates frame by frame): the energy model. Smaller relative energy error, rank correlation not
  better than GFN-FF (n=10 seeds).
- **ConfSearch**: the documented switches are `-method`, `-md_method` and `-opt_method`
  ([CONFSEARCH_DUAL_METHOD.md](CONFSEARCH_DUAL_METHOD.md)); `-opt_method` does both the per-cycle
  re-optimisation and the final ranking, there is no separate ranking-only method. So
  `-opt_method deltaml:<main>.deltaml.json` optimises and ranks on the joint model; to rank with
  the energy model, re-evaluate the final ensemble afterwards with
  `curcuma -sp <final ensemble>.xyz -method deltaml:<basename>.deltaml.energy.json`. Not tested
  end to end.

## Measured on a real system

Single split, M = 2000, old ACSF grid (8 shells, eta 8). The energy gain in this table is not
confirmed by the multi-seed held-out measurement above.

353 conformers of a 90-atom C/H/N molecule (`test_cases/confs/conf.xyz` -- these are GFN2
minima, so "distance from the reference minimum" is directly measurable), base GFN-FF,
reference GFN2, 282 training / 71 held-out, M = 2000. Held-out unless noted:

| lambda | energy MAE | gradient RMSE | Spearman | mean dist. to GFN2 min | NVE drift |
|---|---|---|---|---|---|
| base GFN-FF | 0.824 | 9.53e-3 | 0.387 | 0.635 A | 1.56 |
| 1e-8 | 0.738 | **1.95e-3** | 0.490 | 0.587 A | 1.59 |
| **1e-6** | 0.714 | 2.06e-3 | **0.517** | 0.567 A | **1.11** |
| 1e-4 | **0.709** | 3.83e-3 | 0.350 | -- | -- |
| 1e-2 | 0.790 | 8.09e-3 | 0.303 | **0.474 A** | -- |

Energies in kcal/mol, gradients in Eh/A, drift in kcal/mol/ps (GFN2 itself: 0.10). The
geometry column is the mean over ten conformers optimised on each surface; a single
structure is far too noisy for this -- individual values scatter between 0.34 and 0.68 A
with no clear trend.

Read honestly: **the gradients improve by a factor of four to five, the geometries by 7 to
25 %, and the relative conformer energies by about 13 %.** The ranking improves only at the
weaker regularisations. There is no single best lambda -- 1e-6 wins on ranking and energy
conservation, 1e-2 on geometry -- so the default is 1e-6 and the choice belongs to the
application.

Steady-state cost per force evaluation in molecular dynamics (same molecule, 500 steps,
the number that matters because the model is loaded once and evaluated thousands of times):

| Method | ms/step | relative to GFN2 |
|---|---|---|
| GFN-FF | 0.91 | 88x faster |
| GFN-FF + Delta-ML (M = 2000) | 4.83 | **16.6x faster** |
| GFN2 | 80.1 | -- |

The correction costs about four times a GFN-FF evaluation while remaining an order of
magnitude cheaper than the reference. This is why the model is sparse: a dense kernel over
all 25380 training environments would have been slower than the reference it replaces. How
far M can be pushed is set by the reference/base cost ratio, 88x here and four to five
orders of magnitude for DFT.

## With a DFT reference the case is different

Everything above uses GFN2 as the reference, where the cost ratio to GFN-FF is only 88x.
That is the wrong regime to judge the method in. Repeated with ORCA r2SCAN-3c on caffeine
(24 atoms, 96 frames from a 300 K GFN-FF trajectory, M = 400, lambda = 1e-6):

| | GFN-FF | Delta-ML | factor |
|---|---|---|---|
| relative energy MAE | 2.998 | 1.300 kcal/mol | **2.3x** |
| gradient RMSE | 2.17e-2 | 1.05e-2 Eh/A | 2.1x |
| Spearman | 0.825 | 0.914 | -- |
| NVE drift | 0.41 | 0.49 kcal/mol/ps | comparable |

The energies improve 2.3x here against 1.16x with GFN2 as reference. That is the more
interesting number: GFN-FF and GFN2 are built on related approximations, so their difference
is small and largely unsystematic, while the difference to a density functional carries more
structure and is correspondingly more learnable.

The economics change even more. r2SCAN-3c costs 31.4 s per gradient on this molecule against
11.8 ms for GFN-FF -- a ratio of 2660x, thirty times the GFN2 case -- and evaluating the
correction costs 0.33 ms per step (0.48 against 0.15 for the bare base method), i.e. about
one hundredth of one percent of a reference calculation.

Building the model cost one reference calculation per training structure, so it breaks even
at exactly the training-set size and everything past that is close to free:

| structures evaluated | pure DFT | Delta-ML | |
|---|---|---|---|
| 200 | 1.74 h | 0.84 h | 2.1x |
| 1000 | 8.72 h | 0.84 h | 10.4x |
| 5000 | 43.6 h | 0.85 h | 51x |

**Caveats, and they are not small.** One rigid molecule, 96 structures, a held-out split of
19, and M = 400 -- all of these numbers carry real uncertainty. A residual of 1.3 kcal/mol is
still above thermal energy at room temperature, so this does not yet rank conformers
reliably; it makes GFN-FF roughly twice as good for a tenth of the cost, which is a different
and more modest claim.

## Two use cases, and they do not behave alike

Measured on a second system -- a 66-atom trisaccharide, GFN-FF to GFN2, M = 2000,
lambda = 1e-6 -- with three training sets differing only in where the structures came from:

| training set | `Delta` span | energy MAE (model / base) | gradient RMSE (model / base) | mean dist. to GFN2 min |
|---|---|---|---|---|
| 302 raw 400 K dynamics frames | 282 kcal/mol | 15.86 / 27.42 | 2.56e-2 / 4.62e-2 | 0.894 A |
| 150 GFN-FF-optimised structures | 126 kcal/mol | 10.37 / 12.01 | 1.08e-2 / 2.58e-2 | 0.793 A |
| 60 GFN2-optimised structures | 430 kcal/mol | 12.08 / 64.76 | 1.75e-2 / 5.80e-2 | 1.012 A |
| bare GFN-FF | -- | -- | -- | **0.546 A** |

Every model improves the energies and the forces. **Every model also moves the optimised
minimum further from the reference than the uncorrected base method does** -- and the
training set with the largest energy gain (5.4x) produces the worst geometry.

That is not a contradiction. A lower RMS force error says nothing about *where* the gradient
vanishes: a correction can be more accurate on average and still displace the stationary
point, and a floppy sugar with shallow minima amplifies exactly that. The 90-atom conformer
molecule above, which is stiffer, improves by 7 to 25 % on the same metric.

So the two applications need separate treatment and separate parameters:

- **Conformer work** (training data near minima, used for ranking and pre-optimisation):
  measure the rank correlation and the distance to the reference minimum. `lambda = 1e-6`
  for ranking, `1e-2` for geometry, on the system measured here.
- **Trajectory work** (training data from dynamics, used for forces): measure the gradient
  RMSE and the NVE drift. Geometry is not the relevant criterion, and on a floppy system it
  will look bad no matter what.

**Do not assume the geometric gain transfers between systems.** It did not here, and the
direction of the effect reversed. Measure it for your system before relying on it.

## Does it survive molecular dynamics?

Energy conservation in NVE is the test that exposes a rippled learned surface: the defect is
invisible in the energy error and accumulates linearly along a trajectory. **Start each run
from a structure optimised on that same surface**, otherwise the comparison measures how far
each surface's minimum is from the shared starting point rather than its conservation.

90-atom conformer molecule, 2 ps, 0.5 fs step, thermostat off, each run started from its own
optimised minimum:

| surface | drift [kcal/mol/ps] | RMS fluctuation | mean T |
|---|---|---|---|
| GFN2 | 0.10 | 0.19 | 154 K |
| GFN-FF | 1.56 | 4.91 | 153 K |
| Delta-ML, lambda = 1e-8 | 5.84 | 5.79 | 157 K |
| Delta-ML, lambda = 1e-4 | 2.71 | 4.77 | 156 K |
| Delta-ML, lambda = 1e-2 | 1.58 | 4.74 | 153 K |

The learned surface conserves energy as well as the base method at every regularisation
tested, and better than it at lambda = 1e-6 (1.11 against 1.56). An earlier version of this
document reported a steep trade-off here, with only strong regularisation giving stability;
that was a unit defect in the gradient, not physics, and it is gone.

Skipping the pre-optimisation is still not a detail: starting from an unoptimised geometry,
the test measures how far each surface's minimum sits from the shared starting point rather
than its conservation.

## Where the improvement shows up, and where it does not

Same 90-atom conformer model, distance of the optimised minimum from the GFN2 minimum
(Kabsch RMSD over 90 atoms):

| surface | RMSD to the reference minimum |
|---|---|
| GFN-FF | 0.377 A |
| Delta-ML, lambda = 1e-4 | **0.119 A** |

Pre-optimisation on the corrected surface lands 3.2x closer to the reference minimum, which
is the point of the exercise. But this is system dependent: on the trisaccharide the same
comparison gave 0.288 A against 0.273 A, a 5 % gain, while its *energies* improved 4.8x. A
sugar's geometry is dominated by soft dihedrals and hydroxyl rotations that a small energy
correction barely moves. Do not assume the geometric gain transfers -- measure it.

## GPU

`-gpu cuda` offloads the fit. The training cost is dominated by one symmetric rank-k
update per structure -- `M^2 * 3N`, which is 2.2 GFlop at M = 2000 and N = 90 against
roughly 0.16 GFlop for building the design block itself -- so that single cuBLAS DSYRK is
what gets moved, and everything else stays on the CPU.

353 conformers, 90 atoms, RTX 5080:

| M | CPU | CUDA | speedup | held-out MAE (both) |
|---|---|---|---|---|
| 2000 | 26.2 s | 16.6 s | 1.6x | 0.2992 kcal/mol |
| 4000 | 78.6 s | 33.9 s | 2.3x | 0.2903 kcal/mol |
| 8000 | 246.9 s | 86.6 s | 2.9x | 0.2775 kcal/mol |

The speedup grows with M because the rank-k update's share of the work grows; the ceiling
is set by the descriptor, Jacobian and design-block construction that remain on the host.
Accuracy is identical to every printed digit -- cuBLAS sums in a different order, so the
systems agree to rounding rather than bit-for-bit, and that difference does not survive
into the reported errors.

The row-major design block is handed to cuBLAS unchanged: read with leading dimension M it
is already the column-major transpose the call wants, so the two storage conventions cancel
and nothing is repacked.

Prediction moves too: centres, weights and scaling are uploaded once, and each step
exchanges only the standardised descriptor in and the feature-space adjoint out -- a few
tens of kilobytes, against the ~5 MB the compact Jacobian blocks would cost. Those stay on
the host and are contracted there. 90 atoms, M = 2000, per dynamics step:

| | CPU | CUDA |
|---|---|---|
| GFN-FF alone | 0.31 ms | 0.68 ms |
| Delta-ML total | 5.74 ms | 3.16 ms |
| correction alone | 5.43 ms | **2.48 ms** |

The correction is 2.2x faster; the correction value itself is identical to all ten printed
digits. Note GFN-FF is *slower* on the GPU at this size, which eats part of the gain --
`-gpu` is a single global switch, so the best placement per component is not expressible.

## ConfSearch: not established, and the roles do not line up

Untestable on this branch -- ConfSearch currently segfaults here with `gfnff` and `uff` as
well, so the failure is pre-existing and says nothing about Delta-ML. What can be said from
the measurements is that the three roles ConfSearch needs a method for do not play to this
method's strengths:

| role | what it needs | Delta-ML today |
|---|---|---|
| metadynamics sampling | forces | strong: 4.6x better gradients, NVE stable |
| optimising the bias structures | minimum positions | weak and system dependent, went the wrong way on one of two systems |
| ranking in ConfScan | relative energies | joint model not better than GFN-FF held-out; energy model smaller error, Spearman not better (see "Held-out accuracy") |

There is also a structural tension: metadynamics exists to drive the system into regions it
has not visited, which is exactly where a kernel model knows nothing and its correction
decays to zero. `uncertainty_abort` is the right countermeasure but would end the biased
runs continuously.

And the economics do not favour it either. ConfSearch already samples with GFN-FF, which is
cheap; Delta-ML costs four to five times more per step than bare GFN-FF. It would only pay
by replacing the reference in the optimisation and ranking steps -- the two roles above where
it is weakest.

## Open: the anchor freezes more than it should

The topology anchor exists to guarantee that GFN-FF uses one identical bond graph
everywhere. Measurement says it currently freezes more than that, and only on the CPU.

Same structure, same model, base method only:

| | anchor == input | anchor != input |
|---|---|---|
| pure GFN-FF, CPU | -16.28842478 | -16.29445505 |
| pure GFN-FF, CUDA | -16.28842476 | -16.29445503 |
| Delta-ML base, CPU | -16.28842478 | **-16.29389806** |
| Delta-ML base, CUDA | -16.28842476 | -16.29444061 |

With the anchor equal to the input all four agree to 2e-8. As soon as the structure moves
away from the anchor, the CPU path departs by 5.6e-4 Eh (0.35 kcal/mol) while the GPU stays
within 1.4e-5 of the geometry-consistent value. The GPU refreshes the coordination numbers
and the lists derived from them; the CPU keeps them pinned at the anchor.

The GPU is the correct behaviour here. A coordination number is a smooth function of the
coordinates and belongs to the dynamic state by design; pinning it makes the energy wrong
for the actual geometry and, worse, inconsistent with the gradient -- which is exactly the
signature of the NVE drift reported above. Fixing the anchor to freeze the bond graph only
is the obvious next step, and the drift numbers should be remeasured afterwards.

**Measured over a whole data set (Oct 2026):** the `energy_base` stored by `-activelearn` for the 353 conformers is
not what `-sp -method gfnff` returns for the same geometry: fresh minus stored = -1.21 +- 0.70 kcal/mol (min -3.05,
max 0.00), as large as the GFN2 spread of the set (2.2 kcal/mol). The model learns `E_ref - E_stored` and is applied
on top of the fresh energy by `-sp`; the true GFN2 minimum has GFN-FF rank 7 with the stored and 31 with the fresh
energy. The labelling and the single-point path must use the same base energy before the held-out numbers can be
transferred to deployment (details and the effect on the ranking test:
[DELTA_ML_USECASE_TEST_2026-10.md](DELTA_ML_USECASE_TEST_2026-10.md)).

**Cause and fix (Oct 2026):** not the anchor. `-activelearn` reused one `EnergyCalculator` per
thread and called `setMolecule()` for each structure. A re-initialised `GFNFF` keeps
`EEQSolver::m_external_srab` pointing at the distance array of the *previous* structure, so the
Phase-2 EEQ solve during topology setup uses stale distances, and `Calculation()` then skips its
own Phase-2 solve because the geometry has not changed since setup. Only the Coulomb term moved
(charges off by up to 0.07 e); clearing that pointer in a test build, or a new calculator per
structure, removed the deviation (40 trisaccharide frames: -4.26 +- 1.87 kcal/mol to 1e-14 Eh),
while dropping the bond pinning or `cache_topology false` changed nothing. Labelling now builds
fresh base and reference calculators per structure (as `-sp` and `DeltaMLMethod` do), keeps the
anchor bonds (which `DeltaMLMethod` applies identically: deployed base energy = stored label to
the 1e-10 print precision, 7 structures), forces one thread per calculator also in the `global`
scope (`-threads 4` used to give every GFN-FF instance 4 threads) and writes the dataset in input
order (formerly thread-chunk order, so a given `-activelearn.seed` now selects a different
train/test split than before). Result: all 401 frames of both use-case sets (66 and 90 atoms, GFN2
reference) equal `-sp -dump_gradient` exactly (max deviation 0 for all four labels) and `-threads 1`
and `-threads 4` write byte-identical datasets. Regression test: `cli_deltaml_04_label_consistency`.
A general fix belongs in `GFNFF::InitialiseMolecule` (reset the external distance pointer); any
other caller that re-initialises one GFNFF object for a new geometry is still affected.

## A numerical trap worth knowing about

The fit accumulates normal equations `(Phi^T W Phi + lambda K_MM) beta = Phi^T W y` rather
than forming the full design matrix, because `Phi` has one row per energy plus 3N rows per
structure and would not fit in memory. The price is that normal equations **square the
condition number**, and `K_MM` is itself near-singular once many centres crowd together in
descriptor space.

Measured on the 353-conformer set at M = 4000, with and without the conditioning ridge:

| `jitter_relative` | training MAE | held-out MAE |
|---|---|---|
| 0 | 0.145 kcal/mol | **0.554 kcal/mol** |
| 1e-8 | 0.099 kcal/mol | **0.167 kcal/mol** |

Note the shape of the failure: the *training* error stayed fine while the held-out error
blew up. The factorisation reported success and returned a wildly oscillating `beta` that
still fit the training data. The visible symptom is accuracy getting **worse** as
`n_centers` grows, which is otherwise impossible for a regularised fit -- if you ever see
that, suspect conditioning before suspecting the model.

## Verification status

Machine-tested only. What was checked:

- Descriptor Jacobian against finite differences: clean h^2 convergence down to the roundoff
  floor; translation invariance 2e-16, rotation invariance 1e-15 (both exact analytic
  identities, therefore sharper evidence than the FD comparison itself)
- Model gradient against finite differences: 1.2e-11 Eh/Angstrom
- Energy from the gradient path and from the plain path agree bit-for-bit
- JSON round-trip reproduces energies and gradients exactly
- Learning curve descends monotonically on a synthetic target with a known analytic form
- End-to-end fit on a real 90-atom conformer set (C/H/N), see the table above

What was **not** checked: other molecule classes, charged or open-shell systems, transition
metals, transferability across molecules, long-time MD stability on the corrected surface,
and any DFT reference. Treat results as provisional until tested on real problems.

## Limitations

- Fixed topology only, one molecule per model
- Training cost grows as M^2; centres are capped at 12000 with an explicit error
- Model files store all centre descriptors as JSON and grow with M (a few MB at M = 1000)
- Active learning runs but its benefit is unproven; the offline path is the validated one
- Geometric improvement is system dependent (3.2x on one system, 1.05x on another)
- Joint model: relative energies not better than the base held-out (one molecule, n=5-10 seeds); energy model: smaller error, ranking not better
- MD on the corrected surface needs `lambda` around 1e-2, which costs most of the accuracy
- No GPU path for the kernel or descriptor evaluation yet
