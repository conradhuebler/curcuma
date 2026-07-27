# Delta-ML: cheap method, reference-level ranking

**Status**: 🤖 AI-generated, ⚙️ machine-tested. Human production testing pending.


Learns the difference between a cheap base method and an expensive reference:

```
E_ref(R) ~ E_base(R) + Delta(R)
```

The gain is statistical. `Var(Delta)` is far smaller than `Var(E_ref)` and `Delta` is a
smoother function of geometry, so a few hundred reference points suffice where direct
learning would need tens of thousands. The purpose here is conformer work: run conformer
search and pre-optimisation at GFN-FF cost while ranking conformers at the reference level.

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

| Parameter | Default | Effect |
|---|---|---|
| `-ref_method` | `gfn2` | the level being learned; must match the base method's dispersion treatment |
| `-n_centers` | 1000 | kernel centres M; the accuracy/cost lever |
| `-kernel_sigma` | 0 (auto) | Gaussian width; strong lever, too narrow cannot interpolate, too wide oversmooths |
| `-activelearn.lambda` | 1e-4 | regulariser relative to the data term; raise to 1e-2 if NVE drifts. The flat `-lambda` does not route -- use the dotted form |
| `-jitter_relative` | 1e-8 | conditioning ridge on the normal equations; see below |
| `-use_forces` | true | train on gradients as well as energies |
| `-test_fraction` | 0.2 | held-out share, split randomly because neighbouring frames are correlated |
| `-center_sweep` | -- | comma-separated M values to fit and compare in one run |

## Measured on a real system

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
- MD on the corrected surface needs `lambda` around 1e-2, which costs most of the accuracy
- No GPU path for the kernel or descriptor evaluation yet
