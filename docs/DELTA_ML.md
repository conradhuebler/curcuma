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

353 conformers of a 90-atom C/H/N molecule (`test_cases/confs/conf.xyz`), base GFN-FF,
reference GFN2, 282 training / 71 held-out structures, split randomly. `Delta` spans
4.9 kcal/mol across the set. Held-out numbers:

| M | rel. energy MAE | gradient RMSE | Spearman |
|---|---|---|---|
| base GFN-FF | 0.647 kcal/mol | 5.04e-3 Eh/A | 0.408 |
| 200 | 0.460 | 2.10e-3 | 0.521 |
| 500 | 0.365 | 1.54e-3 | 0.750 |
| 1000 | 0.243 | 1.19e-3 | 0.849 |
| 2000 | 0.194 | 1.02e-3 | 0.878 |
| 4000 | 0.167 | 0.90e-3 | 0.893 |

Monotone and saturating, as a regularised fit should be. At M = 4000 the relative conformer
energies are 3.9x better and the gradients 5.6x better than bare GFN-FF, with the rank
correlation rising from 0.41 to 0.89.

Steady-state cost per force evaluation in molecular dynamics (same molecule, 500 steps, the
number that matters because the model is loaded once and evaluated thousands of times):

| Method | ms/step | relative to GFN2 |
|---|---|---|
| GFN-FF | 0.91 | 88x faster |
| GFN-FF + Delta-ML (M = 2000) | 4.83 | **16.6x faster** |
| GFN2 | 80.1 | -- |

The correction costs about four times a GFN-FF evaluation, so it dominates the base method
while remaining an order of magnitude cheaper than the reference. This is precisely why the
model is sparse: a dense kernel over all 25380 training environments would have been slower
than the reference it replaces. How far M can be pushed is set by the reference/base cost
ratio, which is 88x here and four to five orders of magnitude for DFT.

Single-point runs additionally pay for parsing the model file (2.8 MB at M = 2000), which is
why `-sp` timings look worse than the per-step cost above.

## Train on optimised structures, not raw dynamics frames

This is the single biggest lever on how much the model is worth, larger than any
hyperparameter. Measured on a 66-atom trisaccharide, GFN-FF to GFN2, same molecule, same
M = 2000, only the training set differs:

| training set | `Delta` span | energy MAE (model / base) | gain |
|---|---|---|---|
| 302 raw frames from 30 ps at 400 K | 371 kcal/mol | 27.5 / 12.4 kcal/mol | **1.3x** |
| 150 GFN-FF-optimised structures | 91 kcal/mol | 2.6 / 12.4 kcal/mol | **4.8x** |

Raw high-temperature frames teach the model about a region of configuration space that is
both far rougher than, and irrelevant to, the minima you actually optimise into. Worse, a
model trained that way made the optimised geometry *worse* than bare GFN-FF (0.405 A from
the reference minimum against 0.288 A). Optimise the sampled structures before labelling
them -- `curcuma -opt frames.xyz -method gfnff` -- and train on the result.

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

At `lambda = 1e-2` the learned surface conserves energy as well as the base method it
corrects. On the trisaccharide with an optimised training set it does better than the base
(0.96 against 1.22). Regularisation is what buys this, and it costs accuracy -- pick the
point on the curve your application needs, and check it rather than assuming it.

Skipping the pre-optimisation is not a detail: the same trisaccharide model started from an
unoptimised geometry drifted 29.9 kcal/mol/ps and heated from 300 K to 544 K.

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
