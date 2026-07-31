# NEB-MD: NEB-like Molecular Dynamics

> Status: AI-generated, machine-tested (CLI test `cli_nebmd_01_gfnff_ethane` passes).
> NOT ✅ TESTED or ✅ APPROVED. Validate against a known reaction coordinate before research use.

NEB-MD (`curcuma -nebmd`) drives a chain of `N` replicas (images) along a reaction
path between a start and an end structure with **molecular dynamics** instead of the
classic LBFGS band optimisation. Each interior image runs its own MD (independent
velocities, thermostat, kinetic energy); the images are coupled by spring forces to
their neighbours so they move constrained along the path instead of freely relaxing
to the nearest minimum.

## Read this first

**What works.** The band machinery, the CI-NEB optimiser (FIRE + climbing image), the
path collective variables `s`/`z` with FD-verified analytic gradients, reading an
external MEP, and the WHAM plumbing. On a clean test case (ethane rotation) the CI-NEB
reproduces the experimental barrier: GFN2 2.88 vs ~2.9 kcal/mol.

**What does not work yet.** Free-energy barriers. On helicene, against references of
37.8 (r2SCAN-3c OptTS), 30.4 (GFN2 OptTS) and 40.8 (ORCA GFN2 NEB) kcal/mol, this code
gives 58.4 (mean force) and 13-19 (umbrella+WHAM). **Do not quote free energies from
`-nebmd.pmf` or `-nebmd.umbrella`.**

**The two things that bite, in order.**
1. *The generated path can be the wrong mechanism.* The internal NEB drove helicene
   through a planar structure with **zero imaginary modes** - not a saddle at all. Always
   verify a NEB saddle with a frequency calculation (exactly one imaginary mode), and
   prefer `-nebmd.path_file` with an externally converged MEP.
2. *The free-energy estimator is wrong - not the sampling.* A 100 ps band run shows the
   image drift is stationary (1.11 A at 400 fs, 1.13 A at 100 ps) and its mean (0.75 A)
   equals the free molecule's own thermal motion (0.755 A): the springs hold a proper
   equilibrium. But the mean-force barrier *falls* with better statistics (58.4 -> 4.69
   kcal/mol), because subtracting the spring bias makes the two terms cancel. The naive
   bias subtraction needs replacing by a constrained-ensemble estimator including the
   metric term.

Sections below are ordered method-first, then validation, then the record of what was
measured and which earlier conclusions it overturned.

## Concept

- `N` images are placed along an initial path: the end structure is RMSD-reordered
  (atom order) and Kabsch-aligned onto the start frame, then linearly interpolated.
- Each **interior** image is a `SimpleMD` instance driven through its stepwise API
  (`prepareRun` / `step` / `applyExternalForces`).
- The coordinator (`NebMD::stepBand`) reads every image's current positions + true
  gradient, computes the tangent and spring force, injects an external force into
  each image, and advances all interior images by one MD step (lockstep).
- Endpoints are clamped by default (`-nebmd.fixed_endpoints true`): not integrated,
  their positions stay fixed; their single-point energies feed the improved tangent.

## Force model

For interior image *i*, atom *a* (positions R in Angstrom, true gradient
`g = dE/dR` in Eh/Angstrom from `image.gradient()`, unit tangent `tau`):

```
F_spring = k_a * (R[i+1] + R[i-1] - 2 R[i])          # full harmonic neighbour spring
```

**Sign convention (easy to get wrong)**: `SimpleMD::applyExternalForces()` adds its
argument to `m_eigen_gradient` (i.e. to `dE/dR`), and the integrator then applies
`force = -gradient`. So to apply a desired extra FORCE you must pass its
**negative**. `SimpleMD` already applies `F_true = -g`.

- **Nudged** (`-nebmd.nudged true`, default): desired total force is
  `F_true_perp + F_spring_parallel = -g + (g.tau) tau + F_spring_par`, so the extra
  force on top of what SimpleMD already applies is `(g.tau) tau + F_spring_par` and
  the injected gradient contribution is
  ```
  F_ext = -[ (g . tau) tau + (F_spring . tau) tau ]
  ```
  (NEB nudging: true force perpendicular to the path, spring force parallel.)
  This prevents corner-cutting and the band sliding downhill.

- **Non-nudged** (`-nebmd.nudged false`): `F_ext = -F_spring` (plain elastic band;
  full true force + full spring). Simpler but the band can slide and corner-cut.

The tangent is the improved (energy-weighted) Henkelman-Jonsson 2000 tangent.

### Timing caveat (resolved)

`SimpleMD::Verlet()`/`Rattle()` re-evaluate the gradient at the *new* position
mid-step (`Energy()` overwrites `m_eigen_gradient`). The external force injected
via `applyExternalForces` is now re-applied after `WallPotential()` so it acts on
**both** velocity half-kicks (and the position update), making the scheme
time-symmetric for a constant external force over the step. Plain `-md` is
unchanged (the buffer is empty when no external force is queued).

## Spring constant and timestep

`SimpleMD` works in internal units of Angstrom / fs / amu / Eh, so `k_spring` is
in **Eh/Angstrom^2** (or **Eh/(amu*Angstrom^2)** when mass-weighted).

- `-nebmd.mass_weighted_spring true` (default): `k_a = k * m_a`, so the spring
  frequency `omega = sqrt(k)` is mass-independent. This keeps hydrogen atoms
  stable (otherwise H would need a tiny timestep).
- Default `k_spring = 0.005`, `time_step = 0.5 fs` (passed via `-simplemd.*`):
  `omega ~= 0.07 / fs`, spring period ~ 90 fs, ~180 steps per period - stable.
- Raise `k` for a stiffer band (closer to classic NEB spacing) at the cost of a
  smaller `dt`; lower `k` for looser coupling.

## Usage

```
curcuma -nebmd start.xyz end.xyz \
    -nebmd.nimages 12 \
    -simplemd.method gfnff -simplemd.time_step 0.5 -simplemd.max_time 200 \
    -simplemd.temperature 300 -simplemd.thermostat csvr \
    -nebmd.k_spring 0.005 -nebmd.nudged true -nebmd.fixed_endpoints true \
    -nebmd.dump_frequency 20
```

Recommended workflow when a barrier is the goal (rather than just a relaxed band):

```
# 1. get a real path - either an external MEP ...
curcuma -nebmd start.xyz end.xyz -nebmd.path_file neb_MEP_trj.xyz \
    -md.method gfn2 -md.max_time 200
# ... or optimise one here (then VERIFY the saddle with a frequency calculation)
curcuma -nebmd start.xyz end.xyz -nebmd.optimize true -nebmd.idpp false \
    -nebmd.opt_fmax 0.02 -md.method gfn2 -md.max_time 1
```

Note `-nebmd.idpp false` with `-nebmd.optimize true`: IDPP improves the starting
energies but was measured to seed the band in the wrong valley on helicene.

### Parameter reference (added during the free-energy work)

| flag | meaning |
|------|---------|
| `-nebmd.path_file` | read the reference path from a multi-XYZ (e.g. ORCA `neb_MEP_trj.xyz`); frame count sets `nimages` |
| `-nebmd.optimize` | run a classical CI-NEB (FIRE) before the MD |
| `-nebmd.opt_fmax` | force convergence threshold, Eh/A (0.02 is realistic; 0.01 stalls on strained systems) |
| `-nebmd.climbing` | climbing image, latched after pre-relaxation |
| `-nebmd.umbrella` | umbrella sampling on the path CV `s` instead of neighbour springs |
| `-nebmd.umbrella_k` | restraint on `s`; **0 = auto** (targets sigma ~ half the window spacing) |
| `-nebmd.umbrella_kz` | tube restraint on the perpendicular CV `z`; 0 = off |
| `-nebmd.umbrella_z_tolerance` | free tube radius before `z` bites (see the free-MD scale: ~1 A, not 0.05) |
| `-nebmd.umbrella_keep_springs` | keep neighbour coupling in umbrella mode (diagnostic - the bias is NOT removed by WHAM) |
| `-nebmd.reparametrize` | redistribute images to equal arclength every N steps |
| `-nebmd.pmf` | mean-force estimator (springs mode only; disabled automatically under umbrella) |

### Parameter routing

- **Band parameters**: `-nebmd.*` (`nimages`, `nudged`, `fixed_endpoints`,
  `k_spring`, `mass_weighted_spring`, `dump_frequency`, `write_initial_path`,
  `proton_transfer`, `force_reorder`, `threads`).
- **The energy method is shared** by the band optimisation, the endpoint single points
  and the MD images. It is taken from `-simplemd.method` / `-md.method` if present,
  otherwise from `-nebmd.method` / the global `-method` (default `gfnff`).
  (Bug fixed Jul 2026: `m_method` was read from the `nebmd` scope only, so
  `-md.method gfn2` silently optimised the band with the `gfnff` default - GFN2 and
  GFN-FF returned a bit-identical barrier, which is how it was caught.)
- **MD parameters for the images** can be passed as `-simplemd.*` (preferred),
  `-md.*`, or bare single-owner flags (`-time_step`, `-max_time`, `-thermostat`,
  `-method`). `executeNebMD` forwards `controller["md"]` and ambiguous bare
  flags (e.g. `-temperature`, multi-owner casino+simplemd) into
  `controller["simplemd"]` so the images receive them; `-simplemd.*` wins on
  conflict. The energy method default is `gfnff`.
- `threads` defaults to 1 per image: the gfnff `CitationRegistry` has a known
  race at `threads > 1` (see `src/capabilities/CLAUDE.md`). Keep 1 for gfnff.

### Endpoints

- Same atom order in both files: only Kabsch alignment is applied (atom order
  preserved).
- Different atom order: a Munkress reorder (`-nebmd.force_reorder true`) is
  attempted; if the result still does not match the start atom order, NEB-MD
  errors out (pre-align with `curcuma -nebprep start.xyz end.xyz`).

## Outputs (BMT directory, copied back to CWD)

- `<base>.neb.init.xyz`  - initial interpolated path (multi-XYZ, N frames).
- `<base>.neb.path.xyz`  - band snapshot every `dump_frequency` steps (multi-XYZ,
  N frames per snapshot, comment `path step=.. image=.. E=.. T=..`).
- `<base>.neb.energy.csv` - `step,image,Epot,Ekin,T` per row.
- `<base>.neb.final.xyz` - final band (N frames).

Use BMT (the default) so each run is isolated. With `-no_bmt`, GFN-FF parameter
caches (`*.param.json`) and topology caches (`*.topo.json`) persist in the CWD.
Both are now geometry-validated on load (`param.json` via a geometry signature in
`ForceField::tryLoadAutoParameters`; `topo.json` via the existing bond-list
fingerprint) and regenerated when the structure changes, so stale caches no
longer yield wrong energies. NEB-MD images also disable per-image restart
loading (`no_restart=true`, `setRestart(false)`) so stale snapshot files cannot
poison a fresh run.

## Topology of interpolated images (GFN-FF caveat)

GFN-FF derives its bond topology (bonds, rings, hybridisation) from the geometry.
Linear interpolation between two endpoints can stretch a bond past the detection
threshold at an intermediate image, or bring two atoms into bonding range that
are not bonded at the endpoints. That image then gets a **different topology**,
so the GFN-FF energy/gradient is discontinuous along the band and the
"minimum-energy path" is meaningless.

> **Caveat (found later, see the CI-NEB validation below): IDPP improves the initial
> energies but can seed the band in the WRONG valley.** On [6]helicene the optimised
> IDPP band ends up with a double maximum at 275 kcal/mol, while plain linear
> interpolation optimises to a clean 56 kcal/mol single maximum. When you combine the
> path with `-nebmd.optimize true`, try `-nebmd.idpp false` first.

**IDPP path (default on, `-nebmd.idpp`)** — Instead of using
the raw cartesian interpolation, the initial path is refined with the
Image-Dependent Pair Potential (Smidstrup, Pedersen, Stokbro, Jonsson,
*J. Chem. Phys.* **140**, 214106 (2014)):

```
S_i(R) = sum_{a<b} ( d_target,i(a,b) - d(a,b) )^2 / d(a,b)^4
```

The target distance matrix `d_target,i` is the linear interpolation of the two
endpoint *distance matrices* (a far better interpolation of a molecular structure
than cartesian averaging), and the `1/d^4` weight makes short contacts dominate —
which is exactly the "take the repulsion into account" part: a pair closer than
its target is penalised enormously, so the optimiser pushes those atoms apart.

The refinement runs as a NEB on the IDPP surface (`-nebmd.idpp_iterations`,
default 200): each image feels the perpendicular IDPP force **plus a parallel
spring to its neighbours** with the bisecting tangent, so the images stay evenly
spaced and the path stays connected while the clashes are resolved — i.e. the
path accounts for the neighbouring images, not just for each image in isolation.

Measured on the [6]helicene A->B path (42 atoms, 8 images, GFN-FF), image
energies at step 0 in Eh:

| image | linear | linear + clash repair | **IDPP** |
|-------|--------|-----------------------|----------|
| 1     | -9.187 | -9.187                | **-9.216** |
| 2     | -8.462 | -8.456                | **-8.971** |
| 3     | -7.781 | -8.751                | **-8.499** |
| 4     | -7.765 | -8.744                | **-8.499** |

IDPP objective 2.77e+02 -> 2.31e+00. The IDPP path is smooth and symmetric
(images 3/4 degenerate, 1<->6 and 2<->5 paired), as an inversion path should be.
More importantly, after 200 MD steps the IDPP band stays **above** the endpoints
(-9.29 ... -8.72 Eh) whereas the plain interpolated band collapses **below** them
(-9.71 Eh at images 3/4) — i.e. without IDPP the "barrier" is a clash artefact.

**Overlap repair (default on, `-nebmd.repair_overlaps`)**: after IDPP (or instead
of it, when `-nebmd.idpp false`), and before the topology check,
NEB-MD relaxes the interpolated images with the same soft clash potential that
`PolymerBuild::resolveOverlaps` uses — steepest descent on `(r_cut/r)^4` over
non-bonded pairs (`r_cut = 1.15 * (r_cov,i + r_cov,j)`, capped 0.05 A/step,
bonded pairs taken from image 0 are exempt, endpoints untouched). This removes
the atom-atom interpenetration that linear interpolation produces for
non-trivial conformational changes. Measured on a [6]helicene A->B path (42
atoms, 8 images): 4 images repaired, the two worst images going from -7.78 /
-7.76 Eh to -8.75 / -8.74 Eh at step 0 (i.e. ~630 kcal/mol of pure clash
artefact removed).

NEB-MD then checks the topology (`NebMD::checkPathTopology`):
it compares the bond-connectivity signature (from `Molecule::DistanceMatrix`) of
every image against image 0 and emits a `[WARN]` if any image differs, naming
the image and the bond counts. The check is advisory (the run continues) because
a borderline bond at the threshold can trigger a false positive; treat a
warning as a signal that the path is not a clean conformational change.

For **bond-breaking/forming reactions** (where the topology genuinely must
change along the path) use a QM method (`-simplemd.method gfn2`) instead of
GFN-FF, or supply a better initial path (e.g. one prepared with
`curcuma -nebprep`). For conformational/rotational paths the topology is
invariant and no warning appears.

## Band optimisation - classical CI-NEB (`-nebmd.optimize`, default OFF)

Before (or instead of) the MD, the band can be relaxed onto the minimum energy path
with a classical nudged elastic band:

```
curcuma -nebmd start.xyz end.xyz -nebmd.nimages 12 -nebmd.optimize true \
    -nebmd.idpp false -nebmd.opt_fmax 0.01 -md.method gfnff
```

- Optimiser: **FIRE** (Bitzek et al., PRL 97, 170201 (2006)). L-BFGS is not usable
  here: the nudged force is not the gradient of any potential (the parallel true
  force is projected out), so a line search has no energy to minimise along.
- **Climbing image** (`-nebmd.climbing`, default on; Henkelman, Uberuaga & Jonsson,
  JCP 113, 9901 (2000)): the highest image inverts its parallel force and converges
  onto the saddle point instead of being held back by the springs.
- `-nebmd.opt_k` is the optimisation spring constant (independent of the MD
  `k_spring`), `-nebmd.opt_fmax` the force convergence threshold in Eh/Angstrom.
- Output: `<base>.neb.opt.xyz` (path with per-image `dE` in the comment) plus the
  forward/reverse barrier and reaction energy on the console.

### Validation

**Ethane rotation (staggered -> staggered, 120 deg), GFN-FF, 9 images.** Barrier must
sit in the middle of the band and forward must equal reverse:

| quantity | value |
|----------|-------|
| dE(forward), GFN-FF | 2.32-2.67 kcal/mol (image 4, the middle) |
| dE(reverse) | same as forward |
| reaction dE | -0.00 kcal/mol (identical conformers - correct) |
| **dE(forward), GFN2** | **2.88 kcal/mol** |
| rigid-scan reference (GFN-FF SP at ideal geometries) | 2.72 kcal/mol |
| experiment | ~2.9 kcal/mol |

GFN2 reproduces the experimental rotation barrier almost exactly (2.88 vs ~2.9) and is
only ~1.8x more expensive per single point than GFN-FF for this system (34 ms vs 19 ms,
42 atoms), so for barriers it is the better default whenever the system size allows.

The NEB value is *below* the rigid scan, which is the correct direction (the path is
allowed to relax at the saddle). FIRE demonstrably works: with `opt_fmax 0.0005` the
force drops 0.0123 -> 0.0018 Eh/A over 43 iterations and the barrier falls 2.42 -> 2.32.

**IMPORTANT - IDPP hurts here.** On [6]helicene the IDPP path optimises to a *double*
maximum (275 kcal/mol at images 4 and 7 with a 49 kcal/mol valley between) and a band
that folds back on itself, whereas plain linear interpolation optimises to a clean
symmetric single-maximum profile (56 kcal/mol). IDPP improves the *starting* energies
but can seed the band in the wrong valley. **Use `-nebmd.idpp false` together with
`-nebmd.optimize true`** unless you have checked that IDPP helps for your system.

### Root cause: the band converges to the WRONG MECHANISM (verified against a DFT OptTS)

A proper r2SCAN-3c transition-state optimisation (`! r2SCAN-3c OptTS Freq`, converged in
27 cycles) settles this definitively:

| | barrier |
|---|---------|
| **r2SCAN-3c, true optimised TS** | **37.8 kcal/mol** |
| literature [6]helicene | ~36 kcal/mol |
| r2SCAN-3c on our NEB TS geometry | 72.1 kcal/mol |
| GFN2 on our NEB TS geometry | 69.6 kcal/mol |

The DFT reference reproduces the experimental barrier, so **the NEB TS sits 34.3 kcal/mol
above the real saddle point** - on the same level of theory. This is not a PES question.

Comparing the two structures shows why, and it **corrects an earlier claim in this
document** that "the mechanism is right, only the geometry within it is wrong":

| | NEB TS | true OptTS |
|---|--------|------------|
| C-frame planarity | 0.0425 | **0.3594** |
| shortest H...H | 1.64, 1.94, 2.00 A | 1.83, 2.31, 2.39 A |
| RMSD between them | \multicolumn{2}{c}{1.35 A} | |

The real transition state is **not planar** (0.359, nearly as helical as the minimum at
0.399). The helicene *twists past* itself; it does not flatten. Our band flattens it
(0.043) and squeezes the terminal hydrogens together - the largest displacements between
the two structures are terminal H atoms moving 2.2-2.9 A.

So the defect is that linear interpolation between the two enantiomers runs straight
through the planar structure, and the optimiser cannot leave that basin: it converges to
a *different, higher* mechanism. That is a path-generation problem, not an optimiser
tolerance problem - which is why tightening `opt_fmax` or adding images never helped.

**Two independent confirmations** (both with the standalone `xtb` GFN2, `/opt/bin/xtb`):

1. *GFN2 is fine on the right geometry.* Single points on the **DFT-optimised** TS give
   **32.5 kcal/mol** (vs 37.8 DFT, ~36 experiment) - but **69.6 kcal/mol** on our NEB TS.
   Same method, same code: the geometry is what differs.
2. *Our "TS" is not a transition state at all.* `xtb --ohess` on the NEB TS reports
   **`# imaginary freq. = 0`**. A genuine first-order saddle has exactly one imaginary
   mode; ours has none, so it is simply a point on a hillside that the band could not
   descend.

Practical consequence: **verify any NEB saddle with a frequency calculation** before
trusting the barrier. One imaginary mode = transition state; zero = the band has not
found a saddle.

### Barrier vs image count, both methods (after all optimiser fixes)

`-nebmd.optimize true -nebmd.idpp false -nebmd.opt_fmax 0.02`, helicene:

| N | GFN-FF | status | GFN2 | status |
|---|--------|--------|------|--------|
| 8  | 56.69 | CONV (95 it)  | 58.65  | stalled |
| 12 | 59.62 | CONV (120 it) | 135.28 | stalled |
| 16 | 66.71 | stalled       | 148.00 | stalled |
| 20 | 81.51 | stalled       | 177.19 | stalled |

Forward and reverse agree to 0.01 kcal/mol in every run (the symmetric path is
reproduced correctly), but **the barrier still rises monotonically with N and most runs
stall** - and GFN2 diverges far worse than GFN-FF. None of these numbers is converged;
they are a diagnostic of the wrong-mechanism problem above, not method comparisons.

Reference values for the same system:

| | barrier |
|---|---------|
| r2SCAN-3c, verified TS (1 imaginary mode, -42.42 cm-1) | **37.8** |
| GFN2 on that verified TS geometry | 32.5 |
| literature | ~36 |
| our NEB TS (0 imaginary modes) | 69.6-72.1 |

### Optimiser bugs found and fixed along the way

r2SCAN-3c single points on **exactly the geometries the NEB produced** (ORCA, the same
level at which the endpoint structures were optimised) settle the question of whether
the semiempirical PES is to blame:

| method (identical geometries) | barrier |
|-------------------------------|---------|
| GFN-FF | 55.4 kcal/mol |
| GFN2 | 69.6 kcal/mol |
| **r2SCAN-3c (DFT)** | **72.1 kcal/mol** |
| literature [6]helicene | ~36 kcal/mol |

**DFT confirms GFN2 (72.1 vs 69.6).** So the semiempirical surfaces are fine - the
error is in the *geometry* the band converges to, and all three methods agree on it.

The defect is visible in the structure: the TS has an **H...H contact at 1.64 A**,
whereas the real planar helicene TS has ~1.9-2.0 A between the terminal hydrogens
(the vdW contact is ~2.4 A). The band squeezes the helix ends into each other instead
of letting them slide past one another. That single over-compressed contact accounts
for the ~36 kcal/mol excess, and it is genuine Pauli repulsion - which is why every
method reproduces it.

The C-frame planarity at the TS is 0.043, i.e. the mechanism (flattening the helix) is
right; it is the *terminal-ring geometry within* that mechanism that is wrong. This is
an optimiser/path problem: the band is not relaxing the perpendicular degrees of freedom
in the saddle region.

Consistent with that, the barrier depends on the *numerical* spring constant, which it
must not - in a converged CI-NEB the climbing image sits on the saddle with its spring
term switched off, so `opt_k` cannot influence the result.

| images | `opt_k=0.1` | `opt_k=0.5` |
|--------|-------------|-------------|
| 8  | **111.1 kcal/mol** | 48.4 kcal/mol |
| 20 | 65.3 kcal/mol | 70.5 kcal/mol |

All four runs **stalled** (none reached `opt_fmax`). A purely numerical parameter moving
the barrier by a factor of 2.3 means the reported values are unconverged snapshots and
the 48-111 kcal/mol spread is optimiser noise. The apparent "GFN-FF 56 vs GFN2 69.6"
difference and the monotone rise with image count are artefacts of the same thing - the
PES is not what varies when only `nimages` or `opt_k` change.

For reference, the evidence that the underlying methods are fine: GFN2 reproduces the
ethane rotation barrier (2.88 vs ~2.9 kcal/mol experimental), and independent spring-free
MDs of the two enantiomers give dG = -0.14 kcal/mol (correct, see the control experiment
above). The defect is in the band optimiser, and fixing it is the prerequisite for any
helicene number.

Known concrete defects: FIRE stalls at ~2x the force threshold with the timestep
collapsing (mitigated by a `dt` floor + stall detection, but the stall itself remains);
the step is unpreconditioned; `opt_k` is not scaled with the image spacing, so the band
gets effectively stiffer as `nimages` grows.

### Historical record: the image-count dependence that exposed this

[6]helicene (C26H16, verified enantiomers: signed volume +2.40 / -2.40, identical atom
order, 2.50 A RMSD), GFN-FF, linear interpolation, CI-NEB:

| images | result | barrier |
|--------|--------|---------|
| 8  | converged (328 it) | 44.3 kcal/mol |
| 12 | converged (344 it) | 56.0 kcal/mol |
| 16 | not converged (1000 it) | 60.8 kcal/mol |
| 20 | not converged (3000 it) | 65.2 kcal/mol |

The barrier **rises monotonically with the image count** and stops converging - the
opposite of the expected behaviour, where a denser band resolves the saddle better and
the barrier should settle. The literature value is ~36 kcal/mol (Martin & Marchant
1974). So the profiles are internally consistent (forward = reverse, reaction dE = 0
to 0.01 kcal/mol, symmetric) but the absolute barrier is **not converged and too high**.

**What "not converged" means concretely.** The criterion is `fmax < opt_fmax` (largest
residual force on any atom, Eh/Angstrom). On helicene/20 images the trace is: fmax falls
cleanly 0.998 -> 0.021 in ~100 iterations, then **freezes at 0.0205-0.0212** for the
remaining 1400. It does not diverge or oscillate - it stalls at ~2x the threshold with
the band essentially relaxed. The FIRE timestep collapsed to `dt = 0.000` (constant
`P = v.F < 0` events), i.e. the optimiser was alive but taking no steps. Fixed Jul 2026
by flooring `dt` at 1e-3 and adding stall detection (stops after 200 non-improving
iterations with an explicit warning instead of burning through `opt_iterations`). The
fix makes the behaviour diagnosable; it does **not** improve the barrier (still 65.3).

**GFN2 gives 69.6 kcal/mol** (12 images), also without converging in 800 iterations.
At the time this looked like "GFN2 is worse than GFN-FF (56.0)"; the `opt_k` control
experiment above shows both numbers are unconverged snapshots, so they cannot be
compared at all. What it *did* rule out is the original hypothesis that GFN-FF simply
lacks parameters for the strained inversion geometry.

What was checked and is **not** the cause:

- *Wrong mechanism?* No - but the geometry within it is wrong. The C-frame planarity
  (smallest/largest singular value of the centred carbon coordinates) falls monotonically
  0.399 -> 0.043 at the saddle and rises symmetrically again, i.e. the band does go
  through the flattened helix, which is the accepted inversion mechanism. The defect is
  the over-compressed terminal contact within that structure (H...H 1.64 A, see above).
- *NEB bookkeeping error?* No. Independent single points on the extracted TS geometry
  reproduce the reported barriers exactly (GFN-FF 55.4, GFN2 69.6 kcal/mol) - the
  energies really are that high at those geometries.
- *Endpoint not a minimum?* No. Optimising the input endpoint gains only 2.1 kcal/mol
  (-9.29185 -> -9.29519 Eh with GFN-FF).

So the machinery computes what it claims to compute, and the discrepancy sits in the
saddle-region geometry itself: the band's TS is a *constrained* structure reached by a
path that both methods evaluate as far too strained. Remaining suspects, untested:
insufficient relaxation perpendicular to the path in the saddle region (the spring
constant is uniform, and the FIRE step is unpreconditioned), and the fact that the
climbing image converges last while the surrounding images still move. Until this is
resolved, **do not quote NEB-MD barriers for strained aromatics** - the ethane result
(GFN2 2.88 vs ~2.9 experimental) shows the machinery is correct on a clean case.

## Goal and honest status

The aim is a path-based method that yields **free-energy barriers (dG), not just a
potential-energy saddle**, with (1) a path that reflects the true mechanism and
(2) numbers independent of the number of images.

| requirement | status |
|-------------|--------|
| path CV `s(R)` + analytic gradient | done, validated (FD 1.2e-11) |
| stationary umbrella windows | done (frozen reference path) |
| WHAM estimator | implemented |
| read an external converged MEP (`-nebmd.path_file`) | done |
| tube restraint on `z` (`umbrella_kz`), `dz/dR` FD-verified | done, effect small |
| **usable G(s)** | **not achieved yet** |
| **independence of the image count** | **not achieved yet** |

Current best numbers against a barrier whose reference value is 37.8 (r2SCAN-3c OptTS)
/ 30.4 (GFN2 OptTS) / 40.8 (ORCA GFN2 NEB): PMF gives **58.4** (+45 %), umbrella+WHAM
gives **13-19** (-55 to -70 %). Neither is usable, and they fail in opposite directions.

Two hard obstacles, both measured rather than assumed:

**(a) The generated path can be the wrong mechanism.** For helicene the internal
NEB converges to a *planar* structure (planarity 0.043) with **zero imaginary modes** -
not a saddle at all. A converged ORCA GFN2 NEB on the same system never goes below
planarity 0.24 and gives 40.8 kcal/mol; the true saddles are 30.4 (GFN2 OptTS) and
37.8 (r2SCAN-3c OptTS), both verified by exactly one imaginary mode. Hence
`-nebmd.path_file`: hand in a real MEP instead of trusting the interpolation.

**(b) Umbrella overlap and window holding pull in opposite directions.** The thermal
width of a window is `sigma = sqrt(kT/k)`; WHAM needs `sigma ~ ds/2` (soft k), while
holding an image on a steep flank needs `k*ds/2 > dG/ds` (stiff k). Measured for a
~60 kcal/mol barrier the required stiffness exceeds the overlap limit by a factor of
11 at N=10, 5.3 at N=20, 2.6 at N=40, 1.3 at N=80. **Umbrella windows are not NEB
images** - a usable profile needs far more of them. The code now auto-selects `k`
(`umbrella_k 0`) and reports both failure modes (gaps, drift) instead of silently
returning a number.

### Barriers: where each estimator lands

The target quantity is the **barrier**. Reference values for helicene, all with a
verified saddle (exactly one imaginary mode) or a converged band:

| reference | barrier |
|-----------|---------|
| r2SCAN-3c OptTS (-42.4 cm-1) | **37.8** |
| GFN2 OptTS (-30.6 cm-1) | 30.4 |
| ORCA GFN2 NEB-TS (converged) | 40.8 |
| literature | ~36 |

What this code produces, all sampling the **converged ORCA MEP** with GFN-FF:

| estimator | barrier | error |
|-----------|---------|-------|
| `-nebmd.pmf` (mean force) | **58.4** | +45 % |
| `-nebmd.umbrella` + WHAM | **12.7 - 18.8** | -55 to -70 % |
| NEB HEI, own path | 57 - 82 | not a saddle (0 imaginary modes) |

The two estimators miss in **opposite directions**, and both for a reason that is
directly measurable:

- **PMF too high.** The images drift off the MEP (see below), so they sample *above*
  the path: the reported `<Epot>` barrier is 114.8 kcal/mol where the MEP is 40.8.
  Subtracting the spring bias recovers only part of that.
- **WHAM too low.** The windows slide away from the barrier region (`s0 = 0.571` ->
  `<s> = 0.894`), so the top of the barrier is never visited. A maximum that is never
  sampled does not appear in the histogram, and the profile comes out too flat.

As a cross-check the reaction free energy - which must be exactly 0 for enantiomers -
confirms the same defect: spring-free control **-0.14** (correct), PMF **-46.9** (no
reparametrisation) / **-32.0** (reparam=5), umbrella: windows drift with gaps. Note the
mean-force route gets *worse* on the correct path than on the self-generated one (+4.5),
which rules out "bad path" as the remaining explanation.

The common cause - measuring how far each image moved from its MEP starting point:

```
img  0  0.00   img  4  0.92   img  8  1.10
img  1  0.46   img  5  1.00   img  9  1.08
img  2  0.67   img  6  1.08   img 10  0.99
img  3  0.81   img  7  1.11   img 11  0.80   ... img 13  0.00  (A RMSD)
```

- up to **1.11 A**, largest in the middle of the band where the barrier is. That single
number explains both failures: drifting off the path raises `<Epot>` (PMF too high) and
drifting along it empties the barrier bin (WHAM too low).

So the open problem is not the choice of estimator - both are implemented correctly and
both are undone by the same thing: **the sampling does not stay on the reference path.**
Neither the neighbour springs nor a restraint on `s` alone achieves that for a barrier
this steep.

### 100 ps band MD: the images do NOT leave the path - the estimator is wrong

A full 100 ps band run (14 images on the converged ORCA MEP, GFN2, CSVR, 0.5 fs,
198001 samples per image - by far the best-sampled run in this record) settles where the
error actually is:

| quantity | value | reference |
|----------|-------|-----------|
| PMF barrier | **4.69 kcal/mol** | ~38-41 |
| reaction dG | 2.30 (must be 0) | 0 |
| `<Epot>` barrier | 101.30 | 40.8 |
| "entropic" term | -96.61 | - |

**More sampling makes the mean-force barrier worse, not better**: 58.4 kcal/mol at 400 fs
-> 4.69 at 100 ps. An entropy contribution of -96.6 kcal/mol is unphysical.

Crucially, the drift from the MEP does **not** grow with time:

| run | max drift | mean drift |
|-----|-----------|------------|
| band, 400 fs | 1.11 A | - |
| band, 100 ps (250x longer) | **1.13 A** | **0.75 A** |
| free MD, thermal motion | - | 0.755 A |

The drift is stationary and its mean equals the free molecule's own thermal motion to
three digits. **The images are not running away - the springs hold them in a proper
equilibrium.** That overturns the earlier diagnosis in this document ("the sampling does
not stay on the path"): the sampling is fine, the *estimator* is not.

#### The bias subtraction was wrong (fixed), and what it exposed

`dG/ds = -<F_true.tau> - <F_spring.tau>` treated the spring as an umbrella bias. That is
conceptually wrong: an umbrella bias is an external potential with a **fixed** centre,
while the NEB spring couples to neighbouring images that fluctuate themselves. In a
stationary band the spring is precisely what balances the true force,

```
<F_true . tau> + <F_spring . tau> ~ 0
```

so subtracting one from the other cancels the *signal*, not a bias - which is exactly why
the barrier collapsed as sampling improved (58.4 kcal/mol at 400 fs -> 4.69 at 100 ps).

The estimator now uses `dG/ds = <grad E . tau>` plus the metric term, and reports the
spring balance as a **stationarity diagnostic** instead of folding it in
(den Otter & Briels, JCP 109, 4139 (1998); Fixman, PNAS 71, 3050 (1974)).

Effect on the same 500 fs run on the converged MEP:

| | old (subtracting) | corrected |
|---|-------------------|-----------|
| barrier | 58.40 | 68.02 |
| reaction dG (must be 0) | -32.01 | **-14.94** |

The symmetry error halves. But the new diagnostic immediately shows why the number is
still wrong: `worst |<F.tau>+<Fspring.tau>|/|<F.tau>| = 7.97`, i.e. the spring is a
factor of 8 away from balancing the true force. **The band is not stationary**, so the
premise of the corrected estimator is violated too.

Turning reparametrisation off improves stationarity (7.97 -> 2.01) but lets the arclength
blow up to 414 A as the images bunch together. There is currently **no setting where both
the arclength and the force balance are sound at once** - reparametrisation displaces
images artificially and destroys the balance, while without it the band degenerates. That
is the open problem, and it is a property of the spring-coupled band, not of the
estimator built on top of it.

### Reference scale: what a free 100 ps MD actually explores

Before deciding how tightly to restrain anything, it is worth knowing how much the
molecule moves on its own. A free 100 ps GFN2 MD of helicene at 300 K (CSVR, 0.5 fs,
3962 frames, first 20 % discarded):

| quantity | value |
|----------|-------|
| `<Epot>` | -63.790228 Eh |
| spread | **4.75 kcal/mol** (sd), 34 kcal/mol full range |
| carbon-frame planarity | **0.400 +/- 0.037**, min 0.291 |
| RMSD from the start structure | **0.755 A** mean, 1.03 A max |
| block averages (5 blocks) | spread 0.57 kcal/mol - converged |

Two things follow directly:

**The reaction coordinate is soft.** Thermal fluctuation alone reaches planarity 0.291,
and the true transition state sits at 0.24-0.36. The minimum and the saddle are not
cleanly separated in this coordinate at 300 K.

**Restraint widths must respect the natural motion.** The molecule wanders 0.755 A RMSD
inside its own basin. The image drift measured in the band runs was 1.11 A - only
slightly larger. So part of what earlier sections called "the images leave the path" is
ordinary thermal motion, and a tube narrower than ~0.8 A would freeze out real physics
rather than merely keeping the sampling on the path. The `umbrella_z_tolerance` values
used so far (0.05) were far too tight; a meaningful tube for this system is of order
1 A. This also caps how well any tube restraint can help: the drift that has to be
removed is barely above the thermal noise floor.

(A note on the helicity check: 6 of 3170 frames show the opposite sign with
|helicity| down to 0.63. These are not barrier crossings - at ~36 kcal/mol none are
expected in 100 ps - but grazing passes near the symmetry plane where the sign
criterion becomes numerically unreliable.)

### Tube restraint on `z` and optional neighbour coupling (implemented)

Two additions address exactly this:

- **`-nebmd.umbrella_kz`** - a flat-bottomed restraint on the perpendicular CV:
  `V_z = kz/2 * (z - z_tol)^2` for `z > z_tol`, zero inside. The image may fluctuate
  freely within a tube of radius `-nebmd.umbrella_z_tolerance` around the path and is
  pulled back only outside it. `dz/dR` is analytic and **validated against finite
  differences to 4.8e-14** (`ctest -R test_path_cv`).
- **`-nebmd.umbrella_keep_springs`** - retains the neighbour springs in umbrella mode, so
  an interior image still feels its neighbours instead of being fully independent. The
  endpoints are clamped and exert force without receiving any, which is the intended
  asymmetry.

Assessment of the implementation, stated plainly:

*Correct:* both CV gradients are analytic and FD-verified; the tube is flat-bottomed, so
it does not bias sampling inside the tube; keeping the springs is opt-in and off by
default.

*Compromised:* `umbrella_keep_springs` reintroduces a bias that WHAM does **not** remove -
it trades statistical rigour for staying on the path and is therefore a diagnostic, not a
production setting. The `z` restraint is also not free of consequences: restraining a
coordinate changes the ensemble, and the free energy obtained is then that of the
*tube*, not of the full space. For a narrow tube around a well-defined MEP that is the
intended approximation, but it must be stated rather than hidden.

*Measured effect (converged ORCA MEP, GFN-FF, 400 fs, 14 images):*

| config | WHAM barrier | max window drift in `s` | overlap |
|--------|--------------|-------------------------|---------|
| baseline (`kz=0`) | 16.10 | 0.163 | gaps |
| tube (`kz=5`) | 13.99 | **0.138** | gaps |

The tube restraint reduces the drift by ~15 % and does **not** repair the barrier
(13.99 vs the ~40 reference). Given the free-MD numbers above that is unsurprising: the
drift being removed is close to the thermal noise floor, and `z_tolerance=0.05` was far
below the molecule's natural 0.755 A motion, so the restraint was mostly fighting normal
fluctuation. A physically scaled tube (order 1 A) has not yet been tested.

The perpendicular CV `z(R)` measures
exactly the distance from the path and is already implemented in `PathCV` - it is simply
not wired into the sampling yet. `s` controls progress along the path, `z` would keep the
image on it; only with both is the barrier region actually sampled where it lies.

Measured on the **converged ORCA MEP** (`-nebmd.path_file`, GFN-FF sampling, 8 windows):
the windows still drift by up to 0.323 in `s` at a spacing of 0.143, i.e. more than two
window widths, and one gap remains. So a correct path is *necessary but not sufficient* -
with this few windows the restraint cannot hold them on the flank whatever `k` is. The
next thing to try is 40-80 windows on that MEP, which the auto-`k` formula supports but
which has not been run yet.

Note on what the NEB barrier *is*: it is `E(highest image) - E(first image)`, a
discrete highest-image estimate. With the climbing image it converges onto the saddle,
but only if the band converges - and it can never be strictly independent of the image
count, because it samples the path at discrete points. `G(s)` from WHAM is the
quantity that should be image-count independent; that is the target, not the HEI.

## Free energy / kinetics (EXPERIMENTAL, `-nebmd.pmf`, default OFF)

### Why differentiating and re-integrating E_pot can give a FREE energy

A fair objection: if you take the gradient of the potential energy and integrate it
again, you get the potential energy back - so where does the entropy come from?

The answer is entirely in the **ensemble average**. The free energy is defined by
integrating out all *other* degrees of freedom at fixed arclength,

```
G(s) = -kT ln  INT  delta(s(R) - s) exp(-E(R)/kT) dR
```

and its derivative is

```
dG/ds = < grad E . tau >_s
```

- an average over the hypersurface `s = const`, **not** the gradient at one point.

| what you integrate | what you get |
|--------------------|--------------|
| `grad E . tau` at a single configuration | `E_pot` along the path - no entropy |
| `<grad E . tau>` over the ensemble at fixed s | `G(s)` - includes entropy |

The entropy enters because the image visits *many* configurations at a given `s`. Where
the accessible hypersurface is wide and flat, large gradients cancel in the average and
`<grad E . tau>` is small; where it is narrow, the gradient stays systematically aligned
and the average is large. That width of accessible configuration space *is* the entropic
contribution.

Consistency check in our own data: a rigid trajectory with no fluctuation would give
`<grad E . tau> = grad E . tau` and hence `G(s) = E_pot(s)`, i.e. zero entropy. The short
run (301 samples, images still glued to the starting path) behaves nearly like that,
while the long run separates `<Epot>` and PMF by -58.6 kcal/mol. That the separation is
implausibly large means the average has not converged - but the mechanism is the right
one.

Note also that `<Epot>(s)` and `G(s)` are **not** the same quantity even though both are
ensemble averages: `G = <E> - TS`, and `<grad E . tau>` is not `d<E>/ds`, because moving
`s` also moves the region being integrated over. That change of region is exactly the
entropy. Two points with identical `<Epot>` can have very different `G` if one admits
far more configurations than the other.

### Path collective variables - the general fix for the moving-window problem

The systematic error above comes from `s` being defined by *where the other images
currently are*. A system-specific internal coordinate (a torsion, the helicene
planarity) would fix that but is not general.

The general solution is a **path collective variable** (Branduardi, Gervasio &
Parrinello, *J. Chem. Phys.* **126**, 054103 (2007)), implemented in
`src/capabilities/path_cv.{h,cpp}`. Given a frozen reference path `R_1 ... R_N` it
defines two scalars from the RMSD `d_i = RMSD(R, R_i)` to each frame:

```
s(R) = 1/(N-1) * sum_i (i-1) exp(-lambda d_i^2) / sum_i exp(-lambda d_i^2)
z(R) = -1/lambda * ln sum_i exp(-lambda d_i^2)
```

`s` measures progress along the path (0 -> 1), `z` the distance from it. This is
exactly the "use the RMSD as the metric" idea: it needs a *path*, not a hand-picked
coordinate, so it applies to any system for which a NEB/IDPP path can be produced.

Why it removes the error: the reference frames are **frozen**, so `s` is an explicit
function of `R` alone. The umbrella window no longer migrates during sampling, which
is precisely the defect the spring-free control experiment exposed - and because `s`
is explicit, standard umbrella sampling + WHAM become applicable.

`lambda` is chosen from the mean squared frame spacing (`computeLambdaFromPath()`), the
usual rule of thumb for making neighbouring Gaussians overlap.

Validated by `ctest -R test_path_cv`: `s` rises monotonically along the reference
frames, maps the endpoints near 0 and 1, is deterministic, and the **analytic `ds/dR`
matches finite differences to 1.2e-11**. That gradient check caught a real bug - a
factor of exactly -2, because `RMSDDriver::Gradient()` returns `(ref-tar)/(RMSD*N)`,
i.e. *minus* `d(RMSD)/dR_target`; the sign has to be folded into the chain-rule
prefactor.

**Status**: the CV and its gradient are implemented and validated. Wiring it into a
biased sampling run (umbrella windows on `s` + WHAM) is not done yet.

### Two corrections, and where they stand

**1. The spring bias** - handled in principle, limited in practice. The spring is a known
extra potential, so it is simply subtracted:

```
dG/ds = -<F_true . tau> - <F_spring . tau>
```

which the code does. The formula is exact only for a **stationary** window, though. Our
images migrate (reparametrisation every few steps), so a residual accumulates with
simulation length - precisely what the spring-free control experiment exposed
(-0.14 kcal/mol correct without springs vs +4.52 with them on a long run).

**2. Curvilinear coordinates (Fixman/metric term) - NOT implemented.** Because `s` is a
curvilinear coordinate in 3N space, the volume element is not constant and the estimator
should read

```
dG/ds = <grad E . tau> + kT <d ln sqrt|J| / ds>
```

The second term is missing entirely. It is usually small for cartesian coupling, but the
helicene path is strongly curved, so it is not safely negligible here.

Ranked by what has actually been measured, the error sources are: the over-compressed TS
geometry (~36 kcal/mol, dominant), then the spring bias (~4.5 kcal/mol on long runs),
with the neglected metric term presumably below that.

Motivation: the highest-energy image is a *point* estimate of the potential barrier
at 0 K. A real barrier has a **width** (which sets the TST prefactor and the
tunnelling correction through the curvature at the top) and an **entropic**
contribution from the thermal fluctuations - neither is contained in a single
image. Since NEB-MD already runs N thermostatted replicas, the sampling machinery
for `dG` is in principle present.

What is implemented (opt-in, `-nebmd.pmf true`):

- Mean force along the tangent, accumulated after `-nebmd.pmf_equilibration` steps,
  with the spring (umbrella) contribution subtracted:
  `dG/ds = -<F_true . tau> - <F_spring . tau>`, integrated over the arclength.
- `<base>.neb.pmf.csv` plus an end-of-run evaluation: forward/reverse barrier,
  reaction dG, TS position, FWHM width, curvature at the top, thermal fluctuation,
  and an ASCII profile.
- `-nebmd.mtd_transverse` (default on): when RMSD metadynamics is enabled via
  `-simplemd.rmsd_mtd true`, its bias force is projected perpendicular to the path
  tangent (`SimpleMD::setMTDProjection`), so an image samples its transverse degrees
  of freedom instead of being pushed along the path.

### Control experiment: the springs are the error source

To separate "MD too short" from "springs bias the result", run **two independent,
spring-free MDs** - one from each endpoint - and compare the basins. For enantiomers
every quantity must agree within the statistical error.

Free energy from a plain trajectory: `<Epot>` is just an average, but the **entropy**
depends on the *volume* of phase space visited, not on its mean, so it cannot be read
off directly. For two separate basins the practical estimator is quasi-harmonic
(Schlitter): `S = 0.5 k_B ln det[1 + (k_B T e^2/hbar^2) M sigma]` from the covariance
`sigma` of the mass-weighted, superimposed coordinates - an upper bound, valid within
one roughly harmonic basin (not across a barrier). Then `dG = d<Epot> - T dS`.
(The alternatives - FEP/Zwanzig, thermodynamic integration, umbrella+WHAM - all need a
coupling parameter or fixed windows and do not apply to two plain trajectories.)

Result (helicene, GFN-FF, 6000 steps at 0.5 fs, first third discarded, 402 frames):

| quantity | basin A | basin B | difference |
|----------|---------|---------|------------|
| `<Epot>` | -9.185642 Eh | -9.185129 Eh | -0.32 kcal/mol |
| `T*S_qh` | 31.28 kcal/mol | 31.46 kcal/mol | -0.18 kcal/mol |
| **dG** | | | **-0.14 kcal/mol** |

Against a thermal fluctuation of +/-9.2 kcal/mol that is statistical noise around zero -
correct. Compare with the same system under spring coupling:

| setup | reaction dG (must be 0) |
|-------|--------------------------|
| **independent MD, no springs** | **-0.14 kcal/mol** |
| NEB-MD with springs, short (301 samples) | +0.52 kcal/mol |
| NEB-MD with springs, long (4001 samples) | +4.52 kcal/mol |

So the plain MD is correct and equilibrated, and the error appears **only with the
spring coupling** - growing with simulation length. That is consistent with the
mechanism: the spring bias is subtracted from `dG/ds`, but that subtraction is only
exact for a stationary window. With images migrating (reparametrisation every 5 steps)
the residual accumulates. **The springs, not the MD length and not the energy method,
are the limiting error.**

### What "equilibrated" can and cannot mean here

A *path* has no equilibrium. What can equilibrate is **each image inside its umbrella
window**: it has to sample the degrees of freedom perpendicular to the path at its fixed
arclength position. That is the only condition `<F.tau>` needs. The string itself still
migrates (reparametrisation) - a separate, slower process.

Measured (helicene, 10 images, 6000 steps at 0.5 fs, 2000 discarded), block average over
the first vs the second half of the production frames:

| image | 1st half `<Epot>` | 2nd half | drift |
|-------|-------------------|----------|-------|
| 3 | -9.2192 (sd 0.0034) | -9.2179 (sd 0.0024) | +0.77 kcal/mol |
| 5 | -9.1099 (sd 0.0069) | -9.1120 (sd 0.0007) | -1.35 kcal/mol |

The per-image drift is well below the thermal fluctuation, so **the images are
equilibrated** - short MD is not the limiting factor.

But longer sampling made the *symmetry* worse, not better:

| run | sampling | reaction dG (must be 0) |
|-----|----------|--------------------------|
| short | 300 fs, 150 equilibration, 301 samples | 0.52 kcal/mol |
| long | 3000 fs, 2000 equilibration, 4001 samples | **4.52 kcal/mol** |

That is the important lesson: the good short-run number came from the images still
sticking close to the symmetric starting path, not from convergence. With real sampling
the images drift apart and the asymmetry becomes visible. **A small reaction dG on a
short run is not evidence of a converged profile.** The long run also gives an
implausible entropic term (`<Epot>` barrier 135.5 vs PMF 76.9 kcal/mol, i.e. -58.6
kcal/mol "entropy"), which is a further sign that the mean-force average has not
converged even though the individual images have.

### A pre-optimised path helps, but does NOT fix the estimator (superseded result)

An early 300 fs test suggested that starting from an unconverged CI-NEB was enough:

| starting path | PMF dG(barrier) | reaction dG (must be 0) |
|---------------|-----------------|--------------------------|
| linear interpolation -> MD | 86.2 kcal/mol | 39.5 |
| CI-NEB, not converged -> MD | 71.6 kcal/mol | 0.52 |

**That conclusion did not survive longer sampling.** The same setup run to 3000 fs gives
a reaction dG of **+4.52** instead of 0.52, and on the *converged ORCA MEP* the
mean-force estimator is worse still (**-32 to -47**). The good short-run number came from
the images still sitting on their symmetric starting geometries, not from convergence -
a short run can look symmetric simply because nothing has moved yet.

A pre-optimised path is still worth having (it puts the images in the right valley), but
it does not repair the free-energy estimator. See "Barriers: where each estimator lands".

(The absolute barrier is still too high - see the unresolved helicene/GFN-FF issue
below - but the *sampling* is demonstrably sound only on the pre-optimised path.)

### Two prerequisites (both now available)

1. **Nudging is not conservative.** Projecting out the parallel true force means the
   force field is not the gradient of any potential, so there is no proper NVT
   ensemble - the images heat up (observed: 1857 K at a 300 K setpoint). Free-energy
   estimates therefore require **`-nebmd.nudged false`**, where
   `V = sum_i E(R_i) + (k/2) sum_i |R_{i+1}-R_i|^2` *is* a genuine potential
   (verified: temperatures stay at 300-305 K).
2. **String reparametrisation** (`-nebmd.reparametrize N`, default 0 = off).
   Without it the images drift along the path - they bunch up in the minima and thin
   out over the barrier - so each "umbrella window" keeps moving and the mean force is
   integrated over a shifting grid. `reparametrizeString()` redistributes the images
   to equal arclength every N steps (linear interpolation of the polyline; endpoints
   fixed; velocities rescaled if the temperature ran away).
   Measured on helicene, neighbour spacings along the band:
   - without: `2.7 ... 39.2 A` (path destroyed)
   - with `-nebmd.reparametrize 5`: `1.82 ... 1.98 A` (uniform)

### Validation status (helicene A->B, 10 images, GFN-FF, nudged=false, reparametrize=5)

What is **right** after the fixes: `<Epot>` is symmetric image-by-image
(1<->8: -9.2843/-9.2846, 2<->7: -9.2474/-9.2485, 3<->6: -8.7505/-8.7531,
4<->5: -7.7677/-7.7677) and `dG/ds` is antisymmetric (+0.0042/-0.0041,
+0.0309/-0.0303, ...) - exactly what two enantiomers must give. The barrier top is
found in the middle of the band (image 4 of 1..8), not at an edge.

What is still **wrong**, and why you must not quote these numbers:

- The reaction `dG` comes out as -22 kcal/mol instead of the exact 0 required by
  symmetry. With a perfectly antisymmetric `dG/ds` the integral must vanish; the
  residual is trapezoid error over only 8 interior points on a strongly curved
  integrand. It shrinks with more images.
- The barrier magnitude is unphysical: `<Epot>` rises by ~955 kcal/mol along the
  path, whereas the real helicene inversion barrier is ~35-40 kcal/mol.

  This was checked against two obvious explanations before blaming the path:

  * *Not enough relaxation time?* No. With `k_spring=0.02` the top image equilibrates
    within ~200 steps and then fluctuates around -7.80 +/- 0.01 Eh for 4000 further
    steps - flat, no drift. The high energy is not a transient.
  * *Springs too stiff, image cannot relax sideways?* Partly yes. Loosening to
    `k_spring=0.002` lets the top image settle at -8.13 Eh instead of -7.92 (~130
    kcal/mol lower, stable after ~800 steps). So the restraint does bias the result and
    a soft spring plus a long run matters.

  But even fully relaxed the barrier stays ~730 kcal/mol - an order of magnitude too
  high. The remaining cause is the **path**: interpolating the distance matrix between
  two enantiomers pushes the rings *through* each other, whereas the real inversion is
  a concerted twisting motion. A free-energy profile can only ever be as good as the
  path it is computed on.

So: the machinery (reparametrisation, mean force, spring removal, symmetry) behaves
correctly, but the resulting barrier for this test system is dominated by the poor
path. `-nebmd.pmf` therefore stays **default OFF** and is labelled experimental.
Remaining work for quantitative numbers: a better path (QM method or a
pre-optimised NEB), more images, and block-averaged convergence checks.

## Limitations (v1)

- No climbing image.
- No parallel images (lockstep; `threads = 1` per image for gfnff).
- One shared energy method / thermostat for all images.
- Mild energy drift from the timing caveat above; use a thermostat and short
  relaxation runs, or tune `k` / `dt`.
- Per-image MD restart is not wired (the band state is serialised in
  `curcuma_final.json` but not re-loaded in v1).

## Citations

Registered in the CitationRegistry and emitted automatically by a `-nebmd` run
(see `curcuma_citations.bib`):

- `neb` - Henkelman, G.; Jonsson, H., *J. Chem. Phys.* **113**, 9978 (2000).
  Improved tangent estimate in the NEB method (used for the band tangent).
- `idpp` - Smidstrup, S.; Pedersen, A.; Stokbro, K.; Jonsson, H.,
  *J. Chem. Phys.* **140**, 214106 (2014). IDPP initial path (cited only when
  `-nebmd.idpp` is active, i.e. by default).

The thermostat used by the images cites itself as usual (`csvr`, `berendsen`, ...),
as does the energy method.

## Files

- `src/capabilities/path_cv.h` / `path_cv.cpp` - Branduardi path CVs `s`, `z` with
  analytic gradients; `test_cases/test_path_cv.cpp` validates both against finite
  differences (`ctest -R test_path_cv`).
- `src/capabilities/nebmd.h` / `nebmd.cpp` - the capability.
- `src/capabilities/simplemd.{h,cpp}` - reused stepwise API (not modified in v1).
- `src/capabilities/rmsd.h` - endpoint RMSD reorder + Kabsch alignment.
- `test_cases/cli/nebmd/01_gfnff_ethane/` - CLI test.