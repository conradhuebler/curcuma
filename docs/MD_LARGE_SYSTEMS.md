# MD on a large system (polymer_2x, 7320 atoms)

> 🤖 AI-generated, machine-tested. Re-measured Sep 21, 2026 at commit `4a72e194` on 36 cores /
> 4x RTX A4500. Human production testing pending.

## Why this page exists, and what the problem actually was

`docs/MULTI_GPU.md` recorded that a GFN-FF MD of `test_cases/molecules/larger/polymer_2x.xyz`
heats from 298 K to 20476 K within 40 fs and blamed the structure for not being pre-optimised.
Three rounds of investigation then blamed, in turn, the thermostat, velocity-Verlet truncation
error, and a single collapsing water molecule.

**The actual cause was the clock.** Until commit `ef462fcf` SimpleMD advanced **1.9516 fs per
requested femtosecond**: the step multiplied a velocity in sqrt(Eh/amu) as though it were
Angstrom per femtosecond, so `-dt 1.0` integrated 1.95 fs and `-MaxTime 300` covered 586 fs.
The integrator was internally consistent - energies, forces, temperature and energy conservation
were all correct - which is why nothing caught it for years. See the entry in `AIChangelog.md`
and the ctest `md_time_axis`, which ties the reported time to a Hessian frequency.

The GFN-FF O-H period is 9.17 fs. A 1.95 fs step samples it 4.7 times per period; an
unconstrained X-H bond needs roughly 20. So the page was never about a 1 fs step.

**Everything below was re-measured after the fix. Every time on this page is a real
femtosecond.**

## The headline: 1 fs works

polymer_2x, 7320 atoms, from the GFN-FF-optimised structure, 500 fs, seed 1, 9 threads:

| setting | dE | `<T>` | hottest atom at the end | result |
|---|---:|---:|---:|---|
| `-dt 0.25` NVE | -0.089 Eh | 240 K | - | stable |
| `-dt 0.5` NVE | -0.098 Eh | 233 K | - | stable |
| **`-dt 1.0` NVE** | **+0.049 Eh** | 222 K | 16x | **stable** |
| `-dt 2.0` NVE | - | - | - | diverges, and used to take the process with it (see below) |
| `-dt 0.5` CSVR 300 K | +8.10 Eh | 418 K | **6201x** | **one water collapses at ~480 fs** |
| `-dt 1.0` CSVR 300 K | +6.62 Eh | 360 K | **2121x** | **the same water collapses** |
| `-dt 1.0` CSVR 300 K, `-hydrogen_mass 4` | +4.59 Eh | 296 K | 22x | stable |

"Hottest atom" is the kinetic energy of the hottest atom over the per-atom mean, the local
criterion of `-adaptive_step_local`. A healthy frame of this system sits at 15-23.

Read `dE` only for the NVE rows. **With a thermostat, `dE` is bath work and is supposed to be
nonzero** - the right measure there is whether `<T>` holds the setpoint. NVE settles at 222 K
because the structure starts at a minimum and equipartition puts half the initial kinetic energy
into the potential; the thermostat puts that energy back, which is the +6.6 resp. +4.6 Eh.

So the original goal - **1 fs, unmodified masses, no constraints** - is met in NVE.

**The thermostatted rows are not a step-size problem, and that is the surprise.** Halving the
step makes it *worse* (418 K at 0.5 fs against 360 K at 1.0 fs), and both runs fail at the **same
pair of water molecules** - indices 6522/6523/6524 and 6936/6937/6938. What separates the stable
from the failing runs is not dt but temperature: NVE settles at 222-240 K and holds, the
thermostat keeps 300 K and the site gives way. `-hydrogen_mass 4` holds at 300 K because it
halves the X-H frequency.

That site is a property of *this geometry*, not of the method: water 6522 has **zero oxygen
neighbours within 3.5 A** - it sits free in a cavity, and 342 of the 1904 oxygens in this
structure do. Over 500 fs it crosses 9-11 A and meets another water. A minimised solvated
structure is not an equilibrated one; the cavities are left over from the packing and a local
optimiser cannot close them.

## Using polymer_2x as an MD benchmark

This system exists to measure hardware, not chemistry: how does a GFN-FF or GFN2 MD step scale.
A benchmark runs tens of steps, so none of the long-run behaviour above matters for it - the
first 500 fs are clean in NVE at every step size up to 1.5 fs.

Two structures are kept next to the raw one, both produced with curcuma itself:

| file | method | E | state |
|---|---|---:|---|
| `polymer_2x.xyz` | - | - | raw, as packed |
| `polymer_2x_gfnff_opt.xyz` | GFN-FF | -917.336755 Eh | minimum, max per-atom \|g\| = 0.045 Eh/A |
| `polymer_2x_gfn2_opt.xyz` | GFN2 | -11811.591679 Eh | 400 LBFGS steps, \|grad\| = 0.204, **not converged** |

The recipe. Three step counts and a straight-line fit, so the slope is the pure step cost, the
intercept is the setup and the residual says whether the measurement was worth anything:

```bash
for n in 10 30 50; do
  cp polymer_2x_gfnff_opt.xyz bench.xyz
  time curcuma -md bench.xyz -method gfnff -threads $T -T 300 -dt 1.0 \
       -MaxTime $((n)) -thermostat none -seed 1 -no_bmt -verbosity 0
done
# slope of t(n) = s per step, intercept = setup
```

**Measure on an idle machine.** The first attempt at the table below ran while a GFN2
optimisation held 8 threads and four GPUs, and came out non-monotonic (8 threads slower than 4,
36 slower than 16). That is not a scaling curve, it is contention.

### Baseline: 36 cores (2x Xeon), GFN-FF, 7320 atoms, dt = 1.0 fs

| threads | setup | s/step | speedup | fit residual |
|---:|---:|---:|---:|---:|
| 1 | 22.5 s | 5.05 | 1.00x | +-0.0 s |
| 2 | 13.5 s | 2.95 | 1.71x | +-0.0 s |
| 4 | 8.9 s | 1.93 | 2.62x | +-0.3 s |
| 8 | 5.5 s | 1.45 | 3.48x | +-0.0 s |
| **16** | 6.5 s | **1.25** | **4.04x** | +-0.0 s |
| 24 | 5.8 s | 1.28 | 3.96x | +-1.0 s |
| 36 | 5.0 s | 1.30 | 3.88x | +-0.0 s |

**GFN-FF MD saturates at about 4x around 16 threads and does not improve past it.** By Amdahl
that is a serial fraction of roughly 25 %, and it is the number a GPU has to beat - not the
single-thread time.

### And one GPU does not beat it

| configuration | setup | s/step | against 16 CPU threads |
|---|---:|---:|---:|
| 16 threads, CPU only | 6.5 s | **1.25** | 1.00x |
| 8 threads + 1x RTX A4500 (`-gpu cuda`) | 6.4 s | **1.475** | **0.85x** |
| 16 threads + 1x RTX A4500 | 6.1 s | **1.475** | 0.85x |

One A4500 lands between 8 (1.45) and 16 (1.25) CPU threads, i.e. **the GPU buys nothing for a
GFN-FF MD step at this size** - it is not even breaking even against the CPU saturation point.
Doubling the host threads under the GPU changes the step time by nothing at all (1.475 either
way), so whatever the step is waiting for, it is neither more CPU cores nor more GPU.
The GPU is genuinely used (nvidia-smi shows 100 % utilisation, 1.6 GB), so this is not a silent
fallback; the step simply is not force-kernel bound. `docs/MULTI_GPU.md` reached the same
conclusion from the other side: 1.4 s/step at 7320 atoms against 23 ms/step at 1410 atoms is far
from the ratio the kernels alone would give. The suspects listed there - CN pair-list rebuild
every step, dense EEQ factorisation, per-step synchronisations - are where a faster card would
be wasted.

**What that means for evaluating an H200**: for GFN-FF, the per-step cost has to be attacked in
the code before a faster card can show anything. For GFN2 the picture is the opposite - there
the SCF is dense-linear-algebra bound and memory limited (see below), which is exactly what a
141 GB card changes.

**A trap worth knowing before benchmarking GFN2 here.** A GFN2 optimisation of this system
launched with `-gpu cuda` ran for two hours without completing a single step: it had allocated
19.5 GB on GPU 0, left the GPUs at 0 % utilisation, and was computing on the CPU with 53 GB RSS.
The GPU context for GFN2 at nao = 15444 does not fit on a 20 GB card, the code falls back to the
CPU **silently**, and a CPU GFN2 single point on 7320 atoms takes over an hour. Check
`nvidia-smi` utilisation, not just memory, before trusting a GFN2 GPU timing.

## Where the step limit is, on three atoms and on clusters

A single water, NVE, 5 ps, no constraints. `dE` in Eh:

| start T | dt=0.5 | dt=1.0 | dt=1.5 | dt=2.0 | dt=2.5 |
|---|---:|---:|---:|---:|---:|
| 300 K | -0.0000 | +0.0000 | -0.0001 | -0.0004 | **+1.70** |
| 1200 K | +0.0001 | +0.0002 | -0.0003 | **+1.55** | **+1.07** |
| 1500 K | +0.0002 | +0.0003 | -0.0008 | **+2.24** | **+5.15** |
| 1750 K | +0.0002 | +0.0004 | -0.0009 | **+0.71** | **+12.53** |
| 2000 K | +0.0002 | +0.0004 | -0.0012 | **+1.25** | **+216.42** |
| 2500 K | +0.0001 | +0.0005 | -0.0017 | **+4.55** | **+1.00** |

The limit is between **1.5 and 2.0 fs** and above 1200 K it barely depends on temperature. At
1.0 fs the energy is conserved to 5e-4 Eh even at 2500 K.

The same limit, with the size varied instead of the temperature. Water clusters carved from
polymer_2x, NVE, 200 fs, 300 K, `dE` [Eh] / `<T>` [K]:

| waters | atoms | dt=1.0 | dt=1.5 | dt=2.0 |
|---|---:|---|---|---|
| 40 | 120 | +0.0007 / 191 | -0.0029 / 193 | **+6.30 / 10815** |
| 100 | 300 | +0.0021 / 193 | -0.0097 / 195 | **+11.89 / 7628** |
| 200 | 600 | +0.0040 / 189 | -0.0173 / 192 | **+6.86 / 2324** |
| 400 | 1200 | +0.0074 / 183 | -0.0410 / 186 | **+6742 / 1182349** |
| 800 | 2400 | +0.0055 / 217 | -0.0517 / 215 | **+865 / 74397** |
| 1500 | 4500 | -0.4088 / 222 | -0.5465 / 220 | **+1433 / 65871** |

At 1.0 and 1.5 fs every size conserves; at 2.0 fs every size fails. **There is no size
threshold and no lottery** - an earlier version of this page described a pattern of "40 clean,
100 exploding, 200 clean, 400 and 800 exploding" and read it as a per-oscillator probability.
That pattern was the *marginal* regime at the then-effective 1.95 fs, where the outcome does
depend on the draw. One step further from the limit, it is gone.

What does grow with size is the ordinary accumulation of truncation error, and it is small:
+0.0007 Eh at 120 atoms against +0.0074 at 1200, i.e. 5.8e-6 Eh per atom either way.

**Which limit is this?** The stiffest mode was measured directly (power iteration on the
mass-weighted Hessian, two gradients per iteration) and converges to **3608 cm-1**, an O-H
stretch, period **9.24 fs**. That gives a Verlet *stability* limit `dt < 2/omega` = **2.94 fs**
and the usual *accuracy* rule of thumb `dt <~ T/20` = **0.46 fs**. The measured breakdown sits
at 1.5-2.0 fs, i.e. at about T/5 - conservative against the stability limit and generous
against the rule of thumb. `-hydrogen_mass 4` halves omega and doubles the period, which is why
it buys back the factor of two.

## What breaks when the step is too large

Not the whole system - **one molecule**. In the frame where the `-dt 1.0` CSVR run breaks, **12
of the 7320 atoms hold 75 % of the kinetic energy** and the hottest sits at **2121x** the
per-atom mean, against a steady 15-23x in every healthy frame before it and in the whole NVE and
`hydrogen_mass 4` runs. The atoms are the same ones every time: the waters 6522/6523/6524 and
6936/6937/6938, plus 2100/2101/2102.

At the previously measured collapse the geometry is unambiguous: that water's H-O-H angle has
gone from 103.4° to **3.6°** and its O-H bonds to 0.806/0.859 A. The molecule **alone** is then
worth **+10.96 Eh** against -0.327 Eh in its normal geometry.

**The force field is not at fault.** A rigid H-O-H scan of a single water is smooth and
monotonic across the whole range, including through the 1/sin(theta) point at 180° where the
gradient stays at 0.087 Eh/A: -0.3275 Eh at 103.4°, -0.2847 at 180°, and the other way +0.594 at
20°, +3.311 at 10°, +13.69 at 3.6°, with the gradient rising 0.009 -> 6.6 -> 34.9 -> 261.8 Eh/A.
That is correct H-H repulsion. Reaching 40° already costs 120 kcal/mol, so the geometry is
thermally unreachable and was produced by the integration: at 34.9 Eh/A a 0.5 fs step displaces a
hydrogen by 1.1 A.

**The energy appears in one window, it does not accumulate.** Cutting that water out of every
frame and computing its internal GFN-FF energy gives 0.01 to 0.41 kcal/mol of excitation for the
first 250 fs of the old run - *below* the 0.26 median of twenty other waters - and then 7084
kcal/mol in the next frame.

### A diverging run used to kill the process

At `-dt 2.0` polymer_2x did not merely heat: it crashed inside `SpatialCellList::build()`,
reached from `GFNFF::detectHydrogenBondsNative()`. The cell grid is sized from the bounding box
of the coordinates, that box grows without bound on a diverging trajectory, and the cell count is
its cube - so the allocation cannot be served and the process dies, losing the run together with
its last snapshot.

Fixed in `4a72e194` by giving the grid a budget (cell size doubles until the count fits, in the
worst case down to a single cell). Coarsening is safe: `forEachNeighbor` scans the 3x3x3 block
and filters on the true squared distance, so larger cells only cost time. Reproduced standalone
against both versions of the header, 4 atoms pushed to +-f in all three directions, cutoff 5 A:

| extent | before | after |
|---|---|---|
| 2e3 A | ok | ok |
| **2e4 A** | **abort** (6.4e10 cells requested) | ok |
| 2e5 A | abort | ok |
| 2e6 A | abort | ok |
| NaN | ok | ok |

(The first suspicion was the undefined `static_cast<int>` of a non-finite extent. That turned
out to be harmless by accident - the cast yields a negative number that `std::max(1, ...)`
collapses to one cell - and it is guarded anyway, but it is not what caused the crash.)

## Thermostats

All of them hold the setpoint where the step is safe. 40 waters, `-dt 1.0`, 5 ps, 300 K:

| thermostat | coupling | `<T>` | T at 5 ps | dE |
|---|---:|---:|---:|---:|
| none | - | 281.7 K | 342.4 K | -0.0003 Eh |
| CSVR | 10 fs | 298.2 K | 292.5 K | +0.0201 Eh |
| CSVR | 50 fs | 295.4 K | 274.5 K | +0.0310 Eh |
| CSVR | 200 fs | 300.8 K | 270.5 K | +0.0099 Eh |
| Berendsen | 10 fs | 299.7 K | 278.5 K | +0.0584 Eh |
| Berendsen | 50 fs | 300.5 K | 293.5 K | +0.0110 Eh |
| Nose-Hoover | 10 fs | 299.9 K | 342.1 K | +0.0275 Eh |

So the `-dt 1.0` CSVR failure on polymer_2x is not a thermostat defect. What the thermostat does
is hold the system at 300 K where NVE settles at 222 K, and at 1 fs this system is close enough
to its limit that the difference decides.

**One defect found while checking this, not yet fixed**: `SimpleMD::CSVR()` draws its
chi-squared variate with `m_dof` degrees of freedom where the Bussi-Donadio-Parrinello scheme
uses `m_dof - 1`. The expected kinetic energy then grows by `(1-c)/N_f` per step instead of
staying put. For polymer_2x that is 2.2e-6 per step and cannot explain anything here, but it
scales as 1/N_f, so it is worth 1.6 % per step on a single water. Not measured on a small system
yet.

## The step-rejecting integrator (`-adaptive_step`)

Still useful, but for a different question than before: it does not rescue a correct step, it
rescues a step that is too large. `-dt 2.0` is that regime now.

The design is unchanged and is described in the code (`SimpleMD::IntegratorStep`): measure the
quantity the step is supposed to conserve, and if the step violated it, discard the step and redo
it with a subdivided one. Neither a constraint nor a mass change, and off by default; an explicit
`false` is bit-identical to a binary that never had the feature (`md_adaptive_step` checks that).

Two channels, both calibrated on a running median of the accepted steps rather than on an
absolute energy, because the per-step energy error is a sum over all modes and therefore grows
with the system:

- the **total energy** of the step (`adaptive_step_factor`, default 5), with the kinetic energy
  taken before the thermostat and the bath work subtracted, so the criterion is exact for every
  thermostat;
- the **hottest atom relative to the per-atom mean** (`adaptive_step_local`, on by default,
  `adaptive_step_hot_factor` default 10). This one keeps its contrast at 7320 atoms where the
  global channel loses it, because a violating step stays on a handful of atoms while the
  legitimate fluctuation does not.

Re-measured with the corrected clock, over the first 40 steps of a healthy `dt = 1.0` NVE run
at 300 K:

| system | atoms | median drift | max/median | median hottest-atom ratio |
|---|---:|---:|---:|---:|
| 1 water | 3 | 0.04 kcal/mol | 1.4 | 1.94 |
| 40 waters | 120 | 0.85 | 2.2 | 4.70 |
| 100 waters | 300 | 2.72 | 2.0 | 6.11 |
| 200 waters | 600 | 5.50 | 1.9 | 6.29 |
| polymer | 1410 | 29.02 | 1.4 | 4.53 |
| **polymer_2x** | **7320** | **61.72** | **3.0** | **9.28** |

The median drift spans a factor of 1500 across these systems, which is why an absolute threshold
cannot work. The *relative* spread stays at 1.4-3.0, and the hottest-atom ratio at 1.9-9.3.

On polymer_2x specifically, measured over 60 healthy steps at `dt = 0.5`, the local observable
has the tighter band: drift max/median **3.89** against hottest-atom max/median **1.26**. With
the default factors that puts the global threshold 1.29x above the healthy maximum and the local
one 8x above it - which is the whole argument for the second channel.

> **Not re-measured**: how large the *first* departing step is relative to the healthy band, for
> each system. A first attempt reported the maximum over an entire already-diverged `dt = 2.0`
> run (2.6e6 times the median on polymer_2x), which is not a warning margin and says nothing
> about how much notice a rejection criterion gets. The factors in the code were calibrated
> before the clock fix and are unchanged.

## Two defects found on the way

### 1. The GFN-FF electrostatics cutoff is hard, and it is reachable

`coulomb_r_cut` defaults to 100 Bohr (52.9 A) and the pair loop is a bare
`if (rij > r_cut) continue;` - no switching function. The reference (Fortran `goed_gfnff`) has
**no** cutoff at all; 100 Bohr was picked as "effectively no cutoff" because no pair in any
validation set reaches it. A system wider than ~53 A breaks that assumption.

Measured on a 500-water cluster carved from polymer_2x, moving ONE oxygen by 0.0005 A:

```
dx = +0.0005 A   E = -167.16731610
dx = +0.0010 A   E = -167.16539297     <-- +0.0019 Eh = +5.05 kJ/mol, discontinuous
```

The rest of the scan moves by 5e-6 Eh per step. The jump is 99 % in the Coulomb term, the atom
has 19 partners within +-0.5 A of the cutoff, and `q_O q_H / r` at 100 Bohr is 0.0016 Eh - the
right size. At that geometry the analytic gradient is **0.96 Eh/A wrong** in x and y (finite
difference +0.954 against analytic -0.012). With `-gfnff.coulomb_r_cut 1e9` the scan is smooth
and the same gradient check passes at 1e-6.

`-gfnff.coulomb_r_cut` exposes it; the default stays 100.0, so nothing below ~53 A changes and
every reference set is untouched. Switching it off on polymer_2x costs a single point
3.6 -> 5.4 s. A switching shell (smooth between r_on and r_c) would keep both the performance and
the continuity - not implemented.

**It is not a source of heating.** Re-measured on that 500-water cluster, NVE, 200 fs:

| | dE | `<T>` |
|---|---:|---:|
| `-dt 1.0`, default cutoff | -0.0264 Eh | 223.1 K |
| `-dt 1.0`, `coulomb_r_cut 1e9` | +0.0097 Eh | 223.1 K |
| `-dt 2.0`, default cutoff | +94.2 Eh | 12512 K |
| `-dt 2.0`, `coulomb_r_cut 1e9` | +669.2 Eh | 93227 K |

At a usable step both conserve and the cutoff is worth 0.036 Eh. An earlier version of this page
reported "removing the cutoff takes the drift from +18.8 to +16.0 Eh per 200 fs" - both of those
were measured in the already-diverged regime, where the number says nothing about the cutoff.

### 2. The optimiser can pull a hydrogen into a foreign oxygen

Optimising the same 500-water cluster produced a structure with **H456 0.531 A from O1438** -
half an O-H bond length, and those two waters were **5 A apart** in the input. Two atoms carry
6.5 and 5.7 Eh/A of residual force where the median is 0.0035. An MD started there explodes
immediately.

The potential is not to blame: pushing that hydrogen towards that oxygen with a freshly built
topology raises the energy monotonically, as it must. The optimiser's own energy at its final
geometry (-170.018 Eh) differs from a fresh single point of the same coordinates (-168.695 Eh) by
1.32 Eh - the frozen topology of the iterative path is no longer valid once an atom has changed
its bonding partner. Over 5, 20 and 60 steps the two agree to 0.2 kJ/mol, so this is not a
gradual drift but a cliff the optimiser walks off after a large rearrangement.

**Check the maximum per-atom gradient of any GFN-FF-optimised solvated structure before using
it** - `-dump_gradient` and look at the largest row, not the norm. For reference, the structure
used throughout this page has max |g| = 0.045 Eh/A, RMS 0.0021 and no atom above 0.05.

## Optimising a system of this size

Two convergence controls are **absolute** and therefore unusable at 7320 atoms unless raised:

| Parameter | Default | What happens at 7320 atoms |
|---|---|---|
| `max_energy_rise` | 100 kJ/mol | A single LBFGS step raises the energy by ~260 kJ/mol, so the optimiser aborts immediately with "Energy rise exceeded maximum allowed". |
| `gradient_threshold` | 5e-4 Eh/Bohr | It is a norm over all 3N components. At the relaxed structure the norm is 0.175, i.e. 0.0012 per component - so the criterion can never be met and the run does not terminate. |

The recipe that works:

```bash
curcuma -opt polymer_2x.xyz -method gfnff -threads 36 -optimizer lbfgs \
        -opt.max_energy_rise 20000 -opt.convergence_count 3
# convergence_count 3 = energy + RMSD, dropping the unreachable gradient-norm criterion
```

Restarting the optimiser from its own output (fresh LBFGS memory) is what makes progress - each
round stalls in a line search well before a minimum:

| round | E [Eh] | steps | grad norm |
|---|---|---|---|
| start | -901.94 | - | 2.35 |
| 1 | -913.56 | 27 | 0.79 |
| 2 | -914.20 | 15 | 0.33 |
| 3 | -916.91 | 48 | 0.31 |
| 4 | -916.93 | 5 | 0.17 |
| 5 | -917.34 | 14 | 0.27 |
| 6 | -917.337 | 5 | 0.175 |

Both defaults should arguably scale with system size; that is a code change, not a documentation
matter, and is not made here.

## GFN2

Starting GFN2 from the GFN-FF-optimised structure is worth **14.65 Eh (38500 kJ/mol)**: GFN2
gives -11784.88 Eh on the raw structure and -11799.53 Eh on the GFN-FF-optimised one, before a
single GFN2 step. At ~160 s per step (4x A4500) a GFN2 optimisation of this system is a
multi-hour job, so pre-optimising with GFN-FF is not a convenience but the only affordable route.

A GFN2 optimisation is in progress; 400 iterations took it to -11811.59 Eh with a gradient norm
of 0.204, moving the structure by 1.99 A on average and 7.46 A at most from the GFN-FF minimum -
mostly the 1500 waters rearranging. It has not converged.

> **Watch out for an optimiser result that is the input.** Before the Sep 2026 fix, a
> non-converged single-structure `-opt` wrote `molecules[0]` - the *input* geometry - dressed as
> a result, with a stale energy in the comment line. A 6.8 h GFN2 run was lost that way. It now
> writes `result.final_molecule`; if you have `.opt.xyz` files from before that, check them
> against their input before trusting them.

## What was not tested

- Longer than 500 fs, and no production trajectory - these runs answer "does it stay at the
  setpoint", not "is the sampling correct".
- `-dt 0.75` and `-dt 1.25`, so the limit is bracketed between 1.5 and 2.0 fs and not resolved.
- `-hydrogen_mass` at values other than 4, and its effect on the dynamics (it changes the time
  scale of X-H motion by construction) was not examined.
- GFN2 MD at this size.
- The `adaptive_step` calibration with the corrected clock (see the note in that section).
- Whether the CSVR chi-squared off-by-one matters on a small system.
