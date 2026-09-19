# Stable MD on a large system (polymer_2x, 7320 atoms)

> 🤖 AI-generated, machine-tested. Measurements from Sep 18, 2026 on 36 cores / 4x RTX A4500,
> commit 05d40030. Human production testing pending.

`docs/MULTI_GPU.md` recorded that an MD of `test_cases/molecules/larger/polymer_2x.xyz` heats
from 298 K to 20476 K within 40 fs, and attributed it to the structure not being pre-optimised.
That attribution was wrong. This page is what the measurement says instead.

## The settings that work - over 100 fs, and that is not long enough

| Setting | `<T>` over 100 fs, target 300 K |
|---|---:|
| `-dt 1.0` (the default) | **18527 K** - blows up |
| `-dt 0.5` | 288.7 K |
| **`-dt 0.25`** | **296.1 K** |
| `-dt 1.0 -hydrogen_mass 4` | 273.1 K |

> **Correction (Sep 19, 2026): 100 fs is too short a window, and `-dt 0.5` does not hold.**
> Re-running the same case for **300 fs** (GFN-FF-optimised structure, CSVR at 300 K, seed 1)
> gives `<T>` = **716.6 K**: the running-average potential goes -916.30 -> -913.76 Eh and the
> kinetic 9.55 -> 24.91 Eh, i.e. the total energy rises by about **+8 Eh** *while the thermostat
> is actively trying to remove it*. The potential itself only oscillates (-917.3 at 0 fs, -911.5
> at 100 fs, -918.5 at 300 fs) - it is the kinetic energy that grows monotonically, and the
> 100 fs row above simply ends before the effect is visible. `-dt 0.25` and `-hydrogen_mass 4`
> were **only** measured over 100 fs and are therefore equally unproven at 300 fs. **Do not read
> any row of this table as a validated production setting for this system.**
>
> **And the thermostat is not the problem.** The same case without one - plain NVE, dt = 0.5,
> 300 fs - gains **+67.6 Eh** and reaches `<T>` = **2175 K** (from 273). So energy is created by
> the integration itself, not merely left in by a thermostat that cannot keep up. The potential
> stays flat throughout (-916.3 at the start, -915.2 at the end, never outside -917…-912): all
> of it goes into kinetic energy, i.e. the atoms simply get faster while the structure holds.
>
> **It scales as dt².** NVE gains +191 Eh in 200 fs at dt = 1.0 and +67.6 Eh in 300 fs at
> dt = 0.5, which is +45 Eh per 200 fs - a ratio of **4.24** against the expected 4.
>
> > **Correction (Sep 19, 2026, later the same night): that dt² ratio is real but the conclusion
> > drawn from it - "ordinary velocity-Verlet truncation error" - was wrong.** Truncation error
> > is a property of the whole system and grows with it; this is one molecule. Four measurements
> > settle it, all NVE, dt = 0.5, 200 fs, seed 1, from the same GFN-FF-optimised geometry:
> >
> > | system | atoms | dE |
> > |---|---:|---:|
> > | water cut out of polymer_2x, alone | 4500 | **-0.176 Eh** |
> > | the polymer alone, water removed | 2820 | **+0.021 Eh** |
> > | both together | 7320 | **+45 Eh** |
> > | pure water, 300 / 1200 / 2400 atoms | | +0.002 / +0.009 / +0.165 Eh |
> >
> > Each half conserves the energy, and 4500 atoms of water conserve it as well as 300 do - so it
> > is neither the size nor a generic integration error. **It is a single water molecule
> > collapsing.** In the frame where the run breaks, 12 of the 7320 atoms hold **99.0 %** of the
> > kinetic energy, and the hottest two are the hydrogens of one water (indices 6522/6523/6524);
> > in every healthy frame before it the hottest atom sits at a steady 15-23x the per-atom mean,
> > at the event it reaches **3733x**. That molecule's H-O-H angle has gone from 103.4° to
> > **3.6°** and its O-H bonds to 0.806/0.859 A; the molecule **alone** is then worth **+10.96 Eh**
> > against -0.327 Eh in its normal geometry, which is the entire energy gain of the run.
> >
> > **The force field is not at fault.** A rigid H-O-H scan of a single water (3 atoms, GFN-FF)
> > is smooth and monotonic across the whole range, including through the 1/sin(theta) point at
> > 180° where the gradient stays at 0.087 Eh/A: -0.3275 Eh at 103.4°, -0.2847 at 180°, and
> > going the other way +0.594 at 20°, +3.311 at 10°, +13.69 at 3.6° with the gradient rising
> > 0.009 -> 6.6 -> 34.9 -> 261.8 Eh/A. That is correct H-H repulsion. Reaching 40° already costs
> > **120 kcal/mol** and the observed geometry far more - 250 to 430 kT for a single molecule at
> > the 245 K the run was actually at, i.e. thermally unreachable. So the collapse was driven by
> > the integration, and once the angle is below ~30° every further step is unresolvable: at
> > 34.9 Eh/A a 0.5 fs step displaces a hydrogen by 1.1 A.
> >
> > Two candidate explanations were tested and **eliminated**: `-cleanenergy true`, which rebuilds
> > the energy calculator every step, makes it *worse* (+173.7 Eh against +53.2 on a 300-atom
> > cluster), so it is not stale state between steps; and the frozen setup topology gives bit-
> > identical forces to a freshly perceived one at the breaking geometry (max|g| 176.1013 Eh/A
> > either way), so it is not a missing pair in a cutoff list. On a healthy 120-atom system the
> > MD's own energy and a fresh single point on the written geometry agree to 6e-6 Eh.
> >
> > **And the dt² ratio itself does not survive a matched measurement.** Run the full 7320-atom
> > system at dt = 0.5 for **200 fs** (NVE, seed 1) and it gives **+0.069 Eh, `<T>` = 225.6 K** -
> > conserved. The "+45 Eh per 200 fs" above was not measured; it was the 300 fs number scaled by
> > 200/300, which silently assumes the drift accumulates linearly. It does not: for 250 fs
> > nothing happens, then one molecule collapses and the run gains 67 Eh in the last 50. So the
> > ratio 4.24 compared a real 200 fs run at dt = 1.0 against a scaled number dominated by a
> > single event, and it means nothing. **Never scale a drift to a different window without
> > checking that it is linear in time.**
> >
> > What this changes practically: the step needed is set by the *worst* local event, not by an
> > arithmetic that applies everywhere, so lowering dt globally is the wrong lever. The right one
> > is to reject the individual step - see the local criterion in the next section.

```bash
# what we verified, from the GFN-FF-optimised structure
curcuma -md polymer_2x_gfnff_opt.xyz -method gfnff -threads 36 -T 300 \
        -dt 0.25 -thermostat csvr
# same stability at the full 1 fs step, by making the X-H oscillation slower:
curcuma -md polymer_2x_gfnff_opt.xyz -method gfnff -threads 36 -T 300 \
        -dt 1.0 -hydrogen_mass 4 -thermostat csvr
```

## Two defects found while investigating this (Sep 19, 2026)

The time-step story below is the *practical* recipe. It is not the whole truth: two genuine
defects came out of chasing it, and both matter beyond this system.

### 1. The GFN-FF electrostatics cutoff is hard, and it is reachable

`coulomb_r_cut` defaults to 100 Bohr (52.9 A) and the pair loop is a bare
`if (rij > r_cut) continue;` - no switching function. The reference (Fortran `goed_gfnff`) has
**no** cutoff at all; 100 Bohr was picked as "effective no-cutoff" because no pair in any
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

`-gfnff.coulomb_r_cut` now exposes it; the default stays 100.0, so nothing below ~53 A changes
and every reference set is untouched. The price of switching it off on polymer_2x (7320 atoms)
is a single point 3.6 -> 5.4 s. A switching shell (smooth between r_on and r_c) would keep both
the performance and the continuity - not implemented yet.

**This is not the main cause of the heating**: removing the cutoff took the 500-water NVE drift
from +18.8 to +16.0 Eh per 200 fs. It is a real defect found on the way, not the explanation.

### 2. The optimiser can pull a hydrogen into a foreign oxygen

Optimising the same 500-water cluster produced a structure with **H456 0.531 A from O1438** -
half an O-H bond length, and those two waters were **5 A apart** in the input. Two atoms carry
6.5 and 5.7 Eh/A of residual force where the median is 0.0035. An MD started there explodes
immediately (5e6 K within 200 fs), which is *worse* than starting from the unoptimised cluster.

The potential itself is not to blame: pushing that hydrogen towards that oxygen with a freshly
built topology raises the energy monotonically, as it must. The optimiser's own energy at its
final geometry (-170.018 Eh) differs from a fresh single point of the same coordinates
(-168.695 Eh) by 1.32 Eh - the frozen topology of the iterative path is no longer valid once an
atom has changed its bonding partner. Over 5, 20 and 60 steps the two agree to 0.2 kJ/mol, so
this is not a gradual drift: it is a cliff that the optimiser walks off after a large
rearrangement. **Check the maximum per-atom gradient of any GFN-FF-optimised solvated
structure before using it** - `-dump_gradient` and look at the largest row, not the norm.

## The actual mechanism, on three atoms (Sep 19, 2026)

Everything above was measured on systems of 300 to 7320 atoms. The cause reproduces on **one
water molecule**, and that is where it should be read.

**A single water, unconstrained, dt = 1 fs, NVE over 5 ps:**

| start T | dt = 1.0 fs | dt = 0.5 fs |
|---|---:|---:|
| 1200 K | -0.0024 Eh | +0.0004 Eh |
| **1500 K** | **+4.54 Eh** | +0.0006 Eh |
| 1750 K | **+9.69 Eh** | +0.0007 Eh |
| 2000 K | **+2.24 Eh** | +0.0007 Eh |

There is a **threshold in amplitude**, not a drift: below it the energy is conserved for 10 ps
(1 water +0.0007 Eh, 40 waters -0.013 Eh - the same magnitude at 2 ps and at 10 ps, i.e. the
integrator is a correct symplectic Verlet and pumps nothing), above it the run explodes.

**Why the threshold is at dt = 1 fs.** The GFN-FF O-H bond stiffens steeply under compression.
Curvature of E(r) along the bond, converted to a frequency and the Verlet limit `dt < 2/omega`:

| r(O-H) | E - E0 | omega | dt_stab |
|---|---:|---:|---:|
| 0.70 A | 0.144 Eh | **10758 cm-1** | **0.99 fs** |
| 0.75 A | 0.080 Eh | 8773 cm-1 | 1.21 fs |
| 0.80 A | 0.040 Eh | 7139 cm-1 | 1.49 fs |
| 0.90 A | 0.005 Eh | 4753 cm-1 | 2.23 fs |
| **0.971 A (eq)** | 0 | **3817 cm-1** | **2.78 fs** |
| 1.15 A | 0.020 Eh | 2669 cm-1 | 3.98 fs |

At the equilibrium length dt = 1 fs has a factor 2.8 of margin. Compressed to 0.70 A the margin
is **gone** - and then the feedback closes: more compression -> higher omega -> larger
integration error -> more energy -> more compression.

**Why that makes large systems fail and small ones not.** The per-oscillator probability of
crossing the threshold at 300 K is small but not zero. polymer_2x has **4612 hydrogens**, so
over a few hundred steps the crossing is a certainty; 40 waters (80 O-H) usually get away with
it. The observed pattern - 40 waters clean, 100 waters exploding, 200 clean, 400 and 800
exploding - is exactly that lottery, and it is **not** a size effect in the code. It was
verified step by step in the 100-water case: one O-H oscillation grows over ~10 fs
(1.06 -> 0.81 -> 1.28 -> 0.72 A), the hydrogen then leaves its own oxygen (1.87 A) and hits a
neighbouring one at 0.42 A, where the EEQ charges diverge to -5.45/+4.53 e and produce a force
of order 100 Eh/A. One Verlet step with that force injected 262 Eh into a system whose entire
kinetic energy was 0.43 Eh.

**What was ruled out on the way** (each by measurement, not by argument): the starting structure
(optimised heats *worse* than raw), outlier forces (max |g| 0.046 vs 0.079 in a stable system),
a thread race (gradients agree to 1.7e-14 between 1 and 36 threads), an energy/force
inconsistency (finite differences agree to 1e-6 at the start geometry *and* inside the heating
trajectory), the Coulomb cutoff (removing it makes the drift *worse*: +22.97 vs +11.58 Eh),
water as such (1 and 40 waters conserve energy at dt = 1 fs), stale internal MD state (frozen
topology reproduces the fresh single point to the last digit at the collapsed geometry), and a
CPU/GPU difference (that was the Eh/Bohr unit bug, fixed earlier).

**A guard that was tried and REJECTED.** The diverging EEQ charges are a known pathology of
electronegativity equalisation at short range - the off-diagonal erf(gamma*r)/r rises towards
the diagonal hardness and the 2x2 block becomes near-singular. The reference has the same
property; xtb never meets it because its MD defaults are `hmass=4` and `shake=2` (all bonds),
with a 4 fs step. Rejecting an EEQ solution with |q| above 4 e + |molecular charge| and falling
back to the topology charges (the rule the xTB SCF uses since Known Issue #9) does cut the
damage by a factor of 5 (100 waters at dt = 1 fs: 20953 K / +29.7 Eh -> 4036 K / +5.5 Eh) -
**but it is not in the code, because it breaks the reference sets**: GMTKN55 gfnff went to
MAD 0.198 / max 166 kcal/mol, and the worst case, `IL16/229`, legitimately carries 4.63 e. The
4 e bound is safe for an SCF and is not safe for EEQ topology charges. The patch was reverted
and GMTKN55 verified back to zero differences. It would also have been a safety net rather than
a cure: the step that put two atoms 0.42 A apart was already unphysical.

## The fix: a step-rejecting integrator (`-adaptive_step`, Sep 19, 2026)

Everything above says the same thing: **single steps** violate the integrator's accuracy limit,
and one such step is enough to ruin a trajectory. Lowering `dt` globally or raising the hydrogen
mass both work, but they pay for a handful of bad steps with every step of the run, and the mass
trick changes the dynamics by construction.

Step rejection is the ordinary numerical answer and neither a constraint nor a mass change:
measure the quantity the step is supposed to conserve, and if the step violated it, throw the
step away and redo it with a subdivided time step. The physics is untouched - a smaller step is
still velocity-Verlet on the same potential - and the cost is paid only where it is needed.

```bash
curcuma -md polymer_2x_gfnff_opt.xyz -method gfnff -threads 36 -T 300 \
        -dt 1.0 -thermostat csvr -adaptive_step true
```

**It is off by default** and an explicit `false` is bit-identical to a binary that never had the
feature (`md_adaptive_step` checks exactly that). Full `ctest` after the change: 496 tests, 485
passed, 7 disabled and **4 failed - the same four this repository already had**
(`confscan_dtemplate` flaky, `test_orca_interface`, `xtb_cpscf`,
`cli_curcumaopt_07_opt_multixyz` golden-value drift). No new failure.

### What is measured, and against what

The conserved quantity is `E_pot + E_kin`, with the kinetic energy taken **before** the
thermostat touched the velocities - the thermostat legitimately changes `E_kin` and must not
count as a violation. When a thermostat is active its work over the step is measured and
subtracted, so the criterion is exact for every thermostat, not only CSVR.

### Why the threshold is not a fixed energy

The energy error of velocity-Verlet is not a drift but a bounded oscillation whose amplitude is
the sum over all modes, so it grows with the system. Median per-step `|dE|` at 300 K, dt = 1 fs:

All numbers below are over the **first 40 steps** of the same run, so median and maximum come from
the same window; the last two columns are the worst step of the whole 150-step run.

| system | atoms | median | p90 | max | max/median | worst step | /median |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 water | 3 | 0.35 | 0.52 | 0.58 | **1.7** | 0.6 | 1.7 |
| 40 waters | 120 | 6.27 | 13.7 | 21.4 | **3.4** | 21.4 | 3.4 |
| 100 waters | 300 | 17.2 | 33.3 | 56.6 | **3.3** | 32943.6 | **1919** |
| 200 waters | 600 | 33.5 | 62.7 | 105.3 | **3.2** | 71411.6 | **2134** |
| polymer | 1410 | 195.1 | 255.7 | 326.8 | **1.7** | 326.7 | 1.7 |

(kcal/mol. Measured with `CURCUMA_ADAPTIVE_DEBUG=1`, which prints the drift of every step.) A
fixed number would reject every step of the polymer or no step of the water. What **is**
system-independent is the spread *within* a run: the largest healthy step is 1.7 to 3.4 times the
median across all five systems, while the step that destroys a trajectory is **1919x** resp.
**2134x** it. Almost three orders of magnitude separate the two, and the two systems that survive
the 150 steps never produce anything above their own healthy maximum.

The threshold is therefore `adaptive_step_factor` (default 10) times the **running median of the
steps accepted so far**, capped at `adaptive_step_tol` (default 1.0) times the thermal energy
`N_dof*kB*T/2`. The cap also covers the warm-up, before enough steps exist to form a median; the
healthy maximum measured over all five systems is 0.65 of it, the destructive steps 124 and 134.

**Only a step accepted on the first attempt calibrates the median.** A subdivided step has a
smaller drift than the one it replaced but still a larger one than an ordinary step, so feeding
those back lets a degenerating trajectory raise its own threshold. That was measured, not
argued: on the 100-water cluster the unanchored version left +1.50 Eh with 89 rejections, the
anchored one -0.01 Eh with 8.

### What it does

200 fs of NVE at `dt = 1.0` from the same structure and the same seed, so a row with zero
rejections must reproduce the `off` column exactly - and does. `dE` in Eh, `<T>` in K, the
number in brackets is how many of the 200 steps (120 for the polymer) were redone subdivided.

| system | atoms | off | factor 1.5 | factor 2 | factor 3 | **factor 5** | factor 10 | factor 20 |
|---|---:|---|---|---|---|---|---|---|
| 1 water | 3 | -0.0005 | - | - | - | **-0.0005 (0)** | -0.0005 (0) | -0.0005 (0) |
| 40 waters | 120 | -0.0163 | -0.0166 (74) | -0.0127 (44) | -0.0080 (28) | **-0.0163 (0)** | -0.0163 (0) | -0.0163 (0) |
| polymer | 1410 | +0.0999 | -0.0708 (1) | +0.0999 (0) | +0.0999 (0) | **+0.0999 (0)** | +0.0999 (0) | +0.0999 (0) |
| 100 waters | 300 | **+53.26** / 37303 | +0.051 (143) | +0.222 (129) | +0.507 (105) | **+1.104 (107)** | +2.808 (77) | +3.035 (73) |
| 200 waters | 600 | **+81.37** / 28710 | +0.043 (145) | +0.087 (80) | +0.243 (73) | **+0.456 (75)** | +0.952 (20) | +3.578 (53) |

Read the table as two halves. The upper three systems are healthy at 1 fs: at the default
factor the feature rejects **nothing** there and the trajectory is bit-identical to a run
without it - which is the point, since a rejection on a healthy system is pure cost (40 waters
at factor 1.5 pays 74 rejections for a `dE` that gets *worse*, -0.0166 against -0.0163). The
lower two destroy themselves without it and are brought back to the setpoint with it.

**Why the default is 5**: it is the smallest factor that rejects nothing on any healthy system
measured, while still cutting the two blown-up runs by a factor of 50 to 180. Tightening it to
**2** buys another factor of 5 on a system that still gains energy (100 waters +1.10 -> +0.22 Eh)
and costs 22 % of the steps of a healthy 120-atom system. Loosening it is not useful: 10 and 20
reject *fewer* steps but conserve *worse*, because the steps they let through accumulate.

The false positives are a small-system effect, not a general one: the ratio of the largest
healthy step to the median is 3.4 for the 120-atom cluster but only 1.7 for the 1410-atom
polymer, so the same factor is far more generous on the large system - at factor 1.5 the polymer
rejects exactly one step out of 120.

### The size limit of the GLOBAL criterion (superseded by the local one below)

> **Superseded (Sep 19, 2026).** Everything in this section is still the correct description of
> the **global energy** criterion and of why it fails at 7320 atoms - keep it, it is the
> measurement the local channel was designed from. What is no longer true is the conclusion:
> the local criterion of the following section does rescue polymer_2x
> (+67.59 -> +0.53 Eh, 2175 -> 246 K), and 92 of its 104 rejections are ones the global
> criterion accepted.

On `polymer_2x` itself, 7320 atoms, the **global** criterion does not help - it makes the run
**worse**:

| t | `-dt 1.0` without | `-dt 1.0` with `-adaptive_step` (factor 5) |
|---|---:|---:|
| 0 fs | -917.34 Eh | -917.34 Eh |
| 50 fs | -911.96 | -909.78 |
| 100 fs | **-911.08** | **-894.57** |

(E_pot, 200 fs NVE, same seed.) And it costs about 45 s per step instead of a few, because
roughly half the steps are subdivided.

The reason is measurable and it is not a bug. The criterion works by contrast: the destructive
step has to stand out from the legitimate per-step fluctuation. That fluctuation is a sum over
all modes, so it **grows with the system**, while a local defect - one O-H bond collapsing -
stays local. Healthy phase, first 40 steps, dt = 1 fs, 300 K:

| system | atoms | median drift | max / median (healthy) | destructive step / median |
|---|---:|---:|---:|---:|
| 40 waters | 120 | 0.010 Eh | 3.4 | - |
| 100 waters | 300 | 0.027 Eh | 3.3 | **1919** |
| 200 waters | 600 | 0.053 Eh | 3.2 | **2134** |
| polymer | 1410 | 0.311 Eh | 1.7 | - |
| **polymer_2x** | **7320** | **0.505 Eh** | **5.94** | - |

At 7320 atoms the **healthy** spread alone is 5.94 times the median - wider than the default
factor of 5. The threshold therefore sits *inside* the legitimate distribution: about half the
steps are rejected essentially at random, which changes the trajectory without removing the one
event that matters, and a local collapse worth a few Eh is indistinguishable from a normal
fluctuation of 0.5 to 3 Eh. The contrast is ~600 at 300 atoms and ~1 at 7320.

**So the rule is**: a global energy criterion discriminates while the system is small enough
that one bad bond dominates the total energy error - measured here up to ~1400 atoms. Beyond
that it needs to be **local** (per atom, per bond, or per fragment).

### The local criterion (`-adaptive_step_local`, on by default, Sep 19, 2026)

That local channel now exists. The observable is the plainest one available: the kinetic energy
of the **hottest atom divided by the per-atom mean**. It needs no decomposition of the potential,
it is dimensionless, and it is almost size-independent, because the maximum of N samples of a
chi-squared distribution grows only logarithmically in N. It is calibrated by its own running
median over accepted steps, exactly as the energy drift is, so no absolute number is baked in.

Why it separates where the total energy does not - both observables measured over the same 60
healthy steps of polymer_2x, 7320 atoms, dt = 0.5:

| observable | median | healthy maximum | max / median | default factor | threshold sits at |
|---|---:|---:|---:|---|---|
| total energy drift | 48.96 kcal/mol | 190.45 kcal/mol | **3.89** | 5 | 245, i.e. **1.29x** the healthy maximum |
| hottest atom / mean | 9.70 | 12.19 | **1.26** | 10 | 97, i.e. **8x** the healthy maximum |

The global threshold sits at the edge of the legitimate distribution - that is the failure
documented just above, in one number. The local one has room to spare in both directions, and
the event it has to catch is not marginally above the healthy range but three orders of
magnitude above it: in the frame where the run breaks, the hottest atom reaches **3733x** the
per-atom mean against a steady 15-23x in every frame before it (measured on a 25 fs frame
proxy, so the per-step contrast is larger still).

It is a pure addition: it can only reject a step that the global criterion accepted, never
accept one it rejected. On the deterministic single-water regression case the trajectories with
and without it are **bit-identical to twelve decimals** (dE = -0.007933244988 Eh,
`<T>` = 641.099708218 K either way), and on a 300-atom cluster, where the global criterion
already works, the local one fires **zero** times - it does not misfire on healthy systems.

**What it does on the system this page is about.** polymer_2x, 7320 atoms, 300 fs NVE at
dt = 0.5, same seed, same thread count:

| | dE | `<T>` | rejections |
|---|---:|---:|---|
| without | **+67.59 Eh** | **2175.1 K** | - |
| `-adaptive_step true` | **+0.53 Eh** | **246.2 K** | 104 of 600 steps, **92 of them caught by the local channel alone** |

Those 92 are the measurement that matters: the global energy criterion would have accepted
them. It is also visible in the trajectory - the hottest-atom ratio runs at a flat 15-23 through
frame 12 in both runs (the accepted steps are the same ones), and in the frame where the plain
run reaches 3733 with 12 atoms holding 99.0 % of the kinetic energy, the rejecting run reaches
911 with those atoms holding 24.2 %. So the event is **damped, not eliminated**, which is what
the residual +0.53 Eh is. The run costs about 1.9x the wall time of the plain one.

`-adaptive_step_hot_factor` (default 10) sets the multiple of the running median.
`-adaptive_step_local false` restores the previous behaviour exactly.

**The factor is not the limiting element here, which was checked rather than assumed.**
Tightening it from 10 to 4 on the same polymer_2x run changes nothing measurable: +0.5189 Eh /
245.8 K against +0.5305 Eh / 246.2 K, the same 104 rejections, 85 local-alone instead of 92. So
the residual is structural - the event is damped rather than prevented - and not a threshold
that was set too loose. Do not tune the factor expecting the residual to move.

### What it does not do

- The **global** channel does not scale to arbitrary system size - see the section above. The
  local one was measured to 7320 atoms and nowhere beyond; the observable is only weakly
  size-dependent by construction, but that is an argument, not a measurement.
- It damps a violating event, it does not undo one. On polymer_2x the rejecting run still shows
  the event at 911x the per-atom mean instead of 3733x, and keeps +0.53 Eh of the +67.59.
- It cannot rescue a step that is already unphysical in the potential, only one that is
  unphysical in the integration. A geometry with two atoms 0.42 A apart is wrong either way.
- It costs one extra force evaluation per subdivided step times the number of substeps
  (`adaptive_step_substeps`, default 8), so a run in which most steps are rejected is slower
  than simply halving `dt`. The rejection count is reported at the end of the run; if it is a
  large fraction of the steps, `dt` is wrong for the system and the feature is only papering
  over it.
- It has not been tested with RATTLE, with metadynamics, or across a restart.

## Why it is the time step and not the structure

Three measurements, each of which the "bad structure" explanation fails:

1. **Optimising does not help.** A GFN-FF optimisation takes polymer_2x from -901.94 to
   **-917.34 Eh** and the gradient norm from 2.35 to **0.175**. At `-dt 1.0` that optimised
   structure still heats to `<T>` = 18527 K - **worse** than the raw structure's 2936 K.
2. **The step size does help**, in the same structure and the same thermostat: 0.5 fs and
   0.25 fs both stay at the setpoint.
3. **The NVE drift has a threshold between 0.5 and 1.0 fs.** Over 100 fs without a thermostat:

   | dt | steps per period | dE |
   |---|---:|---:|
   | 0.25 fs | 37 | -0.25 Eh |
   | 0.50 fs | 18 | **-0.05 Eh** |
   | 1.00 fs | 9.2 | **+11.58 Eh** |

**Which limit is this?** The stiffest mode was measured directly - power iteration on the
mass-weighted Hessian, two gradients per iteration, no Hessian needed - and converges to
**3608 cm⁻¹**, a plain X-H (here O-H of water) stretch, period 9.24 fs. That gives

- Verlet **stability** `dt < 2/omega` = **2.94 fs**, so 1 fs has a factor 2.9 of margin, and
- **accuracy** `dt <~ T/20` = **0.46 fs**, which 1 fs misses by a factor 2.

So dt = 1 fs is stable but not accurate here, and the measured table is exactly what that
predicts: clean at 18 steps per period, drifting at 9. `-hydrogen_mass 4` halves omega and
doubles the period, moving `T/20` to 0.92 fs - which is why it restores the 1 fs step.

An earlier version of this page called the 1 fs behaviour "beyond the stability limit". That
was wrong; the stability limit is 2.94 fs and was never reached.

**Size matters, but only as the trigger.** The same NVE at dt = 1.0 on `polymer.xyz` (1410
atoms) drifts **+0.10 Eh**, against +11.58 Eh at 7320 atoms - 22x more per atom for the large
system. So 1 fs is marginal for GFN-FF generally and crosses into instability as the system
grows. Do not read a stable small-molecule MD as evidence that the setting is safe.

## Optimising a system of this size

Two convergence controls are **absolute** and therefore unusable at 7320 atoms unless raised:

| Parameter | Default | What happens at 7320 atoms |
|---|---|---|
| `max_energy_rise` | 100 kJ/mol | A single LBFGS step raises the energy by ~260 kJ/mol, so the optimiser aborts immediately with "Energy rise exceeded maximum allowed". |
| `gradient_threshold` | 5e-4 Eh/Bohr | It is a norm over all 3N components. At the relaxed structure the norm is 0.175, i.e. 0.0012 per component - so the gradient criterion can never be met and the run does not terminate: it keeps taking ~20 s line-search steps that change nothing, up to `max_iterations` (5000). |

The recipe that works, and the numbers it produced:

```bash
curcuma -opt polymer_2x.xyz -method gfnff -threads 36 -optimizer lbfgs \
        -opt.max_energy_rise 20000 -opt.convergence_count 3
# convergence_count 3 = energy + RMSD, dropping the unreachable gradient-norm criterion
```

Restarting the optimiser from its own output (fresh LBFGS memory) is what actually makes
progress - each round stalls in a line search well before a minimum:

| round | E [Eh] | steps | grad norm |
|---|---|---|---|
| start | -901.94 | - | 2.35 |
| 1 | -913.56 | 27 | 0.79 |
| 2 | -914.20 | 15 | 0.33 |
| 3 | -916.91 | 48 | 0.31 |
| 4 | -916.93 | 5 | 0.17 |
| 5 | -917.34 | 14 | 0.27 |
| 6 | -917.337 | 5 | 0.175 |

Both defaults should arguably scale with system size; that is a code change, not a
documentation matter, and is not made here.

## GFN2

Starting GFN2 from the GFN-FF-optimised structure is worth **14.65 Eh (38500 kJ/mol)**: GFN2
gives -11784.88 Eh on the raw structure and -11799.53 Eh on the GFN-FF-optimised one, before a
single GFN2 step. At ~160 s per step (4x A4500) a GFN2 optimisation of this system is a
multi-hour job, so pre-optimising with GFN-FF is not a convenience, it is the only affordable
route.

## What was not tested

- Longer than 200 fs, and no production trajectory - these runs answer "does it stay at the
  setpoint", not "is the sampling correct".
- Thermostats other than CSVR, and `-dt 0.75`, so the stability limit is only bracketed
  between 0.5 and 1.0 fs.
- `-hydrogen_mass` was used at 4 amu only, and its effect on dynamics (it changes the
  time scale of X-H motion by construction) was not examined.
- GFN2 MD at this size.
