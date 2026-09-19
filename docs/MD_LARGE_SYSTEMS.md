# Stable MD on a large system (polymer_2x, 7320 atoms)

> 🤖 AI-generated, machine-tested. Measurements from Sep 18, 2026 on 36 cores / 4x RTX A4500,
> commit 05d40030. Human production testing pending.

`docs/MULTI_GPU.md` recorded that an MD of `test_cases/molecules/larger/polymer_2x.xyz` heats
from 298 K to 20476 K within 40 fs, and attributed it to the structure not being pre-optimised.
That attribution was wrong. This page is what the measurement says instead.

## The settings that work

| Setting | `<T>` over 100 fs, target 300 K |
|---|---:|
| `-dt 1.0` (the default) | **18527 K** - blows up |
| `-dt 0.5` | 288.7 K |
| **`-dt 0.25`** | **296.1 K** |
| `-dt 1.0 -hydrogen_mass 4` | 273.1 K |

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
feature (`md_adaptive_step` checks exactly that).

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

### What it does not do

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
