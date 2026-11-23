# NPT: isotropic barostat for the periodic GFN-FF (CPU)

Status: 🤖 AI-generated, ⚙️ machine-tested (Nov 2026, WP-PERIODIC-NONBONDED N2); no human production testing.
Work-package record with every measurement: [WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md) section N2.

## What it does

`-md` with a unit cell (`unitcell a b c [alpha beta gamma]` in the XYZ comment line; angles 90 when omitted) and the
periodic native CPU GFN-FF energy (`-method gfnff`, `coulomb_ewald` default) can couple the volume to a target
pressure (since Nov 23, 2026 also the GFN-FF CUDA hybrid, `-gpu cuda -gfnff.gpu_periodic true`; WP-PERIODIC-NONBONDED
A1-GPU). Every `pcoupl_interval` steps the cell edges and **every atom position** are scaled by the same factor mu
(atomic scaling, consistent with the atomic virial of N1); the pressure used is the full instantaneous pressure of the
step just finished, `P = (2 Ekin + tr Xi) / (3 V)`. Isotropic only; cubic and orthorhombic cells (scaled uniformly).

| parameter (module `simplemd`) | default | meaning |
|---|---|---|
| `barostat` | `none` | `berendsen` or `crescale` |
| `pressure` | 1.01325 bar | target pressure P0 |
| `tau_p` | 1000 fs | coupling time |
| `compressibility` | 4.5e-5 1/bar | isothermal compressibility beta (water at 300 K); sets the coupling strength, not the result |
| `pcoupl_interval` | 10 | the barostat acts every N steps with the time step N dt |
| `gfnff.pressure_lrc` | true | add the long-range correction of the D4 pairs crossing the evaluation radius to the pressure (below) |
| `gfnff.pressure_hb_crossing` | true | add the crossing term of the HB/XB search radii to the pressure (below; false: pressure before Nov 23, 2026) |
| `gfnff.pressure_crossing_shell` | 0.03 | shell width (relative to the radius) of its estimator |

A barostat switches `compute_pressure` on: the status line gains `P`, `<P>` (bar), `V` (A^3), `rho` (g/cm^3) and `L`
(shortest cell edge, A); the restart files carry `cell` (9 numbers, A), the barostat parameters, the step counter of
the coupling and `<P>`; every trajectory frame carries `unitcell a b c 90.0 90.0 90.0` at the end of its comment line.

```bash
# 2 ps NPT of a water box after a CSVR heat-up (restart file edited: MaxTime)
curcuma -md box432.xyz -method gfnff -threads 4 -restart_file heat.final.json \
        -barostat berendsen -pressure 1.01325 -tau_p 1000 -compressibility 4.5e-5 -pcoupl_interval 10
curcuma -md box432.xyz -method gfnff -thermostat csvr -coupling 100 -T 298.15 -dt 0.5 -max_time 20000 \
        -barostat crescale -tau_p 1000 -seed 7
```

## Pressure: long-range correction of the D4 tail (V-LRC, Nov 9, 2026)

The virial `Xi` of N1 is the strain derivative of the model energy at a fixed pair set; for the homogeneous D4 tail that
gives `tr Xi_tail = 3 E_tail`. The pressure of the untruncated pair potential (g(r) = 1 beyond r_c) needs
`6 E6_tail + 8 E8_tail`; the difference `3 E6_tail + 5 E8_tail` is the impulsive term of the pairs crossing the sharp
radius r_c (classical long-range pressure correction). With `pressure_lrc true` (default) `P`, `<P>` and the barostat
use it; `pressure_lrc false` gives the N1/N2 pressure bitwise. Energies and forces are not affected. Water at the W1
matrix radii: -70.6 / -34.3 / -16.8 bar (432 / 864 / 1728 waters; 3456: about -8 bar from the P1 fit). Derivation and
tests: [WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md) section V-LRC.

**Runs made before the correction** (binaries up to `0351a245`, W1 project: the run matrix, P1, P7, P8, P2, P10) used
the old estimator: their pressure is too high by the amounts above, and their NPT densities are too low by about
0.32 / 0.15 / 0.07 / 0.04 % (432 / 864 / 1728 / 3456 waters; linear estimate with kappa_T = 4.5e-5 /bar from the cutoff
scan, not an NPT rerun). Droplet runs without a cell are not affected.

## Pressure: crossing term of the HB/XB search radii (HB-X, Nov 23, 2026)

A3-P measured the term of the H-bond triples that cross the sharp search radii under a volume change (-31 to -38 bar
on water, missing from the fixed-list virial). With `pressure_hb_crossing true` (default) `P`, `<P>` and the barostat
contain a statistical estimate of it (virial components `hb_crossing`, `xb_crossing`: triples in a shell inside the
radius, energy density extrapolated linearly to the radius); `false` gives the earlier pressure bitwise. Energies and
forces are not affected; CPU and the CUDA hybrid. Water, n = 15 frames: estimate -30.6 bar, residual against the
re-derived-list strain slope -0.15 +- 0.79 bar. **Behaviour change**: the reported P of a water box is about 31 bar
lower than with earlier binaries; NPT runs made before (all of Nov 2026) have densities about 0.2 % too low (rough
linear estimate, kappa_T ~6e-5 /bar, no NPT rerun). Derivation and tests: [WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md)
section HB-X. Still not in P: repulsion and CN crossing parts (below 1.5 bar), the |E_HB| candidate filter (about 0.1 bar).

## Equations (`SimpleMD::applyBarostat`, the code comment carries the derivation)

- Berendsen (Berendsen et al., J. Chem. Phys. 81, 3684 (1984)): `mu^3 = 1 - (N dt / tau_p) beta (P0 - P)`, velocities
  unchanged. Weak coupling: the volume relaxes with `tau_p`, its fluctuations are not those of an ensemble.
- Stochastic cell rescaling, C-rescale (Bernetti and Bussi, J. Chem. Phys. 153, 114107 (2020); the isotropic form
  below is from the work order and from memory, **not checked against the paper's text** -- the one web search
  allowed for the work did not return the equations): `d ln V = -(beta/tau_p)(P0 - P) dt + sqrt(2 kB T beta /
  (V tau_p)) dW`, one Euler-Maruyama step per coupling (dW = sqrt(N dt) xi, xi from the seeded MD generator, drawn after
  the thermostat's numbers of the same step; T = thermostat target), positions scaled by `exp(d ln V / 3)`,
  velocities by `exp(-d ln V / 3)`. Stationary density `exp(-(H + P0 V)/kB T)` in (r, p, V); the ideal-gas limit
  `<V> = N_eff kB T / P0`, relative variance `1/N_eff`, N_eff = dof/3 + 1 (= N + 1 with 3N momentum components; = N
  with curcuma's periodic DOF count 3N - 3) is checked by ctest `pbc_npt_idealgas_*`.
- Placement: kick - drift - force(x_{n+1}) - kick - thermostat - **barostat**. The forces kept for the first half-kick of
  the next step are those of the unscaled x_{n+1} (as in GROMACS; |mu - 1| ~ 1e-4 per coupling at the defaults); no
  extra energy evaluation, and without a barostat the step is bitwise the N1 step.

## Cell update (`EnergyCalculator::setUnitCell` -> `GFNFF::updateCell`)

Between two energy calls the instance changes its cell and keeps everything consistent: cell in A and Bohr, the
workspace cell, the Ewald sums of the Coulomb energy and of the EEQ matrix (new lattice and k vectors at the alpha,
r_cut and k_max fixed at the set-up; n_max per axis follows L), the shared Ewald geometry data (invalidated), the D4
generator cell, the D4 tail (reads V per call), and every pair list and HB/XB detection (forced rebuild at the next
call: the Verlet-skin argument does not cover a changed lattice). Topology, its charges and the EEQ warm start are
kept. Every active cutoff is checked against the new L_min/2 first; a violated check leaves the instance unchanged.

## Frozen cutoffs and the limit abort

The virial is the strain derivative at fixed radii, so when the barostat starts the evaluation radii are frozen at
absolute values (`GFNFF::freezeCutoffsForCellChanges`): the D4 evaluation radius at
`min(current, L_min/2 (1.15)^(-1/3))` (the density may rise by 15 % before it limits; the homogeneous tail covers the
part beyond it), the Ewald r_cut, repulsion (20 Bohr), HB/XB search radii (21.21 Bohr with `hb_accuracy` 0.1) and the
XB radius are absolute already. For the 432-water boxes of the W1 project (L0 = 44.39 Bohr) the D4 radius becomes
21.19 Bohr and the HB/XB search radius stays the limiter: the run stops when `L_min` falls below 42.43 Bohr
(density +12.9 %), with the message `NPT: cell shrank below the cutoff limit L=... Bohr (needs > ... Bohr, <cutoff>)
at step ...; restart written: <basename>.npt_limit.json` holding the last valid state (`StopReason::CellLimit`, exit
after `finalizeRun`). Larger boxes with the same water density are limited the same way (same absolute radii), i.e.
by a density gain of 12.9 % at 432 waters and more at larger boxes.

## Refused

RATTLE (constraint virial not implemented), `cleanenergy`, `-gpu` except the GFN-FF CUDA hybrid (`-gpu cuda
-gfnff.gpu_periodic true`, accepted since Nov 23, 2026, WP-PERIODIC-NONBONDED A1-GPU) and native gfn1/gfn2 on CUDA, wall potentials (they do not scale), a non-periodic energy or no cell, `frag_charge_model ensemble`
with a net charge, an implicit solvent, `dispersion_atm`. `crescale` without a thermostat runs with a warning (NVE +
barostat is not an ensemble).

## Not tested

Real production NPT of a liquid (only 2 ps on the 432-water box, started from a packing: the pressure of such a box
is not equilibrated), anisotropic cells under strain (uniform scaling only), the reactive topology mode with a
barostat (the scan is forced at every cell change, not run), charged boxes, restart files with RATTLE, more than 1296
atoms with a barostat. Human production testing pending.
