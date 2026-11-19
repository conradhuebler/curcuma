# Wall Radius Schedule and Wall Pressure (SimpleMD)

> Status: AI-generated (Nov 2026), machine-tested. Human production testing pending.
> The pressure column is a wall force measure, not a thermodynamic pressure.

A spherical SimpleMD wall (`-wall_type spheric`, `-wall_potential harmonic|logfermi`) can change its
radius during the run. Compressing the wall raises the force the wall exerts on the atoms, which is
reported per step as force per wall area. There is no barostat; the wall is the only pressure handle.

## Parameter

| Parameter | Type | Meaning |
|---|---|---|
| `wall_radius_schedule` | string | `radius:steps:n[;radius:steps:n ...]`, radius in Angstrom |

- Grammar and parser are those of `temp_schedule` (`SimpleMD::ParseSchedule`, `StepRamp`): linear from the
  current radius to `radius` over `n` integration steps, then the next segment; the last radius is held.
- Start radius is `wall_radius`; if it is 0 (default) it is auto-sized from the structure (8.08 Angstrom for the 3 CO + 3 H2O gas phase example) and the first segment starts from that value, so give `-wall_radius` to start at a chosen radius.
- Timing: the radius of step k is set at the top of step k, the same convention as `temp_schedule`; the last CSV row repeats the previous wall values (the final step is printed without a new wall evaluation, existing behaviour).
- Only `steps` segments and radius > 0 are accepted; `reach` or a malformed string disables the schedule with a warning.
- Needs `wall_type spheric`; with another wall type the schedule is disabled with a warning. The `rect` wall has no schedule (not implemented).
- Also acts on `wall_potential pbc` (radius of the container), but no pressure columns are written there.

```bash
curcuma -md gas.xyz -method gfn2 -wall_type spheric -wall_potential logfermi -wall_temp 4000 \
        -wall_radius_schedule "6.0:steps:1000;3.0:steps:3000"
```

## Live control (API)

- `SimpleMD::setWallRadius(double R)`: thread-safe request (`std::atomic`), applied at the top of the next step; cancels the schedule for the rest of the run, like `setTargetTemperature()`.
- `wallRadius()`, `wallForce()` [Eh/Angstrom], `wallPressure()` [GPa], `wallPressureAverage()` [GPa, last 100 steps].
- `setWallTemp()` / `setWallBeta()` change the wall strength live (existing).

## Wall pressure

Per step, in the wall routine: `F = sum_i |F_wall,i|` (magnitude of the radial wall force of every atom,
Eh/Angstrom, the unit of the wall gradient), `P = F / (4 pi R^2)` with the radius of that step,
`1 Eh/Angstrom^3 = 4359.7447 GPa`.

CSV file `<basename>.md.csv` ([MD_CSV_OUTPUT.md](MD_CSV_OUTPUT.md)): three appendix columns after the existing ones
(after `dipole_debye`/`n_unique` if present), written only with a spherical harmonic or logfermi wall, so files of other runs are unchanged:

| Column | Unit | Meaning |
|---|---|---|
| `wall_radius_A` | Angstrom | wall radius of the step |
| `wall_pressure_GPa` | GPa | instantaneous F/(4 pi R^2) |
| `wall_pressure_avg100_GPa` | GPa | mean of the last 100 steps (fewer at the start) |

Limits:

- It is the force of a confining potential per area, not an equilibrium pressure. The value depends on `wall_temp`, `wall_beta` (the wall stiffness), on whether the system is in a quasi-stationary state (a fast radius ramp gives a non-equilibrium transient), and it includes the soft logfermi force of atoms inside the radius.
- Magnitudes are summed; the force of an atom that the wall pushes inwards counts positive. For a gas of non-interacting particles in a hard wall this equals the kinetic pressure; for a soft wall it is only a proxy.
- Console output is unchanged; the columns exist only in the CSV.
- A harmonic wall evaluates `x/|r|`; an atom exactly at the origin gives NaN (existing behaviour, not changed).

## Tests

Machine-tested (`ctest -R cli_simplemd_17`, test directory `test_cases/cli/simplemd/17_wall_radius_schedule`):

- Single He atom outside a harmonic wall: reported force equals `k (r - R)` with `k = wall_temp * kB` (Eh/Angstrom^2), pressure equals `F/(4 pi R^2)` in GPa.
- Radius follows the schedule linearly in the CSV; the pressure of a compressed run rises.
- Without `wall_radius_schedule` the energies of a 200-step GFN2 run were bit-identical to `release/curcuma` (same seed, same compiler flags, columns step to virial_correction_avg of `<basename>.md.csv`, GFN-FF and GFN2, 1 and 4 threads, md5 equal). A build without AVX/march flags differs from `release/curcuma` at the 1e-8 relative level in the temperature, independent of this change.
- End to end: 3 CO + 3 H2O gas phase, GFN2, 2500 K, electronic 1000 K, logfermi wall 4000 K, `6.0:steps:1000;3.0:steps:3000`, 4000 steps, 4 threads, 14 s: radius 8.08 to 3.0 Angstrom, running-mean pressure about 0 to 0.3 GPa above 6 Angstrom and about 3 to 5 GPa at 3.0 to 3.4 Angstrom (noisy, single run, seed 1), no SCF abort.

Not tested: CUDA/ROCm run, GPU-resident MD, `rect` wall, `pbc` container with a schedule, restart files (the radius is not stored in the restart state), qurcuma, periodic systems, comparison of the wall pressure with a virial or equation-of-state pressure.
