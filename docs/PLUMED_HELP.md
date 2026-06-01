# Metadynamics with PLUMED

Enhanced-sampling MD via the PLUMED2 plugin.
Requires building with `-DUSE_Plumed=ON` (off by default).

## CLI flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-mtd` | Bool | `false` | Enable PLUMED metadynamics |
| `-plumed <file>` | String | `plumed.dat` | Path to PLUMED input file (alias: `plumed_file`) |
| `-mtd_dT <int>` | Int | `-1` | Temperature threshold for MTD start (K). Negative = active from step 0 |

## Quick start

```sh
curcuma -md input.xyz -method gfnff -mtd                        # uses default plumed.dat
curcuma -md input.xyz -method gfnff -mtd -plumed plumed.dat     # explicit file
```

A `plumed.dat` file must exist in the working directory (or at the path given by `-plumed`).

Minimal example `plumed.dat` for well-tempered metadynamics on a distance CV:

```
DISTANCE ATOMS=1,2 LABEL=d1
METAD ...
  ARG=d1
  PACE=500 HEIGHT=1.2 SIGMA=0.1
  BIASFACTOR=10 TEMP=298.15
  FILE=HILLS
... METAD
PRINT ARG=d1 STRIDE=100 FILE=COLVAR
```

## How it works

Curcuma embeds PLUMED2 via its C wrapper API. Each MD step:

1. Positions, energy, forces, masses, and virial are passed to PLUMED.
2. PLUMED computes collective variables and bias potentials.
3. PLUMED modifies the force array **in place** — the integrator receives biased forces for the next half-step update.
4. PLUMED writes its own log and output files.

### Unit conversions (Curcuma -> PLUMED)

| Quantity | Curcuma internal | Passed to PLUMED as | PLUMED unit | Factor |
|----------|-----------------|--------------------|-------------|--------|
| Energy   | Hartree          | Hartree            | kJ/mol      | 2625.5 |
| Length   | Angstrom         | **Bohr**           | nm          | 0.0529 |
| Time     | fs              | fs                 | ps          | 1e-3   |
| Mass     | amu             | amu                | amu         | 1      |
| Charge   | e               | e                  | e           | 1      |
| k_B*T    | Hartree          | Hartree            | kJ/mol      | 2625.5 |

Positions are converted from Angstrom to Bohr before passing to PLUMED, so that the
unit system (Bohr, Hartree, Hartree/Bohr) is consistent with the force array.
For periodic systems, box vectors are also converted to Bohr and passed via `setBox`.

These conversions are handled internally; no user action required.

## Thermal equilibration gate

If `-mtd_dT N` (N >= 0) is set, PLUMED calculations are deferred until the average temperature is within N K of the target and at least 10 steps have elapsed. This prevents the bias from destabilising the system during heating.

Default (`-mtd_dT -1`): PLUMED is active from step 0.

## Output files

| File | Written by | Description |
|------|-----------|-------------|
| `plumed_log.out` | Curcuma (via PLUMED `setLogFile`) | PLUMED master log: initialisation info, CV definitions, step summaries |
| `HILLS` | PLUMED (METAD action) | Gaussian hills deposited during metadynamics. Needed for free-energy reconstruction |
| `COLVAR` | PLUMED (PRINT action) | Collective variable values at each STRIDE. Configured in `plumed.dat` |
| Additional files | PLUMED | Depend on `plumed.dat` content: e.g. `fes.dat` (free-energy surface from REWEIGHT_METAD), histogram files, etc. |

File names and write frequency for `HILLS`, `COLVAR`, and all additional outputs are controlled by the `plumed.dat` configuration, not by Curcuma.

## Available PLUMED features (via plumed.dat)

Since Curcuma uses the full PLUMED wrapper, any PLUMED2 collective variable or bias can be defined in the input file:

- **Metadynamics**: standard and well-tempered (METAD)
- **Umbrella sampling**: RESTRAINT, LOWER_WALLS, UPPER_WALLS
- **Steered MD**: MOVINGRESTRAINT
- **Geometric CVs**: distances, angles, dihedrals, torsions
- **Coordination numbers**: COORDINATION
- **Path collective variables**: PATHCV, PATHMSD
- **Multicolvar functions**: GYROVATION, INERTIA, etc.
- **Analysis**: HISTOGRAM, REWEIGHT_METAD
- **Free-energy reconstruction**: sum_hills, reweighting

See <https://www.plumed.org/doc-v2.9/> for the full PLUMED manual.

## Internal RMSD-MTD (alternative)

Curcuma also provides built-in RMSD-based metadynamics (`-rmsd_mtd`) that does not require PLUMED. It biases the simulation away from reference structures using Gaussian hills on the RMSD coordinate. Both systems can coexist but operate independently.

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-rmsd_mtd` | Bool | `false` | Enable internal RMSD-based metadynamics |
| `-rmsd_mtd_k` | Double | `0.1` | Force constant for RMSD bias |
| `-rmsd_mtd_alpha` | Double | `10.0` | Gaussian width for RMSD bias |
| `-rmsd_mtd_pace` | Int | `1` | Bias deposition frequency (steps) |
| `-rmsd_mtd_max_gaussians` | Int | `-1` | Max stored bias structures (-1 = unlimited) |
| `-rmsd_mtd_ref_file` | String | `none` | Reference structures file |
| `-rmsd_mtd_atoms` | String | `-1` | Atom indices for RMSD calculation |
| `-rmsd_mtd_dt` | Double | `1e6` | Bias deposition time |

## Integration details (for developers)

- **Source**: `src/capabilities/simplemd.cpp` lines ~1656-1689 (init), ~2105-2124 (per-step), ~1962-1965 (finalize)
- **Header**: `src/capabilities/simplemd.h` — `#ifdef USE_Plumed` guards
- **CMake**: `CMakeLists.txt` lines 31-32 (option), 152-156 (FetchContent), 828-865 (build + link)
- **API**: PLUMED C wrapper (`plumed_create`, `plumed_cmd`, `plumed_finalize`)
- **Force modification**: `plumed_cmd("performCalc")` modifies `m_eigen_gradient` in place
- **Cleanup**: `plumed_finalize` is called on normal exit, CheckStop, and instability abort