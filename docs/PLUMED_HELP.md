# Metadynamics with PLUMED

Enhanced-sampling MD via the PLUMED2 plugin.
Requires building with `-DUSE_Plumed=ON` (off by default).

> **All PLUMED biases before the fix of Oct 30, 2026 acted with the wrong sign.** PLUMED adds its
> force to the array passed to `setForces`, curcuma passed the gradient, so a `RESTRAINT` repelled
> and `METAD` hills attracted. Results of earlier PLUMED runs with a bias are invalid; plain
> `PRINT` of a CV was unaffected. PLUMED also received kB instead of kB*T, i.e. T = 1 K for every
> action without its own `TEMP`.

**Status of the fix** (AI-generated, machine-tested; human production testing pending).
Tested: test `plumed_bias` (Ar2 `RESTRAINT`, period 0.3971 ps = analytic 0.3971 ps; fails on the
old code with 4e18 nm); ethane torsion, GFN-FF, 300 K: `METAD` 10 ps, 96-112 basin changes (old code
0, torsion pinned to 0.81-1.21 rad); well-tempered `METAD` without `TEMP`, 4 x 1 ns, barrier
2.66 kcal/mol (run range 2.52-2.83) against 2.43 from the histogram of 8 x 1 ns free MD
(run range 2.25-2.57), i.e. about 10 % high with non-converged hills.
Not tested: other systems and CVs, QM methods, restarts, `setVirial`/NPT, GPU together with PLUMED.

## CLI flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-mtd` | Bool | `false` | Enable PLUMED metadynamics |
| `-plumed <file>` | String | `plumed.dat` | Path to PLUMED input file (alias: `plumed_file`) |
| `-mtd_dT <int>` | Int | `-1` | Temperature threshold for MTD start (K). Negative = active from step 0 |
| `-no_plumed_redirect` | Bool | `false` | Disable automatic routing of PLUMED output files into BMT/plumed-files |

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
| Length   | Angstrom         | Angstrom            | nm          | 0.1    |
| Time     | fs              | fs                 | ps          | 1e-3   |
| Mass     | amu             | amu                | amu         | 1      |
| Charge   | e               | e                  | e           | 1      |
| k_B*T    | Hartree          | Hartree            | kJ/mol      | 2625.5 |

`setPositions`/`setForces` receive curcuma's `m_eigen_geometry`/`m_eigen_gradient` arrays
unmodified (Angstrom-based) — there is no in-place Bohr conversion step. PLUMED is told how
to interpret them via a single scalar passed to `setMDLengthUnits` (Angstrom -> nm = 0.1);
the same factor implicitly relates the force units (energyUnits / lengthUnits) since PLUMED
derives those from the two scalars above.

These conversions are handled internally; no user action required.

## Thermal equilibration gate

If `-mtd_dT N` (N >= 0) is set, PLUMED calculations are deferred until the average temperature is within N K of the target and at least 10 steps have elapsed. This prevents the bias from destabilising the system during heating.

Default (`-mtd_dT -1`): PLUMED is active from step 0.

## Output files

When BMT output is active (the default), Curcuma automatically routes all PLUMED output files into a `Basename.plumed_files/` subdirectory inside the BMT directory. This includes files written by PLUMED itself (COLVAR, HILLS, etc.) via automatic `FILE=` path rewriting, and the PLUMED log file.

### BMT output routing (default)

When BMT is active, Curcuma:
1. Creates `Basename.plumed_files/` inside the BMT directory
2. Rewrites `FILE=...` directives in `plumed.dat` to point into this subdirectory
3. Passes the modified content to PLUMED via `readInputLines` (the original `plumed.dat` on disk is **not modified**)
4. Saves the effective `plumed.dat` (with rewritten paths) as `Basename.plumed.dat` in the BMT root for provenance
5. Routes `plumed_log.out` into the subdirectory

Example directory layout:
```
ethanol.md.20260706_143000/        # BMT directory
  ethanol.plumed_files/             # PLUMED output subdirectory
    COLVAR                          # collective variable values
    HILLS                           # Gaussian hills
    plumed_log.out                  # PLUMED log
  ethanol.plumed.dat               # effective plumed.dat (with rewritten paths)
  ethanol.trj.xyz                   # trajectory
  ...
```

### Opt-out: `-no_plumed_redirect`

To disable automatic PLUMED output routing (all PLUMED files stay in CWD):
```sh
curcuma -md input.xyz -method gfnff -mtd -no_plumed_redirect
```

With `-no_plumed_redirect`:
- `plumed_log.out` is **still** routed into BMT (Curcuma controls this directly)
- `COLVAR`, `HILLS`, and other PLUMED output files stay in CWD (no `FILE=` rewriting)
- No provenance `plumed.dat` is saved to BMT

### BMT disabled (`-no_bmt`)

When BMT is disabled entirely, PLUMED output stays in CWD with no subdirectory — the original behavior.

| Scenario | PLUMED log | COLVAR / HILLS / etc. | Provenance plumed.dat |
|----------|-----------|----------------------|---------------------|
| BMT active (default) | `BMT/Basename.plumed_files/` | `BMT/Basename.plumed_files/` | `BMT/Basename.plumed.dat` |
| BMT + `-no_plumed_redirect` | `BMT/Basename.plumed_files/` | CWD (no rewriting) | Not written |
| `-no_bmt` | CWD | CWD | Not written |

### Path rewriting details

Curcuma reads the `plumed.dat` file and prepends the `Basename.plumed_files/` path to each `FILE=` value that:
- Is a bare filename (no `/` in the path)
- Is not an absolute path (does not start with `/`)

Examples:
- `FILE=COLVAR` becomes `FILE=ethanol.plumed_files/COLVAR` (BMT basename = "ethanol")
- `FILE="HILLS"` becomes `FILE="ethanol.plumed_files/HILLS"`
- `FILE=/tmp/test_hills` stays unchanged (absolute path)
- `FILE=subdir/output.dat` stays unchanged (already has a path separator)

### Output file reference

| File | Written by | Description |
|------|-----------|-------------|
| `plumed_log.out` | Curcuma (via PLUMED `setLogFile`) | PLUMED master log: initialisation info, CV definitions, step summaries |
| `HILLS` | PLUMED (METAD action) | Gaussian hills deposited during metadynamics. Needed for free-energy reconstruction |
| `COLVAR` | PLUMED (PRINT action) | Collective variable values at each STRIDE. Configured in `plumed.dat` |
| Additional files | PLUMED | Depend on `plumed.dat` content: e.g. `fes.dat` (free-energy surface from REWEIGHT_METAD), histogram files, etc. |

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
| `-rmsd_econv` | Double | `1e8` | Energy convergence threshold for bias structure addition. Lower values add new reference structures more frequently (more aggressive sampling); higher values make structure addition rarer (more conservative). Default effectively disables the filter. |

### How RMSD-MTD works

RMSD-based metadynamics adds a history-dependent repulsive bias potential during MD, pushing the system away from already-visited conformations. The collective variable (CV) is the RMSD from each stored reference structure.

#### Algorithm per MD step (every `pace` steps)

1. **Geometry extraction** — The current atomic positions of the selected atoms (`-rmsd_mtd_atoms`) are copied into a local geometry.

2. **First-step initialisation** — If no reference structure exists yet, the current geometry is stored as the first reference and written to `Basename.mtd.xyz`.

3. **Bias potential evaluation** (via `BiasThread::execute()`):
   For each stored reference structure *i* (index `i` in `m_biased_structures`):

   a. **RMSD calculation**: Compute `RMSD_i` between current geometry and reference *i* via Kabsch alignment.

   b. **Gaussian hill**: Evaluate the Gaussian contribution:
   ```
   G_i = exp(-RMSD_i^2 * alpha) * dT
   ```
   where `alpha` = `-rmsd_mtd_alpha` (width) and `dT` = `-rmsd_mtd_dt` (well-tempered scaling factor, default 1e6 = effectively standard MTD).

   c. **Well-tempered rescaling** (if `-wtmtd true`):
   ```
   factor_i += exp(-E_accumulated_i / (k_B * T * DT))
   ```
   Standard MTD uses `factor_i = counter_i` instead (each deposition increments the counter).

   d. **Convergence filter** (`rmsd_econv`):
   ```
   if (G_i * rmsd_econv > N_stored):
       counter_i  += 1
       E_accumulated_i += G_i
   ```
   With the default `rmsd_econv = 1e8`, this condition is always satisfied — every reference structure increments its counter at every step. Lower values make the filter more selective: only references where the current RMSD is small enough (i.e. the system is still near that reference) accumulate counts and energy.

   e. **Bias energy accumulation**:
   ```
   V_bias += G_i * factor_i * k
   ```
   `k` = `-rmsd_mtd_k` (force constant, default 0.1).

   f. **Bias gradient**:
   ```
   dV/dR_i = -2 * alpha * k / N_atoms * G_i * factor_i * dT * (dRMSD/dR)
   ```
   The RMSD gradient `dRMSD/dR` comes from the Kabsch alignment driver. Each reference contributes additively to the total bias gradient.

4. **Force modification** — The total bias gradient is added to `m_eigen_gradient` (in-place), so the MD integrator receives biased forces for the next half-step update. Units are Hartree/Bohr throughout.

5. **COLVAR output** — Per-reference and global COLVAR files are appended:
   - Global: `step  RMSD_0  fragment_distances...  V_bias_total`
   - Per-reference: `step  RMSD_i  V_bias_i  counter_i  factor_i`

6. **New reference structure decision**:
   ```
   if (V_bias_total * rmsd_econv < N_stored  AND  rmsd_fix_structure == false):
       add current geometry as new reference structure
       append to Basename.mtd.xyz
   ```
   This is the inverse of the convergence filter in step 3d: when the total bias energy is small (the system is far from all references), a new reference is added. With the default `rmsd_econv = 1e8`, new structures are almost never added (the condition `1e8 * V_bias < N` is virtually never true). Lower `rmsd_econv` values trigger new references more often.

7. **Restart persistence** — On checkpoint, all reference structures (geometry, RMSD, accumulated energy, counter, factor) are saved to `curcuma_restart.json` and restored on restart.

#### Two roles of `rmsd_econv`

| Role | Condition | Effect of low `rmsd_econv` | Effect of high `rmsd_econv` (default 1e8) |
|------|-----------|--------------------------|----------------------------------------|
| Counter update filter (step 3d) | `G_i * econv > N_stored` | Selective: only near references increment | Always true: every reference increments every step |
| New structure gate (step 6) | `V_bias * econv < N_stored` | Frequent: new references added often | Virtually never: fixed set of references from start |

#### Key parameters summary

| Parameter | Symbol | Role |
|-----------|--------|------|
| `-rmsd_mtd_k` | k | Overall bias strength. Higher = stronger repulsion from visited regions |
| `-rmsd_mtd_alpha` | alpha | Gaussian width. Smaller = narrower hills = finer resolution, wider = broader repulsion |
| `-rmsd_mtd_pace` | pace | How often (in steps) the bias is evaluated and applied |
| `-rmsd_mtd_dt` | DT | Well-tempered bias factor. Only used with `-wtmtd true`. Controls adaptive rescaling |
| `-rmsd_econv` | econv | Convergence threshold. Dual role: filters counter updates and gates new reference addition |
| `-wtmtd` | — | Enable well-tempered mode: adaptive Gaussian height via Boltzmann rescaling |
| `-rmsd_mtd_ref_file` | — | File with initial reference structures (multi-structure XYZ) |
| `-rmsd_mtd_atoms` | — | Which atoms to include in RMSD calculation (default: all) |

## Integration details (for developers)

- **RMSD-MTD source**: `src/capabilities/simplemd.cpp` — `BiasThread::execute()` (bias evaluation, lines ~82-142), `ApplyRMSDMTD()` (main loop integration, lines ~2516-2615), init (~767-822), restart I/O (~1247-1274, ~1500-1520)
- **PLUMED source**: `src/capabilities/simplemd.cpp` lines ~2129-2149 (per-step), ~1987-1989 (finalize)
- **Header**: `src/capabilities/simplemd.h` — `BiasThread` class (lines 47-131), `BiasStructure` struct (lines 47-55)
- **CMake**: `CMakeLists.txt` lines 31-32 (option), 152-156 (FetchContent), 828-865 (build + link)
- **API**: PLUMED C wrapper (`plumed_create`, `plumed_cmd`, `plumed_finalize`)
- **Force modification**: `plumed_cmd("performCalc")` modifies `m_eigen_gradient` in place
- **Cleanup**: `plumed_finalize` is called on normal exit, CheckStop, and instability abort