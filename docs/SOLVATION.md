# Solvation Models in Curcuma

Implicit solvation support for GFN and MNDO methods via TBLite and Ulysses interfaces.

## Overview

Curcuma supports implicit solvation models for quantum chemistry calculations:

- **TBLite Interface**: CPCM, GB (Generalized Born), ALPB for GFN methods. ⚠️ **Pending**:
  this path needs `USE_TBLITE=ON`; the reference dev build (`release/`) has TBLite off, so
  the TBLite Interface solvation models (CPCM in particular — not yet native, see below)
  are currently unverified in that build. Native ALPB/GBSA below are unaffected.
- **Ulysses Interface**: GBSA (Generalized Born) for GFN and MNDO methods
- **Native GFN1/GFN2**: ⚠️ self-consistent **ALPB** and **GBSA** (June 2026) —
  AI-implemented, machine-tested vs tblite (≤1e-8 Eh on the validation set), human
  production testing pending. CPCM not yet implemented natively. See below.
- **Native GFN-FF**: ⚠️ self-consistent **ALPB** (June 2026, WP5) — the Born reaction
  field couples into the EEQ charge solve (`A_eeq += B`), so the EEQ charges polarize in
  the solvent. Matches **xtb 6.7.1** (`--gfnff --alpb`) to ≤1e-8 Eh (7 mol × 4 solvents,
  `ctest -L gfnff_solvation`); analytic gradient FD-validated. GFN-FF has no separate GBSA
  model (the reference always uses ALPB, `gfnff_alpb.F90:23-24`), so `-gfnff.solvent_model
  gbsa` maps to ALPB. `-method gfnff -gfnff.solvent water -gfnff.solvent_model alpb`.
  Machine-tested only.

## Native GFN1/GFN2 ALPB / GBSA (no external dependency)

The native `gfn1`/`gfn2` backends include self-consistent ALPB and GBSA solvation
models matching the tblite parameterization (Born electrostatics + CDS surface
tension/H-bond + state shift, see "Reference states"). GFN1 uses CM5 charges, GFN2 Mulliken charges, as in
tblite. The solvent reaction field enters the SCF (the wavefunction polarizes), the
energy is added to the total, and the nuclear gradient is analytic. ALPB uses the
P16 Born kernel; GBSA is ALPB with the analytical-PB shape term off and the
classical Still kernel (exactly tblite's ALPB=11/12 vs GBSA=21/22 distinction).

```bash
# Dotted form required (the flat -solvent is ambiguous across providers):
curcuma -sp mol.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model alpb
curcuma -sp mol.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model gbsa
curcuma -opt mol.xyz -method gfn1 -xtb.solvent dmso                          # defaults to alpb
```

- `-xtb.solvent <name>`: water, dmso, acetone, chloroform, methanol, … (tblite set)
- `-xtb.solvent_model`: `alpb` (default when a solvent is given), `gbsa`, `cpcm` (not yet native), `none`.
  Legacy numeric codes (3=alpb, 2=gbsa, 1=cpcm, 0=none) are still accepted.
- ROCm (`-gpu rocm`): the reaction field v_at += B q_at enters the device SCF potential and the resident loop (`-xtb.rocm_device_solvation`, default true; energy and gradient on the host); see [GPU_TUNING.md](GPU_TUNING.md).

**Validation (machine-tested):** native total energy matches tblite to ≤1e-8 Eh for
7 molecules × {water, dmso, acetone, chloroform} × {gfn1, gfn2} — for **both** models
(ALPB: `ctest -L _solvation`, 56 tests; GBSA: `ctest -L _gbsa`, 56 tests); analytic
gradients FD-validated for both (`ctest -R xtb_solvation_numgrad`). **Not tested:**
other solvents, ions, metals, large systems, long MD/opt stability. See
[SQM_SOLVATION_WP.md](SQM_SOLVATION_WP.md).

## Supported Solvents

Native GFN1/GFN2 (ALPB and GBSA) and native GFN-FF (ALPB) share one name handling (`src/core/solvation/solvent_names.h`):
names are case-insensitive (`Water`, `CHCl3`), spaces and underscores are ignored, and the aliases of xtb 6.7.1 are
accepted. A name that has no parameter set for the chosen method and model is an **error** (exit code 1), never a silent
gas-phase run. `none`, `vacuum` and `gas` mean no solvent.

| Canonical name | Accepted spellings (besides the canonical one) |
|---|---|
| `water` | `h2o` |
| `ch2cl2` | `dichloromethane`, `dichlormethane`, `methylenechloride`, `dcm` |
| `chcl3` | `chloroform`, `trichloromethane`, `tcm` |
| `acetonitrile` | `methylcyanide` |
| `cs2` | `carbondisulfide`, `carbon disulfide` |
| `ether` | `diethylether`, `diethyl ether` |
| `dmf`, `dmso` | `dimethylformamide`, `dimethylsulfoxide` |
| `furane` | `furan` |
| `hexane` | `n-hexane`, `nhexane`, `n-hexan`, `nhexan` |
| `thf` | `tetrahydrofuran` |
| `ethylacetate` | `ethyl acetate` |
| `woctanol` | `octanol wet`, `wet octanol`, `wetoctanol` |

Also available under their own names: `acetone`, `aniline`, `benzaldehyde`, `benzene`, `dioxane`, `ethanol`, `hexadecane`, `methanol`,
`nitromethane`, `octanol`, `phenol`, `toluene` (25 ALPB solvents; the GBSA sets are smaller: GFN1 12, GFN2 14, as in tblite and
xtb, e.g. `dmf` and `hexane` have no GFN1 GBSA parameters). `dichloromethane`, `diethyl ether`, `carbon disulfide`, `ethyl acetate`
and `octanol wet` are curcuma additions. The alias list is that of the xtb source tree in `external/xtb` (`normalizeSolventName`,
newer than the 6.7.1 binary): the 6.7.1 binary accepts only `furane`, `h2o`, `n-hexane` and the canonical names (measured: it rejects `dcm`,
`dichloromethane`, `methylenechloride`, `furan`, `diethylether`, `tetrahydrofuran`, `trichloromethane`, `methylcyanide`), curcuma accepts the longer list.

### Reference states (`-solvent_state`)

One flat argument for GFN1, GFN2 and GFN-FF (`-solvent_state gsolv|bar1M|reference`, case-insensitive; the dotted
`-xtb.solvent_state` and `-gfnff.solvent_state` also work). The state shifts the solvation free energy as `xtb --alpb water [state]`:

| State | Meaning | Shift (298.15 K) |
|---|---|---|
| `gsolv` (default) | gas and solution at the same concentration | 0 |
| `bar1M` | 1 bar ideal gas to 1 mol/L solution | +3.0312e-3 Eh (+1.902 kcal/mol), all solvents |
| `reference` | 1 bar ideal gas to the pure solvent, as in xtb | solvent dependent, e.g. water +1.8138e-3 Eh |

One implementation serves all methods (`solutionStateShift`): ideal-gas molar volume 24.79 L/mol as in xtb and pprcht/gfnff,
CODATA-2018 kB. Against xtb 6.7.1 the shifts agree to 3e-9 Eh (`ctest -R cli_sqm_18`, 3 solvents x 3 methods). tblite derives
24.78960 L/mol from R instead of 24.79 and therefore differs by 1.4e-8 Eh in `bar1M`; it was not ported separately. The temperature is fixed at 298.15 K.

### Not implemented natively

CPCM (tblite `--cpcm`, `solvent_model cpcm` is an error), ddCOSMO (xtb `--cosmo`), CPCM-X (xtb `--cpcmx`), the ionic strength of xtb's GBSA,
a choice of the Born kernel (ALPB always P16, GBSA always Still), the Lebedev grid size (230) and a solvent given only by its dielectric
constant. GFN-FF with `solvent_model gbsa` runs ALPB with a warning, as pprcht/gfnff does; xtb 6.7.1 runs a Still-kernel GBSA with the GFN2 parameters there.

## Usage

### Basic Syntax

```bash
# Single point calculation in water
curcuma -sp molecule.xyz -method gfn2:tblite -solvent water

# Geometry optimization in DMSO
curcuma -opt molecule.xyz -method gfn2:tblite -solvent dmso

# Semi-empirical PM6 in acetone (requires Ulysses)
curcuma -sp molecule.xyz -method pm6:ulysses -solvent acetone
```

### Advanced Options

#### TBLite Solvation Models

```bash
# CPCM (Conductor-like Polarizable Continuum Model)
curcuma -sp molecule.xyz -method gfn2:tblite \\
        -solvent water -solvent_model 1

# GB (Generalized Born)
curcuma -sp molecule.xyz -method gfn2:tblite \\
        -solvent water -solvent_model 2

# ALPB (Analytical Linearized Poisson-Boltzmann)
curcuma -sp molecule.xyz -method gfn2:tblite \\
        -solvent water -solvent_model 3
```

#### Custom Dielectric Constant

```bash
# Use custom epsilon (overrides solvent name)
curcuma -sp molecule.xyz -method gfn2:tblite \\
        -solvent_epsilon 78.4  # Water's dielectric constant
```

## Method Availability

| Method | Provider | Solvation | Models Available |
|--------|----------|-----------|------------------|
| **GFN2** | TBLite | ✅ Yes | CPCM, GB, ALPB |
| **GFN2** | Ulysses | ✅ Yes | GBSA |
| **GFN2** | Native | ❌ Not yet | Planned |
| **GFN1** | TBLite | ✅ Yes | CPCM, GB, ALPB |
| **GFN1** | Ulysses | ✅ Yes | GBSA |
| **GFN1** | Native | ❌ Not yet | Planned |
| **PM6** | Ulysses | ✅ Yes | GBSA |
| **PM6** | Native | ❌ Not yet | Planned |
| **AM1** | Ulysses | ✅ Yes | GBSA |
| **PM3** | Ulysses | ✅ Yes | GBSA |
| **MNDO** | Ulysses | ✅ Yes | GBSA |

## Compilation Requirements

### Enable TBLite Support

```bash
cd build
cmake -DUSE_TBLITE=ON ..
make -j4
```

### Enable Ulysses Support

```bash
cd build
cmake -DUSE_ULYSSES=ON ..
make -j4
```

### Check Compilation Status

```bash
curcuma --version  # Shows enabled features
```

## Output Examples

### TBLite with Water Solvation (Verbosity 1)

```bash
$ curcuma -sp h2o.xyz -method gfn2:tblite -solvent water -verbosity 1

[INFO] Initializing TBLite quantum chemistry method
[PARAM] accuracy: 1
[PARAM] SCF_maxiter: 100
[PARAM] temperature: 300.0 K
[PARAM] solvation_model: GB (Generalized Born)
[PARAM] solvent: water
[INFO] TBLite implicit solvation enabled

[ENERGY] TBLite Energy: -5.123456 Eh
```

### Ulysses with DMSO Solvation (Verbosity 1)

```bash
$ curcuma -sp molecule.xyz -method pm6:ulysses -solvent dmso -verbosity 1

[INFO] Initializing Ulysses quantum chemistry method
[PARAM] method: PM6
[PARAM] solvation_model: GBSA (Generalized Born)
[PARAM] solvent: dmso
[INFO] Ulysses implicit solvation enabled

[ENERGY] Ulysses Energy: -42.789012 Eh
```

## Solvation Models Explained

### CPCM (Conductor-like PCM)
- Treats solvent as conductor
- Fast and robust
- Good for polar solvents
- Reference: Cossi et al., J. Comput. Chem. 24, 669 (2003)

### GB (Generalized Born)
- **Default model for TBLite** (`-solvent_model 2`)
- Pairwise atom-based approximation
- Fast analytical gradients
- Excellent for optimization
- Reference: Still et al., J. Am. Chem. Soc. 112, 6127 (1990)

### ALPB (Analytical Linearized Poisson-Boltzmann)
- Advanced GB variant with parameter fitting
- High accuracy for small molecules
- Reference: Ehlert et al., J. Chem. Theory Comput. 17, 4250 (2021)

### GBSA (Generalized Born + Surface Area)
- **Default model for Ulysses**
- GB model + non-polar cavity term
- Includes atomic surface areas
- Reference: Hawkins et al., J. Phys. Chem. 100, 19824 (1996)

## Energy Decomposition

With `-verbosity 2` or higher:

```bash
[PARAM] electronic_energy: -42.123 Eh
[PARAM] core_repulsion: 15.456 Eh
[PARAM] solvation_energy: -0.234 Eh  # Implicit solvation contribution
[ENERGY] Total Energy: -26.901 Eh
```

## Performance Considerations

- **GB/GBSA**: Fastest (analytical gradients)
- **CPCM**: Medium speed (iterative solver)
- **ALPB**: Similar to GB, highly accurate

For geometry optimization, **GB is recommended** for speed.

## Troubleshooting

### Solvation Not Applied

**Problem**: Solvent parameter ignored
**Solution**: Specify method provider explicitly:
```bash
# Correct (explicit provider)
curcuma -sp mol.xyz -method gfn2:tblite -solvent water

# May fall back to native (no solvation)
curcuma -sp mol.xyz -method gfn2 -solvent water
```

### Compilation Errors

**Problem**: TBLite/Ulysses not found
**Solution**: Check CMake configuration:
```bash
grep "USE_TBLITE\|USE_ULYSSES" build/CMakeCache.txt
```

## References

1. **TBLite Library**: Grimme et al., J. Chem. Phys. 143, 054107 (2015)
2. **Ulysses Package**: Menezes et al., ChemRxiv (2023)
3. **GBSA Model**: Grimme, J. Comput. Chem. 32, 1456 (2011)
4. **ALPB Model**: Ehlert et al., J. Chem. Theory Comput. 17, 4250 (2021)

## See Also

- [Parameter Registry Guide](archive/PARAMETER_MIGRATION_GUIDE.md)
- [TBLite Interface Documentation](../src/core/energy_calculators/qm_methods/tbliteinterface.h)
- [Ulysses Interface Documentation](../src/core/energy_calculators/qm_methods/ulyssesinterface.h)
