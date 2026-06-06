# Solvation Models in Curcuma

Implicit solvation support for GFN and MNDO methods via TBLite and Ulysses interfaces.

## Overview

Curcuma supports implicit solvation models for quantum chemistry calculations:

- **TBLite Interface**: CPCM, GB (Generalized Born), ALPB for GFN methods
- **Ulysses Interface**: GBSA (Generalized Born) for GFN and MNDO methods
- **Native GFN**: Solvation support planned for future implementation

## Supported Solvents

### Common Solvents (TBLite & Ulysses)
- `water` - Most polar solvent
- `methanol`, `ethanol`, `octanol`
- `acetone`, `dmso` (dimethyl sulfoxide), `dmf` (dimethylformamide)
- `acetonitrile`, `nitromethane`
- `benzene`, `toluene`, `aniline`
- `chloroform`, `dichloromethane`
- `thf` (tetrahydrofuran), `dioxane`, `furane`
- `diethyl ether`, `ethyl acetate`
- `hexane`, `hexadecane`
- `carbon disulfide`
- `phenol`, `benzaldehyde`
- `octanol wet` - Wet octanol for partition coefficients

### Special Keywords
- `none` - No solvation (vacuum/gas phase)
- `vacuum` - Equivalent to `none`

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

- [Parameter Registry Guide](PARAMETER_MIGRATION_GUIDE.md)
- [TBLite Interface Documentation](../src/core/energy_calculators/qm_methods/tbliteinterface.h)
- [Ulysses Interface Documentation](../src/core/energy_calculators/qm_methods/ulyssesinterface.h)
