# MNDO Two-Electron Integrals

**Location**: `src/core/energy_calculators/qm_methods/integrals/MNDOIntegrals.hpp`
**Status**: ✅ Fully implemented (November 2025)
**License**: LGPL 2.1+ (extracted from Ulysses)

## Purpose

Standalone implementation of MNDO-type two-electron repulsion integrals using the Dewar-Thiel multipole expansion method for semi-empirical quantum chemistry.

## What are MNDO Integrals?

In semi-empirical methods (MNDO, AM1, PM3, PM6), the computationally expensive four-center electron repulsion integrals are simplified using the **NDDO approximation** (Neglect of Diatomic Differential Overlap):

```
(μν|λσ) ≈ 0  unless  (μ,ν on same atom) AND (λ,σ on same atom)
```

This reduces the problem to calculating **two-center integrals**:

```
γ_AB = <χ_μ(A) χ_ν(A) | 1/r₁₂ | χ_λ(B) χ_σ(B)>
```

where χ are atomic orbitals and r₁₂ is the electron-electron distance.

## Dewar-Thiel Method

Instead of numerical integration, Dewar and Thiel (1977) derived **closed-form analytical formulas** by expanding charge distributions into multipole moments:

- **Monopole** (l=0): s-orbital charge (spherical)
- **Dipole** (l=1): p-orbital charge (directional)
- **Quadrupole** (l=2): d-orbital charge (complex shapes)

All integrals are expressed as sums of terms like:

```
γ = 1/√(R² + f(D_params)²)
```

where:
- **R**: Interatomic distance
- **D_params**: Element-specific orbital expansion parameters
  - D₁: Dipole expansion (p-orbitals)
  - D₂: Quadrupole expansion (d-orbitals)

## Key Features

- ✅ **Analytically exact** within MNDO approximation (no numerical integration)
- ✅ **Computationally efficient** (only sqrt and arithmetic)
- ✅ **Thread-safe** (pure functions, no global state)
- ✅ **Comprehensive documentation** (pedagogical, learning-oriented)
- ✅ **All derivatives** (first and second order for gradients/Hessians)

## API Functions

### Main Integral

```cpp
double curcuma::mndo::mndo_multipole_integral(
    int l1, int m1,        // Quantum numbers for orbital on atom A
    int l2, int m2,        // Quantum numbers for orbital on atom B
    double R_AB,           // Interatomic distance (Bohr)
    double rho_sum,        // Sum of orbital exponents
    const std::vector<double>& D_params  // {D1_A, D2_A, D1_B, D2_B}
);
```

### First Derivative (for forces)

```cpp
double curcuma::mndo::mndo_multipole_integral_dR(...);  // dγ/dR
```

### Second Derivative (for Hessians)

```cpp
double curcuma::mndo::mndo_multipole_integral_dR2(...);  // d²γ/dR²
```

## Quantum Numbers

| l | m | Orbital Type |
|---|---|--------------|
| 0 | 0 | s |
| 1 | -1 | p_y |
| 1 | 0 | p_z |
| 1 | +1 | p_x |
| 2 | -2 | d_yy |
| 2 | -1 | d_yz |
| 2 | 0 | d_z² |
| 2 | +1 | d_xz |
| 2 | +2 | d_xx |
| 2 | 3 | d_xy (special encoding) |

## Example Usage

```cpp
#include "src/core/energy_calculators/qm_methods/integrals/MNDOIntegrals.hpp"

// Calculate (ss|ss) integral
double R_AB = 2.0;  // Bohr
double rho = 1.2 + 1.5;  // ρ_A + ρ_B
std::vector<double> D = {0.0, 0.0, 0.0, 0.0};

double gamma_ss = curcuma::mndo::mndo_multipole_integral(0, 0, 0, 0, R_AB, rho, D);
// Result: 1/√(R² + ρ²) = basic Coulomb repulsion

// Calculate (p_z p_z|p_z p_z) integral
std::vector<double> D_pz = {0.8, 0.0, 0.9, 0.0};  // D1_A, D2_A, D1_B, D2_B
double gamma_pzpz = curcuma::mndo::mndo_multipole_integral(1, 0, 1, 0, R_AB, rho, D_pz);

// Calculate gradient for geometry optimization
double dGamma_dR = curcuma::mndo::mndo_multipole_integral_dR(1, 0, 1, 0, R_AB, rho, D_pz);
```

## Use Cases

1. **Native PM3/PM6 implementation** - Semi-empirical methods without external dependencies
2. **Educational purposes** - Teaching MNDO theory with transparent code
3. **Method development** - Experimenting with new semi-empirical parameterizations
4. **Cross-validation** - Verifying external library calculations (Ulysses, MOPAC)

## References

1. **M. J. S. Dewar, W. Thiel**, *Theor. Chim. Acta*, **46**, 89 (1977)
   DOI: [10.1007/BF00548085](https://doi.org/10.1007/BF00548085)
   → Original MNDO multipole integral formulas

2. **M. J. S. Dewar, W. Thiel**, *J. Am. Chem. Soc.*, **99**, 4899 (1977)
   DOI: [10.1021/ja00457a004](https://doi.org/10.1021/ja00457a004)
   → MNDO parameterization

3. **W. Thiel, A. A. Voityuk**, *Theor. Chim. Acta*, **81**, 391 (1992)
   DOI: [10.1007/BF01134863](https://doi.org/10.1007/BF01134863)
   → Extension to d-orbitals

4. **J. J. P. Stewart**, *J. Mol. Model.*, **13**, 1173 (2007)
   DOI: [10.1007/s00894-007-0233-4](https://doi.org/10.1007/s00894-007-0233-4)
   → PM6 method

## Implementation Notes

- **Source**: Extracted from [Ulysses](https://github.com/Helmholtz-Munich/Ulysses) (LGPL 2.1+)
- **Original file**: `external/ulysses-main/core/src/basissets/2ElectronDewar.hpp`
- **Extraction date**: November 2025
- **Contributors**: Conrad Hübler (extraction), Claude (pedagogical documentation)

## Unit System

All quantities in **atomic units**:
- Distance: Bohr (a₀ = 0.529177 Å)
- Energy: Hartree (Eh = 27.2114 eV)
- Return values: Hartree for integrals, Hartree/Bohr for derivatives
