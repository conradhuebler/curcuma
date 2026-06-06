# GFN2 Parameter Extraction Documentation

## Overview

This document describes how GFN2 parameters were extracted from the TBLite reference implementation and integrated into Curcuma's native GFN2 method.

## Parameter Source

**TBLite Repository**: https://github.com/tblite/tblite
**Source File**: `tblite/src/tblite/xtb/gfn2.f90`
**Extraction Date**: 2025-01-07
**License**: LGPL-3.0-or-later
**Copyright**: (C) 2019-2024 Sebastian Ehlert and contributors

## Original Method Reference

**Publication**: C. Bannwarth, S. Ehlert, S. Grimme
**Journal**: *J. Chem. Theory Comput.* **2019**, *15*, 1652-1671
**DOI**: [10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176)

## Implementation Approach

### Design Philosophy: Data-Driven Parameter Loading

Instead of manually hardcoding each element (which would result in ~4000+ lines of code), we adopted a **compact, data-driven approach**:

1. **Direct Array Extraction**: TBLite Fortran parameter arrays copied to C++ as `constexpr` constants
2. **Loop-Based Population**: Single loop populates all 86 elements from arrays
3. **Easy Verification**: Direct 1:1 mapping to TBLite source makes verification trivial
4. **Maintainable**: Parameter updates only require array data changes, not code restructuring

**Result**: ~420 lines vs ~4000+ for manual element-by-element approach

## Extracted Parameters

### 1. Repulsion Parameters (Lines 64-103 in gfn2.f90)

```fortran
real(wp), parameter :: rep_alpha(86)  ! Exponential decay parameter
real(wp), parameter :: rep_zeff(86)   ! Effective nuclear charge
```

**Physical Meaning**:
- `rep_alpha`: Controls how quickly repulsion decays with distance
- `rep_zeff`: Effective nuclear charge felt by electrons during repulsion

**Usage in GFN2**: Exponential repulsion potential
`V_rep(R_AB) = (Z_eff_A + Z_eff_B) × exp(-α_avg × R_AB) / R_AB`

### 2. Hubbard Parameters (Lines 107-125 in gfn2.f90)

```fortran
real(wp), parameter :: hubbard_parameter(86)  ! Chemical hardness (γ_ss)
```

**Physical Meaning**:
- Base chemical hardness parameter for each element
- Represents the energy cost of adding or removing electrons
- Used in Coulomb energy calculation (ES2 component)

**Usage in GFN2**: Coulomb kernel damping
`γ_AB = 1 / √(R² + 0.5×(γ_AA + γ_BB))`

### 3. Shell Hubbard Corrections (Lines 127-171 in gfn2.f90)

```fortran
real(wp), parameter :: shell_hubbard(0:2, 86)  ! Shell-specific corrections
```

**Structure**: 3D array with corrections for s, p, d shells
**Values**: Stored as `1.0 + correction` (multiplicative factors)

**Physical Meaning**:
- Different atomic shells (s, p, d) have different hardnesses
- Corrections account for shell-specific electron-electron repulsion
- Critical for accurate Coulomb energy in transition metals

**Mapping to C++**:
```cpp
constexpr double SHELL_HUBBARD_CORR[86][3];
// [element_index][shell_index] where shell = 0 (s), 1 (p), 2 (d)
```

**Usage**:
```cpp
elem.gamma_ss = base_hubbard * SHELL_HUBBARD_CORR[idx][0];  // s-shell
elem.gamma_sp = base_hubbard * SHELL_HUBBARD_CORR[idx][1];  // p-shell
elem.gamma_pp = base_hubbard * SHELL_HUBBARD_CORR[idx][2];  // d-shell
```

### 4. Dispersion Parameters (Line 54 in gfn2.f90)

```fortran
real(wp), parameter :: s6 = 1.0_wp, s8 = 2.7_wp, a1 = 0.52_wp, a2 = 5.0_wp
```

**Physical Meaning**: DFT-D4 dispersion correction scaling parameters
**Status**: Documented but **not yet integrated** (D4 integration separate task)
**Expected Impact**: ~2-4 Eh for typical organic molecules

## Implementation Files

### Primary Files
- **Parameter Loader**: `src/core/energy_calculators/qm_methods/gfn2_params_loader.cpp`
- **Header**: `src/core/energy_calculators/qm_methods/gfn2_params_loader.h`
- **GFN2 Implementation**: `src/core/energy_calculators/qm_methods/gfn2.cpp`

### Data Structure

```cpp
struct ElementParams {
    int atomic_number;
    std::string symbol;

    // Repulsion parameters
    double rep_alpha;      // From TBLite rep_alpha(Z)
    double rep_zeff;       // From TBLite rep_zeff(Z)

    // Coulomb parameters (shell-resolved)
    double gamma_ss;       // From hubbard_parameter(Z) × shell_hubbard(0, Z)
    double gamma_sp;       // From hubbard_parameter(Z) × shell_hubbard(1, Z)
    double gamma_pp;       // From hubbard_parameter(Z) × shell_hubbard(2, Z)

    // Dispersion (stub for future D4 integration)
    double c6_base;        // TODO: Extract from D4 database
    double r4_over_r2;     // TODO: Extract from D4 database
};
```

## Verification Procedure

To verify parameters against TBLite source:

1. **Open TBLite source**: `external/tblite/src/tblite/xtb/gfn2.f90`
2. **Compare arrays**: Check that C++ `REP_ALPHA`, `REP_ZEFF`, etc. match Fortran arrays
3. **Element-by-element**: For specific element Z, verify:
   ```cpp
   // C++ (0-based indexing)
   REP_ALPHA[Z-1] == gfn2.f90:rep_alpha(Z)
   ```
4. **Commit hash**: Cross-reference against TBLite commit for reproducibility

## Missing Components (Future Work)

### 1. ~~D4 Dispersion Integration~~ ✅ COMPLETED (March 2026)
**Files**: `src/core/energy_calculators/qm_methods/gfn2.cpp`, `gfn2_xtb_params.hpp`
**Status**: ✅ IMPLEMENTED - D4 dispersion calculation integrated
**Parameters**: `s6=1.0, s8=2.7, a1=0.52, a2=5.0` (from gfn2.f90:54)
**Impact**: ~2-4 Eh contribution now calculated when USE_D4 is defined

### 2. ~~Shell-Resolved Hubbard Parameters~~ ✅ COMPLETED (March 2026)
**Files**: `gfn2_xtb_params.hpp`, `gfn2_params_loader.cpp`
**Status**: ✅ IMPLEMENTED - Exact SHELL_HUBBARD_CORR array extracted from TBLite
**Data Source**: `external/tblite/src/tblite/xtb/gfn2.f90` lines 127-171
**Implementation**: `getShellHubbard()` now computes `γ_shell = HUBBARD[Z] * (1.0 + SHELL_HUBBARD_CORR[Z][shell])`

### 3. Pair-Specific Hamiltonian Scaling (Remaining)
**Data Structure**: `PairParams` in `gfn2_params_loader.h`
**Status**: ⏳ Partially complete - ~15 organic pairs defined
**Future**: Extract complete pair database from TBLite for all 86×86 element pairs

## Quality Assurance

### Parameter Validation

The GFN2 implementation includes runtime validation:

```cpp
if (!m_param_db.hasElement(Z)) {
    CurcumaLogger::warn(fmt::format(
        "Element Z={} missing from GFN2 database - using FALLBACK parameters (INACCURATE)", Z));
}
```

**Expected Behavior**: With complete 86-element database, warnings should **never** appear.

### Energy Validation Test

**Test Case**: Dimethyl ether (C2H6O)
**Reference**: XTB binary with `--gfn2` flag
**Command**:
```bash
# Native GFN2 (without D4)
./curcuma -sp ch3och3.xyz -method gfn2

# XTB reference (with D4)
xtb ch3och3.xyz --gfn2
```

**Expected Results**:
- **Without D4**: Native should be within **1-2 Eh** of XTB reference
- **With D4**: Should match XTB to **<0.01 Eh** (after D4 integration)

## Acknowledgments

**TBLite Authors**: Sebastian Ehlert and contributors
**Original GFN2 Method**: Christoph Bannwarth, Sebastian Ehlert, Stefan Grimme
**Parameter Extraction**: Claude (Anthropic), supervised by Conrad Hübler

## License Compatibility

**TBLite**: LGPL-3.0-or-later
**Curcuma**: GPL-3.0-or-later
**Compatibility**: ✅ GPL-3.0 is compatible with LGPL-3.0 (derivative work allowed)

## References

1. C. Bannwarth, S. Ehlert, S. Grimme, *J. Chem. Theory Comput.* **2019**, *15*, 1652-1671
2. TBLite repository: https://github.com/tblite/tblite
3. GFN2-xTB method documentation: https://xtb-docs.readthedocs.io/

---

*Last Updated*: 2025-01-07
*Claude Generated*: Complete parameter extraction documentation for Phase 1 implementation
