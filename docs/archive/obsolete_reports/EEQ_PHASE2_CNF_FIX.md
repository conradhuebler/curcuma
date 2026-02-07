# EEQ Phase 2 CNF Term Fix - Implementation Notes

## Date: January 4, 2026
## Author: Claude Sonnet 4.5 (Architecture Analysis & Implementation Planning)
## Status: Applied but Effectiveness Under Investigation

## Executive Summary

Applied critical bug fix in EEQ Phase 2 charge calculation where CNF term was incorrectly included in RHS vector. Changed `use_cnf_term=true` to `false` at line 826 in `eeq_solver.cpp`. However, test results show Phase 2 RMS error unchanged at 1.55e-02, indicating the architectural issue may be more complex than initially anticipated.

## Root Cause Analysis

### Fortran Reference (XTB gfnff_ini.f90)

**Phase 1: Topology Charges (line 411)**
```fortran
topo%chieeq(i) = -param%chi(ati)+dxi(i)+param%cnf(ati)*sqrt(dum)
                                        ↑ WITH CNF term
```

**Phase 2: Parameter Preparation (line 715)**
```fortran
topo%chieeq(i) = -param%chi(at(i))+dxi(i)  ! NO CNF term!
                                           ↑ OVERWRITTEN, removes CNF
```

**Phase 3: Geometry Charges (gfnff_engrad.F90:1504)**
```fortran
x(i) = topo%chieeq(i)+param%cnf(at(i))*sqrt(cn(i))
       ↑ from Phase 2 (NO CNF)  + adds NEW CNF with fractional CN
```

### Three Critical Distinctions

| Phase | Chi Formula | CN Type | Purpose |
|-------|------------|---------|---------|
| Phase 1 | -chi + dxi + cnf√(nb_integer) | Integer neighbor count | Topology charges (qa) for parameter generation |
| Phase 2 Prep | -chi + dxi | Integer (unused for RHS) | Parameter base value |
| Phase 3 | -chi + dxi + cnf√(cn_fractional) | Fractional erf-counting | Energy charges (q) for electrostatics |

## Implementation

### Code Change

**File**: `src/core/energy_calculators/ff_methods/eeq_solver.cpp`
**Location**: Line 826 (Phase 2 solve call)

```cpp
// BEFORE (Line 788, WRONG):
current_charges = solveEEQ(A, atoms, cn, dxi, total_charge, topology, true);

// AFTER (Line 826, CORRECT):
current_charges = solveEEQ(A, atoms, cn, dxi, total_charge, topology, false);
```

### Documentation Updates

**Lines 823-825**: Added clear comments explaining Phase 2 behavior
```cpp
// Solve ONCE with corrected parameters (Phase 2: NO CNF term!)
// CRITICAL FIX (Jan 4, 2026): Fortran Phase 2 preparation (gfnff_ini.f90:715)
// overwrites chieeq WITHOUT CNF term. Only Phase 1 (line 411) includes CNF.
```

**Lines 1013-1019**: Added comprehensive solveEEQ documentation
```cpp
// CRITICAL: CNF term handling follows Fortran reference exactly:
// - Phase 1 (goedeckera): chieeq includes CNF with integer nb (gfnff_ini.f90:411)
// - Phase 2 preparation: chieeq OVERWRITTEN without CNF (gfnff_ini.f90:715)
// - Phase 3 (goed_gfnff): RHS adds CNF with fractional CN (gfnff_engrad.F90:1504)
```

**Lines 646-682**: Added complete function documentation for calculateCharges()
```cpp
/**
 * PHASE 1: Topology Charges (goedeckera equivalent)
 * - Uses integer neighbor count (nb) from topology
 * - chieeq = -chi + dxi + cnf*sqrt(nb)  [WITH CNF!]
 * ...
 * PHASE 2: Parameter Preparation + Final Charges
 * - Overwrites chieeq WITHOUT CNF: chieeq = -chi + dxi  [NO CNF!]
 * ...
 * CRITICAL BUG FIX (Jan 4, 2026):
 * - Phase 2 was incorrectly using use_cnf_term=true (line 788)
 * - Changed to use_cnf_term=false to match Fortran line 715
 * - Error reduced from 1.55e-02 e RMS to <1e-5 e (expected)
 */
```

## Validation Results

### Test Setup
**Molecule**: CH₃OCH₃ (dimethyl ether, 9 atoms)
**Test Command**: `./test_cases/test_gfnff_stepwise`
**XTB Reference**: 6.6.1 official GFN-FF

### Test Results (After Fix)

| Metric | Value | Status |
|--------|-------|--------|
| Phase 2 RMS Error | 1.55e-02 e | ⚠️ **UNCHANGED** |
| Carbon charge (Curcuma) | 0.02537 e | 24% too high |
| Carbon charge (XTB Ref) | 0.02055 e | |
| Oxygen charge (Curcuma) | -0.3261 e | 11% too low |
| Oxygen charge (XTB Ref) | -0.3648 e | |

### Key Finding

The Phase 2 RMS error remains at 1.55e-02 after the fix was applied and the code was recompiled. This indicates one of the following:

1. **GFNFF class has embedded EEQ solver**: May not use the standalone EEQSolver class
2. **Different code path**: Charges may come from cached values or alternative calculation
3. **Architecture mismatch**: The two-phase model in EEQSolver may not match how GFNFF uses charges
4. **Multiple CN calculations**: Different CN calculation methods may be used in different phases

## Architectural Questions Requiring Investigation

### Q1: Is GFNFF using EEQSolver or embedded implementation?

**Check**: Search for `calculateCharges` calls in `gfnff_method.cpp`

```cpp
// Does it call EEQSolver::calculateCharges()?
// Or does it have its own implementation?
```

### Q2: Are charges being cached or recomputed?

**Check**: Verify that `gfnff.Calculation()` actually triggers charge recalculation

```cpp
// Are charges cached from previous calls?
// Is there a cache invalidation mechanism?
```

### Q3: Which CNF term is relevant?

**Hypothesis**: The current implementation may be adding CNF correctly in Phase 2 but with wrong CN type (integer instead of fractional)

**Check**:
- Phase 1 uses integer `nb` (topology neighbor count)
- Phase 3 uses fractional `cn` (erf-counting)
- Curcuma may be mixing the two

## References

### Fortran Source Code
- **Phase 1**: `external/gfnff/src/gfnff_ini.f90:411`
- **Phase 2**: `external/gfnff/src/gfnff_ini.f90:715`
- **Phase 3**: `external/gfnff/src/gfnff_engrad.F90:1504`

### Curcuma Implementation
- **EEQSolver**: `src/core/energy_calculators/ff_methods/eeq_solver.cpp:826`
- **GFNFF Class**: `src/core/energy_calculators/ff_methods/gfnff_method.cpp`
- **Test**: `test_cases/test_gfnff_stepwise.cpp:562`

## Next Steps

### For Investigation
1. [ ] Trace execution path from `gfnff.Calculation()` to charge calculation
2. [ ] Check if `gfnff.getCharges()` returns EEQSolver output or cached values
3. [ ] Verify CN calculation (integer vs fractional in different phases)
4. [ ] Add debug output to see which code path is executing

### For Resolution
1. [ ] Identify if GFNFF has embedded EEQ solver that needs updating
2. [ ] Ensure GFNFF delegates to EEQSolver correctly
3. [ ] Validate that CN types match Fortran reference
4. [ ] Re-test after architectural fix

## Conclusion

The critical bug fix has been properly identified and applied to the code. However, the test results indicate the issue is more complex than a simple one-line CNF term fix. The GFNFF class architecture may be using a different charge calculation path or caching mechanism that isn't affected by the EEQSolver Phase 2 modification.

Further investigation is needed to:
1. Understand how GFNFF calculates and caches charges
2. Verify the correct integration between GFNFF and EEQSolver
3. Potentially identify additional code paths that need updating

The fix itself is correct based on Fortran reference analysis, but its effectiveness depends on proper integration with the calling code.
