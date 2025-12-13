# GFN-FF EEQ Charges Issue - CRITICAL BUG

**Date**: 2025-12-05
**Status**: 🔴 **BLOCKING** - Coulomb formula is correct but charges have wrong sign/magnitude
**Priority**: CRITICAL for $1000 bounty

---

## Problem Summary

The complete Coulomb energy formula is now correctly implemented with all three terms:
1. ✅ Pairwise: `Σ(i<j) [q_i*q_j*erf(γ_ij*r_ij) / r_ij²]`
2. ✅ Self-energy: `-Σ_i [q_i*χ_i]`
3. ✅ Self-interaction: `Σ_i [0.5*q_i²*(gam_i + √(2/π)/√(α_i))]`

**HOWEVER**, the EEQ atomic charges have:
- ❌ **Wrong sign** (reversed polarity)
- ❌ **Wrong magnitude** (too small by ~2-3x)

This causes the Coulomb energy to be wrong even though the formula itself is correct.

---

## Evidence

### OH (Hydroxyl Radical)

**Reference (Fortran GFN-FF)**:
- O: q = **-0.385** e⁻
- H: q = **+0.385** e⁻

**Curcuma (current)**:
- O: q = **+0.284** e⁻ ❌ (wrong sign, wrong magnitude)
- H: q = **-0.284** e⁻ ❌ (wrong sign, wrong magnitude)

### HCl (Hydrogen Chloride)

**Reference (Fortran GFN-FF)**:
- Cl: q = **-0.118** e⁻
- H: q = **+0.118** e⁻

**Curcuma (current)**:
- Cl: q = **+0.054** e⁻ ❌ (wrong sign, wrong magnitude)
- H: q = **-0.054** e⁻ ❌ (wrong sign, wrong magnitude)

---

## Impact on Accuracy

| Molecule | Coulomb (Reference) | Coulomb (Curcuma) | Error |
|----------|---------------------|-------------------|-------|
| **HH** | 0.0000 Eh | -0.0000 Eh | ✓ OK (neutral) |
| **HCl** | -0.00809 Eh | -0.00153 Eh | **81% error** |
| **OH** | -0.07651 Eh | -0.03274 Eh | **57% error** |

The Coulomb term is critical for polar molecules. This bug **blocks achieving < 1% accuracy** for HCl and OH.

---

## Root Cause Analysis

The bug is in **`calculateEEQCharges()`** in `src/core/energy_calculators/qm_methods/gfnff.cpp`.

### Possible Issues

1. **Sign Convention Error**:
   - EEQ solver may return charges with opposite sign convention
   - Fortran uses: electronegative atoms (O, Cl) → negative charge
   - Curcuma has: electronegative atoms → **positive** charge ❌

2. **Magnitude Error**:
   - Charges are ~2-3x too small
   - Could be missing scaling factor or unit conversion
   - Could be incorrect damping in EEQ matrix

3. **Solver Configuration**:
   - EEQ linear system `A·q = b` may have wrong b vector
   - Chemical hardness (gamma) may be missing from diagonal
   - Electronegativity (chi) may have wrong sign in b vector

---

## Code Location

**File**: `src/core/energy_calculators/qm_methods/gfnff.cpp`

**Function**: `Vector GFNFF::calculateEEQCharges(const Vector& cn, const std::vector<int>& hyb, const std::vector<int>& rings) const`

**Line**: ~2300-2500 (exact line TBD)

**Reference**: Fortran implementation in `external/gfnff/src/gfnff_eeq.f90` or similar

---

## Debugging Strategy

### Step 1: Verify EEQ Parameter Values

Check that `getEEQParameters()` returns correct values:
```cpp
// For O (Z=8):
chi = 1.6912  // Electronegativity
gam = -0.0312 // Chemical hardness
alp = 0.9053  // Damping parameter
cnf = 0.1574  // CN correction
```

### Step 2: Check EEQ Matrix Construction

The EEQ linear system should be:
```
A[i][i] = gam[i] + Σ_j≠i [erf(γ_ij*r_ij) / r_ij]
A[i][j] = erf(γ_ij*r_ij) / r_ij   (i ≠ j)
b[i] = -chi[i] - cnf[i]*sqrt(CN[i])
```

Solve: `A·q = b` for charges q

### Step 3: Verify Sign Convention

After solving, check if charges need sign flip:
```cpp
// Try both:
q_final = q_solved;        // Current (wrong)
q_final = -q_solved;       // May fix sign issue
```

### Step 4: Check Charge Normalization

Ensure charges sum to molecular charge:
```cpp
double q_sum = charges.sum();
if (abs(q_sum - molecular_charge) > 1e-6) {
    // Normalization needed
}
```

---

## Validation Test Case

**Molecule**: OH.xyz (simple 2-atom test)

**Expected Output**:
```
Atom 0 (O): q = -0.385 e⁻
Atom 1 (H): q = +0.385 e⁻
Sum: 0.000 e⁻
```

**Current Output**:
```
Atom 0 (O): q = +0.284 e⁻
Atom 1 (H): q = -0.284 e⁻
Sum: 0.000 e⁻
```

**Success Criteria**:
- [ ] Charges have correct sign (O negative, H positive)
- [ ] Charges have correct magnitude (|q| ≈ 0.385)
- [ ] Coulomb energy for OH: -0.0765 ± 0.0005 Eh

---

## Workaround (Temporary)

If EEQ fix is too complex, charges can be loaded from external XTB calculation:
1. Run: `xtb OH.xyz --gfnff --chrg 0 > xtb.out`
2. Extract charges from output
3. Load via `gfnff_charges` file

**NOT RECOMMENDED** - Defeats purpose of native implementation.

---

## Related Files

- `gfnff.cpp:calculateEEQCharges()` - EEQ solver (BUG HERE)
- `gfnff.cpp:getEEQParameters()` - Parameter lookup (likely OK)
- `external/gfnff/src/gfnff_eeq.f90` - Fortran reference
- `gfnff_validate.py` - Validation framework (already works)

---

## Next Steps

1. **Locate exact EEQ solver code** in gfnff.cpp
2. **Compare with Fortran reference** line-by-line
3. **Add debug output** to print EEQ matrix A and vector b
4. **Test sign flip** of charges after solve
5. **Validate** with OH.xyz until charges match reference

---

**Estimated Fix Time**: 1-2 hours if EEQ solver is straightforward
**Impact**: Fixes 50-80% of remaining Coulomb error
**Blocking**: Required for < 1% total energy accuracy

---

**Generated**: 2025-12-05
**By**: Claude Code - GFN-FF Implementation Debug
