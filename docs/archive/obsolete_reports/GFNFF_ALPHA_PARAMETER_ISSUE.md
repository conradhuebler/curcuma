# GFN-FF Alpha Parameter Calculation Issue

**Date**: December 31, 2025
**Status**: UNRESOLVED - Systematic error vs XTB 6.6.1
**Impact**: Bond parameter accuracy, affects total energy

## Summary

The alpha (α) parameter calculation in GFN-FF bond parameters shows **systematic errors** compared to XTB 6.6.1 reference:
- **C-H bonds**: 2% too small (need 27% increase)
- **C-O bonds**: 74% too large (need 43% decrease)

## Test Data (CH₃OCH₃)

### Bond Parameters Comparison

| Bond Type | Parameter | Curcuma | XTB 6.6.1 | Error | Error % |
|-----------|-----------|---------|-----------|-------|---------|
| **C-H** (×6) | alp | 0.4723 | 0.4820 | -0.0097 | -2.0% |
| | r0 | 0.9886 Å | 0.9770 Å | +0.0116 Å | +1.2% |
| | fqq | 1.0000 | 1.0060 | -0.0060 | -0.6% |
| **C-O** (×2) | alp | **0.9784** | **0.5620** | **+0.4164** | **+74.1%** ❌ |
| | r0 | 1.2977 Å | 1.2980 Å | -0.0003 Å | -0.02% |
| | fqq | 1.0000 | 1.0470 | -0.0470 | -4.5% |

**Critical Finding**: C-O bonds have **74% error in alpha**, while r0 is nearly perfect (<0.03% error)!

## Current Implementation

### Formula (from code comment line 945)
```
α = srb1 * (1 + fsrb2*ΔEN² + srb3*(bstrength-1))
```

### Actual Code (line 1583-1588)
```cpp
double alpha_term1 = fsrb2 * en_diff * en_diff;
double alpha_term2 = srb3 * bstrength;  // ← DISCREPANCY!
double alpha_sum = 1.0 + alpha_term1 + alpha_term2;
params.alpha = srb1 * alpha_sum;
```

**Discrepancy**: Code uses `srb3 * bstrength`, comment says `srb3 * (bstrength - 1)`

## Attempted Fix: `(bstrength - 1)`

### Hypothesis
Comment formula suggests using `(bstrength - 1)` so that:
- Single bonds (bstrength=1.0): term2 = 0
- Double bonds (bstrength=1.24): term2 = 0.0609
- Triple bonds (bstrength=1.98): term2 = 0.2487

### Test Results

Changed line 1586 from:
```cpp
double alpha_term2 = srb3 * bstrength;
```
To:
```cpp
double alpha_term2 = srb3 * (bstrength - 1.0);
```

**Results**:

| Bond Type | Before Fix | After Fix | XTB Target | Verdict |
|-----------|-----------|-----------|------------|---------|
| **C-H** | 0.4723 | **0.3776** | 0.4820 | ❌ **WORSE** (21.6% error!) |
| **C-O** | 0.9784 | **0.8838** | 0.5620 | ⚠️ Slightly better but still 57% error |

### Analysis

The fix makes C-H bonds **WORSE** while only marginally improving C-O bonds. This proves:

1. **Comment formula is WRONG** (or outdated for XTB 6.6.1)
2. **Original code `srb3 * bstrength` is closer to correct** for C-H
3. **Neither formula is correct** for C-O bonds

## Root Cause Investigation

### Constants Used
```cpp
double srb1 = 0.3731;   // Fortran: gen%srb1
double srb2 = 0.3171;   // Fortran: gen%srb2 (scaled by fsrb2)
double srb3 = 0.2538;   // Fortran: gen%srb3
```

### EN Difference Calculation
```cpp
double en1 = en_gfnff[z1 - 1];  // C: 2.49640820, O: 2.81007174, H: 2.30085633
double en2 = en_gfnff[z2 - 1];
double en_diff = std::abs(en1 - en2);
```

For C-H: ΔEN = |2.496 - 2.301| = 0.195
For C-O: ΔEN = |2.496 - 2.810| = 0.314

### Alpha Calculation Breakdown (Original Formula)

**C-H Bond** (Z=6, Z=1, bstrength=1.0):
```
term1 = 0.3171 * 0.195² = 0.0121
term2 = 0.2538 * 1.0 = 0.2538
sum = 1.0 + 0.0121 + 0.2538 = 1.2659
alpha = 0.3731 * 1.2659 = 0.4723  (XTB: 0.4820, error: -2%)
```

**C-O Bond** (Z=6, Z=8, bstrength=1.0):
```
term1 = 0.3171 * 0.314² = 0.0313
term2 = 0.2538 * 1.0 = 0.2538
sum = 1.0 + 0.0313 + 0.2538 = 1.2851
alpha = 0.3731 * 1.2851 = 0.4794  (XTB: 0.5620, error: -14.7%)
```

**Wait!** Calculated C-O alpha is 0.4794, but test shows 0.9784!

This means something ELSE is wrong - likely the **fsrb2 scaling** or **EN values**.

## Possible Causes (Ranked by Likelihood)

### 1. **fsrb2 Metal Scaling Applied Incorrectly** ⭐⭐⭐⭐⭐
Lines 1569-1579 modify `fsrb2` based on metal type:
```cpp
if (mtyp1 == 4 || mtyp2 == 4) {
    fsrb2 = -srb2 * 0.22;  // TM: negative!
} else if (mtyp1 > 0 || mtyp2 > 0) {
    fsrb2 = srb2 * 0.28;   // Other metals
} else {
    fsrb2 = srb2;  // Non-metals (C, H, O)
}
```

For C-H and C-O (both non-metals), `fsrb2 = srb2 = 0.3171`

**But test shows alp=0.9784 for C-O!** This is ~2× the calculated value.

### 2. **EN Values Different in XTB 6.6.1** ⭐⭐⭐⭐
Current `en_gfnff` values are from older XTB source. XTB 6.6.1 may use different values.

**Check**: Compare `gfnff_rab.f90` from XTB 6.6.1 source vs current values.

### 3. **bstrength Calculated Wrong for C-O** ⭐⭐⭐
C-O should be single bond (bstrength=1.0), but maybe hybridization detection is wrong?

**Check**: Add debug output for `bstrength` calculation.

### 4. **Formula Actually Different in XTB 6.6.1** ⭐⭐
XTB 6.6.1 might use a completely different alpha formula than documented.

**Check**: Read XTB 6.6.1 Fortran source `gfnff_ini.f90` lines ~1200-1300.

### 5. **Unit Conversion Error** ⭐
Alpha is dimensionless, but maybe there's a Bohr/Angstrom issue somewhere?

## Next Steps (Priority Order)

1. **Add Debug Output** to print all intermediate values:
   ```cpp
   - en1, en2, en_diff
   - fsrb2 (before and after metal scaling)
   - bstrength
   - term1, term2, sum
   - Final alpha
   ```

2. **Verify bstrength** for C-O bonds (should be 1.0 for single bonds)

3. **Check XTB 6.6.1 Source Code** for actual alpha formula
   - File: `gfnff_ini.f90` around line 1200-1300
   - Look for variable `alp` or `alpha`

4. **Compare EN Values** between current code and XTB 6.6.1
   - File: `gfnff_rab.f90` lines 62-81

5. **Test with Known Reference** (H₂, HCl where XTB shows 0.000% error)

## Impact on Energy Accuracy

**Current Status** (with 74% alpha error on C-O):
- Bond energy: 7% too large (-1.302 vs -1.216 Eh)
- Total energy: 11.6% error

**Expected with Correct Alpha**:
- Bond energy should improve significantly
- Total energy error should drop below 5%

## References

- **Code**: `src/core/energy_calculators/ff_methods/gfnff_method.cpp` lines 933-1624
- **Test**: `test_cases/test_gfnff_stepwise.cpp` Test 5 (Bond Parameter Validation)
- **Reference Data**: `test_cases/reference_data/ch3och3_reference.json`
- **XTB Source**: `external/gfnff/src/gfnff_ini.f90` (need to check version 6.6.1)

## Conclusion

The alpha parameter calculation has a **systematic error** that differs by bond type:
- C-H: Small error (-2%), acceptable
- C-O: **CRITICAL 74% error** ❌

The formula in the code comment is **WRONG** or **OUTDATED**.
Neither `srb3 * bstrength` nor `srb3 * (bstrength - 1)` produces correct results.

**BLOCKER**: Cannot achieve 1e-6 accuracy without fixing alpha calculation first.

---

**Last Updated**: December 31, 2025
**Requires**: XTB 6.6.1 Fortran source code verification
