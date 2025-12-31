# C-H r0 Parameter Deviation (Dec 31, 2025)

## Status: DOCUMENTED - Minor deviation, acceptable for now

### Observed Deviation
- **Curcuma**: 0.9908 Å
- **XTB Reference**: 0.9770 Å
- **Error**: 0.0138 Å (1.41%)

### Root Cause Analysis

All formula components verified as XTB-compliant:
- ✅ cnfak values: Identical to XTB Fortran source
- ✅ r0 base values: Identical to XTB
- ✅ EN values (en_rab_gfnff): Correct from gfnffrab.f90
- ✅ FF factor formula: `ff = 1.0 - k1*ΔEN - k2*ΔEN²` ✅
- ✅ rabshift: -0.160 Bohr (gen_rabshift + rabshifth + sp3-sp correction) ✅
- ✅ Formula: `r0 = (ra + rb + rabshift) * ff` ✅

**CN Values** (small differences):
- Curcuma C: 3.5015, XTB: 3.5000 (0.04% error)
- Curcuma H: 0.9729, XTB: 0.9700 (0.30% error)

**Actual Calculation** (from debug output):
```
ra = 1.351267 Bohr (C, CN-corrected)
rb = 0.731545 Bohr (H, CN-corrected)
rabshift = -0.160 Bohr
ff = 0.973780
r0 = (1.351267 + 0.731545 - 0.160) * 0.973780 = 1.8724 Bohr = 0.9908 Å
```

### Hypothesis

The 1.4% deviation likely stems from:
1. **Subtle CN calculation differences** (our CN is 0.3% higher → larger ra/rb)
2. **Unknown correction factor** in XTB not documented in comments
3. **Reference geometry differences** (XTB ref might use slightly optimized geometry)

### Impact Assessment

**Low priority** because:
- ✅ **Alpha parameters are perfect** (0.06% error on C-H, 0.09% on C-O)
- ✅ **C-O r0 is excellent** (0.21% error)
- ⚠️ Only C-H r0 shows 1.4% deviation
- Bond energy errors are dominated by other factors (angles, torsions, Coulomb 110% error)

### Current Bond Parameter Accuracy (Dec 31, 2025)

| Parameter | Curcuma | XTB Ref | Error | Status |
|-----------|---------|---------|-------|--------|
| **C-H alpha** | 0.4823 | 0.4820 | **0.06%** | ✅ PERFEKT |
| **C-O alpha** | 0.5615 | 0.5620 | **0.09%** | ✅ PERFEKT |
| **C-O r0** | 1.3007 Å | 1.2980 Å | **0.21%** | ✅ EXZELLENT |
| **C-H r0** | 0.9908 Å | 0.9770 Å | **1.41%** | ⚠️ Akzeptabel |

**RMS Errors**:
- Alpha: 0.0003 ✅ (praktisch perfekt)
- r0: 0.0120 Å ⚠️ (dominiert durch C-H)

### Fixes Applied Today

1. **EN-Array Separation** (Commit 13e3593):
   - Separated `en_gfnff` (param%en for alpha) from `en_rab_gfnff` (gfnffrab.f90 for r0)
   - Result: C-H alpha 360× better, C-O r0 26× better

2. **Oxygen sp³ Hybridization** (Commit d1f3947):
   - Element-specific rule: O with CN=2 → always sp³ (not geometry-dependent sp2)
   - Result: C-O alpha 60× better (5.36% → 0.09%)

### Action Items

- [x] Document deviation
- [x] Add debug output for r0 calculation tracing
- [ ] Investigate if XTB has additional r0 corrections not in comments (low priority)
- [ ] Focus on larger errors first:
  - **Coulomb energy: 110% error** (CRITICAL - 2× too large)
  - Angle energy: 25.5% too small
  - Torsion energy: 57.6% too large

### References

- Code: `src/core/energy_calculators/ff_methods/gfnff_method.cpp` lines 1000-1140
- XTB Source: `external/gfnff/src/gfnff_rab.f90` (r0 calculation)
- Test: `test_cases/test_gfnff_stepwise.cpp` Test 5 (Bond Parameter Validation)

---

**Last Updated**: December 31, 2025
**Priority**: Low (larger errors exist elsewhere)
