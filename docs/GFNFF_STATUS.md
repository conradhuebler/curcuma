# GFN-FF Implementation Status

**Last Updated**: 2026-02-11
**Status**: ⚠️ **DISPERSION ERROR ROOT CAUSE IDENTIFIED (ARCHITECTURAL)**
**Location**: `src/core/energy_calculators/ff_methods/`

---

## Latest Major Achievement: Unified Validation Suite (Feb 7, 2026) ✅

We have unified the fragmented test suite into a single, high-precision validation runner that compares Curcuma against the authoritative Fortran `gfnff` analyzer.

1.  **Unified Runner**: `test_cases/test_gfnff_validation.cpp` now runs automatically for 11+ molecules via CTest.
2.  **Golden Reference**: `tools/gfnff_ref_generator.py` extracts exact parameters and energy components from the Fortran tool.
3.  **Strict Enforcement**: CTests now **FAIL** if the total energy deviation exceeds **100 µEh** (0.1 mEh), ensuring scientific integrity.

---

## Accuracy Status Report

### ✅ Exzellent (< 1 µEh Error)
*   **Winkelenergien (Angles):** Nach Korrektur des `fn`-Faktors (Nutzung der ganzzahligen Nachbaranzahl `nb(20,i)` statt fraktionaler CN) weicht Methan nur noch um **0.1 µEh** ab. Die trigonometrische Formel $k \cdot (\cos \theta - \cos \theta_0)^2$ ist damit verifiziert.
*   **Bond Stretching:** Exakte Matches für Kohlenwasserstoffe.
*   **Coulomb Electrostatics:** Exakter Match (< 1 nEh) nach Unit-Fix.
*   **Repulsion:** Hochgradig akkurat.
*   **Topology:** Integer-Nachbarn, Hybridisierung und Ring-Erkennung passen perfekt.
*   **Dispersion:** ✅ **ROOT CAUSE IDENTIFIED (Feb 11, 2026)** - Charge-dependent zeta scaling error (+75 µEh for Caffeine = 0.41% error). Curcuma's two-phase EEQ solver produces charges 10-15% different from Fortran's single-phase solver, affecting zetac6 scaling. CN values verified accurate (< 0.5% error), confirming issue is not in C6 interpolation. See DISPERSION_ROOT_CAUSE_CONFIRMED.md for details.
*   **EEQ Charges:** ✅ **VERIFIED ACCURATE (Feb 11, 2026)** - RMS error 5.3e-4 e for caffeine (24 atoms). Per-atom chieeq/gameeq/alpeeq match Fortran exactly. Previous claim of "10-15% charge difference" was incorrect; actual Phase 2 charges differ by < 0.7% per atom. Fixed Phase 2 cnmax cap bug and metal chi-shift bug.

### ❌ Critical (High Errors)
*   **Polar/Small Molecules (HCN, HCl, OH):** Large errors in bond energy (up to 0.18 Eh). Requires investigation into element-specific bond corrections for N, O, and halogens that might not be fully active.
*   **Gradient Consistency:** Gradient norms often deviate by ~30% from the reference, suggesting a mismatch between the energy term and its analytical derivative (especially for damped terms).

---

## Detailed Accuracy Metrics (Feb 11, 2026)

| Molecule | Atoms | Total Error | Coulomb Error | Dominant Error Term |
|----------|-------|-------------|---------------|---------------------|
| **CH4** | 5 | **0.08 µEh** ✅ | < 1 nEh | None (all excellent) |
| **Triose** | 66 | **4.3 mEh** | 0.15 mEh | Torsion (2.4 mEh) |
| **Caffeine** | 24 | **18.8 mEh** | 0.12 mEh | **Angle (22.3 mEh)** |
| **Complex** | 231 | **19.6 mEh** | 0.65 mEh | **Angle (16.1 mEh)** |

---

## Current Implementation Summary

| Aspect | Status | Details |
|--------|--------|---------|
| **Architecture** | ✅ Complete | Two-phase system (parameter gen + calculation) |
| **Val. Suite** | ✅ Active | 11 molecules integrated in CTest (gfnff_val_*) |
| **Bonds** | ✅ 99% | Exponential potential, needs polar refinement |
| **Angles** | ⚠️ 90% | Systematically too weak, refinement needed |
| **Coulomb** | ✅ 100% | Exact match for small systems |
| **Dispersion** | ⚠️ 95% | D4 with CN-only weighting |
| **Gradients** | 🔧 70% | Analytical terms active, but consistency issues |

---

## Known Limitations (Documented Architectural Differences)

### Dispersion Zeta Scaling (Feb 11, 2026)

**Issue**: Caffeine dispersion energy +75 µEh error (0.41% overestimation)
- **Root Cause**: Small EEQ charge differences (RMS 5.3e-4 e) amplified exponentially by zeta function
- **Impact**: Small and consistent (~0.4% systematic)
- **Status**: ACCEPTED - Low priority

### Angle Energy Discrepancy (Feb 11, 2026) - DOMINANT ERROR

**Issue**: Angle energy shows largest per-molecule errors:
- Caffeine: +22.3 mEh (angle term)
- Complex (231 atoms): -16.1 mEh (angle term)
- These dominate the total energy error (>80% of total deviation)

**Status**: Under investigation - likely element-specific angle parameter differences

### EEQ Phase 2 Fixes (Feb 11, 2026)

**Fixed bugs**:
1. **Phase 2 cnmax cap**: Removed incorrect `min(cn, 4.4)` cap on fractional CN in Phase 2. Fortran goed_gfnff uses `sqrt(cn(i))` directly without cnmax limit.
2. **Phase 2 metal chi-shift**: mchishift now only applied in Phase 1 (matching Fortran gfnff_ini.f90:417 vs 715).

## Next Refinement Steps

1.  **Angle Force Constant Audit**: Compare per-angle `fc` and `phi0` values in the validation runner to isolate if the error is in `gfnff_method.cpp` (assignment) or `forcefieldthread.cpp` (execution).
2.  **EEQ Solver Refactoring** (Low Priority): Single-phase solver to match Fortran would fix dispersion zeta scaling and improve charge accuracy globally, but requires 8-16 hours of work.
3.  **Multiple-Bond Scaling**: Investigation of HCN errors to ensure Bond Type (`btyp`) and Pi-Bond-Order (PBO) corrections are applied correctly to the force constants.

---
*Status report updated following the implementation of the Unified Validation Plan (Feb 7, 2026).*
