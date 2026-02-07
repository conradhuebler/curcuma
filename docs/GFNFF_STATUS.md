# GFN-FF Implementation Status

**Last Updated**: 2026-02-07
**Status**: ⚠️ **UNIFIED VALIDATION ACTIVE (ACCURACY REFINEMENT ONGOING)**
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
*   **Bond Stretching:** Exact matches for standard hydrocarbons (Methane, CH3OCH3).
*   **Coulomb Electrostatics:** Exact match (< 1 nEh) for small systems after the Jan 29 unit fix.
*   **Repulsion:** Core-core repulsion is highly accurate across all tested organic molecules.
*   **Topology:** Coordination numbers and Hybridization detection match the reference exactly.

### ⚠️ Under Investigation (Abweichung 10-200 µEh)
*   **Angle Bending:** Systematically underestimated (e.g., Methane: 60 µEh vs Ref 69 µEh). This indicates a missing scaling factor in the force constant generation or a subtle difference in the damping function.
*   **Dispersion:** Small deviations (~200 µEh) observed. Likely due to subtle differences in the $a_1, a_2$ BJ-damping parameters or the $r^6$ vs $r^6+R_0^6$ formulation in `forcefieldthread.cpp`.
*   **EEQ Charges:** Minor deviations (0.01-0.03 e) persist in Phase 2 energy charges, which propagates to other terms.

### ❌ Critical (High Errors)
*   **Polar/Small Molecules (HCN, HCl, OH):** Large errors in bond energy (up to 0.18 Eh). Requires investigation into element-specific bond corrections for N, O, and halogens that might not be fully active.
*   **Gradient Consistency:** Gradient norms often deviate by ~30% from the reference, suggesting a mismatch between the energy term and its analytical derivative (especially for damped terms).

---

## Detailed Accuracy Metrics (Current CTest Status)

| Molecule | Atoms | Total Energy Error (Eh) | Status |
|----------|-------|-------------------------|--------|
| **Methane** | 5 | +0.000193 | ❌ (Target: < 0.0001) |
| **CH3OCH3** | 9 | +0.000210 | ❌ (Target: < 0.0001) |
| **Benzene** | 12 | +0.000300 | ❌ (Target: < 0.0001) |
| **HCN** | 3 | +0.212853 | ❌ (Critical) |
| **HH** | 2 | +0.000023 | ✅ (Within limit) |

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

## Next Refinement Steps

1.  **Angle Force Constant Audit**: Compare per-angle `fc` and `phi0` values in the validation runner to isolate if the error is in `gfnff_method.cpp` (assignment) or `forcefieldthread.cpp` (execution).
2.  **BJ-Damping Synchronization**: Verify the GFN-FF specific damping formula $E \propto 1/(r^n + R_0^n)$ vs Curcuma's current implementation.
3.  **Multiple-Bond Scaling**: Investigation of HCN errors to ensure Bond Type (`btyp`) and Pi-Bond-Order (PBO) corrections are applied correctly to the force constants.

---
*Status report updated following the implementation of the Unified Validation Plan (Feb 7, 2026).*
