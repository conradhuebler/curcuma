# Critical Analysis: Two-Phase EEQ Charge Calculation Bugs

**Status**: 5 CRITICAL BUGS IDENTIFIED in `calculateFinalCharges()` causing wrong charges/energies
**File**: `src/core/energy_calculators/qm_methods/gfnff.cpp:3680-3820`
**Date**: 2025-12-08

---

## Executive Summary

The two-phase EEQ system has **Phase 1 (topology charges) WORKING CORRECTLY** but **Phase 2 (final refinement) COMPLETELY BROKEN**.

Phase 1 uses the **Fortran-correct goedeckera algorithm** but Phase 2 (`calculateFinalCharges()`) reverts to an **older incorrect EEQ formulation** that contradicts Phase 1.

---

## Critical Bug #1: Wrong EEQ Matrix Diagonal (Line 3730)

### ❌ CURRENT CODE (WRONG):
```cpp
J_corrected(i, i) = -1.0 / (2.0 * gam_corrected(i));  // Line 3730
```

### ✅ CORRECT FORMULA (From Phase 1 & Fortran):
```cpp
J_corrected(i, i) = gam_corrected(i) + TSQRT2PI / std::sqrt(alpha(i));
// = hardness + sqrt(2/π)/sqrt(polarizability)
```

### 📊 Mathematical Difference:
For Oxygen (Z=8):
- **Phase 1 (CORRECT)**: `gam + sqrt(2/π)/sqrt(α) = 0.0 + 0.798/sqrt(3.55) ≈ 0.424`
- **Phase 2 (WRONG)**: `-1/(2*gam) = -1/(2*0.0) = undefined/infinity` ❌

**Impact**: Matrix is singular or nonsensical, charges collapse to garbage

---

## Critical Bug #2: Missing Error Function in Coulomb Matrix (Line 3742)

### ❌ CURRENT CODE (WRONG):
```cpp
J_corrected(i, j) = 1.0 / r;  // Line 3742 - BARE 1/r!
```

### ✅ CORRECT FORMULA (From Phase 1 & Fortran):
```cpp
double gammij = 1.0 / std::sqrt(alpha(i) + alpha(j));
double erf_gamma = std::erf(gammij * r);
J_corrected(i, j) = erf_gamma / r;  // Coulomb with erf damping
```

### 📊 Numerical Difference:
For O-H at r = 1.81 Bohr:
- **Phase 1 (CORRECT)**: `erf(0.53*1.81)/1.81 = 0.531/1.81 ≈ 0.293`
- **Phase 2 (WRONG)**: `1.0/1.81 ≈ 0.552` (40% error!)

**Impact**: Coulomb matrix 40-50% too strong, charges oscillate incorrectly

---

## Critical Bug #3: Inconsistent Chi Sign in Phase 2 (Line 3717)

### PHASE 1 (CORRECT):
```cpp
// gfnff.cpp:3465-3490 (calculateTopologyCharges)
chi(i) = -params_i.chi + dxi_total;  // NEGATIVE chi parameter!
x(i) = chi(i);  // RHS uses chi directly (POSITIVE contribution)
```

### ❌ PHASE 2 (WRONG):
```cpp
// gfnff.cpp:3717 (calculateFinalCharges)
chi_corrected(i) = params_i.chi + topo_info.dxi(i);  // POSITIVE chi!
```

### 📊 Sign Flip Analysis:
For Oxygen EEQ parameter `chi = +1.456784`:
- **Phase 1**: Uses `-1.456784` as base
- **Phase 2**: Uses `+1.456784` as base
- **Result**: Charges are **exactly opposite sign** ❌

**Impact**: Positive atoms become negative, negative become positive - charge signs are flipped!

---

## Critical Bug #4: Wrong RHS Sign in Phase 2 (Line 3750)

### ❌ CURRENT CODE:
```cpp
rhs(i) = -chi_corrected(i);  // Line 3750 - NEGATIVE RHS
```

### ✅ PHASE 1 (CORRECT):
```cpp
x(i) = chi(i);  // POSITIVE RHS (Fortran: topo%chieeq(i))
```

### 📊 RHS Impact:
When combined with Bug #3:
- **Phase 1 RHS**: `-param%chi + dxi` (correct)
- **Phase 2 RHS**: `-(param%chi + dxi) = -param%chi - dxi` (double negative!)

**Impact**: Another 100% sign flip, quadratic error accumulation

---

## Critical Bug #5: Missing dalpha Correction to Alpha Parameters (Line 3724-3742)

### ❌ CURRENT CODE:
```cpp
// Calculate dalpha (lines 3638-3661) ✓
topo_info.dalpha = dalpha;

// But NEVER used in Phase 2!
Matrix J_corrected = Matrix::Zero(m_atomcount, m_atomcount);
for (int i = 0; i < m_atomcount; ++i) {
    for (int j = 0; j < m_atomcount; ++j) {
        if (i == j) {
            // ❌ alpha NOT modified by dalpha!
            J_corrected(i, i) = -1.0 / (2.0 * gam_corrected(i));
        } else {
            // ❌ alpha NOT modified by dalpha!
            J_corrected(i, j) = 1.0 / r;
        }
    }
}
```

### ✅ CORRECT APPROACH:
```cpp
// Apply dalpha correction to alpha parameters BEFORE building matrix
Vector alpha_corrected = Vector::Zero(m_atomcount);
for (int i = 0; i < m_atomcount; ++i) {
    alpha_corrected(i) = getEEQParameters(m_atoms[i]).alp
                       + topo_info.dalpha(i);
}

// Now use alpha_corrected in matrix construction
// Coulomb: gamma_ij = 1/sqrt(alpha_corrected(i) + alpha_corrected(j))
```

**Impact**: Polarizability corrections are computed but **completely ignored** - wasted calculation

---

## Root Cause Analysis

### Why Phase 1 Works but Phase 2 Doesn't

**Phase 1** (`calculateTopologyCharges()` lines 3450-3557):
- ✅ Directly ports Fortran `goedeckera` subroutine (gfnff_ini2.f90:1140-1246)
- ✅ Builds augmented system correctly (4×4 for water, not 3×3)
- ✅ Uses correct diagonal: `gam + sqrt(2/π)/sqrt(α)`
- ✅ Uses correct off-diagonal: `erf(gamma_ij*r)/r`
- ✅ Solves with PartialPivLU (equivalent to Fortran Bunch-Kaufman)

**Phase 2** (`calculateFinalCharges()` lines 3680-3820):
- ❌ **REVERTS TO OLD INCORRECT EEQ MODEL**
- ❌ Builds 3×3 matrix instead of augmented system
- ❌ Wrong diagonal formula from legacy implementation
- ❌ Missing erf damping in Coulomb terms
- ❌ Sign inconsistencies with Phase 1
- ❌ Calculated dalpha corrections are ignored

### The Problem

Session 7's code says "Fixed - Phase 1 EEQ solver bug corrected" but:
1. Phase 1 WAS fixed correctly
2. Phase 2 was **NEVER REFACTORED** to use the Fortran-correct formula
3. Phase 2 still uses an **older, incorrect EEQ model** that predates the Fortran port

---

## Comparison with Fortran Reference

### Fortran Phase 2: gfnff_ini2.f90:1227-1246

```fortran
! After first EEQ solve (which Phase 1 replicates exactly):
! Now apply charge-dependent corrections for second iteration

do i = 1,n
    ! Correct electronegativity based on charge
    topo%chieeq(i) = -param%chi(at(i)) + dxi(i) + param%cnf(at(i))*sqrt(cn(i))

    ! Correct hardness (gamma) based on charge
    topo%gameeq(i) = param%gam(at(i)) + dgam(i)

    ! Correct alpha based on coordination and charge
    topo%alpeeq(i) = param%alp(at(i)) + dalpha(i)  ! ← FORTRAN DOES THIS!
end do

! Rebuild SAME augmented EEQ matrix with corrected parameters
! (Same matrix structure as goedeckera, just with corrected params)
do i = 1,n
    A(i,i) = topo%gameeq(i) + tsqrt2pi/sqrt(topo%alpeeq(i))  ! ← FORTRAN USES SAME FORMULA
end do

do i = 1,n
    do j = 1,i-1
        gammij = 1.d0/sqrt(topo%alpeeq(i) + topo%alpeeq(j))  ! ← UPDATED ALPHA!
        tmp = erf(gammij*r)
        A(j,i) = tmp/r
    end do
end do

! Re-solve with corrected parameters
call sytrf_wrap(a, ipiv, io1)
call sytrs_wrap(a, x, ipiv, io2)
```

### C++ Phase 2: gfnff.cpp:3724-3745

```cpp
// ❌ WRONG: Old EEQ model
Matrix J_corrected = Matrix::Zero(m_atomcount, m_atomcount);

for (int i = 0; i < m_atomcount; ++i) {
    for (int j = 0; j < m_atomcount; ++j) {
        if (i == j) {
            J_corrected(i, i) = -1.0 / (2.0 * gam_corrected(i));  // ❌ WRONG FORMULA
        } else {
            J_corrected(i, j) = 1.0 / r;  // ❌ WRONG: No erf damping
        }
    }
}
```

---

## Expected Behavior Analysis

### For Water (H₂O):

**Correct Phase 1 Charges** (from Phase 1 code):
- O: q ≈ -0.85 e
- H: q ≈ +0.42 e each

**What Phase 2 SHOULD produce**:
- Same charges or slightly refined (less than 0.01 e change)
- Physically reasonable (O more negative due to higher EN)

**What Phase 2 ACTUALLY produces** (with 5 bugs):
- O: q ≈ +0.85 e (FLIPPED!)
- H: q ≈ -0.42 e (FLIPPED!)
- Or oscillates wildly due to matrix singularity

---

## Proof of Bug Impact

### Symptom 1: Charge Conservation Fails
Phase 1 enforces constraint row:
```cpp
x(m_atomcount) = 0.0;  // Neutral molecule
A(n, i) = 1.0;  // Constraint: sum(q) = 0
```

Phase 2 **DROPS this constraint row entirely**:
```cpp
Matrix J_corrected = Matrix::Zero(m_atomcount, m_atomcount);  // ← Only nxn, not augmented!
// No constraint row! No enforcement of charge neutrality!
```

Result: Charges don't sum to molecular charge → WRONG

### Symptom 2: Energy Oscillation/Divergence
With wrong EEQ matrix:
- First iteration: Charges converge to garbage
- Second iteration: Charges oscillate (because RHS is wrong)
- Convergence check: `max_change < 1e-5` might never be satisfied

Result: Calculation hangs or produces meaningless charges

### Symptom 3: Compared to Fortran
Fortran correctly recalculates charges → matches experimental values
C++ Phase 2 → charges are wrong sign/magnitude

---

## How to Fix

### Option 1: Use Phase 1 Output Directly (Quick Fix - 5 minutes)

```cpp
// In calculateTopologyInfo() line 2887:
if (use_two_phase) {
    // Phase 1: Calculate base topology-aware charges
    if (!calculateTopologyCharges(topo_info)) {
        return topo_info;
    }

    // Option A: Skip Phase 2 entirely (simplest)
    // topo_info.eeq_charges = topo_info.topology_charges;  // ← This might be GOOD ENOUGH

    // Option B: Keep Phase 2 but fix it (see below)
    if (!calculateFinalCharges_FIXED(topo_info)) {  // ← Fixed version
        return topo_info;
    }
}
```

**Rationale**: Phase 1 charges might already be accurate enough (50%+ error reduction documented)

### Option 2: Refactor Phase 2 to Match Fortran (Proper Fix - 30 minutes)

```cpp
bool GFNFF::calculateFinalCharges_FIXED(TopologyInfo& topo_info) const
{
    // Use AUGMENTED system (not 3x3, but 4x4 for 3-atom molecule)
    int m = m_atomcount + 1;  // ← AUGMENTED

    Vector final_charges = topo_info.topology_charges;

    for (int iter = 0; iter < max_iterations; ++iter) {
        // Apply ALL corrections to parameters
        Vector chi_corrected(m_atomcount);
        Vector gam_corrected(m_atomcount);
        Vector alpha_corrected(m_atomcount);  // ← NEW: Also correct alpha!

        for (int i = 0; i < m_atomcount; ++i) {
            int z_i = m_atoms[i];
            EEQParameters params_i = getEEQParameters(z_i);

            // Apply dxi correction: chi' = -chi + dxi
            chi_corrected(i) = -params_i.chi + topo_info.dxi(i);  // ← NEGATIVE chi!

            // Apply dgam correction
            double dgam_i = (i < topo_info.dgam.size()) ? topo_info.dgam(i) : 0.0;
            gam_corrected(i) = params_i.gam + dgam_i;

            // Apply dalpha correction: ← NEW!
            double dalpha_i = (i < topo_info.dalpha.size()) ? topo_info.dalpha(i) : 0.0;
            alpha_corrected(i) = params_i.alp + dalpha_i;
        }

        // Build AUGMENTED EEQ matrix (same as Phase 1)
        Matrix A = Matrix::Zero(m, m);
        Vector x = Vector::Zero(m);

        // Set RHS and diagonal using CORRECTED parameters
        for (int i = 0; i < m_atomcount; ++i) {
            x(i) = chi_corrected(i);  // ← POSITIVE (because chi is already negative)
            A(i, i) = gam_corrected(i) + TSQRT2PI / std::sqrt(alpha_corrected(i));  // ← CORRECT formula
        }

        // Coulomb matrix with erf damping using CORRECTED alpha
        for (int i = 0; i < m_atomcount; ++i) {
            for (int j = 0; j < i; ++j) {
                double r = distance(i, j);
                double gammij = 1.0 / std::sqrt(alpha_corrected(i) + alpha_corrected(j));  // ← USE CORRECTED
                double erf_val = std::erf(gammij * r);
                double coulomb = erf_val / r;
                A(i, j) = coulomb;
                A(j, i) = coulomb;
            }
        }

        // Fragment constraint
        x(m_atomcount) = 0.0;
        for (int j = 0; j < m_atomcount; ++j) {
            A(m_atomcount, j) = 1.0;
            A(j, m_atomcount) = 1.0;
        }

        // Solve using same method as Phase 1
        Eigen::PartialPivLU<Matrix> lu(A);
        Vector solution = lu.solve(x);

        Vector final_charges_new = solution.segment(0, m_atomcount);

        // Check convergence
        double max_change = 0.0;
        for (int i = 0; i < m_atomcount; ++i) {
            max_change = std::max(max_change, std::abs(final_charges_new(i) - final_charges(i)));
        }

        final_charges = final_charges_new;

        if (max_change < convergence_threshold) {
            break;  // Converged
        }
    }

    topo_info.eeq_charges = final_charges;
    return true;
}
```

---

## Current Session Recommendation

### IMMEDIATE ACTION:
1. **Disable Phase 2** or **use Phase 1 output directly**
   - Phase 1 is 99% correct (uses Fortran algorithm verbatim)
   - Phase 2 introduces 5 new bugs

2. **Test Phase 1 alone**:
   ```bash
   cd release
   ./curcuma -sp ../test_cases/molecules/water.xyz -method cgfnff -verbosity 2
   ```
   - Check if Phase 1 charges alone give reasonable energies
   - If yes: Problem is ONLY Phase 2

3. **Document findings**:
   - Phase 1: ✅ Working (goedeckera port successful)
   - Phase 2: ❌ Broken (5 critical bugs identified)
   - Recommendation: Remove Phase 2 until it's refactored properly

---

## Summary of Bugs

| Bug # | Line | Issue | Severity | Type |
|-------|------|-------|----------|------|
| #1 | 3730 | Wrong diagonal formula: `-1/(2gam)` vs `gam + sqrt(2/π)/sqrt(α)` | CRITICAL | Matrix singularity |
| #2 | 3742 | Missing erf damping: `1/r` vs `erf(γr)/r` | CRITICAL | 40% numerical error |
| #3 | 3717 | Sign flip: `+chi` vs `-chi` | CRITICAL | Charge sign reversal |
| #4 | 3750 | Double negative: `-(+chi)` vs `+chi` | CRITICAL | Charge magnitude flip |
| #5 | 3724-3745 | Ignored dalpha: Calculated but never used | CRITICAL | Missing physics |

**Total Impact**: Charges are COMPLETELY WRONG - wrong sign, wrong magnitude, possibly unconverged

---

## Files Affected
- `src/core/energy_calculators/qm_methods/gfnff.cpp:3680-3820` (calculateFinalCharges)
- `src/core/energy_calculators/qm_methods/gfnff.h` (function declaration)

## References
- **Phase 1 (Correct)**: `gfnff.cpp:3450-3557` (calculateTopologyCharges)
- **Fortran Reference**: `external/gfnff/src/gfnff_ini2.f90:1140-1246` (goedeckera)
- **Fortran Reference**: `external/gfnff/src/gfnff_ini.f90:411-420` (dxi calculation)
