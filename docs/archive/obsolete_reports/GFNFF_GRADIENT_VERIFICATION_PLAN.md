# GFN-FF Gradient Verification Plan

**Status**: Phase 1 - Documentation Complete
**Date**: February 1, 2026
**Goal**: Verify and validate all GFN-FF gradient implementations against Fortran reference

---

## Phase 1: Current Implementation Status

### Summary

| Energy Term | Curcuma Method | Fortran Subroutine | Status in execute() | Gradient Implemented |
|-------------|----------------|-------------------|---------------------|---------------------|
| **Bonds** | CalculateGFNFFBondContribution() | egbond (F90:675-721) | ❌ COMMENTED OUT | ✅ Yes |
| **Angles** | CalculateGFNFFAngleContribution() | egbend (F90:857-916) | ❌ COMMENTED OUT | ✅ Yes |
| **Torsions** | CalculateGFNFFDihedralContribution() | egtors (F90:1153-1234) | ❌ COMMENTED OUT | ✅ Yes |
| **Extra Torsions** | CalculateGFNFFExtraTorsionContribution() | egtors (F90:1272-1280) | ❌ COMMENTED OUT | ✅ Yes |
| **Inversions** | CalculateGFNFFInversionContribution() | (in gfnff_ini.f90) | ❌ COMMENTED OUT | ✅ Yes |
| **Dispersion** | CalculateGFNFFDispersionContribution() | gdisp0 (gdisp0.f90) | ❌ COMMENTED OUT | ✅ Yes |
| **Repulsion (bonded)** | CalculateGFNFFBondedRepulsionContribution() | engrad.F90:467-495 | ❌ COMMENTED OUT | ⚠️ Partial |
| **Repulsion (non-bonded)** | CalculateGFNFFNonbondedRepulsionContribution() | engrad.F90:255-276 | ✅ ACTIVE | ⚠️ Partial |
| **Coulomb** | CalculateGFNFFCoulombContribution() | engrad.F90:383-422 | ❌ COMMENTED OUT | ❌ No |
| **Hydrogen Bonds** | CalculateGFNFFHydrogenBondContribution() | abhgfnff_eg1/eg2/eg3 | ✅ ACTIVE | ⚠️ Partial |
| **Halogen Bonds** | CalculateGFNFFHalogenBondContribution() | rbxgfnff_eg | ✅ ACTIVE | ⚠️ Partial |
| **BATM** | CalculateGFNFFBatmContribution() | batmgfnff_eg | ✅ ACTIVE | ⚠️ Energy only |
| **D3/D4 ATM** | CalculateATMContribution() + Gradient() | d3_gradient | ✅ ACTIVE | ✅ Yes |

### Critical Finding: Most GFN-FF Gradients Disabled

In `forcefieldthread.cpp:execute()` (lines 96-180), **most GFN-FF gradient calculations are commented out**:

```cpp
// Lines 116-120: ALL bonded terms disabled
// CalculateGFNFFBondContribution();
// CalculateGFNFFAngleContribution();
// CalculateGFNFFDihedralContribution();
// CalculateGFNFFExtraTorsionContribution();
// CalculateGFNFFInversionContribution();

// Lines 123-131: ALL non-bonded terms disabled
// if (m_dispersion_enabled) {
//     CalculateGFNFFDispersionContribution();
// }
// if (m_repulsion_enabled) {
//     CalculateGFNFFBondedRepulsionContribution();
//     CalculateGFNFFNonbondedRepulsionContribution();  // ONLY THIS IS ACTIVE
// }
// if (m_coulomb_enabled) {
//     CalculateGFNFFCoulombContribution();
// }

// Lines 135-137: HB/XB terms ACTIVE
if (m_hbond_enabled) {
    CalculateGFNFFHydrogenBondContribution();  // ACTIVE
    CalculateGFNFFHalogenBondContribution();   // ACTIVE
}
```

### Active Terms (Currently Calculated)

1. **Non-bonded Repulsion** (line 128) - Partial implementation
2. **Hydrogen Bonds** (line 136) - Partial implementation
3. **Halogen Bonds** (line 137) - Partial implementation
4. **D3/D4 Dispersion** (lines 147, 154) - Native D3/D4 implementation
5. **ATM three-body** (line 162+166) - Complete with gradients
6. **BATM** (line 176) - Energy only, gradients TODO

### Detailed Gradient Implementation Analysis

#### 1. Bond Gradients (CalculateGFNFFBondContribution, line 756)

**Status**: Implemented but DISABLED

**Formula** (from Fortran gfnff_engrad.F90:675-721):
```fortran
! E_bond = k_b * exp(-α * (r - r₀)²)
! g = dE/dr * dr/dx = -2*α*(r-r₀)*E * (dx/r)
```

**Curcuma Implementation** (lines 834-842):
```cpp
if (m_calculate_gradient) {
    double dEdr = -2.0 * alpha * dr * energy;
    m_gradient.row(bond.i) += dEdr * factor * derivate.row(0);
    m_gradient.row(bond.j) += dEdr * factor * derivate.row(1);
}
```

**Issues**:
- Missing CN gradient contribution (dr0/dCN * dCN/dx)
- Missing HB modulation gradient (dhb_cn/dx for egbond_hb)

#### 2. Angle Gradients (CalculateGFNFFAngleContribution, line 846)

**Status**: Implemented but DISABLED

**Formula** (Fortran gfnff_engrad.F90:857-916):
```fortran
! E = k * (cosθ - cosθ₀)² * damp_ij * damp_jk
! g = ∂E/∂θ * ∂θ/∂x * damp + E * ∂damp/∂x
```

**Curcuma Implementation** (lines 993-1034):
- Complete with damping gradients
- Distance-dependent damping included
- Fortran formula matches exactly

#### 3. Torsion Gradients (CalculateGFNFFDihedralContribution, line 1038)

**Status**: Implemented but DISABLED

**Formula** (Fortran gfnff_engrad.F90:1153-1234):
```fortran
! E = V * (1 + cos(n*φ - φ₀)) * damp_ij * damp_jk * damp_kl
! g = ∂E/∂φ * ∂φ/∂x * damp + E * ∂damp/∂x
```

**Curcuma Implementation** (lines 1038-1360):
- Complete with damping gradients
- Cross-center damping formula matches Fortran
- NCI torsion support added (Jan 13, 2026)

**Known Issue**: Energy scaling factor still under investigation (see TORSION_DEBUGGING_SUMMARY.md)

#### 4. Extra Torsion Gradients (CalculateGFNFFExtraTorsionContribution, line 1366)

**Status**: Implemented but DISABLED

**Formula**: Same as primary torsions but without +π phase shift

**Curcuma Implementation**: Complete with full gradient support

#### 5. Repulsion Gradients

##### Non-bonded Repulsion (line 1748) - ACTIVE

**Formula** (Fortran gfnff_engrad.F90:255-276):
```fortran
! E = exp(-α*r^1.5) * repz(i) * repz(j) * repscaln / r
! g = E * (1.5*α*r^1.5 + 1) / r² * dx
```

**Curcuma Implementation**: Check if gradients match Fortran exactly

##### Bonded Repulsion (line 1681) - DISABLED

**Formula** (Fortran gfnff_engrad.F90:467-495):
- Similar to non-bonded but with bonded parameters
- Uses `repscalb` instead of `repscaln`

**Status**: Implementation exists but not validated

#### 6. Coulomb Gradients (CalculateGFNFFCoulombContribution, line 1815)

**Status**: ❌ NOT IMPLEMENTED

**Formula** (Fortran gfnff_engrad.F90:383-422):
```fortran
! E = q_i * q_j * erf(γ_ij * r) / r
! g = ∂E/∂r * ∂r/∂x + ∂E/∂q * ∂q/∂x (chain rule through EEQ)
```

**Curcuma Implementation**: Energy only, gradients missing

**Challenge**: Coulomb gradients require EEQ charge derivatives (dq/dx), which involves solving the EEQ matrix equation for geometry derivatives.

#### 7. HB/XB Gradients (lines 2088, 2356) - ACTIVE

**Status**: Partial implementation

**Formula** (Fortran gfnff_engrad.F90: various abhgfnff_eg subroutines):
- Three-body terms with angle and distance dependence
- Complex damping functions

**Curcuma Implementation**: Needs verification against Fortran

#### 8. BATM Gradients (CalculateGFNFFBatmContribution, line 2955)

**Status**: ❌ NOT IMPLEMENTED (energy only)

**Formula** (Fortran gfnff_engrad.F90:566-603):
```fortran
! batmgfnff_eg subroutine
! Three-body dispersion-like term for 1,4-pairs
```

**Curcuma Implementation**: Energy calculation complete, gradients TODO

---

## Phase 2: Test Framework Requirements

### Test Program Structure

```cpp
// test_cases/test_gfnff_gradients.cpp

class GFNFFGradientTest {
    // 1. Finite-difference gradient calculator
    Vector calculateFiniteDifferenceGradient(Molecule& mol, int atom, int coord);

    // 2. Term-specific gradient tests
    void testBondGradients();
    void testAngleGradients();
    void testTorsionGradients();
    void testRepulsionGradients();
    void testCoulombGradients();
    void testDispersionGradients();
    void testHBondGradients();
    void testXBondGradients();
    void testBatmGradients();

    // 3. Comparison utilities
    double compareGradients(const Vector& analytical, const Vector& numerical);
    bool isWithinTolerance(double error, double tolerance = 1e-6);
};
```

### Test Molecules

1. **H2O** (3 atoms) - Test bonds, angles
2. **CH4** (5 atoms) - Test bonds, angles, non-bonded
3. **CH3OCH3** (9 atoms) - Test all terms including torsions
4. **C2H6** (8 atoms) - Test torsion barriers
5. **H2O dimer** (6 atoms) - Test HB terms

### Success Criteria

- Analytical vs numerical gradient RMS error < 1e-6 for each term
- Maximum component error < 1e-5
- Energy conservation check: |ΔE - g·Δx| < 1e-8

---

## Phase 3: Implementation Priority

### Critical (Blocks Production Use)

1. **Enable all gradient methods** in execute()
2. **Implement Coulomb gradients** (requires EEQ derivatives)
3. **Verify repulsion gradients** against Fortran

### High Priority

4. **BATM gradients** (energy-only currently)
5. **HB/XB gradient validation**
6. **Bond/Angle/Torsion validation** with finite differences

### Medium Priority

7. **CN gradient contributions** (for dynamic r0)
8. **EEQ charge derivatives** (for complete Coulomb gradients)

---

## Fortran Reference Mapping

### Key Fortran Files

| Fortran File | Description | Key Subroutines |
|--------------|-------------|-----------------|
| gfnff_engrad.F90 | Main energy/gradient | gfnff_eg, egbond, egbend, egtors |
| gfnff_gdisp0.f90 | Dispersion | gdisp0, d3_gradient |
| gfnff_rab.f90 | Bond distance functions | gfnffdrab (with gradients) |
| gfnff_hb.f90 | Hydrogen bonding | abhgfnff_eg1, abhgfnff_eg2, abhgfnff_eg3 |
| gfnff_xb.f90 | Halogen bonding | rbxgfnff_eg |
| gfnff_batm.f90 | Bonded ATM | batmgfnff_eg |

### Line Number Reference (gfnff_engrad.F90)

| Term | Energy Lines | Gradient Lines |
|------|--------------|----------------|
| Non-bonded Repulsion | 255-269 | 270-276 |
| EEQ Setup | 330-422 | 383-422 |
| Bonded Repulsion | 467-486 | 487-495 |
| Bonds | 440-453 | (via gfnffdrab) |
| Angles | 508-516 | 512-516 (egbend) |
| Torsions | 525-541 | 534-538 (egtors) |
| BATM | 566-578 | 574-578 (batmgfnff_eg) |
| HB (type 1) | 612-627 | 620-623 (abhgfnff_eg1) |
| HB (type 2) | 629-657 | 651-654 (abhgfnff_eg2/eg3) |
| XB | 687-700 | 696-699 (rbxgfnff_eg) |

---

## Next Steps

1. **Phase 2**: Create test_gfnff_gradients.cpp with finite-difference framework
2. **Phase 3**: Enable and validate each gradient term systematically
3. **Phase 4**: Implement missing Coulomb and BATM gradients
4. **Phase 5**: Full regression test suite
