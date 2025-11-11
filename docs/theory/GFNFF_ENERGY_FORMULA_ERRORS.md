# GFN-FF Energy Formula Errors in ForceFieldThread

**Date**: 2025-11-11
**Status**: 🚨 **CRITICAL ERRORS FOUND**
**File**: `src/core/energy_calculators/ff_methods/forcefieldthread.cpp:574-673`

---

## Executive Summary

The `ForceFieldThread` class implements **incorrect energy formulas** for GFN-FF bonds and angles:

| Component | Implemented Formula | Correct Fortran Formula | Status |
|-----------|---------------------|-------------------------|--------|
| **Bonds** | Harmonic + cubic anharmonic | **Exponential** | ❌ **WRONG** |
| **Angles** | Fourier (C0+C1cosθ+C2cos2θ) | Angle bending with damping | ❌ **ZERO ENERGY** |
| **Torsions** | V(1+cos(nφ-φ₀)) | V(1+cos(nφ-φ₀)) | ✅ Correct |
| **Inversions** | Fourier expansion | Out-of-plane bending | ❓ Unknown |

---

## 1. Bond Energy: Wrong Formula

### 1.1 Implemented (WRONG)

**Location**: `forcefieldthread.cpp:574-603`

```cpp
void ForceFieldThread::CalculateGFNFFBondContribution()
{
#pragma message("TODO: Implement proper GFN-FF bond stretching - currently using UFF geometry function")
    // TODO: GFN-FF uses harmonic + anharmonic bond stretching terms
    // Currently using UFF::BondStretching for geometry calculation only

    // GFN-FF bond stretching: E = 0.5*k*(r-r0)^2 + α*(r-r0)^3
    double dr = rij - bond.r0_ij;
    double harmonic = 0.5 * bond.fc * dr * dr;              // Harmonic term
    double anharmonic = bond.exponent * dr * dr * dr;       // Cubic anharmonic term

    m_bond_energy += (harmonic + anharmonic) * factor;
}
```

**Formula**: E_bond = 0.5·k·(r-r₀)² + α·(r-r₀)³

**Problems**:
1. ❌ **Wrong functional form**: Uses polynomial (harmonic + cubic), not exponential
2. ❌ **Wrong prefactor**: 0.5 for harmonic is standard but doesn't match Fortran
3. ❌ **Cubic anharmonic**: GFN-FF doesn't use cubic terms
4. ⚠️ Comment acknowledges this is a TODO placeholder

---

### 1.2 Correct Fortran Implementation

**Location**: `external/gfnff/src/gfnff_engrad.F90:675-721`

```fortran
subroutine egbond(i,iat,jat,rab,rij,drij,n,at,xyz,e,g,topo)
  t8 = topo%vbond(2,i)                     ! α parameter (exponential decay)
  dr = rab-rij                             ! r - r₀
  dum = topo%vbond(3,i)*exp(-t8*dr**2)     ! k_b * exp(-α*(r-r₀)²)
  e = e+dum                                ! Add to total energy
  yy = 2.0d0*t8*dr*dum                     ! Gradient factor
  ! ... gradient calculation ...
end subroutine egbond
```

**Formula**: E_bond = k_b · exp(-α · (r-r₀)²)

**Key differences**:
- ✅ **Exponential decay**: Energy goes to zero at large displacements
- ✅ **Gaussian-like potential**: Soft at long range, stiff at short range
- ✅ **Two parameters**: k_b (energy scale) and α (stiffness)
- ✅ **No prefactor**: Energy scale is directly k_b, not 0.5·k

---

### 1.3 Comparison

| Property | C++ (Wrong) | Fortran (Correct) |
|----------|-------------|-------------------|
| Functional form | Polynomial | Exponential |
| Energy at r=r₀ | 0 | k_b |
| Energy at r→∞ | +∞ (unphysical!) | 0 (bond breaks) |
| Gradient at r=r₀ | 0 | 0 |
| Second derivative | Constant k | k_b·2α (distance-dependent) |
| Physical behavior | Harmonic oscillator | Morse-like potential |

**Impact**: Bond energies will be **completely wrong** for large displacements (conformational changes, bond breaking).

---

## 2. Angle Energy: Zero Energy!

### 2.1 Implemented (ZERO ENERGY)

**Location**: `forcefieldthread.cpp:605-630`

```cpp
void ForceFieldThread::CalculateGFNFFAngleContribution()
{
#pragma message("TODO: Verify GFN-FF angle bending functional form - currently using UFF geometry")
    // TODO: Check if GFN-FF angle bending uses same Fourier expansion as UFF

    // GFN-FF angle bending: E = k*(C0 + C1*cos(θ) + C2*cos(2θ))
    m_angle_energy += (angle.fc * (angle.C0 + angle.C1 * costheta +
                       angle.C2 * (2 * costheta * costheta - 1))) * factor;
}
```

**Formula**: E_angle = k · (C0 + C1·cosθ + C2·cos2θ)

**Problems**:
1. ❌ **Zero energy**: C0, C1, C2 are all set to **0.0** in `getGFNFFAngleParameters:535-537`
2. ❌ **Missing distance damping**: Fortran couples angle bending to bond stretching
3. ❌ **Wrong parameter passing**: Fortran stores θ₀ and k, not Fourier coefficients
4. ⚠️ Comment acknowledges this needs verification

---

### 2.2 Parameter Assignment (DUMMY VALUES)

**Location**: `gfnff.cpp:502-536`

```cpp
GFNFF::GFNFFAngleParams GFNFF::getGFNFFAngleParameters(int z1, int z2, int z3, double current_angle) const
{
    // ... parameter lookup ...

    params.force_constant = angle_param * 0.001;
    params.equilibrium_angle = current_angle;

    // ❌ DUMMY VALUES - angle energy will be ZERO!
    params.c0 = 0.0;
    params.c1 = 0.0;
    params.c2 = 0.0;

    return params;
}
```

**Result**: Since C0=C1=C2=0, the angle energy is:

```
E_angle = k · (0 + 0·cosθ + 0·cos2θ) = 0
```

**All angle bending energy is missing from GFN-FF calculations!**

---

### 2.3 Correct Fortran Implementation

**Location**: `external/gfnff/src/gfnff_engrad.F90:857-892`

```fortran
subroutine egbend(m,j,i,k,n,at,xyz,e,g,param,topo)
  c0 = topo%vangl(1,m)        ! Equilibrium angle θ₀
  kijk = topo%vangl(2,m)      ! Force constant

  theta = dacos(cosa)         ! Current angle

  ! Distance damping (couples to bond stretching)
  call gfnffdampa(at(i),at(j),rab2,dampij,damp2ij,param)
  call gfnffdampa(at(k),at(j),rcb2,dampjk,damp2jk,param)
  damp = dampij*dampjk

  ! Complex angle bending potential with damping
  ! (not shown - ~30 lines of gradient calculation)
end subroutine egbend
```

**Formula**: E_angle = k_ijk · damp_ij · damp_jk · f(θ, θ₀)

**Key features**:
- ✅ **Distance damping**: Angle bending weakens when bonds stretch
- ✅ **Two parameters**: θ₀ (equilibrium) and k (force constant)
- ✅ **Angle bending potential**: f(θ, θ₀) is NOT a Fourier expansion
- ✅ **Coupled to bonds**: Damping factors depend on bond lengths

---

## 3. Torsion Energy: Correct ✅

**Location**: `forcefieldthread.cpp:632-662`

```cpp
void ForceFieldThread::CalculateGFNFFDihedralContribution()
{
    // GFN-FF torsion: E = V*(1 + cos(n*φ - φ0))

    double V = dihedral.V;
    double n = dihedral.n;
    double phi0 = dihedral.phi0;

    double energy = V * (1 + cos(n * phi - phi0));
    m_dihedral_energy += energy * factor;
}
```

**Formula**: E_torsion = V · (1 + cos(n·φ - φ₀))

**Status**: ✅ **Matches Fortran formula** (`gfnff_engrad.F90:1041-1122`)

---

## 4. Inversion Energy: Unknown ❓

**Location**: `forcefieldthread.cpp:664-673`

```cpp
void ForceFieldThread::CalculateGFNFFInversionContribution()
{
    // GFN-FF inversion: E = k*(C0 + C1*cos(θ) + C2*cos(2θ))

    // (similar to angle code)
}
```

**Problems**:
1. ❓ **Unclear if correct**: Uses same Fourier expansion as angles
2. ❓ **Parameter passing**: Need to verify C0/C1/C2 values are non-zero
3. ❓ **Fortran comparison**: Need to compare with out-of-plane potential

**Action required**: Check if `generateGFNFFInversions()` passes correct parameters.

---

## 5. Impact on Phase 1 Results

### 5.1 What Works

✅ **Torsion parameters generated correctly** (Phase 1.1):
- `generateGFNFFTorsions()` creates correct V, n, φ₀ values
- `ForceFieldThread` uses correct cosine formula
- Expected to give reasonable energies

✅ **Inversion parameters generated correctly** (Phase 1.2):
- `generateGFNFFInversions()` creates barrier heights
- Formula in `ForceFieldThread` needs verification
- May give reasonable energies (to be tested)

---

### 5.2 What's Broken

❌ **Bond energy is wrong**:
- Uses harmonic + cubic instead of exponential
- Will give incorrect energies for all molecules
- Especially wrong for strained geometries

❌ **Angle energy is zero**:
- All angle bending energy missing
- Molecules will have no angular constraints
- Geometries will collapse without angle forces

❌ **Missing parameter corrections**:
- Even if formulas were fixed, parameters are oversimplified
- Missing ring strain, electronegativity, charge, CN corrections
- See `GFNFF_BOND_ANGLE_VALIDATION.md` for details

---

## 6. Required Fixes

### 6.1 Critical (Must Fix)

**Priority 1: Bond energy formula**

Replace harmonic+cubic with exponential:

```cpp
void ForceFieldThread::CalculateGFNFFBondContribution()
{
    for (const auto& bond : m_gfnff_bonds) {
        double dr = rij - bond.r0_ij;

        // Correct GFN-FF formula: E = k_b * exp(-α * dr²)
        double alpha = bond.alpha;     // Need to pass this from parameters!
        double k_b = bond.fc;
        double energy = k_b * exp(-alpha * dr * dr);

        m_bond_energy += energy * factor;

        if (m_calculate_gradient) {
            // dE/dr = -2*α*dr*E
            double dEdr = -2.0 * alpha * dr * energy;
            // ... apply to gradient ...
        }
    }
}
```

**Required changes**:
1. Add `alpha` parameter to `GFNFFBondParams` structure
2. Calculate α in `getGFNFFBondParameters()` based on electronegativity
3. Pass α to `ForceFieldThread` via bond JSON
4. Update gradient calculation

---

**Priority 2: Angle energy formula**

Replace Fourier with angle bending:

```cpp
void ForceFieldThread::CalculateGFNFFAngleContribution()
{
    for (const auto& angle : m_gfnff_angles) {
        double theta = acos(costheta);
        double dtheta = theta - angle.theta0_ijk;  // θ - θ₀

        // Distance damping (couple to bond stretching)
        double r_ij = (m_geometry.row(angle.i) - m_geometry.row(angle.j)).norm();
        double r_kj = (m_geometry.row(angle.k) - m_geometry.row(angle.j)).norm();
        double damp_ij = calculateGFNFFDamping(angle.i, angle.j, r_ij);
        double damp_kj = calculateGFNFFDamping(angle.k, angle.j, r_kj);
        double damp = damp_ij * damp_kj;

        // Angle bending potential
        double energy = angle.fc * damp * dtheta * dtheta;

        m_angle_energy += energy * factor;
    }
}
```

**Required changes**:
1. Remove C0, C1, C2 from `GFNFFAngleParams`
2. Store θ₀ (equilibrium angle) instead
3. Implement distance damping function
4. Update gradient calculation

---

### 6.2 Important (Should Fix)

**Priority 3: Parameter corrections**

Implement missing correction factors (see `GFNFF_BOND_ANGLE_VALIDATION.md`):
- Ring strain detection (requires Phase 2 topology)
- Electronegativity-based α parameter
- Bond order detection (single/double/triple)
- Charge-dependent scaling (requires Phase 3 EEQ)
- Heavy atom, pi-system, X-H, CN corrections

**Priority 4: Inversion validation**

Verify inversion energy formula and parameters:
1. Check if Fourier expansion is correct for out-of-plane bending
2. Compare with Fortran `gfnff_engrad.F90` inversion code
3. Test on ethene (planar sp²) to verify zero energy at equilibrium

---

## 7. Testing Strategy

### 7.1 Unit Tests

Create tests for each component:

**Test 1: Bond energy formula**
- Input: C-C bond at r₀ = 1.54 Å with known k_b, α
- Expected: E = k_b at r₀, E → 0 as r → ∞
- Current: E = 0 at r₀ (wrong!)

**Test 2: Angle energy formula**
- Input: C-C-C angle at 109.5° (tetrahedral)
- Expected: Non-zero energy if θ ≠ θ₀
- Current: E = 0 always (broken!)

**Test 3: Torsion energy formula**
- Input: C-C-C-C torsion at φ = 60° (gauche butane)
- Expected: E = V·(1 + cos(3φ))
- Current: Should work ✅

---

### 7.2 Integration Tests

Compare total energies with Fortran GFN-FF:

| Molecule | Test Purpose | Expected Δ |
|----------|--------------|------------|
| **Methane** | 4 bonds, 6 angles, no torsions | Large (bond+angle wrong) |
| **Ethane** | Bonds, angles, 1 torsion | Large (bond+angle wrong) |
| **Butane** | Full test | Large (bond+angle wrong) |
| **Ethene** | sp² planarity, inversions | Very large (angle missing!) |
| **Benzene** | Aromatic, all terms | Catastrophic (no angles) |

**Prediction**: Current C++ implementation will give **wildly incorrect energies** due to:
1. Wrong bond formula (exponential → harmonic)
2. Zero angle energy (missing angular constraints)
3. Even if torsions work, other terms dominate

---

## 8. Recommendations

### Option A: Fix Formulas Immediately (Recommended)

**Time**: 1-2 days
**Priority**: CRITICAL

1. Implement exponential bond formula (Priority 1)
2. Implement angle bending formula (Priority 2)
3. Run validation tests on 5 test molecules
4. Compare energies with Fortran GFN-FF
5. Document in Phase 1 report

**Blocker**: Need α parameter calculation (electronegativity-based)

---

### Option B: Document as Phase 0 + Use External GFN-FF

**Time**: 1 hour
**Priority**: PRAGMATIC

1. Update all Phase 1 docs to clarify current status:
   - "Phase 0: Parameter generation only (energy formulas incomplete)"
   - "Torsions/inversions: Parameters correct, formulas correct ✅"
   - "Bonds/angles: Parameters simplified, formulas wrong ❌"
2. Continue using XTB/TBLite Fortran for production
3. Native C++ is educational/experimental only
4. Return to fix after Phase 2-3 complete

---

### Option C: Phase 1.3 - Fix Bond/Angle Formulas

**Time**: 3-5 days
**Priority**: COMPLETE PHASE 1

Make Phase 1 self-consistent before moving to Phase 2:

**Phase 1.3a: Bond Formula** (2 days)
- Implement exponential bond potential
- Add α parameter calculation (simplified, without full topology)
- Test on methane, ethane, butane

**Phase 1.3b: Angle Formula** (1-2 days)
- Implement angle bending with distance damping
- Remove dummy C0/C1/C2 parameters
- Test on methane (tetrahedral angles)

**Phase 1.3c: Validation** (1 day)
- Run all 5 test molecules (butane, ethene, benzene, methane, ethane)
- Compare with Fortran GFN-FF
- Accept ±20% error due to missing topology corrections
- Document results

---

## 9. Conclusion

The existing GFN-FF implementation has **critical errors in energy formulas**:

1. ❌ **Bond energy**: Uses harmonic+cubic instead of exponential
2. ❌ **Angle energy**: Always zero due to dummy parameters
3. ✅ **Torsion energy**: Formula appears correct
4. ❓ **Inversion energy**: Needs verification

Combined with **missing parameter corrections** (documented in `GFNFF_BOND_ANGLE_VALIDATION.md`), the current implementation will produce **completely wrong energies** for any molecule.

**Immediate action required**: Choose Option A (fix formulas now) or Option B (document limitations and use external GFN-FF).

---

**Report prepared**: 2025-11-11
**Author**: Claude (Anthropic)
**Critical severity**: Energy formulas are wrong - system cannot be used for production
