# GFN-FF Fortran Functions Reference
## Mapping External Library to Native C++ Implementation

**Purpose**: Quick reference for porting Fortran functions to C++

---

## **Core Energy Calculation Functions**

### 1. Main Energy/Gradient Routine
**File**: `external/gfnff/src/gfnff_engrad.F90`

| Fortran Function | Line | Purpose | C++ Target |
|------------------|------|---------|------------|
| `gfnff_eg()` | 92-671 | Main energy+gradient calculation | `GFNFF::Calculation()` |
| `egbond()` | 675-721 | Bond stretching energy | `GFNFF::calculateBondEnergy()` ✅ |
| `egbond_hb()` | 725-800 | H-bond corrected bonds | `GFNFF::calculateHBondEnergy()` ❌ |
| `egbend()` | 857-917 | Angle bending energy | `GFNFF::calculateAngleEnergy()` ✅ |
| `egbend_nci()` | 978-1037 | NCI-corrected angles | `GFNFF::calculateAngleEnergy_NCI()` ❌ |
| `egtors()` | 1041-1122 | Torsion energy | `GFNFF::calculateTorsionEnergy()` ❌ |
| `egtors_nci()` | 1170-1216 | NCI-corrected torsions | `GFNFF::calculateTorsionEnergy_NCI()` ❌ |
| `goed_gfnff()` | 1274-1400 | EEQ charge calculation | `GFNFF::calculateEEQCharges()` ⚠️ placeholder |
| `dncoord_erf()` | 802-853 | CN derivatives | `GFNFF::calculateCoordinationNumberDerivatives()` ⚠️ zeros |

---

## **Initialization & Topology**

### 2. Main Initialization
**File**: `external/gfnff/src/gfnff_ini.f90` (84,797 lines!)

| Fortran Function | Line | Purpose | C++ Target |
|------------------|------|---------|------------|
| `gfnff_ini()` | 38-2500 | Full topology setup | `GFNFF::initializeForceField()` ⚠️ simplified |
| `gfnff_dlogcoord()` | `gfnff_cn.f90` | CN calculation | `GFNFF::calculateCoordinationNumbers()` ⚠️ basic |
| **Ring Detection** | ~600-900 | Find rings 3-20 | `GFNFF::findSmallestRings()` ❌ returns zeros |
| **Pi-System Setup** | ~1100-1300 | Conjugation detection | `GFNFF::detectPiSystems()` ❌ empty |
| **Hybridization** | ~400-550 | Geometry-based hyb | `GFNFF::determineHybridization()` ⚠️ neighbor count |
| **Torsion Setup** | `gfnff_ini2.f90` | Generate torsion list | `GFNFF::generateGFNFFTorsions()` ❌ |
| **OOP Setup** | `gfnff_ini2.f90` | Generate inversions | `GFNFF::generateGFNFFInversions()` ❌ |

---

## **Parameter Database**

### 3. Parameter Arrays
**File**: `external/gfnff/src/gfnff_param.f90`

| Parameter | Fortran Variable | Line | C++ Location | Status |
|-----------|------------------|------|--------------|--------|
| **Electronegativities** | `chi_angewChem2020(86)` | 87-105 | `gfnff.cpp:814-836` | ✅ Copied |
| **Hardness** | `gam_angewChem2020(86)` | 107-125 | `gfnff.cpp:838-860` | ✅ Copied |
| **CN Corrections** | `cnf_angewChem2020(86)` | 127-145 | Not in C++ yet | ❌ |
| **Polarizabilities** | `alp_angewChem2020(86)` | 147-165 | Hardcoded `5.0` | ❌ |
| **Covalent Radii** | `rad(86)` | Separate file | `gfnff.cpp:430-443` | ✅ Copied |
| **Bond Parameters** | `bond_angewChem2020(86)` | Not visible | `gfnff.cpp:461-482` | ✅ Copied |
| **Angle Parameters** | `angl_angewChem2020(86)` | Not visible | `gfnff.cpp:505-527` | ✅ Copied |

**Note**: Many parameter arrays are in **generated include files** not visible in source

---

## **Helper Functions**

### 4. Geometry & Topology Helpers
**File**: `external/gfnff/src/gfnff_helpers.f90`

| Fortran Function | Purpose | C++ Equivalent |
|------------------|---------|----------------|
| `pbc_gfnff()` | Periodic boundary handling | N/A (no PBC yet) |
| `gfnff_angle()` | Calculate bond angle | Eigen vector math ✅ |
| `gfnff_dihedral()` | Calculate dihedral angle | Need to implement ❌ |
| `gfnff_outofplane()` | Calculate OOP angle | Need to implement ❌ |

---

## **Priority Porting List**

### **Phase 1: Critical Energy Terms** ⚡
1. ✅ **Port Now**: `egtors()` → Torsion energy/gradient
   - Input: 4 atom indices (i,j,k,l), geometry
   - Output: Energy contribution, gradient updates
   - Complexity: Medium (dihedral angle calculation)
   - Lines: ~80 in Fortran

2. ✅ **Port Now**: Torsion setup from `gfnff_ini2.f90`
   - Detect all i-j-k-l sequences from bond list
   - Assign periodicity based on hybridization
   - Lookup barrier from parameter tables

3. ✅ **Port Now**: Out-of-plane energy calculation
   - Input: Central atom + 3 neighbors
   - Output: Inversion energy, gradient
   - Complexity: Low (Wilson angle)
   - Lines: ~50 estimated

### **Phase 2: Topology** 🔍
4. 📝 **Port Next**: Ring detection algorithm
   - Complexity: High (DFS/BFS graph algorithms)
   - Estimate: 200-300 lines C++
   - Critical for: Small ring strain corrections

5. 📝 **Port Next**: Pi-system detection
   - Complexity: Medium (graph connectivity)
   - Estimate: 100-150 lines C++
   - Critical for: Aromatic/conjugated parameter corrections

6. 📝 **Port Next**: Hybridization from geometry
   - Complexity: Medium (angle analysis)
   - Estimate: 150 lines C++
   - Critical for: Accurate parameter assignment

### **Phase 3: Charges** ⚡
7. 🔥 **Complex**: `goed_gfnff()` EEQ solver
   - Complexity: Very High (linear algebra)
   - Estimate: 300-400 lines C++
   - Critical for: Accurate electrostatics
   - **Strategy**: Port incrementally
     1. Matrix setup (without CN corrections)
     2. Linear solver (Eigen LU decomposition)
     3. Add CN corrections
     4. Add fragment constraints

8. 🔥 **Complex**: `dncoord_erf()` CN derivatives
   - Complexity: High (3D tensor derivatives)
   - Estimate: 150 lines C++
   - Critical for: Accurate gradients
   - **Strategy**: Implement after basic EEQ works

---

## **Code Porting Patterns**

### Pattern 1: Energy Term Translation
```fortran
! Fortran (gfnff_engrad.F90:1041)
subroutine egtors(m,i,j,k,l,n,at,xyz,e,g,param,topo)
    real(wp) :: phi, valijkl, rn, fc

    ! Calculate dihedral angle
    call gfnff_dihedral(xyz,i,j,k,l,phi)

    ! Energy: V_n/2 * [1 - cos(n*phi - phi0)]
    valijkl = fc * (1.0_wp - cos(rn*phi - phi0))
    e = e + valijkl

    ! Gradient: Chain rule dE/dphi * dphi/dxyz
    ! ... [gradient code]
end subroutine
```

**C++ Translation**:
```cpp
// C++ (gfnff.cpp - to be added)
double GFNFF::calculateTorsionEnergy(int i, int j, int k, int l,
                                      bool calc_gradient) {
    // 1. Calculate dihedral angle
    double phi = calculateDihedralAngle(i, j, k, l);

    // 2. Get parameters from topology
    auto params = getGFNFFTorsionParameters(i, j, k, l);
    double fc = params.barrier;
    int rn = params.periodicity;
    double phi0 = params.phase;

    // 3. Energy
    double energy = 0.5 * fc * (1.0 - cos(rn * phi - phi0));

    // 4. Gradient (if requested)
    if (calc_gradient) {
        double dE_dphi = 0.5 * fc * rn * sin(rn * phi - phi0);
        Matrix dphi_dxyz = calculateDihedralGradient(i, j, k, l);
        updateGradient(i, j, k, l, dE_dphi * dphi_dxyz);
    }

    return energy;
}
```

### Pattern 2: Parameter Lookup
```fortran
! Fortran: Direct array access with environment corrections
fc = bond_angewChem2020(ati) * bond_angewChem2020(atj)
fc = fc * (1.0 + cn_correction + ring_correction)
```

**C++ Translation**:
```cpp
// C++: Encapsulated with topology awareness
double fc = getBondForceConstant(ati, atj);

// In getBondForceConstant():
double base_fc = sqrt(bond_angewChem2020[ati-1] *
                      bond_angewChem2020[atj-1]);

// Apply topology corrections
TopologyInfo topo = getTopologyInfo(ati, atj);
if (topo.in_small_ring) {
    base_fc *= 1.2;  // Stiffer in small rings
}
if (topo.is_conjugated) {
    base_fc *= 0.9;  // Softer in conjugated systems
}

return base_fc;
```

### Pattern 3: Matrix Operations
```fortran
! Fortran: Manual loops
do i = 1, nat
    do j = 1, nat
        A(i,j) = gam(i) + coulomb(i,j)
    end do
end do
```

**C++ Translation**:
```cpp
// C++: Eigen vectorized operations
Matrix A(nat, nat);

// Diagonal
for (int i = 0; i < nat; ++i) {
    A(i,i) = gam[i];
}

// Off-diagonal (Eigen auto-vectorizes)
for (int i = 0; i < nat; ++i) {
    for (int j = i+1; j < nat; ++j) {
        double coulomb_ij = calculateCoulomb(i, j);
        A(i,j) = A(j,i) = coulomb_ij;
    }
}
```

---

## **Testing Strategy for Each Function**

### Test Template
```cpp
// test_cases/gfnff_native/test_torsions.cpp
TEST(GFNFFNative, TorsionEnergyButane) {
    // 1. Load molecule
    Molecule butane = loadXYZ("butane.xyz");

    // 2. Calculate with native
    GFNFF native_gfnff({{"method", "gfnff"}});
    native_gfnff.InitialiseMolecule(butane);
    double E_native = native_gfnff.Calculation(false);

    // 3. Calculate with external (reference)
    ExternalGFNFFMethod external_gfnff({});
    external_gfnff.setMolecule(butane);
    double E_reference = external_gfnff.calculateEnergy(false);

    // 4. Compare
    EXPECT_NEAR(E_native, E_reference, 0.5 / 627.5);  // 0.5 kcal/mol tolerance
}

TEST(GFNFFNative, TorsionGradientNumerical) {
    // Validate analytical gradient vs numerical
    // ... [finite difference implementation]
}
```

---

## **Quick Reference: File Locations**

### Fortran Source (External)
```
external/gfnff/src/
├── gfnff_engrad.F90          # Main energy/gradient (96k lines)
├── gfnff_ini.f90             # Topology setup (84k lines)
├── gfnff_ini2.f90            # Secondary setup (46k lines)
├── gfnff_param.f90           # Parameters (43k lines)
├── gfnff_cn.f90              # Coordination numbers (6k lines)
├── gfnff_helpers.f90         # Utilities (19k lines)
└── C_interface.f90           # C bindings
```

### C++ Implementation (Native)
```
src/core/energy_calculators/qm_methods/
├── gfnff.h                   # Main class definition (337 lines)
├── gfnff.cpp                 # Implementation (882 lines - needs expansion)
├── gfnff_advanced.h          # Advanced features (179 lines)
├── gfnff_advanced.cpp        # Advanced impl (289 lines - mostly stubs)
├── gfnff_method.h            # ComputationalMethod wrapper (68 lines)
└── gfnff_method.cpp          # Wrapper impl (91 lines)
```

---

## **Estimated Porting Effort**

| Component | Fortran Lines | Estimated C++ | Complexity | Time |
|-----------|---------------|---------------|------------|------|
| Torsion energy | 80 | 120 | Medium | 2 days |
| Inversion energy | 50 | 80 | Low | 1 day |
| Ring detection | 300 | 250 | High | 4 days |
| Pi-system detection | 200 | 180 | Medium | 3 days |
| Hybridization | 150 | 120 | Medium | 2 days |
| EEQ solver | 400 | 350 | Very High | 7 days |
| CN derivatives | 100 | 120 | High | 3 days |
| H-bond detection | 200 | 150 | Medium | 3 days |
| **Total** | **1480** | **1370** | - | **~25 days** |

**Note**: Pure implementation time. Add 50-100% for testing/debugging.

---

## **Next Steps for Developer**

1. **Study Fortran Code**: Start with `gfnff_engrad.F90:1041` (torsions)
2. **Extract Parameters**: Identify all parameter arrays needed for torsions
3. **Write Test First**: Create butane torsion test with reference data
4. **Implement Energy Only**: Get torsion energy working before gradients
5. **Add Gradients**: Analytical derivatives with numerical validation
6. **Integrate**: Add to `generateGFNFFParameters()` and test

**Repeat for each component** until roadmap complete! 🚀
