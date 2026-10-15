# GFN-FF Unit System Analysis
**Date**: December 29, 2025
**Status**: Critical Unit Mismatch Identified

## Problem Summary

The GFN-FF implementation has a **fundamental unit system inconsistency**:

| Component | Unit | Source |
|-----------|------|--------|
| **Molecular Geometry (m_geometry)** | Ångström | XYZ file (no conversion) |
| **Bond Distance (rij)** | Ångström | Calculated from m_geometry |
| **Bond r0_ij Parameter** | **Bohr** | gfnff_method.cpp:1089 |
| **Energy Calculation** | Assumes consistent units |  |

## Evidence

### 1. Parameter Generation (gfnff_method.cpp:1089)

```cpp
// Line 1089-1090: Explicitly states r0_ij is in BOHR
CurcumaLogger::info(fmt::format(
    "RAB_TRANSFORM: r_eq = ({:.8f} + {:.8f} + {:.8f}) * {:.8f} = {:.8f} **Bohr**",
    ra, rb, rabshift, ff, params.equilibrium_distance));
```

**Confirmation**: Code explicitly logs "Bohr"

### 2. Generated Parameters (CH3OCH3.param.json)

```json
{
  "i": 0, "j": 6,        // H-C bond
  "r0_ij": 1.836097      // In Bohr → 0.9716 Å
}
```

**Analysis**:
- If r0_ij is in **Bohr**: 1.836097 × 0.529167 = **0.9716 Å** vs Reference **0.977 Å** → **0.55% error** ✓
- If r0_ij is in **Ångström**: 1.836097 Å vs Reference 0.977 Å → **87.93% error** ✗

**Conclusion**: r0_ij IS in Bohr

### 3. Energy Calculation (forcefieldthread.cpp:720)

```cpp
// Line 716: rij is calculated from m_geometry (Ångström)
double rij = UFF::BondStretching(i, j, derivate, m_calculate_gradient);

// Line 720: Bond distance difference
double dr = rij - bond.r0_ij;  // r(Å) - r(Bohr) = UNIT MISMATCH!

// Line 723: Energy calculation
double energy = k_b * std::exp(-alpha * dr * dr);
```

**Problem**:
- `rij` is in Ångström (from m_geometry, which is from XYZ file)
- `bond.r0_ij` is in Bohr (from parameter generation)
- `dr = rij - bond.r0_ij` mixes units!

### 4. Geometry Loading (molecule.cpp)

```cpp
// Line 250: Direct load from XYZ with NO conversion
m_geometry.row(i) = Eigen::Vector3d(position[0], position[1], position[2]);
// position is in Ångström from XYZ file
// No Bohr conversion anywhere
```

**Conclusion**: m_geometry is in Ångström

## Impact

### Energy Calculation Error

The bond energy formula:
```
E_bond = k_b × exp(-α × (r - r₀)²)
```

With r in Ångström and r₀ in Bohr:
- All bond distances are ~2× smaller than they should be
- Exponential term gets wrong argument
- Calculated energy is completely wrong

**Example**: H-C bond
- r (actual from geometry) = ~1.09 Å
- r₀ (in Bohr) = 1.846 Bohr ≈ 0.9716 Å
- dr = 1.09 Å - 1.846 Bohr = **mixing apples and oranges**

This explains:
- ❌ Test 8: 47.32% error in bond energy even with "reference" parameters
- ❌ Large errors in Coulomb energy (r_ij vs R_ij unit mismatch)
- ❌ Repulsion energy errors

## Solution Options

### Option A: Convert Geometry to Bohr (Global Fix)
**Location**: Molecule class constructor
```cpp
void Molecule::setGeometry(const std::vector<std::array<double, 3>>& geometry) {
    const double angstrom_to_bohr = 1.0 / 0.529167;
    for (size_t i = 0; i < geometry.size(); ++i) {
        m_geometry(i, 0) = geometry[i][0] * angstrom_to_bohr;
        m_geometry(i, 1) = geometry[i][1] * angstrom_to_bohr;
        m_geometry(i, 2) = geometry[i][2] * angstrom_to_bohr;
    }
}
```
**Pros**: Single fix point, all subsequent calculations correct
**Cons**: Affects all code using m_geometry (XYZ output, RMSD, etc.)
**Risk**: High - could break many other features

### Option B: Convert r0_ij to Ångström (Targeted Fix)
**Location**: gfnff_method.cpp line 1083
```cpp
// Convert r0_ij from Bohr to Ångström
params.equilibrium_distance = rtmp * 0.529167;  // Now in Ångström
```
**Pros**: Minimal scope, only affects GFN-FF
**Cons**: Must check ALL parameter generations (angles, torsions, vdW)
**Risk**: Medium - need to verify all parameter types

### Option C: Convert in ForceField::setParameter
**Location**: forcefield.cpp line 457
```cpp
// Convert Bohr→Ångström when loading parameters
b.r0_ij = bond["r0_ij"] * 0.529167;  // Convert to Ångström
```
**Pros**: Single point, protects against future parameter bugs
**Cons**: Assumes all force field methods use Bohr (must verify UFF)
**Risk**: Medium - need external FF validation

## Recommended Approach

**Option B** (Convert in parameter generation) with verification:

1. **Identify all GFN-FF parameter calculations** that output distances:
   - ✅ Bond r0_ij (gfnff_method.cpp:1083)
   - ❓ Angle r0_ij/r0_ik (angles.cpp)
   - ❓ VdW r0_ij (dispersion.cpp)
   - ❓ Repulsion r0_ij (if used)

2. **Apply conversion** at parameter generation
   ```cpp
   // In each parameter generation method
   params.equilibrium_distance = rtmp * 0.529167;  // Bohr → Ångström
   ```

3. **Validate** with Test 5:
   - Bond r0 errors should drop from 88% → <1%
   - Test 8 bond energy should improve from 47% → ?

4. **Run full test suite** to catch other unit issues

## Test Impact

Once fixed:
- **Test 4** (CN): Unchanged (doesn't involve distances)
- **Test 5** (Bond Params): r0_ij errors drop from **88%** to **<1%** ✓
- **Test 8** (Bond Energy): Should improve significantly
  - Current: 47% error with "correct" r0_ij in Bohr
  - Expected: <5% error with r0_ij in correct units

## Related Issues

This unit mismatch may also affect:
- **Coulomb energy** calculation (if using r_ij without conversion)
- **Repulsion energy** (exponential term with mixed units)
- **Angle damping** (uses bond distances)
- **Any other FF method** using Bohr parameters with Ångström coordinates

## Files to Check

1. **gfnff_method.cpp**: Parameter generation
   - Line 1083: r0_ij calculation
   - Search for all `equilibrium_distance = ` assignments

2. **forcefieldthread.cpp**: Energy calculations
   - Line 716: rij calculation (confirm Ångström)
   - Line 720: dr calculation (unit mismatch point)
   - All dispersion/repulsion/coulomb calculations

3. **forcefield.cpp**: Parameter loading
   - Line 457: Bond parameter assignment
   - Line 489: Angle parameter assignment
   - Line 582: VdW parameter assignment

4. **eeq_solver.cpp**: Charge calculation
   - Check if any distances are used
   - Verify units in distance-dependent formulas

## Next Steps

1. Identify ALL distance parameters in GFN-FF
2. Check if conversions are needed in:
   - Parameter generation (gfnff_method.cpp)
   - Parameter loading (forcefield.cpp)
   - Energy calculation (forcefieldthread.cpp)
3. Apply fixes systematically
4. Re-run Test 5 to confirm r0_ij errors drop to <1%
5. Re-run Test 8 to see energy calculation improvement
