# GFN-FF Missing Features and TODOs

**Status**: January 2026
**Last Updated**: After implementing feta, metal shifts, and ring torsion corrections

This document tracks remaining gaps between Curcuma's native GFN-FF implementation (`cgfnff`) and the XTB reference implementation.

---

## Executive Summary

**Current Status**: ~85% complete for organic chemistry, ~70% complete for metal-containing systems

**Recent Completions** (January 2026):
- ✅ feta (metal η-coordination correction)
- ✅ Metal bond shift factors (metal1/2/3_shift, eta_shift)
- ✅ Ring torsion corrections (FR3, FR4, FR5, FR6)

**Critical Remaining Issues**:
1. Angle energy 92% error (force constants 10× too weak)
2. Torsion overcompensation (-542% error from extra sp3-sp3 torsions)
3. Incomplete fheavy metal-ligand corrections

---

## PART 1: CRITICAL (Blocks Accuracy)

### 1.1 Angle Energy 92% Error Investigation ⚠️ URGENT

**Status**: CRITICAL - Under active investigation
**Impact**: Angle forces systematically ~10× too weak
**Current Error**: 74% accuracy (should be >95%)

**Problem**:
- CH₃OCH₃ angle energy: 0.00018 Eh (Curcuma) vs 0.00178 Eh (XTB) = 89.9% error
- Force constants appear correct but total energy is wrong
- Recent fix: fbsmall calculation order (was using uninitialized equilibrium angle)
- Post-fix: Revealed deeper issue with force constant magnitude

**Location**: `gfnff_method.cpp:1660-2060` - `getGFNFFAngleParameters()`

**Investigation Steps Completed**:
- ✅ Fixed fbsmall using uninitialized params.equilibrium_angle (commit fcc00ca)
- ✅ Verified angl/angl2 parameters match Fortran arrays exactly
- ✅ Confirmed formula: fc = fijk × fqq × f2 × fn × fbsmall × feta
- ✅ All individual factors calculated per Fortran reference

**Remaining Hypotheses**:
- ❓ Unit conversion issue (kcal/mol vs Hartree vs atomic units)?
- ❓ Missing scaling factor (base parameters in wrong units)?
- ❓ Energy calculation formula in forcefieldthread.cpp has error?
- ❓ Hybridization detection causing wrong equilibrium angles?

**Next Steps**:
1. Compare step-by-step with XTB verbose output for same molecule
2. Check units of angle_params and angl2_neighbors arrays
3. Verify energy calculation formula in `forcefieldthread.cpp:736-960`
4. Test with simple molecule (H₂O) where all values are known
5. Investigate if equilibrium angle calculation is fundamentally wrong

**Files to Review**:
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:1660-2060`
- `src/core/energy_calculators/ff_methods/forcefieldthread.cpp:736-960`
- `src/core/energy_calculators/ff_methods/gfnff_par.h:339-388`

---

### 1.2 Torsion Overcompensation Fix ⚠️ HIGH PRIORITY

**Status**: Implementation complete but overcompensating
**Impact**: Extra sp³-sp³ torsions too strong (-542% error)
**Current Error**: CH₃OCH₃ torsion: -0.000104 Eh vs +0.000023 Eh (reference)

**Problem**:
- Before extra torsions: +0.000073 Eh (215% error, too large)
- After extra torsions: -0.000104 Eh (542% error in opposite direction!)
- Reference value: +0.000023 Eh

**Implementation**: `gfnff_torsions.cpp:1181-1392` - Extra n=1 torsions for sp3-sp3 bonds

**Mechanism**:
- 6 extra torsions generated for CH₃OCH₃
- Oxygen factor: `ff = -2.00` (too strong)
- Carbon factor: `ff = -0.90`
- Nitrogen factor: `ff = 0.70`

**Required Actions**:
1. Calibrate ff=-2.00 oxygen factor (currently overcompensates)
2. Verify extra torsion count matches XTB 6.6.1 verbose output
3. Check if extra torsions should only apply to specific quartet geometries
4. Test with multiple molecules (ethane, methylamine) to verify heteroatom factors
5. Compare gauche vs anti energy differences with XTB reference

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_torsions.cpp:1181-1392`
- `src/core/energy_calculators/ff_methods/gfnff_par.h:795-797` (torsf_extra constants)

---

### 1.3 Complete fheavy Metal-Ligand Corrections ⚠️ MEDIUM PRIORITY

**Status**: Partially implemented
**Impact**: Transition metal-ligand bond strengths incorrect
**Current**: Basic TM detection works, but incomplete element-specific logic

**XTB Reference**: `gfnff_ini.f90:1237-1248`

**Implemented**:
- ✅ TM-heavy ligand (Z>10): fheavy = 0.65
- ✅ TM-P: fheavy = 1.60 (phosphine ligands)
- ✅ TM-chalcogen (S, Se, Te): fheavy = 0.85
- ✅ TM-halogen (F, Cl, Br, I): fheavy = 1.30

**Missing**:
- ❌ Additional element-specific corrections from gfnff_ini.f90:1212-1245
- ❌ M-CO and M-CN special cases (already implemented but needs validation)
- ❌ Validation with actual metal complexes

**Required Actions**:
1. Review Fortran code for any missing fheavy cases
2. Test with transition metal complexes (Fe-CO, Cu-NH3, etc.)
3. Validate bond strengths against XTB reference

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:1463-1513`

---

## PART 2: HIGH PRIORITY (Improves Accuracy)

### 2.1 Ring Angle Equilibrium Angles

**Status**: Documented but not connected
**Impact**: Ring strain energies incorrect
**Location**: `gfnff_method.cpp:1876-1883`

**XTB Reference**: `gfnff_ini.f90` - Element and ring-specific equilibrium angles

**Current Status**:
```cpp
// Ring-dependent equilibrium angles (Phase 9)
// 3-ring: ~82 degrees (very strained)
// 4-ring: ~96 degrees (strained)
// 5-ring: ~105 degrees (near ideal)
// 6-ring: ~112 degrees (near ideal)
```

**Problem**: Hardcoded angle values noted but **not yet used** in angle parameter generation

**Required Implementation**:
1. Integrate ring size detection into angle equilibrium calculation
2. Override r0_deg based on ring_size for angles in rings
3. Element-specific corrections (C vs Si vs O in rings)
4. Test with small ring molecules (cyclopropane, cyclobutane)

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:1876-1883`

---

### 2.2 Small Ring X-H Bond Corrections (fxh)

**Status**: Placeholder code exists
**Impact**: Minor (3-ring C-H bond accuracy)
**Location**: `gfnff_method.cpp:1347-1352`

**XTB Reference**: `gfnff_ini.f90:1177-1190`

**Current Code**:
```cpp
// Small ring effects (Phase 9)
bool is_3ring = false;  // Placeholder (will be from ring detection in Phase 9)
if (is_3ring) {
    fxh *= 1.05;  // 3-ring CH: +5% (strained C-H bonds)
}
```

**Problem**: `is_3ring` always false, correction **never applied**

**XTB Corrections**:
- 3-ring C-H: fxh = 1.05 (+5% stronger, strained)
- Aldehyde C-H: fxh = 0.95 (-5% weaker, conjugated)
- B-H: fxh = 1.10 (+10%)
- N-H: fxh = 1.06 (+6%)
- O-H: fxh = 0.93 (-7%, stiff)

**Required Implementation**:
1. Connect ring detection to fxh calculation
2. Detect atom membership in 3-rings
3. Implement all element-specific X-H corrections
4. Test with cyclopropane, formaldehyde, ammonia, water

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:1347-1352`

---

### 2.3 Element-Specific Angle f2 Corrections

**Status**: Only water implemented
**Impact**: Element-specific angle accuracy
**Location**: `gfnff_method.cpp:1800-1850`

**XTB Reference**: `gfnff_ini.f90:1463-1599`

**Currently Implemented**:
- ✅ Water (H-O-H): r0=100°, f2=1.20

**Missing Elements**:
- ❌ Nitrogen angle corrections
- ❌ Carbon aromatic angle corrections
- ❌ Silicon angle corrections
- ❌ Phosphorus angle corrections
- ❌ Sulfur angle corrections

**Required Implementation**:
1. Port complete f2 correction logic from Fortran
2. Element-specific equilibrium angles (beyond water)
3. Aromatic system detection for angle corrections
4. Test with diverse molecules (amines, aromatics, organosilicon)

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:1800-1850`
- Reference: `external/gfnff/src/gfnff_ini.f90:1463-1599`

---

### 2.4 Ring Torsion Full Quartet Check

**Status**: Simplified assumption used
**Impact**: Minor (may apply ring torsions to partially-cyclic quartets)
**Location**: `gfnff_torsions.cpp:567`

**Current Implementation**:
```cpp
// For now, assume in_ring means central bond j-k is in ring
// TODO Phase 2: Add quartet ring membership check (requires path finding)
bool all_in_same_ring = in_ring;  // Simplified assumption
```

**XTB Reference**: `gfnff_ini.f90:1819, 1824, 1828` - `ringl == rings4` check

**Problem**:
- Currently checks if central bond j-k is in a ring
- Should check if **all 4 atoms** i-j-k-l are in the **same** ring
- May incorrectly apply ring torsions to edge cases

**Required Implementation**:
1. Implement path-finding to check if all 4 atoms share a common ring
2. Compare `ringstors()` (smallest ring containing quartet) with `ringstorl()` (largest)
3. Only apply ring torsions if quartet is fully contained in ring

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_torsions.cpp:562-568`
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:2398-2463` (ring detection)

---

## PART 3: MEDIUM PRIORITY (Niche Cases)

### 3.1 Aromatic Ring Detection

**Status**: Not implemented
**Impact**: Minor (enhanced dispersion for aromatics)
**Purpose**: Improve D3/D4 accuracy for aromatic systems

**XTB Implementation**: Implicit in π-system detection (piadr, pi bond orders)

**Curcuma Status**:
- Basic π-system checks exist (sp/sp² hybridization)
- No explicit aromatic detection
- D3 already handles aromatics reasonably well without flags

**Required Implementation**:
1. Hückel rule checking (4n+2 π-electrons)
2. Planar ring geometry validation
3. Conjugation verification
4. Enhanced C6 parameters for aromatic systems

**Benefit**: Minor improvement in dispersion energy for benzene, naphthalene, etc.

**Files**: New functionality needed

---

### 3.2 Metal-Specific C6 Parameters

**Status**: Uses generic atomic C6 approximations
**Impact**: Unclear (needs investigation)
**Location**: `gfnff_par.h:505-533`

**XTB Reference**: D3/D4 dispersion with metal-specific C6 values

**Current**:
- Generic atomic C6 approximations used
- May not be accurate for transition metals
- D4 library may provide proper values

**Investigation Needed**:
1. Check if D4 library provides proper metal C6
2. Compare with XTB metal-specific parameters
3. Test with metal-containing molecules
4. Benchmark dispersion energies

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_par.h:505-533`
- `src/core/energy_calculators/ff_methods/d4param_generator.cpp`

---

### 3.3 Hydrogen Bonding Integration

**Status**: Framework exists but not integrated
**Impact**: Minor (D3 captures most HB effects)
**Location**: `gfnff_method.cpp:809`

**XTB Reference**: Full HB energy term with basicity/acidity parameters

**Curcuma Status**:
- ✅ `detectHydrogenBonds()` method exists
- ✅ HB parameters in `gfnff_par.h:628-681` (basicity, acidity arrays)
- ❌ Not integrated into force field energy calculation
- ❌ No HB energy contribution

**Parameters Available**:
- `hb_basicity[86]` - HB acceptor basicity
- `hb_acidity[86]` - HB donor acidity
- Global cutoffs (HB_BACUT, HB_SCUT, HB_ALP)
- Scaling factors (XHACI_GLOBABH, XHACI_COH, XHACI_GLOB)

**Required Implementation**:
1. Add HB energy term to `forcefieldthread.cpp`
2. Calculate HB contributions in energy loop
3. Implement HB gradient calculation
4. Test with hydrogen-bonded systems (water dimer, DNA bases)

**Effort**: HIGH (new energy term implementation)

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:809`
- `src/core/energy_calculators/ff_methods/forcefieldthread.cpp` (energy calculation)
- `src/core/energy_calculators/ff_methods/gfnff_par.h:628-681`

---

### 3.4 Halogen Bonding Integration

**Status**: Framework exists but not integrated
**Impact**: Very minor (specialized chemistry)
**Location**: `gfnff_method.cpp:817`

**XTB Reference**: XB energy term with acidity parameters

**Curcuma Status**:
- ✅ `detectHalogenBonds()` method exists
- ✅ XB parameters in `gfnff_par.h:686-702` (acidity array)
- ❌ Not integrated into force field energy calculation
- ❌ No XB energy contribution

**Parameters Available**:
- `xb_acidity[86]` - Halogen bond acidity
- Global cutoffs (XB_BACUT, XB_SCUT)

**Required Implementation**:
1. Add XB energy term to `forcefieldthread.cpp`
2. Calculate XB contributions for I, Br, Cl
3. Implement XB gradient calculation
4. Test with halogen-bonded systems

**Effort**: HIGH (new energy term implementation)

**Benefit**: Niche (important for drug design, crystal engineering)

**Files**:
- `src/core/energy_calculators/ff_methods/gfnff_method.cpp:817`
- `src/core/energy_calculators/ff_methods/forcefieldthread.cpp`
- `src/core/energy_calculators/ff_methods/gfnff_par.h:686-702`

---

## PART 4: LOW PRIORITY (Refinements)

### 4.1 Hyperconjugation Effects

**Status**: Documented but not implemented
**Impact**: Subtle barrier modulation
**Description**: Fine-tuning of torsional barriers based on orbital overlap

**XTB Implementation**: Implicit in torsion parameter assignments

**Required**: Detailed analysis of Fortran torsion logic for hyperconjugation terms

---

### 4.2 Fragment-Constrained EEQ Charges

**Status**: Not implemented
**Impact**: Multi-fragment systems
**Description**: Independent EEQ for separate molecular fragments

**Use Case**: Host-guest complexes, molecular clusters

**Current**: Treats entire system as single EEQ domain

---

### 4.3 Conjugation Detection for Torsions

**Status**: Partial (notpicon check exists)
**Impact**: π-conjugated system accuracy
**Location**: `gfnff_torsions.cpp:572`

**Current**:
- Basic check: `notpicon = !(hyb_j == 2 && hyb_k == 2)`
- Only checks sp²-sp² hybridization
- Doesn't detect extended conjugation

**Enhancement Needed**:
1. Detect extended conjugated chains
2. Amide system detection (already exists but not fully used)
3. Aromatic system integration
4. Adjust barriers for conjugated vs isolated double bonds

---

## PART 5: VALIDATION STATUS

### Test Molecule Performance (CH₃OCH₃ vs XTB 6.6.1)

| Component | Curcuma (Eh) | XTB Ref (Eh) | Error % | Status |
|-----------|--------------|--------------|---------|--------|
| **Bond**      | -1.225128    | -1.216444    | **+0.71**   | ✅ EXCELLENT |
| **Angle**     | 0.000180     | 0.001780     | **-89.89**  | ❌ CRITICAL |
| **Torsion**   | 0.000074     | 0.000023     | **+215.14** | ⚠️ Overcompensating |
| **Repulsion** | 0.054074     | 0.053865     | **+0.39**   | ✅ EXCELLENT |
| **Coulomb**   | -0.043848    | -0.047825    | **+8.32**   | ✅ GOOD |
| **Dispersion**| -0.001896    | -0.000042    | N/A         | ✅ Working (D4 vs D3) |
| **TOTAL**     | **-1.216546**| **-1.209209**| **+0.61**   | ✅ EXCELLENT |

**Summary**:
- 4/6 components excellent (<1% error)
- 2/6 components have issues (angle, torsion)
- Total energy excellent (0.6% error) due to error cancellation
- Need to fix angle and torsion individually for robust accuracy

---

## PART 6: IMPLEMENTATION ROADMAP

### Phase 1: Critical Fixes (Required for Production)

**Priority**: URGENT
**Timeline**: Immediate

1. **Angle Energy Investigation** (2-3 days)
   - Debug force constant calculation
   - Verify energy formula in forcefieldthread.cpp
   - Test with simple molecules (H₂O, NH₃, CH₄)
   - Target: <5% error

2. **Torsion Calibration** (1-2 days)
   - Adjust ff=-2.00 oxygen factor
   - Test with ethane, methylamine, dimethyl ether
   - Verify extra torsion count matches XTB
   - Target: <10% error

3. **fheavy Completion** (1 day)
   - Review Fortran for missing cases
   - Test with metal complexes
   - Validate bond strengths

### Phase 2: High-Priority Improvements

**Priority**: HIGH
**Timeline**: 1-2 weeks

1. **Ring Angle Equilibrium Angles** (2 days)
   - Connect ring detection to angle calculation
   - Test with cyclopropane, cyclobutane, cyclohexane

2. **Small Ring X-H Corrections** (1 day)
   - Implement is_3ring detection
   - Add all element-specific X-H corrections

3. **Element-Specific f2** (2-3 days)
   - Port complete f2 logic from Fortran
   - Test with diverse molecules

### Phase 3: Medium-Priority Enhancements

**Priority**: MEDIUM
**Timeline**: 2-4 weeks

1. **Aromatic Detection** (3 days)
2. **Metal C6 Investigation** (2 days)
3. **Ring Torsion Quartet Check** (2 days)

### Phase 4: Optional Features

**Priority**: LOW
**Timeline**: As needed

1. **Hydrogen Bonding** (1 week)
2. **Halogen Bonding** (1 week)
3. **Hyperconjugation** (research needed)

---

## PART 7: KNOWN ARCHITECTURAL LIMITATIONS

### 7.1 Parameter Caching

**Current**: Universal caching system works well
**Limitation**: Thread-safety requires manual control
**Solution**: Use `setParameterCaching(false)` for concurrent calculations

### 7.2 Ring Detection Algorithm

**Current**: BFS implementation robust for 3-8 membered rings
**Limitation**: Doesn't detect all ring types (spiro, bridged, fused)
**Impact**: Minor (most common rings work correctly)

### 7.3 EEQ Solver Performance

**Current**: Two-phase system accurate (8.3% error)
**Limitation**: Iterative refinement may be slow for large systems
**Optimization**: Could implement matrix caching for repeated geometries

---

## PART 8: TESTING REQUIREMENTS

### Required Test Cases

**Small Rings**:
- [ ] Cyclopropane (3-ring angles, torsions, strain)
- [ ] Cyclobutane (4-ring puckering)
- [ ] Cyclopentane (5-ring envelope)
- [ ] Cyclohexane (6-ring chair/boat)

**Metal Complexes**:
- [ ] Ferrocene (Fe-C₆H₆, η-coordination)
- [ ] Fe-CO (metal-carbonyl bonding)
- [ ] Cu-NH₃ (metal-amine complexes)
- [ ] Organolithium compounds (Group 1 metals)

**Aromatics**:
- [ ] Benzene (aromatic angles, torsions)
- [ ] Naphthalene (fused rings)
- [ ] Pyridine (aromatic N)

**Hydrogen Bonding**:
- [ ] Water dimer
- [ ] DNA base pairs
- [ ] Protein backbone (amide HB)

**Torsional Barriers**:
- [ ] Ethane (sp³-sp³ baseline)
- [ ] Butane (gauche/anti)
- [ ] Dimethyl ether (O sp³-sp³)
- [ ] Methylamine (N sp³-sp³)

---

## PART 9: CONTACT AND REFERENCES

**Implementation Notes**: See `/home/conrad/.claude/plans/snuggly-skipping-rainbow.md`

**Reference Implementation**:
- XTB: `external/gfnff/src/gfnff_ini.f90`
- Parameters: `external/gfnff/src/gfnff_param.f90`
- Topology: `external/gfnff/src/gfnff_ini2.f90`

**Key Papers**:
- Spicher, S.; Grimme, S. *Angew. Chem. Int. Ed.* **2020** (GFN-FF method)
- Grimme, S. et al. *J. Chem. Phys.* **2010**, 132, 154104 (DFT-D3)
- Caldeweyher, E. et al. *J. Chem. Phys.* **2019**, 150, 154122 (DFT-D4)

**Status Document**: `docs/GFNFF_STATUS.md` (overview, completeness)

**Change Log**: `AIChangelog.md` (development history)

---

**End of Document**
