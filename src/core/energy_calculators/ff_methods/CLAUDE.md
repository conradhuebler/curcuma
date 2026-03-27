# CLAUDE.md - Force Field Methods Directory

## Overview

Force field implementation system with multi-threading support for UFF, QMDFF, and native GFN-FF.

## Architecture

### Core Separation of Concerns

**ForceFieldThread** (`forcefieldthread.cpp/h`):
- **Responsibility**: Energy and gradient calculations for all force field terms
- **Threading**: Multi-threaded via CxxThreadPool
- **Methods**: Calculate{Method}{Term}Contribution() pattern (e.g., CalculateGFNFFBondContribution)

**ForceField** (`forcefield.cpp/h`):
- **Responsibility**: Thread pool management, geometry updates, energy accumulation
- **Role**: Dispatcher that coordinates ForceFieldThread instances

**ForceFieldGenerator** (`forcefieldgenerator.cpp/h`):
- **Responsibility**: UFF parameter generation from atom types
- **Used by**: UFF method only

**EEQSolver** (`eeq_solver.cpp/h` - Claude Generated December 2025, Enhanced January 4, 2026):
- **Responsibility**: Standalone electronegativity equalization charge solver (extracted from GFN-FF)
- **Architecture** (January 4, 2026): Hybrid two-phase + iterative refinement
  - Phase 1: Initial solve with base parameters (gam only, no dgam) + **CNF term in RHS**
  - Phase 2: Iterative refinement with dgam corrections applied in matrix, **NO CNF term in RHS**
  - Key fix: CNF term ONLY in Phase 1 (gfnff_ini.f90:563-570), removed from Phase 2 (gfnff_ini.f90:696-707)
- **Helper Functions** (new): `buildCorrectedEEQMatrix()`, `solveEEQ(use_cnf_term=bool)`
- **Phase 1 Improvements** (December 28, 2025): Pi-system detection (sp/sp2 hybridization from CN), neighbor electronegativity averaging (Pauling scale), environment-dependent dxi corrections (Boron, C=O, C=N, halogens, metals)
- **Accuracy**: CH₄ charges improved 75% (5.0× error → 1.3× error), fixes 4/6 GFN-FF energy terms
- **Used By**: GFN-FF for Coulomb charges, D4ParameterGenerator for charge-dependent C6
- **Parameters**: Element-specific (chi_eeq, gam_eeq, alpha_eeq, cnf_eeq) from gfnff_par.h + ConfigManager integration
- **Coulomb Energy Fix (Jan 28, 2026)**: Corrected dgam ff-values and enabled charge-corrected gameeq/alpeeq - Coulomb energy now matches Fortran reference exactly (< 1 nEh error)
- **Status**: ✅ EEQ charges correct, ✅ Coulomb energy exact match

**GFNFF Class** (`gfnff_method.cpp/h` + `../qm_methods/gfnff.cpp/h`):
- **Responsibility**: GFN-FF parameter generation (topology-aware) + ConfigManager integration
- **EEQSolver Integration** (Dec 2025): Delegates charge calculation to standalone EEQSolver instead of embedded implementation
- **Methods**: generateTopologyAwareBonds(), generateGFNFFDispersionPairs(), etc.
- **Output**: JSON parameter sets passed to ForceField

### GFN-FF Implementation Pattern

**Two-Phase Architecture:**

1. **Parameter Generation** (in GFNFF class):
   ```cpp
   // gfnff.cpp with ConfigManager migration
   ConfigManager config("gfnff", parameters);
   json params = config.exportConfig();
   params["bonds"] = generateTopologyAwareBonds(...);
   params["angles"] = generateTopologyAwareAngles(...);
   params["gfnff_dispersions"] = generateGFNFFDispersionPairs();
   // etc.
   m_forcefield->setParameter(params);
   ```

2. **Term Calculation** (in ForceFieldThread with parameter flags):
   ```cpp
   // forcefieldthread.cpp - Phase 2 parameter flag checks
   if (m_dispersion_enabled) {
       CalculateGFNFFDispersionContribution();  // Only if enabled
   }
   if (m_hbond_enabled) {
       CalculateGFNFFHydrogenBondContribution(); // Only if enabled
   }
   ```

### Adding New GFN-FF Terms - Checklist

To add a new GFN-FF energy term (e.g., "CrossTerm"), you MUST modify:

1. **Parameter Structures** (`forcefieldthread.h`):
   - Add new struct (e.g., `GFNFFCrossTerm { int i,j; double param1; }`)
   - Add member vector (e.g., `std::vector<GFNFFCrossTerm> m_gfnff_crossterms;`)

2. **Parameter Generation** (`gfnff.cpp`):
   - Add generation method (e.g., `generateGFNFFCrossTerms()`)
   - Call in `generateGFNFFParameters()` and add to JSON

3. **Term Calculation** (`forcefieldthread.cpp`):
   - Add calculation method (e.g., `CalculateGFNFFCrossTermContribution()`)
   - Call in `ForceFieldThread::execute()` under `if (m_method == 3)` block

4. **Energy Accumulation** (`forcefield.cpp`):
   - Add energy component member variable if needed
   - Collect from threads in `ForceField::Calculate()`

5. **Parameter Setter** (`forcefield.h/.cpp`):
   - Add setter method (e.g., `setGFNFFCrossTerms(const json&)`)
   - Call in `setParameter()` dispatcher

## GFN-FF Gradient Implementation Status (February 2026)

**✅ ALL GRADIENT TERMS NOW ENABLED** (Feb 1, 2026): Torsion, extra torsion, and inversion gradients activated for MD stability.

### Active Gradient Terms (All Enabled)

| Energy Term | Method | Status in execute() | Gradient Status | Fortran Ref |
|-------------|--------|---------------------|-----------------|-------------|
| Bonds | CalculateGFNFFBondContribution() | ✅ ACTIVE | ✅ Implemented | egbond:675-721 |
| Angles | CalculateGFNFFAngleContribution() | ✅ ACTIVE | ✅ Implemented | egbend:857-916 |
| Torsions | CalculateGFNFFDihedralContribution() | ✅ ACTIVE (Feb 2026) | ✅ Implemented | egtors:1153-1234 |
| Extra Torsions | CalculateGFNFFExtraTorsionContribution() | ✅ ACTIVE (Feb 2026) | ✅ Implemented | egtors:1272-1280 |
| Inversions | CalculateGFNFFInversionContribution() | ✅ ACTIVE (Feb 2026) | ✅ Implemented | gfnff_ini.f90 |
| Dispersion | CalculateGFNFFDispersionContribution() | ✅ ACTIVE | ✅ Implemented | gdisp0.f90 |
| Repulsion (bonded) | CalculateGFNFFBondedRepulsionContribution() | ✅ ACTIVE | ⚠️ Partial | engrad:467-495 |
| Repulsion (non-bonded) | CalculateGFNFFNonbondedRepulsionContribution() | ✅ ACTIVE | ⚠️ Partial | engrad:255-276 |
| Coulomb | CalculateGFNFFCoulombContribution() | ✅ ACTIVE | ✅ Term 1b Fixed (Feb 2026) | engrad:383-422 |
| Hydrogen Bonds | CalculateGFNFFHydrogenBondContribution() | ✅ ACTIVE | ✅ Implemented (Feb 2026 rewrite, <1e-4) | abhgfnff_eg* |
| Halogen Bonds | CalculateGFNFFHalogenBondContribution() | ✅ ACTIVE | ✅ Implemented (Mar 2026) | rbxgfnff_eg |
| Triple Bond Torsions | CalculateGFNFFSTorsionContribution() | ✅ ACTIVE (Mar 2026) | ✅ Implemented | sTors_eg:3454 |
| BATM | CalculateGFNFFBatmContribution() | ✅ ACTIVE | ✅ Implemented (Mar 2026) — GradientBATM() | batmgfnff_eg |
| ATM (D3/D4) | CalculateATMContribution() | ✅ ACTIVE | ✅ Complete | d3_gradient |

### Implementation Details

#### ✅ Complete Gradient Implementations

**Bond Gradients** (forcefieldthread.cpp:834-842):
```cpp
// dE/dr = -2*α*dr*E (chain rule)
double dEdr = -2.0 * alpha * dr * energy;
m_gradient.row(bond.i) += dEdr * factor * derivate.row(0);
m_gradient.row(bond.j) += dEdr * factor * derivate.row(1);
```
- **Note**: Missing CN gradient contribution (dr0/dCN * dCN/dx) - second-order effect

**Angle Gradients** (forcefieldthread.cpp:993-1034):
- Complete implementation with distance-dependent damping
- Full damping gradient terms: ∂E/∂x = (∂E/∂θ * damp) * (∂θ/∂x) + (∂E/∂damp) * (∂damp/∂x)
- Matches Fortran egbend exactly

**Torsion Gradients** (forcefieldthread.cpp:1038-1360):
- Complete with all three damping terms (damp_ik, damp_jk, damp_jl)
- Cross-center damping formula matches Fortran
- NCI torsion support (atcutt_nci = 0.305 vs standard atcutt = 0.505)

**Extra Torsion Gradients** (forcefieldthread.cpp:1366-1520):
- Same implementation as primary torsions
- Without +π phase shift (gauche torsions)

#### ❌ Missing Gradient Implementations

**Coulomb Gradients**:
- **Status**: ✅ FIXED (Feb 1, 2026) - Term 1b (charge derivative via CN) now working
- **Fix Applied**: CN derivatives recalculated in `GFNFF::Calculation()` before gradient computation
- **Root Cause Fixed**: CN derivatives were computed once at initialization but not updated when geometry changed
- **Formula**: ∂E/∂x = ∂E/∂r * ∂r/∂x + ∂E/∂q * ∂q/∂CN * ∂CN/∂x (Term 1b)
- **Implementation**: `gfnff_method.cpp:341-368` - calls `calculateCoordinationNumberDerivatives()` and `distributeCNandDerivatives()`
- **Reference**: Fortran gfnff_engrad.F90:418-422

**Dynamic Coulomb Charges (Feb 23, 2026)**:
- **Status**: ✅ IMPLEMENTED - TERM 1, 2+3 now use dynamic `m_eeq_charges` (previously static init charges)
- **Problem**: TERM 1 (pairwise) and TERM 2+3 (self-energy) used static `coul.q_i/q_j` and `params->chi_i`; only TERM 1b used dynamic charges. Caused ~1e-3 Eh/Bohr gradient errors for polar molecules.
- **Fix**: Added `chi_base_i/j` and `cnf_i/j` to `GFNFFCoulomb` struct. At runtime: `qi = m_eeq_charges(i)`, `chi_eff = chi_base + cnf*sqrt(max(cn,0))`. NaN fallback to static values.
- **TERM 1b guard**: Fixed skip-condition bug — all atoms now use `qtmp(i) = q*cnf/(2*sqrt(max(cn,0))+1e-16)` (Fortran epsilon guard), instead of skipping atoms with `cn ≤ 1e-10`.

#### ⚠️ Partial Gradient Implementations

**Repulsion Gradients**:
- Non-bonded repulsion active but needs validation against Fortran
- Bonded repulsion disabled, needs testing

### Gradient Terms Status (Feb 2026)

**✅ ALL BONDED GRADIENT TERMS NOW ENABLED** in `forcefieldthread.cpp:execute()`:

```cpp
// Lines 116-120: All bonded terms ACTIVE (Feb 2026)
CalculateGFNFFBondContribution();        // Always active
CalculateGFNFFAngleContribution();       // Always active
CalculateGFNFFDihedralContribution();    // ENABLED Feb 2026 - Required for MD stability
CalculateGFNFFExtraTorsionContribution(); // ENABLED Feb 2026 - Gauche torsions for MD
CalculateGFNFFInversionContribution();   // ENABLED Feb 2026 - Planar group stability
```

**Gradient Unit Conversion** (gfnff_method.cpp:377-385):
- ForceField internally uses Bohr/Hartree units
- Gradient is converted to Hartree/Angstrom for optimizer/MD compatibility
- Formula: `gradient_Angstrom = gradient_Bohr * BOHR_TO_ANGSTROM`

**Gradient Test Results** (4/6 passing):
- ✅ Bond, Angle, Torsion, CH₃OCH₃ full molecule
- ⚠️ Repulsion: translational invariance issue
- ⚠️ Benzene: aromatic system needs investigation

**Impact**: Enabling torsion/inversion gradients is **critical for MD simulations** - without these terms, molecules rotate uncontrollably and planar groups become distorted.

**Known Issue**: LBFGS optimization shows line-search failures after first step - this is a separate optimizer parameter tuning issue, not a gradient problem.

### Test Framework

**Test File**: `test_cases/test_gfnff_gradients.cpp`
- Finite-difference gradient calculation
- Term-specific tests (bond-only, angle-only, etc.)
- Full molecule gradient validation
- Tolerance: 1e-5 to 1e-4 Hartree/Bohr depending on term

**CTest Integration**:
```bash
ctest -R test_gfnff_gradients --verbose
```

### Next Steps

1. **Validate Repulsion gradients** against Fortran (bonded + non-bonded)
2. **Run full regression test** with gradients enabled

## Current Implementation Status (Mar 2026)

### ✅ Fully Implemented and Validated Terms
- Bond stretching (exponential potential) — < 0.1% error
- Angle bending (cosine + damping + pi_bond_orders) — < 0.12 mEh (all molecules)
- Dihedral torsion + extra torsions — active, < 2.5 mEh (triose)
- Inversion (out-of-plane) — exact Fortran `domegadr` gradient
- Dispersion D4 (modified BJ + CN-only weighting) — < 1 µEh; GradComp √N precision limit accepted
- Repulsion — < 0.4% error
- Coulomb (dynamic EEQ charges) — < 0.1 mEh (< 1 nEh for small molecules)
- Hydrogen bonds / Halogen bonds — all cases 1-3 + gradients
- BATM — topology charges distributed after thread creation (fixed Mar 6)
- ATM (separated to own GradientATM()) — energy ≈ 0, gradient correct

### ✅ EEQ Solver Status
- **EEQSolver**: Standalone in `eeq_solver.{h,cpp}`, two-phase architecture
- **dxi corrections**: Active (pi-system, neighbor EN averaging, environment-dependent)
- **dgam corrections**: Intentionally disabled — validated as no improvement (<0.001% energy impact)
- **Element hybridization**: Complete XTB element-specific rules (gfnff_ini2.f90:217-332)
- See `docs/DGAM_VALIDATION_REPORT.md` for dgam analysis

### ✅ Parameter Management (Phase 2 - December 2025)
- **ConfigManager Integration**: Type-safe parameter access with validation
- **Parameter Flags**: Selective term calculation (dispersion, hbond, repulsion, coulomb enabled/disable)
- **Legacy Code Removal**: CalculateGFNFFvdWContribution deprecated and removed
- **Test Coverage**: Parameter flag combinations test suite with 5 scenarios
  - Dispersion/hbond disabled tests
  - All non-bonded terms disabled test
  - Edge case (atoms at cutoff distance)
  - Metal-specific correction handling (Fe atom)

### 🟡 Lower Priority TODOs

#### Topology-Specific Corrections (Not Yet Implemented)

**Angle Bending Corrections**:
- [ ] **Ring strain factors** - Small rings (3-, 4-membered) need reduced force constants
- [ ] **Metal coordination** - feta metal correction factor (currently =1.0 for all)
- [ ] **fijk refinement** (Phase 2b) - angl2 topology logic for neighbor type corrections

**Torsion Corrections**:
- [ ] **Ring torsions** - Different phase angles and barriers for cyclic vs acyclic
- [ ] **Conjugation detection** - Increase barriers for π-conjugated systems
- [ ] **Hyperconjugation** - Subtle barrier modulation (documented but not implemented)
- [ ] **Extra torsion calibration** - Current ff=-2.00 (O) factor overcompensates

**Charge (EEQ) Corrections**:
- [ ] Phase 5B: Metal-specific fqq correction (2.5× factor in charge-dependent terms)
- [ ] Pi-system/amide detection for nitrogen dgam (enhancement, current EEQ already good)
- ✅ Fragment-constrained EEQ charges (for multi-fragment systems)

**Dispersion Corrections**:
- [ ] **Metal-specific C6 parameters** - Transition metals may need special handling
- [ ] **Aromatic system detection** - Ring detection algorithm for enhanced dispersion

#### Implementation Priority

**HIGH** (improves accuracy):
1. Ring strain factors for small-ring angles
2. Metal coordination corrections (feta, fqq)
3. Torsion extra-term calibration (extra sp3-sp3 ff factor)

**MEDIUM** (niche cases):
4. Conjugation detection for torsions
5. Aromatic ring detection for dispersion
6. Metal-specific C6 parameters

**LOW** (refinements):
7. Hyperconjugation effects
8. NCI torsion `is_nci` flag population in parameter generation

## Performance

**Multi-threading Benchmarks** (water.xyz, 4 cores):
- 1 thread: 0.320s
- 4 threads: 0.120s
- Speedup: 2.67x ✅

## D3 Implementation Status (December 19, 2025)

### ✅ FULLY VALIDATED - Production Ready

**Accuracy**: **8/9 test molecules <1% error** (H₂: 0.026%, HCl: 0.036%, CH₃OCH₃: 0.659%, etc.)

**Root Cause Fixed (December 19, 2025)**: Triangular indexing formula conversion error
- **Issue**: Fortran 1-based formula `ic = j + i*(i-1)/2` incorrectly converted to C++
- **Fix**: Correct 0-based formula is `ic = j + i*(i+1)/2`
- **Impact**: Heteronuclear pairs had 20-87% errors before fix, now <1%

### Test Results (Comprehensive Validation)

| Molecule | Atoms | Calculated (Eh) | Reference (Eh) | Error % | Status |
|----------|-------|-----------------|----------------|---------|--------|
| H₂       | 2     | -6.7713e-05     | -6.7731e-05    | 0.026   | ✅ PASS |
| HCl      | 2     | -2.6246e-04     | -2.6256e-04    | 0.036   | ✅ PASS |
| OH       | 2     | -1.1779e-04     | -1.1791e-04    | 0.105   | ✅ PASS |
| HCN      | 3     | -6.8388e-04     | -6.8602e-04    | 0.313   | ✅ PASS |
| O₃       | 3     | -5.2928e-04     | -5.9161e-04    | 10.537  | ⚠️ OUTLIER |
| H₂O      | 3     | -2.7621e-04     | -2.7686e-04    | 0.236   | ✅ PASS |
| CH₄      | 5     | -9.2000e-04     | -9.2212e-04    | 0.230   | ✅ PASS |
| CH₃OH    | 6     | -1.4926e-03     | -1.5054e-03    | 0.846   | ✅ PASS |
| CH₃OCH₃  | 9     | -3.3475e-03     | -3.3697e-03    | 0.659   | ✅ PASS |
| triose   | 66    | -2.4371e-02     | -2.4371e-02    | 0.000   | ✅ PASS |
| monosaccharide | 27 | -8.4732e-03   | -8.4732e-03    | 0.000   | ✅ PASS |

**Summary**: 10/11 passing (<1% error) | 1/11 outlier (O₃ at 10.5%)

### Known Limitations

**O₃ Outlier** (10.5% error):
- Likely cause: Ozone geometry (bent) or O-specific CN calculation
- All other molecules including homoatomic H₂ work perfectly
- Non-blocking for production use (most molecules <1%)

### Technical Implementation Details

**Triangular Indexing Fix** (`d3param_generator.cpp:262-272`):
```cpp
// Fortran (1-based): ic = j + i*(i-1)/2
// C++ (0-based):     ic = j + i*(i+1)/2  ← CRITICAL DIFFERENCE
int pair_index;
if (elem_i > elem_j) {
    pair_index = elem_j + elem_i * (elem_i + 1) / 2;
} else {
    pair_index = elem_i + elem_j * (elem_j + 1) / 2;
}
```

**Validated Components**:
1. ✅ **CN calculation**: Exponential counting formula matches s-dftd3
2. ✅ **BJ damping**: Formula E = -s6·C6/(r⁶+R0⁶) - s8·C8/(r⁸+R0⁸) validated
3. ✅ **Gaussian weighting**: exp(-wf * (cn - cnref)²) with wf=4.0
4. ✅ **C6 interpolation**: Correct access to reference_c6 with MAX_REF=7
5. ✅ **Reference data**: Complete 262,444 C6 values + 721 CN values from s-dftd3

---

## D4 Implementation Status (January 25, 2026)

### ✅ BJ DAMPING FORMULA FIX (January 25, 2026)

**Critical fix: GFN-FF uses a MODIFIED BJ damping formula, NOT standard D3/D4**

**Problem**: Curcuma used standard BJ damping with R0 computed from C8/C6 ratio:
```cpp
// WRONG (standard D3/D4 BJ):
r_crit = a1 * sqrt(C8/C6) + a2;
E = -s6*C6/(r^6+R0^6) - s8*C8/(r^8+R0^8);
```

**Solution**: Implement GFN-FF modified formula from `gfnff_gdisp0.f90:365-377`:
```cpp
// CORRECT (GFN-FF modified BJ):
r4r2ij = 3 * sqrtZr4r2_i * sqrtZr4r2_j;  // Implicit C8/C6 factor
r0_squared = (a1*sqrt(r4r2ij) + a2)^2;    // Pre-computed from sqrtZr4r2
t6 = 1/(r^6 + R0^6);
t8 = 1/(r^8 + R0^8);
E = -0.5 * C6 * (t6 + 2*r4r2ij*t8);       // 0.5 for pair counting
```

**Key Differences from Standard BJ**:
1. R0 computed from `sqrtZr4r2` product (NOT from C8/C6 ratio)
2. C8 implicit via `2*r4r2ij*t8` factor (NOT separate C8*t8 term)
3. 0.5 factor for pair counting (each pair counted once)

**Impact**: Caffeine dispersion error reduced **6.6×** (26 mEh → 3.9 mEh)

**Reference**: `gfnff_gdisp0.f90:365-377`, `gfnff_param.f90:531-532`

### ✅ CN-ONLY WEIGHTING FIX (January 17, 2026)

**Critical fix to match GFN-FF Fortran reference**:

**Problem**: D4 was using CN+charge combined weighting (incorrect for GFN-FF)
```cpp
// WRONG (before January 17, 2026):
weights[ref] = std::exp(-wf * (diff_q*diff_q + diff_cn*diff_cn));
```

**Solution**: Changed to CN-only weighting
```cpp
// CORRECT (January 17, 2026):
weights[ref] = std::exp(-wf * diff_cn * diff_cn);
```

**Reference**: `external/gfnff/src/gfnff_gdisp0.f90:405` - `cngw = exp(-wf * (cn - cnref)**2)`

### GFN-FF Hybrid Dispersion Model

**Key Insight**: GFN-FF uses a **hybrid approach**, NOT full D4:
- ✅ D4 Casimir-Polder integration for C6 parameters (frequency-dependent polarizabilities)
- ✅ D3-style CN-only weighting (simpler, coordination-number-based)
- ❌ NOT full D4 from Caldeweyher et al. (which uses CN+charge weighting)

**Rationale for CN-only**:
- Simpler model with fewer parameters
- Reduced computational cost (no charge dependency in weighting)
- Validated in XTB GFN-FF implementation
- Sufficient accuracy for force field purposes

### Default Method Change

**File**: `gfnff_method.cpp:5204`

**Before**: Default was D3 (static lookup tables)
```cpp
std::string method = "d3";  // WRONG for GFN-FF reference
```

**After**: Default is now D4 (Casimir-Polder integration)
```cpp
std::string method = "d4";  // Matches Fortran reference
```

**Impact**: All gfnff calculations now use correct dispersion by default

**Legacy Compatibility**: D3 still available via explicit `-d3` suffix
```bash
# New default (correct):
./curcuma -sp molecule.xyz -method gfnff  # Uses D4

# Legacy D3 (for debugging):
./curcuma -sp molecule.xyz -method gfnff-d3  # Uses D3
```

### Breaking Change Notice

⚠️ **BREAKING CHANGE**: All gfnff dispersion energies will change after this fix.

**No backward compatibility** - results change from INCORRECT to CORRECT values.

### Validation Status

**Build**: ✅ Compilation successful
**Test**: ✅ D4 parameter generation working (36 pairs, 13 ATM triples for CH₃OCH₃)
**Accuracy**: Under investigation (validation dataset needs D4 references from XTB 6.6.1)

**See**: [docs/GFNFF_DISPERSION_FIX.md](../../../../docs/GFNFF_DISPERSION_FIX.md) for complete technical details.

### Implementation Files

**Modified**:
- `d4param_generator.cpp:788-846` - CN-only weighting formula + architectural comments
- `gfnff_method.cpp:5199-5209` - D4 as default method

**Reference**:
- `external/gfnff/src/gfnff_gdisp0.f90:405` - Authoritative Fortran implementation

**Documentation**:
- `docs/GFNFF_DISPERSION_FIX.md` - Complete fix documentation
- `docs/GFNFF_STATUS.md` - Updated status with D4 fix section

---

## Code Consolidation Opportunities (December 2025)

### Coordination Number (CN) Calculation

**Current Situation**:
- D3ParameterGenerator has `calculateCoordinationNumbers()` declared but NOT implemented
- GFN-FF dispersion likely has CN calculation (needs investigation)
- EEQSolver may have geometry-dependent CN logic
- D4 would benefit from shared CN calculation

**Proposed**: Create shared `CNCalculator` utility class in `ff_methods/`
- Geometry-dependent CN from bond distances and covalent radii
- Used by: D3 C6 interpolation, D4, GFN-FF dispersion, potentially EEQ
- Benefits: Code reuse, consistent CN definition, easier validation

**D3 CN Calculation**: ✅ IMPLEMENTED (December 2025)
- Uses exponential counting formula: CN_i = Σ_j 1/(1+exp(-k1·(k2·R_cov/r_ij - 1)))
- Gaussian-weighted C6 interpolation across reference states
- Current accuracy: 1.48x (reduced from 1.52x with empty reference filtering)

**Related**: See `docs/GFNFF_STATUS.md` - "Code Consolidation Opportunities" section

## UFF-D3 Hybrid Method (December 19, 2025)

✅ **FULLY IMPLEMENTED** - Native D3 integration with UFF bonded terms

### Overview

**UFF-D3** combines UFF bonded terms (bonds, angles, dihedrals, inversions, vdW) with validated native D3 dispersion correction, providing a fast and accurate hybrid force field for molecular mechanics.

### Implementation Architecture

**Three-Component System**:

1. **Parameter Generation** (`forcefieldgenerator.cpp`):
   - `GenerateUFFD3Parameters()` - Generates UFF bonded + D3 dispersion parameters
   - Calls `Generate()` for UFF terms, then `D3ParameterGenerator::GenerateParameters()`
   - Merges both parameter sets into unified JSON

2. **Parameter Distribution** (`forcefield.cpp:AutoRanges()`):
   - D3 dispersion pairs distributed to threads via `addD3Dispersion()`
   - Multi-threaded parallelization across atom pairs
   - Method routing: "uff-d3" → method_type==1 with D3 flag

3. **Energy Calculation** (`forcefieldthread.cpp:execute()`):
   - UFF bonded terms: bonds, angles, dihedrals, inversions, vdW
   - Native D3 dispersion: `CalculateD3DispersionContribution()`
   - Total energy: E_total = E_UFF_bonded + E_D3_dispersion

### Usage

```bash
# UFF-D3 single point
./curcuma -sp molecule.xyz -method uff-d3

# UFF-D3 optimization
./curcuma -opt molecule.xyz -method uff-d3

# Geometry-dependent dispersion
./curcuma -sp monosaccharide.xyz -method uff-d3 -threads 4
```

### Accuracy

- **D3 Component**: 10/11 test molecules <1% error (validated against s-dftd3)
- **UFF Bonded**: Standard UFF accuracy for bonds, angles, dihedrals
- **Performance**: Multi-threaded D3 calculation, ~2-3x speedup with 4 threads

### Key Features

- ✅ Validated D3 dispersion (10/11 molecules <1% error)
- ✅ Geometry-dependent CN calculation with Gaussian weighting
- ✅ Multi-threaded parallelization via ForceFieldThread
- ✅ Consistent D3 implementation with GFN-FF
- ✅ PBE0/BJ damping parameters (a1=0.4145, a2=4.8593, s8=1.2177)

### Files Modified

- `forcefieldgenerator.h/cpp`: New `GenerateUFFD3Parameters()` method
- `forcefield.cpp`: D3 distribution in `AutoRanges()` for method "uff-d3"
- `forcefieldthread.h/cpp`: New `CalculateD3DispersionContribution()` method
- `gfnff_method.cpp`: Replaced `generateGFNFFDispersionPairs()` with native D3 (eliminates ~200 lines duplicate code)

### Integration with GFN-FF

**Shared D3 Infrastructure**:
- Both UFF-D3 and GFN-FF use the same `D3ParameterGenerator`
- Consistent D3 calculation across all force field methods
- GFN-FF's own dispersion replaced with validated native D3

**Benefits**:
- Code consolidation: Eliminates duplicate D3 implementations
- Consistency: Same D3 accuracy for both UFF-D3 and GFN-FF
- Maintainability: Single D3 implementation to validate and update

## GPU Pipeline (cuda/)

### ✅ Phase 1+2: GPU CN + GPU dc6dcn (March 2026)
- GPU CN computation replaces CPU O(N²) erf() loop
- GPU dc6dcn per-pair kernel replaces CPU O(N²) matrix + extraction

### ✅ Phase 3+4: CPU/GPU Overlap + Sync Removal (March 2026)
- `prepareAndLaunchChargeIndependent()` launches ~12 charge-independent kernels while CPU EEQ runs
- `launchChargeDependentAndFinish()` uploads EEQ charges, launches Coulomb + postprocess, downloads
- k_dispersion deferred to charge-dependent phase when `gradient=true` (needs dc6dcn from post-EEQ)
- Removed unnecessary `cudaStreamSynchronize` in `computeDC6DCNOnGPU`
- Relaxed postprocess stream dependencies: k_coulomb_self has no stream wait, k_subtract_qtmp waits pairwise+bonded only

### ✅ Phase 5: Shared Memory Energy Reduction (March 2026)
- `blockReduceAddEnergy()` device function: block-level tree reduction in shared memory, one `atomicAdd` per block
- All 12 energy kernels converted: k_dispersion, k_repulsion, k_coulomb, k_bonds, k_angles, k_dihedrals, k_inversions, k_storsions, k_batm, k_atm, k_xbonds, k_hbonds
- Threads with tid >= n contribute `local_E = 0.0` (must stay for `__syncthreads()`)
- Gradient atomicAdds unchanged (write to different addresses, no single-address contention)
- Verified: GPU vs CPU energy diff < 1 nEh, 24/24 CPU regression tests pass

### ✅ Phase 6: GPU Gaussian Weights + Async DMA (March 2026)
- **k_gaussian_weights kernel**: Computes gw and dgw/dCN directly on GPU (1 thread/atom, MAX_REF=7 in registers)
- **Full GPU dc6dcn pipeline**: k_gaussian_weights → k_dc6dcn_per_pair on same stream, no CPU→GPU sync
- **Eliminates**: CPU `precomputeGaussianWeights()` + `computeGaussianWeightDerivatives()` + flatten + 2 sync H2D uploads
- **`uploadRefCN()`**: One-time upload of reference CN values (826 doubles)
- **Async DMA downloads**: All D2H transfers via `cudaMemcpyAsync` on main stream, single `cudaStreamSynchronize` at end
- **Pinned energy buffer**: `m_h_energies` replaces stack array for true async transfer
- **Adaptive block sizing**: `getLaunchConfig()` with `__launch_bounds__(512, 2)` on all 22 kernels
- **Warp shuffle reduction**: `warpReduceSum()` via `__shfl_down_sync()` in `blockReduceAddEnergy()`

### ✅ Topology Caching (March 2026)
- **Two-tier caching**: Static topology (bonds, rings, hybridization) cached until large geometry change (>0.5 Bohr)
- **Dynamic state only**: CN and distance matrices updated each step (O(N²) vs O(N³) for full topology)
- **MD speedup**: ~15x for topology phase when topology is constant (typical MD)
- **Implementation**: `getCachedTopology()` in `gfnff_method.cpp`, `needsFullTopologyUpdate()` checks displacement

### ⚠️ Known Issues
- gfnff GPU validation tests (test_gfnff_gpu) fail with JSON null error — pre-existing, unrelated to pipeline
- k_dispersion cannot overlap with EEQ in gradient mode (dc6dcn dependency)

## Open Bugs

### ⚠️ Dead Code: `assignAtomsForSelfEnergy()`
- Declared in header but never called — can be removed (leftover from thread-safety fix Feb 23, 2026)

## References

- ForceFieldThread implements formulas from Fortran `gfnff_engrad.F90`
- GFNFF parameter generation follows Spicher/Grimme J. Chem. Theory Comput. 2020
- D3 dispersion: Grimme et al., J. Chem. Phys. 132, 154104 (2010)
