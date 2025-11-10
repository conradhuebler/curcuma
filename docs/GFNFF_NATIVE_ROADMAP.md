# GFN-FF Native C++ Implementation Roadmap
## Replacing Fortran Library with cgfnff as Default

**Goal**: Replace external Fortran GFN-FF library with native C++ implementation (`cgfnff`) for better control, maintainability, and educational value.

**Current Status**: ~50% complete - Bonds/Angles functional, critical terms missing

**Total Estimated Effort**: 8 phases, ~6-8 weeks full-time development

---

## **Architecture Vision**

```
┌─────────────────────────────────────────────────────────────┐
│                    NATIVE GFN-FF (cgfnff)                    │
├─────────────────────────────────────────────────────────────┤
│  ✅ Bonds (anharmonic)      ❌ Torsions (phase 1)           │
│  ✅ Angles (Fourier)        ❌ Inversions (phase 1)         │
│  ⚠️  Topology (basic)       ❌ EEQ Charges (phase 3)        │
│  ✅ Parameters Z=1-86       ❌ Ring Detection (phase 2)     │
│  ⚠️  D3/D4 (exists)         ❌ Non-bonded (phase 4)         │
├─────────────────────────────────────────────────────────────┤
│              ForceField Backend (reuse existing)             │
│         CurcumaLogger | ConfigManager | MethodFactory        │
└─────────────────────────────────────────────────────────────┘
```

**Design Principles**:
- **Educational First**: Code clarity > complex optimizations
- **Literature References**: Every formula documented with Paper equations
- **Incremental Validation**: Each phase tested against Fortran reference
- **Zero Fortran Dependencies**: Pure C++ with Eigen

---

## **Phase 1: Critical Energy Terms** ⚡ HIGHEST PRIORITY

**Duration**: 1-2 weeks
**Goal**: Implement missing torsion and inversion energy/gradient calculations

### 1.1 Torsion Implementation
**Reference**: `gfnff_engrad.F90:1041-1122` (`egtors` subroutine)

**Tasks**:
- [ ] Extract torsion parameters from Fortran code
- [ ] Implement `calculateTorsionEnergy(i,j,k,l)` in `gfnff.cpp`
- [ ] Implement `calculateTorsionGradient(i,j,k,l)` with analytical derivatives
- [ ] Add torsion parameter generation in `generateGFNFFParameters()`
- [ ] Integrate with `ForceFieldThread` (already has type=3 support)

**Formula** (from Paper):
```
E_tors = V_n/2 * [1 - cos(n*φ - φ0)] + topology_corrections
```

**Code Structure**:
```cpp
// In gfnff.cpp
struct GFNFFTorsionParams {
    double barrier;        // V_n from parameter table
    int periodicity;       // n (1,2,3)
    double phase;          // φ0
    double topology_fc;    // ring/conjugation corrections
};

json GFNFF::generateGFNFFTorsions() const {
    // 1. Detect all i-j-k-l torsion sequences
    // 2. For each: determine periodicity from hybridization
    // 3. Look up barrier from element types
    // 4. Apply ring/conjugation corrections
    // 5. Return JSON array for ForceField
}
```

**Validation**:
- Test against Fortran `egtors()` for ethane, butane, benzene
- Compare energies at φ = 0°, 60°, 120°, 180°
- Verify gradient with numerical differentiation

### 1.2 Inversion/Out-of-Plane Implementation
**Reference**: `gfnff_ini2.f90` (OOP setup) + `gfnff_engrad.F90` (energy)

**Tasks**:
- [ ] Implement `calculateInversionEnergy(center, i, j, k)`
- [ ] Implement `calculateInversionGradient()` with Wilson angle
- [ ] Add inversion detection (sp2 centers only)
- [ ] Parameter generation for carbonyls, imines, amides

**Formula**:
```
E_inv = k_inv * (θ_oop - θ0)^2
θ_oop = angle between bond and plane of other three atoms
```

**Code Structure**:
```cpp
struct GFNFFInversionParams {
    double force_constant;  // k_inv from hybridization
    double equilibrium;     // θ0 (usually 0° for planar)
};

json GFNFF::generateGFNFFInversions() const {
    // 1. Find all sp2 centers (C=O, C=N, N-planar, etc.)
    // 2. For each: identify three neighbors
    // 3. Calculate out-of-plane angle
    // 4. Assign force constant based on atom type
}
```

**Critical Cases**:
- Carbonyl groups (C=O)
- Amide groups (special treatment in angewChem2020_1)
- Aromatic nitrogens

### 1.3 ForceField Integration
**Files**: `src/core/energy_calculators/ff_methods/forcefield.cpp`

**Tasks**:
- [ ] Add `addGFNFFTorsion()` method
- [ ] Add `addGFNFFInversion()` method
- [ ] Ensure `setParameter()` routes torsions/inversions correctly
- [ ] Test with simple molecules (ethane, formaldehyde)

**Deliverables**:
- ✅ Torsion energy/gradient working for hydrocarbons
- ✅ Inversion energy/gradient working for carbonyls
- ✅ Tests pass against Fortran reference (±0.1 kcal/mol)

---

## **Phase 2: Topology Algorithms** 🔍

**Duration**: 1.5-2 weeks
**Goal**: Implement advanced topology detection for parameter assignment

### 2.1 Ring Detection
**Reference**: `gfnff_ini.f90:600-900` (ring detection algorithms)

**Tasks**:
- [ ] Implement `findSmallestRings()` - DFS/BFS based
- [ ] Detect 3-membered (strained), 4-membered, 5-membered (aromatic), 6-membered
- [ ] Store ring sizes per atom for parameter corrections
- [ ] Implement ring strain corrections for small rings

**Algorithm** (simplified):
```cpp
std::vector<int> GFNFF::findSmallestRings() const {
    // 1. Build adjacency graph from bonds
    // 2. For each atom, BFS to find shortest path back
    // 3. Store smallest ring containing each atom
    // 4. Special handling for fused rings

    std::vector<int> ring_sizes(m_atomcount, 0);
    // ... DFS/BFS implementation
    return ring_sizes;
}
```

**Validation**:
- Cyclohexane: all atoms in 6-ring
- Benzene: all atoms in 6-ring + aromatic flag
- Cubane: all atoms in 4-ring (strained)
- Adamantane: complex fused ring system

### 2.2 Pi-System Detection
**Reference**: `gfnff_ini.f90:1100-1300` (pi-system setup)

**Tasks**:
- [ ] Implement `detectPiSystems()` - conjugated fragment finder
- [ ] Identify sp2/sp hybridized atoms
- [ ] Connect into conjugated chains/rings
- [ ] Mark aromatic systems (Hückel 4n+2 rule)

**Code Structure**:
```cpp
std::vector<int> GFNFF::detectPiSystems(const std::vector<int>& hyb) const {
    std::vector<int> pi_fragments(m_atomcount, 0);

    // 1. Find all sp2/sp atoms
    // 2. Connect bonded sp2 atoms into fragments
    // 3. Check aromaticity for rings
    // 4. Return fragment IDs

    return pi_fragments;
}
```

**Critical Cases**:
- Benzene: single aromatic pi-system
- Butadiene: conjugated chain
- Pyridine: aromatic with heteroatom

### 2.3 Enhanced Hybridization Detection
**Reference**: `gfnff_ini.f90:400-550` (hybridization assignment)

**Current**: Simple neighbor counting
**Target**: Geometry-based detection

**Tasks**:
- [ ] Calculate bond angles around each atom
- [ ] Classify: linear (sp), trigonal (sp2), tetrahedral (sp3)
- [ ] Handle special cases: metals, hypervalent atoms
- [ ] Store hybridization for parameter lookup

**Algorithm**:
```cpp
std::vector<int> GFNFF::determineHybridization() const {
    // For each atom:
    // - If 2 neighbors + angle ~180°: sp
    // - If 3 neighbors + planar: sp2
    // - If 4 neighbors + tetrahedral: sp3
    // - Special cases: metals, P, S with expanded valence
}
```

**Deliverables**:
- ✅ Ring detection working for common ring sizes (3-8)
- ✅ Pi-system detection for conjugated molecules
- ✅ Hybridization accurate for 95% of organic molecules
- ✅ Topology info used in parameter assignment

---

## **Phase 3: EEQ Charge Calculation** ⚡ CRITICAL

**Duration**: 2-3 weeks
**Goal**: Replace placeholder charges with proper EEQ matrix solver

### 3.1 EEQ Theory Implementation
**Reference**: `gfnff_engrad.F90:1274-1400` (`goed_gfnff` subroutine)

**Background**: Extended Electronegativity Equalization
- Solves linear system: **A·q = b** for atomic charges
- A matrix contains: hardness (diagonal) + Coulomb terms (off-diagonal)
- Constraint: Σq = total_charge

**Tasks**:
- [ ] Implement `buildEEQMatrix()` - construct A matrix
- [ ] Implement constraint handling (Lagrange multiplier)
- [ ] Solve linear system with Eigen solvers
- [ ] Apply fragment-wise charge constraints
- [ ] Calculate EEQ energy contribution

**Matrix Structure**:
```cpp
Vector GFNFF::calculateEEQCharges(const Vector& cn,
                                   const std::vector<int>& hyb,
                                   const std::vector<int>& rings) const {
    int n = m_atomcount;
    Matrix A = Matrix::Zero(n+1, n+1);  // +1 for constraint
    Vector b = Vector::Zero(n+1);

    // Build EEQ matrix
    for (int i = 0; i < n; ++i) {
        // Diagonal: chemical hardness + CN correction
        double gam_i = getEEQParameters(i, {...}).gam;
        double cn_corr = cnf_angewChem2020[m_atoms[i]-1] * cn[i];
        A(i,i) = gam_i + cn_corr;

        // Off-diagonal: Coulomb interaction
        for (int j = i+1; j < n; ++j) {
            double r_ij = (m_geometry.row(i) - m_geometry.row(j)).norm();
            double gamma_ij = 1.0 / sqrt(r_ij*r_ij + /* damping */);
            A(i,j) = A(j,i) = gamma_ij;
        }

        // RHS: electronegativity + CN correction
        double chi_i = getEEQParameters(i, {...}).chi;
        b[i] = -chi_i - cn_corr;
    }

    // Charge constraint: Σq = total_charge
    for (int i = 0; i < n; ++i) {
        A(n, i) = A(i, n) = 1.0;
    }
    b[n] = m_charge;

    // Solve A·q = b
    Vector q_extended = A.colPivHouseholderQr().solve(b);
    return q_extended.head(n);  // Return charges without Lagrange multiplier
}
```

### 3.2 Coordination Number Derivatives
**Reference**: `gfnff_engrad.F90:802-853` (`dncoord_erf`)

**Why needed**: CN appears in many energy terms → need ∂E/∂CN for gradients

**Tasks**:
- [ ] Implement `calculateCoordinationNumberDerivatives()`
- [ ] Store as 3D tensor: dcn[xyz][i][j]
- [ ] Use in EEQ gradient calculation
- [ ] Use in bond/angle parameter derivatives

**Formula**:
```cpp
std::vector<Matrix> GFNFF::calculateCoordinationNumberDerivatives(
    const Vector& cn, double threshold) const {

    std::vector<Matrix> dcn(3, Matrix::Zero(m_atomcount, m_atomcount));

    for (int i = 0; i < m_atomcount; ++i) {
        for (int j = i+1; j < m_atomcount; ++j) {
            Vector r_ij = m_geometry.row(j) - m_geometry.row(i);
            double r = r_ij.norm();

            // dCN/dr from exponential decay function
            double rcov_sum = getCovalentRadius(m_atoms[i]) +
                             getCovalentRadius(m_atoms[j]);
            double arg = -16.0 * (r/rcov_sum - 1.0);
            double exp_val = exp(arg);
            double dcn_dr = 16.0 / rcov_sum * exp_val / pow(1.0 + exp_val, 2);

            // Chain rule: dCN/dx = dCN/dr * dr/dx
            for (int dim = 0; dim < 3; ++dim) {
                double dr_dx = r_ij[dim] / r;
                dcn[dim](i,j) = dcn_dr * dr_dx;
                dcn[dim](j,i) = -dcn[dim](i,j);  // Antisymmetric
            }
        }
    }
    return dcn;
}
```

### 3.3 EEQ Gradient Implementation
**Tasks**:
- [ ] Calculate ∂E_EEQ/∂q (charge derivatives)
- [ ] Calculate ∂q/∂xyz (geometry derivatives of charges)
- [ ] Combine with bond/angle gradient contributions
- [ ] Verify total gradient with numerical differentiation

**Deliverables**:
- ✅ EEQ charges agree with Fortran (±0.01 e)
- ✅ Dipole moments correct (±0.1 Debye)
- ✅ Charge gradients correct (numerical test)
- ✅ Tests: H2O, CH3OH, formamide, benzene

---

## **Phase 4: Non-Bonded Interactions** 🌐

**Duration**: 1-2 weeks
**Goal**: Integrate D3/D4 dispersion and repulsion terms

### 4.1 D4 Dispersion Integration
**Current**: D4 interface exists but not coupled to GFN-FF
**Reference**: `gfnff_engrad.F90` uses built-in D4 parameters

**Tasks**:
- [ ] Add D4 calculation call in `GFNFF::Calculation()`
- [ ] Use GFN-FF specific D4 parameters (not generic D4)
- [ ] Combine D4 gradient with bonded gradients
- [ ] Handle CN-dependent D4 coefficients

**Code Integration**:
```cpp
// In gfnff.cpp::Calculation()
double GFNFF::Calculation(bool gradient) {
    // 1. Calculate bonded terms (already working)
    double E_bonded = m_forcefield->Calculate(gradient);

    // 2. Calculate D4 dispersion
    double E_d4 = 0.0;
    if (m_parameters["dispersion"]) {
        // Use existing D4Interface with GFN-FF parameters
        json d4_config = {
            {"method", "gfnff"},
            {"cn", std::vector<double>(m_charges.data(),
                                       m_charges.data() + m_atomcount)}
        };
        // Call D4Thread with topology-aware coefficients
        E_d4 = calculateD4Energy(d4_config);
    }

    // 3. Calculate repulsion (short-range correction)
    double E_rep = calculateRepulsion();

    m_energy_total = convertToHartree(E_bonded + E_d4 + E_rep);
    return m_energy_total;
}
```

### 4.2 Repulsion Term
**Reference**: `gfnff_engrad.F90` (repulsion damping)

**Formula**:
```
E_rep = Σ_ij Z_i*Z_j / r_ij * damp(r_ij)
```

**Tasks**:
- [ ] Implement short-range repulsion
- [ ] Add damping function (avoid singularity)
- [ ] Calculate repulsion gradient
- [ ] Test on ionic systems (NaCl, MgO)

### 4.3 Hydrogen Bond Correction
**Reference**: `gfnff_engrad.F90:725-800` (`egbond_hb`)

**Tasks**:
- [ ] Detect A-H...B hydrogen bonds
- [ ] Calculate H-bond energy with angular dependence
- [ ] Add H-bond gradient
- [ ] Test on water dimer, formamide dimer

**Deliverables**:
- ✅ D4 energy matches Fortran (±0.1 kcal/mol)
- ✅ Repulsion prevents unrealistic geometries
- ✅ H-bonds correct for water clusters
- ✅ Total non-bonded energy validated

---

## **Phase 5: Parameter Corrections** 🔧

**Duration**: 1 week
**Goal**: Fix problematic parameter assignment strategies

### 5.1 Equilibrium Values from Tables
**Current Problem**:
```cpp
// WRONG: Uses current geometry
params.equilibrium_distance = distance;
```

**Solution**: Look up from GFN-FF reference structures

**Tasks**:
- [ ] Extract equilibrium bond lengths from Fortran parameters
- [ ] Create `getEquilibriumBondLength(z1, z2, bond_order)` lookup
- [ ] Create `getEquilibriumAngle(z1, z2, z3, hybridization)` lookup
- [ ] Apply topology corrections (ring strain, conjugation)

**New Code**:
```cpp
GFNFF::GFNFFBondParams GFNFF::getGFNFFBondParameters(
    int z1, int z2, double distance) const {

    GFNFFBondParams params;

    // Force constant (already correct)
    double bond_param_1 = bond_angewChem2020[z1-1];
    double bond_param_2 = bond_angewChem2020[z2-1];
    params.force_constant = sqrt(bond_param_1 * bond_param_2);

    // Equilibrium distance from TABULATED VALUES (not current!)
    params.equilibrium_distance = getEquilibriumBondLength(z1, z2);

    // Apply hybridization corrections
    int hyb1 = m_hybridization[...];  // from topology
    if (hyb1 == 2) {  // sp2: shorter bond
        params.equilibrium_distance *= 0.95;
    }

    // Apply ring strain corrections
    int ring_size = m_ring_sizes[...];
    if (ring_size == 3) {
        params.force_constant *= 1.2;  // Stiffer in small rings
    }

    params.anharmonic_factor = getAnharmonicity(z1, z2);
    return params;
}
```

### 5.2 Scaling Factor Validation
**Tasks**:
- [ ] Verify angle force constant scaling (currently `* 0.001`)
- [ ] Compare against Fortran `egbend()` force constants
- [ ] Document why each scaling factor is needed
- [ ] Add unit tests for parameter consistency

### 5.3 Topology-Aware Parameters
**Reference**: `gfnff_ini2.f90` (parameter modification based on topology)

**Tasks**:
- [ ] Implement `applyTopologyCorrections()` for bonds
- [ ] Implement `applyTopologyCorrections()` for angles
- [ ] Use CN, hybridization, ring info for corrections
- [ ] Test on strained systems (cubane, cyclopropane)

**Deliverables**:
- ✅ Bond lengths match experimental ±0.05 Å
- ✅ Angles match experimental ±5°
- ✅ Force constants reproduce vibrational frequencies
- ✅ Tests: C-C, C=C, C≡C bonds in different environments

---

## **Phase 6: Validation & Testing** ✅

**Duration**: 1-2 weeks
**Goal**: Comprehensive validation against Fortran reference

### 6.1 Test Suite Development
**Files**: Create `test_cases/gfnff_validation/`

**Test Categories**:

#### Energy Tests
```bash
test_cases/gfnff_validation/
├── 01_hydrocarbons/
│   ├── methane.xyz         # Simple sp3
│   ├── ethene.xyz          # sp2, double bond
│   ├── ethyne.xyz          # sp, triple bond
│   ├── butane.xyz          # Torsion test
│   ├── cyclohexane.xyz     # Ring test
│   └── expected_energies.json
├── 02_heteroatoms/
│   ├── water.xyz
│   ├── ammonia.xyz
│   ├── methanol.xyz
│   └── expected_energies.json
├── 03_conjugated/
│   ├── benzene.xyz         # Aromatic
│   ├── naphthalene.xyz     # Fused rings
│   ├── butadiene.xyz       # Conjugated chain
│   └── expected_energies.json
├── 04_strained/
│   ├── cyclopropane.xyz    # Ring strain
│   ├── cubane.xyz          # Extreme strain
│   ├── norbornane.xyz      # Bridged
│   └── expected_energies.json
└── 05_functional_groups/
    ├── formaldehyde.xyz    # C=O
    ├── formamide.xyz       # Amide
    ├── acetic_acid.xyz     # Carboxylic acid
    └── expected_energies.json
```

#### Gradient Tests
- [ ] Numerical vs. analytical gradient comparison
- [ ] Geometry optimization convergence
- [ ] Force constant matrix (Hessian) validation

#### Property Tests
- [ ] Atomic charges vs. Fortran EEQ
- [ ] Dipole moments vs. experimental
- [ ] Bond orders (Wiberg) validation

### 6.2 Reference Data Generation
**Script**: `scripts/generate_gfnff_reference.sh`

```bash
#!/bin/bash
# Generate reference energies with external GFN-FF library

for molecule in test_cases/gfnff_validation/*/*.xyz; do
    echo "Processing $molecule"

    # External Fortran GFN-FF
    ./curcuma -sp $molecule -method gfnff > ${molecule%.xyz}_extern.out

    # Native C++ cgfnff
    ./curcuma -sp $molecule -method cgfnff > ${molecule%.xyz}_native.out

    # Extract energies and compare
    python scripts/compare_energies.py ${molecule%.xyz}_extern.out \
                                       ${molecule%.xyz}_native.out
done
```

### 6.3 Accuracy Benchmarks
**Targets**:
- Energy: ±0.5 kcal/mol vs. Fortran
- Gradient: ±0.001 Hartree/Bohr vs. Fortran
- Charges: ±0.02 e vs. Fortran
- Geometry optimization: Same minimum structure

**Validation Report**:
```
# test_cases/gfnff_validation/VALIDATION_REPORT.md

## Native cgfnff vs. External gfnff Comparison

| Test Set          | Molecules | ΔE (kcal/mol) | ΔGrad (%) | Pass |
|-------------------|-----------|---------------|-----------|------|
| Hydrocarbons      | 5         | 0.32 ± 0.15   | 1.2       | ✅   |
| Heteroatoms       | 3         | 0.45 ± 0.22   | 2.1       | ✅   |
| Conjugated        | 3         | 0.68 ± 0.31   | 3.5       | ⚠️   |
| Strained Rings    | 3         | 1.12 ± 0.55   | 5.8       | ❌   |
| Functional Groups | 3         | 0.51 ± 0.19   | 2.4       | ✅   |

**Overall Pass Rate**: 80% (16/20 molecules within tolerance)
```

**Deliverables**:
- ✅ 20+ test molecules with reference data
- ✅ Automated validation script
- ✅ CI/CD integration (ctest)
- ✅ Validation report documenting accuracy

---

## **Phase 7: Performance Optimization** ⚡

**Duration**: 1 week
**Goal**: Ensure native implementation is competitive with Fortran

### 7.1 Profiling
**Tools**: `gprof`, `perf`, `valgrind --tool=callgrind`

**Tasks**:
- [ ] Profile on large molecule (500+ atoms)
- [ ] Identify bottlenecks (likely EEQ matrix solve)
- [ ] Compare timing vs. Fortran library
- [ ] Set performance target: ≤2x slower than Fortran

### 7.2 Optimizations
**Strategies**:
- [ ] **Neighbor lists**: Avoid O(N²) distance calculations
- [ ] **Matrix operations**: Use Eigen's optimized routines
- [ ] **Threading**: Parallelize bond/angle calculations
- [ ] **Caching**: Store topology (only recalculate on bond breaking)

**Example Optimization**:
```cpp
// Before: O(N²) bond search every step
for (int i = 0; i < N; ++i) {
    for (int j = i+1; j < N; ++j) {
        if (distance(i,j) < threshold) {
            // Add bond
        }
    }
}

// After: Neighbor list (O(N) update)
class NeighborList {
    std::vector<std::vector<int>> neighbors;
    double cutoff = 10.0;  // Å

    void update(const Matrix& geometry) {
        // Only check atoms within cutoff
        // Use cell lists for O(N) scaling
    }
};
```

### 7.3 Memory Optimization
**Tasks**:
- [ ] Avoid unnecessary matrix copies
- [ ] Use `Eigen::Ref<>` for zero-copy function arguments
- [ ] Profile memory usage (valgrind massif)
- [ ] Set memory target: ≤1.5x Fortran memory

**Deliverables**:
- ✅ Benchmark suite (small/medium/large molecules)
- ✅ Performance report comparing Fortran vs. C++
- ✅ Identified optimizations yield ≥30% speedup
- ✅ Memory usage documented

---

## **Phase 8: Default Integration & Documentation** 📚

**Duration**: 1 week
**Goal**: Make `cgfnff` the default, deprecate Fortran library

### 8.1 MethodFactory Update
**File**: `src/core/energy_calculators/method_factory.cpp`

**Change Priority**:
```cpp
// OLD: External > XTB > Native
"gfnff": External GFN-FF > XTB > Native cgfnff

// NEW: Native > External > XTB
"gfnff": Native cgfnff > External GFN-FF > XTB
```

**Tasks**:
- [ ] Update `MethodFactory` priority list
- [ ] Add deprecation warning for external library
- [ ] Test fallback mechanism still works
- [ ] Update CMake default: `USE_GFNFF=OFF` (native always available)

### 8.2 Documentation
**Files to Update**:

#### CLAUDE.md
```markdown
### Native GFN-FF (cgfnff) - DEFAULT METHOD ✅
- **Status**: Production-ready native C++ implementation
- **Features**: Complete GFN-FF functionality including:
  - ✅ Bonds, Angles, Torsions, Inversions
  - ✅ EEQ Charge Calculation
  - ✅ D3/D4 Dispersion
  - ✅ Hydrogen Bonding
  - ✅ Topology-aware parameters
- **Performance**: ~1.5x slower than Fortran (acceptable tradeoff)
- **Validation**: 95% agreement with reference implementation
- **Advantages**:
  - No Fortran dependencies
  - Educational code clarity
  - Full integration with Curcuma ecosystem
  - Easy debugging and extension
```

#### docs/GFNFF_NATIVE.md (NEW)
```markdown
# Native GFN-FF Implementation (cgfnff)

## Overview
Curcuma's native C++ implementation of the GFN-FF force field
(Spicher & Grimme, Angew. Chem. 2020).

## Usage
```bash
# Default method (uses native implementation)
./curcuma -sp molecule.xyz -method gfnff

# Explicit native call
./curcuma -opt molecule.xyz -method cgfnff
```

## Architecture
[Detailed architecture documentation]

## Implementation Details
[Reference to code sections with explanations]

## Validation
[Link to test suite and accuracy reports]

## Literature References
[Complete citation list]
```

#### README.md
```markdown
### Quantum Mechanical Methods
- **GFN-FF (default)**: Native C++ implementation - No dependencies
- **GFN2-xTB**: Via TBLite > Ulysses > XTB (priority-based)
- **Extended Hückel Theory**: Native implementation
- **UFF**: Universal Force Field

### Force Field Methods
- **GFN-FF (cgfnff)**: Geometry/Frequency/Noncovalent Force Field
- **UFF**: Universal Force Field
- **QMDFF**: Quantum Mechanically Derived Force Fields
```

### 8.3 Build System
**CMakeLists.txt Changes**:
```cmake
# Native GFN-FF is always compiled (no flag needed)
set(curcuma_gfnff_native_SRC
    src/core/energy_calculators/qm_methods/gfnff.cpp
    src/core/energy_calculators/qm_methods/gfnff_advanced.cpp
    src/core/energy_calculators/qm_methods/gfnff_method.cpp
)
add_library(curcuma_gfnff_native ${curcuma_gfnff_native_SRC})
target_link_libraries(curcuma_core curcuma_gfnff_native)

# External Fortran library (optional, deprecated)
option(USE_GFNFF_EXTERNAL "Use external Fortran GFN-FF library (deprecated)" OFF)
if(USE_GFNFF_EXTERNAL)
    message(WARNING "External GFN-FF library is deprecated. Native cgfnff is now default.")
    add_subdirectory(external/gfnff)
    target_compile_definitions(curcuma_core PRIVATE USE_GFNFF_EXTERNAL)
endif()
```

### 8.4 Testing
**CI/CD Integration**:
```yaml
# .github/workflows/test.yml
- name: Test Native GFN-FF
  run: |
    cd build
    ctest -R "gfnff_native_" --output-on-failure

- name: Validate vs External (if available)
  run: |
    if [ -d external/gfnff/src ]; then
      python scripts/validate_gfnff.py
    fi
```

**Deliverables**:
- ✅ `cgfnff` is default for `gfnff` method name
- ✅ External library optional (backward compatibility)
- ✅ Complete documentation updated
- ✅ CI tests pass with native implementation

---

## **Maintenance & Future Work** 🔮

### Ongoing Maintenance
- [ ] Monitor for GFN-FF parameter updates (check xtb repo quarterly)
- [ ] Add new elements if coverage expands beyond Z=86
- [ ] Performance tuning based on user feedback

### Future Enhancements
- [ ] **ALPB/GBSA Solvation**: Port solvent models
- [ ] **Periodic Boundaries**: For crystals/surfaces
- [ ] **Analytical Hessian**: Second derivatives for frequencies
- [ ] **GPU Acceleration**: CUDA/OpenCL for large systems
- [ ] **Machine Learning**: Parameter refinement from QM data

---

## **Success Criteria** ✅

**Phase 1-4 (Functional)**:
- ✅ Energy calculation complete (all terms)
- ✅ Gradient calculation complete (analytical)
- ✅ Tests pass for 20+ molecules

**Phase 5-6 (Accurate)**:
- ✅ Energy within 0.5 kcal/mol of Fortran
- ✅ Gradients within 1% of Fortran
- ✅ Optimized geometries identical

**Phase 7 (Performant)**:
- ✅ Runtime ≤2x Fortran library
- ✅ Memory ≤1.5x Fortran library
- ✅ Scales to 1000+ atom systems

**Phase 8 (Production)**:
- ✅ Default method in MethodFactory
- ✅ Documentation complete
- ✅ CI/CD tests passing
- ✅ No Fortran dependencies required

---

## **Risk Mitigation**

### Technical Risks
| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| EEQ solver numerical instability | Medium | High | Use robust Eigen solvers, add regularization |
| Gradient validation fails | Medium | High | Extensive numerical differentiation tests |
| Performance unacceptable | Low | Medium | Profile early, optimize incrementally |
| Parameter incompatibility | Low | High | Validate each parameter against Fortran |

### Schedule Risks
| Risk | Mitigation |
|------|------------|
| Phase overruns | Prioritize phases 1-4, defer optimizations |
| Validation failures | Allocate 50% time buffer in Phase 6 |
| Scope creep | Freeze features after Phase 4, document future work |

---

## **Resources & References**

### Primary Literature
1. **Spicher & Grimme** (2020). "Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems" *Angew. Chem. Int. Ed.* 59, 15665
2. **Grimme et al.** (2017). "A robust and accurate tight-binding quantum chemical method for structures, vibrational frequencies, and noncovalent interactions" *J. Chem. Theory Comput.* 13, 1989

### Code References
- **External Implementation**: `external/gfnff/src/` (Fortran reference)
- **XTB Source**: https://github.com/grimme-lab/xtb (original implementation)
- **Current Native**: `src/core/energy_calculators/qm_methods/gfnff.cpp`

### Tools
- **Eigen**: Matrix operations and solvers
- **CurcumaLogger**: Verbosity control
- **ConfigManager**: Parameter management
- **ForceField**: Energy/gradient calculation backend

---

## **Timeline Summary**

```
Week 1-2:   Phase 1 - Torsions & Inversions
Week 3-4:   Phase 2 - Topology Algorithms
Week 5-7:   Phase 3 - EEQ Charges (complex!)
Week 8-9:   Phase 4 - Non-bonded Terms
Week 10:    Phase 5 - Parameter Fixes
Week 11-12: Phase 6 - Validation
Week 13:    Phase 7 - Optimization
Week 14:    Phase 8 - Integration & Docs

Total: 14 weeks (3.5 months) conservative estimate
```

**Fast Track** (if focused): 6-8 weeks by parallelizing phases and accepting initial performance overhead

---

## **Getting Started**

### Immediate Next Steps
1. **Read External Code**: Study `gfnff_engrad.F90` to understand energy calculations
2. **Test Current State**: Run existing `cgfnff` on simple molecules, document failures
3. **Start Phase 1**: Implement torsions (highest value/effort ratio)
4. **Set Up Validation**: Generate reference data from external library

### Development Workflow
```bash
# 1. Implement feature in gfnff.cpp
vim src/core/energy_calculators/qm_methods/gfnff.cpp

# 2. Build
cd build && make -j4

# 3. Test against reference
./test_gfnff_feature.sh

# 4. Validate
python scripts/compare_to_fortran.py

# 5. Document
vim docs/GFNFF_NATIVE_ROADMAP.md  # Update progress
```

---

**End of Roadmap** - Ready to become the reference GFN-FF implementation in C++! 🚀
