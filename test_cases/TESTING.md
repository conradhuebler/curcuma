# Curcuma Testing Framework Guide

Session-ready documentation for testing all computational methods in Curcuma.

---

## Quick Start

### Run All Tests
```bash
cd /path/to/curcuma/build/test_cases
ctest --output-on-failure
```

### Run Specific Method Tests
```bash
# GFN-FF tests
ctest -R "gfnff" --output-on-failure

# Energy method tests
ctest -R "energy" --output-on-failure

# RMSD tests
ctest -R "rmsd" --output-on-failure
```

### Running External Verification
```bash
# Unit Tests use hardcoded molecules (no paths needed)
./test_gfnff_unified

# External tools need XYZ files
xtb test_cases/molecules/larger/CH4.xyz --gfnff
./build/external/simple-dftd3/release/app/s-dftd3 input.xyz
```

---

## Two-Tier Molecule Strategy

| Use Case | Source | Path Required? |
|----------|--------|----------------|
| **Unit Tests (C++)** | TestMoleculeRegistry (hardcoded) | ❌ No - self-contained, path-independent |
| **External Programs** | XYZ files in molecules/ | ✅ Yes - for xtb, s-dftd3, etc. |

**Why Hardcoded?** Unit tests run anywhere without file system dependencies. Tests are portable and execute in any build directory. External CLI tools (xtb, dftd3) require XYZ file input.

---

## Universal Molecule Registry

**Location**: `test_cases/core/test_molecule_registry.h`

### Available Molecules

#### Dimers (2 atoms)
| Name | Registry Key | XYZ File | Atoms | Tests |
|------|--------------|----------|-------|-------|
| Hydrogen | `"H2"` | dimers/HH.xyz | H-H | Bond energy, all methods |
| Hydrogen Chloride | `"HCl"` | dimers/HCl.xyz | Cl-H | Heteronuclear bond |
| Hydroxyl Radical | `"OH"` | dimers/OH.xyz | O-H | Radical, charge |

#### Trimers (3 atoms)
| Name | Registry Key | XYZ File | Atoms | Tests |
|------|--------------|----------|-------|-------|
| Water Trimer | `"H2O"` | trimers/water.xyz | O-H-H | Angle energy, electrostatic |
| Hydrogen Cyanide | `"HCN"` | trimers/HCN.xyz | H-C-N | Linear molecule |
| Ozone | `"O3"` | trimers/O3.xyz | O-O-O | Ozone test |

#### Larger (5+ atoms)
| Name | Registry Key | XYZ File | Atoms | Tests |
|------|--------------|----------|-------|-------|
| Methane | `"CH4"` | larger/CH4.xyz | C + 4H | Multiple equivalent bonds, all methods |
| Methanol | `"CH3OH"` | larger/CH3OH.xyz | C₁O₁H₆ | All energy terms + torsion |
| Dimethyl Ether | `"CH3OCH3"` | larger/CH3OCH3.xyz | C₂O₁H₆ | Complex torsions |
| Benzene | `"C6H6"` / `"benzene"` | larger/C6H6.xyz | C₆H₆ | Aromatic ring, conjugated systems |
| Monosaccharide | `"monosaccharide"` | larger/monosaccharide.xyz | 27 atoms | Medium-large system |
| Triose | `"triose"` | larger/triose.xyz | 66 atoms | Large system, performance test |

### Using TestMoleculeRegistry

```cpp
#include "test_cases/core/test_molecule_registry.h"
using namespace TestMolecules;

// Create molecule - path-independent! (Unit Tests)
curcuma::Molecule mol = TestMoleculeRegistry::createMolecule("CH4");

// Get XYZ path for external verification
std::string xyz_path = TestMoleculeRegistry::getXyzPath("CH4");
// Returns: "molecules/larger/CH4.xyz"
```

### Check Available Molecules

```cpp
// Get all molecules
auto all_mols = TestMoleculeRegistry::getAllMoleculeNames();
// Returns: {"H2", "HCl", "OH", "CH4", "CH3OH", "CH3OCH3", "C6H6", "HCN", "H2O", "O3", "monosaccharide", "triose"}

// Get by category
auto dimers = TestMoleculeRegistry::getDimers();
auto trimers = TestMoleculeRegistry::getTrimers();
auto larger = TestMoleculeRegistry::getLargerMolecules();
```

---

## Adding New Molecules

### Step 1: Create XYZ File
```bash
# Create file in appropriate category
vim test_cases/molecules/<category>/<name>.xyz
```

XYZ format:
```
<number of atoms>
<comment line (element count, description)>
<element> <x> <y> <z>
...
```

### Step 2: Generate Reference Energies

**For QM Methods (GFN1, GFN2, etc.):**
```bash
# Using XTB (requires xtb in PATH)
xtb molecule.xyz --gfn1 > molecule_gfn1.out
xtb molecule.xyz --gfn2 > molecule_gfn2.out
```

**For Force Fields (UFF, GFN-FF, QMDFF):**
```bash
# Using Curcuma
./curcuma -sp molecule.xyz -method uff > molecule_uff.out
./curcuma -sp molecule.xyz -method cgfnff > molecule_gfnff.out
```

**For Dispersion (D3, D4):**
```bash
# Using s-dftd3 (if available)
./build/external/simple-dftd3/release/app/s-dftd3 molecule.xyz
```

### Step 3: Add to test_molecule_registry.cpp

Edit template (add to `s_molecule_registry` map):
```cpp
{
    "NewMol", {
        .name = "NewMol",
        .description = "Description from XYZ file",
        .category = "dimers|trimers|larger",
        .atoms = {
            {atomic_number, Eigen::Vector3d(x, y, z)},
            // ... all atoms from XYZ file
        },
        .reference_energies = {
            {"gfn1", total_energy},      // XTB GFN1
            {"gfn2", total_energy},      // XTB GFN2
            {"uff", total_energy},      // Curcuma UFF
            {"cgfnff", total_energy},   // Curcuma native GFN-FF
            {"d3", d3_energy},          // D3 dispersion only
        },
        .tolerances = {
            {"gfn1", 1e-6},            // Hartree tolerance
            {"gfn2", 1e-6},
            {"uff", 1e-5},            // Force field tolerance
            {"cgfnff", 1e-6},
            {"d3", 1e-8},             // D3 tolerance
        },
        .atom_count = N
    }
},
```

### Step 4: Add XYZ Path Mapping

Edit `s_xyz_paths` map in test_molecule_registry.cpp:
```cpp
{"NewMol", "molecules/<category>/NewMol.xyz"},
```

### Step 5: Update Test Suite

Add reference data to appropriate test file:
- `test_gfnff_unified.cpp` for GFN-FF energy components
- `test_energy_methods.cpp` for general method validation
- Or create new test file

### Step 6: Register Test in CMakeLists.txt

```cmake
add_executable(test_newmethod test_newmethod.cpp)
target_link_libraries(test_newmethod curcuma_core test_molecule_registry)
add_test(NAME test_newmethod COMMAND test_newmethod)
```

### Step 7: Run Tests
```bash
make test_newmethod
ctest -R "newmethod" --verbose
```

---

## Test Suite Architecture

```
test_cases/
├── TESTING.md                              ← This file (general blueprint)
│
├── core/
│   └── test_molecule_registry.h/.cpp        ← UNIVERSAL molecule library
│       - 12 molecules (dimers, trimers, larger)
│       - Hardcoded coordinates (path-independent)
│       - XYZ path mapping (for external tools)
│       - Reference data storage
│
├── molecules/                               ← XYZ files for external programs
│   ├── dimers/
│   │   ├── HH.xyz, HCl.xyz, OH.xyz
│   ├── trimers/
│   │   ├── HCN.xyz, water.xyz, O3.xyz
│   └── larger/
│       ├── CH4.xyz, CH3OH.xyz, CH3OCH3.xyz, C6H6.xyz
│       ├── monosaccharide.xyz, triose.xyz
│
├── test_gfnff_unified.cpp                   ← GFN-FF energy components
├── test_gfnff_d3.cpp                        ← D3 integration validation
├── test_energy_methods.cpp                  ← General method validation
├── test_molecule.cpp                        ← Molecule data structure tests
│
├── AAAbGal.cpp                              ← RMSD integration test
├── confscan.cpp                             ← ConfScan integration test
├── rmsd/                                    ← RMSD unit tests
│
└── cli/                                     ← End-to-end tests
    ├── rmsd/               ← RMSD CLI tests (6/6 passing)
    ├── confscan/           ← ConfScan CLI tests (7/7 passing)
    ├── curcumaopt/         ← Optimization CLI tests (6/6 passing)
    └── simplemd/           ← MD CLI tests (0/7 - JSON issue)
```

---

## Creating New Test Files

### 1. Unit Test Template

```cpp
/**
 * Test Name
 * Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * Purpose: Brief description of what this test validates
 */

#include <cassert>
#include <iostream>
#include <cmath>

#include "test_cases/core/test_molecule_registry.h"
#include "src/core/molecule.h"
#include "src/core/energycalculator.h"
#include "src/core/config_manager.h"

using namespace TestMolecules;

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Test: [Test Name]" << std::endl;
    std::cout << "========================================" << std::endl;

    // Create molecule using TestMoleculeRegistry
    Molecule mol = TestMoleculeRegistry::createMolecule("CH4");

    // Initialize calculator
    EnergyCalculator calc;
    calc.setMolecule(mol);

    // Run test
    double tolerance = 1e-6;
    double energy = calc.calculateEnergy();

    // Validate
    const MoleculeData& data = TestMoleculeRegistry::getMolecule("CH4");
    double reference = data.reference_energies.at("your_method");

    if (std::abs(energy - reference) < tolerance) {
        std::cout << "✅ Test PASSED" << std::endl;
        return 0;
    } else {
        std::cout << "❌ Test FAILED" << std::endl;
        std::cout << "   Expected: " << reference << " Eh" << std::endl;
        std::cout << "   Got:      " << energy << " Eh" << std::endl;
        return 1;
    }
}
```

### 2. Register in CMakeLists.txt

```cmake
add_executable(test_newmethod test_newmethod.cpp)
target_link_libraries(test_newmethod curcuma_core test_molecule_registry)

if(USE_PCH)
    target_precompile_headers(test_newmethod REUSE_FROM curcuma_core)
endif()

add_test(NAME test_newmethod COMMAND test_newmethod)
set_tests_properties(test_newmethod PROPERTIES TIMEOUT 60)
```

### 3. Run Test

```bash
make test_newmethod
ctest -R "test_newmethod" --verbose
```

---

## Method-Specific Testing

### GFN-FF (Native Force Field)

**Test Files:**
- `test_gfnff_unified.cpp` - Energy component validation
- `test_gfnff_d3.cpp` - D3 dispersion integration

**Energy Components Tested:**
- Bond stretching
- Angle bending
- Torsion
- Repulsion (bonded + non-bonded)
- Electrostatic (EEQ charges)
- Dispersion (D4)

**Usage:**
```cpp
Molecule mol = TestMoleculeRegistry::createMolecule("CH4");

GFNFF gfnff;
gfnff.setParameter("verbosity", 1);  // Control output
gfnff.Initialise(mol);
double energy = gfnff.CalculateEnergy(true, true);  // Energy, gradient
```

**Validation Reference:** XTB GFN-FF (xtb --gfnff)

### Energy Methods (GFN1, GFN2, etc.)

**Test File:**
- `test_energy_methods.cpp` - General method validation

**Available Methods:**
- QM: gfn1, gfn2, gfn0, ipea1, pm3, pm6, eht
- FF: uff, qmdff, cgfnff
- Providers: tblite, ulysses, xtb

**Usage:**
```cpp
Molecule mol = TestMoleculeRegistry::createMolecule("CH4");

json config = {{"method", "gfn2"}};
EnergyCalculator calc(config);
calc.setMolecule(mol);
double energy = calc.calculateEnergy();
```

**Validation Reference:** XTB, TBLite, Ulysses (respectively)

### UFF (Universal Force Field)

**Usage:**
```cpp
Molecule mol = TestMoleculeRegistry::createMolecule("CH4");

json config = {{"method", "uff"}};
EnergyCalculator calc(config);
calc.setMolecule(mol);
double energy = calc.calculateEnergy();
```

**Validation Reference:** UFF literature values, XTB UFF

### D3/D4 Dispersion

**Test Files:**
- `test_gfnff_d3.cpp` - D3 integration
- `test_dispersion.cpp` - General dispersion testing

**Usage:**
```cpp
Molecule mol = TestMoleculeRegistry::createMolecule("CH4");

D3ParameterGenerator d3_gen(config);
d3_gen.GenerateParameters(atoms, geometry);
double d3_energy = d3_gen.getTotalEnergy();
```

**Validation Reference:** s-dftd3 standalone (simple-dftd3)

---

## Debugging Failed Tests

### 1. Check Molecule Coordinates

```cpp
Molecule mol = TestMoleculeRegistry::createMolecule("CH4");
mol.print();  // Verify coordinates
```

### 2. Verify Against External Reference

```bash
# Get XYZ path
std::string xyz_path = TestMoleculeRegistry::getXyzPath("CH4");

# Run external reference
cd build/test_cases
xtb $xyz_path --gfn1    # For QM method
xtb $xyz_path --gfnff   # For force field
```

### 3. Extract Energy Components

```cpp
// Enable verbose output
json config = {{"method", "cgfnff"}, {"verbosity", 3}};
EnergyCalculator calc(config);
calc.setMolecule(mol);
calc.calculateEnergy();
// Look for: BondEnergy, AngleEnergy, Repulsion, etc.
```

### 4. Compare Reference Energies

```cpp
const MoleculeData& data = TestMoleculeRegistry::getMolecule("CH4");
if (data.hasReferenceEnergy("gfn2")) {
    double ref = data.reference_energies.at("gfn2");
    double tol = data.tolerances.at("gfn2");
    // Validate calculated value
    ASSERT_NEAR(calc_energy, ref, tol);
}
```

### 5. Gradient Validation

```cpp
EnergyCalculator calc;
calc.setMolecule(mol);
calc.setCalculateGradient(true);
calc.calculateEnergy();

Matrix grad = calc.getGradient();
// Compare with numerical gradient
```

---

## Reference Data Sources

| Method | Reference Tool | Generation Date |
|--------|----------------|-----------------|
| GFN1 | XTB 6.6.1 | Current |
| GFN2 | XTB 6.6.1 | Current |
| GFN-FF | XTB 6.6.1 (--gfnff) | Current |
| UFF | Curcuma (--method uff) | Current |
| D3 | s-dftd3 standalone | Dec 2025 |
| D4 | XTB (--gfn2) + D4 | Dec 2025 |

**Updating References:**
```bash
# For each method, generate reference energies
xtb molecule.xyz --gfn1 > molecule_gfn1_ref.out
xtb molecule.xyz --gfn2 > molecule_gfn2_ref.out
./curcuma -sp molecule.xyz -method uff > molecule_uff_ref.out
./curcuma -sp molecule.xyz -method cgfnff > molecule_cgfnff_ref.out
```

---

## Test Tolerances

| Method Type | Typical Tolerance | Notes |
|-------------|-------------------|-------|
| **QM Methods** | 1e-6 Eh | High precision (semi-empirical) |
| **Force Fields** | 1e-5 Eh | Lower precision (empirical) |
| **Dispersion (D3)** | 1e-8 Eh | Very precise (single term) |
| **Dispersion (D4)** | 1e-8 Eh | Very precise (single term) |
| **Energy Components** | 1e-6 Eh | Individual terms in GFN-FF |
| **Gradients** | 1e-6 Eh/Bohr | Analytical vs numerical |

---

## CLI End-to-End Tests

**Location**: `test_cases/cli/`

**Architecture**: Tests execute from BUILD TREE, keeping SOURCE TREE clean

**Test Categories**:
- `rmsd/` - RMSD calculations (6/6 passing ✅)
- `confscan/` - Conformational searching (7/7 passing ✅)
- `curcumaopt/` - Geometry optimization (6/6 passing ✅)
- `simplemd/` - Molecular dynamics (0/7 - JSON issue)

**Running CLI Tests**:
```bash
# All CLI tests
ctest -R "cli_" --output-on-failure

# Specific category
ctest -R "cli_rmsd_" --output-on-failure

# Individual test
ctest -R "cli_rmsd_01" --verbose
```

**CLI Test Pattern**:
```bash
#!/bin/bash
source ../../test_utils.sh

run_test() {
    $CURCUMA -capability input.xyz > stdout.log 2> stderr.log
    assert_exit_code $? 0
}

validate_results() {
    local value=$(extract_value_from_output stdout.log)
    assert_scientific_value "2.87214" "$value" "0.0001" "RMSD"
}

main() {
    test_header "Test Name"
    cleanup_before
    run_test && validate_results
    print_test_summary
}

main "$@"
```

---

## Common Test Patterns

### Pattern 1: Energy Validation
```cpp
auto mol = TestMoleculeRegistry::createMolecule("CH4");
EnergyCalculator calc({{"method", "gfn2"}});
calc.setMolecule(mol);
double energy = calc.calculateEnergy();

const auto& data = TestMoleculeRegistry::getMolecule("CH4");
ASSERT_NEAR(energy, data.reference_energies.at("gfn2"), 1e-6);
```

### Pattern 2: Component Breakdown (GFN-FF)
```cpp
GFNFF gfnff;
gfnff.setParameter("verbosity", 3);
gfnff.Initialise(mol);
gfnff.CalculateEnergy(true, true);
// Check individual components in output
```

### Pattern 3: Parameter Validation
```cpp
// Test specific parameter combinations
json config = {{"method", "cgfnff"}, {"dispersion", true}, {"hbond", false}};
EnergyCalculator calc(config);
calc.setMolecule(mol);
double energy = calc.calculateEnergy();
```

### Pattern 4: Verification Against External Tool
```cpp
// Get XYZ path and build command
std::string xyz = TestMoleculeRegistry::getXyzPath("CH4");
std::string cmd = "cd build/test_cases && xtb " + xyz + " --gfn2";

// Execute and parse results
FILE* pipe = popen(cmd.c_str(), "r");
// ... read reference output
```

---

## FAQ

**Q: Why not just load XYZ files everywhere?**
A: Unit tests without file dependencies are portable and run in any build directory. Only external CLI tools (xtb, dftd3) need file input.

**Q: Can I add a molecule without XYZ file?**
A: No, all molecules must have both hardcoded coordinates (for unit tests) and an XYZ file (for external verification). This ensures consistency.

**Q: Which tolerance should I use?**
A: QM methods: 1e-6 Eh; Force fields: 1e-5 Eh; Dispersion: 1e-8 Eh. Adjust based on molecule size and term sensitivity.

**Q: How do I test gradients?**
A: Set calculateGradient=true, then compare analytical gradient with numerical gradient using finite differences.

**Q: Where should I put new test files?**
A: Unit tests in `testCases/`. CLI tests in `testCases/cli/<category>/`. Follow existing patterns in each location.

---

## Maintenance Guidelines

### When Adding New Energy Terms
1. Update reference data in test_molecule_registry.cpp
2. Add tolerance values
3. Update test suite to validate new term
4. Document in TESTING.md

### When External Tools Update
1. Regenerate all reference energies from XYZ files
2. Update reference_energies in test_molecule_registry.cpp
3. Run all tests and adjust tolerances if needed

### When Molecule Coordinates Change
1. Update both hardcoded coordinates (test_molecule_registry.cpp) and XYZ file
2. Regenerate reference energies
3. Run all affected tests

---

## Test Status (Current)

| Category | Tests | Passing | Notes |
|----------|-------|---------|-------|
| **Energy Methods** | test_energy_methods.cpp | ✅ | All methods |
| **GFN-FF Components** | test_gfnff_unified.cpp | ✅ | All 7 terms |
| **GFN-FF D3** | test_gfnff_d3.cpp | ✅ | Integration |
| **Molecule** | test_molecule.cpp | ✅ | Data structure |
| **CLI - RMSD** | 6 tests | ✅ 100% | 6/6 passing |
| **CLI - ConfScan** | 7 tests | ✅ 100% | 7/7 passing |
| **CLI - Curcumaopt** | 6 tests | ✅ 100% | 6/6 passing |
| **CLI - SimpleMD** | 7 tests | ❌ | JSON null crash |

---

## Architecture Notes

### Two-Phase Design (EnergyCalculator)

**Phase 1: Parameter Generation**
- Method-specific parameters (GFN-FF topology, UFF bonds, etc.)
- Generated from molecular geometry
- Output: JSON parameter set

**Phase 2: Energy/Gradient Calculation**
- Parameter set used for computation
- Multi-threading support
- All energy terms computed in parallel

### MethodFactory Resolution

Priority-based method fallback:
```cpp
"gfn2" → TBLite → Ulysses → XTB
"gfn1" → TBLite → Ulysses → XTB
"cgfnff" → Native (explicit)
"uff" → ForceField → parameter generation
```

---

*Generated: December 24, 2025*
*Use this blueprint for all testing development in Curcuma*