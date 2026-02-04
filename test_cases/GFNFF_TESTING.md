# GFN-FF Testing Guide

Session-ready documentation for GFN-FF testing framework.

---

## Quick Start

### Run Unit Tests (Hardcoded Molecules)
```bash
cd /path/to/curcuma/build/test_cases
ctest -R "gfnff" --output-on-failure
```

### Run External Verification (XYZ Files)
```bash
# Get XYZ path programmatically
test_cases/molecules/larger/CH4.xyz

# Run external program
xtb test_cases/molecules/larger/CH4.xyz --gfnff
```

---

## Two-Tier Molecule Strategy

| Use Case | Source | Path Required? |
|----------|--------|----------------|
| **Unit Tests (C++)** | TestMoleculeRegistry (hardcoded) | ❌ No - self-contained |
| **External Programs** | XYZ files in molecules/ | ✅ Yes - for xtb/dftd3 |

**Why Hardcoded?** Unit tests run anywhere without file system dependencies. Tests are path-independent and can execute in any build directory.

---

## Available Molecules

### Dimers (2 atoms)

| Name | Registry Key | XYZ File | Atoms | Tests |
|------|--------------|----------|-------|-------|
| Hydrogen | `"H2"` | dimers/HH.xyz | H-H | Bond energy |
| Hydrogen Chloride | `"HCl"` | dimers/HCl.xyz | Cl-H | Heteronuclear bond |
| Hydroxyl Radical | `"OH"` | dimers/OH.xyz | O-H | Radical, charge |

### Trimers (3 atoms)

| Name | Registry Key | XYZ File | Atoms | Tests |
|------|--------------|----------|-------|-------|
| Water Trimer | `"H2O"` | trimers/water.xyz | O-H-H | Angle energy |
| Hydrogen Cyanide | `"HCN"` | trimers/HCN.xyz | H-C-N | Linear molecule |
| Ozone | `"O3"` | trimers/O3.xyz | O-O-O | Ozone test |

### Larger (5+ atoms)

| Name | Registry Key | XYZ File | Atoms | Tests |
|------|--------------|----------|-------|-------|
| Methane | `"CH4"` | larger/CH4.xyz | C + 4H | Multiple equivalent bonds |
| Methanol | `"CH3OH"` | larger/CH3OH.xyz | C₁O₁H₆ | All energy terms + torsion |
| Dimethyl Ether | `"CH3OCH3"` | larger/CH3OCH3.xyz | C₂O₁H₆ | Complex torsions |
| Benzene | `"C6H6"` / `"benzene"` | larger/C6H6.xyz | C₆H₆ | Aromatic ring |
| Monosaccharide | `"monosaccharide"` | larger/monosaccharide.xyz | 27 atoms | Medium-large system |
| Triose | `"triose"` | larger/triose.xyz | 66 atoms | Large system |

---

## Using TestMoleculeRegistry in Tests

### C++ Unit Tests (No Paths Required)

```cpp
#include "test_cases/core/test_molecule_registry.h"
using namespace TestMolecules;

// Create molecule - path-independent!
curcuma::Molecule mol = TestMoleculeRegistry::createMolecule("CH4");

// Run test
GFNFF gfnff;
gfnff.setParameter("charges", ...);
double energy = gfnff.calculateEnergy(mol);
```

### Get XYZ Path for External Verification

```cpp
// Get XYZ path programmatically (relative to build/test_cases/)
std::string xyz_path = TestMoleculeRegistry::getXyzPath("CH4");
// Returns: "molecules/larger/CH4.xyz"

// Run external program
std::string cmd = "xtb " + xyz_path + " --gfnff";
system(cmd.c_str());
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

## Adding New Molecules to Registry

### Step 1: Create XYZ File

```bash
# Create file in appropriate category
vim test_cases/molecules/<category>/<name>.xyz
```

XYZ format example:
```
<number of atoms>
<comment line>
<element> <x> <y> <z>
...
```

### Step 2: Generate Reference Energies

```bash
cd /tmp  # or temporary directory
xtb molecule.xyz --gfnff
```

Extract reference energies from XTB output:
```
bond energy:       -0.123456 Eh
angle energy:        0.000123 Eh
torsion energy:       0.000000 Eh
repulsion energy:    0.087654 Eh
electrostat energy: -0.001234 Eh
dispersion energy:   -0.002345 Eh
```

### Step 3: Add to test_molecule_registry.cppEdit template:
```cpp
{
    "NewMol", {
        .name = "NewMol",
        .description = "Description from XYZ file",
        .category = "dimers|trimers|larger",
        .atoms = {
            {atomic_number, Eigen::Vector3d(x, y, z)},
            // ... all atoms
        },
        .reference_energies = {
            {"gfnff_xtb", total_energy},      // if available
            {"d3", d3_energy},                // if available
        },
        .tolerances = {
            {"gfnff_xtb", 1e-6},              // Hartree tolerance
            {"d3", 1e-8},                    // D3 tolerance
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

### Step 5: Update Documentation

Add molecule entry to the appropriate table in this file.

### Step 6: Run Tests

```bash
make test_gfnff_unified
ctest -R "gfnff" --verbose
```

---

## Debugging Failed Tests

### 1. Check Molecule Coordinates

```cpp
Molecule mol = TestMoleculeRegistry::createMolecule("CH4");
mol.print();  // Verify coordinates match reference
```

### 2. Verify Against XYZ File (External Reference)

```bash
# Get XYZ path
std::string xyz_path = TestMoleculeRegistry::getXyzPath("CH4");

# Run XTB for reference
cd build/test_cases
xtb $xyz_path --gfnff

# Compare component energies with test output
```

### 3. Extract Energy Components in Test

```cpp
GFNFF gfnff;
gfnff.setParameter("verbosity", 3);  // Enable component output
double energy = gfnff.calculateEnergy(mol);
// Look for: BondEnergy, AngleEnergy, Repulsion, etc.
```

### 4. Compare Reference Energies

```cpp
const MoleculeData& data = TestMoleculeRegistry::getMolecule("CH4");
if (data.hasReferenceEnergy("d3")) {
    double ref_energy = data.reference_energies.at("d3");
    double tolerance = data.tolerances.at("d3");
    // Compare with calculated value
}
```

---

## Directory Structure

```
test_cases/
├── GFNFF_TESTING.md                      ← This file
│
├── core/
│   └── test_molecule_registry.h/.cpp     ← SHARED hardcoded library
│       - H2, HCl, OH (dimers)
│       - HCN, H2O, O3 (trimers)
│       - CH4, CH3OH, CH3OCH3, C6H6 (larger)
│       - monosaccharide (27 atoms)
│       - triose (66 atoms)
│
├── test_gfnff_unified.cpp                ← Uses TestMoleculeRegistry
├── test_gfnff_d3.cpp                     ← Uses TestMoleculeRegistry
│
└── molecules/                            ← XYZ for external programs
    ├── dimers/
    │   ├── HH.xyz
    │   ├── HCl.xyz
    │   └── OH.xyz
    ├── trimers/
    │   ├── HCN.xyz
    │   ├── water.xyz
    │   └── O3.xyz
    └── larger/
        ├── CH4.xyz
        ├── CH3OH.xyz
        ├── CH3OCH3.xyz
        ├── C6H6.xyz
        ├── monosaccharide.xyz
        └── triose.xyz
```

---

## Test Suite Reference

### test_gfnff_unified.cpp

**Purpose**: Main regression test suite for energy components.

**Tests**:
- Energy component breakdown (bond, angle, torsion, repulsion, electrostatic, dispersion)
- Parameter accuracy validation (vbond, vangl, vtor)
- Comprehensive validation (CH3OH - all terms)
- Parameter flag combinations

**Tolerance**:
- Total energy: 1e-6 Eh
- Components: 1e-6 Eh
- Parameters: 1e-8

### test_gfnff_d3.cpp

**Purpose**: D3 dispersion integration validation.

**Tests**:
- GFN-FF D3 vs standalone D3 consistency
- D3 energy accuracy for all molecules
- Dispersion-only validation

**Tolerance**:
- D3 energies: 1e-4 Eh (stricter for some molecules)

---

## Reference Data Sources

| Source | Method | Generation Date |
|--------|--------|-----------------|
| XTB 6.6.1 (8d0f1dd) | GFN-FF energy components | Dec 2025 |
| DFT-D3 standalone | Dispersion-only energies | Dec 2025 |

For updating reference energies:
```bash
xtb molecule.xyz --gfnff  # For full GFN-FF
./dftd3 molecule.xyz     # For D3-only
```

---

## Common Patterns

### Creating Unit Test Molecule

```cpp
// Standard pattern - path-independent
using namespace TestMolecules;
Molecule mol = TestMoleculeRegistry::createMolecule("CH4");
```

### Creating External Verification Command

```cpp
// Get XYZ path and build command
std::string xyz_path = TestMoleculeRegistry::getXyzPath("CH4");
std::string build_dir = std::filesystem::current_path().string();
std::string cmd = "cd " + build_dir + "/test_cases && xtb " + xyz_path + " --gfnff";

// Execute and capture output
FILE* pipe = popen(cmd.c_str(), "r");
// ... read results
```

### Using Reference Energies

```cpp
const MoleculeData& data = TestMoleculeRegistry::getMolecule("CH4");

// Check if reference exists
if (data.hasReferenceEnergy("d3")) {
    double ref = data.reference_energies.at("d3");
    double tol = data.tolerances.at("d3");
    // Validate calculated value
    ASSERT_NEAR(calc_energy, ref, tol);
}
```

---

## FAQ

**Q: Why not just load XYZ files everywhere?**
A: Unit tests without file dependencies are more portable and run in any build directory. XYZ files are only needed for external programs (xtb, dftd3) that require file input.

**Q: Can I add a molecule without XYZ file?**
A: No, all molecules must have both hardcoded coordinates (for unit tests) and an XYZ file (for external verification). This ensures consistency between test and reference.

**Q: What if my molecule name differs from XYZ filename?**
A: Add both to the `s_xyz_paths` map (see test_molecule_registry.cpp). Example: `{"C6H6", "molecules/larger/C6H6.xyz"}, {"benzene", "molecules/larger/C6H6.xyz"}`

**Q: How do I update reference energies?**
A: Run XTB/D3 externally using the XYZ file, then update the `reference_energies` map in test_molecule_registry.cpp.

**Q: What's the tolerance for energy accuracy?**
A: Standard: 1e-6 Eh for total energy, 1e-8 for D3 dispersion. Adjust based on molecule size and term sensitivity.

---

## Maintenance Notes

**When adding new energy terms**:
1. Update reference data in test_molecule_registry.cpp
2. Add tolerance values for new term
3. Update test suite to validate new term

**When XTB updates**: regenerate all reference energies from XYZ files using new XTB version.

**When molecule coordinates change**: update both hardcoded (test_molecule_registry.cpp) and XYZ file simultaneously.

---

*Generated: December 24, 2025*
*Use this guide to quickly start GFN-FF development in any new session*