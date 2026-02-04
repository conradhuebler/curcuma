# D3/D4 Dispersion Test Reference Value Guide

## Overview

This guide explains how to generate and fill in the reference values for the D3/D4 dispersion energy test suite (`test_dispersion.cpp`). The test suite uses `TODO_REFERENCE` placeholders (value = 0.0) that you need to fill in with calculated reference values.

## Quick Start

1. **Generate reference values** using XTB 6.6.1 or ORCA 5.0.3 (see below)
2. **Update test_dispersion.cpp** with the calculated energies
3. **Run the tests** to validate against your reference values

## How to Generate Reference Values

### Method 1: Using XTB 6.6.1 (Recommended)

XTB is freely available and provides both D3 and D4 dispersion corrections.

#### Installation

```bash
conda install -c conda-forge xtb
# or download from: https://github.com/grimme-lab/xtb
```

#### D3 Calculations

**Basic D3 with PBE0/BJ damping:**
```bash
xtb test_cases/validation/benzene.xyz --gfn 2 --d3 --func pbe0
```

Expected output includes:
```
Dispersion correction:
[stuff here]
dispersion energy: -0.012345678 Eh
```

**With three-body ATM term:**
```bash
xtb test_cases/validation/benzene.xyz --gfn 2 --d3 --func pbe0 --atm
```

**Different functional (B3LYP):**
```bash
xtb test_cases/validation/benzene.xyz --gfn 2 --d3 --func b3lyp
```

**Other damping functions:**
```bash
# Zero damping
xtb test_cases/validation/benzene.xyz --gfn 2 --d3 --func pbe0 --zero

# Bjm damping (modified Becke-Johnson)
xtb test_cases/validation/benzene.xyz --gfn 2 --d3 --func pbe0 --bjm
```

#### D4 Calculations

**Basic D4 with PBE0:**
```bash
xtb test_cases/validation/benzene.xyz --gfn 2 --d4 --func pbe0
```

D4 automatically includes three-body corrections. Expected output:
```
dispersion energy: -0.013456789 Eh
```

**Different functional:**
```bash
xtb test_cases/validation/benzene.xyz --gfn 2 --d4 --func b3lyp
```

### Method 2: Using ORCA 5.0.3

ORCA is a commercial quantum chemistry package but offers educational licenses.

#### D3 Calculations

Create `benzene_d3.inp`:
```
! B3LYP D3BJ def2-TZVP
* xyzfile 0 1 benzene.xyz
```

Run and extract:
```bash
orca benzene_d3.inp > benzene_d3.out
grep "Dispersion correction" benzene_d3.out
```

#### D4 Calculations

Create `benzene_d4.inp`:
```
! B3LYP D4 def2-TZVP
* xyzfile 0 1 benzene.xyz
```

Run:
```bash
orca benzene_d4.inp > benzene_d4.out
grep "Dispersion correction" benzene_d4.out
```

### Method 3: Script to Generate All Reference Values

```bash
#!/bin/bash
# generate_dispersion_references.sh

MOLECULES=("HH" "benzene" "ethene" "butane" "HCl" "OH")
FUNCTIONALS=("pbe0" "b3lyp")

echo "Generating D3 reference values..."
for mol in "${MOLECULES[@]}"; do
    for func in "${FUNCTIONALS[@]}"; do
        echo "D3: $mol with $func"
        xtb test_cases/validation/${mol}.xyz --gfn 2 --d3 --func ${func} > ${mol}_${func}_d3.out 2>&1
        grep "dispersion energy" ${mol}_${func}_d3.out | head -1
    done
done

echo -e "\nGenerating D4 reference values..."
for mol in "${MOLECULES[@]}"; do
    for func in "${FUNCTIONALS[@]}"; do
        echo "D4: $mol with $func"
        xtb test_cases/validation/${mol}.xyz --gfn 2 --d4 --func ${func} > ${mol}_${func}_d4.out 2>&1
        grep "dispersion energy" ${mol}_${func}_d4.out | head -1
    done
done
```

## Expected Value Ranges

To sanity-check your reference values:

### D3 Dispersion Energies

| Molecule | Size | Typical D3 Energy | Notes |
|----------|------|-------------------|-------|
| **H2** | 2 atoms | ~-5e-6 to -1e-5 Eh | Very small, numerical noise possible |
| **HCl** | 2 atoms | ~-1e-4 to -3e-4 Eh | Minimal dispersion |
| **OH** | 2 atoms | ~-1e-4 to -3e-4 Eh | Radical system |
| **Ethene (C2H4)** | 6 atoms | ~-3e-3 to -5e-3 Eh | Small organic molecule |
| **Benzene (C6H6)** | 12 atoms | ~-1e-2 to -2e-2 Eh | Aromatic π-π interactions |
| **Butane (C4H10)** | 14 atoms | ~-5e-3 to -1e-2 Eh | Alkane van der Waals |

### D4 Dispersion Energies

D4 is generally **more negative** (stronger dispersion) than D3:

| Molecule | Size | Typical D4 Energy | vs D3 |
|----------|------|-------------------|--------|
| **H2** | 2 atoms | ~-1e-5 to -2e-5 Eh | Similar magnitude |
| **HCl** | 2 atoms | ~-2e-4 to -5e-4 Eh | ~20-50% larger |
| **Ethene** | 6 atoms | ~-4e-3 to -7e-3 Eh | ~20-30% larger |
| **Benzene** | 12 atoms | ~-1.5e-2 to -2.5e-2 Eh | ~20-30% larger |
| **Butane** | 14 atoms | ~-7e-3 to -1.2e-2 Eh | ~20-30% larger |

### Three-Body Corrections (ATM)

For D3 with ATM enabled, expect:
- **Small molecules (H2, HCl)**: ATM ~0% effect
- **Medium molecules (ethene)**: ATM ~1-3% effect
- **Large aromatic (benzene)**: ATM ~5-10% effect
- **Alkanes**: ATM ~2-5% effect

## How to Update test_dispersion.cpp

### Locate the Test Case

In `test_dispersion.cpp`, find the reference you want to update. For example:

```cpp
// D3-2: Benzene - PBE0/BJ - Aromatic π-π
m_references.push_back({
    "test_cases/validation/benzene.xyz", "Benzene", 12,
    "D3", "pbe0", "bj", false,
    -1.0, -1.0, -1.0,
    TODO_REFERENCE,  // TODO: Fill from XTB
    DEFAULT_TOLERANCE_D3,
    // ...
});
```

### Calculate the Energy

```bash
xtb test_cases/validation/benzene.xyz --gfn 2 --d3 --func pbe0
# Output: dispersion energy: -0.012345678 Eh
```

### Update the Reference

Replace `TODO_REFERENCE` with the actual value:

```cpp
// OLD:
    TODO_REFERENCE,  // TODO: Fill from XTB

// NEW:
    -0.012345678,    // XTB 6.6.1 reference value
```

### Complete Example

```cpp
// D3-2: Benzene - PBE0/BJ - Aromatic π-π
m_references.push_back({
    "test_cases/validation/benzene.xyz", "Benzene", 12,
    "D3", "pbe0", "bj", false,
    -1.0, -1.0, -1.0,
    -0.012345678,    // XTB 6.6.1: xtb benzene.xyz --gfn 2 --d3 --func pbe0
    DEFAULT_TOLERANCE_D3,
    "XTB 6.6.1: xtb benzene.xyz --gfn 2 --d3 --func pbe0",
    "Aromatic system with significant π-dispersion"
});
```

## Verification

After updating reference values:

### 1. Run the Test Suite

```bash
cd release  # or your build directory
./test_dispersion --verbose
```

Expected output:
```
DFT-D3/D4 Dispersion Energy Test Suite
================================================================================
Total test cases: 36

[ OK ] D3_Benzene_pbe0_bj     Energy: -0.012345678 Eh (ref: -0.012345678) Err: 0.000000e+00
[ OK ] D3_HH_pbe0_bj           Energy: -0.000012345 Eh (ref: -0.000012345) Err: 0.000000e+00
...

================================================================================
TEST SUMMARY
================================================================================
Total tests:   36
Passed:        36
Failed:        0
Skipped:       0 (method not compiled)
Success rate:  100.0%

CSV report saved to: dispersion_test_results.csv
```

### 2. Generate CSV Report

```bash
./test_dispersion
# Creates: dispersion_test_results.csv
```

View the report:
```bash
cat dispersion_test_results.csv | head -10
```

### 3. Cross-Check Values

Compare your values with literature or other implementations:

| System | Our D3 | XTB Reference | Difference |
|--------|--------|---------------|-----------|
| Benzene | -0.01234 | -0.01235 | 0.01% ✓ |
| Butane | -0.00567 | -0.00568 | 0.18% ✓ |

## Troubleshooting

### Issue: Very Small Values for H2

H2 has minimal dispersion. If you get:
```
Error mismatch: 0.000000123 > 1e-7
```

This is likely numerical noise. Valid options:
1. **Increase tolerance**: Change `1e-7` to `1e-6` for H2 tests
2. **Accept as is**: The test infrastructure marks small placeholders as TODO

### Issue: D4 Not Compiled

If you see `[SKIP] D4_...` in output, D4 is not compiled. Either:
1. Compile with `-DUSE_D4=ON` (requires proper LAPACKE installation)
2. Focus on D3 tests only

### Issue: XTB Version Mismatch

Different XTB versions may give slightly different D3/D4 values:
- **XTB 6.4.x**: Slightly different parameters
- **XTB 6.5.x+**: Current version, recommended

Document your XTB version in the test reference source field.

## Reference Data Template

Here's a template for organizing your calculated values:

```
Molecule: benzene
File: test_cases/validation/benzene.xyz
Atoms: 12
Date: 2025-12-15
XTB Version: 6.6.1

D3 Tests:
- PBE0/BJ: -0.012345678 Eh
- B3LYP/BJ: -0.013456789 Eh
- PBE0/BJ/ATM: -0.012987654 Eh
- PBE0/Zero: -0.012111111 Eh

D4 Tests:
- PBE0: -0.013456789 Eh
- B3LYP: -0.014567890 Eh
- Custom (s8=1.5): -0.015678901 Eh
```

## Quality Assurance

### Consistency Checks

1. **D4 > D3**: D4 energies should be ~20-30% more negative than D3
2. **Three-body < 10%**: ATM correction typically <10% of total dispersion
3. **Functional dependence**: B3LYP typically gives larger dispersion than PBE0

### Tolerance Guidelines

- **Clean values** (expected energy ~-0.01): Tolerance 1e-6
- **Small values** (expected energy ~-1e-4): Tolerance 1e-7
- **Large systems** (>50 atoms): Tolerance 1e-5 acceptable

## References

- **XTB**: [GitHub - grimme-lab/xtb](https://github.com/grimme-lab/xtb)
- **D3**: Grimme et al., J. Chem. Phys. 132, 154104 (2010)
- **D4**: Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019)
- **GFN-FF**: Spicher & Grimme, Angew. Chem. Int. Ed. 59, 15665 (2020)

---

*For questions about specific reference values, consult the XTB documentation or the dispersion correction papers listed above.*
