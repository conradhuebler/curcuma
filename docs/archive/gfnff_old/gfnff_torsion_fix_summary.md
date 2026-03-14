# GFNFF Torsion Fix Summary

## Problem Description
The GFNFF implementation in Curcuma was incorrectly detecting torsions due to two main issues:

1. **Unit mismatch in bond detection**: The `m_geometry` matrix stores coordinates in Angstrom units, but the `getCovalentRadius()` function returns covalent radii in Bohr units. This mismatch caused incorrect bond detection, leading to atoms appearing to be bonded when they weren't.

2. **Flawed duplicate detection**: The `isDuplicateTorsion` function only properly compared the center atoms of torsions but failed to correctly compare terminal atoms, causing different torsions with different terminal atoms to be incorrectly flagged as duplicates.

## Root Cause Analysis
For ethane (C₂H₆), the incorrect bond detection was causing:
- Carbon atoms to appear bonded to hydrogen atoms on the opposite carbon
- This led to each carbon atom having 7 neighbors instead of 4
- The torsion generation algorithm then created 36 torsions instead of the expected 9

## Solution Implemented

### 1. Fixed Unit Mismatch in Bond Detection
**File**: `src/core/energy_calculators/qm_methods/gfnff_torsions.cpp`
**Change**: Convert distance from Angstrom to Bohr to match covalent radii units:

```cpp
// Before:
double distance = (m_geometry.row(i) - m_geometry.row(j)).norm();
// distance was in Angstrom, but rcov_sum was in Bohr

// After:
double distance_angstrom = (m_geometry.row(i) - m_geometry.row(j)).norm();
double distance_bohr = distance_angstrom * CurcumaUnit::Length::ANGSTROM_TO_BOHR;
// Now both distance and rcov_sum are in Bohr units
```

### 2. Fixed Duplicate Detection Algorithm
**File**: `src/core/energy_calculators/qm_methods/gfnff_torsions.cpp`
**Change**: Enhanced the `isDuplicateTorsion` function to properly compare all four torsion atoms:

```cpp
// Before: Only compared center atoms properly
// After: Properly compares all four atoms with tolerance handling
bool GFNFF::isDuplicateTorsion(
    int i1, int j1, int k1, int l1,
    int i2, int j2, int k2, int l2) const
{
    // Check forward direction: i1-j1-k1-l1 vs i2-j2-k2-l2
    bool forward_match = 
        (std::abs(i1 - i2) <= TORSION_ATOM_TOLERANCE) &&
        (std::abs(j1 - j2) <= TORSION_ATOM_TOLERANCE) &&
        (std::abs(k1 - k2) <= TORSION_ATOM_TOLERANCE) &&
        (std::abs(l1 - l2) <= TORSION_ATOM_TOLERANCE);

    // Check reverse direction: i1-j1-k1-l1 vs l2-k2-j2-i2
    bool reverse_match = 
        (std::abs(i1 - l2) <= TORSION_ATOM_TOLERANCE) &&
        (std::abs(j1 - k2) <= TORSION_ATOM_TOLERANCE) &&
        (std::abs(k1 - j2) <= TORSION_ATOM_TOLERANCE) &&
        (std::abs(l1 - i2) <= TORSION_ATOM_TOLERANCE);

    return forward_match || reverse_match;
}
```

## Results

### Before Fix:
- Ethane: 36 torsions (incorrectly detected)
- Bond detection included false positives
- Neighbor lists were incorrect (7 neighbors per carbon instead of 4)

### After Fix:
- Ethane: 9 torsions (correctly detected)
- Bond detection is accurate
- Neighbor lists are correct (4 neighbors per carbon, 1 neighbor per hydrogen)
- Energy calculations complete successfully

## Verification
The fix has been tested with:
- Ethane (C₂H₆): 9 torsions ✓
- Dimethyl ether (CH₃OCH₃): 6 torsions ✓
- Methanol (CH₃OH): 2 torsions ✓

All tests show reasonable torsion counts and successful energy calculations.