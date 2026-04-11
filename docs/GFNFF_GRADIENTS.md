# GFN-FF Analytical Gradients

## Implementation Status
**Status**: ✅ ALL BONDED TERMS ENABLED (Feb 2026)

All analytical gradients for bonded interactions (Bonds, Angles, Torsions, Inversions) and Non-bonded (Dispersion, Repulsion, Coulomb) are implemented and integrated into the `ForceFieldThread` execution loop.

## Gradient Terms Summary

| Term | Status | Reference | Notes |
|------|--------|-----------|-------|
| **Bonds** | ✅ Active | `egbond` | Exact match achieved |
| **Angles** | ✅ Active | `egbend` | Includes distance damping derivatives |
| **Torsions** | ✅ Active | `egtors` | Properly handles NCI damping factors |
| **Extra Torsions** | ✅ Active | `egtors` | sp3-sp3 gauche corrections active |
| **Inversions** | ✅ Active | `gfnff_ini` | Planar stability for aromatics |
| **Dispersion** | ✅ Active | `gdisp0` | D4 analytical derivatives |
| **Coulomb** | ✅ Active | `engrad` | **Fixed Feb 2026**: Includes Term 1b (CN derivatives) |
| **BATM** | 📋 Planned | `batmgfnff` | Triple terms energy active, grad TODO |

## Key Technical Solutions

### Coulomb Gradient: Term 1b Fix ✅
The GFN-FF Coulomb gradient depends on coordination numbers ($\partial q / \partial CN \cdot \partial CN / \partial x$). This "Term 1b" was previously missing or using stale CN values.
- **Solution**: Force recalculation of CN derivatives in `GFNFF::Calculation()` before every gradient pass.

### NaN Stability in Torsions ✅
Resolved issues where small dihedral angles or collinear atom triplets caused NaNs in the gradient.
- **Fix**: Implemented safe normalization and epsilon-guards in `calculateDihedralGradient`.

### Unit Consistency ✅
ForceField internally calculates gradients in `Hartree/Bohr`. Curcuma optimizers expect `Hartree/Angstrom`.
- **Conversion**: $\text{grad}_{\text{Angstrom}} = \text{grad}_{\text{Bohr}} \cdot 0.529177$

## Validation Results (Status Feb 2026)
Validation against numerical finite differences (using `test_cases/test_gfnff_gradients.cpp`):
- ✅ **Standard Hydrocarbons**: Excellent consistency (error < 1e-5).
- ⚠️ **Benzene**: Sign-flip issues in some components identified, investigation ongoing.
- ⚠️ **Repulsion**: Translational invariance warnings in large systems.

## Usage
To run gradient validation for a molecule:
```bash
./release/test_cases/test_gfnff_gradients <molecule.xyz>
```
