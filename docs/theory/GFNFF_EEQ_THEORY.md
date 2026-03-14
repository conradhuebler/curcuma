# GFN-FF Electronegativity Equalization (EEQ) Theory

## Overview
The EEQ method in GFN-FF is used to calculate geometry-dependent atomic partial charges. It is based on the principle that the electronegativities of all atoms in a molecule must equalize.

## Mathematical Formulation
The EEQ system solves a linear system of equations:
$$ \mathbf{J} \cdot \mathbf{q} = \mathbf{x} $$
where:
- $\mathbf{J}$ is the hardness/interaction matrix.
- $\mathbf{q}$ are the atomic charges.
- $\mathbf{x}$ is the electronegativity vector.

### Matrix Elements
- **Diagonal (Self-interaction):** $J_{ii} = \eta_i + \frac{\sqrt{2/\pi}}{\sqrt{\alpha_i}}$
- **Off-diagonal (Coulomb):** $J_{ij} = \frac{\text{erf}(\gamma_{ij} r_{ij})}{r_{ij}}$ with $\gamma_{ij} = \frac{1}{\sqrt{\alpha_i + \alpha_j}}$

## Two-Phase Implementation
To achieve high accuracy matching the Fortran reference, Curcuma uses a two-phase approach:

### Phase 1: Topology Charges ($q_a$)
- **Purpose**: Calculate base charges based on topology and coordination numbers.
- **Formulation**: Uses base parameters and applies the CNF (Coordination Number Factor) term in the RHS vector.
- **Reference**: Fortran `goedeckera` subroutine.

### Phase 2: Refined Energy Charges ($q$)
- **Purpose**: Fine-tune charges based on the real geometric environment.
- **Formulation**: Re-solves the EEQ system using parameters corrected by Phase 1 results ($d\chi$, $d\gamma$, $d\alpha$).
- **Key Discovery**: Phase 2 uses a modified damping and does NOT include the CNF term in the RHS (it is already captured in the $d\chi$ correction).

## Critical Improvements & Bug Fixes

### 1. The Unit Bug (Jan 2026) ✅
Covalent radii were previously passed in Bohr instead of Angstrom during topological distance calculation, leading to ~2x errors in $q_a$. Fixed by standardizing on Angstrom for `rcov` inputs.

### 2. CNF Term Synchronization ✅
Resolved inconsistency where the CNF term was applied twice or omitted. It is now correctly applied ONLY in Phase 1, following the Fortran reference `gfnff_ini.f90:563-570`.

### 3. dgam (Hardness) Corrections ✅
Charge-dependent gamma corrections were implemented. Although valid, experimental validation (`DGAM_VALIDATION_REPORT.md`) showed that they provide negligible accuracy gains (< 0.001%) for standard organic systems and remain disabled for Phase 2 by default to reduce noise.

### 4. Fragment Constraints ✅
For multi-molecule systems, the EEQ system is augmented with additional rows/columns to enforce charge neutrality for each molecular fragment independently, preventing unphysical charge transfer across vacuum.

## References
- Spicher, S.; Grimme, S. "Robust Atomistic Modeling..." *Angew. Chem. Int. Ed.* 2020.
- `src/core/energy_calculators/ff_methods/eeq_solver.cpp`
