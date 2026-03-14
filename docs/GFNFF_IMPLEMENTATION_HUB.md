# GFN-FF Implementation Documentation Hub

**Last Updated**: 2026-02-07
**Status**: ⚠️ **UNIFIED VALIDATION ACTIVE (ACCURACY REFINEMENT ONGOING)**

This hub provides a central entry point for all documentation related to the native GFN-FF implementation in Curcuma.

---

## 🗺️ Documentation Map

| Category | Document | Description |
| :--- | :--- | :--- |
| **Status** | **[GFNFF_STATUS.md](GFNFF_STATUS.md)** | **Current project status and accuracy metrics (start here!)** |
| **Theory** | [GFNFF_EEQ_THEORY.md](theory/GFNFF_EEQ_THEORY.md) | Electronegativity Equalization (EEQ) and charge model |
| | [GFNFF_TORSION_THEORY.md](theory/GFNFF_TORSION_THEORY.md) | Torsion and extra-torsion potential formulation |
| | [GFNFF_INVERSION_THEORY.md](theory/GFNFF_INVERSION_THEORY.md) | Out-of-plane and planarity stabilization |
| **Mechanics** | **[GFNFF_GRADIENTS.md](GFNFF_GRADIENTS.md)** | Analytical gradients and derivative consistency |
| | [GFNFF_PHYSICS_AUDIT.md](GFNFF_PHYSICS_AUDIT.md) | In-depth audit of physical parameters and damping |
| | [GFNFF_DISPERSION_FIX.md](GFNFF_DISPERSION_FIX.md) | Details on D4 integration and BJ-damping formulas |
| **Infras.** | [PARAMETER_FLOW_ARCHITECTURE.md](PARAMETER_FLOW_ARCHITECTURE.md) | How parameters flow from JSON to FF kernels |

---

## 🛠️ Implementation Architecture

Curcuma uses a **Two-Phase Design Pattern** for GFN-FF to ensure scientific accuracy while maintaining high performance:

1.  **Phase 1: Parameter Generation** (`GFNFF` class)
    *   Topology detection (Coordination Numbers, Hybridization, Pi-systems).
    *   Dynamic parameter assignment based on the chemical environment.
    *   Output: Standardized JSON parameter arrays.
2.  **Phase 2: Energy Calculation** (`ForceFieldThread`)
    *   High-performance, multi-threaded execution kernels.
    *   Analytical gradients for all 7 energy terms.

---

## 🧪 Validation & Testing

We maintain a strict validation suite to reach the goal of **< 1 µEh accuracy** compared to the Fortran reference.

*   **Unified Runner**: `release/test_cases/test_gfnff_validation`
*   **Golden References**: Standardized JSON data in `test_cases/reference_data/`
*   **Command**: `ctest -R "gfnff_val_"`

---

## 📜 Historical Archive

For historical context, design decisions, and session summaries, refer to the archive:
*   [Development History](archive/history/)
*   [Obsolete Reports & Bug Analyses](archive/obsolete_reports/)

---
*Educational Note: The GFN-FF implementation prioritizes algorithmic transparency. The source code in `src/core/energy_calculators/ff_methods/` is extensively commented with references to the primary literature.*
