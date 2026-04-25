# Native xTB Implementation Status (neue modulare Architektur)

**Last Updated**: 2026-04-25
**Implementation**: AI-generated, machine-tested — **human production testing pending**
**Location**: `src/core/energy_calculators/qm_methods/`

---

## Overview

Curcuma has **two native xTB implementations**:

1. **Old monolithic** (`gfn1.cpp`, `gfn2.cpp`) — Status: 0/7 (GFN2) and 2/7 (GFN1) tests vs. TBLite. **Deprecated.**
2. **New modular** (`xtb_native.cpp` + kernel modules) — TBLite-architecture port. **Active development.**

This document tracks the **new modular implementation** only. For the migration plan, see [`NATIVE_XTB_ROADMAP.md`](NATIVE_XTB_ROADMAP.md).

---

## Architecture

```
qm_methods/
├── xtb_native.h/.cpp      — SCF driver, public API, basis builder
├── xtb_h0.cpp            — Overlap, H0, CN, repulsion
├── xtb_coulomb.cpp       — Isotropic shell-resolved Coulomb (gamma matrix)
├── xtb_thirdorder.cpp    — Third-order onsite corrections
├── xtb_multipole.cpp     — GFN2 multipole (dipole + quadrupole)
├── xtb_scf.cpp           — Fock builder, eigensolver, Mulliken populations
└── parameters/
    ├── gfn1_params.hpp   — Auto-extracted from TBLite
    ├── gfn2_params.hpp   — Auto-extracted from TBLite
    └── xtb_params_extra.hpp — Pauling EN, radii, CN functions, kshell/kpair
```

**Key design**: Mirrors TBLite's module split (`xtb/h0.f90`, `coulomb/multipole.f90`, `scf/potential.f90`, etc.). Parameters are auto-extracted by `scripts/extract_xtb_params.py` directly from TBLite Fortran source — no manual parameter tables.

---

## Implementation Status

### Implemented

| Component | Status | Notes |
|-----------|--------|-------|
| **Basis (STO→CGTO)** | ✅ | Gram-Schmidt orthogonalization, s/p/d shells |
| **Overlap matrix** | ✅ | Analytical CGTO overlap vs. TBLite dump |
| **H0 (self-energy + CN shift)** | ✅ | `getSelfEnergies()`, `getHamiltonianH0()` |
| **H0 off-diagonal scaling** | ✅ | `hscale` computed on-the-fly in `getHamiltonianH0()` |
| **Coordination numbers** | ✅ | `cn_exp()` (GFN1), `cn_gfn()` (GFN2) |
| **Isotropic Coulomb (ES2)** | ✅ | Shell-resolved Klopman-Ohno gamma matrix |
| **Third-order (ES3)** | ✅ | Atom-resolved (GFN1), shell-resolved (GFN2) |
| **Multipole (AES2, GFN2)** | ✅ | Dipole + quadrupole integrals, damping, interaction matrices |
| **SCF loop** | ✅ | Linear damping, closed-shell, integer occupation |
| **Mulliken populations** | ✅ | Shell + atomic charges; GFN2: atomic dipoles/quadrupoles |
| **Repulsion energy** | ✅ | Exponential pairwise model |
| **Energy decomposition** | ✅ | Electronic + coulomb + third-order + multipole + repulsion |

### Missing / Stubs

| Component | Status | Notes |
|-----------|--------|-------|
| **Analytical gradients** | ❌ Not implemented | Required for optimization/MD. See Roadmap AP 4. |
| **D4 dispersion** | ❌ Stub | `m_E_dispersion = 0.0`. DFT-D4 interface exists but not wired. |
| **DIIS / advanced mixer** | ❌ Not implemented | Only linear damping (0.4). `diis_accelerator.h` exists but unused. |
| **Thermal smearing** | ❌ Declared but unused | `m_electronic_temp = 300 K` never used. Integer occupation only. |
| **Open-shell / unrestricted** | ❌ Not implemented | Test driver rejects odd-electron systems. |
| **d-shell multipole integrals** | ⚠️ Partial | `ao_to_type()` returns -1 for `ang >= 2`. d-functions skipped in multipole. |

---

## Test Status

### Kernel-level tests (vs. TBLite dumps)

| Test | Validates | Status |
|------|-----------|--------|
| `test_xtb_overlap` (3.1) | CGTO overlap vs. TBLite | ✅ max |ΔS| < 1e-8 |
| `test_xtb_h0` (3.2) | Bare Hamiltonian H0 vs. TBLite | ✅ max |ΔH| < 1e-7 |
| `test_xtb_coulomb` (3.3) | Gamma matrix structural checks | ✅ Symmetry, diagonal, limits |
| `test_xtb_scf_snapshot` (3.6) | Full SCF from dumped S/H0 | 🔄 In progress (Δε_max < 1e-4 Eh gate) |

**Test molecules**: H2, He2, LiH, H2O, CH4, NH3, C6H6

### End-to-end integration tests

| Method | Old implementation | New implementation |
|--------|-------------------|-------------------|
| Native GFN2 (`ngfn2`) | 0/7 vs. TBLite | **Wired via GFN2Method wrapper** (AP2 done). Energies returned, but ~35–60 mEh off TBLite. Crash in `addMultipolePotential` fixed (uninit `dp_at`). |
| Native GFN1 (`ngfn1`) | 2/7 vs. TBLite | **Wired via GFN1Method wrapper** (AP2 done). Energies returned, ~5–70 mEh off TBLite. |

**Note**: The new `XTB` class is now accessible via `MethodFactory` through `GFN2Method`/`GFN1Method` wrappers (AP2 complete). Default `gfn2`/`gfn1` still routes through TBLite priority chain (AP3 pending).

---

## Differences: Old vs. New Implementation

| Aspect | Old (`gfn2.cpp`) | New (`xtb_native.cpp`) |
|--------|-----------------|----------------------|
| **Architecture** | Monolithic (~1236 lines) | Modular (6 files, TBLite-mirroring) |
| **Parameters** | Manual tables, incomplete | Auto-extracted from TBLite Fortran |
| **hscale** | Pre-built matrix (buggy) | On-the-fly in H0 builder (correct) |
| **Multipole** | Simplified dipole-dipole only | Full dipole + quadrupole with damping |
| **SCF mixer** | DIIS (8 vectors) | Linear damping only |
| **Gradients** | Implemented (in old GFN2) | **Not yet implemented** |
| **MethodFactory** | Wired as `ngfn2` / `ngfn1` | **Not wired** |

---

## Known Issues

1. ~~**`hscale` comment stale**~~ — **FIXED in AP1** (b0dbfc2)
2. ~~**`m_h0.rad` not filled**~~ — **FIXED in AP1** (b0dbfc2)
3. ~~**`UpdateMolecule()` not overridden**~~ — **FIXED in AP1** (b0dbfc2): `UpdateMolecule()` overridden with cache invalidation.
4. **`Calculation(bool gradient)` ignores `gradient`**: The `gradient` parameter is accepted but gradients are not computed. Returns energy only. **Gradient stub added in AP2** — returns zero matrix with warning. (Will break optimization/MD until AP 4.)
5. **`ngfn2` crash in `addMultipolePotential()`** — **FIXED in AP2**: `m_wfn.dp_at`/`qp_at` were uninitialized in `buildReferenceOccupations()`. Crash occurred on first SCF iteration for any GFN2 molecule.

---

## References

- **Roadmap**: [`NATIVE_XTB_ROADMAP.md`](NATIVE_XTB_ROADMAP.md)
- **Multipole details**: [`GFN2_MULTIPOLE_IMPLEMENTATION.md`](GFN2_MULTIPOLE_IMPLEMENTATION.md)
- **Old implementation status**: `docs/archive/NATIVE_QM_IMPLEMENTATION_STATUS.md`
- **Old implementation guide**: `docs/archive/NATIVE_QM_METHODS_IMPLEMENTATION.md`
- **TBLite source**: `external/tblite/src/tblite/xtb/`

---

*Status updated: 2026-04-25*
*Next update: After MethodFactory wiring + gradient implementation*
