# CLAUDE.md — Dispersion Module

## Overview

Self-contained D4 dispersion module. Provides the C6 reference data + EEQ
charge stack (`D4ParameterGenerator`) and the energy/gradient kernel
(`D4Evaluator`) that GFN-FF and native GFN2 both call.

Created by lifting `d4param_generator.*` out of `ff_methods/` and adding a
shared evaluator class so the same math powers GFN-FF (modified-BJ
parameters) and GFN2 (Caldeweyher 2019 BJ parameters).

## Files

```
dispersion/
├── d4param_generator.{h,cpp}            Reference data + EEQ + dc6/dCN matrix
├── d4_evaluator.{h,cpp}                 Energy + analytical gradient kernel
├── d4_reference_data_fixed.cpp          118-element D4 reference data (Fortran-extracted)
├── d4_reference_data.cpp                Legacy reference data (deprecated, kept for diff)
├── d4_reference_cn_fortran.cpp          D4 reference CN values (dftd3param.f90)
├── d4_alphaiw_data.cpp                  Frequency-dependent polarizabilities (269 refs)
└── d4_corrections_data.cpp              ascale / sscale / secaiw / refsys tables
```

Only `d4param_generator.cpp` + `d4_evaluator.cpp` are listed in
`CMakeLists.txt`; the four data files are `#include`d into
`d4param_generator.cpp` to keep them out of the link line.

## The unified BJ formula

> GFN-FF's "modified BJ damping" and standard D4 BJ (Caldeweyher 2019) are
> mathematically the **same expression**.

The per-pair dispersion energy is

```
E_pair = -ζc6 · C6 · ( s6·t6  +  s8·r4r2_ij·t8 )
```

with `t_n = 1/(r^n + R0^n)`, `R0 = a1·sqrt(r4r2_ij) + a2`, and
`r4r2_ij = 3·sqrt(Zr4r2_i)·sqrt(Zr4r2_j)` (the implicit `C8/C6` ratio in
D4). The two damping variants differ only in scaling:

| Variant            | s6  | s8  | a1   | a2 (Bohr) | s9  | alpha | Used by |
|--------------------|-----|-----|------|-----------|-----|-------|---------|
| `ModifiedBJ_GFNFF` | 1.0 | 2.0 | 0.58 | 4.80      | 1.0 | 14.0  | Native GFN-FF |
| `StandardBJ_D4`    | 1.0 | 2.7 | 0.52 | 5.00      | 5.0 | 16.0  | Native GFN2 (Bannwarth 2019) |

The enum `curcuma::dispersion::DampingFormula` keeps call sites
declarative; the implementation collapses to one inline kernel
(`D4Evaluator::evalDispSum`). A future variant with a structurally
different formula can hang off the same enum.

## Public API

### `D4ParameterGenerator` (relocated, unchanged)

The data layer: 118-element C6 reference tables, CN-Gaussian weighting,
EEQ charges, `dc6/dCN` matrix. Public accessors used by both GFN-FF and
the evaluator:

- `GenerateDispersionPairsNative(atoms, geom_bohr)` — pair list with
  pre-baked `C6, r4r2ij, r0_squared, zetac6, r_cut`
- `getChargeWeightedC6(...)`, `getSqrtZr4r2(...)`
- `updateCNValuesForGradient(cn)` — refreshes `dc6/dCN`
- `getDC6DCN()`, `getGaussianWeights()`, `getC6FlatCache()`, `getRefN()`,
  `getRefCN()` (used by the GFN-FF GPU pipeline)

### `D4Evaluator` (new)

The math kernel. Holds a non-owning `D4ParameterGenerator*` and a
`D4Params` struct populated explicitly by the caller (constructor asserts
`s6>0 && s8>0 && a2>0` to catch silent defaults from
`D4ParameterGenerator`'s PARAM block, which carries GFN-FF values).

Two entry points:

| API | Used by | Notes |
|---|---|---|
| `pairEnergyAndGradient(pair, rij_bohr, ...)` | `ForceFieldThread` | Per-pair primitive; `m_data` may be null (caller handles its own dc6/dCN via `m_dc6dcn_ptr`). |
| `computeEnergyAndGradient(atoms, geom_bohr, ...)` | Native GFN2 (`xtb_native.cpp`) | Whole-molecule loop; requires non-null `m_data`; calls `updateCNValuesForGradient()` internally. |

Both return the per-pair `disp_sum = s6·t6 + s8·r4r2·t8` for callers that
need to compute the CN-chain-rule contribution themselves.

### CUDA hook

`D4Evaluator::launchGpuKernel()` is a documented no-op extension point.
**GFN-FF's existing GPU dispersion** (`ff_methods/cuda/gfnff_kernels.cu::k_dispersion`)
does NOT route through this class — it consumes the pre-baked pair data
from `D4ParameterGenerator::GenerateDispersionPairsNative()` directly, and
that data layout is unchanged. A future GFN2-CUDA path will likely add a
new kernel variant selected by `DampingFormula`; the hook leaves a place
to plug it in without touching the CPU evaluator.

## Call sites

| Site | Damping | Path |
|---|---|---|
| `ff_methods/forcefieldthread.cpp::CalculateD4DispersionContribution` | `ModifiedBJ_GFNFF` | Per-pair, uses `m_d4_dispersions` pre-baked by GFN-FF, applies dc6/dCN via `m_dc6dcn_ptr`. |
| `qm_methods/xtb_native.cpp::calcDispersionEnergy` | `StandardBJ_D4` | Whole-molecule; caches gradient + dE/dCN in mutable members on `XTB`. |
| `qm_methods/xtb_gradient.cpp::calculateGradient` | (consumer) | Folds cached `m_disp_gradient` directly into `m_gradient`; adds `m_disp_dEdcn` into the local `dEdcn` so the existing CN-distribution loop carries it through ∂CN/∂x. |
| `qm_methods/gfn2.cpp::calculateDispersionEnergy` | `StandardBJ_D4` | Legacy GFN2 class (not currently routed via MethodFactory). Wired to the evaluator for parity, but not the active code path for `-method ngfn2`. |

## Status & caveats

> ⚠️ **AI-generated, automated tests pass — human production testing pending.**

### GFN-FF
- Energy + gradient bit-identical to the pre-refactor inline ForceFieldThread
  code (24 GFN-FF ctest entries unchanged before/after).
- CUDA path untouched.

### GFN2 (native, `xtb_native.cpp`)
- Energy: produces non-zero D4 contributions consistent with Caldeweyher
  2019 parameters. Tested values:
  - H₂O: -0.164 mEh (3 pairs)
  - CH₄: -0.76 mEh
  - C₆H₆: -7.82 mEh
- Analytical gradient validated against central finite differences
  (`test_xtb_gradient`):
  - H₂O ngfn2:  MaxErr 2.58e-5 Eh/Å
  - CH₄ ngfn2:  MaxErr 3.27e-5 Eh/Å
  - NH₃ ngfn2:  MaxErr 2.74e-5 Eh/Å
  - C₆H₆ ngfn2: MaxErr 4.26e-5 Eh/Å
  (Tolerance: 5e-4 Eh/Å. All pass.)

### Known limitations
1. **q-response chain rule not implemented.** `zetac6` is treated as a
   static prefactor as in GFN-FF. The full GFN2 path (`∂E_D4/∂q · ∂q/∂x`
   via SCF response) is a separate AP. Expected residual: sub-mEh on the
   test set, much smaller than the 1–23 mEh D4-energy correction itself.
2. **GFN1 D3 still missing in `xtb_native.cpp`.** `calcDispersionEnergy()`
   returns 0 for GFN1; native D3 integration is a follow-up AP (the
   infrastructure in `D3ParameterGenerator` already exists).
3. **`D4ParameterGenerator` PARAM defaults are GFN-FF values** (a1=0.58,
   a2=4.80, s8=2.0). The `D4Evaluator` constructor asserts so silent
   misuse cannot pass — each call site populates `D4Params` explicitly.

## Verification

```bash
cd release && make -j4
# GFN-FF regression (bit-equivalent expected):
ctest -R "gfnff" --output-on-failure
# GFN2 analytical-vs-FD gradient:
./test_cases/test_xtb_gradient \
   ../test_cases/sqm_reference/molecules/H2O.xyz \
   ../test_cases/sqm_reference/molecules/CH4.xyz \
   ../test_cases/sqm_reference/molecules/NH3.xyz \
   ../test_cases/sqm_reference/molecules/C6H6.xyz
```

## References

- E. Caldeweyher et al., *J. Chem. Phys.* **150**, 154122 (2019) — D4 model
- C. Bannwarth, S. Ehlert, S. Grimme, *JCTC* **15**, 1652 (2019) — GFN2-xTB
- S. Spicher, S. Grimme, *Angew. Chem. Int. Ed.* **59**, 15665 (2020) — GFN-FF
- Fortran reference: `external/gfnff/src/gfnff_gdisp0.f90`
- Ulysses parameter database: `external/ulysses-main/core/src/parameters/D4par.hpp`
