# D4 Dispersion Fix for Heavy Elements (I, 4d/5d metals)

Date: 2026-07-16
Status: 🤖 AI-generated, ⚙️ machine-verified against tblite reference (human production testing pending)

## Symptom

Native GFN2/GFN1 D4 dispersion under-bound heavy elements (I Z=53, all 4d/5d
metals such as Ru/Rh/Pd/W/Ir/Pt), deviating from the tblite reference by
~1e-4 to 2e-3 Eh per heavy atom. All other energy terms were already
bit-exact vs tblite; the residual was entirely the D4 dispersion term.
3d metals (Z<=30) and light main-group elements were already exact.

## Root cause

File: `src/core/energy_calculators/dispersion/d4param_generator.cpp`,
`D4ParameterGenerator::initializeParameterData()` (the `m_r4_over_r2`
initializer, previously lines ~145-154).

The per-element `r4Overr2` table (the `<r^4>/<r^2>` moment ratio that sets the
C8/C6 ratio and the BJ damping radius `R0 = a1*sqrt(r4r2ij) + a2`) was only
populated for **H through Kr (Z<=36)**. The line

```cpp
if (m_r4_over_r2.size() < MAX_ELEM) m_r4_over_r2.resize(MAX_ELEM, 10.0);
```

silently filled **every element heavier than Kr with a placeholder 10.0**.

For iodine the true value is 8.4110 but 10.0 was used. A too-large r4Overr2
inflates `R0` (the damping radius), which cuts off more of the attractive
dispersion tail -> under-binding. This is exactly why heavy elements were wrong
while Z<=36 (correctly tabulated) matched tblite.

Consumption path: `D4Evaluator::computeEnergyAndGradient` (native GFN2, called
from `qm_methods/xtb_native.cpp::calcDispersionEnergy`) reads
`m_data->getSqrtZr4r2(Z)`, which is `sqrt(0.5 * m_r4_over_r2[Z-1] * sqrt(Z))`.

## Fix

Replaced the truncated 36-element `m_r4_over_r2` initializer with the complete
118-element table copied verbatim from the Fortran D4 reference
`external/gfnff/src/dftd4param.f90:134-157` (`r4Overr2(max_elem)`, identical to
the dftd4/tblite data). Only Z>36 entries change; H-Kr values are byte-identical
to before, so light and 3d-metal results are untouched. The
`resize(..., 10.0)` guard is now a dead no-op (kept as defensive padding).

Verified the new array has exactly 118 entries with I(53)=8.4110,
Rh(45)=9.5414, Ir(77)=8.3549 matching the Fortran source.

## Verification (release_tblite, USE_TBLITE=ON)

Native `-method gfn2` vs reference `-method tblite-gfn2`, `-charge 0 -no_bmt`.

### Heavy elements (the fix) — before -> after

| System | Native total (before) | Native total (after) | tblite total | diff after |
|---|---|---|---|---|
| I2 (2 I) | -7.65209556 | **-7.65233253** | -7.65233253 | 0 (8 dp) |
| PR22 (2 Ir) | -73.0496...(off 2.26e-3) | **-73.05222509** | -73.05222509 | 0 (8 dp) |

I2 dispersion component: -0.00106257 (before) -> **-0.00129953** (after);
xtb/tblite target -0.00129954.
PR22 dispersion component: -0.06841881 (before) -> **-0.07067917** (after).

Additional heavy-metal spot checks (native == tblite to 8 dp after fix):
PR01 -39.15884721, PR08 -113.05583223, ED10 -39.26983160.

### Regression controls (must stay exact) — native == tblite, unchanged

| System | native | tblite |
|---|---|---|
| H2O (light) | -5.07036982 | -5.07036982 |
| CH4 (light) | -4.17507458 | -4.17507458 |
| C6H6 (light) | -15.87853118 | -15.87853118 |
| ED02 (Fe, 3d) | -27.85077571 | -27.85077571 |
| ED03 (Ni, 3d) | -23.34102052 | -23.34102052 |

All bit-identical; no regression.

## Notes / caveats

- The same `m_r4_over_r2` table also feeds native GFN-FF dispersion
  (`GenerateDispersionPairsNative`). Since only Z>36 entries changed, the
  light/neutral GFN-FF validation set (all Z<=36) is unaffected; GFN-FF results
  for heavy elements are now corrected too (previously also using the 10.0
  placeholder), but were not separately re-validated here.
- Transition-metal D4 reference data (alphaiw, refcn) beyond r4Overr2 was
  already correct (iodine alphaiw matched the Fortran source exactly); the sole
  discrepancy was the truncated r4Overr2 table.
- The `ctest -L d4_diag` entries report "Not Run" in release_tblite because the
  `diag_curcuma_d4_potential` helper is not built in that directory — a
  pre-existing harness gap, not a regression.
