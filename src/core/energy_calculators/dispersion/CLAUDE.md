# CLAUDE.md - dispersion/

D4 dispersion: reference data, EEQ charge model, C6 weighting (`D4ParameterGenerator`) and the
energy/gradient kernel (`curcuma::dispersion::D4Evaluator`). Native D3 lives in
`ff_methods/d3param_generator.*`, not here.

> ⚠️ **AI-generated, automated tests pass - human production testing pending.**

## Files

- `d4param_generator.{h,cpp}` - reference data, CN-Gaussian weights, charge scaling, `dc6/dCN`, GFN-FF pair list
- `d4_evaluator.{h,cpp}` - `D4Params`, `DampingFormula`, two-body kernel + `computeATM()` (three-body)
- `d4_charge_model.{h,cpp}` - single-shot dftd4 EEQ with analytic `dq/dx` (`D4ChargeModel`)
- `d4_ncoord.{h,cpp}` - dftd4 EN-weighted covalent CN (GFN2 path); `d4_charge_scaling.h` - shared dftd4 zeta
- Data, `#include`d into `d4param_generator.cpp` (not compiled on their own): `d4_reference_data_fixed.cpp`, `d4_reference_cn_fortran.cpp`, `d4_alphaiw_data.cpp`, `d4_corrections_data.cpp`
- Compiled units (CMakeLists.txt `curcuma_core_SRC`): `d4param_generator`, `d4_evaluator`, `d4_charge_model`, `d4_ncoord`

## Formula and users

`E_pair = -zeta_c6 * C6 * (s6*t6 + s8*r4r2_ij*t8)`, `t_n = 1/(r^n + R0^n)`, `R0 = a1*sqrt(r4r2_ij) + a2`.

| User | Code path | s6 | s8 | a1 | a2 (Bohr) | three-body |
|------|-----------|----|----|----|-----------|-----------|
| GFN-FF | `GenerateDispersionPairsNative()`, then `FFWorkspace::calcD4Dispersion` (own inline kernel; GPU `k_dispersion`) | 1.0 | 2.0 | 0.58 | 4.80 | off (`-gfnff.dispersion_atm`) |
| Native GFN2 | `XTB::calcDispersionEnergy` calls `D4Evaluator::computeEnergyAndGradient` + `computeATM` | 1.0 | 2.7 | 0.52 | 5.00 | s9 5.0, alpha 16, cutoff `-xtb.d4_atm_cutoff` |
| Native GFN1 | `D3ParameterGenerator::createForGFN1()` (D3, not D4) | 1.0 | 2.4 | 0.63 | 5.00 | none |

- GFN2 couples D4 self-consistently: `XTB::addDispersionPotential` adds `dE_D4/dq` each SCF iteration; the reference build runs once per geometry (`m_d4_prepared`, reset in `XTB::Calculation`)
- GFN2 GPU D4 runs in the backend contexts (`xtb_gpu_context.cu`, `xtb_hip_context.hip`, Vulkan), not through `D4Evaluator` (its empty `launchGpuKernel()` hook was removed on 2026-10-01)
- `-d4_charge_source` (GFN2, `xtb` scope): `mulliken` default (variational response), `eeq` (single-shot EEQ), `cpscf` (explicit Z-vector); see [docs/D4_Q_RESPONSE.md](../../../../docs/D4_Q_RESPONSE.md)
- GFN-FF uses `zetac6` from topology charges as a fixed per-pair prefactor (no `dq/dx` term)

## Invariants and traps

- `D4ParameterGenerator`'s PARAM defaults are GFN-FF values (a1 0.58, a2 4.80, s8 2.0); every caller must fill `D4Params` explicitly
- `D4Evaluator`'s guard on unset params is an `assert`, compiled out in `release/` (`-O3 -DNDEBUG`)
- The C6 reference matrix is built lazily for the elements present; `c6CacheCoversAtoms()` rebuilds it when a reused generator meets a new element
- GFN-FF pair list: built at `PAIR_BUILD_CUTOFF_BOHR` 60, evaluated at `PAIR_EVAL_CUTOFF_BOHR` 50; per-pair C6 must be refreshed after geometry changes (`-gfnff.dispersion_c6_update`)

## Open items

- The `d4_a1` PARAM help text says "GFN2-xTB: 0.63"; native GFN2 uses 0.52 (`xtb_native.cpp`, `D4Params` blocks)
- `GenerateParameters()` carries two in-code fallbacks for `d4_a1`/`d4_a2`: 0.58/4.80 for the pair `R0`, 0.44/4.60 for the exported `d4_damping` JSON; whether the latter is ever read is not checked
- GFN2 D4 residuals vs tblite: [docs/GFN2_D4_STATUS.md](../../../../docs/GFN2_D4_STATUS.md); `ctest -L d4_diag` is registered only when `release_tblite/dumps/*_gfn2.json` exist (none in `release/`)

### Open: native D3/D4 as a general add-on correction (note of 2026-06-26, full text in the archive)

- `DispersionMethod` (`-method d3|d4`) still routes to the external s-dftd3/cpp-d4 (`USE_D3`/`USE_D4`)
- `uff-d3` already evaluates native D3 pairs (`forcefield.cpp:180,354`); D3 presets: `createForGFN1/GFNFF/UFFD3/PBE0/BLYP/B3LYP/TPSS/PBE/BP86`, `createForMethod`
- D4 has no functional presets; only native GFN2 fills `D4Params`
- Potential plan:
  1. Add a native-backed add-on: re-point `DispersionMethod` at the native generators, or a `NativeDispersionMethod` adding `E_disp + grad E_disp` onto the host
  2. D4: functional presets (PBE0-D4, B3LYP-D4, ...) or a dftd4-style lookup table; decide the charge source (EEQ default, host Mulliken opt-in)
  3. A CLI flag on UFF/HF/DFT hosts (e.g. `-disp d3bj -disp_functional b3lyp` or `-d4`) folding the add-on into `calculateEnergy()`/`getGradient()`
  4. Validate against the legacy s-dftd3/cpp-d4 paths ([docs/TECHNICAL_DEBT.md](../../../../docs/TECHNICAL_DEBT.md) section 5) before deprecating them

## Verification

- `ctest -R "gfnff"` (GFN-FF D4 via FFWorkspace), `ctest -R d4_dedq` (`dE/dq`, EEQ `dq/dx`), `ctest -R xtb_gradient` (GFN2 analytic vs FD: H2O, CH4, NH3)

## References

- E. Caldeweyher et al., *J. Chem. Phys.* **150**, 154122 (2019) - D4 model
- C. Bannwarth, S. Ehlert, S. Grimme, *JCTC* **15**, 1652 (2019) - GFN2-xTB
- S. Spicher, S. Grimme, *Angew. Chem. Int. Ed.* **59**, 15665 (2020) - GFN-FF
- Fortran reference: `external/gfnff/src/gfnff_gdisp0.f90`; Ulysses parameters: `external/ulysses-main/core/src/parameters/D4par.hpp`

---

Previous version (measurements, per-phase status, refactoring history, removed 2026-10-01):
[docs/archive/DISPERSION_NOTES_2026-10.md](../../../../docs/archive/DISPERSION_NOTES_2026-10.md)
