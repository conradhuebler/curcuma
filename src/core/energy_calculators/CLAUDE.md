# CLAUDE.md - src/core/energy_calculators/

Every energy method behind one interface. `EnergyCalculator` (`src/core/energycalculator.*`) owns one
`std::unique_ptr<ComputationalMethod>` that `MethodFactory::create(name, json)` builds.

## Layout

- `computational_method.h` - the abstract interface
- `method_factory.{h,cpp}` - table-driven registry (`methodTable()`, one `MethodDescriptor` per family)
- `gpu_plugin.{h,cpp}` - runtime `dlopen` loader for `libcurcuma_{cuda,rocm,vulkan}.so` ([docs/GPU_PLUGIN_STARTUP.md](../../../docs/GPU_PLUGIN_STARTUP.md))
- `qm_methods/` - native xTB/EHT/NDDO and external QM wrappers (own CLAUDE.md)
- `ff_methods/` - FFWorkspace engine, GFN-FF, UFF/QMDFF/CG (own CLAUDE.md)
- `dispersion/` - D4 data, charge model and evaluator (own CLAUDE.md)
- `ENERGY_SYSTEM_OVERVIEW.md` - CLI-to-library parameter flow (Oct 2025, not re-checked)

## ComputationalMethod contract

- Pure virtuals include `setMolecule`, `updateGeometry`, `calculateEnergy(bool gradient)`, `getGradient`, `getCharges`, `hasGradient`, `getEnergyDecomposition`
- `getGradient()` returns **Eh/Angstrom** (the header only says "appropriate units"); guarded by the `gradient_unit_contract` ctest
- Optional virtuals (`getCN`, orbital data, per-term energies, warm start) have no-op defaults

## Method resolution (`methodTable()`, first matching name wins)

- `gfn1`, `gfn2`: `NativeXtbMethod`; `-gpu cuda|rocm|vulkan` goes through the plugins
- `eht`: `EHTMethod`; `pm3`, `am1`, `mndo`, `pm6`: native `NDDOMethod` (their rows precede the Ulysses row)
- `gfnff`, `gfnff-fast`: native GFN-FF; `uff`, `uff-d3`, `qmdff`, `cg`, `cg-lj`: `ForceFieldMethod`
- Explicit external providers: `xtb-gfn1/2` (TBLite, else xtb binary), `tblite-gfn1/2`, `ipea1` (TBLite), `xtb-gfnff`, Ulysses names (`ugfn2`, `rm1`, `*-d3h4x`, ...), `d3`/`d4` (s-dftd3/cpp-d4), ORCA composites
- Unknown names raise `MethodCreationException` with "Did you mean" suggestions; `curcuma -methods` prints the table

## Adding a method

1. A `ComputationalMethod` subclass with its `PARAM` block
2. One row in `MethodFactory::methodTable()`
3. Its JSON sub-scope name in `MethodFactory::methodParameterScopes()`, the list EnergyCalculator, opt, MD, Hessian and ConfSearch forward

## Parameter flow and traps

- `XTBMethod`, `TBLiteMethod`, `UlyssesMethod` and `ForceFieldMethod` (for `ForceFieldGenerator`) hand a `ConfigManager` to the code they wrap
- Native wrappers merge JSON onto in-code defaults (`NativeXtbMethod::getDefaultConfig()`, `value(key, fallback)` in EHT, GFN-FF)
- For GFN-FF a PARAM default change alone has no runtime effect (root CLAUDE.md); for the other native wrappers not checked
- `EnergyCalculator`'s JSON constructors re-attach only the scopes in `methodParameterScopes()`; its `ConfigManager` constructors re-attach none
- Verbosity: wrappers map CurcumaLogger levels 0-3 onto the libraries (`XTB_VERBOSITY_*`, `tblite_set_context_verbosity` 0/1/3)

---

Previous version (ConfigManager migration phases, code examples, status tables, removed 2026-10-01):
[docs/archive/ENERGY_CALCULATORS_NOTES_2026-10.md](../../../docs/archive/ENERGY_CALCULATORS_NOTES_2026-10.md)
