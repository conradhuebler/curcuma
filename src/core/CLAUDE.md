# CLAUDE.md - src/core/

Core data structures and shared infrastructure. All energy methods live in
`energy_calculators/` (own CLAUDE.md, with `qm_methods/`, `ff_methods/`, `dispersion/`).

## Ownership

| Area | Files |
|------|-------|
| Molecule | `molecule.{h,cpp}` (atoms, geometry, charges, fragments, XYZ/JSON export), `xyz_comment_parser.*`, `fileiterator.*` (multi-XYZ iteration) |
| Energy dispatch | `energycalculator.{h,cpp}`: one `std::unique_ptr<ComputationalMethod>` built by `MethodFactory::create()` |
| Parameters | `parameter_macros.h` (PARAM blocks), `parameter_registry.*`, `config_manager.*`, `parameter_validation.h` |
| Logging, citations | `curcuma_logger.*`, `citation_database.*`, `citation_registry.*` |
| Units, elements | `units.h` (`CurcumaUnit`, CODATA 2018), `elements.h`, `periodic_table.*` |
| Numerics | `curcuma_eigen_config.h`, `blas_threads.h`, `math_compat.h`, `portable_{erf,exp,log,acos}.h`, `charge_extrapolation.h`, `intra_parallel_context.h` |
| GPU infrastructure | `gpu_device_pool.*` (batch workers over GPUs), `gpu_fallback.*` (counts every CPU fallback) |
| Other | `solvation/` (GBSA, solvent tables), `functional_groups.*`, `topology.h`, `hbonds.h`, `form_factors.h`, `pseudoff.*`, `imagewriter.hpp` |

## Parameter system

- PARAM blocks in headers are scanned at build time into `generated/parameter_registry.h` (target `GenerateParams`)
- `ConfigManager(module, json)`: `get<T>(key[, default])`, case-insensitive, aliases via `ParameterRegistry::resolveAlias()`, dot keys (`"topological.save_image"`)
- New capabilities use the registry (root CLAUDE.md, "Parameter Definition Standards")

## Invariants and traps

- The PARAM scan globs `src/*.h` at configure time: a new header, or a PARAM in a `.hpp`/`.cpp`, is not seen until it is a `.h` and CMake has re-run
- `ConfigManager("energycalculator", controller)` drops method sub-scopes; the JSON constructor of `EnergyCalculator` re-merges those listed in `MethodFactory::methodParameterScopes()`
- `curcuma_eigen_config.h` must precede every Eigen include (pulled in by `pch_base.h` and `global.h`)
- Unit constants only from `units.h`; `-DUSE_LEGACY_UNIT_CONSTANTS=ON` restores the old per-site values ([docs/UNIT_CONSTANTS.md](../../docs/UNIT_CONSTANTS.md))
- `ForceField` keeps an auto parameter file `<input>.param.json` (`setParameterCaching()`); GFN-FF writes `<input>.topo.json`. Both are caches that can outlive a code change (root CLAUDE.md traps)

## Open items

- Memory use for large systems (>1000 atoms) listed as open; no measurement recorded here
- Molecule refactoring roadmap: [REFACTORING_ROADMAP.md](REFACTORING_ROADMAP.md), comment formats that must not break: [XYZ_COMMENT_FORMATS.md](XYZ_COMMENT_FORMATS.md)
- Older TODO notes (2025, not re-checked): [REFACTORING_TODO.md](REFACTORING_TODO.md), [UNIFIED_INTERFACE_TODO.md](UNIFIED_INTERFACE_TODO.md)

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

---

Previous version (history, completed items, performance notes, removed 2026-10-01): [docs/archive/CORE_NOTES_2026-10.md](../../docs/archive/CORE_NOTES_2026-10.md)
