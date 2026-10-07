# CLAUDE.md - test_cases/

All tests run by `ctest`: C++ unit/integration tests, CLI end-to-end scripts and reference-value gates.
List them with `ctest --test-dir release -N`; select with `-R <regex>` or `-L <label>`.
Longer guide (Feb 2026, not re-checked): [TESTING.md](TESTING.md).

## Where tests are registered

- Root `CMakeLists.txt` (`add_subdirectory(test_cases)` onward): `molecule_comprehensive`, `AAAbGlc_*`, `confscan_*`, `energy_methods`, `parameter_io`, `confstat_*`
- `test_cases/CMakeLists.txt`: the `test_*.cpp` executables (`gfnff_*`, `xtb_gradient_*`, `d4_dedq`, `gradient_unit_contract`, `md_*`, ...), `gfnff_val_*` (one per `reference_data/*.ref.json`)
- `test_cases/cli/CMakeLists.txt`: CLI tests, `add_cli_test(CATEGORY NAME)` registers `cli_<category>_<name>`; GPU categories only when the backend is built
- `test_cases/sqm_reference/CMakeLists.txt`: `sqm_val_*` (1e-8 Eh energy gates for gfn1/gfn2), `sqm_solv_*`, `sqm_gbsa_*`, `gfnff_solv_*`, `d4_diag_*`

## CLI tests (`cli/<category>/<NN_name>/run_test.sh`)

- `add_cli_test` copies the test directory into the build tree at configure time and runs it there, so the source tree stays clean
- Start from `cli/template_test.sh`: source `test_utils.sh` via `SCRIPT_DIR`, then `run_test()` and `validate_results()`
- Validate the science (values with tolerance, structure counts, drift), not only the exit code; output files via `find_output_file` (BMT-aware)
- Reference values and their origin: `cli/GOLDEN_REFERENCES.md`; documented bugs: `cli/KNOWN_BUGS.md`

## Unit and integration tests

- Test molecules come from the registry (rule below); `test_energy_methods.cpp` uses `AAA-bGlc/A.xyz` (117 atoms; host plus methyl beta-D-glucopyranoside, directory renamed from AAA-bGal on 2026-10-01) instead
- New unit test `test_<name>.cpp`: `add_executable`, `target_link_libraries(... curcuma_core test_molecule_registry)`, `add_test` with `TIMEOUT` and `LABELS` in `test_cases/CMakeLists.txt`
- Document the tolerance and where each reference value comes from (program, version, settings)

## Molecule Registry - MANDATORY Rule

**NEVER hardcode molecule geometry in test files.** All test molecules live in `core/test_molecule_registry.cpp`:
```cpp
#include "core/test_molecule_registry.h"
curcuma::Molecule mol = TestMolecules::TestMoleculeRegistry::createMolecule("CH4", false); // false = keep Angstrom
```
- Available: `H2`, `HCl`, `OH`, `Cl2`, `HCN`, `H2O`, `H2O_dimer`, `NH3`, `O3`, `CH4`, `CH3OH`, `CH3OCH3`, `C6H6`, `monosaccharide`, `triose`
- To add one: edit `core/test_molecule_registry.cpp` (atoms in Angstrom) and link `test_molecule_registry` to the test target
- Why: hardcoded geometry gives geometry-dependent pass/fail, duplicates data and is hard to audit

## Structure library

- `test_cases/structures/` holds every test structure with provenance (program, method, level, or the literature source); rules and naming in its `README.md`, validation with `python3 scripts/structlib.py check`
- New structures only through `structlib.py add`; a changed geometry is a new id; no program output in the library
- Migration state: CLI test directories with a `structures.txt` copy their structures from the library (`add_cli_test`); the others still read local copies (`legacy_paths` in the manifest map them); all consumers read the library (CLI tests via `structures.txt`, the others through staged copies in the build tree, see its README); no structure file may be added anywhere else

## Traps

- `createMolecule(name)` defaults to `scale_coordinates = true` and returns **Bohr**; pass `false` for Angstrom
- CLI scripts take `$CURCUMA`, else the first of `release/`, `debug/`, `build/`, `release_rocm/`, ... in the project root, whatever build tree ctest runs in
- `cli/errors/*` hardcode `release/curcuma`; `sqm_val_*` use the build's own binary (`$<TARGET_FILE:curcuma>`)
- Configure-time copies and globs: an edited `run_test.sh` or a new `reference_data/*.ref.json` is picked up only after CMake re-runs
- `add_cli_test` injects `PROJECT_ROOT` after the line `#!/bin/bash`; a script with another shebang does not get it
- `WILL_FAIL` marks `sqm_val*` molecules not yet at 1e-8 (`_GFN1_XFAIL`/`_GFN2_XFAIL`: `complex`, plus `He2` for gfn1); they pass while the gap persists
- `d4_diag_*` and `confscan_molalign` are registered only when their inputs (`release_tblite/dumps/`, `molalign` binary) exist
- `*/03_invalid_*` CLI tests check graceful fallback; `curcumaopt/03_invalid_method` runs a valid gfnff optimisation. Error paths are covered by `cli/errors/`
- `AAAbGlc incr` exists in `AAAbGlc.cpp` but its ctest entry is commented out
- `energy_methods` reference comments still name TBLite / external GFN-FF although `gfn1`, `gfn2`, `gfnff` now resolve to native code; not checked whether the test passes

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

### Future Development
- [ ] Fix curcumaopt JSON null-Fehler (HÖCHSTE PRIORITÄT)
- [ ] Erweitere wissenschaftliche Validierung für alle Tests
- [ ] Implementiere "expected failure" Pattern für invalid_method Tests
- [ ] Füge Performance-Benchmarks hinzu (Regression Tests)
- [ ] Erweitere test_molecule.cpp für geplantes SOA/AOS Refactoring

### Testing Philosophy
- **Test-Driven Development**: Schreibe Tests vor Refactorings
- **Scientific Accuracy**: Validiere Physik, nicht nur Code
- **Regression Prevention**: Capture current behavior als Golden Reference
- **Documentation**: Tests sind lebende Dokumentation der Features

---

Previous version (status counts, per-test descriptions, completed issues, removed 2026-10-01):
[docs/archive/TEST_CASES_NOTES_2026-10.md](../docs/archive/TEST_CASES_NOTES_2026-10.md)
