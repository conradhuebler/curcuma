# CLAUDE.md - src/helpers/

Standalone programs outside the `curcuma` binary: interface smoke tests, micro-benchmarks and
Python re-implementations of GFN-FF terms. None of them is registered in `ctest`.

## Files and build status (root CMakeLists.txt)

| File | Target | Built when |
|------|--------|------------|
| `main.cpp` | `curcuma_helper` | always; modes `-statistic file.dat [-bins N] [-force]` (histogram), `-allxyz` (prints a note only, conversion commented out) |
| `parallel_scf.cpp` | `parallel_scf` | always; times `ParallelEigenSolver` on a 4000x4000 matrix |
| `ulysses_helper.cpp` | `ulysses_helper` | `USE_ULYSSES` (no `HELPERS` gate) |
| `xtb_helper.cpp`, `tblite_helper.cpp`, `gfnff_helper.cpp`, `dftd3_helper.cpp`, `dftd4_helper.cpp` | same name | `HELPERS` plus `USE_XTB` / `USE_TBLITE` / `USE_GFNFF` / `USE_D3` / `USE_D4` |
| `cli_test.cpp` | `cli_helper` | never: target commented out, includes the missing `src/tools/cli_parser.h` |
| `gfnff_test.cpp` | none | never: includes the missing `qm_methods/gfnff.h` |
| `imagewrite.cpp`, `storage_bench.cpp`, `polymer_topo.cpp`, `gfnff_term_validator.cpp` | none | never; standalone `main()` programs (Eigen image writer, parameter-lookup storage benchmark, matrix topology analysis, GFN-FF term check) |
| `gfnff_reference_validator.py`, `gfnff_term_validator.py`, `validate_ch3oh.py` | - | Python scripts re-implementing GFN-FF formulas from the Fortran reference |

## Traps

- `HELPERS` is not a declared `option()`; pass `-DHELPERS=ON` explicitly to get the interface helpers
- `curcuma_helper` prints its usage as `curcuma_tools ...`; the built binary is `curcuma_helper`
- The reference-set runners for GFN-FF/xTB validation live in `scripts/` (root CLAUDE.md), not here

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Testing priorities, benchmarking requirements, and development tool specifications to be defined by operator/programmer*

---

Previous version (status notes, removed 2026-10-01): [docs/archive/HELPERS_NOTES_2026-10.md](../../docs/archive/HELPERS_NOTES_2026-10.md)
