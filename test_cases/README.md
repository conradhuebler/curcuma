# test_cases

Everything `ctest` runs, plus the structures and reference data behind it. Rules for contributors and for AI
sessions are in [CLAUDE.md](CLAUDE.md); the longer guides [TESTING.md](TESTING.md) and
[GFNFF_TESTING.md](GFNFF_TESTING.md) are older and not re-checked.

## Running

```sh
cd release
ctest -N                           # list all tests
ctest -R cli_confscan_ --output-on-failure
ctest -j4                          # everything (145 s on the development machine)
```

Build and test in `release/`. State of 2026-10-09: 337 registered tests, 332 pass, two fail since before this
cleanup (`test_orca_interface`, `xtb_cpscf`), three are disabled (`cli_curcumaopt_02_gfn2_single_point`,
`cli_sqm_10_gfn2_optimization_smoke`, `cli_sqm_11_gfn2_provider_check`). Passing tests show that the code does what the tests
expect, not that a method is correct (root `CLAUDE.md`).

## Layout

| Directory | Contents |
|---|---|
| `structures/` | The structure library: every molecule a test reads, with charge, spin and provenance (`structures/README.md`) |
| `unit/` | C++ unit and integration tests (one executable each) and the molecule registry that reads the library |
| `cli/` | End-to-end tests of the command line, one directory per scenario (`cli/README.md`) |
| `scripts/` | Script tests registered in ctest (`check_*.py`, `large_system_modes.sh`, `test_parameter_io.sh`) |
| `sqm_reference/` | Gates of the native GFN1/GFN2 and D3/D4 code against TBLite and dftd4 reference values, with their generators |
| `reference_data/` | Reference JSON of the GFN-FF validation (`*.ref.json`) and the xtb outputs behind them |
| `reorder/`, `rmsd/` | Sources of the two small RMSD test programs |
| `MOR41-testset/`, `GMTKN55-testset/`, `s30lci_test_set/` | Benchmark sets; MOR41 and GMTKN55 are fetched on demand (`scripts/fetch_testset.py`) |

## Adding a test

- A new structure goes into the library first: `python3 scripts/structlib.py add ...`, then
  `python3 scripts/structlib.py check`. No structure file anywhere else.
- Unit test: source in `unit/`, target and `add_test` in `CMakeLists.txt`, molecule from the registry.
- CLI test: directory `cli/<category>/<NN_name>/` with `run_test.sh` and `structures.txt`, registered with
  `add_cli_test(<category> <NN_name>)` in `cli/CMakeLists.txt`.
- Document the tolerance and where a reference value comes from (program, version, settings).
