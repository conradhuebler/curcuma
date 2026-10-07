# Structure library

One place for every molecular structure that a test uses as input or reference. `manifest.json` describes each
structure; `scripts/structlib.py` validates it (`python3 scripts/structlib.py check`, also run by the pre-commit
hook). Large benchmark sets are not stored here, they are fetched (see `sets` in the manifest).

## Status of the migration

- **Done**: the library, the manifest, the checker, and the import of the structures that were scattered over `test_cases/`
  (67 distinct structures from 187 files; the other 120 were byte-identical copies). Every entry lists its `legacy_paths`.
- **Phase 2 (tests read the library)**: done for the CLI tests and for every other consumer in `test_cases/`, `scripts/` and the
  registry. CLI category notes follow. A migrated test directory has a `structures.txt` and no
  local copy of the structure; `python3 scripts/structlib.py usage` lists which tests use which structure. All CLI test directories
  are migrated (51 directories, 122 local copies replaced); in 49 of them the inputs in the build tree were compared byte for
  byte before and after in three build configurations, `cli/cg/01_single_point` and `cli/gfnff_gpu/01_gfnff_gpu_singlepoint`
  are not registered in any of them and are only covered by the hash check of `structlib.py migrate`. Left in `cli/`: 18
  tracked program outputs (`*.centered.xyz`, `*.reordered.xyz`, `optimized.xyz`, an optimiser trajectory) and the deliberately
  invalid `errors/02_invalid_xyz_format/invalid.xyz`.
- **Other consumers**: `test_cases/CMakeLists.txt` and `sqm_reference/CMakeLists.txt` stage the library at configure time into the
  build tree under the former file names (`curcuma_stage_structures` in `structures.cmake`, for example
  `build/test_cases/molecules/larger/caffeine.xyz`); the tests read the staged copies, so caches such as `*.topo.json` do not end up in
  the library. The `geometry_file` field of 35 reference JSON files, the paths in `TestMoleculeRegistry`, five Python helpers and
  three benchmark scripts point to the library; scripts can use `structlib.py path <id>` and `structlib.py stage`. All 64 remaining
  legacy copies were removed (`structlib.py retire`).
- **Not done**: the `TestMoleculeRegistry` still carries hard-coded atoms; five tracked program outputs outside `cli/` remain.
  The full test suite gives the same result as before the migration (327 pass, 3 known failures, 7 disabled). 23 tracked files that look like program
  output and one file that cannot be parsed (the invalid-format input of `cli/errors/02_invalid_xyz_format`) are listed
  under `not_imported` in the manifest and are not part of the library.
- **Provenance of the legacy entries** was only derived from what the repository records (comment lines, xtb `.out`
  files, `reference_data/*.json`, two README files). Everything else is `unknown`, which is a statement of ignorance,
  not a guess. `python3 scripts/structlib.py report` prints the current numbers.

## Layout

`<class>/<id>.xyz`, one directory per class, ids unique over the whole library:

| class | contents |
|---|---|
| `atoms` | a single atom or ion |
| `small` | one fragment, up to 12 atoms |
| `medium` | one fragment, 13 to 60 atoms |
| `large` | one fragment, more than 60 atoms |
| `clusters` | two or more non-bonded fragments (covalent-radius criterion), dimers, solvent clusters |
| `metals` | contains a transition metal |
| `bulk` | more than 500 atoms: polymer chains, solvent boxes |
| `ensembles` | multi-frame files (conformer sets, trajectories used as input) |
| `cg` | coarse-grained systems (`.vtf` or element 226 beads) |

## Rules

1. **Format**: plain `.xyz` (single frame unless class `ensembles`), coordinates in **Angstrom**, charge and spin in the
   manifest. Other formats (`.vtf`) only where the test needs that format.
2. **The comment line (line 2) carries the provenance** for every new structure, machine readable:
   `id=<id> charge=<q> spin=<unpaired electrons> level=<program-version/method> source=<kind>[:<reference>]`.
   `structlib.py add` writes it. Legacy files keep their old comment line; the manifest is authoritative for them.
3. **Provenance is mandatory** and says what is actually known. `kind` is one of
   `optimized` (program and method required), `program_output` (program), `literature` / `database` / `experimental`
   (reference: DOI, data set and entry, license), `constructed` (description of how the geometry was built) or
   `unknown`. **`unknown` is only allowed for legacy entries**; new structures must state their origin.
4. **Optimised geometries record the level**: program, version, method, basis or parameter set, solvent, convergence
   criterion, and energy and gradient norm at that level when available. Evidence (the output file) is referenced, not
   copied into the library. A structure from elsewhere (literature, database, another group's reference geometry)
   says so in `reference` and names its license.
5. **Reference calculations at a geometry are not the provenance of the geometry.** An xtb single point on a structure
   (`*.out`) or a reference JSON goes into `reference_calculations`, the optimisation level into `provenance`.
6. **A geometry never changes.** The `sha256` is checked. A changed geometry gets a new id with `supersedes`, because
   reference values depend on the exact coordinates. The same molecule in another orientation is a new id with
   `variant_of`, it is not merged automatically.
7. **No program output** in the library (`.out`, `.log`, optimiser trajectories, `*.reordered.xyz`, result tables).
   Output that a test needs as an expected result belongs to that test.
8. **Size**: at most 1.5 MB per file (legacy files above that give a warning). Larger sets are fetched with
   `scripts/fetch_testset.py` and listed under `sets`.
9. **Names say what the structure is**, lowercase `[a-z0-9_.+-]`, built as `<system>[.<variant>]`:
    - `<system>` is the chemical name or a short descriptor, with explicit counts instead of codes:
      `peo201x2-water1500` (two PEO chains of 201 units in 1500 water molecules), `urea400-water1000`, `water8-cluster`,
      `aaa-bgal` (host-guest: `<host>-<guest>`), `hydronium-water2`. Never `2x`, `_final`, `new`, `test`, `input`.
    - Ensembles end in `-frames<N>` (what the frames are is stated in `description`).
    - `.<variant>` says how a variant differs from the plain system and is only added when variants exist:
      the level (`.gfn2-relaxed`, `.xtb-gfn2-opt`, `.uff-analytic-grad`), an idealised geometry (`.ideal-c2v`,
      `.ideal-d6h`), a measured parameter (`.r0.741` for a bond length in Angstrom), a conformer number
      (`.conf76`, as in the source ensemble), or the purpose (`.d-shell-validation`, `.opt-start`).
    - A name that rests on inference (composition, a file or directory name, a comment line) says so in
      `description`; the old names stay in `aliases` so that a search for `polymer_2x` finds the new entry.
    - `.alt<k>` marks a legacy variant whose difference to the main variant is not recorded anywhere; it keeps
      `needs_name` until someone names it.
10. **Roles**: `equilibrium`, `non-equilibrium`, `transition-state`, `ensemble`, `stress` (deliberately difficult), `invalid-input`
    (deliberately malformed, for error-path tests), `unspecified` (legacy only).

## Using a structure in a CLI test

Put a `structures.txt` into the test directory, one line per structure: `<id>` (file name `<id>.xyz`) or
`<id> as <file name>` when the script expects another name. `add_cli_test` copies the library file byte for byte into
the build directory of the test at configure time. Do not keep a local copy of a library structure next to it.
`python3 scripts/structlib.py migrate <test dir>` converts an existing test directory (it refuses when a local file
differs from the library file) and `check` verifies the lists. The legacy import
(`scripts/structlib_import_legacy.py`) cannot be repeated once tests are migrated.

## Adding a structure

```sh
python3 scripts/structlib.py add path/to/mol.xyz --id my_mol --class small --charge 0 --spin 0 \
    --role equilibrium --kind optimized --program xtb --version 6.7.1 --method gfn2 \
    --convergence "opt tight" --energy -5.070 --evidence docs/...
python3 scripts/structlib.py check
```

## Manifest fields

`id`, `class`, `file`, `format`, `formula` (Hill), `natoms`, `frames`, `charge`, `spin`, `role`, `size_bytes`,
`sha256`, `provenance{kind, program, version, method, basis, solvent, convergence, energy_eh, gnorm, reference,
license, description, evidence, level_note, comment}`, `reference_calculations[]`, `derived_from`, `supersedes`,
`description`, `analysis` (derived properties, e.g. the sugar configuration of a guest, with the method that produced it), `aliases`, `variant_of`, `needs_name`, `legacy`, `legacy_paths[]`. Top level: `sets[]` (fetched benchmark sets), `not_imported`.

## Commands

`check` (errors and warnings, exit code 1 on errors), `report` (counts per class and provenance, share of
`unknown`; `--list-unknown` names them), `usage [ID...]` (which tests read which structure, derived from the
legacy paths and from CMake and script text), `add`.
