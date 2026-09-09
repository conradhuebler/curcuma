# WP: A programmatic tool API — results, log sinks, and a richer parameter registry

## Status (Sep 2026): planned, nothing implemented yet

Every work package below is open. Nothing in this document describes existing behaviour
unless it is explicitly marked as such under "Starting point". Line numbers were measured
on this branch (`llm-core`, based on `origin/master` at `a991adcb`) and will drift.

## Why

qurcuma is growing an LLM tool layer: a registry of named operations with typed parameters
and structured results, driven from a chat dock, the command palette and headless tests.
That layer talks to curcuma **in process** (qurcuma links `curcuma_core` + `curcuma_cap`);
the `QProcess`-on-the-CLI path stays for ORCA and xtb and is not extended for curcuma.

In process, four things a subprocess would have provided for free have to be built here:
stdout is not captured by a pipe, the working directory is shared, global state is shared
across threads, and there is no `kill`. That is what this WP is about.

None of it is LLM-specific. A result contract, a log sink and a self-describing parameter
registry equally serve the CLI's own help output, `-export_config`, the Python bindings
(`origin/claude/add-python-interface-…`, whose history records the same problem: "Replace
broken bindings with honest placeholders for CLI-only features"), and GUI form generation.

## Starting point (measured on this branch)

What already exists and works:

- `ParameterRegistry` (`src/core/parameter_registry.h:29`) with `getForModule()` and
  `getDefaultJson()`, fed by `PARAM(...)` macros scraped at build time by
  `scripts/param_parser/main.cpp` (wired in `CMakeLists.txt:1435-1462`).
  30 `BEGIN_PARAMETER_DEFINITION` blocks (one of them, `module`, is the documentation
  example in `parameter_macros.h`), 684 `PARAM` lines across the headers.
- `ConfigManager` (`src/core/config_manager.h:60`), merging registry defaults with user JSON.
- `CAPABILITY_REGISTRY` (`src/main.cpp:2398`), 36 commands with description, category and
  accepted formats.
- The `-import_config` / `-export_run` round trip (`docs/CLI_ROUND_TRIP.md`).
- `GeometryTools` (`src/tools/geometry.h`) and `TrajectoryStatistics`
  (`src/capabilities/trajectory_statistics.h:45`, mean/stddev/min/max/median/variance plus
  `exportStatistics() -> json`).

What is missing, with the evidence:

1. **No result contract.** `CurcumaMethod` (`src/capabilities/curcumamethod.h:35`) has
   `start()` but no `Results()` (`grep -c "virtual json Results"` → 0). Result documents are
   assembled in `main.cpp`, e.g. the `.rmsd.json` block. Only two drivers expose results at
   all, with different types: `TrajectoryAnalysis::getResults() -> const json&`
   (`trajectoryanalysis.h:100`) and `UnifiedAnalysis::AnalysisThread::getResults() ->
   const std::vector<json>&` (`analysis.h:198`). In process there is no file to read back,
   so without `Results()` the result is simply unreachable.

2. **No log sink.** `CurcumaLogger` (`src/core/curcuma_logger.h:34`) is all-static and writes
   through `fmt::print` to stdout. Beyond it, 831 raw `std::cout` sites in `src/`
   (`main.cpp` 163, `analysis.cpp` 144, `simplemd.cpp` 81).

3. **Global verbosity state assumes single-threaded LIFO nesting.** `CurcumaMethod` saves
   `CurcumaLogger::get_verbosity()` in its constructors (`curcumamethod.cpp:49` and `:102`)
   and restores it in the destructor. Two concurrent instances on different threads overwrite
   each other, and one destructor restores a level captured on another thread. Every
   destructor also calls `printCitations()` to stdout.

4. **Cancellation goes through a file in the current directory.** `CurcumaMethod::CheckStop()`
   (`curcumamethod.cpp:263`) tests for a file named `stop` in the CWD. That is process-global:
   one stop aborts every concurrent run, and a leftover file aborts the next one at step 0.

5. **BMT writes into the caller's directory.** Each run creates
   `Basename.Keyword.TIMESTAMP/` unless suppressed.

6. **The parameter registry does not know enough for a schema.** Five gaps:
   - *No exposure tier.* Filtering on category `"Basic"` does not work: exactly 27 parameters
     carry it, spread over six files (`simplemd.h` 9, `curcumaopt.h` 6, `casino.h` 5,
     `ancopt_optimizer.h` 4, `gfnff.h` 2, `eeq_solver.h` 1). For `confscan`, `rmsd`,
     `docking`, `hessian`, `analysis` and `confsearch` the set is **empty**. Categories are
     help-output groupings (`Walls`, `SCF`, `Cost Matrix`, …), not a visibility level; both
     concepts share one field.
   - *No allowed values.* `thermostat` carries `berendsen|andersen|nosehover|csvr|none` in its
     help text; so does `wall_type` (`none|spheric|rect`).
   - *No units.* K, fs, Å, Hartree, Å⁻¹ all live in prose, inconsistently.
   - *No dependencies.* `wall_temp` only matters when `wall_type != none`, `rmsd_mtd_k` only
     when `rmsd_mtd`. qurcuma reimplements that gating by hand.
   - *No structured parameters.* `temp_regions` is a JSON array of objects, read at
     `simplemd.cpp:2206-2216`, but it is **not registered as a `PARAM`**. It never appears in
     `-export_config`, has no `ConfigManager` default, and cannot be described by any schema.
     `ParamType` knows only `String|Int|Double|Bool`, so selection expressions
     (FragString grammar) and paths are all "String".

7. **No module-level metadata.** `-list_modules` prints "simplemd (78 parameters)" with no
   description, and the command↔module mapping is recorded nowhere: only 10 of 36 command
   names match a module; `md`→`simplemd`, `dock`→`docking`, `traj`→`trajectoryanalysis`; 22
   commands have no module at all. A trap: the *command* `orca` runs an existing ORCA input
   file, while the *module* `orca` holds ORCA method parameters.

8. **Six geometry commands duplicate the same 60 lines.** `executeBond` (`main.cpp:181`),
   `executeAngle` (`:334`), `executeTorsion` (`:475`), `executeGyration` (`:1301`),
   `executeCentroid` (`:1499`), `executeDistance` (`:1539`) each redo argv parsing, existence
   checks, `FileIterator` iteration, index validation, the per-frame loop,
   `TrajectoryStatistics`, and three output formats by hand. None has a PARAM block, a schema
   or a `Results()`.

## Work packages

One WP is one buildable commit. The 29 existing ctests stay green after each.

### WP1 — Extend the parameter registry (additive)

`ParameterDefinition` (`parameter_registry.h:19`) gains `tier`, `allowed`, `unit`,
`min`/`max`, `requires`, `deprecated`/`replaced_by`. `ParamType` gains `StringList`, `Json`
and the semantic markers `Selection` and `Path`. A new
`ModuleDefinition { name, description, category, commands[] }` records the command mapping
from gap 7.

The `PARAM` macro is **already variadic** and expands to nothing
(`parameter_macros.h:41-43`); only the extractor is rigid. Its regex
(`scripts/param_parser/main.cpp:161`) demands exactly six arguments and ends on `\}\s*\)`.
It gains one **optional** trailing group holding a flat annotation string:

```cpp
// unchanged, all existing PARAM lines keep working:
PARAM(temperature, Double, 298.15, "Target temperature.", "Basic", {"T"})

// new, only where it pays:
PARAM(thermostat, String, "csvr", "Thermostat type.", "Thermostat", {},
      "tier=primary; enum=csvr|berendsen|andersen|nosehover|none")
PARAM(wall_temp, Double, 298.15, "Energy scale of the wall potential.", "Walls", {},
      "tier=advanced; unit=K; min=0; requires=wall_type!=none")
```

A flat `key=value; …` keeps the regex extractor simple. If the grammar grows past a handful
of keys, the 249-line scraper should be replaced by a real tokenizer — its own comment
history records a fixed regex typo, and it runs at build time over every header, so a
mistake there breaks the build or silently drops parameters.

**Done when:** every existing `PARAM` line is untouched and produces the same definition; a
test compares the per-module definition count before and after (`simplemd` 78,
`confsearch` 71, `analysis` 50, `confscan` 44, `gfnff` 35, `rmsd` 33); `-list_modules` shows
module descriptions.

**A bug to fix while in there (found 09.09.2026).** The extractor silently drops
multi-line PARAMs whose help text is written as adjacent string literals. Measured: the
generator reports 662 definitions and warns once about `gfnff.h` around line 310; `solvent` and
`solvent_model` (`gfnff.h:298` and `:304`) are **not** in the generated registry at all, so
`ConfigManager` has no defaults for them and `-export_config gfnff` does not list them.

The cause is the accumulate-and-match loop: it warns and *clears the accumulator* as soon as the
text so far contains "PARAM" and any `)`, and these help texts contain one
(`'chloroform'). 'none' (default)`) long before the PARAM's own closing paren. The same
parameter written on one 436-character line (`xtbinterface.h:42`) goes through fine, which is why
the source carries the note "PARAMs stay single-line".

Fix: only attempt the match, and only warn, once the PARAM's opening parenthesis is balanced,
counting parentheses outside string literals. The count then has to go 662 → 664.

### WP2 — Annotate what is exposed

The ~60 parameters of the commands a GUI or tool layer exposes: `tier`, `enum` for
`thermostat` / `wall_type` / `wall_potential` / `optimizer` / `topology_mode` / `gpu`,
`unit`, and `requires` for the wall, RMSD-MTD and temperature-ramp blocks. Register
`temp_regions` as a `Json` parameter (gap 6).

**Done when:** `-export_config simplemd` contains `temp_regions`; `-help_module simplemd`
shows units and allowed values.

### WP3 — Thread-local log sink

`CurcumaLogger` gains a sink registered per scope (RAII `push_sink`/`pop_sink`), records
carrying a scope id, so a caller can capture the output of one run while another runs
concurrently on a different thread. Then the loud `std::cout` paths in `simplemd.cpp` (81)
and `analysis.cpp` (144) move onto the logger. `main.cpp`'s 163 sites stay: those commands
run as processes, where a pipe captures them anyway.

**Done when:** a registered sink sees level and text; CLI output is unchanged in plain mode.

### WP4 — Per-run state instead of global state

Thread-local verbosity instead of the global save/restore (gap 3). A cancellation **handle
per run** instead of the CWD `stop` file (gap 4), with the file kept as a fallback for CLI
use. `setOutputDir()` / `no_bmt` so an embedded run writes nothing into the caller's
directory (gap 5).

**Done when:** two concurrent `CurcumaMethod` instances do not disturb each other; an
embedded run leaves no files behind.

### WP5 — `Results()`

`virtual json Results() const` on `CurcumaMethod`, implemented first on `RMSDDriver`, then
`sp` / `opt` / `md`. **Pure: it returns JSON and never touches disk.** Artifact writing stays
in the CLI layer, so an embedded caller does not litter the working directory.

**Done when:** the emitted `.rmsd.json` is byte-identical to before, but assembled in the
driver.

### WP6 — A measurement capability

`src/capabilities/measurement.{h,cpp}`: a `Kind` enum
(`Distance, Angle, Dihedral, Gyration, Centroid, RmsdToReference`), a PARAM block
(atoms, kind, unit, window, frame range), operating on a `Molecule` or a `FileIterator`
trajectory, using `GeometryTools` and `TrajectoryStatistics`, returning `Results()`.
The six commands from gap 8 become thin dispatches onto it.

**Dependency, and the reason this WP is not free on `master`:** `GeometryTools` on this
branch has `Distance` (`geometry.h:33`) and `Centroid` (`:40`) but **no `Angle`, no
`Dihedral`, no `GyrationRadius`**. Those were added on `reactff2` by commit `09a451c4`
("Add angle, dihedral and radius of gyration; fix the planar alignment guard"), whose own
message makes the same argument this WP makes. That commit also carries an
`RMSDFunctions::BestFitRotation` fix (a degeneracy guard that misfired on every planar
structure). Bringing it to `master` is the natural prerequisite, but it touches RMSD
numerics that other work is currently active in, so it is an explicit decision, not a
side effect of this WP.

**Done when:** `curcuma -angle`, `-torsion`, `-distance` produce the same numbers as before
on both a single structure and a trajectory; `-export_config measurement` yields a schema.

### WP7 — Export the capability table

`src/core/capability_table.{h,cpp}`: **metadata only** (name, description, category,
formats) plus a JSON exporter and `curcuma -capabilities-json`. `main.cpp` keeps its handlers
and references the table. The handlers take `argc/argv`, print to stdout and return `int`;
they are most of `main.cpp`'s 2960 lines and moving them buys nothing here.

Also: the banner and status lines of the introspection commands move to stderr, so
`-export_config <module>` writes pure JSON to stdout. Today the JSON follows a 20-line
ANSI-coloured banner on stdout, which any parser has to cut away.

**Done when:** a test asserts every dispatched command appears in the table;
`-export_config rmsd | jq .` succeeds.

### WP8 — Directed external potentials

The motivating case is agentic docking: a model loads a receptor and a guest and has to get the
guest into the cavity. To do that it needs an actuator -- it has to be able to *pull* on a set of
atoms and then watch what happens.

**What exists today.** `SimpleMD::applyExternalForces(const Geometry& forces)`
(`simplemd.h:269`) takes an additive per-atom force matrix, but `m_external_forces` is documented
as "cleared after use" (`:503`) and applies between two `step()` calls. It is an injection, not a
potential: the caller has to re-apply it every step, nothing is recorded in the controller, and a
run cannot be reproduced from its configuration. The only declarative potentials in SimpleMD are
the walls (`wall_*`, 11 PARAMs).

**What is missing.** A configured, persistent set of external potentials, evaluated inside the
force loop the way the walls are, described by parameters so that a run carries its own
definition and `-export_run` reproduces it. Three forms cover the cases in sight:

| Form | Meaning | Docking use |
|---|---|---|
| `constant_force` | fixed force on an atom set along a direction | steer the guest |
| `centroid_harmonic` | harmonic restraint of an atom set's centroid to a point | hold the guest at the cavity centre |
| `distance_harmonic` | harmonic restraint between two atom sets' centroids | draw two fragments together, or keep them apart |

Atom sets are named with the existing selection grammar (`FragString2Indicies`), so `"F1"` is the
guest and no new syntax appears.

**Two requirements that are easy to miss.**

*It has to be a structured parameter.* A list of potentials is an array of objects, and
`ParamType` currently knows only `String|Int|Double|Bool` -- this is the same gap that leaves
`temp_regions` unregistered (see the starting point above). WP8 therefore depends on WP1's `Json`
type, and is the second consumer that proves it.

*It has to be settable while the run is going.* An agent adjusts the pull and watches the
response; a value that can only be given at startup is useless for that. The pattern already
exists in SimpleMD for the thermostat setpoint and the wall parameters, which qurcuma changes
live through thread-safe setters -- external potentials need the same, not just a controller
entry read once in `Initialise()`.

**Done when:** a run started with two potentials in its controller reproduces from
`-export_run`; `-export_config simplemd` describes them; a potential added or changed mid-run
takes effect on the next step; and the forces show up in the energy/gradient bookkeeping rather
than being added behind its back.

**Deliberately not in scope:** finding the cavity. That is a separate question and probably not
curcuma's, at least not first -- see the guiding-scenario section in qurcuma's
`docs/WP-llm-tool-layer.md`.

## Branch discipline

This branch (`llm-core`) is based on `origin/master` and stays free of `reactff2`
dependencies, so its commits can go upstream cleanly. Each WP lands in **new files** where
possible, which keeps a later cherry-pick mechanical.

`reactff2-llm` is the integration branch: `reactff2` merged with `llm-core`. That is what
qurcuma builds against, because qurcuma calls `reactff2`-only APIs (`reactiveBonds()`,
`setWallTemp()`, `stopReason()`, the gyration radius). Merge direction is
`llm-core` → integration, never back.

`main.cpp` is the collision hotspot (WP6, WP7) and several other efforts touch it, so those
two WPs come last.
