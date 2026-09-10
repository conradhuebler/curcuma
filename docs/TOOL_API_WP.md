# WP: A programmatic tool API — results, log sinks, and a richer parameter registry

## Status (09.09.2026)

| WP | | |
|---|---|---|
| WP1 | Extend the parameter registry | **done** |
| WP2 | Annotate what is exposed | **done** — 13 modules, 65 enum values, 15 relevance conditions |
| WP3 | Thread-local log sink | **done** — sink plus 148 migrated prints |
| WP4 | Per-run state instead of global | **done** — verbosity, stop, output directory |
| WP5 | `Results()` | **partly** — the contract, `RMSDDriver`, `SimpleMD`; `opt`/`sp` and the rest are open |
| WP6 | A measurement capability | **done** — the capability and its test; the six CLI commands are not yet dispatched onto it |
| WP7 | Export the capability table | open |
| WP8 | Directed external potentials | **done** — three forms, live-settable, each carrying its accumulated work |
| WP9 | Free-energy protocols (work, PMF, FEP) | open |

Anything without a "done" marker is planned and does not exist. Nothing in this document describes existing behaviour
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

1. **No result contract.** *(WP5: addressed for RMSDDriver and SimpleMD.)* `CurcumaMethod` (`src/capabilities/curcumamethod.h:35`) has
   `start()` but no `Results()` (`grep -c "virtual json Results"` → 0). Result documents are
   assembled in `main.cpp`, e.g. the `.rmsd.json` block. Only two drivers expose results at
   all, with different types: `TrajectoryAnalysis::getResults() -> const json&`
   (`trajectoryanalysis.h:100`) and `UnifiedAnalysis::AnalysisThread::getResults() ->
   const std::vector<json>&` (`analysis.h:198`). In process there is no file to read back,
   so without `Results()` the result is simply unreachable.

2. **No log sink.** *(WP3: done.)* `CurcumaLogger` (`src/core/curcuma_logger.h:34`) is all-static and writes
   through `fmt::print` to stdout. Beyond it, 831 raw `std::cout` sites in `src/`
   (`main.cpp` 163, `analysis.cpp` 144, `simplemd.cpp` 81).

3. **Global verbosity state assumes single-threaded LIFO nesting.** *(WP4: done.)* `CurcumaMethod` saves
   `CurcumaLogger::get_verbosity()` in its constructors (`curcumamethod.cpp:49` and `:102`)
   and restores it in the destructor. Two concurrent instances on different threads overwrite
   each other, and one destructor restores a level captured on another thread. Every
   destructor also calls `printCitations()` to stdout.

4. **Cancellation goes through a file in the current directory.** *(WP4: done, `requestStop()`.)* `CurcumaMethod::CheckStop()`
   (`curcumamethod.cpp:263`) tests for a file named `stop` in the CWD. That is process-global:
   one stop aborts every concurrent run, and a leftover file aborts the next one at step 0.

5. **BMT writes into the caller's directory.** *(WP4: done, `pinOutputDir()`.)* Each run creates
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

**A bug found and fixed on the way in (09.09.2026) — larger than it first looked.** The
extractor silently dropped multi-line PARAMs whose help text is written as adjacent string
literals. The single warning it emitted pointed at `gfnff.h`, which suggested two lost
parameters; the true number was **22**, the rest disappearing without any warning at all.

The cause was the accumulate-and-match loop: it warned and *cleared the accumulator* as soon as
the text so far contained "PARAM" and any `)`, and these help texts contain one
(`'chloroform'). 'none' (default)`) long before the PARAM's own closing paren. The same parameter
written on one 436-character line (`xtbinterface.h:42`) went through fine, which is why the
source carries the note "PARAMs stay single-line".

Fixed by waiting for the PARAM's own parenthesis to close, counting only parentheses outside
string literals, and by skipping comment lines — `gfnff.h` keeps a removed parameter as
`// PARAM(eeq_distance_cutoff, ...) - REMOVED` for the record, which was reported as malformed on
every build. Count 662 → **684**, nothing lost. What came back: 18 `eeq_solver` parameters
(the whole PCG/extrapolation block), `gfnff.solvent`, `gfnff.solvent_model` and two
`eeq_refactor_*`. None of them had a `ConfigManager` default or appeared in `-export_config`
before.

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

**The capability is done** (`d06b153e`), with `test_measurement` checking it against geometry
worked out by hand rather than against a previous run of the same code. The dependency above
was resolved by taking **only** the `GeometryTools` half of `09a451c4` (`62a8bb34`): those are
pure additions and change no existing number, while the `BestFitRotation` fix in the same
commit is RMSD numerics another line of work is active in and stays where it is.

Three decisions the implementation settled:

- **How many atoms a kind needs is a property of the kind.** `requiredAtoms()` is the only
  place that knows, and an angle over two atoms is refused with the count in the message
  rather than measured as something else.
- **A centroid is reported as a position**, not reduced to a scalar with statistics over it.
  The mean of three coordinates is a number nobody should be able to quote.
- **Selections are resolved per frame.** `"F2"` on frame 900 is not the same index set as on
  frame 0, and for a measurement over a trajectory that is the correct reading, not a bug.

**Still open: the six commands in `main.cpp` are not dispatched onto it yet.** They keep their
own argv parsing and frame loops, so for now the capability is a second way to compute the same
thing rather than the only one — which is the situation this WP exists to end. It was left out
deliberately: `main.cpp` is the collision hotspot named under Branch discipline, and the
rewiring is mechanical once that file is quiet.

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

**Done** (`619b43f8`), in `src/capabilities/external_potentials.{h,cpp}`. What the
implementation settled that the plan above left open:

- **Work is accumulated per potential**, the sum over steps of F·dr. WP9 stage 1 asked for
  exactly this and it belongs to the potential rather than to a separate pass. For
  `constant_force` the energy alone is origin-dependent (E = −ΣF·r), so the work is the
  quantity that means anything, and the header says so rather than leaving it to be found out.
- **A selection that matches nothing is refused at parse time.** A bias that silently acts on
  nothing is worse than one that will not start.
- **The evaluation has no reference to SimpleMD**: it takes a geometry and a gradient. That is
  deliberate, because the optimiser wants the same potentials — see below.

**Open, and the natural next step: the same potentials in a geometry optimisation.** Every
optimiser already has an external-force hook, but it takes a flat force vector, i.e. the same
transient injection. Restrained optimisation is the more useful case of the two:

| Form | In MD | In an optimisation |
|---|---|---|
| `centroid_harmonic` | hold a fragment while it is heated | place a guest in a cavity, relax everything else around it |
| `distance_harmonic` | draw two fragments together | a **relaxed scan**: step `r0` outward, minimise at each value, and the binding curve falls out |
| `constant_force` | steer, and integrate the work | **questionable**: E = −F·r is unbounded, so the minimiser translates the set along the force and there may be no minimum at all |

The accumulated work has no meaning for an optimiser either: there is no trajectory, only
whatever path the minimiser took. So the wiring is not symmetric — the two harmonic forms
belong in the optimiser, `constant_force` should be refused there, and `work` reported only
for MD.

**Attempted and reverted, with what it showed.** A first pass wired all three backends:
`OptimizerDriver` gained the potentials plus a bias evaluated in Angstrom units, and each
backend converted at its own gradient site — ANCOpt hands out coordinates in Angstrom while
LBFGSpp and the native optimisers work in Bohr. It built. It did not work, and it was taken
back out rather than left in place, because a restraint that is wired and silently does
nothing is the exact failure mode this layer keeps running into.

What is known, so the next attempt does not start over:

- **The plumbing is reached.** A probe in the bias showed it called with the right shapes
  (18 coordinates for 6 atoms, one potential) from the backends.
- **Nothing moved.** With `k = 5 Eh/Å²` pulling two waters from 6.0 Å towards `r0 = 4.0 Å`,
  the centroid distance stayed at exactly 6.000 Å on all three, and the reported energy
  carried none of the bias (0.0038 Eh, the plain UFF value, where the restraint alone would
  be several Eh).
- **The backends fail differently**: `native_lbfgs` ran its 400 iterations and returned six
  atoms; `lbfgspp` and `ancopt` returned zero iterations and an empty molecule. All three
  reported `success = false`. That an unconverged run hands back the *input* geometry would
  explain the unchanged distance on all three at once, and is the first thing to check.
- **The measurement is the right one.** A wrong Angstrom/Bohr factor does not crash, it moves
  the restrained minimum, so a restrained distance landing on `r0` is what proves each
  conversion. The optimiser half of `test_external_potentials.cpp` was written around exactly
  that and can be restored from this branch's history.

**Deliberately not in scope:** finding the cavity. That is a separate question and probably not
curcuma's, at least not first -- see the guiding-scenario section in qurcuma's
`docs/WP-llm-tool-layer.md`.

### WP9 — Free-energy protocols (work, PMF, FEP)

Asked for a binding *free* energy, a model can today compute a potential-energy difference and
say so. Everything beyond that is out of reach, and the gap is not one tool but three different
depths of change. Listing them together because the cheapest is nearly free once WP8 exists and
the deepest touches the energy calculators.

**What already exists and is easy to overlook.** `rmsd_mtd` is a real enhanced-sampling bias, not
a toy: RMSD-space metadynamics with an exact gradient of the bias potential, well-tempered
reporting, a deposition stride and a cap on stored hills (`simplemd.h`, "RMSD-MTD", 13 PARAMs).
The walls are a second declarative biasing form. Restart files carry a run's state
(`curcuma_restart.json`). So biasing infrastructure is present; what is missing is a *collective
variable* other than RMSD, the bookkeeping around a biased run, and the statistics on top.

**Stage 1 — the work along a steered pull.** With WP8's `constant_force` in place, a pull is
already a steered MD. What is missing is that nobody adds up what it did: the accumulated work
`∫F·dr` of the external potentials over the trajectory, reported per step and in `Results()`.
That single number turns a pull into a Jarzynski estimate once several pulls are averaged, and
into Crooks with the reverse direction. Cheap, and it makes the actuator quantitative instead of
merely visible.

**Stage 2 — a restrained collective variable, and histograms.** Umbrella sampling along a
distance coordinate is the usual route to a host–guest PMF and needs **no** Hamiltonian scaling
at all — only a harmonic restraint on a CV and the CV's value per step:

| Piece | Meaning |
|---|---|
| CV definition | distance between two atom-set centroids first; the same selection grammar as WP8 |
| `cv_harmonic` | restraint of that CV to a target with a force constant, both settable mid-run |
| CV trace | the value each step, in `Results()`, so windows can be histogrammed |
| WHAM/MBAR | the estimator over a set of windows |

The restraint form is WP8's `distance_harmonic` with a target that moves, so stages 1 and 2 share
their machinery. The estimator itself is arithmetic over collected histograms and does not have to
live in the force loop; it could equally be a separate capability that reads window files.

**Stage 3 — λ-coupling and soft-core (FEP/TI).** Scaling the non-bonded interactions between two
atom groups by a coupling parameter, with soft-core to keep the potential finite as atoms vanish,
plus `∂V/∂λ` per step for TI and `ΔU` between neighbouring λ for BAR/MBAR. This is the only stage
that reaches into the energy calculators rather than into `SimpleMD`, and it has to be done per
method: GFN-FF's non-bonded terms, and the semiempirical methods separately, where "decoupling"
is not even well defined for the SCF part. That is a research question as much as an
implementation, and it is the reason this stage is named last rather than first.

**Done when:** stage 1 — a steered run reports its accumulated external work in `Results()`, and
two runs in opposite directions can be combined; stage 2 — a set of umbrella windows on a
centroid distance produces a PMF whose barrier is reproducible from the run configurations;
stage 3 — not scoped yet, deliberately.

**Deliberately not in scope:** treating a potential-energy difference as a free energy anywhere in
the output. If ΔG is not what was computed, the result says ΔE.

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
