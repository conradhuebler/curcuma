# Technical debt: parameter routing and configuration (curcuma)

Commit `c3f570f` (2026-10-31), clean tree. Read-only analysis, no repo file changed.

## How the claims were checked

- **Registry contents**: there is no `release/` build in this checkout. I compiled `scripts/param_parser/main.cpp` as is and
  ran it on all 191 `src/**/*.h` (the same glob CMake uses). Result: 774 PARAM definitions in 29 modules, plus 14
  `MODULE_INFO` blocks. A Python script parsed the generated header. Scratch files (not in the repo):
  `.../scratchpad/pp/parameter_registry.h`, `defs.json`, `analysis.txt`, `vals.txt`.
- **Runtime routing**: I built a harness from the unmodified `config_manager.cpp`, `parameter_registry.cpp`, the
  generated registry and the verbatim text of `CLI2Json` (`main.cpp:646-1124`). `CurcumaLogger` and the three
  `MethodFactory` scope functions were replaced by stubs. The harness runs real invocations (`.../scratchpad/harness/`).
  Findings marked **[run]** were reproduced this way. **[read]** means code reading only. **[unverified]** means a
  plausible consequence that was not executed.
- A full `curcuma` binary was not built and no ctest was run.

---

## 1. Findings ranked by severity

### Severity HIGH: wrong or silently ignored user input

**H1. `-traj` and `-dock` flags miss their module, because the keyword-to-module map is hand-kept and incomplete.** [run]
- `CLI2Json` maps CLI keywords to modules with a hard-coded 9-entry map (`main.cpp:727-737`). The registry already
  records the same mapping through `MODULE_INFO` (`traj`->`trajectoryanalysis`, `dock`->`docking`, `distance/angle/...`
  ->`measurement`). `ParameterRegistry::modulesForCommand()` (`parameter_registry.cpp:328`) exists but has **no caller**.
- `curcuma -traj t.xyz -stride 5` gives `{"traj":{"stride":5}}` and a warning ("registered in multiple modules ...
  Keeping flag in command module"). `ConfigManager("trajectoryanalysis", ...)` (`trajectoryanalysis.cpp:32`) skips
  object-valued keys, so stride stays at the registry default 1. **The flag is ignored.**
- `curcuma -dock -host A.xyz -guest B.xyz -step_x 4` prints three warnings and then works by accident (auto-routing
  moves the flags to `controller["docking"]`). But `executeDocking` reads `controller["dock"]["host"]`
  (`main.cpp:2246`), which is now empty. So `initializeBMT` is skipped and docking output goes to the CWD, against the
  mandatory BMT rule. That part is [read].

**H2. A global `-method` overwrites `rmsd.method` inside ConfScan.** [run]
- `curcuma -confscan c.xyz -method gfn2` gives `rmsd.method = "gfn2"`.
- Cause: `CLI2Json` copies global values into `controller["global"]` (`main.cpp:1119-1122`). The multi-module
  `ConfigManager` then fills **every** module from `global` for each key the module does not set itself
  (`config_manager.cpp:245-256`). `rmsd.method` is the RMSD alignment method; the global `method` is the energy method.
- ConfScan then warns "Unknown RMSD method 'gfn2'" and falls back to subspace (`confscan.cpp:401-404`). The user's
  energy method is not applied either, because ConfScan reads `energy_method` (`confscan.cpp:254`).
- The single-module constructor has the same global fallback (`config_manager.cpp:59-68`). In general, any module
  parameter whose name matches a global (`method`, `threads`, `verbosity`, `charge`, `spin`) is overwritten by the
  global value.

**H3. The CLI number detector accepts numeric prefixes, which turns file names into numbers or crashes.** [run]
- `std::stod(next)` (`main.cpp:858`) succeeds on "2026run.json", "1abc.log" or "3.xyz", because it parses only a
  prefix. Results:
  - `-export_run 2026run.json` gives `export_run: 2026.0`, then `std::string = controller["export_run"]`
    (`main.cpp:3071`) throws `type_error.302`.
  - `-bak 1abc.log` gives `bak: 1.0`.
  - `-reference_xyz 3.xyz` gives 3.0.
- Only flags registered as `String` are protected (the guard at `main.cpp:877-894`). `import_config`, `export_run`,
  `bak` and the other 21 hard-coded globals (except the 5 also registered in modules) are not.
- `-charge 1e999` throws an uncaught `std::out_of_range` from `stod`. Only `invalid_argument` is caught
  (`main.cpp:860`).
- The final `else` branch (`main.cpp:917-927`) is unreachable.

**H4. `-opt ... -write_trajectory 0`, the form shown in the `-opt` help text, aborts the optimisation.** [read + nlohmann behaviour run]
- The help prints `-write_trajectory <0|1>` (`main.cpp:1999`). The CLI stores 0.0, and `MergeJson` passes it through.
- `OptimizationDispatcher::optimizeStructure` calls `merged_config.value("write_trajectory", true)`
  (`optimizer_factory.cpp:295`). nlohmann throws `type_error.302` "type must be boolean, but is number" (confirmed in
  the harness). The exception is caught at `main.cpp:2156` as "optimizer threw exception".
- Bool coercion is written by hand in at least 3 places with different rules: `main.cpp:1868-1875`,
  `main.cpp:3058-3062`, `optimizer_driver.cpp:166-170`. There is no central coercion.

**H5. Optimizer defaults have three layers, and "is this a user override?" is guessed by comparing against defaults.** [read]
- The three default sources:
  - the registry module `opt` (`curcumaopt.h:67-107`)
  - `OptimizerInterfaceJson` (`optimizer_interface.cpp:28-77`): `lbfgs_min_step` 1e-20, `method` "lbfgspp"; here
    `method` means the optimizer, while in `opt` it means the energy method
  - `LBFGSppOptimizer::GetDefaultConfiguration` (`lbfgspp_optimizer.cpp:306-320`): `lbfgs_max_line_search` 2, while
    the registry says 20
- `executeOptimization` merges the registry defaults into the config (`main.cpp:2030-2031`). The dispatcher then
  treats every value that differs from `GetDefaultConfiguration()` as an explicit user override
  (`optimizer_factory.cpp:277-286`).
- Consequences:
  - On `-opt`, the registry default (20) replaces the LBFGSpp default (2) as if the user had set it.
  - A user value that happens to equal the default is dropped when a `-convergence_preset` is active.
  - Direct callers that skip the registry merge (SimpleMD, RMSDTraj, PolymerBuild, Docking: `simplemd.cpp:791`,
    `rmsdtraj.cpp:609`, `polymerbuild.cpp:2609/2784/3199`, `docking.cpp:360/368/404`) get the other set.
    [unverified which values each caller ends up with]

**H6. Hard-coded fallback defaults diverge from the PARAM defaults (generalises trap #31).** [read]
- 736 `value(key, default)` / `get<T>(key, default)` calls in `src/`; 665 have a literal default. 67 of those differ
  from every registry default of the same name. Most are deliberate "not given" sentinels (confsearch `-1.0`). The
  substantive ones:
  - `FragmentScfDriver::configureXtb` defaults `d4_charge_source` to `"eeq"` (`xtb_fragment_scf.cpp:54,58`). The PARAM
    (`native_xtb_method.h:32`) and the main path (`native_xtb_method.cpp:86,105`) use `"mulliken"`, and the comment says
    "Mirror NativeXtbMethod::applyConfig". So gradients from `-large_system_mode fragments` use a different D4
    q-response than the normal path unless the user sets it. [read, not run]
  - `D3ParameterGenerator` reads with fallbacks `d3_s8 1.0, a1 0.4, a2 4.0` (`d3param_generator.cpp:43-46`, and 7 more
    sites with variants `0.40/4.20`, `0.44/4.60`). The PARAMs are `2.0/0.58/4.80` (`d3param_generator.h:36-43`). Which
    one wins depends on whether the caller's `ConfigManager` was built for module `d3param`. GFN-FF itself uses
    `2.85/0.80/4.60` (`gfnff_method.cpp:12281-12283`); forcefieldgenerator uses `2.7/0.45/4.0`
    (`forcefieldgenerator.cpp:49-53`). [read; the effective value per path was not run]
  - `tbliteinterface.cpp:64`: `solvent_model` fallback 2, PARAM 0 (deliberate auto-enable, but undocumented in the
    PARAM).
  - `simplemd.cpp:907`: `topology_mode` fallback `"auto"`, the GFN-FF PARAM and constructor say `"constant"`
    (`gfnff_method.cpp:582,605`). This only gates the RATTLE refusal, which is therefore skipped when the scope is
    missing.
- GFN-FF: `gfnff_method.cpp` has 143 fallback reads of 96 keys. 47 keys duplicate a `gfnff` PARAM default by hand.
  49 keys are not registered at all, among them `dispersion`, `hbond`, `coulomb`, `repulsion`, `use_simplified_pbo`,
  `hb_cell_list_min_atoms`, `cache_topology`, `print_timing`, `threads`. Those have no help, no flat routing and no
  `-export_run` entry. The constructor merges a private default table (`gfnff_method.cpp:574-587`), not the registry.

**H7. Capabilities that receive only their own subdocument lose method scopes and globals.** [read]
- `executeCasino` passes `MergeJson(registry casino, controller["casino"])` (`main.cpp:2314-2317`). Casino builds its
  `EnergyCalculator` from that (`casino.cpp:137-138`). So `-gfnff.*`, flat-routed `-cn_cutoff_bohr`, `-xtb.*` and
  `-scf.*` never reach Casino's method.
- Casino also always finds `verbosity` (registry default 2) and calls `CurcumaLogger::set_verbosity(2)`
  (`casino.cpp:59-75`). The `verbose` branch is dead and the global `-verbosity` default 1 is overridden.
- `executeRMSDTraj` does the same (`main.cpp:2353-2358`): only `controller["rmsdtraj"]`, no `rmsd` module, no globals.

### Severity MEDIUM: inconsistency that will cause the next H-class bug

**M1. Unknown or misspelled keys are accepted silently at every layer.** [run]
- `-sp w.xyz -temprature 400` gives `{"opt":{"temprature":400}}` without a message.
- `ConfigManager` adds unknown keys as new entries (`config_manager.cpp:127-133`, `300-306`); `has("temprature")` is
  true.
- `MergeJson` adds them lowercased (`global.h:213ff`).
- Only the `scf` scope warns about unknown keys (`method_factory.cpp:631-634`). TODO.md:20 lists the open-shell flags
  `-multi`, `-uhf` that are "silently ignored".
- Positional tokens after a flag value are dropped silently: `-threads 4 5` loses the 5 [run].

**M2. The single-module `ConfigManager` takes every top-level scalar of the controller into the module.** [run]
- `config_manager.cpp:71-81`: `ConfigManager("opt", {"opt":{},"temperature":900})` gets `has("temperature") == true`.
- Combined with M1, any capability sees keys meant for others. The ConfigManager for `eeq_solver` is built from the
  whole GFN-FF parameter set (`gfnff_method.cpp:661-671`) and loads all 110 `gfnff` defaults through
  `canonicalModule` [run].

**M3. Global parameters exist in up to four copies, and the list is hand-kept outside the registry.**
- 21 names in `global_params` (`main.cpp:695-712`). 14 of them (all `gpu*`, `no_bmt`, `bak`, `noprogress`,
  `export_run`, `import_config`) are registered nowhere: no help text, no type, no export default.
- The other 5 (`verbosity`, `threads`, `method`, `charge`, `spin`) are also registered as ordinary PARAMs in 5-9
  modules each, with different defaults (`verbosity`: casino 2, opt 1, trajectoryanalysis 2, native_lbfgs 2, gfnff 0)
  and different types (`spin` is Int in opt/simplemd/confsearch, Double in xtb/tblite).
- The `spin` help says "Spin multiplicity" in opt (`curcumaopt.h:75`, `main.cpp:1982`) and "0=singlet, 1=doublet" in
  simplemd (`simplemd.h:813`). These are two different semantics in user-facing text.
- Each global ends up in `controller[param]`, `controller["global"][param]`, `controller[module][param]`, and for
  `-sp`/`-md` also `controller[keyword][param]` [run: `-sp -spin 1` yields four copies].
  `confsearch.cpp:105` and `simplemd.cpp:874` copy `global` once more by hand.
- After `-import_config`, an edited top-level value and a stale `global` copy can disagree.

**M4. `controller[keyword]` duplicates are written but almost never read.**
- `main.cpp:1109-1111` stores the command's parameters a second time under `"sp"` or `"md"`. The only reader of a
  keyword-named scope is the docking case above (`main.cpp:2246`); `controller["md"]`/`["sp"]` have no reader [read,
  grep].
- `-export_run` writes these stale duplicates into the replay file.

**M5. Name ambiguity is structural.**
- 91 canonical names exist in more than one module. 6 have conflicting types (`spin`, `lbfgs_line_search` Int/String,
  `damping` Double/String, `fragments` String/Bool, `accuracy` Int/String, `solvent_model` String/Int). 35 have
  conflicting defaults.
- 128 flat spellings (name or alias, lowercased) have more than one owner, and 20 aliases equal another module's
  canonical name (for example `scf` aliases `max_iterations`/`diis_start` collide with `opt`). This is the root of
  TODO.md:14: `-sp -max_iterations 3` stays in `opt` and never reaches the SCF [run].
- The resolution order lives in imperative code, not data: command module, then method scope, then `scf`, then a
  single owner, else warn-and-keep (`main.cpp:949-1011`).
- The type-conflict check in `validateRegistry` is compiled out (`#ifdef DEBUG_PARAMETER_REGISTRY`,
  `parameter_registry.cpp:373`).
- `validateRegistry()` failure only prints a warning (`main.cpp:2882-2884`).

**M6. Static JSON defaults instead of the registry (violates the "Parameter Definition Standards").**
- `CasinoJson` (`casino.h:99-125`, used at `casino.cpp:33,44`) duplicates the registered `casino` module. Its comment
  `casino.h:34` claims it was replaced.
- `ConfStatJson` (`confstat.h:35`, keys `Cutoff/Temp/Threshold/Method`, not snake_case, no module).
- `AnalyseNCIPlotJson` (`analysenciplot.h:24`).
- `TDAEngineJson` (`tda_engine.h:34`, duplicates `analysis.topological_*`).
- `OptimizerInterfaceJson` (H5).
- `OptJsonPrivate` and the `CurcumaOptJson` macro (`curcumaopt.h:40-56`): no users, dead.
- Manual merges in `main.cpp` that the doc's own "Golden Rule" forbids (`PARAMETER_FLOW_ARCHITECTURE.md`:
  "execute* functions: Just pass controller"): opt `main.cpp:2030`, casino `main.cpp:2314`, rmsdtraj `main.cpp:2353`.
- `Json2KeyWord`: 123 calls in 12 capability files. It throws `int -1`, not a `std::exception`, and does no alias
  resolution, so the registered aliases of Casino, Docking, ConfSearch and others work only where a `ConfigManager` is
  used.

**M7. Key normalisation is re-implemented about 30 times, with different rules.**
- `::tolower` transforms: 20 in `config_manager.cpp`, 2 in `config_manager.h`, 4 in `parameter_registry.cpp`, 4 in
  `global.h`, 8 in `method_factory.cpp`. Hyphen removal happens at only 17 of those sites.
- `ConfigManager::getFromModule` (`config_manager.h:249-261`), the nested and dot-notation paths of `findKey`
  (`config_manager.cpp:450-569`) and `findOwnerModules` behave differently.
- The comment "canonicalize hyphens to underscores" is wrong: the code *erases* `-`, so a JSON key `max-iterations`
  becomes `maxiterations` and does not match `max_iterations` (the CLI path is safe because `main.cpp:778` replaces
  `-` with `_`).
- `findKey` has four fallback strategies, among them a `module.key` prefix and a double-nested `config[module][key]`.
  Its warnings go to `std::cerr` under an inverted `#ifndef DEBUG_CONFIG_MANAGER` (`config_manager.cpp:436,461`).

**M8. `ConfigManager::get<T>(key, default)` swallows every exception.**
- `config_manager.h:268-274`. A wrong JSON type (for example a bool stored as 0.0) is reported as "not found" and
  replaced by the literal default [run: `E throws: ... not found`].

**M9. The sub-scope lists are maintained in three places and contain dead entries.**
- `MethodFactory::methodParameterScopes()` (`method_factory.cpp:568-573`) lists `uff`, `qmdff`, `eht`, `d3`, `d4`;
  no registry module has those names.
- `scope_modules` in `CLI2Json` (`main.cpp:721-723`) adds `modern_optimizer`, which does not exist.
- `parameterScopeForMethod("eht")` returns the unregistered `eht` (`method_factory.cpp:601`).
- `ConfigManager("energycalculator", ...)` is used in 4 places but `energycalculator` has no PARAM block. Its
  re-attachment of stripped scopes is technical-debt item D-15/D-1.

### Severity LOW: hygiene

- **L1 Naming.** 3 non-snake canonical names: `confsearch.startT/endT/deltaT` (`confsearch.h:264-266`). 76 non-snake
  aliases out of 241 aliases on 228 parameters. 23 aliases differ from their canonical name only by case
  (`Charge`, `Step_X`, `sLE`, ...) and are redundant because all matching is case-insensitive. Alias
  `update-rotation` (`rmsd`) cannot be reached from the CLI. No alias carries the `deprecated` annotation (0 of 774),
  so there is no path to remove old names.
- **L2 Unreachable flat flags.** `-verbose`, `-v`, `-silent`, `-quiet` are intercepted before routing
  (`main.cpp:781-805`). The canonical PARAM `verbose` of `analysis`, `casino`, `polymerbuild` and `trajectoryanalysis`
  can be set only in dotted form.
- **L3 Dead or misnamed PARAMs** [read].
  - `ancopt.opt_level` (`ancopt_optimizer.h:377`): the code reads `optimization_level` (`ancopt_optimizer.cpp:454`),
    which is unregistered.
  - `ancopt.anc_lbfgs_history` (`:383`) is never read; the member is fixed at 12 (`ancopt_optimizer.h:172`).
  - `d3param.d3_rs6/d3_rs8` are never read.
  - `polymerbuild.ljlm_lm_max_iter/ljlm_lm_tolerance/ljlm_k_bond` have no string-literal reader [likely dead; a
    constructed key was not ruled out].
  - The 15 `analysis.topological_*` hits are false positives: they are read as `"topological.x"` and depend on
    `findKey`'s dot-to-underscore fallback (`config_manager.cpp:520-538`), another hidden mapping.
- **L4 PARAM annotations are not enforced.** 18 PARAMs declare `enum=`, 16 `min/max`, 15 `unit`. They are used only by
  `printHelp`. `parameter_validation.h` (`ParamValidation::*`) has no include and no caller; it is dead code, yet it is
  listed in `src/core/CLAUDE.md:12`.
- **L5 Dead code in `CLI2Json`.** `setNestedJsonValue` always splits dotted keys, so `key` never holds a flat
  `"a.b"` name. The `is_flat_dotted` branch (`main.cpp:1034-1062`) and the `pname.find('.')` check (`:946`) are
  unreachable [read].
- **L6 Dotted versus flat precedence.** `-gfnff.cn_cutoff_bohr 5.5 -cn_cutoff_bohr 6` gives 6 regardless of order: the
  flat reroute is applied after the dotted value (`main.cpp:1014-1019`) [run]. The more specific form should win.
- **L7 `-export_run`** detects touched modules with `getForModule(it.key())`, without `canonicalModule`
  (`main.cpp:3107`), so `eeq_solver` gets no defaults (documented in PARAMETER_SYSTEM.md). It also exports `global`
  and `controller[keyword]` duplicates.
- **L8 Stale documentation.**
  - The `PARAMETER_SYSTEM.md` "Migration Status" table (`:413-421`) says casino/simplemd/opt/confscan/hessian/rmsd are
    "Pending"; all are registered.
  - `-export-config` and `-import-config` are spelled with hyphens there.
  - `PARAMETER_FLOW_ARCHITECTURE.md` gives line numbers from Oct 2025 (`CLI2Json` "184-342", now 675-1124).
  - `config.get<int>("MaxIterations")` is said to find `max_iterations` (`:292`). It does, but only via the alias
    `MaxIter`, not by case folding.
  - `ENERGY_SYSTEM_OVERVIEW.md` is marked "not re-checked".
- **L9** The hand-written `-opt` help (`main.cpp:1966-2015`) duplicates the registry help. Today the printed defaults
  agree with `curcumaopt.h`, but nothing keeps them in step.

Already listed in `docs/TECHNICAL_DEBT.md`; I confirmed they still apply but did not re-investigate them: D-15/D-1 (energycalculator
scope stripping), D-16 (stringly config keys), Q-37 (`setParameters` bypasses the registry), Q-54 (EHT has no PARAMs),
Q-57 (mutable ConfigManager), X-M8 (scoped-config fallback written 5 times), F-Q1 (unregistered `use_two_phase_eeq`).

---

## 2. Proposal: unified parameter routing

### Target architecture

One resolution pass in `main`, one typed result, and consumers that only read.

```
argv --> Tokenizer --> Resolver ------------------------> RunConfig --> ConfigManager views --> code
          (typed by    (CommandSpec + registry:             (per module: defaults
           registry)    flag -> (module, canonical name))    + user values, validated;
                                                             one "global" scope)
```

1. **The registry is the only schema.**
   - Every parameter that code reads is a PARAM, including the 21 globals. Add a module `global` with
     `MODULE_INFO(... commands = all)`; `verbosity`, `threads`, `method`, `charge`, `spin` and `gpu*` live there once.
     Remove their duplicates from opt/simplemd/confsearch/... or keep them only as declared "inherits global"
     annotations (`inherit=global`).
   - Add a `scope=` annotation for sub-scopes (`scope=method` for gfnff/xtb/tblite/..., `scope=submodule` for rmsd
     inside confscan).
2. **A `CommandSpec` derived from `MODULE_INFO`.**
   - Fields: command, primary module, sub-modules (for example confscan -> {confscan, rmsd}), and whether method
     scopes apply.
   - It replaces `keyword_to_module` (`main.cpp:727`), `scope_modules` (`:721-723`),
     `MethodFactory::methodParameterScopes()` and the hand-written module lists in the `ConfigManager` constructors.
3. **Tokenizer.**
   - Value typing comes from the resolved definition: String/Path stays a string, Int/Double need a full-match parse
     (`std::from_chars`), Bool accepts `true/false/1/0/yes/no`.
   - Unregistered flags keep today's heuristics but with a full-match numeric test.
4. **Resolver.** A fixed precedence, written as data and documented in one table:
   1. dotted `-mod.name`
   2. primary module of the command
   3. sub-modules of the command
   4. the scope of the selected `-method` (then `scf`)
   5. `global`
   6. a single other owner, with a warning
   7. otherwise an error with "did you mean" suggestions (`StringUtils::find_closest_matches` already exists,
      `method_factory.cpp:744`)

   `-allow_unknown_flags` keeps the legacy pass-through for a transition period.
5. **RunConfig.**
   - `{global:{...}, <module>:{...}}`, with registry defaults filled for every module the command and the method touch.
   - Enum/min/max annotations are validated once.
   - No copies at the top level and no `controller[keyword]`.
   - `-export_run` writes it verbatim; `-import_config` feeds the Resolver (dotted keys), so imported values go through
     the same validation.
6. **ConfigManager becomes a read-only view** over one module plus `global`, for names declared `inherit=global` only.
   - No flat top-level fallback, no prefix or double-nesting guesses.
   - `get<T>(key)` fails on a type error instead of reporting "not found".
   - One `normalizeKey()` function shared with the registry.
7. **Methods and capabilities read only through the view.**
   - `GFNFF::GFNFF(const json&)`, `NativeXtbMethod`, `FragmentScfDriver`, `D3/D4ParameterGenerator`, Casino, ConfStat,
     NCI and TDA get their defaults from the registry. The literal fallbacks are deleted, which ends trap #31 by
     construction.

### Migration order (lowest risk first)

| Step | Change | Risk | Guard |
|---|---|---|---|
| 0 | Move `CLI2Json`/`setNestedJsonValue`/`mergeJsonRecursive` from `main.cpp` into `src/core/cli_parser.{h,cpp}` without changing them; add a ctest with golden JSON for about 30 invocations (the probes in this report) | none | the golden test is the acceptance criterion for every later step |
| 1 | Registry lint script (naming, type conflicts, dead PARAMs, literal-fallback vs PARAM mismatches) in `scripts/`, run from `check_docs.py` | none | report only |
| 2 | Safe first steps (section 3) | low | golden test plus `ctest` |
| 3 | Warn on unknown keys: in the Resolver (CLI), in `ConfigManager` for keys that are neither PARAM nor declared global; warn only | low | there are no numeric changes |
| 4 | Remove duplicate defaults in the numeric code: GFN-FF constructor merges `getDefaultJson("gfnff")` first, then delete the literal fallbacks one by one; the same for xtb, fragment SCF, D3/D4 generators | medium-high (energies) | `scripts/refset_regression.py` MOR41/GMTKN55 before and after, bitwise for unchanged defaults; delete `.topo.json` caches |
| 5 | `CommandSpec` from `MODULE_INFO`; remove `keyword_to_module`, `scope_modules`, `methodParameterScopes()` duplicates | medium | golden test |
| 6 | Module `global`; `ConfigManager` loses its top-level-flat and blanket-global fallbacks (fixes H2, M2) | high (many capabilities read globals through those fallbacks) | golden test, full ctest, `-no_bmt` and BMT runs of every `cli_*` group |
| 7 | Single copy of globals; delete `controller[keyword]` and top-level copies; readers use `ConfigManager("global")` | high | grep for `controller["threads"]`, `controller.value("verbosity"...)` (15 call sites) and similar first |
| 8 | Migrate Casino/ConfStat/NCI/TDA/RMSDTraj/Optimizer defaults to the registry; drop `Json2KeyWord` (123 calls) and `MergeJson` in capabilities; the optimizer gets one default layer | medium per capability | per-capability ctests |
| 9 | Enforce the PARAM annotations (enum/min/max), add `deprecated` on old aliases, remove them after one release | low-medium | lint |

---

## 3. Safe first steps (mechanical, small commits)

Each is independent. "Behaviour" says what changes for a user.

1. **Strict numeric detection in `CLI2Json`** (`main.cpp:855-864`).
   - Replace `std::stod(next)` with `size_t pos; std::stod(next, &pos); if (pos != next.size()) isNumber = false;` and
     catch `std::out_of_range` (treat as a string, or error).
   - Delete the unreachable `else` (`:917-927`).
   - Behaviour: `-export_run 2026run.json`, `-bak 1abc.log` and `-charge 1e999` stop being numbers or crashing.
     Numbers like `1e-6`, `-3` and `.5` are unchanged.
2. **Add `dock`->`docking` and `traj`->`trajectoryanalysis`** to `keyword_to_module` (`main.cpp:727-737`), and read
   `controller["docking"]["host"]` in `executeDocking` (`main.cpp:2246`).
   - `controller[keyword]` stays as the compatibility copy (`main.cpp:1109`).
   - Behaviour: `-traj -stride 5` takes effect; docking stops warning and gets its BMT directory.
   - Leave the `measurement` commands for step 5 of the migration (wider effect).
3. **Delete dead code.**
   - `src/core/parameter_validation.h` (0 includes) and its mention in `src/core/CLAUDE.md:12`.
   - `OptJsonPrivate`, `getCurcumaOptJson` and the `CurcumaOptJson` macro (`curcumaopt.h:40-56`, no users; re-grep
     before deleting).
   - `"modern_optimizer"` in `scope_modules` (`main.cpp:723`).
   - The unreachable `is_flat_dotted` branch (`main.cpp:1034-1062`) and the `pname.find('.')` check (`:946`).
4. **Fix the dead or misnamed ancopt PARAMs.**
   - Make `ancopt_optimizer.cpp:454` read `opt_level` (alias `optimization_level`), or rename the PARAM.
   - Read `anc_lbfgs_history` into `m_anc_lbfgs_history` in the ancopt config loader.
   - Default values are unchanged (0 and 12), so there is no numeric change.
5. **Make `FragmentScfDriver::configureXtb`'s `d4_charge_source` fallback `"mulliken"`** (`xtb_fragment_scf.cpp:54,58`),
   as its comment says.
   - Small, but it changes fragment-mode gradients: run the fragment SCF ctests and one gradient check
     (`scripts/gradient_compare.py`).
6. **Remove the 23 case-only aliases** (list in `.../scratchpad/pp/analysis.txt`; files `curcumaopt.h`, `rmsd.h`,
   `docking.h`, `confscan.h`, `simplemd.h`, `confsearch.h`).
   - `resolveAlias` and `ConfigManager` match case-insensitively, and `findOwnerModules` lowercases. Run
     `make GenerateParams` and the golden test.
7. **Rename `confsearch.startT/endT/deltaT`** to `start_temperature/end_temperature/delta_temperature`, with the old
   names as aliases (`confsearch.h:264-266`).
   - Change the readers in `confsearch.cpp` to the new names. First check that ConfSearch reads them through a
     path that resolves aliases (`ConfigManager`), not `Json2KeyWord`.
8. **`ConfigManager` hygiene without behaviour change.**
   - Route the two `std::cerr` warnings in `findKey` (`config_manager.cpp:436-440,461-465`) through
     `CurcumaLogger::warn` and fix the inverted `#ifndef`.
   - Correct the "hyphens to underscores" comments (5 sites) to say "hyphens removed".
9. **`-export_run` touched-module detection** uses `ParameterRegistry::canonicalModule()` (`main.cpp:3107`) and skips
   `global` and `controller[keyword]` duplicates (`main.cpp:3074-3123`).
10. **Turn the debug-only type-conflict check in `validateRegistry`** (`parameter_registry.cpp:373-391`) into a
    build-time warning in `scripts/param_parser/main.cpp`. It runs at `make GenerateParams`, prints nothing at runtime,
    and the 6 conflicts become visible to developers.
11. **Documentation corrections** (no code): the `PARAMETER_SYSTEM.md` migration table and hyphenated command names;
    the `PARAMETER_FLOW_ARCHITECTURE.md` line numbers; the `-spin` help text in `curcumaopt.h:75` and `main.cpp:1982`
    should match the actual 2S convention used by `Molecule::setSpin`. Check the convention before writing it down
    (TECHNICAL_DEBT Q-2 states `m_spin = multiplicity - 1`).

Not "safe" despite being small, because they change user-visible results: H2 (removing the global fallback), H5
(optimizer default layers), H6 (literal fallbacks in numeric code), H7 (Casino/RMSDTraj receiving the full controller).
They belong to migration steps 4, 6 and 8.

## Not checked

- No full binary was built and no ctest was run. The harness stubs `MethodFactory::parameterScopeForMethod` and
  `readsSharedScfScope` with a simplified gfn1/gfn2 rule.
- The effective D3 parameters per call path (H6), the optimizer values each direct caller ends up with (H5), and the
  Casino energy configuration (H7) were derived by reading only.
- The literal-fallback scan is regex-based: calls spanning lines and keys built at runtime are missed, and the 67
  "mismatches" include sentinel values that are intentional.
