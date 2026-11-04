# Technical debt: build system, tests, documentation, repo hygiene, architecture seams

Scope: read-only audit of `/home/user/curcuma` at `c3f570f` (shallow clone, 391 commits visible, so `git log` dates
before 2026-06-28 are the graft point and not real ages). No repository file was changed. Status: 🤖 AI-generated audit,
not human-reviewed. Complements `docs/TECHNICAL_DEBT.md`, which covers method code (GFN-FF, xTB, QM interfaces,
dispatch) and has almost nothing on build, tests or docs.

Not done / not verifiable here: no build tree exists (`release/`, `build/` absent) and configuring would run
FetchContent into `external/` inside the checkout, so **no cmake, make or ctest was run**. Everything below is from
reading files, `git grep`, and small Python scans. Items that would need a build to confirm are marked **(unverified)**.

Severity: high = can hide wrong results or failures; medium = maintenance/reproducibility cost; low = cosmetic.
Effort: S (< 1 h, one file), M (half day), L (multi-day).

---

## 1. Findings ranked

### 1.1 High severity

| # | Finding | Evidence | Effort |
|---|---|---|---|
| H1 | **CI smoke tests cannot fail on a crashing binary.** Steps pipe curcuma into `tail` with no `shell: bash` / `pipefail`; GitHub's default Linux/macOS shell is `bash -e {0}`, so the step's status is that of `tail`. CI never runs `ctest`. | `.github/workflows/ci-variants.yml:43-44,76,92`; `linux-gpu.yml:83-84` (`... -method gfnff $g 2>&1 \| tail -n 12`); only `ccpp.yml:52` sets `shell: bash`; `grep ctest .github/workflows scripts/ci` = 0 hits. GitHub default-shell semantics taken from the docs, not tested in a run. | S |
| H2 | **CI never runs the docs/structure checks.** `check_docs.py` and `structlib.py check` only run in an opt-in pre-commit hook. | `scripts/git-hooks/pre-commit`; no workflow mentions either script. | S |
| H3 | **Numerical path depends on what is installed, not on options.** `find_package(BLAS/LAPACK)` and `find_library(lapacke)` run unconditionally (`CMakeLists.txt:1483-1517`) and set `HAVE_BLAS`/`HAVE_LAPACKE` PUBLIC even with `USE_BLAS=OFF`; `eeq_solver.cpp:30` turns on `EIGEN_USE_BLAS` (GFN-FF EEQ charges), `rf_solver.cpp`/`lbfgs.cpp` turn on `EIGEN_USE_LAPACKE`. Two machines with the same cmake line can produce different EEQ/optimiser arithmetic; this matters for the project's bit-identity claims. | `CMakeLists.txt:1483-1497,1505-1508`; `git grep HAVE_BLAS -- src` (5 files). | M |
| H4 | **WILL_FAIL tests pass on a crash.** 3 xfail cases (`gfn2 complex`, `gfn1 He2`, `gfn1 complex`) are registered with `WILL_FAIL TRUE` for CPU, CUDA, ROCm and Vulkan-DC variants; ctest inverts any non-zero exit, so a segfault or a missing reference file also reads as "expected failure". | `test_cases/sqm_reference/CMakeLists.txt:314-315,431-432,492-493,524-525`, lists at `:292-293`. Better: `PASS_REGULAR_EXPRESSION` on the residual line, or have `validate_sqm.py` exit 0 within a documented band. | S-M |
| H5 | **`-dipole` throws away its result.** `executeDipole` computes `linear_vector` and `nonlinear_vector` and then prints only "Dipole optimization completed."; the line before is an AI truncation artefact: `// Output results... (truncated for brevity, keep existing logic)`. `OptimiseDipoleScaling.h` prints nothing (only commented couts). | `src/main.cpp:1396-1400`; `src/capabilities/optimiser/OptimiseDipoleScaling.h:79-160`. Only occurrence of such an artefact in `src/` (grep for "truncated for brevity" etc.). | S |
| H6 | **External-library builds are not built anywhere in CI.** All CI configurations pass `-DUSE_XTB=OFF -DUSE_TBLITE=OFF -DUSE_D3=OFF -DUSE_D4=OFF` (Ulysses/GFNFF default OFF); Vulkan is never built. The verified bugs Q-1..Q-7 in `docs/TECHNICAL_DEBT.md` sit exactly in code CI never compiles. | `scripts/ci/linux_container_build.sh:34`, `linux_host_gpu_build.sh:65`, `build_windows.bat:21`, `ci-variants.yml:71`. | M (CI minutes, Fortran deps) |

### 1.2 Medium severity

| # | Finding | Evidence | Effort |
|---|---|---|---|
| M1 | **37 capability TUs (27,926 lines) compiled twice.** `curcuma_cap_SRC` is a strict subset of `curcuma_core_SRC` (37 of 38 entries overlap, the 38th is a duplicate line). Only two tests and the external qurcuma link `curcuma_cap`; `curcuma` links only `curcuma_core`. `src/capabilities/CLAUDE.md:19` codifies "add the .cpp to `curcuma_core_SRC` and `curcuma_cap_SRC`". Same symbols in two static libs = latent ODR/link-order trap. | `CMakeLists.txt:492-616` vs `618-656`; `:1545`; `test_cases/CMakeLists.txt:165,589`; `docs/TOOL_API_WP.md:25`. Fix: `add_library(curcuma_cap INTERFACE)` + `target_link_libraries(curcuma_cap INTERFACE curcuma_core)`, keeping the name for qurcuma. | S-M |
| M2 | **`HELPERS` option is undeclared, `Helper` is declared and unused.** Five `if(HELPERS)` blocks; `option(Helper ...)` read nowhere. The helpers CLAUDE.md already lists this as a trap instead of fixing it. | `CMakeLists.txt:45-46` vs `:799,843,875,932,977`; `src/helpers/CLAUDE.md` Traps. | S |
| M3 | **`C17=OFF` is a dead configuration.** OFF sets `CXX_STANDARD 14` on 7 targets, but 14 files include `<filesystem>` and 11 use structured bindings unconditionally. `CMAKE_CXX_STANDARD 17` is also set at the top. | `CMakeLists.txt:6-8,42,770-778,804-808,849-853,881-885,981-985,1591-1595`. That C++14 fails to compile: inferred, **(unverified)**. | S |
| M4 | **Version embedding: unused global -D flags, dead `_UNIX` branch.** `add_definitions(-DGIT_COMMIT_HASH/BRANCH/DATE/TAG)` reach every TU (incl. nvcc) but no source uses the macros (`version.h.in` uses `@VAR@` substitution). `version.h.in` tests `#ifdef _UNIX`, but `_UNIX` is only a CMake variable (`CMakeLists.txt:321`), never `#cmakedefine`d, so the commit date is always "Forever young." and `_DEBUG` is never set by CMake on GCC ("RELEASE Mode" always). `version.h` is included by `global.h` (fan-in 47 headers), so every configure after a commit changes `version.h` and recompiles most of the tree. | `CMakeLists.txt:313-318`; `src/version.h.in:10,19,24`; `git grep -w _UNIX` (2 hits, both above); include-graph scan. | S |
| M5 | **Ulysses is fetched from an unpinned branch, even with `USE_ULYSSES=OFF`.** `GIT_TAG main`; with OFF it is still populated "for the include path", but no TU outside `USE_ULYSSES` includes Ulysses headers (`ulyssesinterface.cpp:25` guards the include; `src/pch/pch_ulysses.h` is included by nothing). Default build depends on gitlab.com and on whatever `main` is today. | `CMakeLists.txt:461-481,277`; `git grep pch_ulysses` (0 includes). That OFF builds without it: **(unverified)**. | S |
| M6 | **Two assertion tests are built but never run.** `rmsd_test` (expects RMSD 2.87214) and `reorder_test` (2.05127) have `add_executable` but no `add_test`. Also never run: `test_gfnff_grad_traj`, `test_gpu_numgrad`, 7 `diag_*`/`dump_*` tools (those are diagnostics, fine). | scan of all `add_executable` vs `add_test` in the four CMake files; `test_cases/rmsd/main.cpp:36-53`; `src/capabilities/rmsd/CLAUDE.md:39` notes it. Whether they still pass: **(unverified)**. | S |
| M7 | **Test registration split and labels missing.** 19 tests are registered in the root `CMakeLists.txt` (`:1604-1645`) for executables defined in `test_cases/CMakeLists.txt`, none with LABELS; 12 more in `test_cases/CMakeLists.txt` lack LABELS (`cg_potentials`, `external_potentials`, `geometry_tools`, `test_periodic_table`, `test_fast_exp`, `test_eeq_two_phase`, `test_gfnff_topology_reentrancy`, `test_polymerbuild_*`, `cif`, `measurement_capability`, `gfnff_react_fd_gradient`). `ctest -L` cannot select them. Also `target_sources(Continuous PRIVATE qmdff_par.h)` attaches a header to CTest's dashboard target (`:1600-1603`), purpose unclear. | Python scan of `add_test`/`set_tests_properties`. | S |
| M8 | **Disabled tests waiting for a decision.** `cli_curcumaopt_02_gfn2_single_point` (runs uff, duplicates 01), `cli_sqm_10_gfn2_optimization_smoke` (wrong log markers), `cli_cg_01_single_point` commented out (`cli/CMakeLists.txt:258`), `AAAbGlc_incremental` commented out (`CMakeLists.txt:1609`). Each has a reason comment; none has a TODO.md entry I could find by grep. | `test_cases/cli/CMakeLists.txt:89-91,327-329,258`. | S each |
| M9 | **GPU backends: three hand-maintained copies.** CUDA (23.1k lines `.cu/.cuh`), ROCm (21.4k lines, hand-ported HIP, not hipified), Vulkan (31.7k lines incl. 47 generated `.spv.inc`). CMake plugin boilerplate (mirror `curcuma_core` COMPILE_DEFINITIONS, PIC, output dir, `add_dependencies(curcuma ...)`) is written three times; ROCm hardcodes `/opt/rocm` 9 times. `scripts/gpu_parity_check.py` exists as mitigation. Committed SPIR-V `.spv.inc` are regenerated only by hand (`vulkan/shaders/compile_shaders.sh`), no check that they match the 47 `.comp` sources. | `CMakeLists.txt:1052-1464`; `wc -l` per backend. | L (kernels), S (a `curcuma_add_gpu_plugin()` function) |
| M10 | **Documentation volume and drift.** `docs/` holds 156 top-level `.md` (288 tracked files under `docs/`): 35 `*WP*`, 13 `*PLAN*`, 8 `*ROADMAP*`, 4 `*AUDIT*`. 12 have no inbound reference anywhere (`ANALYSIS_PARALLELIZATION`, `D4_HEAVY_ELEMENT_FIX`, `GFNFF_EEQ_METAL_CHARGE_FIX`, `GFNFF_METAL_BOND_FIX`, `GFN_XTB_DOCUMENTATION_AUDIT`, `GPU_OPTIMIZATION_ROADMAP`, `NATIVE_GBSA_IMPLEMENTATION`, `PARAMETER_CACHING`, `PR40_D4_RESIDUAL_FIX`, `SCATTERING_ANALYSIS_GUIDE`, `TBLITE_SOLVATION_IMPLEMENTATION`, `VERBOSITY_SYSTEM_ANALYSIS`). Applying `check_docs.py`'s backtick-path rule to non-archive `docs/*.md`: **240 paths in 73 of 181 files match no tracked file** (mostly deleted `gfn1.cpp`, `gfn2.cpp`, `forcefieldthread.*`, `gfnff_advanced.cpp`). Example: `docs/BUILD_SYSTEM_ROADMAP.md` (Nov 2025) says builds 1,3,4,5 are "BLOCKED" by `gfnff.cpp`/`forcefieldthread.h`, both deleted. | scan script (scratchpad `pathscan.py`, same regex as `check_docs.py`). | M (archive pass) |
| M11 | **Stale architecture docs inside `src/`** (not CLAUDE.md, so outside `check_docs.py`): `src/core/REFACTORING_ROADMAP.md` (221 l), `REFACTORING_TODO.md` (85), `UNIFIED_INTERFACE_TODO.md` (60), `src/capabilities/OPTIMIZATION_MODERNIZATION.md` (59, 0 inbound refs, describes an `OptimizerStrategy` class that does not exist), `energy_calculators/ENERGY_SYSTEM_OVERVIEW.md` (633) and `qm_methods/QM_ARCHITECTURE.md` (505). `ENERGY_SYSTEM_OVERVIEW.md:45,355` says `gfn2: TBLite -> Ulysses -> XTB`, contradicting root CLAUDE.md (native canonical since AP3); both docs cite `src/capabilities/curcumaopt.cpp` (deleted 2026-10-01) and `docs/PARAMETER_MIGRATION_GUIDE.md` (moved to `docs/archive/`). | `git grep -c 'class OptimizerStrategy' -- src` = 0. | S-M |
| M12 | **Root CLAUDE.md contradicts itself and the code.** (a) `:145` "None is ✅ TESTED" vs `:202` "gfnff ... ✅ **COMPLETE**"; (b) `:387` "Thread Safety: Zero overhead at Level 0 ✅" and `:365` "✅ FULLY IMPLEMENTED" vs open item "Verbosity cannot be scoped across CxxThreadPool workers (#3)" and TECHNICAL_DEBT Q-10/Q-49; (c) `:191` `MethodFactory::createMethod(...)` does not exist (it is `MethodFactory::create`, `energycalculator.cpp:216`); (d) `:225` QMInterface example shows `Calculation(bool gradient, bool verbose)` and pure-virtual `InitialiseMolecule`, real header has `Calculation(bool gradient = false)` (`abstract_interface.h:123`) and a non-pure default (`:68`); (e) `:198` "pm3: Native only (H, C, N, O supported)" vs `nddo_params.h:90` "PM3 parameter set (30 elements: H through Bi)"; (f) `:177` "COMPLETELY REFACTORED (January 2025)" is a dated status line the "Where Things Go" rules forbid (not caught: the check only matches headings); (g) Build section per-category counts "(6/6 passing)" are labelled historical but still printed. | line numbers from `grep -n`. | S |
| M13 | **AIChangelog.md is not an index.** Rule: one line per fact, long entries in `docs/changelog/`. 29 lines exceed 600 characters, one is 2,500 (`AIChangelog.md:40`). `TODO.md` has 3 lines > 600 characters, and two open entries sit after section 8 ("removed from the old list"), i.e. in the wrong section. `check_docs.py` checks neither file's line length. | `awk 'length>600'`; `TODO.md` tail. | S-M |

### 1.3 Low severity

| # | Finding | Evidence | Effort |
|---|---|---|---|
| L1 | `USE_D3` help text says "Compile cpp-D4". | `CMakeLists.txt:36-37` | S |
| L2 | `analysenciplot.cpp` listed twice in both source lists. | `CMakeLists.txt` core and cap lists | S |
| L3 | Header files in `curcuma_core_SRC` (`formats.h`, `geometry.h`, `general.h`, `forcefieldfunctions.h`, `orcainterface.h`, `orca_method.h`) - harmless, inconsistent. | same | S |
| L4 | `find_package(BLAS/LAPACK)` appears three times (`:86-87`, `:954-955`, `:1483-1484`); BLAS linked twice when `USE_BLAS=ON` (`:1489-1490` and `:1522-1528`). | | S |
| L5 | `-g` is in `CMAKE_CXX_FLAGS` for every build type (debug info in Release; only `curcuma`/`curcuma_helper` are stripped with `-s`, plugins and libs are not). `-Wsign-compare ... -Wno-sign-compare` cancel; `-Wno-narrowing` hides narrowing; `-Wimport` is obsolete. Warning count of a real build: **(unverified)**. | `CMakeLists.txt:347-362` | S |
| L6 | `CMAKE_MODULE_PATH` points at `cmake/`, which does not exist. `cmake_policy(SET CMP0169 OLD)` keeps the deprecated direct `FetchContent_Populate` (will break on a future CMake). | `CMakeLists.txt:249-251,259` | S / M |
| L7 | `target_precompile_headers(curcuma REUSE_FROM ...)` called twice (core, then cap); the second overwrites the first. | `CMakeLists.txt:765-766` | S |
| L8 | `km` (C Hungarian, `src/capabilities/c_code/`) comment says "nothing in the codebase includes its public interface"; false: `rmsd.cpp:32` and `rmsd/rmsd_assignment.cpp:27` include it and call `assign()`, but only on dead branches (`bool use_c = false;` at `rmsd.cpp:1320`; `MunkresAssignmentSolver` has no caller, already a TODO deletion candidate). The library is still built and linked. | `CMakeLists.txt:723-736` | S (with the TODO §7 decision) |
| L9 | `src/pch/pch_ulysses.h` included by nothing. 4 helper `.cpp` with `main()` and no target (already in TODO §7). | `git grep pch_ulysses` | S |
| L10 | Only include cycle among 209 headers: `src/tools/formats.h` <-> `src/tools/general.h` (both header-only, 692 and 528 lines). No `src/core` file includes `src/capabilities` (layering holds at include level, but the `curcuma_core` library itself contains all capabilities, see M1). | scratchpad `inc.py` (Tarjan SCC) | M |
| L11 | 54 `fmt::print`/`cout` lines with non-ASCII characters (rule "No UTF symbols in terminal output"), e.g. em dash in `ancopt_optimizer.cpp:502,507,513`; 788 raw `std::cout/cerr` in `src/*.cpp`. | `git grep -P` | S per file |
| L12 | 15 of 26 `src/` directories have no CLAUDE.md although the root rule says each source dir has one (e.g. `capabilities/handlers`, `capabilities/optimiser`, `core/solvation`, `ff_methods/cuda`, `ff_methods/rocm`, `qm_methods/vulkan`). Some are deliberately covered by the parent. | dir scan | S |

### 1.4 Repository hygiene (stray root items)

| Item | What it is | Referenced by | Proposal |
|---|---|---|---|
| `STIL_Funding.jpg` (928 KB) | funding logo | `README.md:249` (raw GitHub URL) | keep or move to `misc/`, update README URL; operator decision |
| `misc/` (6 files, 14 MB) | logos, `.desktop` | `curcuma_I.png` + `curcuma.desktop`: AppImage step `linux-gpu.yml:98-100`; `curcuma_II.png`: README logo; `curcuma_III.png` (4.6 MB) and `curcuma_WH.svg`: no reference except Inkscape export metadata inside the SVG | `curcuma_III.png` is a deletion candidate; the three PNGs are 26% of the 51 MB tracked tree |
| `patches/tblite/` | patch for the tblite reference dumper | `docs/SQM_VALIDATION.md`, `test_cases/sqm_reference/dump_tblite_reference.cpp` | keep (in use) |
| `plans/analytical_d3_gradient_plan.md` | one plan | `docs/SQM_VALIDATION_ROADMAP.md` | move to `docs/` (or archive); one-file top-level dir |
| `tools/gfnff_ref_generator.py` | GFN-FF reference generator | `docs/GFNFF_STATUS.md` | move to `scripts/` (48 scripts already there); update the one link |
| `python_scripts/distance_histogram.py` | standalone script | nothing (`git grep python_scripts` = 0) | delete or move to `scripts/`; operator decision |
| `test_mol/OH.xyz` | 2 lines, no atom-count/comment header (not valid XYZ) | nothing (`git grep 'test_mol/'` = 0; other "OH.xyz" mentions are archive docs) | delete; library already has `small/hydroxyl.xyz` |
| structure files outside `test_cases/structures/` | `docs/validation/{benzene,cyclopropane}.xyz` (cited by `docs/theory/PHASE2_TOPOLOGY_DETECTION.md`), `src/molecules/{pdmaema,peo,ppo}.xyz` (no reference found by grep; possibly user-facing polymer fragments, **unverified**), `test_cases/cli/simplemd/15_gfnff_react_water_noevent/input.xyz` (added Sep 2026 after the library migration, not byte-identical to any library entry), `cli/errors/02_invalid_xyz_format/invalid.xyz` (deliberate, listed in manifest) | | import the cli/simplemd/15 water hexamer into the library (`structlib.py add`) and add `structures.txt`; decide on `docs/validation` and `src/molecules` |

### 1.5 Architecture seams (build-visible)

- Only `EnergyCalculator` calls `MethodFactory::create` (`energycalculator.cpp:216`); no capability constructs a method class directly (grep for `new GFNFF`, `XTBInterface`, `TBLiteInterface`, ... in `src/capabilities`, `src/main.cpp`: 0). Capabilities reach into the factory only for `methodParameterScopes()` (`confsearch.cpp:1107`, `hessian.cpp:229`, `simplemd.cpp:887`, `main.cpp:721,1325`). Seam is clean.
- Exceptions: `main.cpp:1409` uses `OrcaInterface` directly for `-orca` (pass-through of an existing input; acceptable), `qmdfffit.cpp` includes `forcefieldgenerator.h`/`qmdff_par.h`/`uff_par.h`, `polymerbuild.cpp` includes `cn_calculator.h`.
- `src/main.cpp`: 3,161 lines, 37 `execute*` handlers; the CLI dispatcher is the largest seam and where H5 hid.
- Library layering is nominal: `curcuma_core` contains every capability (M1), so "core" cannot be linked without capabilities; GPU plugins resolve core symbols against the `-rdynamic` executable (documented in CMake).
- Still-dead `MethodFactory` lists (TECHNICAL_DEBT D-7 partly open): `m_ff_methods`, `m_xtb_methods`, `m_d3_methods`, `m_d4_methods` have 2 occurrences each (declaration + definition), no reader.

### 1.6 `scripts/check_docs.py` output

```
$ python3 scripts/check_docs.py
check_docs: 0 error(s), 0 warning(s)          (exit 0)

$ python3 scripts/check_docs.py --all-docs
check_docs: 15 error(s), 0 warning(s)
```
The 15 `--all-docs` errors: 7 are real broken links, 8 are false positives from math notation parsed as links
(`[...](μ,ν)`, `[...](1,0)` in `docs/AP6b_multipole_F_discrepancy_WP.md:127,163` and `docs/GFN2_MULTIPOLE_IMPLEMENTATION.md:197,198`).
Real ones: `docs/CONFIG_MANAGER.md:412`, `docs/PARAMETER_SYSTEM.md:467`, `docs/SOLVATION.md:263` -> `PARAMETER_MIGRATION_GUIDE.md` (now in `docs/archive/`);
`docs/archive/PARAMETER_MIGRATION_GUIDE.md:479` -> `PARAMETER_SYSTEM.md` (needs `../`); `docs/GFNFF_IMPLEMENTATION_HUB.md:16,17` -> `theory/GFNFF_{TORSION,INVERSION}_THEORY.md` (missing);
`docs/theory/GFNFF_COMPLETE_GUIDE.md:395` -> `GFNFF_IMPLEMENTATION_HUB.md` (needs `../`); `docs/wp4/WP-EEQ-many-fragment-diagnosis.md:9` -> `WP-EEQ-P2-profile-plan.md` (missing).
Gaps in the checker itself: it does not check AIChangelog.md/TODO.md line length, code blocks in CLAUDE.md, markdown outside CLAUDE/README/TODO, or bold dated status lines.
`python3 scripts/structlib.py check`: 80 structures, 0 errors, 6 warnings (generic `.alt1` names, `dihydrogen.r0.47` distance).

---

## 2. Safe first steps (small, low-risk commits)

Each is independent. "Verify" = what to run before committing (in `release/`, per root CLAUDE.md).

1. **CI: make smoke tests fail on failure** - `.github/workflows/ci-variants.yml`, `linux-gpu.yml`: add `shell: bash` (gives `-eo pipefail`) or `set -o pipefail` to the steps with `| tail`. Verify: push to a `fix/ci-smoke-pipefail` branch; temporarily feed a missing input file and see the job go red.
2. **CI: docs and structure check job** - new job in `ccpp.yml` (ubuntu, python3 only): `python3 scripts/check_docs.py && python3 scripts/structlib.py check`. Verify: both exit 0 today (shown above).
3. **Declare `HELPERS`, drop `Helper`** - `CMakeLists.txt:45-46`: `option(HELPERS "Build interface helper programs" OFF)`; update the trap line in `src/helpers/CLAUDE.md`. Verify: `cmake -LH .. | grep HELPERS`; `make` unchanged with default OFF.
4. **Fix `USE_D3` description, remove duplicate `analysenciplot.cpp` lines** - `CMakeLists.txt:36-37` and both source lists. Verify: `make -j4` exit status 0, `ctest -N` count unchanged.
5. **Register `rmsd_test` and `reorder_test`** - `test_cases/CMakeLists.txt`: `add_test(NAME rmsd_basic COMMAND rmsd_test WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})` (inputs `input_aa.xyz`/`input_ab.xyz` are staged there at `:3-5`), `LABELS "rmsd;unit"`. Verify: `ctest -R 'rmsd_basic|reorder_basic' --output-on-failure`; if the hard-coded 2.87214/2.05127 no longer match, report rather than change the numbers.
6. **Add LABELS to the 31 unlabelled tests** and move the 19 root-registered `add_test` lines (`CMakeLists.txt:1604-1645`) into `test_cases/CMakeLists.txt` next to their executables. Verify: `ctest -N` before/after gives the same names and count; `ctest -L unit -N` lists them.
7. **Remove unused `-DGIT_*` definitions and fix `_UNIX`** - delete `add_definitions(-DGIT_...)` at `CMakeLists.txt:313-318` (keep `configure_file` for `version.h`); in `src/version.h.in` replace `#ifdef _UNIX` with `#ifdef __unix__` (or add `#cmakedefine _UNIX` to `global_config.h.in`). Verify: `git grep -w -E 'GIT_COMMIT_HASH|GIT_BRANCH|GIT_COMMIT_DATE|GIT_COMMIT_TAG' -- src ':!src/version.h.in'` = 0 hits; build; curcuma's banner now shows the commit date.
8. **Drop the dead C++14 branch** - replace the seven `if(C17) ... else() ... CXX_STANDARD 14` blocks with `CXX_STANDARD 17`; keep `#cmakedefine C17` as always-on or remove the 8 `#ifdef C17` in `curcumamethod.cpp`, `simplemd.cpp`, `main.cpp` in a second commit. Verify: full build, `ctest`.
9. **Pin Ulysses** - `CMakeLists.txt:465`: `GIT_TAG <commit sha>` instead of `main`. Verify: fresh build dir configures; `USE_ULYSSES=ON` build if the operator has it.
10. **Delete unreferenced files** - `test_mol/OH.xyz`, `src/pch/pch_ulysses.h`; ask the operator about `python_scripts/distance_histogram.py`, `misc/curcuma_III.png`. Verify: `git grep` for each basename returns nothing outside archive docs; build.
11. **Fix the 7 real broken links** (list in 1.6) and teach `check_docs.py` to skip link targets containing `,` or non-ASCII. Verify: `python3 scripts/check_docs.py --all-docs` -> 0 errors.
12. **Root CLAUDE.md corrections** (M12 a-g): `createMethod` -> `create`, QMInterface snippet from the real header, PM3 element range, drop "✅ COMPLETE"/"Thread Safety ✅"/"COMPLETELY REFACTORED (January 2025)". Verify: `python3 scripts/check_docs.py`; root stays under 500 lines.
13. **Archive stale `src/` architecture docs** - `git mv` `src/core/REFACTORING_ROADMAP.md`, `REFACTORING_TODO.md`, `UNIFIED_INTERFACE_TODO.md`, `src/capabilities/OPTIMIZATION_MODERNIZATION.md` to `docs/archive/`; fix `ENERGY_SYSTEM_OVERVIEW.md:45,355` priority line. Verify: `git grep` for each basename and update the 1-2 inbound links; `check_docs.py`.
14. **Restore `-dipole` output** (behaviour change, small): print `linear_vector`/`nonlinear_vector` via `CurcumaLogger::result` and/or write them through `BMTUtils::outputPath()`; remove the truncation comment. Needs operator choice of format. Verify: run `curcuma -dipole` on a small conformer file with and without `-no_bmt`.

Larger follow-ups (not "first steps"): collapse `curcuma_cap` into an INTERFACE target (M1, check qurcuma first); make BLAS/LAPACKE linking follow `USE_BLAS` only, or at least print the decision in the run banner (H3; changes numerics on machines that have liblapacke, so needs a regression run per `docs/REGRESSION_CHECK.md`); replace `WILL_FAIL` with a residual-band check (H4); a CI job with `USE_TBLITE=ON`/`USE_D3/D4=ON` (H6); a `curcuma_add_gpu_plugin()` CMake function (M9); one archive pass over `docs/` WP/PLAN files whose work is done (M10).
