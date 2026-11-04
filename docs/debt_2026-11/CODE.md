# Code-level technical debt in `src/` (curcuma, HEAD c3f570f, 2026-11-04)

Read-only audit; no repository file was changed. AI-generated, not human-reviewed.

**Method and limits**
- Scope: 413 source files under `src/` (`.cpp/.h/.hpp/.cu/.cuh/.hip/.hiph/.inc/.comp`), 251 709 lines.
- No build tree exists in the checkout (`release/` holds only `main.cpp`, no `generated/`), so **nothing below is
  compile-verified**. Every "dead" claim is grep-verified: the identifier occurs only at its declaration and
  definition across `src/`, `test_cases/`, `scripts/` and all `CMakeLists.txt`, with comments stripped. Macro
  token-pasting and calls by name from outside the repo can defeat this; such cases are flagged.
- Helper scripts used (scratchpad): `orphans.sh` (files referenced nowhere), `dead.py` (member functions with
  <= 2 tokens), `decl.py` (declared, never defined). Raw output: `dead_raw.txt`, `dead_filtered.txt`,
  `blockcode.txt`, `unicode.txt`, `env.txt` in the same scratchpad directory.
- `docs/TECHNICAL_DEBT.md` is stale for `gfnff_method.cpp`: it says 11 738 lines (now 13 060); its line numbers are
  off by about 1 300; F-E3 (`m_forcefield` forks) is resolved (0 hits of `m_forcefield` in `gfnff_method.cpp`); F-Q3
  (`gfnff_advanced.cpp`) and D-4 (`energycalculator_enums.h`) files no longer exist. Still true: F-Q1, F-Q2, F-Q5,
  F-E6, F-E7, F-E9, F-W1, X-L2.

---

## 1. Findings ranked by severity and effort

Severity: H = wrong output possible / mandated rule broken in a hot path; M = maintainability, divergence, build
cost; L = cosmetic or dead code. Effort: S < 1 h, M = a day, L = several days.

### H / S-M

**H1. Unicode in terminal output (project forbids it): 80 lines with box drawing, arrows, check marks or emoji;
229 output lines with any non-ASCII character in 36 files.**
- Printed on every run at program end: `src/core/curcuma_logger.cpp:230-232` (citation block, `═` x63).
- GFN-FF energy report at verbosity 2: `src/core/energy_calculators/ff_methods/gfnff_method.cpp:148,183,230,278`
  (`─`/`═`, 13 lines in the file).
- Check marks and emoji: `src/capabilities/analysis.cpp:1049` (`✓`, 34 lines in file),
  `ff_methods/forcefield_method.cpp:545` (`✅`), `d3param_generator.cpp:426` (`└─`), `eeq_solver.cpp:683` (`→`),
  `simplemd.cpp:686` (`→`), `rmsd.cpp:404` (`↔`), `casino.cpp:912` (`→`), `scnp_parser.cpp:470`.
- Characters by count: `─` 197, `═` 195, `Å` 54, `•` 33, `σ` 18, `—` 17, `→` 15, `⁻¹²³` 36, `°` 7.
- Verify after fix: the Python scan in `unicode.txt`'s header (non-ASCII inside string literals on output lines)
  returns 0 outside `src/helpers/`.

**H2. File output that bypasses `outputPath()`/`BMTUtils::outputPath()` (mandated).**
- `ff_methods/gfnff_method.cpp:11099` `std::ofstream diag_file("gfnff_diag_charges.json")` (verbosity >= 3; was
  F-Q5, still open).
- `dispersion/d4param_generator.cpp:635` `"gfnff_diag_d4.json"` (verbosity >= 3).
- `ff_methods/eeq_solver.cpp:4529-4530` `charges_dump_N{}.json` (config `dump_charges`).
- `capabilities/confscan.cpp:1229` `confscan_saved_<n>.json` (only with `WriteMoreInfo`).
- Intentional CWD copies, check they are documented: `main.cpp:1724` (`rmsd_out_cwd`), `confsearch.cpp:1611`,
  `main.cpp:134` (`stop` signal file).
- These also make file I/O depend on the print level (root CLAUDE.md "no numerical path may depend on the print
  level" is about numerics, but these writes collide between parallel runs in one CWD).

**H3. Two GFN-FF EEQ code paths, one unregistered.** `gfnff_method.cpp:10715-10718` reads JSON key
`use_two_phase_eeq` (not a PARAM, default true); the `false` branch at `:11111` calls the legacy
`calculateEEQCharges` (`:7785`, 244 lines). Untested divergent path reachable by a typo-free JSON key. Was F-Q1.
Effort M: delete branch + function, or register as PARAM with a test.

**H4. Function-local mutable `static` counters in compute code (thread-unsafe under `CxxThreadPool`, state leaks
across molecules).** `gfnff_method.cpp:1240` (`call_no`), `:5383` (`debug_once`), `:5402` (`debug_count`), `:5800`
(`bond_count`); `xtb_scf.cpp:384,410`, `xtb_native.cpp:568` (`warned*`), `rocm/gfnff_rocm.hip:5711`
(`uploaded_dev/val`). 71 function-local mutable statics overall. The GFN-FF ones gate debug printing only (data race,
no numeric effect); `gfnff_rocm.hip:5711` caches a device upload key and should be checked for multi-GPU
correctness. Was F-E9.

### M / S (build and configuration)

**M1. `curcuma_cap` is a second compile of 38 capability TUs that `curcuma_core` already compiles.**
`CMakeLists.txt:618-660` vs `:493-617`: every file in `curcuma_cap_SRC` is also in `curcuma_core_SRC`
(`comm -12` of the two lists = 38 files, including `simplemd.cpp` 5158 lines, `confscan.cpp`, `analysis.cpp`).
`curcuma` links only `curcuma_core` (`:1545`); `curcuma_cap` is linked only by `test_external_potentials` and
`test_measurement` (`test_cases/CMakeLists.txt:165,589`), which link `curcuma_core` too. Also
`target_precompile_headers(curcuma REUSE_FROM ...)` is given twice (`:765-766`).

**M2. `src/capabilities/analysenciplot.cpp` listed twice in `curcuma_core_SRC`** (`CMakeLists.txt:504,511`).

**M3. Unused compile definitions force full rebuilds.** `CMakeLists.txt:313,317,318` add
`-DGIT_COMMIT_HASH/-DGIT_COMMIT_DATE/-DGIT_COMMIT_TAG` to every TU; no source uses these macros (`grep -rnw
GIT_COMMIT_HASH src` hits only `version.h.in`'s `@...@` placeholders). In addition `src/core/global.h` includes
`src/version.h` (git hash) and is included by 85 files, while only `tools/info.h:34-35` and `molecule.cpp:1510`
read the strings.

**M4. `-g` hardcoded into `CMAKE_CXX_FLAGS` for GCC in all build types** (`CMakeLists.txt:363`), plus
`-Wsign-compare ... -Wno-sign-compare` in the same line. Release binaries carry debug info; check intent.

**M5. PCH configuration is partly fiction.** `src/pch/pch_external.h` contains only `CxxThreadPool.hpp` and
empty `#ifdef USE_XTB/USE_TBLITE` stubs; `CMakeLists.txt:713-719` print "Ulysses PCH module enabled" and "DFT
dispersion correction PCH modules enabled" although neither is precompiled; `src/pch/pch_ulysses.h` is included
nowhere. `pch_core.h` pulls `uff_par.h`/`qmdff_par.h` (namespace-scope `const json` tables, 16 statics) into every
core TU; `QMDFFParameterJson` has no reader.

**M6. External dependencies not pinned.** `CMakeLists.txt:466` ulysses `GIT_TAG main`; `:793,833` mctc-lib
`GIT_TAG "HEAD"`; nlohmann json v3.9.1 (2020) by `file(DOWNLOAD)` at `:437`. Everything else is pinned by commit.
`external/` is not vendored and has no submodules (`.gitmodules` is a comment; FetchContent populates `external/`,
gitignored lines 234-246).

**M7. `#pragma message` in headers: 13 sites, 3 in widely included headers** - `tools/general.h:385,398`
(included by 40-45 files), `tools/formats.h:51`. Each prints once per TU per build. Others:
`simplemd.cpp:4845`, `rmsd.cpp:487`, `forcefieldgenerator.cpp:546,656,815,823`, `forcefield.cpp:88,98`,
`forcefieldfunctions.h:279` (inside a dead function, see S8).

**M8. Misleading CMake comment about `km`.** `CMakeLists.txt:722-727` says nothing calls `assign()`; in fact
`rmsd.cpp:1330` calls it inside `if (use_c)` with `bool use_c = false;` (`:1320`), and
`rmsd/rmsd_assignment.cpp:67` calls it inside the never-instantiated `MunkresAssignmentSolver`. On MSVC `km` is not
built, so the reference only links because the optimiser drops the constant-false branch.

### M / M-L (structure and duplication)

**M9. God functions** (start-end by column-0 brace, verified):

| Lines | Function | Location |
|---|---|---|
| 2064 | `XTB::Calculation` (SCF loop, modes, GPU dispatch, extrapolation, diagnostics) | `qm_methods/xtb_native.cpp:496-2559` |
| 1368 | `GFNFFGpuMethodImpl::calculateEnergy` | `qm_methods/gfnff_gpu_method_impl.h:982` |
| 1163 | `GFNFF::getGFNFFBondParameters` | `ff_methods/gfnff_method.cpp:4671-5833` |
| 1070 | `GFNFF::calculateTopologyInfoOnce` | `ff_methods/gfnff_method.cpp:10496-11565` |
| 945 | `EEQSolver::calculateFinalCharges` | `ff_methods/eeq_solver.cpp:3690` |
| 879 | `ConfSearch::start` | `capabilities/confsearch.cpp:87` |
| 834 | `PolymerBuild::buildSubchain` | `capabilities/polymerbuild.cpp:2181` |
| 811 | `GFNFF::detectHydrogenBondsNative` | `ff_methods/gfnff_method.cpp:9129` |
| 789/788 | `getGFNFFTorsionParameters` / `getGFNFFAngleParameters` | `gfnff_torsions.cpp:511`, `gfnff_method.cpp:5835` |
| 779 | `XTB::calculateGradient` | `qm_methods/xtb_gradient.cpp:51` |
| 582/572 | `launchChargeDependentAndFinish` CUDA / ROCm twins | `cuda/ff_workspace_gpu.cu:3531`, `rocm/gfnff_rocm.hip:9522` |
| 450 | `CLI2Json` | `main.cpp:675` |

**M10. Top 15 files by line count and what they mix**

| Lines | File | Mixes |
|---|---|---|
| 32 838 | `ff_methods/d3_reference_c6.cpp` (3.1 MB) | generated D3 C6 table only (data, compile cost) |
| 13 060 | `ff_methods/gfnff_method.cpp` | topology, all parameter generation, EEQ dispatch, HB/XB detection, CN derivatives, JSON (de)serialisation, numgrad, verbosity-3 dumps, dead wrappers (S5) |
| 10 675 | `ff_methods/rocm/gfnff_rocm.hip` | ROCm kernels + FFWorkspaceHip + EEQ glue; 59 % of the CUDA lines appear verbatim (M11) |
| 8 186 | `qm_methods/cuda/xtb_gpu_context.cu` | SCF, eigensolver variants, f32 probes, integrals, memory estimates; 59 `catch (...)` |
| 6 072 | `qm_methods/rocm/xtb_hip_context.hip` | ROCm twin of the above (37 % verbatim) |
| 5 486 | `core/solvation/tblite_solvation_params.h` | data tables in a header (1 includer) |
| 5 249 | `ff_methods/eeq_solver.cpp` | 6 solver back-ends, topology/final charges, legacy "Phase 1" (`:3118` marked DEPRECATED), dumps, 63 `std::cout` |
| 5 158 | `capabilities/simplemd.cpp` | integrators, thermostats, MTD/RMSD-MTD (`EvaluateBias` 449 lines), walls, restart I/O, printing, empty Rattle stubs |
| 4 550 | `ff_methods/cuda/gfnff_kernels.cu` | all GFN-FF CUDA kernels |
| 4 479 | `ff_methods/cuda/ff_workspace_gpu.cu` | buffers, launches, host orchestration |
| 3 937 | `ff_methods/cuda/eeq_solver_gpu.cu` | GPU EEQ solvers |
| 3 816 | `qm_methods/xtb_native.cpp` | one 2064-line `Calculation()` plus setup |
| 3 311 | `capabilities/polymerbuild.cpp` | building, CN, docking, optimisation glue |
| 3 161 | `main.cpp` | CLI parsing, 37 `execute*` handlers, help registry, 238 `std::cout/cerr`, 4 `exit()` |
| 3 153 | `ff_methods/gfnff.h` | the whole GFNFF class, nested structs, GPU seams |

**M11. CUDA/ROCm copy-paste.** Normalised verbatim line overlap (lines >= 25 chars, cuda/hip tokens neutralised):
GFN-FF `gfnff_kernels.cu + ff_workspace_gpu.cu + eeq_solver_gpu.cu` vs `gfnff_rocm.hip + eeq_solver_hip.hiph`:
3 771 of 6 413 lines (59 %). xTB `xtb_gpu_context.cu + xtb_gpu_integrals_device.cuh` vs `xtb_hip_context.hip +
xtb_hip_integrals.hiph`: 1 851 of 5 040 (37 %). `xtb_gpu_method.cpp` vs `xtb_hip_method.cpp`: 19 %. The GFN-FF
method wrapper is already shared (`gfnff_gpu_method_impl.h`, both `gfnff_{gpu,hip}_method.cpp` are about 50 lines)
and is the model to follow. Same `use_mixed_precision = false; ///< TODO: Restore` in both twins
(`cuda/ff_workspace_gpu.cu:782`, `rocm/gfnff_rocm.hip:6623`).

**M12. Duplicated numeric routines on the CPU side.**
- SCC pack/unpack lambdas in three places (`xtb_native.cpp:1019,1034`, `xtb_dc.cpp:255,265`,
  `xtb_sparse.cpp:208,218`) next to member `XTB::packSccState/unpackSccState` (`xtb_native.cpp:413,433`) - X-L2.
- Coordination numbers: `CNCalculator` (`cn_calculator.cpp:44,120,256,357`) plus own routines in
  `EEQSolver::calculateCoordinationNumbers` (`eeq_solver.cpp:5060`), `XTB::computeCoordinationNumbers`
  (`xtb_h0.cpp:32`), `computeD4CovalentCN` (`d4_ncoord.cpp:104`), `polymerbuild.cpp:1782`.
- Numerical gradients: 5 implementations (`energycalculator.cpp:555`, `computational_method.h:172`,
  `gfnff_method.cpp:12940` and `:12991`, `forcefield.cpp:575`); two of them are dead (S6).
- Covalent-radius tables: at least 7 (`elements.h:460,685`, `gfnff_par.h:296,351`, `cn_calculator.cpp:27`,
  `d4param_generator.cpp:704`, dead `hbonds.h:77`). Values differ by method on purpose; do not merge blindly.

**M13. Two abstract bases.** 8 classes implement `QMInterface` (incl. `XTB`, `QMDriver` -> `EHT`, `NDDO`) and 15
implement `ComputationalMethod`; every external backend has a `*_method.cpp` adapter plus a `*interface.cpp`
(D-6, unchanged). `forcefield.h:66` / `forcefield.cpp:65` "TODO: merge QMInterface and ForceField" is the same item.

**M14. Console output bypassing `CurcumaLogger`.** `std::cout/cerr`: 809 non-comment uses; library code
`src/core` 249 (`eeq_solver.cpp` 63, `parameter_registry.cpp` 49, `ParallelEigenSolver.hpp` 37,
`config_manager.cpp` 32, `molecule.cpp` 16), `src/capabilities` 175, `main.cpp` 238, `helpers` 149.
`printf` 55, `fmt::print` 306 (allowed by `src/CLAUDE.md`). Several `std::cout` sites are not verbosity-gated
(e.g. `ParallelEigenSolver.hpp:117,143,454,471` timing prints) - check against Level 0 silence.

**M15. 107 distinct `getenv("...")` switches** (`xtb_native.cpp` 24, `xtb_gpu_context.cu` 15, `gfnff_rocm.hip` 13,
`xtb_hip_context.hip` 10). 10 are mentioned nowhere in `docs/` or `README.md`: `CURCUMA_VK_TEST_SECULAR`,
`CURCUMA_SCF_JAC_LINTEST`, `CURCUMA_SCF_JAC_LINSCALE`, `CURCUMA_PLAIN`, `CURCUMA_NO_COLOR`,
`CURCUMA_LOBPCG_DEBUG`, `CURCUMA_ADAPTIVE_DEBUG`, `ORCA_ROOT`, `ORCA_PATH`, `FORCE_COLOR`.

**M16. `exit()` inside library code.** `capabilities/confscan.cpp:339` (inconsistent step lengths),
`capabilities/rmsd.cpp:399` (atom-list length mismatch), `capabilities/curcumamethod.cpp:267` (`exit(0)` on
`-help` from any capability), dead `core/hbonds.h:224`. `main.cpp` has 4 more (acceptable at top level). 132
`throw`, 95 `catch (...)` (59 in `xtb_gpu_context.cu`, mostly allocation guards that return false; spot-checked
`:461,481,643,744,881,3576`), so error handling mixes three styles.

**M17. Static JSON defaults still in headers despite the Parameter Registry rule.** `casino.h:99` (`CasinoJson`,
including restart options whose implementation is never called, see S10), `confstat.h:35`, `tda_engine.h:34`,
`analysenciplot.h:24`, `forcefieldgenerator.h:38-99` (six non-const `static json` objects, one copy per TU).

### L (dead code, stale comments)

- **L1. Stale TODO comments: 54** (no FIXME/HACK/XXX). Clusters: GFN-FF 11 in `gfnff_method.cpp` (several obsolete:
  `:4492` "Implement GFN-FF specific topology calculation" on a function that only returns true; `:11591,11598`
  inside a dead overload; `:12261`), "make pure virtual" x5 on `start()` (`curcumamethod.h:71,82`,
  `confscan.h:267`, `analysenciplot.h:40`, `curcumaopt.h:262`; `CurcumaMethod::start` already is pure virtual),
  Ulysses solvation tables x3 (`solvation_parameters.h:155,186,203`), unimplemented output handlers
  (`handlers/topology_handler.cpp:55`, `handlers/rdf_handler.cpp:54`), `native_xtb_method.h:148`
  (`getDipole()` returns zero), `qm_methods/gfnff_method.h:35-36` (asks whether the file is obsolete; it is the
  live adapter - F-W1 naming clash with `ff_methods/gfnff_method.cpp`).
- **L2. Commented-out code**: 27 `/* ... */` blocks with >= 4 code lines (353 lines; largest
  `hbonds.h:860` 41, `forcefieldgenerator.cpp:553` 38, `molecule.cpp:172` 33, `pseudoff.h:62` 23,
  `formats.h:672` 19, `forcefieldgenerator.cpp:730` 17); about 216 `//`-commented code lines (top: `simplemd.cpp`
  27, `molecule.cpp` 20, `qmdfffit.cpp` 15, `rmsdtraj.cpp` 14, `confscan.cpp` 14, `topology.h` 12). Full list:
  `blockcode.txt`.
- **L3. "DEPRECATED" markers with live code**: `eeq_solver.cpp:3118` (Phase 1 wrapper still called from
  `gfnff_method.cpp:10890-10945`), `getDeviceCoordsPtr()` (S7), `gfnff_method.cpp:12489,12811` (comment-only
  tombstones after removed code - delete the comments with S5).
- **L4. No `#if 0`, no `INACTIVE (date)` markers, no `.backup/.orig/.bak/.rej` files** (tracked or untracked).
  `DEBUG_CONFIG_MANAGER` (29 `#ifdef` blocks in `config_manager.cpp`) and `DEBUG_PARAMETER_REGISTRY` are never
  defined by CMake - debug-only code that never compiles in CI.
- **L5. Comment-only remnants of removed engines**: `ForceFieldThread` 40 mentions in 11 files, all comments
  (e.g. `ff_workspace_gfnff.cpp` 14, `ff_workspace_uff.cpp` 8, `d4param_generator.cpp:561`), `SwitchMethod` 2
  (`method_factory.h:37,140`). No code left.
- **L6. Raw owning `new`: 116 sites** (`confsearch.cpp` 13 incl. `new Molecule` into `m_in_stack`, `ulysses.cpp`
  12, `main.cpp` 10 for capability objects, `rmsdtraj.cpp` 8), 70 `delete`. Thread objects handed to
  `CxxThreadPool` are owned by the pool; the `Molecule*` stacks and `main.cpp` capability objects are not.
- **L7. Test structures inside `src/`**: `src/molecules/{pdmaema,peo,ppo}.xyz`; root CLAUDE.md says no structure
  file outside `test_cases/structures/`. Only reader: orphan `helpers/gfnff_term_validator.cpp:51`.
- **L8. Doc claims checked against code and found wrong**: `ff_methods/CLAUDE.md:55` says `uff-d3` uses
  `ForceFieldGenerator::GenerateUFFD3Parameters()` - that function has no caller (S9); `uff-d3` sets
  `m_parameter["d3"]=1` in `forcefieldgenerator.cpp:253-256`. `core/CLAUDE.md:12` lists `parameter_validation.h`
  and `:17` lists `hbonds.h`, `pseudoff.*` as owned infrastructure; all three are unused (S1). `tools/CLAUDE.md:28,34`
  lists `trajectory_helpers.h` (unused duplicate, S1).

---

## 2. Safe deletions and cleanups (small commits, no behaviour change)

General verification for every item: (a) the grep given returns only the lines being removed; (b) full build in
`release/` (`make -j4`, check the exit status, not a grep for "error"); (c) `ctest --output-on-failure` gives the
same pass/fail set as before (TODO.md records three known failures). For PARAM-bearing headers also diff
`release/generated/parameter_registry.h` before/after `make GenerateParams`. None of this was run here (no build tree).

### Files with no includer and no CMake entry

**S1. Delete orphan headers / TUs** (each: `grep -rn "<basename>" src test_cases CMakeLists.txt scripts` shows
only the file itself and CLAUDE.md/doc lines):
- `src/core/hbonds.h` (923 lines; contains a commented `int main()`, `exit(1)`, `printf`, an own covalent-radius
  table; class `H4Correction` referenced nowhere else)
- `src/core/parameter_validation.h` (153 lines; `validate_*` functions have no caller) - also fix
  `src/core/CLAUDE.md:12`
- `src/capabilities/optimization_validator.h` (149 lines; `OptimizationValidator` only in itself)
- `src/tools/trajectory_helpers.h` (48 lines; duplicates `TrajectoryWriter::convertToTrajectoryJSON/
  createTrajectoryJSON`, `trajectory_writer.h:88,95`) - also fix `src/tools/CLAUDE.md:28,34`
- `src/pch/pch_ulysses.h` (27 lines)
- `src/core/pseudoff.cpp` (37 lines, not in CMake) together with the declaration
  `static double LennardJones(...)` at `src/core/pseudoff.h:32-33` (`grep -rn "LennardJones(" src` -> only these).
  `pseudoff.h` itself stays (inline `LJPotential`, `DistancePenalty` are used by `LevMar*.h`, `polymerbuild.cpp`).
- `src/core/energy_calculators/qm_methods/def2-SV(P).dat` (131 KB; `grep -rl "def2-SV" .` -> none)
- `src/core/energy_calculators/qm_methods/basissetparser.hpp` (652 lines; namespace `BasisSetParser`, `ElementBasis`,
  `BasisShell`, `exampleUsage` used nowhere) plus its `#include` at `qm_methods/qm_driver.cpp:35`. Needs a compile.
- Already on TODO.md 7: `src/helpers/{imagewrite,storage_bench,polymer_topo,gfnff_term_validator}.cpp` (1 708 lines,
  standalone `main()`, no CMake target, re-verified: zero references). `gfnff_term_validator.cpp` is the only reader
  of `src/molecules/` (L7).
- Unreferenced Python next to them: `src/helpers/{gfnff_reference_validator,gfnff_term_validator,validate_ch3oh}.py`
  (only mentioned in `src/helpers/CLAUDE.md` and an archive note) - operator decision, not code.

### Dead functions inside live files (grep: identifier occurs only at declaration + definition)

**S2. `src/capabilities/curcumaopt.h`: keep only the PARAM block.** Delete `OptJsonPrivate` (`:48-56`),
`getCurcumaOptJson`/`#define CurcumaOptJson` (`:39-45`), `LBFGSInterface`, `SPThread`, `OptThread`, `OptMThread`,
`CurcumaOpt` (`:114-315`) and the now-unneeded includes. `grep -rnw "CurcumaOptJson\|getCurcumaOptJson\|
LBFGSInterface\|OptMThread\|OptJsonPrivate" src` -> only this file; the header has no includer. Note `OptThread`
also exists in `confsearch.h:46` (ODR hazard if both were ever included). Verify: generated registry diff for
module `opt` is empty.

**S3. `forcefieldfunctions.h` dead helpers** (`src/core/energy_calculators/ff_methods/forcefieldfunctions.h`):
`UFF::Dihedral` (`:95-110`), `UFF::Torsion` (`:112-175`), `UFF::OutOfPlane` (`:177-235`, with the incomplete
derivative "TODO" at `:225-232`), whole `namespace QMDFF` (`:238-292`, incl. `#pragma message("check")`).
`grep -rnE "(UFF|QMDFF)::(Dihedral|Torsion|OutOfPlane|LJStretchEnergy|AngleDamping|AngleBend|LinearAngleBend)" src`
-> none; no `using namespace UFF/QMDFF` anywhere. Keep `BondStretching`, `AngleBending`, `NormalVector` (used in
`ff_workspace_uff.cpp`, `ff_workspace_gfnff.cpp`).

**S4. `rmsd.cpp` constant-false branch.** Remove `bool use_c = false;` and the `if (use_c) {...}` block
(`src/capabilities/rmsd.cpp:1320-1335`, raw `new[]`/`delete[]` included). Then correct the CMake comment
(`CMakeLists.txt:722-727`). Removing `km`/`c_code/` and `extern "C" #include interface.h` (`rmsd.cpp:31-33`)
additionally requires deciding the fate of `rmsd/rmsd_assignment.cpp` (TODO.md 7, operator decision).

**S5. GFN-FF dead members** (`src/core/energy_calculators/ff_methods/`):
- `gfnff_method.cpp:12400-12487` `bool GFNFF::calculateTopologyCharges(TopologyInfo&) const` + decl `gfnff.h:2554`
  (all `calculateTopologyCharges(` calls are `m_eeq_solver->...`)
- `gfnff_method.cpp:12717-12813` `bool GFNFF::calculateFinalCharges(TopologyInfo&, int, double)` + decl
  `gfnff.h:2675` + tombstone comments `:12489-12491`, `:12811-12813` (was F-Q2)
- declared, never defined, never called: `gfnff.h:2344` (`Vector calculateTopologyCharges(const Vector& cn, ...)`),
  `gfnff.h:2401` (`Vector calculateFinalCharges(const TopologyInfo&) const`), `gfnff.h:1660`
  (`generateBatmTriplesNative`)
- `gfnff_method.cpp:11567-11601` overload `getEEQParameters(int, const TopologyInfo&)` + decl `gfnff.h:2320` (every
  call site passes one argument: `:7817,7831,8003,10697,12841`)
- `gfnff_method.cpp:12940-12977` `NumGradFixedCharges` + decl `gfnff.h:702`
- `gfnff_method.cpp:12258-12311` `extractDispersionConfig` + decl `gfnff.h:1673`
- `gfnff_method.cpp:773-780` `getEEQDistanceCutoff` + decl `gfnff.h:990`
- `gfnff_method.cpp:13051-13060` `canSkipD4GaussianWeightsUpdate`, `recordD4CNValues` + decls `gfnff.h:1395-1396`
  (check the member variables they touch for other readers before removing those)
- `gfnff_torsions.cpp:1381-1524` `calculateDihedralGradient` (decl `gfnff.h:1902`), `gfnff_torsions.cpp:2255-2269`
  `generateGFNFFSTorsions` (decl `:1552`)
- `gfnff_inversions.cpp:339-464` `calculateInversionGradient` (decl `:1975`), `:537-655`
  `getGFNFFInversionParameters` (decl `:1996`), `:889-906` `generateGFNFFInversions` (decl `:1565`)
- Not dead but a no-op: `calculateTopology()` (`gfnff_method.cpp:4490-4502`, always `return true`, called once at
  `:4024`) - inline away together with the stale log text (F-E7).
- About 1 000 lines in total; verify additionally with `ctest -L structinv` and one GMTKN55/MOR41 regression run
  (`scripts/refset_regression.py`) since these files are the GFN-FF core.

**S6. Force-field numerical-gradient chain.** `ForceFieldMethod::calculateNumericalGradient`
(`ff_methods/forcefield_method.cpp:420-431`) has no caller; it is the only caller of `ForceField::NumGrad`
(`forcefield.cpp:575`, via `forcefield_method.cpp:426`). Also `ForceField::setESPs` (`forcefield.cpp:556-573`) and
`ForceFieldMethod::isParameterCachingEnabled` (`:343`).

**S7. GPU no-op accessors.** `getDeviceCoordsPtr()` returns `nullptr` and has no caller:
`cuda/ff_workspace_gpu.h:579-581`, `cuda/ff_workspace_gpu.cu:4324-4327`, `rocm/gfnff_rocm.hip:10338-10341`.
`FFWorkspaceGPU::shouldSkipDc6dcn()` (`cuda/ff_workspace_gpu.h:647`, TODO.md 7). Verify with a CUDA and a ROCm build.

**S8. `ForceFieldGenerator::GenerateUFFD3Parameters()`** (`forcefieldgenerator.cpp:1157-1256`, decl
`forcefieldgenerator.h:170`): no caller. Delete and fix `ff_methods/CLAUDE.md:55` in the same commit (L8).

**S9. Empty stubs.** `SimpleMD::Rattle_Verlet_First/Rattle_Constrain_First/Rattle_Verlet_Second`
(`simplemd.cpp:4384-4395`, decls `simplemd.h:448` ff.) - empty bodies, no callers.

**S10. Unwired Casino restart code** (`casino.cpp:1128-1168 performCheckpoint`, `:1200-1240 detectAndLoadRestart`,
`:1282-1295 loadFromRestartFile`, `:1296-1331 handleTrajectoryRestart`): no callers, yet `casino.h:118-124`
advertises `auto_restart`, `checkpoint_frequency`, `checkpoint_on_signal`. Either wire it or delete code and
options together (deleting the options changes `-help` output; treat as a behaviour-visible change).

**S11. Remaining grep-dead members** (121 functions, about 2 370 lines in total, list in `dead_filtered.txt`).
Low-risk examples: `EnergyCalculator::initializeCommon/OrbitalOccuptations/setGeometryFile/PrintAvailableMethods/
SaveToFile/GetMethodInfo/resetVerbosity/checkForNaN/getNonBondedEnergy` (`energycalculator.cpp:170-912`),
`EEQSolver::solveWithPCG` (`eeq_solver.cpp:2737-2828`), `XTB::checkConvergence_impl` (`xtb_scf.cpp:1063-1109`,
was X-S2), `TrajectoryWriter::writeScatteringPerFrameFiles` (`trajectory_writer.cpp:1427-1527`),
`D4ParameterGenerator::calculateFrequencyDependentPolarizabilities/getEffectiveC6/getR4OverR2/
getAtomicPolarizability`, `MethodFactory::checkCompilationFlag/isUlyssesMethod/getMethodInfo`,
`OrcaInterface::createInputFile/executeOrcaProcess/runOrca`, `Molecule::getConnectivtiy/AnalyseIntermoleculeDistance/
ElementsRemoved/hasMixedSystem/getAtomicAtoms/WriteJsonFile`, `CurcumaLogger::addCitation/active_sink_count/
get_elapsed_time`, `ParameterRegistry::findModule/modulesForCommand`. Some are setters of public classes
(`OptimizerDriver::set*`, `TrajectoryStatistics::get*`, `LBFGS::*`): delete only after operator review. Exclude
`CurcumaLogger::debug_var/debug_timing` (used through macros) and `ParameterRegistry::addDefinition/addModule`
(called by generated code from `scripts/param_parser/main.cpp:216,281`).

### Build and text cleanups

**S12. `CMakeLists.txt`**: remove duplicate `analysenciplot.cpp` (`:511`); remove `add_definitions(-DGIT_COMMIT_*)`
(`:313,317,318`; verify `grep -rnw "GIT_COMMIT_HASH\|GIT_COMMIT_DATE\|GIT_COMMIT_TAG" src` -> only `version.h.in`);
drop the misleading PCH status messages (`:712-719`); keep one `REUSE_FROM` (`:765-766`). Larger but still
mechanical: delete `curcuma_cap` and link the two tests to `curcuma_core` only (M1; verify `ctest -R
"external_potentials|measurement"`).

**S13. `main.cpp:723`**: drop `"modern_optimizer"` from `scope_modules` (no `BEGIN_PARAMETER_DEFINITION
(modern_optimizer)` exists; TODO.md 7). Verify `ctest -R cli_` unchanged (a flag aimed at that module would now warn
instead of being silently accepted - check no test relies on it).

**S14. Delete the 27 commented-out code blocks and the `//` code lines** (L2), the 13 `#pragma message` lines (M7)
and the comment-only tombstones (L3, L5). Verification: compile only; `git diff --stat` touches comments only
(`git diff -U0 | grep '^[-+]' | grep -v '^[-+]\s*//\|^[-+]\s*/\*\|^[-+]\s*\*'` empty apart from pragma lines).

**S15. ASCII-only output** (H1): replace `─/═/•/→/✓/✅/↔/└─` in string literals and `Å/σ/°/±` in terminal output
with ASCII (`-`, `=`, `*`, `->` is also forbidden by root CLAUDE.md so use `to`, `OK`, `Angstrom`, `sigma`, `deg`,
`+/-`). Not strictly "no behaviour change": output text changes, so `cli_*` tests that compare output may need new
references. Do it per file and run `ctest -R cli_` after each.
