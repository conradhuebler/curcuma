# CLAUDE.md - Capabilities Directory

## Overview

The capabilities directory contains high-level molecular modeling applications and computational chemistry tasks. These modules use the core computational engines to perform complex multi-step calculations and analyses.

## Structure

```
capabilities/
├── confscan.cpp/h         # Conformational scanning along reaction coordinates
├── confsearch.cpp/h       # Systematic conformational searching
├── curcumaopt.cpp/h       # Geometry optimization algorithms
├── simplemd.cpp/h         # Molecular dynamics simulation
├── rmsd.cpp/h            # Structure comparison and alignment
├── rmsdtraj.cpp/h        # Trajectory RMSD analysis
├── hessian.cpp/h         # Second derivative calculations
├── persistentdiagram.cpp/h # Topological data analysis
├── optimiser/            # Optimization algorithms
│   ├── lbfgs.cpp/h       # LBFGS optimization
│   └── LevMar*.h         # Levenberg-Marquardt variants
└── c_code/               # C interface code (Hungarian algorithm)
```

## Key Capabilities

### Conformational Analysis
- **ConfScan**: Systematic scanning of conformational space along defined coordinates
- **ConfSearch**: Automated conformational searching with energy filtering
  - Dual-method (Jun 2026): `md_method` (explore + pre-opt) vs `opt_method` (per-cycle Phase 3b re-opt + final ranking); both empty -> `method`. Phase 3b skipped when equal. See [docs/CONFSEARCH_DUAL_METHOD.md](../../docs/CONFSEARCH_DUAL_METHOD.md)
  - Restart (Jun 2026): `-restart` writes a self-contained checkpoint (bias pool + cumulative + seeds + energies + schedule) after every MD/cycle to CWD + BMT; resume skips pre-opt, restores bias pool (`SharedBiasPool::restoreStructures`), continues from `next_T`. See [docs/CONFSEARCH_RESTART.md](../../docs/CONFSEARCH_RESTART.md)
  - **Registry-backed (Jul 2026)**: ConfSearch owns 67 PARAMs (`confsearch.h`), so `-opt_method`/`-thermostat`/`-restart`/`-wall_*` are no longer stolen by the auto-router. Reads go through `ConfigManager` (resolves aliases). `-T` is not a ConfSearch parameter — use `-startT`/`-endT`. 11 SimpleMD-legacy keys (`unique`, `respa`, `dipole`, `impuls_scaling`, `MaxTopoDiff`, `rescue`, `cleanenergy`, `wall_xl/yl/zl`, `rattle_tolerance`, `printOutput`) are deliberately unregistered: no module owns them, so claiming them would steal them from `-md`
  - **`ChildConfig()`/`FilterConfig()` (Jul 2026)**: one source of truth for every child computation (MD, 4 opt sites, 2 ConfScan filters) — carries charge/spin/gpu/verbosity + method sub-scopes (`gfnff`/`xtb`/…). Replaces ad-hoc `{method,threads,gpu}` JSONs that dropped charge and solvation
  - **Phase 3c torsion recombination (Jul 2026)**: `-confgen_phase true` runs ConfGen on the cycle's deduplicated minima, appends the chemically valid new conformers to `<base>.bias.opt.accepted.xyz` and lets every downstream phase treat them like metadynamics hits (no special case). Off by default -- costs one optimisation per proposal. See [docs/CONFSEARCH_PROPOSALS.md](../../docs/CONFSEARCH_PROPOSALS.md)
  - **Reporting + per-cycle output (Jul 2026)**: reference energies are taken after the INITIAL optimisation (not after cycle 1's MTD) and reported per level of theory without the meaningless cross-method delta; every phase reports its ensemble's energies (`-ensemble_report N`); each cycle writes `<base>.cycleNN_TxxxK.<method>.xyz` + `<base>.best_per_cycle.<method>.xyz` on BOTH levels (`-cycle_output`); `CurcumaLogger::section()` gives visible block structure at verbosity 1. See [docs/CONFSEARCH_REPORTING.md](../../docs/CONFSEARCH_REPORTING.md)
  - **`-topology_factor` (Jul 2026, default 1.3)**: every topology comparison used `Molecule::DistanceMatrix()` at factor 1.5, which counts a compressed 1-3 contact as a bond (107-atom peptide: a 1-3 C...C pair at 2.25 A becomes "bonded" at 2.28 A, i.e. on a 2-degree angle fluctuation). Measured: 92 % of MD frames counted as topology changes at 1.5, 0 % at 1.3 — Phase 4 was discarding most conformers. One cycle at 500 K: 3 -> 34 conformers, lowest gfnff -18.7594 -> -18.7800. Same value ConfGen already used
  - **Topology gate (Jul 2026, `-snapshot_topology_gate` default ON)**: MD snapshots whose bond topology differs from the reference are dropped BEFORE Phase 2 — they are rejected by the Phase 4 filter anyway (after a full optimisation) and are the geometries that make GFN-FF return a finite energy with a NaN gradient (1/r in `hb`/`atm`/`batm`). Measured: butane 3000 K, 15 of 974 kept; 800 K, nothing rejected. **`-snapshot_clash_ratio` (0.55)** screens raw distances too — an already-bonded pair that collapses forms no new bond, so the topology test misses it while hb/atm/batm 1/r still blow up. **Failed structures are saved** to `<stage>.failed.xyz` (input geometry, named `failed_NNN_<reason>`; `OptThread::methodReportedNaN()` supplies the signal, since a NaN gradient reaches OptimizationResult as plain non-convergence). **NaN guards**: comparisons with NaN are false, so a diverged optimisation passed the energy-rise limit and every convergence test — the driver now aborts on non-finite energy/coordinates, and no non-finite geometry is ever written into the pool. Empty result → cycle skipped with a warning. Related: `LBFGSppObjectiveFunction` now treats a non-finite GRADIENT as an error (it only checked the energy). **Repair (opt-in, `-repair_snapshots`)**: instead of discarding, restrain the offending pairs back (missing bond -> reference length, spurious contact -> outside bonding range), optimise, release, optimise freely, keep only if the topology matches. Bounded: `-repair_max` 20/cycle, `-repair_max_bonds` 2. Measured butane 3000 K: 0/20 and 3/20 repaired
  - **File naming (Jul 2026)**: `<base>.<purpose>.<method>[.opt][.accepted].xyz` before the temperature loop, `<base>.cycleNN_TxxxK.<purpose>.<method>[...].xyz` inside it (`stageBase()` / `cycleStage()`). Purposes: initial, explore, bias_pool, refine_crude, refine, ensemble. Per-cycle files are no longer overwritten. Final result keeps `<base>.cumulative.opt.accepted.xyz`
  - **Restart fix (Jul 2026)**: checkpoint energies are read null-tolerantly — a non-finite energy (e.g. `initial_energy_opt` in a single-method run) is written by nlohmann as JSON `null`, and `value<double>()` throws on null, which aborted every single-method `-restart` resume
  - **Batch progress (Jul 2026)**: `PerformOptimisation` was silent during AND after its batch (OptThread leaves the shared logger level at 0; it was restored only at function end, swallowing the per-structure table and the summary). Level now restored right after `StartAndWait()`; an `OptThread::setOnComplete` hook drives `CurcumaLogger::progress_bar` (stdout, `i/N`, honours `-noprogress`/TTY) at verbosity 1, per-structure lines at >= 2. Same for the MD batch. CxxThreadPool bars stay at verbosity >= 3 — they are stderr/percent-only and are the per-permutation flood in ConfScan (pool reused per candidate); in ConfSearch one pool == one batch
  - **Two-stage Phase 3b (Jul 2026, opt-in)**: `-phase3b_two_stage true` = crude opt of everything -> dedup at that level -> accurate opt of the survivors. The intermediate filter has its energy cutoff disabled (dedup only; crude energies must not rank). Flags: `-phase3b_preopt_preset`, `-phase3b_preopt_max_iter`, `-phase3b_preset`, `-phase3b_filter`
  - **RMSD-aware seeding (Jul 2026)**: `SelectSeeds()` — `-seed_selection diverse` (default) walks the energy ranking but requires `-seed_min_rmsd` (default `2 x -rmsd`) permutation-aware spacing between seeds, relaxing the spacing only if too few basins exist; `-seed_selection energy` = old pure ranking. No effect at the default `-seed_rank 1`. See [docs/CONFSEARCH_SEEDING.md](../../docs/CONFSEARCH_SEEDING.md)
- **Geometry restraints** (`optimisation/geometry_restraints.h`, Jul 2026): harmonic dihedral `E = 1/2 k (phi - phi_0)^2` AND distance `E = 1/2 k (r - r_0)^2` (config keys `dihedral_restraints` / `distance_restraints`) added at the two E/G choke points (`LBFGSppObjectiveFunction::operator()`, `OptimizerDriver::evaluateEnergyAndGradient`) via `config["dihedral_restraints"]`. Energy AND gradient (a line search accepts by energy); geometry Å, gradient Eh/Bohr (×0.529177); deviation wrapped to (−π, π]. FD-verified 7.9e-12 (dihedral) / 3.6e-11 (distance) Eh/Bohr — `test_geometry_restraints`
- **ConfGen** (`confgen.cpp/h`, Jul 2026): ensemble analysis for targeted conformer proposals. **Restrained build (P0, Jul 2026, `-restrained_build` default ON)**: when the rigid torsion setting clashes, the proposal is built by driving the torsions from the clash-free template with dihedral restraints instead — same clash/topology gates afterwards, free optimisation for the reported energy, proposals whose torsions miss their target by >30° dropped. Uses `torsion_space.{h,cpp}` (rotatable-bond detection via BFS ring test, rotamer-state clustering, `setDihedral` with numerically verified rotation sign) + matched-pair per-term energy differences. Found and fixed a missing `Repulsion` term in `getEnergyDecomposition()` (was off by 2781 kJ/mol for GFN-FF). `-generate true` builds recombined state vectors, optimises them and judges them with the force field (mandatory topology check with an explicit `-topology_factor`, like-for-like energy reference from re-optimised templates). Tests: `test_torsion_space` (32 assertions), `cli_confgen_01_matched_pairs`, `cli_confgen_02_generate`. See [docs/CONFSEARCH_PROPOSALS.md](../../docs/CONFSEARCH_PROPOSALS.md)
- **RMSD Analysis**: Structure comparison, alignment, and trajectory analysis
  - JSON output: `<target>.rmsd.json` with rmsd, rmsd_raw, permutation, reference_xyz, reorder_xyz, file provenance (always generated, even no-reorder mode)

### Optimization

> **All non-external optimizers below are 🤖 AI-generated. None are ✅ TESTED or ✅ APPROVED.**

- **CurcumaOpt**: Geometry optimization dispatcher — legacy system, human-tested
- **`-opt` multi-XYZ** (`main.cpp`/`optimizer_factory.cpp`):
  - ⚙️ Machine-tested — all frames are optimised and written to `.opt.xyz` in input order.
  - Parallel dispatch with `-threads N`; workers are independent, but step-table output is suppressed during the batch to avoid interleaved stdout.
  - After the batch finishes, an ordered per-frame summary is printed (index, status, iterations, final energy).
  - Live `CxxThreadPool` progress bar for parallel batches is pending an update of `external/CxxThreadPool` (see `docs/OPT_MULTIXYZ_PARALLELISM_WP.md`).
- **LBFGSpp**: Wrapper around external LBFGSpp library — external code, wrapper is AI-generated
- **ANCOPT** (`ancopt_optimizer.cpp/h`): 🤖 AI-generated port of XTB's AncOpt (Grimme)
  - ⚙️ Machine-tested: CH4/UFF converges (4 steps); Tier L path runs on 1410-atom polymer+UFF
  - Large-system enhancements (Apr 2026, all 🤖 AI-generated, not ✅ TESTED):
    - **Tier L** (600–2000 atoms): Truncated Lanczos ANC (`generateANCLanczos`, top-k modes), implicit T/R projection (O(N²) vs O(N³)), `detrotra8` gated behind `n3<=1800`
    - **Tier XL** (>2000 atoms): L-BFGS in ANC subspace (`calculateLBFGSStepInternal`), drops dense nvar×nvar Hessian
    - **Shared RF solver**: `RFSolver::lanczosLowestEigenpair` + `calculateRFStep` in `optimisation/rf_solver.h/.cpp`
    - **Structured advisory**: tier, algorithm path, per-phase timing at verbosity 2
    - **EIGEN_USE_LAPACKE**: per-file in rf_solver.cpp (safe, no `I` variable conflict)
  - Not tested: QM gradients, transition states, linear molecules, XL tier (>2000 atoms)
- **OptimizerDriver** (`optimizer_driver.cpp/h`): 🤖 AI-generated base class (Template Method)
- **OptimizerFactory** / **OptimizationDispatcher**: 🤖 AI-generated
- **Native L-BFGS / DIIS / RFO** (`optimisation/lbfgs.cpp`): 🤖 AI-generated — see `optimisation/CLAUDE.md`
- **Constrained optimization**: Fix atomic positions by setting gradient = 0

### Molecular Dynamics
- **SimpleMD**: Basic molecular dynamics with various thermostats
  - ✅ Coarse-graining support with automatic system detection
  - ✅ PBC wrapping for periodic boundary conditions
  - ✅ 10x timestep scaling for pure CG systems
  - ✅ VTF trajectory output for CG systems
  - ✅ Orientational dynamics infrastructure (prepared for Phase 6 ellipsoids)
- **NEB Docking**: Nudged elastic band for transition state searches
- **Trajectory Analysis**: Analysis tools for MD trajectories

### Advanced Analysis
- **Hessian**: Second derivative calculations for normal modes
- **Persistent Diagrams**: Topological data analysis for molecular structures
- **Enhanced TDA (dMatrix replacement)**: Complete topological data analysis with TDAEngine
- **Pairmapper**: Advanced structure matching algorithms

### BMT Output Directory Integration (CurcumaMethod)
- **Default**: `CurcumaMethod::createBMTDir(keyword)` creates `Basename.Keyword.YYYYMMDD_HHMMSS/` and sets `m_output_dir`
- **`initializeBMT()`** (main.cpp): Helper that calls `setFile()`, `createBMTDir()`, and registers `-bak` files
- **`addBakFile()`** / **`processBakFiles()`**: Register and copy files back to CWD after calculation
- **`outputPath()`**: Route all output through BMT directory when set; returns bare filename when BMT is disabled
- **Commands using BMT**: md, opt, hessian, qmdfffit, confsearch, confscan, confstat, dock, analysis, rmsd
- **Standalone BMT**: `BMTUtils::` functions used directly for analysis/rmsd (non-CurcumaMethod handlers)
- **Status**: 🤖 AI-generated, machine-tested — human production testing pending

## Development Guidelines

### Interface Design
- All capabilities use the core EnergyCalculator for energy/gradient evaluations
- Consistent parameter handling through JSON-based configuration
- Progress reporting and result persistence for long calculations

### Performance Considerations
- Multi-threading support where applicable
- Memory-efficient handling of large trajectory data
- Automatic checkpointing for resumable calculations

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Future development tasks and visions to be defined by operator/programmer*

## Variable Section

### Parameter System - ConfigManager Layer
✅ **Status**: Production-ready across 4+ modules (analysis, rmsd, simplemd, confscan)
- ConfigManager provides type-safe parameter access with hierarchical dot notation
- Multi-module parameter routing fixed (Oct 26, 2025)
- Migration ongoing: Replace Json2KeyWord calls with `config.get<T>("param")`
- Reference: `src/core/config_manager.h`, example: `analysis.cpp`

### Current Development
- Enhanced conformational search algorithms
- Improved trajectory analysis tools
- Better integration with quantum chemical methods
- ⚙️ **IMPLEMENTED** (AI): Parameter Registry System, dMatrix Integration, RMSD Code Restructuring
- ⚙️ **IMPLEMENTED** (AI, Oct 29, 2025): Help System Dynamic Generation
- ⚙️ **IMPLEMENTED** (AI, Oct 30, 2025): SimpleMD CG Integration Phase 1-4
- ⚙️ **IMPLEMENTED** (AI, Nov 2025): SimpleMD CG Integration Phase 5 - VTF trajectory output
- ⚙️ **IMPLEMENTED** (AI, Jan 2026): Analysis Output Refactoring - Registry-based handler architecture
- Pending: Unit system migration, RMSD Strategy pattern (Phase 3), CG Phase 6 (ellipsoidal extensions)

### New Analysis Output Architecture
- **New Analysis Output Architecture**: Handler-based system with registry pattern
- **File Naming Schema**: basename.general.csv, basename.NNN.type.csv, basename.type_statistics.csv
- **Extensible Design**: IAnalysisOutputHandler interface for new analysis types
- **Benefits**: Eliminates duplication, single point of change, automatic file generation

### Known Issues
- Memory optimization needed for large systems (>1000 atoms)
- ConfScan output: at verbosity 1 each pass shows a live progress bar (on by default; `-confscan.progress false` / global `-noprogress` / non-TTY disable it) and a clean per-pass summary; per-structure `Accept/Reject` detail is at verbosity >=2. The internal per-batch `CxxThreadPool` bars (reorder/check pools, and ConfSearch's own MD/opt pools) are separate detail output gated at **verbosity >= 3** — below that they are forced to `ProgressBarType::None`, so a default ConfSearch run shows only the overall progress + phase lines, not one pool bar per permutation batch (Jul 2026).
- **ConfScan reorder geometry truncation (⚠️ open, found June 2026)**: `ConfScan::Reorder` calls `mol1->ApplyReorderRule(t->ReorderRule())` (`confscan.cpp:1439`) with a rule sized to the *compared* atom set. Heavy-only (`-rmsd.protons false`) → reordered structures written with heavy-atom count (e.g. 55 vs 114); `get_rmsd` → empty rule → reordered/rejected structures written with 0 atoms. Filtering counts are correct; only the written geometry of reordered structures is reduced. Fix: rule must map all atoms, or `ApplyReorderRule` must preserve unmapped atoms. CLI tests 02/06 surfaced this.
- **ConfScan threaded reorder stall (⚠️ open, found June 2026)**: the reorder path with `threads>1` can intermittently deadlock/stall (observed: heavy-only `-rmsd.protons false`, and restart-resume double-run; at high and low load). CLI tests 06 (heavy) and 07 (restart) pin `threads=1`.
- **GFN-FF topology cache made optimisation energies drift — fixed for ConfSearch (Jul 2026)**: with `topology_mode=auto` every mid-run topology re-derivation shifts the energy scale; optimising a hot MD snapshot then follows those jumps (measured: "converged" to -9.168083 Eh, written geometry really worth -8.668213 Eh, same 52 bonds). One such structure in the cumulative pool becomes the reference and collapses the result to "1 unique conformer of 482". `ChildConfig()` now sets `gfnff.topology_mode = constant` for every child (explicit user setting still wins)
- **ConfScan IS energy-sorted — verified (Jul 2026)**: `m_ordered_list` is a `std::multimap<double,int>` keyed on the energy, so no `std::sort` is needed or wanted. Proven on a deliberately shuffled 44-conformer ensemble (accepted file strictly ascending, same 14 structures as unshuffled). The earlier claim "ConfScan sorts nowhere" was wrong; `Reorder` and the accepted write-out now re-establish the order locally through a multimap so the property is checkable where it is relied on. The `break`->`continue` change of the energy cutoff stays (robust, matches the write-out path); the real cause of `Processed: 1 / 200` was an artefact structure with a corrupt energy acting as `m_lowest_energy`
- **ripser was not re-entrant across instances — patched at configure time (Jul 2026)**: `external/ripser/ripser.h` declares its simplex enumerators as function-local `static`s holding `const ripser&` + `const binomial_coeff_table&`; from the SECOND ripser instance in a process on, those reference a destroyed object. curcuma builds one instance per structure (ConfScan descriptor), so every descriptor after the first was read from freed memory — usually silent (same heap address), a hard segfault as soon as the allocation pattern changed (reproduced with two ConfScan passes in one ConfSearch cycle). CMake strips the `static` idempotently; the durable fix belongs in the conradhuebler/ripser fork
- **Bias structures carried the reference molecule's energy — fixed (Jul 2026)**: `BiasStructure::energy` was never set at deposit and the ConfSearch export copies a reference molecule; snapshots therefore showed a wrong but plausible energy. SimpleMD stores `m_Epot`, ConfSearch sets it on all export paths
- **Dual-method final statistics mixed both PES — fixed (Jul 2026)**: the closing "search lowered the energy by X kJ/mol" line subtracted the md_method initial energy from the accepted pool's minimum, which is an opt_method energy in dual mode (`gfnff: -9.005 -> -85.017 Eh`, "199567 kJ/mol"). Now: ranking gain on the opt PES only, plus a separate `exploration (<md_method>) lowered its own minimum by ...` line (md initial -> md running best). `initial_energy_opt` added to the checkpoint so a `-restart` resume keeps the opt-PES anchor
- **Legacy RMSD-MTD restart aborted the run — fixed (Jul 2026)**: `restart["rmsd_ref_file"]` stored a bare `<basename>.mtd.xyz` although every writer routes through `outputPath()` (BMT dir), and the reload fed it to `FileIterator`, which throws → uncaught exception, exit 1. Path now recorded via `outputPath()`, reload resolves a candidate list and degrades to a warning + empty bias pool, `-rmsd_mtd_ref_file` misses are reported instead of thrown, `m_bias_json` index bounds-checked. Trigger was `BMTUtils::createBMTDir()` handing two same-second runs the *same* directory (now `_2`, `_3`, … suffix), which also let run 2 resume run 1's `curcuma_restart.json`. Final hills dump now truncates exactly once (`first_hill_dump`) instead of relying on `i == j && i == 0`; shared-pool `.mtd.xyz` writes routed through `outputPath()`. `cli_simplemd_13_rmsd_mtd_legacy_ab` is green again
- **RATTLE froze inter-fragment motion — fixed (Jul 2026)**: `SimpleMD::Rattle()` called `RemoveRotations()` unconditionally on every constrained step, zeroing every *fragment's* linear+angular momentum each step. Any collective force (RMSD-MTD bias, walls) could then no longer move fragments relative to each other, and the removed `6*nfrag` DOF were not in the DOF bookkeeping (reported T too low → thermostat overheats). ConfSearch turns RATTLE on automatically for every cycle at `T >= rattle_threshold_temp` (400 K), so multi-fragment ConfSearch runs looked like "the bias does not grip". Removed; COM/rotation removal is the main loop's job (`remove_com_motion`/`remove_com_mode`, applied for both integrators). A/B acetic-acid dimer (3 ps, gfnff, T=500 K, seed 42, `-rattle 2`): fragment-COM span 0.43 → 2.17 A, hills 44 → 127; single-molecule runs unchanged (butane 55→54, 90-atom 249→249 hills). Pinned by `cli_simplemd_15_rattle_mtd_fragments`
- SimpleMD wall potential: ✅ **Sign errors fixed (June 2026)** — `ApplyRectHarmonicWalls` (min-wall term was subtracted instead of added → atoms below `r_min` pushed further out) and `ApplySphericLogFermiWalls` (`gradient -= dV/dr` flipped the radial force outward). Both now add `+dV/dr` to `m_eigen_gradient` (= dE/dr; force = -gradient), matching `ApplySphericHarmonicWalls` and `ApplyRectLogFermiWalls` which were already correct.
- **ConfSearch efficiency/robustness (Phase A-C) — roadmap & open TODOs**: see [docs/CONFSEARCH_ROADMAP.md](../../docs/CONFSEARCH_ROADMAP.md). Big items: (1) verbosity-ownership rework (global CurcumaLogger level is leaked/clamped by sub-objects → ConfSearch logs hidden, RATTLE report forced to std::cout — `FIXME` in `SimpleMD::InitConstrainedBonds`); (2) `CitationRegistry::cite` thread race → crash at gfnff `threads>1` (workaround: threads=1); (3) Phase C `cluster`/`weighted` calibration is experimental/unvalidated. Cross-run bias heating bounded by **defaults ON**: `rmsd_mtd_freeze_inherited`+`temp_abort` (`rmsd_mtd_max_height` opt-in); bare `-startT 500` no longer blows up (TODO #4; intra-run wide-hill blow-up still open). Verbosity-ownership (1) is now largely resolved (CurcumaMethod base RAII); as of Jul 2026 the ~30 unconditional `std::cout`/`fmt::print` sites in `simplemd.cpp` (incl. the RATTLE report) are gated on `m_verbosity >= 1`, so ConfSearch runs its child MD at verbosity 0 (silent) at ConfSearch verbosity 1 and prints a compact `MD runs: X/N done | ~R running | last <s> | elapsed <s>` counter instead; MD-module output returns at verbosity 2. Standalone `-md` at verbosity 1 is unchanged.

### Unported Features from Old CurcumaOpt (TODO)

The following features existed in `curcumaopt.cpp` but are **not yet ported** to the new `OptimizerDriver`/`OptimizationDispatcher` system:

1. **Parallel batch optimization** (`curcumaopt.cpp:276`) — `ProcessMolecules()` used `CxxThreadPool` (SPThread/OptThread) to optimize multiple molecules in parallel. The new system processes batches serially (loop over single-molecule calls). Implement `OptimizationBatchRunner` using CxxThreadPool.

2. **Hydrogen-only optimization** (`opt_h`, `curcumaopt.cpp:521`) — per-atom constraints derived from element type (element==1 → movable). Useful for optimizing H positions while heavy atoms are fixed. Needs integration into `OptimizerDriver` constraint system.

3. **Hessian after optimization** (`hessian=1/2`, `curcumaopt.cpp:237,323`) — computes normal modes using `Hessian` class after optimization, saves `hessian.json` and `scf.json`. Port by calling `Hessian` class inside `executeOptimization()` in `main.cpp` when `hessian` parameter is set.

4. **Molecular orbital diagram** (`mo_scheme`, `curcumaopt.cpp:404`) — `WriteMO()` generates TikZ LaTeX orbital diagram, `WriteMOAscii()` ASCII variant. Parameters: `mo_homo`, `mo_lumo`, `mo_scale`. Needs EnergyCalculator `Energies()` / `NumElectrons()` results plumbed through result struct.

5. **"stop" file interrupt** (`curcumaopt.cpp:665`) — checks for `./stop` file on disk during optimization loop to allow early termination. Simple to add to `OptimizerDriver::Optimize()` loop.

6. **`fusion` mode** (`curcumaopt.cpp:688`) — skips `Molecule::Check()` validity gate, needed for unusual bonding / fusion compounds. Add `bool fusion` flag to `OptimizationContext` and check in driver loop.

7. **Dipole moment output for GFN2** (`curcumaopt.cpp:395`) — printed via `EnergyCalculator::Dipole()` after SP/opt when method is gfn2. Not available in current `OptimizationResult` struct.

**Priority**: Items 1 (parallel batch) and 2 (opt_h constraints) are most commonly used. Items 3-7 are lower priority.

---

*This documentation covers all molecular modeling capabilities and applications*