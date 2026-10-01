# CLAUDE.md - Capabilities Directory

User-facing commands built on `src/core`: optimisation, MD, conformer search and scan, RMSD, trajectory analysis,
Hessian, docking. Energies and gradients always come from `EnergyCalculator`. Docs: [CONFSEARCH_DUAL_METHOD](../../docs/CONFSEARCH_DUAL_METHOD.md),
[CONFSEARCH_RESTART](../../docs/CONFSEARCH_RESTART.md), [CONFSEARCH_ROADMAP](../../docs/CONFSEARCH_ROADMAP.md), BMT in [../tools/CLAUDE.md](../tools/CLAUDE.md).

## Layout and ownership

- `curcumamethod.*`: `CurcumaMethod` base (scoped logger verbosity, BMT dir, `outputPath()`, `-bak` files); 15 built capability classes derive from it
- `optimizer_driver.*`, `optimizer_factory.*`, `optimizer_interface.*`: `OptimizerDriver` loop, `OptimizerFactory`, `OptimizationDispatcher`; the only `-opt` path
- `lbfgspp_optimizer.*`, `ancopt_optimizer.*`, `native_optimizer_adapters.*`: concrete drivers; native algorithms in [optimisation/](optimisation/CLAUDE.md)
- `rmsd.*` (`RMSDDriver`, strategies in [rmsd/](rmsd/CLAUDE.md)), `confscan.*`, `confsearch.*`, `shared_bias_pool.*`, `simplemd.*`, `analysis*`, `handlers/`, `trajectory*`
- Others: `hessian`, `docking`, `nebdocking`, `qmdfffit`, `persistentdiagram`, `tda_engine`, `pairmapper`, `casino`, `polymerbuild`, `confstat`, `rmsdtraj`
- `optimiser/`: LevMar headers (docking, NEB, qmdfffit), `OptimiseDipoleScaling.h` (main.cpp); `c_code/`: C Hungarian solver used by `rmsd`
- Dead: `curcumaopt.cpp` (legacy `CurcumaOpt`), `native_lbfgs_optimizer.cpp`, `optimisation/modern_optimizer_simple.cpp` are not compiled;
  `munkress_2.h`, `optimiser/{LBFGSppInterface,LevMarNEBPseudoFF,Proton,XTBDocking}.h` are included nowhere (grep)

## Checklists

- **New capability**: derive from `CurcumaMethod`, `PARAM` block in the header, read via `ConfigManager`; add the `.cpp` to `curcuma_core_SRC` and `curcuma_cap_SRC`
- ... dispatch in `main.cpp`, call `initializeBMT()` before the run, write every file through `outputPath()`; test with and without `-no_bmt`
- ... register every flag the command reads: a flat flag that another module registers is routed to that module
- **New optimizer**: `OptimizerType` + names in `parseOptimizerType()`/`optimizerTypeToString()`, creator and description in `optimizer_factory.cpp`,
  `-opt` help in `main.cpp` `executeOptimization`, `optimizer` PARAM help in `curcumaopt.h`; implement the pure virtuals (`optimizer_driver.h:149-159`, `:226-230`)
- **New analysis output**: implement `IAnalysisOutputHandler` (4 virtuals), `registerHandler()` in the `AnalysisOutputDispatcher` constructor

## Invariants and traps

- `-opt` PARAMs live in `curcumaopt.h` although `CurcumaOpt` is not built: the registry is generated from every header (`GLOB_RECURSE`), built or not
- `-optimizer auto` (default) always means LBFGSpp; `selectOptimalOptimizer()` has no caller. `lbfgs` is the native L-BFGS; unknown names fall back to LBFGSpp
- `OptimizerDriver` owns loop, convergence and output; every abort returns the last accepted structure, which `-opt` writes even if not converged
- `stall_steps`/`stall_rmsd` (20 / 1e-6 A) end a frozen run as "No progress"; a zero step counts as converged only if the driver criteria hold
- Multi-XYZ `-opt -threads N`: one `CxxThreadPool` worker and `EnergyCalculator` per frame, logger silenced, then an input-order summary (`optimizeBatch`)
- ConfSearch registers its flags (71 PARAMs) so the auto-router leaves them alone; `-T` is not one (`-startT`/`-endT`); unregistered SimpleMD keys: `confsearch.h:241-244`
- ConfSearch children (MD, four opt sites, ConfScan filter) take their config only from `ChildConfig()`/`FilterConfig()` (charge, spin, gpu, verbosity, method scopes)
- Verbosity 1: ConfSearch child MD silent plus an `MD runs:` counter; ConfScan/ConfSearch pool bars only at verbosity >= 3; ConfScan pass bar obeys `-confscan.progress`/`-noprogress`
- `RMSDDriver` writes `<target>.rmsd.json` (rmsd, rmsd_raw, permutation, reference_xyz, reorder_xyz, file names) to BMT and CWD, also without reordering
- BMT via `initializeBMT()`: md, hessian, qmdfffit, confsearch, confscan, confstat, dock, rmsd, polymerbuild; via `BMTUtils::`: sp, opt, analysis
- Analysis files: `basename.general.csv`, `basename.NNN.<type>.csv`, `basename.<type>_statistics.csv`
- Pure-CG MD (`Molecule::isCGSystem()`, all atoms CG): timestep x10, PBC wrapping, VTF trajectory (`cg_write_vtf`)

## Status

- All non-external optimizers are 🤖 AI-generated and not human-tested; LBFGSpp is external code behind an AI-generated wrapper
- BMT output: 🤖 AI-generated, machine-tested, human production testing pending. `-opt` multi-XYZ: ⚙️ machine-tested (`cli_curcumaopt_07_opt_multixyz`)
- ANCOpt (port of xtb AncOpt): 🤖 AI-generated, ⚙️ machine-tested in manual runs (archive); no ctest passes `-optimizer`. Not tested: QM gradients, TS, linear molecules
- ANCOpt tiers: n3 > `anc_lanczos_threshold` (1800) uses Lanczos ANC capped at `anc_lanczos_k` (500), so L-BFGS-in-ANC (nvar > `anc_lbfgs_threshold` 2000) is unreachable by default

## Open items

- `CurcumaOpt` features without a driver counterpart: `opt_h` (first priority in the former note), Hessian after opt, `mo_scheme`, `./stop` file, GFN2 dipole,
  `fusion` (needs a `Molecule::Check()` gate, which `OptimizerDriver` does not have)
- Parallel `-opt` batch has no live progress bar: `CxxThreadPool` updates it only in legacy mode (work-package note OPT_MULTIXYZ_PARALLELISM_WP in docs/, not committed to git yet)
- ConfScan `Reorder`: a thread disabled for one candidate can stay disabled for the next (exclude-list `continue`, re-enable only inside
  `if (reorder && keep_molecule)`), and a skipped thread keeps its old `m_keep_molecule`; a fix changes accepted counts, needs its own review
- ConfScan `force_reorder` baseline kept one structure more than the tiered pipeline (44-structure ensemble), not root-caused ([CONFSCAN_REORDER_TIMING_WP](../../docs/CONFSCAN_REORDER_TIMING_WP.md))
- ConfSearch: Phase C `cluster`/`weighted` calibration experimental; wide-hill MTD blow-up open (roadmap items 3, 4)
- Pure CG: `cg_timestep_scaling`/`cg_timestep_factor` have no effect, `SimpleMD::Initialise()` hard-sets 10.0 after reading them (code reading, not run)
- CG ellipsoids/rotation not implemented (`m_cg_enable_rotation = false`); `Json2KeyWord` still in 10 built capability files (migrate to `ConfigManager`)

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Future development tasks and visions to be defined by operator/programmer*

---

Previous version (dated fixes, measurements, status logs, removed 2026-10-01): [docs/archive/CAPABILITIES_NOTES_2026-10.md](../../docs/archive/CAPABILITIES_NOTES_2026-10.md)
