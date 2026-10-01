# CLAUDE.md - Curcuma Development Guide

## Overview

**Curcuma** is a molecular modelling and simulation toolkit for computational chemistry, force fields, and quantum calculations.

**Educational Purpose**: Teaching and research platform prioritizing pedagogical clarity over complex software engineering. Goal: learn, understand, and implement computational chemistry methods without getting lost in C++ abstractions.

## Very General Instructions for AI Coding
- Avoid flattery, compliments, or positive language. Be clear and concise. Do not use agreeable language to deceive.
- Do comprehensive verification before claiming completion
- Show me proof of completion, don’t just assert it
- Prioritize thoroughness over speed
- If I correct you, adapt your method for the rest of the task
- No completion claims until you can demonstrate zero remaining instances
- Dont use git -A to blindly add files

## AI-Generated Content and Validation Policy

### Status Labels — Definitions
These labels are used throughout CLAUDE.md files and documentation. Only the human operator may assign ✅ TESTED or ✅ APPROVED.

| Label | Meaning | Who sets it |
|-------|---------|-------------|
| 🤖 AI-generated | Code written by AI, not reviewed by human | AI |
| ⚙️ Machine-tested | Passes automated tests (CI, ctest) | AI |
| 👁️ Human-reviewed | Human has read and understood the code | Human only |
| ✅ TESTED | Human has run it on real problems and it behaves correctly | **Human only** |
| ✅ APPROVED | Human confirms correctness, ready for production | **Human only** |

**The AI must never write ✅ TESTED or ✅ APPROVED on its own work.**

### Conservative Self-Assessment Rules for AI
When documenting implemented features, the AI must apply these rules:

1. **Automated tests pass ≠ correct** — tests only cover what was anticipated. Unknown failure modes exist.
2. **Agreement with reference on test molecules ≠ general correctness** — the reference comparison is only as broad as the test set.
3. **No gaps visible ≠ no gaps exist** — absence of a known bug is not the same as correctness. Especially for AI-generated scientific code: the most dangerous bugs are those that produce plausible but wrong results.
4. **"Implemented" means the code compiles and runs** — it does not imply physical correctness, numerical stability across all inputs, or completeness relative to the reference method.
5. **When in doubt, add a caveat** — a caveat that turns out to be unnecessary is harmless. A missing caveat on wrong code causes user errors.

### Required Documentation for New AI-Generated Features
Every new method or capability added by AI must include in its CLAUDE.md:
- What was tested (which molecules, which conditions)
- What was **not** tested (system classes, edge cases, conditions)
- What is **not implemented** relative to the reference method
- A note that human production testing is pending until the human removes it

## General Instructions

- Each source code dir has a CLAUDE.md with basic information of the code and logic
- **Keep CLAUDE.md files FOCUSED and CONCISE** - ONE clear idea per bullet, max 1-2 lines
  - ❌ DON'T: Multi-paragraph explanations, code examples, historical details
  - ✅ DO: Brief statements with links to detailed docs if needed
  - ✅ DO: "✅ **Feature name** - Brief description" for completed items
- **Instructions blocks** contain operator-defined future tasks and visions for code development
- Only include information important for ALL subdirectories in main CLAUDE.md
- Preserve new knowledge from conversations but keep it brief
- Always suggest improvements to existing code
- **Keep entries concise and focused to save tokens**
- **Keep git commits concise and focused**
- **Rule of thumb**: If a CLAUDE.md section exceeds 20 lines, consider if it's better placed elsewhere
- Newly added features need a precise and short documentation under docs/, a link to the documentation from claude.md and a note in the readme

## Development Guidelines

### Code Organization
- Each `src/` subdirectory contains detailed CLAUDE.md documentation
- Variable sections updated regularly with short-term information
- Preserved sections contain permanent knowledge and patterns
- Instructions blocks contain operator-defined future tasks and visions

### Implementation Standards

#### Educational-First Design Principles
- **Core functionality visibility**: Always provide clear, direct access to the computational chemistry implementation
- **Minimal abstraction layers**: Avoid unnecessary templates, inheritance hierarchies, or design patterns that obscure the scientific content
- **Algorithm transparency**: The actual mathematical/physical implementation should be easily locatable and readable
- **Documentation focus**: Emphasize *what* the code does scientifically, not just *how* it's structured
- **Learning-oriented comments**: Include references to equations, papers, and theoretical background in code comments
- **Method implementation clarity**: Each computational method should have a clear entry point with minimal indirection
- **Accuracy:** 100 % with respect to referenz implementation for any scientific method

#### Code Organization for Learning
- **Flat over hierarchical**: Prefer simple, direct implementations over complex class hierarchies
- **Self-contained modules**: Each computational method should be understandable without deep knowledge of the entire system
- **Clear naming**: Function and variable names should reflect their scientific meaning
- **Minimal templates**: Only use templates when absolutely necessary; prefer explicit types for clarity
- **Direct implementations**: Avoid hiding core algorithms behind layers of abstractions

#### Standard Development Practices
- Mark new functions as "Claude Generated" for traceability
- Document new functions briefly (doxygen ready) with scientific context
- Document existing undocumented functions if appearing regularly (briefly and doxygen ready)
- Remove TODO Hashtags and text if done and approved
- Implement comprehensive error handling and logging
- Maintain backward compatibility where possible
- **Always check and consider instructions blocks** in relevant CLAUDE.md files before implementing
- Reformulate and clarify task and vision entries if not already marked as CLAUDE formatted
- In case of compiler warning for deprecated functions, replace the old function call with the new one
- Implement timing analysis for complex functions
- Keep track of significant improvements in AIChangelog.md, one line per fact (the file is an index; long entries live in docs/changelog/, open tasks in TODO.md, old TODO text in docs/TODO_ARCHIVE_2026-10.md)
- **Complex Architecture Documentation**: Factory patterns, dispatchers, and multi-step workflows require comprehensive inline documentation following docs/ARCHITECTURE_DOCUMENTATION.md standards
- **BMT output compatibility (MANDATORY)**: Every capability that writes output files MUST route them through `outputPath()` (CurcumaMethod subclasses) or `BMTUtils::outputPath()` (standalone handlers). Hardcoded CWD paths are not permitted. Verify with `-no_bmt` (legacy) and default BMT mode before merging.
- **No UTF symbols in terminal output**: Do not use Unicode box-drawing characters, emoji, arrows (->), checkmarks, or any non-ASCII symbols in fmt::print/std::cout output. Use plain ASCII only. Reason: breaks output in many terminal emulators, log files, and remote shells. CurcumaLogger colored output is exempt (uses ANSI codes, not Unicode).

#### Parameter Definition Standards (MANDATORY for new capabilities)
- **ALL new capabilities MUST use Parameter Registry System** - no static JSON configurations
- **Definition Location**: Define parameters in capability header using PARAM macros within `BEGIN_PARAMETER_DEFINITION(module)` block
- **Naming Convention**: Use **snake_case** exclusively (`max_iterations`, not `MaxIterations` or `maxIterations`)
- **Include Required**: Add `#include "src/core/parameter_macros.h"` to capability header
- **Help Text**: Provide comprehensive, user-facing descriptions for each parameter
- **Categories**: Group parameters logically (Basic, Algorithm, Output, Advanced)
- **Aliases**: Add old parameter names as aliases for backward compatibility during migration
- **Type Safety**: Use correct ParamType (String, Int, Double, Bool) matching C++ type
- **Constructor**: Use `ParameterRegistry::getInstance().getDefaultJson("module")` instead of static JSON
- **Build Verification**: Run `make GenerateParams` and check for validation warnings
- **Documentation**: See reference implementation in `src/capabilities/analysis.h`
- **Migration Guide**: Follow [docs/archive/PARAMETER_MIGRATION_GUIDE.md](docs/archive/PARAMETER_MIGRATION_GUIDE.md) for existing capabilities

#### Copyright and File Headers
- **Copyright ownership**: All copyright remains with Conrad Hübler as the project owner and AI instructor
- **Year updates**: Always update copyright year to current year when modifying files
- **Claude contributions**: Mark Claude-generated code sections but copyright stays with Conrad
- **Format**: `Copyright (C) 2019 - 2026 Conrad Hübler <Conrad.Huebler@gmx.net>`
- **AI acknowledgment**: Add Claude contribution notes in code comments, not copyright headers

#### Code Structure Guidelines
- **Main computational functions**: Should be easily findable and readable without deep C++ knowledge
- **Algorithm documentation**: Include mathematical formulations and literature references
- **Parameter transparency**: Make method parameters and their physical meaning obvious
- **Debugging accessibility**: Provide easy ways to inspect intermediate results and algorithm steps
- **Educational examples**: Include well-commented example usage in documentation

## Current Capabilities

Long form with all measurements: [docs/CAPABILITY_NOTES_2026.md](docs/CAPABILITY_NOTES_2026.md). Campaign history and
numbers: [docs/KNOWN_ISSUES_ARCHIVE.md](docs/KNOWN_ISSUES_ARCHIVE.md).

> **Default method is `gfnff` for every capability** (`-sp`, `-opt`, `-md`, `-hessian`, `-confsearch`, `-casino`): fast,
> general purpose. **`gfn2` is the accurate one** (native GFN2-xTB, ~100x slower, with electronic structure). `uff` is
> no capability's default; it stays the default only of `ForceFieldGenerator`'s own `method` parameter, which selects
> the force field that class builds (uff | uff-d3 | qmdff | cg).

> ⚠️ **All native QM/FF methods are 🤖 AI-implemented and ⚙️ machine-tested only. None is ✅ TESTED.** Validate against an
> external reference before research use.

| Method | Status / what was checked | Docs |
|---|---|---|
| `gfn2` native | 11/12 `sqm_reference` at 1e-8 vs tblite (231-atom `complex` open); GMTKN55 vs xtb 6.7.1 within 0.017 kcal/mol (D4 ATM cutoff split, `-xtb.d4_atm_cutoff 40.0` removes it); open shell, d-shell, MOR41 95/95 < 1e-6 Eh vs tblite | SQM_VALIDATION.md, SCF_MODES.md, SQM_DSHELL_WP.md |
| `gfn1` native | 14/16 `sqm_reference` at 1e-8 (He2, `complex` open); halogen-bond term; residual 0.011 kcal/mol is the xtb/tblite STO-6G split (`-xtb.sto6g_legacy_4sp true`) | SQM_WP2_gfn1_accuracy.md |
| `gfnff` native | port of GFN-FF; vs pprcht/gfnff over every structure: MOR41 and GMTKN55 max 0.07 kcal/mol (Known Issue #25); S30L-CI max 0.58, from the deliberate sTors deviation (0.05 with `-gfnff.storsion_reference_loop_bug true`, #12); vs xtb differs where xtb and pprcht differ (#6, #21) | GFNFF_STATUS.md, S30L_GFNNF_VALIDATION.md |
| `pm3`/`am1`/`mndo` native | 21/21 tests vs Ulysses (< 4 µEh) | MNDO_INTEGRALS.md |
| `eht` native | machine-tested | |
| External | TBLite (`tblite-*`, `ipea1`), XTB (`xtb-*`), Ulysses (`ugfn2`, PM6, RM1 ...), ORCA, DFT-D3/D4, `xtb-gfnff` | |
| Force fields | UFF, QMDFF, coarse-grained `cg` (`-load_ff_json`), parameter caching for all | CLEANUP_2026_09.md |

- **GPU** (`-gpu cuda|rocm|vulkan|auto`): runtime `dlopen` plugins, CPU-only runs need none. CUDA (run on H200, RTX 5090, RTX 5080, RTX A4500, GTX 1660) and ROCm (Radeon 890M gfx1150 only) for gfn1/gfn2/gfnff, Vulkan opt-in. Multi-GPU, tuning flags, ROCm build unverified: [SQM_GPU.md](docs/SQM_GPU.md), [MULTI_GPU.md](docs/MULTI_GPU.md), [GPU_TUNING.md](docs/GPU_TUNING.md), [GPU_PLUGIN_STARTUP.md](docs/GPU_PLUGIN_STARTUP.md). GPU fallbacks are counted; `-gpu_strict true` exits at the first.
- **SCF** (native xTB): Broyden default; `-scf_mode`, `-scf_guess h0|eeq|fragments`; NaN or unconverged SCF is an error (`-scf_allow_unconverged true`); `-scf_extrapolation aspc|gauss`; `-large_system_mode fragments|dc|sparse`; threads `-threads N`: [SCF_MODES.md](docs/SCF_MODES.md), [SQM_LARGE_SYSTEMS.md](docs/SQM_LARGE_SYSTEMS.md), [SQM_THREADING.md](docs/SQM_THREADING.md).
- **Solvation**: TBLite (CPCM/GB/ALPB), Ulysses GBSA, native gfn1/gfn2 ALPB+GBSA (≤1e-8 Eh vs tblite), native GFN-FF ALPB (≤1e-8 Eh vs xtb 6.7.1); native CPCM pending: [SOLVATION.md](docs/SOLVATION.md), [SQM_SOLVATION_WP.md](docs/SQM_SOLVATION_WP.md).
- **Optimisation**: LBFGS and others, constraints; default `auto` optimiser can abort on a high-symmetry system whose only active mode is totally symmetric (use `-opt.optimizer lbfgs`).
- **ConfSearch / ConfScan / RMSD**: ConfScan filter protocol for structures with topological symmetry (energy, rotational constants, Vietoris-Rips barcodes, permutation reuse): preprint [ChemRxiv 10.26434/chemrxiv.15009180/v1](https://chemrxiv.org/doi/full/10.26434/chemrxiv.15009180/v1); dual-method (`-md_method`, `-opt_method`), restartable (`-restart`), MTD screen: [CONFSEARCH_DUAL_METHOD.md](docs/CONFSEARCH_DUAL_METHOD.md), [CONFSEARCH_RESTART.md](docs/CONFSEARCH_RESTART.md), [CONFSEARCH_MTD_SCREEN.md](docs/CONFSEARCH_MTD_SCREEN.md).
- **MD**: SimpleMD, temperature ramps/regions ([TEMPERATURE_RAMP.md](docs/TEMPERATURE_RAMP.md)), step-rejecting integrator `-adaptive_step` (off by default), PLUMED `-mtd`. **MD times recorded before `ef462fcf` (Sep 2026) must be multiplied by 1.9516.** Large-system MD: [MD_LARGE_SYSTEMS.md](docs/MD_LARGE_SYSTEMS.md).
- **Other**: `-interaction` (E(AB)-E(A)-E(B), S30L modes), NEB, trajectory analysis (parallel, TrajectoryWriter), scattering, Hessian, persistent diagram, orbitals.
- **Output**: default BMT directory (Basename.Method.Timestamp), `-bak` copies files back, `-no_bmt` writes to CWD (see `src/tools/CLAUDE.md`).
- **Core**: MNDO integrals ([MNDO_INTEGRALS.md](docs/MNDO_INTEGRALS.md)); one force-field engine `FFWorkspace` (legacy `ForceField`/`ForceFieldThread` GFN-FF path removed Sep 2026, [CLEANUP_2026_09.md](docs/CLEANUP_2026_09.md)).

Older (2025) work: `AIChangelog.md` and git history.

## Architecture

### Core Components

#### Energy Calculator (`src/core/energycalculator.h/cpp`)
**COMPLETELY REFACTORED (January 2025)** - New unified polymorphic architecture:

##### **New Architecture**
- **Polymorphic Design**: Single `ComputationalMethod` interface for all QM/MM methods
- **MethodFactory**: Priority-based method resolution with hierarchical fallbacks
- **Unified Interface**: `calculateEnergy()`, `getGradient()`, consistent API across all methods
- **Method Priority System**: `gfn2`/`gfn1` → Native xTB (AP3); `ipea1` → TBLite; `ugfn2` → Ulysses
- **Thread-Safe**: Full multi-threading support maintained
- **Universal Verbosity**: Consistent output control across all computational methods

##### **Method Resolution (New)**
```cpp
// New MethodFactory system (replaces old SwitchMethod)
std::unique_ptr<ComputationalMethod> method = 
    MethodFactory::createMethod("gfn2", config);
double energy = method->calculateEnergy();
```

##### **Supported Method Hierarchies** (AP3, April 2026)
- **gfn1/gfn2**: Native curcuma xTB (canonical); `ipea1`/`ugfn2` for other providers
- **xtb-gfn1/xtb-gfn2**: External GFN — TBLite (USE_TBLITE) → XTB binary (USE_XTB), like `xtb-gfnff`
- **tblite-gfn1/tblite-gfn2**: TBLite explicitly (forces that backend)
- **eht**: Native only (always available, no dependencies)
- **pm3**: Native only (H, C, N, O supported, no dependencies)
- **uff/qmdff**: ForceField wrapper with parameter generation
- **gfnff**: Native C++ GFN-FF (always available, ✅ **COMPLETE**)
- **xtb-gfnff**: Fortran/XTB GFN-FF — ExternalGFNFF (USE_GFNFF) → XTB (USE_XTB)

#### Force Field System (`src/core/forcefield.h/cpp`)
Modern force field engine with:
- **Threading** via the single `FFWorkspace` engine (the legacy `ForceFieldThread` engine was removed Sep 2026)
- **Universal parameter caching** - automatic save/load as JSON
- **Method-aware loading** - validates parameter compatibility
- **Multi-threading safety** - controllable caching for concurrent calculations

#### QM Interface (`src/core/energy_calculators/qm_methods/interface/abstract_interface.h`)
Unified interface for all quantum mechanical methods:
```cpp
class QMInterface {
    virtual bool InitialiseMolecule() = 0;
    virtual double Calculation(bool gradient, bool verbose) = 0;
    virtual bool hasGradient() const = 0;
    virtual Vector Charges() const = 0;
    virtual Vector BondOrders() const = 0;
};
```

### File Organization

```
curcuma/
├── src/
│   ├── capabilities/          # High-level molecular modeling tasks
│   │   ├── confscan.cpp      # Conformational scanning
│   │   ├── confsearch.cpp    # Conformational searching  
│   │   ├── optimizer_driver.cpp # Geometry optimization (OptimizerDriver)
│   │   ├── simplemd.cpp      # Molecular dynamics
│   │   └── rmsd.cpp          # Structure analysis
│   ├── core/                 # Core computational engines
│   │   ├── energycalculator.cpp      # NEW: Unified polymorphic dispatcher
│   │   ├── molecule.cpp              # Molecular data structures
│   │   ├── curcuma_logger.cpp        # Universal logging system
│   │   ├── energy_calculators/       # NEW: All computational methods organized here
│   │   │   ├── computational_method.h     # Base interface for all methods
│   │   │   ├── method_factory.cpp         # Priority-based method creation
│   │   │   ├── qm_methods/                # QM method implementations & wrappers
│   │   │   │   ├── eht.cpp                # Extended Hückel Theory + verbosity
│   │   │   │   ├── xtbinterface.cpp       # XTB interface + verbosity
│   │   │   │   ├── tbliteinterface.cpp    # TBLite interface + verbosity
│   │   │   │   ├── ulyssesinterface.cpp   # Ulysses interface + verbosity
│   │   │   │   ├── gfnff_method.cpp         # ComputationalMethod wrapper
│   │   │   │   ├── orcainterface.cpp      # ORCA interface
│   │   │   │   ├── dftd3interface.cpp     # DFT-D3 dispersion corrections
│   │   │   │   ├── dftd4interface.cpp     # DFT-D4 dispersion corrections
│   │   │   │   ├── *_method.cpp           # Polymorphic method wrappers
│   │   │   │   └── interface/             # Abstract interfaces
│   │   │   └── ff_methods/                # Force field implementations
│   │   │       ├── forcefield.cpp         # Force field engine + verbosity  
│   │   │       ├── forcefieldgenerator.cpp # Parameter generation + verbosity
│   │   │       ├── ff_workspace*.cpp      # Partitioned energy/gradient engine (all FF methods)
│   │   │       ├── gfnff_method.cpp      # Native GFN-FF implementation
│   │   │       ├── gfnff.h               # GFN-FF class interface
│   │   │       ├── gfnff_inversions.cpp   # GFN-FF inversion terms
│   │   │       ├── gfnff_torsions.cpp     # GFN-FF torsion terms
│   │   │       ├── qmdff.cpp              # QMDFF implementation
│   │   │       ├── eigen_uff.cpp          # UFF implementation
│   │   │       └── *_par.h                # Parameter databases
│   ├── tools/                # Utilities and file I/O
│   │   ├── formats.h         # File format handling (XYZ, MOL2, SDF)
│   │   └── geometry.h        # Geometric calculations
│   └── helpers/              # Development and testing tools
├── test_cases/               # Validation and benchmark molecules
├── external/                 # Third-party dependencies
└── CMakeLists.txt           # Build configuration
```

## Past Developments

Moved to [docs/CAPABILITY_NOTES_2026.md](docs/CAPABILITY_NOTES_2026.md) (GFN-FF/GFN1/GFN2 cleanup Sep 2026, `-interaction`, GFN-FF ring torsions, GPU HB freeze) and `AIChangelog.md`.

## Build and Test Commands

**ALWAYS build and test in the `release/` directory.** This is the canonical build for regression comparisons. The `build/` directory is for development experiments only.

```bash
# Build — always use release/
cd release
make -j4

# Run all tests
ctest --output-on-failure

# Run specific test categories
ctest -R "cli_rmsd_" --output-on-failure      # RMSD CLI tests (6/6 passing)
ctest -R "cli_confscan_" --output-on-failure  # ConfScan CLI tests (7/7 passing)
ctest -R "cli_simplemd_" --output-on-failure  # SimpleMD CLI tests (7/7 passing)
ctest -R "cli_curcumaopt_" --output-on-failure # Opt CLI tests (6/6 passing)

# Run individual CLI test with verbose output
ctest -R "cli_rmsd_01" --verbose

# Legacy: Manual curcuma execution
./curcuma -sp input.xyz -method uff           # UFF single point
./curcuma -rmsd ref.xyz target.xyz            # RMSD calculation
./curcuma -opt input.xyz -method gfn2         # GFN2 optimization
```

**Test count**: `ctest -N` in `release/` lists 333 tests (Oct 1, 2026); run `ctest` for the current result. The per-category counts above are historical.

## Project Management

- **Knowledge store (Obsidian vault)**: `~/Nextcloud/Obsidan/Wissen/`, German, operator-owned. Read its `00 Regeln.md` before writing anything into it.
  - `Projekte/curcuma *.md`: project status; `Labor/curcuma *.md`: lab journal for computational campaigns (append-only, mandatory fields incl. commit, dirty state, diff); `Offene Fragen/`: scientific questions without an owner; `Wissen/curcuma.md` and related notes: reusable method knowledge.
  - Code, changelog and bug history stay in this repository; method knowledge, campaign records and project status go to the vault. Cross-project agent rules: `~/.claude/CLAUDE.md`.
- **Prioritized TODO list**: [TODO.md](TODO.md). **Module docs**: each `src/` subdirectory has a CLAUDE.md.
- **Status and validation**: [docs/GFNFF_STATUS.md](docs/GFNFF_STATUS.md), [docs/SQM_VALIDATION.md](docs/SQM_VALIDATION.md), [docs/GMTKN55_VALIDATION.md](docs/GMTKN55_VALIDATION.md), [docs/MOR41_VALIDATION.md](docs/MOR41_VALIDATION.md), [docs/S30L_GFNNF_VALIDATION.md](docs/S30L_GFNNF_VALIDATION.md), [docs/GRADIENT_VALIDATION.md](docs/GRADIENT_VALIDATION.md); technical debt: [docs/TECHNICAL_DEBT.md](docs/TECHNICAL_DEBT.md).
- **rev-gfnff backlog**: places where GFN-FF could be better than the reference, each checked against an external reference (r2SCAN-3c, GFN2, DLPNO, experiment): [docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md). Port fidelity stays the default.
- **Regression check before/after a change**: `scripts/refset_regression.py` runs MOR41/GMTKN55 through two curcuma binaries (reference built from a worktree) and reports which energies moved; `scripts/scan_convergence.py` separates real changes from non-converged SCF. See [docs/REGRESSION_CHECK.md](docs/REGRESSION_CHECK.md).
- **Benchmark sets**: `scripts/fetch_testset.py` fetches MOR41 and GMTKN55 on demand (S30L is manual/paywalled); `scripts/gmtkn55_compare.py` caches energies in `_run/energies.json` and reuses them unless `--recompute` is passed (see Validation Traps); `scripts/testset_perf.py` benchmarks wall-clock. See [docs/TESTSET_RETRIEVAL.md](docs/TESTSET_RETRIEVAL.md), [docs/GMTKN55_VALIDATION.md](docs/GMTKN55_VALIDATION.md).
- **S30L-CI** (manually supplied, not in `fetch_testset.py`): 30 host-guest complexes with explicit counterions, `test_cases/s30lci_test_set/` (per-structure dirs gitignored, `reference_s30lci`/`README` tracked), run via `scripts/s30lci_gfnff_compare.py`.

## Workflow States
- **ADD**: Features to be added
- **WIP**: Currently being worked on
- **ADDED**: Basically implemented
- **TESTED**: Works (by operator feedback)
- **APPROVED**: Move to changelog, remove from CLAUDE.md

## Where Things Go

One home per kind of information; copies drift apart. `scripts/check_docs.py` enforces the mechanical parts: run it
before committing, or enable the hook with `git config core.hooksPath scripts/git-hooks`.

| Kind | Home |
|---|---|
| Rules, invariants, traps, layout of a directory | the `CLAUDE.md` of that directory (root 500 lines, others 120, no line over 600 characters) |
| Open code tasks and defects | `TODO.md`, open items only, at most 3 lines each |
| Bug report, measurement, campaign | one dated document in `docs/` plus one index line in `AIChangelog.md` |
| Method status and validation | `docs/<METHOD>_STATUS.md` or `docs/*_VALIDATION.md` |
| Test structures (geometry, charge, spin, provenance) | `test_cases/structures/` (manifest, naming and provenance rules in its README); no structure file anywhere else |
| Scientific questions without an owner, project status | the vault (`Offene Fragen/`, `Projekte/`); `TODO.md` links to the note by name, no copy |
| History | git log and the `AIChangelog.md` index |

- No dated status sections, "Completed" lists or fix narratives in a CLAUDE.md. No new markdown file in the repository root (README, CLAUDE, TODO, AIChangelog only) and no completion or summary reports; the changelog line is the record.
- A statement about the code ("X exists", "default is Y") is checked against the code when it is written; a path in backticks must exist. What was not checked says so.
- Closing a task is one change: remove its TODO entry, add the changelog line, close or strike the linked vault note.
- Dead code: code without a caller is deleted in the same change that removes its last caller (git is the archive); no `.backup`/`.orig` files. Code that stays off on purpose carries the comment `INACTIVE (date): reason, switch` and a TODO entry. Known candidates: `TODO.md`, section "Löschkandidaten".
- Open questions: code-bound ones in `TODO.md`, section "Entscheidungen und Fragen"; scientific ones in the vault. Each side holds links, not copies.

## Git Best Practices
- **Only commit source files**: Use `git add <file>` for specific files, never `git add -A` without review
- **Review before committing**: Always check `git diff` and `git status` to avoid accidental commits
- **Build before commit**: Ensure `make -j4` succeeds and no compiler warnings/errors exist
- **Commit message format**: Start with action verb (Fix, Add, Improve, Refactor), follow with brief description
- **Include Co-Author info**: All commits include Claude contribution notes with proper attribution
- **Test artifacts stay local**: Build outputs and temporary test files are ignored by .gitignore
- **Branch names**: `feature/<topic>` for new capabilities, `fix/<topic>` for bug fixes; `<topic>` is 2-4 lowercase ASCII words in kebab-case naming the subject, no dates or issue numbers (e.g. `feature/gfnff-solvation`, `fix/bmt-dir-collision`)
- **Why the scheme matters**: `.github/workflows/ccpp.yml` builds `feature/**` and `fix/**` automatically and publishes each as its own `ci-<branch>` prerelease — a branch outside the scheme has to be added to the workflow by hand

## Standards

### Universal Logging System (`src/core/curcuma_logger.h/.cpp`)
**✅ FULLY IMPLEMENTED** across all computational methods

#### **Verbosity Levels**
- **Level 0**: **Silent Mode** - Zero output (critical for optimization/MD)
- **Level 1**: **Minimal Results** - Final energies, convergence status
- **Level 2**: **Scientific Analysis** - HOMO/LUMO, energy decomposition, molecular properties  
- **Level 3**: **Complete Debug** - Full orbital listings, timing, algorithm details

#### **Color Scheme & Functions**
- `CurcumaLogger::error()` - Always visible (red)
- `CurcumaLogger::warn()` - Level ≥1 (orange)
- `CurcumaLogger::success()` - Level ≥1 (green)
- `CurcumaLogger::result()` - Level ≥1 (white, neutral reporting of scientific results)
- `CurcumaLogger::info()` - Level ≥2 (default)
- `CurcumaLogger::param()` - Level ≥2 (blue, structured output)
- `CurcumaLogger::energy_abs()` - Energy output with units
- `CurcumaLogger::citation()` - registers a reference, printed at program end at any level (see docs/LOGGING_SYSTEM.md)

#### **Implementation Status**
- **QM Methods**: EHT, XTB, TBLite, Ulysses ✅
- **Force Fields**: ForceField, ForceFieldGenerator ✅  
- **Native Libraries**: XTB/TBLite verbosity synchronized ✅
- **Thread Safety**: Zero overhead at Level 0 ✅

### Unit System (`src/core/units.h`)
- **Centralized**: All constants in `CurcumaUnit` namespace with CODATA-2018 values
- **Internal**: Atomic units (Hartree, Bohr, atomic time) 
- **Output**: Auto-select user-friendly units (kJ/mol, Å, fs)
- **Educational**: Clear naming and comprehensive documentation
- **Migration**: Replace scattered constants with centralized functions

### JSON Controller System (`src/main.cpp`, CLI2Json)
- **Consistent Parameter Passing**: All methods use `controller["methodname"]` subdocuments
- **Examples**: `Hessian(controller["hessian"])`, `QMDFFFit(controller["qmdfffit"])`, `ModernOptimizer(..., controller["opt"])`
- **CLI Arguments**: Automatically split into controller subdocuments via `CLI2Json()`
- **Structure**: `controller[keyword][parameter]` - e.g. `controller["opt"]["verbosity"]`, `controller["hessian"]["MaxIter"]`
- **Global Parameters**: `verbosity`, `threads`, `method`, `gpu` are additionally duplicated at top-level
- **Flat-flag auto-routing (2026)**: Any registered PARAM is reachable by its flat name (`-cn_cutoff_bohr 5.5` routes to `controller["gfnff"]["cn_cutoff_bohr"]` because the registry records ownership). Same-name in the active command's module wins; truly ambiguous names (multiple owners, none matching) warn and stay in the command module. Dotted form `-<module>.<param>` always works for disambiguation. Unregistered/legacy flags stay in the command module (unchanged).
- **JSON round-trip (2026)**: `-export_run file.json` writes the resolved controller plus `_command`, `_input`, and full registry defaults for every touched module. `-import_config file.json` performs a recursive deep merge (CLI wins at every depth). Invoking `curcuma -import_config run.json` reads `_command`/`_input` from the JSON, so the file alone is enough to replay a run. See [docs/CLI_ROUND_TRIP.md](docs/CLI_ROUND_TRIP.md).

## Planned Development

### TRAJECTORY ANALYSIS CONSOLIDATION
**Status**: ✅ Phases 1-3 complete (Jan 2026) — see [docs/ANALYSIS_CONSOLIDATION_PLAN.md](docs/ANALYSIS_CONSOLIDATION_PLAN.md)
- ✅ TrajectoryWriter, analysis.cpp migration, TrajectoryStatistics extended
- ⏳ Phase 4: Migrate `trajectoryanalysis.cpp` + `rmsdtraj.cpp` (optional)
- ⏳ Phase 5-6: Cleanup geometry commands, ProgressTracker (optional)

---

### Breaking Changes (Test-Driven)
- **Molecule data structure refactoring**: Hybrid SOA/AOS design for better performance
  - **PHASE 1**: ✅ Comprehensive test suite with refactoring-specific validation
    - ctest `molecule_comprehensive` (target `test_molecule`, source `test_cases/test_molecule.cpp`): 15 test categories
    - `src/core/REFACTORING_ROADMAP.md`: Detailed phase-by-phase plan
    - Tests include current behavior AND validation for planned improvements
    - Specific tests for: XYZ parser unification, cache granularity, fragment O(1) lookup, type safety
  - **PHASE 2**: XYZ Comment Parser unification (eliminate 10 duplicate functions)
    - **CRITICAL**: Production comment formats must not break (ORCA, XTB, simple energy)
    - See `src/core/XYZ_COMMENT_FORMATS.md` for required format compatibility
  - **PHASE 3**: Granular cache system (replace single m_dirty flag)
  - **PHASE 4**: Fragment system O(1) lookups (replace std::map)
  - **PHASE 5**: Type-safe ElementType enum (replace int elements)
  - **PHASE 6**: Unified atom structure with zero-copy geometry access
  - **CRITICAL**: All existing functionality must remain API-compatible

## Validation Method and Traps

Learned during the Jul - Sep 2026 GFN-FF / xTB validation (evidence per item in
[docs/KNOWN_ISSUES_ARCHIVE.md](docs/KNOWN_ISSUES_ARCHIVE.md)). Reference sets and runners: `scripts/refset_regression.py`,
`scripts/gmtkn55_compare.py`, `scripts/mor41_validation.py`, `scripts/s30lci_gfnff_compare.py`, `scripts/gradient_compare.py`.

- **Arbitrate with a third opinion.** For GFN-FF run pprcht/gfnff (`external/gfnff`, `cmake -Dbuild_exe=ON`): curcuma != pprcht == xtb is a curcuma bug; curcuma == pprcht != xtb is a reference split (#13, #21). tblite, not xtb, is the reference for gfn1/gfn2 (#27, #30).
- **Term diff before theory.** Compare bond/angle/torsion/Coulomb/dispersion/HB/XB separately (`-verbosity 2`); a constant error per structural motif, doubling with the motif, is one wrong parameter (#12, #16).
- **Characterise a set by every structure, never a sample** (#23, #24).
- **Caches fake results.** `gmtkn55_compare.py` reuses `_run/energies.json` (drop the keys or `--recompute`); GFN-FF writes `<basename>.topo.json` next to the input, so delete it before re-measuring and use a fresh directory per structure when file names repeat (#11, #28). Cache keys must carry every input the content depends on (#21c).
- **Check that the build happened.** Test `make`'s exit status, not a grep for "error"; a changed default argument in a header may not trigger a rebuild (#15, #30).
- **Hardcoded fallbacks are the defaults** for `GFNFF::GFNFF(const json&)`: changing a PARAM default in `gfnff.h` alone has no runtime effect (#31).
- **Element-indexed tables: check the length first.** Three 15-entry-short tables were undefined behaviour for Z >= 72 (#22).
- **No numerical path may depend on the print level** (#38) or on atom order (#34); per-structure benchmarks cannot see reused-calculator bugs (#32).
- **Gradients**: validate against a central finite difference of curcuma's own energy with Angstrom displacements (`gradient_unit_contract`); the `getGradient()` contract is Eh/Angstrom (#28). A converged SCF criterion must cover everything the mixer mixes (#29).
- **Deliberate deviations from a reference** are opt-in/out switches, listed with defaults in the archive's switch table (`-xtb.sto6g_legacy_4sp`, `-xtb.d4_atm_cutoff`, `-gfnff.storsion_reference_loop_bug`, `-gfnff.frag_charge_model reference`, ...).

## Open Items

Complete table with sources: [docs/KNOWN_ISSUES_ARCHIVE.md](docs/KNOWN_ISSUES_ARCHIVE.md#open-items). New bug reports go to the archive
and `AIChangelog.md`, not into this file; only open items and traps stay here.

- GFN-FF energy depends on atom numbering in pprcht and curcuma (911/2557 GMTKN55+MOR41 structures; REV_GFNFF_TODO #12) (#34)
- ROCm GFN-FF: no Coulomb term by default (F-1), `-scf_mixed_precision false` ignored (G2-13), wave32 assumption (F-17) (#37)
- Bit-identical results only at a fixed CPU thread count; GPU not run-to-run reproducible (#36)
- `shouldUpdateHBXB()` RMSD trigger off by `sqrt(natoms)` (ported from the reference, left); D4 ATM C6 frozen at setup (#33)
- `DC13/c20bowl` Coulomb -6.29 kcal/mol, not root-caused (#14); gfn1/gfn2 gradients vs xtb at exactly coincident coordinates (#28)
- Charge placement (#31) untested for periodic systems, GPU runtime, many fragments; `-DUSE_PORTABLE_MATH` unconfirmed on Windows (#7)
- Verbosity cannot be scoped across `CxxThreadPool` workers (#3); ConfSearch wide-hill MTD blow-up ([CONFSEARCH_ROADMAP.md](docs/CONFSEARCH_ROADMAP.md))
- GFN-FF is a force field: 60-70 kcal/mol from DLPNO on MOR41 reactions (#6); accuracy work in [docs/REV_GFNFF_TODO.md](docs/REV_GFNFF_TODO.md)
