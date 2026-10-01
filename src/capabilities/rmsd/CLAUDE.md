# CLAUDE.md - RMSD Module

Alignment and atom-reordering strategies of `RMSDDriver` (`../rmsd.h/.cpp`); the driver is used by `-rmsd`, ConfScan, ConfSearch,
rmsdtraj, SimpleMD, NEB docking and `OptimizerDriver`. Status: 🤖 AI-generated refactor (Jan 2025), coverage under Tests.

## Strategies (`rmsd_strategies.*`, factory `AlignmentStrategyFactory::createStrategy`)

| ID | `-rmsd.method` | Class | Approach (enum comment) |
|----|----------------|-------|-------------------------|
| 1 | `incr` | `IncrementalAlignmentStrategy` | legacy incremental, threaded |
| 2 | `template` | `TemplateAlignmentStrategy` | fragment template + Kuhn-Munkres |
| 3 | `hybrid0` | `HeavyTemplateStrategy` | heavy-atom template (legacy) |
| 4 | `subspace` (default), `hybrid` | `AtomTemplateStrategy` | selected-element template + Kuhn-Munkres |
| 5 | `inertia`, `free` | `InertiaAlignmentStrategy` | inertia frame + Kuhn-Munkres |
| 6 | `molalign` | `MolAlignStrategy` | external `molalign` binary (`molalign_bin`) |
| 7 | `dtemplate` | `DistanceTemplateStrategy` | distance template + Kuhn-Munkres |
| 10 | `predefined` | `PredefinedOrderStrategy` | given order, no reordering |

- Name map: `method_map` in `../rmsd.cpp`; unknown names warn and fall back to `subspace`
- Heavy, Atom and Distance templates call `IncrementalAlignmentStrategy` internally (`createStrategy(1)`)

## Other files

- `rmsd_functions.h`: `RMSDFunctions` Kabsch best fit (`BestFitRotation`, weighted variants), `getRMSD`
- `munkres.h` (local) and `../c_code/` (C Hungarian): assignment solvers behind `RMSDDriver::SolveCostMatrix()`
- `rmsd_costmatrix.*` (`CostMatrixCalculator`, 4 overloads, 6 cost functions) and `rmsd_assignment.*` (`MunkresAssignmentSolver`): compiled, no caller outside their own files

## Architecture and invariants

- `AlignmentStrategy::align(RMSDDriver*, const AlignmentConfig&)` returns `AlignmentResult`; strategies reach private driver state via `friend` (`../rmsd.h:102-109`)
- `AlignmentConfig` is defined in `../rmsd.h` and built by `CreateAlignmentConfig()`; `InitializeAlignmentStrategy()` runs inside `LoadControlJson()`
- `LoadControlJson()` delegates to `LoadFragmentAndThreading-`, `LoadAlignmentMethod-`, `LoadElementTemplate-`, `LoadCostMatrix-`, `LoadAtomSelection-`,
  `LoadFileOrderParameters()` and `DisplayConfigurationSummary()`
- Strategies build cost matrices with `RMSDDriver::MakeCostMatrix()` and solve them with `RMSDDriver::SolveCostMatrix()`; Kuhn-Munkres is O(n^3)

## Tests

- `cli_rmsd_01..06` run the default method only; ctests `confscan_subspace/free/dtemplate/template/molalign` reach strategies 4, 5, 7, 2, 6 (molalign only if the binary is found)
- `hybrid0` (3) and `predefined` (10) have no test found by grep; `rmsd_test` (`test_cases/CMakeLists.txt`) is built but not registered with `add_test`

## Open items

- Duplicate code: `CostMatrixCalculator`/`MunkresAssignmentSolver` were extracted from `MakeCostMatrix`/`SolveCostMatrix` but nothing calls them; wire in or remove (operator decision)
- Uncalled leftovers in `RMSDDriver`: `MolAlignLib()`, `PrepareHeavyTemplate()`, both `PrepareAtomTemplate()`, `PrepareDistanceTemplate()`; the strategies hold copies
- `-rmsd.method` help (`../rmsd.h`) lists `hungarian`, which is not in `method_map` (falls back to `subspace`), and omits the aliases `hybrid`, `hybrid0`, `free`

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Future development tasks and visions to be defined by operator/programmer*

---

Previous version (refactoring phases and completion claims, removed 2026-10-01): [docs/archive/RMSD_NOTES_2026-10.md](../../../docs/archive/RMSD_NOTES_2026-10.md)
