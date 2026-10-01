# CLAUDE.md - Force Field Methods Directory

UFF, UFF-D3, QMDFF, coarse-grained `cg` and the native GFN-FF, all on one energy/gradient engine.
Status: 🤖 AI-implemented, ⚙️ machine-tested, no human production testing (see root `CLAUDE.md`, validation in
[docs/GFNFF_STATUS.md](../../../../docs/GFNFF_STATUS.md)). The former long version of this file (gradient status of Feb 2026,
D3/D4 fix records, GPU phase log, TODO checkboxes) is kept in
[docs/archive/FF_METHODS_NOTES_2026-10.md](../../../../docs/archive/FF_METHODS_NOTES_2026-10.md), with a list of the statements in it that are out of date.

## Architecture (one engine)

- **`FFWorkspace`** (`ff_workspace.h/.cpp`, `ff_workspace_gfnff.cpp`, `ff_workspace_uff.cpp`, `ff_workspace_cg.cpp`): the only energy/gradient engine. Partitions every interaction list over `CxxThreadPool` workers, reduces, then applies the CN chain rule and the Coulomb self-energy in `postProcess()`. Kernels are `calc*` methods (bonds, angles, dihedrals, inversions, STorsions, dispersion, repulsion, Coulomb, HB/XB, ATM, BATM, UFF/QMDFF, CG pairs). Input is a `GFNFFParameterSet` moved in with `setInteractionLists()`.
- **`ff_terms.h`**: term structs shared by the workspace, `ForceField`, `GFNFF` and the GPU SoA headers.
- **`ForceField`** (`forcefield.cpp/h`): UFF, UFF-D3, QMDFF, CG. Parameter generation in `ForceFieldGenerator`, JSON parameter caching, then hands the lists to its own `FFWorkspace`. `-method cg` takes `-load_ff_json FILE` (`cg_default`, `cg_per_atom`, `pair_interactions`, `bonds`); sphere gradient analytic, ellipsoids by finite differences.
- **`GFNFF`** (`gfnff_method.cpp/h`, `gfnff_torsions.cpp`, `gfnff_inversions.cpp`, `gfnff_frag_charge.cpp`): topology, Phase-1 EEQ, FT-HMO pi-bond orders, native parameter generation (`generateGFNFFParameterSet()`), per-step CN and Phase-2 EEQ (`prepareCNAndEEQ()`), H/X-bond re-detection, ALPB, then `m_workspace->calculate()`. Owns the `CxxThreadPool`. `NumGrad()` and `NumGradFixedCharges()` run on the workspace.
- **`EEQSolver`** (`eeq_solver.cpp/h`): Phase 1 topological, Phase 2 geometric (dxi/dgam/alpha corrections). Default Schur-Cholesky; projected PCG (`solveWithProjectedPCG`) from `eeq_ppcg_min_atoms` / `eeq_ppcg_min_nfrag` (0 = exact). The cached-factor mode is off (`eeq_refactor_eps_bohr` 0).
- **No periodic boundary conditions** on the CPU path (documented, not re-checked); the GPU path has its own PBC handling.
- **GPU** (`cuda/`, `rocm/`): one host wrapper template `GFNFFGpuMethodImpl<Backend>` in `qm_methods/gfnff_gpu_method_impl.h`; headers shared by CUDA and HIP through `gpu_rt.h`, `rocm/*_hip.h` are shims. ROCm was run on one Radeon 890M; the Sep 2026 header unification of the HIP side was not compiled (no ROCm SDK on the machine that made it). Details: [docs/GPU_TUNING.md](../../../../docs/GPU_TUNING.md), [docs/MULTI_GPU.md](../../../../docs/MULTI_GPU.md).

## Adding a GFN-FF term (workspace engine)

1. Struct in `gfnff_parameters.h` (or `ff_terms.h` if shared with UFF/QMDFF) plus a vector in `GFNFFParameterSet`.
2. `generateXxxNative()` in `GFNFF`, called from `generateGFNFFParameterSet()`.
3. Store it in `FFWorkspace::setInteractionLists()`, add its range to the partition struct in `partition()`, and a `FFEnergyComponents` / `FFTermTimings` field.
4. `FFWorkspace::calcXxx(int partition)` in `ff_workspace_gfnff.cpp` (energy into `acc.energy`, gradient into `acc.gradient`, CN chain rule into `acc.dEdcn`), call it from `executeGFNFF()`, add the reduction in `reduce()`.
5. Extend `GFNFFEnergyReport` and the verbosity-2 table in `GFNFF::Calculation()`.
6. Mirror the kernel in `cuda/gfnff_kernels.cu` (and `rocm/gfnff_rocm.hip`, SoA upload in `cuda/gfnff_soa.h`), or gate the term CPU-only.
7. Add the term to `test_cases/test_gfnff_validation.cpp` against the Fortran reference.

## Invariants and traps

Each item was a real defect; evidence in [docs/KNOWN_ISSUES_ARCHIVE.md](../../../../docs/KNOWN_ISSUES_ARCHIVE.md).
- **Units**: `GFNFF` works in Bohr; `ComputationalMethod::getGradient()` is Eh/Angstrom (#28). Constants come from `src/core/units.h`.
- **Port fidelity first**: the Fortran sources (`external/gfnff`, pprcht) decide; every deliberate deviation is an opt-in or opt-out switch. Several reference arrays have two variants (pre- and post-Hueckel pi membership, `imetal` versus the raw metal flag, `bpair` from `nbondmat` versus a BFS distance): use the one the reference uses (#15, #21, #23, #24, #25).
- **Energy-only calls on a reused instance** must see the current CN and D4 C6 (#32). Pair lists that depend on geometry (non-bonded repulsion, D4 pairs, HB/XB, explicit Coulomb) must be refreshed during MD (#33).
- **Defaults live in code, not only in PARAMs**: `GFNFF::GFNFF(const json&)` reads `m_parameters.value(key, FALLBACK)`; change the fallback as well as the PARAM (#31).
- **Topology cache** (`<basename>.topo.json`): delete it before re-measuring after any topology or charge change.
- **No numerical path may depend on the verbosity level** (#38) or on thread arrival order (#36).
- **Element-indexed tables**: check the length against the highest Z used (#22).

## Dispersion

- Shared D3 kernel (`d3param_generator.cpp`): validated only through the native GFN1 path against tblite (10 of 12 `sqm_reference` molecules at 1e-8, 2026-05-31). `uff-d3` and standalone D3 are not validated.
- GFN-FF dispersion follows the reference: D4 C6 from Casimir-Polder integration with CN-only Gaussian weights and the modified BJ damping of `gfnff_gdisp0.f90`, not standard D3/D4 (`../dispersion/d4param_generator.cpp`). The three-body ATM term is not in the reference and is off by default (`dispersion_atm`).
- `uff-d3`: `ForceFieldGenerator::GenerateUFFD3Parameters()` merges UFF terms with the D3 pair list.
- A shared CN utility exists (`cn_calculator.h`); `D3ParameterGenerator` and `EEQSolver` still have their own CN routines.

## Open items

- Not implemented relative to the reference (notes of Apr 2026, not re-checked): hyperconjugation torsion modulation, metal-specific C6.
- EEQ on GPU: FMM matrix-vector product (O(N log N)) and ROCm block-Jacobi mirror are open; ROCm many-fragment systems go to the exact CPU PCG above `eeq_rocm_cpu_fragment_threshold`.
- `k_dispersion` cannot overlap with the EEQ solve in gradient mode (needs `dc6dcn` from the post-EEQ charges).
- ROCm GFN-FF correctness gaps and GFN-FF open items: root `CLAUDE.md` "Open Items", [docs/MULTI_GPU_GAPS.md](../../../../docs/MULTI_GPU_GAPS.md), [docs/REV_GFNFF_TODO.md](../../../../docs/REV_GFNFF_TODO.md).

## References

FFWorkspace implements the formulas of Fortran `gfnff_engrad.F90`; parameter generation follows Spicher and Grimme, Angew. Chem. Int. Ed. 59 (2020) 15665; D3: Grimme et al., J. Chem. Phys. 132, 154104 (2010).
