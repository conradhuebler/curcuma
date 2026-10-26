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
- **`EEQSolver`** (`eeq_solver.cpp/h`): Phase 1 topological, Phase 2 geometric (dxi/dgam/alpha corrections). Default `auto` (`solve_method`): Schur-Cholesky below the thresholds, projected PCG (`solveWithProjectedPCG`) from `eeq_ppcg_min_atoms` / `eeq_ppcg_min_nfrag` (0 = exact), `localinv` preconditioner. The cached-factor mode is off (`eeq_refactor_eps_bohr` 0). Its PARAMs are in the `gfnff` block of `gfnff.h` (category EEQ); it still reads the JSON scope `eeq_solver`, so `-eeq_solver.<name>` keeps working and wins over `-gfnff.<name>`.
- **No periodic boundary conditions** on the CPU path (documented, not re-checked); the GPU path has its own PBC handling.
- **GPU** (`cuda/`, `rocm/`): one host wrapper template `GFNFFGpuMethodImpl<Backend>` in `qm_methods/gfnff_gpu_method_impl.h`; headers shared by CUDA and HIP through `gpu_rt.h`, `rocm/*_hip.h` are shims. ROCm was run on a Radeon 890M (gfx1150) and an RX 9060 XT (gfx1200, Oct 2026: gfnff energies equal CPU at 8 decimals, 24-7320 atoms); the Sep 2026 header unification of the HIP side was not compiled (no ROCm SDK on the machine that made it). Details: [docs/GPU_TUNING.md](../../../../docs/GPU_TUNING.md), [docs/MULTI_GPU.md](../../../../docs/MULTI_GPU.md).

## Adding a GFN-FF term (workspace engine)

1. Struct in `gfnff_parameters.h` (or `ff_terms.h` if shared with UFF/QMDFF) plus a vector in `GFNFFParameterSet`.
2. `generateXxxNative()` in `GFNFF`, called from `generateGFNFFParameterSet()`.
3. Store it in `FFWorkspace::setInteractionLists()`, add its range to the partition struct in `partition()`, and a `FFEnergyComponents` / `FFTermTimings` field.
4. `FFWorkspace::calcXxx(int partition)` in `ff_workspace_gfnff.cpp` (energy into `acc.energy`, gradient into `acc.gradient`, CN chain rule into `acc.dEdcn`), call it from `executeGFNFF()`, add the reduction in `reduce()`.
5. Extend `GFNFFEnergyReport` and the verbosity-2 table in `GFNFF::Calculation()`.
6. Mirror the kernel in `cuda/gfnff_kernels.cu` (and `rocm/gfnff_rocm.hip`, SoA upload in `cuda/gfnff_soa.h`), or gate the term CPU-only.
7. Add the term to `test_cases/unit/test_gfnff_validation.cpp` against the Fortran reference.

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
- ROCm (Oct 2026): FP64 `atomicAdd` is a compare-and-swap loop on gfx1200, so pair kernels gather per atom instead (dispersion, Coulomb tiles, non-bonded repulsion `CURCUMA_ROCM_REP_GATHER=0` restores the old path, HB per-triple slots + `k_slot_gather` and triples grouped by case type, `CURCUMA_ROCM_HB_GATHER=0` restores the old path). HB/XB refresh on ROCm is lean (`setSkipHostHBWorkspaceUpdate`: no copy into the CPU workspace, parallel pack and slot sort); the repulsion lists go through one pinned block with a parallel CSR build (`gpu_staged_list_upload`).
- ROCm EEQ PPCG matvec runs in FP32 (`rocblas_ssymv` on an FP32 matrix copy) with FP64 vectors and true-residual updates; `CURCUMA_ROCM_EEQ_FP32=0` restores FP64 (`CURCUMA_ROCM_EEQ_FP32_DELTA` sets the refresh ratio).
- ROCm below 500 atoms runs single-stream with lean step bookkeeping, a one-workgroup device EEQ for N <= 128 and <= 8 fragments, Coulomb tile shape by size (`CURCUMA_ROCM_SINGLE_STREAM`, `_LEAN_STEP`, `_EEQ_SMALL`, `_COUL_TILE`); `has_stream_overlap` is true on ROCm (early dispersion only).
- GPU (Oct 2026): `-gfnff.gpu_hb_detect_on_device` (default true; trait `has_device_hb_detect`, CUDA and ROCm) builds the H-bond triples on the device, `gpu_hb_detect_verify` compares them with the host list; ROCm orders them by rows (case, A) plus a per-row sort, no key sort library. Host inputs: `GFNFF::prepareHBDetectInputs()` is a copy of the O(N) part of `detectHydrogenBondsNative()` and must stay in step with it ([docs/GPU_STEP_PAUSES_2026-10-14.md](../../../../docs/GPU_STEP_PAUSES_2026-10-14.md)).
- ROCm (WP9, Oct 2026): `gpu_rep_pairs_on_device` also on ROCm (ordered device rebuild with the gather CSR, host fallback `CURCUMA_ROCM_FAIL_REPDEV`), `rocm_dc6dcn_tables` (CUDA tabulation, default off by measurement), `rocm_device_schur` (Schur step of the dense EEQ on the device, default on); notes in [docs/ROCM_RESIDENT_LOOP_WP.md](../../../../docs/ROCM_RESIDENT_LOOP_WP.md).
- ROCm (WP10, Oct 2026): `rocm_fixed_point` (default true) accumulates forces, dE/dCN, energies and CN sums as fixed-point integers (CUDA scales, `fix_add`, `d_gfnff_fix`); a non-finite or out-of-range value makes the step NaN; the 128 x 128 Coulomb tile walks in lock step. Bitwise run-to-run identical on gfx1200 (AI, machine-tested).
- `gpu_disp_pairs_on_device` (default true on CUDA and ROCm) builds the D4 list there in the host order; ROCm builds the gather CSR on the device too (`k_disp_rows_*`).
- GPU GFN-FF with net charge and several fragments (frag_charge_model ensemble) runs on the CPU method: `requestsCpuFallback()` in the wrapper, `EnergyCalculator::setMolecule()` swaps, counted GPU fallback; `-gfnff.frag_charge_model reference` keeps the GPU.
- GPU GFN-FF with an implicit solvent (`-gfnff.solvent`, ALPB) also runs on the CPU method (counted fallback, `cpuFallbackReason()`); before Oct 17, 2026 the GPU path returned the gas-phase energy without a message.
- ROCm (Oct 2026): the gradient-mode dispersion is the fused gather `k_dispersion_gather_uv` (C6 and dC6/dCN per edge from per-atom u/v vectors, no `k_dc6dcn_per_pair`); `CURCUMA_ROCM_DISP_FUSED=0` selects the unfused path.
- GPU `k_dispersion` reads no EEQ charge (`zetac6` holds the Phase-1 topology charges, C6/dC6dCN come from the CN); since K2 (Oct 2026) it runs right after `dc6dcn`, parallel to the EEQ solve.
- GFN-FF dispersion follows the reference: D4 C6 from Casimir-Polder integration with CN-only Gaussian weights and the modified BJ damping of `gfnff_gdisp0.f90`, not standard D3/D4 (`../dispersion/d4param_generator.cpp`). The three-body ATM term is not in the reference and is off by default (`dispersion_atm`).
- `uff-d3`: `ForceFieldGenerator::GenerateUFFD3Parameters()` merges UFF terms with the D3 pair list.
- A shared CN utility exists (`cn_calculator.h`); `D3ParameterGenerator` and `EEQSolver` still have their own CN routines.

## Open items

- Not implemented relative to the reference (notes of Apr 2026, not re-checked): hyperconjugation torsion modulation, metal-specific C6.
- EEQ on GPU: FMM matrix-vector product (O(N log N)) is open; ROCm runs the device projected PCG (Jacobi or `localinv`, built on the device), many-fragment systems use the exact CPU EEQ only when it stalls or `eeq_rocm_cpu_fragment_threshold` is given explicitly.
- ROCm GFN-FF correctness gaps and GFN-FF open items: root `CLAUDE.md` "Open Items", [docs/MULTI_GPU_GAPS.md](../../../../docs/MULTI_GPU_GAPS.md), [docs/REV_GFNFF_TODO.md](../../../../docs/REV_GFNFF_TODO.md).

## Reactive topology mode (`-gfnff.react`)

Detail: [docs/GFNFF_REACT_TOPOLOGY.md](../../../../docs/GFNFF_REACT_TOPOLOGY.md). Not re-checked after the merge of master (Oct 2026).
- NVT-only (dE_jump at rebuild events, logged); rebuild == fresh-init bit-identical (ctest `gfnff_react_fd_gradient`, `cli_simplemd_13/14/15`, `gfnff_react_filters`).
- Public topology/event API on `GFNFF`: `reactiveBonds()`, `reactiveBondOrders()`, `reactiveRebuildCount()`, `consumeReactEvents()`; reached through `ComputationalMethod::gfnffInstance()` (a vtable hook, because the GPU method lives in a dlopen plugin).
- The react set owns the bond list even when empty (no silent geometric re-detection); refractory blocks N scans; exchange resolution decrements sigma counts and skips same-scan formations.
- `rattle` is refused with react (constraints are frozen at initialisation).

## References

FFWorkspace implements the formulas of Fortran `gfnff_engrad.F90`; parameter generation follows Spicher and Grimme, Angew. Chem. Int. Ed. 59 (2020) 15665; D3: Grimme et al., J. Chem. Phys. 132, 154104 (2010).
