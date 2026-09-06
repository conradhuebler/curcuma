# GFN-FF / GFN1 / GFN2 cleanup and speedup (September 2026)

> 🤖 AI-generated, machine-tested. Every change below was verified to leave the energies of
> 86 single points (gfnff / gfn1 / gfn2, 3 to 3000 atoms, incl. charged, transition-metal
> and multi-fragment systems) identical to the last printed digit (12 decimals) and the
> GFN-FF gradients bit-identical; a 100-step polymer MD ends with identical energies. Human
> production testing is pending.

## Goal

Review the native GFN-FF, GFN1 and GFN2 implementations, remove dead and duplicated code,
make them faster on CPU (and keep the GPU paths working) without changing a single number,
and leave a structure that makes adding further QM or force-field methods cheap.

## What was removed (about 20 000 lines)

| Area | Removed | Why it was safe |
|---|---|---|
| qm_methods | `gfn1_params.hpp`, `gfn2_xtb_params.hpp`, `gfn2-xtb_param.hpp` (duplicate tables), `integrals.h`, `test_parser.cpp`, `vulkan/prototype/` | zero includers; the live tables are `parameters/gfn{1,2}_params.hpp` |
| qm_methods | `am1/pm3/pm6/mndo.{h,cpp}` + their `*_method` wrappers | not in CMake since the unified `NDDOMethod` (commit 85570d24); `am1`/`pm3`/`mndo`/`pm6` energies unchanged |
| qm_methods | `XTB::MakeOverlap/MakeH` stubs, `QMDriver` base of the native xTB | only existed to satisfy `QMDriver`; drops `STOIntegrals`/`LofthusOverlap` from the xTB include graph |
| ff_methods | the legacy `ForceField`/`ForceFieldThread` GFN-FF engine (`forcefieldthread.{h,cpp}` 4670 lines, GFN-FF/JSON parts of `forcefield.cpp` 3185 -> 887), `gfnff_advanced.*`, `forcefieldderivaties.h` | `GFNFF` had evaluated only through `FFWorkspace` since March 2026 but still fed the dead engine every step; `-method cg`/`d3` never reached it (not registered in the factory) |
| ff_methods/gfnff_method.cpp | legacy sync, JSON parameter path (`generateGFNFFParameters` + per-term JSON generators), diagnostics that only ran on the dead engine, 8 uncalled helpers (13350 -> 10690 lines) | dead by construction; validation-test JSON exports now read the native parameter set |
| ff_methods/eeq_solver | `EEQSolverCache` + `buildSmartEEQMatrix`, dense Floyd-Warshall, `buildNeighborLists`, `calculateEEQEnergy`, emoji self-test, write-only N x N distance cache | zero callers |
| dispersion | `d4_reference_data.cpp` (unreferenced twin of `_fixed`), `qmdff_terms.h` | zero includers |

Shared term structs now live in `ff_methods/ff_terms.h`.

## Exact (bit-identical) speedups

| Change | File | Effect |
|---|---|---|
| Persistent LAPACK `dsyevd`/`ssyevd` scratch instead of a 2·nao² allocation + memset per SCF iteration | `xtb_scf.cpp` | GFN1 polymer/1410 191 -> 132 s, GFN2 113 -> 99 s |
| `as_cgto_shell()` hoisted out of the gradient / CPSCF shell-pair loops (2 heap allocations per shell pair) | `xtb_gradient.cpp`, `xtb_response.cpp` | part of the gfn2 gradient gain |
| No per-step feed of the dead engine (parameter sync incl. a full CN-derivative rebuild, CN/charge distribution, shared-distance pointers, H-/X-bond JSON) | `gfnff_method.cpp` | GFN-FF SP 1.35-1.9x |
| H-bond bond-gradient: CSR index of `m_hb_grad_entries` by H atom (was a full linear scan per H-bonded bond) | `ff_workspace_gfnff.cpp` | water box bond term 141 ms -> few ms per step |
| `getenv` hoisted out of the H-bond triple loop (1.4 M calls per step on the water box) | `ff_workspace_gfnff.cpp` | |
| EEQ Phase 2: zero only the constraint rows/cols, `Eigen::Ref` views instead of N x N copies of `A_nn` and the fragment matrix | `eeq_solver.cpp` | |
| EEQ `TopologyInput` cached per topology version (was a deep copy of all neighbour lists per MD step) | `gfnff_method.cpp` | |
| Parameter set moved into the workspace; external copy holds bonded terms only unless a GPU wrapper calls `setKeepFullParameterSet(true)` (two 600 MB copies on the water box before) | `gfnff_method.cpp` | |

Final measurement (release/, -O3, AVX2, OpenBLAS, `OMP_NUM_THREADS=1`; baseline = db0c760f,
final = 215b3af1 incl. the projected-PCG EEQ below; run-to-run variance on this box ~10 %):

| System | Run | before | after | factor |
|---|---|---|---|---|
| polymer 1410 | gfnff SP cold / warm (8 thr) | 774 / 560 ms | 430 / 429 ms | 1.8 / 1.3 |
| polymer 1410 | gfnff SP cold / warm (1 thr) | 765 / 721 ms | 495 / 439 ms | 1.55 / 1.6 |
| water box 3000 (999 fragments) | gfnff SP cold / warm (8 thr) | 7.07 / 6.47 s | 2.09 / 1.99 s | 3.4 / 3.3 |
| water box 3000 | gfnff MD step (8 thr) | 1092 ms | 309 ms | 3.5 |
| polymer 1410 | gfnff MD 100 fs (8 thr) | 12.5 s | 7.5 s (round 2, PPCG default) | 1.65 |
| polymer 1410 | gfn1 SP (8 thr) | 190.9 s | 131.6-143.6 s | 1.33-1.45 |
| polymer 1410 | gfn2 SP (8 thr) | 113.1 s | 82.4-98.7 s | 1.15-1.37 |
| complex 231 | gfnff SP (1 thr) | 76 ms | 45 ms | 1.7 |
| triose 66 | gfn2 SP (1 thr) | 134 ms | 111 ms | 1.2 |

Energies: identical to the last digit for every run below 500 atoms; systems on the
projected-PCG default (polymer, water box) differ by <= 2e-12 Eh and <= 1e-12 Eh/Bohr from the
baseline. All MD runs end with identical energies. All timings were taken while four unrelated
`-confsearch` jobs of the operator (~24 cores) were running on the same box, so absolute times
are inflated and ratios carry ~10 % noise; both binaries were measured under the same load.

Small molecules (< 100 atoms) are dominated by process start-up (20-50 ms) and did not change.

## Many-fragment EEQ: projected PCG (approximate by tolerance, operator-approved)

The Schur-Cholesky solve needs `nfrag + 1` triangular solves (`A_nn^{-1} C^T`), i.e. O(N² nfrag):
on a 3000-atom water box (999 fragments) that was 887 ms per solve, 927 ms of a 1092 ms MD step.
`EEQSolver::solveWithProjectedPCG` runs ONE conjugate-gradient iteration sequence on the
constraint tangent space (projection = subtract the per-fragment mean, O(N); Jacobi
preconditioner projected the same way), warm-started from the previous step's charges. It is
selected automatically for `nfrag >= 32` and `N >= 500` (`eeq_ppcg_min_nfrag`, `eeq_ppcg_min_atoms`;
0 disables, `solve_method ppcg` forces) and converges to `eeq_ppcg_tol` (1e-10 relative).

| water box 3000 / 999 fragments, 8 threads | exact Schur-Cholesky | projected PCG |
|---|---|---|
| Phase-2 solve | 887 ms | 59 ms (14 iterations) |
| EEQ per MD step | 927 ms | 143 ms |
| MD step total | 1092 ms | 309 ms |
| energy | -328.184291307772 | -328.184291307772 |
| max gradient difference | — | 3e-10 Eh/Bohr |

Single-fragment and small systems are untouched (polymer/1410 identical). The constrained
problem is very well conditioned: the per-fragment constraint removes the long-range
charge-transfer modes, so 2-14 iterations reach 1e-13..1e-9 residuals.

## Structure for future methods

- `MethodFactory::methodTable()`: one `MethodDescriptor` row per method family (names, family,
  description, availability probe, providers, creator). `create()`, the `-methods` listing,
  `getMethodInfo()` and the availability checks read the table; unknown names get suggestions,
  unavailable providers a build hint. Adding a method = one row + the `ComputationalMethod`
  subclass with its `PARAM` block.
- `MethodFactory::methodParameterScopes()`: the JSON sub-scope names (`gfnff`, `xtb`, ...)
  that EnergyCalculator, opt/sp, SimpleMD, ConfSearch and the CLI router forward. Previously
  five slightly different private copies (the source of several "flag never arrived" bugs).
- `FFWorkspace` is the only force-field engine; the checklist for a new GFN-FF term is in
  `src/core/energy_calculators/ff_methods/CLAUDE.md`.
- Native xTB implements `QMInterface` directly; `QMDriver` is only for the STO/GTO EHT/NDDO path.
- `test_portable_math` pins the vendored fdlibm `erf/acos/exp/log` (the `CURCUMA_PORTABLE_MATH`
  determinism layer) to the host libm within 2 ulp plus special values.
- `scripts/regression_bench.sh` compares two binaries over a molecule ladder (energy to 12
  digits, max |dgrad|, wall time) — the tool used for every claim in this document.

## Verified GPU

`release_cuda` (RTX 5080) rebuilt after the refactor: `-gpu cuda` GFN-FF and GFN2 single
points on the 231-atom complex agree with the CPU (GFN-FF 1e-7 Eh as documented before, GFN2
to 8 decimals), GFN-FF GPU MD runs. ROCm/HIP could not be compiled here (no SDK); the HIP
wrapper received the same two-line change as the CUDA wrapper (`setKeepFullParameterSet`).

## Round 2 (same day): the open items

- **Load-dependent MD tests fixed at the root.** Phase-1 EEQ (`calculateTopologyChargesMultiRHS`,
  `solveEEQ`) ran its LAPACK solve without the `ScopedBlasThreads` guard, so OpenBLAS used all
  cores and its threaded kernels gave topology charges that differed by ~1e-13 with machine
  load; the chaotic 10 ps trajectories turned that into pass/fail coin flips under `ctest -j`.
  With the guard both drift tests pass 3/3 rounds under `-j8`. The MD test scripts also remove
  stale `input.snapshots/` and `input.topo.json` before running.
- **SimpleMD fails loud** when the energy method cannot be created (was a segfault in
  `FastEnergy()`): "MD setup failed: <factory reason>", and `prepareRun()` says why nothing runs.
- **`-method cg` restored** on the workspace engine (`ff_workspace_cg.cpp`, `-load_ff_json FILE`
  with `cg_default` / `cg_per_atom` / `pair_interactions` / `bonds`): same pair energy as the
  removed engine, now with an analytic gradient for spheres (finite differences for ellipsoids);
  `Elements::String2Element` accepts numeric atomic numbers ("226" = CG bead). Two-sphere
  check reproduces the closed form (E = 0.40328774 Eh, |g| = 1.058157e-1);
  `cli_simplemd_08_cg_spheres` is a registered ctest now.
- **Projected PCG is the default Phase-2 EEQ solve above 500 atoms for any fragment count**
  (`eeq_ppcg_tol` 1e-12): polymer/1410 solve 44 -> 16 ms (62 iterations cold, fewer warm),
  energy/gradient deviation 1e-12, MD 100 fs 10.4 -> 7.3 s with an identical trajectory.
  `eeq_ppcg_min_nfrag 0` restores the exact dense path; systems below 500 atoms stay exact.
- **xTB:** the 17 declared-but-never-defined `XTBMethod` members (X-M1) are gone;
  `as_cgto_shell`/`ao_to_type` were already shared in `xtb_ao_utils.hpp` (X-I5 was stale).
- **GPU wrappers:** the duplicated citation block and the stale "CPU residual" header
  description were removed in both the CUDA and the HIP wrapper.

## Round 3 (same day): GPU architecture — every backend is a plugin

Operator-approved plan (three tiers; tier 3 "unify kernel bodies" deliberately NOT done blind):

- **Plugin symmetry.** ROCm and Vulkan were compiled *into* `curcuma_core` while CUDA was a
  dlopen plugin, so `method_factory.cpp` / `energycalculator.cpp` carried a
  `#if USE_CUDA … #elif USE_ROCM …` ladder and every backend choice was a different core
  binary. Now `libcurcuma_{cuda,rocm,vulkan}.so` follow one recipe (entry files
  `qm_methods/{cuda,rocm,vulkan}/gpu_plugin_entry_*.cpp`, one `add_library(... SHARED)`
  block each, compile definitions mirrored from the core). The core has **zero** backend
  `#ifdef`s: `resolveGpuMode()` probes `gpu_plugin::available(backend)` at runtime, `-gpu auto`
  takes the first plugin found (cuda, rocm, vulkan), a missing plugin warns *"plugin
  libcurcuma_<b>.so is not present … cmake -DUSE_<B>=ON"* and falls back to CPU, `-methods`
  lists the plugins next to the binary. Details: [GPU_PLUGIN_STARTUP.md](GPU_PLUGIN_STARTUP.md).
- **Verified** (RTX 5080): CPU build — `ldd` shows no CUDA/Vulkan/HIP library, energies
  unchanged (caffeine gfn2 -42.14723025, complex gfnff -37.24064873 Eh), fallback warning
  shown. Vulkan plugin — caffeine/complex × gfn1/gfn2 identical to the pre-plugin binary to
  the printed digit (-42.14723025 / -44.50985543 / -329.52714784 / -343.17980354 Eh),
  35-step caffeine `-opt` on the device, `-gpu vulkan` for GFN-FF reports "shaders not
  ported" and uses CPU. CUDA plugin — same four xTB numbers plus complex gfnff -37.24064863
  and water box -328.18429136 Eh, `test_gfnff_gpu` 4/4. ROCm — CMake block mirrors CUDA
  line by line, **unverified** (no SDK here).
- **Pre-existing Vulkan test failures, unchanged by this round** (fail identically with the
  pre-plugin binary of the same commit): `cli_gpu_gradient_01/02_vulkan_*` expect a
  "Gradient norm" line that `-sp -verbosity 2` no longer prints (the CPU reference run fails
  the same way; the baseline db0c760f binary does not print it either), and 4 of the 20
  opt-in `sqm_val_vkdc_*` divide-and-conquer eigensolve tests (C6H6/acetic-acid-dimer/
  caffeine gfn1, caffeine gfn2) give NaN on this NVIDIA card with `CURCUMA_VK_TRIDIAG_SOLVE=dc`
  (the D&C was validated on an AMD RADV device, see SQM_VULKAN.md). The default Vulkan
  eigensolve path is fine.
- `cli_simplemd_08_cg_spheres` failed only because a stale 0-byte `input.snapshots/input.trj.xyz`
  from a killed run was picked up by `find_output_file`; the script now starts clean.

## Still open (deliberately)

- **GFN-FF PBC on the CPU path** (unit cell ignored by the workspace) — postponed by the operator.
- **ROCm plugin build** needs a machine with hipcc/rocSOLVER to confirm the converted CMake block.
- `cli_simplemd_13_rmsd_mtd_legacy_ab` fails only when run in parallel with its sibling
  (shared output names), `cli_sqm_11` needs TBLite, `cli_curcumaopt_07` has a pre-existing
  golden-value drift, `test_bmt_utils` / `test_orca_interface` / `xtb_cpscf` / `confscan_dtemplate`
  fail identically on the baseline.
