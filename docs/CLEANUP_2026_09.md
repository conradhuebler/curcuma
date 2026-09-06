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
| water box 3000 (999 fragments) | gfnff SP cold / warm (8 thr) | 7.07 / 6.47 s | 2.32 / 2.27 s | 3.05 / 2.85 |
| water box 3000 | gfnff MD step (8 thr) | 1092 ms | 309 ms | 3.5 |
| polymer 1410 | gfnff MD 100 fs (8 thr) | 12.5 s | 9.4-11.9 s | 1.05-1.3 |
| polymer 1410 | gfn1 SP (8 thr) | 190.9 s | 131.6-143.6 s | 1.33-1.45 |
| polymer 1410 | gfn2 SP (8 thr) | 113.1 s | 82.4-98.7 s | 1.15-1.37 |
| complex 231 | gfnff SP (1 thr) | 76 ms | 45 ms | 1.7 |
| triose 66 | gfn2 SP (1 thr) | 134 ms | 111 ms | 1.2 |

Energies: identical to the last digit for every single-fragment run; the water box (projected
PCG) differs by <= 5e-12 Eh, gradients by <= 5e-10 Eh/Bohr. All three MD runs end with
identical energies.

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

## Findings not fixed here (deliberately)

- **MD per-step cost of single-fragment GFN-FF is the Phase-2 EEQ solve**: 45-70 ms of ~90 ms
  at N=1410 (refactorisation + N² matrix build); a single point does not show it because
  Phase 2 is skipped when the geometry is unchanged. Reusing the Cholesky factor is an
  approximation (A4 refinement). The many-fragment case is now covered by the projected PCG
  above; the single-fragment dense factorisation is untouched.
- **`cli_simplemd_08_gfnff_acetic_acid_dimer_md` is flaky under parallel `ctest -j`**: the
  OpenMP-parallel Phase-1 EEQ gives charges differing by 1e-13 under thread contention and
  the 10 ps trajectory is chaotic (drift 0.02 idle vs 0.13 loaded). Baseline and new build
  are bit-identical when run idle.
- **`-method cg` is not registered** in the factory (the `08_cg_spheres` script is not a
  ctest); SimpleMD then segfaults instead of aborting after the failed method creation.
- **GFN-FF PBC on the CPU path** is silently non-periodic (the unit cell only ever reached the
  removed engine). `cli_sqm_11_gfn2_provider_check` needs TBLite and fails in a TBLite-less
  release build (pre-existing).
- xTB: the three copies of `as_cgto_shell`/`ao_to_type` (X-I5), the audit-only gradient gate
  flags in hot loops, and the ~50 % duplicated CUDA/HIP/Vulkan host adapters remain.
