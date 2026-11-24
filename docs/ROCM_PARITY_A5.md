# ROCm parity A5: remaining GFN1/GFN2 gaps (Nov 24, 2026)

Item A5 of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md). AI-generated, machine-tested on one RX 9060 XT (gfx1200,
16 GB) only; human production testing pending. Base for every comparison: the worktree HEAD (`feature/rocm-parity` plus the
plan patch) built once before the changes (`a5_base`), same build flags.

## Result in one paragraph

(1) `-scf_precond frag` on the ROCm resident loop was already ported (WP5, Oct 24, 2026, `-xtb.rocm_resident_precond`, default
true); the plan entry "ignored with a warning" was out of date. It was re-measured on the current tree (with the WP10 fixed-point
gradient) and got a ctest (`cli_sqm_23`). (2) `-large_system_mode` has no device path on CUDA either: the drivers run on the CPU
on every backend. The CUDA screened (sparse) integral storage (`-gpu_sparse_integrals`) is not ported: it is not a mechanical
mirror (reason below). What changed instead: three silent CPU paths now say so, two of them as counted GPU fallbacks.

## Changes

| File | Function | Change |
|---|---|---|
| `src/core/energy_calculators/qm_methods/xtb_gpu_adapter.h` | `XtbGpuAdapter::setMolecule` | after the CPU pipeline: if a `-large_system_mode` driver owns the calculation and the device is live, one warning and `reportGpuFallback` (`-gpu_strict true` stops, exit 3). Shared by CUDA, ROCm, Vulkan |
| `src/core/energy_calculators/qm_methods/native_xtb_method.h` | `largeSystemDriverActive()` | new inline accessor (`m_c1_driver != nullptr`) |
| `src/core/energy_calculators/qm_methods/xtb_hip_method.cpp` | `HipScfBackend::beginBasis`, new `beginComputed`, `warnDeviceFallback` | a failed device basis setup or integral build is a warning plus a counted fallback, as `XtbGpuScfBackend` does on CUDA (`xtb_gpu_method.cpp`). Before: a verbosity-2 info line ("GPU resident GFN2 begin() failed; running CPU SCF") and `-gpu_strict true` did not stop |
| `src/core/energy_calculators/qm_methods/xtb_hip_method.cpp` | `XtbHipComputationalMethod` constructor | `-gpu_sparse_integrals on` prints one warning (screened storage is CUDA only); `auto` (default) and `off` stay quiet. Not a counted fallback: nothing moves to the CPU |
| `test_cases/cli/sqm/23_gfn2_scf_precond_rocm/`, `test_cases/cli/CMakeLists.txt` | new ctest `cli_sqm_23_gfn2_scf_precond_rocm` (`USE_ROCM`, label `rocm`, skip 77) | cli_sqm_13 with `-gpu rocm` plus a run-to-run md5 check of the frag gradient |
| `docs/SCF_PRECOND.md` | | stale "Not implemented: GFN1 on CUDA" removed (done Oct 25, 2026); ROCm re-check line |
| `docs/SQM_LARGE_SYSTEMS.md` | | paragraph "GPU": drivers are CPU only, counted fallback |

Plugin ABI change: no. `ComputationalMethod`, the plugin entry points and every class shared across the `dlopen` boundary are
unchanged; the new members (`m_large_system_reported` in `XtbGpuAdapter`, `m_fallback_warned` in `HipScfBackend`) live in classes
that are instantiated and used inside the plugin only, the accessor is inline and reads an existing member.

## Measurements

Conditions: `-sp -threads 4 -verbosity 2 -dump_gradient`, n = 1 per cell, GPU runs under `flock`. Peo201-chain segmented
(`-scf_precond_segment_atoms 3 -scf_precond_max_fragment_atoms 20`). dE and dG from the full-precision header and the Eh/A rows of
the gradient files.

SCF iterations none / frag, default threshold (1e-5):

| System (atoms) | GFN2 CPU | GFN2 ROCm | GFN1 CPU | GFN1 ROCm |
|---|---|---|---|---|
| urea3-water8 (48) | 10 / 6 | 10 / 6 | 9 / 6 | 9 / 6 |
| urea-water 197 unrelaxed | 23 / 33 (guard) | 23 / 33 (guard) | 13 / 8 | 13 / 8 |
| caffeine (24) | 12 / 4 | 12 / 4 | 10 / 4 | 10 / 4 |
| acetic acid dimer (16) | 12 / 6 | 12 / 6 | 9 / 5 | 9 / 5 |
| peo201-chain (1410) | 11 / 7 | 11 / 7 | 9 / 9 | 9 / 9 |

- ROCm = CPU iteration count in 20 of 20 cells; device M uploaded in 10 of 10 ROCm frag runs, no "ignored", no fallback.
  CUDA reference (SCF_PRECOND.md, other card and thread count): GFN1 urea3-water8 11 / 8, GFN2 CPU 11 / 6, unrelaxed 197 24 / 33.
- Energies: ROCm frag vs CPU frag max 3.4e-10 Eh (peo201-chain GFN2), ROCm frag vs ROCm none max 1.5e-9 Eh (same); every
  other pair <= 8.8e-11 Eh.
- Gradients at the default threshold: ROCm frag vs CPU frag up to 6.7e-6 Eh/A, but CPU none vs CPU frag also up to 4.7e-6 Eh/A:
  the SCF threshold, not the device. With `-scf_threshold 1e-10` (8 systems: the four small ones, GFN1 and GFN2): ROCm frag vs CPU
  frag max dE 5.1e-13 Eh, max dG 4.9e-11 Eh/A; CPU none vs CPU frag max dG 8.3e-11 Eh/A; ROCm none vs CPU none max dG 5.3e-11 Eh/A.
  Iterations at 1e-10 equal between CPU and ROCm in 16 of 16 cells.
- Second ROCm frag run: gradient file md5-identical in 10 of 10 (fixed point default).
- Without the new switches, final binary against the base: gradient file md5 identical in 24 of 24 runs (GPU none, GPU frag, CPU
  none; GFN1/GFN2 on urea3-water8, caffeine, acetic acid dimer, unrelaxed 197).
- Negative control: `-xtb.rocm_resident_precond false` gives the "ignored" warning, no device M, 10 iterations, so `cli_sqm_23`
  (b) would fail on that path.
- `-large_system_mode fragments -gpu rocm` (urea3-water8): energy identical to the CPU run (-82.84237071 Eh, gradient md5
  equal), one warning, one counted fallback; with `-gpu_strict true` exit 3. Before: same energy, no message.
- Dense device memory limit: urea400-water1000 (6200 atoms) fails the ROCm basis setup (GFN2 nao 14000, 1.57 GB per dense
  buffer; GFN1 nao 17600, 2.48 GB), now with a warning and exit 3 under `-gpu_strict true` (2 s); before, the run continued on the
  CPU after a verbosity-2 line. peo201-water500 (2910 atoms) runs on the device (WP numbers). The limit between 2910 and 6200 atoms
  was not bracketed.
- ctest (final binary, `-R "cli_sqm|large_system|rocm" -j1`): see the section Tests below.

## Not ported: screened (sparse) integral storage

CUDA (`cuda/xtb_gpu_context.cu`, Sep/Oct 2026): `buildScreenedPairs` (atom-pair cutoff from the basis, pattern row/col/colptr/perm,
auto when < 50 % of the AO pairs survive), `k_sp_gather`/`k_sp_scatter`, `k_multipole_ints_sp`, `k_build_fock_sp`, `k_pop_ao_sp(P)`,
`k_multipole_moments_sp(P)`, `k_density_sp_sym` (density on the pattern only), `k_grad_h0_pulay_sp(P)`, `k_grad_cn_onsite_spP`,
and about 220 references to the sparse state spread over the integral build, Fock build, density, moments, gradient, downloads,
the distributed eigensolver and the memory estimate. Reasons it was not mirrored here:

- The ROCm context (`rocm/xtb_hip_context.hip`) is a separate implementation with its own buffer set (`dS0`, `dH0`, `dP`, `dDpInt`,
  `dQpInt`, used in at least 30 places across the host-driven loop, the resident loop, pseudo-diagonalisation, the gradient and the
  lazy host copies of WP8b). The kernels can be copied, but without rewiring every consumer they would be dead code; the rewiring
  is a new design of the ROCm data flow, not a mirror.
- The CUDA default is `auto`: a port that follows CUDA switches large ROCm systems to the screened path by default, which changes
  the bits of every large-system run (the gate "bitwise equal without a new switch" then needs a ROCm default of `off`).
- What it would buy on ROCm: memory, not speed of the eigensolve. S stays dense (Cholesky), F/C stay dense; the 10 nao^2
  arrays H0, dp(3), qp(6) and P shrink to nnz. At 6200 atoms (GFN2, nao 14000) these ten arrays are 10 x 1.57 GB = 15.7 GB dense (computed, not measured). The
  pattern-only density saved 26 % of the CUDA SCF time on polymer_2x (CUDA comment), not measured on ROCm.
- Estimate: several hundred lines of HIP plus host changes, each consumer to be gated against the dense path bit for bit; a work
  package of its own (proposal below), not part of A5.

`-large_system_mode` (fragments, dc, sparse): CPU drivers on every backend (`xtb_gpu_method.cpp` comment "GPU +
large_system_mode is not yet wired and runs on the CPU fragment driver"). Nothing to port; now reported as a fallback.

## Not tested

- CUDA and Vulkan with the adapter change (`XtbGpuAdapter::setMolecule`): compiled only for ROCm here (no CUDA toolchain, Vulkan
  not built). On CUDA the same warning and counted fallback are expected; a CUDA run of `-large_system_mode fragments -gpu cuda
  -gpu_strict true` should exit 3.
- The skip branch of `cli_sqm_23` (no ROCm device) was not run.
- `-scf_precond frag` on ROCm: open shell and charged systems (WP5 measured them), Hessian and MD on the GPU, guard variants
  other than the default, systems above 1410 atoms, mixed precision off.
- Fault injection of the new ROCm counted fallback other than the real memory limit (no injection hook for `beginBasis`).

## Proposed central-file lines (for the coordinator)

- TODO.md: "ROCm: screened integral storage (`-gpu_sparse_integrals`, CUDA `buildScreenedPairs` and the `*_sp` kernels) not
  ported; ROCm GFN1/GFN2 stops at the dense device memory (6200 atoms fail on 16 GB); own work package, ROCm default off (docs/ROCM_PARITY_A5.md)."
- AIChangelog.md: "Nov 24, 2026: ROCm parity A5: -scf_precond frag on the ROCm resident loop re-measured (iterations = CPU 20/20,
  dG 4.9e-11 Eh/A at 1e-10, md5 stable), ctest cli_sqm_23; -large_system_mode with -gpu and a failed ROCm basis/integral build are
  counted GPU fallbacks; -gpu_sparse_integrals on warns on ROCm (docs/ROCM_PARITY_A5.md)."
- GPU_PARITY_CUDA_ROCM.md table A row "Multi-GPU eigensolver (cuSOLVERMg/Mp), sparse integrals": ROCm column "no; sparse integrals:
  `-gpu_sparse_integrals on` warns, dense storage; a failed dense allocation is a counted fallback since Nov 24, 2026 (A5)".
  Table A row `supportsResidentLoop`: unchanged (precond already listed).
- `scripts/gpu_parity_known.json`: none needed by this change (`beginComputed`/`beginBasis` are not `supports*` overrides).
- Plugin ABI: no change.

## Tests

`ctest -R "cli_sqm|large_system|rocm" -j1` with the final binary (only the `curcuma` target and the ROCm plugin built): 137 selected,
129 passed (74 `structinv_rocm_*`/`structcomp_rocm_*`, 18 `cli_sqm_*` including the new `cli_sqm_23`, the `large_system` test), 2
disabled (`cli_sqm_10`, `cli_sqm_11`), 5 not run (executables `gfnff_rocm_energy_only_geometry_*`, `hbond_energies_rocm` not built
in this target set), 1 failed: `cli_sqm_16_scf_criterion_xtb` ("gfn1: tight accuracy iterations 7, normal 7"), CPU only, fails
identically with the base binary, so not caused by this change.
