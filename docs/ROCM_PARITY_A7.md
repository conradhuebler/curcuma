# ROCm parity A7: screened (sparse) integral storage (Nov 24, 2026)

Item A7 of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md), the work package proposed in
[ROCM_PARITY_A5.md](ROCM_PARITY_A5.md) ("Not ported: screened (sparse) integral storage"). 🤖 AI-generated, ⚙️ machine-tested
on one RX 9060 XT (gfx1200, 16 GB) only; human production testing pending. Base for every bitwise comparison: the worktree
HEAD (`feature/rocm-parity`, `ebd1b357`) built once before the changes, same flags (`-DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON
-DUSE_AVX2=ON -DCMAKE_BUILD_TYPE=Release`); the `curcuma` executable of base and final build is md5-identical, only the plugin
`libcurcuma_rocm.so` differs.

## Result in one paragraph

`-gpu_sparse_integrals off|auto|on` now works on ROCm with the CUDA values and thresholds (integrals below 1e-20 dropped,
element-pair cutoff from the basis, `auto` = screened when fewer than 50 % of the AO pairs survive). Only the default differs:
without the key ROCm stores dense (CUDA default `auto`), so no existing ROCm run changes. A relaxed 4848-atom system (GFN2 nao 10500) that
fails dense on the 16 GB card runs screened on the device (10.5 % of the AO pairs stored) and agrees with the CPU to 4.1e-11 Eh
and 1.6e-11 Eh/A; the full urea400-water1000 (6200 atoms, nao 14000) fits screened (peak 14.5 GB device memory) but its
unrelaxed geometry does not converge, on the CPU either (sphere cut, see below). Screening error against the dense device
run at the CUDA threshold: at most 4.1e-12 Eh and 1.1e-12 Eh/A (10 cells). Without the switch the gradient files are
md5-identical to the base build in 16 of 16 runs.

## What was ported (files, functions)

| File | What |
|---|---|
| `src/core/energy_calculators/qm_methods/rocm/xtb_hip_sparse.hiph` (new) | kernels `k_sp_gather`, `k_multipole_ints_sp`, `k_build_fock_sp`, `k_density_sp_sym`, `k_popband_spP`, `k_multipole_moments_spP`, `k_res_mom_part_spP`, `k_grad_cn_onsite_spP`, `k_grad_h0_pulay_spP`; host `hipScreenRadius` (verbatim CUDA `screenRadius`). Twins of the dense ROCm kernels, same expressions; atomics through `xacc_add` (fixed point with `-xtb.rocm_fixed_point`) |
| `rocm/xtb_hip_sparse_impl.hiph` (new, included inside `XtbHipContext::Impl`) | state (`sp_mode`, `sp_on`, pattern and value buffers), `spCaptureBasis`, `spBuildPattern` (port of CUDA `buildScreenedPairs` without the dense-build screen), `spPrepare` (storage decision + allocation per geometry), `spEnsureDense`, `spFreeDense`, `spFreeArrays`, consumers `spFock`, `spDensity`, `spPopBand`, `spScatterHost`, `spDownloadCholesky` |
| `rocm/xtb_hip_context.hip` | one `if (sp_on)` branch per consumer, the dense statement unchanged in the `else`: `ensureResident`/`begBasis` (dense nao^2 buffers deferred only when `sp_mode != 0`), `resBegin` (host-uploaded integrals: dense, as CUDA), `computeIntegrals` (`spPrepare`, H0 scratch in `dCw`, gather, Cholesky probe in `dS`), `computeMultipole`, `downloadMP`, `resSolve`, `resSolveMP`, `resSolveMPResident`, `fockAndSolve`, `resDensity`, `resFinalize`, `multipoleMomentsImpl`, `resStep` (density, populations, moments), `gradient` (W on the pattern, pseudo-step W from G_oo, screened H0/Pulay and CN kernels), spin functions (refuse), `downloadH0`, `downloadCholesky`; public `setSparseIntegrals`, `sparseIntegrals`, `sparseFraction`, `sparseCutoffBohr`, `sparseError` |
| `rocm/xtb_hip_context.h` | declarations of the five public functions |
| `xtb_hip_method.cpp` | the A5 warning "ignored on ROCm" replaced by `setSparseIntegrals` (CUDA mapping); `HipScfBackend::beginComputed` prints the CUDA storage line (verbosity 2) and passes `sparseError()` into the fallback warning |
| `CMakeLists.txt` | the two new `.hiph` in the `DEPENDS` list of the hipcc command |

Storage with screened storage active (GFN2): dense `dS0` (S, read by every eigensolver path, the factor probe and the
pseudo-diagonalisation), `dF`, `dS`, `dCw` (and the FP32 copies of the mixed-precision solve); on the pattern: S, H0, P, W and
the nine multipole components. Not allocated: dense H0, P, W, L, dp_int, qp_int (13 nao^2 doubles less than the dense GFN2 set).

## Measurements

Conditions: `-sp -threads 4 -verbosity 2 -dump_gradient`, n = 1 per cell, GPU runs under `flock`; dE and max |dG| from the
full-precision header and the Eh/A rows of the gradient files (script: pure Python). "pattern" = stored fraction of the AO
pairs reported by the run (GFN2 / GFN1; largest element-pair cutoff 29-31 Bohr GFN2, 34-36 Bohr GFN1).

Gate table, `-scf_threshold 1e-10`, screened = `-gpu_sparse_integrals on`, dense = same final binary without the switch:

| System (atoms) | pattern | screened vs dense, GFN2 dE / dG | GFN1 dE / dG | screened vs CPU, GFN2 dE / dG | GFN1 dE / dG |
|---|---|---|---|---|---|
| hydrogen-chloride.d-shell-validation (2) | 100 / 100 % | 0 / 0 | 0 / 0 | 0 / 8.8e-12 | 8.9e-16 / 2.9e-11 |
| caffeine (24) | 100 / 100 % | 0 / 0 | 0 / 0 | 0 / 4.2e-12 | 2.1e-14 / 5.3e-11 |
| pep20-peptide.gfnff-min (309) | 73.4 / 84.9 % | 5.7e-14 / 6.4e-15 | 1.1e-13 / 2.1e-14 | 1.6e-12 / 5.5e-13 | 1.3e-12 / 1.1e-11 |
| water216-box.sc-s1 (648) | 59.0 / 78.0 % | 0 / 1.6e-14 | 4.6e-13 / 1.1e-14 | 4.8e-12 / 4.0e-12 | 9.1e-12 / 1.6e-11 |
| peo201-chain (1410) | 41.4 / 54.5 % | 1.4e-12 / 3.0e-14 | 4.1e-12 / 1.1e-12 | 1.8e-12 / 1.3e-12 | 1.2e-10 / 1.4e-11 |

- Screening error (screened vs dense, the CUDA rule of documenting it): max dE 4.1e-12 Eh, max dG 1.1e-12 Eh/A over the 10
  cells; it contains the dropped integrals (< 1e-20 each) and the different summation order of the pattern density
  (`k_density_sp_sym` against the rocBLAS GEMM). With nothing dropped (100 %) the gradient file is md5-identical to dense.
- Against the CPU every cell is within 1.2e-10 Eh / 5.3e-11 Eh/A; the 1.2e-10 Eh of peo201-chain GFN1 is the dense device
  path as well (dense vs CPU 1.1e-10 Eh), not the screen. Gate 1e-8 met in 10 of 10 cells.
- Second screened run: gradient file md5-identical in 10 of 10 (fixed-point accumulation, default) plus the large run below.
- Without the switch against the base build: md5-identical in 10 of 10 at `-scf_threshold 1e-10` and 4 of 4 at the default
  threshold with the final plugin (caffeine, peo201-chain, GFN1/GFN2); `-gpu_sparse_integrals off` (peo201-chain GFN2) and
  `auto` on water216 GFN2 (59 % -> dense through the deferred allocation `spEnsureDense`) are md5-identical to dense; `auto` on
  peo201-chain GFN2 (41 %) is md5-identical to `on`.
- Other paths, peo201-chain, screened vs dense / vs CPU (dE Eh, dG Eh/A): host-driven loop (`-xtb.rocm_resident_loop false`)
  GFN2 9.1e-13, 1.6e-14 / 0, 1.3e-12 and GFN1 1.4e-12, 1.1e-12 / 0, 1.4e-11; pseudo-diagonalisation with FP64 steps (20 steps
  taken, gradient W from G_oo) 0, 7.0e-15 / 4.1e-12, 5.3e-12; open shell (charge 1, doublet, resident loop) 1.8e-12, 6.3e-13 /
  1.5e-11, 3.7e-11. Spin polarisation with `on`: one warning, dense storage, gradient md5-identical to the dense spin run.
- Finite differences of the assembled surface (pep20, central, h = 1e-4 A, `-scf_threshold 1e-11`, atoms 0/150/300, 9
  components): screened GFN2 max |FD - analytic| 1.68e-8 Eh/A, dense GFN2 at the same coordinates 1.79e-8, screened GFN1
  2.25e-8 (dense GFN1 not run). The FD error of the method itself is of this size (energy noise / 2h); screened and dense agree.

Large systems (GFN2, `-gpu_strict true`):

| System | atoms, nao | dense (base and final) | screened |
|---|---|---|---|
| peo201x2-water1500.gfn2-relaxed, both chains + 676 nearest waters (whole molecules, cut by connectivity) | 4848, 10500 | exit 1: multipole matrices and the FP32 eigensolve do not fit, the device loop fails at iteration 0 | 10.5 % of the pairs; device run, 23 SCF iterations (CPU 23), SCF 654 s (CPU 4 threads 2794 s, shared machine); vs CPU dE 4.1e-11 Eh, dG 1.6e-11 Eh/A (`-scf_threshold 1e-10`); second run md5-identical |
| urea400-water1000 | 6200, 14000 | basis setup fails (A5) | 11.5 % of the pairs, cutoff 30.3 Bohr; 5 device iterations of about 30 s with `-scf_max_iterations 4`, peak device memory 14.49e9 B (sampled every 2 s while this run held the GPU); the SCF diverges (max dq 4.8) |

- The urea400-water1000 geometry is unrelaxed: a 4626-atom sphere cut of it diverges in the same way on the CPU (iteration 3:
  +4.5e4 Eh) and on the device, so no converged energy of this system exists to compare; the relaxed PEO system above is
  the accuracy test of the large-system path.
- Timing of the small systems was not compared: the GPU lock was shared with other agents.

## Not tested

- Any GPU other than gfx1200; peak memory of the 4848-atom run (only the 6200-atom run was sampled); whether 6200 atoms fit
  with `-scf_mixed_precision false` or with `-scf_pseudo_diag` (extra factor buffers).
- Geometry optimisation and MD with screened storage (a new pattern per geometry; only single points were run), Hessian.
- Solvation (ALPB/GBSA/ddCOSMO) with screened storage, `-scf_precond frag` with screened storage, `-xtb.rocm_cached_cholesky`
  / `-xtb.rocm_explicit_linv` / `-xtb.rocm_density_syrk` with screened storage (the syrk switch is ignored while screened).
- Elements beyond H, C, N, O, Cl; the host copies (`downloadH0`, `downloadMultipoleInts`, `downloadCholesky` scatter
  paths) are compiled but no consumer asked for them in the runs above (lazy host copies).
- A converged SCF of urea400-water1000 itself.

## Not implemented relative to CUDA

- Dense S stays on the device (CUDA keeps only the factor L dense); the atom-pair screen of the dense S/H0 build
  (CUDA `dAtKind`/`dScreenCut2`, `CURCUMA_GPU_FULL_DENSE_S`) is not ported: the dense S/H0 build computes every pair.
- No memory estimate before allocation (CUDA `estimateStorageBytes`/`gpu_memory_check`): a failed `hipMalloc` fails the integral
  build, which is the counted fallback of A5 (`-gpu_strict true` exits 3).
- Spin polarisation (`-xtb.spin_polarization`) with screened storage is refused on the device (CUDA supports it): the spin
  functions return false and the caller's fallback runs.
- Distributed (multi-GPU) density on the pattern: not applicable (one GPU, no ROCm multi-GPU density).
- Periodic cells: the host-uploaded integrals are always dense (same as CUDA).

## Proposed central-file lines (for the coordinator)

- TODO.md: replace the A5 line "ROCm: screened integral storage ... not ported" by "ROCm screened storage (A7): spin
  polarisation refused (dense), no memory estimate before allocation, dense S kept; opt/MD/solvation with it untested
  (docs/ROCM_PARITY_A7.md)."
- AIChangelog.md: "Nov 24, 2026: ROCm parity A7: -gpu_sparse_integrals on|auto on ROCm (CUDA thresholds, ROCm default dense):
  screened vs dense <= 4.1e-12 Eh / 1.1e-12 Eh/A, vs CPU <= 1.2e-10 Eh; 4848-atom GFN2 system (dense does not fit 16 GB) on the
  device, 4.1e-11 Eh vs CPU; default md5 = base 16/16; ctest cli_sqm_24 (docs/ROCM_PARITY_A7.md)."
- GPU_PARITY_CUDA_ROCM.md table A row "sparse integrals": ROCm column "yes (A7, Nov 24, 2026), default off (CUDA auto); spin
  polarisation with it refused (dense); no memory estimate".
- `scripts/gpu_parity_known.json`: no entry needed. `scripts/gpu_parity_check.py` on the final tree: 0 unexplained differences;
  the CUDA-only dense-P twins (`k_grad_h0_pulay_sp`, `k_multipole_moments_sp`, `k_pop_ao_sp`, `k_pop_ao_spP`, `k_sp_scatter`)
  and the new ROCm names (`k_popband_spP`, `k_res_mom_part_spP`) are listed as info only.
- Plugin ABI: no change. New public functions of `XtbHipContext` and new members of `HipScfBackend`/`Impl` live inside
  `libcurcuma_rocm.so`; `ComputationalMethod`, the C entry points and the shared classes are unchanged (the `curcuma` executable
  of base and final build is md5-identical).
- docs/SQM_LARGE_SYSTEMS.md / GPU_TUNING.md (not central, left unchanged here): the `-gpu_sparse_integrals` row should say "ROCm:
  default off".

## Tests

New ctest `cli_sqm_24_gfn2_sparse_integrals_rocm` (`test_cases/cli/sqm/24_gfn2_sparse_integrals_rocm/`, `USE_ROCM`, labels
`sqm;gpu;rocm;cli`, skip 77): pep20 GFN2 and GFN1 screened against dense (1e-8 Eh, 1e-8 Eh/A), storage line, `-gpu_strict true`
exit 0, md5 of a second run; passed (12.9 s). The skip branch was not run.

`flock ... ctest -R "cli_sqm|large_system|rocm" -j1` with the final binary (targets `curcuma` and `curcuma_rocm` only): 150
selected, 133 passed (incl. `cli_sqm_24`, `cli_sqm_23_*`, the `structinv_rocm_*`/`structcomp_rocm_*` set, `large_system`),
2 disabled (`cli_sqm_10`, `cli_sqm_11`), 13 not run (test executables `gfnff_rocm_pbc_*`, `gfnff_rocm_energy_only_geometry_*`,
`pbc_xtb_gpu_*_rocm`, `hbond_energies_rocm` not built in this target set), 2 failed and unrelated to this change:
`gfnff_rocm_charge_groups` (`ModuleNotFoundError: numpy` in its Python checker) and `cli_sqm_16_scf_criterion_xtb` (CPU only,
fails identically with the base binary, ROCM_PARITY_A5.md).
