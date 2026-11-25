# ROCm parity SPE: spin polarisation with screened storage; A1r, G1c, E1 on the host-driven loop (Nov 25, 2026)

Two items of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md). 🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT
(gfx1200, 16 GB) only; human production testing pending. Base for every bitwise comparison: the worktree HEAD (`dd4812c4`) built
once before the changes, same flags (`-DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON -DUSE_AVX2=ON -DCMAKE_BUILD_TYPE=Release`). The
`curcuma` executable of base and final build is md5-identical (`c2098167...`); only `libcurcuma_rocm.so` differs. Kernel ISA of the
device object (`a2/isa.py`): 145 of 145 base kernels identical, two new kernels (`k_spin_w_combine_sp`, `k_twice_lower`).

## Result

1. Spin polarisation with `-gpu_sparse_integrals on|auto` runs on the device. Before, ROCm printed one warning and kept dense
   storage. Unlike CUDA (dense P, Pm, W next to screened S/H0), P = Pa + Pb, Pm = Pa - Pb and W' live on the screened pattern,
   like the closed-shell A7 path, so the later consumers are the existing screened kernels. Screened vs dense at most
   1.4e-12 Eh and 2.0e-13 Eh/A; vs CPU at most 3.2e-12 Eh and 1.2e-10 Eh/A.
2. A1r: not ported (measured: slower without FP64 GEMM emulation, which rocBLAS lacks). E1 ported, default on as on CUDA
   (SCF 1.7 % / 6 % faster on 216 / 432 H2O). G1c ported; it acts only with `-scf_pseudo_diag true`, which stays opt-in on ROCm
   (WP7), and then makes the SCF 37 % / 11 % faster than E1 alone.

## Item 1: what was ported

| File | What | CUDA reference |
|---|---|---|
| `rocm/xtb_hip_sparse_impl.hiph` | `spSpinEnsureSp`, `spSpinFreeSp`, `resSpinKeepAlphaSp`, `spSpinChannelSp` (k_scale_cols + `k_density_sp_sym` per spin), `resDensitySpinSp` (`k_spin_combine` element-wise on the nnz arrays, `k_popband_spP` for P and Pm), `resFinalizeSpinSp` (host P, Pm zero outside the pattern, as `resFinalize` with screened storage), `spinWSp`; state `spPmSp`, `spWbSp`, `pmsp_nnz` | `XtbGpuContext::residentSpinKeepAlpha/residentDensitySpin/residentFinalizeSpin/spinWeightedDensity` |
| `rocm/xtb_hip_sparse.hiph` | kernel `k_spin_w_combine_sp`: W' = (Wa + Wb) - 1/2 Pm (u_row + u_col) on the pattern | `k_spin_w_shift` |
| `rocm/xtb_hip_context.hip` | one `if (sp_on)` line in `resSpinKeepAlpha`, `resDensitySpin`, `resFinalizeSpin`, `spinW`; `freeSpin` calls `spSpinFreeSp`; the A7 refusal in `gradient` removed (with `w_spin` the screened kernels read dSpW) | - |
| `xtb_hip_method.cpp` | the refusal warning (spin + screened -> dense) removed | - |
| `test_cases/cli/sqm/24_gfn2_sparse_integrals_rocm/run_test.sh` | case (e): pep20 triplet GFN2 screened against dense, no refusal warning | - |

Memory with screened storage and spin: dense spCa (alpha vectors) and dF (beta), pattern spPmSp and spWbSp (2 nnz doubles); the
dense spPm (n^2) is not allocated.

## Item 1: measurements

`-sp -gradient -dump_gradient -threads 4 -verbosity 2 -scf_threshold 1e-10`, `-spin 2 -xtb.spin_polarization true`, GPU runs with
`-gpu_strict true` (no fallback, exit 0) under the GPU lock; screened = `-gpu_sparse_integrals on`, dense = same binary without it;
n = 1 per cell; dE from the gradient-file header, dG over the Eh/A rows.

| System (atoms, pattern GFN2 / GFN1) | method | screened vs dense dE / dG | screened vs CPU dE / dG | dense vs CPU dE / dG | 2nd screened run |
|---|---|---|---|---|---|
| dioxygen.bent-triplet (2, 100 %) | GFN2 | 0 / 0 (md5 equal) | 0 / 3.3e-12 | - | md5 equal |
| same | GFN1 | 0 / 0 (md5 equal) | 5.3e-15 / 3.5e-14 | - | md5 equal |
| pep20-peptide.gfnff-min (309, 73.4 / 84.9 %) | GFN2 | 0 / 1.4e-13 | 2.8e-13 / 4.7e-11 | 2.8e-13 / 4.7e-11 | md5 equal |
| same | GFN1 | 0 / 2.0e-13 | 0 / 1.2e-10 | 0 / 1.2e-10 | md5 equal |
| peo201-chain (1410, 41.4 %) | GFN2 | 1.4e-12 / 9.2e-14 | 3.2e-12 / 4.0e-11 | 4.6e-12 / 4.0e-11 | not run |

- Screening error (screened vs dense): at most 1.4e-12 Eh, 2.0e-13 Eh/A; CPU gate 1e-8 met in 5 of 5 cells. The GPU-CPU
  difference is that of the dense device path (same size with and without screening).
- Finite differences of the assembled surface: pep20 triplet GFN2 screened, central, h = 1e-4 A, `-scf_threshold 1e-11`, atoms
  0/150/300, 9 components: max |FD - analytic| 1.47e-8 Eh/A (A7 closed-shell pep20: screened 1.68e-8, dense 1.79e-8, the FD noise
  of the method).
- SCF iterations: pep20 GFN2 65 (dense, screened, CPU equal); peo201-chain 50 screened, 50 dense, 52 CPU. Device SCF times of the
  host-driven spin loop: pep20 GFN2 12.4 s screened, 12.5 s dense (CPU 4 threads 6.8 s); peo201-chain 123.6 s screened, 129.4 s
  dense (CPU 205.8 s); one run each on a shared machine.
- Negative control: the base binary with `-gpu_sparse_integrals on` on the pep20 triplet prints the refusal warning and gives the
  gradient file of the dense run (md5 equal).

## Item 2: A1r, G1c, E1 (periodic host-driven loop)

All three are periodic-only on CUDA (`xtb_native.cpp` sets them only with a cell). Component timing (standalone rocBLAS/rocSOLVER
program `spe/bench/bench.hip` in the session scratchpad, not in the repository; random symmetric F, SPD S; median of 5 (n = 1296)
and 3 (n = 2592) calls, ms):

| n | dsygvd | ssygvd | ssyevd | dgemm n^3 | sgemm n^3 | density dgemm (nocc = 2n/3) | complement dsyrk (nv = n/3) | dpotrf + dpotri |
|---|---|---|---|---|---|---|---|---|
| 1296 (216 H2O GFN2) | 192 | 104 | 46 | 11.0 | 0.7 | 7.3 | 2.1 | 21.6 |
| 2592 (432 H2O GFN2) | 759 | 266 | 143 | 82.3 | 4.8 | 56.1 | 17.4 | 112.1 |

### A1r: not ported

The CUDA refinement pays only through cuBLAS FP64 fixed-point emulation (dgemm 41.7 -> 6.7 ms on the RTX 5080); its own `auto`
picks `fp64_direct` on every device without emulation ("refinement with native FP64 GEMMs is slower"). rocBLAS has no FP64 GEMM
emulation, so on ROCm the strategy would always be `fp64_direct` and the port would be dead code. Cost of the refinement with
native FP64 GEMMs from the table: FP32 start (ssyevd + two strsm) plus about 13 n^3 dgemm (4 per pass for two updates and the
acceptance pass, one orthonormalisation) = 143 + 13 x 82 = about 1210 ms at n = 2592 against 759 ms for dsygvd, and
46 + 13 x 11 = 189 ms at n = 1296 against 192 ms. No gain at either size, so no port.

### E1: ported, default on (CUDA default)

- Files: `rocm/xtb_hip_hostloop.hiph` (new, kernel `k_twice_lower`, verbatim CUDA), `rocm/xtb_hip_hostloop_impl.hiph` (new, inside
  `Impl`: `hlComplementWanted`, `hlDensityComplement`, `hlFree`), `rocm/xtb_hip_context.hip` (one branch in `resDensity` after
  the unchanged Cw build, include lines, `hlFree` in `freeResident`, public `setHostLoopDensityComplement`,
  `hostLoopDensityComplementCount`), `rocm/xtb_hip_context.h`, `xtb_hip_method.cpp` (`HipScfBackend::setHostLoopDensityComplement`
  override, a verbosity-2 count line), `CMakeLists.txt` (DEPENDS of the hipcc command). CUDA reference
  `XtbGpuContext::Impl::densityComplement`.
- P = 2 (S^-1 - C_v C_v^T) when every occupied column has occ = 2 exactly, dense storage only; S^-1 by dpotrf + dpotri on a copy
  of dS0, once per S (`s_version`; CUDA reuses its resident factor). Any failure runs the GEMM. `CURCUMA_PBC_GPU_HOST=dcompl`
  keeps the GEMM (the core token).

### G1c: ported, acts only with `-scf_pseudo_diag true` (ROCm default off since WP7)

- Files: `rocm/xtb_hip_hostloop_impl.hiph` (`hlEigensolve`), `rocm/xtb_hip_context.hip` (`resSolve`, `resSolveMP`,
  `resSolveMPResident` call `hlEigensolve`, which is the old full solve unless the host step is armed; the gate of
  `eigensolveStep` accepts a host step; public `setHostLoopPseudo`, which writes the `loop_*` occupation fields and
  `ps_last_dq` as CUDA does), `xtb_hip_method.cpp` (`HipScfBackend::setHostLoopPseudo` override). CUDA reference
  `XtbGpuContext::setHostLoopPseudo` and the `host_pseudo` gate of `eigensolveResidentFock`.
- Difference from CUDA: FP32 steps only. The FP64 pseudo step (`-scf_pseudo_diag_fp64`) stays with the resident loop, whose
  gradient rebuilds W from G_oo; the host-loop gradient builds W from the eigenvalues.
- With a cell the verbosity-2 line "pseudo-diagonalisation (last resident SCF)" counts the host-loop steps since the start of the
  process (no resident loop begin resets it).

### E1/G1c measurements

gfn2, `-sp -gradient -threads 4 -gpu rocm -gpu_strict true`, water216-box.sc-s1 with its cell (L = 18.6258 A) and a 432-H2O cell
made of two copies stacked along z (18.6258 x 18.6258 x 37.2516 A, scratch file, not in the library).

Accuracy, `-scf_threshold 1e-9`, n = 1 per cell (dE Eh / max dG Eh/A):

| comparison | 216 H2O | 432 H2O |
|---|---|---|
| E1 (new default) vs `dcompl` (GEMM) | 6.8e-13 / 4.0e-12 | 1.8e-12 / 8.6e-12 |
| `dcompl` vs base build | md5 equal | md5 equal |
| new default vs CPU | 1.6e-12 / 7.4e-12 | 3.2e-12 / 2.3e-11 |
| G1c (`-scf_pseudo_diag true`) vs new default | 0 / 9.9e-12 | 1.4e-12 / 1.2e-11 |
| `-scf_pseudo_diag true` with `CURCUMA_PBC_GPU_HOST=pseudo` vs new default | md5 equal | md5 equal |
| second run (default; G1c) | md5 equal; md5 equal | md5 equal; md5 equal |

17 SCF iterations in every case at 1e-9; E1 built 17 of 17 densities; G1c took 5 (216) and 4 (432) pseudo steps, none rejected.
GFN1 216 H2O (n = 1): E1 vs `dcompl` 2.1e-12 Eh / 5.5e-12 Eh/A, vs CPU 9.1e-13 / 1.5e-11, G1c vs E1 6.8e-13 / 1.1e-11; 14
iterations, 14 of 14 densities by E1, 2 pseudo steps.

Time, default threshold (9 SCF iterations everywhere), three alternating rounds base / `dcompl` / new / G1c, each under the GPU lock
(other agents shared the CPU), SCF time from the log (ms) and wall time of the whole `-sp -gradient` run (s), median of 3:

| variant | 216 SCF | 216 wall | 432 SCF | 432 wall |
|---|---|---|---|---|
| base build | 1305 (1274-1315) | 2.28 | 4424 (4365-4435) | 7.79 |
| new, `dcompl` (E1 off) | 1312 (1308-1322) | 2.30 | 4432 (4418-4435) | 7.70 |
| new default (E1) | 1283 (1273-1289) | 2.23 | 4164 (4149-4165) | 7.39 |
| new + `-scf_pseudo_diag true` (E1 + G1c) | 812 (802-816) | 1.81 | 3721 (3714-3723) | 6.99 |

- E1 saves 22 ms (1.7 %) and 260 ms (5.9 %) per SCF; the ranges of 432 do not overlap, those of 216 overlap (the 216 gain is within the spread). G1c on top saves
  471 ms (37 %) and 443 ms (11 %). Whether `-scf_pseudo_diag` should default on for ROCm is a decision outside this item (WP7 left
  it off because of extra iterations with `-scf_precond frag`); not changed here.

## Gates without the switches

- New binary against the base, `-gpu rocm` at the default threshold: 14 of 14 gradient files md5-identical (caffeine GPU and CPU,
  O2 spin 2, Pt7 spin 8 at 2000 K, pep20 closed shell dense and screened (`-gpu_sparse_integrals on`), pep20 triplet dense; each
  GFN2 and GFN1). Screened closed-shell runs are untouched by the change (md5 equal to the base).
- ctest (final build, under the GPU lock, `-j1`): `-R "cli_sqm_2[34]|pbc_xtb_gpu_parity_rocm|sqm_val_rocm|structinv_rocm|cli_sqm_16_spin"`
  79 of 79 passed, including `pbc_xtb_gpu_parity_rocm` (now with E1 in its periodic runs), `cli_sqm_23_spin_polarization_rocm` and
  `cli_sqm_24` with the new case (e) (screened spin vs dense 0 Eh, 9.1e-14 Eh/A at its 1e-9 threshold). Not run: the full ctest suite.

## Not tested

- Spin polarisation with screened storage: geometry optimisation and MD (a new pattern per geometry), solvation, `-gpu_sparse_integrals
  auto` (only `on`), d-shell elements beyond Pt7 at 100 % pattern, systems where the dense spin run does not fit the 16 GB card (memory
  gain not measured), charged open shells.
- E1/G1c: open-shell or charged cells (the core does not arm them), MD with a cell (one geometry per run only), GFN1 timing, cells
  other than water, `-xtb.rocm_cached_cholesky` / `-xtb.rocm_density_syrk` together with E1, failure paths of `hlDensityComplement`
  (no fault injection), any GPU other than gfx1200.
- A1r: only the component benchmark, no implementation.

## Not implemented relative to CUDA

- A1r (reason above). E1 builds S^-1 from its own factor instead of the resident L. G1c without FP64 pseudo steps on the host loop.
- Spin + screened storage keeps P, Pm and W on the pattern; CUDA keeps them dense (a different layout, not a missing feature).

## Proposed central-file lines (for the coordinator)

- TODO.md: drop "spin polarisation refused (dense)" from the A7 line; add "ROCm: `-scf_pseudo_diag` default for ROCm (G1c makes
  the periodic SCF 37 % / 11 % faster on 216 / 432 H2O, WP7 keeps it off), decide (docs/ROCM_PARITY_SPE.md)."
- AIChangelog.md: "Nov 25, 2026: ROCm parity SPE: spin polarisation with -gpu_sparse_integrals on the device (P, Pm, W' on the
  pattern; screened vs dense <= 1.4e-12 Eh / 2.0e-13 Eh/A, vs CPU <= 3.2e-12 Eh); E1 density complement and G1c host-loop pseudo
  step ported (periodic SCF -1.7 % / -5.9 %, with -scf_pseudo_diag -37 % / -11 % on 216 / 432 H2O); A1r not ported (no FP64 GEMM
  emulation in rocBLAS, refinement 1.6x slower than dsygvd at n = 2592) (docs/ROCM_PARITY_SPE.md)."
- GPU_PARITY_CUDA_ROCM.md: row "sparse integrals": remove "spin polarisation with it refused (dense)"; rows A1r / G1c / E1: ROCm
  "no (measured, no FP64 emulation)" / "yes, with -scf_pseudo_diag (ROCm default off)" / "yes (default on)".
- GPU_PERIODIC_GAPS.md (ROCm part): A1r not planned on ROCm, G1c and E1 closed.
- qm_methods/CLAUDE.md: in the "ROCm screened integral storage" bullet replace "spin polarisation keeps dense storage" by "spin
  polarisation on the pattern (SPE)"; in the "Periodic cell on ROCm" bullet replace "not ported: ..., A1r, G1c, E1" by "not ported:
  ..., A1r (no FP64 emulation); G1c, E1 since SPE".
- `scripts/gpu_parity_known.json`: `densityComplement`, `setHostLoopPseudo`, `setHostLoopDensityComplement`, `k_twice_lower` now
  exist on ROCm (names `hlDensityComplement`, `k_twice_lower`); `refineEigenFP64`, `selectEigenStrategy`, `setEigenRefinement` and
  the `k_rf_*` kernels stay CUDA-only (reason: no FP64 emulation). Not checked whether the script lists them.
- Plugin ABI: no change. The `GpuScfBackend` virtuals (`setHostLoopPseudo`, `setHostLoopDensityComplement`) already existed;
  the new `XtbHipContext` functions live inside `libcurcuma_rocm.so`; the `curcuma` executable is md5-identical to the base.
