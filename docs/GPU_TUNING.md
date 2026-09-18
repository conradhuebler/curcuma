# GPU and Large-System Tuning Options

> 🤖 AI-generated options and measurements, machine-tested only. Human production testing pending.

Reference for manual tweaking of multi-GPU runs and large systems (thousands of atoms). Every option below lists its default, what it changes, when to touch it, and what was measured. Measurements: 4x RTX A4500 (20 GB, Ampere sm_86, PCIe), 36 cores, CUDA 13.3, Sep 2026, unless noted. Background and history: [MULTI_GPU.md](MULTI_GPU.md).

CLI flags are global keys (`-gpu_device 2`) unless written with a scope (`-gfnff.gpu_coulomb_implicit false`). Environment variables are read at start-up.

## 1. Device selection and batch distribution

| Option | Default | Effect |
|---|---|---|
| `-gpu cuda\|rocm\|vulkan\|auto\|none` | `none` | Backend. `auto` = first available plugin. |
| `-gpu_device N` | -1 (runtime default, device 0) | Pins a single calculation to device N. Indices are after `CUDA_VISIBLE_DEVICES` / `ROCR_VISIBLE_DEVICES` (SLURM sets these). Invalid index: error + CPU fallback. |
| `-gpu_devices 0,2` | all visible | Device pool for batch workers. |
| `-gpu_workers_per_device N` | 1 | Batch workers sharing one device. |

- **Batch workers** (`-sp`/`-opt` on multi-XYZ, ConfSearch, Hessian) lease a device slot each; surplus workers wait. The pool is active when more than one device is visible, or when `-gpu_devices`/`-gpu_workers_per_device` is given.
- **Host threads per GPU worker** = cores / GPU slots (CPU batch workers stay single-threaded).
- **When to change `-gpu_workers_per_device`**: small/medium molecules where one context does not saturate the card. Measured, 8x polymer (1410 atoms) GFN2: 4 GPUs 28.5 s, 8 slots (2 per GPU) 25 s. Watch device memory: each worker holds its own context.
- **When GPU batching does not pay**: small molecules. 40 conformers of 90 atoms: 4 GPUs 2.0 s vs CPU with 36 workers 0.58 s.
- `curcuma -methods` lists the devices a plugin sees.

## 2. Native GFN1/GFN2 on the GPU (CUDA)

| Option | Default | Effect / when to change |
|---|---|---|
| `-gpu_memory_check true\|false` | true | Refuses a basis whose estimated device memory exceeds free memory (warning + CPU run). Set false only if the estimate is too pessimistic for a card that just fits. |
| `-gpu_sparse_integrals auto\|on\|off` | auto | Screened pair storage for S, H0 and the GFN2 multipole integrals. `auto` = when fewer than 50 % of AO pairs survive the distance screen (integrals < 1e-20 dropped). `on` forces it (used for validation: identical energies on small molecules), `off` forces dense nao^2 storage. polymer_2x: 6 % of pairs, device memory ~50 -> 17 GB. |
| `-scf_fp32_threshold X` | 1e-5 on the GPU, 1e-3 on the CPU | Mixed precision: FP32 eigensolves until max\|dq\| < X, then FP64. Larger X = more FP64 iterations (slower on FP64-weak consumer cards, closer CPU gradients). polymer (1410 atoms): 1e-3 -> 4 FP64 solves, SCF 6.7 s, gradient vs CPU 2.6e-6; 1e-5 -> 1 FP64 solve, 4.5 s, 7.0e-6 Eh/A. Below ~3e-6 the FP32 noise floor adds iterations. |
| `-scf_mixed_precision true\|false` | **false on a device with full-rate FP64** (A100, H100/H200/GH200, B200, V100 - compute capability 6.0/7.0/8.0/9.x/10.x), true elsewhere | Mixed precision is a consumer-GPU trick: it pays where FP64 runs at 1/32 to 1/64 of FP32 (A4500: FP64 iteration 38 s vs FP32 9.4 s) and hurts where FP64 is half of FP32. Measured on an H200 NVL, polymer_2x: **FP32 iteration 5.29 s, FP64 3.36 s** - and the FP32 phase converged to a fixed point **1.1 kcal/mol off** (it reported max\|dq\| 7.1e-6 while the true residual was 9.4e-3), which the SCF only discovered through its FP64 check, at the cost of 42 iterations instead of 15. |
| `-scf_threshold 1e-8` | 1e-5 | Tight SCF. Needed for GPU gradients that must agree with the CPU beyond ~1e-4 Eh/A (the loose default lands the device Broyden at a different point inside the tolerance). |

### Multi-GPU eigensolve (one large molecule on several GPUs)

The per-iteration full-spectrum eigensolve of the device-resident SCF can be spread over several GPUs; integrals, Fock build, density, mixing and gradient stay on the calculation's own device (`-gpu_device`, default 0). Off by default.

| Option | Default | Effect / when to change |
|---|---|---|
| `-gpu_eigensolver_devices all\|0,1,2,3\|none` | **all visible devices** when more than one is visible and this calculation does not run inside a batch worker, else off | Devices for the eigensolve (at least two); may include the calculation's own device. `none` switches it off. Batch workers keep it off: each of them already owns one device. |
| `-gpu_eigensolver_backend auto\|mp\|mg` | auto | `mp` = cuSOLVERMp over NCCL, `mg` = cusolverMg (CUDA toolkit, deprecated in CUDA 13). `auto` tries mp, then mg. |
| `-gpu_eigensolver_min_nao N` | 4000 | Below N basis functions the single-GPU solver is used: the column scatter/gather and NCCL collectives cost more than they save for small matrices. This gate is what keeps the default harmless for everyday molecules. |
| `-gpu_eigensolver_block N` | 128 | Column block size of the 1 x ndev distribution. Standalone n = 15444, 4 GPUs: nb 64/128/256 -> FP64 21.9/20.7/20.9 s, FP32 5.50/5.24/5.36 s. |
| `-gpu_eigensolver_fp32 true\|false` | true | Also distribute the FP32 iterations of the mixed-precision SCF (`mp` only; `mg` never gets FP32 solves, see below). Set false to keep FP32 on the calculation's own GPU. |

- On a failure the SCF continues with the single-GPU solver; at verbosity >= 1 a warning names the reason, verbosity 2 reports backend, GPUs and number of solves after each SCF.
- The solvers live in `libcurcuma_cuda_mgpu.so` next to `libcurcuma_cuda.so` and are loaded on first use, so a host without cuSOLVERMp/NCCL still runs all single-GPU calculations.
- What is distributed: FP64 iterations run the whole generalized solve on the GPUs (reduction `sygst`, `syevd`, back-transform `trsm`; the Cholesky factor is distributed once per geometry). FP32 iterations distribute only `syevd`: the cuSOLVERMp FP32 reduction is not faster than one GPU (0.56 vs 0.57 s) and the cuBLASMp FP32 back-transform is slower (1.74 vs 0.38 s).
- Memory: the eigensolver workspaces and the per-GPU solver buffers are released as soon as the SCF has converged, before the post-SCF phase and the gradient.
- **Every distributed solve is verified** (`-gpu_eigensolver_verify`, default true). A random vector is pushed through the matrix twice - once directly and once through the returned eigenpairs (for the generalized path with the metric S = L L^T) - and the eigenvalues are checked for a genuinely permuted spectrum. A backend that fails is dropped for that precision, the input is restored and the solve returns to the calculation's own device, with a warning naming the residual. Measured cost: polymer_2x 183.3 s with verification vs 184.9 s without, i.e. inside the noise; the gradients agree to 1.2e-14.
- Why it is not a one-off check: **cusolverMg passes the first solve and degrades afterwards** (polymer, nao 3222: solve 1 residual 1.4e-6, solve 2 8.7e-2 with buffer reuse, 1.2e-3 with everything rebuilt per call). A single gate at the start let a run converge to an energy **1.65 kcal/mol wrong**; with per-solve verification the same run is rejected and lands on the correct energy. cusolverMg is therefore useful only where it verifies: at nao 558 it passes every solve, at 3222 its FP32 is rejected and only its FP64 is used.
- **What FP32 costs on a datacenter GPU** (operator, H200 NVL, Sep 18, 2026): besides being slower
  per iteration there, the FP32 phase can converge to the wrong fixed point at this size. curcuma
  never accepts convergence on an FP32 step, so the result stayed correct
  (-11784.87804452 Eh), but the run needed 42 iterations. Two things now prevent that: mixed
  precision defaults to OFF on such a device (see above), and if an FP64 iteration finds a residual
  more than 10x above what FP32 last claimed, the SCF says so and stays in FP64 for the rest.
- **Second machine, same verdict** (operator, Sep 18, 2026, 2x RTX PRO 5000 Blackwell, a build
  where only cusolverMg was found): `cusolverMg: its FP32 eigenpairs failed verification here
  (relative residual 0.003513, 0 eigenvalues out of ascending order); FP32 iterations stay on
  device 0 and only FP64 is distributed`. The run then converged in 12 iterations to
  -11784.87804452 Eh, i.e. exactly the reference energy, in 91 s (FP32 iteration ~4.0 s on one
  Blackwell, the single FP64 one 25.3 s). So Mg's FP32 defect is not specific to the A4500 box,
  and the fallback does what it is supposed to do on hardware we cannot test here.
- `CURCUMA_GPU_EIG_VERIFY_ALWAYS=1` prints the residual of every verification, `CURCUMA_GPU_EIG_CORRUPT=1` deliberately damages the eigenvectors once to prove the check reacts (it reports residual 0.43 and falls back, and the run still gives the correct energy).
- Results: energies identical to the single-GPU run at the printed precision; with `-scf_threshold 1e-9` energies and gradients agree to <= 1.2e-9 (complex, 231 atoms, gfn1 and gfn2, mp and mg). At the loose default threshold the gradients differ by up to 1e-5 Eh/A, because the FP32 iterations take a slightly different path inside the tolerance (same as GPU vs CPU).

Measured, polymer_2x GFN2 single point + gradient (7320 atoms, nao 15444), 4x RTX A4500 (PCIe, P2P), SCF on device 0, `CURCUMA_GPU_PROFILE=1`, energy -11784.87804452 Eh in every row:

| Configuration | Wall | SCF iterations (FP32+FP64) | per FP32 iteration | per FP64 iteration | peak device 0 | peak other GPUs |
|---|---|---|---|---|---|---|
| single GPU | 419 s | 14+2 | 12.2 s | 85 s | 17.9 GB | - |
| `-gpu_eigensolver_devices 0,1` (syevd only, before the FP64 generalized path) | 327 s | 15+1 | 11.7 s | ~67 s | 17.7 GB | 5.4 GB |
| `-gpu_eigensolver_devices all` (syevd only, same build) | 260 s | 13+1 | ~10.6 s | ~49 s | 15.4 GB | 3.2 GB |
| `-gpu_eigensolver_devices all` (current: FP64 generalized, early release) | 242 s | 13+1 | 9.4 s | 38.2 s | 13.5 GB | 5.0 GB |
| `-gpu_eigensolver_devices 1,2,3` (device 0 not in the solver) | 244 s | 11+1 | 10.3 s | 48.5 s | 12.3 GB | 7.0 GB |

- Iteration counts differ between rows because the FP32 path lands at slightly different points; compare the per-iteration columns.
- Both features are ON by default from 4000 basis functions up when several GPUs are visible (operator decision, Sep 2026); the rows above with fewer devices need the explicit option.
- Measure on idle GPUs: the `mem` column of the profiler is device-wide, so a concurrent job (e.g. `ctest -L gpu`) inflates it and slows the run - one comparison in this file was repeated for that reason.
- FP64 iteration split (all 4 GPUs): sygst 3.8 s (single GPU: two trsm 21 s), syevd 19.5 s (49 s), trsm back-transform 10.9 s (10.5 s). FP32: syevd 4.6 s (7.5 s).
- Including device 0 is ~10 % faster per iteration; excluding it lowers device 0 by ~1.2 GB and raises the others by ~2 GB. Use `1,2,3` when device 0 is the one that does not fit.
- The pattern density is distributed separately, see the next section.
- Gradients at the default `scf_threshold` differ from the single-GPU run by up to 4.6e-4 Eh/A (different FP32 path inside the loose tolerance, same as GPU vs CPU); with `-scf_threshold 1e-9` they agree to <= 1.2e-9 (complex).
### Checklist: making ONE calculation use several GPUs

What is distributed inside a single point / optimisation step:

| part | device | share of an A4500 run (polymer_2x, 4 GPUs) |
|---|---|---|
| eigensolve (FP32 and FP64 iterations) | all listed GPUs | 86 of 194 s |
| screened-pattern density | all listed GPUs | 10 s |
| integrals, Fock build, potential, charges, gradient, post-SCF | the calculation's own device only | ~98 s |

So the ceiling for one calculation is set by the part that stays on one device - on the A4500 box the whole single point went 419 -> 194 s (2.16x on 4 GPUs), not 4x, and that is the honest expectation.

**Prerequisites, in the order they bite:**

1. **cuSOLVERMp, cuBLASMp and NCCL at build time.** Without them only `cusolverMg` (CUDA toolkit) is available, and curcuma gives Mg **FP64 solves only** - its FP32 eigenvectors are wrong at large n (measured: polymer_2x diverged from the first iteration). Since mixed precision is the default, that means roughly one of twenty iterations is distributed and the run looks single-GPU. Verify at configure time:
   `cmake .. -DCUSOLVERMP_ROOT=... -DCUBLASMP_ROOT=... -DNCCL_ROOT=...` must print
   `-- Multi-GPU eigensolver: cuSOLVERMp <path>, NCCL <path>`.
2. **The devices must be listed as INDICES.** `-gpu_devices 2` means "device number 2", not "two GPUs" - use `-gpu_devices all` or `0,1`. Same for `-gpu_eigensolver_devices` / `-gpu_density_devices` (`all`, `solver` or a list).
3. **Check what actually ran**, at `-verbosity 2`, after the SCF:
   `GPU multi-GPU eigensolver: cuSOLVERMp on 2 GPUs, 20 solves` - the backend name and the solve count are the test. `cusolverMg on 2 GPUs, 1 solves` means prerequisite 1 is missing.
   `GPU distributed density: pattern density on devices [0,1], 21 steps`.
4. **The libraries must be found at run time too**: their directories are baked into the RPATH of `libcurcuma_cuda_mgpu.so`, so a compute node needs the same paths (shared file system) or `LD_LIBRARY_PATH`.

Note that the status lines for the density and for a successful eigensolve verification are
`info` (verbosity >= 2); only a failed verification is a warning. A run at the default verbosity
that prints nothing about the density did not necessarily skip it - use `-verbosity 2`.

**Measuring it on a new machine** (nothing here is measured on H200 yet - PCIe A4500 numbers do not transfer to NVLink):

```
CURCUMA_GPU_PROFILE=1 curcuma -sp big.xyz -method gfn2 -gpu cuda -verbosity 2 \
    -gpu_eigensolver_devices none -gpu_density_devices none      # baseline, one GPU
CURCUMA_GPU_PROFILE=1 curcuma -sp big.xyz -method gfn2 -gpu cuda -verbosity 2   # defaults
```

The profile's `eig FP32 / eig FP64 / density P` rows say how much of the run is distributable at all; if they are not the majority, more GPUs cannot help much and the answer is a bigger share on the device instead (integrals, Fock, gradient are all single-device today).

### Distributed pattern density

`P(r,c) = sum_k Cw(r,k) C(c,k)` is a sum over the occupied columns, so each GPU can evaluate the whole screened pattern over its own slice of columns and the partials are added on the calculation's device. This is exact, not an approximation. Only a column slice of C (n x kn) and one partial pattern array travel per SCF step; the pattern indices are uploaded once per geometry.

| Option | Default | Effect / when to change |
|---|---|---|
| `-gpu_density_devices all\|0,1,2\|solver\|none` | **all visible devices** under the same condition as the eigensolve above, else off | Helper devices for the density (the calculation's own device is always one of the workers and is skipped in the list). `solver` reuses `-gpu_eigensolver_devices`, `none` switches it off. Needs the screened storage (`-gpu_sparse_integrals`, automatic for large systems). |
| `-gpu_density_min_nao N` | 4000 | Below N basis functions the single-device kernel is used (polymer, nao 2975: the split is 10 % slower). |
| `-gpu_multipole_otf auto\|on\|off` | auto | GFN2 multipole interaction matrices: stored, or rebuilt per iteration in the kernel. `auto` rebuilds above ~1 GB of matrices (~2700 atoms), which is what makes 7320 atoms fit a 20 GB card (18 nat^2 doubles = 7.7 GB, plus the same again as host upload copies). Costs ~50 ms/iteration at 1410 atoms, where storing is the better choice (measured: 9.13 s stored vs 9.87 s on the fly). Energies identical either way. Replaces `CURCUMA_GPU_MP_OTF`, which still works and wins over the flag. |

- Measured, polymer_2x GFN2 on 4x RTX A4500 (identical with the defaults and with explicit `-gpu_eigensolver_devices all -gpu_density_devices solver`, both 194 s): density 3.35 -> **0.85 s per SCF step** (46.9 -> 10.2 s total), FP32 iteration 9.4 -> 6.9 s, wall **242 -> 194 s**, energy identical. Device 0 peak 13.5 -> 13.9 GB, helpers 5.0 -> 5.8 GB.
- Too small to pay: polymer (1410 atoms, nao 2975) goes 50 -> 55 ms per step, i.e. the transfers cost more than the split saves. The split needs at least 16 columns per worker and falls back to the single device otherwise.
- On any failure it falls back to the single-device kernel for the rest of the run and warns (at verbosity >= 1).
- Rejected alternatives for the same kernel: a row-major variant (contiguous rows, occupations applied in the kernel) was **8x slower** (4.69 vs 0.60 s over 12 steps on polymer); a banded GEMM does not help because the atom order of polymer_2x is not spatially local - a column block spans 93-99 % of the rows although only 7.9 % of the atom pairs are inside the cutoff.

- Not measured: `mg` end to end on polymer_2x in FP64 (standalone FP64 26.4 s vs `mp` 20.7 s on 4 GPUs), GFN1 on polymer_2x, NVLink systems (H200).

Environment variables:

| Variable | Effect |
|---|---|
| `CURCUMA_GPU_PROFILE=1` | Stream-synchronised per-phase device timings (integrals, potential, Fock, reduce/syevd/back-transform per precision, density, charges, Broyden) plus the post-SCF host phases. Adds synchronisation cost; for diagnosis only. |
| `CURCUMA_GPU_MP_OTF=0\|1` | The old form of `-gpu_multipole_otf`. Still honoured and still wins over the flag, for reproducing older runs. |

Automatic large-system behaviour (no option; triggered above ~5000 basis functions on the CUDA device-resident path):
- the host does not build the dense dipole/quadrupole integrals or the post-SCF host Fock matrix (host fallbacks build them on demand; disabled for `-d4_charge_source cpscf`);
- a single point without gradient keeps P and C on the device and takes the SCC energies from the last device step (`MolecularOrbitals()` is empty in that case).

Reference numbers, polymer_2x GFN2 (7320 atoms, nao 15444), one A4500: 375 s, 16.1 GB device, 15.6 GB host, E identical to the CPU reference (4468 s). One H200 (operator, before these changes): 110 s.

## 3. GFN-FF on the GPU (CUDA)

| Option | Default | Effect / when to change |
|---|---|---|
| `-gfnff.coulomb_implicit true\|false` | true | **CPU**: evaluate the N^2/2 Coulomb pairs on the fly from the per-atom EEQ charges and alpeeq instead of building and storing the pair list. The stored list is 128 bytes per pair - 3.4 GB and ~0.5 s of pure write bandwidth at 7320 atoms, and threading that loop changes nothing (measured). polymer_2x single point, 36 threads: wall 4.49 -> 3.93 s, host RSS 6.66 -> 3.72 GB, pair setup 496 -> 0.1 ms, and the Coulomb evaluation itself 1733 -> 1096 ms (summed over threads) because it no longer streams the list. MD (polymer, 100 steps) is unchanged at 7.7 s, small molecules unchanged. Energies identical, gradients within 4e-16. Set false for the stored list. Not used with `eeq_distance_cutoff > 0`. |
| `-gfnff.gpu_coulomb_implicit true\|false` | true | The device enumerates all Coulomb pairs itself instead of a host-built N^2/2 list. polymer_2x: -1.5 s setup, -6 GB host memory, energy identical. Ignored with `eeq_distance_cutoff > 0`. Set false only to compare against the pair-list path. |
| `-gfnff.gpu_disp_pairs_on_device true\|false` | false | Builds the D4 dispersion pair list on the device. polymer_2x: dispersion pairs 630 -> 32 ms, wall 6.8 -> 5.6 s, energy identical, gradient 5e-16; device peak higher (list lives on the device). Kept off by default until tested on more systems - a good first switch for large systems. |
| `-gfnff.eeq_mixed_precision true` | false | FP32 EEQ factorisation + FP64 refinement (few-fragment paths). For FP64-weak cards; measure per card. |
| `-gfnff.eeq_rmsd_threshold X` (Bohr) | 0 | MD/opt: reuse the EEQ Cholesky factor while the geometry moved less than X per atom (0 = refactor every step, exact reference behaviour). |
| `-gfnff.gpu_block_size N` | 0 (auto) | Kernel launch block size. |
| `-gfnff.gpu_cn_pair_regen*`, `gpu_cn_pair_cutoff_factor` | see `curcuma -help gfnff` | CN-derivative pair list refresh during MD/opt. |

Environment variables:

| Variable | Effect |
|---|---|
| `CURCUMA_GFNFF_PROFILE=1` | Setup phase timings summed over both q-loop topology passes (the verbosity-2 report shows only the last pass) and the GPU upload phases. |

Threading note: GFN-FF topology and parameter generation use OpenMP with the GFN-FF thread budget (`-threads N`); inside batch workers they use the batch's per-worker share. Before Sep 2026 they silently ran on one thread.

Reference numbers, polymer_2x GFN-FF single point (7320 atoms, 1502 fragments), one A4500: 6.8 s (5.6 s with `gpu_disp_pairs_on_device`), 4.8 GB host. CPU with 36 threads: ~6 s.

## 4. CPU knobs of the native GFN1/GFN2 SCF

These apply with and without a GPU (the eigensolve and the reduction run on the host
whenever the device-resident path is not active). Defaults are the measured optimum on the
36-core development box; the optimum is machine-dependent, which is what section 5 is for.

| Option | Default | Effect / when to change |
|---|---|---|
| `-threads N` | 1 | Intra-molecule threads, including the eigensolve. The dense divide-and-conquer eigensolve is memory-bandwidth-bound, so more threads than memory channels can lose: polymer (nao 3222) on 36 cores measured 8 threads 27.1 s, **16 threads 23.8 s**, 24 threads 26.5 s, 36 threads 28.6 s. |
| `-eigensolver_max_threads N` | 0 (follow `-threads`) | Caps the eigensolve alone. Use it when the integrals and the gradient want all cores but the eigensolve peaks earlier. Replaces `CURCUMA_EIG_MAX_THREADS`, which still works and still wins. |
| `-scf_reduce auto\|sygst\|trsm` | auto | Reduction of F to standard form with the cached Cholesky factor of S. `sygst` does half the flops but threads poorly, `trsm` does twice the flops as two BLAS3 solves and keeps scaling; `auto` takes `trsm` from `-scf_reduce_threads` up. Measured n=3222 on 36 cores (OpenMP OpenBLAS), dsygst vs 2x dtrsm: 620/900 ms at 1 thread, 195/192 at 8, 190/**123** at 16, 308/**156** at 36. The two agree to 8e-15 elementwise. |
| `-scf_reduce_threads N` | 8 | Where `auto` switches. Measure with `sygst` vs `trsm` at your thread count before changing it. |
| `-scf_mixed_precision true\|false` | true on CPU and on consumer GPUs, **false on full-rate-FP64 GPUs** (compute capability 6.0/7.0/8.0/9.x/10.x, e.g. H100/H200/Blackwell) | FP32 eigensolves far from convergence, FP64 near it. On an RTX A4500 turning it off costs a factor 1.9 (polymer, 9.2 -> 17.7 s); on an H200 it *is* the slower path (5.29 s vs 3.36 s per iteration) and can converge to a false fixed point. An explicit flag always wins over the device default. |
| `-scf_fp32_threshold X` | 1e-3 (CPU), 1e-5 (GPU) | Where the FP32 phase ends. Smaller = longer in FP32 = faster and less safe. Must stay above `-scf_threshold`. |
| `-scf_fp32_stall_patience N` | 3 | FP32 iterations without real progress (max\|dq\| not improving by 30 %) before the rest of the SCF runs in FP64. FP32 eigenvectors carry ~1e-7 noise, which puts a floor under dq on a large system; without the guard the SCF hovers there. 0 disables. |
| `-scf_fp32_false_fixpoint_factor X` | 10 | If an FP64 iteration reports a residual more than X times what FP32 last claimed, FP32 converged to a false fixed point and the rest runs in FP64. Measured on an H200: FP32 claimed 7.1e-6 at an energy 1.1 kcal/mol off while the truth was 9.4e-3. 0 disables. |
| `-eigensolver mkl\|native\|purify\|lobpcg` | mkl | `mkl` (LAPACK dsyevd) is the fastest on CPU; the others are the GPU-portable / research paths. |

What is deliberately NOT a flag, and why:

| Variable | Why it stays an environment variable |
|---|---|
| `CURCUMA_EIG_PROFILE`, `CURCUMA_GPU_PROFILE`, `CURCUMA_XTB_REDUCE_PROBE`, `CURCUMA_GPU_EIG_VERIFY_ALWAYS`, `CURCUMA_GPU_EIG_CORRUPT` | Diagnostics. They print or inject, they do not tune. |
| `CURCUMA_EIG_TRED2=eigen\|blocked\|scalar` | Only reachable with `-eigensolver native`, itself an opt-in research path. All three are bit-identical; see [SQM_EIGENSOLVE_GPU.md](SQM_EIGENSOLVE_GPU.md). |
| `CURCUMA_GPU_EIG_MG_REUSE=1` | Restores a known-wrong behaviour (cusolverMg returns bad eigenvectors on a reused descriptor) for testing. Making it a documented flag would advertise a setting nobody should use. |
| `CURCUMA_GFNFF_GPU_RESIDENT_HBQ=1` | A deliberate deviation from the reference (live instead of frozen H-bond charges), i.e. a physics experiment, not a performance knob. |

## 5. Measuring it on your own machine

`scripts/tuning_sweep.py` scans these knobs on the hardware it runs on and prints a command
line. Every value it tries is a CLI flag, so nothing has to be rebuilt:

```bash
python scripts/tuning_sweep.py test_cases/molecules/larger/polymer.xyz --method gfn2
python scripts/tuning_sweep.py big.xyz --method gfn2 --gpu cuda --repeats 3 --json sweep.json
python scripts/tuning_sweep.py big.xyz --gpu cuda --knobs gpu_eigensolver_devices,threads
python scripts/tuning_sweep.py big.xyz --gpu cuda --dry-run     # what would be measured
```

`--dry-run` lists the settings and the number of runs without starting any (15 knobs / 92
runs for gfn2 on a 4-GPU 36-core box), and after the baseline the script prints a rough
total time, so an unattended cluster run can be sized before it is submitted.

It runs a baseline, then each knob on its own, then the winners together (a combination can
be slower than its parts, so that is measured rather than assumed), and it **checks every
run's energy against the baseline** - a setting that moves the energy by more than
`--energy-tol` kcal/mol is printed as SUSPECT and never recommended. `--json` keeps the full
table including a hardware note, so a result from the cluster can be read here.

### A full sweep, measured

polymer_2x (7320 atoms, nao 15444), gfn2, 4x RTX A4500, 45 runs (`--repeats 1`), about 3 h.
Baseline 143.95 s, 12 SCF iterations, E = -11784.87804452 Eh. **All 43 runs that completed
returned that same energy to the last printed digit**, including the distributed eigensolve
and the distributed density.

| Knob | values (wall time) |
|---|---|
| `-threads` | 1 / 2 / 4 / 8 / 16 / 32 / 36: **144 s throughout** |
| `-eigensolver_max_threads` | 0 / 4 / 8 / 16: 144 s |
| `-scf_reduce` | auto 143, sygst 144, trsm 144 |
| `-scf_mixed_precision` | true 145, **false 420 (0.34x)** |
| `-scf_fp32_threshold` | 1e-3 **223**, 1e-4 **166**, 1e-5 **144** |
| `-scf_guess` | eeq 144, h0 157 (14 instead of 12 iterations) |
| `-gpu_eigensolver_devices` | all 144, **none 220 (1.53x for the split)**, 0,1 207 |
| `-gpu_eigensolver_backend` | auto 144, mp 144, **mg 2226 (0.06x), 37 iterations** |
| `-gpu_eigensolver_block` | 64 183, 128 145, 256 **140**, 512 148 |
| `-gpu_eigensolver_fp32` | true 145, false 176 |
| `-gpu_eigensolver_verify` | true 144.26, false 144.32 - **free** |
| `-gpu_density_devices` | all 144, **none 193 (1.34x for the split)** |
| `-gpu_multipole_otf` | auto 144, on 144, off: does not fit the card (see below) |
| `-gpu_sparse_integrals` | auto 144, on 144, off: no result within 40 min |

What this says beyond the individual numbers:

- **Both multi-GPU defaults earn their keep at this size** (1.53x eigensolve, 1.34x density)
  and only at this size - at 1410 atoms the same knobs measured flat, which is what the
  4000-basis-function gate is for.
- **`-threads` is irrelevant on the GPU path here.** The integrals are built on the device
  and the host multipole integrals are never built at all ("deferred, 17.2 GB"), so the whole
  setup phase measures 20.63 s at 1 thread and 20.53 s at 36, inside a run whose SCF alone is
  113 s on the device. This says nothing against threading on the CPU path, where the same method on
  polymer measured 44.2 -> 18.8 s.
- **Mixed precision is the one large lever on a consumer/workstation card** (2.9x), and
  `-scf_fp32_threshold` is a real dial rather than a safety switch: 1e-3 costs 55 % over the
  GPU default of 1e-5.
- **`mg` (cusolverMg) is a memory fallback, never a performance option**: 15x slower and 37
  instead of 12 iterations, because the per-solve verification keeps rejecting its
  eigenvectors. The energy is still exact, which is the verification doing its job.
- **Two `off` settings simply cannot run this system**, which is the empirical case for both
  `auto` defaults. `-gpu_multipole_otf off` is instructive about where the limit sits: the 7.2 GB
  of stored matrices still fit the 20 GB card, and the resident SCF loop's own buffers then do
  not, so the failure arrives one step later. It now says so ("the stored GFN2 multipole
  interaction matrices hold 7.2 GB of device memory - rebuild them per iteration instead")
  instead of the bare "GPU resident SCF step failed at iteration 0".

Two things to know before trusting a sweep:
- Run it the way the production job runs (same allocation, same `CUDA_VISIBLE_DEVICES`). The
  multi-GPU knobs in particular only mean something for what the process actually sees.
- The result is for that structure, method and mode. A knob that wins on a 7000-atom single
  point can lose on a 200-atom optimisation - the multi-GPU paths for instance do nothing
  below 4000 basis functions by design, and a sweep on a small molecule will say so.

## 6. Build

### Did this build get cuSOLVERMp, or only the fallback?

The distributed eigensolver is **optional**, and a build without it is easy to miss - that is
how the Sep 2026 H200 run ended up measuring the cusolverMg fallback (15x slower than a single
GPU on polymer_2x) while looking like "multi-GPU does not help". Three ways to check, in the
order they are cheap:

```bash
cmake .. ... 2>&1 | grep "=== curcuma multi-GPU eigensolver"   # one grep-able summary line
ldd release/libcurcuma_cuda_mgpu.so | grep -Ei "cusolvermp|cublasmp|nccl"   # the built artifact
curcuma -methods            # "multi-GPU eigensolver: mp (cuSOLVERMp)" vs "mg (cusolverMg fallback)"
```

The configure summary is one of
`cuSOLVERMp + NCCL`, `cusolverMg only - optional cuSOLVERMp missing`, or
`none - one molecule stays on one GPU (everything else works)`, and the "not found" branch names
each missing variable individually plus the pip/HPC-SDK command that supplies it.

**What is lost without them**: only the split of ONE molecule's eigensolve. Single-GPU runs,
batch distribution over several GPUs (`-gpu_devices`) and the distributed density
(`-gpu_density_devices`) are unaffected - and where only `mg` is available,
`-gpu_eigensolver_devices none` is faster than using it.

`-DCURCUMA_REQUIRE_MULTI_GPU_EIGENSOLVER=ON` turns the "not found" case into a configure error,
for machines where the distributed eigensolve is the point of the build. Default is OFF.

- Multi-GPU eigensolver (`CURCUMA_MULTI_GPU_EIGENSOLVER`, default ON): cusolverMg is found in the CUDA toolkit. cuSOLVERMp, cuBLASMp and NCCL are not part of the toolkit; pass their locations, e.g. from the pip wheels `nvidia-cusolvermp-cu13`, `nvidia-cublasmp-cu13`, `nvidia-nccl-cu13`:
  `cmake .. -DCUSOLVERMP_ROOT=<site-packages>/nvidia/cu13 -DCUBLASMP_ROOT=<site-packages>/nvidia/cublasmp/cu13 -DNCCL_ROOT=<nccl prefix>`
  (or an NVIDIA HPC SDK `math_libs` / `comm_libs` directory). The library directories are baked into the RPATH of `libcurcuma_cuda_mgpu.so`, so the run host must see the same paths (shared file system) or have them on `LD_LIBRARY_PATH`. `cmake` prints which backends were found.

- `CMAKE_CUDA_ARCHITECTURES` defaults to `75;80;86;89;90;100;120` with CUDA >= 13 (Turing to Blackwell, incl. H100/H200 sm_90). For fast local builds pass only your card, e.g. `-DCMAKE_CUDA_ARCHITECTURES=86`. A binary without the card's architecture falls back to JIT or fails the kernels.
