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

- Measured, polymer_2x GFN2 on 4x RTX A4500 (identical with the defaults and with explicit `-gpu_eigensolver_devices all -gpu_density_devices solver`, both 194 s): density 3.35 -> **0.85 s per SCF step** (46.9 -> 10.2 s total), FP32 iteration 9.4 -> 6.9 s, wall **242 -> 194 s**, energy identical. Device 0 peak 13.5 -> 13.9 GB, helpers 5.0 -> 5.8 GB.
- Too small to pay: polymer (1410 atoms, nao 2975) goes 50 -> 55 ms per step, i.e. the transfers cost more than the split saves. The split needs at least 16 columns per worker and falls back to the single device otherwise.
- On any failure it falls back to the single-device kernel for the rest of the run and warns (at verbosity >= 1).
- Rejected alternatives for the same kernel: a row-major variant (contiguous rows, occupations applied in the kernel) was **8x slower** (4.69 vs 0.60 s over 12 steps on polymer); a banded GEMM does not help because the atom order of polymer_2x is not spatially local - a column block spans 93-99 % of the rows although only 7.9 % of the atom pairs are inside the cutoff.

- Not measured: `mg` end to end on polymer_2x in FP64 (standalone FP64 26.4 s vs `mp` 20.7 s on 4 GPUs), GFN1 on polymer_2x, NVLink systems (H200).

Environment variables:

| Variable | Effect |
|---|---|
| `CURCUMA_GPU_PROFILE=1` | Stream-synchronised per-phase device timings (integrals, potential, Fock, reduce/syevd/back-transform per precision, density, charges, Broyden) plus the post-SCF host phases. Adds synchronisation cost; for diagnosis only. |
| `CURCUMA_GPU_MP_OTF=0\|1` | Force the GFN2 multipole interaction matrices to be stored (0) or rebuilt per iteration on the device (1). Default: rebuilt above ~1 GB of matrices (~2700 atoms). On-the-fly costs ~50 ms/iteration at 1410 atoms. |

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

## 4. Build

- Multi-GPU eigensolver (`CURCUMA_MULTI_GPU_EIGENSOLVER`, default ON): cusolverMg is found in the CUDA toolkit. cuSOLVERMp, cuBLASMp and NCCL are not part of the toolkit; pass their locations, e.g. from the pip wheels `nvidia-cusolvermp-cu13`, `nvidia-cublasmp-cu13`, `nvidia-nccl-cu13`:
  `cmake .. -DCUSOLVERMP_ROOT=<site-packages>/nvidia/cu13 -DCUBLASMP_ROOT=<site-packages>/nvidia/cublasmp/cu13 -DNCCL_ROOT=<nccl prefix>`
  (or an NVIDIA HPC SDK `math_libs` / `comm_libs` directory). The library directories are baked into the RPATH of `libcurcuma_cuda_mgpu.so`, so the run host must see the same paths (shared file system) or have them on `LD_LIBRARY_PATH`. `cmake` prints which backends were found.

- `CMAKE_CUDA_ARCHITECTURES` defaults to `75;80;86;89;90;100;120` with CUDA >= 13 (Turing to Blackwell, incl. H100/H200 sm_90). For fast local builds pass only your card, e.g. `-DCMAKE_CUDA_ARCHITECTURES=86`. A binary without the card's architecture falls back to JIT or fails the kernels.
