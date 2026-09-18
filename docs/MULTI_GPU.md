# Multi-GPU (work in progress, branch `feature/multi-gpu`)

> 🤖 AI-generated, not human-tested. Numbers below are measurements, not guarantees.
>
> All user-facing options and environment variables with defaults and when to change them: [GPU_TUNING.md](GPU_TUNING.md).

Plan: `~/.claude/plans/neuer-branch-wir-wollen-happy-marshmallow.md` (phases 0-6).
Two cases: (A) batch of structures distributed over GPUs, (B) one large system on several GPUs.

## Phase 0 baseline (Sep 17, 2026, master c4cde7f2)

Hardware: 4x RTX A4500 (sm_86, 20 GB, PCIe, P2P ok), 36 cores, 125 GB RAM, CUDA 13.3.
Build: `release/`, `CMAKE_CUDA_ARCHITECTURES=75`. Logs: `scripts/gpu_bench/baseline_logs/`,
runner `scripts/gpu_bench/bench.sh` (wall, peak RSS, peak GPU memory via `tm.py`).
Caveat: the polymer (1410 atoms) runs had ollama occupying 9-12 GB per GPU; the polymer_2x runs had free GPUs.

### Single point

| System | Method | CPU 36 thr | 1x A4500 | Energy CPU vs GPU |
|---|---|---|---|---|
| polymer (1410) | gfn2 | 41.8 s (SCF 12 it 39.2 s) | 12.6 s (SCF 11 it 7.6 s), GPU +2.5 GB | identical to 1e-8 |
| polymer (1410) | gfn1 | 50.6 s (SCF 48.7 s) | 12.9 s (SCF 10.8 s) | identical to 1e-8 |
| polymer (1410) | gfnff | 1.1 s | 2.8 s | 3e-7 Eh |
| polymer_2x (7320) | gfnff | 9.6 s wall, energy call 163 ms | 13.2 s (-threads 36), 20.5 s (-threads 1), energy call 137 ms | 5.7e-7 Eh |

- GFN-FF single point is >98 % CPU setup (topology + parameters); the GPU energy call is 137 ms. Gate G0e: passed -> multi-GPU cannot speed up a GFN-FF SP.
- GPU gfn2/gfn1 wall minus reported TOTAL: ~0.7-1 s untimed (context/plugin init, uploads).
- GFN-FF CPU vs GPU energy differs by 5.7e-7 Eh on polymer_2x (not investigated; GPU EEQ cutoff suspected).

### GFN-FF MD, polymer_2x, 50 steps, dt 1 fs, -threads 36 (GPU run used the wrong gradient unit, see below)

| | Wall | Per step (approx.) | Peak |
|---|---|---|---|
| CPU | 355 s | ~6.9 s | 9.7 GB RSS |
| 1x A4500 | 83 s | ~1.4 s | 10.9 GB RSS, 2.2 GB GPU |

- CPU MD heats up (T 298 -> 20476 K at 40 fs) while the GPU run stays near 300 K. polymer_2x is **not pre-optimized** (operator note), so large initial forces make divergent trajectories plausible; not treated as a bug. Use an optimized structure for MD timings.
- ~1.4 s/step on GPU for 7320 atoms is far above the 23 ms/step measured at 1410 atoms -> per-step cost does not come from the force kernels alone (suspects: CN pair-list rebuild each step, dense EEQ potrf, syncs; see plan Phase 6).

### Dense eigensolver, n = 15444 (GFN2 nao of polymer_2x), random symmetric matrix

Tool: `test_cases/cuda/bench_syevd_mg.cpp` (standalone, build line in its header). Solve time only; data already on the device(s). Eigenvalues Dn vs Mg agree to 1e-13.

| Solver | GPUs | FP64 | FP32 | memory per GPU (FP64) |
|---|---|---|---|---|
| cusolverDnXsyevd | 1 | 50.2 s | **7.6 s** | 9.1 GB |
| cusolverMgSyevd | 1 | 61.4 s | 23.7 s | 8.0 GB |
| cusolverMgSyevd | 2 (PIX pair 0,1) | 37.5 s | 15.8 s | 4.9 GB |
| cusolverMgSyevd | 2 (NODE pair 0,2) | 37.5 s | 15.8 s | 4.9 GB |
| cusolverMgSyevd | 4 | 26.4 s | 12.6 s | 4.0 GB |

- Gate G0d (Mg 4 GPU >= 1.3x Dn 1 GPU): **FP64 passes (1.9x)**, **FP32 fails** - single-GPU Dn FP32 (7.6 s) beats every Mg configuration.
- On the A4500 FP32 is 6.6x faster than FP64. Multi-GPU therefore only helps the FP64 iterations (the final polishing steps of the mixed-precision SCF) and as a memory enabler (4.0 vs 9.1 GB per device).
- PCIe topology (PIX vs NODE) makes no measurable difference here.
- `cusolverMg` is marked **deprecated** in CUDA 13.3 (successor: cuSOLVERMp on NCCL/CAL). Integrated as FP64-only fallback; the eigenvalue-sum check used here did not catch its unusable FP32 eigenvectors at this size (see step 3).
- cuSOLVERMp 0.9.1 (same matrix): FP64 1 GPU 50.2 s, 2 GPUs 36.3 s, 4 GPUs 20.7 s (nb 128); FP32 4 GPUs 5.24 s, 2 GPUs 7.36 s; 1.5-2.7 GB per GPU. Tool: `test_cases/cuda/bench_syevd_mp.cpp`.
- Cusolver Mg layout pitfall: the column blocks are dealt out **cyclically** (block b -> device b % ndev). A contiguous layout returns wrong eigenvalues without any error.

### Not yet measured
- gfn2 polymer_2x CPU reference: -11784.87804452 Eh, 11 SCF iterations at ~320-430 s each, 4468 s wall (24 threads, shared machine), 49 GB RSS. gfn1 polymer_2x not yet measured.

## Device selection (Phase 1, implemented)

- `-gpu_device N` (global CLI key): pins native gfn1/gfn2 (CUDA/ROCm/Vulkan) and GFN-FF (CUDA/ROCm) to device N. Indices are the runtime's, i.e. after `CUDA_VISIBLE_DEVICES` / `ROCR_VISIBLE_DEVICES` (SLURM sets these). Default -1 = unchanged behaviour (device 0 / current device).
- Invalid index: error message and CPU fallback.
- `curcuma -methods` lists the devices each plugin sees.
- Verified on 4x A4500: gfn2 `complex` on device 3 and GFN-FF `polymer` on device 2 (nvidia-smi PCI bus), energies bit-identical to device 0.
- The xTB `__constant__` tables are now uploaded per device (was once per process -> uninitialised tables on any second device).
- ROCm and Vulkan changes are compile-unverified (no SDK on the dev box).

## H200 report (operator, Sep 17, 2026)
- GFN2 on the ~7k-atom molecule: **110 s on one H200** (so it did run on the GPU; the missing sm_90 architecture was NOT the problem there).
- 1 -> 2 GPUs gave no speedup: expected, nothing used a second device (no device selection existed).

## Large-system GPU memory (case B, step 1a implemented)
- `XtbGpuContext::estimateResidentBytes()` sums every resident buffer with cuSOLVER-queried workspace sizes; `beginBasis()` refuses a basis that does not fit and frees partial allocations (previously: silent partial allocation, leak, CPU fallback visible only at verbosity 2).
- The refusal is a warning at default verbosity, e.g. polymer_2x GFN2 on an A4500: `needs about 49.9 GB ... only 19.4 of 19.6 GB are free on device 2; this calculation runs on the CPU`. `-gpu_memory_check false` disables the check.
- Estimate vs nvidia-smi peak: complex 0.31 GB vs 0.49 GB, polymer 2.36 GB vs 2.51 GB (peak includes ~0.2-0.4 GB CUDA context).
- Which matrices are really needed per SCF iteration: dense only C (holds F), L, the syevd workspace and a transient P; H0, S and the 9 multipole integrals can share one screened AO-pair list (16 % of pairs within 40 Bohr on polymer_2x). The host currently also builds dense copies of the multipole integrals (17 GB) and S/H0/L/gamma on the CUDA path. Next steps.

## Large-system GPU memory, steps 1b/1c (implemented, Sep 17, 2026)

Result: **GFN2 on polymer_2x (7320 atoms, nao 15444) runs on ONE RTX A4500 (20 GB)**: device peak 17.0 GB (was ~50 GB), 619 s total, E = -11784.87804452 Eh, **identical to the CPU reference** (-11784.87804452 Eh, 24 threads, 4468 s wall while other jobs ran). Host peak RSS 39 GB (CPU run: 49 GB).

What changed (all CUDA, `xtb_gpu_context.cu`):
- **Screened pair storage** for S, H0 and the 9 GFN2 multipole integrals (`-gpu_sparse_integrals auto|on|off`, auto = when < 50 % of AO pairs survive). Atom-pair cutoff from each element's smallest primitive exponent and largest coefficient, integrals < 1e-20 dropped (~30-36 Bohr). Kernels `k_multipole_ints_sp`, `k_build_fock_sp`, `k_pop_ao_sp`, `k_multipole_moments_sp`, `k_grad_h0_pulay_sp` perform the same arithmetic as the dense ones.
- **On-the-fly multipole interaction**: the 18 nat^2 matrices (7.7 GB at 7320 atoms) are rebuilt per pair in `k_multipole_potential_otf` / `k_energy_multipole_otf` above 1 GB of matrices (`CURCUMA_GPU_MP_OTF=0/1` forces it). Also removes the host upload copy.
- **Exclusive eigensolver workspaces**: FP32 and FP64 syevd buffers never coexist; the gradient frees both before allocating W. Cw is n x ncol.
- Memory estimate follows the peak (max of the two workspaces, screened sizes).

Validation (A4500, `-gradient`):

| system | method | dense vs screened E | max abs gradient diff | screened fraction |
|---|---|---|---|---|
| complex | gfn2 / gfn1 | identical (8 dp) | 1.0e-15 / 1.1e-15 | 93.8 % / 98.4 % (forced) |
| polymer | gfn2 / gfn1 | identical (8 dp) | 7.4e-16 / 7.6e-16 | 41.4 % / 54.5 % |
| complex / polymer | gfn2 OTF on vs off | identical (8 dp) | 6.0e-15 / 9.4e-15 | - |

GPU peak polymer GFN2 2509 -> 1813 MiB. OTF costs ~50 ms per iteration at 1410 atoms (hence only for large systems). 205/205 GPU + sqm validation ctests pass.

polymer_2x GFN2 profile on one A4500 (11 iterations, 4 in FP64):

| phase | time | share |
|---|---|---|
| eig FP64 syevd (4) | 200 s | 36 % |
| density P + populations (11) | 143 s | 26 % |
| eig FP64 reduce (4) | 86 s | 16 % |
| eig FP32 syevd (7) | 53 s | 10 % |
| eig FP64 back-transform (4) | 43 s | 8 % |
| host setup / post-SCF | 47 s / 30 s | - |

### Pattern-only density, deferred host integrals, GPU FP32 threshold 1e-5 (Sep 17, 2026)
- **`k_density_sp`**: with screened storage the SCF loop computes P only at stored pairs (populations, moments and band energy read nothing else); dense P is rebuilt on demand for the gradient and the host download. complex/polymer: energies identical to 12 digits, gradients <= 1e-15. polymer density+populations 1357 -> 744 ms (41 % pattern).
- **Deferred host multipole integrals**: for nao above ~5000 (> 2 GB of integrals) on the CUDA resident path the host no longer builds the 9 dense nao^2 dipole/quadrupole matrices, nor the post-SCF host Fock matrix (only read by debug dumps; `getFock()` is empty then). Host fallbacks (host SCF loop, host gradient) build them on demand; deferral is off for `d4_charge_source=cpscf`.
- **GPU default `scf_fp32_threshold` = 1e-5** (operator decision; CPU stays 1e-3; an explicit value wins).

polymer_2x GFN2, one A4500, before -> after these three:

| | before | after |
|---|---|---|
| total wall | 619 s | **410 s** |
| host peak RSS | 39 GB | **22 GB** |
| device peak | 17.0 GB | 17.1 GB |
| host setup | 47 s | 21 s |
| SCF iterations (FP64) | 11 (4) | 16 (2) |
| density P + populations | 143 s | 53 s |
| eig FP32 syevd | 53 s (7) | 106 s (14) |
| eig FP64 syevd / reduce / back-transform | 200 / 86 / 43 s | 100 / 43 / 21 s |
| post-SCF host energies | 30 s | 40 s |

Energy identical to the CPU reference (-11784.87804452 Eh). Remaining host hot spots: post-SCF energies (40 s) and setup (21 s). ctest: 207/208 GPU/SQM/CPSCF/gradient tests pass; `xtb_cpscf` fails identically without these changes (pre-existing, CLAUDE.md known failures).

### Post-SCF host phase (Sep 17, 2026)
Breakdown printed with `CURCUMA_GPU_PROFILE=1` or `-verbosity 3`. Changes:
- **D4 on the device (CUDA)**: 2-body energy/gradient/dEdcn/dEdq (`k_d4_grad`, port of the ROCm kernel) and ATM 3-body (`k_d4_atm_nl`, ROCm kernel + host neighbour list within the 25 Bohr cutoff; for alp = 16 the damping power uses q^5 cbrt(q) instead of pow). polymer: D4 2-body 173 -> 7 ms, ATM 2130 ms (host, 8 threads) -> 1840 ms (device; the A4500 is FP64-weak). Energies and gradients unchanged (max 3e-16).
- **Large single point keeps the wavefunction on the device**: with deferred host multipole integrals and no gradient requested, P and C are not rebuilt/downloaded (16 s), the host potential is not rebuilt (5 s) and the Coulomb/third-order/multipole/band energies are taken from the last device-resident step (5 s of O(nat^2) host work). `ensureHostWavefunction()` downloads them on demand (host gradient fallback). `MolecularOrbitals()` is empty in that case.

polymer_2x GFN2, one A4500:

| | 3 steps ago | now |
|---|---|---|
| post-SCF | 40 s | **4.9 s** |
| total wall | 410 s | **375 s** |
| host peak RSS | 22 GB | **15.6 GB** |
| device peak | 17.1 GB | 16.1 GB |

All energy components identical to 8 dp (Electronic, Coulomb, third-order, multipole, repulsion, dispersion; total -11784.87804452 Eh). 208/208 GPU/SQM/D4/gradient ctests pass. Remaining post-SCF: D4 4.1 s, repulsion 0.8 s.

### FP32 threshold (decided: 1e-5 on the GPU)
polymer GFN2, one A4500, `-scf_fp32_threshold`:

| threshold | iterations | FP64 syevd calls | SCF | energy | max abs gradient diff vs CPU |
|---|---|---|---|---|---|
| 1e-3 (default) | 11 | 4 | 6.71 s | -2088.25340678 | 2.6e-6 |
| 1e-4 | 11 | 2 | 5.05 s | same | - |
| 1e-5 | 12 | 1 | 4.50 s | same (1e-12) | 7.0e-6 |
| 3e-6 | 15 | 1 | 5.41 s | same | - |

complex: 3.5e-4 Eh/A GPU-vs-CPU gradient difference at both thresholds (known: the loose default `scf_threshold` lands the device Broyden at a different point; use `-scf_threshold 1e-8` for gradients).

## Step 2: GFN-FF setup (Sep 17, 2026, in progress)

`CURCUMA_GFNFF_PROFILE=1` prints every setup phase summed over both q-loop passes (the verbosity-2 report shows only the last pass) and the GPU upload phases.

Findings on polymer_2x (7320 atoms, 1502 fragments):
- **OpenMP was pinned to one thread**: `CxxThreadPool` calls `omp_set_num_threads(1)` process-wide, so every `omp parallel for` in the topology (distance matrix, CN, Dijkstra) ran serially even with `-threads 36`. Topology and parameter generation now open the GFN-FF thread budget with `ScopedBlasThreads` (inside a batch worker the batch's intra budget). All affected loops write row/atom-local results; energy identical to 12 digits.
- **nb_hc / nb_nometal** neighbour lists (2 x ~1 s per pass) parallelised over rows: 2.2 s -> 0.29 s. GFN-FF CN 0.33 -> 0.02 s.
- **Implicit Coulomb pairs on the GPU** (`-gpu_coulomb_implicit`, default true): the device gathers over all atom pairs (`k_coulomb_implicit`, gamma_ij = 1/sqrt(alp_i + alp_j) from per-atom alpeeq) instead of reading the host-built N^2/2 list (24.5 M pairs, 2.7 GB, kept twice on the host). Not used with `eeq_distance_cutoff > 0`.

| polymer_2x GFN-FF single point | before | now |
|---|---|---|
| CPU (-threads 36), wall | 9.5 s | 6.9 s |
| GPU (1x A4500), wall | 13.2 s | **7.4 s** |
| GPU, host peak RSS | 11.0 GB | **4.8 GB** |

Energies: CPU -898.805534472558 identical; GPU pair list vs implicit -898.805533905291 / ...292, gradients <= 2e-15 (complex, polymer, polymer_2x). ctest gfnff/GPU/MD/opt: 215/216 (only the known `cli_curcumaopt_07` golden drift).

Second round (same day), all row-parallel with per-row buffers appended in row order, results identical:
- EEQ phase 1 Coulomb matrix fill (26.8 M erf): 1.61 -> 0.99 s (both passes)
- bond list: 342 -> 43 ms; topo distances (BFS) + nbondmat + BATM list: 522 -> 292 ms; the verbosity-3-only 1,4-pair debug count no longer runs at lower verbosity
- GPU workspace upload 1427 -> 402 ms (no Coulomb pair list)

| polymer_2x GFN-FF single point, 1x A4500 | start of step 2 | now |
|---|---|---|
| wall | 13.2 s | **6.8 s** |
| with `-gpu_disp_pairs_on_device true` (existing option, default off) | - | **5.6 s** (energy identical, gradient 5e-16) |
| host peak RSS | 11.0 GB | 4.8 GB (4.1 GB with device dispersion pairs) |

Remaining: EEQ solves ~2.6 s (projected PCG, 32 + 71 iterations per pass, matvec already threaded BLAS and memory-bandwidth bound), distance matrix 0.16 s, torsions 0.19 s.

## GPU phase profiler
`CURCUMA_GPU_PROFILE=1 curcuma -sp mol.xyz -method gfn2 -gpu cuda -verbosity 1` prints stream-synchronised per-phase device timings (integrals, potential, Fock, reduce / syevd / back-transform for FP32 and FP64, density, charges, energy, Broyden). Off by default (no synchronisation cost). `-verbosity 2` additionally shows a `pre-SCF` line (device uploads, EEQ guess) that was previously only inside TOTAL.

polymer (1410 atoms, nao 3222), GFN2, one A4500, before -> after the occupied-column density:

| phase | before | after |
|---|---|---|
| eig FP64 syevd (4 calls) | 2665 ms | same |
| density P + populations (11 calls) | 2138 ms | **1227 ms** |
| eig FP32 syevd (7 calls) | 1026 ms | same |
| eig FP64 reduce (4 calls) | 891 ms | same |
| SCF total | 7670 ms | **6765 ms** |

Energies unchanged at the printed precision (complex -329.52714784, polymer -2088.25340678). The density GEMM now uses the occupied columns only (last occ > 1e-12), as the CPU path always did. Observation: 4 of 11 iterations run in FP64 and cost 43 % of the SCF on this consumer card.

## Case A: batch of structures over several GPUs (implemented)

- **Pool**: `src/core/gpu_device_pool.{h,cpp}`. `main()` configures it from `-gpu`, `-gpu_devices` ("0,2", default all visible), `-gpu_workers_per_device` (default 1). Batch workers take a `GpuDeviceLease` (blocks until a slot is free, least-loaded device); `EnergyCalculator::createMethod` passes the leased device as `gpu_device`.
- **Wired into**: `-sp` multi-structure batch (new), `-opt` multi-structure batch, `CurcumaOpt` SP/Opt threads, ConfSearch MD + optimisation workers, numerical Hessian workers.
- **Default kept**: a single visible GPU without `-gpu_devices`/`-gpu_workers_per_device` does not activate the pool (workers share device 0 as before).
- **ConfSearch**: no longer switches the GPU off for `-threads > 1` when a pool is active (surplus workers wait for a slot).
- **Host threads**: a GPU batch worker gets `cores / GPU slots` intra-molecule threads (`intraThreadBudget`, `intra_parallel_context.h`); CPU batch workers stay serial as before.
- **`-sp` on a multi-XYZ file** (new): all frames, energies in input order, `<basename>.sp.xyz` with energies in the comment line (BMT-aware). Workers = max(`-threads`, GPU slots).

### Measurements (4x RTX A4500, GFN2, energies bit-identical in every configuration)

8 perturbed copies of `polymer` (1410 atoms). A CPU GFN2 run on 24 cores was running in parallel during these measurements.

| GPU slots | before host-thread budget | with budget |
|---|---|---|
| 1 (`-gpu_devices 0`) | 175 s | 100 s |
| 2 | 82 s | - |
| 4 | 42 s | **28.5 s** |
| 8 (`-gpu_workers_per_device 2`) | 31 s | 25 s |

40 conformers of AAA-bGal (90 atoms): 1 GPU 5.2 s, 4 GPUs 2.0 s, CPU 36 workers **0.58 s** -> for molecules of this size the CPU batch is faster; GPU batch pays off from a few hundred atoms per structure.

## Bug found on the way: GPU GFN-FF gradient units
- `GFNFFGpuMethodImpl::getGradient()/copyGradientTo()` returned Eh/Bohr; the `ComputationalMethod` contract is Eh/Angstrom (Known Issue #28 fixed only the CPU wrapper). GPU GFN-FF MD forces were 1.89x too small, 18/19 `gfnff_gpu*` ctests failed on master with a gradient norm ratio of 0.529. Fixed (CUDA + ROCm share the template); 19/19 pass.
- This also explains why the GPU MD on the unoptimized polymer_2x looked stable while the CPU MD heated up.

### Pre-existing build issue
- `test_cases/sqm_reference/test_xtb_cuda_*` do not compile (`curcuma::xtb::gpu` undeclared) since `USE_CUDA` was removed from the core definitions; main binary and `libcurcuma_cuda.so` build fine.

## Step 3: multi-GPU eigensolve (implemented, Sep 17, 2026)

- `DistributedEigensolver` (`qm_methods/cuda/xtb_distributed_eigensolver.{h,cpp}`) in its own library `libcurcuma_cuda_mgpu.so`, loaded on first use by `libcurcuma_cuda.so` (`xtb_distributed_eigensolver_loader.cpp`), so a host without cuSOLVERMp/NCCL keeps single-GPU runs.
- Backends: `mp` = cuSOLVERMp + cuBLASMp over NCCL (one host thread per GPU, communicators from `ncclCommInitAll`); `mg` = cusolverMg, FP64 only.
- Hook: `XtbGpuContext::eigensolveResidentFock`. FP64: scatter F (and L once per geometry), `sygst` + `syevd` + `trsm` on the GPUs, gather C. FP32: `syevd` only. Falls back to the single-GPU cuSOLVER path on any failure, with a warning.
- Pitfalls found: cuBLASMp grid destruction is collective (a sequential teardown deadlocks at exit); stdout redirected to a file is block-buffered, so a long run looks hung (use `stdbuf -oL`).
- Also fixed on the way: the gradient allocated W before releasing the eigensolver workspaces (raised the device peak); workspaces are now released when the SCF ends.
- Validation: complex (gfn1, gfn2) at `-scf_threshold 1e-9` vs single GPU, energy and gradient <= 1.2e-9 for mp (4 GPUs, 2 GPUs not including device 0) and mg (FP64); polymer_2x energy identical at the printed 8 decimals; 200/200 `ctest -L gpu`.
- Numbers and options: [GPU_TUNING.md](GPU_TUNING.md#multi-gpu-eigensolve-one-large-molecule-on-several-gpus).

## Step 3b: distributed pattern density (implemented, Sep 17, 2026)

- `-gpu_density_devices all|list|solver`: the screened-pattern density is a sum over the occupied columns, so every device evaluates the full pattern over its own column slice and the partials are added with a daxpy on the calculation's device (exact).
- polymer_2x GFN2: density 46.9 -> 10.2 s, wall 242 -> **194 s** (single GPU 419 s), energy identical; complex at `-scf_threshold 1e-9` agrees with the single-GPU run to 1.2e-9; `ctest -L gpu` 200/200.
- Not worth it below ~3000 basis functions (polymer, nao 2975: 50 -> 55 ms per step).
- Trap found here: `cudaDeviceEnablePeerAccess` leaves `cudaErrorPeerAccessAlreadyEnabled` pending, and the next `cudaGetLastError()` reported it as a kernel-launch failure - the peer calls now consume it.

## Defaults (operator decision, Sep 17, 2026)

Both multi-GPU paths for one large molecule are ON when more than one device is visible and the
calculation does not run inside a batch worker (`leasedGpuDevice() < 0`), gated at 4000 basis
functions (`gpu_eigensolver_min_nao`, `gpu_density_min_nao`). `-gpu_eigensolver_devices none` /
`-gpu_density_devices none` switch them off. polymer_2x with no flags at all: 194.2 s, device-0
peak 13.9 GB, energy identical to the single-GPU run - the same numbers as with the options set
explicitly. Below the gate the status line says so ("not used (nao below ...)"), and
`ctest -L gpu` (200 tests, all small molecules) is unaffected: 200/200.

## What to measure on the H200 node (Sep 18, 2026)

Ordered by what curcuma cannot answer here. Everything below runs from the repository; nothing
needs a code change. Bring back the `.json` files and the profile output.

**0. Check the build first - without this the rest measures the wrong thing.**
```bash
bash scripts/find_mgpu_libs.sh        # what this node has (try `module load nvhpc` first)
cmake .. -DCMAKE_CUDA_ARCHITECTURES=90 -DCUSOLVERMP_ROOT=... -DCUBLASMP_ROOT=... -DNCCL_ROOT=... \
      -DCURCUMA_REQUIRE_MULTI_GPU_EIGENSOLVER=ON        # fails the configure if they are missing
cmake .. ... 2>&1 | grep "=== curcuma multi-GPU eigensolver"   # or check the summary by hand
curcuma -methods            # "multi-GPU eigensolver: mp (cuSOLVERMp)", plus the H200s
```
The Sep 17 H200 run fell back to cusolverMg because cuSOLVERMp/NCCL were missing, and Mg is 15x
slower than the single-GPU path on our own measurement - a build without them looks like
"multi-GPU does not help" when nothing multi-GPU ran.

**1. The sweep, for the record.** `python scripts/tuning_sweep.py
test_cases/molecules/larger/polymer_2x.xyz --method gfn2 --gpu cuda --repeats 2 --json h200.json`
(~2 h). It checks every energy against the baseline, so a divergence shows up as SUSPECT rather
than as a number to be believed.

**2. Is mixed precision really the wrong default there?** The current default is FP64-only on
full-rate-FP64 cards, decided from one operator log (5.29 s FP32 vs 3.36 s FP64 per iteration,
and an FP32 phase converging 1.1 kcal/mol off). Measure `-scf_mixed_precision true|false` x
`-scf_fp32_threshold 1e-5|1e-6` and watch the ITERATION COUNT as much as the wall time: if
`true` now converges in the same number of iterations as `false`, the false-fixed-point guard
has removed the reason for the default and it should be reconsidered.

**3. The per-phase profile - the one measurement that is still missing everywhere.**
`CURCUMA_GPU_PROFILE=1 curcuma -sp polymer_2x.xyz -method gfn2 -gpu cuda -verbosity 3`. It says
how much of the run is eigensolve, density, Fock, integrals. That decides whether distributing
the integrals and the Fock build (column-block ownership) is worth writing at all: if the
eigensolve plus density is already 80 % of the run on that hardware, the answer is no.

**4. Does the split scale past 2 GPUs on NVLink?** `-gpu_eigensolver_devices 0 | 0,1 | 0,1,2,3`
(and the same for `-gpu_density_devices`). On PCIe we measured 1.53x and 1.34x at 4 devices with
7320 atoms. NVLink should do better per device - or the far faster single card may leave nothing
to win, which is equally worth knowing.

**5. What the 141 GB card allows that ours does not.** `-gpu_multipole_otf off` and
`-gpu_sparse_integrals off` both fail or time out on a 20 GB A4500. If they run there, compare
them against the defaults: storing the multipole matrices may well beat rebuilding them when
memory is free, and that would justify making the `auto` threshold device-memory-aware instead
of a fixed ~2700 atoms.

**6. How large can a system get now?** polymer_2x needs about 14 GB on one device. With 141 GB
the interesting question is where the next wall is (nao^2 buffers, the int32 indexing noted in
the size-guard work), so: one run on the largest structure you have.

**7. Batch throughput**, if the node has 4 or 8 cards: `-sp` on a multi-XYZ file with
`-gpu_devices 0,1,2,3 -gpu_workers_per_device 1|2`. Independent of everything above and the
easiest real-world win.

## Operator runs on other hardware (Sep 17/18, 2026)

- **2x H200 NVL, build without cuSOLVERMp/NCCL**: the eigensolve fell back to cusolverMg, which
  curcuma then only used for FP64 - one distributed solve out of 21, so the run looked
  single-GPU. polymer_2x gfn2: 15 iterations / 52 s before the distributed density, 21 / 66 s
  with it; the extra iterations are the FP32 noise-floor effect that the stagnation guard now
  addresses (`scf_fp32_threshold` == `scf_threshold` == 1e-5), not a defect of the density path.
  What that machine needs is a build WITH cuSOLVERMp + cuBLASMp + NCCL.
- **2x RTX PRO 5000 Blackwell, Mg-only build, current code**: the per-solve verification rejected
  Mg's FP32 eigenpairs (relative residual 3.5e-3) and kept FP64 distributed. Both runs below
  converged in 12 iterations to -11784.87804452 Eh, the same energy as our 4x A4500 runs:

  | polymer_2x gfn2, same machine | one GPU | `-gpu_density_devices all` |
  |---|---:|---:|
  | FP32 iteration | 4.62 s | **4.00 s** |
  | FP64 iteration | 32.5 s | **25.3 s** |
  | total | 100 s | **91 s** |

  This is the first multi-GPU measurement on hardware other than the A4500 box, and it is worth
  reading carefully: the 9 % came WITHOUT a distributed FP32 eigensolve (Mg was rejected), i.e.
  from the distributed density plus the one distributed FP64 solve. On 2 GPUs the FP32 eigensolve
  is worth little anyway (A4500: 12.2 -> 11.7 s per iteration on two devices), so the way to more
  on that machine is more devices or a build with cuSOLVERMp, not the current backend.
- **Full knob sweep on the A4500 box** (Sep 18, 2026): all 15 knobs on polymer_2x, 45 runs, one
  energy for all of them; the numbers and what they mean are in
  [GPU_TUNING.md](GPU_TUNING.md#a-full-sweep-measured). The same sweep is the thing to run on the
  H200 node: `python scripts/tuning_sweep.py test_cases/molecules/larger/polymer_2x.xyz --method
  gfn2 --gpu cuda --repeats 2 --json h200.json`.
- Nothing here is a curcuma measurement on NVLink hardware: the per-phase profile
  (`CURCUMA_GPU_PROFILE=1`) has not been taken on either machine, so how much of those runs is
  distributable at all is still unknown.
