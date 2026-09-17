# Multi-GPU (work in progress, branch `feature/multi-gpu`)

> 🤖 AI-generated, not human-tested. Numbers below are measurements, not guarantees.

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
- `cusolverMg` is marked **deprecated** in CUDA 13.3 (successor: cuSOLVERMp on NCCL/CAL). Mg integration only as a fallback; NCCL/cuSOLVERMp path to be evaluated once NCCL is installed.
- Cusolver Mg layout pitfall: the column blocks are dealt out **cyclically** (block b -> device b % ndev). A contiguous layout returns wrong eigenvalues without any error.

### Not yet measured
- gfn2/gfn1 polymer_2x on CPU (~45 GB RAM, estimated > 1 h; running) and on GPU (needs ~45 GB -> silent CPU fallback expected on 20 GB).

## Device selection (Phase 1, implemented)

- `-gpu_device N` (global CLI key): pins native gfn1/gfn2 (CUDA/ROCm/Vulkan) and GFN-FF (CUDA/ROCm) to device N. Indices are the runtime's, i.e. after `CUDA_VISIBLE_DEVICES` / `ROCR_VISIBLE_DEVICES` (SLURM sets these). Default -1 = unchanged behaviour (device 0 / current device).
- Invalid index: error message and CPU fallback.
- `curcuma -methods` lists the devices each plugin sees.
- Verified on 4x A4500: gfn2 `complex` on device 3 and GFN-FF `polymer` on device 2 (nvidia-smi PCI bus), energies bit-identical to device 0.
- The xTB `__constant__` tables are now uploaded per device (was once per process -> uninitialised tables on any second device).
- ROCm and Vulkan changes are compile-unverified (no SDK on the dev box).

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
