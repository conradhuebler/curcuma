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
| `-scf_mixed_precision false` | true | FP64 eigensolves only. Try on datacenter cards with strong FP64 (H100/H200) - not yet measured there. |
| `-scf_threshold 1e-8` | 1e-5 | Tight SCF. Needed for GPU gradients that must agree with the CPU beyond ~1e-4 Eh/A (the loose default lands the device Broyden at a different point inside the tolerance). |

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

- `CMAKE_CUDA_ARCHITECTURES` defaults to `75;80;86;89;90;100;120` with CUDA >= 13 (Turing to Blackwell, incl. H100/H200 sm_90). For fast local builds pass only your card, e.g. `-DCMAKE_CUDA_ARCHITECTURES=86`. A binary without the card's architecture falls back to JIT or fails the kernels.
