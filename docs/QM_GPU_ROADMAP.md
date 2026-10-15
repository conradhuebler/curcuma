# Native ab-initio QM (HF / HF-3c / KS-DFT): performance and GPU roadmap

> Status: 🤖 AI-authored. Sections 1-2 are **implemented and machine-tested** (CPU only,
> Sep 2026). Section 3 is a plan; stages Q-G1, Q-G2 and Q-G3 are implemented (Oct 2026,
> AI-generated, machine-tested on one RTX 5090 only), Q-G4 onward are proposals.

The native xTB GPU work ([SQM_GPU.md](SQM_GPU.md), [SQM_GPU_ROADMAP.md](SQM_GPU_ROADMAP.md),
[SQM_PERFORMANCE.md](SQM_PERFORMANCE.md), [MULTI_GPU.md](MULTI_GPU.md)) produced a set of
lessons that transfer directly to the Gaussian-basis HF engine (`QMEngine`,
`src/core/energy_calculators/qm_methods/qm_*.cpp`). The profile is different, though:
in xTB the eigensolve dominates; in HF the **4-centre ERIs and the J/K (Fock) builds**
dominate, and the eigensolve is small.

## 1. Done (Sep 2026): CPU restructuring in the shape a GPU needs

| xTB lesson | applied to the QM engine |
|---|---|
| Shell-pair-blocked integral kernels (per-pair data once, not per AO component) | `buildERI` groups the AO list into shells, precomputes per **shell pair** the primitive-pair table (product exponent/centre + the three Hermite tables up to the shell maxima), and per **primitive quartet** evaluates Boys + R once for all component quartets. Was: everything recomputed per AO component quartet (6^4 = 1296x for a d quartet). |
| Screened integrals (7k-atom GFN2 on one card) | Schwarz screening at shell-pair level, `-qm.eri_screening` (default 1e-12; 0 = exact). |
| On-the-fly instead of stored intermediates (`-gpu_multipole_otf`) | cartesian -> spherical 6d->5d transform applied **per shell quartet**; the cartesian n^4 tensor is never formed (was: full cartesian tensor + an O(n^8) direct 4-index sum). |
| Thread the hot loops, gate tiny systems | OpenMP over bra shell pairs (ERI) and over rows (J, K); `-qm.threads N`; J/K stay serial below 32 functions. |
| Flat data + BLAS-shaped contractions | J = ERI(n^2 x n^2) * vec(P) (GEMV per row); K as one GEMV per (mu, lam) block -- both map 1:1 onto cuBLAS. |

Measured on this container (4 cores, no BLAS), `bench_qm_integrals` (old kernels
re-implemented inside the benchmark, results compared element-wise):

| system (def2-SVP) | ERI old | ERI blocked 1 thread | 4 threads | max diff |
|---|---:|---:|---:|---:|
| H2O (25 cart.) | 460 ms | 39 ms (11.7x) | 16 ms (28x) | 2.2e-16 |
| formaldehyde (40) | 2837 ms | 178 ms (15.9x) | 62 ms (45x) | 2.2e-16 |
| water dimer (50) | 7039 ms | 419 ms (16.8x) | 145 ms (49x) | 2.2e-16 |

Spherical transform, formaldehyde: 480 ms (old O(n^8)) -> folded into the ERI build
(ERI + transform together 47 ms on 4 threads). Full run, benzene HF/def2-SVP (114
functions, 15 SCF iterations, 4 threads): 14.3 s with a separate transform pass ->
**5.9 s** on the fly (ERI 3.8 s, SCF 1.9 s of which J/K 1.6 s); energy -230.53597729 Eh
= PySCF to 8 decimals. The pre-Sep-2026 engine was not timed on benzene; extrapolated
from the kernel ratios it would take minutes (ERI) to tens of minutes (transform).

Validation of the new kernels: `ctest -L qm_2e` (full ERI tensor vs an independent
Python McMurchie-Davidson witness, 1e-10) and `ctest -L qm_1e`, `-L qm_hf3c` all pass.

**Analytic gradient (WP8, Sep 2026)**: same shell-blocked machinery, derivative tables built
one angular-momentum step higher on the differentiated centre, two-step McMurchie-Davidson
contraction (ket Hermite sums folded into R once per ket component), R^0 built on the simplex
t+u+v <= L instead of the (L+1)^4 box, no per-quartet heap allocation. Benzene HF-3c: 2e gradient
8.5 s (1 thread) / 2.2 s (4 threads), ~7x the ERI build because it runs every ORDERED bra pair
against the canonical ket pairs. Next steps there: canonical quartets with translational
invariance (derivatives on A, B, C; D from their sum) -- **done Sep 2026**: 8.5 -> 4.5 s
(1 thread), 2.2 -> 1.24 s (4 threads), now ~3.6x the ERI build -- and density-weighted
screening (open).

**What was the next wall**: the stored ERI tensor (n^4 doubles: 1.35 GB at 114 functions,
12.8 GB at 200). Since Sep 2026 the SCF switches to integral-direct J/K above
`-qm.eri_max_memory_mb` (section 2), so the limit is now time, not memory.

## 2. Integral-direct SCF (done Sep 2026, CPU)

Recompute screened shell quartets every iteration and digest them straight into J and K
(no stored tensor), with density-weighted screening `Q_ab Q_cd max|P|` and incremental
Fock builds on dP (Almlöf, Faegri, Korsell, J. Comput. Chem. 3, 385 (1982); Häser,
Ahlrichs, J. Comput. Chem. 10, 104 (1989)). Memory becomes O(n^2); cost per iteration
grows but shrinks with screening as dP -> 0. This is also the only form that makes sense
on a GPU (device memory is the binding limit there, as the 7k-atom GFN2 work showed).

**Implemented** (`qmint::DirectJK`, `-qm.scf_direct on|off|auto`): incremental builds on dP
with a full rebuild every 8 iterations, density-weighted screening, per-thread J/K
accumulators. Benzene/def2-SVP peak memory 1320 -> 34 MB at 5.6x the time (28.4 vs 5.0 s);
naphthalene/def2-SVP (180 functions, stored tensor 8.4 GB) runs in 195 s and matches PySCF
to 8 decimals (print precision). `auto` switches to direct above `-qm.eri_max_memory_mb` (4000). The cost per
build is now the shell-quartet kernel, which is what stage Q-G2/Q-G3 below would move to
the device.
Since then the kernel itself is ~3x faster (tabulated Boys function, ket-primitive-first
contraction, primitive screening, larger pair on the bra side): benzene direct 28.4 -> 8.6 s,
naphthalene 195 -> 58 s; details in [NATIVE_QM_IMPLEMENTATION.md](NATIVE_QM_IMPLEMENTATION.md).

## 3. Proposed GPU stages (mirroring the xTB stages)

Each stage keeps the rules that made the xTB port trustworthy: runtime `dlopen` plugin
(`libcurcuma_cuda.so`, [GPU_PLUGIN_STARTUP.md](GPU_PLUGIN_STARTUP.md)), CPU fallback on
every hook, and a bit-for-bit (or documented-ulp) comparison against the CPU path before
anything is enabled by default.

| Stage | content | xTB precedent |
|---|---|---|
| Q-G1 | **done Oct 2026** (see below): J/K from the stored tensor on the device, host SCF unchanged | Stage 1 (single offloaded hot spot) |
| Q-G2 | **done Oct 2026** (see below): device ERI kernel, one work item per shell quartet, fed by the flat `ShellPair`/`PrimPair` tables of section 1 | Stage 3 (integrals on device) |
| Q-G3 | **done Oct 2026** (see below): integral-direct J/K on the device (section 2), density-weighted screening, FP64 accumulation | `-gpu_multipole_otf`, screened integrals |
| Q-G4 | device-resident SCF loop: Fock -> Löwdin-reduced eigensolve (reuse the xTB `cusolverDnDsyevd` path) -> density -> DIIS; host polls O(1) scalars per iteration | Stage 6 (fully resident loop) |
| Q-G5 | mixed precision far from convergence (FP32 ERI/J/K with FP64 correction on dP), never accepting convergence on an FP32 step | X-AP3 (FP32 eigensolve early in the SCF) |
| Q-G6 | multi-GPU: distribute shell-quartet batches over the device pool | [MULTI_GPU.md](MULTI_GPU.md) |

**Q-G1 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: `-gpu cuda|auto` uploads the stored
ERI tensor once per geometry (`cuda/qm_gpu_jk.cu`, plugin entry `curcuma_cuda_create_qm_jk`, seam
`qm_gpu_backend.h`, loader `gpu_plugin::createQmJK`); each Fock build moves only P in and J, K out.
Upload failure / missing plugin falls back to the CPU and is counted (`reportGpuFallback("qm_jk")`).
Benzene HF/def2-SVP (114 functions, 12 iterations, 4 CPU threads; re-measured with the Q-G3 build,
`-gpu cuda -qm.gpu_eri false -qm.scf_direct off` vs `-gpu none -qm.scf_direct off`): J/K builds
679 -> 10.8 ms (0.9 ms per build), SCF loop 757 -> 30 ms, plus a one-off 165 ms upload; whole single
point 1.44 -> 1.15 s; energy -230.5359772878 identical to 10 decimals (= PySCF). The value
"686 -> 179 ms" first recorded here was wrong (cause not reconstructed; the 10.8 ms were already seen
in the Q-G2 session). The ERI build (617 ms) is the dominant cost, which is what Q-G2 addresses.
`ctest -L qm_gpu` (CPU vs GPU, 1e-9 Eh; skipped without a device). Not tested: other GPUs,
ROCm/Vulkan (no QM entry point). The tensor must fit in device memory (n <= ~250 on a 32 GB card);
beyond that Q-G3 takes over.

**Q-G2 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: with `-gpu cuda|auto` the stored
tensor is computed on the device (`cuda/qm_gpu_eri.cu`, entry `curcuma_cuda_create_qm_eri_jk`) and handed
to the Q-G1 J/K builder without a host copy. The host flattens the same shell/primitive-pair tables,
Hermite coefficients, Schwarz factors and Boys grid the CPU `buildERI` uses (`qmint::buildEriGpuTables`,
seam struct `QMGpuEriTables`), so both screens decide identically. Kernels: one thread per shell quartet
for the classes with <= 27 cartesian components (template per (LA,LB,LC,LD), R^n kept in one in-place
simplex), one block per quartet for the larger d classes (R^n in shared memory, one component per
thread, no atomics); 6d->5d per quartet; each tensor element written once (canonical orbit
representative) then copied to its 8 permutations; classes on 8 streams. FP64 throughout.
`-qm.gpu_eri false` keeps the Q-G1 upload path. In `scf_direct auto` a device-built tensor is used even
above `-qm.eri_max_memory_mb` (it needs no host memory) as long as it fits the card. Fallbacks (counted as
`qm_eri`): f functions or higher (the GTO basis type has no f anyway), tensor larger than free device
memory, any CUDA error -> CPU ERI + Q-G1 upload, or CPU direct SCF when the host may not hold the tensor.

Measured (HF/def2-SVP, 4 CPU threads, `-gpu cuda`, warm JIT cache; the plugin is built for sm_75 and
JIT-compiled for sm_120 on the first run after each rebuild, which took 3.6 s once):

| system (spherical n) | CPU ERI build | device ERI kernels | GPU ERI path incl. host tables + context | max abs tensor diff | E diff |
|---|---:|---:|---:|---:|---:|
| H2O (24) | 2.8 ms | 18.1 ms | 98 ms (context+module load 79 ms) | 8.9e-16 | 0 (8 dp) |
| formaldehyde (38) | 13.1 ms | 33.5 ms | 114 ms (context 79 ms) | 5.8e-15 | 0 (8 dp) |
| benzene (114) | 612 ms | 71.9 ms | 170 ms (host tables 16 ms, context 81 ms) | 5.8e-15 | 0 (8 dp) |
| naphthalene (180) | 3.1 s (stored, `scf_direct off`) | 251 ms | 389 ms (host tables 50 ms, context 79 ms) | 5.8e-15 | 0 (8 dp) |

Tensor diffs from `CURCUMA_QM_GPU_ERI_CHECK=1` (download + element-wise compare with the CPU tensor, same
screening; also with `-qm.eri_screening 0` and `-qm.cartesian_d true` on H2O/formaldehyde/benzene, all
<= 5.8e-15, 0 elements above 1e-12). Whole single points (wall): benzene CPU 1.46 s -> Q-G1 1.06 s -> Q-G2
0.45 s; naphthalene CPU (auto = integral-direct) 25.6 s -> Q-G2 0.95 s, E = -383.07275442 Eh = PySCF
-383.0727544204 for the same hand-built D2h geometry. In this session the Q-G1 benzene J/K builds took
10.8 ms in total (12 iterations); the Q-G1 paragraph is corrected accordingly. Small molecules gain nothing: the ~80 ms context + module load dominates, the ERI kernels
themselves are latency-bound (one quartet's primitives run serially). `ctest -L qm_gpu`: `qm_gpu_eri_*`
(tensor <= 1e-12, energy <= 1e-9 vs CPU, incl. cartesian, unscreened, HF-3c/MINIX) and `qm_gpu_jk_*`.
`-opt -method hf-3c -gpu cuda` (formaldehyde) reaches the same minimum energy as the CPU (6 dp printed).

Not tested / not implemented: other GPUs and compute capabilities, ROCm/Vulkan (no QM entry point),
f functions (CPU fallback), systems beyond naphthalene (the tensor must fit into device memory; ~260
functions on a 32 GB card), multi-GPU, the analytic gradient (2e part stays on the CPU), timings with an
sm_120-native build. The block kernel for d classes is a first, unoptimised version (benzene: ~2/3 of the
kernel time); the thread kernel for classes up to 81 components did not compile in reasonable time.

**Q-G3 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: integral-direct J/K on the
device, no stored tensor (`cuda/qm_gpu_direct.cu`, entry `curcuma_cuda_create_qm_direct_jk`, loader
`gpu_plugin::createQmDirectJK`, hook `QMEngine::buildJK` in `qm_scf.cpp`). The integral kernels are the
Q-G2 ones, moved to `cuda/qm_gpu_eri_kernels.cuh` and instantiated with `DIRECT = true` in their own
TU (nvcc needs ~9 min for the thread-kernel instances of one mode; two TUs compile in parallel). Per
Fock build: a kernel forms Pmax(A,B) = max|P| per shell block; every canonical quartet of the pair
classes (X, Y) is enumerated on the device as the rectangle pairs(X) x pairs(Y) (non-canonical half
returns at once; no O(npair^2) quartet list, `buildEriGpuTables(..., list_quartets = false)`); the
density-weighted screen is the expression of `qmint::DirectJK` in the same order, so device and CPU
compute the same quartets; each quartet block is contracted at once into Jt, Kt (weight deg/8, each
target summed over two indices, then one FP64 `atomicAdd`); the host symmetrises. The incremental
builds on dP and the full rebuild every 8 iterations are `runSCF()`'s, unchanged. Device memory:
4 n^2 + nshell^2 doubles plus the pair tables (measured 6.3 MB for benzene, 12.6 MB for pyrene,
`cudaMemGetInfo` difference, 2 MB granularity).

Selection (`-gpu cuda|auto`): `-qm.gpu_direct auto` (default) runs every integral-direct SCF on the
device, i.e. `scf_direct on`, or `scf_direct auto` when the stored tensor fits neither
`-qm.eri_max_memory_mb` on the host nor the device; a stored tensor that fits the device is still
preferred (faster, below). `on` forces the device-direct path for every `-gpu` run, `off` keeps
integral-direct SCF on the CPU. `-qm.gpu_max_memory_mb` caps the device memory assumed by the
stored-tensor fit test (for testing the not-fitting case on small molecules). Fallback to the CPU
`DirectJK` on any error, counted as `qm_direct`; a stored tensor that does not fit is no longer
counted as a `qm_eri` fallback when the device-direct path follows.

Validation (HF/def2-SVP; J/K from `CURCUMA_QM_GPU_JK_CHECK=1`, first build = full SAD density;
energies printed to 10 decimals; PySCF 2.14 `RHF`, def2-svp, `conv_tol 1e-12`). Naphthalene,
anthracene, pyrene are planar hexagonal-lattice geometries (r(CC) 1.40, r(CH) 1.09 A), not library
structures:

| system (n) | quartets dev = CPU | max dJ, dK vs CPU direct | vs stored device tensor | E (CPU stored/direct = device stored/direct) | PySCF |
|---|---:|---:|---:|---:|---:|
| H2O (24) | 3081 | 4.4e-15, 4.0e-15 | 3.6e-15, 1.8e-15 | -75.9534444233 | -75.9534444233 |
| formaldehyde (38) | 13734 | 8.9e-15, 1.5e-14 | 5.3e-15, 1.3e-12 | -113.7703303142 | -113.7703303143 |
| benzene (114) | 507551 | 1.8e-14, 1.5e-14 | 3.0e-12, 3.1e-12 | -230.5359772878 | -230.5359772878 |
| naphthalene (180) | not run | not run | not run | -383.0665753163 | -383.0665753163 |
| anthracene (246) | 4035829 | 4.3e-14, 1.7e-14 | 1.4e-11, 5.7e-12 | -535.5888176085 | -535.5888176085 |
| pyrene (274, tensor 45 GB) | 5658915 | 5.3e-14, 2.0e-14 | does not fit | -611.2999868478 | -611.2999868478 |

E_elec of the device-direct run equals the CPU run to all 10 printed decimals in every case above and
in the ctests (cartesian 6d, unscreened, HF-3c/MINIX included). The 1e-12 to 1e-11 differences to the
stored tensor are the screening difference (Schwarz only vs density-weighted); CPU direct vs CPU
stored shows the same 3.0e-12 on benzene. `-opt -method hf-3c -gpu cuda -qm.gpu_direct on`
(formaldehyde) reaches -113.115469 Eh, the reference minimum -113.1154688146.

Timing (RTX 5090, `-qm.threads 4` unless noted, warm JIT cache; per build = J/K time / builds, a mix
of full and incremental builds that is the same for all paths; device-direct = kernels plus
transfers as reported by the builder, excluding the ~70-130 ms setup):

| system (n) | iter | CPU stored | CPU direct | Q-G2 device stored | Q-G3 device direct | whole run: CPU direct / Q-G2 / Q-G3 |
|---|---:|---:|---:|---:|---:|---:|
| H2O (24) | 12 | 0.1 ms | 1.5 ms | < 0.1 ms | 6.9 ms | 0.04 / 0.30 / 0.35 s |
| formaldehyde (38) | 16 | 0.2 ms | 7.3 ms | < 0.1 ms | 11.2 ms | 0.15 / 0.31 / 0.45 s |
| benzene (114) | 12 | 56.6 ms | 348 ms | 0.9 ms | 49.0 ms | 4.64 / 0.44 / 0.97 s |
| naphthalene (180) | 16 | host limit | 1496 ms | 6.7 ms | 198 ms | 26.0 / 0.94 / 3.92 s |
| anthracene (246) | 19 | host limit | 3763 ms | 28.3 ms | 494 ms | 76.6 / 2.23 / 10.8 s |
| pyrene (274) | 21 | host limit | 5652 ms (24 threads: 1487) | does not fit | 765 ms | 126.5 (24 threads: 34.0) / - / 18.0 s |

Pyrene with plain `-gpu cuda` (auto) selects the device-direct path by itself (18.2 s, no fallback
reported); with `-qm.scf_direct off` it falls back (counted) to a 45 GB host tensor (58.6 s). Where
the tensor fits, device-direct is 17-30x slower per build than the stored device tensor (one ERI
evaluation per build instead of once per geometry), so `auto` keeps the stored path. Against the CPU
direct path it is ~7x faster per build at 4 CPU threads and 1.9x at 24 threads (pyrene); below ~40
functions it is slower than the CPU (81 class launches per build). The cost is the Q-G2 quartet
kernel itself: naphthalene device-direct computes 71.7 M quartets in 3.37 s over 17 builds, the Q-G2
tensor build 5.2 M quartets in 250 ms, both ~21 quartets per microsecond, so the digest adds little.

Not tested: other GPUs and compute capabilities, ROCm/Vulkan (no entry point), f functions (CPU
fallback), systems beyond 274 functions, other basis sets than def2-SVP and MINIX, multi-GPU, the
gradient (2e part stays on the CPU), an sm_120-native build. Not reproducible bit for bit: the
atomics make the last bits of J/K schedule-dependent; the computed quartet count of the incremental
builds differed by 13 of 352461446 between two pyrene runs (the CPU direct path differs too between 4
and 24 threads: 352461552 vs 352462036). Not implemented: compaction of the screened quartets (every
rectangle element is launched), device-side symmetrisation, the resident SCF (Q-G4), FP32 (Q-G5); the
`auto` choice between stored and direct looks only at memory, not at speed.

Caveats carried over: speed-ups are hardware-specific and must be measured (the Vulkan
FP32 path gave no win); FMA contraction changes the last ulp, so tests compare energies
at a stated tolerance rather than bit patterns; and a device kernel that silently falls
back to the CPU must say so at verbosity >= 1.

An algorithmic alternative that also reduces cost on the CPU is **RI-J / density fitting**
(3-index integrals, O(n^3) memory); it changes the energy at the 1e-5 Eh level and is
therefore an opt-in method variant, not a replacement of the exact path.
