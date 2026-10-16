# Native ab-initio QM (HF / HF-3c / KS-DFT): performance and GPU roadmap

> Status: 🤖 AI-authored. Sections 1-2 are **implemented and machine-tested** (CPU only,
> Sep 2026). Section 3 is a plan; stages Q-G1 to Q-G4, Q-G7 (2e gradient) and Q-G8 (1e gradient)
> are implemented (Oct 2026, AI-generated, machine-tested on one RTX 5090 only), Q-G5 and Q-G6 are
> proposals.

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
| Q-G4 | **done Oct 2026** (see below): device-resident SCF loop: Fock -> Löwdin-reduced eigensolve (`cusolverDnDsyevd`) -> density -> DIIS; host polls O(1) scalars per iteration | Stage 6 (fully resident loop) |
| Q-G5 | mixed precision far from convergence (FP32 ERI/J/K with FP64 correction on dP), never accepting convergence on an FP32 step | X-AP3 (FP32 eigensolve early in the SCF) |
| Q-G6 | multi-GPU: distribute shell-quartet batches over the device pool | [MULTI_GPU.md](MULTI_GPU.md) |
| Q-G7 | **done Oct 2026** (see below): two-electron part of the analytic RHF gradient on the device, integral-direct | Stage 3 kernels reused for the gradient |
| Q-G8 | **done Oct 2026** (see below): one-electron part of the analytic RHF gradient on the device, from the Q-G7 tables | Stage 3 kernels reused for the gradient |

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
functions on a 32 GB card), multi-GPU, the analytic gradient (2e part on the CPU at the time; on the device since Q-G7), timings with an
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
rectangle element is launched), device-side symmetrisation, the resident SCF (Q-G4, added later, below), FP32 (Q-G5); the
`auto` choice between stored and direct looks only at memory, not at speed.

**Q-G4 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: device-resident RHF SCF loop
(`cuda/qm_gpu_scf.cu`, entry `curcuma_cuda_create_qm_scf`, seam `qm_gpu_scf_backend.h` -- a separate header
so that the ~9 min kernel TUs need not recompile when it changes -- loader `gpu_plugin::createQmScf`, driver
`QMEngine::runSCFDevice` in `qm_scf.cpp`). The loop control stays on the host and is the host loop's: same
Fock cadence (full, or incremental on P - P_ref with a full rebuild every 8 builds in direct mode), same
convergence test ||P_new - P||_F < `scf_threshold`, same damping, same guess/warm start, same final full J/K
build (`finishConvergedSCF`, shared). The DIIS/ADIIS class was split into subspace logic
(`QMScfAccelerator`: slots, coefficients from the m x m Gram/ADIIS matrices) and storage (`HostDiisStore`
on the host, the device history in `QMGpuSCF`), so both loops solve the same small problems with the same
code. Per iteration on the device: F = H + J - 1/2 K with J/K from the existing builders through the new
`QMGpuJK::buildDevice` (P, J, K device pointers; stored tensor Q-G1/Q-G2 or device-direct Q-G3, the latter
symmetrised on the device with the host expression); e = X^T (FPS - SPF) X and max|e| (cuBLAS); the Gram
matrix <e_i, e_j> and the ADIIS model reduced on the device (deterministic two-stage sums), solved on the
host (m <= 8); F_use = sum c_i F_i; F' = X^T F_use X, `cusolverDnDsyevd`, C = X C'; P_new = 2 C_occ C_occ^T
(`cublasDsyrk` + mirror, exactly symmetric); ||P_new - P|| and E = 1/2 Tr(P_new (H + F)). Transfers per
iteration: three or four scalars, the m x m matrices, m coefficients; only the converged P, C and eps come
back (the final F from the shared full build). X is built on the host once per geometry (Löwdin S^-1/2, as
before). The state object (matrices, cuBLAS/cuSOLVER handles, warmed up once) survives geometry changes and
is rebuilt only when n changes.

Selection: `-qm.gpu_scf auto` (default) uses the device loop when J/K run on the device and the basis has at
least `-qm.gpu_scf_min_nbf` (200) functions; `on` for any size (no device J/K builder = counted fallback);
`off` keeps the host loop. Any device error falls back to the host loop from the same start density, counted
as `qm_scf`. Checked: `-gpu cuda -qm.gpu_scf on -qm.gpu_eri false -qm.scf_direct on -qm.gpu_direct off`
(J/K on the CPU) reports `17 x qm_scf` over an `-opt` and reaches the same minimum. No numerical path depends
on the print level; `CURCUMA_QM_SCF_TRACE=1` prints one line per converged SCF (loop, n, iterations, start,
E_elec) where `-opt` / batch `-sp` silence the method.

Validation (HF/def2-SVP, `-qm.threads 4`; E_tot = E_elec + E_nn from the 10-decimal printout; ORCA 6.1
`! HF def2-SVP TightSCF NoAutoStart`; CPU = `-gpu none -qm.threads 24`). Device loop (`-qm.gpu_scf on`) and
host loop (`off`) use the same J/K builder (stored device tensor; pyrene device-direct):

| system (n) | iter dev / host / CPU | E_tot device loop | = host loop, = CPU | ORCA 6.1 | dev - ORCA |
|---|---:|---:|---:|---:|---:|
| H2O (24) | 12 / 12 / 12 | -75.9534444233 | yes, yes (10 dp) | -75.953444424 | 7e-10 |
| formaldehyde (38) | 16 / 16 / 16 | -113.7703303142 | yes, yes | -113.770330312 | -2.2e-9 |
| benzene (114) | 12 / 12 / 12 | -230.5359772878 | yes, yes | -230.535977292 | 4.2e-9 |
| naphthalene (180) | 16 / 16 / 16 | -383.0665753163 | yes, yes | -383.066575319 | 2.7e-9 |
| anthracene (246) | 19 / 19 / 19 | -535.5888176085 | yes, yes | -535.588817613 | 4.5e-9 |
| pyrene (274) | 21 / 21 / 21 | -611.2999868478 | yes, yes | -611.299986845 | -2.8e-9 |

ctests `qm_gpu_scf_*` (E_elec device loop vs host loop vs CPU <= 1e-9, equal iteration counts): H2O,
formaldehyde, benzene, benzene device-direct (incremental builds, 12241450 quartets), H2O cartesian 6d,
H2O HF-3c/MINIX, formaldehyde `-qm.scf_mode adiis` (17 iterations), H2O `plain` damping (43 iterations); all
E_elec identical to 10 decimals. `-opt -method hf-3c -gpu cuda -qm.gpu_scf on` (formaldehyde): 17 SCFs, all
in the device loop, 16 warm-started, iteration counts identical to host loop and CPU, one E_elec of 17
differs by 1e-10 (last printed digit), final energy -113.11546881 Eh (reference minimum -113.1154688146).
An unconverged SCF (benzene, `-qm.scf_max_iterations 4`) is reported as an error by both loops with the same
diagnostic energy (-230.50392362 Eh). A 4-frame batch `-sp` (naphthalene, perturbed naphthalene, anthracene, perturbed anthracene; state rebuilt at
the n change): identical E_elec and iterations, device vs host loop.

Timing (RTX 5090, warm JIT cache, ms; "rest" = SCF loop minus Fock builds, per iteration; host eigensolve =
Eigen/LAPACK on 4 threads, device = X^T F X + syevd; setup = SCF state incl. cuSOLVER/cuBLAS warm-up, once
per process and basis size, not in the loop time):

| system (n) | iter | loop host / device | Fock host / device | eigensolve per iter host / device | rest per iter host / device | device setup |
|---|---:|---:|---:|---:|---:|---:|
| H2O (24) | 12 | 0.6 / 7.1 | 0.2 / 0.2 | 0.02 / 0.48 | 0.03 / 0.57 | 72 |
| formaldehyde (38) | 16 | 2.2 / 17.7 | 0.5 / 0.3 | 0.09 / 0.97 | 0.11 / 1.09 | 71 |
| benzene (114) | 12 | 28.5 / 51.4 | 10.7 / 9.9 | 1.20 / 3.28 | 1.48 / 3.46 | 72 |
| naphthalene (180) | 16 | 170 / 176 | 108 / 106 | 2.79 / 4.14 | 3.84 / 4.38 | 78 |
| anthracene (246) | 19 | 698 / 635 | 537 / 531 | 5.71 / 5.14 | 8.45 / 5.45 | 79 |
| pyrene (274, direct) | 21 | 16206 / 16061 | 15956 / 15930 | 7.73 / 5.92 | 11.90 / 6.25 | 75 |

Whole single points are unchanged within noise (benzene 0.45 / 0.54 s, anthracene 2.23 / 2.24 s, pyrene
18.13 / 18.13 s, two pyrene repeats agree within 0.1 s). The device loop is slower below ~200 functions:
`cusolverDnDsyevd` is latency-bound (standalone, same card: 0.44 ms at n = 24, 3.1 ms at 114, 5.8 ms at 274;
`syevj` was 2-4x slower and is not used), and the state setup costs ~75 ms once. Above that it saves 3-6 ms
per iteration (host DIIS matrix products and eigensolve), which is < 1 % of an SCF whose J/K build is
integral-direct; the stage matters only if J/K become much cheaper (Q-G5) or n grows well beyond 274.

Not tested: other GPUs (a consumer card has a low FP64 rate; a data-centre card may change the balance,
not measured), ROCm/Vulkan (no entry point), n > 274, other basis sets than def2-SVP and MINIX, open shell (not
implemented anywhere), linear-dependent bases (X built on the host as before, no canonical orthogonalisation
with dropped vectors), multi-GPU, the device loop inside MD.
Not implemented: the DIIS subspace solve on the device (m x m, host), device-side gradient (CPU at the time; 2e part on the device since Q-G7),
a speed-based auto choice (threshold only, measured on this card with 4 CPU threads; with more CPU threads
the host eigensolve gets faster and the crossover moves up).

**Q-G7 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: two-electron part of the analytic
RHF gradient on the device (`cuda/qm_gpu_grad.cu`, entry `curcuma_cuda_create_qm_grad2e`, seam
`qm_gpu_grad_backend.h` -- separate so that the ~9 min kernel TUs need not recompile -- loader
`gpu_plugin::createQmGrad2e`, hook `QMEngine::computeGradient` / `gpuGradientTwoElectron` in `qm_engine.cpp`).
Same quantity as `qmint::gradientTwoElectron`: g = 1/2 sum D (mn|ls)^x with D = P P - 1/4 (P P + P P) of the
cartesian density, canonical shell quartets weighted 1/2 deg, derivatives on A, B, C, the one on D from
translational invariance, one-atom quartets skipped. Integral-direct (no derivative-integral or ERI
tensor; device memory P, Pmax plus pair tables, 25 MB for pyrene). Tables: `qmint::buildEriGradGpuTables` =
the Q-G3 tables in the cartesian basis with the Hermite tables rebuilt one step higher on both shells
(the centre derivative 2a g(l+1) - l g(l-1) then is a modified 1D Hermite coefficient; d-shell
derivatives reach f-type orders through these tables, so s, p, d run on the device, nothing falls back
for def2-SVP/MINIX). Contraction: per ket primitive pair Dket_ab(h) = sum_cd D c c E E E, per primitive
quartet K_ab(tuv) = sum_h Dket_ab(h) R(tuv + h) contracted with the derivative bra coefficients (d/dA,
d/dB); d/dC from the same routine on (CD|AB). One warp per shell quartet; "primitive mode" (lanes over
primitive quartets, R in thread-local memory) when the bra has < 32 Hermite entries; per-block shared
accumulation, FP64 atomics (last bits schedule-dependent). Screen: Q_AB Q_CD (Pmax_AB Pmax_CD + 1/4 (...)) <
`eri_screening`/100 (the derivative integrals are not bounded by the Schwarz factors; at 1e-12 instead
of 1e-14 benzene was off by 2.6e-10 vs 6.2e-12 Eh/Bohr, for 5 % fewer quartets); `-qm.eri_screening 0`
evaluates every quartet like the CPU. The 1e part stayed on the CPU (on the device since Q-G8).

Selection: `-qm.gpu_grad auto` (default) and `on` both use the device whenever `-gpu cuda|auto` is active
(the device was faster than 4 CPU threads at every measured size, H2O included); `off` keeps the CPU. Any
failure falls back to the CPU, counted as `qm_grad`. Independent of the SCF path: checked after the stored
device tensor (Q-G2, default), the uploaded CPU tensor (Q-G1, `-qm.gpu_eri false`), device-direct J/K
(Q-G3) and the device-resident loop (Q-G4); cartesian 6d, unscreened, HF-3c/MINIX. Diagnostics:
`CURCUMA_QM_GPU_GRAD_CHECK=1` also runs the CPU 2e gradient from the same density and prints the
difference; `CURCUMA_QM_SCF_TRACE=1` prints one line per gradient (`2e=device|cpu`, times), visible in
`-opt`; `CURCUMA_QM_GPU_GRAD_PROFILE=1` serialises the class launches and times them.

Validation (HF/def2-SVP, `scf_threshold 1e-10`, 24 CPU threads; "2e dev - CPU" from the same density via
the check switch; totals from `-sp -dump_gradient`, Eh/Bohr; PySCF 2.14 `RHF` conv_tol 1e-12,
`nuc_grad_method()`; ORCA 6.1 `! HF def2-SVP TightSCF NoAutoStart EnGrad`). CPU reference runs used the
stored host tensor (`-qm.scf_direct off`) for naphthalene to pyrene; pyrene's device run had to go
device-direct (tensor does not fit the card) and was converged to 1e-8 only:

| system (natoms) | 2e dev - CPU | sum of dev 2e | total dev run - CPU run | dev - PySCF | dev - ORCA |
|---|---:|---:|---:|---:|---:|
| H2O (3) | 3.1e-15 | 4.7e-15 | 3.6e-15 | 1.0e-9 | 1.9e-7 |
| formaldehyde (4) | 1.2e-14 | 8.0e-15 | 1.2e-14 | 1.6e-9 | 1.4e-6 (VeryTightSCF: 2.8e-7) |
| benzene (12) | 6.2e-12 | 1.3e-13 | 6.2e-12 | 5.1e-11 | 1.7e-7 |
| naphthalene (18) | 5.0e-12 | 2.5e-13 | 5.1e-12 | 3.5e-10 | not run |
| anthracene (24) | 3.5e-12 | 6.0e-14 | 3.4e-12 | 7.2e-10 (*) | 8.6e-8 |
| pyrene (26) | 1.5e-11 | 2.4e-13 | 3.8e-10 (SCF 1e-8) | 2.1e-9 | not run |

(*) PySCF flags anthracene as not converged against `conv_tol_grad 1e-7` (also after `newton()`); its energy
equals curcuma's to 1e-12. The CPU gradient has the same distance to PySCF and ORCA as the device one in
every row; the ORCA distance shrinks with ORCA's SCF threshold (formaldehyde 1.4e-6 -> 2.8e-7), so it is
ORCA's convergence, not the gradient. Sum of the total device gradient <= 5e-12. Finite difference of the
device-run energy (H2O, Angstrom displacements, central, 15-decimal energies): 2.1e-8 Eh/Angstrom at h = 1e-4,
5.2e-7 at 5e-4 (h^2 scaling, truncation). ctests `qm_gpu_grad_*` (9): H2O (+FD), formaldehyde, benzene,
benzene device-direct, benzene device-resident SCF, H2O Q-G1 upload, H2O cartesian, formaldehyde unscreened
(<= 1e-12), H2O HF-3c.

Timing (single gradient, warm JIT cache, ms; device = whole 2e call incl. host tables and setup; the 1e part
is on the CPU in both columns and listed separately):

| system (n cart.) | 2e CPU 4 thr | 2e CPU 24 thr | 2e device | of which kernels | 1e CPU 4 thr | speed-up vs 24 thr |
|---|---:|---:|---:|---:|---:|---:|
| H2O | 12.0 | 5.8 | 7.9 | 6.9 | 2.7 | 0.7x |
| formaldehyde | 77.2 | 28.5 | 19.3 | 17.2 | 7.5 | 1.5x |
| benzene | 5301 | 1671 | 188 | 172 | 121 | 8.9x |
| naphthalene | 34508 | 10863 | 828 | 790 | 415 | 13x |
| anthracene | 120956 | 35818 | 2032 | 1957 | 1049 | 18x |
| pyrene | 192100 | 56277 | 3134 | 3038 | 1421 | 18x |

Whole `-opt` (wall, same start geometry, identical optimiser path and final geometry to the 6 printed
decimals): formaldehyde HF-3c 17 gradients, GPU 0.47-0.57 s vs CPU 0.24 s (4 threads; small system, plugin
and context dominate), final -113.115469 Eh = reference minimum; distorted benzene HF/def2-SVP 25 steps, 29
gradients, GPU 13.4 s vs CPU 203.5 s (4 threads) / 76.5 s (24 threads), final -230.536045 Eh in all three. In
the benzene GPU run the 2e gradients took 5.4 s and the CPU 1e part 3.3 s: the 1e part is now the next
cost (it is O(n^2 natoms) with no screening).

Not tested: other GPUs (consumer card, low FP64 rate; a data-centre card will shift the numbers), ROCm/Vulkan
(no entry point), f functions (none in the basis files; the table builder refuses L > 2 and the CPU runs),
more than 26 atoms / 274 functions, basis sets other than def2-SVP and MINIX, open shell (not implemented
anywhere), MD with the device gradient, `-gpu_strict`. Not implemented: 1e gradient on the device (added as
Q-G8), primitive screening, batching of several primitive quartets per warp in entry mode (the warp runs one at a time;
Boys on one lane), reuse of the device tables between SCF and gradient (rebuilt per gradient, 15-95 ms).
Observed while validating (not caused by this stage): the integral-direct SCF (CPU `DirectJK` and Q-G3 alike)
stalls near ||dP|| ~ 1e-8 for naphthalene and pyrene and does not reach `scf_threshold` 1e-9 / 1e-10 within
100 iterations, so tight gradients of direct runs are limited by the SCF, not by the gradient.

**Q-G8 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: one-electron part of the analytic
RHF gradient on the device (`cuda/qm_gpu_grad1e.cu`, entry `curcuma_cuda_create_qm_grad1e`, seam
`qm_gpu_grad1e_backend.h`, loader `gpu_plugin::createQmGrad1e`, hook `QMEngine::computeGradient` /
`gpuGradientOneElectron`). Same quantity as `qmint::gradientOneElectron`: sum P dT + sum P dV (basis part on the
function centres, operator/Hellmann-Feynman part on the moving nucleus via d/dC = -(d/dA + d/dB)) - sum W dS,
Eh/Bohr; the nuclear repulsion gradient stays on the host (O(natoms^2), as before). Tables: the Q-G7 ones
(`qmint::buildEriGradGpuTables`), now built once per gradient and handed to both device parts; the 1e kernel
adds only shell centres, nuclear positions and charges. Integrals by McMurchie-Davidson from the raised
Hermite tables: overlap s(i,j) = E^ij_0, kinetic from the identity of `primitiveKinetic` (its A derivative
needs s(LA+2, j), taken from E^{i+1,j}_0 = X_PA E^ij_0 + E^ij_1), attraction from R^0_tuv(p, P-C) with the
device Boys table (`boysDev`, unchanged). One thread per (canonical shell pair, nucleus) plus one slot per pair
for S and T, each contracting with P and W at once; shared-memory accumulation, FP64 atomics (last bits
schedule-dependent); no screening (the CPU has none). s, p, d shells on the device (the derivative orders reach
f through the raised tables); L > 2 is refused by the table builder and the CPU runs (counted fallback).

Selection: `-qm.gpu_grad1e auto` (default) and `on` use the device whenever `-gpu cuda|auto` is active (it
was not slower than 4 or 24 CPU threads at any measured size, H2O included); `off` keeps the CPU. Independent of
`-qm.gpu_grad` (with `gpu_grad off` the tables are built for the 1e part alone). Any failure falls back to the
CPU, counted as `qm_grad1e` (checked: new binary with the Q-G7 plugin, which lacks the entry point, reports `1
x qm_grad1e` and gives the CPU gradient, 1.8e-14 Eh/Angstrom from the CPU run). `CURCUMA_QM_GPU_GRAD_CHECK=1`
also runs the CPU 1e gradient from the same P, W and prints `QM GPU grad1e check: ...`; the trace line of
`CURCUMA_QM_SCF_TRACE=1` now reads `1e=device|cpu 2e=device|cpu tables_ms 2e_ms 1e_ms`.
Not done: fusing the 1e and 2e uploads into one call (would mean recompiling `qm_gpu_grad.cu`); each part
uploads its own copy of the tables (21 MB for pyrene) and P (plus W for 1e), measured transfers inside the 1e
call <= 0.13 ms.

Validation (HF/def2-SVP, `scf_threshold 1e-10` (pyrene 1e-8, device-direct SCF), 24 CPU threads; "1e dev - CPU"
from the same P, W via the check switch; totals Eh/Bohr from `-sp -dump_gradient`; PySCF 2.14 as in Q-G7; ORCA 6.1
`! HF def2-SVP ExtremeSCF NoAutoStart EnGrad`, last orbital gradient <= 1e-9 in every run). CPU reference runs as
in Q-G7 (stored host tensor for naphthalene to pyrene):

| system (natoms) | 1e dev - CPU | sum of dev 1e | total dev run - CPU run | dev - PySCF | dev - ORCA ExtremeSCF |
|---|---:|---:|---:|---:|---:|
| H2O (3) | 5.0e-14 | 2.7e-15 | 5.0e-14 | 1.0e-9 | 7.0e-9 |
| formaldehyde (4) | 9.1e-13 | 2.2e-15 | 9.0e-13 | 1.6e-9 | 8.7e-9 |
| benzene (12) | 1.4e-14 | 1.4e-14 | 6.2e-12 | 5.1e-11 | 4.1e-9 |
| naphthalene (18) | 1.6e-12 | 1.8e-14 | 5.2e-12 | 3.5e-10 | 4.2e-9 |
| anthracene (24) | 2.2e-12 | 2.8e-14 | 4.9e-12 | 7.2e-10 (*) | 4.3e-9 |
| pyrene (26) | 1.5e-12 | 8.9e-15 | 3.7e-10 (SCF 1e-8) | 2.1e-9 | 3.9e-9 |

(*) as in Q-G7. The "total dev - CPU" column is now dominated by the 2e part and the SCF, not the 1e part. The
CPU 1e gradient itself is translationally invariant only to 5e-12 (sum_cpu); the device skips the pairs whose
contributions cancel exactly (two-centre terms on one atom, attraction of a one-atom pair by its own nucleus).
**ORCA residual explained**: the 1e-7 to 1e-6 differences to ORCA quoted in Q-G7 came from ORCA's SCF orbital
gradient, not from integral thresholds: `VeryTightSCF` stops at TolG 2e-6 (formaldehyde: last orbital gradient
1.3e-6, distance 2.76e-7), `ExtremeSCF` (TolG 1e-9) gives 8.74e-9, and adding `Thresh 1e-14 TCut 1e-16` to a tighter
SCF changes nothing (8.72e-9). The remaining 4-9e-9 Eh/Bohr to ORCA (PySCF: 1e-9 to 2e-9) is not explained.
Finite difference of the device-run energy (1e and 2e on the device; H2O, central, Angstrom displacements,
15-decimal energies): 4.4e-9 / 2.1e-8 / 8.2e-8 Eh/Angstrom at h = 5e-5 / 1e-4 / 2e-4 (h^2 scaling, truncation).
ctests `qm_gpu_grad1e_*` (6, 1e dev - CPU <= 1e-10 required): H2O (+FD), formaldehyde with the 2e part on the CPU
(`-qm.gpu_grad off`), benzene, benzene device-direct SCF, H2O cartesian 6d, H2O HF-3c/MINIX; the 9 `qm_gpu_grad_*`
tests check the 1e line too when it is present (default). `ctest -L "qm_1e|qm_2e|qm_hf3c|qm_gpu|qm_grad"`: 97 of 99
pass; the two `qm_opt_hf3c_*` failures are the missing numpy of the system python (both pass with the venv python,
E(opt) - E(ref) 2e-12 / -4e-12 Eh).

Timing (single gradient, RTX 5090, warm JIT cache, ms; "1e device" = whole 1e call incl. object setup, uploads
and download; "tables" = host build of the shared device tables, previously counted inside the 2e call; 1e CPU
from runs with `-qm.gpu_grad1e off`):

| system (n cart.) | tables | 1e device | of which kernel | 1e CPU 4 thr | 1e CPU 24 thr | 2e device | speed-up 1e vs 24 thr |
|---|---:|---:|---:|---:|---:|---:|---:|
| H2O | 0.8 | 2.2 | 1.98 | 3.1 | 2.8 | 6.9 | 1.3x |
| formaldehyde | 1.9 | 2.3 | 2.08 | 7.7 | 4.1 | 17.3 | 1.8x |
| benzene | 14.9 | 2.8 | 2.14 | 117 | 33.6 | 175 | 12x |
| naphthalene | 37.2 | 4.0 | 2.85 | 414 | 114 | 789 | 28x |
| anthracene | 73.5 | 5.5 | 3.81 | 1046 | 252 | 1958 | 46x |
| pyrene | 93.2 | 6.3 | 4.09 | 1411 | 339 | 3042 | 54x |

The kernel has a floor of ~2 ms (H2O and benzene alike): one thread walks all primitive pairs and component pairs
of its (pair, nucleus) serially, so the slowest thread (a contracted s or a d pair) sets the time. Whole `-opt`
(distorted benzene HF/def2-SVP, `-qm.threads 4`, 25 steps, 29 gradients, identical step table to the 6 printed
decimals, final -230.536045 Eh): 13.4 s before (Q-G7 binary, 1e on the CPU 3.29 s) -> 10.2 s (1e on the device
42 ms in total, tables 0.41 s, 2e 5.06 s); `-qm.gpu_grad1e off` with the new binary 13.3 s. HF-3c `-opt`
(formaldehyde, 17 gradients, 16 warm-started SCFs): same final geometry as the CPU run to the printed 6 decimals,
-113.115469 Eh.

Not tested: other GPUs, ROCm/Vulkan (no entry point), f functions (none in the basis files; refused, CPU),
basis sets other than def2-SVP and MINIX, more than 26 atoms, open shell (not implemented anywhere), point charges
or effective nuclear charges other than Z (the operator uses the element number, as the CPU), more than 400 atoms
(global-atomics branch of the kernel never ran), MD with the device gradient, `-gpu_strict`. Not implemented:
splitting a (pair, nucleus) over several threads (would lower the ~2 ms floor), fused 1e/2e upload, reuse of the
device tables between gradients or between SCF and gradient.

Caveats carried over: speed-ups are hardware-specific and must be measured (the Vulkan
FP32 path gave no win); FMA contraction changes the last ulp, so tests compare energies
at a stated tolerance rather than bit patterns; and a device kernel that silently falls
back to the CPU must say so at verbosity >= 1.

An algorithmic alternative that also reduces cost on the CPU is **RI-J / density fitting**
(3-index integrals, O(n^3) memory); it changes the energy at the 1e-5 Eh level and is
therefore an opt-in method variant, not a replacement of the exact path.
