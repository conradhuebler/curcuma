# Native ab-initio QM (HF / HF-3c / KS-DFT): performance and GPU roadmap

> Status: 🤖 AI-authored. Sections 1-2 are **implemented and machine-tested** (CPU only,
> Sep 2026). Section 3 is a **proposed** plan, not implemented. No GPU code for
> the QM engine exists yet.

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
| Q-G1 | J/K from the stored tensor on the device (cuBLAS GEMV/GEMM), host SCF unchanged | Stage 1 (single offloaded hot spot) |
| Q-G2 | device ERI kernel: one work item per shell quartet (or per primitive quartet for high L), fed by the flat `ShellPair`/`PrimPair` tables of section 1 | Stage 3 (integrals on device) |
| Q-G3 | integral-direct J/K on the device (section 2), density-weighted screening, FP64 accumulation | `-gpu_multipole_otf`, screened integrals |
| Q-G4 | device-resident SCF loop: Fock -> Löwdin-reduced eigensolve (reuse the xTB `cusolverDnDsyevd` path) -> density -> DIIS; host polls O(1) scalars per iteration | Stage 6 (fully resident loop) |
| Q-G5 | mixed precision far from convergence (FP32 ERI/J/K with FP64 correction on dP), never accepting convergence on an FP32 step | X-AP3 (FP32 eigensolve early in the SCF) |
| Q-G6 | multi-GPU: distribute shell-quartet batches over the device pool | [MULTI_GPU.md](MULTI_GPU.md) |

Caveats carried over: speed-ups are hardware-specific and must be measured (the Vulkan
FP32 path gave no win); FMA contraction changes the last ulp, so tests compare energies
at a stated tolerance rather than bit patterns; and a device kernel that silently falls
back to the CPU must say so at verbosity >= 1.

An algorithmic alternative that also reduces cost on the CPU is **RI-J / density fitting**
(3-index integrals, O(n^3) memory); it changes the energy at the 1e-5 Eh level and is
therefore an opt-in method variant, not a replacement of the exact path.
