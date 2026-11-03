# Native ab-initio QM (HF / HF-3c / KS-DFT): performance and GPU roadmap

> Status: 🤖 AI-authored. Sections 1-2 are **implemented and machine-tested** (CPU only,
> Sep 2026). Section 3 is a plan; stages Q-G1 to Q-G4, Q-G7 (2e gradient) and Q-G8 (1e gradient)
> and Q-G5 (FP32 phase of the direct SCF and mixed-precision RI-JK) are implemented (Oct 2026, AI-generated,
> machine-tested on one RTX 5090 only), Q-G6 (multi-GPU) is a proposal. Section 4: Kohn-Sham DFT on the GPU
> (Oct 30 2026, AI-generated, machine-tested on one RTX 5090 only); 4.4 open shell (stage 6a), 4.5 f and g shells
> (stage 6b, Nov 2026).

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
Oct 25, 2026: Hermite-matrix CPU kernel (`mdQuartetT`, contraction and 5d folded into per-primitive-pair
matrices, compile-time class sizes) and a Schwarz-sorted quartet loop: one direct build 4-5x faster on one
thread, CPU HF C80H162 712 -> 231 s on 24 threads; energies unchanged to 10 decimals
([NATIVE_QM_IMPLEMENTATION.md](NATIVE_QM_IMPLEMENTATION.md#cpu-eri-kernel-hermite-matrix-mcmurchie-davidson-oct-2026),
[QM_TIMING_ORCA_PYSCF_2026-10.md](QM_TIMING_ORCA_PYSCF_2026-10.md)).
Same day: with `-gpu cuda` and `-qm.gpu_direct auto`, `ensureGpuJK` decides from the card's total memory
whether the Q-G2 tensor can fit before building `buildEriGpuTables` (was: tables first, ~8 s at C40H82 and
~160 s at C80H162 for nothing); HF C80H162 197 -> 30 s with the default flags.

**Incremental-build floor (Oct 16, 2026; AI-generated, machine-tested).** An incremental build on
dP is screened with the same absolute threshold as a full one (`Q_ab Q_cd max|dP| < eri_screening`),
so each build drops up to `eri_screening` per quartet however small dP is; the accumulated Fock matrix
carries a path-dependent error that does not shrink with dP. Measured (HF/def2-SVP, CPU, 16 threads,
diagnostic build comparing each Fock matrix with a full build and an unscreened one): naphthalene full
build vs unscreened 5.0e-11 Eh (fixed), accumulated incremental vs full 2e-11 to 1.2e-10 Eh (varies),
DIIS error max|e| stuck near 6e-10. ||dP||_F (AO basis) then floats at about 1e-8 (benzene), 3e-8
(naphthalene), 1e-7 to 3e-7 (anthracene, pyrene) and never reaches 1e-9 or 1e-10 (CPU direct and
device-direct alike). The floor moves with the screening threshold (naphthalene: 1e-8 to 7e-8 at 1e-12,
1e-10 to 1e-9 at 1e-14, below 1e-10 at 1e-16) and disappears with a full build every iteration (23 iterations,
the same as the stored tensor), so it comes from the incremental screening and not from DIIS or the
convergence test. The full build is a fixed function of P (its screen sees only P, the quartet count stays
constant near convergence), so full builds converge like the stored tensor.
Fix (`DirectFockCadence` in `qm_scf.cpp`, host and device loop): once ||dP|| (below 1e-4) has not
dropped below half its smallest earlier value for 2 iterations in a row, the rest of that SCF uses full
builds. Convergence test and `eri_screening` are unchanged; verbosity 2 prints one line when it fires.

| HF/def2-SVP, `-qm.scf_direct on` | before (60 it max) | after: iterations (switch at) | stored tensor | abs(E after - ORCA 6.1) |
|---|---|---|---|---|
| benzene (library `benzene.validation`), 1e-10, 4 threads | not converged | 24 (17) | - | PySCF: 1.1e-11 |
| naphthalene 1e-10 | not converged | 29 (21) | 23 | 2.7e-9 |
| anthracene 1e-10 | not converged | 30 (22) | - | 4.4e-9 |
| pyrene 1e-10, CPU | not converged | 34 (24) | - | 2.8e-9 |
| pyrene 1e-10, `-gpu cuda -qm.gpu_direct on` (device-resident loop) | not converged | 36 (24) | - | 2.8e-9 |

At the default `scf_threshold` 1e-6 the rule never fired (CPU: benzene, naphthalene, anthracene, pyrene;
device-direct: pyrene): same iteration counts and energies to 10 decimals; SCF loop times within the
run-to-run spread (pyrene CPU 44.4-44.6 s before, 44.6-44.8 s after; device 16.05 s both). The quartet
counts of repeated runs differ by up to ~700 in both versions (threaded summation order in J/K). Cost
when it fires: 2-8 iterations more than the stored tensor or a full-build-only run, those with full builds. The device 2e gradient
(screen `eri_screening`/100) is one build with the full density and has no such floor; device vs CPU
gradient after a direct SCF at 1e-10: 5.1e-12 (naphthalene), 1.5e-11 Eh/Bohr (pyrene).
`-qm.scf_mode adiis`, naphthalene 1e-10: not converged before, 27 iterations (switch at 21) after.
Not tested: systems larger than pyrene, plain mode, warm-started geometry sequences at tight
thresholds, other GPUs than the RTX 5090. The floor grows with size (benzene to pyrene: 1e-8 to 3e-7), so somewhat larger
systems may stall at the default 1e-6 (extrapolated, not measured); the rule applies there as well. Regression check:
`qm_direct_jk` part (4).

**DIIS history at the switch (Oct 16, 2026; AI-generated, machine-tested).** At the switch the DIIS/ADIIS history
now drops every entry whose Fock matrix came from an incremental build (`QMScfAccelerator::dropIncrementalEntries()`,
host and device loop); entries from the periodic full rebuilds stay. Variants compared with a temporary switch
(HF/def2-SVP, CPU direct, 24 threads, 4 repeats; device: pyrene, 5-8 repeats): clear the whole history (A) gives
nearly the same counts but up to +3 on device pyrene/ADIIS 1e-10; keeping the newest 1-3 entries (B) gains little or
loses (methanol 1e-10 keep3: 32 vs 26). Iterations before -> after (switch iteration in brackets):

| system (switch) | 1e-10 diis | 1e-10 adiis | 1e-8 diis | 1e-8 adiis |
|---|---|---|---|---|
| methanol (24) | 26 -> 25 | 21 (no switch) | 17 (no switch) | 17 (no switch) |
| benzene (17-18) | 26 -> 21 | 25 -> 20 | 22 -> 19 | 20 -> 18 |
| naphthalene (21) | 29 -> 25 | 26-27 -> 25-26 | 24 -> 22 | 23 -> 22 |
| anthracene (22-23) | 29-30 -> 29-30 | 31-32 -> 31-32 | 26 -> 24 | 28 -> 25 |
| pyrene CPU (24) | 34 -> 33 | 34 -> 34-35 | 32 -> 26 | 29 -> 28 |
| pyrene device loop / gpu_scf off (24) | 35-37 -> 34-36 | 35-36 -> 36-37 | 32 -> 26 | 29 -> 28 |

Methanol 1e-11: diis 28 -> 26-27, adiis 30-32 -> 27. The only case that got worse is pyrene ADIIS at 1e-10 (median
+1, single runs up to +2 on the device); J/K summation order makes repeated runs differ by about 1 iteration.
E_elec unchanged to 10 decimals in every run; default threshold 1e-6 unchanged (naphthalene 16/18, anthracene 19/19
iterations, diis/adiis, rule does not fire). Not tested: plain mode, systems larger than pyrene, warm-started
geometry sequences, other GPUs than the RTX 5090. The gap to the stored tensor remains (naphthalene 25 vs 23).

**Large systems, n ~ 500-2900 (Oct 19, 2026; AI-generated, machine-tested, RTX 5090 only).** Profile of the exact
GPU-direct SCF (`-qm.ri off -gpu cuda -qm.gpu_direct on -qm.gpu_scf on`, all-anti alkanes from
`scripts/local_mp2/lmp2_common.py alkane(n)`, HF/def2-SVP, `-qm.threads 1`): the J/K builds were 99 % of the SCF loop
(eigensolve n=1930: 98 ms per iteration, DIIS 42 ms); per build two costs grew faster than the quartet count: (a) the
launch of the full pairs(X) x pairs(Y) rectangle per class, O(npair^2) threads (C80H162: 467k shell pairs, the first
build with 182 M quartets took 29.9 s, a full one with 1.28 G quartets 48 s); (b) the block kernel of the classes with
more than 27 cartesian components (d shells, (pp|pp)): 13 % of the quartets but 83 % of the kernel time (C20H42,
`CURCUMA_QM_GPU_DIRECT_PROFILE=1`, 1-18 quartets/us against 160-860 for the thread kernel). On top, the old "sad" guess
(each atom's block of the molecular core Hamiltonian) started far off: ||dP|| above 50 for 9 (C40) to 14 (C80)
iterations, C80 27 iterations. Changes:

- quartet enumeration (`cuda/qm_gpu_direct.cu`): pairs of each class sorted by Schwarz factor; per build the host walks
  the bra list and keeps per bra the ket prefix with Q_i Q_j max|P| >= thr ("ragged rows", 25 ms per C80 build); only
  class pairs X <= Y are launched; the exact screen per quartet is unchanged, so quartet counts equal the CPU `DirectJK`.
- two-step McMurchie-Davidson kernel for the large classes (`eriHermiteKernel`): ket Hermite sums once per ket
  component pair and bra Hermite index, then the bra sums; 1.4x on those classes. Tried and dropped: one warp per
  quartet (no gain) and warps over primitive quartets with shared atomics (slower).
- `-qm.scf_guess atom` (new default): superposition of spherically averaged, fractionally occupied atomic RHF densities.
- `-qm.scf_early_screening 1e-10` (new default, direct SCF, cold start only): looser screening until ||dP|| <
  `-qm.scf_early_switch` (1e-3), then a full build with `eri_screening`; convergence only with `eri_screening`.
- `-qm.eri_prim_density true` (new default): primitive quartets screened against (eri_screening 1e-3)/pq (CPU and device).
- setup: nuclear-attraction matrix builds the Hermite tables once per primitive pair instead of per (pair, nucleus);
  bit-identical, C40 setup with 1 thread 24.4 -> 8.5 s.

| system (n) | old: iter, SCF loop, wall | new: iter, SCF loop, wall | quartets per full build | E_elec new - old |
|---|---|---|---|---|
| C20H42 (490) | 18, 22.5 s, 27 s | 13, 5.9 s, 7.9 s | 70 M | 0 (10 dp) |
| C40H82 (970) | 23, 152 s, 185 s (load 20 at start) | 13, 23.8 s, 31.6 s | 305 M | 0 |
| C60H122 (1450) | 18, 328 s, 406 s | 13, 53.6 s, 72.0 s | 704 M | 0 |
| C80H162 (1930) | 27, 1172 s, 1381 s | 13, 95.5 s, 129.6 s | 1.27 G | 0 |
| C120H242 (2890) | SCF diverges/oscillates, did not converge (killed after 20 iterations, ~1 h; earlier session 3.9 h) | 13, 216 s, 306 s | not printed | no reference |

Single runs, machine load 1-5 (another agent active). Wall-time exponents over n: old C40->C80 2.9 (C60->C80 4.3), new
2.05 (C80->C120 2.13); the J/K quartet count grows as n^2 (exact Coulomb, no multipole far field). C80 new: 10.0 G
quartets in 93.7 s of builds (107 quartets/us, old 24/us), final build 13 s, setup ~20 s (1e integrals with 1 thread;
`-qm.threads 24` cuts it). Step by step on C20 (SCF loop): 22.5 s -> 14.8 (enumeration + Hermite kernel) -> 7.4 (atom
guess, 13 iterations) -> 6.4 (early screening) -> 5.9 s (primitive density screen); C40: 45.5 -> 30.1 -> 23.8 s. Old code on
C120 (rerun this session): ||dP|| 4e2 to 6e3 over 20 iterations with no downward trend (iteration 5: 2.3e3, 13: 6.3e3,
20: 1.7e3), ~170 s per Fock build of which ~120 s rectangle launch (the first build, 0.4 G quartets, took 128 s);
killed after ~1 h. The cause of the earlier 3.9 h run is therefore the oscillating SCF from the old guess times the
slow builds. Energies: E_elec of benzene, naphthalene, anthracene (ideal acenes),
C8H18, C20-C80 equal the old code to all 10 printed decimals on the GPU-direct path and the CPU-direct path (CPU: 2.2-2.6x
faster from the guess alone); `ctest -L "qm_gpu|qm_2e|qm_1e|qm_hf3c|qm_grad|qm_mp2|qm_ri|qm_lmp2|qm_ccsd"` 259/259 (the two
`qm_opt_hf3c_*` with the venv python). Not tested: other elements than H, C, O at this size, 3D-compact systems (more
quartets per pair), other GPUs (a data-centre FP64 card shifts the balance), basis sets other than def2-SVP/MINIX,
warm-started geometry sequences of large systems, `-qm.scf_threshold` below 1e-6 at this size, ROCm/Vulkan (no entry
point). Not implemented: multipole far field for J (the n^2 wall), FP32 for small quartets (Q-G5), a faster kernel for
the large classes (still ~80 % of the build time). (The last two followed on Oct 20, below.)

**Thread-per-quartet McMurchie-Davidson kernels for every class (Oct 20, 2026; AI-generated, machine-tested, RTX 5090
only).** Profile of the Oct 19 code (C20H42, `CURCUMA_QM_GPU_DIRECT_PROFILE=1`): the block kernels of the classes with
more than 27 cartesian components ran at 2-24 quartets/us because one block per quartet synchronised several times per
primitive quartet (Boys by one thread, R level by level, ket sums, bra sums) with integer index decoding and global
loads in the dependency chain; the classes themselves need little FP64 work per quartet. New
`cuda/qm_gpu_direct_md.cu` (instances) and `cuda/qm_gpu_direct_md_impl.cuh` (templates): one thread per quartet (or a
compile-time slice of the ket components, chosen per class from the register state; the slice is `blockIdx.y`, so it
is uniform per block), no barrier; per primitive pair a setup table F of the Hermite expansion of every active
component pair with contraction coefficients and the 6d->5d transform folded in; per bra primitive the ket sums are
contracted over the ket primitives, then one bra sum; Hermite indices whose E coefficient vanishes for the component
pair are dropped at compile time (compile-time loops; a `#pragma unroll` form needed ~25 min of nvcc for 8 classes, the
compile-time form 4 min for all 45). Same screening and quartet counts, FP64 atomics as before. The MD kernels are
also faster than the old thread kernel on every small class (1.1-1.6x), so they run for all classes; the old kernels
remain for cartesian d shells (`-qm.cartesian_d`) and as diagnostics (`CURCUMA_QM_GPU_DIRECT_MD=0` old kernels,
`=large` old thread kernel for the small classes). Block size 32 (C40 SCF loop with the FP32 phase: 32: 4.54 s, 64: 4.59 s, 128: 4.88 s, 256: 4.84 s; FP64 only: 32:
5.79 s, 128: 5.88-6.21 s).

Per class, C20H42 (quartets/us, serialised profile launches, all-FP64 SCF, old = `CURCUMA_QM_GPU_DIRECT_MD=0`):
(pp|pp) 12.2 -> 101, (ps|dp) 13.7 -> 115, (pp|ds) 18.7 -> 113, (ps|pd) 10.5 -> 77, (sp|dp) 15.9 -> 129, (pp|dp) 8.5 ->
72, (pp|pd) 6.5 -> 50, (ss|dd) 17.4 -> 153, (ds|dp) 10.9 -> 88, (pp|dd) 4.3 -> 18, (pd|pd) 3.8 -> 12, (dd|dd) 1.9 ->
3.8 (the 28 large classes 2.0-8.8x); small classes 1.1-1.6x ((ss|ss) 875 -> 998, (ss|pp) 360 -> 444, (ps|pp) 187 ->
269). Kernel time of one SCF 6551 -> 1782 ms, of which the large classes 5293 -> 788 ms (81 % -> 44 %). The classes
with a split over 15-25 lanes ((pd|pd), (dd|dd), (dp|dp), (pd|dp)) gain least: every lane recomputes Boys and R.
Whole SCF (HF/def2-SVP, `-gpu cuda -qm.gpu_direct on -qm.gpu_scf on -qm.ri off`, `-qm.gpu_eri_precision fp64`; old =
same binary with `CURCUMA_QM_GPU_DIRECT_MD=0`, which is the Oct 19 kernel set; single runs, load ~1-2):

| system (n) | old: SCF loop, wall, quartets/us | new: SCF loop, wall, quartets/us | iterations | E_elec |
|---|---|---|---|---|
| benzene (114) | 289 ms, 0.8 s, 33 | 149 ms, 0.7 s, 75 | 10 | equal (10 dp) |
| naphthalene (180) | 1238 ms, 1.8 s, 39 | 379 ms, 0.9 s, 145 | 13 | equal |
| anthracene (246) | 3087 ms, 3.9 s, 41 | 806 ms, 1.5 s, 169 | 15 | equal |
| C20H42 (490) | 5.90 s, 7.9 s, 89 | 1.64 s, 3.1 s, 345 | 13 | equal |
| C40H82 (970) | 23.8 s, 31.6 s, 92 | 5.79 s, 11.5 s, 401 | 13 | equal |
| C60H122 (1450) | 53.6 s, 72.0 s, 93 | 12.4 s, 25.3 s, 429 | 13 | equal |
| C80H162 (1930) | 95.5 s, 130.1 s, 94 | 21.5 s, 45.8 s, 446 | 13 | equal |
| C120H242 (2890) | 216.3 s, 304.8 s, 94 | 51.1 s, 116.1 s, 434 | 13 | equal |

quartets/us = shell quartets computed in the SCF / sum of the build kernel times (all builds incl. the final one). The
wall time of the large alkanes is now about half setup (one-electron integrals and tables with the default thread
count; not touched here). Device memory of the F tables: 1 MB (benzene), 19 MB (C20), 139 MB (C80), 221 MB (C120). Tests: `ctest -L "qm_gpu|qm_2e|qm_grad"` 94/94 (the two
`qm_opt_hf3c_*` with the venv python). Not tested: other GPUs (block size and mode thresholds were chosen on the RTX
5090), elements beyond H/C/O at this size, 3D-compact systems, other basis sets than def2-SVP/MINIX with d shells.

**FP32 integrals in the early SCF phase, Q-G5 (Oct 20, 2026; AI-generated, machine-tested, RTX 5090 only).**
`-qm.gpu_eri_precision fp64|fp32|mixed|auto` (default auto = mixed on a device with FP64 at most 1/16 of FP32, from 150
basis functions; else fp64) and `-qm.gpu_fp32_switch` (1e-3). `cuda/qm_gpu_direct_fp32.cu` instantiates the MD kernels
with `float`: R, Boys (FP32 copy of the same grid/Taylor table, `__expf`), Hermite tables F, density, contraction and
digest sums in FP32; the product-centre difference P - Q is formed in FP64 before the cast (long chains: coordinates up
to ~300 Bohr); Schwarz/density screening unchanged (FP64), atomics into J/K FP64. `DirectFockCadence` (host and device
loop) leaves the FP32 phase when ||dP|| < `gpu_fp32_switch` or after 3 FP32 iterations without a new smallest ||dP||,
with a full FP64 rebuild, drops the DIIS entries of the FP32 builds, and accepts convergence only with the FP32 phase
over; the final J/K build is FP64. Cold starts only. FP32 error per class (`CURCUMA_QM_GPU_FP32_CHECK=1`, C20H42, 8
FP32 builds of the SCF): max|J_FP32 - J_FP64| / max|J| 3e-8 to 1.7e-7, for K 5e-8 to 4.6e-7 ((ss|dp)), absolute at most
4.6e-6 ((ss|ss) J); no class stood out. The FP32 Boys Taylor offset i*dT - T uses a split dT (the plain i*0.05f - T is
off by up to ~2e-6 near T = 36); this did not change the measured maxima. FP32 vs FP64 rate per class (C20, profile): 1.4-1.6x for the low-L classes ((ss|ss),
(ss|pp), (ss|ds), (ss|dp), (ss|dd)), 2-5x for the rest ((sd|pd) 5.1x, (pp|pd) 4.4x, (pp|pp) 3.4x): the kernels are not
FP64-bound alone (screening, FP64 atomics, table loads), so the 64x hardware ratio does not show.

| system (n) | fp64: SCF loop, iter | mixed: SCF loop, iter | FP32 builds / all | first FP32 builds fp64 -> fp32 | E_elec mixed - fp64 |
|---|---|---|---|---|---|
| benzene (114) | 149 ms, 10 | 173 ms, 11 (auto: not used) | 7/12 | - | 0 (10 dp) |
| naphthalene (180) | 379 ms, 13 | 343 ms, 15 | 9/16 | - | 0 |
| anthracene (246) | 806 ms, 15 | 618 ms, 16 | 9/17 | 467 -> 143 ms (3.3x) | 0 |
| C20H42 (490) | 1.64 s, 13 | 1.26 s, 14 | 8/15 | 961 -> 369 ms (2.6x) | 0 |
| C40H82 (970) | 5.79 s, 13 | 4.54 s, 14 | 8/15 | 3.43 -> 1.58 s (2.2x) | 0 |
| C60H122 (1450) | 12.4 s, 13 | 10.1 s, 14 | 8/15 | - | 0 |
| C80H162 (1930) | 21.5 s, 13 | 17.8 s, 14 | 8/15 | 12.5 -> 6.7 s (1.9x) | -1e-10 |
| C120H242 (2890) | 51.1 s, 13 | 43.3 s, 15 (rerun after the Boys-offset fix: 41.4 s, 14) | 8/16 | 29.2 -> 15.8 s (1.85x) | -2e-10 |

Single runs (the rerun of naphthalene..C80 after the Boys-offset fix reproduced the table within 1 %), load ~1. Against
the Oct 19 baseline (all FP64, old kernels) the C120H242 SCF loop goes 216 -> 41-43 s, the wall time 305 -> 107-112 s. The gain of the FP32 phase is 10-23 % of the SCF loop: after the switch the SCF still needs a
full FP64 build and 5-6 FP64 iterations (the DIIS restart costs 1-2 iterations), and the FP32 speedup falls with chain
length because the low-L classes dominate. Measured and rejected: keeping the FP32 Fock matrices in the DIIS history
(anthracene/C20/C40: 1-2 iterations more than dropping them), later switches 3e-4..3e-5 and `fp32` (FP32 down to 1e-5):
equal or 1-3 iterations more, no gain. Benzene: one extra iteration outweighs the cheaper builds, hence the 150-function
floor of auto. Test: `qm_gpu_scf_benzene_fp32` (mixed forced, device loop = host loop = CPU direct within 1e-9 Eh).
Not tested: other GPUs (in particular whether auto picks fp64 on an HPC card), heavier elements (larger exponents:
(-2 rho)^L in FP32 is the range risk), tighter `-qm.scf_threshold`, warm starts (FP32 not used there). Not implemented:
FP32 for the incremental builds of the converging phase (their error would scale with ||dP||; not done because the
convergence rule here requires FP64 builds in that phase), warp-level reduction of the J_ab atomics.

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
| Q-G5 | mixed precision far from convergence (FP32 ERI/J/K with FP64 correction on dP), never accepting convergence on an FP32 step. **Done for RI-JK only (Oct 19 2026)**: `-qm.ri_mixed_precision` (FP32 tensor, SGEMM/SSYRK, then FP64 arithmetic plus exact batched corrections), and the RI tensor streamed from host memory (`-qm.gpu_ri_stream`); exact path still FP64, see [WP10 "RI-JK ressourcensparend"](DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md). **Exact device-direct path (Oct 20 2026)**: `-qm.gpu_eri_precision` (FP32 integrals until ||dP|| < `-qm.gpu_fp32_switch`, then FP64; default auto on FP64-poor cards), section 2 | X-AP3 (FP32 eigensolve early in the SCF) |
| Q-G6 | multi-GPU: distribute shell-quartet batches over the device pool | [MULTI_GPU.md](MULTI_GPU.md) |
| Q-G7 | **done Oct 2026** (see below): two-electron part of the analytic RHF gradient on the device, integral-direct | Stage 3 kernels reused for the gradient |
| Q-G8 | **done Oct 2026** (see below): one-electron part of the analytic RHF gradient on the device, from the Q-G7 tables | Stage 3 kernels reused for the gradient |

**Q-G1 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: `-gpu cuda` uploads the stored
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

**Q-G2 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: with `-gpu cuda` the stored
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

Selection (`-gpu cuda`): `-qm.gpu_direct auto` (default) runs every integral-direct SCF on the
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
and 24 threads: 352461552 vs 352462036). Not implemented at the time: compaction of the screened quartets (every
rectangle element is launched; replaced by sorted ragged rows on Oct 19 2026, section 2), device-side symmetrisation, the resident SCF (Q-G4, added later, below), FP32 (Q-G5); the
`auto` choice between stored and direct looks only at memory, not at speed.

**Q-G4 (Oct 2026, 🤖 AI-generated, ⚙️ machine-tested; RTX 5090 only)**: device-resident RHF SCF loop
(`cuda/qm_gpu_scf.cu`, entry `curcuma_cuda_create_qm_scf`, seam `qm_gpu_scf_backend.h` -- a separate header
so that the ~9 min kernel TUs need not recompile when it changes -- loader `gpu_plugin::createQmScf`, driver
`QMEngine::runSCFDevice` in `qm_scf.cpp`). The loop control stays on the host and is the host loop's: same
Fock cadence (full, or incremental on P - P_ref with a full rebuild every 8 builds in direct mode, full only after a stall, section 2), same
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

Selection: `-qm.gpu_grad auto` (default) and `on` both use the device whenever `-gpu cuda` is active
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

Selection: `-qm.gpu_grad1e auto` (default) and `on` use the device whenever `-gpu cuda` is active (it
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

**RI-JK scaling (Oct 2026, AI-generated, machine-tested):** the setup solve B = L^-1 T runs as threaded level-3
calls on the host (heptacene n=510: 16.6 -> 2.1 s), optionally as streamed cublasDtrsm (`-qm.gpu_ri_fit on`,
`cuda/qm_gpu_ri_fit.cu`, default off: 1.8 s on the RTX 5090, slower below n ~ 450, 15.3 vs 23.9 s at n=970); the
device K build (GEMM + SYRK) is at about half the FP64 peak of the card. GPU wall heptacene 31.9 -> 13.8 s, CPU
65.7 -> 23.0 s. Since then the 3-index tensor is stored packed (lower triangle, n(n+1)/2 columns, 8 naux n(n+1)/2
bytes; device K unpacks per row block with the `unpackPairs` kernel in `cuda/qm_gpu_ri.cu`, plugin entry
`curcuma_cuda_create_qm_ri_packed`): heptacene GPU 13.8 -> 12.4 s, CPU 23.2 -> 18.7 s, upload 444 -> 250 ms; memory
limit on the 32 GB card n <= ~1155-1195 (estimated), C40H82 (n=970, 16.9 GB) now runs on the card, 146 -> 62 s.
Tables, energies, what was not tested:
[WP10 doc, "Skalierung RI-JK" and "Gepackter 3-Index-Tensor"](DFT_ROADMAP/WP10_korrelation_ri_mpn_ccsd.md#skalierung-ri-jk-okt-2026).

## 4. Kohn-Sham DFT on the GPU (Oct 2026)

> 🤖 AI-generated, ⚙️ machine-tested on one RTX 5090 (sm_120, plugin built for sm_75 and JIT-compiled from PTX on
> the first run after a rebuild). Human production testing pending.

`-method lda|pbe|b3lyp -gpu cuda` (`-qm.gpu_xc auto|on|off`, default auto). Code: `cuda/qm_gpu_xc.cu`, seam
`qm_gpu_xc_backend.h`, plugin entry `curcuma_cuda_create_qm_xc`, engine side `QMEngine::ensureGpuXC`, `xcFor`,
`xcGradient` (`qm_scf.cpp`). Tests: `ctest -L qm_gpu_dft` (34, reference = CPU run of the same binary).

- **One source of the formulas**: the XC point kernels, the Becke/Stratmann polynomials and one AO component live in
  `qm_xc_kernels.h` as `inline` / `__host__ __device__` functions; `qm_xc.cpp`, `qm_grid.cpp` and the .cu call them.
  After the move the CPU results are bitwise unchanged (`-dump_gradient` energies and gradients of hf/lda (both VWN)/
  pbe/b3lyp (both VWN) on H2O, NH3, benzene, 18 runs).
- **Grid**: built on the host (`buildMolecularGrid`: benzene 3.7 ms, caffeine 20-28 ms, against SCF loops of 0.15 and
  1.8 s on the GPU) and uploaded once per geometry; the device object dies with the geometry (`invalidateXC`).
- **Per block of B points** (B from the memory budget, max 32768; `-qm.gpu_xc_max_memory_mb`,
  `-qm.gpu_xc_block_points`): AO values with first/second derivatives and the spherical transform (one thread per
  point and shell), X = phi P (DGEMM), rho and grad rho, the functional, W += phi^T Y (DGEMM; `cublasDsyr2k` was
  40 ms per call for n = 24 and was replaced), V = W + W^T. E_xc and N from per-point arrays summed in fixed 256-point
  chunks. No floating-point atomics in the XC part.
- **SCF**: the device-resident loop takes V_xc from the device object (`QMGpuSCF::setKohnSham`, F = H + J + V_xc -
  (a_x/2) K, E with -1/2 Tr(P V_xc) + E_xc); RI-J/RI-JK keep the host loop (as HF) with device V_xc. Host and device
  loop share the driver (`runSCF`/`runSCFDevice`), so DIIS, cadence and convergence test are the same.
- **Gradient**: XC basis part (AO deriv 2 for a GGA), point part of the grid response and Becke weight derivative
  (one thread block per point, O(natoms^2) scratch per block, natoms <= 700) on the device; the Q-G7 2e kernel takes the
  exchange weight kx (`compute(P, thr, kx, g)`; kx = 1 the RHF branch, kx = 0 no exchange products and no exchange
  screening term). Plugin entries renamed with the changed vtables: `_create_qm_scf_ks`, `_create_qm_grad2e_kx`.
- **Fallbacks** (counted as `qm_xc`, `-gpu_strict true` exits with code 3): no plugin, shells above f, more than 32
  primitives, not even 32 points fit the memory budget, any device error; V_xc and the gradient then run on the CPU.

G1, V_xc and E_xc on the same density (every SCF build, `CURCUMA_QM_GPU_XC_CHECK=1`), lda vwn5/vwn_rpa, pbe, b3lyp
vwn5/vwn_rpa, def2-SVP, sg1 grid:

| system (n, points) | max rel max|dV| | max rel |dE_xc| | max |dN| | dev / CPU (4 thr) per build |
|---|---:|---:|---:|---:|
| H2O (24, 29100) | 8.9e-16 | 9.4e-16 | 8.9e-15 | 0.2-0.3 / 2.4-4.0 ms |
| NH3 (29, 38800) | 1.1e-15 | 1.5e-15 | 1.4e-14 | 0.3-0.4 / 3.6-6.5 ms |
| benzene (114, 116400) | 4.4e-15 | 3.3e-15 | 1.1e-13 | 6.6-7.7 / 49-69 ms |

The deviation is rounding order (GEMM vs Eigen sums, device exp/cbrt/log within an ulp of glibc); the points above
the density threshold are the same set. Caffeine b3lyp (n = 246, 232800 points): 50 vs 191 ms per V_xc build (8 CPU
threads); the FP64 GEMMs dominate on this FP64-poor card (benzene: X = phi P 7.4 ms, W 4.5 ms, AO 1.5 ms per build).

G2, SCF energies and times (def2-SVP, `-qm.scf_threshold 1e-8`, CPU 8 threads `-gpu none`, GPU `-gpu cuda -qm.gpu_scf
on`; direct: CPU `-qm.scf_direct on`, GPU `-qm.gpu_direct on`):

| system | method | variant | |dE| GPU-CPU | iterations CPU/GPU | SCF loop CPU / GPU | wall CPU / GPU | GPU loop |
|---|---|---|---:|---:|---:|---:|---|
| benzene | pbe | exact | 0 | 12/12 | 718 / 152 ms | 1.06 / 0.73 s | device |
| benzene | b3lyp | exact | 0 | 11/11 | 925 / 140 ms | 1.29 / 0.61 s | device |
| benzene | pbe | RI-J | 0 | 12/12 | 442 / 182 ms | 0.78 / 0.60 s | host |
| benzene | b3lyp | RI-JK | 0 | 11/11 | 421 / 275 ms | 0.54 / 0.56 s | host |
| benzene | b3lyp | direct | 0 | 13/13 | 856 / 282 ms | 0.96 / 0.74 s | device |
| caffeine | pbe | exact | 1e-10 | 24/24 | 13.3 / 1.9 s | 14.1 / 3.3 s | device |
| caffeine | b3lyp | exact | 0 | 22/22 | 11.8 / 1.8 s | 12.6 / 3.2 s | device |
| caffeine | pbe | RI-J | 0 | 24/24 | 5.5 / 1.4 s | 6.0 / 2.8 s | host |
| caffeine | b3lyp | RI-JK | 0 | 22/22 | 5.5 / 1.8 s | 6.2 / 2.6 s | host |
| caffeine | b3lyp | direct | 0 | 22/24 | 11.8 / 2.2 s | 12.6 / 3.0 s | device |

|dE| at the 10 printed decimals. Caffeine direct: the 2 extra iterations come from the Q-G5 FP32 phase
(`-qm.gpu_eri_precision auto` switches it on from 150 functions on this card); with `-qm.gpu_eri_precision fp64` the
device loop needs 22 iterations as the CPU, same energy. The ctests compare at `scf_threshold 1e-10` (7 cases, all
|dE| <= 1e-10 Eh, iterations equal except benzene b3lyp direct 21/20).

G3, gradient (CPU 8 threads; same-density comparisons from `ctest -L qm_gpu_dft`: XC basis/point/weight parts
<= 8.6e-15 Eh/Bohr, 2e part <= 4.5e-12 Eh/Bohr):

| system | method | max |g(GPU run) - g(CPU run)| | XC gradient CPU / GPU | 2e CPU / GPU | wall CPU / GPU |
|---|---|---:|---:|---:|---:|
| benzene | pbe | 1.4e-11 Eh/Bohr | 79 / 13 ms | 2.9 / 0.17 s | 4.1 / 0.9 s |
| benzene | b3lyp | 4.5e-12 | 77 / 13 ms | 2.9 / 0.17 s | 4.4 / 0.9 s |
| caffeine | pbe | 5.8e-10 | 500 / 72 ms | 68.2 / 1.9 s | 83.7 / 5.4 s |
| caffeine | b3lyp | 5.6e-10 | 506 / 72 ms | 68.2 / 1.9 s | 82.1 / 5.2 s |

(two independent SCFs at 1e-8; the caffeine difference is the SCF convergence, the same-density parts agree to 1e-14.)
`-opt` (tight, CPU run 4 threads): H2O, formaldehyde, benzene with pbe and b3lyp reach the same minimum as the CPU
`-opt` (energy of the optimised structure identical at 10 decimals, coordinates identical at the 6 written
decimals); wall GPU / CPU: H2O 1.0/1.4 s (pbe), 1.4/2.6 s (b3lyp); formaldehyde 2.9/9.7, 3.0/10.0 s; benzene 5.8/96.8,
6.6/115.7 s.

Reproducibility: with the stored device tensor the SCF energy and the XC gradient are bitwise repeatable (benzene b3lyp,
3 runs, `-qm.gpu_grad off -qm.gpu_grad1e off`); with the device 1e/2e gradient (atomics, Q-G7/Q-G8) and the direct J/K
(atomics, Q-G3) the last bits change from run to run (gradient norm in the 13th digit, direct energy 4e-13 Eh). Not
enforced (the gfn* rules of GPU_REPRODUCIBILITY.md do not cover the native QM).

Not tested: other GPUs (only the RTX 5090), ROCm and Vulkan (no entry point; `qm_xc_kernels.h` is written to compile
with hipcc, neither built nor tested), multi-GPU, systems above caffeine (24 atoms, 246 functions), grids other than
sg1 and the forced block sizes 64/96 points, the Stratmann partition on the device, cartesian d (`-qm.cartesian_d`),
basis sets other than def2-SVP, MD. Not implemented: a J-only device 2e gradient with its own screen for LDA/PBE
(kx = 0 only drops the exchange products), FP32 phases of the XC part. (AO screening and the weight derivative
above 700 atoms: section 4.1.)

### 4.1 AO screening per block and weight derivative without atom limit (Oct 30 2026)

> 🤖 AI-generated, ⚙️ machine-tested (CPU 24 cores, RTX 5090). Human production testing pending.

`-qm.xc_ao_threshold` (default 1e-12, 0 = the code path and results of before, bitwise). `QMEngine::ensureXCGrid`
sorts the points of every atom into compact blocks of <= 128 points (`qmgrid::sortGridBlocks`, recursive median
split along the longest box axis) and lists per block the shells whose reach touches the block's bounding box
(`qmgrid::buildBlockScreen`). The reach `AOEvaluator::shellRadii` is the one formula for CPU and GPU: the largest r
at which a rigorous bound of |phi|, |grad phi| or |d2 phi| (sum_k |c_k| r^L e^(-a_k r^2) and the derivative
analogues, times the largest |Q| column sum) reaches the threshold, capped at the old cutoff a_min r^2 = 80. CPU
(`buildVxc`/`buildXCGradient` with the screen): per block phi_S, P_SS, X_SS scattered into the matrix of its chunk
(union of the chunk's functions; 64 chunks, independent of the thread count, 1 vs 4 threads bitwise). GPU: consecutive
blocks merged into groups of <= 2048 points (`-qm.gpu_xc_block_points`, at least one host block) with the union of
their shells; the AO kernel writes zeros for a shell that the point's own block does not keep (per-block mask), so
every point sees exactly the functions of the CPU and the two agree to rounding; P_SS gathered, W_SS = phi_S^T Y_S
scattered group by group on one stream (no atomics), gradient columns scattered.

Gate (a): threshold 0 against the binary before the change, `-dump_gradient` files bitwise identical: CPU 48/48 (H2O,
NH3, benzene x hf, lda vwn5/vwn_rpa, pbe, b3lyp vwn5/vwn_rpa, pbe without grid response, b3lyp RI-JK x 1 and 4
threads); GPU 12/12 (H2O, benzene, caffeine x lda, pbe, b3lyp, pbe without response; `-qm.gpu_grad off
-qm.gpu_grad1e off`, since the device 1e/2e gradient has atomics and differs run to run even for HF).

Gate (b): default threshold against 0, def2-SVP, sg1, `-qm.scf_threshold 1e-10` (alkanes 1e-7), complete runs:

| system (atoms, functions) | method | |dE| Eh | max |dg| Eh/Bohr | functions kept per block |
|---|---|---:|---:|---:|
| H2O (3, 24) | pbe / b3lyp | 0 / 4.3e-14 | 1.7e-13 / 2.0e-14 | 99.5 % |
| benzene (12, 114) | pbe / b3lyp | 2.3e-13 / 1.7e-13 | 4.5e-14 / 4.3e-14 | 89.6 % |
| caffeine (24, 246) | pbe RI-J (CPU) / pbe (GPU) | 2.3e-12 / 1.4e-12 | 5.2e-13 / - | 67.6 % |
| trisaccharide C18H32O16 (66, 636) | pbe RI-J | 4.5e-12 | 1.4e-12 | 41.5 % |
| C40H82 (122, 970) | pbe RI-J CPU / GPU | 6.4e-12 / 2.7e-12 | 4.5e-13 / 3.6e-13 | 20.2 % |

The differences are at the level of the SCF convergence noise. Even 1e-8 gives |dE| <= 3.3e-11 Eh and
max |dg| <= 5.1e-11 Eh/Bohr (caffeine, trisaccharide), 1e-10 keeps 3-4 % fewer functions than 1e-12; the default
1e-12 keeps four orders of magnitude of margin to that, for basis sets that were not tested (diffuse functions).

Gate (c), XC time per V_xc build and per XC gradient (basis + point part; the Becke weight derivative separately),
pbe RI-J, def2-SVP, sg1, CPU 8 threads, RTX 5090, linear alkanes generated in the scratch directory (not in the
repository):

| alkane (atoms, functions, points) | kept | V_xc CPU off / on | V_xc GPU off / on | XC grad CPU off / on | XC grad GPU off / on | weight deriv CPU / GPU |
|---|---:|---:|---:|---:|---:|---:|
| C10H22 (32, 250, 310400) | 60.4 % | 278 / 170 ms | 66 / 66 ms | 586 / 282 ms | 82 / 101 ms | 0.21 / 0.03 s |
| C20H42 (62, 490, 601400) | 36.7 % | 2018 / 412 ms | 411 / 200 ms | 4478 / 763 ms | 464 / 305 ms | 2.0 / 0.19 s |
| C40H82 (122, 970, 1183400) | 20.2 % | 23196 / 919 ms | 2981 / 489 ms | 27950 / 1765 ms | 3044 / 746 ms | 16.2 / 1.4 s |
| C80H162 (242, 1930, 2347400) | 10.6 % | - / 2031 ms | 21391 / 1122 ms | - / 3967 ms | 21881 / 1707 ms | 130 / 11.6 s |

With screening the V_xc build grows by 2.4x, 2.2x, 2.2x (CPU) and 3.0x, 2.4x, 2.3x (GPU) per doubling of the chain:
close to linear above C20; without it 7-12x (CPU) and 6-7x (GPU). CPU C80 unscreened not run (estimated > 2 min per
build). GPU stages with screening (C40, `CURCUMA_QM_GPU_XC_PROFILE=1`): the two DGEMMs X = phi P and W = phi^T Y are
84 % of the time on this FP64-poor card; groups of 1024/2048/4096/8192 points: 529/489/544/605 ms per build, hence
2048. The device groups carry 33 % of the functions (union over 16 blocks) where the 128-point CPU blocks keep 20 %
(the GEMMs run over the union). The table was taken before the per-block zero mask was added to the AO kernel;
re-measured with it: C40 485 ms / 730 ms, C80 1110 ms / 1671 ms (V_xc / XC gradient), the same within noise.
The Becke weight derivative (O(npts natoms^2), not screened) now dominates every XC gradient from ~60 atoms on.

Weight derivative without atom limit (`weightGradKernel`, `QMGpuXC::weightGradient`, entry
`curcuma_cuda_create_qm_xc_screen`, since WP11 stage 3 `_create_qm_xc_mgga`): the old kernel held 8 natoms doubles in shared memory (48 KB, hence 700 atoms) and
4 natoms^2 in global scratch per block. Now everything per point is in global scratch, mu is recomputed (the same
expression as the stored value), S and dS are one triangle of h and h' (2 natoms^2 + 8 natoms doubles per block), and a
fixed number of 4 virtual blocks per SM is launched in waves that fit 50 % of the free memory, so the result does not
depend on the free memory. At 1-700 atoms the result is bitwise that of the old kernel (gate (a) GPU above includes
it). The Becke partition is not local (no cutoff without changing the result), so no neighbour list.
`ctest -L qm_gpu_dft` (`test_qm_gpu_wgrad`, synthetic f on the grid, rel = max|g_dev - g_cpu| / max|g_cpu|):

| structure (atoms) | grid (points, f != 0) | rel | CPU (threads) | GPU | device scratch |
|---|---|---:|---:|---:|---:|
| caffeine (24) | sg0 (93840, 80434) | 5.0e-15 | 0.07 s (4) | 0.007 s | 7.6 MB |
| peo201-chain (1410), library | 4x6 (33840, 29006) | 7.0e-15 | 178 s (8) | 7.9 s | 15.8 GB (2 waves) |
| peo201-chain (1410), stratmann | 4x6 (33840, 29006) | 1.9e-15 | 176 s (8) | 65 s * | 0.8 GB (29 waves) * |
| C250H502 (752), scratch | sg0 (2940320, 2520274) | 2.5e-13 | 2966 s (24) | 134 s | 6.2 GB (1 wave) |

\* run while another job held most of the device memory: 29 waves instead of 1-2, hence slower; the result does not
depend on the waves. The ctests use every 16th point for peo201 (2115 points) to keep the CPU reference short. C250:
grid 85 s (host), host peak 1.1 GB, device peak ~7.4 GB (sampled with nvidia-smi). Memory per wave block
8 (2 natoms^2 + 8 natoms) bytes: 9 MB at 752, 32 MB at 1410 atoms; time O(npts natoms^2).

C250H502 complete KS gradient (def2-SVP, 6010 functions, sg0, pbe, `-gpu cuda`, `-qm.scf_threshold 1e-6`): not
finished within the 45 min limit of the attempt. Measured: grid 98 s, AO screening keeps 3.6 % of the functions per
block, device-direct SCF iterations 8-33 s each (Fock build dominates, up to 1.3e10 quartets), iteration 14 after
18 min; what ran in the remaining 27 min (tight SCF continuation, device 2e gradient, XC gradient) is not known
(output buffered when the run was stopped). The weight derivative alone takes 134 s on the device there (table
above), so the XC part is not the limit; the 2e gradient and the SCF at 6010 functions were not timed.

ctest after the change: `-L "qm_1e|qm_2e|qm_grad|qm_hf3c|qm_grid|qm_lda|qm_pbe|qm_b3lyp|qm_dft_grad|qm_ri|qm_mp2"` (CPU)
green except the two numpy-dependent `qm_opt_hf3c_*`; `-L qm_gpu` 103/104: `qm_gpu_dft_scf_benzene_pbe_rij` fails
on the iteration count (GPU 14, CPU 16, tolerance 1; energies equal at the 10 printed decimals). Cause: the CPU
RI-J SCF of benzene sits at the dP = 1e-10 border (threshold 0: dP 9.9e-11 at iteration 14); the screened V_xc
(1e-13-level change) moves it to 1.8e-10 and two more iterations. Not fixed; the test tolerance was not changed.

Thread independence: the screened XC part keeps fixed chunks (benzene b3lyp exact and the `qm_dft_grad_threads_*`
cases bitwise 1 vs 4 threads). Whole runs with RI-J (caffeine) or the CPU integral-direct J/K (caffeine b3lyp, also
with threshold 0) differ between 1 and 4 threads in the last bits (|dE| 1.6e-12 Eh, |dg| 4e-13 Eh/Angstrom); that
comes from those J/K paths, not from the screening, and was there before.

Not tested: basis sets other than def2-SVP (diffuse functions: the bound stays rigorous, the gain smaller), grids
other than sg0/sg1, cartesian d, the Stratmann partition with screening, compact 3D molecules above 66 atoms (only
chains were timed), a transform that mixes shells (screening then switches itself off with a warning; not
exercised), other GPUs than the RTX 5090, ROCm. Not implemented: screening of the Becke weight derivative (it is
O(npts natoms^2) on CPU and GPU and dominates every XC gradient from ~60 atoms on), point screening by weight.


### 4.2 Meta-GGA TPSS and SCAN on the device (Oct 30 2026, WP11 stage 3)

> 🤖 AI-generated, ⚙️ machine-tested (RTX 5090 only). Human production testing pending.

`-method tpss|scan -gpu cuda`: V_xc/E_xc, the device-resident KS SCF loop and the XC gradient (basis part with the
tau term, point part, Becke weight derivative) run in `cuda/qm_gpu_xc.cu`; the former counted `qm_xc` fallback for
meta-GGAs is gone, `-gpu_strict true` runs through. Setup flag `QMGpuXCSetup::mgga` (appended), plugin entry renamed to
`curcuma_cuda_create_qm_xc_mgga`. The point functional is `evaluatePointMGGA` of `qm_xc_kernels.h` (dual numbers, the
CPU source). Per block: one strided-batched DGEMM [phi, d_x phi, d_y phi, d_z phi] P, tau per point, four DGEMMs for
W = phi^T Y + sum_k d_k phi^T diag(w v_tau / 4) d_k phi (V = W + W^T), gradient term (w v_tau / 2) sum_j d_k d_j phi o
(d_j phi P); AO screening groups as for LDA/GGA, no FP atomics, nothing depends on the free memory.

Same density vs the CPU (H2O, NH3, benzene, caffeine; sg1 and 150x974): max rel |dV| <= 1.7e-15, rel |dE_xc| <=
7.4e-15, XC gradient parts <= 4.1e-13 Eh/Bohr; SCF energies equal at 10 decimals (caffeine TPSS 1e-10). Per V_xc build,
CPU 8 threads / GPU: benzene sg1 67 / 34 ms, 150x974 896 / 468 ms, caffeine sg1 310 / 180 ms, 150x974 3982 / 2120 ms;
the four times larger FP64 GEMM work dominates on this card (benzene 150x974: X 152, W 237 of 490 ms per build, PBE
49 and 65 ms). Wall benzene tpss `-opt` CPU 35.1 s / GPU 5.8 s (same minimum). LDA/PBE/B3LYP: CPU dumps bitwise
unchanged (51/51 against the 7842f654 binary), GPU 15/18 bitwise and the three RI-JK cases within the run-to-run noise
of the reference. Tables, the rejected batched-W variant and what was not tested:
[WP11, "Stufe 3"](DFT_ROADMAP/WP11_meta_gga.md#stufe-3-meta-gga-auf-cuda-30-okt-2026).

### 4.3 Double hybrids on the device (Oct 31 2026, WP12 stage 3)

> 🤖 AI-generated, ⚙️ machine-tested (RTX 5090 only). Human production testing pending.

`-method b2plyp|b2gp-plyp|dsd-blyp|dsd-pbep86|revdsd-pbep86` (+ `-d3bj` names) `-gpu cuda`. **Energies**: P86 (PZ81 and
VWN5 local part) runs through the GGA path, B95 through the meta-GGA path of `cuda/qm_gpu_xc.cu`; both reach the device
via `evaluatePoint`/`evaluatePointMGGA` of `qm_xc_kernels.h` (no new .cu code, the former CPU-only gate in
`QMEngine::ensureGpuXC` is gone), so V_xc/E_xc, the device-resident KS SCF loop and the RI-MP2 sum run without a counted
fallback and `-gpu_strict true` runs through. **Gradient** (`qm_dh_grad_gpu.cpp`, each part with its own counted CPU
fallback):

| part | device | how |
|---|---|---|
| KS gradient (1e, 2e with a_x, XC basis/points/Becke weights) | yes | existing Q-G7/Q-G8 and XC gradient paths (P86/B95 now included) |
| F_xc[D] in every Z-vector iteration | yes | `QMGpuXC::fxc`: two densities per block, `qmxc::fxc::evaluatePointFxc` (jets of `qm_xc_fxc_kernels.h`, host and device), then the V_xc assembly |
| P2-directional XC gradient | yes | `QMGpuXC::responseGradient`: basis and point parts, Becke weight part from its energy density |
| PT2 1e derivative integrals | yes | Q-G8 with (P2, W) |
| PT2 2e, separable (P2 x P, a_x exchange) | yes | polarisation of the Q-G7 kernel: 1/2 [G(P + P2) - G(P - P2)], two passes, no new kernel |
| PT2 2e, canonical Gamma | yes | MP2 2e kernel (`cuda/qm_gpu_mp2_grad.cu`) with P = P2 = 0 |
| RI-PT2 3-index part | yes | `RIJK::threeIndexGradientDevice` (`cuda/qm_gpu_ri_grad.cu`); the 2-index metric part stays on the host as in the RI-JK gradient |
| PT2 intermediates (RI/canonical amplitudes, unrelaxed density, Lagrangian, Z-vector algebra) | no | host code, counted `qm_dh_grad_pt2`; `-gpu_strict true` with a gradient exits 3 |

New device seam: `QMGpuXCSetup::fxc` (appended; set only for double hybrids, so the block sizes of every other
functional are unchanged) and the virtuals `QMGpuXC::fxc`/`responseGradient`; plugin entry renamed to
`curcuma_cuda_create_qm_xc_fxc`. New kernels `densityPairKernel`, `respRowKernel`, `respColKernel`; no FP atomics in the
XC part, block sizes and waves independent of the free memory.

GPU vs CPU, same binary, same grid (`check_qm_gpu_dft.py`; "same density" = `CURCUMA_QM_GPU_XC_CHECK=1` /
`CURCUMA_QM_GPU_GRAD_CHECK=1` / `CURCUMA_QM_GPU_RI_GRAD_CHECK=1`, every device build also on the CPU; SCF `-qm.scf_threshold
1e-10`, 4 CPU threads; gradients in Eh/Bohr):

| molecule | grid | functionals | max rel\|dV\| | rel\|dE_xc\| | \|dE_KS\| / \|dE_total\| | iterations dev/CPU |
|---|---|---|---:|---:|---:|---|
| H2O | sg1 / 150x974 | all 10 / 4 | 1.0e-15 / 1.5e-15 | 6.6e-16 / 4.5e-15 | 0 / 0 | equal |
| NH3 | sg1 / 150x974 | all 10 / 4 | 7.3e-16 / 1.6e-15 | 7.7e-16 / 3.2e-15 | 0 / 0 | equal |
| benzene | sg1 / 150x974 | all 10 / 4 | 8.6e-16 / 2.1e-15 | 7.4e-16 / 6.4e-15 | 0 / <= 1.0e-12 | equal / +-1 |
| caffeine | sg1 | 3 | 1.0e-15 | 8.6e-16 | 0 / <= 1.5e-11 | 32/29, 26/29, 27/29 |

V_xc columns: dsd-pbep86 (VWN5 and PZ81 local part) and dsd-pbeb95-d3bj (b2plyp for caffeine); SCF: all 10 table
functionals on sg1 (b2plyp, dsd-pbep86-d3bj, dsd-pbeb95-d3bj, revdsd-pbep86 on 150x974), device-resident loop with
`-gpu_strict true`. Caffeine at 1e-10: 2-3 iterations apart at equal energies, the slow-convergence pattern of section
4.2 (TPSS caffeine 30/32); not a ctest case.

| molecule | grid | F_xc rel\|dV\| | XC response grad | PT2 2e | RI 3-index | KS XC grad (basis/points/weights) | KS 2e | \|g(GPU run) - g(CPU run)\| |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| H2O | sg1 / 150x974 | 7.5e-16 / 1.4e-15 | 1.0e-17 / 6.1e-17 | 3.5e-15 / 2.7e-15 | 1.5e-16 / 3.5e-17 | <= 3.5e-14 / <= 8.2e-14 | 8.4e-15 / 3.6e-15 | 6.5e-14 / 1.2e-13 |
| NH3 | sg1 / 150x974 | 1.1e-15 / 2.1e-15 | 1.6e-17 / 4.9e-17 | 3.1e-15 / 4.0e-15 | 8.3e-17 / 4.9e-17 | <= 3.6e-14 / <= 5.3e-14 | 6.2e-15 / 5.3e-15 | 3.7e-14 / 6.6e-14 |
| benzene | sg1 / 150x974 | 1.2e-15 / 3.8e-15 | 3.3e-17 / 1.5e-16 | 1.2e-12 / 1.2e-12 | 4.0e-16 / 3.4e-16 | <= 4.5e-14 / <= 1.1e-13 | 8.4e-12 / 7.5e-12 | 8.4e-12 / 8.1e-12 |
| caffeine | sg1 | 2.8e-15 | 9.5e-17 | 7.9e-12 | 2.6e-15 | <= 3.9e-14 | 2.3e-11 | 1.6e-10 (SCF end) |

Gradient functionals: b2plyp, dsd-pbep86-d3bj, dsd-pbeb95-d3bj, revdsd-pbep86 on sg1, dsd-pbep86-d3bj and dsd-pbeb95-d3bj
on 150x974, caffeine dsd-pbep86-d3bj and dsd-pbeb95-d3bj (RI-PT2, frozen core, grid response). The 1e-12 of benzene in the 2e columns is the quartet screen and the
atomics of the Q-G7 kernel (section 4, G3: 2e part <= 4.5e-12 on the same density).

Times (CPU 8 threads `-gpu none` vs GPU `-gpu cuda -qm.gpu_scf on`, def2-SVP, sg1, `-qm.scf_threshold 1e-8`, frozen core
and P86 on VWN5 (defaults), three runs each, median of the wall time; SCF loop and V_xc per build = median of the
printed values, V_xc per build = total of all builds / builds; quiet machine, 24 cores, RTX 5090):

| system | functional | \|dE\| Eh | iterations CPU/GPU | SCF loop CPU / GPU | V_xc per build CPU / GPU | wall energy CPU / GPU |
|---|---|---:|---:|---:|---:|---:|
| benzene | dsd-pbep86-d3bj | 0 | 12/12 | 1091 / 197 ms | 34.7 / 12.0 ms | 1.53 / 0.76 s |
| benzene | dsd-pbeb95-d3bj | 0 | 12/12 | 1469 / 455 ms | 65.9 / 33.4 ms | 1.93 / 1.04 s |
| caffeine | dsd-pbep86-d3bj | 1.8e-11 | 23/23 | 11.30 / 2.05 s | 141 / 56 ms | 12.65 / 3.98 s |
| caffeine | dsd-pbeb95-d3bj | 4.7e-11 | 22/22 | 14.82 / 4.66 s | 319 / 178 ms | 16.40 / 6.71 s |

| system | functional | wall energy + gradient CPU / GPU | Z-vector GPU (iterations) | PT2 2e + RI GPU | Z-vector CPU (iterations) | 2e pass CPU |
|---|---|---:|---:|---:|---:|---:|
| benzene | dsd-pbep86-d3bj | 5.83 / 1.74 s | 197 ms (9) | 393 ms | 1038 ms (9) | 2857 ms |
| benzene | b2plyp | 5.80 / 1.71 s | 189 ms | 391 ms | - | - |
| benzene | dsd-pbeb95-d3bj | 6.97 / 2.39 s | 517 ms (9) | 393 ms | 1525 ms (9) | 3053 ms |
| caffeine | dsd-pbep86-d3bj | 121.0 / 14.6 s | 2262 ms (18) | 4116 ms | 33458 ms (100, stalled) | 68883 ms |
| caffeine | b2plyp | 120.1 / 14.4 s | 2121 ms | 4104 ms | - | - |
| caffeine | dsd-pbeb95-d3bj | 161.0 / 21.1 s | 5680 ms (18) | 4116 ms | 71169 ms (100, stalled) | 68744 ms |

GPU parts from `CURCUMA_QM_SCF_TRACE=1` of the timed runs; CPU parts from one `test_qm_dh_grad` run each with the same
settings (`--fc --p86 vwn5 --threads 8`). Not the same scope in the last two columns: the CPU 2e pass contains the KS
and the PT2 4-index part (one pass, stage 2) plus the RI 3-index part; the GPU column is the PT2 part (two Q-G7 passes) plus
RI, the KS 2e part runs earlier on the device inside the KS gradient. **CPU Z-vector of caffeine:** stops at the
100-iteration limit with max residual 3.4e-10 (dsd-pbep86-d3bj) and 7.6e-9 (dsd-pbeb95-d3bj) instead of 1e-10 (accepted
below 1e-6 with a warning); on the GPU the same equations converge in 18 iterations. Not investigated; the likely cause is
the screening noise of the host integral-direct J/K (246 functions, the stored tensor would need 29 GB), the effect the
MP2 gradient avoids by rebuilding the tensor. Part of the caffeine speed-up (8.3x, 7.6x) comes from that stall; benzene:
3.4x, 2.9x. `-opt` benzene (sg1, default criteria): dsd-pbep86-d3bj CPU 69.5 s / GPU 15.7 s, b2plyp 80.5 / 18.0 s, final
energies of the optimised structures (`-sp` on the CPU) equal at 12 decimals (-231.549022022659, -231.818464545902 Eh).

Not tested: other GPUs than the RTX 5090, ROCm/Vulkan (no entry point), molecules above caffeine, `-md`, canonical PT2
gradient above benzene on the device (n^4 Gamma upload), Stratmann partition, cartesian d, `-qm.xc_ao_threshold 0` with the
new kernels, multi-GPU. Not implemented: the PT2 intermediates on the device, tau from the occupied orbitals (as for the
meta-GGAs), FP32 phases.

### 4.4 Open shell, UHF and UKS on the device (Nov 2026, stage 6a)

> 🤖 AI-generated, ⚙️ machine-tested (RTX 5090 only, release build sm_75 PTX run through the JIT). Human production
> testing pending.

`-method hf|hf-3c|lda|pbe|b3lyp -spin <n> -gpu cuda` (UHF/UKS of [stage 5b](NATIVE_QM_IMPLEMENTATION.md#stage-5b-open-shell-uhf-and-uks-nov-2-2026)).
Until stage 5b an open shell ran every part on the CPU (counted `qm_uhf`); `QMEngine::gpuRequested()` no longer excludes
it and the counted key is gone.

| part | device | how |
|---|---|---|
| J[Pa + Pb], K[Pa], K[Pb], stored tensor (Q-G1/Q-G2) | yes | `QMGpuJK::buildSpin`, `jkSpinKernel` (`cuda/qm_gpu_jk.cu`): every tensor element read once for the three contractions |
| same, integral-direct (Q-G3) | yes | `QmGpuDirectJkCuda::buildSpin`: each screened quartet computed once and digested into J (Pa + Pb) and both K (`DevTables::P2/K2`); screen with max(\|Pa\|, \|Pb\|) per shell block (`pmaxSpinKernel`), the expression of the new CPU `qmint::DirectJK::buildSpin` (UHF on the CPU uses it too) |
| RI-J / RI-JK | yes | two `buildJK` calls as in stage 5b, contractions on the device (`-qm.gpu_ri`); no two-density RI kernel |
| V_xc,a, V_xc,b, E_xc | yes | `QMGpuXC::vxcSpin` (`densitySpinKernel`); point kernels ONCE in `qm_xc_spin_kernels.h` (host and device; moved unchanged from `qm_xc_spin.cpp`, CPU UKS bitwise as before); no AO screening (as the CPU UKS) |
| 1e gradient (Q-G8) | yes | existing kernel with P = Pa + Pb and W = W_a + W_b |
| 2e gradient (Q-G7) | yes | `QMGpuGrad2e::computeUHF`: D = Pt Pt - kx/2 sum_s (Ps Ps + Ps Ps), the expression and order of `qmint::gradientTwoElectronUHF`; screen bound Pt_AB Pt_CD + kx (S_AC S_BD + S_AD S_BC), S = max(\|Pa\|, \|Pb\|); kx = 0 (lda, pbe) goes through `compute(Pt, 0)` |
| XC gradient incl. grid response | yes | `QMGpuXC::gradientSpin`: basis part per spin, point part of both spins (`gradRowSpinKernel`), Becke weight derivative of e_x + e_c |
| SCF loop (DIIS/ADIIS, diagonalisation, densities, convergence test) | no | host loop of `runUSCF` with the device builds, counted `qm_uhf_scf` once per engine unless `-qm.gpu_scf off` (then chosen, not a fallback); `-gpu_strict true` exits 3 otherwise |
| ECP gradient, gCP/SRB/D3 of hf-3c, RI 2-index metric part | no | host, as for RHF (not counted, as for RHF) |
| f/g shells (def2-TZVP and larger) | yes (stage 6b) | until stage 6b counted `qm_eri`/`qm_grad`; section 4.5 |

Why no device-resident UHF loop (stage 2 of the task): for the medium case below the host part of the loop (eigensolve
of two Fock matrices, DIIS, densities) is 0.54 s of 2.39 s (HF) and 0.33 s of 4.81 s (B3LYP); the J/K and V_xc builds
dominate and run on the device. A resident loop would need a two-spin version of `QMScfAccelerator`'s device inner
products and of `qm_gpu_scf.cu`; not done.

New seam members (appended, entries renamed with the vtables so an older plugin is refused, not called through a short
vtable): `QMGpuJK::buildSpin` (`_create_qm_jk_spin`, `_create_qm_eri_jk_spin`, `_create_qm_direct_jk_spin`),
`QMGpuXCSetup::spin` and `QMGpuXC::vxcSpin/gradientSpin` (`_create_qm_xc_spin`), `QMGpuGrad2e::computeUHF`
(`_create_qm_grad2e_uhf`). **Trap:** the beta density and exchange pointers of the direct digest share the slots of the
Q-G2 quartet list in a union. Two appended `DevTables` fields changed the NVVM code (FMA count) of 17 unchanged Q-G2
`eriThreadKernel` instances (PTX diff of `qm_gpu_eri.cu`) and with it the closed-shell device tensor by 1.7e-13 Eh on
benzene; with the union the closed-shell energies are bitwise again.

GPU vs CPU (ctest `-L qm_gpu_uhf`, 18 tests, 9 s wall at -j4, 34 s test time; `dump_qm_energy --gpu cuda` against the CPU
run of the same binary, `-qm.scf_threshold 1e-10`; gradients in Eh/Bohr; only `qm_uhf_scf` may be counted):

| case | \|dE\| GPU - CPU | \|d<S^2>\| | max \|dg\| | GPU vs PySCF E / g |
|---|---:|---:|---:|---:|
| UHF CH2 triplet (stored / direct / RI-JK) | 1.4e-14 / 0 / 1.4e-14 | <= 4.4e-16 | <= 3.1e-15 | 2.8e-14 / 1.2e-12 |
| UHF CH3 (distorted), O2 triplet, K atom | <= 5.7e-14 | <= 8.0e-15 | 4.4e-14 | <= 4.6e-13 / 7.9e-13 |
| UKS pbe CH2, b3lyp CH3 (stored / direct), lda H2O+ | <= 2.8e-14 | <= 5.3e-15 | <= 3.6e-14 | <= 1.6e-11 / 4.6e-12 |
| UKS b3lyp NH2 RI-JK, hf-3c CH3 | <= 2.8e-14 | <= 2.2e-15 | 2.6e-14 | 1.5e-11 / - |

FD gate of the device path (4-point central differences of `-gpu cuda` energies, 0.005 Bohr, 3 coordinates each): UHF
CH2 3.8e-11, UKS pbe CH2 2.4e-11, UKS b3lyp CH3 1.7e-10 Eh/Bohr.

Medium case: caffeine radical cation (24 atoms, charge +1, `-spin 1`, def2-SVP, 246 functions, sg1 for B3LYP,
`-qm.scf_threshold 1e-8`, CPU 8 threads `-gpu none` (integral-direct, the host does not hold the 29 GB tensor) vs
`-gpu cuda` (the stored tensor fits the 32 GB card; `-qm.scf_direct on` for the device direct rows); one run each, wall
time of the process; RTX 5090, 24-core host):

| method | path | iterations | J/K(+V_xc) builds | wall energy | wall energy + gradient (2e / 1e / XC part) |
|---|---|---:|---:|---:|---:|
| UHF | CPU direct | 60 | 31.3 s | 32.7 s | 102.1 s (68.4 / 0.58 / - s) |
| UHF | GPU stored, one pass | 60 | 1.60 s | 3.6 s | 5.6 s (1.91 / 0.076 / - s) |
| UHF | GPU stored, two passes | 60 | 3.39 s | 5.4 s | - |
| UHF | GPU direct, one pass | 60 | 5.45 s | 6.4 s | - |
| UHF | GPU direct, two passes | 60 | 9.99 s | 11.1 s | - |
| UKS b3lyp | CPU direct | 36 | 32.1 s | 33.6 s | 102.6 s (68.8 / 0.58 / 0.73 s) |
| UKS b3lyp | GPU stored, one pass | 36 | 4.21 s | 6.1 s | 8.1 s (1.88 / 0.076 / 0.115 s) |
| UKS b3lyp | GPU stored, two passes | 36 | 5.26 s | 7.0 s | - |
| UKS b3lyp | GPU direct, one pass | 36 | 6.54 s | - | 9.5 s (1.88 / 0.075 / 0.115 s) |
| UKS b3lyp | GPU direct, two passes | 36 | 9.31 s | 10.3 s | - |

"Two passes" = `CURCUMA_QM_UHF_TWO_PASS=1` (diagnostic: one single-density `buildJK` per spin, the stage 5b scheme). The
one-pass build is 2.1x (stored) and 1.8x (direct) faster for the J/K part of UHF; the direct one-pass build computes
half the quartets (774.9e6 vs 1547.4e6 over the SCF). Differences on caffeine+: GPU stored vs CPU direct |dE| 9.3e-12
(UHF) / 1.9e-11 (UKS) Eh, |d<S^2>| 3.9e-12 / 5.8e-13, max |dg| 2.3e-11 / 2.9e-11 Eh/Bohr (the CPU direct screen);
GPU direct vs CPU direct UHF |dE| 0 in this run (bits schedule-dependent), UKS 1.4e-12 Eh, max |dg| 2.5e-11; one vs two
passes <= 1.1e-12 Eh. Total speed-up energy + gradient: 18x (UHF), 13x (UKS b3lyp).

Closed shell unchanged: 20 dumps (hf/b3lyp/pbe/lda/hf-3c/mp2, stored, direct, RI-JK, water/methanol/benzene, energies and
gradients, CPU and `-gpu cuda`, 4 threads) against the binary of a4a23b5e: every CPU dump md5-identical, every GPU energy
bit-identical; GPU gradients differ by <= 4.4e-14 Eh/Bohr, the size of their own run-to-run spread (atomics; 2.8e-14 on
benzene with the old binary); integral-direct energies at 4 threads differ run to run in both binaries (host
`omp critical` order, device atomics), CPU direct at 1 thread is md5-identical. `ctest -L qm_` 978/978 passed.

Not tested: other GPUs than the RTX 5090, multi-GPU, ROCm and Vulkan (no open-shell entry points, untouched), molecules
above caffeine+, `-md` with `-gpu cuda` on an open shell (`-opt` only one smoke run: b3lyp CH3 radical, sg1, converged in
13 steps, last gradient norm 1.5e-5 Eh/Bohr; no ctest), f/g bases on
the device (fallback), ECP open shells on the device (the ECP gradient stays host), `-qm.gpu_eri_precision` FP32 phases
(the open-shell builds are FP64 only), `-qm.ri_mixed_precision` and the low-memory RI modes with UHF, the Stratmann
partition. Not implemented: device-resident UHF/UKS loop, two-density RI-JK kernel, AO screening on the UKS path (CPU and
device), incremental (dP) direct builds for UHF.

### 4.5 f and g shells on the device (Nov 2026, stage 6b)

> 🤖 AI-generated, ⚙️ machine-tested (RTX 5090 only, release build sm_75 PTX through the JIT). Human production
> testing pending.

Until stage 6b `qmint::buildEriGpuTables` refused L > 2, so with def2-TZVP and larger every ERI, J/K and gradient part
ran on the CPU. Survey at 592ba35c (`-gpu cuda`, `dump_qm_energy`, H2O/NH3/CH4 TZVP/QZVPP, caffeine TZVP, CH3 UHF TZVP,
CH3I TZVP): counted fallbacks `qm_eri` (stored tensor, every SCF), `qm_direct` (`-qm.scf_direct on`), `qm_grad` and
`qm_grad1e` (every gradient), `qm_mp2_grad` (mp2 gradient), `qm_ri_grad` (RI-JK/RI-J gradient and the RI-PT2 gradient of
the double hybrids); RI-J/RI-JK contractions and XC already ran on the device (the RI tensor comes from the CPU, the
device XC takes up to g). After stage 6b none of them is counted on these cases; what remains counted is not f/g
specific: `qm_uhf_scf` (host UHF loop, stage 6a) and `qm_dh_grad_pt2` (host PT2 intermediates, WP12).

| part | how (files) |
|---|---|
| tables | `buildEriGpuTables` takes L <= 4; pair classes `LA*cb + LB` with cb = 3 when every shell is s/p/d (the layout of before, unchanged) and cb = 5 with f/g shells; the plugin derives cb from `sh_L` (`qmgpu::pairClassBase`), so `QMGpuEriTables` and the entry names are unchanged |
| ERI, stored tensor (Q-G2) and integral-direct J/K (Q-G3) | every class with an f or g shell: `cuda/qm_gpu_eri_generic.cu`, one runtime-L McMurchie-Davidson kernel in Hermite-matrix form (F tables of the active component pairs per primitive pair, contraction and spherical transform folded in, as `mdInit`; R^0 level by level in shared memory, degree-ordered simplex up to (gg\|gg)); one warp per quartet for the small classes, one block of 128 threads for the large ones; the pair with the higher L is put on the bra side (the ket Hermite step per primitive quartet stays small); ket sums and the quartet block in a per-group slot of global memory, so all f/g launches go to one stream; same `quartetOf` screen and primitive break as the s/p/d kernels; FP64 also in an FP32 phase (`-qm.gpu_eri_precision`: only the s/p/d classes run FP32) |
| s/p/d classes in an f/g basis | unchanged kernels; `mdInit` skips f/g pairs and now accepts d transform entries <= 1e-14 outside the 5d pattern (the S-metric orthonormalisation leaves up to 9.6e-17 in def2-QZVPP, 1.7e-30 in caffeine def2-TZVP; exactly zero in def2-SVP); without that the integral-direct build fell back to the old hermite kernels (H2O QZVPP direct 8.7 s instead of 4.9 s) |
| 2e gradient (Q-G7, RHF/KS/UHF) and MP2 2e gradient | the warp kernels were already runtime-L; changed: `sF` layout `sfSize(SFB)` (Boys + powers for L <= 17), f/g classes with much work per quartet run `grad2eBlockKernel` (one block per quartet, scratch in global memory); `qm_gpu_mp2_grad.cu` got the same edits (it is a copy) |
| 2e gradient screen | `qmint::gpuGradScreen`: eri_screening / 100 for s/p/d (unchanged), / 1e4 with f/g shells: at / 100 CH4/def2-QZVP was 1.0e-10 Eh/Bohr from the unscreened CPU gradient (same density), at / 1e4 4.5e-13, for 0.8 % more quartets |
| 1e gradient (Q-G8) | template on the Hermite order: `MAXL = 5` (s/p/d, code of before), `MAXL = 9` for f/g bases |
| RI 3-index gradient | second instance of `riGrad3cKernel` (Ltot <= 13, orbital g raised + auxiliary g, one warp per block for shared memory); its own Boys table (the s/p/d one bit for bit as before); auxiliary h still refused |
| ECP matrix and gradient | stay on the CPU (decision below) |

Checks (`dump_qm_energy` GPU vs CPU of the same binary, 8 threads, `-qm.scf_threshold 1e-10`; gradients Eh/Bohr):

| case | \|dE\| | max \|dg\| |
|---|---:|---:|
| HF H2O, NH3, CH4, HF with def2-TZVP, TZVPP, QZVP, QZVPP (16) | <= 1.9e-13 | <= 2.1e-12 (CH4 QZVP; the others <= 5.6e-13) |
| pbe/b3lyp H2O TZVP, b3lyp H2O QZVP, pbe NH3 QZVPP | <= 1.4e-13 | <= 1.8e-13 |
| integral-direct HF H2O TZVP/QZVP, b3lyp NH3 QZVP | <= 7.1e-14 | <= 4.6e-13 |
| RI-JK HF H2O TZVP/QZVP, RI-J pbe H2O QZVP | <= 2.3e-13 | <= 1.1e-13 |
| mp2 gradient H2O TZVP/QZVP, ri-mp2 H2O QZVPP, b2plyp gradient H2O TZVP | <= 1.9e-13 | <= 2.3e-13 |
| UHF CH3, NH2, CH2 triplet TZVP, direct UHF CH2 QZVP, UKS b3lyp CH3 TZVP | <= 2.8e-14 | <= 5.0e-14 |
| ECP + f: HF CH3I, HI, HgCl2 TZVP, b3lyp HI TZVP | <= 9.1e-13 | <= 7.9e-12 (HgCl2, \|g\| up to 24 in the parts) |
| caffeine def2-TZVP HF / b3lyp (sg1), `-qm.scf_threshold 1e-8`, CPU direct vs GPU direct | 2.3e-13 / 2.5e-12 | 2.2e-10 / 2.2e-10 (SCF threshold) |

Same density on both sides (`CURCUMA_QM_GPU_ERI_CHECK`, `_JK_CHECK`, `_GRAD_CHECK`, `_RI_GRAD_CHECK`): device tensor vs CPU
`buildERI` element-wise <= 2.0e-14 (H2O QZVPP, CH4 QZVP, HgCl2 TZVP, H2O QZVP cartesian); direct J/K vs CPU `DirectJK`
<= 2.1e-14 with the same quartet count; 2e gradient <= 6.8e-14 (H2O QZVP/QZVPP, cartesian, UHF CH2 TZVP), 4.4e-13 (CH4
QZVP), 3.8e-12 (HgCl2, |g| up to 24); 1e gradient <= 7.9e-12 (HgCl2, |g| up to 47), else <= 6e-15; MP2 2e 7.3e-15; RI
3-index 3.9e-14. The device tensor has exact zeros where the CPU has rounding-level values and vice versa
(`zero_pattern_mismatch` 29932 of 1.9e8 elements for H2O QZVPP spherical, 0 in the cartesian basis); max difference
1.7e-14. The 2e gradient numbers are with the f/g screen of `gpuGradScreen`; caffeine GPU 2e gradient 28 s in a run
with the md kernels off and 42 s in the run of the table (same kernel code, not investigated; the CPU was busy with
another 8-thread job in both). FD gate on the device path (4-point, 0.005 Bohr, ctest): HF H2O TZVP 2.5e-10, b3lyp H2O QZVP (g shells)
1.8e-10, UHF CH2 TZVP 3.7e-11 Eh/Bohr; against PySCF: HF H2O QZVPP 5.4e-13 Eh / 5.6e-13 Eh/Bohr, b3lyp H2O QZVP 2.6e-12 Eh
/ 9.2e-12 Eh/Bohr.

Closed shell s/p/d unchanged: the 20 dumps of stage 6a (hf/b3lyp/pbe/lda/hf-3c/mp2, stored, direct, RI-JK, water/
methanol/benzene, CPU and `-gpu cuda`, 4 threads) against the binary of 592ba35c (same build options): every CPU dump
md5-identical except integral-direct at 4 threads (differs run to run in both binaries, 2.3e-13 Eh), every stored-path
GPU energy bit-identical, GPU gradients within their own run-to-run spread (<= 1.4e-14 vs 1.4e-14 old vs old). Open
shell (8 dumps of stage 6a): identical energies except O2/TZVP, which now runs its ERI on the device (5.7e-14 Eh). PTX
is NOT a usable criterion here: compiling the unchanged 592ba35c `qm_gpu_eri.cu` twice gave 10 of 63 kernels with
different PTX (3 with a different FP instruction sequence); the measured energies are the criterion.

Timing (RTX 5090, 24-core host, wall time of the process incl. start-up):

| case | CPU 1 thread | CPU 8 threads | GPU |
|---|---:|---:|---:|
| H2O def2-QZVPP (117 functions) HF energy, stored tensor | 3.85 s | 1.80 s | 1.28 s |
| same, energy + gradient | 18.6 s | 3.92 s | 2.42 s (2e gradient 1.10 s) |
| same, `-qm.scf_direct on` on the device | - | - | 4.9 s (generic f/g classes ~75 % of the direct kernels) |
| caffeine def2-TZVP (494 functions) HF, integral-direct, SCF / 2e gradient / process | - | 129 s / 1147 s / 1287 s | 39 s / 42 s / 85 s |
| caffeine def2-TZVP b3lyp (sg1), same | - | 157 s / 1132 s / 1302 s | 60 s / 34 s / 99 s |

ECP on the device (task step 4, not done): CPU ECP matrix 98 ms (HgCl2 TZVP), 85 ms (CH3I), 38 ms (HI) on 1 thread; ECP
gradient (1e part on the device minus its 10 ms) ~110 ms (HgCl2) and ~80 ms (CH3I) on 1 thread. The GPU runs of these
molecules take 2.0-2.9 s in total with the 2e gradient 0.27-0.99 s on the device, so a device ECP would save at most ~9 %
(HgCl2, 1 thread) and less with threads; not worth a device port of the projection code. Kept on the CPU.

CPU f/g path (task step 5, not done): H2O def2-QZVPP HF energy 3.85 s and energy + gradient 18.6 s on 1 thread (generic
`mdQuartet` for l_a + l_b > 4 and the general gradient path); no contained, bitwise-verifiable change was attempted.

Not tested: other GPUs than the RTX 5090, multi-GPU, ROCm and Vulkan (no f/g kernels there, untouched), h orbital shells
(the parser rejects them), auxiliary h in the RI gradient (refused, counted `qm_ri_grad`), FP32 phases of f/g classes
(they run FP64), molecules larger than caffeine, `-opt`/`-md` with f/g bases on the GPU, the low-memory RI modes with
f/g. Compile time (nvcc -ptx per TU, eight TUs in parallel, before / after): `qm_gpu_eri.cu` 832 / 825 s,
`qm_gpu_direct.cu` 836 / 841 s, `qm_gpu_direct_md.cu` 271 / 280 s, `qm_gpu_direct_fp32.cu` 260 / 261 s, the gradient
TUs below 2.1 s, new `qm_gpu_eri_generic.cu` 0.6 s: the generic kernel is runtime-L (no unrolled templates), the plugin
build time is unchanged within noise.
