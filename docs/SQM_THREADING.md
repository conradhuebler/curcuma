# Native GFN1/GFN2 SCF — intra-molecule multi-threading

> ⚠️ 🤖 AI-generated / ⚙️ machine-tested. Energies and gradients are bit-identical
> between serial and threaded runs on the validation set; not yet human production
> tested. See the caveats at the end.

## What this is

The native `curcuma::xtb::XTB` SCF (GFN1/GFN2) can use multiple cores **within a
single calculation**. Before, the whole SCF ran serial (MKL pinned to one thread) so
that Curcuma's coarse, *molecule-level* parallelism — one independent SCF per
`CxxThreadPool` worker in MD/conformer/batch runs — was never oversubscribed.

Intra-molecule threading helps the case the molecule-level model cannot: a **single
large molecule** (`-sp`/`-opt` of e.g. the 231-atom `complex`). It uses the project
`CxxThreadPool` (the GFN-FF idiom), not OpenMP.

## How to use it

Pass `-threads N`. A single-molecule `-sp`/`-opt`/MD run then fans the SCF out over up
to `N` cores:

```bash
curcuma -sp complex.xyz -method gfn2 -threads 8
```

Default (no `-threads`) is serial — identical to before. There is nothing else to set.

## Auto-gate (why it is safe)

Intra-molecule threading is honoured **only** when it will not oversubscribe:

- **Molecule-level batches stay coarse.** The concurrent batch workers
  (`CurcumaOpt::SPThread`/`OptThread`, the numerical-Hessian `HessianThread`) raise
  `curcuma::intraParallelSuppressed()` (RAII, `src/core/intra_parallel_context.h`) for
  the duration of their task; the SCF inside them runs serial. So N molecules × N
  intra-threads never happens. Single `-sp`/`-opt` and MD run on the main thread where
  the flag is false → full thread budget. Conformer search optimises sequentially on
  the main thread, so each conformer's SCF is free to thread.
- **Small work stays serial.** A size guard (`effectiveIntraThreads`,
  `kMinWorkPerThread`) keeps tiny systems serial — the pool dispatch would cost more
  than the work. Only molecules with enough shells/atoms/AOs parallelise.

The gate lives in tracked curcuma code, not the vendored `CxxThreadPool` header.

## What is parallelised

| Region | Where | Mechanism |
|---|---|---|
| Overlap + H0 integrals | `xtb_h0.cpp` `getHamiltonianH0` | `CxxThreadPool`, stripe over shells (disjoint rows) |
| GFN2 multipole integrals | `xtb_multipole.cpp` `setupMultipole` | stripe over AOs / source atoms |
| Gradient (dominant: overlap & multipole integral derivatives) | `xtb_gradient.cpp` | stripe over atoms, thread-local grad/dEdcn + reduce |
| Fock build | `xtb_scf.cpp` `buildFock` | stripe over AO rows (disjoint) |
| Eigensolve (dsygst/dsyevd/dtrsm) | `xtb_scf.cpp` `solveEigen` | MKL threads (`MklThreadScope`), not the pool — **capped at 8**, see below |

`buildFock`–`solveEigen` are the *only* regions threaded inside the SCF loop (others
are too small to amortise a per-iteration dispatch). The pool is persistent (workers
reused across iterations and geometry steps).

## Eigensolve thread cap (2026-07)

The per-iteration MKL eigensolve is capped at `min(intra_threads, 8)`. The D&C
`dsyevd` is memory-bandwidth-bound and **regresses** past ~8 threads (complex/231:
215 ms @8 vs 290 ms @16 on a 16-core Ryzen), while the hand-threaded
Fock/gradient/setup regions still use the full `-threads N`. Override with
`CURCUMA_EIG_MAX_THREADS`. Runs at ≤8 threads are byte-unchanged; at
`-threads 16` vs gxtb cold-start the ratio went gfn1 1.18→0.86×, gfn2 1.75→1.54×.

## Performance (complex, 231 atoms, nao=558, single-point E+grad)

> ⚠️ **The per-phase table below is from 2026-06 and its t1 column is now stale.**
> The 2026-07 shell-pair-blocked integral kernels cut the single-core setup from
> 209 to 91 ms and the gfn1 gradient from 132 to 73 ms
> ([docs/SQM_PERFORMANCE.md](SQM_PERFORMANCE.md#integral-setup-2026-07)), so the
> t1 numbers here — and therefore the speedup ratios, which had more serial work
> to amortise — no longer apply. **Re-measure before quoting.** The t8 numbers
> and the qualitative picture (eigensolve-bound) still hold.

Per-phase, gfn2, min-representative (`-verbosity 3`), **as measured 2026-06**:

| Phase | t1 (2026-06) | t8 | speedup |
|---|---|---|---|
| overlap + H0 | 88 ms | 17 ms | 5.2× |
| multipole setup | 116 ms | 30 ms | 3.8× |
| **setup total** | **206 ms** | **49 ms** | **4.2×** |
| build Fock | 92 ms | 28 ms | 3.3× |
| solve eigen | 614 ms | 215 ms | **2.85×** (after WP1) |
| gradient | 218 ms | 59 ms | 3.7× |
| **TOTAL** | **1550 ms** | **679 ms** | **2.3×** (after WP1) |

Current single-core reference for the same buckets: overlap+H0 28 ms, multipole
setup 56 ms, setup total 91 ms, gradient 175 ms (gfn2) / 73 ms (gfn1).

Total min-of-5 scaling, complex/231 (after WP1+WP2+WP2b — threaded MKL + D4 cache/threading):

| method | t1 | t2 | t4 | t8 |
|---|---|---|---|---|
| gfn1 | 1149 | 738 | 478 | **372 (3.09×)** |
| gfn2 | 1371 | 896 | 631 | **513 (2.67×)** |

(After WP1 alone it was gfn1 375 / gfn2 650 at t8; WP2/2b cut the D4 in-SCF cost.)

**Single-core unchanged** (taskset -c 0, default, no `-threads`): gfn1 1160 ms,
gfn2 1336 ms (≈ the pre-WP1 1200/1362; slightly faster from the `W`-as-gemm).

## The eigensolve — was the wall, now threads (WP1, fixed)

Originally the per-iteration eigensolve did **not** thread (dsyevd unchanged t1→t8 at
nao=558 *and* 4028, even with `MKL_NUM_THREADS=8`). Root cause: the build linked
**`libmkl_sequential`** alongside the threaded layer (CMake's `find_package(BLAS)`
resolved to the sequential MKL; it loaded first and won), so MKL's BLAS/LAPACK was
single-threaded regardless of any setting — `MklThreadScope` was a no-op. The old
single-core note "MKL threading the eigensolve: no effect at 231 atoms" was a
build-config artifact, not a `dsyevd` limitation.

**WP1 removed the sequential MKL** (guard `find_package(BLAS/LAPACK)` linkage with
`if(NOT USE_MKL)` so only the explicit threaded `mkl_gnu_thread` links; `main.cpp` sets
`MKL_Set_Num_Threads(1)` as the global default so existing BLAS users and
molecule-level batches are unchanged, and only `MklThreadScope(eig_threads)` around
`solveEigen` bumps it). The eigensolve now threads:

| gfn2 complex/231 | t1 | t8 | speedup |
|---|---|---|---|
| dsyevd | 425 ms | 133 ms | 3.2× |
| back-transform | 63 ms | 11 ms | 5.6× |
| solve eigen (total) | 614 ms | 215 ms | 2.85× |

This is why the totals above scale 2.3–3.1× rather than the pre-WP1 1.3–1.4×. gfn1
scales best (3.07×); gfn2 is gated by the serial D4 in-SCF potential (~180 ms),
addressed by WP2. Further headroom (custom CPU divide-and-conquer eigensolver,
purification, GPU) is tracked in [SQM_THREADING_WP.md](SQM_THREADING_WP.md).

## Tested / not tested

- **Tested:** energy and gradient bit-identical between t1 and t8 on `complex` (gfn1
  and gfn2); native-xTB ctest suite green at the default (serial) and unchanged with
  threading; small-molecule gradient tests run via the serial size-guard path.
- **Not tested:** correctness of the auto-gate under a real oversubscribing batch
  (verified by construction/guards, not by a stress benchmark); systems > 1410 atoms;
  thread counts beyond the host core count; Windows/macOS thread behaviour.
- **Not implemented:** threading of the D4 in-SCF potential (`addDispersionPotential`,
  a separate dispersion subsystem) and of the eigensolve beyond MKL.
