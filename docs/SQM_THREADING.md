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
| Eigensolve (dsygst/dsyevd/dtrsm) | `xtb_scf.cpp` `solveEigen` | MKL threads (`MklThreadScope`), not the pool |

`buildFock`–`solveEigen` are the *only* regions threaded inside the SCF loop (others
are too small to amortise a per-iteration dispatch). The pool is persistent (workers
reused across iterations and geometry steps).

## Performance (complex, 231 atoms, nao=558, single-point E+grad)

Per-phase, gfn2, min-representative (`-verbosity 3`):

| Phase | t1 | t8 | speedup |
|---|---|---|---|
| overlap + H0 | 88 ms | 17 ms | 5.2× |
| multipole setup | 116 ms | 30 ms | 3.8× |
| **setup total** | **206 ms** | **49 ms** | **4.2×** |
| build Fock | 92 ms | 28 ms | 3.3× |
| solve eigen | 605 ms | 607 ms | **1.0×** |
| gradient | 218 ms | 59 ms | 3.7× |
| **TOTAL** | **1593 ms** | **1072 ms** | **1.5×** |

Total min-of-5 scaling (unpinned, 32-core host):

| method | t1 | t2 | t4 | t8 |
|---|---|---|---|---|
| gfn1 | 1192 | 1080 | 976 | 922 (1.29×) |
| gfn2 | 1509 | 1244 | 1122 | 1062 (1.42×) |

**Single-core unchanged** (taskset -c 0, MKL=1, default): gfn1 1200 ms (was 1221),
gfn2 1362 ms (was 1364).

## The eigensolve wall — and why (re-validated, root cause found)

The embarrassingly-parallel parts scale 3–5×, but they are ~30 % of the run, so the
overall speedup is ~1.3–1.4×, plateauing by t4. The bottleneck is the per-iteration
eigensolve, which does **not** thread:

| nao | dsyevd t1 | dsyevd t8 |
|---|---|---|
| 558 (231 atoms) | 417 ms/it | 419 ms/it |
| 4028 (1410 atoms) | 9442 ms/it | 9378 ms/it |

It is unchanged at *both* sizes, and unchanged even with `MKL_NUM_THREADS=8` forced.
**Root cause: this build links `libmkl_sequential` (the CMake MKL discovery pulls in
both `mkl_sequential` and `mkl_gnu_thread`; the sequential layer loads first and
wins).** So MKL's BLAS/LAPACK is single-threaded regardless of any thread setting, and
`MklThreadScope` around `solveEigen` is currently a harmless no-op. This *explains* the
earlier single-core note "MKL threading the eigensolve: no effect at 231 atoms" — it
was a build-config artifact, not a `dsyevd` limitation. The CxxThreadPool wins (setup,
gradient, Fock) are independent of MKL and scale at *both* sizes (polymer setup
10065→2469 ms, 4.1×).

**Two ways to unlock the eigensolve (follow-ups, not in this change):**
1. Link threaded MKL only (drop `mkl_sequential`, keep `mkl_gnu_thread` + libgomp).
   Then `MklThreadScope` takes effect immediately — but it must stay gated by the
   intra-thread budget so molecule-level batches are not oversubscribed.
2. A custom, GPU-ready divide-and-conquer symmetric solver (deferred track) replacing
   only the `dsyevd` core and keeping the dsygst/dtrsm reduction — avoids the MKL
   threading-layer dependency entirely. See the implementation plan.

## Tested / not tested

- **Tested:** energy and gradient bit-identical between t1 and t8 on `complex` (gfn1
  and gfn2); native-xTB ctest suite green at the default (serial) and unchanged with
  threading; small-molecule gradient tests run via the serial size-guard path.
- **Not tested:** correctness of the auto-gate under a real oversubscribing batch
  (verified by construction/guards, not by a stress benchmark); systems > 1410 atoms;
  thread counts beyond the host core count; Windows/macOS thread behaviour.
- **Not implemented:** threading of the D4 in-SCF potential (`addDispersionPotential`,
  a separate dispersion subsystem) and of the eigensolve beyond MKL.
