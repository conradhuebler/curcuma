# Native GFN1/GFN2 SCF — single-molecule performance work packages

> ⚠️ 🤖 AI-generated roadmap. Each WP is energy-neutral unless stated; validate
> bit-identity (or documented tolerance) against the serial baseline.

Follows the intra-molecule threading already landed (see [SQM_THREADING.md](SQM_THREADING.md)).
Baseline at that point — **complex/231, gfn2, t8 = 1072 ms**, where the SCF is dominated by:

| region | t8 | note |
|---|---|---|
| eigensolve | 605 ms (57%) | dsyevd 418 + BLAS-3 (reduce 76 / back 62 / density 47 = 185) |
| D4 in-SCF potential | ~180 ms (17%) | serial; full D4 gradient recomputed every iter |
| setup / gradient / Fock / pop | ~160 ms | already threaded (CxxThreadPool) |

The eigensolve does not thread because the build links **sequential** MKL. The WPs
attack the eigensolve and the D4 potential, cheapest first.

---

## WP1 — Unlock threaded MKL (+ free BLAS-3 wins)  ✅ DONE (machine-tested)
**Result:** sequential MKL removed (`ldd`/`readelf` show only `mkl_gnu_thread`); the
eigensolve now threads — **dsyevd 425→133 ms (3.2×)**, solve-eigen 614→215 (2.85×).
complex/231 TOTAL **t8: gfn1 375 ms (3.07×), gfn2 650 ms (2.26×)** vs pre-WP1 1.29/1.42×.
Single-core unchanged/slightly faster (gfn1 1160, gfn2 1336 ms). Energy bit-identical;
native-xTB + GFN-FF ctests green (only the pre-existing `gfnff_numgrad_fixed_charges`).

**Goal:** make the linked MKL actually thread, so the eigensolve's BLAS-3 parts
(reduce/back-transform/density, ~185 ms) and `dsyevd` use cores; gated so molecule-level
batches stay serial.

**Root cause:** `curcuma_core` links the threaded layer (`mkl_gnu_thread`, CMakeLists
L976) *and* the sequential layer pulled in by the unconditional `find_package(BLAS/LAPACK)`
(L933–942, also L962–972); sequential loads first and wins (confirmed by `ldd` + the
runtime `MKL_THREADING_LAYER=GNU` no-op).

**Approach:**
1. Guard the `find_package(BLAS/LAPACK)` link block (L933–942) and the `USE_BLAS` block
   (L962–972) with `if(NOT USE_MKL)` so only the explicit threaded MKL (L976) is linked.
2. Keep MKL **serial by default** globally (`mkl_set_num_threads(1)` once at startup in
   `main.cpp`, under `USE_MKL`) so existing GFN-FF/UFF/optimizer BLAS calls are unchanged
   and molecule-level batches are not oversubscribed. Only the existing
   `MklThreadScope(eig_threads)` around `solveEigen` bumps it (already gated).
3. **Free win:** convert the gradient's energy-weighted density `W` from a rank-1 loop
   (`xtb_gradient.cpp:~92`) to one `gemm` (`C_occ·diag(2ε)·C_occᵀ`) → BLAS-3, threads.

**Files:** `CMakeLists.txt`, `src/main.cpp`, `xtb_gradient.cpp`.
**Risk:** medium (global link change touches all BLAS users — must verify no oversubscribe
in batch, no numerical change). **Exit:** `ldd` shows no `mkl_sequential`; eigensolve
BLAS-3 parts scale with `-threads`; energies bit-identical; native-xTB ctests green;
batch (conformer/Hessian) wall unchanged.

## WP2 — Trim the D4 in-SCF potential  ⏳
**Goal:** stop recomputing the full D4 Cartesian gradient every SCF iteration.
**Approach:** `addDispersionPotential` (`xtb_native.cpp`) calls
`computeEnergyAndGradient(with_gradient=true)` only to read `dEdq`, discarding the
gradient. Add a `dEdq`-only path to `D4Evaluator` (no per-pair Cartesian accumulation).
Optionally then parallelize the D4 evaluator's O(nat²) loops (helps GFN-FF too).
**Files:** `dispersion/d4_evaluator.*`, `xtb_native.cpp`. **Risk:** low (algorithmic,
energy-neutral). **Exit:** potential-build per-iter drops; energy+gradient bit-identical.

## WP3 — Seeded / iterative eigensolver  ⏳
**Goal:** exploit that near convergence the occupied subspace barely rotates — refine the
previous eigenvectors instead of a full `dsyevd` each iteration.
**Approach:** subspace iteration / block-Davidson seeded from the previous `C` (and across
geometry steps in opt/MD). Keep `ε` (needed for Fermi smearing + the gradient `W`). Fall
back to a full solve on the first iteration / when the residual is large.
**Files:** `xtb_scf.cpp` (`solveEigen`), `xtb_native.h`. **Risk:** medium (convergence /
robustness; must not change the converged fixed point). **Exit:** fewer full solves; energy
bit-identical at convergence; SCF iteration count unchanged.

## WP4 — Custom divide-and-conquer eigensolver (CPU first, GPU-ready)  ⏳
**Goal:** replace MKL `dsyevd` with our own **threaded CPU** divide-and-conquer symmetric
eigensolver, designed so the same kernels port to GPU later. Removes the dependence on the
MKL threading layer and is the long-term path to a scalable, portable eigensolve.
**Approach:** keep the `dsygst`/`dtrsm` reduction + back-transform (WP1 threads them);
replace only the tridiagonal-eigensolve core with a blocked Cuppen divide-and-conquer
(reduce to tridiagonal, recursively split, rank-1 update / secular equation, deflation),
parallelised over subproblems via `CxxThreadPool` and written with tile/batch structure
suitable for a future cuSOLVER/custom-CUDA backend. Assess
`ParallelEigenSolver.hpp` (already in `qm_methods/`) as a seed or remove it if dead.
Gate behind a config switch so MKL stays default until parity is proven.
**Files:** new `qm_methods/dc_eigensolver.*`, `xtb_scf.cpp`. **Risk:** high (correctness of
a from-scratch eigensolver; eigenvalue accuracy feeds energy/gradient). **Exit:** matches
MKL eigenpairs to ~1e-10; energy bit-identical; competitive wall at N≥500 single molecule.

## WP5 — O(N) / strategic (very large single molecules)  ⏳
**Goal:** break the O(N³) eigensolve wall for >1000-atom single molecules.
**Approach (menu, pick per need):**
- **Divide-and-conquer DFT (spatial fragmentation):** overlapping fragments, each solved
  independently (embarrassingly parallel + GPU), density stitched → O(N). Approximate.
- **Density-matrix purification** (McWeeny / TC2): no diagonalisation, only dense matmuls
  + traces → threads/GPU perfectly; needs an `ε`-free energy-weighted density (via FPS) and
  a finite-T treatment (or T=0).
- **Mixed precision:** FP32 early SCF iterations, FP64 near convergence (~2× CPU, ~8–16× GPU).
- **GPU eigensolve** (cuSOLVER `syevd` / MAGMA) under `USE_CUDA` — the per-iter offload.
**Files:** new modules; `USE_CUDA` path. **Risk:** high / research. **Exit:** O(N) scaling
demonstrated on the 1410-atom polymer with bounded energy error.

---

**Order:** WP1 → WP2 (cheap, measurable, de-risk the eigensolve ceiling) → WP3 → WP4
(CPU D&C) → WP5. WP4/WP5 share the GPU-portable kernel design.
