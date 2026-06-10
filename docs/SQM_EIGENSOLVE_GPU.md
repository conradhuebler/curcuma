# SQM Eigensolve — MKL-free, GPU-portable kernels

> ⚠️ AI-implemented, machine-tested only. Validated bit-identical / to-tolerance vs MKL & Eigen on
> the molecules below; not human production tested. **MKL `dsyevd` stays the default**; everything
> here is opt-in.

The native GFN1/GFN2 SCF is eigensolve-bound. MKL `dsyevd` already wins on CPU at ~200 atoms, so the
goal of these levers is **GPU-readiness**: MKL-free, GEMM-only kernels that are *correct*, as the CPU
foundation for a future GPU port. CPU parity with MKL is not the bar. All three are selected via
`-eigensolver <backend>` (or `CURCUMA_EIG_TRED2` for the tridiagonalization backend).

## 1. WY-block-Q blocked tridiagonalization (`CURCUMA_EIG_TRED2=blocked`)

The native eigensolver's own blocked Householder reduction `A = Q·T·Qᵀ` (`blockedTridiag` in
`native_eigensolver.cpp`), the MKL-free alternative to Eigen's `Tridiagonalization`. A panel of `nb`
columns is reduced, the trailing block updated once by a BLAS-3 rank-2k update, and **the orthogonal
`Q` is accumulated once per panel by a compact-WY block apply** `Q ← Q·(I − V·T·Vᵀ) = Q − ((Q·V)·T)·Vᵀ`
(`T` built columnwise like LAPACK `dlarft`). This replaced a per-column BLAS-2 rank-1 update, so the
kernel is now **GEMM-only** (bar the inherent panel `symv`, as in LAPACK `dsytrd`) → GPU-portable.

- Backends: `eigen` (default, Eigen vectorized — routes BLAS-2 via MKL), `blocked` (this), `scalar` (baseline).
- Measured (complex/231, nao=558): blocked `tred2` **t8 22.2 → 10.20 ms** (Eigen 9.99, on-par), t1 16.4 ms.
- Bit-identical eigenpairs (ctest `native_eigensolver{,_blocked,_scalar}`, recon ~2e-14 to n=558) and SCF
  energy (caffeine −42.14723025, complex −329.52707823 Eh).

## 2. Density-matrix purification (`-eigensolver purify`)

Builds the 0 K idempotent density **directly from F, with no diagonalization** — trace-correcting TC2
(Niklasson 2002, `purifyDensity`): in the orthonormal basis `Ã = L⁻¹·F·L⁻ᵀ` (already formed), Gershgorin
init to [0,1], then `P ← P²` (trace too high) / `2P − P²` (too low); idempotency monitored by
`Tr(P)−Tr(P²)` (reuses the `P²` already formed) → **one GEMM per iteration**. Density `P_AO = 2·L⁻ᵀ·P̃·L⁻¹`;
the gradient's energy-weighted density is supplied GEMM-only as `W = 2·L⁻ᵀ·(P̃·Ã·P̃)·L⁻¹` (no eigenpairs).
The pure-GEMM/trace GPU bridge.

- **0 K only:** requires `-electronic_temperature 0` and a HOMO-LUMO gap. With finite-T smearing or on
  non-convergence it **warns and falls back to the eigensolver**. No HOMO/LUMO or `mulliken`-CPSCF in this mode.
- Validated: ctest `purification` (idempotency / Tr / projector / W vs Eigen ~1e-12 to n=558); SCF energy
  bit-identical to MKL@0K (caffeine, complex gfn2 −329.52707823 / gfn1 −343.17980352); H2O `-opt` energies
  AND gradient norms identical step-by-step (the `W̃` route is correct).
- Deferred: finite-T (grand-canonical) purification.

## 3. Seeded block LOBPCG (`-eigensolver lobpcg`) — EXPERIMENTAL

Iterative lowest-`k` eigensolver (`lobpcgLowest`): block LOBPCG (Knyazev) for the lowest `nocc + buffer`
eigenpairs, **subspace-recycled across SCF iterations** (the previous iteration's vectors seed the next,
`m_eig_seed`). Diagonal preconditioner; a robust generalized Rayleigh-Ritz filters the `[X|W|P]` subspace
through its Gram matrix (drops near-null directions). Guard vectors: only the lowest `nConverge` states
must reach `tol`; the extras absorb the slow gapless boundary eigenvector. GEMM/SpMM-based → GPU-friendly.

- **Net-loss on dense GFN (by design):** at ~50% occupancy `k ≈ n/2`, the ~3k Rayleigh-Ritz subspace costs
  more than one `dsyevd` (complex single-core: mkl 1.32 s vs lobpcg 3.75 s, ~2.8×). LOBPCG pays only with
  `k ≪ n` or sparse `A` — i.e. coupled with C1 (fragmentation/sparsity, not yet implemented). Falls back to
  `dsyevd` on non-convergence. Kept as the validated GPU-portable iterative kernel.
- Validated: ctest `lobpcg` (eigenvalues ~1e-13, subspace ~1e-10 vs Eigen, incl. warm-start); SCF energy
  bit-identical to MKL (caffeine, complex @300 K); H2O `-opt` energies + gradient norms identical.

## Status / scope

| Lever | Flag | Role | Validation |
|-------|------|------|------------|
| WY-block-Q | `CURCUMA_EIG_TRED2=blocked` | GEMM-only tridiagonalization | bit-identical, on-par w/ Eigen |
| Purification | `-eigensolver purify` | GEMM-only 0 K density (no diag) | bit-identical @0 K |
| LOBPCG | `-eigensolver lobpcg` | iterative, recycled (experimental) | bit-identical, net-loss on dense |

**Not implemented (deferred):** C1 sparsity/fragmentation for >1000 atoms (a further approximation);
finite-T grand-canonical purification. See memory `sqm-scf-eigensolve-is-the-wall` and
`AIChangelog.md`.
