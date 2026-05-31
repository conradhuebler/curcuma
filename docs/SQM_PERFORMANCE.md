# Native GFN1/GFN2 Performance

> ⚠️ AI-generated, machine-tested (ctest). Not human production-tested.
> Speedups measured on one Linux/MKL machine; numbers are indicative.

Optimizations to the native `curcuma::xtb::XTB` SCF path (2026-05-31). All are
energy-neutral within the 1e-8 validation tolerances — see verification below.

## What was changed

1. **Eigensolver: factor the overlap once, not every iteration**
   (`xtb_scf.cpp`, `xtb_native.{h,cpp}`)
   `solveEigen` previously built a fresh `GeneralizedSelfAdjointEigenSolver(F, S)`
   every SCF iteration, re-doing the Cholesky of the geometry-constant overlap `S`
   each time with Eigen's built-in (QR-based, non-divide-and-conquer) path.
   Now the Löwdin orthonormalizer `X = S^{-1/2}` is built **once per geometry**
   (`buildOrthonormalizer`), and each iteration solves the cheap **standard**
   problem `(X·F·X)·C~ = C~·ε` via MKL's divide-and-conquer `dsyevd`
   (declared `extern "C"`; Eigen `SelfAdjointEigenSolver` fallback when MKL is
   absent). Back-transform `C = X·C~` preserves `CᵀSC = I`, so the gradient and
   CPSCF paths are unchanged. The density is built with one weighted GEMM
   (`C·diag(occ)·Cᵀ`) instead of per-orbital rank-1 updates.

2. **MKL pinned to one thread for the SCF** (`MklSerialScope`, `xtb_native.h`)
   MKL spawns an OpenMP thread team per BLAS/LAPACK call. For the serial SCF over
   small/medium dense matrices this per-call overhead is measured *slower* than
   serial up to at least 231 atoms. The whole `Calculation()` runs MKL-serial via
   a thread-local guard (`MKL_Set_Num_Threads_Local(1)`). This is also the right
   foundation for parallel workloads: run many calculations through the project's
   CxxThreadPool and each stays serial — no MKL oversubscription.

3. **GFN2 skips the unused D4 JSON pair list** (`d4param_generator.{h,cpp}`)
   `D4ParameterGenerator::GenerateParameters` builds a `nlohmann::json` object per
   atom pair (~12 string-keyed inserts each) under OpenMP. Native GFN2 reads the
   C6 reference cache directly through `D4Evaluator` and never consumes that list,
   so it now opts out via `setBuildPairLists(false)`. GFN-FF still builds it (its
   `ForceFieldThread` consumes `d4_dispersion_pairs`) and is unaffected.

## Measured (min of 3, `release/`, MKL + AVX2)

| Molecule (atoms) | Method | Baseline SCF | Optimized SCF | Speedup |
|------------------|--------|-------------:|--------------:|--------:|
| complex (231)    | gfn1   | 4251 ms      | 2435 ms       | 1.75×   |
| complex (231)    | gfn2   | 3157 ms      | 2078 ms       | 1.52×   |
| triose (66)      | gfn1   | 42 ms        | 34 ms         | 1.24×   |
| caffeine (24)    | gfn1   | 5.9 ms       | 4.1 ms        | 1.44×   |

Large systems are eigensolver-bound, so the eigensolver change dominates there.
For small/medium GFN2 the single-point `SCF_ms` is dominated by **one-time**
MKL/OpenMP runtime init + D4 setup charged to iteration 0 (e.g. caffeine gfn2:
iter0 ≈ 25 ms, iter1+ ≈ 0.6 ms each). That one-time cost amortizes across the many
geometries of an `-opt`/MD run.

## Verification

```bash
cd release && make -j4 curcuma
# Native GFN correctness (energies, gradients, CPSCF, D4): 12/12 pass
ctest -R "xtb_gradient|xtb_cpscf|d4_dedq|d4_diag|ngfn" --output-on-failure
# GFN-FF unaffected by the D4 pair-list opt-out: 17 gfnff_val pass
ctest -R "gfnff_val" --output-on-failure
```

Energies are bit-stable (e.g. gfn2 H2O = −5.07036982 Eh, equals the recorded
reference). Pre-existing unrelated failures (confirmed by a stash-baseline rebuild):
`gfnff_numgrad_fixed_charges`, `cli_gfnff_01/02` (expect an `*.opt.xyz` from a
single point; energies pass), and `cli_sqm_11` (`ipea1` needs TBLite, not built).

## Not done / further headroom

- Fock isotropic loop + Mulliken triple-loop vectorization — marginal: steady-state
  per-iteration is already ~0.6 ms and large systems are eigensolver-bound.
- Explicit multi-threaded eigensolve — not beneficial at ≤231 atoms (MKL threading
  loses to serial there); would only help much larger systems.
- The `gfn1_method` / `gfn2_method` wrapper layer still duplicates the
  `ComputationalMethod` adapter; GFN2's wrapper additionally carries an
  error-handling/orbital-analysis layer GFN1 lacks. Consolidating them is a
  separate refactor (kept out of this performance change to keep the diff reviewable).
