# Native GFN1/GFN2 SCF Performance

> ⚠️ AI-generated, machine-tested (ctest). Not human production-tested.
> Numbers measured on one machine (AMD Ryzen 9 9950X3D, Zen5; MKL serial;
> GCC `-march=native -O3`); indicative, not guaranteed.

Single authoritative performance record for the native `curcuma::xtb::XTB`
GFN1/GFN2 path (the canonical `gfn1`/`gfn2` backend). Supersedes the scattered
per-WP perf notes. Methodology, measured breakdown, what worked, what did not,
and the residual gap to the Fortran references.

> **Single-core** record below. For **intra-molecule multi-threading** (`-threads N`
> on one large molecule) see [SQM_THREADING.md](SQM_THREADING.md). Key cross-finding:
> the "MKL threading the eigensolve — no effect" note here is a build artifact — this
> build links `libmkl_sequential`, so MKL BLAS/LAPACK never threads; the multi-thread
> wins (setup/gradient/Fock 3–5×) come from the CxxThreadPool, not MKL.

## Methodology — reproduce with `scripts/sqm_bench.sh`

```bash
scripts/sqm_bench.sh [N_REPEATS=3] [caffeine triose complex]
```

- Single core: `taskset -c 0`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`.
- Workload: energy **+ gradient** (`curcuma -sp` always computes the gradient;
  tblite/xtb asked for `--grad`). min-of-N wall time.
- References: `tblite` (Fortran + MKL) and `xtb` (Fortran, static MKL). **xtb is
  run `--norestart`** — otherwise it reads back `xtbrestart` and "converges" in
  ~1 iteration, a warm-start artifact that under-reports its cold cost ~2×.
- curcuma is profiled with the `-verbosity 2` timing breakdown (setup / SCF /
  post-SCF / gradient) and a `-verbosity 3` per-SCF-iteration breakdown
  (potential / Fock / eigensolve sub-phases / populations).

## Headline result (complex, 231 atoms, single core, energy+gradient)

| method | curcuma (orig) | curcuma (now) | tblite | xtb (cold) |
|--------|---------------:|--------------:|-------:|-----------:|
| gfn1   | 2982 ms        | **1221 ms**   | 2562   | 1367       |
| gfn2   | 2944 ms        | **1364 ms**   | 1567   | 979        |

- **gfn1: −59%. Now beats both tblite and xtb.**
- **gfn2: −54%. Now beats tblite; 1.39× xtb (was ~3×).**
- Same on the smaller systems: triose gfn1 67 ms (tblite 74, xtb 97);
  triose gfn2 88 ms; caffeine in the 15–30 ms range (setup/one-time-init bound).

All energies bit-identical to the pre-optimization values and to the tblite
reference within tolerance — full native-GFN suite 45/45 green
(`sqm_val_*`, `ecomp_*`, `d4_diag_*`, `d4_dedq`, `xtb_gradient_*`, `xtb_cpscf`,
`gfn{1,2}_align`).

## What was found (deep timing analysis)

The original bottleneck story ("eigensolver-bound, near-optimal") was incomplete.
Per-iteration the eigensolve was on par with tblite; the real costs were
elsewhere. Measured contributors and the fixes:

1. **Initial guess (biggest single win).** Default was bare-H0 (`scf_guess=h0`):
   the SCF spent ~6 early iterations recovering from a guess far from the
   solution. The already-implemented **EEQ guess** (single-shot dftd4 EEQ shell
   charges) starts in the right basin. Made default. complex SCF iterations
   gfn1 **35→16**, gfn2 **34→22** — energy bit-identical (the SCF fixed point is
   guess-independent). `scf_guess=h0` still available.

2. **Eigensolve transform overhead.** `solveEigen` reduced the generalized
   problem `F C = S C ε` by forming a **dense** `S^{-1/2}` and doing
   `X·F·X` + a dense back-transform = ~6·nao³ flops/iter on top of `dsyevd`.
   Replaced with the textbook reduction: cache the **Cholesky factor L** of the
   constant overlap (`buildOrthonormalizer`, Eigen LLT), per iteration reduce
   with LAPACK **`dsygst`** and back-transform with one **triangular solve**
   (`dtrsm` via Eigen) — ~2·nao³, and a cheaper one-time setup (LLT vs a full
   eigendecomposition of S). complex transform+back: gfn1 393→190 ms, gfn2
   304→160 ms. (Trap fixed: the cached factor must be **column-major**
   `Eigen::MatrixXd` for the Fortran `dsygst`; the project `Matrix` is row-major.)

3. **Density build.** `P = C·diag(occ)·Cᵀ` used all nao columns. With Fermi
   smearing only the occupied + a few fractional columns carry weight; restrict
   the GEMM to those (`leftCols(ncol)`, occ > 1e-12). complex density GEMM
   gfn1 129→58 ms, gfn2 101→56 ms (energy-neutral, dropped columns weigh <1e-12).

4. **GFN2 D4 charge-response (biggest gfn2 win).** The post-SCF D4 gradient ran
   an unconditional **Mulliken CPSCF / Z-vector** charge response —
   `computeMullikenChargeResponse`, **574 ms = 88% of the gfn2 D4 cost**, even
   though the documented default `d4_charge_source="eeq"` specifies the cheap
   analytic single-shot-EEQ response (`D4ChargeModel::addChargeResponseGradient`,
   an adjoint solve on the small EEQ linear system, ~ms). The routing was never
   wired. Wired it: `eeq` (default) uses the analytic response, `mulliken` keeps
   the CPSCF. The **energy weighting stays on the SCF Mulliken charges** either
   way, so the energy is unchanged; only the gradient path differs.
   complex gfn2 post-SCF **653→83 ms**, gradient still within the <5e-5 FD target
   (`xtb_gradient_*` green).

5. **Over-tight convergence default.** Binding criterion was `max|dq_shell| <
   1e-6`; the total energy reaches <1e-8 Eh several iterations earlier. Loosened
   the default to **`1e-5`** (energy bit-identical across the whole molecule set;
   complex gfn2 22→19 it, triose 13→11, caffeine 15→12). Still tighter than xtb's
   production default. MD/opt can tighten via `-scf_threshold`.

### Per-iteration breakdown after the fixes (complex, `-verbosity 3`)

gfn1 (nao 680): solve-eigen dominates (`dsyevd` ~37 ms/it; reduce ~6.5,
back ~5.4, density ~3.6); potential/Fock/populations < 0.3 ms/it.
gfn2 (nao 558): `dsyevd` ~22 ms/it; reduce ~3.9, back ~3.4, density ~2.5;
potential build ~10 ms/it (multipole + D4 dE/dq), Fock ~2 ms/it.

## What did NOT help / not pursued

- **Forcing MKL instruction set on AMD** (`MKL_ENABLE_INSTRUCTIONS=AVX2/AVX512`):
  no effect. The legacy `MKL_DEBUG_CPU_TYPE` trick is gone in modern MKL. xtb,
  tblite and curcuma all link the **same MKL**, so the eigensolve floor is shared
  — the BLAS library is not the differentiator.
- **Partial diagonalization** (`dsyevr`/`dsygvx`, lowest nocc+buffer vectors):
  not pursued. The gradient needs only occupied orbitals, so it is feasible for
  the default path, but at 44–54 % occupation the saving on `dsyevd` is modest,
  and a truncated `m_wfn.C` breaks the `mulliken` CPSCF path and orbital
  analysis. Earlier measurements (and this occupancy) suggest a net loss; left as
  documented headroom.

## Residual gap to xtb (gfn2 complex, 1364 vs 979 ms)

xtb runs the same MKL, so its per-iteration eigensolve is comparable; its
advantage is (a) **fewer SCF iterations** (15 vs 19 — partly because xtb's
default convergence is looser than curcuma's even after the 1e-5 change), and
(b) it folds the gradient into the SCC step whereas curcuma's setup (179 ms,
incl. 107 ms multipole-integral build) and gradient (224 ms) are separate. gfn1
already beats xtb because its larger basis makes the eigensolve dominate, where
curcuma's leaner setup/gradient wins. Closing the gfn2 gap further would need
either a looser convergence (`-scf_threshold 5e-5` → 16 it, energy-identical) or
setup/gradient work — both lower-value than the wins above.

## Verification

```bash
cd release && make -j16 curcuma
ctest -R "sqm_val_|ecomp_|d4_diag|d4_dedq|xtb_gradient|xtb_cpscf|gfn1_align|gfn2_align"
# 45/45 pass; energies bit-stable (e.g. gfn2 H2O = -5.07036982 Eh,
# gfn2 complex = -329.52707823 Eh, gfn1 complex = -343.17980352 Eh)
scripts/sqm_bench.sh 3
```

Pre-existing, unrelated ctest failures (not touched here): `AAAbGal_dtemplate`,
`gfnff_numgrad_fixed_charges`, `cli_gfnff_01/02`, `cli_sqm_11` (ipea1 needs
TBLite, not built).
