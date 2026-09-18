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

| method | curcuma (orig) | curcuma (2026-06) | curcuma (2026-07) | tblite | xtb (cold) |
|--------|---------------:|------------------:|------------------:|-------:|-----------:|
| gfn1   | 2982 ms        | 1221 ms           | **791 ms**        | 2562   | 1297       |
| gfn2   | 2944 ms        | 1364 ms           | **998 ms**        | 1567   | 938        |

- **gfn1: −73%. 1.6× faster than xtb.**
- **gfn2: −66%. At parity with xtb (1.06×; was ~3×, then 1.39×).**
- Idle-machine medians of 7, same run as the cur/gx table below.
- The 2026-07 column combines the [shell-pair-blocked integral
  kernels](#integral-setup-2026-07) with the [FP32 mixed-precision
  eigensolve](#fp32-mixed-precision-eigensolve--now-on-by-default-2026-07),
  which is now on by default.
- Same on the smaller systems: triose gfn1 67 ms (tblite 74, xtb 97);
  triose gfn2 88 ms; caffeine in the 15–30 ms range (setup/one-time-init bound).

Energies are **not** bit-identical to the pre-2026-07 values any more, by
construction: the blocked kernels shift the last ulp via FMA contraction and the
mixed-precision eigensolve moves 3 of 14 reference energies in the 12th decimal.
Both are quantified in their sections below and stay 10000x inside the 1e-8
tblite gate; the native-GFN suite is green
(`sqm_val_*`, `xtb_gradient_*`, `xtb_cpscf`, `gfn{1,2}_align`, 143/143 with
large-system, numgrad, solvation and extrapolation included).

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

## Integral setup (2026-07)

Everything above targeted the SCF. The integral kernels had never been touched,
and they held a straightforward ~2.9x redundancy.

### The defect

`getHamiltonianH0()` already looped shell pairs, but called `cgto_overlap()`
**once per AO-COMPONENT pair**. The GFN2 multipole build looped AO pairs
outright, and the gradient called `cgto_overlap_grad()` /
`cgto_multipole_grad_transformed()` per component too. So for a p-p shell pair
the full `nprim_a x nprim_b` primitive loop ran nine times, and each leaf
recomputed the Gaussian product centre, `gamma`, and
`S00 = pow(pi/gamma,1.5)*exp(-ai*aj/gamma*R2)` — a libm `pow` **per primitive
pair per component**. None of that depends on the cartesian powers; only the
1-D moments do, and they need just `(la, lb)` per axis.

The stale comment above the overlap loop already claimed the integral was
"computed once and broadcast to all AO pairs". Now it actually is.

### Sizing (complex/231, GFN2: nao 558, nsh 340)

Primitive-pair evaluations, `(sum over AO of nprim)^2` vs
`(sum over shells of nprim)^2`:

| granularity | evaluations |
|---|---:|
| per AO component (before) | `(122*3 + 436*4)^2 = 2110^2 = 4.45e6` |
| per shell pair (after)    | `(122*3 + 109*8)^2 = 1238^2 = 1.53e6` |

**2.90x** — the ceiling for the primitive-setup part of each bucket.

### Measured (complex/231, single core)

| bucket | before | after |
|---|---:|---:|
| overlap + H0 | 82.1 ms | **28.0 ms** |
| multipole setup | 119.7 ms | **56.3 ms** |
| — AO dipole/quadrupole integrals | 109.4 ms | **42.7 ms** |
| — origin shift + traceless | 8.6 ms | 11.3 ms |
| setup total | 209 ms | **91 ms** |
| gradient (gfn1) | 132 ms | **73 ms** |
| gradient (gfn2) | 195 ms | 175 ms |

gfn2 total 1199 -> 1083 ms; gfn1 total 1142 -> 1017 ms.

The gfn2 gradient improves less because only its overlap half is blocked;
`cgto_multipole_grad_transformed()` is still per-component (~63 ms) — the
clearest remaining lever.

### Numerics: why this is NOT bit-identical, and why that is acceptable

The blocked kernels are **algebraically exact**. Compiled with
`-ffp-contract=off` they agree with the per-component path to **exactly
0.000e+00** (verified by a temporary per-element cross-check over all AO pairs
of triose). The deviation comes solely from GCC contracting `a*b+c` into FMA
differently in the restructured code — the build enables `-mfma` and GCC
defaults to `-ffp-contract=fast`.

Measured effect on every reference molecule: **energies bit-identical to all 12
printed digits**; gradients differ by at most **1.7e-14 Eh/Bohr** (2.7e-11
relative), against the project's 1e-8 validation gate.

Forcing `-ffp-contract=off` in these headers was **tried and rejected**: it makes
the kernels exactly reproducible, but the same code serves the gradient, whose
time then goes 190 -> 272 ms. That costs more than the optimisation gains.

Worth knowing when reading any "bit-identical" claim in this project: with
`-mfma` + `-ffp-contract=fast` the golden values are already specific to this
compiler and CPU, not a portable property.

### Design rules for anyone touching these kernels

- The **primitive contraction order is part of the numerical contract**. Keep the
  accumulation i-major/j-minor.
- Keep `pow(M_PI/gamma, 1.5)` verbatim; `t*sqrt(t)` rounds differently.
- `cgto_overlap` forms the product centre as `(ai*xa + aj*xb) / gamma`;
  `cgto_overlap_grad` uses `* invg`. These round differently and are kept
  separate **on purpose** — do not unify them.
- **Triangular S/H0 + mirror is not bit-identical** and was rejected: `S(mu,nu)`
  and `S(nu,mu)` are computed independently today, and the transposed call sums
  the same terms in transposed order. IEEE addition is not associative, so
  mirroring moves the last ulp. (The individual factors *are* swap-invariant:
  `gamma = ai+aj` and `ai*aj` are exact under commutativity.) A bit-identical
  variant exists — emit both blocks from one primitive pass — but after the
  blocking above the only remaining transcendental is `exp`, so it is worth
  ~15-20 ms while breaking the disjoint-row striping invariant. Documented
  headroom, not implemented.

### Blocking the multipole GRADIENT: tried, measured, reverted

The same treatment was applied to `cgto_multipole_grad_transformed()` (the last
per-component kernel, ~100 ms of the gfn2 gradient). It is **slower** and was
reverted. Two measurements, both on complex/231:

| variant | gfn2 gradient |
|---|---:|
| per-component (kept) | **173 ms** |
| blocked, full 2x2x3 moment table | 239 ms |
| blocked, table bounded by angular momentum | 202 ms |

Why it loses, where the overlap kernels win:

1. **s-heavy pairs dominate by count.** Building the full moment table costs 36
   `moment1d` calls per primitive pair, but an s-s pair only ever needs 9. Bounding
   the table by `ang` recovers part of that (239 -> 202 ms) but not enough.
2. **The hoisted work is cheap here.** The overlap kernels hoist a libm `pow`;
   this one hoists `moment1d`, whose transcendental is a `sqrt`. Meanwhile the
   per-component assembly is ~95 lines, so the moments are a small share of the
   kernel — there is little to amortise.
3. **The per-pair state is large.** `MpGradPair` is 44 doubles vs 13 for
   `OverlapPrimPair`, so the table must be written to memory and re-read per
   component instead of staying in registers.

Conclusion: the remaining gfn2 gradient cost is the **assembly**, not the moment
evaluation. Blocking is the wrong tool for it; a cheaper assembly (or fewer
components via symmetry) would be the lever.

### Remaining headroom

- Memoizing the primitive-pair invariants across shell pairs by shell *type*
  (all atoms of an element yield byte-equal `alpha`/`coeff`), which would remove
  the remaining `pow` entirely rather than just de-duplicating it per shell pair.
- The `origin shift` pass grew slightly (8.6 -> 11.3 ms) and is memory-bound: 9
  `nao x nao` temporaries, ~67 MB of traffic. Fusing it into the blocked kernel
  must keep reading `m_S(mu,nu)` (not the kernel's own `S`, a different
  expression tree).

## curcuma vs gxtb (cross-method, 2026-07-25)

Produced by `scripts/bench_vs_gxtb.py`, which is the only harness that measures
all three native methods against the same reference binary. It gives gxtb a
fresh temp dir and `--norestart` per run (reusing `xtbrestart` warm-starts gxtb
to ~3 iterations and understates its time), and pins curcuma with `-threads K`
plus `OMP/MKL_NUM_THREADS=K`.

```
scripts/bench_vs_gxtb.py test_cases/sqm_reference/molecules/complex.xyz \
    gfnff,gfn1,gfn2 1,16 5 1
```

complex/231, energy+gradient, median of 5, gxtb 6.7.1 cold-start.
**`cur/gx < 1.0` means curcuma is faster.**

Median of **7**, measured on an **idle** machine (load 0.24):

| method | curcuma K=1 | gxtb K=1 | cur/gx K=1 | curcuma K=16 | gxtb K=16 | cur/gx K=16 |
|---|---:|---:|---:|---:|---:|---:|
| gfnff | 36 ms | 34 ms | **1.05** | 35 ms | 29 ms | 1.19 |
| gfn1  | 791 ms | 1297 ms | **0.61** | 326 ms | 379 ms | **0.86** |
| gfn2  | 998 ms | 938 ms | **1.06** | 406 ms | 301 ms | 1.35 |

Before this round the same table read gfnff 1.26/1.59, gfn1 0.89/0.91,
gfn2 1.32/1.38.

**Single core: gfn1 is 1.6x faster than gxtb; gfn2 and gfnff are at parity
(within 6%). At 16 threads gfn1 stays ahead but gfn2 and gfnff fall behind** —
gxtb scales better, and the remaining limit is analysed in
[SQM_THREADING.md](SQM_THREADING.md#where-the-scaling-limit-actually-is-2026-07).

> ⚠️ This table is **invalidated by any change to the integral, SCF or setup
> path**. Re-measure and record the commit hash when quoting it. Measured at
> commit `c42e2feb`, single socket, MKL, `-mfma`.

> ⚠️ **Measure on an idle machine.** An earlier revision of this table
> (gfnff 1.27/1.37, gfn1 0.78/0.89, gfn2 1.20/1.36) was taken while an unrelated
> 4-5 core job was running and overstated every ratio; gfnff was worst hit
> (1.27 vs the true 1.05). Check `uptime` first.

Caveat on gfnff: at ~36 ms total the run is dominated by process start plus
one-time setup, so the ratio is sensitive to noise and to whether the
`.topo.json` cache is present.

## CPU vs GPU is a question of SIZE, not of method quality (2026-07)

Measured on an idle machine, RTX 5080, energy+gradient, against the CPU at
`-threads 16`. Energies identical in every pair.

**complex, 231 atoms (nao 558) — the GPU loses:**

| method | CPU t16 | CUDA | |
|---|---:|---:|---|
| gfn1 | 349 ms | 621 ms | 1.8x slower |
| gfn2 | 519 ms | 718 ms | 1.4x slower |
| gfnff | 38 ms (t1) | 298 ms | 7.8x slower |

**polymer, 1410 atoms (nao 3222) — gfn2 on the GPU wins decisively:**

| phase | CPU t16 | CUDA | |
|---|---:|---:|---|
| setup | 1198 ms | 1050 ms | 1.14x |
| **SCF (11 it)** | 11999 ms | **3407 ms** | **3.5x** |
| post-SCF | 706 ms | 915 ms | 0.77x (GPU slower) |
| **gradient** | 1374 ms | **363 ms** | **3.8x** |
| **TOTAL** | **15368 ms** | **5922 ms** | **2.6x faster** |

(E = -2088.25340678 on both.)

So there is a **crossover between 231 and 1410 atoms for gfn2**, and it lands
exactly where the hardware argument predicts: the GPU wins the O(nao³) eigensolve
(SCF 3.5x) and the gradient (3.8x), while setup and post-SCF — small, serial,
transfer-bound — stay flat or regress.

### For gfnff, a single point is the WRONG benchmark

A gfnff single point pays the full topology build and parameter generation for
*one* energy evaluation — the worst possible case for a device whose pipeline
targets the per-step path. Measured that way, gfnff/polymer looks like a 1.65x
GPU **loss** (820 ms vs 497 ms).

Measured as MD, which is what gfnff is actually used for, the verdict inverts.
polymer/1410, `-threads 16`, two step counts so the slope separates setup from
per-step cost:

| | 50 fs | 200 fs | **per step** | setup (intercept) |
|---|---:|---:|---:|---:|
| CPU t16 | 2475 ms | 8375 ms | **39.33 ms** | ~508 ms |
| CUDA | 1279 ms | 2103 ms | **5.49 ms** | ~1004 ms |

**7.2x faster per MD step.** Trajectory identical (Epot -203.361820,
Etot -201.406573 on both).

The GPU setup is ~2x the CPU's, which is exactly why the single point misleads:
the two lines cross at **~15 MD steps**, and everything beyond that is GPU
territory.

Practical guidance:
- **gfn1/gfn2, large systems** → `-gpu cuda` (2.6x at 1410 atoms).
- **gfnff, any MD or optimisation past ~15 steps** → `-gpu cuda` (7.2x per step).
- **Single points, and anything small** → CPU.

Two earlier revisions of this section were wrong and are corrected above: first
"the GPU is slower than the CPU" (generalised from 231 atoms only), then "gfnff
does not cross over" (generalised from a single point). Both errors came from
benchmarking a workload that does not represent how the method is used.

Note the Fortran references have **no GPU path at all**, so the GPU is not a
comparison against them; it is a comparison against curcuma's own CPU path.

> Aside, observed while collecting this: a cold gfnff single point on
> `mixture.xyz` (6200 atoms, **1400 fragments**) did not finish within 7 minutes,
> against ~21 s documented in
> [GFNFF_PERFORMANCE_LEVERS.md](GFNFF_PERFORMANCE_LEVERS.md). The many-fragment
> EEQ dispatch is known not to be tuned yet (operator confirmed); the exact
> solver is preferred by default for fragments in contact
> (`eeq_contact_prefer_exact`), which for 1400 fragments means a dense
> (natoms+nfrag)² factorisation. Not investigated further here.

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

## FP32 mixed-precision eigensolve — now ON by default (2026-07)

Ported from the ROCm work, where FP64 is 1/32-1/64 of FP32. The CPU/MKL
implementation (`ssygst` + `ssyevd`, `xtb_scf.cpp`) had existed since then but was
opt-in, documented as "~10% faster". That figure was measured when setup+gradient
were ~40% of runtime; after the blocked integrals the eigensolve is **58%**
(631 of 1080 ms for gfn2), so the same trick is worth much more:

| complex/231, 1 core | FP64 | mixed | eigensolve |
|---|---:|---:|---|
| gfn1 | 1027 ms | **793 ms** | 771 -> 561 ms |
| gfn2 | 1077 ms | **937 ms** | 630 -> 490 ms |

Convergence is never accepted on an FP32 step, so the fixed point is reached in
FP64. Cost over the 14-molecule reference set: energies agree to **1e-12 Eh**
(11/14 bit-identical, 3 move in the 12th decimal), gradients to **6e-7 Eh/Bohr**.
For scale, the 1e-8 tblite gate is 10000x looser than the energy shift, and the
default `scf_threshold=1e-5` already costs 1.1e-6 Eh/Bohr on its own — four times
more. Use `-scf_mixed_precision false` when maximum gradient precision matters.

## Mixed precision: reduce in FP64, diagonalise in FP32 (2026-09)

The mixed-precision branch above did the WHOLE generalized reduction in FP32
(`ssygst`) as well. On a larger basis and with threads that is a severe
pessimization, because this machine's OpenMP OpenBLAS has no optimised
single-precision triangular kernels:

| n = 3222, 8 threads (polymer, in-process probe) | FP64 | FP32 |
|---|---:|---:|
| `?sygst` (reduce) | **193 ms** | 2430 ms |
| two `?trsm` (same reduction) | **203 ms** | 2558 ms |
| `?syevd` (eigensolve) | 1341 ms | **920 ms** |

So FP32 is the right choice for the eigensolve and the wrong one for the
reduction, by more than a factor of ten. The branch now reduces with `dsygst`,
casts the reduced matrix once, and calls `ssyevd`:

| polymer (1410 atoms, nao 3222), gfn2, 36 threads | wall | reduce/it | syevd/it |
|---|---:|---:|---:|
| before, `-scf_mixed_precision true` (default) | 44.2 s | 1523 ms | 921 ms |
| before, `-scf_mixed_precision false` | 33.8 s | 221 ms | 1346 ms |
| **after, default** | **29.3 s** | 230 ms | 920 ms |
| after, `-scf_mixed_precision false` | 34.3 s | 223 ms | 1341 ms |

Energies identical to the printed 8 decimals in all four; the FP64 path is
bit-identical to before (gradient diff 0.0), and the mixed-precision gradient
moved *closer* to it (1.9e-7 -> 6.3e-8 Eh/Angstrom). `ctest`: the 65 validation
tests pass, full suite shows only the four documented pre-existing failures.

`CURCUMA_XTB_REDUCE_PROBE=1` prints the table above for the running system
(sizes and thread counts as used) - worth doing once per machine, since the
verdict depends entirely on the BLAS build.

## `-threads` now decides the eigensolve too (2026-09)

The eigensolve used to be capped at 8 threads regardless of `-threads`, because
the FP32 reduction dominated and regressed past that. With the FP64 reduction the
optimum moved, and it is machine-dependent, so the cap is gone and `-threads`
decides (`CURCUMA_EIG_MAX_THREADS` still caps it independently if wanted).
polymer, gfn2, this 36-core box: `-threads 8` 28.8 s, **`-threads 16` 24.9 s**,
`-threads 24` 27.0 s, `-threads 36` 30.2 s; complex/231 is unaffected (1.11 s,
the size gate keeps it serial). The D&C eigensolve is memory-bandwidth-bound, so
more threads than memory channels still lose - that is now a user choice.

## The three memory-bound O(nat^2)/O(nao^2) SCF loops (2026-09)

Profiling polymer (1410 atoms, nao 3222) instead of complex (231) changed the
picture: three loops that were "too small to thread" at 231 atoms dominated
everything except the eigensolve. All three sweep many separate nat x nat or
nao x nao matrices at the same index pair, so they are memory bound, not compute
bound (the multipole ones touch 18 matrices, i.e. ~286 MB per call at nat = 1410).

| polymer, gfn2, -threads 16 | before | after |
|---|---:|---:|
| potential build (`addMultipolePotential`) | 172 ms/it | **55 ms/it** |
| `energyMultipole` (inside "energy/mix") | 158 ms/it | **18 ms/it** |
| populations (GFN2 multipole moments) | 214 ms/it | **77 ms/it** |
| reduce (`dsygst` -> two `dtrsm` above 8 threads) | 231 ms/it | **169 ms/it** |
| **wall** | 24.9 s | **21.0 s** |

The potential build is threaded over the target atom with disjoint writes, so it
stays bit-identical. The energy and the populations use per-thread partial sums
added in a fixed thread order, which reassociates the outer sum: over polymer the
energy stays identical to 12 decimals and the gradient moves by 1.5e-14, and the
converged FP64 result was in fact unchanged (diff 0.0). Switching the reduction
to triangular solves is a rounding-level change (3.9e-14 on the gradient).

## What is left on the CPU (polymer, gfn2, -threads 16, 21.0 s)

setup 1.5 s, SCF 15.3 s, post-SCF 1.6 s. Inside the SCF per iteration: `dsyevd`
694 ms (54 %), reduce 169, Fock 100, populations 77, density 64, potential 55,
back-transform 50, energies 34. So the eigensolve IS the CPU floor now, and the
alternatives were measured rather than assumed:

| n = 3222, 16 threads | time |
|---|---:|
| `dsyevd` (all vectors, what curcuma uses) | **1021 ms** |
| `dsyevr` (all vectors) | 2122 ms |
| `dsyevr` (lowest 1699 = occupied + 5 %) | 2228 ms |
| curcuma's own D&C (`-eigensolver native`, in-run) | 5144 ms/it |

Partial diagonalisation does not pay even when only the occupied block is
requested, which is the same conclusion the GPU reached (AP1) for a different
reason. No MKL is installed on this machine; with MKL the floor would likely be
lower.

## D4 ATM: neighbour lists, and why threading looked broken (2026-09)

The post-SCF D4 three-body (ATM) term scaled badly with threads - polymer, 1410 atoms:
3904 ms at 1 thread, 1700 at 16, 770 at 36 (5.1x on 18 cores / 36 hardware threads).
Measurements ruled out the obvious suspects one at a time: the thread count in force is
the requested one (printed at verbosity 3), the machine is a single NUMA node, and the
stride partitioning over the outer atom balances a triangular loop. What is left is the
loop itself: it is O(N^3) over ALL triples with three distance tests inside, i.e. the
work is dominated by tests that fail.

The CUDA kernel (`k_d4_atm_nl`) has used a neighbour list for this since Jun 2026; the
CPU now does the same. `j` and `k` are taken from the neighbours of `i` instead of from
all atoms below it, the lists are ascending and scanned in the dense loop's order, so the
triple set and the accumulation order are unchanged - gradients are **bit-identical**
(polymer and complex, gfn1 and gfn2).

| polymer, ATM phase | dense | neighbour list |
|---|---:|---:|
| 16 threads | 1719 / 1564 ms | 1524 / 1543 ms |
| 36 threads | 956 / 896 ms | 802 / 821 ms |

Only ~10 %, because polymer is compact: at the 25 Bohr ATM cutoff **34 %** of all atom
pairs are neighbours (482 per atom), so the dense loop was not wasting much. The gain is
a function of that fraction, and it falls fast with size - polymer_2x (7320 atoms) has
**4.8 %** (352 neighbours per atom), where the dense loop visits ~433x more triples than
contribute. That case is not timed here: a gfn2 single point on 7320 atoms takes over an
hour on this CPU.

`CURCUMA_D4_ATM_DENSE=1` restores the dense scan for an A/B measurement.

A second idea did NOT pay and was kept only because it is free: factorising
`r0_ij = a1*sqrt(3*r4r2_i*r4r2_j) + a2` into per-atom square roots removes three `sqrt`
per triple, and changed nothing measurable - the loop is not sqrt-bound.

## The iteration-count gap is a criterion artifact, not slower iterations

gxtb converges `complex` in **15** iterations, curcuma in **19** — but per
iteration curcuma is already the faster of the two:

| | iterations | SCF time | per iteration |
|---|---:|---:|---:|
| gxtb 6.7.1 | 15 | 0.595 s | 39.7 ms |
| curcuma | 19 | 0.733 s | **38.6 ms** |

The difference is what "converged" means. curcuma tests `max|dq_sh|`
(`xtb_native.cpp`, `.cwiseAbs().maxCoeff()`); xtb tests **RMSdq**. Over 340 shells
a max-norm is systematically stricter at the same numeric threshold, so at
`scf_threshold=1e-5` curcuma converges *further* than xtb does — which the
gradient data confirms (1.1e-6 Eh/Bohr vs fully-converged).

Comparing at equal convergence instead of equal threshold:

| `-scf_threshold` | iterations | total | dE vs 1e-8 | max\|dgrad\| |
|---|---:|---:|---:|---:|
| 1e-5 (default) | 19 | 1077 ms | 2.8e-11 | 1.1e-06 |
| 5e-5 | 16 | **935 ms** | 1.2e-09 | 9.8e-06 |
| 1e-4 | 15 | 898 ms | 8.4e-09 | 2.3e-05 |

At 5e-5 — roughly xtb's effective convergence level — curcuma matches gxtb's
932 ms. The default is deliberately NOT loosened: the energy stays excellent but
the gradient degrades an order of magnitude, which matters for `-opt`/MD.

## Residual gap to xtb (gfn2 complex, 1083 vs 979 ms)

xtb runs the same MKL, so its per-iteration eigensolve is comparable; its
advantage is (a) **fewer SCF iterations** (15 vs 19 — partly because xtb's
default convergence is looser than curcuma's even after the 1e-5 change), and
(b) it folds the gradient into the SCC step whereas curcuma's setup and gradient
are separate passes. gfn1 already beats xtb because its larger basis makes the
eigensolve dominate, where curcuma's leaner setup/gradient wins.

After the 2026-07 integral work the setup is 91 ms and the gradient 175 ms, so
**the SCF is now 70% of the gfn2 runtime** (762 of 1083 ms) and the iteration
count is the dominant remaining lever — not the integrals. A looser convergence
(`-scf_threshold 5e-5` → 16 it, energy-identical) closes most of what is left.

## Verification

```bash
cd release && make -j16 curcuma
ctest -R "sqm_val_|ecomp_|d4_diag|d4_dedq|xtb_gradient|xtb_cpscf|gfn1_align|gfn2_align"
# 45/45 pass; energies bit-stable. Current values (2026-07):
#   gfn2 H2O     = -5.070369821862 Eh
#   gfn2 complex = -329.527147840995 Eh
#   gfn1 complex = -343.179803543151 Eh
# NOTE the older -329.52707823 / -343.17980352 in earlier revisions of this doc
# predate the electronic free-energy (Fermi entropy) term and the F5 variational
# D4 fix. Use `-dump_gradient` (12-dp energy, 14-digit gradient) as the gate;
# the 8-dp "Single Point Energy" print is far too coarse at -329 Eh.
scripts/sqm_bench.sh 3
```

Pre-existing, unrelated ctest failures (not touched here): `AAAbGal_dtemplate`,
`gfnff_numgrad_fixed_charges`, `cli_gfnff_01/02`, `cli_sqm_11` (ipea1 needs
TBLite, not built).
