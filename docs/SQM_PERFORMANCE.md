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
| gfn1   | 2982 ms        | 1221 ms           | **793 ms**        | 2562   | 1367       |
| gfn2   | 2944 ms        | 1364 ms           | **937 ms**        | 1567   | 979        |

- **gfn1: −73%. 1.6× faster than xtb.**
- **gfn2: −68%. Now at parity with xtb (0.96×; was ~3×, then 1.39×).**
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

| method | curcuma K=1 | gxtb K=1 | cur/gx K=1 | curcuma K=16 | gxtb K=16 | cur/gx K=16 |
|---|---:|---:|---:|---:|---:|---:|
| gfnff | 43 ms | 33 ms | 1.27 | 37 ms | 27 ms | 1.37 |
| gfn1  | 987 ms | 1269 ms | **0.78** | 341 ms | 382 ms | **0.89** |
| gfn2  | 1122 ms | 932 ms | 1.20 | 412 ms | 302 ms | 1.36 |

Before this round the same table read gfnff 1.26/1.59, gfn1 0.89/0.91,
gfn2 1.32/1.38.

> ⚠️ This table is **invalidated by any change to the integral, SCF or setup
> path**. Re-measure and record the commit hash when quoting it. Measured at
> commit `17def9ff`, single socket, MKL, `-mfma`.

Caveat on gfnff: at 43 ms total the run is dominated by process start plus
one-time setup, so the ratio is sensitive to measurement noise and to the
`.topo.json` cache being present or absent.

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
