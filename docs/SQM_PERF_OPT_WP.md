# Work Package: Native GFN1/GFN2 performance — CPU + GPU

> 🤖 AI-prepared work package (2026-06-04), not yet implemented. Companion to
> [docs/SQM_GPU.md](SQM_GPU.md), [docs/SQM_PERFORMANCE.md](SQM_PERFORMANCE.md) and
> the threading/eigensolve notes. Goal: make `-method gfn1|gfn2` faster on **both**
> the threaded CPU and the CUDA GPU path, without touching the validated energies
> (1e-8 vs tblite) or gradients (~1e-15 vs the CPU reference).
> Per CLAUDE.md AI-policy: every item below is a proposal; nothing here is
> implemented or human-tested. Mark new code "Claude Generated"; do **not**
> self-assign ✅ TESTED.

## 0. Measured baseline (complex / 231 atoms, `-sp`, `-threads 8`)

Phase breakdown from `-verbosity 3` (after the D3 + mixed-precision + auto-thread
work of 2026-06-04):

| phase | GFN2 GPU | GFN2 CPU | GFN1 GPU | note |
|---|---|---|---|---|
| setup (integrals) | ~50 ms | ~50 ms | ~55 ms | overlap+H0, γ, GFN2 multipole AO ints |
| **SCF (19 it)** | **~358 ms** | **~354 ms** | ~254 ms | **eigensolve ~280 ms + host potential ~68 ms** |
| post-SCF E | ~75 ms | ~75 ms | ~24 ms | GFN2 **D4** dispersion / GFN1 D3 (now fast) |
| gradient | ~16 ms | ~44 ms | ~15 ms | device vs host |
| **TOTAL** | **~519 ms** | **~556 ms** | **~345 ms** | + one-time ~300 ms GPU-context init (wall) |

Eigensolve detail (GFN1, per-iter sum): reduce 49 ms / **dsyevd 164 ms** / back 11 ms.
CPU eigensolve scales ~2.2× (t1 608 → t8 272 ms — sub-linear, as dense D&C does).
GPU eigensolve is mixed-precision (FP32 early / FP64 polish): 519→278 ms already.

**The SCF (eigensolve + host potential) is ~68 % of the run. It is the target.**

## 1. Tier-1 — biggest, cross-cutting (CPU **and** GPU)

### AP1 — Partial diagonalization (only occupied + buffer eigenpairs)

- **Why:** every SCF iteration solves the full `nao×nao` spectrum, but the density
  only needs ~`nocc` vectors; Fermi smearing (default 300 K) needs a thin window a
  few `kT` above HOMO. At ~50 % occupancy, solving `nocc(+buffer)` instead of `nao`
  is ~**1.5–2× on the dsyevd** (the single largest cost).
- **CPU** (`xtb_scf.cpp::solveEigen`): replace LAPACK `dsyevd` with `dsyevr`/`dsyevx`,
  range `il=1 … iu=nocc+buffer` (e.g. buffer = max(8, 0.05·nao)). The reduction
  (cached `L` + trsm) and back-transform are unchanged.
- **GPU** (`cuda/xtb_gpu_context.cu::eigensolveResidentFock`): replace
  `cusolverDnDsyevd` with `cusolverDnDsyevdx` (expert API, eigenvalue/index range);
  keep the FP32/FP64 split.
- **Risk:** the occupation bisection (`occupationsFromEps`) needs eigenvalues
  spanning the Fermi window — validate that `iu` always covers it (warn + widen, or
  fall back to full solve, if the highest computed eigenvalue is within a few `kT`
  of the Fermi level). Metallic / tiny-gap systems may need a larger buffer.
- **Impact:** ~30–40 % of the eigensolve → ~100 ms on complex, both paths.
- **Effort/risk:** medium / medium. **Validation:** `gfn{1,2}_validation` +
  `gpu_gfn{1,2}_validation` @1e-8 stay green; a tiny-gap molecule check.

### AP2 — D4 dispersion hot loop (GFN2 post-SCF ~75 ms)

- **Why:** same class as the D3 win (commit 4b41562). `d4_evaluator.cpp` calls
  `D4ParameterGenerator::weightedC6Gfn2` per pair (the AP6b exact path,
  `per_reference_charge=true`), which re-walks the 7×7 reference C6 block per pair.
- **Approach:** hoist the element-pair reference block out of the inner ref loop
  (mirror `D3ParameterGenerator::refC6Block`, commit 4b41562); cache per distinct
  element pair. The pair list is already geometry-fixed and the energy loop is
  already threaded (`d4_evaluator.cpp:203`), so the win is the per-pair block walk.
- **Impact:** ~75 → ~40 ms (GFN2), both paths. **Effort/risk:** low / low (energy
  bit-identical; gate `gfn2_validation` + `test_d4_*`).

## 2. Tier-2 — GPU-specific

### AP3 — Isotropic potential build on the device (the "option b" refactor)

- **Why:** the per-iteration `v_sh = γ·q_sh (+ third-order + in-SCF D4 scalar)` still
  runs on the host (~68 ms threaded) and crosses the bus each iteration. The device
  already has γ resident and computes the shell charges.
- **Approach:** keep `q_sh`/`v_sh` on the device; a `k_coulomb_vsh` (γ·q_sh gemv via
  cuBLAS) + the third-order shift fold into `residentSolve`. The host keeps only the
  pieces that genuinely need host data (D4 charge coupling can be downloaded as a
  length-`nat` vector, as today). Removes the per-iter v_ao round-trip.
- **Refs:** host side `addCoulombShellPotential`/`addThirdOrderPotential`
  (`xtb_coulomb.cpp`/`xtb_thirdorder.cpp`), the v_ao expansion + `solve(v_ao)` in
  `xtb_native.cpp`.
- **Impact:** ~68 ms/iter region + removes the per-iter transfer; fully on-device SCF.
- **Effort/risk:** medium-high / medium (restructures the potential assembly).

### AP4 — Skip the redundant host integral build (download from device)

- **Why:** on the GPU-resident path the host *still* runs
  `getHamiltonianH0`/`buildGammaMatrix`/`setupMultipole`(AO)/`buildOrthonormalizer`
  (~47 ms threaded) even though `computeIntegrals` rebuilds them on the device.
- **Approach:** when the device build will be used, skip those host calls and
  download `dS/dH0/dL/dGamma/dDpInt/dQpInt` into `m_S/m_H0/m_X/m_gamma/m_dp_int/m_qp_int`
  (the `download*` methods already exist on `XtbGpuContext`); keep the cheap host
  part of `setupMultipole` (mrad/amat, O(nat²)). Requires moving the device
  `beginBasis`+`computeIntegrals` ahead of the host integral build in
  `XTB::Calculation` and a guard so the CPU fallback path is unchanged.
- **Impact:** ~47 ms/step on `-opt`/`-md`. **Effort/risk:** medium / medium-high
  (setup-ordering restructure). **Validation:** needs a small **multi-step** test
  (a full `-opt` is too slow to iterate on) — e.g. a 2–3-cycle `-opt` GPU-vs-CPU
  energy check, plus `gpu_gfn{1,2}_validation`.

### AP5 — (Big bet) Device density-matrix purification — no eigensolve

- **Why:** the FP64 eigensolve is the GPU wall (GeForce FP64 ≈ 1/64 FP32).
  Palser–Manolopoulos purification builds P directly via **GEMM** (no
  diagonalization) at T=0 with a HOMO–LUMO gap. On a GeForce, FP32 GEMM is the
  fastest op; for gapped organic molecules (gap ≫ kT ⇒ T=0 ≈ 300 K energy) this can
  replace the eigensolve with ~10–20 FP32 GEMM iterations of O(N³).
- **Refs:** the CPU `purify` path (`xtb_sparse.cpp`, `-eigensolver purify`) is the
  blueprint; it already supplies `W` directly for the gradient (no eigenpairs).
- **Impact:** potentially the **largest GPU win** on gapped systems; could make the
  GPU pull decisively ahead of the CPU.
- **Effort/risk:** high / high (device purification kernels + the T=0/gap gate +
  Fermi-occupation handling; keep the eigensolve as the default/fallback).

## 3. Tier-3 — smaller / investigate

- **AP6 — CPU eigensolve threading audit.** It threads (~2.2× at t8) but
  `CMakeLists.txt:964` notes a `libmkl_sequential`-loads-first conflict; verify the
  threaded layer (`mkl_gnu_thread`, `CMakeLists.txt:1009`) actually drives `dsyevd`
  (check `OMP_NUM_THREADS`/`mkl_set_num_threads` reaches it). Possible free CPU win.
- **AP7 — CUDA stream overlap.** Run the host potential build concurrently with the
  device kernels (separate stream / events) to hide latency.
- **AP8 — Reuse resident dP/dC in the GPU gradient.** `computeGradient` re-uploads
  P/C each step; the device already has them resident after the SCF (~2 ms/step).
- **AP9 — More aggressive FP32 (GPU).** Raise `scf_fp32_threshold` so fewer
  iterations are FP64. Easy but accuracy-risky (the `!m_eig_fp32` convergence guard
  must still hold the 1e-8 gate).

## 4. Suggested order

1. **AP1 (partial diagonalization)** — attacks the 68 % on both CPU and GPU;
   compounds with the GPU mixed precision already in place.
2. **AP2 (D4 hoist)** — quick, safe GFN2 win.
3. **AP4 (skip redundant host build)** + **AP3 (potential on device)** — the GPU
   per-step `-opt`/`-md` levers (validate with a small multi-step test, not a full
   optimization).
4. **AP5 (device purification)** — the high-ceiling GPU bet for gapped systems.

## 5. Validation (all items)

- Energies: `gfn{1,2}_validation` + `gpu_gfn{1,2}_validation` @1e-8 vs tblite stay
  green with the documented xfails (He2, complex).
- Gradients: `ctest -L gpu_gradient` (device vs CPU) stay green; a 2–3-cycle `-opt`
  GPU-vs-CPU final-energy check for the multi-step (AP3/AP4) items.
- No-CUDA `release/`: `gfn{1,2}_validation` unchanged (host changes are `#ifdef`-free).
- `compute-sanitizer` on any new device kernel. Bench harness: `scripts/sqm_bench.sh`.
