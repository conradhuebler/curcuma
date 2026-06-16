# Work Package: Vulkan GPU eigensolver — closing the gap to rocSOLVER

**Status (2026-06-16): COMPLETE for the planned items.** V-AP5 replaced the cyclic Jacobi
with a Householder tridiagonalization eigensolver (`solveSymTridiag` in
`src/core/energy_calculators/qm_methods/vulkan/xtb_vulkan_context.cpp`); this WP then made it
fully GPU-resident (EIG-1), parallelised the back-transform (EIG-2A), blocked it into BLAS-3
GEMMs via compact-WY (EIG-2B), and added an opt-in FP32 mixed-precision path (EIG-4). The
eigensolve is **no longer the overwhelming bottleneck** it was — it is now ~2× ROCm rather than
~40×, and the residual Vulkan↔ROCm gap is spread across the eigensolve tail (tql2 +
back-transform) and the non-eigensolve SCF/setup. EIG-3 (host tql2) is deferred as a
diminishing return.

> 🤖 AI-generated. Done items are ⚙️ machine-tested (correctness validated); open items are
> proposals. Numbers are on a Radeon 890M (RADV). Only ✅ TESTED/APPROVED by the operator are
> validated. ⚠️ The CPU/GPU clocks were raised mid-development and the box runs other jobs —
> **only trust fresh, same-session, back-to-back numbers**; cross-session ms are not comparable.

---

## Progress

| item | status | result |
|------|--------|--------|
| EIG-0 profiling hook | ✅ | `CURCUMA_VK_EIG_PROFILE=1` (+ `-FP32` line); committed |
| EIG-1 fully-GPU tridiagonalization | ✅ | `tri_house`+`tri_kw` fold the host ‖x‖/pᵀv reductions onto the GPU; tridiag fence syncs **1112 → 3** (n=558) |
| EIG-2A workgroup-per-column back-transform | ✅ | `tri_applyl` 256-lane shared reduction; submits → 3 |
| EIG-2B WY-blocked GEMM back-transform | ✅ | `tri_vfull`+`wy_buildt`+`gemm_g` (compact-WY panels, b=32); back-transform **27.7 → ~20 ms (1.4×)**, bit-identical |
| EIG-3 host tql2 | deferred | sequential host O(n³), now ~29 ms (smallest phase); not a clean win — see re-scope |
| EIG-4 FP32 mixed precision | ✅ opt-in | `solveSymTridiag32`; helps the tridiagonalization (1.5×) only — back-transform/tql2 don't; net ~5%, `-scf_mixed_precision true` |

**Fresh same-clock (`complex`/231 GFN2 `-sp`, shared box ⇒ ±noise):**

| eigensolve / iter | total | vs ROCm |
|---|---|---|
| cyclic Jacobi (the old default) | 3016 ms | 61 s | — |
| tridiag FP64 (EIG-1+2A) | **~129 ms** *(solve-eigen timer 161 ms incl. X·F·X / X·C̃ GEMMs)* | **~4.2 s** | 2.0× / 2.5× |
| tridiag FP32 mixed (EIG-4, opt-in) | ~106 ms | ~4.0 s | 1.6× / 2.4× |
| ROCm `dsygvd`/`ssygvd` | 66 ms | 1.67 s | 1.0× |

Per-phase (`EIGPROF`, fresh): FP64 tridiag **~66**, tql2 **~30**, back-transform **~20** ms
(EIG-2B WY; was 30 with EIG-2A). **The remaining ~2× to ROCm is now algorithmic (tridiag +
tql2 + back-transform tail) and spread into the rest of the SCF — not a bug.**

---

## 0. Original baseline (EIG-0 — DONE; numbers below predate the clock bump)

Per-phase timing is wired behind `CURCUMA_VK_EIG_PROFILE=1` (prints one `[EIGPROF]` line
per `solveSymTridiag` call). `complex`/231, GFN2 `-sp`, n=558, steady-state SCF iteration:

| phase | time/iter | what it is | why it costs |
|-------|-----------|------------|--------------|
| **tridiagonalization** | **166 ms** | GPU `tri_matvec` (p=A·v) + `tri_rank2` (A−=vwᵀ+wvᵀ) per column; host does ‖x‖, pᵀv, v/β/w | **1112 host↔GPU fence syncs** (2·(n−2)): the host reductions sit in the per-column loop, so every step is a `submit()`+`vkWaitForFences` round-trip |
| **host tql2** | **73 ms** | EISPACK `cpuTriEig` on the n×n tridiagonal → eps + Z | single-core O(n³) eigenvector (Givens) accumulation |
| **back-transform** | **100 ms** | GPU `tri_applyl`: C = Q·Z (apply n−3 reflectors) | one thread **per Z column**, each doing a serial length-m `vᵀZ` reduction + rank-1 update → underutilized + serial |
| **total** | **~337 ms** | | ROCm reference ≈ **70 ms** |

Take-away: all three phases are improvable, and none is intrinsically O(20·n³) anymore (that
was the Jacobi). The tridiagonalization is **latency-bound** (fence syncs, not FLOPs); the
back-transform is **occupancy-bound**; tql2 is **single-core-bound**.

---

## 1. Work items

### EIG-1 — Fully-GPU tridiagonalization (kill the 1112 fence syncs) — ✅ DONE (2026-06-16)
Implemented as designed: `tri_house` + `tri_kw` (one-workgroup shared-memory reductions)
fold the host ‖x‖/pᵀv reductions onto the GPU; scalars live in the combined `tScal`
buffer `[diag|off|β]`; `tri_house→tri_matvec→tri_kw→tri_rank2` record into chunked
submits. Fence syncs **1112 → 3** (n=558). Correctness unchanged (residual 3.2e-14).

**Problem:** the hybrid does the O(n) reductions on the host (via the mapped coherent
buffers), so each Householder step needs the GPU result back before the host can build `w`,
forcing 2 blocking submits/step.

**Approach:** move the per-step scalar work onto the GPU so the whole reduction records into
a few chunked command buffers (exactly like the already-batched back-transform / chunked
Jacobi):
- `tri_house` (one workgroup): reduce ‖x‖² over column k (shared-memory tree), then build
  `v` in-place and `β`, `e[k]` — needs only `‖x‖` and `x₀` (closed form already derived in
  `solveSymTridiag`).
- extend `tri_matvec` (or a `tri_kw` kernel) to also reduce `pᵀv` and write `w = β·p − K·v`
  on the GPU (one workgroup reduction + elementwise).
- record `tri_house → tri_matvec → tri_kw → tri_rank2` for all k into chunked submits
  (`jacobiRoundsPerSubmit`-style budget); the host stays out of the loop.

**Submits:** ~1112 → ~O(n / rps) ≈ a handful. **Expected: 166 → ~25–40 ms.**
**Effort:** medium (2 small reduction kernels + reorchestration). **Risk:** low — FP64
reduction order changes ~1e-13; validate with `CURCUMA_VK_TEST_TRIDIAG`.

### EIG-2A — Parallelize the reflector back-transform — ✅ DONE (2026-06-16)
Option A implemented: `tri_applyl` is now one workgroup (256 lanes) per Z column with a
shared-memory tree reduction of `vᵀZ` + parallel rank-1 update; submits batch to 3.
(Option B / WY-blocked GEMM back-transform later landed as **EIG-2B** below.)

### EIG-2B — WY-blocked GEMM back-transform — ✅ DONE (2026-06-16)
Option B implemented. The n−3 per-reflector passes over Z are grouped into panels of b=32
reflectors; each panel's product H_{k0}·…·H_{k0+b-1} is the compact-WY form `Q = I − V·T·Vᵀ`
and applied as three rectangular FP64 GEMMs: `Y = Vᵀ·Z` (b×n), `Y2 = T·Y` (b×n), `Z −= V·Y2`
(n×n) — cutting Z memory passes from ~n to ~3·(n/b) and reusing V tiles across all columns.
New shaders: **`tri_vfull`** (scatter the compact panel reflectors to a full n×b layout,
β=0 columns zeroed so degenerate reflectors drop out), **`wy_buildt`** (build the b×b T from
Gram(Vfull)+β by the LAPACK forward/columnwise recurrence, one workgroup), **`gemm_g`** (a
general rectangular M×N×K tiled GEMM with transA + subtract-accumulate — the existing
`gemm.comp` is square-only, so WY needs this). All panels record into one chunked submit (no
host readback between them). **Result (fresh `complex`/231 GFN2 `-sp`): back-transform
27.7 → ~20 ms (~1.4×), energies bit-identical** to EIG-2A/CPU (gfn1+gfn2 over the sqm set
incl. the 231-atom `complex`, dE = 0; standalone eigensolver resid 3.4e-14 / ortho 8.2e-15 to
n=558; gfn2 `-opt` identical; `ctest -L gpu` 2/2). The win is modest (the b=32 tall-skinny
GEMMs are still occupancy-limited on the iGPU, and tql2 noise swamps the *total*), but it is a
real, isolated back-transform speed-up with no downside, so WY is the default
(`CURCUMA_VK_BACKXF=reflector` reverts to EIG-2A).

**Problem:** `tri_applyl` is one-thread-per-column with a serial length-m inner loop (same
underutilization the gradient had before V-AP6).

**Approach (pick one):**
- **A (low-risk):** one **workgroup per Z column** with a shared-memory reduction of `vᵀZ`,
  then a parallel rank-1 update — the V-AP6 pattern applied to the back-transform.
- **B (faster, more work):** block the reflectors via the **compact WY representation**
  (`Q = I − V·T·Vᵀ` for a panel of `b` reflectors) and apply them as **two GEMMs** using the
  already-tiled `gemm.comp`. Turns the back-transform into BLAS-3 the iGPU runs efficiently.

**Expected: 100 → ~25–40 ms.** **Effort:** A small / B medium. **Risk:** A low; B medium
(panel/triangular-T bookkeeping).

### EIG-3 — Speed up the tridiagonal solve (`cpuTriEig`/tql2)
**Problem:** single-core O(n³) eigenvector accumulation (73 ms).

**Options:**
- **(a)** thread the Givens-rotation loop over Z rows on the project `CxxThreadPool`
  (the `for k in 0..n` z-update at each QL step is embarrassingly parallel) — modest, low risk.
- **(b)** eigenvalues-only on host (O(n²), bisection/`sterf`) + eigenvectors on the GPU via
  Cuppen **divide-and-conquer** (the existing CPU `tridiagDC`/`rank1Eigen` in
  `native_eigensolver.cpp` is a portable reference) or inverse iteration — bigger, and
  inverse iteration needs cluster reorthogonalization (robustness risk on degenerate spectra).
- **(c)** leave it — it is the smallest phase.

**Expected: 73 → ~30 ms (a) / ~20 ms (b).** **Effort:** (a) low / (b) high. **Risk:** (a)
low / (b) medium-high (correctness on clustered eigenvalues).

> **Re-scope (2026-06-16):** fresh-measured tql2 is now only **~29 ms** (the smallest of the
> three phases). Option (a) threads a per-QL-iteration O(n) Givens loop — fine-grained, so the
> pool overhead likely eats the gain; (b) is high-effort/high-risk. **Defer EIG-3** unless EIG-4
> lands and tql2 becomes the critical path. The biggest remaining lever is the **GPU** work
> (tridiag 69 + back 30 = 99 ms) → **EIG-4**.

### EIG-4 — FP32 mixed precision (parallels ROCm's lever) — ✅ DONE (opt-in, 2026-06-16)
Implemented `solveSymTridiag32`: cast Ã→FP32, FP32 `tri_*_f32` tridiagonalization +
back-transform on FP32 mirror buffers (`tA32`…`tZ32`), FP64 `tql2` on the cast-to-double
tridiagonal (host), cast eps/vectors back. Routed via `solveAtilJacobi`'s `fp32` branch
(`-scf_mixed_precision true`; FP64 corrector once `max|dq| < scf_fp32_threshold` → energy FP64).

**Result (fresh, `complex` GFN2 `-sp`, shared box):** only the tridiagonalization benefits —
**76→51 ms/it (1.5×)**, arithmetic/bandwidth-bound as predicted; the **back-transform
(28→26 ms) and host tql2 do NOT** (latency/host-bound — the *same* lesson as the X-AP3 FP32
Jacobi). Net SCF ~3700→~3560 ms (**~5%**), energy bit-identical (10/10 sqm). **Kept opt-in**
(not defaulted): modest gain + the FP32 fixpoint shifts the gradient ~1e-7 at the loose
default `scf_threshold`. A future EIG-2B (FP32 *and* WY-blocked GEMM back-transform) could
extend the win to the back-transform, but the back-transform isn't FP64-arithmetic-bound, so
the lever there is occupancy/algorithm, not precision.

---

## 2. Outcome (delivered order: EIG-1 → EIG-2A → EIG-4)

Fresh same-clock (`complex` GFN2 `-sp`, shared box ⇒ ±noise):

| | Jacobi (old) | EIG-1+2A (FP64) | + EIG-4 (FP32, opt-in) | ROCm |
|--|----|--------------|---------|------|
| eigensolve / iter | 3016 ms | **~129 ms** | **~106 ms** | 66 ms |
| `complex` GFN2 `-sp` | 61 s | **~4.2 s** | **~4.0 s** | 1.67 s |

EIG-4's gain landed lower than the ~1.5–2× hoped for the *GPU phases*: FP32 sped up only the
tridiagonalization (arithmetic-bound, 1.5×); the back-transform and host tql2 did not (the
X-AP3 lesson — they are latency/host-bound), so the net is ~5% and it stays opt-in. The
residual ~2× to ROCm beyond the eigensolve is the host integral build / γ / D4 and setup —
tracked separately in `docs/SQM_PERF_OPT_WP.md`.

---

## 3. Validation (all items)
- `CURCUMA_VK_TEST_TRIDIAG=1`: standalone eigensolver vs the Jacobi/Eigen — residual
  ‖Av−λv‖ ≤ ~1e-12, orthonormality ≤ ~1e-14, to n=558.
- Energies **bit-identical** to ROCm/CPU on the 12-molecule `sqm_reference` set incl. the
  231-atom `complex`, GFN1 **and** GFN2 (`-sp`, 8 dp).
- GFN1/GFN2 `-opt` (NH3, triose) converge to ROCm's minimum.
- `CURCUMA_VK_EIG_PROFILE=1` per-phase numbers before/after each item.
- `ctest -L gpu` (Vulkan) green; FP32 path (EIG-4) keeps the converged energy FP64.

---

## 4. Closing notes / what's left
- **Done & committed:** `be159a1` (V-AP5 eigensolver + device-lost fix), `b67834d`
  (EIG-1 + EIG-2A), `b0e5161` (EIG-4 FP32). The `CURCUMA_VK_EIG_PROFILE` hook (EIG-0) and the
  `CURCUMA_VK_TEST_TRIDIAG` standalone check are committed.
- **EIG-2B done (2026-06-16):** WY-blocked GEMM back-transform — `tri_vfull`+`wy_buildt`+
  `gemm_g`, back-transform 27.7 → ~20 ms (~1.4×), bit-identical. Default ON
  (`CURCUMA_VK_BACKXF=reflector` reverts to EIG-2A).
- **Deferred (diminishing returns, do on a quiet box where the small gains are measurable):**
  - **EIG-3** — host tql2 (~29 ms): GPU Cuppen divide-and-conquer (`tridiagDC`/`rank1Eigen`
    in `native_eigensolver.cpp` is a portable reference) or threaded Givens. High effort / risk.
- **Beyond the eigensolve:** the rest of the Vulkan↔ROCm gap is the host integral build / γ /
  D4 and setup — tracked separately in `docs/SQM_PERF_OPT_WP.md`.
- **Naming collision (open):** the perf commits used `V-AP5`/`V-AP6`, which clash with
  `docs/SQM_GPU_ROADMAP.md`'s planned `V-AP5` (device potential build) / `V-AP6` (GFN-FF).
  Reconcile — renumber the roadmap stages or retag the perf line `V-PERF-*`.
