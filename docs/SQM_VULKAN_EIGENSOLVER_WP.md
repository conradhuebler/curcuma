# Work Package: Vulkan GPU eigensolver — closing the gap to rocSOLVER

**Status (2026-06-15):** V-AP5 replaced the cyclic Jacobi with a Householder
tridiagonalization eigensolver (`solveSymTridiag` in
`src/core/energy_calculators/qm_methods/vulkan/xtb_vulkan_context.cpp`), cutting the
`complex`/231 GFN2 eigensolve from **2629 → 367 ms/iter (7.2×)**. The eigensolve is now
the single dominant cost of a Vulkan GFN2 `-sp` (≈7 s of 8.7 s on `complex`). This WP
plans the remaining work to approach ROCm's **~70 ms/iter** (`rocsolver_dsygvd`/`ssygvd`,
FP32-mixed).

> 🤖 AI-generated plan. Items are proposals with estimated gains; numbers are measured on a
> Radeon 890M (RADV) unless noted. Only ✅ TESTED/APPROVED by the operator are validated.
> ⚠️ The CPU/GPU clocks were raised mid-development and the box runs other jobs — **only
> trust fresh, same-session, back-to-back numbers**; cross-session ms are not comparable.

---

## Progress (2026-06-16)

| item | status | note |
|------|--------|------|
| EIG-0 profiling hook | ✅ | `CURCUMA_VK_EIG_PROFILE=1` |
| EIG-1 fully-GPU tridiagonalization | ✅ DONE | `tri_house`+`tri_kw` fold the host ‖x‖/pᵀv reductions onto the GPU; tridiag fence syncs **1112 → 3** for n=558 |
| EIG-2A workgroup-per-column back-transform | ✅ DONE | `tri_applyl` 256-lane shared reduction; submits → 3 |
| EIG-3 tql2 | open | sequential host O(n³); not a clean win (see below) |
| EIG-4 FP32 mixed precision | open | the real remaining lever; gain uncertain on this iGPU |

**Fresh same-clock (`complex`/231 GFN2 `-sp`, shared box ⇒ ±noise), eigensolve/iter & total:**
cyclic Jacobi 3016 ms / 61 s → tridiag (EIG-1+2A) **161 ms / 4.35 s** (18.7×) → ROCm 66 ms / 1.67 s.
Per-phase (`EIGPROF`, fresh): tridiag **69 ms**, host tql2 **29 ms**, back-transform **30 ms**
(= 128 ms standard eigensolve; the `solve eigen` timer's extra ~33 ms is the X·F·X / X·C̃ GEMMs).
**Remaining gap to ROCm is ~2.4×, now algorithmic (FP32 + tql2), not a bug.**

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
(Option B / WY-blocked GEMM back-transform remains a possible future refinement.)

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

### EIG-4 — FP32 mixed precision (parallels ROCm's lever) — do after EIG-1/2
Once the GPU phases are batched, add FP32 variants of `tri_house`/`tri_matvec`/`tri_rank2`/
`tri_applyl` (+ FP32 tql2) used while `max|dq| > scf_fp32_threshold`, reverting to FP64 near
convergence (the `-scf_mixed_precision` plumbing + `m_eig_fp32` flag already exist; the
X-AP3 note found FP32 did **not** help the *Jacobi* because it was dispatch-bound — but the
tridiagonalization is arithmetic/bandwidth-bound, so FP32 should pay here).

**Expected: ~1.5–2× on the GPU phases.** **Effort:** low-medium (FP32 shader variants reuse
ploJ). **Risk:** FP32 Householder accuracy far from convergence (bounded — FP64 corrector
near convergence, energy stays FP64).

---

## 2. Suggested order & target

1. ~~**EIG-1**~~ ✅ — biggest, prerequisite for EIG-4.
2. ~~**EIG-2A**~~ ✅ — reused the V-AP6 reduction pattern.
3. **EIG-4** (FP32) — the remaining lever (targets the 99 ms of GPU phases).
4. **EIG-3** — deferred (tql2 now only ~29 ms; see re-scope).

Fresh same-clock (`complex` GFN2 `-sp`, shared box):

| | Jacobi (old) | EIG-1+2A (now) | + EIG-4 (est.) | ROCm |
|--|----|--------------|---------|------|
| eigensolve / iter | 3016 ms | **161 ms** | ~100–120 ms | 66 ms |
| `complex` GFN2 `-sp` | 61 s | **4.35 s** | ~3.5–4 s | 1.67 s |

(The residual gap to ROCm beyond the eigensolve is the host integral build / γ / D4 —
separate, see `docs/SQM_PERF_OPT_WP.md`. EIG-4's gain is genuinely uncertain on this iGPU —
the X-AP3 FP32 Jacobi was net-neutral because it was dispatch-bound; the tridiag is
memory/ALU-bound so FP32 *should* help, but **validate on a quiet box** before trusting ms.)

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

## 4. Housekeeping
- **Naming collision:** the perf commits used `V-AP5` (eigensolver) / `V-AP6` (gradient),
  which clash with `docs/SQM_GPU_ROADMAP.md`'s planned `V-AP5` (device potential build) /
  `V-AP6` (GFN-FF on Vulkan). Reconcile — e.g. renumber the roadmap's feature stages or
  retag the perf line `V-PERF-*`. Update `SQM_GPU_ROADMAP.md` either way.
- The `CURCUMA_VK_EIG_PROFILE` hook is currently an uncommitted local addition in
  `xtb_vulkan_context.cpp` (EIG-0 deliverable); commit it with EIG-1.
