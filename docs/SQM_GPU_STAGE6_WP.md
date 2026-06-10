# Work Package: Stage 6 — fully device-resident GFN1/GFN2 SCF iteration

> ✅ **IMPLEMENTED + machine-validated for GFN2 (2026-06-05)** — see the
> [Stage 6 section in docs/SQM_GPU.md](SQM_GPU.md#stage-6--fully-device-resident-gfn2-scf-iteration-2026-06-05-done--validated)
> for what landed, the validation, and the honest performance note (residency
> milestone, not a measured `-sp` speed-up; the GFN1 fused loop was deferred —
> GFN1 keeps its Stage-2a resident path). The WP below is the original proposal.
>
> 🤖 AI-prepared work package (2026-06-04). Companion to
> [docs/SQM_GPU.md](SQM_GPU.md) and [docs/SQM_GPU_STAGE5_WP.md](SQM_GPU_STAGE5_WP.md).
> Stage 5 (committed: EEQ charge model, device `q_at`, in-SCF D4 `dE/dq`, and the
> **whole** GFN2 potential build — `γ·q_sh` + shell third-order + multipole
> `v_dp`/`v_qp`/`v_at` scalar shift) made the per-iteration *potential* and the
> integrals/eigensolve/density/gradient device-resident. This WP specifies the
> **remaining host work inside the SCF loop**: occupation, populations, charge
> mixing, energy and convergence — the pieces that force the per-iteration
> `eps`/`pop_ao`/`q_sh` round-trips. Per CLAUDE.md AI-policy: every item is a
> proposal; mark new code "Claude Generated"; do **not** self-assign ✅ TESTED.

## 0. Why this is the real lever

Measured on complex/231, GFN2, `-gpu cuda` (Stage 5 in place):

| SCF phase | share | bound by |
|---|---|---|
| **eigensolve + per-iter round-trips** | **~81 %** | **latency** (the `eps`↓ / `occ`↑ / `pop_ao`↓ / `q_sh`↑ syncs), *not* eigensolve compute |
| potential build | now device (Stage 5) | — |
| populations / energy / mix / convergence | small host compute | — |

Stage 5's honest value note already flagged this: FP32≈FP64≈partial-diag all land
at ~330 ms SCF, i.e. the SCF is **latency-bound by the per-iteration host
round-trips**, not by eigensolve FLOPs. Stage 5 swapped the `v_ao` upload for a
`q_sh`/`dp_at`/`qp_at` upload (≈ neutral transfer) and moved the potential compute
to the device — but the host still drives the loop, so **every iteration still
crosses the bus** for `eps` (down), `occ` (up), `pop_ao` (down), the multipole
moments (down), and the mixed SCC (up). This WP removes those.

**Target:** per SCF iteration nothing of size > O(1) crosses the bus — only a
scalar `max|Δq|` (+ the running energy) comes down for the host convergence check.
The eigensolve stays on the device (already there). The expected win is the bulk
of the ~81 % latency region on iterative runs (`-opt`/`-md`), where the context is
warm and the round-trip count dominates.

## 1. What still runs on the host each iteration (the targets)

From `XTB::Calculation` (xtb_native.cpp) per iteration, post-Stage-5:
1. `occupationsFromEps(eps, occ, ncol)` — Fermi/integer occupations from `eps`
   (needs `eps` ↓, returns `occ` ↑).
2. `updatePopulationsFromPopAo(pop_ao)` + `multipoleMoments` — `q_sh`/`q_at`
   (device B1 already can) and `dp_at`/`qp_at` from the resident density (needs
   `pop_ao` ↓ / moments ↓).
3. `energyCoulombShell()` + `energyThirdOrder()` + `energyMultipole()` + `gpu_band`
   — the SCC energy (host, from charges + `γ`/`amat`).
4. Broyden mixing of `[q_sh; dp_at; qp_at]` (`broyden.update`, host) → next input
   (`q_sh`/`dp_at`/`qp_at` ↑).
5. `checkConvergence` on `max|Δq_sh|` and `ΔE` (host scalars).

## 2. Work items

### WP-S6.1 — device occupation (Fermi) + band/electronic scalars
- Port `occupationsFromEps` (Fermi smearing at `electronic_temperature`, integer
  fill at 0 K, the Nelec constraint + the chemical-potential bisection) to a device
  kernel that consumes the resident `eps` and writes a resident `occ`. The density
  build already takes `occ`; keep it resident → no `eps`↓/`occ`↑ round-trip.
- Reductions (electron count, band energy `Σ occ·eps`) via cuBLAS/`cub`.
- **Caveat:** the bisection is a small serial loop — fine on one block; validate the
  chemical potential vs the host to ~1e-12.

### WP-S6.2 — device populations resident (q_sh / dp_at / qp_at)
- B1 already gives resident `q_at`; add resident `q_sh` (shell reduction of
  `pop_ao`, the `n0_sh − n_sh` form) and keep `dp_at`/`qp_at` resident (the
  multipole-moments kernel already computes them — stop downloading). No `pop_ao`↓.

### WP-S6.3 — device SCC energy
- `E_coulomb = ½ q_shᵀ γ q_sh` (cuBLAS `Dsymv`+`Ddot` on resident `dGamma`/`q_sh`),
  shell third-order `Σ q_sh²·Γ_s` (reduction), multipole energy (the `amat`
  contraction, mirrors `energyMultipole`; the `amat` are already resident from
  Stage 5 B3), band `Σ occ·eps`. Emit one scalar `E_scc` ↓ for the host display +
  convergence. (D4/repulsion/halogen are post-SCF, unchanged.)

### WP-S6.4 — device Broyden mixer (the crux)
- Port `BroydenMixer` (modified Broyden / Johnson 1988, `broyden_mixer.h`): the SCC
  vector `[q_sh; dp_at; qp_at]` (length `nsh+9·nat`), the stored differences and the
  small `m×m` history dot-product system. The vector ops + the (tiny) history solve
  run on the device; the mixed vector stays resident as the next iteration's input.
  No `q_sh`/`dp_at`/`qp_at` ↑. **Validate** bit-faithfully vs the host mixer on a
  fixed input sequence (the history makes this the highest-risk item).

### WP-S6.5 — device convergence + the loop seam
- `max|Δq_sh|` (reduction) + `|ΔE|` → one scalar (or a packed 2-vector) ↓ per
  iteration; the host only reads it to decide break/continue. New backend entry
  `residentScfStep(...)` (or a device-driven inner loop with a host poll every k
  iterations) that fuses solve → occ → density → populations → energy → mix →
  convergence on the device. Keep the existing per-iteration host path as the
  fallback (gate exactly as Stages 2–5).

## 3. Validation
- `gpu_gfn{1,2}_validation` @1e-8 vs tblite stays green (xfails unchanged).
- Component: device occupation/energy/mixer each vs the host at a frozen state
  (~1e-12), mirroring the `sqm_cuda_*` pattern.
- Full `-sp`/`-opt` GPU energy bit-stable vs the host loop (the Stage-5 bar:
  complex `−329.52707823`, CH4 `-opt` `−4.17521844`).
- `compute-sanitizer` clean; no-CUDA `release/` unchanged (`#ifdef`-free host,
  inert `GpuScfBackend` defaults). Bench `scripts/sqm_bench.sh` (pass `-threads N`).

## 4. Effort / risk / sequencing
- **Effort:** high — five device kernels (Fermi, populations, energy, Broyden,
  convergence) + a fused step entry. **Risk:** high — it owns the SCF fixed point;
  the Broyden history port (S6.4) is the delicate part.
- **Sequencing:** S6.1→S6.2→S6.3 are independent and individually validatable;
  S6.4 (mixer) is the gate to actually removing the `q_sh`↑ traffic; S6.5 fuses.
- **Value:** the only Stage that targets the measured ~81 % latency region. Weigh
  against the eigensolve-compute path (already FP32-accelerated, compute is *not*
  the bottleneck). Do this only if the GPU SCF is being pushed for `-opt`/`-md`
  throughput; for `-sp` the one-time context init dominates a single point anyway.
