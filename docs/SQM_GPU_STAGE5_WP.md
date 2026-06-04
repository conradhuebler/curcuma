# Work Package: Stage 5 — GFN2 in-SCF D4/EEQ potential on the GPU

> 🤖 AI-prepared work package (2026-06-04), not yet implemented. Companion to
> [docs/SQM_GPU.md](SQM_GPU.md), [docs/SQM_PERF_OPT_WP.md](SQM_PERF_OPT_WP.md)
> (AP3/AP4) and [docs/SQM_GPU_STAGE3_WP.md](SQM_GPU_STAGE3_WP.md). Stages 2b/3/4
> made the GFN2 SCF + integrals + nuclear gradient device-resident; AP4
> (`27a39fe`) removed the redundant host integral build. This WP specifies the
> remaining host work in the SCF loop: the **per-iteration isotropic potential
> build** (the AP3 target, expanded to include the D4 piece the user identified).
> Per CLAUDE.md AI-policy: every item below is a proposal; nothing here is
> implemented or human-tested. Mark new code "Claude Generated"; do **not**
> self-assign ✅ TESTED.

## 0. Where we are / what Stage 5 is

After Stages 2b/3/4 + AP4, a GFN2 `-sp`/`-opt` GPU step keeps H0/S/L/γ, the
density P, MO coefficients C, and the multipole integrals resident; the
eigensolve, density, populations, multipole moments and the nuclear gradient run
on the device. **The one thing still computed on the host every SCF iteration is
the isotropic AO potential** `v_ao`, assembled in `XTB::Calculation` and uploaded
to `residentSolve`:

```
v_ao(μ) = v_sh(ao2sh[μ]) + v_at(ao2at[μ])
  v_sh  = γ·q_sh                 (addCoulombShellPotential)   — host gemv
        + q_sh²·Γ_sh             (addThirdOrderPotential, GFN2 shell-resolved)
  v_at += dE_D4/dq_A             (addDispersionPotential)      — host D4 pair loop
        + multipole scalar shift
```

Stage 5 moves this build onto the device so the SCF iteration crosses the bus
only for the values that genuinely originate on the host (the Broyden-mixed shell
charges), removing the per-iteration `v_ao` round-trip.

### Honest value assessment (read before prioritising)

Measured on complex/231, GFN2, `-gpu cuda`, `-verbosity 3` (this session, after
AP2/AP4):

| SCF phase | time | per-it | note |
|---|---|---|---|
| **solve eigen (+ round-trips)** | **275.7 ms** | 14.5 | **81 % of SCF — latency-bound** |
| potential build | 46.4 ms | 2.44 | **of which D4 33.6 ms (1.77/it)**; non-D4 ~13 ms |
| populations | 8.0 ms | 0.42 | |
| energy/mix | 8.8 ms | 0.46 | |

- **The eigensolve/round-trip (276 ms) is the dominant cost and Stage 5 does NOT
  touch it.** Measured this session: FP32≈FP64≈partial-diag all ~330 ms SCF, i.e.
  the SCF is **latency-bound by the per-iteration `eps`/`pop_ao` download+sync
  round-trips**, not by eigensolve compute. The larger GPU lever is a separate WP
  (device-side occupation → no per-iter `eps`/`pop` round-trip).
- **Stage 5 targets ~46 ms/run** (the whole host potential build), of which the
  D4 piece (~34 ms) is the real content — the non-D4 AP3 part is only ~13 ms and
  the γ·q_sh gemv itself is sub-millisecond.
- **The D4 part is genuinely per-iteration and NOT hoistable:**
  `addDispersionPotential` feeds the current **SCF Mulliken charges**
  (`m_wfn.q_at`, via `setTopologyCharges`) into the per-reference C6 scaling each
  iteration, so `dE_D4/dq` changes every SCF step. Moving it off the host means
  computing it on the device, not caching it.
- **Conclusion:** Stage 5 is a real ~46 ms/it-region win **on top of** banking
  the eigensolve latency separately. It is a multi-hour CUDA port with a new
  kernel; do it only if the GPU SCF is being pushed hard. AP2 (committed, ~10 %)
  and AP4 (committed, ~16–20 ms/step) are the cheap wins already in hand.

## 1. Precedent (what already exists on the GPU)

- **EEQ on the GPU:** `ff_methods/cuda/eeq_solver_gpu.{cu,h}` (GFN-FF path).
- **D4 dispersion on the GPU:** `ff_methods/cuda/gfnff_kernels.cu::k_dispersion`
  (energy + gradient) — BUT it consumes GFN-FF's **pre-baked single-`zetac6`**
  pair C6, not the GFN2 **per-reference charge-weighted** C6. So Stage 5 needs a
  **new** kernel, not a reuse of `k_dispersion`.
- **D4 reference data is already accessor-exposed for GPU upload:**
  `D4ParameterGenerator::getC6FlatCache()` (118·118·7·7 flat C6 block),
  `getRefN()`, `getRefCN()`, `getGaussianWeights()` — geometry-fixed, upload once.
- **Device-resident GFN2 scaffolding:** `XtbGpuContext::residentSolve` /
  `residentSolveMultipole`, `dGamma` resident, `dPop`/density on device.

## 2. Work items

### WP-S5.1 — device `q_at` (atomic SCF charges on the device)
- **Why:** the D4 potential needs the per-atom SCF charges `q_at`; they are
  host-only today (`updatePopulationsFromPopAo` runs on the host from the
  downloaded `pop_ao`).
- **Approach:** the device already computes `pop_ao` (`residentDensity`,
  `k_pop_ao`). Add a device reduction `q_at(A) = z_A − Σ_{μ∈A} pop_ao(μ)` (the
  atom part of `updatePopulations`), keeping `q_at`/`q_sh` resident. Needs the
  AO→atom / AO→shell maps (already uploaded for the multipole path: `dAo2at`).
- **Caveat:** Broyden mixes the **shell** SCC vector on the host, so the
  *mixed-input* `q_sh` for the next iteration still originates on the host and is
  uploaded (length-nsh). Stage 5 uploads `q_sh` instead of `v_ao` (≈ same size:
  nsh+nat ≈ nao — **no net transfer win**, the win is removing host compute).

### WP-S5.2 — device GFN2 D4 potential kernel (the hard part)
- **Why:** the 33.6 ms/run in-SCF `dE_D4/dq`.
- **Port:** the AP6b exact per-reference path —
  `D4ParameterGenerator::buildAtomRefW` (per-atom CN-Gaussian × charge-zeta
  weights, `d4_zeta`) + `contractC6Gfn2` (7×7 block over `m_c6_flat_cache`) +
  `D4Evaluator::computeEnergyAndGradient`'s BJ-damping pair loop + `computeATM`,
  emitting `dEdq(A)`. Inputs (upload once per geometry): C6 flat cache, refN,
  refCN/refCovCN, refq, the geometry-fixed pair list / CN. Per iteration: the
  device `q_at` (WP-S5.1) drives the per-reference zeta scaling.
- **Reuse the AP2 split:** `buildAtomRefW`/`contractC6Gfn2` are the device-kernel
  blueprint (per-atom weights once, then the 7×7 contraction) — the same hoist
  shape that AP2 introduced on the CPU.
- **Scope notes:** energy-coupling only needs `dE/dq` (not the Cartesian
  gradient — that stays the Stage-4 device gradient, which already has its own D4
  path). Mind the ATM `s9=5` term. `d4_charge_source=mulliken` is the default-off
  CPSCF case — Stage 5 covers the in-SCF coupling which uses `q_at` regardless.

### WP-S5.3 — fold AP3 (γ·q_sh + third-order) into the device build
- Once `q_sh`/`v_at` are device-resident: `v_sh = γ·q_sh` (cuBLAS `Dgemv` on the
  resident `dGamma`) + the GFN2 shell third-order `q_sh²·Γ_sh`; expand
  `v_ao(μ)=v_sh(ao2sh)+v_at(ao2at)` in a kernel. This is the literal AP3, now
  cheap because the data is already on the device.

### WP-S5.4 — `residentSolve` variant + remove the host round-trip
- New entry `residentSolvePotential(q_sh, …)` that builds `v_ao` on the device
  (WP-S5.2 + S5.3) and folds straight into the existing Fock build + eigensolve.
- The host loop uploads only the Broyden-mixed `q_sh` (length-nsh) and the
  multipole scalar inputs; no `v_ao` upload. The multipole anisotropic potential
  (`v_dp`/`v_qp`) path is unchanged.

## 3. Validation
- **Energy:** `gpu_gfn2_validation` @1e-8 vs tblite stays green (xfail `complex`);
  complex GFN2 `-sp` energy bit-stable vs the host-potential path.
- **Component:** new device `dEdq(A)` matches the CPU `addDispersionPotential`
  output elementwise (~1e-12) at a frozen set of `q_at` (mirror the
  `sqm_cuda_*` integral component tests).
- **Multi-step:** 2–3-cycle `-opt` GPU-vs-CPU final energy (a full `-opt` is too
  slow to iterate on), as for AP4.
- **No-CUDA `release/`:** unchanged (host changes `#ifdef`-free).
- **`compute-sanitizer`** clean on the new kernel. Bench: `scripts/sqm_bench.sh`.

## 4. Effort / risk
- **Effort:** high (a new GFN2-D4 device kernel + device `q_at` + a resident-solve
  variant). **Risk:** medium-high — it lands in the validated per-iteration hot
  path (currently 90/90 GPU tests). Gate behind the device-resident path with a
  host fallback, exactly as Stages 2b/3/4.
- **Sequencing:** independent of the eigensolve-latency WP; that one (device-side
  occupation, no per-iter `eps`/`pop` round-trip) is the larger lever and should
  be weighed first if only one is pursued.
