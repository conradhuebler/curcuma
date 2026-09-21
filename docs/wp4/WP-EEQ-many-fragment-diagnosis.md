# WP-EEQ — why the factor cache never engages on many-fragment systems

**Status:** ⚙️ diagnosis complete and measured; the fix is NOT implemented.
**Date:** 2026-07-26
**System:** `mixture.xyz` — 6200 atoms, **1400 fragments**, the only many-fragment
case in the test set.

This supersedes the hypothesis in
[WP-EEQ-P2-profile-plan.md](WP-EEQ-P2-profile-plan.md), which assumed the
per-step EEQ cost is the Cholesky factorisation and that loosening
`eeq_refactor_eps_bohr` would amortise it. The threshold is not the binding
constraint — see below.

## 1. The EEQ solve dominates the MD step

`prepareCNAndEEQ` is **82–91 % of per-step calculation time** at ~850 ms/step.
Everything else — force-field terms, integrator, I/O — shares the remaining
~10 %.

Single-point parameter generation, for orientation (after the A1c fix):

| phase | ms |
|---|---:|
| **EEQ Phase 1** | **3307** |
| Dispersion pairs | 1155 |
| **EEQ Phase 2** | **867** |
| HB/XB detection | 628 |
| CN + hybridisation + rings | 193 |
| Distance matrix | 90 |
| Topo distances + BATM | 65 |
| FT-HMO | 33 |
| Pi charges (ipis, 400 pi-systems) | 0.30 |

Phase 1 internals (`CURCUMA_EEQ_PROFILE` instrumentation, since removed):

| step | ms |
|---|---:|
| dxi + chi/gam/alpha | 0.2 |
| allocate A (7600², 462 MB) | 78 |
| Dijkstra topological distances | 85 |
| Coulomb fill (19.2M `erf`) | 226 |
| **dispatchSolve** | **2904** |

The build is cheap; **88 % is the solve**, and at ~1.9e11 flops in 2.9 s
(≈65 GFLOP/s) that solve runs at hardware speed. The cost is the *algorithm*:
SchurCholesky is O(N³/3 + N²·nfrag), and at nfrag=1400 the Schur term
(1.08e11) exceeds the Cholesky itself (7.9e10).

## 2. The alternative solvers are worse, not better

10-step MD, mixture, `-threads 16`:

| `solve_method` | Phase 1 | Phase 2 | total SP | energy |
|---|---:|---:|---:|---|
| **cholesky** (default) | 3290 ms | 869 ms | 12470 ms | −857.18483727 |
| pcg | 7342 ms | 13152 ms | 28784 ms | −857.18483859 |
| auto | 3301 ms | 886 ms | 12491 ms | −857.18483727 |
| batched | 488 ms | 320 ms | 9080 ms | **−856.10168551** |

- **PCG is 2–15× slower.** 1400 fragments in contact couple strongly and the
  block-Jacobi preconditioner does not capture that, so the iteration count
  explodes. The O(N²·k) argument fails because k is large.
- **Batched is 7–10× faster but wrong by 1.08 Eh (680 kcal/mol)** — it drops the
  cross-fragment Coulomb by construction, which is exactly why
  `eeq_contact_prefer_exact` routes contacting fragments away from it.

**The current default is the best of the three.** There is no better solver to
switch to; the win has to come from *not solving from scratch every step*.

## 3. The factor cache cannot engage — and the threshold is not why

Measured cache decisions over 10 MD steps (`-eeq_solver.verbosity 3`; note the
solver reads its own `verbosity`, the CLI flag does not reach it):

| `eeq_refactor_eps_bohr` | refine | cache hits | refactorisations |
|---|---|---:|---:|
| 0.05 (default) | 1 | **0** | 6 |
| 0.50 | 1 | 1 | — |
| 100.0 (geometry trigger effectively off) | 1 | **1** | — |

Even with the geometric trigger disabled the factor is rebuilt almost every
step. Wall time is flat across all settings (35–39 s, inside noise).

### The actual chain

1. `needsFullTopologyUpdate()` (`gfnff_method.cpp`) compares
   `max|Δr|` over **all 3N components** against a fixed **0.5 Bohr**. Over 18600
   components in a 300 K box, some atom exceeds that within a few steps, so it
   fires almost every step.
2. → topology recompute → **Phase 1 EEQ** runs.
3. → Phase 1 clears `m_pending_geometry` (`eeq_solver.cpp:2769`) to re-arm its
   documented contract.
4. → Phase 1's solve therefore takes the "local solve" branch, which sets
   **`m_chol_cache.valid = false`** (`eeq_solver.cpp:2189`).
5. → the next Phase-2 call has no factor and must refactorise.

**Three coupled invalidation mechanisms.** Gating only the explicit
`invalidateCholeskyCache()` in the topology-rebuild path (tried: 2 hits instead
of 0, wall time unchanged) does not help, because steps 3–4 invalidate anyway.

The size-scaling max-norm appears **twice** in this code path — here and in
`eeq_refactor_eps_bohr` — and in both places it means "always refactor" once the
system is large enough. It is the same defect the threading work hit in a
different guise: a criterion tuned on small molecules that silently inverts on
large ones.

## 4. Iterative refinement does work

The A4 refinement (`eeq_refine_iters`, default 1) was documented as having no
measurable benefit — that was measured on **polymer** (1410 atoms, 1 fragment),
where the EEQ is not the bottleneck. Here it matters:

| eps | refine | final Epot after 10 steps |
|---|---|---|
| 100.0 | 1 | −876.756919 |
| 100.0 | **0** | **−876.818559** |

Reusing a stale factor without refinement drifts **0.062 Eh = 39 kcal/mol in 10
steps**. With refinement the trajectory is stable. So the safety mechanism for
factor reuse is in place and validated on the system that needs it — only the
reuse itself never happens.

## 5. What a fix would have to do

Not implemented; listed so the next attempt starts from the measurement rather
than the hypothesis.

1. **Make the topology trigger size-aware.** A max over 3N components against a
   fixed bound cannot work across two orders of magnitude in N. Note the
   criterion is not *wrong* for its stated purpose (any single atom moving
   0.5 Bohr can change connectivity) — the problem is what it costs downstream.
2. **Decouple "topology was recomputed" from "the factor is worthless".** If the
   recomputed connectivity is unchanged, dxi/dgam are unchanged and the factor is
   as good as it was. Comparing the new bond/neighbour lists against the previous
   ones is cheap and exact, and would let the cache survive the common case.
3. **Let Phase 1 stop clobbering the Phase-2 cache.** Step 4 above is a side
   effect of Phase 1's local-solve branch, not a deliberate policy.

Expected payoff if the cache does engage: the cache-hit path is one triangular
solve plus the nfrag×nfrag Schur solve (~1e9 flops) against a full
~1.9e11-flop rebuild — roughly a 30× cut on 82–91 % of the MD step.

## 6. Related, already fixed

`A1c` (commit `754246d3`) removed a pathological O(N³)-per-fragment loop in the
ipis block: mixture went from **not finishing in 7 minutes to 12.3 s**. That was
a different bug in the same area and is not what this document is about.

## Reproduce

```bash
# per-step share and cache decisions (the solver has its OWN verbosity)
curcuma -md mixture.xyz -method gfnff -threads 16 -maxtime 1e1 -print 100 \
        -eeq_solver.verbosity 3 -verbosity 2 -no_bmt | grep -E "EEQ-Cache|prepareCNAndEEQ"

# solver comparison
for m in cholesky pcg auto batched; do
  curcuma -sp mixture.xyz -method gfnff -threads 16 \
          -eeq_solver.solve_method $m -verbosity 2 -no_bmt
done
```

> Measure on an idle machine — multi-thread wall times on this box move 20–30 %
> under an unrelated 4–5 core job.
