# RMSD-MTD bias speedup: Gaussian-cutoff screen + pool cap

Status: 🤖 AI-generated, ⚙️ machine-tested (Jul 2026). Not ✅ TESTED/APPROVED.

## Problem

The RMSD-metadynamics bias is `V(x) = Σ_i k·counter_i · Σ_p exp(−α·RMSD(x, x_{i,p})²)`
summed over every stored hill `x_i` (p = identity + symmetry images). Each MD step iterates
**every** hill and runs a full Kabsch/SVD alignment + analytic gradient for it
(`SimpleMD::ApplyRMSDMTD`, `src/capabilities/simplemd.cpp`). In ConfSearch the shared pool
grows across temperature cycles, so per-step cost is `O(N_hills · N_perms)` SVD alignments and
grows linearly with the pool — the search slows down exactly when it has found many structures.

## Stage 1 — Gaussian-cutoff screen (physics-preserving), default ON

Adapts the standard PLUMED "cutoff + neighbour list" idea to RMSD space, using a **rigorous RMSD
lower bound** so far hills are skipped before the expensive Kabsch, with zero effect on energy,
force, or the visited/deposition bookkeeping.

- **Lower bound (Mirsky's inequality).** For geometric-centered subset coordinate matrices X
  (walker) and Y (hill), `RMSD_min(X,Y) ≥ ‖σ(X) − σ(Y)‖₂ / √N_sub`, where σ are the three singular
  values of the centered `N_sub×3` matrix = principal radii of gyration (√ eigenvalues of the 3×3
  gyration tensor XᵀX). A rotation leaves σ unchanged, so this is a true lower bound. It is also
  **permutation-invariant**, so one check bounds all symmetry images of a hill at once.
- **Cutoff.** Skip a hill when `n_images · exp(−α·L²) < eps_step`, with `L` the lower bound and
  `eps_step = min(rmsd_mtd_cutoff_tol, global_count/rmsd_econv)`. Tying `eps_step` to the visited
  threshold keeps the visited/deposition gate exact. Because `L ≤ RMSD`, a skipped hill's true
  Gaussian is even smaller — no near hill is ever missed.
- **Caching + kernel micro-opt.** Each hill's centered subset + σ are cached (lazily, keyed by the
  stable `BiasStructure::index`, cleared each MD run). The walker is centered once per step. Survivors
  use `RMSDDriver::BestFitRMSDCentered()`, which skips the two per-call re-centerings of `BestFitRMSD`.
  This makes screen-ON cost ≤ legacy **even at 0 % screening**.
- Disabled automatically when RMSF weights are active (the unweighted σ-bound is invalid there).

### Effectiveness (honest scope)

The screen is **always correct**, but its *skip rate* depends on how much the pool spans in shape:
- The σ (radius-of-gyration) bound is **loose for conformers of similar overall shape** — two
  distinct conformers (RMSD ~1–3 Å) can have nearly identical radii of gyration, so the bound stays
  below the cutoff and the hill is not skipped.
- Measured: **triose (compact sugar, 44–56 hill pool) ≈ 0.1 % skipped**; butane ≈ 0 %. The screen
  pays off for pools with **large shape diversity** (compact ↔ extended, unfolding, long chains).
- For dense conformer pools, the reliable lever is the **pool cap** below; the screen still helps via
  the unconditional centered-kernel micro-opt and never hurts correctness.

### Parameters (registered in both `simplemd` and `confsearch`)

| flag | default | meaning |
|------|---------|---------|
| `-rmsd_mtd_screen` | `true` | enable the screen; `false` = evaluate every hill (legacy) |
| `-rmsd_mtd_cutoff_tol` | `1e-8` | Gaussian tolerance; smaller = more conservative |
| `-rmsd_mtd_screen_margin` | `0.0` | extra safety radius (RMSD units) added to the cutoff |

At verbosity ≥1 an RMSD-MTD run prints, once at the end,
`RMSD-MTD screen: <computed> hills computed, <screened> screened (<pct>% Kabsch fits skipped); MTD wall time <ms> ms`.

## Stage 2 — Pool cap (`rmsd_mtd_max_gaussians`), default off

Previously `rmsd_mtd_max_gaussians` was read but never enforced. It now caps the shared pool
**between temperature cycles** via `SharedBiasPool::capToSize()`: keep every persistent (fed-back
optimised) minimum plus the highest-counter non-persistent snapshots up to the cap, dropping the
rarely-visited rest. This bounds `N_hills` — and therefore the per-step bias cost — hard. It changes
the sampling slightly (a coarser bias), so it is opt-in (`-rmsd_mtd_max_gaussians N`, `-1` = unbounded).
This is the practical lever when the pool grows to "really many" structures.

## Verification (machine-tested)

- **Local path equivalence**: butane `-md -rmsd_mtd`, gfnff, screen off vs on, same seed →
  **bit-identical COLVAR over 3000 steps** (max|Δbias| = 0, max|Δrmsd_ref| = 0), identical pool (4=4).
- **Shared path equivalence**: butane ConfSearch, screen off vs on → identical accepted conformer
  (−1.959357 Eh) and identical bias-pool trace.
- **Pool cap**: butane ConfSearch `-rmsd_mtd_max_gaussians 3` → "Bias pool capped to 3 … dropped N".
- **Regression**: `ctest -R "cli_simplemd_|cli_confscan_|routing_03"` → 20/20 pass.

## Not tested / limitations

- Effectiveness on genuinely large shape-diverse pools (unfolding/long chains) is expected but not
  yet benchmarked; MTD wall-time was below the ms timer resolution on the small pools tried.
- The σ-bound is deliberately conservative; a tighter (still rigorous) bound would need O(N²) work
  (distance-matrix RMSD), which costs more than the Kabsch it would save — not pursued.
- Parallelising the survivor loop was **not** done: in ConfSearch the cores are already saturated by
  walker-level parallelism, so it yields ~0 there while adding thread-safety risk to a hot path.
