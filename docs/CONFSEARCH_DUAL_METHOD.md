# ConfSearch Dual-Method Workflow (`md_method` / `opt_method`)

> 🤖 AI-generated, machine-tested (Jun 2026). Human production testing pending.

ConfSearch can use a **cheap method for exploration** and a **more accurate
method for refinement/ranking**. Example: explore with `gfnff`, then optimize
and rank the discovered conformers with `gfn2` (or an ORCA method).

## Parameters

| Flag | Default | Role |
|------|---------|------|
| `-method` | `gfnff` | Single method used everywhere (legacy/back-compat). |
| `-md_method` | *(empty -> `method`)* | MD exploration **and** the pre-optimization. |
| `-opt_method` | *(empty -> `method`)* | Accurate per-cycle re-optimization **and** final ranking. |

Both new flags fall back to `-method` when left empty, so an existing
`-method gfnff` run is **byte-for-byte unchanged**. Flat (`-md_method ...`) and
dotted (`-confsearch.md_method ...`) forms both work.

> **Fixed July 2026.** Until ConfSearch was migrated into the ParameterRegistry,
> `md_method`/`opt_method` were registered only by the **polymerbuild** module, so
> the CLI auto-router silently moved both flags to `controller["polymerbuild"]` and
> ConfSearch never saw them — it fell back to `method` and logged
> `Single-method mode`. The dotted form did not help either (the `confsearch.`
> prefix is stripped before the router runs), which left `-import_config` as the
> only working route. ConfSearch now owns these names and both forms work.

## Per-cycle pipeline

```
MD exploration            (md_method)
  -> Phase 2  fast opt          (md_method)
  -> Phase 3  RMSD/energy filter (md_method)        <- "filter between"
  -> Phase 3b accurate re-opt    (opt_method) [NEW] <- only if opt_method != md_method
  -> Phase 4  EXPLORATION (md_method): topo + seed select + global min + bias
              REFINEMENT  (opt_method): cumulative pool + bias    [dual only]
...repeat over the temperature schedule...
Final deduplication / ranking   (opt_method)
```

The md-level optimize+filter (Phases 2-3) reduces the per-cycle set before the
accurate method runs, so `opt_method` (incl. expensive ORCA) only sees the
deduplicated survivors.

### Two PES, never compared

`md_method` and `opt_method` live on **different potential energy surfaces**.
ConfSearch keeps the two energy worlds strictly separate:

- **Exploration decisions stay on `md_method`.** Seed selection for the next
  cycle, the running exploration global minimum, the per-cycle "new best"
  progression, and the seed energy window all read the **md_method** minima. A
  basin discovered by `md_method` is **never discarded because `opt_method`
  ranks it higher** — e.g. a gfnff-stable conformer that `r2scan` finds less
  stable is still carried forward.
- **`opt_method` energies only feed the refinement/ranking.** The re-optimized
  structures fill the cumulative pool (its window is relative to that cycle's
  lowest **opt_method** energy) and the final ranking. They are never subtracted
  from / compared to `md_method` energies.
- **Both optimized geometries enter the bias pool.** The `md_method` minimum
  (which drives the next `md_method` MD) and the `opt_method` minimum are both
  deposited as persistent bias structures. This is safe because the RMSD-MTD
  bias force is purely geometry-based (`W = k·counter`); the stored `energy` is
  metadata only and never enters the force, so no cross-PES comparison occurs.

**Skip rule:** when `opt_method == md_method` (the default single-method case)
Phase 3b and the separate refinement step are skipped entirely — Phase 4 runs
its single-PES path and one minimum per conformer is fed back, exactly as before.

## Example

```bash
curcuma -confsearch input.xyz -md_method gfnff -opt_method gfn2
```

Per cycle this produces `*.bias.opt.accepted.xyz` (md-level, filtered; used for
exploration/seeds/bias) and, when the methods differ,
`*.bias.opt.accepted.opt.xyz` (opt-level, re-optimized; used for the cumulative
pool + bias). The final `*.cumulative.opt.accepted.xyz` is ranked at `opt_method`.

## Caveats (not yet human-tested)

- The topology reference is set from the `md_method` pre-optimized geometry. The
  exploration side checks `md_method` geometries against it; the refinement side
  checks the `opt_method` geometries against the same reference (0/1 bond
  connectivity). A method that systematically pushes a bond length across the
  connectivity cutoff could change topo-reject counts on the refinement side.
- Next-cycle MD seeds are the `md_method`-optimized geometries (exploration stays
  on the cheap PES); the next MD also runs at `md_method`.
- ORCA/r2scan as `opt_method` across many temperature cycles is feasible only
  because the per-cycle set is filtered first; still, total cost scales with the
  number of cycles. Use a short temperature schedule for expensive `opt_method`.
- The `md_method` and `opt_method` minima are usually in the same basin, so the
  two deposited bias geometries sit close together (≈ one slightly taller hill).
