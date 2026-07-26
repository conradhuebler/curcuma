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

- **The final statistics report each PES against itself** (fixed Jul 2026). The
  accepted conformers are ranked on `opt_method`, so the closing "search lowered
  the energy by X kJ/mol" line uses the **`opt_method`** energy of the input
  structure. The `md_method` progression is reported separately as
  `exploration (<md_method>) lowered its own minimum by ... -- separate PES`,
  comparing the `md_method` initial energy to the `md_method` running best.
  Before the fix the `md_method` initial energy was subtracted from the
  `opt_method` global minimum, which produced absurd numbers
  (`gfnff: -9.005 -> -85.017 Eh`, "lowered by 199567 kJ/mol"). The `opt_method`
  reference energy is part of the checkpoint (`initial_energy_opt`), so a
  `-restart` resume keeps the ranking-PES comparison; checkpoints written before
  Jul 2026 simply omit that line.

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

## One topology reference per PES (fixed Jul 26, 2026)

The caveat that used to stand here — "the refinement side checks the `opt_method` geometries against
the `md_method` reference; a method that systematically pushes a bond length across the connectivity
cutoff could change topo-reject counts" — turned out to be not a caveat but a defect. Observed on a
real run:

```
ConfSearch: opt_method (gfn2) refinement: 0 structures -> cumulative + bias, 38 rejected (topo).
ConfSearch: gfn2 best: -161.600500 Eh (+0.00 kJ/mol vs. initial -161.600500 Eh)
```

**Every** re-optimised structure of the cycle was discarded while the exploration side kept 10 of the
same structures — the whole ranking side of that temperature was lost silently, and the final
ranking thinned out accordingly. The two methods simply disagree about one contact near the
covalent-radius cutoff; that is a difference in bond-length prediction, not a chemical reaction.

Now each PES is checked against its own reference:

| side | geometries | reference |
|---|---|---|
| exploration | `md_method` minima | input structure optimised at `md_method` |
| refinement | `opt_method` minima | input structure optimised at `opt_method` |

Both references are taken at the same point of the run (the initial optimisation), both are stored in
the restart checkpoint (`topo_ref`, `topo_ref_opt`), and a resume from an older checkpoint falls back
to the single reference. Three diagnostics were added because the failure was silent:

- a warning when the two methods disagree about the input structure's topology at all,
- the first differing atom pair on the refinement side (the exploration side already had this),
- a loud warning when *all* structures of a cycle are rejected on either side.

## Other caveats (not yet human-tested)

- Next-cycle MD seeds are the `md_method`-optimized geometries (exploration stays
  on the cheap PES); the next MD also runs at `md_method`.
- ORCA/r2scan as `opt_method` across many temperature cycles is feasible only
  because the per-cycle set is filtered first; still, total cost scales with the
  number of cycles. Use a short temperature schedule for expensive `opt_method`.
- The `md_method` and `opt_method` minima are usually in the same basin, so the
  two deposited bias geometries sit close together (≈ one slightly taller hill).
