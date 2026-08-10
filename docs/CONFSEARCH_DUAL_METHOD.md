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

One temperature stage, repeated `repeat` times (Aug 2026):

```
  Step 1 EXPLORE    MD with RMSD bias, one run per seed   (md_method)
  Step 2 RELAX      geometry optimisation                 (md_method)
  Step 3 REDUCE     RMSD deduplication                    (md_method)
  Step 4 RECOMBINE  ConfGen proposals, optional           (md_method)
  -- re-seed from the cross-cycle pool, next repetition --
  Step 5 REFINE     accurate re-optimisation              (opt_method)  once per temperature
  Step 6 SELECT     energy window + topology + seeding
...next temperature...
Final deduplication / ranking                             (opt_method)
```

Steps 1-4 are chained: each repetition re-seeds from the best structures of ALL
cycles, so repetition 2 starts from what repetition 1 found. Every repetition adds
its minima to a per-stage pool; that pool is deduplicated once and is what step 5
receives, so the accurate method sees everything the stage produced and not just
the last batch. Intermediate repetitions seed on the md_method surface even under
`-seed_pes opt`, because no opt_method structures exist before step 5 has run.
Step 5 runs once per temperature, at the end of the last repetition -- it accounts for ~61 % of the
total runtime, so chaining it would roughly triple the cost of `repeat 4` without
adding a single trajectory. Steps 2-3 reduce the per-cycle set before the accurate
method runs, so `opt_method` (incl. expensive ORCA) only sees the deduplicated
survivors.

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

### How far the two surfaces agree — measured, and it is worse than "weakly" (Aug 2026)

> 🤖 AI-measured on ONE system (WEKLQ, a 107-atom peptide, gfnff/gfn2). Treat the numbers as a
> property of that system, not as a general statement about the two methods.

Across a whole 141-structure ensemble the two rankings correlate at r = 0.40, which is the number
this document carried until now. **Within one cycle — which is the situation that actually decides
anything — they are anti-correlated:**

| quantity | cycle A (13 structures) | cycle B (14 structures) |
|---|---|---|
| Pearson r (gfnff vs gfn2, optimised on both) | **−0.32** | **−0.46** |
| Spearman rho | −0.13 | −0.16 |
| gfn2 rank of the deepest gfnff structure | 9 of 13 | 10 of 14 (+68 kJ/mol) |

The structures are chemically valid (0 topology changes), so this is not an artefact of broken
geometries. The practical consequence was measured directly: a run that reached **121 kJ/mol deeper
on gfnff** produced a **67 kJ/mol worse** gfn2 result. Exploring deeper on the cheap surface does not
make the accurate answer better — which is the reason the search now selects on the ranking surface
before anything is deduplicated (`-phase3b_two_stage`, default since Aug 2026).

**Would a single point be enough?** No. A gfn2 single point on the gfnff geometry ranks better than
the gfnff energy (r = +0.40, rho = +0.35, best pick at rank 2 of 14) but still not well, because the
relaxation from the gfnff geometry to the gfn2 minimum spans **149–368 kJ/mol** — several times the
conformer differences themselves (~80 kJ/mol). The term that has to be removed is exactly the one a
single point leaves in. Only a crude optimisation on the ranking surface removes it, which is what
the two-stage mode does.

## Example

```bash
curcuma -confsearch input.xyz -md_method gfnff -opt_method gfn2
```

Per cycle this produces `*.<cycle>[_rN].s3_reduce.<md_method>.xyz` (md-level,
filtered; used for exploration/seeds/bias) and, when the methods differ,
`*.<cycle>.s5_refine.<opt_method>.opt.xyz` (opt-level, re-optimised; used for the
cumulative pool + bias). Every file carries the step that produced it. The final
`*.cumulative.opt.accepted.xyz` is ranked at `opt_method` and keeps its name -- it
is the deliverable, not an intermediate.

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

## When the two surfaces pull in different directions (`-surface_feedback`, Aug 2026)

A dual-method run assumes the cheap surface is a rough version of the accurate one. Measured on the
107-atom peptide it is not: within a cycle the two rank the same structures at **r = −0.32 … −0.46**,
and the disagreement is not noise but has a direction. GFN-FF rewards every additional hydrogen bond
(its HB term carries the largest share of the energy spread), so the cold cycles drifted into
10–11-bridge folds — GFN-FF put them **230 kJ/mol below everything else it had found**, GFN2 put the
same structures **52–95 kJ/mol above its own best**, and it removed 3–4 of their bridges during
optimisation (10–11 → 6–7). The exploration was walking away from the answer, on purpose, guided by
its own energy.

`-surface_feedback true` measures this instead of assuming it. After each cycle every ranked
structure also gets a single point on the *exploration* surface, so both energies exist for the same
geometry, and

```
delta_i = (E_rank,i - min E_rank) - (E_explore,i - min E_explore)
```

is positive exactly for the structures the cheap method over-rewards. `delta` is then correlated
with four **generic** shape descriptors — hydrogen bonds, close contacts, radius of gyration,
buriedness of the polar atoms — and the strongest one, if it passes `-surface_feedback_min_r` (0.4),
decides which structures get a `-surface_feedback_strength` (3×, clamped to 10) taller hill in the
shared bias pool. The exploration then leaves that direction instead of being drawn back into it.

The coordinate is deliberately **not** hard-coded. "Too many hydrogen bonds" is true for this peptide
(whose reference conformer has 7) and false for a polyol host; a fixed threshold would be an
algorithm for one molecule. Validated on the measured case: over 60 ranked structures of one cycle
the detector gives `r(delta, hydrogen bonds) = +0.44` while the other three descriptors sit at
−0.05, −0.02 and −0.03 — it finds the right coordinate without being told, and stays silent when
nothing correlates.

Requirements and limits: dual-method only (`md_method != opt_method`), at least
`-surface_feedback_min_structures` (20) structures with both energies, so it starts working from the
second cycle. The extra cost is one cheap single point per ranked structure per cycle.

### Measured end to end — it works, and it does not help (Aug 10, 2026)

Two otherwise identical seven-stage runs on the 107-atom peptide, one with the feedback, one
without. It does what it is built for: the share of over-bridged structures on the exploration side
falls from **10.6 % to 5.6 %** on average from cycle 2 on — six cycles in a row, with the two ranges
not overlapping (8.7–13.3 % against 2.6–6.8 %), while cycle 1, where the feedback is still inactive,
is identical in both (1.3 %).

The result does not move:

| | baseline | with feedback |
|---|---|---|
| structures | 1060 | 813 |
| deepest vs the reference | **+28.4** | **+28.6** kJ/mol |
| smallest RMSD | 2.40 A | **2.25 A** |
| smallest H-bond Hamming | **4** | 5 |
| smallest torsion distance | 7 of 29 | **6 of 29** |

Two numbers marginally better, two marginally worse, the energy equal to 0.2 kJ/mol. The reason is
in the same data: on the RANKING surface the over-bridging was never a mass phenomenon (0.4–3.7 %
per cycle) because the accurate re-optimisation discards or repairs those structures anyway — it
removes 3–4 of the excess bridges by itself. The feedback therefore saves time in a dead end; it
does not supply what is missing, namely a move that crosses the 6–7 torsions to the reference in one
step.

**Therefore off by default**, and worth switching on only when the exploration cost in that dead end
matters (it is one extra single point per ranked structure per cycle). The general lesson is worth
more than the feature: an intervention that demonstrably improves its own target metric does not
thereby improve the result — measuring only the target metric would have called this a success.

For a system whose behaviour is already known there is the manual alternative `-hbond_excess_max N`
(off by default, and deliberately without a sensible default): a structure with more than
`reference + N` hydrogen bonds is not used as a seed for the next cycle, and with
`-hbond_excess_reject true` it also stays out of the ensemble. The soft form is the recommended one —
those structures are real minima, they simply are not worth exploring further.

## Other caveats (not yet human-tested)

- Next-cycle MD seeds are the `md_method`-optimized geometries (exploration stays
  on the cheap PES); the next MD also runs at `md_method`.
- ORCA/r2scan as `opt_method` across many temperature cycles is feasible only
  because the per-cycle set is filtered first; still, total cost scales with the
  number of cycles. Use a short temperature schedule for expensive `opt_method`.
- The `md_method` and `opt_method` minima are usually in the same basin, so the
  two deposited bias geometries sit close together (≈ one slightly taller hill).
