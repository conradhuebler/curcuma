# ConfSearch: reference energies, per-cycle output, two-stage refinement

Status: 🤖 AI-generated, ⚙️ machine-tested (Jul 2026). Not ✅ TESTED / ✅ APPROVED.

What a conformer search reports, what it writes and how its files are named — plus two bugs
these changes uncovered (a non-re-entrant ripser, and optimisation batches that were silent).

## 1. The reference energy is the optimised INPUT structure

The search reports how much it gained. That number is only meaningful against a fixed reference,
and the reference has to be measured *before* any MD ran.

It used to be taken from the lowest structure of the **first metadynamics cycle** — after Phase 1-3
of cycle 1. Everything that cycle 1 found was therefore absorbed into the baseline and reported as
zero gain. The reference is now taken directly after the initial geometry optimisation, on each
level of theory separately:

```
ConfSearch: Reference energies (optimised input structure, one per level of theory):
ConfSearch:   gfnff      (exploration) = -1.958615 Eh
ConfSearch:   gfn2       (ranking)     = -13.664177 Eh
ConfSearch:   the two are different potential-energy surfaces -- no difference between them is reported
```

The previous line printed a **difference between the two methods**:

```
ConfSearch: Initial energies: gfnff=-18.739260 Eh, gfn2=-161.600419 Eh (delta=375081.97 kJ/mol)
```

375081.97 kJ/mol is not a physical quantity — it is the different zero of energy of the two methods.
It is gone; every energy difference in the output now compares one method against itself.

If the initial optimisation produces no usable energy at all, the first cycle's lowest structure
becomes the reference, with a warning. That is the only remaining path to the old behaviour.

## 2. Every phase reports the energies of the ensemble it produced

A structure count says nothing about what a phase produced. Each phase now reports the ensemble it
wrote, read back from the file:

```
ConfSearch: RMSD filtering complete. 36 structures accepted.
ConfSearch: REDUCE -- accepted ensemble [gfnff]: 36 structure(s), lowest -18.742131 Eh, span 41.28 kJ/mol
ConfSearch:     #1   -18.742131 Eh  (+0.00 kJ/mol)
ConfSearch:     #2   -18.741062 Eh  (+2.81 kJ/mol)
ConfSearch:     #3   -18.74055  Eh  (+4.15 kJ/mol)
ConfSearch:     ... and 33 more, up to +41.28 kJ/mol
```

`-ensemble_report N` sets how many conformers are listed (0 = summary line only).

## 3. Every file says which cycle, which purpose, which method

The old names described neither the stage nor the method that produced the file, and the per-cycle
working files overwrote each other, so what a given temperature did was gone as soon as the next one
started. `input.input.opt.xyz` also just read like a typo.

```
<base>.<purpose>.<method>[.opt][.accepted].xyz               (before the temperature loop)
<base>.cycleNN_TxxxK.<purpose>.<method>[.opt][.accepted].xyz (inside a cycle)
```

The cycle tag comes first so a directory listing groups a cycle together:

| file | what it is |
|---|---|
| `input.initial.gfnff.xyz` / `.opt.xyz` | input structures, before / after the pre-optimisation |
| `input.initial.gfn2.xyz` / `.opt.xyz` | the same at the ranking level (dual-method) |
| `input.cycle01_T500K.explore.gfnff.xyz` | MD snapshots of that cycle |
| `input.cycle01_T500K.explore.gfnff.opt.xyz` | after Phase 2 |
| `input.cycle01_T500K.explore.gfnff.opt.accepted.xyz` | after the Phase 3 dedup |
| `input.cycle01_T500K.bias_pool.gfnff.xyz` | the full RMSD-MTD bias pool (was `.mtd.xyz`) |
| `input.cycle01_T500K.refine_crude.gfn2.*` | Phase 3b stage 1 + its dedup (two-stage mode only) |
| `input.cycle01_T500K.refine.gfn2.xyz` / `.opt.xyz` | Phase 3b accurate stage |
| `input.cycle01_T500K.ensemble.gfnff.xyz` | the cycle's result on the exploration level |
| `input.cycle01_T500K.ensemble.gfn2.xyz` | the cycle's result on the ranking level |
| `input.best_per_cycle.<method>.xyz` | one frame per cycle: its most stable structure |
| `input.cumulative.opt.xyz` / `.cumulative.opt.accepted.xyz` | **unchanged** — the final deliverable keeps its established name |

Two consequences worth knowing:

- Each Phase-3b stage starts from a **copy** of its input rather than chaining another suffix onto
  the previous file. One extra write per cycle buys a name that states the method — the old chain
  ended in `<base>.bias.opt.accepted.opt.accepted.opt.xyz`.
- Nothing is overwritten between cycles any more, so a long run keeps every intermediate. That is
  deliberate (it is the run's provenance, and BMT gives each run its own directory), but it does
  mean a 6-cycle run on a large system writes ~6× more intermediate XYZ than before.

## 4. Per-cycle ensembles and most stable structures, on both levels

Besides the stage files of section 3, every cycle writes its RESULT — the topology-valid structures
that survived Phase 4 — for **each level of theory**:

| File | Content |
|---|---|
| `<base>.cycleNN_TxxxK.ensemble.<method>.xyz` | the cycle's topology-valid ensemble, **energy-sorted** |
| `<base>.best_per_cycle.<method>.xyz` | one frame per cycle: its most stable structure (a trajectory of the search's progress) |

```
ConfSearch: cycle01_T500K ensemble [gfnff]: 36 structure(s), most stable -18.742131 Eh, span 41.28 kJ/mol -> <bmt>/input.cycle01_T500K.ensemble.gfnff.xyz
ConfSearch:   most stable structure of this cycle appended to <bmt>/input.best_per_cycle.gfnff.xyz
```

Switch off with `-cycle_output false`. Ordering is done by a `std::multimap` keyed on the energy —
the container sorts, there is no explicit sort call (see section 6).

## 5. Two-stage high-level re-optimisation (`-phase3b_two_stage`, **default since Aug 2026**)

Step 5 (REFINE) re-optimises the cycle's deduplicated minima at `opt_method` and is the most
expensive step of a dual-method run. Its input was deduplicated on a **different** potential-energy
surface (`md_method`), so structures that survive as distinct there can collapse onto one another as
soon as the accurate method relaxes them — and each of those collapses was paid for at full price.

**Why it became the default.** Deduplicating and selecting on the exploration surface was measured
to be actively wrong for this class of system, not merely suboptimal: on a 107-atom peptide the
gfnff and gfn2 energies *within one cycle* correlate at **r = −0.32 / −0.46** (rank −0.13 / −0.16),
i.e. the structure gfnff likes best sits at gfn2 rank 9–10 of 13–14, +68 kJ/mol above the gfn2
minimum. A gfn2 single point on the gfnff geometry is better (r = +0.40, best pick at rank 2 of 14)
but still weak, because the relaxation from the gfnff geometry to the gfn2 minimum spans
**149–368 kJ/mol** — several times the conformer differences themselves. Only a crude optimisation
on the ranking surface removes that term, and that is what stage 1 does. Details and the full table:
[CONFSEARCH_DUAL_METHOD.md](CONFSEARCH_DUAL_METHOD.md). Set `-phase3b_two_stage false` for the old
single-stage behaviour.

```
-phase3b_two_stage true:
  stage 1  crude optimisation of EVERY structure   (preset phase3b_preopt_preset, default 'loose')
  filter   ConfScan dedup at that level            (phase3b_filter, default true)
  stage 2  accurate optimisation of the survivors  (preset phase3b_preset, default 'normal')
```

```
ConfSearch: Phase 3b stage 1/2 -- crude pre-optimisation at gfn2 (preset 'loose')
ConfSearch: Phase 3b stage 1 (crude) [gfn2]: 3 structure(s), lowest -13.665128 Eh, span 276.08 kJ/mol
ConfSearch: Phase 3b filter at gfn2: 2 of 3 crude structures are distinct minima (1 accurate optimisation(s) saved)
ConfSearch: Phase 3b stage 2/2 -- accurate optimisation at gfn2 (preset 'normal')
```

The crude structures are **thrown away** — only the selection they produce is used, so the accuracy
of the final ensemble is that of the accurate stage alone.

**The intermediate filter deduplicates, it does not rank.** Its energy cutoff is explicitly disabled
(`max_energy = -1`): a structure that sits high after a loose relaxation can still fall inside the
window once optimised properly, and discarding it on the crude number would lose a conformer for
good. The energy window is applied in Phase 4, on fully optimised energies. (Measured on the butane
smoke case: with the cutoff active the filter kept 1 of 3, without it 2 of 3.)

| Flag | Default | Meaning |
|------|---------|---------|
| `-phase3b_two_stage` | `true` | run crude -> filter -> accurate instead of one accurate optimisation |
| `-phase3b_preopt_preset` | `loose` | convergence preset of the crude stage |
| `-phase3b_preopt_max_iter` | `0` | step cap for the crude stage; 0 keeps the preset value |
| `-phase3b_preset` | `normal` | convergence preset of the accurate stage (also of single-stage Phase 3b) |
| `-phase3b_filter` | `true` | run the dedup between the stages |

How much it pays depends on how much the two surfaces disagree; every run states which way the flag
is set. It only exists in dual-method runs (step 5 is skipped when `opt_method == md_method`).

## 6. ConfScan sorting: verified, no `std::sort` needed

Claim under test: *ConfScan should sort automatically; the container should do it.*

**It does.** `m_ordered_list` is a `std::multimap<double,int>` keyed on the energy (`openFile`), and
every pass inherits that order. Verified by deliberately shuffling a 44-conformer ensemble and
running the scan: `*.accepted.xyz` comes out strictly ascending in energy, and the accepted set is
identical to the unshuffled run (14 structures).

An earlier comment in `confscan.cpp` claimed the opposite ("ConfScan sorts nowhere: structures are
processed in file order"). That was wrong and is corrected in place. Two spots now re-establish the
ordering locally instead of inheriting it, so the property is checkable where it is relied on
(`m_lowest_energy` is taken from the first structure of a pass; the energy cutoff and `maxrank` are
applied while walking the list):

- `ConfScan::Reorder` orders its input list through a multimap,
- the write-out of `*.accepted.xyz` does the same.

Both are no-ops on correctly ordered input and cost one multimap build per pass.

## The bug this uncovered: ripser is not re-entrant across instances

Adding a second ConfScan pass inside one temperature cycle segfaulted immediately in
`ripser::compute_dim_0_pairs`. Cause, in the fetched `external/ripser/ripser.h`:

```cpp
static simplex_boundary_enumerator facets(0, *this);
static simplex_coboundary_enumerator cofacets(*this);
```

Both enumerators store `const ripser& parent` and `const binomial_coeff_table&`. A **function-local
static** is constructed on the first call and never again — so from the second ripser instance in a
process onwards, those references point at a **destroyed object**. curcuma builds one ripser
instance per structure (ConfScan's persistent-diagram descriptor), so this is hit constantly. It
normally goes unnoticed because a fresh instance often lands on the same heap address; it crashes as
soon as the allocation pattern changes between two ConfScan runs.

Removing `static` fixes it (verified: 0 crashes where it previously crashed every time). Since
`external/ripser` is FetchContent-populated (and git-ignored), CMake applies the edit at configure
time; the step is idempotent and prints

```
-- ripser: removed function-local static simplex enumerators (dangling reference across instances)
```

**The proper fix belongs upstream in https://github.com/conradhuebler/ripser.** Once it is pushed
there and the `GIT_TAG` in `CMakeLists.txt` is bumped, the configure step silently does nothing.

Note the second consequence: before this fix, every persistent-diagram descriptor computed after the
first ripser object was destroyed read `parent.n` and the binomial table from freed memory. Those
descriptors are ConfScan's pre-filter for skipping RMSD comparisons, so silently wrong values were
possible — not only crashes. ConfScan's own test set is unchanged by the fix (`confscan_free`,
`confscan_hybrid`, `confscan_molalign`, all `cli_confscan_*` pass; `confscan_dtemplate` fails before
and after — a pre-existing, thread-count dependent golden-value test).

## 7. What counts as "the same molecule": `-topology_factor` (default 1.3)

Every topology comparison in the search used `Molecule::DistanceMatrix()`, whose covalent-radius
factor is **1.5**. That is far too generous for this decision, and the consequence is not subtle.

Measured on a 107-atom peptide (C34H55N11O7), 10 ps of gfnff MD at 500 K:

| covalent factor | "bonds" in the reference | MD frames counted as topology changes |
|---|---|---|
| **1.5** (old) | 110 | **92 %** |
| 1.4 | 108 | 0 % |
| **1.3** (now) | 108 | **0 %** |
| 1.2 | 108 | 8 % |

The molecule never reacted. The two extra "bonds" at 1.5 are **compressed 1-3 contacts**: atoms 5
and 12, both bonded to atom 4, sit 2.25 Å apart, and the 1.5 criterion calls them bonded at 2.28 Å —
i.e. as soon as their angle closes from 107° to 109°. A tetrahedral angle is 109.47°, so the "bond"
appears and disappears on a two-degree thermal fluctuation.

Phase 4 therefore discarded most of the conformers the search had just produced. Same molecule, one
ConfSearch cycle at 500 K, 4 × 10 ps of MD:

| RMSD-MTD setting | conformers @1.5 | conformers @1.3 | lowest gfnff @1.5 | @1.3 |
|---|---|---|---|---|
| defaults, 2 ps | 3 | **34** | −18.748979 | **−18.776946** |
| defaults, 10 ps | 1 | **54** | −18.759449 | −18.777362 |
| `bias_calibration couple` | 2 | **25** | −18.755169 | −18.778346 |
| couple + `rmsd_mtd_k 0.05` | 3 | **33** | −18.756473 | **−18.780005** |
| couple + `rmsd_mtd_k 0.1` | 1 | **54** | −18.735526 | −18.774007 |

**10 to 50 times more conformers, and a minimum ~55 kJ/mol deeper.** ConfGen had already hit this
and defined its own `-topology_factor` (1.3); ConfSearch kept using 1.5. It now has the same
parameter, used for the reference topologies, the snapshot gate, both Phase 4 filters and the
repair.

## 8. Topology gate before the optimisation (`-snapshot_topology_gate`, default ON)

A conformer search explores ONE molecule. A snapshot that formed or broke a bond is not a conformer
of it, and the Phase 4 filter rejects it — but only *after* a full optimisation has been paid for it.
Those same geometries are the ones GFN-FF cannot differentiate:

```
[ERROR] GFN-FF: combined gradient contains NaN/Inf - scanning per-term contributions
[ERROR]   [NaN] term=hb first_atom=8 axis=0
[ERROR]   [NaN] term=batm first_atom=8 axis=0
[ERROR] NaN/Inf in gfnff gradient - energy -2.13026119 Eh is finite and returned, but the gradient
        is invalid (optimization/MD should abort or restart)
```

`hb` and the three-body `atm`/`batm` terms carry 1/r factors, so two atoms on top of each other give
a finite energy with an infinite derivative. (The per-term scan is a diagnostic and a line search
re-enters the same geometry up to `max_line_search` times, so one bad structure produced ~20 of these
nine-line blocks. GFN-FF now prints the full scan for the first occurrence and every 100th; raise
`-verbosity` to 2 for all of them.) The snapshots are therefore screened against the
reference topology *before* Phase 2:

```
ConfSearch: topology gate: 959 of 974 MD snapshots failed the check (42 formed a bond (collision),
            917 broke one), none repaired -- 15 snapshots go into the optimisation
```

(butane at 3000 K, deliberately destructive — 959 optimisations saved in one cycle). Formed and
broken bonds are counted separately because they mean different things: a formed bond is a collision
(the NaN source), a broken one is fragmentation.

The criterion is the same covalent-radius rule as Phase 4 (factor 1.5), which is generous enough that
thermal stretching does not trip it: a C-C bond must exceed 2.28 Å to count as broken. Measured at
800 K on butane, *all* snapshots pass, so a normal run is unaffected. When nothing survives at all,
the cycle is skipped with an explicit warning rather than running the phases on an empty file.

**A topology test is not enough.** A pair that is ALREADY BONDED forms no new bond when it
collapses, so the topology check sees nothing wrong — while the 1/r factors in the hydrogen-bond and
three-body terms still blow up. That is exactly what was observed: snapshots passed the gate and
their optimisation then reported NaN in `hb`/`atm`/`batm` with a finite energy. `-snapshot_clash_ratio`
(default 0.55) therefore screens the raw distances as well: any pair closer than 55 % of the sum of
its covalent radii is a collapsed geometry, bonded or not (for C-H that is 0.59 Å against a real bond
of 1.09 Å — deeply pathological, never thermal). 0 disables it.

**The structures that still fail are written out.** Whatever survives the screens and then breaks
the method during the optimisation is saved for analysis instead of only being counted:

```
ConfSearch: 3 structure(s) the method could not handle written to
            <bmt>/input.cycle02_T550K.explore.gfnff.failed.xyz (input geometries, named with the
            reason) -- reproduce with: curcuma -sp <structure> -method gfnff
```

What is saved is the **input** geometry — that is the thing that makes the force field produce NaN,
and it is what you reload, inspect and reduce to a bug report; the optimised geometry is the diverged
one and useless. Each frame is named `failed_007_method_nan [gfnff]`, with the reason being
`method_nan` (the method reported NaN/Inf), `no_geometry` or `non_finite`. Plain non-convergence is
deliberately **not** included — ConfSearch accepts unconverged structures on purpose, so that would
dump most of the batch. The file is per stage and per cycle, like every other output.

Verified end to end: a butane with one H placed exactly on a C is caught, written, and the saved
frame reproduces the failure on its own (`curcuma -sp` on it fails in the GFN-FF charge setup).

**And NaN defeats numeric guards.** Every comparison with NaN is false, so a diverged optimisation
sailed past the energy-rise limit and every convergence test, and kept re-evaluating a dead geometry
(20 line-search evaluations per step). Three places now check explicitly rather than by comparison:
the optimiser aborts the structure on a non-finite energy or coordinate and keeps the last valid
geometry; ConfSearch never writes a non-finite geometry into the pool (it falls back to the input
structure, or drops the entry and says so); and the gate screens non-finite snapshots up front.

Two related knobs: `-snapshot_topology_gate false` restores the old behaviour (optimise everything),
and `-topo_check true` makes the MD itself abort as soon as the molecule fragments — the gate stops
you paying for broken structures, `topo_check` stops the dynamics from producing them.

**Independently fixed:** the LBFGSpp objective only checked the ENERGY for NaN, so a finite energy
with a NaN gradient went straight into the line search. A non-finite gradient is now an error, which
makes the driver take a zero step and stop on that structure instead of following garbage.

## 9. Live progress for the optimisation and MD batches

Two different mechanisms exist, and mixing them up is what caused the silence:

| | writes to | shows | honours `-noprogress` / non-TTY |
|---|---|---|---|
| `CxxThreadPool::setProgressBar` | **stderr** | percent only, no label, no `i of N` | no (only the `CxxThreadBar` env var) |
| `CurcumaLogger::progress_bar` | **stdout** | label + `i/N` + percent | yes |

The pool bars were switched off below verbosity 3 in Jul 2026 for a good reason — but only for
**ConfScan**, where the pool is `Reset()`+`StartAndWait()` **per candidate structure**, so a scan
produced one 100-character bar per permutation batch. That gating was applied to ConfSearch's MD and
optimisation pools at the same time, where it is not needed (one pool = one batch = one bar) and
where it left the batch completely silent.

Worse, `ConfSearch::PerformOptimisation` was silent **after** the batch too: `OptThread::execute()`
sets the shared static logger level to the worker level (0), and the level was only restored at the
end of the function — after the per-structure table and the batch summary. At verbosity 1 a
43-structure batch printed `Optimizing 43 structures using 1 thread(s) [gfnff]` and then nothing at
all. The level is now restored immediately after `StartAndWait()`.

What is shown now:

| verbosity | during the batch | after |
|---|---|---|
| 1, terminal | one live bar: `Optimising [gfnff] [########----] 34% 17/50` (MD: `MD runs [gfnff] …`) | one summary line with count + wall time |
| 1, redirected / `-noprogress` | nothing for opt; the MD counter lines are kept so a log still records progress | same summary |
| >= 2 | one line per finished structure (index, steps, energy, per-structure and elapsed time) | the ordered per-structure table |
| >= 3 | additionally the CxxThreadPool stderr bar | as above |

ConfScan's per-permutation pools stay at `None` below verbosity 3 — unchanged. ConfScan's own
overall bar (`updateProgress`, stdout) and these batch bars can never be live at the same time,
because the ConfSearch phases run strictly sequentially.

`curcuma -opt` on a multi-frame file gets a plain `--- Structure 3/17 (114 atoms) ---` header per
frame instead of a bar: on that path the optimiser prints its per-iteration table at verbosity 1,
which would shred an in-place bar.

## 10. Cheap pre-filter between RELAX and REDUCE (Aug 2026)

The deduplication in step 3 is quadratic — 133 structures need ~170 s, 10 000 would need days — and
it is fed by a batch of relaxed MD snapshots that mostly sit in the *same* minima. Two exact, cheap
criteria therefore run first:

| Flag | Default | Meaning |
|------|---------|---------|
| `-reduce_prefilter_window` | `100.0` | drop relaxed structures more than this many kJ/mol above the cycle minimum (0 = off) |
| `-reduce_prefilter_energy_tol` | `1e-6` | two relaxed energies agreeing to this many Hartree are the SAME minimum; keep one (0 = off) |

```
ConfSearch: pre-filter before REDUCE: 27 of 37 structure(s) kept (0 outside the 100 kJ/mol window,
            10 already at an identical minimum)
```

Both are decisions a human would not argue with: a structure 400 kJ/mol above the minimum is not a
conformer anyone is looking for, and two optimisations that ended at the same energy to 1e-6 Eh
ended in the same well. The collapse costs O(N log N) (sort by energy, compare neighbours) instead
of O(N²) Kabsch fits. The RMSD pass still decides everything else — the tolerance is far too tight
to merge two genuinely different conformers.

Measured motivation: with `-bias_reset deposits` the topology gate suddenly let 10 751 snapshots
through in one cycle, all of which would have gone into that quadratic pass. It also removes the
trigger of a ConfScan pathology, without fixing it: ConfScan derives its adaptive energy threshold
as a `std::max` over all pairs closer than `getrmsd_thresh`, so **one** outlier pair set ΔE to
190 kJ/mol and 7480 of 7484 structures were rejected as "duplicates" (4 conformers left). The
`std::max` derivation is still fragile for other callers.

## 10b. Keeping the wrong species in the bias (`-bias_rejected`, opt-in)

A structure that fails the topology check is the wrong *species*, so it must not enter the ensemble
— but the RMSD-MTD bias is purely geometric (the stored energy is metadata and never enters the
force), so the region it occupies can still be marked as visited. `-bias_rejected true` deposits
those structures as persistent hills; they stay out of the cumulative pool and out of the seed
selection.

```
ConfSearch: 32 topology-rejected structure(s) deposited as persistent hills -- they stay out of the
            ensemble, but the next trajectory is pushed away from them
```

Measured on the 107-atom peptide (one 600 K cycle, 32 rejects / 64 valid): the rejects sit
2.24–3.54 Å (median 2.82) from the nearest valid conformer, while the valid conformers sit
2.20–4.80 Å from *each other*, and **none** of the 32 is within the 1.25 Å deduplication threshold
of a valid structure. A hill on a reject is therefore no closer to a good conformer than the hills
the search already places on its own minima.

What it does not do: the proton transfer that produces those rejects happens during the
**optimisation** of a snapshot, so the trajectory never visits the rejected geometry itself. Section
11 addresses the cause — and with `-hold_polar_h` those structures survive as neutral conformers and
enter the bias on the regular path anyway.

## 11. Holding the polar X-H bonds during every optimisation (`-hold_polar_h`, opt-in)

A proton transfer turns a conformer into a **tautomer**, which the topology filter then discards —
after the optimisation has been paid for. Measured on the 107-atom peptide, one 600 K cycle: **32 of
96** structures were rejected for the identical transfer, H107 from O57 to N25, and those zwitterions
were the *deepest structures of the cycle* — all ten lowest and 29 of the 32 lowest, the best of them
72.5 kJ/mol below the best neutral conformer and 17 kJ/mol below the (neutral) GOAT reference. The
transfer happens **during the relaxation of the snapshot** (step 2), not in the snapshot itself, so
the snapshot topology gate cannot see it.

Recovering them afterwards does not work: rebuilding the neutral tautomer from all 32 and
re-optimising sends 30 of 31 straight back into the zwitterion, and with an O–H restraint the neutral
form sits ~93 kJ/mol higher and returns to the original zwitterion energy to six decimals the moment
the restraint is released. In those geometries the neutral tautomer is not a minimum at all —
preventing the transfer is the only useful handling.

`-hold_polar_h true` (old name `-refine_hold_polar_h` still accepted) puts a harmonic distance
restraint (`-hold_polar_h_force`, default 5 Eh/Å²) on every H bound to N/O/F/S **in the reference
structure**, at the reference bond length, in **every** optimisation of the search — step 2 RELAX and
both stages of step 5 REFINE, each against its own surface's reference:

```
ConfSearch: holding 7 polar X-H bond(s) at their reference length in every optimisation
            (k = 5 Eh/A^2) -- prevents proton transfer into a tautomer
```

A bond sitting at its equilibrium feels no force from the restraint, so the reported energies are
essentially unperturbed; only a proton that starts to leave feels anything. C-H bonds are left
alone deliberately — they do not migrate, and restraining them would only stiffen the molecule.

## Formatting

`CurcumaLogger::header()` is gated at verbosity >= 2, so at the default level a run printed one
undifferentiated stream of `[RESULT]` lines. `CurcumaLogger::section(title, major)` was added
(verbosity >= 1): a blank line plus `--- title ---`, or a full-width `=` frame for cycle boundaries.
Temperature cycles and all phases use it, and every cycle ends with a blank line.
