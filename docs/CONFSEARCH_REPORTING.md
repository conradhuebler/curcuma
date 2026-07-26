# ConfSearch: reference energies, per-cycle output, two-stage refinement

Status: 🤖 AI-generated, ⚙️ machine-tested (Jul 2026). Not ✅ TESTED / ✅ APPROVED.

Four changes to what a conformer search reports and produces, plus one bug they uncovered.

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
ConfSearch: Phase 3 accepted ensemble [gfnff]: 36 structure(s), lowest -18.742131 Eh, span 41.28 kJ/mol
ConfSearch:     #1   -18.742131 Eh  (+0.00 kJ/mol)
ConfSearch:     #2   -18.741062 Eh  (+2.81 kJ/mol)
ConfSearch:     #3   -18.74055  Eh  (+4.15 kJ/mol)
ConfSearch:     ... and 33 more, up to +41.28 kJ/mol
```

`-ensemble_report N` sets how many conformers are listed (0 = summary line only).

## 3. Per-cycle ensembles and most stable structures, on both levels

The working files of a cycle (`*.bias.opt.accepted.xyz` and its `opt_method` counterpart) are
overwritten by the next cycle, so what a given temperature produced was not recoverable afterwards.
Every cycle now writes, for **each level of theory**:

| File | Content |
|---|---|
| `<base>.cycleNN_TxxxK.<method>.xyz` | the cycle's topology-valid ensemble, **energy-sorted** |
| `<base>.best_per_cycle.<method>.xyz` | one frame per cycle: its most stable structure (a trajectory of the search's progress) |

```
ConfSearch: cycle01_T500K ensemble [gfnff]: 36 structure(s), most stable -18.742131 Eh, span 41.28 kJ/mol -> <bmt>/input.cycle01_T500K.gfnff.xyz
ConfSearch:   most stable structure of this cycle appended to <bmt>/input.best_per_cycle.gfnff.xyz
```

Switch off with `-cycle_output false`. Ordering is done by a `std::multimap` keyed on the energy —
the container sorts, there is no explicit sort call (see section 5).

## 4. Two-stage high-level re-optimisation (`-phase3b_two_stage`, opt-in)

Phase 3b re-optimises the cycle's deduplicated minima at `opt_method` and is the most expensive step
of a dual-method run. Its input was deduplicated on a **different** potential-energy surface
(`md_method`), so structures that survive as distinct there can collapse onto one another as soon as
the accurate method relaxes them — and each of those collapses was paid for at full price.

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
| `-phase3b_two_stage` | `false` | run crude -> filter -> accurate instead of one accurate optimisation |
| `-phase3b_preopt_preset` | `loose` | convergence preset of the crude stage |
| `-phase3b_preopt_max_iter` | `0` | step cap for the crude stage; 0 keeps the preset value |
| `-phase3b_preset` | `normal` | convergence preset of the accurate stage (also of single-stage Phase 3b) |
| `-phase3b_filter` | `true` | run the dedup between the stages |

Whether it pays depends on how much the two surfaces disagree — it is off by default and every run
states which way it is set. It only exists in dual-method runs (Phase 3b is skipped when
`opt_method == md_method`).

## 5. ConfScan sorting: verified, no `std::sort` needed

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

## 6. Live progress for the optimisation and MD batches

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

## Formatting

`CurcumaLogger::header()` is gated at verbosity >= 2, so at the default level a run printed one
undifferentiated stream of `[RESULT]` lines. `CurcumaLogger::section(title, major)` was added
(verbosity >= 1): a blank line plus `--- title ---`, or a full-width `=` frame for cycle boundaries.
Temperature cycles and all phases use it, and every cycle ends with a blank line.
