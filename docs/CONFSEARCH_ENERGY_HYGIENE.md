# ConfSearch: two passes that fix the energies a decision is made on

**Status: AI-generated, machine-tested. Not human production tested.**

A conformer search decides everything on numbers it produced itself: which structures fall inside the
energy window, which one represents a basin after deduplication, which ones seed the next cycle, and
what the run finally reports. One systematic offset sits on those numbers when `-hold_polar_h` is on --
the restraint energy -- and section 2 removes it. Section 1 documents a second offset that was claimed,
measured wrongly, and retracted.

| pass | when | flag | what it removes |
|---|---|---|---|
| RELAX polish | after every RELAX batch, before window/dedup/seeding | `-relax_polish_rank` (**0**, retracted) | nothing measurable -- see section 1 |
| polar-H release | once, after the final deduplication | `-release_polar_h_rank` (100) | the restraint energy of `-hold_polar_h` |

---

## 1. RETRACTED: "the optimiser stops early" — `-relax_polish_rank` (now default 0)

**Status (2026-09-05): the measurement below was wrong, and the pass is off by default.** It is kept
as an opt-in for a method whose optimiser really does stop early; nothing measured so far does.

### The original measurement, and its error

40 structures were drawn from `s2_relax_kept` across all seven temperature stages of a finished
production run (`20260823_211047`) and re-optimised with the same convergence preset that had produced
them:

```
median 3.80   mean 3.84   p90 6.10   max 6.74 kJ/mol   (n = 40)
```

The first version of this document read that as one premature optimiser stop ("pass 2 recovers
everything, passes 3 and 4 recover 0.00") and stated that the run had used no restraints in RELAX. It
had: the run's command line (fish history) and `input_md_params.json` carry `-hold_polar_h true`. The
re-optimisation was done **without** the restraints. So the 3.80 kJ/mol were the restraint energy of
16 polar X-H bonds, released — exactly the "reading trap" section 2 warns about, committed one section
earlier.

### The control measurement

Six RELAX products of the same run, two each from 550 K, 450 K and 350 K, re-optimised twice with the
same preset: once free, once with the same 16 X-H restraints RELAX had used (k = 5 Eh/A^2, targets
from `input.s0_initial.gfn2.opt.xyz`):

| structure | E_relax (Eh) | free: gain (kJ/mol) / steps | restrained: gain / steps |
|---|---|---|---|
| cycle02_T550K_r3 #1 | -161.656185 | 3.33 / 54 | **0.00 / 1** |
| cycle02_T550K_r3 #2 | -161.648600 | 2.69 / 62 | 0.00 / 1 |
| cycle04_T450K_r2 #1 (parent of the record) | -161.663829 | 4.05 / 72 | 0.00 / 1 |
| cycle04_T450K_r2 #2 | -161.657581 | 4.16 / 56 | 0.00 / 4 |
| cycle06_T350K_r2 #1 | -161.658078 | 3.49 / 54 | 0.00 / 1 |
| cycle06_T350K_r2 #2 | -161.652563 | 6.31 / 56 | 0.00 / 2 |

Free: median 3.8 kJ/mol, the number of the original measurement. Restrained: 0.00 in every case, the
optimiser confirms convergence within one to four steps. The RELAX minima were converged. There is no
premature stop.

### What the record actually was

The chain still reads as before -- seed -14.55, RELAX product -15.73, "record" -19.79 at 0.02 A from
its parent -- but the 4.06 kJ/mol are the restraint energy of that conformer. The recombination step
re-optimised its templates **without** the restraints (the polar-H pass-through to ConfGen came later,
see [CONFSEARCH_PROPOSALS.md](CONFSEARCH_PROPOSALS.md)), so that run's pool mixed two energy scales:
restrained RELAX minima and free re-scored ones. Every "recombination product" that was really a
re-optimised template sat 0.4-5.6 kJ/mol (median 2.85, n = 226) below its own parent at identical
geometry. That is the artefact behind the "RECOMBINE is the dominant seed supplier" reading, and it is
gone in runs where both sides are restrained (second run: +0.01 kJ/mol median, n = 114).

The correct treatment of the offset is section 2: keep the guard during the search, release it once at
the end, report free minima.

### Where the pass stands now

- `-relax_polish_rank` default **0**. Under `-hold_polar_h` it is skipped even when requested (on a
  restrained run it measured 0.000 kJ/mol median over 200 structures, max 0.985, and diverged in 2 of
  10 restarts).
- The earlier claim that a second run's reported best (-8.93) was "2.34 kJ/mol worse than a structure
  the run already held (-11.27)" compared a restrained energy with a free one and is withdrawn.

## 2. `-hold_polar_h` is a constraint — `-release_polar_h_rank`

`-hold_polar_h` holds every reference N/O/F/S-H bond at its reference length in every optimisation of
the search. It has to: without it the deepest structures of a cycle are tautomers (measured: 32 of 96
in one 600 K cycle, and they were all ten of the lowest). But a restraint is a constraint, so every
energy computed under it sits above the free minimum of the same conformer.

**Measured**: 8 structures of a running production search, re-optimised without the restraints —

```
gains 0.41 ... 3.37 kJ/mol,  median 1.92
topology unchanged in all three checked   (largest displacement 0.06 A, on one hydrogen)
```

An **8-fold spread** across eight structures. That offset does not cancel in a ranking: two
conformers 1 kJ/mol apart can swap places purely on how much restraint energy each carries. And
because `REFINE` also runs restrained (`confsearch.cpp` — `opt_accurate` / `opt_crude`), the reported
final energies are restrained minima too.

> **Reading trap this resolves.** The same 8 structures, re-optimised *externally without*
> restraints, look exactly like the premature-stop effect of section 1 — median 1.92 kJ/mol. They are
> not: the in-run polish pass, which correctly keeps the restraints, finds 0.003 kJ/mol on the same
> structures. Any measurement of "how far from the minimum is this pool" must run under the same
> constraints as the pool, or it measures the constraint.

### What the pass does

Once, at the very end, after the final deduplication (deliberately after: the restrained energies
decided which structures are distinct, which is the comparison they are consistent for). The
`-release_polar_h_rank` lowest conformers of the accepted ensemble are re-optimised without the
restraints, and:

- topology unchanged and lower → the free geometry and energy replace the restrained ones;
- **proton moved** → the structure keeps its **restrained** geometry and is counted in a warning. It
  is a real conformer of the right species that only exists under the constraint; dropping it would
  lose a conformer the search legitimately found, and releasing it would report a tautomer.

The ensemble is re-sorted afterwards, because the gain differs per structure. Bounded by rank because
an accepted ensemble runs into the hundreds (940 in one measured run) while everything a reader ranks
sits at its top. `0` = off, `-1` = all. No effect when `-hold_polar_h` is off.

Files: `<base>.s7_release.<opt_method>.xyz` (input) and `.opt.xyz` (result), so the pass is auditable.

---

## What these passes do NOT do

They do not find conformers. Everything above is about the numbers attached to structures the search
already has. The measurement that motivated section 1 also showed the opposite case — a recombination
proposal that *did* find something: `proposal_from_1_d3_repaired`, 16.8 kJ/mol below the best RELAX
structure of its repetition, **2.72 A** from the nearest known minimum, topology intact. That is
`-confgen_phase`'s contribution, and it is a different mechanism from the second optimisation pass
that the same step happened to provide. See [CONFSEARCH_PROPOSALS.md](CONFSEARCH_PROPOSALS.md).

## Open

- Section 1's open questions are closed: there was no premature stop. Whether an unrestrained run
  would show one has not been measured, and nothing in the data suggests it.
