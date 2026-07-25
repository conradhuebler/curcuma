# ConfSearch seed selection: energy + RMSD diversity

Status: 🤖 AI-generated, ⚙️ machine-tested (Jul 2026). Not ✅ TESTED/APPROVED.

## What a "seed" is

After every temperature cycle ConfSearch has a set of optimised, topology-valid minima. Those that
fall inside `seed_energy_window` (kJ/mol above the running global minimum) are candidates for the
**seeds** of the next, colder cycle: each seed starts `repeat` independent MD runs. Seeds are the
scarcest resource in the search — `repeat * n_seeds` MD runs are the whole cost of a cycle.

## Problem with pure energy ranking

The old rule was `seed_rank` lowest energies, full stop. If the three most stable minima are small
variations of the *same* basin (0.3–0.5 Å apart) and the first structurally different conformer is
2 kJ/mol higher, all MD time goes into re-exploring one region while the neighbouring basin is
dropped. The RMSD-MTD bias partly compensates (it pushes away from what it already saw), but the
starting points still all sit in the same funnel.

## Rule since Jul 2026 (`seed_selection`, default `diverse`)

Walk the same energy ranking, but require geometric spacing:

1. Sort the window survivors by energy. Seed #1 is always the most stable structure.
2. Accept the next-lowest candidate only if its **permutation-aware best-fit RMSD**
   (`ConfSearch::PermRMSD` — Kabsch + the symmetry reorder rules ConfScan collected) to *every*
   already accepted seed is `>= r_min`.
3. Stop at `seed_rank` seeds (`seed_rank = 0` means: no count limit, spacing still applies).
4. If fewer than `seed_rank` seeds were found, relax `r_min -> r_min/2 -> r_min/4 ...` and repeat.
   The search is never starved of seeds; in the worst case the result equals the old energy ranking.

`r_min = seed_min_rmsd`, or `seed_diversity_factor * rmsd` when `seed_min_rmsd = 0` (the default:
2 x the dedup threshold, i.e. seeds sit at least one dedup radius apart from each other).

`seed_selection energy` restores the pre-Jul-2026 behaviour exactly.

## Parameters

| Flag | Default | Meaning |
|------|---------|---------|
| `-seed_selection` | `diverse` | `diverse` (energy + RMSD spacing) or `energy` (lowest N only) |
| `-seed_rank` | `1` | Max seeds per cycle. `0` = every candidate that passes window + spacing |
| `-seed_min_rmsd` | `0.0` | Spacing in Angstrom. `0` derives it from `rmsd` |
| `-seed_diversity_factor` | `2.0` | Multiplier on `-rmsd` when `seed_min_rmsd = 0` |
| `-seed_energy_window` | `50.0` | kJ/mol window vs the global minimum that defines the candidates |

**Note:** with the default `seed_rank 1` only one seed is kept, so the diversity rule has nothing to
choose and the run is bit-identical to before. Diversity only matters with `-seed_rank N > 1` or
`-seed_rank 0`.

## Output

At verbosity >= 1 each cycle reports the selection and the actual spacing achieved:

```
ConfSearch: seed selection (diverse): 4 of 26 structure(s) kept, spacing >= 2.50 A
ConfSearch:   seed  1: dE =   +0.00 kJ/mol, RMSD to closest previous seed = --
ConfSearch:   seed  2: dE =   +1.83 kJ/mol, RMSD to closest previous seed = 3.11 A
ConfSearch:   seed  3: dE =   +4.02 kJ/mol, RMSD to closest previous seed = 2.74 A
```

A `(relaxed to X.XX A -- not enough distinct basins)` suffix means the requested spacing could not
be met and was lowered; that is a hint that `seed_rank` is larger than the number of genuinely
distinct low-energy basins found so far.

## Cost

`O(n_kept * n_candidates)` Kabsch fits per cycle on a handful of structures — negligible next to a
single MD step.

## Not tested / not implemented

- No test on systems where the permutation cache is large (many equivalent groups); the RMSD is then
  the min over all cached images, which is correct but slower.
- The selection is greedy, not an optimal max-min (farthest-point) cover. A candidate that would
  have enabled a better overall spread can be blocked by an earlier, lower-energy pick.
- Energy and diversity are not traded off continuously (no score function) — energy order is strict,
  diversity is a hard filter on top.
- No genetic/crossover step exists in ConfSearch (CREST's "GC" phase has no counterpart here).
