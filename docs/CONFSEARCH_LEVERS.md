# ConfSearch: the levers, ranked by what the chains of 13 runs say

**Status: AI-generated analysis (Sep 2026) over the provenance chains of 13 finished or aborted
WEKLQ runs (107-atom peptide, gfn2 ranking). Counterfactual, not prospective: for every event of a
run's chain of running minima (structure that lowered the best energy) we ask which setting would
have cut it. Percentages are the share of the total descent (kJ/mol) that such events carried.**

Total descent over the 13 chains: 1806 kJ/mol in 205 events. Most of it is the initial fall from the
raw input (+150 to +80 kJ/mol) in the first repetition; the part that decides a search is the tail
below +20 kJ/mol, and the record run is the only chain that went far below the reference.

| lever | current default | share of descent it carried | where it decided | evidence class |
|---|---|---|---|---|
| **re-seeding cadence** (repetitions 3+ of a stage) | `-repeat` 2 in most runs, 5 in the record run | 10.7 % overall, **26 %** in both repeat-5 runs | **all four links of the record chain below the reference are in repetitions 2-5**, three of them in r3-r5; with repeat 2 that chain does not exist | 13 runs; the record itself n = 1 |
| **trajectory length** > 1000 fs | `-explore_md_time` 2000 | 7.7 % overall, **17-29 %** in the three 2000-fs runs | events at 1140, 1340, 1520, 1790, 1890 fs; 1000-fs runs lose them by construction | 3 runs |
| **proposals** (RECOMBINE) | `-confgen_phase` off | 12.1 % overall, 11-35 % where on | almost only at 600 K, from the raw start; the novelty gate discarded the record's skeleton once (fixed, `-new_energy_gain`) | 5 runs |
| **densification MD** | `-refine_md` off | 54-56 % in the two hybrid v2 runs (includes the early descent) | the depth mechanism inside a found region; per gfn2 optimisation it yields as many top-50 hits as exploration | 2 runs |
| **seed fan** (seeds ranked 5-10) | `-seed_rank` 10 | 7.7 % overall, 5-24 % in 8 of 13 runs | the "bad" seeds; any cut of seed_rank removes these events | 8 runs |
| **cold stages** (<= 400 K) | ladder to 300 K | 3.0 % overall, 12 % in two 1000-fs runs, 0 in every 2000-fs run | `-endT 450` saves ~40 % of the budget of a full ladder; the saving is worth more as repetitions at 600-450 K than as cold stages | 13 runs |
| **compactness** | none | -- | the RMSD bias leaves the seed's Rg upwards in 90 % of deposits; compact products only inherited | 17343 deposits |
| **bias metric** | all-atom RMSD | -- | 75 % of the hill-to-hill displacement is hydrogen motion | 17044 pairs |

What the table does NOT say: it counts what was found, not what a different setting would have found
instead (the cut events would partly have appeared elsewhere). It ranks the levers by how much of the
measured success depended on them.

## The levers as switches (Sep 2026)

| lever | switch | state |
|---|---|---|
| re-seeding cadence without exploration cost | `-refine_md true -refine_md_chain true` (`-refine_md_chain_max` 4) | built, default off, A/B open |
| energy-aware novelty | `-new_energy_gain 1.0 -new_pattern_min 1` (ConfGen) | built, default ON, retrospectively justified (61 + 38 cases) |
| compactness | `-rg_flood true` (+ pool hand-over in ConfSearch) | built, default off; three-arm test (Sep 5, n = 1 per arm): deposit-Rg distribution unchanged (median 5.46-5.51 A in all four arms), energies within the control spread -- rejected |
| heavy-atom bias metric | `-rmsd_mtd_heavy_only true` | built, default off, untested |
| mixed energy scales in one pool | `-hold_polar_h` reaches every optimisation; `-release_polar_h_rank` at the end | fixed; `-relax_polish_rank` retracted |
| redundant re-scoring of templates | see [CONFSEARCH_PROPOSALS.md](CONFSEARCH_PROPOSALS.md), "re-scoring" | 72-78 % of RECOMBINE's optimisations, 25-38 % of a run |

## The experiment that turns the ranking into a number

Equal budget (gfn2 optimisations + MD steps/100, logged), three seeds per arm, `-startT 600 -endT 450`,
2000 fs exploration, `-seed_rank 10`, `-hold_polar_h`:

- A: `-repeat 2` (baseline);
- B: `-repeat 2 -refine_md true -refine_md_chain true` (cadence lever without exploration cost);
- C: `-repeat 5` cut at the budget of A+B (cadence lever the expensive way);
- D: B + `-confgen_phase true` (proposals with the new novelty rule).

Read-outs: best on the free scale, conformers <= 20 kJ/mol, chain length and the repetition index of
every sub-reference link, run-to-run spread within each arm. The lever is confirmed when B beats A in
all three seeds or by more than A's own spread. Second molecule for the mechanics:
`test_cases/molecules/larger/triose.xyz`.
