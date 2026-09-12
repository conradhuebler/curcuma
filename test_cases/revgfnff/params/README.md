# rev-gfnff parameter sets (over-coordination term E_over)

Override documents for `-gfnff.param_file FILE` (sparse, deep-merged over the built-in tables;
also the input format of `scripts/revgfnff_fit.py`). Both sets are selectable without a file
through `-gfnff.rev_over_preset stage1a | fit2026-09-12`; the fit is the built-in default since
2026-09-12 (operator decision).

| file | p_over H / C / N / O [Eh] | over_shift | valence N / O | origin |
|---|---|---|---|---|
| `rev_over_stage1a.json` | 0.3 / 0.3 / 0.3 / 0.3 | 0.5 | 3 / 2 | stage-1a placeholder (Sep 11, 2026) |
| `rev_over_fit_2026-09-12.json` | 0 / 0.209 / 0.339 / 0.994 | 0.870 | 2.536 / 2.531 | LM fit against class C (hyper-coordination curves) + GMTKN55 barriers BH76/WCPT18/PX13/BHPERI/BHDIV10 with the class-D guard (`docs/REV_GFNFF_ROADMAP.md` WP3) |

Barrier MAD in kcal/mol (gfnff / stage1a / fit): BH76 56.4 / 59.4 / 54.8, BHDIV10 39.9 / 47.8 / 33.2,
BHPERI 35.4 / 33.9 / 33.3, PX13 143 / 251 / 152, WCPT18 31.8 / 38.2 / 26.8; class-C rms(dE)
103 -> 42. React-MD smoothness (`scripts/revgfnff_jump_stats.py --dt 0.25`) is the same with both.

Compare by hand:

    release/curcuma -sp struc.xyz -method revgfnff -gfnff.rev_over_preset stage1a
    release/curcuma -sp struc.xyz -method revgfnff -gfnff.param_file test_cases/revgfnff/params/rev_over_stage1a.json
    python3 scripts/revgfnff_barrier_terms.py --method revgfnff --recompute [--extra -gfnff.param_file ...]
