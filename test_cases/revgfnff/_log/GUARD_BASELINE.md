# Guard baseline: gfnff / revgfnff vs the published references (2026-09-13)

## Provenance

| field | value |
|---|---|
| binary | `release/curcuma`, copied unmodified into the session scratchpad |
| md5 / size / mtime | `58512a18c523456e8db7ca826ad84897` / 23712744 B / 2026-09-12 17:41:57 +0200 |
| embedded version | `0.0.263-66-g4ef8d30c`, branch `reactff2-llm` (`curcuma -version`) |
| = commit | **4ef8d30c** "Fix the stale total energy in the final MD table row" (2026-09-12 12:50) |
| vs HEAD `29146ad0` | 4 `src/` commits behind: **pre-`c_ij` AND pre-`e36d9925`** (stage 3a (i) r0 fix) |
| r0 fix absent, 2nd proof | ch4 C-H at 1.4 r_eq: revgfnff error **+15.8** kcal/mol, the recorded PRE-fix value (`R0_FIX_STATUS.md`, post-fix +2.6); h2o O-H +16.2 (post-fix -1.4). Equilibrium rev-gfnff offset +0.0004 kcal/mol, a fixed build is ~-0.02/bond |
| working tree | the in-flight `c_ij` edits (`ff_workspace*`, `gfnff*.*`, mtimes 2026-09-13 21:41-22:51) postdate the 09-12 link and are NOT in this binary |
| threads | 2 per process, 4 processes = 8 |
| topology hygiene | fresh directory + unique basename per structure; batch runs WITHOUT `-batch_reuse_topology` (cache_topology off); no `<basename>.topo.json` reused |
| invocation | `-method gfnff|revgfnff`, `-charge`/`-spin` from `.CHRG`/`.UHF`, `-verbosity 0 -no_bmt` |

The r0 fix is rev-only, so the **gfnff column is a valid current baseline**; the **revgfnff column is as of 4ef8d30c** (see section 5).

## 1 + 2. Conformer sets and S66 (kcal/mol) - reactions, vs the published reference

Scored with `scripts/gmtkn55_reactions.py`'s own `.res` parse / `evaluate()` / `score()`, i.e. the same quantity behind `test_cases/GMTKN55-testset/_results/reactions_summary.md`. Reaction, not structure: a conformer pair, or dimer minus monomers. max = signed largest-|error| reaction.

| set | gfnff MAD | gfnff max | revgfnff MAD | revgfnff max | n |
|---|---:|---:|---:|---:|---:|
| ACONF | 0.155 | -0.358 | 0.155 | -0.358 | 15 |
| ICONF | 3.310 | -20.157 | 3.310 | -20.157 | 17 |
| MCONF | 0.589 | +1.766 | 0.589 | +1.766 | 51 |
| PCONF21 | 1.648 | +5.253 | 1.648 | +5.253 | 18 |
| S66 | 0.825 | -2.649 | 0.825 | -2.649 | 66 |
| **conformers, all 8 sets** | **1.492** | -20.157 | 1.492 | -20.157 | 285 |
| conformers, all 8 sets (RMSD) | 2.303 | - | 2.303 | - | 285 |
| conformers, the 4 sets above only | 1.171 | -20.157 | 1.171 | -20.157 | 101 |

revgfnff == gfnff here: every set agrees to 3 decimals; per reaction the two differ by at most 2.0e-4 kcal/mol (351 reactions), per structure by at most 3.1e-3 (UPU23/1e; mean 1.3e-3 - the near-constant r0 offset cancels in a reaction).

## 3. Class D - MD frames vs r2SCAN-3c (20 systems x 25 frames = 500 points)

dE_k = E_k - E_0 of that system (the method's energy zero cancels), MAD/RMS over the 25 points; grad_RMS = sqrt(mean_k(|g_cur - g_ref|^2 / 3N)) in kcal/mol/Angstrom, both gradients Eh/Angstrom (reference Eh/Bohr / BOHR, as `scripts/revgfnff_data.py`) - the convention of `revgfnff_fit.FitContext.evaluate`. Per-frame perception: each frame gets its own topology, which is what a user running frames one at a time gets.

| system | gfnff dE_MAD | revgfnff dE_MAD | gfnff grad_RMS | revgfnff grad_RMS | n |
|---|---:|---:|---:|---:|---:|
| c2h4_1000K | 5.890 | 5.890 | 9.931 | 9.931 | 25 |
| c2h4_2000K | 11.751 | 11.751 | 14.253 | 14.254 | 25 |
| c2h6_1000K | 1.303 | 1.303 | 5.108 | 5.108 | 25 |
| c2h6_2000K | 2.991 | 2.991 | 7.221 | 7.221 | 25 |
| ch3cl_1000K | 2.592 | 2.592 | 11.783 | 11.783 | 25 |
| ch3cl_2000K | 3.840 | 3.840 | 14.533 | 14.533 | 25 |
| ch3f_1000K | 5.633 | 5.633 | 25.327 | 25.327 | 25 |
| ch3f_2000K | 8.920 | 8.920 | 27.301 | 27.301 | 25 |
| ch3nh2_1000K | 2.557 | 2.557 | 8.405 | 8.405 | 25 |
| ch3nh2_2000K | 4.478 | 4.477 | 12.320 | 12.321 | 25 |
| ch3oh_1000K | 3.578 | 3.578 | 16.241 | 16.241 | 25 |
| ch3oh_2000K | 4.986 | 4.986 | 16.571 | 16.571 | 25 |
| ch4_1000K | 1.939 | 1.939 | 5.291 | 5.291 | 25 |
| ch4_2000K | 6.241 | 6.241 | 10.868 | 10.868 | 25 |
| h2co_1000K | 10.526 | 10.526 | 30.183 | 30.183 | 25 |
| h2co_2000K | 14.206 | 14.206 | 32.312 | 32.312 | 25 |
| h2o_1000K | 0.819 | 0.819 | 8.295 | 8.295 | 25 |
| h2o_2000K | 1.285 | 1.285 | 10.609 | 10.610 | 25 |
| nh3_1000K | 0.814 | 0.814 | 7.721 | 7.721 | 25 |
| nh3_2000K | 1.151 | 1.151 | 9.866 | 9.866 | 25 |
| **pooled** | **4.775** | 4.775 | **16.305** | 16.305 | 500 |
| pooled, dE_RMS | 7.391 | 7.391 | - | - | 500 |
| pooled, max frame dE | 39.40 | 39.40 | - | - | 500 |
| pooled, `-batch_reuse_topology true` | 4.807 | 4.807 | 16.466 | 16.467 | 500 |

Both topology modes are given because the fitter's class-D guard runs with reuse=true while a normal run re-perceives: they differ per system by up to 4.5 % (h2co_2000K 14.21 vs 14.85), elsewhere <0.3 %. On every row revgfnff == gfnff to <=1.3e-4 kcal/mol (dE_MAD) and <=9.0e-4 (grad_RMS) - the r0 offset cancels in a relative energy.

## 4. Recorded numbers: cited vs re-derived

| number | where recorded | status |
|---|---|---|
| gfnff 285 conformer reactions MAD **1.49** | `REV_GFNFF_ROADMAP.md` validation matrix; `REV_GFNFF_DATA_BASIS.md`; `REV_GFNFF_TODO.md` entry 9; `_results/reactions_summary.md` | **cited, independently reproduced** = 1.4924 over the full 285 (all 8 subsets are present locally, so nothing was extrapolated) |
| gfnff S66 MAD **0.83** (n=66) | `FABLE_ROADMAP_REVIEW.md:225`; `reactions_summary.md` S66 = 0.825 | **cited, independently reproduced** = 0.8252 |
| gfnff ACONF/ICONF/MCONF/PCONF21 0.155 / 3.310 / 0.589 / 1.648 | `reactions_summary.md` per-subset table | **cited, reproduced identically** (0.1553 / 3.3099 / 0.5888 / 1.6482) including every max (+/-0.3578, -20.1572, +1.7665, +5.2535) |
| class D (500 points) | nowhere - the review calls it "never used as an accuracy target" | **not recorded; measured here for the first time** (section 3) |
| revgfnff on any of these sets | nowhere (`reactions_summary.md` holds gfnff/gfn1/gfn2/xtb/pbeh3c only; `_results/` has no revgfnff reaction file) | **first measurement here** |

Not re-derived (out of scope, same validation matrix): MOR41 vs pprcht 0.00052, GMTKN55 NCI 7.76 / charged NCI 38.9, barriers 42.9, S30L-CI 0.386.

## 5. What the revgfnff column does and does not cover

`release/curcuma` predates stage 3a (i), so revgfnff here is **pre-r0-fix**. Per `R0_FIX_STATUS.md` that fix preserves relative energies to 0.001 kcal/mol but moves the stretched-bond region by up to 13 kcal/mol (ch4 C-H at 1.4 r_eq +15.8 -> +2.6; h2o O-H +16.2 -> -1.4):

- section 1+2 (conformers, S66): should move by <=0.001 kcal/mol per reaction - expected to survive.
- section 3 (class D): dE_MAD/dE_RMS and grad_RMS will **change** on the elongated frames; the values above are a pre-fix floor, not a post-fix number.
- gfnff (first column everywhere): unaffected; the fix is rev-only and `gfnff` was verified bit-identical (`aef7aaf0`).

Re-measure all three after the next `release/` build; the revgfnff class-D floor is still open.
