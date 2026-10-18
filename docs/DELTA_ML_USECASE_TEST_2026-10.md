# Delta-ML application test: ranking, geometry, cost (Oct 2026)

Status: 🤖 AI-generated scripts, ⚙️ machine-tested. One molecule. Not human-reviewed. Measurements by an agent
(Sonnet), the label/deployment analysis below by the orchestrator. Scripts: `scripts/deltaml_usecase/`
(absolute paths to the original working directory are hard-coded; adapt `common.py`).

## Question

Does Delta-ML (GFN-FF base, GFN2 reference) improve conformer ranking and pre-optimisation compared with GFN-FF alone,
and at what cost? Earlier measurements only covered energy MAE and gradient RMSE on random held-out splits
([DELTA_ML.md](DELTA_ML.md)).

## Setup

- Data: 353 conformers of one 90-atom C/H/N molecule from a trajectory (dataset JSON of `-activelearn`).
- Models: A = default (1000 centres), B = best joint configuration of the parameter scan (4000 centres, lambda 1e-7,
  rcut_radial 7, zeta 1,2,4,8,16), C = the energy-only model (`.deltaml.energy.json`, lambda 1e-4).
- Test 1: 5-fold CV x 2 repeats (n = 10 folds, about 70 candidates each), predictions through `-sp -method deltaml:<model>`.
  Test 2: 20 random conformers, models trained on the other 333, optimised with `-opt.optimizer lbfgs` against
  GFN2 optimisations. Test 3: ms per MD step by run-time difference, 5 repeats, 8 threads.

## Test 0: how hard is the ranking?

The relative GFN2 energy of the 353 conformers spans only 2.2 kcal/mol (std 0.56; 25, 91, 288 conformers within 0.5, 1,
2 kcal/mol of the minimum). GFN-FF spans 6.85 kcal/mol (std 1.23), Pearson 0.568, Spearman 0.490 against GFN2, RMSE after a
linear fit 0.464 kcal/mol. The model errors found earlier (0.6 to 1.1 kcal/mol) are of the size of the whole GFN2 spread.

## Label/deployment mismatch (found while evaluating, affects Tests 1 and 2)

The base energy stored in the dataset (`energy_base`, computed by a reused calculator with the pinned bond graph of
the anchor) is **not** the energy `-sp -method gfnff` returns for the same geometry: fresh minus stored =
-1.21 +- 0.70 kcal/mol over all 353 structures (min -3.05, max 0.00, no trend with the index). The model is trained on
`E_ref - E_stored` but `-sp -method deltaml:<model>` adds the correction to the fresh energy. The cause is not
diagnosed; frozen quantities of a reused GFN-FF calculator (open items in the root CLAUDE.md: D4 ATM C6 frozen at
setup, `shouldUpdateHBXB()` trigger) are candidates. Consequences measured:

- True GFN2 minimum: rank 7 of 353 with the stored base energy, rank 31 with fresh GFN-FF. The baseline of Test 1 as
  first reported (stored energies) was therefore more favourable than what `-sp -method gfnff` gives.
- Test 1 evaluated in two frames from the same 706 predictions (`scripts/deltaml_usecase/frames.py`; n = 10 folds):

| frame | method | top-1 | recall@5 | regret [kcal/mol] | Spearman | paired Spearman vs base |
|---|---|---|---|---|---|---|
| as trained (stored base + correction) | base | 0.50 | 0.72 | 0.12 | 0.48 | - |
| | A | 0.10 | 0.28 | 0.78 | 0.33 | -0.146 +- 0.095 (0/10 better) |
| | B | 0.30 | 0.38 | 0.88 | 0.38 | -0.099 +- 0.112 (3/10) |
| | C | 0.60 | 0.64 | 0.08 | 0.47 | -0.009 +- 0.085 (7/10) |
| as deployed (fresh GFN-FF + correction) | base | 0.30 | 0.30 | 0.63 | 0.52 | - |
| | A | 0.10 | 0.20 | 1.19 | 0.34 | -0.179 +- 0.136 (2/10) |
| | B | 0.10 | 0.12 | 0.89 | 0.43 | -0.091 +- 0.138 (3/10) |
| | C | 0.20 | 0.40 | 0.64 | 0.67 | +0.155 +- 0.078 (9/10) |

Reading: A and B are worse than the base in both frames (Spearman and recall@5). C is level with the base in the
as-trained frame; its Spearman gain exists only in the mixed frame, so it cannot be credited to the model.
The mismatch itself is a result: GFN-FF energies of the same geometry differ by 1.2 +- 0.7 kcal/mol depending on how
the calculator was set up.

## Test 2: pre-optimisation (n = 20, lbfgs, all runs exit 0)

| method | RMSD to GFN2 minimum [A] | GFN2 energy at the result minus GFN2 minimum [kcal/mol] | steps (median) |
|---|---|---|---|
| GFN-FF | 0.83 +- 0.45 | 10.45 +- 0.93 | 636 |
| A | 0.95 +- 0.32 | 69.65 +- 11.47 | 624 |
| B | 0.92 +- 0.24 | 85.58 +- 9.85 | 724 |
| C | 0.59 +- 0.34 | 9.83 +- 0.80 | 487 |

Paired against GFN-FF: A +59 +- 11 kcal/mol (0 of 20 better), B +75 +- 10 (0 of 20); C -0.62 +- 0.84 (14 of 20, within the
spread), RMSD -0.25 +- 0.36 A (18 of 20). At the geometries found by A and B the model's own energy is 100 to 165 kcal/mol
below GFN2: spurious minima. Cause not diagnosed (candidates: rough fit with small lambda away from the training data,
see the lbfgs failure of the old default grid in DELTA_ML.md; the stored/fresh mismatch is too small to explain it).
The starts are not exact GFN2 minima (RMSD start to GFN2 minimum 0.24 +- 0.28 A, max 1.06 A). C was run as a control; its
forces are no better than GFN-FF by construction.

## Test 3: cost (ms per MD step, 8 threads, n = 5)

GFN-FF 0.32 +- 0.05, model A 3.56 +- 0.08 (11x GFN-FF, 3.9 % of GFN2), model B 8.06 +- 0.07 (25x, 8.9 %), GFN2 90.2 +- 3.3.
Model A does not scale with the thread count (1 thread 3.61 ms). Fit time A 12 s, B 145 s, energy model 1.5 s and 9.5 s.
Delta-ML therefore costs 11 to 25 times GFN-FF, not "GFN-FF cost".

## Conclusion

For this molecule, with these models, Delta-ML does not improve conformer ranking or pre-optimisation relative to GFN-FF:
the joint models are worse in ranking and produce spurious minima, the energy model is level with the base in the
as-trained frame. Not tested: other molecules or energy ranges (this set is a trajectory within 2.2 kcal/mol), a DFT
reference, charged or open-shell systems, MD stability of the new defaults, a labelling that is consistent with
deployment (rerun after the mismatch is understood), the 8000-centre limit.

Caveats: the conformers are trajectory neighbours and the CV folds share training data, so the spread is
underestimated; the two repeats use the same 353 structures; an agent run overwrote the untracked file
`gfnff_diag_d4.json` in the repository root.
