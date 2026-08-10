# QMDFF@GFN2 as a conformer-search surface — measured

**Status: 🤖 AI-generated, ⚙️ machine-measured. Two ensembles, one reference method.**

**Short answer: it does not work on the system it was meant for.** On the 107-atom
peptide QMDFF@GFN2 ranks conformers *worse* than GFN-FF, and the choice of training
conformer moves the result by as much as the method change does. The details below
matter, because the failure is informative and the fit itself is not at fault.

## Why this was tried

A conformer search needs a cheap surface to **explore** on and an accurate one to
**rank** on, and the two disagree. Measured on this peptide (Obsidian note
*Beschreibung und Rangfolge von Konformeren*), GFN-FF and GFN2 correlate at
**r = −0.32 … −0.46** within a search cycle. QMDFF@GFN2 — a force field fitted to the
ranking method itself — was the proposed structural fix: one surface for both jobs,
and a term decomposition of the *right* surface.

Curcuma can now build it: `-qmdfffit mol.xyz -method gfn2` fits the full QMDFF energy
expression to a GFN2 Hessian ([QMDFF.md](QMDFF.md)).

## Ensembles

| | A: ConfScan test set | B: WEKLQ peptide |
|---|---|---|
| file | `test_cases/confscan/input.xyz` | `…/WEKLQ/input.confsearch.20260727_102718/input.cumulative.opt.accepted.xyz` |
| conformers | 44 | 139 |
| atoms | 114 | 107 |
| GFN2 spread | 38.7 kJ/mol | 93.2 kJ/mol |

B is the system from the note. The QMDFF parameter set is fitted **once**, at one
conformer, then applied unchanged to the whole ensemble — what a conformer search
would do.

## Result 1 — ranking

| ensemble | method | Pearson r | Spearman ρ | MAD [kJ/mol] | spread | rank of the GFN2 minimum |
|---|---|---|---|---|---|---|
| A (44) | GFN-FF | 0.481 | 0.452 | **11.2** | 44.6 | 22 / 44 |
| A (44) | QMDFF@GFN2 | **0.772** | **0.774** | 14.4 | 69.0 | 5 / 44 |
| B (139) | GFN-FF | **0.332** | **0.288** | **28.2** | 169.8 | 1 / 139 |
| B (139) | QMDFF@GFN2 | 0.217 | 0.115 | 66.0 | 219.3 | 1 / 139 * |

\* **This entry is self-fulfilling and must not be read as a success.** On B the fit
was placed at the GFN2 minimum, and QMDFF takes r₀/θ₀ from its training structure, so
every bonded term is exactly zero there by construction. The force field is *defined*
to put its own training conformer first. On A the fit sat at conformer 0 while the
minimum was 17, so the 5/44 there is a real result.

**On the peptide QMDFF@GFN2 loses to GFN-FF on every non-trivial measure**: rank
correlation collapses from 0.288 to 0.115 (i.e. essentially none), the MAD more than
doubles, and the spread is over-estimated by a factor 2.4.

Ensemble A was simply an easy case. Note that GFN-FF already correlates at r = +0.48
there, whereas the peptide is the regime the note measured at r = −0.32 … −0.46.
Generalising from A alone would have been wrong.

## Result 2 — the fit is not what fails

The obvious suspects were checked and cleared:

| check | B (WEKLQ) |
|---|---|
| Hessian fit R² | **0.982** (residual 0.134) |
| geometry is a GFN2 stationary point | yes — `tr_content` = 1.4e-4 |
| force constants at the non-negativity floor | 17 of 364 |
| conformers with a different H-attachment pattern than the training structure | **0 of 139** |

That last row rules out the note's own protomer worry: no proton has migrated, so
every conformer shares the training structure's bond list. The parameter set is a
good local model — it reproduces 98 % of the GFN2 Hessian's variance — and it still
cannot rank an ensemble spanning 93 kJ/mol.

This is the note's first structural caveat, confirmed quantitatively: **QMDFF is a
curvature expansion about one structure, and a conformer search operates exactly where
that expansion is no longer valid.** A better fit cannot fix it; only refitting per
basin could.

## Result 3 — the training conformer matters as much as the method

`scripts/qmdff_reference_sensitivity.py` repeats the whole procedure on ensemble A with
the fit placed at four conformers spanning the GFN2 energy range:

| training conformer | its own GFN2 rank | Spearman ρ | MAD | rank of the GFN2 minimum |
|---|---|---|---|---|
| 17 | 1 / 44 | 0.703 | 16.3 | 1 / 44 * |
| 19 | 15 / 44 | **0.751** | 12.1 | 8 / 44 |
| 37 | 30 / 44 | 0.619 | 7.1 | **28 / 44** |
| 10 | 44 / 44 | 0.508 | 17.9 | 5 / 44 |

\* again self-fulfilling — the fit is at the minimum.

- **ρ ranges 0.508 … 0.751 — a spread of 0.24.** The entire claimed GFN-FF → QMDFF
  improvement on this ensemble (0.45 → 0.77) is *within the range produced by moving
  the training conformer*. The method change and the arbitrary choice of where to fit
  are of comparable size, which means ensemble A's headline number was not a property
  of the method.
- **Fitting at the minimum is not optimal.** Conformer 19 (rank 15) beats the minimum
  (0.751 vs 0.703). There is no rule available in practice anyway — a search does not
  know the minimum in advance.
- **The worst case is actively harmful.** Fitting at conformer 37 puts the true
  minimum at rank 28/44, i.e. *worse* than GFN-FF's 22/44.

A near/far split by heavy-atom RMSD from the training structure gave 9.7 vs 6.1 kJ/mol
— nominally better far away, which is the opposite of the expectation. That metric
removes a per-subset energy offset, so a subset that is uniformly shifted looks
accurate; it therefore cannot separate "well described" from "wrong by a constant" and
**no conclusion is drawn from it**. A proper test needs an offset-free error measure.

## Result 4 — interpretability of the terms

Exact decomposition Var(E) = Σᵢ Cov(Eᵢ, E):

| term | A: GFN-FF | A: QMDFF@GFN2 | B: GFN-FF | B: QMDFF@GFN2 |
|---|---|---|---|---|
| hydrogen bonds | 46.9 % | 15.7 % | 21.9 % | −6.7 % |
| electrostatics | 3.5 % | 16.5 % | **37.3 %** | 25.4 % |
| dispersion | 2.6 % | −3.7 % | 21.3 % | 3.0 % |
| repulsion | 24.3 % | 26.8 % | 8.1 % | **28.8 %** |
| angle bend | 16.2 % | 22.1 % | 0.6 % | 18.6 % |
| torsion | 9.7 % | 20.6 % | 7.3 % | 22.6 % |
| bond stretch | 0.6 % | 2.0 % | 4.6 % | 7.7 % |

The **GFN-FF column for B reproduces the note's peptide finding** independently:
non-covalent terms (electrostatics + HB + dispersion) carry **80.5 %**, torsions
**7.3 %**. That is a useful consistency check on both the note and this pipeline.

The QMDFF columns look more balanced — torsions roughly triple, to ~21–23 %. It is
tempting to conclude that "torsions describe the wrong physics" is method-dependent.
**That conclusion is not supported for B**, because a term decomposition is only
meaningful if the surface being decomposed is the one you care about, and on B the
QMDFF surface correlates with GFN2 at ρ = 0.115. Decomposing it describes QMDFF, not
GFN2. On A, where ρ = 0.77, the statement has some support — but see Result 3 for how
much of A's ρ is an artefact of the training choice.

Independently of accuracy, a decomposition attributes energy to whichever *functional
form* can absorb it. QMDFF has no cooperative hydrogen-bond machinery, so a GFN2-fitted
QMDFF pushes what GFN-FF books as "hydrogen bond" into its Born–Mayer repulsion and
point-charge electrostatics — visibly, in both ensembles. On B the QMDFF HB term even
comes out **anti-correlated** (−6.7 %).

## What this means for the use case

The structural argument in the note is sound: exploring and ranking on one surface is
the right goal. This measurement says the **single-fit** realisation does not reach it
for a flexible peptide.

Options, roughly in order of cost:

1. **Refit per basin** — the note's own proposal, and now the only remaining candidate.
   Cost: one GFN2 Hessian (6N = 642 gradient calls here) per fitting point. Untested.
   The sensitivity study suggests the *number* of fits matters more than their placement.
2. **Use QMDFF@GFN2 only locally** — as an MD/relaxation surface within a basin
   (`-md_method`), never as the ranking or dedup criterion. This is consistent with
   both the note's "select and dedup on the ranking surface" rule and with the present
   data.
3. **Accept GFN-FF for exploration** and put the effort into cheaper GFN2 single points
   on the ranking side instead.

Option 1 is the only one that tests the original idea. It is a well-defined next
experiment: fit at *k* well-separated conformers, assign each conformer to the nearest
fit, and re-measure ρ as a function of *k*.

## Caveats

- Two ensembles, one molecule each, one reference method. B is the relevant system;
  A is now known to be unrepresentative.
- The fitted force field is not a stationary point of itself; the small test case
  reports 14 imaginary modes. Frequencies from it are not meaningful without
  re-optimisation.
- Parameters are frozen to one geometry and one method (charges, C6, r₀, θ₀). Inherent
  to QMDFF.
- With n = 139, the 95 % CI on ρ = 0.115 is roughly [−0.05, 0.28] and on ρ = 0.288
  roughly [0.13, 0.43]; the two overlap slightly, so "QMDFF is worse" is better stated
  as "QMDFF is not better, and its MAD is unambiguously worse (66 vs 28 kJ/mol)".

## Reproducing

```bash
# ranking + term decomposition
scripts/qmdff_conformer_evaluation.py <ensemble>.xyz --curcuma release/curcuma \
    --workdir eval --threads 8

# sensitivity to the training conformer
scripts/qmdff_reference_sensitivity.py <ensemble>.xyz --curcuma release/curcuma \
    --workdir refscan --n-refs 4 --threads 4
```

Both write JSON next to their work directory. `--skip-fit` reuses an existing fit and
reports which conformer it came from.
