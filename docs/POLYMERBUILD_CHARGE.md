# polymerbuild: charge of the chain (Oct 2026)

Status: AI-generated, machine-tested; one chain (PDMAEMA, `(AB)5`, n=1 per number) measured, no human review.

## What polymerbuild does with a charge

- Fragment charge: `** Charge = N **` in the XYZ comment (native format, the comment must contain the word `Curcuma`).
  `src/molecules/pdmaema_h.xyz` (unit with N-H) carries `Charge = 1`.
- The formal charge of the result is the sum over the sequence; it goes into the header of `polymer_A_final.xyz`.
- `-charge_model gfnff|gfn1|gfn2` writes `polymer_A_final.charges` (atom index, element, charge, sum = formal charge).
- Every optimisation step of the build still runs with charge 0, a warning says so.

## Measured: where the charge sits ((AB)5, 267 atoms, +5, built structure, `-sp -charge 5 -dump_atomic_charges`)

| model | protonated side chain | neutral side chain | rest (backbone, methyls) |
|---|---|---|---|
| gfn2 | +0.84 e | -0.04 e | +1.03 e |
| gfn1 | +0.84 e | -0.01 e | +0.85 e |
| gfnff (EEQ) | +0.19 e | +0.15 e | +3.34 e |

GFN-FF treats the connected chain as one fragment and spreads the net charge by electronegativity equalisation; it cannot
localise it on the ammonium groups. Per-atom correlation over 267 atoms: gfnff-gfn1 r=0.92, gfnff-gfn2 r=0.94, gfn1-gfn2
r=0.93 (Spearman 0.68 / 0.76 / 0.89). Per-side-chain sums (10 chains): gfn1-gfn2 r=1.00, gfnff against either r=0.19-0.20.
The atom-wise agreement hides that the group sums disagree. Script: `scripts/polymer_charge_compare.py`.

## Does the charge change the built geometry? ((AB)3, 161 atoms, optimised from the built structure)

| | N+ C-N / N-H (A) | heavy-atom RMSD to the built structure |
|---|---|---|
| built (gfnff, charge 0) | 1.506 / 1.049 | - |
| gfnff, charge 0 | 1.506 / 1.049 | 0.26 A |
| gfnff, charge 3 | 1.509 / 1.048 | 1.61 A |
| gfn2, charge 3 | 1.479 / 1.027 | 1.05 A |

Locally the charge changes nothing in GFN-FF; the gap to GFN2 is the model, not the charge. The global conformation does
move with the charge (1.6 A), so a chain built neutral is a different conformer than one built charged.

## Not done

- Charge in the optimisation steps of the build (would need the partial charge of each growing chain at every step and every join).
- Localised charges for GFN-FF: `GFNFF::setCharges()` exists in C++, no CLI access; not tried.
- Other charge models than the three above, other polyelectrolytes, solvated chains.
