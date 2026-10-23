# polymerbuild: charge of the chain (Oct 2026)

Status: AI-generated, machine-tested; one chain (PDMAEMA, `(AB)5`, n=1 per number) measured, no human review.

## What polymerbuild does with a charge

- Fragment charge: `** Charge = N **` in the XYZ comment (native format, the comment must contain the word `Curcuma`).
  `src/molecules/pdmaema_h.xyz` (unit with N-H) carries `Charge = 1`.
- The formal charge of the result is the sum over the sequence; it goes into the header of `polymer_A_final.xyz`.
- `-charge_model gfnff|gfn1|gfn2` writes `polymer_A_final.charges` (atom index, element, charge, sum = formal charge).
- Every optimisation step of the build still runs with charge 0, a warning says so.
- `polymer_A_final.groups.json`: one EEQ charge group per charged monomer (see below), for `-gfnff.charge_groups`.

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

## EEQ charge groups: `-gfnff.charge_groups FILE`

GFN-FF fixes the charge per covalently connected fragment (one constraint row in the EEQ matrix, Lagrange multiplier).
A chain is one fragment, so the net charge is spread by electronegativity. With `charge_groups` every user group gets its own
row `sum(q_i) = Q_g`; the atoms outside all groups form one more group with the remaining charge. File:
`{"groups": [{"atoms": [0-based indices], "charge": 1.0}, ...]}`. It reuses the override of the `frag_charge_model`
variants (`GFNFF::setFragmentOverride`), so Phase 1 and Phase 2 and the gradient use the same constraint matrix.
Only inside one connected molecule; groups that overlap, fall outside the molecule or leave an inconsistent total are
rejected with an error and the run continues without groups. Topology cache is bypassed, the GPU method falls back to the CPU
(`requestsCpuFallback`; the GPU edit is not compiled in the CPU build).

Measured ((AB)5, 267 atoms, +5, one group = one protonated unit of 27 atoms, charge +1):

| | unit charge sum (mean [min, max]) | neutral units (sum) | per-atom RMSD to gfn2 (e) |
|---|---|---|---|
| gfn2 | 0.944 [0.882, 0.983] | +0.28 | - |
| gfn1 | 0.933 [0.882, 0.987] | +0.33 | 0.078 |
| gfnff, no groups | 0.505 [0.366, 0.679] | +2.48 | 0.062 |
| gfnff, unit groups | 1.000 (constraint) | 0.00 | 0.062 |

The group sum is set by construction; the distribution inside a group is not improved (RMSD unchanged). The result depends
on the group definition: for the 13-atom ammonium head group (N, H, three C, their H) the protonated side chain gets 0.48 e,
for the whole side chain (19 atoms) at +1.0 / +0.84 e it gets 1.00 / 0.84 e and the rest 0.71 / 1.23 e (gfn2: 0.84 / 1.03 e).
The formal +1 is 6% above the QM unit charge 0.944; edit the group charge in the file for a different value.

Checks: gradient against a central finite difference of the own energy (perturbed (AB)3, 18 components): max 1.2e-5 Eh/A
with groups, 1.1e-5 without, largest component 0.47 Eh/A. `-opt` converged (754 steps), NVE MD over 1 ps: Etot drift
+0.0022 Eh with groups, +0.0031 Eh without (n=1 each). Test `gfnff_charge_groups` (script `check_charge_groups.py`).
`cli_simplemd_14_gfnff_react_h_recombination` fails with and without this change (checked by stashing it, n=1).

## Not done

- Charge in the optimisation steps of the build (would need the partial charge of each growing chain at every step and every join).
- Groups in the optimisation steps of the build (they would need the groups of each growing chain), several fragments, GPU.
- Other charge models than the three above, other polyelectrolytes, solvated chains.
