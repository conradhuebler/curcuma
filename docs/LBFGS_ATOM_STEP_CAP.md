# `-lbfgs_max_atom_step`: per-atom step cap for the lbfgspp optimiser

AI-generated (Claude, Oct 30, 2026), machine-tested only; human production testing pending.

## Problem

`-opt` (default `auto` = LBFGSpp) tore a C=O bond in a solvated system of 5747 atoms ((AB)5 PDMAEMA + 1825 H2O + 5 Cl-,
replica L1): 1.21 -> 1.99 A in iteration 4, 2.41 A after 200 iterations. Diagnosis (n = 1 structure, deterministic):
the initial inverse-Hessian scale H0 = (s'y/y'y) I is learned from the soft water modes and is about 10 times too
large for C=O; the Strong-Wolfe line search looks at the total energy only, the water gains hide the local break
(+594 kJ/mol for the two atoms, -2740 kJ/mol for the rest), and the bounded GFN-FF bond term has no restoring force from
about 1.9 A on. The same overshoot occurs in replica L2 (1.87 A) and the 3071-atom probe (1.73 A) but they recover.
Full report: `~/curcuma_runs/probe_large/diag_L1_bond/ANALYSIS_2026-10-29.md` (not in the repository).

## Switch

| | |
|---|---|
| `-lbfgs_max_atom_step <A>` | largest displacement of any single atom in one L-BFGS step, Angstrom; **0 = off (default)** |
| Where | `LineSearchBacktrackingCapped` in `src/capabilities/lbfgspp_capped_linesearch.h`, a copy of LBFGSpp's `LineSearchBacktracking` (`external/LBFGSpp` is a fetched, untracked dependency) |
| Rule | trial step length <= cap / max_i |d_i|, also for the expansion of the Wolfe branch; at the cap a step that wants to grow is accepted as is |
| Cap off | same arithmetic as the original; coordinates bit-identical (below) |
| State | `thread_local` value set by `LBFGSppOptimizer` before every step |

Native `lbfgs`/`rfo` use a global trust radius on the step norm, ANCopt clips each ANC component (1 Bohr); both are other mechanisms.

## What was tested

- L1 droplet (5747 atoms), 30 iterations: cap 0 and the pre-change binary give bit-identical final coordinates (max |dx| = 0) and the same energies (-661.392110 Eh); cap 0.2 A keeps C241-O242 at 1.214 A (cap 0: 2.392 A), energy -661.399470 Eh.
- Earlier, with the diagnostic patch (initial-step cap only): 200 iterations, c = 0.2/0.1 A, L1 -671.293/-671.292 Eh against -671.139 Eh uncapped; L2 within 0.005 Eh of uncapped (n = 1 each).
- `ctest -R opt_lbfgs_atom_step_cap`: acetic-acid dimer, coordinates perturbed by 0.12 A: cap 0 == no flag, uncapped step 0.266 A, capped <= 0.05 A, energies differ by 2.1e-4 Eh (tolerance 1e-3 set after that run).
- `ctest -R cli_curcumaopt_` passes (6 of 6 enabled).

## Not tested

Other structures and cap values, MOR41/GMTKN55 optimisations with the cap on, the cap with the Armijo/Wolfe variants (`-lbfgs_line_search 1/2`), `topology_mode auto`, GPU, a full-length 5747-atom optimisation to convergence. The default stays 0; 0.2 A is suggested for large solvated systems only.

## Related findings (not fixed)

`alpha=` printed at verbosity 2 is always 1.0; the `lbfgs_line_search` help lists "4=Backtracking" but LBFGSpp accepts 1-3; `lbfgs_max_line_search` has two defaults (2 in the driver, 20 in the registry; the registry value reaches the solver); the "Gradient Norm [Eh/Bohr]" column is in Eh/A.
