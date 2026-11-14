# SimpleMD + GFN-FF: step time against system size (Oct 29, 2026)

AI-generated measurement (Claude), not reviewed by a human; full report with all tables and raw data outside the repository:
`~/curcuma_runs/probe_large/steptime_diag/ANALYSIS_2026-10-29.md`.

Setup: release/curcuma at 106552c3, CPU, 8 threads on cores 0-7 (cores 16-23, the SMT siblings sharing the 96 MB L3, idle), flags of the
solvated-polymer production (`-gfnff.charge_groups`, `cn_zero_guard true`, `topology_mode constant`), droplets carved from the equilibrated
L1 frame, 200 steps, n = 3 per size.

| atoms (chain + water) | 1501 | 3003 | 4505 | 5747 |
|---|---:|---:|---:|---:|
| ms/step | 23.5 | 66.2 | 121.8 | 192.7 |

- Log-log exponent 1.55 +- 0.02 (n = 12); pure water 1.53 +- 0.03 (n = 12), so chain and charge groups do not change it. T = aN + bN^2 fits within 4.3 %; the N^2 part is 60 % of the step at 3000 and 74 % at 5750 atoms.
- 3003 -> 5747 atoms is x2.91, not the x3.2 quoted earlier (that compared runs under different load). Load matters: runs with the SMT siblings busy were 1.5-1.8x slower.
- Cause: all-pairs work without cutoff (as in the reference): EEQ Phase-2 matrix and PPCG matvecs (dense), Coulomb, distance array; the D4 pair list (60 Bohr) still holds 78 % of all pairs for a 52-54 A droplet. Share of the excess over linear (n = 1 instrumented run per size): D4 weights and C6 refresh 25 %, PPCG iterations 20 %, D4 energy 10 %, zeroing the constraint rows of A 9 %, EEQ matrix build 8 %, localinv rebuild 8 %, Coulomb 6 %.
- Cache: the D4 list (90 -> 344 MB) and the half of A read by the PPCG (81 -> 132 MB) leave the L3 between the sizes (supported by the data, n = 1 per size).
- The constraint rows give no dense Schur complement on the auto path; they cost only zeroing (11 ms at 5747) and 196 MB. `cn_zero_guard` and the groups do not change the scaling. Solver variants: `eeq_ppcg_precond diag` 17 % faster at 3003 atoms, equal at 5747; cholesky 3x and pcg 39x slower at 5747.
- Production stayed flat (299-308 ms/step in every tenth of 10.5 h): periodic rebuilds are about 8 % of the step.

Fix ideas, none implemented: constraint block only for Schur/LU (-6 %, 0.5-1 day); pAp > 0 check in the PPCG instead of the serial definiteness check of the localinv rebuild (-5 %, 1 day); fuse the D4 C6 refresh into the energy loop (-10 to -15 %, 2-4 days, accuracy against the reference to be shown again); linear scaling only with cutoffs/Ewald/FMM, which changes energies against the reference (decision of the operator).

Not tested: GPU, other T, more than 5747 atoms, the 32 MB L3 of the second CCD.
