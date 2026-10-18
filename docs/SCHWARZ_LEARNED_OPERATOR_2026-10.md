# Learned two-level Schwarz operator for the inverse RI metric: Python experiment, Oct 2026

Status: 🤖 AI-generated, ⚙️ machine-tested (numpy/PySCF prototype, nothing in the curcuma engine). Learned runs: 3 seeds
each; everything else one run (n=1), one machine (24 cores, shared with other jobs, so wall times are upper bounds).
Human review pending. Follows [LOCAL_MP2_FEASIBILITY_2026-10.md](LOCAL_MP2_FEASIBILITY_2026-10.md) (Q4: the inverse
metric is not local) and [LOCAL_CCSD_PROTOTYPE_2026-10.md](LOCAL_CCSD_PROTOTYPE_2026-10.md).

## Question

Köhler and Rheinbach (TU Freiberg, arXiv:2606.18525v2, read from the arXiv PDF) learn a softmax-free operator
`M^-1 = Phi Q0 K0^T Phi^T + sum_i R_i^T D_i^1/2 Q_i K_i^T D_i^1/2 R_i` (local low-rank blocks on overlapping
subdomains plus a coarse low-rank block, their eq. 17) that approximates the inverse 1D Poisson matrix, and report
lower loss than a global `Q K^T` at 8x fewer parameters, with the global model degrading at large learning rates
(their Table 1). Does this structure help for a quantum-chemistry operator? Target: `G ~ J^-1`, `J` the Coulomb metric
of def2-universal-JKFIT (curcuma's file).

## Setup and decisions

- Scripts: `scripts/schwarz_learn/` (`sl_common.py` operators/evaluation, `sl_learn.py` models/Adam, `check_grad.py`,
  `run_free.py`, `run_transfer.py`, `run_lrsweep.py`, `rank_needed.py`, `tables.py`). Geometries, aux parser and cached
  RI-JK RHF/def2-SVP from `scripts/local_mp2/lmp2_common.py` (alkanes C4-C20, water clusters, acenes; C6H14 added).
- **How G is used:** only as "apply an inverse", fitting coefficients `c = G v` with `v = (P|rho)` (RI-J) or
  `v = (P|ia)` (RI-MP2), then the robust (Dunlap) expressions `E_J = (rho|c) - 1/2 c^T J c` and
  `(ia|jb) = A_ia c_jb + c_ia A_jb - c_ia^T J c_jb`. Their error is second order, `-1/2 dc^T J dc`, so no symmetric
  factor `J^-1/2` is learned. Non-robust ("plain") expressions with the symmetrised G are reported as well.
- **Loss** (analogue of the paper's sample-weighted wMSE): mean relative J-norm fit error
  `dc^T J dc / (v^T J^-1 v)` over right-hand sides = twice the relative robust RI-J error per right-hand side.
  Pool per molecule: `(P|rho)`, `(P|rho_A)` per atom, all `(P|ia)` (frozen core); 20 % of the `(P|ia)` held out
  ("val."). `(P|rho)` is in the training pool of `run_free.py` (1 of ~10^3 columns), not in the transfer tests.
- **Learning:** Coulomb-normalised variables (`Jn = S^-1 J S^-1`, `S = diag(J)^1/2`), hand-written Adam (no torch in the
  venv, nothing installed), as in the paper: beta 0.9/0.999, no weight decay, ReduceLROnPlateau (0.5, patience 200),
  batch 256; lr 1e-3 (the paper's default) and 3000 steps (paper 2000). Initial scale `0.01/sqrt(n)` for the product
  (the paper's h-based scale has no analogue). Gradients checked against central differences (`check_grad.py`,
  max relative deviation 3e-7).
- **Subdomains:** atom-centred, all aux functions on atoms within `R_sd` = 4 or 6 Bohr; partition of unity
  `D_i = diag(1/m_P)` as in the paper. Coarse space: one column per atom, the charge integrals of its functions
  (multipole coarse space; the paper's interface hats have no analogue). Coarse level either learned `Q0 K0^T`
  (rank = natm, as the paper's `r0 = N-1`) or fixed Galerkin `Phi (Phi^T J Phi)^-1 Phi^T`.
- **Comparators:** exact `J^-1`; global rank k (truncated eigen of `J^-1`, largest eigen of `J`, data-optimal rank k =
  closed-form reduced-rank regression for the loss, the best any global rank-k operator can do on the training data);
  classical two-level Schwarz with exact local inverses (`asm` plain additive, `asm-pu` the paper's weighting,
  `ras-core` the core rows of each local inverse; coarse none/additive/balancing); closed-form local reduced-rank
  blocks (`local RRR`, same structure and budget as the learned one); learned global `Q K^T` at the Schwarz rank bound
  (natm*r_l + n0, "rank-matched", as the paper) and at the same parameter count ("param-matched"); identity and Jacobi
  for the preconditioner columns. PCG (symmetrised G) and right-preconditioned GMRES for `J x = (P|rho)` to 1e-10.
- Parameter counts: trainable parameters, or stored numbers for closed forms, divided by the dense Cholesky factor
  `naux(naux+1)/2`; flops per applied vector divided by `2 naux^2`.

## Result 1: no structured or learned approximate inverse is accurate enough for one-shot fitting

Target: `|dE| < 1e-6` Eh/atom (vault falsifier for local fits). Reached by no learned and no Schwarz operator in this study;
only the exact inverse and closed-form global data-optimal operators reach it, at 93 % of full rank when `(P|rho)` is in
the fit set (Result 2) and 97-99.9 % when it is not (table below).

Rank a global operator needs (`rank_needed.py`, best possible global rank-k operator, fitted on all AO-pair right-hand
sides, tested on `(P|rho)`; val. = held-out `(P|ia)` with P_k fitted on `(P|ia)`):

| system | naux | k needed, dE_J < 1e-6 Eh/atom | k/naux | val. err (ia) at k=0.5 naux | at 0.9 naux | at 0.98 naux | dE_J/atom at 0.9 naux |
|---|---|---|---|---|---|---|---|
| alkane4 | 480 | 475 | 0.990 | 2.3e-03 | 2.7e-06 | 2.8e-08 | -6.2e-04 |
| alkane8 | 924 | 919 | 0.995 | 1.8e-03 | 1.5e-06 | 1.1e-08 | -1.3e-03 |
| alkane12 | 1368 | 1359 | 0.993 | 1.8e-03 | 1.4e-06 | 1.1e-08 | -3.5e-03 |
| water4 | 452 | 451 | 0.998 | 8.9e-04 | 3.0e-07 | 6.3e-09 | -2.7e-03 |
| water8 | 904 | 903 | 0.999 | 8.1e-04 | 4.0e-07 | 7.5e-09 | -2.6e-02 |
| acene2 | 894 | 867 | 0.970 | 2.5e-03 | 8.2e-07 | 1.2e-08 | -2.5e-04 |
| acene3 | 1230 | 1193 | 0.970 | 2.3e-03 | 6.2e-07 | 9.2e-09 | -3.2e-04 |

Any operator with rank bound below this (the Schwarz structure has rank <= natm*r_l + n0) is limited by this curve on
the same data. **The fitting operator is essentially full rank** on chemically relevant right-hand sides: the
Frobenius-optimal truncation of `J^-1` (largest eigenvalues of `J^-1`, i.e. the near-null space of `J`) is the worst
possible choice, 6e-5 Frobenius error with -20 Eh/atom RI-J error on C8H18.

## Result 2: per-molecule comparison (free parameters per subdomain, as in the paper)

Seed spread in parentheses (min..max over 3 seeds); "indef." = PCG breakdown (symmetrised G indefinite or singular);
"n.c." = not converged in 2000 (PCG) / 400 (GMRES) iterations.

#### alkane8: natm 26, naux 924, E_J 527.940 Eh, E_corr -1.159797 Eh, 3000 Adam steps, lr 0.001, Cholesky of J 2.6 ms; local density fit R=8 (outside the operator class): dE_J -8.7e-07 Eh

| operator | params / dense | flops / dense | val. rel. J-norm err | Frob. err (norm.) | dE_J robust (Eh/atom) | dE_MP2 robust (Eh/atom) | dE_MP2 plain (Eh/atom) | PCG it. | GMRES it. | train s |
|---|---|---|---|---|---|---|---|---|---|---|
| exact J^-1 (Cholesky) | 1.000 | 1.000 | 0.0e+00 | 0.0e+00 | 0.0e+00 | -1.0e-16 | 9.0e-13 | 2 | 2 | - |
| none (identity) | 0.000 | 0.000 | 9.6e+03 | 1.0e+00 | -3.1e+07 | - | - | n.c. | n.c. | - |
| Jacobi diag(J)^-1 | 0.002 | 0.001 | 2.5e+01 | 1.0e+00 | -2.5e+04 | - | - | n.c. | n.c. | - |
| global trunc-eig J^-1, rank 234 | 0.506 | 0.506 | 1.0e+00 | 6.1e-05 | -2.0e+01 | 4.5e-02 | 4.5e-02 | indef. | n.c. | - |
| global large-eig J, rank 234 | 0.506 | 0.506 | 2.0e-01 | 1.0e+00 | -5.2e-01 | 2.5e-03 | 2.5e-03 | indef. | n.c. | - |
| global data-optimal, rank 234 | 0.506 | 0.506 | 4.7e-02 | 1.0e+00 | -5.3e-01 | -6.4e-05 | -6.4e-05 | indef. | n.c. | - |
| global trunc-eig J^-1, rank 858 | 1.857 | 1.857 | 3.5e-01 | 2.2e-07 | -1.9e+01 | 3.5e-02 | 3.5e-02 | indef. | 144 | - |
| global large-eig J, rank 858 | 1.857 | 1.857 | 5.6e-04 | 1.0e+00 | -1.6e-06 | 4.0e-06 | 4.0e-06 | indef. | 153 | - |
| global data-optimal, rank 858 | 1.855 | 1.857 | 4.6e-07 | 4.3e-01 | -7.0e-08 | 4.1e-10 | 4.1e-10 | indef. | 141 | - |
| global data-optimal, rank 35 | 0.076 | 0.076 | 6.6e-01 | 1.0e+00 | -1.4e+01 | 5.9e-03 | 5.9e-03 | indef. | n.c. | - |
| global data-optimal, rank 141 | 0.305 | 0.305 | 1.8e-01 | 1.0e+00 | -1.1e+01 | -8.7e-04 | -8.7e-04 | indef. | n.c. | - |
| global data-optimal, rank 104 | 0.225 | 0.225 | 2.9e-01 | 1.0e+00 | -1.2e+01 | -1.1e-03 | -1.1e-03 | indef. | n.c. | - |
| global data-optimal, rank 416 | 0.899 | 0.900 | 3.7e-03 | 1.0e+00 | -1.1e-01 | -3.0e-06 | -3.0e-06 | indef. | n.c. | - |
| Schwarz exact asm (R=4.0) | 0.870 | 0.866 | 4.4e+01 | 1.0e+00 | -7.3e+03 | -4.9e+02 | -4.3e+00 | 142 | 116 | - |
| Schwarz exact asm-pu (R=4.0) | 0.870 | 0.866 | 1.6e+04 | 7.1e-01 | -1.1e+08 | -2.9e+09 | -1.4e+03 | 1203 | 225 | - |
| Schwarz exact asm-pu + coarse bal (R=4.0) | 0.929 | 0.931 | 1.8e+01 | 7.1e-01 | -1.8e+01 | -2.5e+02 | -1.2e+00 | 287 | 148 | - |
| Schwarz exact ras-core + coarse bal (R=4.0) | 0.929 | 0.283 | 1.7e+01 | 5.7e-01 | -3.3e+00 | -2.5e+02 | -5.7e-01 | indef. | 239 | - |
| Schwarz exact asm (R=6.0) | 6.837 | 6.830 | 2.1e+02 | 8.0e+00 | -1.0e+04 | -4.2e+03 | -1.4e+01 | 27 | 26 | - |
| Schwarz exact asm-pu (R=6.0) | 6.837 | 6.830 | 1.8e+05 | 4.3e-01 | -8.3e+08 | -3.5e+11 | -1.1e+04 | 360 | 92 | - |
| Schwarz exact asm-pu + coarse bal (R=6.0) | 6.896 | 6.895 | 5.9e+01 | 4.3e-01 | -8.1e+02 | -1.8e+04 | -5.5e+00 | 112 | 65 | - |
| Schwarz exact ras-core + coarse bal (R=6.0) | 6.896 | 0.575 | 1.5e+01 | 3.6e-01 | -8.0e+00 | -3.9e+02 | -6.8e-01 | indef. | 89 | - |
| local RRR r=8 + Galerkin (R=4.0) | 0.153 | 0.000 | 1.6e+04 | 9.8e-01 | -1.1e+08 | -2.9e+09 | -1.4e+03 | n.c. | n.c. | - |
| learned Schwarz R=4.0 r_l=8 coarse=learned (3 seeds) | 0.156 | 0.080 | 6.4e-02 (6.3e-02..6.4e-02) | 1.0e+00 | -9.5e-01 | 1.9e-04 | 8.5e-06 | indef. | n.c. | 57 |
| learned Schwarz R=4.0 r_l=8 coarse=galerkin (3 seeds) | 0.153 | 0.080 | 5.7e-02 (5.6e-02..5.7e-02) | 1.0e+00 | -1.1e+00 | 7.3e-04 | 4.5e-04 | indef. | n.c. | 47 |
| learned Schwarz R=4.0 r_l=32 coarse=learned (3 seeds) | 0.614 | 0.310 | 1.3e-02 (1.3e-02..1.3e-02) | 1.0e+00 | -4.0e-01 | 7.2e-05 | 2.8e-05 | indef. | n.c. | 200 |
| learned Schwarz R=4.0 r_l=32 coarse=galerkin (3 seeds) | 0.611 | 0.310 | 1.3e-02 (1.3e-02..1.3e-02) | 1.0e+00 | -3.3e-01 | 7.3e-05 | -1.5e-05 | indef. | n.c. | 106 |
| learned Schwarz R=6.0 r_l=8 coarse=learned (3 seeds) | 0.453 | 0.229 | 5.8e-02 (5.8e-02..5.9e-02) | 1.0e+00 | -8.1e-01 | 9.2e-05 | 1.6e-04 | indef. | n.c. | 117 |
| learned Schwarz R=6.0 r_l=8 coarse=galerkin (3 seeds) | 0.450 | 0.229 | 5.2e-02 (5.1e-02..5.2e-02) | 1.0e+00 | -9.2e-01 | 5.5e-04 | 6.1e-04 | indef. | n.c. | 101 |
| learned Schwarz R=6.0 r_l=32 coarse=learned (3 seeds) | 1.803 | 0.905 | 1.4e-02 (1.4e-02..1.4e-02) | 1.0e+00 | -3.2e-01 | 5.1e-05 | 6.1e-05 | indef. | n.c. | 108 |
| learned Schwarz R=6.0 r_l=32 coarse=galerkin (3 seeds) | 1.800 | 0.905 | 1.3e-02 (1.2e-02..1.3e-02) | 1.0e+00 | -2.9e-01 | 4.6e-05 | 7.1e-05 | indef. | n.c. | 104 |
| learned global rank 234 (rank-matched, 3 seeds) | 1.012 | 0.506 | 4.9e-02 (4.9e-02..4.9e-02) | 1.0e+00 | -4.8e-01 | 2.7e-06 | 2.5e-05 | indef. | n.c. | 102 |
| learned global rank 35 (param-matched, 3 seeds) | 0.151 | 0.076 | 6.7e-01 (6.7e-01..6.7e-01) | 1.0e+00 | -1.4e+01 | 5.4e-03 | 5.5e-03 | indef. | n.c. | 61 |
| learned global rank 858 (rank-matched, 3 seeds) | 3.710 | 1.857 | 1.1e-02 (1.0e-02..1.1e-02) | 1.0e+00 | -2.9e-01 | 7.3e-05 | 1.0e-04 | indef. | n.c. | 111 |
| learned global rank 141 (param-matched, 3 seeds) | 0.610 | 0.305 | 1.8e-01 (1.8e-01..1.8e-01) | 1.0e+00 | -1.1e+01 | -7.5e-04 | -6.7e-04 | indef. | n.c. | 24 |
| learned global rank 104 (param-matched, 3 seeds) | 0.450 | 0.225 | 2.9e-01 (2.9e-01..2.9e-01) | 1.0e+00 | -1.2e+01 | -9.2e-04 | -8.5e-04 | indef. | n.c. | 32 |
| learned global rank 416 (param-matched, 3 seeds) | 1.799 | 0.900 | 1.2e-02 (1.2e-02..1.3e-02) | 1.0e+00 | -3.4e-01 | 6.6e-05 | 9.9e-05 | indef. | n.c. | 46 |

#### water4: natm 12, naux 452, E_J 271.979 Eh, E_corr -0.814691 Eh, 3000 Adam steps, lr 0.001, Cholesky of J 0.7 ms; local density fit R=8 (outside the operator class): dE_J -1.6e-07 Eh

| operator | params / dense | flops / dense | val. rel. J-norm err | Frob. err (norm.) | dE_J robust (Eh/atom) | dE_MP2 robust (Eh/atom) | dE_MP2 plain (Eh/atom) | PCG it. | GMRES it. | train s |
|---|---|---|---|---|---|---|---|---|---|---|
| exact J^-1 (Cholesky) | 1.000 | 1.000 | 0.0e+00 | 0.0e+00 | 1.4e-14 | 4.6e-17 | 1.9e-13 | 2 | 2 | - |
| none (identity) | 0.000 | 0.000 | 2.5e+03 | 1.0e+00 | -6.4e+06 | - | - | n.c. | n.c. | - |
| Jacobi diag(J)^-1 | 0.004 | 0.001 | 1.3e+01 | 1.0e+00 | -6.9e+03 | - | - | n.c. | n.c. | - |
| global trunc-eig J^-1, rank 108 | 0.478 | 0.478 | 9.9e-01 | 2.2e-04 | -2.2e+01 | 6.8e-02 | 6.8e-02 | indef. | n.c. | - |
| global large-eig J, rank 108 | 0.478 | 0.478 | 2.9e-01 | 1.0e+00 | -1.0e+00 | 1.5e-03 | 1.5e-03 | indef. | n.c. | - |
| global data-optimal, rank 108 | 0.477 | 0.478 | 8.3e-02 | 1.0e+00 | -8.5e-01 | -4.4e-04 | -4.4e-04 | indef. | n.c. | - |
| global trunc-eig J^-1, rank 396 | 1.752 | 1.752 | 4.3e-01 | 2.0e-06 | -2.0e+01 | 5.7e-02 | 5.7e-02 | indef. | 121 | - |
| global large-eig J, rank 396 | 1.752 | 1.752 | 1.3e-03 | 1.0e+00 | -3.4e-02 | 7.0e-05 | 7.0e-05 | indef. | 128 | - |
| global data-optimal, rank 396 | 1.748 | 1.752 | 6.0e-07 | 7.4e-01 | -2.9e-09 | 2.1e-09 | 2.1e-09 | indef. | 117 | - |
| global data-optimal, rank 29 | 0.128 | 0.128 | 5.3e-01 | 1.0e+00 | -1.4e+01 | -1.8e-03 | -1.8e-03 | indef. | n.c. | - |
| global data-optimal, rank 116 | 0.512 | 0.513 | 6.2e-02 | 1.0e+00 | -6.7e-01 | -1.3e-04 | -1.3e-04 | indef. | n.c. | - |
| global data-optimal, rank 41 | 0.181 | 0.181 | 4.2e-01 | 1.0e+00 | -1.4e+01 | -2.7e-03 | -2.7e-03 | indef. | n.c. | - |
| global data-optimal, rank 164 | 0.724 | 0.726 | 1.3e-02 | 1.0e+00 | -3.7e-01 | 3.5e-05 | 3.5e-05 | indef. | n.c. | - |
| Schwarz exact asm (R=4.0) | 1.163 | 1.157 | 1.1e+01 | 1.8e+00 | -1.3e+03 | -1.5e+01 | -1.3e+00 | 43 | 40 | - |
| Schwarz exact asm-pu (R=4.0) | 1.163 | 1.157 | 4.4e+03 | 3.6e-01 | -8.4e+06 | -1.7e+08 | -2.1e+02 | 206 | 89 | - |
| Schwarz exact asm-pu + coarse bal (R=4.0) | 1.221 | 1.226 | 9.2e+00 | 3.6e-01 | -3.5e+02 | -5.1e+02 | -9.6e-01 | 78 | 55 | - |
| Schwarz exact ras-core + coarse bal (R=4.0) | 1.221 | 0.359 | 2.9e+00 | 4.0e-01 | -1.4e+01 | -5.4e+00 | -7.2e-02 | indef. | 76 | - |
| Schwarz exact asm (R=6.0) | 2.479 | 2.473 | 2.8e+01 | 3.4e+00 | -1.7e+03 | -1.1e+02 | -2.8e+00 | 25 | 24 | - |
| Schwarz exact asm-pu (R=6.0) | 2.479 | 2.473 | 1.0e+05 | 2.5e-01 | -1.6e+08 | -1.0e+11 | -5.1e+03 | 223 | 87 | - |
| Schwarz exact asm-pu + coarse bal (R=6.0) | 2.537 | 2.541 | 2.5e+02 | 2.4e-01 | -1.0e+04 | -4.7e+05 | -2.4e+01 | 87 | 55 | - |
| Schwarz exact ras-core + coarse bal (R=6.0) | 2.537 | 0.552 | 2.3e+00 | 1.6e-01 | -7.5e+00 | -7.0e+00 | -7.5e-02 | indef. | 47 | - |
| local RRR r=8 + Galerkin (R=4.0) | 0.256 | 0.000 | 4.4e+03 | 9.7e-01 | -8.4e+06 | -1.7e+08 | -2.1e+02 | indef. | n.c. | - |
| learned Schwarz R=4.0 r_l=8 coarse=learned (3 seeds) | 0.259 | 0.134 | 1.1e-01 (1.0e-01..1.1e-01) | 1.0e+00 | -1.4e+00 | -4.5e-04 | -2.5e-04 | indef. | n.c. | 17 |
| learned Schwarz R=4.0 r_l=8 coarse=galerkin (3 seeds) | 0.256 | 0.134 | 9.4e-02 (9.4e-02..9.6e-02) | 1.0e+00 | -1.5e+00 | 5.1e-04 | 4.4e-04 | indef. | n.c. | 23 |
| learned Schwarz R=4.0 r_l=32 coarse=learned (3 seeds) | 1.029 | 0.520 | 7.5e-03 (7.5e-03..7.5e-03) | 1.0e+00 | -2.9e-01 | 6.3e-05 | 1.8e-04 | indef. | n.c. | 27 |
| learned Schwarz R=4.0 r_l=32 coarse=galerkin (3 seeds) | 1.026 | 0.520 | 7.0e-03 (7.0e-03..7.1e-03) | 1.0e+00 | -2.5e-01 | 9.6e-05 | 6.9e-05 | indef. | n.c. | 26 |
| learned Schwarz R=6.0 r_l=8 coarse=learned (3 seeds) | 0.365 | 0.187 | 1.0e-01 (1.0e-01..1.1e-01) | 1.0e+00 | -1.4e+00 | -4.3e-04 | -2.1e-04 | indef. | n.c. | 52 |
| learned Schwarz R=6.0 r_l=8 coarse=galerkin (3 seeds) | 0.362 | 0.187 | 9.2e-02 (9.2e-02..9.2e-02) | 1.0e+00 | -1.5e+00 | 4.0e-04 | 5.2e-04 | indef. | n.c. | 41 |
| learned Schwarz R=6.0 r_l=32 coarse=learned (3 seeds) | 1.452 | 0.732 | 7.7e-03 (7.7e-03..7.7e-03) | 1.0e+00 | -3.2e-01 | 5.3e-05 | 1.3e-04 | indef. | n.c. | 110 |
| learned Schwarz R=6.0 r_l=32 coarse=galerkin (3 seeds) | 1.450 | 0.732 | 7.1e-03 (7.1e-03..7.1e-03) | 1.0e+00 | -2.5e-01 | 5.4e-05 | 8.9e-05 | indef. | n.c. | 62 |
| learned global rank 108 (rank-matched, 3 seeds) | 0.954 | 0.478 | 8.4e-02 (8.4e-02..8.4e-02) | 1.0e+00 | -9.8e-01 | -3.7e-04 | -3.5e-04 | indef. | n.c. | 12 |
| learned global rank 29 (param-matched, 3 seeds) | 0.256 | 0.128 | 5.3e-01 (5.3e-01..5.3e-01) | 1.0e+00 | -1.4e+01 | -1.6e-03 | -1.4e-03 | indef. | n.c. | 7 |
| learned global rank 396 (rank-matched, 3 seeds) | 3.497 | 1.752 | 4.0e-03 (4.0e-03..4.4e-03) | 1.0e+00 | -8.4e-02 | 3.7e-05 | -3.8e-05 | indef. | n.c. | 70 |
| learned global rank 116 (param-matched, 3 seeds) | 1.024 | 0.513 | 6.2e-02 (6.2e-02..6.2e-02) | 1.0e+00 | -6.8e-01 | -9.1e-05 | -4.9e-05 | indef. | n.c. | 50 |
| learned global rank 41 (param-matched, 3 seeds) | 0.362 | 0.181 | 4.2e-01 (4.2e-01..4.2e-01) | 1.0e+00 | -1.4e+01 | -2.6e-03 | -2.5e-03 | indef. | n.c. | 35 |
| learned global rank 164 (param-matched, 3 seeds) | 1.448 | 0.726 | 1.4e-02 (1.4e-02..1.4e-02) | 1.0e+00 | -5.0e-01 | 1.2e-04 | 1.8e-05 | indef. | n.c. | 37 |

#### acene2: natm 18, naux 894, E_J 640.251 Eh, E_corr -1.293203 Eh, 3000 Adam steps, lr 0.001, Cholesky of J 11.6 ms; local density fit R=8 (outside the operator class): dE_J -5.6e-07 Eh

| operator | params / dense | flops / dense | val. rel. J-norm err | Frob. err (norm.) | dE_J robust (Eh/atom) | dE_MP2 robust (Eh/atom) | dE_MP2 plain (Eh/atom) | PCG it. | GMRES it. | train s |
|---|---|---|---|---|---|---|---|---|---|---|
| exact J^-1 (Cholesky) | 1.000 | 1.000 | 0.0e+00 | 0.0e+00 | -3.2e-14 | 7.2e-16 | 2.7e-11 | 2 | 2 | - |
| none (identity) | 0.000 | 0.000 | 1.5e+04 | 1.0e+00 | -7.3e+07 | - | - | n.c. | n.c. | - |
| Jacobi diag(J)^-1 | 0.002 | 0.001 | 3.2e+01 | 1.0e+00 | -5.6e+04 | - | - | n.c. | n.c. | - |
| global trunc-eig J^-1, rank 162 | 0.362 | 0.362 | 1.0e+00 | 1.5e-04 | -3.6e+01 | 7.2e-02 | 7.2e-02 | indef. | n.c. | - |
| global large-eig J, rank 162 | 0.362 | 0.362 | 2.9e-01 | 1.0e+00 | -1.4e+00 | 1.1e-03 | 1.1e-03 | indef. | n.c. | - |
| global data-optimal, rank 162 | 0.362 | 0.362 | 1.5e-01 | 1.0e+00 | -1.2e+01 | -3.0e-04 | -3.0e-04 | indef. | n.c. | - |
| global trunc-eig J^-1, rank 594 | 1.329 | 1.329 | 8.9e-01 | 4.7e-07 | -3.5e+01 | 7.2e-02 | 7.2e-02 | indef. | n.c. | - |
| global large-eig J, rank 594 | 1.329 | 1.329 | 9.5e-03 | 1.0e+00 | -1.1e-01 | 2.6e-04 | 2.6e-04 | indef. | n.c. | - |
| global data-optimal, rank 594 | 1.327 | 1.329 | 2.2e-04 | 9.8e-01 | -8.9e-05 | 7.5e-07 | 7.5e-07 | indef. | n.c. | - |
| global data-optimal, rank 29 | 0.065 | 0.065 | 6.9e-01 | 1.0e+00 | -2.6e+01 | 1.2e-02 | 1.2e-02 | indef. | n.c. | - |
| global data-optimal, rank 118 | 0.264 | 0.264 | 2.6e-01 | 1.0e+00 | -1.9e+01 | -1.0e-04 | -1.0e-04 | indef. | n.c. | - |
| global data-optimal, rank 75 | 0.168 | 0.168 | 4.2e-01 | 1.0e+00 | -2.2e+01 | 4.0e-04 | 4.0e-04 | indef. | n.c. | - |
| global data-optimal, rank 302 | 0.675 | 0.676 | 2.0e-02 | 1.0e+00 | -2.0e-01 | 3.2e-05 | 3.2e-05 | indef. | n.c. | - |
| Schwarz exact asm (R=4.0) | 0.906 | 0.903 | 3.3e+01 | 9.2e-01 | -8.1e+03 | -6.8e+02 | -6.6e+00 | 93 | 81 | - |
| Schwarz exact asm-pu (R=4.0) | 0.906 | 0.903 | 2.8e+04 | 9.7e-01 | -2.7e+08 | -6.7e+09 | -5.3e+03 | 515 | 130 | - |
| Schwarz exact asm-pu + coarse bal (R=4.0) | 0.949 | 0.951 | 2.1e+01 | 9.7e-01 | -1.2e+02 | -6.0e+02 | -2.2e+00 | 192 | 110 | - |
| Schwarz exact ras-core + coarse bal (R=4.0) | 0.949 | 0.304 | 3.0e+01 | 9.5e-01 | -3.6e+02 | -9.0e+02 | -1.9e+00 | indef. | 191 | - |
| Schwarz exact asm (R=6.0) | 5.516 | 5.512 | 1.3e+02 | 3.3e+00 | -9.3e+03 | -2.8e+03 | -1.5e+01 | 24 | 22 | - |
| Schwarz exact asm-pu (R=6.0) | 5.516 | 5.512 | 2.7e+06 | 6.1e-01 | -3.4e+09 | -6.0e+14 | -1.3e+06 | 975 | 137 | - |
| Schwarz exact asm-pu + coarse bal (R=6.0) | 5.559 | 5.560 | 2.5e+03 | 6.1e-01 | -2.6e+03 | -5.6e+07 | -1.8e+02 | 237 | 91 | - |
| Schwarz exact ras-core + coarse bal (R=6.0) | 5.559 | 0.659 | 1.2e+02 | 4.1e-01 | -9.9e+02 | -1.7e+04 | -5.7e+00 | indef. | 121 | - |
| local RRR r=8 + Galerkin (R=4.0) | 0.131 | 0.000 | 2.8e+04 | 1.0e+00 | -2.7e+08 | -6.7e+09 | -5.3e+03 | n.c. | n.c. | - |
| learned Schwarz R=4.0 r_l=8 coarse=learned (3 seeds) | 0.133 | 0.069 | 1.6e-01 (1.6e-01..1.6e-01) | 1.0e+00 | -2.0e+00 | -4.3e-04 | -4.2e-04 | indef. | n.c. | 50 |
| learned Schwarz R=4.0 r_l=8 coarse=galerkin (3 seeds) | 0.131 | 0.069 | 1.5e-01 (1.5e-01..1.5e-01) | 1.0e+00 | -1.9e+00 | 6.6e-04 | 3.8e-04 | indef. | n.c. | 53 |
| learned Schwarz R=4.0 r_l=32 coarse=learned (3 seeds) | 0.528 | 0.266 | 1.7e-02 (1.7e-02..1.8e-02) | 1.0e+00 | -4.0e-01 | 1.3e-04 | 8.3e-05 | indef. | n.c. | 224 |
| learned Schwarz R=4.0 r_l=32 coarse=galerkin (3 seeds) | 0.526 | 0.266 | 1.7e-02 (1.7e-02..1.7e-02) | 1.0e+00 | -3.1e-01 | 1.2e-04 | -8.9e-05 | indef. | n.c. | 116 |
| learned Schwarz R=6.0 r_l=8 coarse=learned (3 seeds) | 0.339 | 0.172 | 1.6e-01 (1.6e-01..1.6e-01) | 1.0e+00 | -1.8e+00 | -5.1e-04 | -4.3e-04 | indef. | n.c. | 69 |
| learned Schwarz R=6.0 r_l=8 coarse=galerkin (3 seeds) | 0.337 | 0.172 | 1.5e-01 (1.5e-01..1.5e-01) | 1.0e+00 | -1.8e+00 | 4.5e-04 | 4.5e-04 | indef. | n.c. | 91 |
| learned Schwarz R=6.0 r_l=32 coarse=learned (3 seeds) | 1.351 | 0.679 | 1.9e-02 (1.9e-02..1.9e-02) | 1.0e+00 | -2.8e-01 | 1.6e-04 | -2.0e-05 | indef. | n.c. | 80 |
| learned Schwarz R=6.0 r_l=32 coarse=galerkin (3 seeds) | 1.350 | 0.679 | 1.8e-02 (1.8e-02..1.8e-02) | 1.0e+00 | -2.7e-01 | 1.1e-04 | 3.7e-05 | indef. | n.c. | 86 |
| learned global rank 162 (rank-matched, 3 seeds) | 0.724 | 0.362 | 1.5e-01 (1.5e-01..1.5e-01) | 1.0e+00 | -8.1e+00 | -2.8e-04 | -2.1e-04 | indef. | n.c. | 50 |
| learned global rank 29 (param-matched, 3 seeds) | 0.130 | 0.065 | 6.9e-01 (6.9e-01..6.9e-01) | 1.0e+00 | -2.6e+01 | 1.1e-02 | 1.2e-02 | indef. | n.c. | 68 |
| learned global rank 594 (rank-matched, 3 seeds) | 2.655 | 1.329 | 1.4e-02 (1.4e-02..1.5e-02) | 1.0e+00 | -2.4e-01 | 1.1e-04 | 1.2e-04 | indef. | n.c. | 89 |
| learned global rank 118 (param-matched, 3 seeds) | 0.527 | 0.264 | 2.5e-01 (2.5e-01..2.6e-01) | 1.0e+00 | -1.9e+01 | -1.1e-04 | -5.9e-05 | indef. | n.c. | 29 |
| learned global rank 75 (param-matched, 3 seeds) | 0.335 | 0.168 | 4.2e-01 (4.2e-01..4.2e-01) | 1.0e+00 | -2.2e+01 | 4.4e-04 | 5.4e-04 | indef. | n.c. | 28 |
| learned global rank 302 (param-matched, 3 seeds) | 1.350 | 0.676 | 2.5e-02 (2.5e-02..2.5e-02) | 1.0e+00 | -3.3e-01 | 1.3e-04 | 1.3e-04 | indef. | n.c. | 40 |

Global truncated-eigen rows at the parameter-matched ranks and local RRR rows other than R=4, r_l=8 are omitted here (same pattern); `tables.py` prints all.

Reading:
- **At equal parameter count the learned Schwarz operator beats the learned global one** by 3-15x in loss when the
  budget is below the dense size (C8H18, R=4: 5.7e-2 vs 6.7e-1 at r_l=8, 1.3e-2 vs 1.8e-1 at r_l=32; same pattern for
  (H2O)4 and naphthalene). This is the paper's parameter-efficiency observation, reproduced. At 1.3-1.8x the dense size
  (R=6, r_l=32) the gap shrinks to 1-2x (C8H18: 1.3e-2 vs 1.2e-2).
- **The gain is capacity per parameter, not optimisability:** the learned param-matched global model reaches the
  closed-form optimum of its rank (C8H18 rank 35: 0.67 vs 0.66; rank 141: 0.18 vs 0.18). Adam has no difficulty there.
- **Where Adam does fall short is high accuracy:** at rank 858 (C8H18) the closed-form optimum is 4.6e-7, learned global
  1.1e-2; at rank 396 ((H2O)4) 6.0e-7 vs 4.0e-3. All learned models (Schwarz and global) stall at 4e-3..2e-2 after
  3000 steps; a link to cond(J) = 7e7-2.5e9 is plausible but not tested.
- Downstream, the best learned operators give 0.25-0.4 Eh/atom RI-J and 5e-5..1.6e-4 Eh/atom RI-MP2 error (robust);
  the target is 1e-6. The plain RI-MP2 expression is sometimes smaller than the robust one (error cancellation, not
  better fits).
- Classical Schwarz with exact local inverses is worse as a one-shot inverse (3..1e9 Eh/atom) but is a working
  **preconditioner** (plain additive, R=6: 24-27 PCG iterations; Jacobi and identity do not converge in 2000). The
  learned operators are rank deficient and indefinite and cannot precondition at all. For these molecule sizes the
  overlapping subdomains make the classical operator larger than the dense inverse (0.9-6.9x).
- The closed-form local reduced-rank blocks fail like `asm-pu` (overlap double counting of near-null directions);
  the learned blocks are therefore doing something non-trivial, they just do not reach useful accuracy.
- Reference outside the operator class: fitting each atomic density `rho_A` with the aux functions within 8 Bohr
  (`q4_metric.py` method) gives 1.6e-7..8.7e-7 Eh in total. Locality lives in the right-hand side partition, not in
  the operator.

## Result 3: transfer with parameters shared across equivalent subdomains

#### Transfer alkanes: train alkane6, alkane8; test alkane12, alkane16, alkane20; 3000 steps

- R, r_l = 4.0_8: 14 classes for 46 training subdomains, max |Jn_ii - Jn_jj| within a class 5.9e-15
- R, r_l = 4.0_32: 14 classes for 46 training subdomains, max |Jn_ii - Jn_jj| within a class 5.9e-15
- R, r_l = 6.0_8: 20 classes for 46 training subdomains, max |Jn_ii - Jn_jj| within a class 5.9e-15
- R, r_l = 6.0_32: 20 classes for 46 training subdomains, max |Jn_ii - Jn_jj| within a class 5.9e-15

| system | role | operator | params | val. rel. J-norm err | dE_J robust (Eh/atom) | dE_MP2 robust (Eh/atom) | GMRES it. | unseen subdomains |
|---|---|---|---|---|---|---|---|---|
| alkane6 | train | Jacobi diag(J)^-1 | 702 | 2.8e+01 | -1.7e+04 | - | n.c. | - |
| alkane6 | train | Schwarz exact asm (R=6.0) | 1995954 | 1.8e+02 | -5.8e+03 | -2.5e+03 | 19 | - |
| alkane6 | train | Schwarz exact asm-pu + coarse bal (R=6.0) | 2010906 | 7.6e+01 | -9.7e+02 | -2.4e+04 | 59 | - |
| alkane6 | train | Schwarz exact ras-core + coarse bal (R=6.0) | 2010906 | 8.6e+00 | -5.4e+00 | -1.0e+02 | 58 | - |
| alkane6 | train | global data-optimal, rank 660 | 463320 | 1.3e-07 | -5.2e-08 | 1.1e-10 | 91 | - |
| alkane6 | train | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 34368 | 5.3e-02 (5.3e-02..5.3e-02) | -9.7e-01 | 7.2e-04 | n.c. | 0 (dev 6e-15) |
| alkane6 | train | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 137472 | 1.1e-02 (1.1e-02..1.1e-02) | -2.0e-01 | 6.7e-05 | n.c. | 0 (dev 6e-15) |
| alkane6 | train | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 139008 | 4.7e-02 (4.7e-02..4.8e-02) | -8.9e-01 | 5.7e-04 | n.c. | 0 (dev 0e+00) |
| alkane6 | train | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 556032 | 1.1e-02 (1.1e-02..1.1e-02) | -2.2e-01 | 5.4e-05 | n.c. | 0 (dev 0e+00) |
| alkane8 | train | Jacobi diag(J)^-1 | 924 | 3.1e+01 | -2.5e+04 | - | n.c. | - |
| alkane8 | train | Schwarz exact asm (R=6.0) | 2921694 | 2.1e+02 | -1.0e+04 | -4.2e+03 | 26 | - |
| alkane8 | train | Schwarz exact asm-pu + coarse bal (R=6.0) | 2946993 | 5.4e+01 | -8.1e+02 | -1.8e+04 | 65 | - |
| alkane8 | train | Schwarz exact ras-core + coarse bal (R=6.0) | 2946993 | 1.4e+01 | -8.0e+00 | -3.9e+02 | 90 | - |
| alkane8 | train | global data-optimal, rank 858 | 792792 | 2.5e-07 | -1.7e-07 | 1.6e-10 | 141 | - |
| alkane8 | train | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 34368 | 5.6e-02 (5.6e-02..5.6e-02) | -1.1e+00 | 7.6e-04 | n.c. | 0 (dev 6e-15) |
| alkane8 | train | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 137472 | 1.2e-02 (1.2e-02..1.2e-02) | -2.8e-01 | 7.4e-05 | n.c. | 0 (dev 6e-15) |
| alkane8 | train | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 139008 | 5.2e-02 (5.1e-02..5.2e-02) | -9.5e-01 | 5.9e-04 | n.c. | 0 (dev 6e-15) |
| alkane8 | train | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 556032 | 1.2e-02 (1.2e-02..1.2e-02) | -2.3e-01 | 5.5e-05 | n.c. | 0 (dev 6e-15) |
| alkane12 | test | Jacobi diag(J)^-1 | 1368 | 3.5e+01 | -4.1e+04 | - | n.c. | - |
| alkane12 | test | Schwarz exact asm (R=6.0) | 4773174 | 2.5e+02 | -2.0e+04 | -6.4e+03 | 37 | - |
| alkane12 | test | Schwarz exact asm-pu + coarse bal (R=6.0) | 4827267 | 3.8e+01 | -5.5e+02 | -1.2e+04 | 75 | - |
| alkane12 | test | Schwarz exact ras-core + coarse bal (R=6.0) | 4827267 | 2.0e+01 | -1.4e+01 | -8.7e+02 | 133 | - |
| alkane12 | test | global data-optimal, rank 1254 | 1715472 | 4.9e-07 | -4.9e-07 | 4.6e-10 | 241 | - |
| alkane12 | test | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 34368 | 5.9e-02 (5.9e-02..6.0e-02) | -1.5e+00 | 8.0e-04 | n.c. | 0 (dev 6e-15) |
| alkane12 | test | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 137472 | 1.2e-02 (1.2e-02..1.2e-02) | -5.3e-01 | 8.3e-05 | n.c. | 0 (dev 6e-15) |
| alkane12 | test | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 139008 | 5.6e-02 (5.6e-02..5.6e-02) | -1.1e+00 | 6.3e-04 | n.c. | 0 (dev 6e-15) |
| alkane12 | test | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 556032 | 1.3e-02 (1.3e-02..1.3e-02) | -2.9e-01 | 6.3e-05 | n.c. | 0 (dev 6e-15) |
| alkane16 | test | Jacobi diag(J)^-1 | 1812 | 4.0e+01 | -5.5e+04 | - | n.c. | - |
| alkane16 | test | Schwarz exact asm (R=6.0) | 6624654 | 2.7e+02 | -3.1e+04 | -7.6e+03 | 44 | - |
| alkane16 | test | Schwarz exact asm-pu + coarse bal (R=6.0) | 6718341 | 2.9e+01 | -4.2e+02 | -9.5e+03 | 79 | - |
| alkane16 | test | Schwarz exact ras-core + coarse bal (R=6.0) | 6718341 | 2.4e+01 | -1.6e+01 | -1.2e+03 | 187 | - |
| alkane16 | test | global data-optimal, rank 1650 | 2989800 | 7.0e-07 | -9.5e-07 | 6.9e-10 | 341 | - |
| alkane16 | test | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 34368 | 6.2e-02 (6.1e-02..6.2e-02) | -2.0e+00 | 8.3e-04 | n.c. | 0 (dev 6e-15) |
| alkane16 | test | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 137472 | 1.3e-02 (1.3e-02..1.3e-02) | -8.7e-01 | 8.7e-05 | n.c. | 0 (dev 6e-15) |
| alkane16 | test | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 139008 | 6.0e-02 (5.9e-02..6.0e-02) | -1.2e+00 | 6.6e-04 | n.c. | 0 (dev 6e-15) |
| alkane16 | test | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 556032 | 1.3e-02 (1.3e-02..1.3e-02) | -3.6e-01 | 6.7e-05 | n.c. | 0 (dev 6e-15) |
| alkane20 | test | Jacobi diag(J)^-1 | 2256 | 4.5e+01 | -6.9e+04 | - | n.c. | - |
| alkane20 | test | Schwarz exact asm (R=6.0) | 8476134 | 2.8e+02 | -4.1e+04 | -8.4e+03 | 51 | - |
| alkane20 | test | Schwarz exact asm-pu + coarse bal (R=6.0) | 8620215 | 2.3e+01 | -3.4e+02 | -7.7e+03 | 88 | - |
| alkane20 | test | Schwarz exact ras-core + coarse bal (R=6.0) | 8620215 | 2.4e+01 | -1.8e+01 | -1.3e+03 | 239 | - |
| alkane20 | test | global data-optimal, rank 2046 | 4615776 | 8.2e-07 | -1.5e-06 | 7.0e-10 | n.c. | - |
| alkane20 | test | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 34368 | 6.1e-02 (6.1e-02..6.2e-02) | -2.5e+00 | 8.4e-04 | n.c. | 0 (dev 7e-15) |
| alkane20 | test | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 137472 | 1.3e-02 (1.3e-02..1.3e-02) | -1.3e+00 | 9.0e-05 | n.c. | 0 (dev 7e-15) |
| alkane20 | test | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 139008 | 6.1e-02 (6.0e-02..6.1e-02) | -1.3e+00 | 6.7e-04 | n.c. | 0 (dev 7e-15) |
| alkane20 | test | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 556032 | 1.4e-02 (1.4e-02..1.4e-02) | -4.4e-01 | 7.0e-05 | n.c. | 0 (dev 7e-15) |

#### Transfer acenes: train acene2, acene3; test acene4, acene5; 3000 steps

- R, r_l = 4.0_8: 14 classes for 42 training subdomains, max |Jn_ii - Jn_jj| within a class 6.4e-15
- R, r_l = 4.0_32: 14 classes for 42 training subdomains, max |Jn_ii - Jn_jj| within a class 6.4e-15
- R, r_l = 6.0_8: 22 classes for 42 training subdomains, max |Jn_ii - Jn_jj| within a class 6.4e-15
- R, r_l = 6.0_32: 22 classes for 42 training subdomains, max |Jn_ii - Jn_jj| within a class 6.4e-15

| system | role | operator | params | val. rel. J-norm err | dE_J robust (Eh/atom) | dE_MP2 robust (Eh/atom) | GMRES it. | unseen subdomains |
|---|---|---|---|---|---|---|---|---|
| acene2 | train | Jacobi diag(J)^-1 | 894 | 3.5e+01 | -5.6e+04 | - | n.c. | - |
| acene2 | train | Schwarz exact asm (R=6.0) | 2206842 | 1.3e+02 | -9.3e+03 | -2.8e+03 | 22 | - |
| acene2 | train | Schwarz exact asm-pu + coarse bal (R=6.0) | 2223999 | 2.1e+03 | -2.6e+03 | -5.6e+07 | 91 | - |
| acene2 | train | Schwarz exact ras-core + coarse bal (R=6.0) | 2223999 | 1.1e+02 | -9.9e+02 | -1.7e+04 | 123 | - |
| acene2 | train | global data-optimal, rank 594 | 531036 | 1.0e-04 | -3.7e-04 | 7.5e-07 | n.c. | - |
| acene2 | train | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 41856 | 1.4e-01 (1.4e-01..1.5e-01) | -2.0e+00 | 6.0e-04 | n.c. | 0 (dev 6e-15) |
| acene2 | train | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 167424 | 1.5e-02 (1.5e-02..1.5e-02) | -5.3e-01 | 1.5e-04 | n.c. | 0 (dev 6e-15) |
| acene2 | train | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 169632 | 1.4e-01 (1.4e-01..1.4e-01) | -1.8e+00 | 4.0e-04 | n.c. | 0 (dev 0e+00) |
| acene2 | train | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 678528 | 1.6e-02 (1.6e-02..1.6e-02) | -4.8e-01 | 1.5e-04 | n.c. | 0 (dev 0e+00) |
| acene3 | train | Jacobi diag(J)^-1 | 1230 | 4.1e+01 | -9.2e+04 | - | n.c. | - |
| acene3 | train | Schwarz exact asm (R=6.0) | 3438264 | 1.6e+02 | -1.9e+04 | -5.8e+03 | 30 | - |
| acene3 | train | Schwarz exact asm-pu + coarse bal (R=6.0) | 3469314 | 2.0e+03 | -2.0e+03 | -2.4e+07 | 116 | - |
| acene3 | train | Schwarz exact ras-core + coarse bal (R=6.0) | 3469314 | 1.8e+02 | -2.4e+03 | -4.9e+04 | 165 | - |
| acene3 | train | global data-optimal, rank 792 | 974160 | 1.6e-04 | -2.0e-03 | 1.0e-06 | n.c. | - |
| acene3 | train | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 41856 | 1.6e-01 (1.6e-01..1.6e-01) | -2.3e+00 | 5.3e-04 | n.c. | 0 (dev 6e-15) |
| acene3 | train | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 167424 | 1.6e-02 (1.6e-02..1.6e-02) | -6.7e-01 | 1.7e-04 | n.c. | 0 (dev 6e-15) |
| acene3 | train | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 169632 | 1.5e-01 (1.5e-01..1.5e-01) | -1.9e+00 | 3.2e-04 | n.c. | 0 (dev 6e-15) |
| acene3 | train | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 678528 | 1.8e-02 (1.8e-02..1.8e-02) | -5.3e-01 | 1.5e-04 | n.c. | 0 (dev 6e-15) |
| acene4 | test | Jacobi diag(J)^-1 | 1566 | 4.5e+01 | -1.3e+05 | - | n.c. | - |
| acene4 | test | Schwarz exact asm (R=6.0) | 4669686 | 1.7e+02 | -3.1e+04 | -8.7e+03 | 40 | - |
| acene4 | test | Schwarz exact asm-pu + coarse bal (R=6.0) | 4718697 | 1.3e+03 | -1.8e+03 | -1.9e+07 | 144 | - |
| acene4 | test | Schwarz exact ras-core + coarse bal (R=6.0) | 4718697 | 2.2e+02 | -3.6e+03 | -8.6e+04 | 243 | - |
| acene4 | test | global data-optimal, rank 990 | 1550340 | 2.1e-04 | -5.0e-03 | 1.2e-06 | n.c. | - |
| acene4 | test | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 41856 | 1.6e-01 (1.6e-01..1.7e-01) | -2.7e+00 | 5.0e-04 | n.c. | 0 (dev 6e-15) |
| acene4 | test | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 167424 | 1.8e-02 (1.8e-02..1.8e-02) | -8.4e-01 | 1.8e-04 | n.c. | 0 (dev 6e-15) |
| acene4 | test | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 169632 | 1.6e-01 (1.6e-01..1.6e-01) | -2.0e+00 | 3.0e-04 | n.c. | 0 (dev 6e-15) |
| acene4 | test | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 678528 | 2.1e-02 (2.1e-02..2.1e-02) | -6.9e-01 | 1.8e-04 | n.c. | 0 (dev 6e-15) |
| acene5 | test | Jacobi diag(J)^-1 | 1902 | 4.9e+01 | -1.7e+05 | - | n.c. | - |
| acene5 | test | Schwarz exact asm (R=6.0) | 5901108 | 1.9e+02 | -4.3e+04 | -1.1e+04 | 47 | - |
| acene5 | test | Schwarz exact asm-pu + coarse bal (R=6.0) | 5972148 | 1.1e+03 | -1.5e+03 | -1.6e+07 | 166 | - |
| acene5 | test | Schwarz exact ras-core + coarse bal (R=6.0) | 5972148 | 2.5e+02 | -4.5e+03 | -1.2e+05 | 309 | - |
| acene5 | test | global data-optimal, rank 1188 | 2259576 | 2.4e-04 | -9.0e-03 | 1.5e-06 | n.c. | - |
| acene5 | test | learned shared R=4.0 r_l=8 + Galerkin (3 seeds) | 41856 | 1.7e-01 (1.7e-01..1.7e-01) | -3.1e+00 | 4.7e-04 | n.c. | 0 (dev 6e-15) |
| acene5 | test | learned shared R=4.0 r_l=32 + Galerkin (3 seeds) | 167424 | 1.8e-02 (1.8e-02..1.8e-02) | -1.1e+00 | 1.9e-04 | n.c. | 0 (dev 6e-15) |
| acene5 | test | learned shared R=6.0 r_l=8 + Galerkin (3 seeds) | 169632 | 1.6e-01 (1.6e-01..1.7e-01) | -2.1e+00 | 3.0e-04 | n.c. | 0 (dev 6e-15) |
| acene5 | test | learned shared R=6.0 r_l=32 + Galerkin (3 seeds) | 678528 | 2.1e-02 (2.1e-02..2.1e-02) | -8.1e-01 | 2.1e-04 | n.c. | 0 (dev 6e-15) |

Reading:
- Sharing is exact here: the generated chains are perfectly periodic, equal keys give identical local metric blocks
  (deviation <= 7e-15), and every subdomain of the test molecules was seen in training (0 unseen).
- **Transfer works in the relative sense:** loss on unseen C12-C20 is the training loss +5..30 % (6.1e-2 vs 5.6e-2 at
  R=4, r_l=8; 1.4e-2 vs 1.2e-2 at R=6, r_l=32), tetracene/pentacene likewise +5..30 %; the parameter count is constant
  (34k-680k). This is what per-subdomain free parameters (the paper) cannot do.
- RI-MP2 error per atom stays constant with size (6-9e-5 Eh/atom at r_l=32, alkanes), RI-J error per atom grows
  (-0.2 to -1.3 Eh/atom); both far above 1e-6.
- The oracle (global data-optimal fitted on the test molecule itself at the Schwarz rank bound, 62-94 % of naux)
  reaches 1e-7..8e-7 loss for alkanes but only 1e-4..2.4e-4 for acenes (62-66 % rank).
- Classical Schwarz iterations grow with size (asm R=6: GMRES 19 -> 51 from C6 to C20). The plain additive variant has
  no coarse level; the balanced variants with the atomic-charge coarse space need more iterations (59 -> 88), so this
  coarse space does not remove the growth.

## Result 4: learning-rate sweep (the paper's Table 1 in a QC setting)

#### Learning-rate sweep, alkane4 (naux 480), 2000 steps, 3 seeds; Schwarz params 34760, ranks {'rank-matched': 126, 'param-matched': 36}

- closed-form optimum, global rank 126 (rank-matched): train 3.31e-02, val 4.96e-02
- closed-form optimum, global rank 36 (param-matched): train 4.19e-01, val 4.69e-01

| lr | model | params | final train loss median (min..max) | val loss median | Frob. err (norm.) | diverged |
|---|---|---|---|---|---|---|
| 1e-04 | schwarz | 34760 | 2.4e-01 (2.4e-01..2.4e-01) | 2.6e-01 | 1.0e+00 | 0/3 |
| 1e-04 | global rank-matched | 120960 | 8.0e-02 (7.9e-02..8.0e-02) | 1.0e-01 | 1.0e+00 | 0/3 |
| 1e-04 | global param-matched | 34560 | 4.3e-01 (4.3e-01..4.3e-01) | 4.8e-01 | 1.0e+00 | 0/3 |
| 3e-04 | schwarz | 34760 | 1.0e-01 (1.0e-01..1.0e-01) | 1.2e-01 | 1.0e+00 | 0/3 |
| 3e-04 | global rank-matched | 120960 | 4.2e-02 (4.2e-02..4.2e-02) | 6.0e-02 | 1.0e+00 | 0/3 |
| 3e-04 | global param-matched | 34560 | 4.2e-01 (4.2e-01..4.3e-01) | 4.8e-01 | 1.0e+00 | 0/3 |
| 1e-03 | schwarz | 34760 | 6.0e-02 (5.9e-02..6.0e-02) | 7.1e-02 | 1.0e+00 | 0/3 |
| 1e-03 | global rank-matched | 120960 | 3.7e-02 (3.7e-02..3.7e-02) | 5.4e-02 | 1.0e+00 | 0/3 |
| 1e-03 | global param-matched | 34560 | 4.2e-01 (4.2e-01..4.2e-01) | 4.7e-01 | 1.0e+00 | 0/3 |
| 3e-03 | schwarz | 34760 | 4.8e-02 (4.7e-02..4.8e-02) | 5.8e-02 | 1.0e+00 | 0/3 |
| 3e-03 | global rank-matched | 120960 | 3.6e-02 (3.5e-02..3.6e-02) | 5.2e-02 | 1.0e+00 | 0/3 |
| 3e-03 | global param-matched | 34560 | 4.2e-01 (4.2e-01..4.2e-01) | 4.7e-01 | 1.0e+00 | 0/3 |
| 1e-02 | schwarz | 34760 | 4.4e-02 (4.4e-02..4.5e-02) | 5.4e-02 | 1.0e+00 | 0/3 |
| 1e-02 | global rank-matched | 120960 | 3.9e-02 (3.9e-02..3.9e-02) | 5.6e-02 | 1.0e+00 | 0/3 |
| 1e-02 | global param-matched | 34560 | 4.2e-01 (4.2e-01..4.2e-01) | 4.7e-01 | 1.0e+00 | 0/3 |
| 3e-02 | schwarz | 34760 | 4.5e-02 (4.4e-02..4.5e-02) | 5.4e-02 | 1.0e+00 | 0/3 |
| 3e-02 | global rank-matched | 120960 | 5.7e-02 (5.6e-02..6.0e-02) | 7.8e-02 | 1.0e+00 | 0/3 |
| 3e-02 | global param-matched | 34560 | 4.2e-01 (4.2e-01..4.2e-01) | 4.7e-01 | 1.0e+00 | 0/3 |

Reading: in this QC setting the learned Schwarz operator is better than the rank-matched global one only at the
largest learning rate (3e-2: val 5.4e-2 vs 7.8e-2), worse at every lr <= 1e-2 (e.g. 1e-3: 7.1e-2 vs 5.4e-2), and nothing
diverges. The global model's 1100x collapse at large lr in the paper's Table 1 does not appear. The parameter-matched
global model sits at its closed-form optimum (train 0.42 = optimum 0.42) at every lr: no optimisation difficulty, only
too little capacity.

## Verdict per question

| question | result | verdict |
|---|---|---|
| Can a learned (or classical) structured G replace `J^-1` in one-shot fitting? | best learned: 0.25-0.4 Eh/atom RI-J, 5e-5..2e-4 Eh/atom RI-MP2 (robust); classical Schwarz 3..1e9 Eh/atom | **no**, 5-6 orders (RI-J) and 2 orders (RI-MP2) above 1e-6 Eh/atom |
| Is there a low-rank structure to exploit? | the best global operator needs 97-99.9 % of naux for 1e-6 Eh/atom RI-J (7 systems) | **no**, the fitting operator is full rank on relevant right-hand sides |
| Parameter efficiency of the Schwarz structure (paper's claim) | 3-15x lower loss than the global model at equal parameters below the dense size; 1-2x at 1.3-1.8x dense size | **reproduced**, but at loss 1e-2..1e-1, useless accuracy |
| Optimisability (paper's Table 1) | global model reaches its closed-form optimum at low rank; Schwarz better than rank-matched global only at lr 3e-2 (5.4e-2 vs 7.8e-2), worse at lr <= 1e-2; no divergence | **not reproduced** as a general advantage; at high accuracy every learned model stalls far above the optimum (1e-2 vs 5e-7) |
| Transfer / size extensivity with shared parameters | loss on unseen C12-C20 and tetracene/pentacene = training loss (+0..30 %), constant parameter count, 0 unseen subdomains | **works** for exactly periodic chains; absolute RI-J error per atom grows with size; water clusters: no equivalent subdomains, no transfer |
| Learned G as preconditioner | rank deficient, indefinite: PCG breaks down, GMRES does not converge | **no**; classical additive Schwarz (R=6) gives 19-51 GMRES/PCG iterations where Jacobi fails |
| Cost | applying the learned R=4, r_l=8 operator: 7-13 % of dense flops; training 17-224 s per molecule (numpy, shared machine) vs 0.7-12 ms for the Cholesky of J | training cost only pays with transfer; accuracy rules it out anyway |

**Recommendation (agent):** do not pursue learned Schwarz operators for the RI metric in curcuma. The metric is not a
bottleneck (Q4: global Cholesky 0.3-2 s up to naux ~10^4), the operator has no low-rank structure to exploit, and
the accuracy needed for RI (relative fit error ~1e-7) is 5 orders below what the learned structure reaches. If
Schwarz ideas are kept, the classical additive Schwarz operator as a **preconditioner** for an iterative metric solve
is the only variant here that works (19-51 iterations at R=6 Bohr on up to 2256 aux functions), and local fitting
of the right-hand side (q4 method, <= 1e-6 Eh total at R=8) is the working local alternative. A learned component
would have to be trained for preconditioning (spectrum of `GJ`), not for approximation; this was not tested.

## Not tested

- Secondary target (operators in the amplitude equations): not done.
- Larger training budgets (>3000 steps), other optimisers (L-BFGS, second order), symmetric parametrisation
  `Q_i Q_i^T`, a loss for preconditioning (e.g. condition number of `GJ`) instead of approximation.
- Other aux sets (def2-SVP/C, -JFIT), other elements than H/C/O, distorted geometries (the transfer test uses ideal,
  exactly periodic chains where shared classes are exactly equivalent), water transfer (impossible with this sharing).
- Coarse spaces beyond atomic charges (dipoles: code path exists, not run), coarse rank sweeps, `R_sd` > 6 Bohr.
- Compiled timings; all times are numpy on a shared machine.

## Literature (one short search, not exhaustive)

Learned components of Schwarz methods exist for PDEs: learned coarse spaces (arXiv:2607.06261; graph-based meshfree
coarse spaces arXiv:2603.26149), learned interface conditions for optimised Schwarz (Taghibakhshi et al., NeurIPS 2022,
arXiv:2205.09833). For RI the established route to locality is a different metric (attenuated Coulomb metric, Jung,
Sodt, Gill, Head-Gordon 2005) or local fit domains; ML work on density-fitting coefficients predicts the coefficients,
not the metric inverse (e.g. arXiv:2104.12457). No work learning an approximate inverse RI metric was found; this is
not a novelty claim.

## Reproduce

```
export LMP2_WORK=<dir> LMP2_THREADS=8; P=~/.venvs/qmref/bin/python
$P -I scripts/schwarz_learn/check_grad.py
$P -I scripts/schwarz_learn/rank_needed.py
$P -I scripts/schwarz_learn/run_free.py alkane8 3000 3     # also water4, acene2
$P -I scripts/schwarz_learn/run_transfer.py alkanes 3000 3 # also acenes
$P -I scripts/schwarz_learn/run_lrsweep.py alkane4 2000
$P -I scripts/schwarz_learn/tables.py
```
