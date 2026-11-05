# Seminumerical exchange: COSX / RIJCOSX (Nov 2026)

> 🤖 AI-generated, ⚙️ machine-tested. Not ✅ TESTED. COSX is an **approximation** of the exchange; it is opt-in and
> energy-only. Human production testing is pending.

Code: `src/core/energy_calculators/qm_methods/qm_cosx.{h,cpp}` (builder, engine glue), the potential integrals
`qmint::GridPotentialIntegrals` in `qm_integrals.cpp`. Tests: `ctest -L qm_cosx` (`test_cases/qm_cosx/`).
Reference script: `scripts/cosx_pyscf_reference.py`.

## Usage

```
curcuma -sp mol.xyz -method hf    -qm.exchange_method cosx                 # COSX: exact J, seminumerical K
curcuma -sp mol.xyz -method b3lyp -qm.exchange_method cosx -qm.ri j        # RIJCOSX: RI-J + seminumerical K
```

| Parameter | Default | Meaning |
|---|---|---|
| `-qm.exchange_method` | `exact` | `exact` (analytic K, or the RI-K of `-qm.ri jk`) or `cosx` |
| `-qm.cosx_grid` | `35x110` | radial x Lebedev points per atom, unpruned (`qmgrid::parseGridLevel`) |
| `-qm.cosx_fit_overlap` | `true` | overlap fitting (below) |
| `-qm.cosx_threshold` | `1e-11` | pair/primitive skip: bound x max\|F\| below it |

Radial scheme, partition and size adjustment are those of the DFT grid (`-qm.grid_radial`, `-qm.grid_partition`,
`-qm.grid_radii_adjust`); AO values per block are screened with `-qm.xc_ao_threshold`.

## Equations

Closed-shell exchange matrix (spin-summed density P, Fock matrix F = H + J - K/2):

    K_mn = sum_ls P_ls (ml|ns) = int phi_m(r) sum_ls phi_l(r) P_ls A_ns(r) d^3r,
    A_ns(r) = int phi_n(r') phi_s(r') / |r - r'| d^3r'.

COSX (Neese, Wennmohs, Hansen, Becker, Chem. Phys. 356, 98 (2009); idea Friesner 1985) integrates the first electron
on a grid {r_g, w_g} and the second analytically:

    F_gs = w_g sum_l phi_l(g) P'_ls          (points x n)
    G_gm = sum_s A_ms(g) F_gs                (A_ms(g): nuclear-attraction integral of a unit charge at r_g)
    K    = sym( sum_g phi(g) G(g)^T ),       sym(X) = (X + X^T)/2

A_ab(C) = sum_k (2 pi/p_k) sum_h E_k[ab][h] R_h(p_k, P_k - C) uses the Hermite tables E_k of the ERI kernel
(`makeMDPair`: contraction coefficients and the spherical transform folded in) and the Boys function; for
p_k |P_k - C|^2 > 40 + 4L the Boys function is replaced by its asymptote (relative error < 1e-14), and a one-centre
pair whose primitives are all in that range is evaluated once as derivatives of 1/r with the merged table
sum_k (pi/p_k)^(3/2) E_k. Per block of 128 points and per shell pair the primitive loop is a small GEMM
R (points x Hermite) * E_k^T.

**Overlap fitting** (Neese 2009, section 2.3; PySCF `sgx` `fit_ovlp`): P' = S~^-1 S P with the numerical overlap
S~ = sum_g w_g phi(g) phi(g)^T, so that the quadrature integrates every product phi_m phi_l against a constant like
the analytic S. The fit acts on one index only; K is symmetrised afterwards. Consequence: the COSX operator is not
the derivative of the energy expression -1/4 Tr(P K), so the SCF energy is not variational with respect to the grid
error (one reason there is no analytic gradient). Setup refuses a grid whose S~ is (near) singular in the basis
(smallest generalised eigenvalue of (S~, S) <= 1e-8).

**Screening**: a shell pair (A, B) is skipped at a point when b_AB max(|F_pA|, |F_pB|) < `cosx_threshold`, with the
heuristic pair bound b_AB = sum_k (2 pi/p_k) max_ab sum_h |E_k[ab][h]| (2 p_k)^(|h|/2) (not rigorous); the same test
per primitive pair with its share of b_AB. Measured effect at the default 1e-11: 1.5e-9 Eh (H2O/def2-SVP).

**Determinism**: blocks are summed in 16 fixed chunks, so 1 and 4 threads give bitwise the same energy
(`qm_cosx_threads_h2o`).

J comes from RI-J (`-qm.ri j`, RIJCOSX) or from the exact CPU path: the stored tensor, or integral-direct, where the
existing kernel digests K as well (no J-only direct kernel exists, so COSX with direct J saves nothing).

## Validation of the implementation

`scripts/cosx_pyscf_reference.py` runs PySCF 2.14 RHF/RKS with K replaced by the same COSX expression on curcuma's own
grid (dumped with `CURCUMA_QM_COSX_DUMP=<file>`), AO values from `mol.eval_gto` and the potential integrals from
libcint `int1e_grids` (an independent integral code):

| Case | Grid | curcuma - PySCF (Eh) |
|---|---|---|
| H2O HF/def2-SVP, fit | 35x110 | -1.4e-9 |
| NH3 B3LYP/def2-SVP, fit | 25x86 | +3.2e-10 |
| H2CO (distorted) HF/def2-SVP, no fit | 20x50 | +2.6e-10 |
| CH4 HF/def2-TZVP, fit | 25x86 | +1.4e-9 |
| H2O (distorted) HF/def2-SVP, fit and no fit, threshold 0 | 35x110 | <= 1e-10 |

The remaining 1e-9 is the default screen (`-qm.cosx_threshold 0` gives 1e-10). No comparison with ORCA's RIJCOSX was
made (different, pruned grids; ORCA not available in this session).

## Error against the exact exchange

E(COSX) - E(exact), micro-Eh, same SCF settings (`-qm.scf_threshold 1e-8`), J exact. Molecules: H2O, CH4, C2H6
(`test_cases/structures/small`), NH3 (`qm_1e`), H2CO distorted (`qm_grad`), benzene (`qm_hf3c`), C4 = n-butane,
idealised all-trans geometry (generated, not in the structure library).

HF / def2-SVP

| molecule | 20x50 | 25x86 | 35x110 | 50x194 | 35x110 no fit |
|---|---|---|---|---|---|
| H2O | +26.9 | +11.8 | +2.3 | -0.7 | +33.5 |
| NH3 | -174.1 | -14.9 | +7.0 | +1.1 | +84.3 |
| CH4 | +995.3 | -3.1 | -3.9 | -3.8 | -15.9 |
| H2CO | +10.7 | -128.2 | -32.5 | +1.9 | +14.5 |
| C2H6 | +697.4 | +361.6 | +3.7 | -23.5 | -1.9 |
| C4H10 | -1599.6 | -946.2 | -76.4 | +26.7 | -287.6 |
| C6H6 | +1812.6 | +85.7 | -39.2 | -5.4 | -60.5 |

HF / def2-TZVP

| molecule | 20x50 | 25x86 | 35x110 | 50x194 | 35x110 no fit |
|---|---|---|---|---|---|
| H2O | +339.9 | -47.2 | -7.1 | -0.4 | +16.7 |
| NH3 | -213.0 | +1.3 | -2.1 | +0.6 | +113.3 |
| CH4 | -28.7 | -2.0 | -1.9 | -0.7 | -28.9 |
| H2CO | +263.1 | -78.3 | +13.0 | +1.0 | -7.8 |
| C2H6 | -371.5 | -44.0 | -2.6 | -2.0 | +22.0 |
| C4H10 | -649.3 | +107.3 | +7.0 | -4.3 | -330.6 |
| C6H6 | -1619.9 | -107.7 | -26.2 | +13.4 | -177.7 |

B3LYP / def2-SVP (exact exchange weight 0.2, so about a fifth of the HF error)

| molecule | 20x50 | 25x86 | 35x110 | 50x194 | 35x110 no fit |
|---|---|---|---|---|---|
| H2O | +3.7 | +2.7 | +0.4 | -0.1 | +6.6 |
| NH3 | -37.4 | -3.3 | +1.0 | +0.2 | +16.2 |
| CH4 | +182.2 | -0.3 | -0.6 | -0.7 | -3.3 |
| H2CO | -3.3 | -24.9 | -6.2 | +0.4 | +2.0 |
| C2H6 | +130.5 | +69.6 | +0.6 | -4.6 | +0.9 |
| C4H10 | -312.8 | -182.8 | -16.4 | +5.5 | -58.8 |
| C6H6 | +383.2 | +28.8 | -1.5 | -1.4 | -16.3 |

B3LYP / def2-TZVP

| molecule | 20x50 | 25x86 | 35x110 | 50x194 | 35x110 no fit |
|---|---|---|---|---|---|
| H2O | +65.6 | -10.2 | -1.4 | -0.1 | +3.4 |
| NH3 | -42.7 | +0.6 | -0.1 | +0.2 | +22.1 |
| CH4 | -14.3 | -0.3 | -0.4 | -0.1 | -5.8 |
| H2CO | +42.2 | -16.5 | +2.9 | +0.2 | -2.3 |
| C2H6 | -75.4 | -6.4 | -0.6 | -0.4 | +5.6 |
| C4H10 | -121.3 | +17.0 | +1.7 | -0.8 | -73.9 |
| C6H6 | -312.9 | -20.4 | -5.7 | +2.1 | -34.8 |

Observations:

- At the default 35x110 the HF error is <= 76 micro-Eh (C4H10/SVP, about 5 micro-Eh per atom), B3LYP <= 16 micro-Eh;
  50x194 brings HF to <= 27 micro-Eh. The error has no sign and does not decrease monotonically with the grid
  (H2CO/SVP: 20x50 +10.7, 25x86 -128): it is a quadrature error of an unpruned grid, not a smooth function of size.
- Overlap fitting lowered the error at 35x110 in 23 of 28 cases; in 5 (H2CO in all four settings, C2H6 HF/SVP
  +3.7 vs -1.9 without) it did not. On 15x26 (H2CO HF/SVP, not in the tables) it took the error from 2.8e-2 to 4.1e-4 Eh.
- On a fine grid without fit the error goes to zero (H2O HF/SVP 75x302: -1.6e-7 Eh, ctest `qm_cosx_fine_h2o`).
- RIJCOSX: the COSX error with RI-J equals the one with exact J to 5e-10 Eh (H2O, `qm_cosx_rij_h2o`); the RI-J error
  (-1.07e-4 Eh for H2O/SVP) is larger than the COSX error there.
- The default 35x110 was chosen from these tables as the smallest tested grid whose HF error stayed below 1e-4 Eh on
  every molecule; it was not tuned per molecule. ORCA's COSX grids are pruned and smaller; not compared.

## Timings (4 vCPU)

HF, `-qm.threads 4` on the 4 vCPU of this machine (idle apart from the run), `-qm.scf_threshold 1e-6`, def2-universal-JFIT
for RI-J, def2-universal-JKFIT for RI-JK. Seconds per Fock build (J and K), in brackets the wall time of the whole run
(s). "stored (auto)" is `-qm.scf_direct auto -qm.eri_max_memory_mb 6000`: the stored tensor for benzene/SVP and
C4H10/SVP, integral-direct for all larger cases. "RI-J + exact K" = existing `-qm.ri j` (the direct digest forms K and
J). C_n = n-alkane C_nH_(2n+2), idealised all-trans geometries generated for this measurement (not in the structure
library). Raw logs were not kept in the repository.

| system | n_bf | stored (auto) | direct | RI-JK | RI-J + exact K | RIJCOSX 35x110 | RIJCOSX 25x86 | 35x110 / direct | 25x86 / direct |
|---|---|---|---|---|---|---|---|---|---|
| C6H6/SVP | 114 | 0.0551 (2) | 0.156 (2) | 0.0136 (0.5) | 0.145 (2) | 1.59 (18) | 0.844 (10) | 10.2 | 5.4 |
| C4/SVP | 106 | 0.0412 (2) | 0.107 (2) | 0.00708 (0.3) | 0.0941 (2) | 1.4 (19) | 0.865 (11) | 13.1 | 8.1 |
| C8/SVP | 202 | 0.65 (10) | 0.664 (10) | 0.103 (3) | 0.626 (10) | 5.8 (77) | 3.54 (47) | 8.7 | 5.3 |
| C12/SVP | 298 | 1.64 (27) | 1.58 (26) | 0.444 (11) | 1.69 (28) | 12.8 (183) | 7.66 (110) | 8.1 | 4.9 |
| C16/SVP | 394 | 2.92 (48) | 2.95 (48) | 1.28 (29) | 3.04 (52) | 22.5 (321) | 13.7 (198) | 7.6 | 4.7 |
| C20/SVP | 490 | 5.14 (84) | 4.92 (80) | 2.88 (64) | 5.07 (87) | 34.2 (489) | 20.0 (290) | 7.0 | 4.1 |
| C6H6/TZVP | 222 | 1.75 (24) | 1.73 (23) | 0.0446 (2) | 1.76 (24) | 4.85 (59) | 2.84 (38) | 2.8 | 1.6 |
| C4/TZVP | 184 | 0.926 (14) | 0.907 (13) | 0.0265 (1) | 0.904 (14) | 4.06 (53) | 2.40 (32) | 4.5 | 2.6 |
| C8/TZVP | 356 | 5.78 (91) | 5.89 (93) | 0.356 (10) | 5.57 (89) | 17.2 (244) | 10.7 (153) | 2.9 | 1.8 |
| C12/TZVP | 528 | - | 14.3 (227) | 1.35 (33) | - | 39.9 (568) | 23.8 (343) | 2.8 | 1.7 |
| C16/TZVP | 700 | - | 27.5 (441) | 3.92 (93) | - | 69.1 (987) | 39.9 (577) | 2.5 | 1.5 |

Error of the same runs against the exact energy (micro-Eh; `-qm.scf_threshold 1e-6`, so the last digit is not
significant): the COSX part (RIJCOSX - RI-J) grows with the chain, the RI-J part is larger at 35x110:

| system | RI-JK | RI-J + exact K | RIJCOSX 35x110 | RIJCOSX 25x86 | COSX part at 35x110 |
|---|---|---|---|---|---|
| C6H6/SVP | +116 | -303 | -343 | -218 | -39 |
| C4/SVP | +30 | -387 | -464 | -1333 | -76 |
| C8/SVP | +81 | -618 | -555 | -2810 | +63 |
| C12/SVP | +129 | -860 | -667 | -4280 | +194 |
| C16/SVP | +177 | -1103 | -779 | -5750 | +324 |
| C20/SVP | +226 | -1346 | -892 | -7220 | +454 |
| C6H6/TZVP | +130 | -282 | -309 | -390 | -26 |
| C4/TZVP | +28 | -284 | -277 | -176 | +7 |
| C8/TZVP | +78 | -454 | -538 | -259 | -84 |
| C12/TZVP | +128 | - | -809 | -347 | - |
| C16/TZVP | +178 | - | -1081 | -434 | - |

Findings:

- **No speedup in any measured case.** RIJCOSX is slower than the exact integral-direct J/K at every size: 4.1-13x
  per build with def2-SVP, 1.5-4.5x with def2-TZVP; RI-JK is faster than both everywhere (C16H34/TZVP 3.9 s vs 27.5 s
  direct vs 39.9 s RIJCOSX 25x86 per build) and more accurate than RIJCOSX at either grid.
- The ratio falls with system size and basis (TZVP, 25x86: 2.6, 1.8, 1.7, 1.5 from C4 to C16). A crossover against
  the exact direct path lies beyond C16H34/def2-TZVP (700 functions) if the trend continues; it was not reached, and
  it would not be one against RI-JK. At 25x86 the COSX error is already several mEh for C12-C20/SVP.
- Scaling from C12 to C16 (SVP): direct x1.87, RIJCOSX x1.75 for x1.33 atoms (about N^2.2 and N^1.9); the screen
  does not yet make COSX linear at this length: C12H26/SVP 25x86 evaluated 51 % of all (shell pair, point)
  combinations at `-qm.cosx_threshold 1e-11`, 41 % at 1e-9 (6.6 vs 7.8 s per build, energy moved 7e-10 Eh); at 1e-7
  the SCF did not converge in 100 iterations (the screen depends on the density and adds noise).
- Profile (benzene/SVP, 20x50, 1 thread, instrumented build not kept): 93 % of the primitive evaluations come from
  two-centre shell pairs; Boys function plus Hermite recursion about 60 % of the kernel time, the GEMM over the
  Hermite index about 30 %. The release build has no AVX (`USE_AVX2` off), so a point-vectorised Boys function would
  gain little here. Faster variants (pruned grids, a multipole far field for two-centre pairs, batched Boys values,
  ORCA's incremental screening) are not implemented.
- Published COSX speedups (Neese 2009) use pruned grids with fewer points per atom, further screening and larger
  bases with higher angular momentum, where four-centre integrals are more expensive; none of this was measured or
  reproduced here.

## Limits

- Not implemented: analytic gradient (requests are an error naming the option; `-opt`/`-md` fail, never a gradient of
  another energy), open shell (UHF/UKS), ECP atoms (refused, not validated), `-qm.ri jk` together with COSX,
  `-qm.ri_mixed_precision`, stability analysis, the analytic Hessian (falls back to the numerical Hessian, which
  then fails on the gradient), correlated methods and double hybrids (refused after the SCF), GPU (counted fallback
  `qm_cosx`, COSX and its J build run on the CPU), grid pruning, multipole far field for two-centre pairs, the
  incremental-density screening of ORCA, hf-3c (ignores the option with a warning), pure functionals (no exact
  exchange: option ignored with a warning).
- Not tested: elements beyond Ne (no ECP; Na-Kr all-electron untested), f/g in large systems (def2-TZVP f on C/N/O
  only), charged molecules, cartesian d (`-qm.cartesian_d`), grids other than the five in the tables, SCF
  convergence behaviour on difficult cases, systems above C20H42.
- The pair bound used for screening is a heuristic, not a rigorous bound.
